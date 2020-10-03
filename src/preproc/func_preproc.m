
function BFP_outfile = func_preproc(BFPPATH,t1,fmri,func_dir,TR,config)

% SCRIPT TO PREPROCESS THE FUNCTIONAL SCAN
% authors: ajoshi and choisoyo at USC
% for more information go to http://brainsuite.org/bfp/
%     this script has been modified from www.nitrc.org/projects/fcon_1000
%% set up parameters and directories
% default parameters
v = {'SimRef',1;...
    'FSLRigid',1;...
    'FWHM',6;...
    'HIGHPASS',0.005;...
    'LOWPASS',0.1;...    
    };
for i = 1:size(v,1)
if isfield(config, v{i,1})
%     config.(v{i,1})=str2double(config.(v{i,1}));
    if isnan(config.(v{i,1}))
        config.(v{i,1}) = v{i,2};
    end
else
    config.(v{i,1})=v{i,2};
end
end

FSLRigidReg = config.FSLRigid;
FWHM = config.FWHM;
hp = config.HIGHPASS;
lp = config.LOWPASS;
SimRef = config.SimRef;

nuisance_template=fullfile(BFPPATH,'supp_data','nuisance.fsf');
[~,n_vols]=unix(['3dinfo -nv ',fmri,'.nii.gz']); n_vols=str2double(n_vols);
TRstart=0;
TRend=n_vols-1;
sigma=str2double(FWHM)/2.3548;

[~,bname,~]=fileparts(fmri);example=[bname,'.example'];

disp('---------------------------------------');
disp('BFP fMRI PREPROCESSING !');
disp('---------------------------------------');

cwd = pwd;
cd(func_dir);
%% write out log file
ver_file = fullfile(BFPPATH, 'bfp_version.txt');

if ~exist(ver_file, 'file')
    ver_file = fullfile(BFPPATH, '/src/preproc', 'bfp_version.txt');
    if ~exist(ver_file, 'file')
        error('BFP directory: %s \n: Directory does not exist\n',BFPPATH);
    end
end
fid = fopen(ver_file, 'r');
[ver,~] = fscanf(fid, '%s',Inf);
fclose(fid);

logfname=[fmri,'.log.txt'];
fp=fopen(logfname,'a+');
t = now;
d = datetime(t,'ConvertFrom','datenum');
fprintf(fp,'\n%s\n',d);
fprintf(fp, 'func_preproc BFP version: %s\n', ver);
fprintf(fp,'fMRI: %s \nT1: %s \n', fmri,t1);
%% Deoblique
disp('Deobliquing');
if ~exist([fmri,'.dr.nii.gz'],'file')
    unix(['3dcalc -a ',fmri,'.nii.gz[',num2str(TRstart),'..',num2str(TRend),'] -expr ''a'' -prefix ',fmri,'.dr.nii.gz']);
    unix(['3drefit -deoblique ',fmri,'.dr.nii.gz']);
    fprintf(fp, '--Deoblique \n');
else
    disp('file found. skipping step')
    fprintf(fp, 'Deoblique file found. skipping step. \n');
end
%% Reorient into fsl friendly space (what AFNI calls RPI)
disp('Reorienting');
if ~exist([fmri,'.ro.nii.gz'],'file')
    unix(['3dresample -orient RPI -inset ',fmri,'.dr.nii.gz -prefix ',fmri,'.ro.nii.gz']);
    fprintf(fp, '--Reorient \n');
else
    disp('file found. skipping step')
    fprintf(fp, 'Reorient file found. skipping step. \n');
end
%% Get reference image used for motion correction and registration
disp('Getting image for motion correction')
outfile = [fmri,'.ro.mean.nii.gz'];
if ~exist(outfile,'file')
    if SimRef
        disp('Using SimRef...this may take a while if data is large...')
        orig = load_untouch_nii_gz([fmri,'.ro.nii.gz']);
        v_ref = fMRI_findRefv(orig,10,0);
        new = orig;
        new.img = squeeze(orig.img(:,:,:,v_ref));
        new.hdr.dime.dim(1)=3;
        new.hdr.dime.dim(5)=1;
        save_untouch_nii_gz(new,outfile);
        csvwrite([fmri,'.ssim.vref.txt'],v_ref);
        disp(['Reference volume computed using timepoint #',num2str(v_ref)])
        fprintf(fp,['--Option SimRef: Reference volume is timpoint #',num2str(v_ref),'\n']);
        clear new orig
    else
        unix(['3dTstat -mean -prefix ',outfile,' ',fmri,'.ro.nii.gz']);
        disp('Reference volume computed by averaging all volumes.')
        fprintf(fp,'--Reference volume computed by averaging all volumes.\n');
    end
else
    disp('file found. skipping step')
end
%% Motion correct to average of timeseries
disp('Motion correcting');
infile = outfile; clear outfile %[fmri,'.ro.mean.nii.gz'];
outfile = [fmri,'.mc.nii.gz'];
if ~exist(outfile,'file')
    unix(['3dvolreg -verbose -Fourier -twopass -base ',infile,' -zpad 4 -prefix ',fmri,'.mc.nii.gz -1Dfile ',fmri,'.mc.1D ',fmri,'.ro.nii.gz']);
    fprintf(fp, '--Motion Correction \n');
else
    disp('file found. skipping step')
    fprintf(fp, 'Motion Correction file found. skipping step. \n');
end

if ~exist([fmri,'.mco.txt'],'file')
    unix(['fsl_motion_outliers -i ',outfile,' -o ',fmri,'.mco -s ',fmri,'.mco.txt -p ',fmri,'.mco.png --dvars --nomoco -v']);
end
if ~exist([fmri,'.mc.ssim.txt'],'file')
    s(:,2) = motionEval(outfile, infile);
    s(:,1) = motionEval([fmri,'.ro.nii.gz'],infile);
    p = figure('visible','off');
    plot(s(:,1)); hold; plot(s(:,2));
    if exist([fmri,'.ssim.vref.txt'],'file')
        v_ref = importdata([fmri,'.ssim.vref.txt']);
        line([v_ref v_ref], [0 1],'Color','black','LineStyle','--');
        legend({'original','motion corrected','refence volume'},'Location','southeast');
    else
        legend({'original','motion corrected'},'Location','southeast');
    end
    ylim([0,1.1]);
    ylabel('SSIM');xlabel('vol no.');
    saveas(p,[fmri,'.mc.ssim.png']);
    csvwrite([fmri,'.mc.ssim.txt'],s);
end
    
clear infile
%% Get image for use in registration
disp('Getting image for coregistration...')
if ~exist([example,'.func.nii.gz'],'file')
    if SimRef
        disp('Using SimRef...')
        unix(['cp ',fmri,'.ro.mean.nii.gz ',example,'.func.nii.gz'])
    else
        disp('Using eigth image')
        unix(['3dcalc -a ',fmri,'.mc.nii.gz[7] -expr ''a'' -prefix ',example,'.func.nii.gz']);
    end
else
    disp('file found. skipping step')
end
%% FUNC->T1
disp('Performing registration to T1')
if ~exist([example,'.func2t1.nii.gz'],'file')
    if FSLRigidReg > 0
        disp('Using FSL rigid registration');
        unix(['flirt -ref ',t1,'.bfc.nii.gz -in ',example,'.func.nii.gz -out ',example,'.func2t1.nii.gz -omat ',example,'.func2t1.mat -cost corratio -dof 6 -interp trilinear']);
        unix(['convert_xfm -inverse -omat t12',example,'.func.mat ',example,'.func2t1.mat']);
        fprintf(fp, '--Registration to T1: Option FSL\n');
    else
        disp('Using USC rigid registration');
        moving_filename = [example,'.func.nii.gz'];
        static_filename = [t1,'.bfc.nii.gz'];
        output_filename = [example,'.func2t1.nii.gz'];
        usc_rigid_reg(moving_filename, static_filename, output_filename, 'inversion');
        fprintf(fp, '--Registration to T1: Option USC rigid registration\n');
    end
else
    disp('file found. skipping step')
    fprintf(fp, 'Registration to T1 files found. skipping step. \n');
end
%% Remove skull/edge detect
disp('Skull stripping');
outfile = [fmri,'.ss.nii.gz'];
if ~exist(outfile,'file')
    if str2double(config.T1mask) > 0
        if FSLRigidReg > 0
            unix(['flirt -ref ',example,'.func.nii.gz -in ',t1,'.mask.nii.gz -out ',fmri,'.mask.nii.gz -applyxfm -init t12',example,'.func.mat -interp nearestneighbour']);
            fprintf(fp, '--Option T1: mask fMRI\n');
        else
            transform_data_affine([t1,'.mask.nii.gz'], 's', [fmri,'.mask.temp.nii.gz'], [example,'.func.nii.gz'], [t1,'.bfc.nii.gz'], [fmri,'.example.func2t1.rigid.registration.result.mat'], 'nearest');
            msk = load_untouch_nii_gz([fmri,'.mask.temp.nii.gz']);
            msk.img(msk.img==255)=1;
            save_untouch_nii_gz(msk,[fmri,'.mask.nii.gz']);
            unix(['rm ',fmri,'.mask.temp.nii.gz'])
            fprintf(fp, '--Option T1: Skull Strip fMRI\n');
        end
    else
        unix(['3dAutomask -prefix ',fmri,'.mask.nii.gz -dilate 1 ',example,'.func.nii.gz']);
        fprintf(fp, '--Option autothreshold: Skull Strip fMRI\n');
    end
    unix(['3dcalc -a ',fmri,'.mc.nii.gz -b ',fmri,'.mask.nii.gz -expr ''a*b'' -prefix ',fmri,'.ss.nii.gz']);
else
    disp('file found. skipping step')
    fprintf(fp, 'Skull Strip fMRI file found. skipping step. \n');
end
%% skull strip image used for registration
disp( 'Skull stripping reference image')
orig = load_untouch_nii_gz([example,'.func.nii.gz']);
msk = load_untouch_nii_gz([fmri,'.mask.nii.gz']);
orig.img(msk.img==0)=0;
save_untouch_nii_gz(orig,[example,'.func.nii.gz']);
%% Spatial smoothing
disp('Spatial Smoothing');
infile = outfile; clear outfile % ss or mco 
outfile = [fmri,'.sm.nii.gz'];
if ~exist(outfile,'file')
    unix(['fslmaths ', infile,' -kernel gauss ',num2str(sigma),' -fmean -mas ',fmri,'.mask.nii.gz ',outfile]);
    fprintf(fp, '--Spatial Smoothing: sigma %s\n',num2str(sigma));
else
    disp('file found. skipping step')
    fprintf(fp, 'Spatial Smoothing file found. skipping step. \n');
end
%% Grandmean scaling
disp('Grand-mean scaling');
infile = outfile; clear outfile %[fmri,'.sm.nii.gz']
outfile = [fmri,'.gms.nii.gz'];
if ~exist(outfile,'file')
    unix(['fslmaths ',infile,' -ing 10000 ',outfile,' -odt float']);
    fprintf(fp, '--Grand Mean Scaling \n');
else
    disp('file found. skipping step')
    fprintf(fp, 'Grand Mean Scaling file found. skipping step. \n');
end
%% Temporal filtering
disp('Band-pass filtering');
gmsfile = [fmri,'.gms.nii.gz']; %%%%remove
infile = outfile; clear outfile %_gms
outfile = [fmri,'.filt.nii.gz'];
[~,n_vols]=unix(['3dinfo -nv ',infile]); n_vols=str2double(n_vols);
if ~exist(outfile,'file')
    %check if have min number of volumes for high pass filter
    minvol = 1/(str2double(TR)*str2double(hp));
    if minvol <= n_vols
        unix(['3dFourier -lowpass ',num2str(lp),' -highpass ',num2str(hp),' -retrend -prefix ',outfile,' ',infile]);
        fprintf(fp, '--Temporal Filter\n');
    else
        unix(['3dFourier -lowpass ',num2str(lp),' -retrend -prefix ',outfile,' ',infile]);
        fprintf(fp, '--Temporal Filter: too few volumes. Only lowpass filter was run.\n');
    end
else
    disp('file found. skipping step')
    fprintf(fp, 'Temporal Filter file found. skipping step.\n');
end
%% Detrending
if str2double(config.RunDetrend) > 0
    disp('Removing linear and quadratic trends');
    infile = outfile; clear outfile; %_filt
    outfile = [fmri,'.pp.nii.gz'];
    detrendfile = [fmri,'.pp.nii.gz']; %%%%remove
    if ~exist(outfile,'file')
        unix(['3dTstat -mean -prefix ',fmri,'.filt.mean.nii.gz ',infile]);
        unix(['3dDetrend -polort 2 -prefix ',fmri,'.dt.nii.gz ',infile]);
        unix(['3dcalc -a ',fmri,'.filt.mean.nii.gz -b ',fmri,'.dt.nii.gz -expr ''a+b'' -prefix ',outfile]);
        fprintf(fp, '--Option Detrend: Yes\n');
    else
        disp('file found. skipping step')
        fprintf(fp, 'Detrend file found. skipping step.\n');
    end
else
    fprintf(fp, '--Option Detrend: No\n');
end
%% FUNC->standard (3mm)
disp('Performing registration to standard space')
if ~exist([example,'.func2standard.nii.gz'],'file')
    if FSLRigidReg > 0
        unix(['flirt -ref standard.nii.gz -in ',example,'.func.nii.gz -out ',example,'.func2standard.nii.gz -omat ',example,'.func2standard.mat -cost corratio -dof 6 -interp trilinear']);
        % # Create mat file for conversion from subject's anatomical to functional
        unix(['convert_xfm -inverse -omat standard2',example,'.func.mat ',example,'.func2standard.mat']);
        fprintf(fp, '--Registration to Standard: Option FSL\n');
    else
        disp('Using USC rigid registration');
        moving_filename = [example,'.func.nii.gz'];%fullfile(funcDir,sprintf('%s_%s_bold_example_func.nii.gz',subid,sessionid{ind}));
        static_filename = 'standard.nii.gz';
        output_filename = [example,'.func2standard.nii.gz'];
        usc_rigid_reg(moving_filename, static_filename, output_filename, 'inversion');
        fprintf(fp, '--Registration to Standard: Option USC rigid registration\n');
    end
else
    disp('file found. skipping step')
    fprintf(fp, 'Registration to Standard file found. skipping step\n');
end
%% Nuisance Signal Regression
if str2double(config.RunNSR) > 0
    % ## 13.
    [~,fmribase,~] = fileparts(fmri);
    nuisance_dir=[func_dir,'/nuisance_',fmribase];
    %
    disp(' --------------------------------------------');
    disp(' !!!! RUNNING NUISANCE SIGNAL REGRESSION !!!!');
    disp(' --------------------------------------------');
    %
    %
    % ## 14. make nuisance directory
    if exist(nuisance_dir,'dir')
        rmdir(nuisance_dir,'s')
    end
    mkdir(nuisance_dir)
    %
    % # 15. Seperate motion parameters into seperate files
    disp('Splitting up subject motion parameters');
    fmc = [fmri,'.mc.1D'];
    unix(['awk ''{print $1}'' ',fmc,' > ',nuisance_dir,'/mc1.1D']);
    unix(['awk ''{print $2}'' ',fmc,' > ',nuisance_dir,'/mc2.1D']);
    unix(['awk ''{print $3}'' ',fmc,' > ',nuisance_dir,'/mc3.1D']);
    unix(['awk ''{print $4}'' ',fmc,' > ',nuisance_dir,'/mc4.1D']);
    unix(['awk ''{print $5}'' ',fmc,' > ',nuisance_dir,'/mc5.1D']);
    unix(['awk ''{print $6}'' ',fmc,' > ',nuisance_dir,'/mc6.1D']);
    
    %
    % # Extract signal for global, csf, and wm
    % ## 16. Global
    disp('Extracting global signal for subject');
    unix(['3dmaskave -mask ',fmri,'.mask.nii.gz -quiet ',outfile,' > ',nuisance_dir,'/global.1D'])
    %
    % ## 17. csf
    if FSLRigidReg > 0
        unix(['flirt -ref ',example,'.func.nii.gz -in ',t1,'.pvc.label.nii.gz -out ',fmri,'.pvc.label.nii.gz -applyxfm -init t12',example,'.func.mat -interp nearestneighbour']);
    else
        transform_data_affine([t1,'.pvc.label.nii.gz'], 's', [fmri,'.pvc.label.nii.gz'], [example,'.func.nii.gz'], [t1,'.bfc.nii.gz'], [fmri,'.example.func2t1.rigid.registration.result.mat'], 'nearest');
    end
    %
    pvclbl = load_untouch_nii_gz([fmri,'.pvc.label.nii.gz']);
    csfmsk = pvclbl;
    csfmsk.img = zeros(size(pvclbl.img));
    csfmsk.img(pvclbl.img==1 | pvclbl.img==6)=1;
    save_untouch_nii_gz(csfmsk,[fmri,'.csf.mask.nii.gz'],2);
%     unix(['fslmaths ',fmri,'.pvc.label.nii.gz -thr 5.5 -bin ',fmri,'.csf.mask.nii.gz']);
    unix(['fslmaths ',fmri,'.pvc.label.nii.gz -thr 2.5 -uthr 3.5 -bin ',fmri,'.wm.mask.nii.gz -odt char']);
    %
    % echo "Extracting signal from csf"
    unix(['3dmaskave -mask ',fmri,'.csf.mask.nii.gz -quiet ',outfile,' > ',nuisance_dir,'/csf.1D'])
    %
    % ## 18. wm
    disp('Extracting signal from white matter for subject');
    unix(['3dmaskave -mask ',fmri,'.wm.mask.nii.gz -quiet ',outfile,' > ',nuisance_dir,'/wm.1D'])
    %
    % ## 6. Generate mat file (for use later)
    % ## create fsf file
    %
    [~,n_vols]=unix(['3dinfo -nv ',fmri,'.pp.nii.gz']); n_vols=str2double(n_vols);
    disp('Modifying model file');
    unix(['sed -e s:nuisance_dir:"',nuisance_dir,'":g <',nuisance_template,' >',nuisance_dir,'/temp1']);
    unix(['sed -e s:nuisance_model_outputdir:"',nuisance_dir,'/residuals.feat":g <',nuisance_dir,'/temp1 >',nuisance_dir,'/temp2']);
    unix(['sed -e s:nuisance_model_TR:"',num2str(TR),'":g <',nuisance_dir,'/temp2 >',nuisance_dir,'/temp3']);
    unix(['sed -e s:nuisance_model_numTRs:"',num2str(n_vols),'":g <',nuisance_dir,'/temp3 >',nuisance_dir,'/temp4']);
    unix(['sed -e s:nuisance_model_input_data:"',func_dir,'/',outfile,'":g <',nuisance_dir,'/temp4 >',nuisance_dir,'/nuisance.fsf']);
    %
    % #rm ${nuisance_dir}/temp*
    %
    disp('Running feat model');
    unix(['feat_model ',nuisance_dir,'/nuisance']);
    %
    
    pp = load_untouch_nii_gz(outfile);
    msk = load_untouch_nii_gz([fmri,'.mask.nii.gz']);
    ppimg = pp.img(logical(msk.img));
    minVal = num2str(min(ppimg));
    
    %[~,minVal]=unix(['3dBrickStat -min -mask ',fmri,'.mask.nii.gz ',fmri,'.pp.nii.gz']);
    %minVal=str2double(minVal);
    %
    % ## 7. Get residuals
    disp('Running film to get residuals');
    unix(['film_gls --rn=',nuisance_dir,'/stats --noest --sa --ms=5 --in=',outfile,' --pd=',nuisance_dir,'/nuisance.mat --thr=',minVal]);
    %
    % ## 8. Demeaning residuals and ADDING 100
    unix(['3dTstat -mean -prefix ',nuisance_dir,'/stats/res4d.mean.nii.gz ',nuisance_dir,'/stats/res4d.nii.gz']);
    if exist([fmri,'.res.nii.gz'],'file')
        delete([fmri,'.res.nii.gz'])
    end
    unix(['3dcalc -a ',nuisance_dir,'/stats/res4d.nii.gz -b ',nuisance_dir,'/stats/res4d.mean.nii.gz -expr ''(a-b)+100'' -prefix ',fmri,'.res.nii.gz']);
    
    RS_infile = [fmri,'.res.nii.gz'];
    BFP_outfile = [fmri,'.res2standard.nii.gz'];
    fprintf(fp, '--Option Nuissance Signal Regression: Yes\n');
else
    if str2double(config.RunDetrend) > 0
        RS_infile = detrendfile;
        BFP_outfile = [fmri,'.pp2standard.nii.gz'];
        fprintf(fp, '--Option Nuissance Signal Regression: No\n');
    else
        RS_infile = gmsfile;
        BFP_outfile = [fmri,'.gms2standard.nii.gz'];
        fprintf(fp, '--Option Nuissance Signal Regression: No\n');
    end
end
%% Resampling residuals to standard MNI space
if FSLRigidReg > 0
    unix(['flirt -ref ',func_dir,'/standard -in ',RS_infile,' -out ',BFP_outfile,' -applyxfm -init ',func_dir,'/',example,'.func2standard.mat -interp trilinear']);
else
        transform_data_affine(RS_infile, 'm', BFP_outfile, [example,'.func.nii.gz'], 'standard.nii.gz', [fmri,'.example.func2standard.rigid.registration.result.mat'], 'linear');
end
fprintf(fp, '--Transform fMRI data to standard space\n');
fclose(fp);
cd(cwd);


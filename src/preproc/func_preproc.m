
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

usc_rigid_reg_bin = fullfile(BFPPATH, 'usc_rigid_reg.sh');
transform_data_affine_bin = fullfile(BFPPATH, 'transform_data_affine.sh');

[pth,bname,ext]=fileparts(fmri);example=[bname,'_example'];

disp('---------------------------------------');
disp('BFP fMRI PREPROCESSING !');
disp('---------------------------------------');

cwd = pwd;
cd(func_dir);
%% Dropping first # TRS & Deoblique
disp('Dropping first TRs')
disp('Deobliquing');
if ~exist([fmri,'_dr.nii.gz'],'file')
    unix(['3dcalc -a ',fmri,'.nii.gz[',num2str(TRstart),'..',num2str(TRend),'] -expr ''a'' -prefix ',fmri,'_dr.nii.gz']);
    unix(['3drefit -deoblique ',fmri,'_dr.nii.gz']);
else
    disp('file found. skipping step')
end
%% Reorient into fsl friendly space (what AFNI calls RPI)
disp('Reorienting');
if ~exist([fmri,'_ro.nii.gz'],'file')
    unix(['3dresample -orient RPI -inset ',fmri,'_dr.nii.gz -prefix ',fmri,'_ro.nii.gz']);
else
    disp('file found. skipping step')
end
%% Get reference image used for motion correction and registration
disp('Getting image for motion correction')
if ~exist([fmri,'_ro_mean.nii.gz'],'file')
    if SimRef
        disp('Using SimRef...')
        orig = load_untouch_nii_gz([fmri,'_ro.nii.gz']);
        [x,y,z,nvol] = size(orig.img);
        k = round(nvol/2);
        s= zeros(nvol,2);
        disp('Calculating similarity measures...')
        for i = 1:nvol
            data1 = orig.img(round(x/4):round(3*x/4),round(y/4):round(3*y/4),round(z/4):round(3*z/4),k);
            data2 = orig.img(round(x/4):round(3*x/4),round(y/4):round(3*y/4),round(z/4):round(3*z/4),i);
            s(i,2) = ssim(data1,data2);
        end
        %
        data = orig.img(:,:,:,s(:,2)>0.9);
        mdata = mean(data,4);
        new = orig;
        new.img = mdata;
        new.hdr.dime.dim(1)=3;
        new.hdr.dime.dim(5)=1;
        save_untouch_nii_gz(new,[fmri,'_ro_mean.nii.gz']);
        disp(['Reference volume computed using ',num2str(sum(s(:,1)>0.9)),' volumes'])
        clear new orig data mdata data1 data2 x y z nvol
    else
        unix(['3dTstat -mean -prefix ',fmri,'_ro_mean.nii.gz ',fmri,'_ro.nii.gz']);
        disp('Reference volume computed by averaging all volumes.')
    end
else
    disp('file found. skipping step')
end
%% Motion correct to average of timeseries
disp('Motion correcting');
if ~exist([fmri,'_mc.nii.gz'],'file')
    unix(['3dvolreg -verbose -Fourier -twopass -base ',fmri,'_ro_mean.nii.gz -zpad 4 -prefix ',fmri,'_mc.nii.gz -1Dfile ',fmri,'_mc.1D ',fmri,'_ro.nii.gz']);
    orig = load_untouch_nii_gz([fmri,'_mc.nii.gz']);
    [x,y,z,nvol] = size(orig.img);
    k = round(nvol/2);
    disp('Calculating similarity measures to evaluate motion correction...')
    for i = 1:nvol
        data1 = orig.img(round(x/4):round(3*x/4),round(y/4):round(3*y/4),round(z/4):round(3*z/4),k);
        data2 = orig.img(round(x/4):round(3*x/4),round(y/4):round(3*y/4),round(z/4):round(3*z/4),i);
        s(i,1) = ssim(data1,data2);
    end
%     p = figure;
%     if size(s,2)==2
%         plot(s(:,1)); hold; plot(s(:,2));
%         line([k k], [0 1],'Color','black','LineStyle','--'); 
%         legend({'motion corrected','original','refence volume'},'Location','southeast');
%     else
%         plot(s(:,1));
%         line([k k], [0 1],'Color','red','LineStyle','--'); 
%         legend({'motion corrected','refence volume'},'Location','southeast');
%     end
%     ylim([0,1.1]);
%     ylabel('SSIM');xlabel('vol no.');
%     saveas(p,[fmri,'_mc_ssim.png']);
%     close all
else
    disp('file found. skipping step')
end
%% Get image for use in registration
disp('Getting image for coregistration...')
if ~exist([example,'_func.nii.gz'],'file')
    if SimRef
        disp('Using SimRef...')
        unix(['cp ',fmri,'_ro_mean.nii.gz ',example,'_func.nii.gz'])
    else
        disp('Using eigth image')
        unix(['3dcalc -a ',fmri,'_mc.nii.gz[7] -expr ''a'' -prefix ',example,'_func.nii.gz']);
    end
else
    disp('file found. skipping step')
end
%% FUNC->T1
disp('Performing registration to T1')
if ~exist([example,'_func2t1.nii.gz'],'file')
    if FSLRigidReg > 0
        disp('Using FSL rigid registration');
        unix(['flirt -ref ',t1,'.bfc.nii.gz -in ',example,'_func.nii.gz -out ',example,'_func2t1.nii.gz -omat ',example,'_func2t1.mat -cost corratio -dof 6 -interp trilinear']);
        %     # Create mat file for conversion from subject's anatomical to functional
        unix(['convert_xfm -inverse -omat t12',example,'_func.mat ',example,'_func2t1.mat']);
    else
        disp('Using USC rigid registration');
        moving_filename = [example,'_func.nii.gz'];
        static_filename = [t1,'.bfc.nii.gz'];
        output_filename = [example,'_func2t1.nii.gz'];
%         if isdeployed
%             cmd = sprintf('%s %s %s %s %s %s', usc_rigid_reg_bin, moving_filename, static_filename, output_filename, 'inversion');
%             unix(cmd);
%         else
            usc_rigid_reg(moving_filename, static_filename, output_filename, 'inversion');
%         end
    end
else
    disp('file found. skipping step')
end
%% Remove skull/edge detect
disp('Skull stripping');
if ~exist([fmri,'_ss.nii.gz'],'file')
    if str2double(config.T1mask) > 0
        if FSLRigidReg > 0
            unix(['flirt -ref ',example,'_func.nii.gz -in ',t1,'.mask.nii.gz -out ',fmri,'_mask.nii.gz -applyxfm -init t12',example,'_func.mat -interp nearestneighbour']);
        else
%             if isdeployed
%                 cmd = sprintf('%s %s %s %s %s %s %s %s', transform_data_affine_bin, [t1,'.mask.nii.gz'], 's', [fmri,'_mask.nii.gz'], [example,'_func.nii.gz'], [t1,'.bfc.nii.gz'], [fmri,'_example_func2t1.rigid_registration_result.mat'], 'nearest');
%                 unix(cmd);
%             else
                transform_data_affine([t1,'.mask.nii.gz'], 's', [fmri,'_mask.temp.nii.gz'], [example,'_func.nii.gz'], [t1,'.bfc.nii.gz'], [fmri,'_example_func2t1.rigid_registration_result.mat'], 'nearest');
%                 unix(['3dcalc -a ',fmri,'_mask.temp.nii.gz -expr ''a/255'' -prefix ',fmri,'_mask.nii.gz']);
                msk = load_untouch_nii_gz([fmri,'_mask.temp.nii.gz']);
                msk.img(msk.img==255)=1;
                save_untouch_nii_gz(msk,[fmri,'_mask.nii.gz']);
                unix(['rm ',fmri,'_mask.temp.nii.gz'])
%             end
        end
    else
        unix(['3dAutomask -prefix ',fmri,'_mask.nii.gz -dilate 1 ',example,'_func.nii.gz']);
    end
    unix(['3dcalc -a ',fmri,'_mc.nii.gz -b ',fmri,'_mask.nii.gz -expr ''a*b'' -prefix ',fmri,'_ss.nii.gz']);
else
    disp('file found. skipping step')
end
%% skull strip image used for registration
disp( 'Skull stripping reference image')
orig = load_untouch_nii_gz([example,'_func.nii.gz']);
msk = load_untouch_nii_gz([fmri,'_mask.nii.gz']);
orig.img(msk.img==0)=0;
save_untouch_nii_gz(orig,[example,'_func.nii.gz']);
%% Spatial smoothing
disp('Spatial Smoothing');
if ~exist([fmri,'_sm.nii.gz'],'file')
    unix(['fslmaths ',fmri,'_ss.nii.gz -kernel gauss ',num2str(sigma),' -fmean -mas ',fmri,'_mask.nii.gz ',fmri,'_sm.nii.gz']);
else
    disp('file found. skipping step')
end
%% Grandmean scaling
disp('Grand-mean scaling');
if ~exist([fmri,'_gms.nii.gz'],'file')
    unix(['fslmaths ',fmri,'_sm.nii.gz -ing 10000 ',fmri,'_gms.nii.gz -odt float']);
else
    disp('file found. skipping step')
end
%% Temporal filtering
disp('Band-pass filtering');
gmsfile = [fmri,'_gms.nii.gz'];
if ~exist([fmri,'_filt.nii.gz'],'file')
    unix(['3dFourier -lowpass ',num2str(lp),' -highpass ',num2str(hp),' -retrend -prefix ',fmri,'_filt.nii.gz ',gmsfile]);
else
    disp('file found. skipping step')
end
%% Detrending
if str2double(config.RunDetrend) > 0
    disp('Removing linear and quadratic trends');
    detrendfile = [fmri,'_pp.nii.gz'];
    maskstep_infile = detrendfile;
    if ~exist(detrendfile,'file')
        unix(['3dTstat -mean -prefix ',fmri,'_filt_mean.nii.gz ',fmri,'_filt.nii.gz']);
        unix(['3dDetrend -polort 2 -prefix ',fmri,'_dt.nii.gz ',fmri,'_filt.nii.gz']);
        unix(['3dcalc -a ',fmri,'_filt_mean.nii.gz -b ',fmri,'_dt.nii.gz -expr ''a+b'' -prefix ',detrendfile]);
    else
        disp('file found. skipping step')
    end
else
    maskstep_infile = gmsfile;
end
%% FUNC->standard (3mm)
disp('Performing registration to standard space')
if ~exist([example,'_func2standard.nii.gz'],'file')
    if FSLRigidReg > 0
        unix(['flirt -ref standard.nii.gz -in ',example,'_func.nii.gz -out ',example,'_func2standard.nii.gz -omat ',example,'_func2standard.mat -cost corratio -dof 6 -interp trilinear']);
        % # Create mat file for conversion from subject's anatomical to functional
        unix(['convert_xfm -inverse -omat standard2',example,'_func.mat ',example,'_func2standard.mat']);
    else
        disp('Using USC rigid registration');
        moving_filename = [example,'_func.nii.gz'];%fullfile(funcDir,sprintf('%s_%s_bold_example_func.nii.gz',subid,sessionid{ind}));
        static_filename = 'standard.nii.gz';
        output_filename = [example,'_func2standard.nii.gz'];
%         if isdeployed
%             cmd = sprintf('%s %s %s %s %s %s', usc_rigid_reg_bin, moving_filename, static_filename, output_filename, 'inversion');
%             unix(cmd);
%         else
            usc_rigid_reg(moving_filename, static_filename, output_filename, 'inversion');
%         end
    end
else
    disp('file found. skipping step')
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
    %unix(['mkdir -p ',nuisance_dir]);
    if exist(nuisance_dir,'dir')
        rmdir(nuisance_dir,'s')
    end
    mkdir(nuisance_dir)
    %
    % # 15. Seperate motion parameters into seperate files
    disp('Splitting up subject motion parameters');
    unix(['awk ''{print $1}'' ',fmri,'_mc.1D > ',nuisance_dir,'/mc1.1D']);
    unix(['awk ''{print $2}'' ',fmri,'_mc.1D > ',nuisance_dir,'/mc2.1D']);
    unix(['awk ''{print $3}'' ',fmri,'_mc.1D > ',nuisance_dir,'/mc3.1D']);
    unix(['awk ''{print $4}'' ',fmri,'_mc.1D > ',nuisance_dir,'/mc4.1D']);
    unix(['awk ''{print $5}'' ',fmri,'_mc.1D > ',nuisance_dir,'/mc5.1D']);
    unix(['awk ''{print $6}'' ',fmri,'_mc.1D > ',nuisance_dir,'/mc6.1D']);
    %
    % # Extract signal for global, csf, and wm
    % ## 16. Global
    disp('Extracting global signal for subject');
    unix(['3dmaskave -mask ',fmri,'_mask.nii.gz -quiet ',fmri,'_pp.nii.gz > ',nuisance_dir,'/global.1D'])
    %
    % ## 17. csf
    if FSLRigidReg > 0
        unix(['flirt -ref ',example,'_func.nii.gz -in ',t1,'.pvc.label.nii.gz -out ',fmri,'.pvc.label.nii.gz -applyxfm -init t12',example,'_func.mat -interp nearestneighbour']);
    else
        
%         if isdeployed
%             cmd = sprintf('%s %s %s %s %s %s %s %s', transform_data_affine_bin, [t1,'.pvc.label.nii.gz'], 's', [fmri,'.pvc.label.nii.gz'], [example,'_func.nii.gz'], [t1,'.bfc.nii.gz'], [fmri,'_example_func2t1.rigid_registration_result.mat'], 'nearest');
%             unix(cmd);
%         else
            transform_data_affine([t1,'.pvc.label.nii.gz'], 's', [fmri,'.pvc.label.nii.gz'], [example,'_func.nii.gz'], [t1,'.bfc.nii.gz'], [fmri,'_example_func2t1.rigid_registration_result.mat'], 'nearest');
%         end
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
    unix(['3dmaskave -mask ',fmri,'.csf.mask.nii.gz -quiet ',fmri,'_pp.nii.gz > ',nuisance_dir,'/csf.1D'])
    %
    % ## 18. wm
    disp('Extracting signal from white matter for subject');
    unix(['3dmaskave -mask ',fmri,'.wm.mask.nii.gz -quiet ',fmri,'_pp.nii.gz > ',nuisance_dir,'/wm.1D'])
    %
    % ## 6. Generate mat file (for use later)
    % ## create fsf file
    %
    disp('Modifying model file');
    unix(['sed -e s:nuisance_dir:"',nuisance_dir,'":g <',nuisance_template,' >',nuisance_dir,'/temp1']);
    unix(['sed -e s:nuisance_model_outputdir:"',nuisance_dir,'/residuals.feat":g <',nuisance_dir,'/temp1 >',nuisance_dir,'/temp2']);
    unix(['sed -e s:nuisance_model_TR:"',num2str(TR),'":g <',nuisance_dir,'/temp2 >',nuisance_dir,'/temp3']);
    unix(['sed -e s:nuisance_model_numTRs:"',num2str(n_vols),'":g <',nuisance_dir,'/temp3 >',nuisance_dir,'/temp4']);
    unix(['sed -e s:nuisance_model_input_data:"',func_dir,'/',fmri,'_pp.nii.gz":g <',nuisance_dir,'/temp4 >',nuisance_dir,'/nuisance.fsf']);
    %
    % #rm ${nuisance_dir}/temp*
    %
    disp('Running feat model');
    unix(['feat_model ',nuisance_dir,'/nuisance']);
    %
    
    pp = load_untouch_nii_gz([fmri,'_pp.nii.gz']);
    msk = load_untouch_nii_gz([fmri,'_mask.nii.gz']);
    ppimg = pp.img(logical(msk.img));
    minVal = num2str(min(ppimg));
    
    %[~,minVal]=unix(['3dBrickStat -min -mask ',fmri,'_mask.nii.gz ',fmri,'_pp.nii.gz']);
    %minVal=str2double(minVal);
    %
    % ## 7. Get residuals
    disp('Running film to get residuals');
    unix(['film_gls --rn=',nuisance_dir,'/stats --noest --sa --ms=5 --in=',fmri,'_pp.nii.gz --pd=',nuisance_dir,'/nuisance.mat --thr=',minVal]);
    %
    % ## 8. Demeaning residuals and ADDING 100
    unix(['3dTstat -mean -prefix ',nuisance_dir,'/stats/res4d_mean.nii.gz ',nuisance_dir,'/stats/res4d.nii.gz']);
    if exist([fmri,'_res.nii.gz'],'file')
        delete([fmri,'_res.nii.gz'])
    end
    unix(['3dcalc -a ',nuisance_dir,'/stats/res4d.nii.gz -b ',nuisance_dir,'/stats/res4d_mean.nii.gz -expr ''(a-b)+100'' -prefix ',fmri,'_res.nii.gz']);
    
    RS_infile = [fmri,'_res.nii.gz'];
    BFP_outfile = [fmri,'_res2standard.nii.gz'];
else
    if str2double(config.RunDetrend) > 0
        RS_infile = detrendfile;
        BFP_outfile = [fmri,'_pp2standard.nii.gz'];
    else
        RS_infile = gmsfile;
        BFP_outfile = [fmri,'_gms2standard.nii.gz'];
    end
end
%% Resampling residuals to standard MNI space
if FSLRigidReg > 0
    unix(['flirt -ref ',func_dir,'/standard -in ',RS_infile,' -out ',BFP_outfile,' -applyxfm -init ',func_dir,'/',example,'_func2standard.mat -interp trilinear']);
else
%     if isdeployed
%         cmd = sprintf('%s %s %s %s %s %s %s %s', transform_data_affine_bin, RS_infile, 'm', BFP_outfile, [example,'_func.nii.gz'], 'standard.nii.gz', [fmri,'_example_func2standard.rigid_registration_result.mat'], 'linear');
%         unix(cmd)
%     else
        transform_data_affine(RS_infile, 'm', BFP_outfile, [example,'_func.nii.gz'], 'standard.nii.gz', [fmri,'_example_func2standard.rigid_registration_result.mat'], 'linear');
%     end
end
cd(cwd);


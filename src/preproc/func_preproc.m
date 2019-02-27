
function BFP_outfile = func_preproc(BFPPATH,t1,fmri,func_dir,TR,nuisance_template, FWHM, hp,lp,FSLRigidReg, config)

%##########################################################################################################################
%## SCRIPT TO PREPROCESS THE FUNCTIONAL SCAN
%## parameters are passed from 0_preprocess.sh
%## This is based on batch_process.sh script
%## Written by the Underpants Gnomes (a.k.a Clare Kelly, Zarrar Shehzad, Maarten Mennes & Michael Milham)
%## for more information see www.nitrc.org/projects/fcon_1000
%##
%##########################################################################################################################


[~,n_vols]=unix(['3dinfo -nv ',fmri,'.nii.gz']); n_vols=str2double(n_vols);
TRstart=0;
TRend=n_vols-1;
sigma=str2double(FWHM)/2.3548;

usc_rigid_reg_bin = fullfile(BFPPATH, 'usc_rigid_reg.sh');
transform_data_affine_bin = fullfile(BFPPATH, 'transform_data_affine.sh');

[pth,bname,ext]=fileparts(fmri);example=[bname,'_example'];
%
% ## directory setup
% ##########################################################################################################################
% ##---START OF SCRIPT----------------------------------------------------------------------------------------------------##
% ##########################################################################################################################
%
disp('---------------------------------------');
disp('BFP fMRI PREPROCESSING !');
disp('---------------------------------------');

cwd = pwd;
cd(func_dir);



% ## 1. Dropping first # TRS
disp('Dropping first TRs');
unix(['3dcalc -a ',fmri,'.nii.gz[',num2str(TRstart),'..',num2str(TRend),'] -expr ''a'' -prefix ',fmri,'_dr.nii.gz']);
%
% ##2. Deoblique
disp('Deobliquing');
unix(['3drefit -deoblique ',fmri,'_dr.nii.gz']);
%
% ##3. Reorient into fsl friendly space (what AFNI calls RPI)
disp('Reorienting');
unix(['3dresample -orient RPI -inset ',fmri,'_dr.nii.gz -prefix ',fmri,'_ro.nii.gz']);
%
% ##4. Motion correct to average of timeseries
disp('Motion correcting');
unix(['3dTstat -mean -prefix ',fmri,'_ro_mean.nii.gz ',fmri,'_ro.nii.gz']);
unix(['3dvolreg -Fourier -twopass -base ',fmri,'_ro_mean.nii.gz -zpad 4 -prefix ',fmri,'_mc.nii.gz -1Dfile ',fmri,'_mc.1D ',fmri,'_ro.nii.gz']);
%
% ##6. Get eighth image for use in registration
unix(['3dcalc -a ',fmri,'_mc.nii.gz[7] -expr ''a'' -prefix ',example,'_func.nii.gz']);

% ## 12.FUNC->T1
% ## You may want to change some of the options
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
    if isdeployed
        cmd = sprintf('%s %s %s %s %s %s', usc_rigid_reg_bin, moving_filename, static_filename, output_filename, 'inversion');
        unix(cmd);
    else
        usc_rigid_reg(moving_filename, static_filename, output_filename, 'inversion');
    end
end
% ##5. Remove skull/edge detect
disp('Skull stripping');
if str2double(config.T1mask) > 0
    if FSLRigidReg > 0
        unix(['flirt -ref ',example,'_func.nii.gz -in ',t1,'.mask.nii.gz -out ',fmri,'_mask.nii.gz -applyxfm -init t12',example,'_func.mat -interp nearestneighbour']);
    else
        if isdeployed
            cmd = sprintf('%s %s %s %s %s %s %s %s', transform_data_affine_bin, [t1,'.mask.nii.gz'], 's', [fmri,'_mask.nii.gz'], [example,'_func.nii.gz'], [t1,'.bfc.nii.gz'], [fmri,'_example_func2t1.rigid_registration_result.mat'], 'nearest');
            unix(cmd);
        else
            transform_data_affine([t1,'.mask.nii.gz'], 's', [fmri,'_mask.temp.nii.gz'], [example,'_func.nii.gz'], [t1,'.bfc.nii.gz'], [fmri,'_example_func2t1.rigid_registration_result.mat'], 'nearest');
            unix(['3dcalc -a ',fmri,'_mask.temp.nii.gz -expr ''a/255'' -prefix ',fmri,'_mask.nii.gz']);
            unix(['rm ',fmri,'_mask.temp.nii.gz'])
        end
    end
else
    unix(['3dAutomask -prefix ',fmri,'_mask.nii.gz -dilate 1 ',example,'_func.nii.gz']);
end
unix(['3dcalc -a ',fmri,'_mc.nii.gz -b ',fmri,'_mask.nii.gz -expr ''a*b'' -prefix ',fmri,'_ss.nii.gz']);

%
% ##6. Get eighth image for use in registration
disp( 'Getting example_func for registration')
unix(['rm ',example,'_func.nii.gz'])
unix(['3dcalc -a ',fmri,'_ss.nii.gz[7] -expr ''a'' -prefix ',example,'_func.nii.gz']);

%
% ##7. Spatial smoothing
disp('Spatial Smoothing');
unix(['fslmaths ',fmri,'_ss.nii.gz -kernel gauss ',num2str(sigma),' -fmean -mas ',fmri,'_mask.nii.gz ',fmri,'_sm.nii.gz']);
%
% ##8. Grandmean scaling
disp('Grand-mean scaling');
unix(['fslmaths ',fmri,'_sm.nii.gz -ing 10000 ',fmri,'_gms.nii.gz -odt float']);
%

% ##9. Temporal filtering
disp('Band-pass filtering');
gmsfile = [fmri,'_gms.nii.gz'];
unix(['3dFourier -lowpass ',num2str(lp),' -highpass ',num2str(hp),' -retrend -prefix ',fmri,'_filt.nii.gz ',gmsfile]);
%

if str2double(config.RunDetrend) > 0
    % ##10.Detrending
    disp('Removing linear and quadratic trends');
    unix(['3dTstat -mean -prefix ',fmri,'_filt_mean.nii.gz ',fmri,'_filt.nii.gz']);
    unix(['3dDetrend -polort 2 -prefix ',fmri,'_dt.nii.gz ',fmri,'_filt.nii.gz']);
    detrendfile = [fmri,'_pp.nii.gz'];
    unix(['3dcalc -a ',fmri,'_filt_mean.nii.gz -b ',fmri,'_dt.nii.gz -expr ''a+b'' -prefix ',detrendfile]);
    maskstep_infile = detrendfile;
else
    maskstep_infile = gmsfile;
end
%
% ##11.Create Mask
disp('Generating mask of preprocessed data');
unix(['fslmaths ',maskstep_infile,' -Tmin -bin ',fmri,'_pp_mask.nii.gz -odt char']);
%

% ## 12.FUNC->standard (3mm)
% ## You may want to change some of the options
if FSLRigidReg > 0
    unix(['flirt -ref standard.nii.gz -in ',example,'_func.nii.gz -out ',example,'_func2standard.nii.gz -omat ',example,'_func2standard.mat -cost corratio -dof 6 -interp trilinear']);
    % # Create mat file for conversion from subject's anatomical to functional
    unix(['convert_xfm -inverse -omat standard2',example,'_func.mat ',example,'_func2standard.mat']);
else
    disp('Using USC rigid registration');
    moving_filename = [example,'_func.nii.gz'];%fullfile(funcDir,sprintf('%s_%s_bold_example_func.nii.gz',subid,sessionid{ind}));
    static_filename = 'standard.nii.gz';
    output_filename = [example,'_func2standard.nii.gz'];
    if isdeployed
        cmd = sprintf('%s %s %s %s %s %s', usc_rigid_reg_bin, moving_filename, static_filename, output_filename, 'inversion');
        unix(cmd);
    else
        usc_rigid_reg(moving_filename, static_filename, output_filename, 'inversion');
    end
end

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
    unix(['mkdir -p ',nuisance_dir]);
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
    disp(['Extracting global signal for subject']);
    unix(['3dmaskave -mask ',fmri,'_pp_mask.nii.gz -quiet ',fmri,'_pp.nii.gz > ',nuisance_dir,'/global.1D'])
    %
    % ## 17. csf
    if FSLRigidReg > 0
        unix(['flirt -ref ',example,'_func.nii.gz -in ',t1,'.pvc.label.nii.gz -out ',t1,'.func.pvc.label.nii.gz -applyxfm -init t12',example,'_func.mat -interp nearestneighbour']);
    else
        
        if isdeployed
            cmd = sprintf('%s %s %s %s %s %s %s %s', transform_data_affine_bin, [t1,'.pvc.label.nii.gz'], 's', [t1,'.func.pvc.label.nii.gz'], [example,'_func.nii.gz'], [t1,'.bfc.nii.gz'], [fmri,'_example_func2t1.rigid_registration_result.mat'], 'nearest');
            unix(cmd);
        else
            transform_data_affine([t1,'.pvc.label.nii.gz'], 's', [t1,'.func.pvc.label.nii.gz'], [example,'_func.nii.gz'], [t1,'.bfc.nii.gz'], [fmri,'_example_func2t1.rigid_registration_result.mat'], 'nearest');
        end
    end
    %
    unix(['fslmaths ',t1,'.func.pvc.label.nii.gz -thr 5.5 -bin ',t1,'.func.csf.mask.nii.gz']);
    unix(['fslmaths ',t1,'.func.pvc.label.nii.gz -thr 2.5 -uthr 3.5 -bin ',t1,'.func.wm.mask.nii.gz']);
    %
    % echo "Extracting signal from csf"
    unix(['3dmaskave -mask ',t1,'.func.csf.mask.nii.gz -quiet ',fmri,'_pp.nii.gz > ',nuisance_dir,'/csf.1D'])
    %
    % ## 18. wm
    disp('Extracting signal from white matter for subject');
    unix(['3dmaskave -mask ',t1,'.func.wm.mask.nii.gz -quiet ',fmri,'_pp.nii.gz > ',nuisance_dir,'/wm.1D'])
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
    [~,minVal]=unix(['3dBrickStat -min -mask ',fmri,'_pp_mask.nii.gz ',fmri,'_pp.nii.gz']);
    %minVal=str2double(minVal);
    %
    % ## 7. Get residuals
    disp('Running film to get residuals');
    unix(['film_gls --rn=',nuisance_dir,'/stats --noest --sa --ms=5 --in=',fmri,'_pp.nii.gz --pd=',nuisance_dir,'/nuisance.mat --thr=',minVal]);
    %
    % ## 8. Demeaning residuals and ADDING 100
    unix(['3dTstat -mean -prefix ',nuisance_dir,'/stats/res4d_mean.nii.gz ',nuisance_dir,'/stats/res4d.nii.gz']);
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

%
% ## 9. Resampling residuals to MNI space
if FSLRigidReg > 0
    unix(['flirt -ref ',func_dir,'/standard -in ',RS_infile,' -out ',BFP_outfile,' -applyxfm -init ',func_dir,'/',example,'_func2standard.mat -interp trilinear']);
else
    if isdeployed
        cmd = sprintf('%s %s %s %s %s %s %s %s', transform_data_affine_bin, RS_infile, 'm', BFP_outfile, [example,'_func.nii.gz'], 'standard.nii.gz', [fmri,'_example_func2standard.rigid_registration_result.mat'], 'linear');
        unix(cmd)
    else
        transform_data_affine(RS_infile, 'm', BFP_outfile, [example,'_func.nii.gz'], 'standard.nii.gz', [fmri,'_example_func2standard.rigid_registration_result.mat'], 'linear');
    end
end
%
cd(cwd);


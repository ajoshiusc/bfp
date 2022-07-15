clear all;close all;clc;
addpath(genpath('/home/ajoshi/projects/bfp/src'))
fmri = '/home/ajoshi/sub08001/func/sub08001_rest_bold';
subbasename = '/home/ajoshi/sub08001/anat/sub08001_T1w';
config.FSLPATH = '/home/ajoshi/webware/fsl'
config.FSLOUTPUTTYPE='NIFTI_GZ';
config.AFNIPATH = '/home/ajoshi/abin';
config.FSLRigidReg=0;
config.MultiThreading=1;
BFPPATH='/home/ajoshi/projects/bfp';
[~,subid] = fileparts(subbasename);
sessionid = 'rest';



GordSize=96854;
setenv('PATH', [getenv('PATH'),':',config.FSLPATH,':',config.FSLPATH,'/bin']);
setenv('FSLDIR', config.FSLPATH);
setenv('FSLOUTPUTTYPE',config.FSLOUTPUTTYPE);
% some newer afni versions throw warnings for non-float data
% creating parsing errors. This takes care of that.
setenv('AFNI_NIFTI_TYPE_WARN','NO');
setenv('PATH', [getenv('PATH'),':',config.AFNIPATH,':',config.AFNIPATH,'/bin']);

alff_ext={'ALFF','fALFF','ALFF_Z'};

func_dir = fileparts(fmri);
infile = [fmri,'_gms.nii.gz'];
outbase = fmri;

GOrdSurfIndFile=fullfile(BFPPATH,'supp_data','bci_grayordinates_surf_ind.mat');

[~,b] = unix(['3dinfo -tr ',[fmri,'.nii.gz']]); % get TR
TR = num2str(str2double(b));
LP = num2str(0.01);
HP = num2str(0.1);
mask_file = [fmri,'_mask.nii.gz'];

cmd = [fullfile(BFPPATH,'dev/alff/createALFF.sh '), outbase, ' ', infile, ' ',mask_file, ' ', TR, ' ', LP,',',HP];
unix(cmd);

example = [fmri,'.example'];


% Resampling ALFF_images to standard space
for i = 1:length(alff_ext)

    if config.FSLRigidReg > 0
        unix(['flirt -ref ',func_dir,'/standard.nii.gz -in ',fmri,'_',alff_ext{i},'.nii.gz -out ',fmri,'_',alff_ext{i},'2standard.nii.gz -applyxfm -init ',fmri,'_example_func2standard.mat -interp trilinear']);
    else
        transform_data_affine([fmri,'_',alff_ext{i},'.nii.gz'], 'm', [fmri,'_',alff_ext{i},'2standard.nii.gz'], [example,'.func.nii.gz'], fullfile(func_dir,'standard.nii.gz'), [fmri,'.example.func2standard.rigid_registration_result.mat'], 'linear');
    end

    fmri2surfFile=fullfile(func_dir,sprintf('%s_%s_bold2surf.mat',subid,sessionid));
    GOrdFile=fullfile(func_dir,sprintf('%s_%s_bold.%s.GOrd.mat',subid,sessionid,alff_ext{i}));
    resample2surf(subbasename,[fmri,'_',alff_ext{i},'2standard.nii.gz'],fmri2surfFile,config.MultiThreading);

    load(GOrdSurfIndFile,'ind_left','ind_right');

    a=load(fmri2surfFile);
    left_go_data=a.datal_atlas(ind_left);
    right_go_data=a.datar_atlas(ind_left);

    lsz=size(left_go_data,1);
    rsz=size(right_go_data,1);
    data=zeros(GordSize,1);
    data(1:lsz,:)=left_go_data;
    data(1+lsz:lsz+rsz,:)=right_go_data;

    fprintf('Saving file: %s\n', GOrdFile);

    save(GOrdFile,'data');

end


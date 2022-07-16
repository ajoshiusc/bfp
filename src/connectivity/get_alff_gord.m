function get_alff_gord(config, fmribase,anatbase)

GordSize=96854;
setenv('PATH', [getenv('PATH'),':',config.FSLPATH,':',config.FSLPATH,'/bin']);
setenv('FSLDIR', config.FSLPATH);
setenv('FSLOUTPUTTYPE',config.FSLOUTPUTTYPE);
% some newer afni versions throw warnings for non-float data
% creating parsing errors. This takes care of that.
setenv('AFNI_NIFTI_TYPE_WARN','NO');
setenv('PATH', [getenv('PATH'),':',config.AFNIPATH,':',config.AFNIPATH,'/bin']);

alff_ext={'ALFF','fALFF','ALFF_Z'};

func_dir = fileparts(fmribase);
infile = [fmribase,'_gms.nii.gz'];
outbase = fmribase;

GOrdSurfIndFile=fullfile(config.BFPPATH,'supp_data','bci_grayordinates_surf_ind.mat');

[~,b] = unix(['3dinfo -tr ',[fmribase,'.nii.gz']]); % get TR
TR = num2str(str2double(b));
LP = num2str(0.01);
HP = num2str(0.1);
mask_file = [fmribase,'_mask.nii.gz'];

cmd = [fullfile(config.BFPPATH,'src/connectivity/createALFF.sh '), outbase, ' ', infile, ' ',mask_file, ' ', TR, ' ', LP,',',HP];
unix(cmd);

example = [fmribase,'.example'];


% Resampling ALFF_images to standard space
for i = 1:length(alff_ext)

    if config.FSLRigidReg > 0
        unix(['flirt -ref ',func_dir,'/standard.nii.gz -in ',fmribase,'_',alff_ext{i},'.nii.gz -out ',fmribase,'_',alff_ext{i},'2standard.nii.gz -applyxfm -init ',fmribase,'_example_func2standard.mat -interp trilinear']);
    else
        transform_data_affine([fmribase,'_',alff_ext{i},'.nii.gz'], 'm', [fmribase,'_',alff_ext{i},'2standard.nii.gz'], [example,'.func.nii.gz'], fullfile(func_dir,'standard.nii.gz'), [fmribase,'.example.func2standard.rigid_registration_result.mat'], 'linear');
    end

    fmribase2surfFile=fullfile(func_dir,sprintf('data2surf.mat'));
    GOrdFile=fullfile(sprintf('%s_bold.%s.GOrd.mat',fmribase,alff_ext{i}));
    resample2surf(anatbase,[fmribase,'_',alff_ext{i},'2standard.nii.gz'],fmribase2surfFile,config.MultiThreading);

    load(GOrdSurfIndFile,'ind_left','ind_right');

    a=load(fmribase2surfFile);
    left_go_data=a.datal_atlas(ind_left);
    right_go_data=a.datar_atlas(ind_right);

    lsz=size(left_go_data,1);
    rsz=size(right_go_data,1);
    data=zeros(GordSize,1);
    data(1:lsz,:)=left_go_data;
    data(1+lsz:lsz+rsz,:)=right_go_data;

    fprintf('Saving file: %s\n', GOrdFile);

    save(GOrdFile,'data');
    delete(fmribase2surfFile);

end
end
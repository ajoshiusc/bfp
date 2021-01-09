function fMRI_epicorr_t1based(epi_filename, ref_filename, fmri,t1,opts)
%   epi_filename: input epi file
%   ref_filename: reference epi volume that will be used for correction
%   fmri: fmri basename (BFP BIDS format)
%   t1: t1 basename (BrainSuite BIDS format)
%   opts: see preprocessing config file from BFP 
%%
[example, ext] = remove_extension(ref_filename);
struct_filename = [t1,'.bfc.nii.gz']; %t1 bfc file
%%  directory for intermediate files
[~,fmribase,~] = fileparts(fmri);
b0dir=[fmri,'_b0'];
fmri_b0 = fullfile(b0dir,fmribase);
mkdir(b0dir);
%%  log file
logfname=[fmri_b0,'.log.txt'];
fp=fopen(logfname,'a+');
t = now;
d = datetime(t,'ConvertFrom','datenum');
fprintf(fp,'\n%s\n',d);
fprintf(fp, 'Running T1-based distortion correction \n');
%%  initial rigid registration
regmat = [fmri,'.example.func2t1.rigid_registration_result.mat'];
if ~exist(regmat,'file')
    msg = ['Initial registration file not found: ', regmat];
    disp(msg);
    fprintf(fp, '%s\n',msg);
    msg = 'Running initial rigid registration using INVERSION';
    disp(msg);
    fprintf(fp, '%s\n',msg);
    moving_filename = ref_filename;
    static_filename = [t1,'.bfc.nii.gz'];
    output_filename = [example,'.func2t1.nii.gz'];
    opts.similarity = 'inversion';
    register_files_affine(moving_filename, static_filename, output_filename, opts);
    clear opts
end
%%  reorient reference/example file
refRAS = [remove_extension(ref_filename),'.RAS.nii.gz'];
if ~exist(refRAS,'file')
    if exist(ref_filename,'file')
        reorient_nifti_sform(ref_filename,refRAS);
        msg = ['reference file reoriented to RAS format: ', refRAS];
        disp(msg);
        fprintf(fp, '%s\n',msg);
    else
        msg = ['error! reference file not found: ', ref_filename];
        disp(msg)
        fprintf(fp, '%s\n', msg);
        error(msg);

    end
else
    msg = ['reference file already reoriented'];
    disp(msg)
    fprintf(fp, '%s\n', msg);
end
%%  reorient mc/epi file to RAS format
epiRAS = [remove_extension(epi_filename),'.RAS.nii.gz'];
if ~exist(refRAS,'file')
    if exist(epi_filename,'file')
        reorient_nifti_sform(epi_filename,epiRAS);
        msg = ['input epi file reoriented to RAS format: ', epiRAS];
        disp(msg);
        fprintf(fp, '%s\n',msg);
    else
        msg = ['error! input epi file not found: ', epi_filename];
        disp(msg)
        fprintf(fp, '%s\n', msg);
        error(msg);

    end
else
    msg = ['input epi file already reoriented'];
    disp(msg)
    fprintf(fp, '%s\n', msg);
end
%%  mask epi
if opts.maskepi
    transform_data_affine([t1,'.mask.nii.gz'], 's', [fmri_b0,'.mask.temp.nii.gz'], ...
        [example,'.func.nii.gz'], [t1,'.bfc.nii.gz'], [fmri,'.example.func2t1.rigid_registration_result.mat'], 'nearest');
    fprintf(fp, '--Option T1: Skull Strip fMRI\n');
    
    msk = load_untouch_nii_gz([fmri,'.mask.temp.nii.gz']);
    msk.img(msk.img>0)=1;
    save_untouch_nii_gz(msk,[fmri,'.mask.nii.gz']);
    unix(['rm ',fmri,'.mask.temp.nii.gz'])
    
    fprintf(fp, '--Masking method: BDP liberal \n');
    output_mask_file = [fmri_b0,'.ro.mean.mask.nii.gz'];
    if_aggresive = 0;
    maskHeadPseudoHist(ref_filename, output_mask_file, if_aggresive)
else
    fprintf(fp, '--Masking method: 3dAutomask dilate 1 \n');
    output_mask_file = [fmri_b0,'.ro.mean.mask.nii.gz'];
    unix(['3dAutomask -prefix ', output_mask_file,' -dilate 1 ',fmri,'.ro.mean.RAS.nii.gz']);
end
%%  create temp volume with ref as 1st volume
outfile = [fmri_b0,'.corr.ref.nii.gz'];

disp('creating temp reference volume for t1 based distortion correction');
fprintf(fp, 'input file <.corr.ref.nii.gz> concatenates reference volume with motion corrected data \n');
fprintf(fp, '\t reference volume will be removed in final output \n');

epi = load_untouch_nii_gz(epi_filename);
ref = load_untouch_nii_gz(ref_filename);
epi_temp = epi;
sz = size(epi.img);
sz(4) = sz(4) + 1;
epi_temp.img = zeros(sz);
epi_temp.img(:,:,:,1) = ref.img;
epi_temp.img(:,:,:,2:end) = epi.img;
epi_temp.hdr.dime.dim(2:5) = sz;
if sz(4)>1
    epi_temp.hdr.dime.dim(1) = 4;
else
    epi_temp.hdr.dime.dim(1) = 3;
end
save_untouch_nii_gz(epi_temp,outfile);

%%  T1-based distortion correction
fprintf(fp, 'Running T1-based distortion correction \n');
infile = outfile;
outfile = [fmri_b0,'.corr.refout.nii.gz'];
struct_output_filename = [fmri_b0,'.t1.corr.nii.gz'];

opts.epi_mask = output_mask_file;
opts.struct_mask = [t1,'.mask.nii.gz'];
opts.rigid_reg_mat = regmat;

opts.num_threads = 60;
fprintf(fp, '--number of threads: %d \n', opts.num_threads);

opts.non_uniformity_correction = false;
fprintf(fp, '--bias field correction: %s \n', mat2str(opts.non_uniformity_correction));
% opts.similarity = 'inversion';

% opts.rigid_similarity = 'inversion';
% fprintf(fp, '--rigid_similarity: %s \n', opts.rigid_similarity);

% opts.reg_res = 1.8;
% fprintf(fp, '--reg_res: %s \n', opts.reg_res);

EPI_correct_files_registration_INVERSION(infile, struct_filename, ...
    outfile, struct_output_filename, opts);
fprintf(fp, 'Successful! Corrected file is %s \n', outfile);
clear opts
%%  remove reference volume
infile = outfile;
outfile = [fmri,'.corr.nii.gz'];
epi_temp = load_untouch_nii_gz(infile);
epi_corr = epi_temp;
epi_corr.img(:,:,:,1) = [];
epi_corr.hdr.dime.dim(5) =epi_corr.hdr.dime.dim(5)-1;
if sz(4)>1
    epi_corr.hdr.dime.dim(1) = 4;
else
    epi_corr.hdr.dime.dim(1) = 3;
end
save_untouch_nii_gz(epi_corr,outfile);
fprintf(fp, 'Final Output: %s', outfile);
%% for testing: reregister to t1
v = lung([fmri_b0,'.corr.refout.nii.gz']);
v.img(:,:,:,2:end) = [];
v.hdr.dime.dim(1) = 3;
v.hdr.dime.dim(5) = 1;
outfile = [fmri_b0,'.corr.refvol.nii.gz'];
save_untouch_nii_gz(v,outfile)

moving_filename = outfile;
static_filename = [t1,'.bfc.nii.gz'];
output_filename = [fmri_b0,'.corr.func2t1.nii.gz'];
clear opts
opts.similarity = 'inversion';
opts.nthreads = 60;
register_files_affine(moving_filename, static_filename, output_filename, opts);

transform_data_affine([t1,'.pvc.label.nii.gz'], 's', [fmri_b0,'.pvc.label.nii.gz'], ...
    ref_filename, [t1,'.bfc.nii.gz'], [fmri_b0,'.corr.func2t1.rigid_registration_result.mat'], 'nearest');
%%
transform_data_affine(outfile, 'm', [fmri_b0,'.corr.func2t1.origreg.nii.gz'], ...
      [fmri_b0,'.corr.ref.0_diffusion.nii.gz'], [t1,'.bfc.nii.gz'], regmat, 'spline');

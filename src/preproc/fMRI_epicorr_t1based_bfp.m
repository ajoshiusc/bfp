function fMRI_epicorr_t1based_bfp(epi_filename, ref_filename, fmri,t1,config)
%   epi_filename: input epi file
%   ref_filename: reference epi volume that will be used for correction
%   fmri: fmri basename (BFP BIDS format)
%   t1: t1 basename (BrainSuite BIDS format)
%   config: see preprocessing config file from BFP 
%%
[example, ext] = remove_extension(ref_filename);
struct_filename = [t1,'.bfc.nii.gz']; %t1 bfc file
%%  directory for intermediate files
[~,fmribase,~] = fileparts(fmri);
b0dir=[fmri,'_b0'];
fmri_b0 = fullfile(b0dir,fmribase);
mkdir(b0dir);
%%  log file
if ~exist('config.log_filename','var')
    logfname = [fmri_b0,'.log.txt'];
else
    logfname = config.log_filename;
end
fp=fopen(logfname,'a+');
t = now;
d = datetime(t,'ConvertFrom','datenum');
fprintf(fp,'\n%s\n',d);
fprintf(fp, 'Running T1-based distortion correction \n');
%%  reorient reference/example file
% refRAS = [fmri_b0,'.ref.RAS.nii.gz'];
% if ~exist(refRAS,'file')
%     if exist(ref_filename,'file')
%         reorient_nifti_sform(ref_filename,refRAS);
%         msg = ['reference file reoriented to RAS format: ', refRAS];
%         disp(msg);
%         fprintf(fp, '%s\n',msg);
%     else
%         msg = ['error! reference file not found: ', ref_filename];
%         disp(msg)
%         fprintf(fp, '%s\n', msg);
%         error(msg);
% 
%     end
% else
%     msg = ['reference file already reoriented'];
%     disp(msg)
%     fprintf(fp, '%s\n', msg);
% end
%%  reorient mc/epi file to RAS format
% epiRAS = [remove_extension(epi_filename),'.RAS.nii.gz'];
% if ~exist(epiRAS,'file')
%     if exist(epi_filename,'file')
%         reorient_nifti_sform(epi_filename,epiRAS);
%         msg = ['input epi file reoriented to RAS format: ', epiRAS];
%         disp(msg);
%         fprintf(fp, '%s\n',msg);
%     else
%         msg = ['error! input epi file not found: ', epi_filename];
%         disp(msg)
%         fprintf(fp, '%s\n', msg);
%         error(msg);
% 
%     end
% else
%     msg = ['input epi file already reoriented'];
%     disp(msg)
%     fprintf(fp, '%s\n', msg);
% end
%%  create temp volume with ref as 1st volume
infile = [fmri_b0,'.infile.nii.gz'];

if ~exist(infile,'file')
    disp('creating temp reference volume for t1 based distortion correction');
    fprintf(fp, 'input file concatenates reference volume as first volume with motion corrected data: %s\n',infile);
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
    save_untouch_nii_gz(epi_temp,infile);
else
    disp('File found: temp reference volume for t1 based distortion correction');
    fprintf(fp, 'input file concatenates reference volume as first volume with motion corrected data: %s\n',infile);
end
%% reorient input file
infileRAS = [fmri_b0,'.infile.RAS.nii.gz'];
if ~exist(infileRAS,'file')
    reorient_nifti_sform(infile,infileRAS);
end
% t1RAS = [fmri_b0,'.t1.bfc.RAS.nii.gz'];
% reorient_nifti_sform(struct_filename,t1RAS);
%%  mask epi
epimask_filename = [fmri_b0,'.mask.nii.gz'];
if config.epit1corr_mask == '0'
    fprintf(fp, '--Masking method: t1-based mask \n');
    copyfile([fmri,'.mask.nii.gz'],epimask_filename);    
elseif config.epit1corr_mask=='1'
    fprintf(fp, '--Masking method: BDP liberal \n');
    if_aggresive = 0;
    maskHeadPseudoHist(infileRAS, epimask_filename, if_aggresive);
elseif config.epit1corr_mask=='2'
    fprintf(fp, '--Masking method: BDP aggressive \n');
    if_aggresive = 1;
    maskHeadPseudoHist(ref_filename, epimask_filename, if_aggresive);
elseif config.epit1corr_mask=='3'
    if exist(epimask_filename,'file')
        delete(epimask_filename);
    end
    fprintf(fp, '--Masking method: AFNI 3dAutomask dilate 2 \n');
    unix(['3dAutomask -prefix ', epimask_filename,' -dilate 2 ',ref_filename]);
end
reorient_nifti_sform(epimask_filename,epimask_filename);
%%  initial rigid registration
output_filename = [fmri_b0,'.ref2t1.nii.gz'];
regmat = [remove_extension(output_filename),'.rigid_registration_result.mat'];
if ~exist(regmat,'file')
    msg = 'Running initial rigid registration';
    disp(msg);
    fprintf(fp, '%s\n',msg);
    msg = 'Running initial rigid registration using INVERSION';
    disp(msg);
    fprintf(fp, '%s\n',msg);
    moving_filename = ref_filename;
    static_filename = [t1,'.bfc.nii.gz'];
    opts.moving_mask = epimask_filename;   
    opts.static_mask = [t1,'.mask.nii.gz'];
    
    opts.similarity = 'inversion';
    opts.nthreads = str2num(config.epit1corr_numthreads);
    register_files_affine(moving_filename, static_filename, output_filename, opts);
    clear opts
end
%%  T1-based distortion correction
fprintf(fp, 'Running T1-based distortion correction \n');
outfile_corr = [fmri_b0,'.infile.corr.nii.gz'];
struct_output_filename = [fmri_b0,'.t1.corr.nii.gz'];

opts.epi_mask = epimask_filename;
opts.struct_mask = [t1,'.mask.nii.gz'];
% opts.rigid_reg_mat = regmat;

opts.num_threads = str2num(config.epit1corr_numthreads);
fprintf(fp, '--number of threads: %d \n', opts.num_threads);

% opts.non_uniformity_correction = false;
% fprintf(fp, '--bias field correction: %s \n', mat2str(opts.non_uniformity_correction));
% opts.similarity = 'inversion';

% opts.rigid_similarity = 'inversion';
% fprintf(fp, '--rigid_similarity: %s \n', opts.rigid_similarity);

% opts.reg_res = 1.8;
% fprintf(fp, '--reg_res: %s \n', opts.reg_res);

% opts.debug = true;

EPI_correct_files_registration_INVERSION(infileRAS, struct_filename, ...
    outfile_corr, struct_output_filename, opts);
fprintf(fp, 'Successful! Corrected file is %s \n', outfile_corr);
clear opts
%%  remove reference volume
outfile_final = [fmri,'.corr.nii.gz'];
epi_temp = load_untouch_nii_gz(outfile_corr);
epi_corr = epi_temp;
epi_corr.img(:,:,:,1) = [];
epi_corr.hdr.dime.dim(5) =epi_corr.hdr.dime.dim(5)-1;
if sz(4)>1
    epi_corr.hdr.dime.dim(1) = 4;
else
    epi_corr.hdr.dime.dim(1) = 3;
end
save_untouch_nii_gz(epi_corr,outfile_final);
fprintf(fp, 'Final Output: %s', outfile_final);
%% for testing: reregister to t1
v = lung(outfile_corr);
v.img(:,:,:,2:end) = [];
v.hdr.dime.dim(1) = 3;
v.hdr.dime.dim(5) = 1;
outfile = [fmri_b0,'.infile.corr.refvol.nii.gz'];
save_untouch_nii_gz(v,outfile)

moving_filename = outfile;
opts.moving_mask = epimask_filename;  
static_filename = [t1,'.bfc.nii.gz'];
output_filename = [fmri_b0,'.corr.func2t1.nii.gz'];
opts.nthreads = str2num(config.epit1corr_numthreads);

clear opts
opts.similarity = 'inversion';
register_files_affine(moving_filename, static_filename, output_filename, opts);

transform_data_affine([t1,'.pvc.label.nii.gz'], 's', [fmri_b0,'.pvc.label.nii.gz'], ...
    ref_filename, [t1,'.bfc.nii.gz'], [fmri_b0,'.corr.func2t1.rigid_registration_result.mat'], 'nearest');
%%
transform_data_affine(outfile, 'm', [fmri_b0,'.corr.func2t1.origreg.nii.gz'], ...
      [fmri_b0,'.corr.ref.0_diffusion.nii.gz'], [t1,'.bfc.nii.gz'], regmat, 'spline');

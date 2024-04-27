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
    save_untouch_nii_gz(epi_temp,infile,64);
else
    disp('File found: temp reference volume for t1 based distortion correction');
    fprintf(fp, 'input file concatenates reference volume as first volume with motion corrected data: %s\n',infile);
end
%% reorient input file
infileRAS = [fmri_b0,'.infile.RAS.nii.gz'];
if ~exist(infileRAS,'file')
    reorient_nifti_sform(infile,infileRAS);
end
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
% reorient_nifti_sform(epimask_filename,epimask_filename);
%%  T1-based distortion correction
fprintf(fp, 'Running T1-based distortion correction \n');
outfile_corr = [fmri_b0,'.infile.corr.nii.gz'];
struct_output_filename = [fmri_b0,'.t1.corr.nii.gz'];

opts.epi_mask = epimask_filename;
opts.struct_mask = [t1,'.mask.nii.gz'];

opts.num_threads = str2num(config.epit1corr_numthreads);
fprintf(fp, '--number of threads: %d \n', opts.num_threads);

fprintf(fp, '--initialize rigid registration method: %s \n', config.epit1corr_rigidsim);
opts.rigid_similarity = config.epit1corr_rigidsim;

fprintf(fp, '--bias field correction: %s \n', mat2str(config.epit1corr_bias));
if config.epit1corr_bias == '0'
    opts.non_uniformity_correction = false;
elseif config.epit1corr_bias == '1'
    opts.non_uniformity_correction = true;
end

EPI_correct_files_registration_INVERSION(infileRAS, struct_filename, ...
    outfile_corr, struct_output_filename, opts);
fprintf(fp, 'Successful! Corrected file is %s \n', outfile_corr);
clear opts
%%  remove reference volume
outfile_final = [fmri,'.epicorr.nii.gz'];
epi_temp = load_untouch_nii_gz(outfile_corr);
epi_corr = epi_temp;
epi_corr.img(:,:,:,1) = [];
epi_corr.hdr.dime.dim(5) =epi_corr.hdr.dime.dim(5)-1;
sz = size(epi_corr.img);
if sz(4)>1
    epi_corr.hdr.dime.dim(1) = 4;
else
    epi_corr.hdr.dime.dim(1) = 3;
end
save_untouch_nii_gz(epi_corr,outfile_final);
fprintf(fp, 'Final Output: %s', outfile_final);
end

function [out_nii, file_out] = fixBSheader(hdr_src_file, infile, output_file)
% A hack for (inconsistent) headers of BrainSuite/SVReg files. If input files, hdr_src_file and
% infile, overlay correctly in BrainSuite then this function modifies header of infile such that it
% is consistent with the hdr_src_file and still overlays correctly in BrainSuite. Following steps
% are followed: 
%  1. check hdr_src_file for correct sform and cannonical sform header
%  2. reorient image of infile to make sform header cannonical 
%  3. Perform some basic sanity checks on dimension etc.. 
%
% Note that this function in NOT completely generalized in nature because of following assumptions
% about hdr_src_file. These assumptions are same as that in nifti files written by BrainSuite (&
% hence hdr_src_file should ideally be a file written by BrainSuite). 
%   1. hdr_src_file must have valid sform matrix. 
%   2. sform matrix must be in canonical form 
%

workDir = tempname();
mkdir(workDir);

flg = check_nifti_file(hdr_src_file);
if flg == 2
   error('BDP:InvalidBSheader','%s must have valid sform matrix.', escape_filename(hdr_src_file))
end

[~, M] = reorient_nifti_sform(hdr_src_file, fullfile(workDir, [randstr(10) '.nii.gz']));
if ~isequal(diag(abs(diag(M))),M)
   error('BDP:InvalidBSheader','%s must have cannonical sform header.', escape_filename(hdr_src_file))
end
hdr_src_nii = load_untouch_nii_gz(hdr_src_file);

[~, outFile] = check_nifti_file(infile, workDir);
out_nii = reorient_nifti_sform(outFile, fullfile(workDir, [randstr(10) '.nii.gz']));

out_nii.hdr.hist.sform_code = hdr_src_nii.hdr.hist.sform_code;
out_nii.hdr.hist.srow_x = hdr_src_nii.hdr.hist.srow_x;
out_nii.hdr.hist.srow_y = hdr_src_nii.hdr.hist.srow_y;
out_nii.hdr.hist.srow_z = hdr_src_nii.hdr.hist.srow_z;
out_nii.hdr.hist.qform_code = 0; % unset qform_code

% throw away scl_slope & scl_inter as BrainSuite does not use it
out_nii.hdr.dime.scl_slope = 0;
out_nii.hdr.dime.scl_inter = 0;

% check for possible inconsistency
if norm(out_nii.hdr.dime.dim(2:4)-hdr_src_nii.hdr.dime.dim(2:4)) > 1e-3   
   disp(['out.hdr.dime.dim: ' num2str(out_nii.hdr.dime.dim)])
   disp(['header_source.hdr.dime.dim: ' num2str(hdr_src_nii.hdr.dime.dim)])
   error('BDP:InvalidBSheader','hdr.dime.dim does not match for %s and %s', ...
      escape_filename(hdr_src_file), escape_filename(infile));
end

if norm(out_nii.hdr.dime.pixdim(2:4)-hdr_src_nii.hdr.dime.pixdim(2:4)) > 1e-3   
   disp(['out.hdr.dime.pixdim: ' num2str(out_nii.hdr.dime.pixdim)])
   disp(['header_source.hdr.pixdime.dim: ' num2str(hdr_src_nii.hdr.dime.pixdim)])
   error('BDP:InvalidBSheader','hdr.dime.pixdim does not match for %s and %s',...
      escape_filename(hdr_src_file), escape_filename(infile));
end

if ~isequal(size(out_nii.img(:,:,:,1)), size(hdr_src_nii.img(:,:,:,1)))
   disp(['size(out.img) = ' num2str(size(out_nii.img))])
   disp(['size(hdr_src_nii.img) = ' num2str(size(hdr_src_nii.img))])
   error('BDP:InvalidBSheader','Image size does not match for %s and %s', ...
      escape_filename(hdr_src_file), escape_filename(infile));
end

if nargin==3
   file_out = save_untouch_nii_gz(out_nii, output_file, workDir);
end

rmdir(workDir, 's');
end

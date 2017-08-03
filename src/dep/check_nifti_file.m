% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2017 The Regents of the University of California and the University of Southern California
% Created by Anand A. Joshi, Chitresh Bhushan, David W. Shattuck, Richard M. Leahy 
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; version 2.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
% USA.


function [fileFlag, outFile] = check_nifti_file(nifti_file, file_workdir)
% Check input files for following: 
%   - If it exists on disk
%   - does it has the right extension: .nii/.nii.gz/.img/.hdr
%   - Defined header transforms
%      - when only qform present, then add sform (if file_workdir is defined)
%
%  fileFlag = 1 when valid sform matrix 
%           = 2 when only qform matrix (may NOT be checked for validity)
%           = 3 when both qform and sform matrix exist, & sform matrix is valid
%
% outFile is either empty or a filename with vaild sform matrix.
%

outFile = nifti_file;

if ~ischar(nifti_file)
   error('BDP:CheckNiftiFile:WrongInput', 'nifti_file must be character input, was given \n%s\n', escape_filename(nifti_file));
end

% filename tests
if ~exist(nifti_file, 'file')
   error('BDP:FileDoesNotExist','NIfTI file does not exist: %s', escape_filename(nifti_file));
   
elseif length(nifti_file)<4
   error('BDP:CheckNiftiFile:WrongFileName', ['NIfTI file must end in .nii/.nii.gz/.img/.hdr \n' ...
      'Input file %s is too short for that.'], escape_filename(nifti_file))
   
elseif length(nifti_file)<7 && ~strcmp(nifti_file(end-3:end), '.nii') && ~strcmp(nifti_file(end-3:end), '.img') ...
      && ~strcmp(nifti_file(end-3:end), '.hdr')
   error('BDP:CheckNiftiFile:WrongFileName', ['NIfTI file must end in .nii/.nii.gz/.img/.hdr \n' ...
      'Input file does not: %s'], escape_filename(nifti_file))

elseif ~strcmp(nifti_file(end-6:end), '.nii.gz') && ~strcmp(nifti_file(end-3:end), '.nii') ...
      && ~strcmp(nifti_file(end-3:end), '.img') && ~strcmp(nifti_file(end-3:end), '.hdr')
   error('BDP:CheckNiftiFile:WrongFileName', ['NIfTI file must end in .nii/.nii.gz/.img/.hdr \n' ...
      'Input file does not: %s'], escape_filename(nifti_file))
else
   
   nii = load_untouch_nii_gz(nifti_file);
end


if (nii.hdr.hist.qform_code <= 0) && (nii.hdr.hist.sform_code <= 0)
   %error('BDP:CheckNiftiFile:HeaderError', 'None of sform_code or qform_code are set in the header of file: %s', ...
    nii=add_sform(nii, 1);  
    fileFlag=1;   
%escape_filename(nifti_file));
   
elseif (nii.hdr.hist.qform_code <= 0) && (nii.hdr.hist.sform_code > 0)
   fileFlag = 1;
   
elseif (nii.hdr.hist.qform_code > 0) && (nii.hdr.hist.sform_code <= 0)
   if exist('file_workdir', 'var')
      [~, name, ~] = fileparts(nifti_file);
      outFile = [file_workdir filesep remove_extension(name) '.sform_added.nii.gz'];
      add_sform(nifti_file, outFile);
   else
      fprintf('\nOnly qform matrix is defined in file %s. You can run add_sform to add sform matrix to the nifti file.', nifti_file);
      outFile = [];
   end
   fileFlag = 2;

else % both qform & sform +ve
   fileFlag = 3;
   
end


% sform checks
if nii.hdr.hist.sform_code > 0
   sform(1,:) = nii.hdr.hist.srow_x;
   sform(2,:) = nii.hdr.hist.srow_y;
   sform(3,:) = nii.hdr.hist.srow_z;
   
   Mrest = sform(1:3, 1:3)*diag(1./abs(nii.hdr.dime.pixdim(2:4)));

   if isequal(round(sform), zeros(3,4))
      warning('BDP:CheckNiftiFile:HeaderError', 'sform_code is set but sform matrix is a zero in file: %s', escape_filename(nifti_file));
      
   elseif abs(abs(det(Mrest))-1)>1e-2  % check if Mrest is rotation matrix
      sform
      error('BDP:CheckNiftiFile:HeaderError', ['Absolute determinant of sform matrix is not 1. Header seems to be' ...
         'wrong (not a rotation matrix) in file: %s'], escape_filename(nifti_file))
   end
end

end


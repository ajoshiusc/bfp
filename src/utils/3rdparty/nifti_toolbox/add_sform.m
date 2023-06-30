function [outnii, outFile] = add_sform(dataIn, outFile, force_flag)
% Adds sform parameters using the exisiting qform parameters, only if sform is unset or zero
% matrix in the standard NIFTI file. Saves the data in .gz file. It does not alter image matrix
% or the qform parameters.
%
%   dataFile - .nii or .nii.gz filename with full path
%   outFile -  output file name.
%   force_flag - [optional] when true, tries to write sform parameters even if sform is set
%                and/or NOT zero matrix. 
%
% Similar to public domain nifti1_io.c in implementation - 
%    http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1_io.c
%

if nargin==1
   force_flag = false;
elseif nargin==2 && ischar(outFile)
   force_flag = false;
elseif nargin==2 && ~ischar(outFile)
   force_flag = outFile;
   clear outFile
end

if ischar(dataIn)
   nii = load_untouch_nii_gz(dataIn);
elseif ~isfield(dataIn,'untouch') || dataIn.untouch ~= 1
   error('BDP:addQform:incorrectNiftiLoading','Please use ''load_untouch_nii.m'' to load the structure.');
else
   nii = dataIn;
   clear dataIn
end

sform(1,:) = nii.hdr.hist.srow_x;
sform(2,:) = nii.hdr.hist.srow_y;
sform(3,:) = nii.hdr.hist.srow_z;

if ~force_flag && (nii.hdr.hist.qform_code == 0)
   outnii = nii;
   fprintf('\n qform is NOT set in header. Terminating without any changes.\n');
   
elseif ~force_flag && (nii.hdr.hist.sform_code > 0) && ~isequal(sform, zeros(3,4))
   outnii = nii;
   fprintf('\n sform is already set & is not all zeros. Terminating without any changes.\n');
   
else
   if force_flag
      disp('Force mode. Output may be completely wrong. Make sure of what you are doing.')
   end
   
   % get qform parameters
   b = nii.hdr.hist.quatern_b;
   c = nii.hdr.hist.quatern_c;
   d = nii.hdr.hist.quatern_d;
   a = sqrt(1-(b*b + c*c + d*d));
   
   if ~isreal(a) && abs(a)>1e-2
      error('quatern_b^2 + quatern_c^2 + quatern_d^2 must be <= 1. Something is wrong with the qform header.')
   else
      a = abs(a);
   end
   
   R = [a*a+b*b-c*c-d*d,  2*b*c-2*a*d,       2*b*d+2*a*c;
      2*b*c+2*a*d,      a*a+c*c-b*b-d*d,   2*c*d-2*a*b;
      2*b*d-2*a*c,      2*c*d+2*a*b,       a*a+d*d-c*c-b*b;];
   
   
   if abs(det(R)-1)>1e-2 && force_flag
      warning('matrix R is not proper, ie determinant is not 1. Will generate sform parameters anyway..');
      
   elseif abs(det(R)-1)>1e-2 && ~force_flag
      error('matrix R is not proper, ie determinant is not 1. Check the header. det(R) = %f', det(R));
   end
   
   
   if nii.hdr.dime.pixdim(1)==0
      qfac = 1;
      warning('qfac or pixdim(1) is set to 0 in header. Will assume qfac=1 according to NIfTI-1 guidelines.');
   elseif nii.hdr.dime.pixdim(1) < 0
      qfac = -1;
   else
      qfac = 1;
   end
   
   % modify sform parameters
   outnii = nii; 
   nii.img = [];
   outnii.hdr.hist.sform_code = 1; % NIFTI_XFORM_SCANNER_ANAT
   
   R(:,1) = R(:,1)*nii.hdr.dime.pixdim(2);
   R(:,2) = R(:,2)*nii.hdr.dime.pixdim(3);
   R(:,3) = R(:,3)*nii.hdr.dime.pixdim(4)*qfac;
   
   outnii.hdr.hist.srow_x(1:3) = R(1,:);
   outnii.hdr.hist.srow_y(1:3) = R(2,:);
   outnii.hdr.hist.srow_z(1:3) = R(3,:);
   
   outnii.hdr.hist.srow_x(4) = nii.hdr.hist.qoffset_x;
   outnii.hdr.hist.srow_y(4) = nii.hdr.hist.qoffset_y;
   outnii.hdr.hist.srow_z(4) = nii.hdr.hist.qoffset_z;
   
   if exist('outFile', 'var')
      outFile = save_untouch_nii_gz(outnii, outFile);
   end
end

end


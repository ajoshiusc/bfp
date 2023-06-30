function [outnii, outFile] = add_qform(dataIn, outFile, force_flag)
% Adds qform (Quaternion) parameters to nifti header. The parameters are defined to match sform
% matrix. This function is intended to be used with NIFTI1 files/structs with NIFTI1_XFORM_CODE set
% to NIFTI_XFORM_SCANNER_ANAT (sform_code = 1; qform_code = 1). This implies that sform affine
% matrix should not have any shear component. This function does not alter image matrix or sform
% parameters. 
% Similar to public domain nifti1_io.c in implementation - 
%    http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1_io.c
%
% Usage: 
%   outnii = add_qform(dataIn, force_flag)
%   [outnii, outFile] = add_qform(dataIn, outFile)
%   [outnii, outFile] = add_qform(dataIn, outFile, force_flag)
%   add_qform(dataIn, ...)
%
%   dataIn - Nifti filename with full path or nifti matlab struct
%   outFile - Output filename
%   force_flag - When true, tries to set qform parameters even if qform is set (and skips many
%                sanity checks)
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

if (nii.hdr.hist.qform_code == 0 && nii.hdr.hist.sform_code == 0) && ~force_flag
   outnii = nii;
   warning('BDP:addQform:AnalyzeNii', 'Analyze format found. Nothing to do.');
   
elseif nii.hdr.hist.sform_code==0 && ~force_flag
   outnii = nii;
   warning('BDP:addQform:NoSform', 'sform not found. Nothing to do.');
   
elseif nii.hdr.hist.qform_code>0 && ~force_flag
   outnii = nii;
   warning('BDP:addQform:QformSet', 'qform is already set. Nothing to do.');
   
else
   outnii = nii;
   nii.img = [];
   
   % read sform matrix
   R = zeros(3);
   R(1,1:3) = nii.hdr.hist.srow_x(1:3);
   R(2,1:3) = nii.hdr.hist.srow_y(1:3);
   R(3,1:3) = nii.hdr.hist.srow_z(1:3);
   
   % Compute pixdim
   pixscl = sqrt(sum(R.^2, 1));
   if any(pixscl==0)
      fprintf('Some columns of sform seems to be wrong. Assuming unit scaling for those dimensions.\n');
      temp = [0 0 0];
      temp(pixscl==0) = 1;
      R = R + diag(temp);
      pixscl = pixscl + temp;
   end
   
   % remove scaling from matrix
   R = R./pixscl([1 1 1], :);
   
   % sanity checks on R
   if (abs(det(R))-1)>1e-2
      R
      error('BDP:addQform','Matrix R is not proper, i.e. determinant is not 1. Check the header. det(R) = %f', det(R));
   end
   
   if abs(det(R)-1)<1e-2 % det(R)==1
      outnii.hdr.dime.pixdim(1) = 1;
      
   else
      outnii.hdr.dime.pixdim(1) = -1;
      R(:,3) = -1 * R(:,3);
   end
   
   a = R(1,1) + R(2,2) + R(3,3) + 1;
   
   if a>0.5
      a = 0.5 * sqrt(a);
      b = 0.25 * (R(3,2)-R(2,3))/a;
      c = 0.25 * (R(1,3)-R(3,1))/a;
      d = 0.25 * (R(2,1)-R(1,2))/a;
      
   else   % trickier cases
      xd = 1 + R(1,1) - (R(2,2)+R(3,3));  % 4*b*b
      yd = 1 + R(2,2) - (R(1,1)+R(3,3));  % 4*c*c
      zd = 1 + R(3,3) - (R(1,1)+R(2,2));  % 4*d*d
      
      if xd>1
         b = 0.5 * sqrt(xd) ;
         c = 0.25 * (R(1,2)+R(2,1)) /b;
         d = 0.25 * (R(1,3)+R(3,1)) /b;
         a = 0.25 * (R(3,2)-R(2,3)) /b;
      elseif yd>1
         c = 0.5 * sqrt(yd) ;
         b = 0.25 * (R(1,2)+R(2,1)) /c;
         d = 0.25 * (R(2,3)+R(3,2)) /c;
         a = 0.25 * (R(1,3)-R(3,1)) /c;
      else
         d = 0.5 * sqrt(zd) ;
         b = 0.25 * (R(1,3)+R(3,1)) /d;
         c = 0.25 * (R(2,3)+R(3,2)) /d;
         a = 0.25 * (R(2,1)-R(1,2)) /d;
      end
      
      if a<0
         b=-b;
         c=-c;
         d=-d;
         a=-a;
      end
   end
   
   outnii.hdr.hist.qoffset_x = nii.hdr.hist.srow_x(4);
   outnii.hdr.hist.qoffset_y = nii.hdr.hist.srow_y(4);
   outnii.hdr.hist.qoffset_z = nii.hdr.hist.srow_z(4);
   
   outnii.hdr.hist.quatern_b = b;
   outnii.hdr.hist.quatern_c = c;
   outnii.hdr.hist.quatern_d = d;
   outnii.hdr.hist.qform_code = 1; % NIFTI_XFORM_SCANNER_ANAT
   
   if isequal(outnii.hdr.dime.pixdim(2:4), [0 0 0])
      outnii.hdr.dime.pixdim(2:4) = pixscl;
   elseif norm(outnii.hdr.dime.pixdim(2:4)-pixscl)>0.1
      outnii.hdr.dime.pixdim(2:4) = pixscl;
      
      nii.hdr.dime.pixdim(2:4)
      pixscl
      warning('BDP:addQform:InconsistentPixdim', bdp_linewrap(['Preset pixdim was not consistent with' ...
         'computed pixdim (pixscl). pixdim is overwritten by pixscl.']));
   end
   
end

if exist('outFile', 'var')
   outFile = save_untouch_nii_gz(outnii, outFile);
end

end


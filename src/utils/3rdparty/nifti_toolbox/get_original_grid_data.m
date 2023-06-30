function [data, X_vol, Y_vol, Z_vol, res, Tvol, vol] = get_original_grid_data(vol)
% Returns the actual sampling grid points (in real-world coordinate) and the data at those grid
% points. This function only uses information from sform matrix of nifti header. qform matrix
% is ignored, even if present. No reslicing or interpolation.
%
% Returned data points may be NOT be of 'ndgrid' form. For ndgrid type data points see
%  get_grid_data.m
%

min_res = 0.001; % minimum expected resolution, in mm, along any dim
zero_thresh = 0.05; % max abs value which is considered as good as zero

if ischar(vol)
   vol = load_untouch_nii_gz(vol, true);
elseif ~isfield(vol,'untouch') || vol.untouch ~= 1
   error('Please use ''load_untouch_nii.m'' to load the structure.');
end

if vol.hdr.hist.sform_code<=0
   error(['sform_code is not set in file %s. This file is not supported. You can try using '...
      'add_sform.m first.'], [fileBaseName(vol.fileprefix) '.nii']);
end

Tvol = eye(4);
Tvol(1,:)= vol.hdr.hist.srow_x(1:4);
Tvol(2,:)= vol.hdr.hist.srow_y(1:4);
Tvol(3,:)= vol.hdr.hist.srow_z(1:4);

% Check if Tvol is a ~zero matrix
if sum(abs(vect(Tvol(1:3,1:3))) < zero_thresh) == 9 
   error(['sform seems to be approximately a zero matrix in file %s. This is not supported - please recheck the header. You can try using '...
      'add_sform.m first, if the file has valid q-form matrix.'], [fileBaseName(vol.fileprefix) '.nii']);
end


% Check if even one of the resoultion is ~0
res = abs(vol.hdr.dime.pixdim(2:4)); % resolution. Some software uses negative pixdims to represent a spatial flip
if sum(res < min_res) > 0
   error(['The resolution information in pixdim seems to be close to zero (atleast in one of the dimension) in file %s. This is not unlikely and not supported - please recheck the header. You can try using '...
      'add_sform.m first, if the file has valid q-form matrix.'], [fileBaseName(vol.fileprefix) '.nii']);  
end

Mrest = Tvol(1:3, 1:3)*diag(1./res); % Get only rotation part of the matrix (remove effect of resolution); Mrest can be proper/improper rotation matrix
eval = eig(Mrest);
if sum(abs(eval(:)) < zero_thresh) > 0  % eval can be -ve or complex, but ignoring sign they should be roots of unity if Mrest is rotation matrix 
   error(['The resolution information in sform seems to be close to zero (atleast in one of the dimension) in file %s. This is not supported- please recheck the header. You can try using '...
      'add_sform.m first, if the file has valid q-form matrix.'], [fileBaseName(vol.fileprefix) '.nii']);    
end;

% Skip this check - To allow sform to be an arbitary affine matrix as per
% official nifti-1 guidelines (when sform_code>1). See: https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
% if abs(abs(det(Tvol)) - prod(res)) > zero_thresh % determinant can be -ve (eg. improper rotation matrix)
%    error(['Either the resolution information in sform does not match the header pixdim field in file %s OR the sform matrix has large shear/skew component. Please recheck the header. You can try using '...
%       'add_sform.m first, if the file has valid q-form matrix.'], [fileBaseName(vol.fileprefix) '.nii']);    
% end

data = vol.img;
vol_size = size(data);

% Voxel indexing starts from 0
[X_vol, Y_vol, Z_vol] = ndgrid(0:(vol_size(1)-1), 0:(vol_size(2)-1),  0:(vol_size(3)-1));

c(1,:) = X_vol(:); clear X_vol;
c(2,:) = Y_vol(:); clear Y_vol;
c(3,:) = Z_vol(:);
c(4,:) = 1;
c = Tvol*c;
X_vol = reshape(c(1,:), size(Z_vol));
Y_vol = reshape(c(2,:), size(Z_vol));
Z_vol = reshape(c(3,:), size(Z_vol));
clear c


end

function [data, X_vol, Y_vol, Z_vol, res, Tvol, vol] = get_original_grid_data(vol)
% Returns the actual sampling grid points (in real-world coordinate) and the data at those grid
% points. This function only uses information from sform matrix of nifti header. qform matrix
% is ignored, even if present. No reslicing or interpolation.
%
% Returned data points may be NOT be of 'ndgrid' form. For ndgrid type data points see
%  get_grid_data.m
%

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

if isequal(round(Tvol(1:3,1:3)), zeros(3,3))
   error(['sform seems to a zero matrix in file %s. This is not supported. You can try using '...
      'add_sform.m first, if the file has valid q-form matrix.'], [fileBaseName(vol.fileprefix) '.nii']);
end

res = abs(vol.hdr.dime.pixdim(2:4)); % resolution. Some software uses negative pixdims to represent a spatial flip
data = double(vol.img); % convert to double

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

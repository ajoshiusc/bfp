function [vol_reslice, nii] = myreslice_nii(vol, method, X_grid, Y_grid, Z_grid, fileOut, gz)
% Reslice the nifti image (VOL) to a grid in world-coordinates i.e. applies the sform transformation
% on the image matrix. When [X_grid, Y_grid, Z_grid] is not defined, the image is resliced on a grid
% with same resolution as input image. When [X_grid, Y_grid, Z_grid] are defined then VOL is
% resliced on the specified grid. No successive interpolation is used. Note that resolution is NOT
% matched when [X_grid, Y_grid, Z_grid] are defined. 
%
%  Usage: 
%     vol_reslice = myreslice_nii(vol, method)
%     vol_reslice = myreslice_nii(vol, method, fileOut)
%     vol_reslice = myreslice_nii(vol, method, fileOut, gz)
%     vol_reslice = myreslice_nii(vol, method, X_grid, Y_grid, Z_grid)
%     vol_reslice = myreslice_nii(vol, method, X_grid, Y_grid, Z_grid, fileOut)
%     vol_reslice = myreslice_nii(vol, method, X_grid, Y_grid, Z_grid, fileOut, gz)
%     [vol_reslice, nii] = myreslice_nii(...)
%
% Inputs: 
%    vol    - nifti file name or nii stucture
%    method - interpolation type; set to 'isotropic' for isotropic output with cubic 
%             interpolation (skip X_grid, Y_grid, Z_grid) 
%    fileOut- Output nifti filename
%    X_grid, Y_grid, Z_grid - ndgrid points for interpolation. If this grid point is NOT like
%                             ndgrid, then it will generate INCORRECT nifti files, BUT correct
%                             vol_reslice.
%    gz - boolean flag. When set to false it does not zip the output file. 
%
% Also see:
%   interp3_nii.m - for interpolation from one nifti space to another.
%   myreslice_nii_match_res.m - for matching resolution with [X_grid, Y_grid, Z_grid]
%   reslice_nii.m - in nifti toolbox for better nifti checks and qform support. 
%

[data_orig, X_vol, Y_vol, Z_vol, res, Tvol, vol] = get_original_grid_data(vol);

if strcmpi(method, 'isotropic')
   method = 'cubic';
   res = min(res)*[1 1 1];
end

if nargin<=4
   if exist('X_grid', 'var')
      fileOut = X_grid; 
   end
   if exist('Y_grid', 'var')
      gz = Y_grid; 
   end
   [X_grid, Y_grid, Z_grid] = ndgrid(min(X_vol(:)):res(1):max(X_vol(:)), ...
                                  min(Y_vol(:)):res(2):max(Y_vol(:)), ...
                                  min(Z_vol(:)):res(3):max(Z_vol(:)));
else
   res = [X_grid(2,1,1)-X_grid(1,1,1)  Y_grid(1,2,1)-Y_grid(1,1,1)  Z_grid(1,1,2)-Z_grid(1,1,1)];
end

origin = -1*([X_grid(1,1,1) Y_grid(1,1,1) Z_grid(1,1,1)]./res) + 1;

% get location in original grid of vol 
c(1,:) = X_grid(:); clear X_grid;
c(2,:) = Y_grid(:); clear Y_grid;
c(3,:) = Z_grid(:);
c(4,:) = 1;
c = (inv(Tvol))*c;
X_grid = reshape(c(1,:), size(Z_grid));
Y_grid = reshape(c(2,:), size(Z_grid));
Z_grid = reshape(c(3,:), size(Z_grid));
clear c

% Voxel indexing starts from 0
vol_size = size(data_orig);
[X_vol, Y_vol, Z_vol] = ndgrid(0:vol_size(1)-1, 0:vol_size(2)-1,  0:vol_size(3)-1);


if ndims(data_orig) == 3
   vol_reslice = interpnNNboundary(X_vol, Y_vol, Z_vol, double(data_orig), X_grid, Y_grid, Z_grid, method, 0);
   
elseif ndims(data_orig) == 4
   vol_reslice = zeros([size(X_grid) size(data_orig, 4)]);
   for k = 1:size(data_orig, 4);
      vol_reslice(:,:,:,k) = interpnNNboundary(X_vol, Y_vol, Z_vol, double(data_orig(:,:,:,k)), X_grid, Y_grid, Z_grid, method, 0);
   end
   
else
   error('vol has unsupported dimensions.');
end

vol_reslice(~isfinite(vol_reslice))=0;


if exist('fileOut', 'var') || nargout == 2
   if strcmpi(method, 'nearest') % same data type can support the resliced volume
      nii = make_nii(vol_reslice, res, origin, vol.hdr.dime.datatype, 'resliced volume');
   else      
      nii = make_nii(vol_reslice, res, origin, 64, 'resliced volume');
   end
end

if exist('fileOut', 'var') 
   if exist('gz', 'var') && ~gz
      save_nii_wrapper(nii, fileOut);
   else
      save_nii_gz(nii, fileOut);
   end
end
end


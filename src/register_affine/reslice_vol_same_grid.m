function [vol1_new, vol2_new, X_grid, Y_grid, Z_grid, res, G_new] = reslice_vol_same_grid(v1, v2, res, method, outFile1, outFile2, gz)
% Reslice/Resamples v1 and v2 on the same grid point. (v1, v2 should be registered together)
%  v1, v2 - nii structure or nifti file
%  res - resolution in mm of target output image
%  outFile1, outFile2 - (optional) output file name
%  gz - (optional) set to false to disable compressesd output.
%
%  vol1_new, vol2_new - resliced volume
%
% This function should give good results (No successive interpolation) with all nifti files -
% even with oblique nifti-headers. Also, this function does NOT match the resolution of the images
% before re-sampling!
%

% Blank pixels around the image
padVox = 3;

if length(res)==1
   res = res*[1 1 1];
end

[vol1, X_vol1, Y_vol1, Z_vol1, res1, T1] = get_original_grid_data(v1);
[vol2, X_vol2, Y_vol2, Z_vol2, res2, T2] = get_original_grid_data(v2);

% Find the bounding box
[x_grid, y_grid, z_grid] = ndgrid( min([X_vol1(:); X_vol2(:)]) : max([X_vol1(:); X_vol2(:)]), ...
                                   min([Y_vol1(:); Y_vol2(:)]) : max([Y_vol1(:); Y_vol2(:)]), ...
                                   min([Z_vol1(:); Z_vol2(:)]) : max([Z_vol1(:); Z_vol2(:)]));
clear X_vol1 Y_vol1 Z_vol1 X_vol2 Y_vol2 Z_vol2

% vol 1
[x_vox, y_vox, z_vox, x_orig, y_orig, z_orig] = get_voxel_grid(x_grid, y_grid, z_grid, T1, size(vol1));
tempv = interpnNNboundary(x_orig, y_orig, z_orig, squeeze(vol1(:,:,:,1)), x_vox, y_vox, z_vox, 'nearest', 0);
tempv(isnan(tempv))=0;
clear x_vox y_vox z_vox x_orig y_orig z_orig

[bb_1, bb_2, bb_3] = find_bounding_box(tempv>0);
x = [x_grid(bb_1(1),bb_2(1),bb_3(1)) x_grid(bb_1(end),bb_2(end),bb_3(end))];
y = [y_grid(bb_1(1),bb_2(1),bb_3(1)) y_grid(bb_1(end),bb_2(end),bb_3(end))];
z = [z_grid(bb_1(1),bb_2(1),bb_3(1)) z_grid(bb_1(end),bb_2(end),bb_3(end))];
clear tempv s bb_1 bb_2 bb_3


% vol 2
[x_vox, y_vox, z_vox, x_orig, y_orig, z_orig] = get_voxel_grid(x_grid, y_grid, z_grid, T2, size(vol2));
tempv = interpn(x_orig, y_orig, z_orig, squeeze(vol2(:,:,:,1)), x_vox, y_vox, z_vox, 'nearest', 0);
tempv(isnan(tempv))=0;
clear x_vox y_vox z_vox x_orig y_orig z_orig

[bb_1, bb_2, bb_3] = find_bounding_box(tempv>0);
x = [x x_grid(bb_1(1),bb_2(1),bb_3(1)) x_grid(bb_1(end),bb_2(end),bb_3(end))];
y = [y y_grid(bb_1(1),bb_2(1),bb_3(1)) y_grid(bb_1(end),bb_2(end),bb_3(end))];
z = [z z_grid(bb_1(1),bb_2(1),bb_3(1)) z_grid(bb_1(end),bb_2(end),bb_3(end))];
clear tempv s bb_1 bb_2 bb_3

bound_box = [min(x) min(y) min(z) max(x) max(y) max(z)];
clear tempv x_grid y_grid z_grid s s1 s2

% pad few voxels
width = bound_box(4:6)-bound_box(1:3);
padd = padVox * res;
bound_box(1:3) = bound_box(1:3) - padd;
bound_box(4:6) = bound_box(4:6) + padd;

G_new{1} = bound_box(1):res(1):bound_box(4);
G_new{2} = bound_box(2):res(2):bound_box(5);
G_new{3} = bound_box(3):res(3):bound_box(6);

[X_grid, Y_grid, Z_grid] = ndgrid(G_new{1}, G_new{2}, G_new{3});
clear X_new Y_new Z_new

[x_vox, y_vox, z_vox, x_orig, y_orig, z_orig] = get_voxel_grid(X_grid, Y_grid, Z_grid, T1, size(vol1));
if ndims(vol1) == 3
   vol1_new = interpnNNboundary(x_orig, y_orig, z_orig, double(vol1), x_vox, y_vox, z_vox, method, 0);
elseif ndims(vol1) == 4
   vol1_new = zeros([size(X_grid) size(vol1, 4)]);
   parfor k = 1:size(vol1, 4);
      vol1_new(:,:,:,k) = interpnNNboundary(x_orig, y_orig, z_orig, double(vol1(:,:,:,k)), x_vox, y_vox, z_vox, method, 0);
   end
else
   error('v1 has unsupported dimensions.');
end
clear vol1

[x_vox, y_vox, z_vox, x_orig, y_orig, z_orig] = get_voxel_grid(X_grid, Y_grid, Z_grid, T2, size(vol2));
if ndims(vol2) == 3
   vol2_new = interpnNNboundary(x_orig, y_orig, z_orig, double(vol2), x_vox, y_vox, z_vox, method, 0);
elseif ndims(vol2) == 4
   vol2_new = zeros([size(X_grid) size(vol2, 4)]);
   parfor k = 1:size(vol2, 4);
      vol2_new(:,:,:,k) = interpnNNboundary(x_orig, y_orig, z_orig, double(vol2(:,:,:,k)), x_vox, y_vox, z_vox, method, 0);
   end
else
   error('v2 has unsupported dimensions.');
end
clear vol2 x_vox y_vox z_vox x_orig y_orig z_orig

vol1_new(isnan(vol1_new))=0;
vol2_new(isnan(vol2_new))=0;
vol1_new(vol1_new<0)=0;
vol2_new(vol2_new<0)=0;

origin = -1*(bound_box(1:3)./res) + 1;

if exist('outFile1','var')
   nii1 = make_nii(vol1_new, res, origin, 64, 'resliced volume');
   if exist('gz', 'var') && ~gz
      save_nii_wrapper(nii1, outFile1);
   else
      save_nii_gz(nii1, outFile1);
   end
end

if exist('outFile2','var')
   nii2 = make_nii(vol2_new, res, origin, 64, 'resliced volume');
   if exist('gz', 'var') && ~gz
      save_nii_wrapper(nii2, outFile2);
   else
      save_nii_gz(nii2, outFile2);
   end
end

end



% get location in original grid of v1
function [x_vox, y_vox, z_vox, x_orig, y_orig, z_orig] = get_voxel_grid(x_grid, y_grid, z_grid, Tvol, vol_size)
c(1,:) = x_grid(:); clear x_grid;
c(2,:) = y_grid(:); clear y_grid;
c(3,:) = z_grid(:);
c(4,:) = 1;
c = (inv(Tvol))*c;
x_vox = reshape(c(1,:), size(z_grid));
y_vox = reshape(c(2,:), size(z_grid));
z_vox = reshape(c(3,:), size(z_grid));
clear c zgrid

[x_orig, y_orig, z_orig] = ndgrid(0:vol_size(1)-1, 0:vol_size(2)-1,  0:vol_size(3)-1);

end

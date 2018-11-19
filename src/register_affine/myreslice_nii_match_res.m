function [vol_reslice, nii] = myreslice_nii_match_res(vol, method, X_grid, Y_grid, Z_grid, fileOut, gz)
% Similar to myreslice_nii.m but it matches the resolution of vol to target resolution before
% interpolation. It is useful in case of downsampling (can not really match resolution while
% upsampling). vol must be 3D volume.
%
% Usage: 
%   vol_reslice = myreslice_nii_match_res(vol, method, X_grid, Y_grid, Z_grid)
%   vol_reslice = myreslice_nii_match_res(vol, method, X_grid, Y_grid, Z_grid, rotate_window)
%   vol_reslice = myreslice_nii_match_res(vol, method, X_grid, Y_grid, Z_grid, fileOut)
%   vol_reslice = myreslice_nii_match_res(vol, method, X_grid, Y_grid, Z_grid, fileOut, gz)
%   [vol_reslice, nii] = myreslice_nii_match_res(...)
% 
% Inputs: 
%   vol - nifti file name or nii stucture
%   method - interpolation type; 
%   X_grid, Y_grid, Z_grid - ndgrid points for interpolation. If this grid point is NOT like ndgrid,
%                            then it will generate INCORRECT nifti files, BUT correct vol_reslice. 
%   fileOut - Output nifti filename
%
% Also see myreslice_nii.m, interp3_nii.m
%


[data_orig, X_vol, Y_vol, Z_vol, ~, Tvol, vol] = get_original_grid_data(vol);

if ndims(data_orig)>3
   error('vol must be a 3D volume');
end

if nargin<5
   error('Atleast 5 inputs are required.');
elseif nargin==6 && islogical(fileOut)
   rotate_window = fileOut;
   clear fileOut
else
   rotate_window = false;
end

[res_vox_target, xyz_ind_target] = findResolution(X_grid, Y_grid, Z_grid);
[res_vox_vol, xyz_ind_vol] = findResolution(X_vol, Y_vol, Z_vol);
origin = -1*([X_grid(1,1,1) Y_grid(1,1,1) Z_grid(1,1,1)]./res_vox_target) + 1;

% match resolution by convolving
res_target = res_vox_target(xyz_ind_target);
res_vol = res_vox_vol(xyz_ind_vol);
res_ratio = res_target./res_vol;

if max(res_ratio)<1
   data_res_matched = data_orig;
   
elseif rotate_window
   % this works accurately only if vol has isotropic resolution
   invInd(xyz_ind_target) = 1:numel(xyz_ind_target);
   res_ratio = res_ratio(invInd); % in target-voxel coordinate
   
   R_grid = findRotGrid(X_grid, Y_grid, Z_grid, Tvol, res_vox_vol, res_vox_target);
   W = rotatedBluringKernal(res_ratio, R_grid);   
   data_res_matched = convnPadded(data_orig, W);
   
else
   % ratio in voxel dimension of data/vol
   invInd(xyz_ind_vol) = 1:numel(xyz_ind_vol);
   res_ratio = res_ratio(invInd); % in vol-voxel coordinate
   
   w1 = gaussian1Dwindow(res_ratio(1));
   w2 = gaussian1Dwindow(res_ratio(2));
   w3 = gaussian1Dwindow(res_ratio(3));
   
   [w1, w2, w3] = ndgrid(w1, w2, w3);
   W = w1 .* w2 .* w3;
   W = W./sum(W(:));
   data_res_matched = convnPadded(data_orig, W);
   
   clear w1 w2 w3 W Y
end
clear data_orig

% interpolate
vol.img = data_res_matched;
vol_reslice = myreslice_nii(vol, method, X_grid, Y_grid, Z_grid);


if exist('fileOut', 'var') || nargout == 2
   if strcmpi(method, 'nearest') % same data type can support the resliced volume
      nii = make_nii(vol_reslice, res_vox_target, origin, vol.hdr.dime.datatype, 'resliced volume');
   else      
      nii = make_nii(vol_reslice, res_vox_target, origin, 64, 'resliced volume');
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

function w = createApodizationWindow(res_ratio, dim_siz)
w = ones(1,dim_siz);

if res_ratio>1 && res_ratio<2
   len_out = floor(dim_siz*(res_ratio-1)/2)*2; % always even size
   len1 = len_out/2;
   h = hamming(len1*2, 'symmetric'); 
   w(1:len1) = h(1:len1);
   w(end-len1+1:end) = h(len1+1:end);
   
elseif res_ratio==2
   w = hamming(dim_siz, 'symmetric');
   
elseif res_ratio>2
   w = zeros(1,dim_siz);
   len_out = floor(dim_siz*(1-(2/res_ratio))/2)*2; % always even size
   len1 = dim_siz-len_out;
   h = hamming(len1, 'symmetric');   
   w(1+(len_out/2):len1+(len_out/2)) = h;
end
end

function w = gaussian1Dwindow(res_ratio)
% res_ration is FWHM

sgm = res_ratio/sqrt(8*log(2));
siz = sgm*6;
x = -ceil(siz/2):ceil(siz/2);
w = exp(-(x.^2/(2*(sgm^2))));
w = w/sum(w(:));

end

function W = rotatedBluringKernal(res_ratio, R_grid)
sgm = res_ratio/sqrt(8*log(2));
siz = sgm*6;
s = ceil(max(siz)/2);
siz = 2*s + 1;

[x, y, z] = ndgrid(-s:s);

Xvol = transpose(inv(R_grid)*([x(:) y(:) z(:)]'));

% Compute: W = mvnpdf(Xvol, [0 0 0], diag(sgm));
sqrt_invsig = transpose(sqrt(1./sgm(:)));
temp = Xvol.*sqrt_invsig(ones(1,size(Xvol,1)), :);
W = exp(-0.5 *sum(temp.^2, 2));
W = reshape(W/sum(W(:)), [siz siz siz]);

end


function [res_vox, xyz_ind] = findResolution(Xg, Yg, Zg)
% res_vox : resolution along different dimension of grid points. 
% xyz_ind : Dimension along which x,y,z are increasing in that order. So the resolution in order of
%           x followed by y and z is res_vox(xyz_ind)
%

% This will fail in 45 degree rotation matrix!
xdiff = abs([Xg(end,1,1) Xg(1,end,1) Xg(1,1,end)] - Xg(1,1,1));
ydiff = abs([Yg(end,1,1) Yg(1,end,1) Yg(1,1,end)] - Yg(1,1,1));
zdiff = abs([Zg(end,1,1) Zg(1,end,1) Zg(1,1,end)] - Zg(1,1,1));

xind = find(xdiff==max(xdiff), 1, 'first');
yind = find(ydiff==max(ydiff), 1, 'first');
zind = find(zdiff==max(zdiff), 1, 'first');

xyz_ind = [xind yind zind];
if length(unique(xyz_ind))~=3 % 45 degree rotation matrix!
   ind_set = [1 2 3];
   [~, IX] = sort(xyz_ind);
   invIX(IX) = 1:numel(IX); % Inverse sort index
   xyz_ind = ind_set(invIX);
end

V0 = [Xg(1,1,1); Yg(1,1,1); Zg(1,1,1)];
V1 = [Xg(2,1,1); Yg(2,1,1); Zg(2,1,1)];
V2 = [Xg(1,2,1); Yg(1,2,1); Zg(1,2,1)];
V3 = [Xg(1,1,2); Yg(1,1,2); Zg(1,1,2)];

res_vox = [norm(V1-V0) norm(V2-V0) norm(V3-V0)];

end


function R_grid = findRotGrid(Xg, Yg, Zg, Tvol, res_vox_vol, res_vox_target)

vol_size = size(Xg);
[X_vol, Y_vol, Z_vol] = ndgrid(0:(vol_size(1)-1), 0:(vol_size(2)-1),  0:(vol_size(3)-1));
Y = [X_vol(:) Y_vol(:) Z_vol(:) ones(numel(Xg), 1)];
X = [Xg(:) Yg(:) Zg(:) ones(numel(Xg), 1)];
T_tar = transpose(Y\X);
clear X Y

R_tar = T_tar*inv(diag([res_vox_target(:); 1]));
R_vol = Tvol*inv(diag([res_vox_vol(:); 1]));

R_tar = R_tar(1:3, 1:3);
R_vol = R_vol(1:3, 1:3);

R_grid = inv(R_vol)*R_tar; % rotates from target-voxel to vol-voxel

end

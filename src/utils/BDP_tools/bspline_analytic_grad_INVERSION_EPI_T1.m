% 
% BDP BrainSuite Diffusion Pipeline
% 
% Copyright (C) 2019 The Regents of the University of California and
% the University of Southern California
% 
% Created by Chitresh Bhushan, Divya Varadarajan, Justin P. Haldar, Anand A. Joshi,
%            David W. Shattuck, and Richard M. Leahy
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
% 


function [X_error, X_error_grad] = bspline_analytic_grad_INVERSION_EPI_T1( ...
                                                    X, Osize1, O_grid, spacing, rigid_scl, ...
                                                    Iepi, Ist, IstINV, res, mask_Iepi, mask_Iepi_less_csf, mask_Ist, opt)                                                 
% Inverts both EPI and T1 (two terms in cost function).
% mask_Iepi_less_csf - Mask using which INVERSION mapping is computed
% 

% Check/set input options
defaultoptions = struct(...
   'penalty', [1 0.00001 0.00001 0.00001],... Required - [alpha kx bigkx ky]
   'intensity_correct', [], ... % This is ignored - distortion is always estimated with correct model. 
   ...                          % This input-option is retained for historical reasons.
   ...                          % Final application of the Jacobian term can be avoided at final application of distortion field.
   'rigid_refine', true, ... 
   'mask_type', 'and',...
   'invEPI_weight', [],... Required - relative weight of error term for inverted-EPI
   'debug', false, ...
   'nthreads', 6);

if ~exist('opt','var')
   opt = defaultoptions;
else
   tags = fieldnames(defaultoptions);
   for i=1:length(tags)
      if(~isfield(opt,tags{i})),  opt.(tags{i})=defaultoptions.(tags{i}); end
   end
end

invEPI_weight = opt.invEPI_weight;

% set parameters
rigid_param = X(1:6)';
O_grid(:,:,:,1) = reshape(X(7:end),Osize1);


% Apply rigid transform to struct image
mode = 1; % linear interpolation and outside pixels set to zero
par = rigid_param.*rigid_scl;
M = par2affineMat(par(1:3),par(4:6));
Ist_rigid = affine_transform_3dvol_double(double(Ist), double(M), double(mode), double(opt.nthreads));
IstINV_rigid = affine_transform_3dvol_double(double(IstINV), double(M), double(mode), double(opt.nthreads));
mask_Ist_rigid = affine_transform_3dvol_double(double(mask_Ist),double(M),double(mode), double(opt.nthreads));


% Transform (& intensity correct) Imv & mask
t_x = bspline_repeated_endpoint_deformation_3d_double_only_x(O_grid(:,:,:,1), size(Iepi), spacing, opt.nthreads);
[x_g, y_g, z_g] = ndgrid(1:size(Iepi,1), 1:size(Iepi,2), 1:size(Iepi,3));
x_g = x_g + t_x;

% opt.intensity_correct
[~, crct_grd] = gradient(t_x);
Iepi_w = interpn(Iepi, x_g, y_g, z_g, 'linear', 0);
Iepi_t = Iepi_w .* (1+crct_grd);

Iepi_t(Iepi_t>1)=1;
Iepi_t(Iepi_t<0)=0;

mask_Iepi_t = interpn(double(mask_Iepi), x_g, y_g, z_g, 'linear', 0);
mask_Iepi_less_csf_t = interpn(double(mask_Iepi_less_csf), x_g, y_g, z_g, 'linear', 0);
clear x_g y_g z_g t_x

% INVERT tranformed EPI 
Iepi_t_inv = opt.histeqInterpolant(1-Iepi_t(:));
Iepi_t_inv = reshape(Iepi_t_inv, size(Iepi_t));
Iepi_t_inv_msked = Iepi_t_inv .* mask_Iepi_t; % mask inverted image


% Calculate the current registration error
if strcmpi(opt.mask_type, 'and')
   mask_err_iEPI = mask_Iepi_t .* mask_Ist_rigid;
   mask_err_iT1 = mask_Iepi_less_csf_t .* mask_Ist_rigid;
else
   mask_err_iEPI = mask_Iepi_t + mask_Ist_rigid - (mask_Iepi_t .* mask_Ist_rigid);
   mask_err_iT1 = mask_Iepi_less_csf_t + mask_Ist_rigid - (mask_Iepi_less_csf_t .* mask_Ist_rigid);
end

diff_img_iEPI = Ist_rigid - Iepi_t_inv_msked;
diff_img_iT1 = IstINV_rigid - Iepi_t;
X_error = ( (invEPI_weight * ((diff_img_iEPI(:).^2)' * mask_err_iEPI(:))) ...
   + ((diff_img_iT1(:).^2)' * mask_err_iT1(:)) ) / numel(diff_img_iEPI);

if opt.debug
   display_volume([Ist_rigid; Iepi_t_inv_msked; IstINV_rigid; Iepi_t; Iepi], [0 1]);
   display_volume(([invEPI_weight*diff_img_iEPI.*mask_err_iEPI; diff_img_iT1.*mask_err_iT1]))
   display_volume([mask_err_iEPI; mask_Ist_rigid; mask_Iepi_t], [0 1]);
   overlay_volumes([(diff_img_iEPI.^2 .*(mask_err_iEPI*10)); Ist_rigid; Iepi_t_inv_msked; Iepi_t], [mask_err_iEPI; mask_err_iEPI; mask_err_iEPI; mask_err_iEPI]/3);
   overlay_volumes([(diff_img_iEPI.^2 .*(mask_err_iEPI*10)); Ist_rigid; Iepi_t_inv; Iepi_t], [mask_err_iEPI; mask_Ist_rigid; mask_Iepi_t; mask_Iepi_t]/3);
end


% gradient
if nargout > 1
   X_error_grad = gradient_analytical_intensity_correct(spacing, rigid_param, rigid_scl, size(O_grid), ...
      Iepi_t, Iepi_w, Iepi_t_inv, Ist_rigid, IstINV_rigid, crct_grd, diff_img_iEPI, diff_img_iT1, ...
      mask_err_iEPI, mask_err_iT1, mask_Ist_rigid, mask_Iepi_t, mask_Iepi_less_csf_t, invEPI_weight, opt);
end

% Add penalty to total error
if opt.penalty(1)>0
   if nargout > 1
      % [SO_error, SO_grad] = bspline_penalty_repeated_endpoint_metal_bending_EPI(O_grid/opt.penalty_op.resize_scl, ...
      %   opt.penalty_op.V, opt.penalty_op.cp_ind, opt.penalty_op.Q);
      
      [SO_error, SO_grad] = bspline_penalty_repeated_endpoint_CP_diff_EPI(O_grid(:,:,:,1), ...
         spacing, res, opt.penalty(2:end));
      
      X_error = X_error + (opt.penalty(1)*SO_error);
      X_error_grad = X_error_grad + (opt.penalty(1)*[zeros(6,1); SO_grad(:)]);
      
   else
      % SO_error = bspline_penalty_repeated_endpoint_metal_bending_EPI(O_grid/opt.penalty_op.resize_scl, ...
      %   opt.penalty_op.V, opt.penalty_op.cp_ind, opt.penalty_op.Q);
      
      SO_error = bspline_penalty_repeated_endpoint_CP_diff_EPI(O_grid(:,:,:,1), spacing, res, opt.penalty(2:end));
      X_error = X_error + (opt.penalty(1)*SO_error);
   end
end


end

function X_grad = gradient_analytical_intensity_correct(spacing, rigid_param, rigid_scl, O_grid_size,...
   Iepi_t, Iepi_w, Iepi_t_inv, Ist_rigid, IstINV_rigid, crct_grd, diff_img_iEPI, diff_img_iT1, ...
   mask_err_iEPI, mask_err_iT1,  mask_Ist_rigid, mask_Iepi_t, mask_Iepi_less_csf_t, invEPI_weight, opt)
% Calculates the analytical gradient of SSD error (with Jacobian intensity correction) w.r.t.
% (bspline ctrl points + rigid params)

g_intsty_map = opt.histeqGradInterpolant(1-Iepi_t(:));
g_intsty_map = reshape(g_intsty_map, size(Iepi_t));

[~, g_Iepi_w] = gradient(Iepi_w);
[~, g_mask_Iepi_t] = gradient(double(mask_Iepi_t));
[~, g_mask_Iepi_less_csf_t] = gradient(double(mask_Iepi_less_csf_t));

temp = -2 *invEPI_weight .* double(mask_err_iEPI) .* diff_img_iEPI;
d_Fmv_t1_iEPI = temp .* (Iepi_t_inv .* g_mask_Iepi_t - double(mask_Iepi_t) .* g_intsty_map .* (1+crct_grd) .* g_Iepi_w);
d_Fmv_t2_iEPI = temp .* (-1 .* double(mask_Iepi_t) .* g_intsty_map .* Iepi_w);

temp = -2 .* double(mask_err_iT1) .* diff_img_iT1;
d_Fmv_t1_iT1 = temp .* (1+crct_grd) .* g_Iepi_w;        
d_Fmv_t2_iT1 = temp .* Iepi_w;

if strcmpi(opt.mask_type, 'and')   
   d_Fmv_t3_iEPI =  invEPI_weight .* (diff_img_iEPI.^2) .* double(mask_Ist_rigid).* g_mask_Iepi_t;
   d_Fmv_t3_iT1 =  (diff_img_iT1.^2) .* double(mask_Ist_rigid).* g_mask_Iepi_less_csf_t;
else
   d_Fmv_t3_iEPI = invEPI_weight .* (diff_img_iEPI.^2) .* (1-double(mask_Ist_rigid)) .* g_mask_Iepi_t;
   d_Fmv_t3_iT1 = (diff_img_iT1.^2) .* (1-double(mask_Ist_rigid)) .* g_mask_Iepi_less_csf_t;
end

d_Fmv = d_Fmv_t1_iEPI + d_Fmv_t3_iEPI + d_Fmv_t1_iT1 + d_Fmv_t3_iT1;
d_Fmv2 = d_Fmv_t2_iEPI + d_Fmv_t2_iT1;
clear d_Fmv_t1_* d_Fmv_t2_* d_Fmv_t3_*

O_grad = zeros(O_grid_size(1:3));
d_P = zeros(numel(mask_err_iEPI)+1, 1); % for fake zero voxel
for zi=1:4
   for zj=1:4
      for zk=1:4
         [xx, yy, zz] = ndgrid(zi:4:O_grid_size(1), zj:4:O_grid_size(2), zk:4:O_grid_size(3));
         CP_ind = sub2ind(O_grid_size(1:3), xx(:), yy(:), zz(:)); % control-points in the set
         
         O_temp = zeros(O_grid_size(1:3));
         O_temp(CP_ind) = 1;
         [bspline_prod, dbspline_prod] = ...
            bspline_repeated_endpoint_deformation_3d_double_only_x_gradient(O_temp, size(mask_err_iEPI), spacing, opt.nthreads);
         dbspline_prod = dbspline_prod./spacing(1);
         
         d_Prod = (d_Fmv .* bspline_prod) + (d_Fmv2 .* dbspline_prod);
         d_P(1:(end-1)) = d_Prod(:);
         
         d_msk = sum(d_P(opt.cp_mask_indx(:, opt.cp_col_index(CP_ind))), 1);
         O_grad(CP_ind) = d_msk(:);
      end
   end
end

% Normalize by number of voxels used for each control point 
O_grad = O_grad ./ prod(spacing*4 + 1);




% analytical Rigid transform gradient
if opt.rigid_refine
   par_current = rigid_param .* rigid_scl;
   
   [gy, gx, gz] = gradient(Ist_rigid);
   temp = 2 * invEPI_weight .* mask_err_iEPI(:) .* diff_img_iEPI(:);
   Gdmsk1_iEPI = temp(:, [1 1 1]) .* [gx(:) gy(:) gz(:)];
   
   [gy, gx, gz] = gradient(IstINV_rigid);
   temp = 2 .* mask_err_iT1(:) .* diff_img_iT1(:);
   Gdmsk1_iT1 = temp(:, [1 1 1]) .* [gx(:) gy(:) gz(:)];
   
   [gy_mask_Ist_rigid, gx_mask_Ist_rigid, gz_mask_Ist_rigid] = gradient(double(mask_Ist_rigid));
   if strcmpi(opt.mask_type, 'and')
      temp = invEPI_weight * (diff_img_iEPI.^2) .* double(mask_Iepi_t);
      temp2 = (diff_img_iT1.^2) .* double(mask_Iepi_less_csf_t);
   else
      temp = invEPI_weight * (diff_img_iEPI.^2) .* double(1-mask_Iepi_t);
      temp2 = (diff_img_iT1.^2) .* double(1-mask_Iepi_less_csf_t);
   end
   temp = temp(:);
   temp2 = temp2(:);
   Gdmsk2_iEPI = temp(:, [1 1 1]) .* [gx_mask_Ist_rigid(:) gy_mask_Ist_rigid(:) gz_mask_Ist_rigid(:)]; % Col has grad along [x y z]
   Gdmsk2_iT1 = temp2(:, [1 1 1]) .* [gx_mask_Ist_rigid(:) gy_mask_Ist_rigid(:) gz_mask_Ist_rigid(:)]; % Col has grad along [x y z]
   
   G_msk = Gdmsk1_iEPI + Gdmsk1_iT1 + Gdmsk2_iEPI + Gdmsk2_iT1;
   clear gx gy gz temp temp2 Gdmsk1_* Gdmsk2_*
   
   
   dM = computedM(par_current);
   % Make center of the image coordinates 0,0
   [xd, yd, zd] = ndgrid(0:size(Iepi_t,1)-1, 0:size(Iepi_t,2)-1, 0:size(Iepi_t,3)-1);
   X_origin = size(Iepi_t)/2;
   S = eye(4);
   S(1:3,end) = -1*X_origin(:);
   Sinv = eye(4);
   Sinv(1:3,end) = X_origin(:);
   
   X = [xd(:) yd(:) zd(:) ones([length(xd(:)) 1])];
   clear xd yd zd
   
   rigid_err_grad = zeros(1,length(rigid_param));
   for par_i=1:6
      % M_curr = squeeze(dM(:,:,par_i));
      M_curr = Sinv * squeeze(dM(:,:,par_i)) * S;
      d_TX = X * transpose(M_curr); % transpose() b/c d_TX is transposed to make next step vectorized
      d_TX(:, 4) = [];    % size of [n_vox 3]
      
      rigid_err_grad(par_i) = transpose(d_TX(:)) * G_msk(:); % dot prod followed by sum over voxels
   end
   
   % Normalize rigid_grad by voxel count
   rigid_err_grad = (rigid_err_grad.*rigid_scl)/numel(mask_Iepi_t);
   
else
   rigid_err_grad = zeros(size(rigid_param));
end


% final grad
X_grad = [rigid_err_grad(:);  O_grad(:)];

end


function dM = computedM(par_current)
rot = par_current(4:6)*(pi/180); % convert to radians

Mt=[1 0 0 par_current(1);
   0 1 0 par_current(2);
   0 0 1 par_current(3);
   0 0 0 1];

Rx=[1 0 0 0;
   0 cos(rot(1)) -sin(rot(1)) 0;
   0 sin(rot(1)) cos(rot(1)) 0;
   0 0 0 1];

Ry=[cos(rot(2)) 0 sin(rot(2)) 0;
   0 1 0 0;
   -sin(rot(2)) 0 cos(rot(2)) 0;
   0 0 0 1];

Rz=[cos(rot(3)) -sin(rot(3)) 0 0;
   sin(rot(3)) cos(rot(3)) 0 0;
   0 0 1 0;
   0 0 0 1];

dRx=[0 0 0 0;
   0 -sin(rot(1)) -cos(rot(1)) 0;
   0 cos(rot(1)) -sin(rot(1)) 0;
   0 0 0 0];

dRy=[-sin(rot(2)) 0 cos(rot(2)) 0;
   0 0 0 0;
   -cos(rot(2)) 0 -sin(rot(2)) 0;
   0 0 0 0];

dRz=[-sin(rot(3)) -cos(rot(3)) 0 0;
   cos(rot(3)) -sin(rot(3)) 0 0;
   0 0 0 0;
   0 0 0 0];

dM = zeros([4 4 6]);

dM(1,4,1) = 1; % gradient due to Tx
dM(2,4,2) = 1; % gradient due to Ty
dM(3,4,3) = 1; % gradient due to Tz
dM(:,:,4) = Mt*dRx*Ry*Rz*(pi/180); % multiplied by (pi/180) to correct for units
dM(:,:,5) = Mt*Rx*dRy*Rz*(pi/180);
dM(:,:,6) = Mt*Rx*Ry*dRz*(pi/180);
end


function X_grad = gradient_analytic(Spacing, O_grid, rigid_param, rigid_scl, O_grid_size, ...
   res, Iepi_t, Iepi_t_inv, Ist_rigid, diff_img, mask_err, mask_Ist_rigid, mask_Iepi_t, opt)
% I made a decision not to support this function -- This function can potentially have bugs. 
% Retained for historical reasons and future reference.
%
% Calculates the analytical gradient of SSD error (WITHOUT Jacobian intensity correction) w.r.t.
% (bspline ctrl points + rigid params)

fprintf('this function is not tested/check yet for formula and implementation!')

g_intsty_map = opt.histeqGradInterpolant(1-Iepi_t(:));
g_intsty_map = reshape(g_intsty_map, size(Iepi_t));

[~, g_Iepi_t] = gradient(Iepi_t);
[~, g_mask_Iepi_t] = gradient(double(mask_Iepi_t));
d_Fmv_1 = -2.* double(mask_err) .* diff_img.* (Iepi_t_inv.*g_mask_Iepi_t - mask_Iepi_t.*g_intsty_map.*g_Iepi_t) ;

if strcmpi(opt.mask_type, 'and')   
   d_Fmv_2 =  (diff_img.^2) .* double(mask_Ist_rigid).* g_mask_Iepi_t;
else
   d_Fmv_2 = (diff_img.^2) .* (1-double(mask_Ist_rigid)) .* g_mask_Iepi_t;
end
d_Fmv = d_Fmv_1 + d_Fmv_2;


% cubic bspline kernal (from basis function & reparameterized coord)
temp1 = 0:(1/Spacing(1)):1;
temp1 = 0.5*(temp1.^3) - (temp1.^2) + 4/6; % b_2(u)
temp2 = (1-(1/Spacing(1))):(-1/Spacing(1)):0;
temp2 = (temp2.^3)/6; % b_0(u)
Bu = [temp1 temp2];

temp1 = 0:(1/Spacing(2)):1;
temp1 = 0.5*(temp1.^3) - (temp1.^2) + 4/6;
temp2 = (1-(1/Spacing(2))):(-1/Spacing(2)):0;
temp2 = (temp2.^3)/6;
Bv = [temp1 temp2];

temp1 = 0:(1/Spacing(3)):1;
temp1 = 0.5*(temp1.^3) - (temp1.^2) + 4/6;
temp2 = (1-(1/Spacing(3))):(-1/Spacing(3)):0;
temp2 = (temp2.^3)/6;
Bw = [temp1 temp2];

xmin = -2*Spacing(1);
xmax = 2*Spacing(1)-1;

ymin = -2*Spacing(2);
ymax = 2*Spacing(2)-1;

zmin = -2*Spacing(3);
zmax = 2*Spacing(3)-1;

[bx, by, bz] = ndgrid((xmin:xmax), (ymin:ymax), (zmin:zmax));
b_prod = Bu(abs(bx)+1) .* Bv(abs(by)+1) .* Bw(abs(bz)+1);
b_prod = repmat(b_prod, floor((O_grid_size(1:3)-1)/4));

O_grad = zeros(O_grid_size(1:3));
for zi=1:4
   for zj=1:4
      for zk=1:4 
         
         bspline_prod = circshift(b_prod, ([zi zj zk]-3).*Spacing);
         pixs_rem = (mod(O_grid_size(1:3)-1, 4).*Spacing) + 1;
         pixs_rem(pixs_rem<0) = 0;
         bspline_prod(end+1:end+pixs_rem(1), :, :) = bspline_prod(1:pixs_rem(1), :, :);
         bspline_prod(:, end+1:end+pixs_rem(2), :) = bspline_prod(:, 1:pixs_rem(2), :);
         bspline_prod(:, :, end+1:end+pixs_rem(3)) = bspline_prod(:, :, 1:pixs_rem(3));

         d_Prod = d_Fmv .* bspline_prod;
         d_Prod = d_Prod(:);                  
         d_Prod(end+1) = 0; % for fake zero voxel
         
         img_ind = opt.grad_img_indx(opt.grad_mask_index(zi,zj,zk):opt.grad_mask_index(zi,zj,zk) ...
            + (opt.grad_mask_count(zi,zj,zk)*opt.grad_ind_sz)-1);
         d_msk = d_Prod(img_ind);
         d_msk = sum(reshape(d_msk, opt.grad_ind_sz, []), 1);
         temp_msk = false(size(O_grad));
         temp_msk(zi:4:size(O_grid,1), zj:4:size(O_grid,2), zk:4:size(O_grid,3)) = true;
         O_grad(temp_msk) = d_msk(:);         
      end
   end
end

% Normalize by number of voxels used for each control point 
O_grad = O_grad ./ prod(Spacing*4 + 1);



% analytical Rigid transform gradient - Same as gradient_analytical_intensity_correct()
par_current = rigid_param .* rigid_scl;

[gy_Ist_rigid, gx_Ist_rigid, gz_Ist_rigid] = gradient(Ist_rigid);
temp = 2 .* mask_err(:) .* diff_img(:);
Gdmsk1 = temp(:, [1 1 1]) .* [gx_Ist_rigid(:) gy_Ist_rigid(:) gz_Ist_rigid(:)];

[gy_mask_Ist_rigid, gx_mask_Ist_rigid, gz_mask_Ist_rigid] = gradient(double(mask_Ist_rigid));
if strcmpi(opt.mask_type, 'and')
   temp = (diff_img.^2) .* double(mask_Iepi_t);
else
   temp = (diff_img.^2) .* double(1-mask_Iepi_t);
end
temp = temp(:);
Gdmsk2 = temp(:, [1 1 1]) .* [gx_mask_Ist_rigid(:) gy_mask_Ist_rigid(:) gz_mask_Ist_rigid(:)]; % Col has grad along [x y z]

G_msk = Gdmsk1 + Gdmsk2;


rot = par_current(4:6)*(pi/180); % convert to radians

Mt=[1 0 0 par_current(1);
   0 1 0 par_current(2);
   0 0 1 par_current(3);
   0 0 0 1];

Rx=[1 0 0 0;
   0 cos(rot(1)) -sin(rot(1)) 0;
   0 sin(rot(1)) cos(rot(1)) 0;
   0 0 0 1];

Ry=[cos(rot(2)) 0 sin(rot(2)) 0;
   0 1 0 0;
   -sin(rot(2)) 0 cos(rot(2)) 0;
   0 0 0 1];

Rz=[cos(rot(3)) -sin(rot(3)) 0 0;
   sin(rot(3)) cos(rot(3)) 0 0;
   0 0 1 0;
   0 0 0 1];

dRx=[0 0 0 0;
   0 -sin(rot(1)) -cos(rot(1)) 0;
   0 cos(rot(1)) -sin(rot(1)) 0;
   0 0 0 0];

dRy=[-sin(rot(2)) 0 cos(rot(2)) 0;
   0 0 0 0;
   -cos(rot(2)) 0 -sin(rot(2)) 0;
   0 0 0 0];

dRz=[-sin(rot(3)) -cos(rot(3)) 0 0;
   cos(rot(3)) -sin(rot(3)) 0 0;
   0 0 0 0;
   0 0 0 0];

dM = zeros([4 4 6]);

dM(1,4,1) = 1; % gradient due to Tx
dM(2,4,2) = 1; % gradient due to Ty
dM(3,4,3) = 1; % gradient due to Tz
dM(:,:,4) = Mt*dRx*Ry*Rz*(pi/180); % multiplied by (pi/180) to correct for units
dM(:,:,5) = Mt*Rx*dRy*Rz*(pi/180);
dM(:,:,6) = Mt*Rx*Ry*dRz*(pi/180);


% Make center of the image coordinates 0,0
[xd, yd, zd] = ndgrid(0:size(Iepi_t,1)-1, 0:size(Iepi_t,2)-1, 0:size(Iepi_t,3)-1);
X_origin = size(Iepi_t)/2;
xd = (xd - X_origin(1));
yd = (yd - X_origin(2));
zd = (zd - X_origin(3));

X = [xd(:) yd(:) zd(:) ones([length(xd(:)) 1])];
clear xd yd zd

rigid_err_grad = zeros(1,length(rigid_param));
for par_i=1:6
   M_curr = squeeze(dM(:,:,par_i));
   d_TX = X * transpose(M_curr); % transpose() b/c d_TX is transposed to make next step vectorized
   d_TX(:, 4) = [];    % size of [n_vox 3]
   
   rigid_err_grad(par_i) = transpose(d_TX(:)) * G_msk(:); % dot prod followed by sum over voxels
end

% Normalize rigid_grad by voxel count
rigid_err_grad = (rigid_err_grad.*rigid_scl)/numel(mask_Iepi_t);


% final grad
X_grad = [rigid_err_grad(:);  O_grad(:)];

end


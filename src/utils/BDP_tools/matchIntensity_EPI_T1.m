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


function [bfield_interp, epi_bfc, bfield_sp] = matchIntensity_EPI_T1(epiImg, sImg, im_res, mask_epi, mask_epi_less_csf, mask_sImg, opts)
% Estimates a smooth multiplicative field (bias-field) for b=0 image, which reduces intensity
% difference after INVERSION of b=0 image. Minimize number of voxels in images to minimize
% computational time. 
% 
% Requirements: 
%    * epiImg & sImg should be matched for resolution & should have intensities in
%      range [0 1]. 
%    * im_res is isotropic resolution in mm (images should be iterpolated on isotropic grid for
%      computation speed)
%    * Images should be coarsly aligned and resampled on same grid. Use after rigid-alignment
%      with INVERSION. 
%    * epi_mask_less_csf should 'ideally' match anatomical features in struct_mask
%    * sImg should be 'bias-field' corrected, for best results. 
%

defaultoptions = struct(...
   'nbins', 128, ... for INVERSION
   'parzen_width', 8, ...
   'spacing_mm', [40 40 40], ... control point spacing in mm
   'verbose', false, ...
   'penalty', 3.5e-2, ...
   'spring_penalty', [3e-3 0.8], ...
   'grad_step', 0.08, ... for numerical gradient
   'centralgrad', false, ...
   'MaxFunEvals', 20, ...
   'num_threads', 10, ...
   'step_size', 500, ...
   'step_size_scale_iter', 0.5, ...
   'mask_type', 'and', ...
   'error_res', 5, ...in mm; compute errors at downsampled resolution >= im_res
   'invert_epi', true, ...
   'invert_T1', true ...
   );

msk_thrs = 0.4;

if(~exist('opts','var')),
   opts = defaultoptions;
else
   tags = fieldnames(defaultoptions);
   for i=1:length(tags)
      if(~isfield(opts,tags{i})),  opts.(tags{i})=defaultoptions.(tags{i}); end
   end
   if(length(tags)~=length(fieldnames(opts))),
      warning('BDP:UnknownOptions', 'Unknown options found.');
   end
end

if isequal(size(epiImg), size(sImg))
   sz_in = size(epiImg);
else
   error('size of epiImg and sImg must be same');
end

if length(im_res)~=1
   error('im_res must be a scalar which is the isotropic resolution of all input images.')
end


% fix image size by padding zeros - make image size to be a perfect fit for control points
spacing = round(opts.spacing_mm./im_res);
num_seg = ceil(sz_in./spacing);
num_seg(num_seg<4) = 4; % Need atleast 4 segments for implemented bspline code

pix_diff = (num_seg.*spacing)+1 - sz_in;
st_pix = floor(pix_diff/2)+1;
end_pix = st_pix + sz_in - 1;

epiImg = padvolWithZeros(epiImg, pix_diff, st_pix, end_pix);
sImg = padvolWithZeros(sImg, pix_diff, st_pix, end_pix);
mask_epi = padvolWithZeros(mask_epi, pix_diff, st_pix, end_pix);
mask_epi_less_csf = padvolWithZeros(mask_epi_less_csf, pix_diff, st_pix, end_pix);
mask_sImg = padvolWithZeros(mask_sImg, pix_diff, st_pix, end_pix);

sizeI = size(epiImg);
O_grid = bspline_grid_generate_deformation(spacing, sizeI);
O_grid = O_grid(:,:,:,1) + 1; % scalar field with unit value everywhere
Osize = size(O_grid);


% Do NOT try multi-resolution.
% Multi-resolution apporach can not work because of implementation of map parameterized by
% b-splines in bspline_phantom_endpoint_deformation_only_x() & others. Further bias-field is
% supposed to be smooth - low resolution difference is enough.


% Compute downsampling parameters
error_res = opts.error_res; % resolution of vol difference in mm >= im_res
sgm = error_res./(im_res.* sqrt(8*log(2))); % sigma
spacing_comp = round(spacing.*(im_res/error_res));
t_sz = spacing_comp.*(Osize-1) + 1; % error computation image-size


if strcmpi(opts.mask_type, 'and')
   mask_err = imerode((mask_epi>msk_thrs)&(mask_sImg>msk_thrs), strel_sphere(1));
else
   mask_err = mask_epi + mask_sImg - (mask_epi .* mask_sImg);
end
indx_epi_less_csf = find(mask_epi_less_csf > msk_thrs);
indx_sImg = find(mask_sImg > msk_thrs);
[~,~,hist_sImg] = histcParzen(sImg(mask_sImg > msk_thrs), [0 1], opts.nbins, opts.parzen_width, opts.num_threads);


fmin_opts = struct(...
   'GradObj','on', ...
   'Display','off', ...
   'MaxFunEvals', opts.MaxFunEvals, ...
   'TolX', 0.01, ...
   'TolFun', 5e-6, ...
   'backtrack_c', 1e-7, ...
   'step_size', opts.step_size, ...
   'step_size_scale_iter', opts.step_size_scale_iter);
if opts.verbose, fmin_opts.Display='iter'; end

INV_opts = opts;
[INV_opts.cp_mask_indx, INV_opts.cp_col_index] = bspline_create_mask_index(Osize, spacing_comp, t_sz);
INV_opts = rmfield(INV_opts, 'penalty');
INV_opts.penalty.const = opts.penalty;
[INV_opts.penalty.V, INV_opts.penalty.cp_ind, INV_opts.penalty.Q] = ...
   bspline_penalty_repeated_endpoint_metal_bending_operator(Osize, spacing_comp, error_res);

% optimize
O_grid = O_grid(:);
O_grid = gradient_descent_backtrack(@(x)bspline_numeric_grad_bfield_INVERSION(x, Osize, spacing_comp, ...
   t_sz, sgm, epiImg, sImg, mask_epi>msk_thrs, mask_sImg>msk_thrs,...
   indx_epi_less_csf, indx_sImg, hist_sImg, mask_err, INV_opts), ...
   O_grid, fmin_opts);
O_grid = reshape(O_grid, Osize);


% compute on original grid using linear interpolation 
bfield_interp = bspline_repeated_endpoint_deformation_3d_double_only_x(O_grid, t_sz, spacing_comp, opts.num_threads);
bfield_interp = my_imresize3d(bfield_interp, [], sizeI, 'linear');
bfield_interp = resize_to_orig(bfield_interp, st_pix, end_pix);
epi_bfc = clip_intensity(resize_to_orig(epiImg, st_pix, end_pix).*bfield_interp);


% compute on original grid using fitted splines
if nargout>2
   bfield_sp = bspline_repeated_endpoint_deformation_3d_double_only_x(O_grid, sizeI, spacing, opts.num_threads);
   bfield_sp = resize_to_orig(bfield_sp, st_pix, end_pix); % removed padded zeros
end


end


function [O_err, O_grad] = bspline_numeric_grad_bfield_INVERSION(O_grid, sz_grid, spacing, t_sz, sgm, ...
   epiImg, sImg, mask_epi, mask_sImg, indx_epi_less_csf, indx_sImg, hist_sImg, mask_err, opt)
% INVERSION is always performed at native resolution.
% [err, grad] is computed after downsamping difference vol to t_sz using gaussian with sgm.

O_grid = reshape(O_grid, sz_grid);

% compute INVERSION map
bfield = bspline_repeated_endpoint_deformation_3d_double_only_x(O_grid, t_sz, spacing, opt.num_threads);
bfield = my_imresize3d(bfield, [], size(epiImg), 'linear');
epi_bfc = clip_intensity(bfield.*epiImg);
[~,~, histeqMap, T_grid] = histeqSmooth(1-epi_bfc(indx_epi_less_csf), ...
   hist_sImg, opt.nbins, opt.parzen_width, opt.num_threads);
invert_epi = griddedInterpolant(T_grid, histeqMap, 'cubic');

[~,~, histeqMap, T_grid] = histeqSmooth(1-sImg(indx_sImg), ...
   epi_bfc(indx_epi_less_csf), opt.nbins, opt.parzen_width);
invert_T1 = griddedInterpolant(T_grid, histeqMap, 'cubic');


% compute error & penalties
diff_img_sq = error_bfield_INVERSION(O_grid, spacing, t_sz, sgm, ...
   epiImg, sImg, mask_epi, mask_sImg, invert_epi, invert_T1, mask_err, opt);

[SP_err, SP_grad] = penaltyCPbfield(O_grid, opt.spring_penalty(2));

if opt.penalty.const(1)>0
   % [SO_error, SO_grad] = bspline_penalty_CP_diff_EPI(O_grid, spacing, t_res, opt.penalty(2:end));
   [SO_error, SO_grad] = bspline_penalty_repeated_endpoint_metal_bending_EPI(O_grid, ...
      opt.penalty.V, opt.penalty.cp_ind, opt.penalty.Q);
else
   SO_error = 0;
   SO_grad = zeros(sz_grid);
end
O_err = mean(diff_img_sq(:)) + opt.penalty.const(1)*SO_error + opt.spring_penalty(1)*SP_err;


% numerical gradient - does not incorporate regularization
temp_E_diff = zeros(prod(t_sz)+1, 1); % for fake zero voxel
Ogrid_size = size(O_grid);
if nargout>1
   O_grad = zeros(sz_grid);
   for zi=1:4
      for zj=1:4
         for zk=1:4
            
            [xx, yy, zz] = ndgrid(zi:4:Ogrid_size(1), zj:4:Ogrid_size(2), zk:4:Ogrid_size(3));
            CP_ind = sub2ind(Ogrid_size(1:3), xx(:), yy(:), zz(:)); % control-points in the set
            
            O_gradpx = O_grid;  % Move grid every fourth grid node
            O_gradpx(CP_ind) = O_gradpx(CP_ind) + opt.grad_step;
            
            E_gradpx = error_bfield_INVERSION(O_gradpx, spacing, t_sz, sgm, ...
               epiImg, sImg, mask_epi, mask_sImg, invert_epi, invert_T1, mask_err, opt);
            
            if opt.centralgrad
               O_gradqx = O_grid;
               O_gradqx(CP_ind) = O_gradqx(CP_ind) - opt.grad_step;
               
               E_gradqx = error_bfield_INVERSION(O_gradqx, spacing, t_sz, sgm, ...
                  epiImg, sImg, mask_epi, mask_sImg, invert_epi, invert_T1, mask_err, opt);
               E_diff = (E_gradpx-E_gradqx)/opt.grad_step/2;
               
            else
               
               E_diff = (E_gradpx-diff_img_sq)/opt.grad_step;
            end
            
            temp_E_diff(1:(end-1)) = E_diff(:);                        
            
            E_msk = sum(temp_E_diff(opt.cp_mask_indx(:, opt.cp_col_index(CP_ind))), 1);
            O_grad(CP_ind) = E_msk(:);
            
         end
      end
   end
   O_grad = O_grad ./ prod(spacing*4 + 1);
   O_grad = O_grad(:) + opt.spring_penalty(1)*SP_grad(:) + opt.penalty.const(1)*SO_grad(:);
end

end



function [diff_img_sq] = error_bfield_INVERSION(O_grid, spacing, t_sz, sgm, ...
   epiImg, sImg, mask_epi, mask_sImg, invert_epi, invert_T1, mask_err, opt)
% error after bias-field correction and inversion

% bias corr at native resolution 
bfield = bspline_repeated_endpoint_deformation_3d_double_only_x(O_grid, t_sz, spacing, opt.num_threads);
bfield = my_imresize3d(bfield, [], size(epiImg), 'linear');
epi_bfc = clip_intensity(bfield.*epiImg);

diff_img_sq = zeros(prod(t_sz), 1);
if opt.invert_epi
   epi_bfc_INV = zeros(size(epiImg));
   epi_bfc_INV(mask_epi) = invert_epi(1-epi_bfc(mask_epi));
   diff_img1 = (sImg-epi_bfc_INV).*mask_err;
   diff_img1 = my_imresize3d(imgaussian(diff_img1, sgm), [], t_sz, 'linear');
   diff_img_sq = diff_img1(:).^2;
end

if opt.invert_T1
   sImg_inv = zeros(size(epiImg));
   sImg_inv(mask_sImg) = invert_T1(1-sImg(mask_sImg));
   diff_img2 = (sImg_inv-epi_bfc).*mask_err;  
   diff_img2 = my_imresize3d(imgaussian(diff_img2, sgm), [], t_sz, 'linear');
   diff_img_sq = diff_img_sq + diff_img2(:).^2;
end
end



function out_vol = padvolWithZeros(vol, pix_diff, st_pix, end_pix)
sz_in = size(vol);
out_vol = zeros(sz_in + pix_diff);
out_vol(st_pix(1):end_pix(1), st_pix(2):end_pix(2), st_pix(3):end_pix(3)) = vol;
end


function vol_out = resize_to_orig(vol, st_pix, end_pix)
vol_out = vol(st_pix(1):end_pix(1), st_pix(2):end_pix(2), st_pix(3):end_pix(3));
end


function img = clip_intensity(img)
img(img<0) = 0;
img(img>1) = 1;
end

function [O_err, O_grad] = penaltyCPbfield(O_grid, thresh)
% computes spring-type penalty for scaling < thresh

osize = size(O_grid);
scl = O_grid; % scaling coefficient at each point

P = zeros(osize);
zmsk = scl<=0;
msk = scl<thresh & ~zmsk;

P(msk) = (1./(2*scl(msk))).^2 - (1./(2*thresh)).^2;
P(zmsk) = 1e10; % high penalty
O_err = mean(P(:));

O_grad = zeros(osize);
O_grad(msk) = -1*((1./(scl(msk).^3))./2)./numel(O_grid);
O_grad(zmsk) = -1e15;

end


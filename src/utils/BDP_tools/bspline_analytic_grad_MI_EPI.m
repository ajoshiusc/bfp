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


function [O_error, O_error_grad] = bspline_analytic_grad_MI_EPI(X, Osize1, O_grid, Spacing, Iepi, Ist, ...
   res, mask_Iepi, mask_Ist, opt)
% Calculates registration error and gradient for 1D b-spline non-rigid registration of two images

% Check/set input options
defaultoptions = struct( ...
   'penalty', [1 0.45 0.45 0.45], ...
   'verbose', false, ...
   'intensity_correct', true, ...
   'moving_mask', false);

if(~exist('opt','var')),
   opt=defaultoptions;
else
   tags = fieldnames(defaultoptions);
   for i=1:length(tags)
      if(~isfield(opt,tags{i})),  opt.(tags{i})=defaultoptions.(tags{i}); end
   end
end

% Check if there is any Mask input
if(~exist('mask_Imv','var') || ~exist('mask_Ist','var'))
   error('EPI_correct_bspline_registration_gradient:No_mask','Masks are required for registation')
end

msk_thresh = 0.4;
O_grid(:,:,:,1) = reshape(X,Osize1);

% Transform Iepi
t_x = bspline_phantom_endpoint_deformation_only_x(O_grid, size(Iepi), Spacing);
[x_g, y_g, z_g] = ndgrid(1:size(Iepi,1), 1:size(Iepi,2), 1:size(Iepi,3));
x_g = x_g + t_x;

if opt.intensity_correct
   [~, crct_grd] = gradient(t_x);
   Iepi_w = interpn(Iepi, x_g, y_g, z_g, 'linear', 0);
   Iepi_t = Iepi_w .* (1+crct_grd);
else
   Iepi_t = interpn(Iepi, x_g, y_g, z_g, 'linear', 0);   
end
mask_Iepi_t = interpn(double(mask_Iepi), x_g, y_g, z_g, 'nearest', 0);
clear x_g y_g z_g t_x

% Calculate the current registration error
if strcmpi(opt.mask_type, 'and')
   mask_err = mask_Iepi_t .* mask_Ist;
else
   mask_err = mask_Iepi_t + mask_Ist - (mask_Iepi_t .* mask_Ist);
end


[O_error, log_hist_mv, log_hist_st, log_hist12, MI_num, MI_den] = registration_error_mutual_info(Iepi_t, Ist, ...
   mask_err>msk_thresh, opt.CFn_opts);


% determine the gradient.
if ( nargout > 1 )
   if (opt.intensity_correct)
      O_error_grad = registration_gradient_mutual_info_intensity_correct(Spacing, res, Ist, Iepi_t, Iepi_w, crct_grd, ...
         size(O_grid), mask_err, log_hist_mv, log_hist12, MI_num, MI_den, opt);
   else
      O_error_grad = registration_gradient_mutual_info(Spacing, res, Ist, Iepi_t, size(O_grid), mask_Ist, mask_Iepi,...
         log_hist_mv, log_hist12, MI_num, MI_den, opt);
   end
   
   if(opt.verbose)
      disp(['Abs Mean Error Gradient' num2str(mean(abs(O_error_grad(:))))]);
      disp(['Abs Max  Error Gradient' num2str(max(abs(O_error_grad(:))))]);
   end
   O_error_grad = O_error_grad(:);
end

% Add penalty to total error
if opt.penalty(1)>0
   SO_error = bspline_penalty_CP_diff_EPI(O_grid(:,:,:,1), Spacing, res, opt.penalty(2:end));
   O_error = O_error + (opt.penalty(1)*SO_error);
end

if opt.verbose
   disp(['Error' num2str(O_error)]);
   disp(['Smoothness Error' num2str(SO_error)]);
end

end


function [err, log_histV, log_histU, log_histVU, MI_num, MI_den] = registration_error_mutual_info(Imv, Ist, mask, opts)
% registration error based on normalized mutual information. Estimates
% joint pdf  using cubic b-spline parzen window.
%
%      err = -1 * (H(V) + H(U)) / H(V,U)
%

Um = Ist(mask>0);
Vm = Imv(mask>0);
clear Ist Imv mask

histVU = mutual_histogram_parzen_variable_size_multithread_double(double(Vm), double(Um), ...
   double(0), double(1), double(opts.nbins), double(opts.win_width), double(opts.nthreads));
histVU = double(histVU);
histV = double(sum(histVU, 1));
histU = double(sum(histVU, 2));
clear Vm Um

if opts.log_lookup % log using lookup table
   
   log_histVU = zeros(size(histVU));
   ind_0 = (histVU < opts.range1(1));
   log_histVU(ind_0) = log(opts.range1(1));
   ind_1 = (histVU <= opts.range1(2)) & ~ind_0;
   ind_2 = (histVU > opts.range1(2));
   ind = histVU(ind_1) .* opts.scale1(ind_1);
   log_histVU(ind_1) = opts.log_lookup1(uint32(ind));
   log_histVU(ind_2) = log(histVU(ind_2)); % outside lookup table
   
   histV_histU = [histV(:); histU(:)];
   log_histV_histU = zeros(size(histV_histU));
   ind_0 = (histV_histU < opts.range2(1));
   log_histV_histU(ind_0) = log(opts.range2(1));
   ind_1 = (histV_histU <= opts.range2(2)) & ~ind_0;
   ind_2 = (histV_histU > opts.range2(2));
   ind = histV_histU(ind_1) .* opts.scale2(ind_1);
   log_histV_histU(ind_1) = opts.log_lookup2(uint32(ind));
   log_histV_histU(ind_2) = log(histV_histU(ind_2)); % outside lookup table
   
   log_histV = log_histV_histU(1:opts.nbins);
   log_histU = log_histV_histU(opts.nbins+1:end);
   
else % masked log computation
   
   c = log(opts.log_thresh);
   log_histVU = ones(size(histVU))*c;
   log_histV = ones(size(histV))*c;
   log_histU = ones(size(histU))*c;
   
   ind = find(histV > opts.log_thresh);
   log_histV(ind) = log(histV(ind));
   
   ind = find(histU > opts.log_thresh);
   log_histU(ind) = log(histU(ind));
   
   ind = find(histVU > opts.log_thresh);
   log_histVU(ind) = log(histVU(ind));
end

p1log = histV(:) .* log_histV(:);
p2log = histU(:) .* log_histU(:);
p12log = histVU .* log_histVU;

% Calculate amount of Information
MI_num = sum(p1log) + sum(p2log);
MI_den = sum(p12log(:));

% Studholme, Normalized mutual information
if(MI_den==0)
   MI_den = eps;
   warning('check this case.')
   keyboard
end

err = -1*(MI_num/MI_den);

end


function O_grad = registration_gradient_mutual_info(Spacing, res, Ist_in, Imv_t_in, O_grid_size, MaskIst, MaskImv_t,...
   log_hist_mv, log_hist12, MI_num, MI_den, opt)
% Calculates the analytical gradient of (-ve) mutual information (w.r.t.
% uniform cubic bspline control points) 
%
% Does NOT include intensity correction terms.
%
% ONLY calculates the gradient for x-coordinated of control points. 
%

global grad_img_indx grad_mask_count grad_mask_index grad_ind_sz;
global bins nthreads img_range;

if isempty(bins), bins = 128; end
if isempty(img_range), img_range = [0 1]; end

mask = (MaskImv_t>0 & MaskIst>0);
clear MaskImv_t MaskIst
[~, Gx] = gradient(Imv_t_in, res(1), res(2), res(3));
clear Gy Gz
Gx = Gx(mask);

Ist = Ist_in(mask)*((bins-1)/(img_range(2)-img_range(1)));
Imv_t = Imv_t_in(mask)*((bins-1)/(img_range(2)-img_range(1)));
clear Ist_in Imv_t_in

h_Ist = zeros([length(Ist) 4]);
d_h_Imv = zeros([length(Ist) 4]);
h_Ist_index = zeros([length(Ist) 4], 'int16');
d_h_Imv_index = zeros([length(Ist) 4], 'int16');

for m = 1:4
   m_bin = m-1;
   
   Ist_msk = logical(mod(floor((Ist-m_bin)/2), 2));
   Ist_p = mod(Ist-m_bin, 2);
   Ist_p(Ist_msk) = Ist_p(Ist_msk)-2;
   Ist_p = abs(Ist_p);
   
   Imv_msk = logical(mod(floor((Imv_t-m_bin)/2), 2));
   Imv_p = mod(Imv_t-m_bin, 2);
   Imv_p(Imv_msk) = Imv_p(Imv_msk)-2;
   Imv_p = -1*Imv_p;
   
   h_Ist_m = zeros(size(Ist_p));
   d_h_Imv_m = zeros(size(Ist_p));
   
   msk = Ist_p<1;
   h_Ist_m(msk) = 0.5*(Ist_p(msk).^3) - Ist_p(msk).^2 + 2/3;
   h_Ist_m(~msk) = Ist_p(~msk).^2 - (Ist_p(~msk).^3)/6 - 2*Ist_p(~msk) + 4/3;
   
   msk = abs(Imv_p)<1;
   d_h_Imv_m(msk) = (3/2)*(abs(Imv_p(msk)).^2) - 2*abs(Imv_p(msk));
   d_h_Imv_m(~msk) = (-0.5)*(abs(Imv_p(~msk)).^2) + 2*abs(Imv_p(~msk)) - 2;
   d_h_Imv_m(Imv_p<0) = -1*d_h_Imv_m(Imv_p<0);
   
   h_Ist_index_m = m_bin + 4*floor((Ist+2-m_bin)/4);
   d_h_Imv_index_m = m_bin + 4*floor((Imv_t+2-m_bin)/4);
   
   msk = (Ist-m_bin)<0 | (Ist-m_bin+2)>=bins;
   h_Ist_m(msk) = 0; % not included 
   h_Ist_index_m(msk) = 1; % dummy, not included 
   
   msk = (Imv_t-m_bin)<0 | (Imv_t-m_bin+2)>=bins;
   d_h_Imv_m(msk) = 0; % not included 
   d_h_Imv_index_m(msk) = 1;  % dummy, not included 
   
   h_Ist(:,m) = h_Ist_m;
   d_h_Imv(:,m) = d_h_Imv_m;
   h_Ist_index(:,m) = h_Ist_index_m + 1; 
   d_h_Imv_index(:,m) = d_h_Imv_index_m + 1;
end
clear Imv_t Ist

% Remove out of index indices, if any
msk = (d_h_Imv_index>bins) | (d_h_Imv_index<0);
d_h_Imv_index(msk) = 1; % dummy, not included 
d_h_Imv(msk) = 0; % not included 

msk = (h_Ist_index>bins) | (h_Ist_index<0);
h_Ist_index(msk) = 1; % dummy, not included 
h_Ist(msk) = 0; % not included 

clear d_h_Imv_m d_h_Imv_index_m h_Ist_m h_Ist_index_m msk Ist_p Imv_p Imv_msk Ist_msk

temp1 = 0:(1/Spacing(1)):1;
temp1 = 0.5*(temp1.^3) - (temp1.^2) + 4/6;
temp2 = (1-(1/Spacing(1))):(-1/Spacing(1)):0;
temp2 = (temp2.^3)/6;
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
clear bx by bz Bu Bv Bw O_grid1 temp1 temp2

b_prod = repmat(b_prod, floor((O_grid_size(1:3)-1)/4));


[p1, p2] = ndgrid(1:4,1:4);
p1 = p1(:);
p2 = p2(:);

% 1st term of gradient
d_P1 = (sum(log_hist_mv(d_h_Imv_index).*d_h_Imv, 2) .* sum(h_Ist, 2)) * MI_den; %/npixels;

% 2nd term of gradient
linearInd = sub2ind([bins bins], d_h_Imv_index(:,p1), h_Ist_index(:,p2));
d_P2 = (sum(log_hist12(linearInd) .* d_h_Imv(:,p1) .* h_Ist(:,p2), 2)) * MI_num; %/npixels;
clear log_hist_mv log_hist12 d_h_Imv_index d_h_Imv h_Ist linearInd h_Ist_index h_Ist p2 p1

mask_size = size(mask);
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
         
         d_Imv_control_point = Gx .* bspline_prod(mask);       
         d_P = zeros(mask_size);
         d_P(mask) = (d_P1-d_P2) .* d_Imv_control_point;         
         d_P = d_P(:);
         d_P(end+1) = 0;
         
         img_ind = opt.grad_img_indx(opt.grad_mask_index(zi,zj,zk):opt.grad_mask_index(zi,zj,zk) ...
            + (opt.grad_mask_count(zi,zj,zk)*opt.grad_ind_sz)-1);
         d_msk = d_P(img_ind);
         d_msk = sum(reshape(d_msk, opt.grad_ind_sz, []), 1);
         temp_msk = false(size(O_grad));
         temp_msk(zi:4:size(O_grad,1), zj:4:size(O_grad,2), zk:4:size(O_grad,3)) = true;
         O_grad(temp_msk) = d_msk(:);
         
      end
   end
end

% Normalize by number of voxels used for each control point 
O_grad = O_grad ./ prod(Spacing*4 + 1);
O_grad = O_grad/(MI_den^2);

% % Normalize O_grad by voxel count
% temp_msk = O_vox_count>0;
% O_grad(temp_msk) = O_grad(temp_msk)./ (O_vox_count(temp_msk));
% O_grad = O_grad/(MI_den^2);

% % old normalization
% O_grad = O_grad.*(bins/(sum(mask(:))*MI_den^2));

end


function O_grad = registration_gradient_mutual_info_intensity_correct(Spacing, res, Ist_in, Iepi_t_in, Iepi_w_in, crct_grd, O_grid_size, ...
   mask_err, log_hist_mv, log_hist12, MI_num, MI_den, opt)
% Calculates the analytical gradient of -ve NMI w.r.t. uniform cubic bspline control points

msk_thresh = 0.4;
mask_err = mask_err>msk_thresh;

[~, Gx] = gradient(Iepi_w_in);
clear Gy Gz
dTx_Gx = (1+crct_grd).*Gx;
dTx_Gx = dTx_Gx(mask_err);
Iepi_w_in = Iepi_w_in(mask_err);
clear Gx crct_grd

Ist = Ist_in(mask_err)*(opt.nbins-1);
Iepi_t = Iepi_t_in(mask_err)*(opt.nbins-1);
clear Ist_in Imv_t_in

% 4 = parzen window width??

h_Ist = zeros([length(Ist) 4]);
d_h_Imv = zeros([length(Ist) 4]);
h_Ist_index = zeros([length(Ist) 4], 'int16');
d_h_Imv_index = zeros([length(Ist) 4], 'int16');

for m = 1:4
   m_bin = m-1;
   
   Ist_msk = logical(mod(floor((Ist-m_bin)/2), 2));
   Ist_p = mod(Ist-m_bin, 2);
   Ist_p(Ist_msk) = Ist_p(Ist_msk)-2;
   Ist_p = abs(Ist_p);
   
   Imv_msk = logical(mod(floor((Iepi_t-m_bin)/2), 2));
   Imv_p = mod(Iepi_t-m_bin, 2);
   Imv_p(Imv_msk) = Imv_p(Imv_msk)-2;
   Imv_p = -1*Imv_p;
   
   h_Ist_m = zeros(size(Ist_p));
   d_h_Imv_m = zeros(size(Ist_p));
   
   msk = Ist_p<1;
   h_Ist_m(msk) = 0.5*(Ist_p(msk).^3) - Ist_p(msk).^2 + 2/3;
   h_Ist_m(~msk) = Ist_p(~msk).^2 - (Ist_p(~msk).^3)/6 - 2*Ist_p(~msk) + 4/3;
   
   msk = abs(Imv_p)<1;
   d_h_Imv_m(msk) = (3/2)*(abs(Imv_p(msk)).^2) - 2*abs(Imv_p(msk));
   d_h_Imv_m(~msk) = (-0.5)*(abs(Imv_p(~msk)).^2) + 2*abs(Imv_p(~msk)) - 2;
   d_h_Imv_m(Imv_p<0) = -1*d_h_Imv_m(Imv_p<0);
   
   h_Ist_index_m = m_bin + 4*floor((Ist+2-m_bin)/4);
   d_h_Imv_index_m = m_bin + 4*floor((Iepi_t+2-m_bin)/4);
   
   msk = (Ist-m_bin)<0 | (Ist-m_bin+2)>=opt.nbins;
   h_Ist_m(msk) = 0; % not included 
   h_Ist_index_m(msk) = 1; % dummy, not included 
   
   msk = (Iepi_t-m_bin)<0 | (Iepi_t-m_bin+2)>=opt.nbins;
   d_h_Imv_m(msk) = 0; % not included 
   d_h_Imv_index_m(msk) = 1;  % dummy, not included 
   
   h_Ist(:,m) = h_Ist_m;
   d_h_Imv(:,m) = d_h_Imv_m;
   h_Ist_index(:,m) = h_Ist_index_m + 1; 
   d_h_Imv_index(:,m) = d_h_Imv_index_m + 1;
end
clear Imv_t Ist
% Remove out of index indices, if any
msk = (d_h_Imv_index>opt.nbins) | (d_h_Imv_index<0);
d_h_Imv_index(msk) = 1; % dummy, not included 
d_h_Imv(msk) = 0; % not included 

msk = (h_Ist_index>opt.nbins) | (h_Ist_index<0);
h_Ist_index(msk) = 1; % dummy, not included 
h_Ist(msk) = 0; % not included 

clear d_h_Imv_m d_h_Imv_index_m h_Ist_m h_Ist_index_m msk Ist_p Imv_p Imv_msk Ist_msk


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

% derivative of cubic bspline kernal (direct form, NOT from basis)
temp1 = 1/Spacing(1):(1/Spacing(1)):1;
temp1 = 1.5*(temp1.^2) - 2*temp1;
temp2 = (1-(1/Spacing(1))):(-1/Spacing(1)):0;
temp2 = -0.5*(temp2.^2);
dBu = [temp1 temp2];
dBu = [-1*dBu(end:-1:1) 0 dBu];

xmin = -2*Spacing(1);
xmax = 2*Spacing(1)-1;

ymin = -2*Spacing(2);
ymax = 2*Spacing(2)-1;

zmin = -2*Spacing(3);
zmax = 2*Spacing(3)-1;

[bx, by, bz] = ndgrid((xmin:xmax), (ymin:ymax), (zmin:zmax));
b_prod = Bu(abs(bx)+1) .* Bv(abs(by)+1) .* Bw(abs(bz)+1);
db_prod = (dBu(bx-xmin+1) .* Bv(abs(by)+1) .* Bw(abs(bz)+1))./Spacing(1);
clear bx by bz Bu Bv Bw O_grid1 temp1 temp2

b_prod = repmat(b_prod, floor((O_grid_size(1:3)-1)/4));
db_prod = repmat(db_prod, floor((O_grid_size(1:3)-1)/4));


[p1, p2] = ndgrid(1:4,1:4);
p1 = p1(:);
p2 = p2(:);

% 1st term of gradient
d_P1 = (sum(log_hist_mv(d_h_Imv_index).*d_h_Imv, 2) .* sum(h_Ist, 2)) * MI_den; %/npixels;

% 2nd term of gradient
linearInd = sub2ind([opt.nbins opt.nbins], d_h_Imv_index(:,p1), h_Ist_index(:,p2));
d_P2 = (sum(log_hist12(linearInd) .* d_h_Imv(:,p1) .* h_Ist(:,p2), 2)) * MI_num; %/npixels;
dP1_dP2_diff = (d_P1-d_P2);

clear log_hist_mv d_h_Imv_index d_h_Imv h_Ist h_Ist_index linearInd log_hist12 p1 p2 d_P1 d_P2

mask_size = size(mask_err);
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
         bspline_prod = bspline_prod(mask_err);
         
         dbspline_prod = circshift(db_prod, ([zi zj zk]-3).*Spacing);
         pixs_rem = (mod(O_grid_size(1:3)-1, 4).*Spacing) + 1;
         pixs_rem(pixs_rem<0) = 0;
         dbspline_prod(end+1:end+pixs_rem(1), :, :) = dbspline_prod(1:pixs_rem(1), :, :);
         dbspline_prod(:, end+1:end+pixs_rem(2), :) = dbspline_prod(:, 1:pixs_rem(2), :);
         dbspline_prod(:, :, end+1:end+pixs_rem(3)) = dbspline_prod(:, :, 1:pixs_rem(3));
         dbspline_prod = dbspline_prod(mask_err);
         
         d_Imv_control_point = (dTx_Gx .* bspline_prod) + (Iepi_w_in .* dbspline_prod);       
         d_P = zeros(mask_size);
         d_P(mask_err) = dP1_dP2_diff .* d_Imv_control_point;
         clear bspline_prod dbspline_prod d_Imv_control_point
         
         d_P = d_P(:);
         d_P(end+1) = 0;
         
         img_ind = opt.grad_img_indx(opt.grad_mask_index(zi,zj,zk):opt.grad_mask_index(zi,zj,zk) ...
            + (opt.grad_mask_count(zi,zj,zk)*opt.grad_ind_sz)-1);
         d_msk = d_P(img_ind);
         d_msk = sum(reshape(d_msk, opt.grad_ind_sz, []), 1);
         temp_msk = false(size(O_grad));
         temp_msk(zi:4:size(O_grad,1), zj:4:size(O_grad,2), zk:4:size(O_grad,3)) = true;
         O_grad(temp_msk) = d_msk(:);        
      end
   end
end

% Normalize by number of voxels used for each control point 
%O_grad = O_grad ./ prod(Spacing*4 + 1);
O_grad = O_grad/(MI_den^2);

% % Normalize O_grad by voxel count
% temp_msk = O_vox_count>0;
% O_grad(temp_msk) = O_grad(temp_msk)./ (O_vox_count(temp_msk));
% O_grad = O_grad/(MI_den^2);

% old normalization
% O_grad = O_grad.*(opt.nbins/(sum(mask(:))*MI_den^2));

end



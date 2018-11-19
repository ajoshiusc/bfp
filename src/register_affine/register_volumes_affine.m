function [M, x, Ireg] = register_volumes_affine(Imoving, Istatic, mask_moving, mask_static, im_res, opts)
% im_res is isotropic voxel resolution in mm. All images are assumed to be (re)sampled on
% isotropic grid. See defaultoptions below for all options.  

% Process inputs & sanity checks
defaultoptions = struct(...
   'dof', 6, ... 1/2/3/6/12
   'axis_translation', [], ...
   'x_scale', [], ...
   'similarity', '', ... 'cr', 'sd', 'mi', 'INVERSION', 'BDP'
   'verbose', true, ...
   'step_size', [], ... step size for optimization
   'grad_step_size', 1e-1', ... step size for numerical gradient
   'centralgrad', true, ... if to use central gradient while computing numerical gradients
   'mask_op', 'and', ... Error mask: and / or of indivdual masks
   'apodize_cost', true, ...
   'max_func_evals', 600, ... max # of func evals at each scale/resolution
   'x_init', [], ... initial estimate of x
   'init_method', 'search', ... none/search/mask -- method of finding initial estimate of x
   'search_range', [60 60 60], ... in degrees, defaults to -60 to 60
   'search_delta', [15 15 15],... search step in degrees
   'search_imres', 5, ... approx isotropic image resolution in mm for search
   'CFn_opts', [], ... structure of cost function options
   'nthreads', 6 ...
   );

interp_method = 1; % linear interpolation and outside pixels set to zero
msk_thrs = 0.25; % threshold to convert a contineous mask to binary


if(~exist('opts','var')),
   opts = defaultoptions;
else
   tags = fieldnames(defaultoptions);
   for i = 1:length(tags)
      if ~isfield(opts, tags{i})
         opts.(tags{i}) = defaultoptions.(tags{i});
      end
   end
   if length(tags)~=length(fieldnames(opts))
      warning('BDP:UnknownOptions','unknown options found');
   end
end

if length(im_res)>1
   error('im_res should be scalar. All images are assumed to be (re)sampled on isotropic grid.')
end

if ~isequal(size(Istatic), size(Imoving))
   error('Images should have same size.')
else
   Istatic = double(Istatic);
   Imoving = double(Imoving);
end

if isempty(opts.similarity)
   error('Please specify the similarity of the images.');
else
   opts.similarity = lower(opts.similarity);
end

if opts.dof<3 && length(opts.axis_translation)~=opts.dof
   error('opts.axis_translation is not set correctly for low dof registration');
end

if opts.verbose
   areg_tic = tic();
end

% set parameter scaling
if isempty(opts.x_scale)
   scale = ones(opts.dof, 1);
elseif ~isequal(length(opts.x_scale), opts.dof)
   error('length of opts.x_scale must be same as opts.dof')
else
   scale = opts.x_scale(:);
end

% set step-size depending on similarity measure
if isempty(opts.step_size)
   switch opts.similarity
      case {'sd', 'inversion', 'bdp'}
         opts.step_size = 1000;
      case {'mi', 'cr'}
         opts.step_size = 35;
      otherwise
         error('Unknown similarity: %s', opts.similarity)
   end
end

% setup cost-function options
if isempty(opts.CFn_opts)
   switch opts.similarity
      case{'inversion', 'bdp'}
         opts.CFn_opts = struct('nbins', 300, 'win_width', 30, 'hybrid_mask_op', false);
         opts.CFn_opts.mi = struct('nbins', 150, 'win_width', 8, 'nthreads', opts.nthreads, 'log_lookup', false, 'log_thresh', 1/(150*150*10));         
      case{'mi'}
         opts.CFn_opts = struct('nbins', 150, 'win_width', 8, 'nthreads', opts.nthreads, 'log_lookup', false, 'log_thresh', 1/(150*150*10) );
      case{'cr'}
         opts.CFn_opts = struct('nbins', 128, 'nthreads', opts.nthreads);
      case{'sd'}
         opts.CFn_opts = struct('hybrid_mask_op', true);
   end
end

% set initial rigid parameters
if isempty(opts.x_init)
   x = zeros(opts.dof, 1);
elseif ~isequal(length(opts.x_init), opts.dof)
   error('length of opts.x_init must be same as opts.dof')
else
   x = opts.x_init(:);
end

% set search parameters
if ~ismember(opts.init_method, {'none', 'search', 'mask'})
   error('Invalid opts.init_method: %s', opts.init_method)
end

if strcmpi(opts.init_method, 'search')
   if opts.dof<6
      opts.search_init = false;
      opts.init_method = 'none'; % search with dof<6 does not make sense
   else
      opts.search_init = true;
      if length(opts.search_range)==1
         opts.search_range = opts.search_range*[1 1 1];
      elseif length(opts.search_range)~=3
         error('opts.search_range must be of length 3 or 1')
      end
      
      if length(opts.search_delta)==1
         opts.search_delta = opts.search_delta*[1 1 1];
      elseif length(opts.search_delta)~=3
         error('opts.search_delta must be of length 3 or 1')
      end
   end
else
   opts.search_init = false;
end

%% Initialization - search or mask based
if ~strcmpi(opts.init_method, 'none')
   if opts.search_init  % search based
      fprintf('\nSearch-based initialization... '); drawnow;
      if opts.verbose
         init_tic = tic();
      end
      
      % INVERSION, if required (at native resolution)
      if strcmpi(opts.similarity, 'inversion') || strcmpi(opts.similarity, 'bdp')
         st_inv = zeros(size(Istatic));
         st_inv(mask_static>msk_thrs) = histeqSmooth(1-Istatic(mask_static>msk_thrs), ...
            Imoving(mask_moving>msk_thrs), opts.CFn_opts.nbins, opts.CFn_opts.win_width);
         static_img = st_inv;
         clear st_inv
      else
         static_img = Istatic;
      end
            
      % smooth & resample image
      fwhm = opts.search_imres/im_res;
      if fwhm>1
         sg = fwhm/sqrt(8*log(2));
         rs_factor = 1.2/fwhm; % little bit smaller than actual resolution
         sm_mov_img = my_imresize3d(imgaussian(Imoving,sg), rs_factor ,[], 'linear');
         sm_stat_img = my_imresize3d(imgaussian(static_img,sg),rs_factor ,[], 'linear');         
         mask_moving_small = my_imresize3d(mask_moving, rs_factor, [], 'linear');
         mask_static_small = my_imresize3d(mask_static, rs_factor, [], 'linear');
      
      else
         sm_mov_img = Imoving;
         sm_stat_img = static_img;
         mask_moving_small = mask_moving;
         mask_static_small = mask_static;
      end
      
      % param for search grid
      th_x = 0:opts.search_delta(1):opts.search_range(1);
      th_x = [-1*th_x(end:-1:2) th_x];      
      th_y = 0:opts.search_delta(2):opts.search_range(2);
      th_y = [-1*th_y(end:-1:2) th_y];      
      th_z = 0:opts.search_delta(3):opts.search_range(3);
      th_z = [-1*th_z(end:-1:2) th_z];      
      [th_x, th_y, th_z] = ndgrid(th_x, th_y, th_z);
      
      par_search = zeros(6, numel(th_x));
      par_search(4:6, :) = [th_x(:) th_y(:) th_z(:)]'; % only rotation, centroid should already be matched
      
      % apodise cost
      sm_stat_img = sm_stat_img .* imgaussian(double(mask_static_small), 1.2);
      sm_mov_img = sm_mov_img .* imgaussian(double(mask_moving_small), 1.2);
      
      % search min cost
      err_msk = true(size(sm_mov_img)); % full volume similarity
      err = zeros(numel(th_x), 1);
      for k = 1:numel(th_x)
         M = par2M(x(1:6) + par_search(:,k), [1 1 1 1 1 1]');
         Imoved = affine_transform_3dvol_double(sm_mov_img, M, 1, opts.nthreads);
         err(k) = error_affine_reg(Imoved, sm_stat_img, err_msk, opts.similarity, opts.CFn_opts);
      end
      [~, ind] = min(err);
      x(1:6) = x(1:6) + par_search(:,ind);
      
      clear Imoved sm_stat_img sm_mov_img mask_moving_small mask_static_small ...
         th_x th_y th_z err err_msk par_search M      
      
   else % mask based init
      fprintf('\nMask-based initialization... '); drawnow;
      if opts.verbose
         init_tic = tic();
      end
      optim = struct(...
         'GradObj','on', ...
         'Display','off', ...
         'MaxFunEvals', 200, ...
         'TolX', 1e-4, ...
         'TolFun', 1e-14, ...
         'step_size', opts.step_size, ...
         'step_size_scale_iter', 0.5);
      
      mask_reg_options = struct( ...
         'grad_step', opts.grad_step_size, ...
         'centralgrad', opts.centralgrad, ...
         'axis_translation', opts.axis_translation, ...
         'nthreads', opts.nthreads ...
         );
      
      if opts.verbose, optim.Display='iter'; end
      sm_mov_mask = imgaussian(double(my_imresize3d(mask_moving,0.5,[],'linear')), 2);
      sm_stat_mask = imgaussian(double(my_imresize3d(mask_static,0.5,[],'linear')), 2);
      x = gradient_descent_adjust_step(@(x)mask_align_error(x, scale, sm_mov_mask, sm_stat_mask, mask_reg_options), x, optim);
      if length(x)<3
         x = x*2;
      else
         x(1:3) = x(1:3)*2;
      end
      clear sm_mov_mask sm_stat_mask mask_reg_options
   end
   
   if opts.verbose      
      disp(['Completed initialization in ' num2str(toc(init_tic)) ' seconds.']); drawnow();
   end
   
end

%% Affine registration

% INVERSION, if required
if strcmpi(opts.similarity, 'inversion') || strcmpi(opts.similarity, 'bdp')
   mask_static_tight = imerode(mask_static, strel_sphere(2));   
   st_inv = zeros(size(Istatic));
   st_inv(mask_static_tight>msk_thrs) = histeqSmooth(1-Istatic(mask_static_tight>msk_thrs), ...
      Imoving(mask_moving>msk_thrs), opts.CFn_opts.nbins, opts.CFn_opts.win_width);
   Istatic_a = st_inv.*mask_static_tight;
   clear st_inv mask_static_tight
else
   Istatic_a = Istatic;
end

optim_opts = struct( ...
   'GradObj','on', ...
   'Display','off', ...
   'MaxFunEvals', opts.max_func_evals, ...
   'TolX', 1e-4, ...
   'TolFun', 1e-14, ...
   'backtrack_c', 1e-5, ...
   'step_size', opts.step_size, ...
   'step_size_scale_iter', 0.5);

if opts.verbose, optim_opts.Display='iter';  end

reg_options = struct( ...
   'grad_step', opts.grad_step_size, ...
   'interp_method', interp_method, ...
   'similarity', opts.similarity, ...
   'centralgrad', opts.centralgrad, ...
   'mask_op', opts.mask_op, ...
   'axis_translation', opts.axis_translation, ...
   'CFn_opts', opts.CFn_opts, ...
   'nthreads', opts.nthreads ...
   );

scale_res    = [10   8    6    4    2.5  im_res]; % in mm
resize_scale = [0.4  0.4  0.5  0.8  1      1   ];
sgma = scale_res./(im_res*sqrt(8*log(2)));

% % remove finer resolution registration than im_res, saves time!
% tmp_msk = scale_res>=im_res;
% if ~any(tmp_msk) % im_res has lower resolution than 10 mm image
%    tmp_msk(1) = true;
%    msg = ['\nRegistration-resolution is set to a value lower than ' num2str(scale_res(1)) 'mm resolution. This function ' ...
%       'can not run registration a coarser resolution than that. So, we will try to register images at ' ...
%       num2str(scale_res(1)) 'mm resolution.\n'];
%    bdp_linewrap(msg);
% end
% if im_res>scale_res(end-1) % fix for equality check (>=) above
%    tmp_msk(end) = false;
% end
% resize_scale = resize_scale(tmp_msk);
% scale_res = scale_res(tmp_msk);
% sgma = sgma(tmp_msk);


for irf=1:length(resize_scale)
   fprintf('\nScale %d of %d running...', irf, length(resize_scale))
   if length(x)<3
      x = x*resize_scale(irf);
   else
      x(1:3) = x(1:3)*resize_scale(irf);
   end
   
   % resample image
   sg = sgma(irf);
   mask_moving_small = my_imresize3d(mask_moving, resize_scale(irf), [], 'linear');
   mask_static_small = my_imresize3d(mask_static, resize_scale(irf), [], 'linear');
   Imoving_small = my_imresize3d(imgaussian(Imoving, sg), resize_scale(irf), [],'linear');
   Istatic_small = my_imresize3d(imgaussian(Istatic_a, sg), resize_scale(irf), [],'linear');
      
   if (strcmpi(opts.similarity, 'inversion') || strcmpi(opts.similarity, 'bdp')) ...
         && opts.CFn_opts.hybrid_mask_op
      if scale_res(irf)>6
         reg_options.mask_op = 'or';
      else
         reg_options.mask_op = 'and';
      end
   end
   
   if opts.apodize_cost
      Imoving_small = Imoving_small .* imgaussian(double(mask_moving_small), 1.2);
      Istatic_small = Istatic_small .* imgaussian(double(mask_static_small), 1.2);
   end

   if opts.verbose
      txt = sprintf('\n\nResize scale = %f, Step size = %f, Res = %f', resize_scale(irf), ...
         optim_opts.step_size, scale_res(irf));
      disp(txt);
   end

   x = gradient_descent_backtrack(@(x) ...
      error_grad_affine_reg(x, scale, Imoving_small, Istatic_small, ...
      mask_moving_small, mask_static_small, reg_options), x, optim_opts);     
   
   if length(x)<3
      x = x/resize_scale(irf);
   else
      x(1:3) = x(1:3)/resize_scale(irf);
   end   
end

%% BDP mode: Refine with MI based measure

if strcmpi(opts.similarity, 'bdp')
   fprintf('\n');
   optim_opts = struct( ...
      'GradObj','on', ...
      'Display','off', ...
      'MaxFunEvals', 30, ...
      'TolX', 1e-3, ...
      'TolFun', 1e-10, ...
      'step_size', 50, ...
      'backtrack_c', 1e-7, ...
      'step_size_scale_iter', 0.5);   
   if(opts.verbose), optim_opts.Display='iter';  end
   
   reg_options = struct( ...
      'grad_step', opts.grad_step_size, ...
      'interp_method', interp_method, ...
      'similarity', 'mi', ...
      'centralgrad', opts.centralgrad, ...
      'mask_op', 'and', ...
      'axis_translation', opts.axis_translation, ...
      'CFn_opts', opts.CFn_opts.mi, ...
      'nthreads', opts.nthreads ...
      );
   
   scale_res    = [3.5  im_res]; % in mm; should be >= im_res
   resize_scale = [0.8    1   ];
   sgma = scale_res./(im_res*sqrt(8*log(2)));
   
   for irf = 1:length(resize_scale)
      fprintf('\nRefine: Scale %d of %d running...', irf, length(resize_scale))
      
      if length(x)<3
         x = x*resize_scale(irf);
      else
         x(1:3) = x(1:3)*resize_scale(irf);
      end
      
      % resample image
      sg = sgma(irf);
      mask_moving_small = my_imresize3d(mask_moving, resize_scale(irf), [], 'linear');
      mask_static_small = my_imresize3d(mask_static, resize_scale(irf), [], 'linear');
      Imoving_small = my_imresize3d(imgaussian(Imoving, sg), resize_scale(irf), [],'linear');
      Istatic_small = my_imresize3d(imgaussian(Istatic, sg), resize_scale(irf), [],'linear');      
      
      if opts.apodize_cost
         Imoving_small = Imoving_small .* imgaussian(double(mask_moving_small), 1.2);
         Istatic_small = Istatic_small .* imgaussian(double(mask_static_small), 1.2);
      end
      
      if opts.verbose
         txt = sprintf('\n\nResize scale = %f, Step size = %f, Res = %f', resize_scale(irf), ...
            optim_opts.step_size, scale_res(irf));
         disp(txt);
      end
      
      x = gradient_descent_backtrack(@(x) ...
         error_grad_affine_reg(x, scale, Imoving_small, Istatic_small, ...
         mask_moving_small, mask_static_small, reg_options), x, optim_opts);
      
      if length(x)<3
         x = x/resize_scale(irf);
      else
         x(1:3) = x(1:3)/resize_scale(irf);
      end      
   end
end


%% Transform the input image
M = par2M(x, scale, opts);
if nargout>2
   Ireg = affine_transform_3dvol_double(double(Imoving), double(M), 3, opts.nthreads);
end

% End time measurement
if opts.verbose
   disp(['Completed register_volumes_affine in ' num2str(toc(areg_tic)) ' seconds.']); drawnow();
end

clearvars -except M x Ireg 
end


function [err, err_grad] = error_grad_affine_reg(par, scale, Imv, Ist, Imv_mask, Ist_mask, opts)

M = par2M(par, scale, opts);

if opts.interp_method == 4  % NN interpolation
   sizeImv = size(Imv);
   [x,y,z] = ndgrid(0:(sizeImv(1)-1), 0:(sizeImv(2)-1), 0:(sizeImv(3)-1));
   xd = x-(sizeImv(1)/2); clear x;
   yd = y-(sizeImv(2)/2); clear y;
   zd = z-(sizeImv(3)/2); clear z;
   
   % Calculate the backwards transformation fields 
   % add 1+(sizeImv/2) to following to get back in original grid starting from index 1
   Bx = M(1,1)*xd + M(1,2)*yd + M(1,3)*zd + M(1,4);
   By = M(2,1)*xd + M(2,2)*yd + M(2,3)*zd + M(2,4);
   Bz = M(3,1)*xd + M(3,2)*yd + M(3,3)*zd + M(3,4);
   
   % consider replacing follwing with griddedInterpolant for efficient processing
   Imoved = interpn(xd, yd, zd, Imv, Bx, By, Bz, 'nearest', 0);
else
   mode = 1; % linear interpolation and outside pixels set to zero
   Imoved = affine_transform_3dvol_double(double(Imv),double(M),double(mode), double(opts.nthreads));
end

% registration error
Imv_t_mask = affine_transform_3dvol_double(double(Imv_mask),double(M),double(1), double(opts.nthreads));

if strcmpi(opts.mask_op, 'and')
   err_msk = Imv_t_mask .* Ist_mask;
elseif strcmpi(opts.mask_op, 'or')
   err_msk = Imv_t_mask + Ist_mask - (Imv_t_mask .* Ist_mask);
else
   error('mask_op must be and/or')
end
err = error_affine_reg(Imoved, Ist, err_msk, opts.similarity, opts.CFn_opts);


% Numerical Gradients
if(nargout>1)
   err_grad = zeros(size(par));
   for i = 1:length(par)
      par2 = par;
      par2(i) = par(i) + opts.grad_step;
      err1 = error_grad_affine_reg(par2, scale, Imv, Ist, Imv_mask, Ist_mask, opts);
      
      if opts.centralgrad
         par2(i) = par(i) - opts.grad_step;
         err2 = error_grad_affine_reg(par2, scale, Imv, Ist, Imv_mask, Ist_mask, opts);
         err_grad(i) = (err1-err2)/opts.grad_step/2;
      else
         err_grad(i) = (err1-err)/opts.grad_step;
      end
   end
end
end


function err = error_affine_reg(Imv, Ist, err_msk, similarity, CFn_opts)
% computes error for different similarity measures

switch similarity
   case {'sd', 'inversion', 'bdp'}
      err = (err_msk(:)' * ((Imv(:)-Ist(:)).^2))/numel(err_msk); % normalized for volume
      
   case 'mi'
      err = my_normalized_mutual_information(Imv, Ist, err_msk>0, CFn_opts);
      
   case 'cr'
      err = neg_correlation_ratio(Imv, Ist, err_msk>0, CFn_opts);
end

end


function err = neg_correlation_ratio(x, y, msk, opts)
% Computes (1-correlation ratio)
% See Roche et al. 1998, Jenkinson & Smith 2001

x = x(msk>0);
y = y(msk>0);

var_y = var(y);
N = numel(y);
num = 0;
bin_edges = linspace(0, 1, opts.nbins);
for k = 2:length(bin_edges)
   if k==length(bin_edges)
      bin_msk = (x>=bin_edges(k-1)) & (x<=bin_edges(k));
   else
      bin_msk = (x>=bin_edges(k-1)) & (x<bin_edges(k));
   end
   
   y_bin = y(bin_msk);
   n = numel(y_bin);
   if n>0
      var_bin = var(y_bin);
   else
      var_bin = 0;
   end
   num = num + n*var_bin;
end

err = num/(N*var_y);
end


function [err, log_histV, log_histU, log_histVU, MI_num, MI_den] = my_normalized_mutual_information(V, U, Mask, opts)
% Returns normalized mutual information. Estimates joint pdf using cubic b-spline parzen window.
%   err = -1 * (H(V) + H(U)) / H(V,U)
%

% Remove unmasked pixels
Vm = V(Mask);
Um = U(Mask);

[histVU] = mutual_histogram_parzen_variable_size_multithread_double( double(Vm), double(Um), ...
   double(0), double(1), double(opts.nbins), ...
   double(opts.win_width), double(opts.nthreads) );

histV = double(sum(histVU, 1));
histU = double(sum(histVU, 2));
histVU = double(histVU);

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


function [err, err_grad] = mask_align_error(par, scale, Imv, Ist, opts)
% compute error between the masks using squared difference as similarity metric

M = par2M(par, scale, opts);

mode = 1; % linear interpolation and outside pixels set to zero
Imoved = affine_transform_3dvol_double(double(Imv),double(M),double(mode), double(opts.nthreads));

err = mean((double(Ist(:))-double(Imoved(:))).^2);

% Numerical Gradients
if(nargout>1)
   err_grad = zeros(size(par));
   for i = 1:length(par)
      par2 = par;
      par2(i) = par(i) + opts.grad_step;
      err1 = mask_align_error(par2, scale, Imv, Ist, opts);
      
      if opts.centralgrad
         par2(i) = par(i) - opts.grad_step;
         err2 = mask_align_error(par2, scale, Imv, Ist, opts);
         err_grad(i) = (err1-err2)/opts.grad_step/2;
      else
         err_grad(i) = (err1-err)/opts.grad_step;
      end
   end
end

end


function M = par2M(x, scale, opt)
% Computes affine matrix from parameter(x) and scale

x = x.*scale;
if length(x) < 3
   t = zeros(1,3);
   t(opt.axis_translation) = x;
   M = par2affineMat(t);
elseif length(x) == 3
   M = par2affineMat(x(1:3));
elseif length(x) == 6
   M = par2affineMat(x(1:3), x(4:6));
elseif length(x) == 12
   M = par2affineMat(x(1:3), x(4:6), [], x(7:12));
else
   error('Length of x must be 1, 2, 3, 6 or 12');
end

end




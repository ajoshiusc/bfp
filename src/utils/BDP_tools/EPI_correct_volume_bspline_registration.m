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


function [Iepi_corr, Istruct_reg, O_trans, spacing, rigid_param, t_x] = EPI_correct_volume_bspline_registration( ...
   Iepi, IepiBFC, Istruct, res, opts)
% Performs 1D non-linear registration for EPI distortion correction. First
% dimension of the input image should be the phase encode dimension.

defaultoptions = struct(...
   'Penalty', [1 0.00001 0.00001], ... [alpha kx bigkx ky], alpha is scaling of overall regularization term
   'similarity', 'inversion-epi-t1', ... bdp / inversion-epi / inversion-t1 / mi / sd
   'MaxRef', 1,... Should be either 1 or 2, Max number of bspline refinements
   'Spacing',[],... space b/w control points in unit of voxels
   'MaskEPI',[],...     tight mask of EPI image containing only brains & CSF
   'MaskEPILessCSF',[],... Tighter mask than MaskEPI, with CSF near boundary removed
   'MaskStruct',[],...     tight mask for structural image containing only brain (no CSF near boundaries)
   'step_size', [], ...    step size for gradient descent
   'step_size_scale_iter', 0.5,... % for gradient descent
   'rigid_scl', 1, ... Rigid scale; scalar or vector of length 6
   'intensity_correct',true,... % This would be largely ignored - distortion is always estimated with correct model.
   ...                          % Its too much work to derive, support and debug an almost-would-never-be-used feature!
   ...                          % This input-option is retained for historical reasons.
   ...                          % If required, final application of the Jacobian term can be avoided at the application of distortion field.
   'MaxFunEvals', 200, ... max number of gradient evaluation
   'verbose', true,... for optimization iterations
   'debug', false, ... for plots
   'nthreads', 6, ...
   'histeqInterpolant', [], ...
   'histeqGradInterpolant', [], ...
   'CFn_opts', [] ...
   );

if ~exist('opts','var')
   opts = defaultoptions;
else
   tags = fieldnames(defaultoptions);
   for i=1:length(tags),
      if ~isfield(opts,tags{i})
         opts.(tags{i}) = defaultoptions.(tags{i});
      end
   end
   if length(tags)~=length(fieldnames(opts))
      warning('BDP:UnknownOptions', 'unknown options found!');
   end
end

if numel(opts.Penalty)~=3 && numel(opts.Penalty)~=4
   error('Penalty length must be 3 or 4.')
end

% set step-size depending on similarity measure
if isempty(opts.step_size)
   switch opts.similarity
      case {'inversion-epi', 'sd'}
         opts.step_size = 2500;
      case {'inversion-t1', 'inversion-epi-t1', 'bdp'}
         opts.step_size = 4500;
         % case {'mi', 'cr'}
         %   opts.step_size = 120;
      otherwise
         error('Unknown similarity: %s', opts.similarity)
   end
end

% setup cost-function options
if isempty(opts.CFn_opts)
   switch opts.similarity
      case {'inversion-epi', 'inversion-t1', 'inversion-epi-t1'}
         opts.CFn_opts = struct('nbins', 600, 'win_width', 50);
      case {'bdp'}
         opts.CFn_opts = struct('nbins', 600, 'win_width', 50);
         opts.CFn_opts.mi = struct('nbins', 150, 'win_width', 8, 'nthreads', opts.nthreads, ...
            'log_lookup', false, 'log_thresh', 1/(150*150*10));
      case{'mi', 'nmi'}
         opts.CFn_opts = struct('nbins', 150, 'win_width', 8, 'nthreads', opts.nthreads, ...
            'log_lookup', false, 'log_thresh', 1/(150*150*10) );
      case{'sd'}
         opts.CFn_opts = struct('hybrid_mask_op', true);
   end
end



% Set parameters
if opts.verbose
   t_corr = tic();
end

[O_trans, spacing, opts.MaxRef] = Make_Initial_Grid(opts, size(Iepi));

Iepi_mask = opts.MaskEPI; 
opts.MaskEPI = []; % free memory
Iepi_mask_less_CSF = opts.MaskEPILessCSF;
opts.MaskEPILessCSF = [];
Istruct_mask = opts.MaskStruct;
opts.MaskStruct = [];

% Rigid setup
rigid_param = [0 0 0 0 0 0];
if length(opts.rigid_scl) == 1
   rigid_scl = ones(size(rigid_param)) * opts.rigid_scl;
elseif length(opts.rigid_scl) == 6
   rigid_scl = opts.rigid_scl;
else
   error('opts.rigid_scl should have length of either 1 or 6.')
end


% EPI correction
if strcmpi(opts.similarity, 'inversion-epi')
   [O_trans, spacing, rigid_param] = EPI_correct_INVERSION_EPI(O_trans, spacing, rigid_param, rigid_scl, opts,...
      IepiBFC, Istruct, res, Iepi_mask, Iepi_mask_less_CSF, Istruct_mask);
   
elseif strcmpi(opts.similarity, 'inversion-t1')
   [O_trans, spacing, rigid_param] = EPI_correct_INVERSION_T1(O_trans, spacing, rigid_param, rigid_scl, opts,...
      IepiBFC, Istruct, res, Iepi_mask, Iepi_mask_less_CSF, Istruct_mask);
   
elseif ismember(opts.similarity, {'inversion-epi-t1', 'bdp'})
   [O_trans, spacing, rigid_param] = EPI_correct_INVERSION_EPI_T1(O_trans, spacing, rigid_param, rigid_scl, opts,...
      IepiBFC, Istruct, res, Iepi_mask, Iepi_mask_less_CSF, Istruct_mask);
else
   error('\nInvalid similarity or not implemented yet: %s\n', opts.similarity)
end

% Mutual information based refinement
% if strcmpi(opts.similarity, 'bdp')
%    par = rigid_param;
%    M = par2affineMat(par(1:3),par(4:6));
%    Ist_rigid = affine_transform_3dvol_double(double(Istruct), double(M), 1, opts.nthreads);
%    mask_Ist_rigid = affine_transform_3dvol_double(double(Istruct_mask),double(M), 1, opts.nthreads);
%    
%    [O_trans, spacing] = EPI_correct_NMI_refine(O_trans, spacing, opts, Iepi, Ist_rigid, res, ...
%       Iepi_mask, mask_Ist_rigid);   
% end


% Transform the input image with estimated transformation
t_x  = bspline_repeated_endpoint_deformation_3d_double_only_x(O_trans, size(Iepi), spacing, opts.nthreads);
[x_g, y_g, z_g] = ndgrid(1:size(Iepi,1), 1:size(Iepi,2), 1:size(Iepi,3));
x_g = x_g + t_x;
Iepi_corr = interpn(Iepi, x_g, y_g, z_g, 'cubic', 0);
if opts.intensity_correct
   [~, crct_grd] = gradient(t_x);
   Iepi_corr = Iepi_corr.*(1+crct_grd);
end

% Apply rigid transform
M = par2affineMat(rigid_param(1:3), rigid_param(4:6));
Istruct_reg = affine_transform_3dvol_double(double(Istruct),double(M),double(3), opts.nthreads);

if opts.verbose
   disp(['Completed EPI_correct_volume_bspline_registration() in ' num2str(toc(t_corr)) ' seconds.']); 
   drawnow();
end
end

function [O_trans, spacing, MaxRef] = Make_Initial_Grid(opts, Imgsize)
% Computes spacing, when not defined, and max number of refine steps before generating bspline
% control point grid.

min_voxel_spacing_bw_control_points = 2;
if isempty(opts.Spacing)
   MaxRef = min(floor(log2(Imgsize/4)));
   spacing = [2^MaxRef 2^MaxRef 2^MaxRef];
   
else
   spacing = round(opts.Spacing);
   t = spacing;
   MaxRef = 0;
   while (nnz(mod(t,2))==0) && (nnz(t<min_voxel_spacing_bw_control_points)==0)
      MaxRef = MaxRef+1;
      t = t/2;
   end
end
O_trans = bspline_grid_generate_deformation(spacing, Imgsize);

if MaxRef<1
   error('Minimum spacing should be 4 voxels, so that MaxRef>=1')
elseif MaxRef>2
   MaxRef = 2;
end

% Limit refinements steps to user input
if ~isempty(opts.MaxRef) && opts.MaxRef<MaxRef
   MaxRef = opts.MaxRef;
end

end

function [O_trans, Spacing] = EPI_correct_NMI_refine(O_trans, Spacing, opts, Iepi, Istruct, im_res, ...
   mask_Iepi, mask_Istruct)

msk_thresh = 0.45;
upsample_factor = 2; % along PED; should be pow of 2

if opts.verbose, disp('Starting MI-based correction...'); drawnow; end

reg_opts = struct(...
   'penalty', opts.Penalty, ...
   'verbose', false, ...
   'intensity_correct', opts.intensity_correct, ...
   'moving_mask', opts.moving_mask);

fmin_opts = struct(...
   'GradObj','on', ...
   'Display','off', ...
   'MaxFunEvals', 30, ...
   'TolX',1e-2, ...
   'TolFun',1e-9, ...
   'step_size', 15, ...
   'step_size_scale_iter', opts.step_size_scale_iter);
if opts.verbose, fmin_opts.Display='iter'; end

scale_res   = [4  im_res]; % in mm; should be >= im_res
resize_scl  = pow2(nextpow2(im_res./scale_res)); % must be power of 2
sgma = scale_res./(im_res*sqrt(8*log(2)));
resize_scl(resize_scl>1 | sgma<0.5) = 1; % donot upsample for super low-res images

for iter = 1:length(scale_res)
   fprintf('\nRefine: Iteration %d of %d running...', iter, length(scale_res));
   
   % resize/smooth images
   Spacing_iter = round(Spacing * resize_scl(iter)); % rounding should not be required? everything is pow of 2
   Spacing_iter(1) = Spacing_iter(1) * upsample_factor; % along PED
   
   temp = size(Iepi)*resize_scl(iter);
   temp(1) = temp(1) * upsample_factor; % along PED
   tsize = round(temp./Spacing_iter).*Spacing_iter + 1;
   scl_frc = tsize./size(Iepi);
   ISres = im_res./scl_frc;
   
   [ISepi, maskIepi, ISstruct, maskIstatic] = ...
      downsampleImagesMasks(Iepi, mask_Iepi, Istruct, mask_Istruct, sgma(iter), tsize, resize_scl(iter));
   
   O_trans = O_trans * resize_scl(iter) * upsample_factor; % along PED
   
   if opts.verbose
      disp(['Current Grid size : ' num2str(size(O_trans,1)) 'x' num2str(size(O_trans,2)) 'x' num2str(size(O_trans,3)) ]); drawnow;
   end
   
   [reg_opts.grad_img_indx, reg_opts.grad_mask_count, reg_opts.grad_mask_index, reg_opts.grad_ind_sz] = ...
      create_mask_index(size(O_trans), Spacing_iter, size(ISepi));
   
   % computer bspline operators
   [reg_opts.cp_mask_indx, reg_opts.cp_col_index] = ...
      bspline_create_mask_index(size(O_trans), Spacing_small, size(ISepi));
   %[reg_opts.penalty_op.V, reg_opts.penalty_op.cp_ind, reg_opts.penalty_op.Q] = ...
   %   bspline_penalty_repeated_endpoint_metal_bending_operator(sizes1, Spacing_small, ISres);
   reg_opts.penalty_op.resize_scl = resize_scl(iter); % to fix for control-point scaling hack!
   
   
   
   % Reshape O_trans from a matrix to a vector.
   O_trans1 = squeeze(O_trans(:,:,:,1));
   sizes1 = size(O_trans1);
   O_trans1 = O_trans1(:);
   
   % Start the b-spline nonrigid registration optimizer
   O_trans1 = gradient_descent_adjust_step(@(x)bspline_analytic_grad_MI_EPI(x, sizes1, O_trans, Spacing_iter, ISepi, ISstruct, ...
      ISres, maskIepi>msk_thresh, maskIstatic>msk_thresh, reg_opts), O_trans1, fmin_opts);
   
   % Reshape O_trans from a vector to a matrix
   O_trans1 = reshape(O_trans1, sizes1);
   O_trans(:,:,:,1) = O_trans1;

   % save current Control points
   if opts.debug
      t_x = bspline_phantom_endpoint_deformation_only_x(O_trans, size(ISepi), Spacing_iter);
      [x_g, y_g, z_g] = ndgrid(1:size(ISepi,1), 1:size(ISepi,2), 1:size(ISepi,3));
      x_g = x_g + t_x;      
      Ireg = interpn(ISepi, x_g, y_g, z_g, 'linear', 0);      
      if opts.intensity_correct
         [~, crct_grd] = gradient(t_x);
         Ireg = Ireg.*(1+crct_grd);
      end      
      overlay_volumes([ISstruct; ISstruct; ISepi], [ISepi; Ireg; Ireg], [0 1]);
      overlay_control_points2png(ISepi, O_trans, Spacing_iter, ['bspline-knots' num2str(iter)])
   end
   
   % control points in original size   
   O_trans(:,:,:,1) = O_trans(:,:,:,1)/upsample_factor;
   O_trans = O_trans./resize_scl(iter);
   
   
end
end

function [O_trans,Spacing] = nonrigid_registration_MI(O_trans, Spacing, options, Imoving, Istatic, res, ...
   MASKmoving, MASKstatic, similarity, MaxItt, x_rigid)
% Non-rigid b-spline grid registration
if(options.verbose>0), disp('Start non-rigid EPI_correct_bspline_registrationgid b-spline grid registration'); drawnow; end

% set registration options.
reg_options = struct('similarity', similarity, ...
   'penaltypercentage', options.Penalty, ...
   'interpolation', options.Interpolation, ...
   'scaling', options.Scaling, ...
   'verbose', false, ...
   'step', 0.1, ... % step for numerical gradient
   'intensity_correct', options.intensity_correct, ...
   'moving_mask', options.moving_mask);

if(strcmpi(similarity,'sd')), reg_options.centralgrad=false; end

% set optimization options (function minimization computation)
fmin_options = struct('GradObj','on', ...
   'GoalsExactAchieve',0, ...
   'StoreN',10, ...
   'HessUpdate','lbfgs',...%,'steepdesc' ...
   'Display','off', ...
   'MaxIter',7, ...
   'DiffMinChange',0.001, ...
   'DiffMaxChange',1, ...
   'MaxFunEvals',options.MaxFunEvals, ...
   'TolX',1e-2, ...
   'TolFun',1e-14, ...
   'GradConstr', true, ...
   'step_size', options.step_size, ...
   'step_size_scale_iter', options.step_size_scale_iter);

if(options.verbose>0), fmin_options.Display='iter'; end

if(options.dense_grid)
   refine_flag = [ 0      1     0];
   blur_flag   = [ 1      1     0];
   blur_param  = [2.5     1     1];
   resize_scl  = [0.5    0.5    1];
else
   refine_flag = [ 0     0     1     0];
   blur_flag   = [ 1     1     1     0];
   blur_param  = [3     2.4   1.7    1];
   resize_scl  = [0.5    1    1     1];
end

for iter = 1:length(blur_flag)
   if blur_flag(iter) == 1
      Hsize = round(2*blur_param(iter))+1;
      ISmoving = imgaussian(Imoving,Hsize/5,[Hsize Hsize Hsize]);
      ISstatic = imgaussian(Istatic,Hsize/5,[Hsize Hsize Hsize]);
      fmin_options.step_size = options.step_size / (1.8^(iter-1));
   else
      ISmoving = Imoving; %imgaussian(Imoving, 0.5, 0.3*[1 1 1]);
      ISstatic = Istatic; %imgaussian(Istatic, 0.5, 0.3*[1 1 1]);
      fmin_options.step_size = 1e4;
      fmin_options.MaxFunEvals = 3;
   end
   
   if refine_flag(iter) == 1
      [O_trans,Spacing] = bspline_refine_grid(O_trans, Spacing, size(Imoving));
      %reg_options.penaltypercentage = reg_options.penaltypercentage*10;
   end
   
   if resize_scl(iter) ~= 1.0
      MASKmovingsmall = my_imresize3d(MASKmoving, resize_scl(iter), [], 'nearest')>0;
      MASKstaticsmall = my_imresize3d(MASKstatic, resize_scl(iter), [], 'nearest')>0;
      ISmoving_small = my_imresize3d(ISmoving, resize_scl(iter), [],'linear');
      ISstatic_small = my_imresize3d(ISstatic, resize_scl(iter), [],'linear');
   else
      MASKmovingsmall = MASKmoving;
      MASKstaticsmall = MASKstatic;
      ISmoving_small = ISmoving;
      ISstatic_small = ISstatic;
   end
   
   Spacing_small = Spacing*resize_scl(iter);
   O_trans = O_trans*resize_scl(iter);
   ISres = res/resize_scl(iter);
   
   if (options.verbose>0),
      disp(['Current Grid size : ' num2str(size(O_trans,1)) 'x' num2str(size(O_trans,2)) 'x' num2str(size(O_trans,3)) ]); drawnow;
   end
   
   
   % Reshape O_trans from a matrix to a vector.
   O_trans1 = squeeze(O_trans(:,:,:,1));
   sizes1 = size(O_trans1);
   O_trans1 = O_trans1(:);
   
   % Start the b-spline nonrigid registration optimizer
   %    O_trans1 = fminlbfgs(@(x)EPI_correct_bspline_registration_gradient(x, sizes1, O_trans, Spacing_small, ISmoving_small, ISstatic_small, ISres, ...
   %                                                             MASKmovingsmall, MASKstaticsmall, reg_options), O_trans1, fmin_options);
   
   O_trans1 = gradient_descent(@(x)EPI_correct_bspline_registration_gradient_MI(x, sizes1, O_trans, Spacing_small, ISmoving_small, ISstatic_small, ISres, ...
      MASKmovingsmall, MASKstaticsmall, reg_options), O_trans1, fmin_options);
   
   % Reshape O_trans from a vector to a matrix
   O_trans1 = reshape(O_trans1, sizes1);
   O_trans(:,:,:,1) = O_trans1;
   
   
   
   % save current control points
   if options.debug
      [t_x, t_y, t_z] = bspline_phantom_endpoint_deformation(O_trans, size(ISmoving_small), Spacing_small);
      clear t_y t_z
      [x_g, y_g, z_g] = ndgrid(1:size(ISmoving_small,1), 1:size(ISmoving_small,2), 1:size(ISmoving_small,3));
      x_g = x_g + t_x;
      
      Ireg = interpn(ISmoving_small, x_g, y_g, z_g, 'linear', 0);
      
      if options.intensity_correct
         [~, crct_grd] = gradient(t_x,  ISres(1), ISres(2), ISres(3));
         Ireg = Ireg.*(1+crct_grd);
      end
      
      overlay_volumes2png(ISstatic_small, Ireg, [0 1], ['bspline' num2str(iter)]);
      clear crct_trsfm x_g y_g z_g t_x t_y t_z
      overlay_control_points2png(ISmoving_small, O_trans, Spacing_small, ['bspline-knots' num2str(iter)])
   end
   
   % control points in original size
   O_trans = O_trans/resize_scl(iter);
   
   fprintf('%d/%d iterations completed.\n', iter, length(blur_flag));
end
end

function [O_trans, Spacing, rigid_param] = EPI_correct_INVERSION_EPI(O_trans, Spacing, rigid_param, rigid_scl, ...
   opts, Iepi, Istruct, im_res, mask_Iepi, mask_Iepi_less_CSF, mask_Istruct)
% inverts EPI image 

fprintf('\nPerforming non-rigid b-spline registration with INVERSION-EPI');
msk_thrs = 0.35; % threshold to convert contineous mask to binary mask
Iepi = normalize_intensity(Iepi, [0 99.9], mask_Iepi>msk_thrs);

reg_opts = struct( ...
   'penalty', opts.Penalty, ...
   'debug', opts.debug, ...
   'intensity_correct', opts.intensity_correct, ...
   'mask_type', 'or', ...
   'histeqInterpolant', opts.histeqInterpolant, ...
   'histeqGradInterpolant', opts.histeqGradInterpolant, ...
   'nthreads', opts.nthreads ...
   );

fmin_opts = struct('GradObj','on', ...
   'Display','off', ...
   'MaxFunEvals',opts.MaxFunEvals, ...
   'TolX',1e-6, ...
   'TolFun',1e-10, ...
   'backtrack_c', 1e-7, ...
   'step_size', opts.step_size, ...
   'step_size_scale_iter', opts.step_size_scale_iter);
if(opts.verbose>0), fmin_opts.Display='iter'; end


CP_refine   = [0    1    0    0    1    0      0  ]; % bspline CP refine flag
scale_res   = [4.5  4.5  4    3    3.5  2.5 im_res]; % in mm; should be >= im_res
resize_scl  = [0.5  0.5  1    1    1    1      1  ]; % must be power of 2

penalty_wt  = [50   40   35   25   30   25    20  ]*25/2; % for CP-diff. 

rigid_refine= [1    1    1    1    1    0      0  ] * (max(rigid_scl)>0); % rigid refine flag
rigid_wt    = [0.2  0.2  0.2  0.2  0.2  0.2    1  ]; % rigid parameter wt, must be >0
sgma = scale_res./(im_res*sqrt(8*log(2)));

for iter = 1:length(scale_res)
   fprintf('\nScale %d of %d running...', iter, length(scale_res));
   
   % update options
   reg_opts.penalty(1) = penalty_wt(iter)*opts.Penalty(1);
   reg_opts.rigid_refine = rigid_refine(iter);
   rigid_param_iter = rigid_param/rigid_wt(iter);
   rigid_scl_iter = rigid_scl*rigid_wt(iter);
   
   if CP_refine(iter) == 1
      [O_trans, Spacing] = bspline_refine_grid_deformation(O_trans, Spacing, size(Iepi));
   end
   
   if opts.verbose>0
      disp(['Current Grid size : ' num2str(size(O_trans,1)) 'x' num2str(size(O_trans,2)) 'x' num2str(size(O_trans,3)) ]); drawnow;
   end
   
   % rounding should not be required? everything is pow of 2
   Spacing_small = round(Spacing * resize_scl(iter)); 
   ISres = im_res.*(Spacing./Spacing_small);
   
   % resample image
   if resize_scl(iter) < 1
      tsize = round(size(Iepi)*resize_scl(iter)./Spacing_small).*Spacing_small + 1;
      scl_frc = tsize./size(Iepi);
      O_trans = O_trans * resize_scl(iter);
      ISres = im_res ./ scl_frc;
      rigid_param_iter(1:3) = rigid_param_iter(1:3) .* scl_frc;
   else
      tsize = size(Iepi);
   end
   
   [ISepi, mask_Iepi_small, ISstruct, mask_Istruct_small] = ...
      downsampleImagesMasks(Iepi, mask_Iepi, Istruct, mask_Istruct, sgma(iter), tsize, resize_scl(iter));
   
   [~, mask_Iepi_less_CSF_small, ~, ~] = ...
      downsampleImagesMasks(Iepi, mask_Iepi_less_CSF, Istruct, mask_Istruct, sgma(iter), tsize, resize_scl(iter));
   
   
   % computer bspline operators
   [reg_opts.cp_mask_indx, reg_opts.cp_col_index] = ...
      bspline_create_mask_index(size(O_trans), Spacing_small, size(ISepi));
   
   % INVERSION intensity map
   [~, reg_opts.histeqInterpolant, reg_opts.histeqGradInterpolant] = unwarpAndInvertImg(...
      ISepi, mask_Iepi_less_CSF_small, ISstruct, mask_Istruct_small, ...
      'inversion-epi', O_trans, Spacing_small, msk_thrs, opts);   
   
   % optimize
   O_trans1 = squeeze(O_trans(:,:,:,1));
   sizes1 = size(O_trans1);
   X = [rigid_param_iter(:); O_trans1(:)];
   X = gradient_descent_backtrack(@(x)bspline_analytic_grad_INVERSION_EPI(x, sizes1, ...
      O_trans, Spacing_small, rigid_scl_iter, ISepi, ISstruct, ISres, ...
      mask_Iepi_small, mask_Iepi_less_CSF_small, mask_Istruct_small, reg_opts), X, fmin_opts);
   rigid_param = X(1:6)' * rigid_wt(iter); % Multiplied by rigid_wt fix scaling
   O_trans1 = X(7:end);
   
   % Reshape O_trans from a vector to a matrix
   O_trans1 = reshape(O_trans1, sizes1);
   O_trans(:,:,:,1) = O_trans1;
   
   if opts.debug
      t_x = bspline_repeated_endpoint_deformation_3d_double_only_x(O_trans(:,:,:,1), size(ISepi), Spacing_small, opts.nthreads);
      [x_g, y_g, z_g] = ndgrid(1:size(ISepi,1), 1:size(ISepi,2), 1:size(ISepi,3));
      x_g = x_g + t_x;
      
      Iepi_corr = interpn(ISepi, x_g, y_g, z_g, 'linear', 0);
      mask_Iepi_corr = interpn(mask_Iepi_small, x_g, y_g, z_g, 'linear', 0);
      % opts.intensity_correct
      [~, crct_grd] = gradient(t_x);
      Iepi_corr = Iepi_corr.*(1+crct_grd);
      
      Iepi_corr_INV = reshape(reg_opts.histeqInterpolant(1-Iepi_corr(:)), size(ISepi)).*mask_Iepi_corr;
      Iepi_INV = reshape(reg_opts.histeqInterpolant(1-ISepi(:)), size(ISepi)).*mask_Iepi_small;
      
      par = rigid_param.*rigid_scl;
      M = par2affineMat(par(1:3),par(4:6));
      Istruct_reg = demon_registration.affine_transform_3d_double(double(ISstruct),double(M),double(1));
      mask_Istruct_reg = demon_registration.affine_transform_3d_double(double(mask_Istruct_small),double(M),double(1));
      
      display_volume([ISstruct Istruct_reg; ISepi Iepi_corr; Iepi_INV Iepi_corr_INV])
      
      err_orig = (Istruct_reg-Iepi_INV).^2;
      err_reg = (Istruct_reg-Iepi_corr_INV).^2;
      display_volume([err_orig; err_reg]); colormap jet; caxis([0 max([err_orig(:); err_reg(:)])])
      
      overlay_volumes([Istruct_reg.*mask_Istruct_reg, Istruct_reg; Istruct_reg.*mask_Istruct_reg, Istruct_reg; Iepi_INV, ISepi], ...
         [Iepi_corr_INV, Iepi_corr; Iepi_INV, ISepi; Iepi_corr_INV, Iepi_corr], [0 1]);
      
      overlay_volumes([Istruct_reg; Istruct_reg; mask_Istruct_reg; mask_Iepi_small], ...
         [ISstruct; mask_Iepi_corr*0.5; mask_Iepi_corr; mask_Iepi_corr], [0 1]);
      
      overlay_control_points2png(ISepi, O_trans, Spacing_small, ['bspline-knots' num2str(iter)])
   end
   
   % control points in original size
   if resize_scl(iter) ~= 1.0
      O_trans = O_trans./resize_scl(iter);
      rigid_param(1:3) = rigid_param(1:3)./scl_frc;
   end
   
   if opts.debug
      t_x = bspline_repeated_endpoint_deformation_3d_double_only_x(O_trans(:,:,:,1), size(Iepi), Spacing, opts.nthreads);
      [x_g, y_g, z_g] = ndgrid(1:size(Iepi,1), 1:size(Iepi,2), 1:size(Iepi,3));
      x_g = x_g + t_x;
      Iepi_corr = interpn(Iepi, x_g, y_g, z_g, 'linear', 0);
      % opts.intensity_correct
      [~, crct_grd] = gradient(t_x);
      Iepi_corr = Iepi_corr.*(1+crct_grd);
      
      par = rigid_param.*rigid_scl;
      M = par2affineMat(par(1:3),par(4:6));
      Istruct_reg = demon_registration.affine_transform_3d_double(double(Istruct),double(M), 1);
      overlay_volumes([Istruct_reg, Istruct, Iepi_corr], [Iepi_corr,Iepi,Iepi], [0 1]);
   end
end

% final rigid parameters
rigid_param = rigid_param.*rigid_scl;

end


function [O_trans, Spacing, rigid_param] = EPI_correct_INVERSION_T1(O_trans, Spacing, rigid_param, rigid_scl, ...
   opts, Iepi, Istruct, im_res, mask_Iepi, mask_Iepi_less_CSF, mask_Istruct)
% inverts T1 image 

fprintf('\nPerforming non-rigid b-spline registration with INVERSION-T1');
msk_thrs = 0.35; % threshold to convert contineous mask to binary mask

Iepi = normalize_intensity(Iepi, [0 99.9], mask_Iepi>msk_thrs);

reg_opts = struct( ...
   'penalty', opts.Penalty, ...
   'debug', opts.debug, ...
   'intensity_correct', opts.intensity_correct, ...
   'mask_type', 'or', ...
   'histeqInterpolant', opts.histeqInterpolant, ...
   'histeqGradInterpolant', opts.histeqGradInterpolant, ...
   'nthreads', opts.nthreads ...
   );

% set optimization options (function minimization computation)
fmin_opts = struct('GradObj','on', ...
   'Display','off', ...
   'MaxFunEvals',opts.MaxFunEvals, ...
   'TolX',1e-6, ...
   'TolFun',1e-10, ...
   'step_size', opts.step_size, ...
   'backtrack_c', 1e-7, ...
   'step_size_scale_iter', opts.step_size_scale_iter);
if(opts.verbose>0), fmin_opts.Display='iter'; end

CP_refine   = [0    1    0    0    1    0      0  ]; % bspline CP refine flag
scale_res   = [4.5  4.5  4    3    3.5  2.5 im_res]; % in mm; should be >= im_res
resize_scl  = [0.5  0.5  1    1    1    1      1  ]; % must be power of 2

penalty_wt  = [50   40   35   25   30   25    20  ]*25/2; % for CP-diff

rigid_refine= [1    1    1    1    1    0      0  ] * (max(rigid_scl)>0); % rigid refine flag
rigid_wt    = [0.2  0.2  0.2  0.2  0.2  0.2    1  ]; % rigid parameter wt, must be >0
sgma = scale_res./(im_res*sqrt(8*log(2)));

for iter = 1:length(scale_res)
   fprintf('\nScale %d of %d running...', iter, length(scale_res));
   
   % update options
   reg_opts.penalty(1) = penalty_wt(iter)*opts.Penalty(1);
   reg_opts.rigid_refine = rigid_refine(iter);
   rigid_param_iter = rigid_param/rigid_wt(iter);
   rigid_scl_iter = rigid_scl*rigid_wt(iter);
   
   if CP_refine(iter) == 1
      [O_trans, Spacing] = bspline_refine_grid_deformation(O_trans, Spacing, size(Iepi));
   end
   
   if opts.verbose>0
      disp(['Current Grid size : ' num2str(size(O_trans,1)) 'x' num2str(size(O_trans,2)) 'x' num2str(size(O_trans,3)) ]); drawnow;
   end
   
   % invert T1 at highest resolution
   Istruct_INV = unwarpAndInvertImg(Iepi, mask_Iepi_less_CSF, Istruct, mask_Istruct, 'inversion-t1', ...
      O_trans, Spacing, msk_thrs, opts);
   
   % rounding should not be required? everything is pow of 2
   Spacing_small = round(Spacing * resize_scl(iter)); 
   ISres = im_res.*(Spacing./Spacing_small);

   % resample image
   if resize_scl(iter) < 1
      tsize = round(size(Iepi)*resize_scl(iter)./Spacing_small).*Spacing_small + 1;
      scl_frc = tsize./size(Iepi);
      O_trans = O_trans * resize_scl(iter);
      ISres = im_res ./ scl_frc;
      rigid_param_iter(1:3) = rigid_param_iter(1:3) .* scl_frc;
   else
      tsize = size(Iepi);
   end
   
   [ISepi, mask_Iepi_small, ISstruct, mask_Istruct_small] = ...
      downsampleImagesMasks(Iepi, mask_Iepi_less_CSF, Istruct_INV, mask_Istruct, sgma(iter), tsize, resize_scl(iter));
   
   % computer bspline operators
   [reg_opts.cp_mask_indx, reg_opts.cp_col_index] = ...
      bspline_create_mask_index(size(O_trans), Spacing_small, size(ISepi));
   
   % optimize
   O_trans1 = squeeze(O_trans(:,:,:,1));
   sizes1 = size(O_trans1);
   X = [rigid_param_iter(:); O_trans1(:)];
   X = gradient_descent_backtrack(@(x)bspline_analytic_grad_SSD_EPI(x, sizes1, ...
      O_trans, Spacing_small, rigid_scl_iter, ISepi, ISstruct, ISres, ...
      mask_Iepi_small, mask_Istruct_small, reg_opts), X, fmin_opts);
   rigid_param = X(1:6)' * rigid_wt(iter); % Multiplied by rigid_wt fix scaling 
   O_trans1 = X(7:end);
   
   % Reshape O_trans from a vector to a matrix
   O_trans1 = reshape(O_trans1, sizes1);
   O_trans(:,:,:,1) = O_trans1;
   
   if opts.debug
      t_x = bspline_repeated_endpoint_deformation_3d_double_only_x(O_trans(:,:,:,1), size(ISepi), Spacing_small, opts.nthreads);
      [x_g, y_g, z_g] = ndgrid(1:size(ISepi,1), 1:size(ISepi,2), 1:size(ISepi,3));
      x_g = x_g + t_x;
      
      Iepi_corr = interpn(ISepi, x_g, y_g, z_g, 'linear', 0);
      mask_Iepi_corr = interpn(mask_Iepi_small, x_g, y_g, z_g, 'linear', 0);
      % opts.intensity_correct
      [~, crct_grd] = gradient(t_x);
      Iepi_corr = Iepi_corr.*(1+crct_grd);
      
      par = rigid_param.*rigid_scl;
      M = par2affineMat(par(1:3),par(4:6));
      Istruct_reg = demon_registration.affine_transform_3d_double(double(ISstruct),double(M),double(1));
      mask_Istruct_reg = demon_registration.affine_transform_3d_double(double(mask_Istruct_small),double(M),double(1));
      
      temp = my_imresize3d(imgaussian(Istruct, sgma(iter)), [], tsize, 'linear');
      Istruct_orig_reg = demon_registration.affine_transform_3d_double(double(temp),double(M),double(1));
      
      display_volume([temp Istruct_orig_reg; ISstruct Istruct_reg; ISepi Iepi_corr])
      
      err_orig = (Istruct_reg-ISepi).^2;
      err_reg = (Istruct_reg-Iepi_corr).^2;
      display_volume([err_orig; err_reg]); colormap jet; caxis([0 max([err_orig(:); err_reg(:)])])
      
      overlay_volumes([Istruct_reg.*mask_Istruct_reg, Istruct_reg; Istruct_reg.*mask_Istruct_reg, Istruct_reg], ...
         [Iepi_corr, Iepi_corr; ISepi, ISepi], [0 0.5]);
      
      overlay_volumes([Istruct_reg; Istruct_reg; mask_Istruct_reg; mask_Iepi_small], ...
         [ISstruct; mask_Iepi_corr*0.5; mask_Iepi_corr; mask_Iepi_corr], [0 1]);
      
      overlay_control_points2png(ISepi, O_trans, Spacing_small, ['bspline-knots' num2str(iter)])
   end
   
   % control points in original size
   if resize_scl(iter) ~= 1.0
      O_trans = O_trans./resize_scl(iter);
      rigid_param(1:3) = rigid_param(1:3)./scl_frc;
   end
   
   if opts.debug
      t_x = bspline_repeated_endpoint_deformation_3d_double_only_x(O_trans(:,:,:,1), size(Iepi), Spacing, opts.nthreads);
      [x_g, y_g, z_g] = ndgrid(1:size(Iepi,1), 1:size(Iepi,2), 1:size(Iepi,3));
      x_g = x_g + t_x;
      Iepi_corr = interpn(Iepi, x_g, y_g, z_g, 'linear', 0);
      % opts.intensity_correct
      [~, crct_grd] = gradient(t_x);
      Iepi_corr = Iepi_corr.*(1+crct_grd);
      
      par = rigid_param.*rigid_scl;
      M = par2affineMat(par(1:3),par(4:6));
      Istruct_reg = demon_registration.affine_transform_3d_double(double(Istruct),double(M), 1);
      Istruct_reg_INV = demon_registration.affine_transform_3d_double(double(Istruct_INV),double(M), 1);
      overlay_volumes([Istruct_reg, Istruct_reg_INV, Istruct_INV], [1.8*Iepi_corr,Iepi_corr,Iepi], [0 1]);
   end
end

% final rigid parameters
rigid_param = rigid_param.*rigid_scl;

end


function [O_trans, Spacing, rigid_param] = EPI_correct_INVERSION_EPI_T1(O_trans, Spacing, rigid_param, rigid_scl, ...
   opts, Iepi, Istruct, im_res, mask_Iepi, mask_Iepi_less_CSF, mask_Istruct)
% inverts both EPI & T1 image 

fprintf('\nPerforming non-rigid b-spline registration with INVERSION-EPI-T1 (INVERSION-BOTH)');
msk_thrs = 0.4; % threshold to convert contineous mask to binary mask

Iepi = normalize_intensity(Iepi, [0 99.9], mask_Iepi>msk_thrs);

reg_opts = struct( ...
   'penalty', opts.Penalty, ...
   'debug', opts.debug, ...
   'intensity_correct', opts.intensity_correct, ...
   'mask_type', 'or', ...
   'histeqInterpolant', opts.histeqInterpolant, ...
   'histeqGradInterpolant', opts.histeqGradInterpolant, ...
   'nthreads', opts.nthreads ...
   );

fmin_opts = struct('GradObj','on', ...
   'Display','off', ...
   'MaxFunEvals',opts.MaxFunEvals, ...
   'TolX',1e-6, ...
   'TolFun',1e-10, ...
   'backtrack_c', 1e-7, ...
   'step_size', opts.step_size, ...
   'step_size_scale_iter', opts.step_size_scale_iter);
if opts.verbose, fmin_opts.Display='iter'; end

CP_refine   = [0    1    0    0    1    0      0  ]; % bspline CP refine flag
scale_res   = [4.5  4.5  4    3    3.5  2.5 im_res]; % in mm; should be >= im_res
resize_scl  = [0.5  0.5  1    1    1    1      1  ]; % must be power of 2

penalty_wt  = [50   40   35   25   30   25    20  ]*25; % for CP-diff

rigid_refine= [1    1    1    1    1    0      0  ] * (max(rigid_scl)>0); % rigid refine flag
rigid_wt    = [0.2  0.2  0.2  0.2  0.2  0.2    1  ]; % rigid parameter wt, must be >0
sgma = scale_res./(im_res*sqrt(8*log(2)));

for iter = 1:length(scale_res)
   fprintf('\nScale %d of %d running...', iter, length(scale_res));
   
   % update options
   reg_opts.penalty(1) = penalty_wt(iter)*opts.Penalty(1);
   reg_opts.rigid_refine = rigid_refine(iter);
   rigid_param_iter = rigid_param/rigid_wt(iter);
   rigid_scl_iter = rigid_scl*rigid_wt(iter);
   reg_opts.invEPI_weight = mean(Iepi(mask_Iepi>msk_thrs))/mean(Istruct(mask_Istruct>msk_thrs));
   
   
   if CP_refine(iter) == 1
      [O_trans, Spacing] = bspline_refine_grid_deformation(O_trans, Spacing, size(Iepi));
   end
   
   if opts.verbose>0
      disp([' Res: ' num2str(scale_res(iter)) 'mm, Grid size: ' num2str(size(O_trans,1)) ...
         'x' num2str(size(O_trans,2)) 'x' num2str(size(O_trans,3)) ]); drawnow();
   end
   
   % invert T1 at highest resolution
   Istruct_INV = unwarpAndInvertImg(Iepi, mask_Iepi_less_CSF, Istruct, mask_Istruct, 'inversion-t1', ...
      O_trans, Spacing, msk_thrs, opts);   
   
   % rounding should not be required? everything is pow of 2
   Spacing_small = round(Spacing * resize_scl(iter)); 
   ISres = im_res.*(Spacing./Spacing_small);
   
   % resample/resize image
   if resize_scl(iter) < 1
      tsize = round(size(Iepi)*resize_scl(iter)./Spacing_small).*Spacing_small + 1;
      scl_frc = tsize./size(Iepi);
      O_trans = O_trans * resize_scl(iter);
      ISres = im_res ./ scl_frc;
      rigid_param_iter(1:3) = rigid_param_iter(1:3) .* scl_frc;
   else
      tsize = size(Iepi);
   end
   [ISepi, mask_Iepi_small, ISstruct, mask_Istruct_small] = ...
      downsampleImagesMasks(Iepi, mask_Iepi, Istruct, mask_Istruct, sgma(iter), tsize, resize_scl(iter));
   
   [~, mask_Iepi_less_CSF_small, ISstructINV, ~] = ...
      downsampleImagesMasks(Iepi, mask_Iepi_less_CSF, Istruct_INV, mask_Istruct, sgma(iter), tsize, resize_scl(iter));
   
   O_trans1 = squeeze(O_trans(:,:,:,1));
   sizes1 = size(O_trans1);
   
   % computer bspline operators
   [reg_opts.cp_mask_indx, reg_opts.cp_col_index] = ...
      bspline_create_mask_index(size(O_trans), Spacing_small, size(ISepi));
   % [reg_opts.penalty_op.V, reg_opts.penalty_op.cp_ind, reg_opts.penalty_op.Q] = ...
   %   bspline_penalty_repeated_endpoint_metal_bending_operator(sizes1, Spacing_small, ISres);
   % reg_opts.penalty_op.resize_scl = resize_scl(iter); % to fix for control-point scaling hack!
   
   % INVERSION intensity map
   [~, reg_opts.histeqInterpolant, reg_opts.histeqGradInterpolant] = unwarpAndInvertImg(...
      ISepi, mask_Iepi_less_CSF_small, ISstruct, mask_Istruct_small, ...
      'inversion-epi', O_trans, Spacing_small, msk_thrs, opts);
   
   % optimize
   X = [rigid_param_iter(:); O_trans1(:)];
   X = gradient_descent_backtrack(@(x)bspline_analytic_grad_INVERSION_EPI_T1(x, sizes1, ...
      O_trans, Spacing_small, rigid_scl_iter, ISepi, ISstruct, ISstructINV, ISres, ...
      mask_Iepi_small, mask_Iepi_less_CSF_small, mask_Istruct_small, reg_opts), X, fmin_opts);
   rigid_param = X(1:6)' * rigid_wt(iter); % Multiplied by rigid_wt fix scaling 
   O_trans1 = X(7:end);
   
   % Reshape O_trans from a vector to a matrix
   O_trans1 = reshape(O_trans1, sizes1);
   O_trans(:,:,:,1) = O_trans1;
   
   if opts.debug
      par = rigid_param.*rigid_scl;
      M = par2affineMat(par(1:3),par(4:6));
      Istruct_reg = demon_registration.affine_transform_3d_double(double(ISstruct),double(M),double(1));
      mask_Istruct_reg = demon_registration.affine_transform_3d_double(double(mask_Istruct_small),double(M),double(1));
      
      Istruct_reg_INV = unwarpAndInvertImg(ISepi, mask_Iepi_less_CSF_small, Istruct_reg, mask_Istruct_reg, 'inversion-t1', ...
         O_trans, Spacing_small, msk_thrs, opts);
      ISstruct_INV = unwarpAndInvertImg(ISepi, mask_Iepi_less_CSF_small, ISstruct, mask_Istruct_small, 'inversion-t1', ...
         O_trans, Spacing_small, msk_thrs, opts);
      [Iepi_corr_INV, ~,~, Iepi_corr, mask_Iepi_corr] = unwarpAndInvertImg(ISepi, mask_Iepi_less_CSF_small, ISstruct, ...
         mask_Istruct_small, 'inversion-epi', O_trans, Spacing_small, msk_thrs, opts);
      Iepi_INV = unwarpAndInvertImg(ISepi, mask_Iepi_less_CSF_small, ISstruct, mask_Istruct_small, 'inversion-epi', ...
         zeros(size(O_trans)), Spacing_small, msk_thrs, opts);
      
      overlay_volumes([Iepi_INV; Iepi_corr_INV],[ISstruct; Istruct_reg])
      overlay_volumes([ISstruct_INV; Istruct_reg_INV],[ISepi; Iepi_corr])
      
      display_volume([ISstruct Istruct_reg; Iepi_INV Iepi_corr_INV; ISepi Iepi_corr; ISstruct_INV Istruct_reg_INV])
      t_x = normalize_intensity(bspline_repeated_endpoint_deformation_3d_double_only_x(...
         O_trans, size(ISepi), Spacing_small, opts.nthreads), [0.1 99.9], (mask_Iepi_corr.*mask_Istruct_reg)>msk_thrs);
      overlay_volumes([ISstruct; Istruct_reg; ISepi; t_x], [ISepi; Iepi_corr; Iepi_corr; t_x])
      
      msk_orig = (mask_Istruct_small + mask_Iepi_small)/2;
      msk_reg = (mask_Istruct_reg + mask_Iepi_corr)/2;
      err_orig1 = (ISstruct-Iepi_INV).^2 .* msk_orig;
      err_reg1 = (Istruct_reg-Iepi_corr_INV).^2 .* msk_reg;
      err_orig2 = (ISstruct_INV-ISepi).^2 .* msk_orig;
      err_reg2 = (Istruct_reg_INV-Iepi_corr).^2 .* msk_reg;
      display_volume([err_orig1, err_reg1; [err_orig1 err_reg1]*reg_opts.invEPI_weight; err_orig2 err_reg2]); 
      colormap jet; % caxis([0 max([err_orig1(:); err_orig2(:)])])
      
      overlay_volumes([Istruct_reg.*mask_Istruct_reg, Istruct_reg; Istruct_reg.*mask_Istruct_reg, Istruct_reg; Iepi_INV, ISepi], ...
         [Iepi_corr_INV, Iepi_corr; Iepi_INV, ISepi; Iepi_corr_INV, Iepi_corr], [0 1]);      
      overlay_volumes([Istruct_reg; Istruct_reg; mask_Istruct_reg; mask_Iepi_small], ...
         [ISstruct; mask_Iepi_corr*0.5; mask_Iepi_corr; mask_Iepi_corr], [0 1]);
      overlay_control_points2png(ISepi, O_trans, Spacing_small, ['bspline-knots' num2str(iter)])
   end
   
   % control points in original size
   if resize_scl(iter) ~= 1.0
      O_trans = O_trans./resize_scl(iter);
      rigid_param(1:3) = rigid_param(1:3)./scl_frc;
   end
   
   if opts.debug
      t_x = bspline_repeated_endpoint_deformation_3d_double_only_x(O_trans(:,:,:,1), size(Iepi), Spacing, opts.nthreads);
      [x_g, y_g, z_g] = ndgrid(1:size(Iepi,1), 1:size(Iepi,2), 1:size(Iepi,3));
      x_g = x_g + t_x;
      Iepi_corr = interpn(Iepi, x_g, y_g, z_g, 'linear', 0);
      % opts.intensity_correct
      [~, crct_grd] = gradient(t_x);
      Iepi_corr = Iepi_corr.*(1+crct_grd);
      
      par = rigid_param.*rigid_scl;
      M = par2affineMat(par(1:3),par(4:6));
      Istruct_reg = demon_registration.affine_transform_3d_double(double(Istruct),double(M), 1);
      overlay_volumes([Istruct_reg, Istruct, Iepi_corr/reg_opts.invEPI_weight/1.5], [Iepi_corr,Iepi,Iepi]/reg_opts.invEPI_weight/1.5, [0 1]);
   end
end

% final rigid parameters
rigid_param = rigid_param.*rigid_scl;

end

function [ISepi, mask_Iepi_small, ISstruct, mask_Istruct_small] = ...
   downsampleImagesMasks(Iepi, mask_Iepi, Istruct, mask_Istruct, sgma, tsize, resize_scl)
% downsamples image and masks such that erosion and blurring on top and bottom EPI slice are
% slightly different (avoid blurring from background). Assumes second dimension of image is readout
% direction.

se = strel_sphere(ceil(1./resize_scl/2));
msk_thrs = 0.2; % threshold to convert contineous mask to binary mask
sz = size(Iepi);

tempM = zeros(sz);
temp = zeros(sz);
M = any(mask_Iepi>msk_thrs, 1);
[m1, m2, m3] = find_bounding_box(M);
tempM(:, m2, m3) = imerode(mask_Iepi(:, m2, m3), se);
temp(:, m2, m3) = imgaussian(Iepi(:, m2, m3), sgma);
if ~isequal(sz, tsize) % when downsampling repeat end slices
   l = max(1, m2(1)-1); h = min(m2(end)+1, sz(2));
   temp(:, l, :) = temp(:, m2(1), :);
   temp(:, h, :) = temp(:, m2(end), :);
   tempM(:, l, :) = tempM(:, m2(1), :);
   tempM(:, h, :) = tempM(:, m2(end), :);
   
   l = max(1, m3(1)-1); h = min(m3(end)+1, sz(3));
   temp(:, :, l) = temp(:, :, m3(1));
   temp(:, :, h) = temp(:, :, m3(end));
   tempM(:, :, l) = tempM(:, :, m3(1));
   tempM(:, :, h) = tempM(:, :, m3(end));
end
ISepi = my_imresize3d(temp, [], tsize, 'linear');
mask_Iepi_small = my_imresize3d(tempM, [], tsize, 'linear');

tempM = zeros(size(Istruct));
temp = zeros(size(Istruct));
M = any(mask_Istruct>msk_thrs, 1);
[m1, m2, m3] = find_bounding_box(M);
tempM(:, m2, m3) = imerode(mask_Istruct(:, m2, m3), se);

temp(:, m2, m3) = imgaussian(Istruct(:, m2, m3), sgma);
if ~isequal(sz, tsize) % when downsampling repeat end slices
   l = max(1, m2(1)-1); h = min(m2(end)+1, sz(2));
   temp(:, l, :) = temp(:, m2(1), :);
   temp(:, h, :) = temp(:, m2(end), :);
   tempM(:, l, :) = tempM(:, m2(1), :);
   tempM(:, h, :) = tempM(:, m2(end), :);
   
   l = max(1, m3(1)-1); h = min(m3(end)+1, sz(3));
   temp(:, :, l) = temp(:, :, m3(1));
   temp(:, :, h) = temp(:, :, m3(end));
   tempM(:, :, l) = tempM(:, :, m3(1));
   tempM(:, :, h) = tempM(:, :, m3(end));
end
ISstruct = my_imresize3d(temp, [], tsize, 'linear');
mask_Istruct_small = my_imresize3d(tempM, [], tsize, 'linear');

end


function [I_inv, histeqInterpolant, histeqGradInterpolant, Iepi_corr, mask_Iepi_corr] = unwarpAndInvertImg(Iepi, msk_Iepi, Ist, msk_Ist, ...
   inv_str, O_trans, Spacing, msk_thrs, opts)
% unwarps epi using currect distortion estimate, compute intensity map and invert images

% unwarp
t_x = bspline_repeated_endpoint_deformation_3d_double_only_x(O_trans, size(Iepi), Spacing, opts.nthreads);
[x_g, y_g, z_g] = ndgrid(1:size(Iepi,1), 1:size(Iepi,2), 1:size(Iepi,3));
x_g = x_g + t_x;
Iepi_corr = interpn(Iepi, x_g, y_g, z_g, 'linear', 0);
mask_Iepi_corr = interpn(msk_Iepi, x_g, y_g, z_g, 'linear', 0);

% opts.intensity_correct
[~, crct_grd] = gradient(t_x);
Iepi_corr = Iepi_corr.*(1+crct_grd);

Iepi_corr = clip_intensity(Iepi_corr);

% intensity map and inversion
if strcmpi(inv_str, 'inversion-epi')
   I_inv = zeros(size(Iepi));
   [I_inv(mask_Iepi_corr>msk_thrs), ~, histeqMap, T_grid] = histeqSmooth(1-Iepi_corr(mask_Iepi_corr>msk_thrs), ...
      Ist(msk_Ist>msk_thrs), opts.CFn_opts.nbins, opts.CFn_opts.win_width);
   
elseif strcmpi(inv_str, 'inversion-t1')
   I_inv = zeros(size(Ist));
   [I_inv(msk_Ist>msk_thrs), ~, histeqMap, T_grid] = histeqSmooth(1-Ist(msk_Ist>msk_thrs), ...
      Iepi_corr(mask_Iepi_corr>msk_thrs), opts.CFn_opts.nbins, opts.CFn_opts.win_width);
end
histeqInterpolant = griddedInterpolant(T_grid, histeqMap, 'cubic');
histeqGradInterpolant = griddedInterpolant(T_grid, gradient(smoothGaussian1D(histeqMap,2), T_grid(2)-T_grid(1)), 'cubic');

end

function img = clip_intensity(img)
img(img<0) = 0;
img(img>1) = 1;
end

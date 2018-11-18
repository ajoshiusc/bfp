function [M_world, ref_loc, x_param] = register_files_affine(moving_filename, static_filename, output_filename, opts)
% Affine registration inteface for nifti files. This function uses nifti headers information. 
%
% Outputs: 
%   M_world - Estimated 4x4 affine transformation matrix for the registration
%   ref_loc - Coordinate location (3x1 vector) which was used as origin for affine transformation.
%             See affine_transform_nii.m for more details.
%   x_param - Transformation parameters
%
%   opts.similarity:
%          'sd' : Sum of Squared difference (same modality)
%          'mi' : Normalized mutual information (Studholme)
%          'cr' : Correlation ratio (Roche et al.)
%          'inversion' : INVERSION. Static image is inverted. 
%          'bdp' : INVERSION followed by NMI
%
%   opts.dof - Specify degree of freedom/number of parameters
%          1  : Only translation along one axis, see notes below
%          2  : Only translation along two axis, see notes below
%          3  : Only translation along three axis, see notes below
%          6  : Rigid, Translation and rotation 
%          12 : Full Affine, Translation, rotation, shear 
%
%      When dof<=3, it is implicitly assumed that translation axis matches the image acquition
%      axis (& NOT world coordinate axis). For dof<=3, resolution and orientation information from
%      header is used. (orientation information is used to re-orient volume, in orthogonal
%      sense). Translation axis is defined by vector opts.axis_translation of appropriate
%      length (length of 2 for dof of 2 & likewise; no need to define for dof=3).
%


% Check/set input options
defaultoptions = struct( ...
   'similarity', '', ...  See comments
   'dof', 6, ...          See comments
   'axis_translation', [],  ... 1/2/3, defines the translation axis for dof<3
   'intensity_norm', '',... histeq/sigmoid/none - intensity normalization prior to registration, when empty sets appropriately based on similarity
   'static_mask', [], ... filename for mask of <static_filename>
   'moving_mask', [], ... filename for mask of <moving_filename>
   'BST_masks', true, ... Interpret mask files as masks saved by BrainSuite (allow inconsistent headers)
   'mask_op', [],     ... Error mask: and / or of indivdual masks
   'pngout', false,   ... write overlay images in png format
   'reg_res', 1.8,    ... Resolution in mm (isotropic) to be used for registration
   'nbins', 128,      ... Number of bins for histeq and MI based registration
   'parzen_width', 8, ... Width of parzen window for histogram estimation
   'log_lookup', false, ... When true, uses log lookup table for computing log()
   'log_thresh', 1/(128*128*2), ... Lower threshold to skip log() when NOT using log-lookup; roughly related to nbins
   'step_size', [],   ... Step size for cost optimization, when empty sets appropriately depending on similarity
   'CFn_opts', struct(), ... cost function options, when empty sets appropriately depending on similarity & other options (overrides other options)
   'nthreads', 6,    ... Number of (possible) parallel threads to use
   'non_uniformity_correction', false, ... Makes sense only if images are corsely-pre-registered
   'init_opts', struct('init_method', 'search', 'search_range',[60 60 60], ...
                       'search_delta',[15 15 15], 'search_imres',5), ...
   'verbose', false,  ... Show optimization/processing details
   'debug', false  ... plot additional images etc. for debugging; independent from verbose
   );

msk_thresh = 0.2;

if(~exist('opts','var')),
   opts = defaultoptions;
else
   tags = fieldnames(defaultoptions);
   for i=1:length(tags)
      if(~isfield(opts,tags{i})),  opts.(tags{i})=defaultoptions.(tags{i}); end
   end
   if(length(tags)~=length(fieldnames(opts))),
      warning('BDP:UnknownOptions','unknown options found');
   end
end

if isempty(opts.similarity)
   error('Please specify the similarity of the images.');
elseif ~ismember(lower(opts.similarity), {'sd', 'mi', 'cr', 'inversion' 'bdp'})
   error('This similarity measure is unknown: %s', opts.similarity)
else
   opts.similarity = lower(opts.similarity);
end

if length(opts.reg_res)>1
   error('options.reg_res should be scalar');
end

if opts.dof<3 && length(opts.axis_translation)~=opts.dof
   error('opts.axis_translation is not set correctly for low dof registration');
end

if ~ismember(opts.dof, [1 2 3 6 12])
   error('opts.dof can not be %d. It must be one of the following: 1, 2, 3, 6 or 12.', opts.dof);
end

if isempty(opts.intensity_norm)
   switch opts.similarity
      case {'mi', 'cr'}
         opts.intensity_norm = 'none';
      case 'sd'
         opts.intensity_norm = 'histeq';
      case {'inversion', 'bdp'}
         opts.intensity_norm = 'sigmoid';
   end
elseif ~ismember(opts.intensity_norm, {'histeq', 'sigmoid', 'none'})
   error('Unknown opts.intensity_norm: %s', opts.intensity_norm)
end

if isempty(opts.mask_op)
   switch opts.similarity
      case {'mi', 'cr'}
         opts.mask_op = 'and';
      case 'sd'
         opts.mask_op = 'or';
      case {'inversion', 'bdp'}
         if opts.dof<=6
            opts.mask_op = 'and'; % default BDP behavior
         else
            opts.mask_op = 'or';
         end
   end
elseif ~ismember(opts.mask_op, {'and', 'or'})
   error('Unknown opts.mask_op: %s', opts.mask_op)
end


% temp_dir
workDir = tempname();
mkdir(workDir);

% define variables
workDir_static_file = fullfile(workDir, 'static_file.nii.gz');
workDir_static_mask = fullfile(workDir, 'static_mask.nii.gz');
workDir_moving_mask = fullfile(workDir, 'moving_mask.nii.gz');
workDir_mov_file_moved = fullfile(workDir, 'moving_file_moved.nii.gz');
workDir_mov_mask_moved = fullfile(workDir, 'moving_mask_moved.nii.gz');

png_centroid_fix = [remove_extension(output_filename) '.approx_align'];
png_out = [remove_extension(output_filename) '.output'];


%% check image dimension
fprintf('\nReading input data...')
mov_nii = load_untouch_nii_gz(moving_filename, workDir);
static_nii = load_untouch_nii_gz(static_filename, workDir);

if ndims(mov_nii.img)==4
   mov_nii.hdr.dime.dim(1) = 3;
   mov_nii.hdr.dime.dim(5) = 1;
   mov_nii.img = squeeze(mov_nii.img(:,:,:,1));
   moving_filename_3D = save_untouch_nii_gz(mov_nii, fullfile(workDir, 'moving_file_3D.nii.gz'));
   fprintf('\nFound 4D moving volume: %s\nBDP will use first volume for registration.\n', escape_filename(moving_filename));
elseif ndims(mov_nii.img)>4
   error('Unsupported dimension')
else
   moving_filename_3D = moving_filename;
end

if ndims(static_nii.img)==4
   static_nii.hdr.dime.dim(1) = 3;
   static_nii.hdr.dime.dim(5) = 1;
   static_nii.img = squeeze(static_nii.img(:,:,:,1));
   static_filename_3D = save_untouch_nii_gz(static_nii, fullfile(workDir, 'static_file_3D.nii.gz'));
   fprintf('\nFound 4D static volume: %s\nBDP will use first volume for registration.\n', escape_filename(static_filename));  
elseif ndims(static_nii.img)>4
   error('Unsupported dimension')
else
   static_filename_3D = static_filename;
end

%% generate/set/transfer masks in appropriate coordinates
fprintf('\nSetting/generating masks...')

% check image headers
[~, static_filename_3D] = check_nifti_file(static_filename_3D, workDir);
[~, moving_filename_3D] = check_nifti_file(moving_filename_3D, workDir);

if opts.dof<3
   w_tempfile = fullfile(workDir, ['moving_file.3D.RAS.' randstr(12) '.nii.gz']);
   reorient_nifti_sform(moving_filename_3D, w_tempfile);
   moving_filename_3D = w_tempfile;
   
   w_tempfile = fullfile(workDir, ['static_file.3D.RAS.' randstr(12) '.nii.gz']);
   [~, rmat] = reorient_nifti_sform(static_filename_3D, w_tempfile);
   if ~isequal(rmat, eye(3))
      static_filename_3D = [remove_extension(static_filename) '.3D.RAS.nii.gz'];
      copyfile(w_tempfile, static_filename_3D);
      fprintf(bdp_linewrap(['\nDOF<3 and hence Static file has been rotated to RAS format and saved. '...
         'The registed output would be written in same coordinate system as RAS file: '...
         escape_filename(static_filename_3D)]));
   end
end

% moving masks
if isempty(opts.moving_mask)
   opts.moving_mask = maskHeadPseudoHist(moving_filename_3D, workDir_moving_mask, true);
   
elseif opts.BST_masks
   w_tempfile = fullfile(workDir, ['moving_file.3D.RAS.' randstr(12) '.nii.gz']);
   reorient_nifti_sform(moving_filename_3D, w_tempfile);
   moving_filename_3D = w_tempfile;
   
   reorient_nifti_sform(opts.moving_mask, workDir_moving_mask);
   [~, workDir_moving_mask] = fixBSheader(moving_filename_3D, workDir_moving_mask, ...
      [remove_extension(workDir_moving_mask) '.BSfix.nii.gz']);
   opts.moving_mask = workDir_moving_mask;
   
else
   interp3_nii(opts.moving_mask, moving_filename_3D, workDir_moving_mask, 'nearest');
   opts.moving_mask = workDir_moving_mask;
end

% static masks
if isempty(opts.static_mask)
   opts.static_mask = maskHeadPseudoHist(static_filename_3D, workDir_static_mask, true);
   
elseif opts.BST_masks
   w_tempfile = fullfile(workDir, ['static_file.3D.RAS.' randstr(12) '.nii.gz']);
   [~, rmat] = reorient_nifti_sform(static_filename_3D, w_tempfile);
   if ~isequal(rmat, eye(3))
      static_filename_3D = [remove_extension(static_filename) '.3D.RAS.nii.gz'];
      copyfile(w_tempfile, static_filename_3D);
      fprintf(bdp_linewrap(['\nStatic file has been rotated to RAS format and saved. The registed output ' ...
         'would be written in same coordinate system as RAS file: ' escape_filename(static_filename_3D)]));
   end
   
   reorient_nifti_sform(opts.static_mask, workDir_static_mask);
   [~, workDir_static_mask] = fixBSheader(static_filename_3D, workDir_static_mask, ...
      [remove_extension(workDir_static_mask) '.BSfix.nii.gz']);
   opts.static_mask = workDir_static_mask;

else
   interp3_nii(opts.static_mask, static_filename_3D, workDir_static_mask, 'nearest');
   opts.static_mask = workDir_static_mask;
end


%% Save masks properly & fix headers
mov_nii = load_untouch_nii_gz(moving_filename_3D, true, workDir);
static_nii = load_untouch_nii_gz(static_filename_3D, true, workDir);
[dm, xm, ym, zm, resm] = get_original_grid_data(opts.moving_mask);
[ds, xs, ys, zs, ress] = get_original_grid_data(opts.static_mask);

% override grid to change it to be along acquisition axis
if opts.dof<=3   
   T_m = diag(abs([mov_nii.hdr.dime.pixdim(2:4) 1]));
   szm = size(dm).*resm;
   T_m(1,end) = -szm(1)/2;
   T_m(2,end) = -szm(2)/2;
   T_m(3,end) = -szm(3)/2;
   mov_nii.hdr.hist.srow_x = T_m(1,:);
   mov_nii.hdr.hist.srow_y = T_m(2,:);
   mov_nii.hdr.hist.srow_z = T_m(3,:);
   [~, xm, ym, zm] = get_original_grid_data(mov_nii);
   
   T_s = diag(abs([static_nii.hdr.dime.pixdim(2:4) 1]));
   szs = size(ds).*ress;
   T_s(1,end) = -szs(1)/2;
   T_s(2,end) = -szs(2)/2;
   T_s(3,end) = -szs(3)/2;
   static_nii.hdr.hist.srow_x = T_s(1,:);
   static_nii.hdr.hist.srow_y = T_s(2,:);
   static_nii.hdr.hist.srow_z = T_s(3,:);
   [~, xs, ys, zs] = get_original_grid_data(static_nii);
end

mov_mask_nii = mov_nii;
mov_mask_nii.img = dm>0;
save_untouch_nii_gz(mov_mask_nii, workDir_moving_mask, 2);

static_mask_nii = static_nii;
static_mask_nii.img = ds>0;
save_untouch_nii_gz(static_mask_nii, workDir_static_mask, 2);

clear static_mask_nii defaultoptions mov_mask_nii

%% Fix for huge translation - approximate head structure alignment
center_disp_M = eye(4);
center_disp = zeros(3,1);

if ~strcmpi(opts.init_opts.init_method, 'none')
   fprintf('\nMatching centroids (approx. align)...')
   xm = xm(dm>0);
   ym = ym(dm>0);
   zm = zm(dm>0);
   % Im = double(mov_nii.img(dm>0));
   % cm = (([xm(:) ym(:) zm(:)]'*Im(:))./sum(Im(:))); % centroid of image 
   cm = [mean(xm(:));  mean(ym(:));  mean(zm(:))]; % centroid of mask
   
   xs = xs(ds>0);
   ys = ys(ds>0);
   zs = zs(ds>0);
   % Is = double(static_nii.img(ds>0));
   % cs = (([xs(:) ys(:) zs(:)]'*Is(:))./sum(Is(:))); % centroid of image 
   cs = [mean(xs(:));  mean(ys(:));  mean(zs(:))]; % centroid of mask 
      
   if opts.dof<3
      center_disp(opts.axis_translation) = cm(opts.axis_translation) - cs(opts.axis_translation);
   else
      center_disp = (cm - cs);
   end
   center_disp_M(1:3,4) = center_disp;
end
clear dm xm ym zm xs ys zs ds cm cs Im Is

affine_transform_nii(mov_nii, inv(center_disp_M), [0 0 0]', workDir_mov_file_moved);
affine_transform_nii(workDir_moving_mask, inv(center_disp_M), [0 0 0]', workDir_mov_mask_moved);


%% mask data & reslice to same grid
fprintf('\nMatching resolution of the images...')
mov_nii = load_untouch_nii_gz(workDir_mov_file_moved, true, workDir);

% first reslice to lower resolution scan
[~, X_grid_mv, Y_grid_mv, Z_grid_mv, res_mv, Tmv] = get_original_grid_data(mov_nii);
[~, X_grid_st, Y_grid_st, Z_grid_st, res_st, Tst] = get_original_grid_data(static_nii);

common_res = max(res_mv, res_st); % apparent resolution of images in mm, used for registration 

if ~isequal(common_res, res_mv)
   xg = min(X_grid_mv(:)):common_res(1):max(X_grid_mv(:)); clear X_grid_mv;
   yg = min(Y_grid_mv(:)):common_res(2):max(Y_grid_mv(:)); clear Y_grid_mv;
   zg = min(Z_grid_mv(:)):common_res(3):max(Z_grid_mv(:)); clear Z_grid_mv;
   [xg, yg, zg] = ndgrid(xg, yg, zg);
   myreslice_nii_match_res(mov_nii, 'linear', xg, yg, zg, workDir_mov_file_moved);
   myreslice_nii_match_res(workDir_mov_mask_moved, 'nearest', xg, yg, zg, workDir_mov_mask_moved);
   clear xg yg zg
end
clear X_grid_mv Y_grid_mv Z_grid_mv

if ~isequal(common_res, res_st)
   xg = min(X_grid_st(:)):common_res(1):max(X_grid_st(:)); clear X_grid_st;
   yg = min(Y_grid_st(:)):common_res(2):max(Y_grid_st(:)); clear Y_grid_st;
   zg = min(Z_grid_st(:)):common_res(3):max(Z_grid_st(:)); clear Z_grid_st;
   [xg, yg, zg] = ndgrid(xg, yg, zg);
   myreslice_nii_match_res(static_nii, 'linear', xg, yg, zg, workDir_static_file);
   myreslice_nii_match_res(workDir_static_mask, 'nearest', xg, yg, zg, workDir_static_mask);
   clear xg yg zg
else
   save_untouch_nii_gz(static_nii, workDir_static_file);
end
clear X_grid_st Y_grid_st Z_grid_st

% reslice data on same grid
[~, ~, X_grid, Y_grid, Z_grid, ~, Gpos] = reslice_vol_same_grid(workDir_mov_mask_moved, workDir_static_mask, ...
   opts.reg_res, 'linear');
m = myreslice_nii_match_res(workDir_mov_file_moved, 'linear', X_grid, Y_grid, Z_grid);
s = myreslice_nii_match_res(workDir_static_file, 'linear', X_grid, Y_grid, Z_grid);
m_mask = myreslice_nii_match_res(workDir_mov_mask_moved, 'linear', X_grid, Y_grid, Z_grid);
s_mask = myreslice_nii_match_res(workDir_static_mask, 'linear', X_grid, Y_grid, Z_grid);

% mask data
m(~(m_mask > msk_thresh)) = 0;
s(~(s_mask > msk_thresh)) = 0;

% compute rigid-reg origin - as used in affine_transform_3d_double.c or affine_transform_3dvol_double.c
ref_loc(1,1) = (Gpos{1}(2) + Gpos{1}(end))/2;
ref_loc(2,1) = (Gpos{2}(2) + Gpos{2}(end))/2;
ref_loc(3,1) = (Gpos{3}(2) + Gpos{3}(end))/2;
clear  X_grid Y_grid Z_grid Gpos

%% intensity correction
fprintf('\nNormalizing intensity of images... ')

if strcmpi(opts.intensity_norm, 'histeq')
   [mImg, m_low, m_high] = normalize_intensity(m, [0.02 99.75], m_mask>msk_thresh);
   [sImg, s_low, s_high] = normalize_intensity(s, [0.02 99.75], s_mask>msk_thresh);
   temp_mv = zeros(size(mImg));
   temp_mv(m_mask>msk_thresh) = histeqSmooth(mImg(m_mask>msk_thresh), sImg(s_mask>msk_thresh), opts.nbins);
   mImg = temp_mv;
   clear temp_mv
   
elseif strcmpi(opts.intensity_norm, 'sigmoid')
   [mImg, sImg, m_high, s_high] = heuristicSigmoidScale(m, s, m_mask>msk_thresh, s_mask>msk_thresh);
   
elseif strcmpi(opts.intensity_norm, 'none')
   [mImg, m_low, m_high] = normalize_intensity(m, [0.02 99.9], m_mask>msk_thresh);
   [sImg, s_low, s_high] = normalize_intensity(s, [0.02 99.9], s_mask>msk_thresh);

else
   error('Unknown opts.intensity_norm: %s', opts.intensity_norm)
end

if opts.debug
   figure;
   [n, b] = histcParzen(s(s_mask>msk_thresh)/s_high, [0 2], 100);
   plot(b, n, 'g'); hold on;
   [n, b] = histcParzen(m(m_mask>msk_thresh)/m_high, [0 2], 100);
   plot(b, n, 'b'); hold on;
   
   [n, b] = histcParzen(sImg(s_mask>msk_thresh), [0 2], 100);
   plot(b, n, '--g', 'LineWidth',2); hold on;
   [n, b] = histcParzen(mImg(m_mask>msk_thresh), [0 2], 100);
   plot(b, n, '--b', 'LineWidth',2); hold off;
   title('Histogram (separate mask)')
   legend('s', 'm','sImg', 'mImg')
   display_volume([s/s_high; sImg; m/m_high; mImg], [0 1]);
end
clear ds m b mov_mask_nii mov_nii static_nii


%% Setup registration options and start affine registration
reg_opt.dof = opts.dof;
reg_opt.axis_translation = opts.axis_translation;
reg_opt.similarity = opts.similarity;
reg_opt.verbose = opts.verbose;
reg_opt.step_size = opts.step_size;
reg_opt.nthreads = opts.nthreads;
reg_opt.mask_op = opts.mask_op;

reg_opt.init_method = opts.init_opts.init_method;
if ~strcmpi(reg_opt.init_method, 'none')
   flds = {'search_range', 'search_delta', 'search_imres'};
   for k = 1:length(flds)
      if isfield(opts.init_opts, flds{k})
         reg_opt.(flds{k}) = opts.init_opts.(flds{k});
      end
   end
end

if strcmpi(reg_opt.similarity, 'mi')
   reg_opt.CFn_opts = catstruct(regOptsMI(opts), opts.CFn_opts);
   
elseif strcmpi(reg_opt.similarity, 'bdp') || strcmpi(reg_opt.similarity, 'inversion')
   reg_opt.CFn_opts = catstruct(struct('nbins', 300, 'win_width', 30, 'hybrid_mask_op', false), opts.CFn_opts);
   reg_opt.CFn_opts.mi = regOptsMI(opts);
end

% non-uniformity correction
if opts.non_uniformity_correction && (strcmpi(reg_opt.similarity, 'bdp') || strcmpi(reg_opt.similarity, 'inversion'))
   fprintf('\nPerforming non-uniformity correction...')
   bfc_opts.num_threads = opts.nthreads;
   bfc_opts.verbose = opts.verbose;
   bfc_opts.penalty = 5e-3;
   bfld = matchIntensity_EPI_T1(sImg, mImg, opts.reg_res, s_mask, s_mask, m_mask, bfc_opts);
   
   sImg = sImg .* bfld;
   sImg(sImg>1) = 1;
   sImg(sImg<0) = 0;
   normalize_intensity(sImg, [0 100], s_mask);
   clear bfld bfc_opts
end

% Visualize input images
if opts.pngout
   overlay_volumes2png(mImg, sImg, [0 1], png_centroid_fix);
   overlay_volumes2png(mImg, sImg, [0 1], [png_centroid_fix '.static.edge'], 'edge');
   overlay_volumes2png(sImg, mImg, [0 1], [png_centroid_fix '.moving.edge'], 'edge');
   fprintf('\nSaved approx. aligned overlay PNG files.')
end

fprintf('\nPerforming affine registration with %d degrees of freedom...', opts.dof)
[~, x, Ireg] = register_volumes_affine(mImg, sImg, m_mask, s_mask, opts.reg_res, reg_opt);

% convert to mm & concatenate centroid shift
fprintf('\nSaving output files...')

if opts.dof<3
   t = zeros(3,1);
   t(opts.axis_translation) = x;
   x = t;
end

x(1:3) = opts.reg_res * x(1:3); % convert to mm

if length(x) == 3
   M_reg = par2affineMat(x(1:3));
elseif length(x) == 6
   M_reg = par2affineMat(x(1:3), x(4:6));
elseif length(x) == 12
   M_reg = par2affineMat(x(1:3), x(4:6), [], x(7:12));
else
   error('Length of x must be 3, 6 or 12');
end

% for returning & saving
x_param = x;
x_param(1:3) = x(1:3) + center_disp;

M_world = center_disp_M * M_reg;
origin = ref_loc;
save([remove_extension(output_filename) '.rigid_registration_result.mat'], ...
   'M_world', 'origin', 'moving_filename', 'static_filename');

affine_transform_nii(moving_filename, inv(M_world), ref_loc, [remove_extension(output_filename) '.headr_transform.nii.gz']);
interp3_nii([remove_extension(output_filename) '.headr_transform.nii.gz'], static_filename_3D, output_filename, 'linear');

% Visualize output file
if opts.pngout   
   mImg = Ireg;
   overlay_volumes2png(mImg, sImg, [0 1], png_out);
   overlay_volumes2png(mImg, sImg, [0 1], [png_out '.static.edge'], 'edge');
   overlay_volumes2png(sImg, mImg, [0 1], [png_out '.moving.edge'], 'edge');
   fprintf('\nSaved output overlay PNG files.')
end

fprintf('\n'); % finish on a new line
rmdir(workDir, 's');
clearvars -except M_world ref_loc x_param
end



function [mImgP, sImgP, m_high, s_high] = heuristicSigmoidScale(mImg, sImg, mask_mImg, maks_sImg)
% Applies a heuristic sigmoid scaling to images - based on register_files_rigid_mprage_0_diffusion()
% Also see sigmoidScale_b0image() for another heuristic Sigmoid Scaling.
%

% intensity correction - MPRAGE (moving file)
upper_prct = 92;
c = 10;
[~, ~, m_high] = normalize_intensity(mImg, [0 upper_prct], mask_mImg);
upper_prct = upper_prct/100;
m_high = m_high/upper_prct;
mImgP = mImg/m_high;
mtemp = mImgP(mImgP>=upper_prct) - upper_prct;
mtemp = erf(mtemp*c)* (1-upper_prct);
mImgP(mImgP>=upper_prct) = mtemp + upper_prct;
mImgP(mImgP>1) = 1;
mImgP(mImgP<0) = 0;

% intensity correction b=0 image (static file)
sigmoid_lvl = 0.5;
upper_prct = 95;
[~, s_low, s_high] = normalize_intensity(sImg, [0.02 upper_prct], maks_sImg);
if s_low<0, s_low=0; end % It is magnitude image after all !
upper_prct = upper_prct/100;
s_high = s_high/upper_prct;
sImgP = (sImg-s_low)/(s_high-s_low);
stemp = sImgP(sImgP>=sigmoid_lvl) - sigmoid_lvl;
stemp = erf(stemp)* (1-sigmoid_lvl); %stemp./sqrt(1+stemp.^2)* (1-sigmoid_lvl);
sImgP(sImgP>=sigmoid_lvl) = stemp + sigmoid_lvl;
sImgP(sImgP>1) = 1;
sImgP(sImgP<0) = 0;

end

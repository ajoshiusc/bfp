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


function [t_epi, M, O_trans, Spacing] = EPI_correct_files_registration_INVERSION(epi_filename, struct_filename, ...
                                                                        epi_output_filename, struct_output_filename, opts)
% Interface for nifti input/output files. Estimates the distortion due to B0 field inhomogenity in
% EPI scans using one dimenstional non-rigid bspline-based registration with various similariy
% measures.
%
% All masks should be in same coordinate as respective files. See maskDWI.m and mask_b0_setup()
% in coregister_diffusion_mprage_pipeline.m  

%% Check/set input options
defaultoptions = struct(...
   'phase_encode_direction', 'y-', ...
   'similarity', 'inversion-epi-t1', ... bdp / inversion-epi / inversion-t1 / inversion-epi-t1 / mi / sd
   'rigid_similarity', 'bdp', ... bdp / inversion / mi / cr
   'struct_mask', [], ... required
   'epi_mask', [], ... required
   'epi_mask_less_csf', [], ... optional
   'intensity_correct', true, ... Jacobian based intensity modulation when performing distortion-correction
   'intensity_norm', 'sigmoid',... sigmoid/histNorm/histeq/none -- intensity normalization prior to registration
   'reg_res', 1.5, ... in mm; isotropic resolution at which data is sampled for registration
   'rigid_refine', 1, ... refine rigid parameters during distortion correction?
   'init_ctrl_spacing', 24, ... in mm; Initial b-spine control point spacing
   'step_size', [], ...  step size for gradient descent, when empty sets appropriately
   'step_size_scale_iter', 0.5,... for backtracking
   'MaxFunEvals', 100, ... max number of gradient evaluation
   'Penalty', [1 0.00000001 0.0000001], ... [alpha kx bigkx ky], alpha is scaling of overall regularization term
   'non_uniformity_correction', true, ... Apply a coarse bias-field correction
   'bfield_penalty', 5e-3, ... 
   'bfield_file', [], ... pre-computed bfield nifti file
   'rigid_reg_mat', [], ... pre-computed rigid-registration result; skips rigid registration.
   'verbose', false, ...  when true, shows optimization iterations
   'debug', false, ...    when true, shows various detailed plots
   'num_threads', 6, ... number of parallel processing threads to use
   'pngout', false, ...
   'overlay_mode', 'redgreen' ...  % 'rview'/'redgreen'/'greenblue'/'yellowblue'
   );

if ~exist('opts','var')
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

if length(opts.reg_res)>1
   error('opts.reg_res should be scalar');
end

msk_thrs = 0.25; % threshold to convert contineous mask to binary mask
workDir = tempname(); 
mkdir(workDir);

overlay_rigid = [remove_extension(epi_output_filename) '.distorted'];
overlay_corrected = [remove_extension(epi_output_filename) '.corrected'];

%% = = = = = = set/load masks = = = = = = = %
disp('Loading data...')

% extract b=0 image from diffusion images
dwi = load_untouch_nii_gz(epi_filename);
if dwi.hdr.dime.dim(1) > 3
   dwi.img = dwi.img(:,:,:,1);
   dwi.hdr.dime.dim(1) = 3;
   dwi.hdr.dime.dim(5) = 1;
   epi_b0_filename = save_untouch_nii_gz(dwi, [remove_extension(epi_filename) '.0_diffusion.nii.gz']);
else
   epi_b0_filename = epi_filename;
end
clear dwi

% check masks
if isempty(opts.struct_mask)
   error('opts.struct_mask must be defined!');
elseif isempty(opts.epi_mask)
   error('opts.epi_mask must be defined!');
elseif isempty(opts.epi_mask_less_csf)
   opts.epi_mask_less_csf = opts.epi_mask;
end

workdir_epi_mask = fullfile(workDir, 'epi_mask.nii.gz');
fixBSheader(epi_b0_filename, opts.epi_mask, workdir_epi_mask);

workdir_epi_mask_less_CSF = fullfile(workDir, 'epi_mask_lesscsf.nii.gz');
fixBSheader(epi_b0_filename, opts.epi_mask_less_csf, workdir_epi_mask_less_CSF);

workdir_struct_mask = fullfile(workDir, 'struct_mask.nii.gz');
fixBSheader(struct_filename, opts.struct_mask, workdir_struct_mask);



%% = = = = Initial rigid registration = = = = %
if isempty(opts.rigid_reg_mat)
   fprintf('\nStarting rigid registration of input images...')
   reg_options = struct( ...
      'moving_mask', workdir_struct_mask, ...
      'static_mask', workdir_epi_mask, ...
      'pngout', opts.pngout, ...
      'verbose', opts.verbose, ...
      'debug', opts.debug, ...
      'similarity', opts.rigid_similarity, ...
      'dof', 6, ...
      'nthreads', opts.num_threads);
   
   [M_world, origin] = register_files_affine(struct_filename, epi_b0_filename, struct_output_filename, reg_options);         
   opts.rigid_reg_mat = [remove_extension(struct_output_filename) '.rigid_registration_result.mat'];
   fprintf('\nRigid registration done.\n')

else
   fprintf('\nUsing registration parameter from mat file: %s\n', opts.rigid_reg_mat);
   load(opts.rigid_reg_mat);
   fileout = [remove_extension(struct_output_filename) '.headr_transform.nii.gz'];
   affine_transform_nii(struct_filename, inv(M_world), origin, fileout);
end

%% = = = = = = = = = EPI correction = = = = = = = = %
fprintf('\nSetting up data for non-rigid registration based distortion correction...')

% load EPI data, mask
epi_mask = load_untouch_nii_gz(workdir_epi_mask);
epi_mask_less_csf = load_untouch_nii_gz(workdir_epi_mask_less_CSF);
epi_data = load_untouch_nii_gz(epi_b0_filename, true, workDir);
epi_mask.img = epi_mask.img>0;
epi_mask_less_csf.img = epi_mask_less_csf.img>0;

struct_fname = [remove_extension(struct_output_filename) '.headr_transform.nii.gz'];
struct_mask = affine_transform_nii(workdir_struct_mask, inv(M_world), origin);
struct_mask.img = struct_mask.img>0;

% Match resolution and interpolate at EPI grid
struct_data = epi_data; 
[~, x_epi, y_epi, z_epi, epi_res, T_epi] = get_original_grid_data(epi_mask);
struct_data.img = myreslice_nii_match_res(struct_fname, 'linear', x_epi, y_epi, z_epi);
struct_mask.img = myreslice_nii(struct_mask, 'linear', x_epi, y_epi, z_epi);
struct_mask.hdr = epi_data.hdr; % copy epi header

% Normalize image intensities
if strcmpi(opts.intensity_norm, 'histeq')
   struct_data.img = normalize_intensity(struct_data.img, [0.5 99.97], struct_mask.img>msk_thrs);
   epi_data.img = normalize_intensity(abs(epi_data.img), [0 99.8], epi_mask.img>msk_thrs);
   temp = zeros(size(struct_data.img));
   temp(struct_mask.img>msk_thrs) = histeqSmooth(struct_data.img(struct_mask.img>msk_thrs), ...
      epi_data.img(epi_mask.img>msk_thrs), opts.nbins);
   struct_data.img = temp;
   clear temp
   
elseif strcmpi(opts.intensity_norm, 'histNorm')
   s1 = normalize_intensity(struct_data.img, [0.2 99.99], struct_mask.img>msk_thrs);
   temp = zeros(size(s1));
   temp(struct_mask.img>msk_thrs) = histeqSmooth(s1(struct_mask.img>msk_thrs), 1200);
   struct_data.img = temp;
   
   % epi_data.img = normalize_intensity(abs(epi_data.img), [0.2 99.5], epi_mask.img>msk_thrs);
   e1 = normalize_intensity(abs(epi_data.img), [0.1 99.5], epi_mask.img>msk_thrs);
   temp = zeros(size(e1));
   temp(epi_mask.img>msk_thrs) = histeqSmooth(e1(epi_mask.img>msk_thrs), 1200);
   epi_data.img = temp;
   clear temp s1 e1
   
elseif strcmpi(opts.intensity_norm, 'none')
   struct_data.img = normalize_intensity(struct_data.img, [0.2 99.99], struct_mask.img>msk_thrs);
   epi_data.img = normalize_intensity(abs(epi_data.img), [0.1 99.5], epi_mask.img>msk_thrs);
   
elseif strcmpi(opts.intensity_norm, 'sigmoid')   
   struct_data.img = normalize_intensity(struct_data.img, [0.2 99.99], struct_mask.img>msk_thrs);
   [~,l,h] = normalize_intensity(abs(epi_data.img), [0.1 70], epi_mask.img>msk_thrs);
   temp = (double(abs(epi_data.img))-l)/(2*h-l); % set 70th percentile intensity to 0.5
   temp(temp>5.5) = 5.5;
   temp(temp<0) = 0;
   epi_data.img = sigmoidScale_b0image(temp);   
   
else
   error('Unknown opts.intensity_norm: %s', opts.intensity_norm)
end

% Correct intensities for bias-field effects in b=0
if opts.non_uniformity_correction
   if isempty(opts.bfield_file)
      fprintf('\nCorrecting intensity non-uniformity in b=0 image...')
      bcorr_res = 2; % in mm
      [epiImg, m_sImg, m_epi_less_csf, sImg, m_epi, Xg, Yg, Zg, Xe, Ye, Ze] = resliceForNonUniformityCorr( ...
         epi_data, epi_mask, epi_mask_less_csf, struct_data, struct_mask, epi_res, bcorr_res);

      bfc_opts.num_threads = opts.num_threads;
      bfc_opts.verbose = opts.verbose;
      bfc_opts.penalty = opts.bfield_penalty;
      bfld = matchIntensity_EPI_T1(epiImg, sImg, bcorr_res, m_epi, m_epi_less_csf, m_sImg, bfc_opts);
      clear epiImg sImg bcorr_res m_epi_less_csf m_sImg m_epi bfc_opts
      
      % Interpolate on epi grid
      bfield = epi_data;     
      bfield.img = interpn(Xg, Yg, Zg, bfld, Xe, Ye, Ze, 'linear', 1); 
      save_untouch_nii_gz(bfield, [remove_extension(epi_output_filename) '.bfield.nii.gz'], 16);
   else
      bfield = load_untouch_nii_gz(opts.bfield_file);
      bfield.img = double(bfield.img);
   end
   epiImg_bfc = epi_data.img .* bfield.img;
   epiImg_bfc(epiImg_bfc>1) = 1;
   epiImg_bfc(epiImg_bfc<0) = 0;
   clear bfield bfld Xg Yg Zg
else
   epiImg_bfc = epi_data.img;
end
clear x_epi y_epi z_epi

% Interpolate data on registration grid (parallel to epi plane)
grid_res = (opts.reg_res)./epi_res; % unit of voxels in original EPI image
reg_epi_res = opts.reg_res*[1 1 1]; % in mm, resolution at which registration is performed
epi_size =  size(epi_mask.img);

reg_grid{1} = -2:grid_res(1):(epi_size(1)+1); % 2 voxel padding on all sides to allow adding some 
reg_grid{2} = -2:grid_res(2):(epi_size(2)+1); % slices & minimize bspline effects near boundary.
reg_grid{3} = -2:grid_res(3):(epi_size(3)+1); % This is further refined with bounding box later.

[X_reg_epi, Y_reg_epi, Z_reg_epi] = ndgrid(reg_grid{1}, reg_grid{2}, reg_grid{3});
[x_epi, y_epi, z_epi] = ndgrid(0:(epi_size(1)-1), 0:(epi_size(2)-1), 0:(epi_size(3)-1));

epi_data = interpn(x_epi, y_epi, z_epi, double(epi_data.img), X_reg_epi, Y_reg_epi, Z_reg_epi, 'linear', 0);
epiImg_bfc = interpn(x_epi, y_epi, z_epi, double(epiImg_bfc), X_reg_epi, Y_reg_epi, Z_reg_epi, 'linear', 0);
epi_mask = interpn(x_epi, y_epi, z_epi, double(epi_mask.img), X_reg_epi, Y_reg_epi, Z_reg_epi, 'linear', 0);
epi_mask_less_csf = interpn(x_epi, y_epi, z_epi, double(epi_mask_less_csf.img), X_reg_epi, Y_reg_epi, Z_reg_epi, 'linear', 0);
struct_data = interpn(x_epi, y_epi, z_epi, double(struct_data.img), X_reg_epi, Y_reg_epi, Z_reg_epi, 'linear', 0);
struct_mask = interpn(x_epi, y_epi, z_epi, double(struct_mask.img), X_reg_epi, Y_reg_epi, Z_reg_epi, 'linear', 0);
clear x_epi y_epi z_epi X_reg_epi Y_reg_epi Z_reg_epi

if opts.debug
   nbins = 200; parzen_width = 16;
   out_epiImg = zeros(size(epi_data));
   [out_epiImg(epi_mask_less_csf>msk_thrs), ~, hMap_epiImg, b_epiImg] = histeqSmooth(1-epi_data(epi_mask_less_csf>msk_thrs), ...
      struct_data(struct_mask>msk_thrs), nbins, parzen_width);
   out_epi_bfc = zeros(size(epi_data));
   [out_epi_bfc(epi_mask_less_csf>msk_thrs), ~, hMap_epi_bfc, T_grid] = histeqSmooth(1-epiImg_bfc(epi_mask_less_csf>msk_thrs),...
      struct_data(struct_mask>msk_thrs), nbins, parzen_width);
   
   figure; plot(b_epiImg, hMap_epiImg, '*-'); hold on;  plot(T_grid, hMap_epi_bfc, 'g*-'); 
   axis equal; xlim([0 1]); ylim([0 1]); grid on; legend('epiImg', 'epiImg_bfc', 'Location', 'North')
   display_volume([struct_data; out_epiImg; out_epi_bfc], [0 1]); colormap jet
   display_volume([struct_data, out_epiImg, out_epi_bfc], [0 1], 1); colormap jet
   
   % Parzen estimate of PDFs
   [~,~, nT1, bT1] = histcParzen(struct_data(struct_mask>msk_thrs), [0 1], nbins, parzen_width);
   [~,~, nEPI_bfc, bEPI_bfc] = histcParzen(epiImg_bfc(epi_mask_less_csf>msk_thrs), [0 1], nbins, parzen_width);
   [~,~, nEPI, bEPI] = histcParzen(epi_data(epi_mask_less_csf>msk_thrs), [0 1], nbins, parzen_width);
   figure; plot(bT1, nT1/sum(nT1), 'g', 'LineWidth',2); hold on;
   plot(bEPI_bfc, nEPI_bfc/sum(nEPI_bfc), '--b'); hold on;
   plot(bEPI, nEPI/sum(nEPI), 'b', 'LineWidth',2); hold off;
   title('Parzen estimate of PDF (separate mask)')
   legend('T1', 'EPI_bfc', 'EPI', 'Location', 'North')
   
   % Gradient plots
   figure; subplot(1,2,1); plot(T_grid, hMap_epi_bfc, 'b', 'LineWidth',2); xlim([-0.2 1.2]);
   title('Intensity Map: Inv EPI -> T1'); xlabel('Inv EPI Intensity'); ylabel('T1 Intensity');
   subplot(1,2,2); plot(T_grid, gradient(hMap_epi_bfc, T_grid(2)-T_grid(1)), 'b', 'LineWidth',2);
   hold on; plot(T_grid, gradient(smoothGaussian1D(hMap_epi_bfc, 2.5), T_grid(2)-T_grid(1)), 'r', 'LineWidth',1);
   title('Gradient of Intensity Map'); xlabel('Gradient'); ylabel('T1 Intensity'); xlim([-0.2 1.2]);
   legend('gradient()', 'gradient(smoothGaussian1D())', 'Location', 'North'); grid on;
end
clear hist_bins histeqT1 histeqT2 out_epi_bfc out_epiImg


% keep data only in bounding box
pad_vox = 1;
size_struct_data = size(struct_data);
[bb_1, bb_2, bb_3] = find_bounding_box((epi_mask>msk_thrs) | (struct_mask>msk_thrs), pad_vox);
[bb_1, bb_2, bb_3] = pad_bbox_PED(bb_1, bb_2, bb_3, size_struct_data, opts.phase_encode_direction, 2);
[bb_1, bb_2, bb_3, mod_msk] = trim_msk_PED(bb_1, bb_2, bb_3, (epi_mask>msk_thrs)&(struct_mask>msk_thrs), ...
   opts.phase_encode_direction, pad_vox);

mod_msk = mod_msk(bb_1, bb_2, bb_3);
epi_mask = epi_mask(bb_1, bb_2, bb_3) .* mod_msk;
epi_mask_less_csf = epi_mask_less_csf(bb_1, bb_2, bb_3) .* mod_msk;
struct_mask = struct_mask(bb_1, bb_2, bb_3) .* mod_msk;
sImg = struct_data(bb_1, bb_2, bb_3);
epiImg = epi_data(bb_1, bb_2, bb_3);
epiImg_bfc = epiImg_bfc(bb_1, bb_2, bb_3);

% track registration grid to obtain rigid reg origin
rigid_reg_grid{1} = reg_grid{1}(bb_1);
rigid_reg_grid{2} = reg_grid{2}(bb_2);
rigid_reg_grid{3} = reg_grid{3}(bb_3);
clear struct_data epi_data mod_msk

% permute data to make 1st dimension the phase encode dim of EPI
switch opts.phase_encode_direction
   case {'x', 'x-'}
      permute_vec = [1 2 3 4];
   case {'y', 'y-'}
      permute_vec = [2 1 3 4];
   case {'z', 'z-'}
      permute_vec = [3 2 1 4];
end
epiImg = permute(epiImg, permute_vec);
epiImg_bfc = permute(epiImg_bfc, permute_vec);
epi_mask = permute(epi_mask, permute_vec);
epi_mask_less_csf = permute(epi_mask_less_csf, permute_vec);
sImg = permute(sImg, permute_vec);
struct_mask = permute(struct_mask, permute_vec);
reg_epi_res = reg_epi_res(permute_vec(1:3));
grid_res_reg = grid_res(permute_vec(1:3));
rigid_reg_grid = rigid_reg_grid(permute_vec(1:3));

% registration
non_rigid_opts = struct(...
   'Penalty', opts.Penalty, ...
   'Spacing', 2.^(round(log2(opts.init_ctrl_spacing./reg_epi_res))), ... % b-spline control point spacing in units of voxels,
   ...                                                 % power of 2 as bspline refinement only reduces
   ...                                                 % size by a factor of 2
   'nthreads', opts.num_threads, ...
   'step_size', opts.step_size, ...
   'step_size_scale_iter', opts.step_size_scale_iter,...
   'intensity_correct', opts.intensity_correct,...
   'MaxFunEvals', opts.MaxFunEvals, ...
   'debug', opts.debug, ...
   'verbose', opts.verbose, ...
   'rigid_scl', double(opts.rigid_refine), ...
   'histeqInterpolant', [], ...
   'histeqGradInterpolant', [], ...
   'similarity', opts.similarity ...
   );

% fix image size by padding zeros - make image size to be a perfect fit for control points
[epiImgP, st_vox, end_vox, szP] = bsplineFixVolSizePad(epiImg, non_rigid_opts.Spacing);
epiImgP_bfc = bsplineFixVolSizePad(epiImg_bfc, non_rigid_opts.Spacing);
sImgP = bsplineFixVolSizePad(sImg, non_rigid_opts.Spacing);
non_rigid_opts.MaskEPI = bsplineFixVolSizePad(epi_mask, non_rigid_opts.Spacing);
non_rigid_opts.MaskEPILessCSF = bsplineFixVolSizePad(epi_mask_less_csf, non_rigid_opts.Spacing);
non_rigid_opts.MaskStruct = bsplineFixVolSizePad(struct_mask, non_rigid_opts.Spacing);
clear epiImg sImg epiImg_bfc epi_mask epi_mask_less_csf struct_mask

% adjust rigid reg grid
for k = 1:length(rigid_reg_grid)
   g_st = [(1-st_vox(k)):-1]*grid_res_reg(k) + rigid_reg_grid{k}(1);
   g_end = [1:(szP(k)-end_vox(k))]*grid_res_reg(k) + rigid_reg_grid{k}(end);
   rigid_reg_grid{k} = [g_st rigid_reg_grid{k} g_end];
end

% Visualize input images
if opts.pngout
   overlay_volumes2png(sImgP, epiImgP_bfc, [0 1], overlay_rigid, opts.overlay_mode, [30 99]);
   overlay_volumes2png(epiImgP_bfc.*non_rigid_opts.MaskEPI, sImgP, [0 1], [overlay_rigid '.edgeM'], 'edge');
   overlay_volumes2png(sImgP, epiImgP_bfc.*non_rigid_opts.MaskEPI, [0 1], [overlay_rigid '.edgeD'], 'edge');
end


fprintf('\nRunning Non-rigid registration...')
% inspect data for debugging
if opts.debug
   overlay_volumes([non_rigid_opts.MaskStruct*0.5, sImgP, epiImgP*1.5; non_rigid_opts.MaskStruct*0.5, sImgP, epiImgP*1.5], ...
      [non_rigid_opts.MaskEPI, non_rigid_opts.MaskEPI, non_rigid_opts.MaskEPI; ...
      non_rigid_opts.MaskEPILessCSF, non_rigid_opts.MaskEPILessCSF, non_rigid_opts.MaskEPILessCSF]*0.5);
end

[~, ~, O_trans, Spacing, rigidX_sImg, t_deform] = EPI_correct_volume_bspline_registration(epiImgP, epiImgP_bfc, sImgP, ...
   reg_epi_res(1), non_rigid_opts);

% Visualize output images
if opts.pngout
   t_x = t_deform;
   [x_g, y_g, z_g] = ndgrid(1:size(epiImgP, 1), 1:size(epiImgP,2), 1:size(epiImgP,3));
   x_g = x_g + t_x;   
   Iepi_corr = interpn(epiImgP_bfc, x_g, y_g, z_g, 'cubic', 0);   
   if opts.intensity_correct
      [~, crct_grd] = gradient(t_x);
      Iepi_corr = Iepi_corr .* (1+crct_grd);
   end
   
   M = par2affineMat(rigidX_sImg(1:3), rigidX_sImg(4:6));
   Istruct_reg = affine_transform_3dvol_double(double(sImgP),double(M),double(3), opts.num_threads);
   
   overlay_volumes2png(Istruct_reg, Iepi_corr, [0 1], overlay_corrected, opts.overlay_mode, [30 99]);
   overlay_volumes2png(Iepi_corr.*non_rigid_opts.MaskEPI, Istruct_reg, [0 1], [overlay_corrected '.edgeM'], 'edge');
   overlay_volumes2png(Istruct_reg, Iepi_corr.*non_rigid_opts.MaskEPI, [0 1], [overlay_corrected '.edgeD'], 'edge');
end

fprintf('\nComposing all deformations...')

% t_x in original grid 
t_x = t_deform(st_vox(1):end_vox(1), st_vox(2):end_vox(2), st_vox(3):end_vox(3),1);
t_x = ipermute(t_x, permute_vec);
clear t_deform
t_x_reg = zeros(size_struct_data);
t_x_reg(bb_1, bb_2, bb_3) = t_x;
clear t_x t_deform

[X_reg_epi, Y_reg_epi, Z_reg_epi] = ndgrid(reg_grid{1}, reg_grid{2}, reg_grid{3});
[x_epi, y_epi, z_epi] = ndgrid(0:(epi_size(1)-1), 0:(epi_size(2)-1), 0:(epi_size(3)-1));
t_epi = interpn(X_reg_epi, Y_reg_epi, Z_reg_epi, double(t_x_reg), x_epi, y_epi, z_epi, 'cubic', 0);
t_epi = grid_res_reg(1)*t_epi; % scale correctly to reflect displacement in original grid
clear t_x X_reg_epi Y_reg_epi Z_reg_epi x_epi y_epi z_epi


% update rigid param
rigidX_sImg(1:3) = rigidX_sImg(1:3).*reg_epi_res; % in mm
rigid_reg_grid(permute_vec(1:3)) = rigid_reg_grid;
origin_vox_non_rigid = ones(4,1);
for k = 1:length(rigid_reg_grid)
   origin_vox_non_rigid(k) = (rigid_reg_grid{k}(2) + rigid_reg_grid{k}(end))/2;
end
origin_mm_non_rigid = T_epi*origin_vox_non_rigid;
[M_final, origin_final] = composeRigidTrans(M_world, origin, rigidX_sImg, origin_mm_non_rigid(1:3), ...
   opts.phase_encode_direction);


fprintf('\nSaving output files...')

% Rigid transform MPRAGE 
fileout = [remove_extension(struct_output_filename) '.headr_transform.nii.gz'];
affine_transform_nii(struct_filename, inv(M_final), origin_final, fileout);
interp3_nii(fileout, epi_b0_filename, struct_output_filename, 'linear');

% Update rigid transformation
static_filename = epi_output_filename;
moving_filename = struct_filename;
M_world = M_final;
origin = origin_final; 
save([remove_extension(struct_output_filename) '.rigid_registration_result.mat'], ...
   'M_world', 'origin', 'moving_filename', 'static_filename');

% transform original image and save
fprintf('\nCorrecting EPI images for distortion...')
epi_orig = load_untouch_nii_gz(epi_filename, true, workDir);
[x_g, y_g, z_g] = ndgrid(1:size(epi_orig.img,1), 1:size(epi_orig.img,2), 1:size(epi_orig.img,3));

switch opts.phase_encode_direction
   case {'x', 'x-'}
      x_g = x_g + t_epi;
      if opts.intensity_correct
         [~, crct_grd] = gradient(t_epi);
      end
   case {'y', 'y-'}
      y_g = y_g + t_epi;
      if opts.intensity_correct
         [crct_grd] = gradient(t_epi);
      end
   case {'z', 'z-'}
      z_g = z_g + t_epi;
      if opts.intensity_correct
         [~, ~, crct_grd] = gradient(t_epi);
      end
   otherwise
      error('options.phase_encode_direction has to be x/x-/y/y-/z/z-');
end

moving_out = epi_orig;
moving_out.img = zeros(size(epi_orig.img));
cpb = ConsoleProgressBar(); % Set progress bar parameters
cpb.setMinimum(1); cpb.setMaximum(size(epi_orig.img,4)); cpb.start();
for k = 1:size(epi_orig.img,4)
   moving_out.img(:,:,:,k) = interpn(double(epi_orig.img(:,:,:,k)), x_g, y_g, z_g, 'cubic', 0);
   if opts.intensity_correct
      moving_out.img(:,:,:,k) = moving_out.img(:,:,:,k) .* (1+crct_grd);
   end
   msg = sprintf('%d/%d volumes done', k, size(epi_orig.img,4));
   cpb.setValue(k); cpb.setText(msg);
end
cpb.stop();


% save corrected data
moving_out.hdr.dime.scl_slope = 0;
moving_out.hdr.dime.scl_inter = 0;
save_untouch_nii_gz(moving_out, epi_output_filename, workDir);

% save distortion map
moving_out.img = t_epi;
moving_out.hdr.dime.dim(1) = 3;
moving_out.hdr.dime.dim(5) = 1;
save_untouch_nii_gz(moving_out, [remove_extension(epi_output_filename) '.distortion.map.nii.gz'], 64, workDir);

rmdir(workDir, 's');

end

function [M_final, ref_loc] = composeRigidTrans(M_rigid, origin_rigid, EPI_rigidX, EPI_origin, PED_dir)
% Combines rigid-registration results with refinements during EPI-correction. 

SO_rigid = eye(4);
SO_rigid(1:3,end) = -1*origin_rigid(:);

SO_EPI = eye(4);
SO_EPI(1:3,end) = -1*EPI_origin(:);

M_EPI = par2affineMat(EPI_rigidX(1:3), EPI_rigidX(4:6));

% permutation matrix corresponding to PED
switch PED_dir
   case {'x', 'x-'}
      P = eye(4);
   case {'y', 'y-'}
      P = [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
   case {'z', 'z-'}
      P = [0 0 1 0; 0 1 0 0; 1 0 0 0; 0 0 0 1];
end

M_all = inv(SO_EPI)*P*inv(M_EPI)*(P')*SO_EPI * inv(SO_rigid)*inv(M_rigid)*SO_rigid;
M_final = inv(M_all);
ref_loc = zeros(3,1);
end

function [bb_1, bb_2, bb_3, mod_msk] = trim_msk_PED(bb_1, bb_2, bb_3, msk_and, phase_encode_direction, pad_num)
% Trims voxels orthogonal to PED, which are not present in both of the masks. The idea is to remove
% slices orthogonal to PED for which no anatomical information is present (otherwise INVERSION will
% try to squish voxels to match background!). 
%
% bb_1, bb_2, bb_3, must be bounding box of OR of masks.

mod_msk = false(size(msk_and));
[b1, b2, b3] = find_bounding_box(msk_and);

switch phase_encode_direction
   case {'x', 'x-'}
      M = any(msk_and>0, 1);
      [m1, m2, m3] = find_bounding_box(M, pad_num);
      m1 = bb_1;
      mod_msk(:, b2, b3) = true;
      
   case {'y', 'y-'}
      M = any(msk_and>0, 2);
      [m1, m2, m3] = find_bounding_box(M, pad_num);
      m2 = bb_2;
      mod_msk(b1, :, b3) = true;
      
   case {'z', 'z-'}
      M = any(msk_and>0, 3);
      M = permute(M, [3 1 2]); % hack to make M 3D vol
      [m3, m1, m2] = find_bounding_box(M, pad_num);
      m3 = bb_3;
      mod_msk(b1, b2, :) = true;
end

if length(m1)<length(bb_1)
   bb_1 = m1;
end

if length(m2)<length(bb_2)
   bb_2 = m2;
end

if length(m3)<length(bb_3)
   bb_3 = m3;
end
end

function [bb_1, bb_2, bb_3] = pad_bbox_PED(bb_1, bb_2, bb_3, sizeImg, phase_encode_direction, pad_num)
% pads voxels along PED

switch phase_encode_direction
   case {'x', 'x-'}
      bb_1 = [bb_1(1)-pad_num:bb_1(1)-1 bb_1 bb_1(end)+1:bb_1(end)+pad_num];
      bb_1(bb_1<=0) = [];
      bb_1(bb_1>sizeImg(1)) = [];
      
   case {'y', 'y-'}
      bb_2 = [bb_2(1)-pad_num:bb_2(1)-1 bb_2 bb_2(end)+1:bb_2(end)+pad_num];
      bb_2(bb_2<=0) = [];
      bb_2(bb_2>sizeImg(2)) = [];
      
   case {'z', 'z-'}
      bb_3 = [bb_3(1)-pad_num:bb_3(1)-1 bb_3 bb_3(end)+1:bb_3(end)+pad_num];
      bb_3(bb_3<=0) = [];
      bb_3(bb_3>sizeImg(3)) = [];
      
   otherwise
      error('Unknown phase_encode_direction')
end

end

function [epiImg, m_sImg, m_epi_less_csf, sImg, m_epi, Xg, Yg, Zg, Xe, Ye, Ze] = resliceForNonUniformityCorr( ...
   epi_data, epi_mask, epi_mask_less_csf, struct_data, struct_mask, epi_res, bcorr_res)
% Simple reslicing for epi_data. Assumes epi_data is in RAS.

T = diag([epi_res(:); 1]);
enii = epi_data;
enii.hdr.hist.srow_x = T(1,:);
enii.hdr.hist.srow_y = T(2,:);
enii.hdr.hist.srow_z = T(3,:);
[~, Xe, Ye, Ze] = get_original_grid_data(enii);


enii2 = enii;
enii.img = struct_mask.img; 
enii2.img = epi_mask.img; 
[m_sImg, m_epi, Xg, Yg, Zg] = reslice_vol_same_grid(enii, enii2, bcorr_res, 'linear');

enii.img = epi_data.img;
epiImg = myreslice_nii_match_res(enii, 'linear', Xg, Yg, Zg);

enii.img = struct_data.img;
sImg = myreslice_nii_match_res(enii, 'linear', Xg, Yg, Zg);

enii.img = epi_mask_less_csf.img;
m_epi_less_csf = myreslice_nii(enii, 'linear', Xg, Yg, Zg);
end

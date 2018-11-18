function transform_data_affine(data_file, data_coord, output_file, ...
   reg_moving_file, reg_static_file, reg_mat_file, method)
% General function to transform nii/dfs using affine matrix estimated using
% register_files_affine.m. No sophisticated checks are performed for
% file-coordinate checks. 
%
%  data_file - filename of the data file to be transformed. It could be
%              a volume (.nii/.nii.gz) or surface (.dfs/.dfs.gz)
%  data_coord - Possible values: mov/m/static/s. Indicates the
%               coordinate-space of the data_file. 'moving' or 'm' represents
%               that the data_file is in the same coordinate space as the
%               reg_moving_file. Similarly, 'static' or 's' indicates same
%               coordinate space as reg_static_file
%  output_file - Target output filename
%  reg_moving_file - Moving file used during registration
%  reg_static_file - Static file used during registration
%  reg_mat_file - Mat file saved as output after affine regitration
%  method - (Optional) interpolation method. 'linear' when not specified.
%           Ignored when used with surface files. 
%

workDir = tempname();
mkdir(workDir);

% check image headers
[~, reg_moving_file] = check_nifti_file(reg_moving_file, workDir);
[~, reg_static_file] = check_nifti_file(reg_static_file, workDir);


% detect file type
[~, ext] = remove_extension(data_file);
if ismember(ext, {'.nii', '.nii.gz'})
   filetype = 'vol';
   [~, data_file] = check_nifti_file(data_file, workDir);
   
elseif ismember(ext, {'.dfs', '.dfs.gz'})
   filetype = 'surf';
else
   error('Unrecognized or unsupported format of data file %s\n', ext)
end

if ~exist('method', 'var')
   method = 'linear';
end

% set transformation direction
param = load(reg_mat_file);
if lower(data_coord(1))=='m' % moving to static
   M_mat = inv(param.M_world);
   ref_loc = param.origin;
   native_file = reg_moving_file;
   target_file = reg_static_file;
   
elseif lower(data_coord(1))=='s' % static to moving
   M_mat = param.M_world;
   ref_loc = param.origin;
   native_file = reg_static_file;
   target_file = reg_moving_file;
   
else
   error('Unrecognized data_coord: %s\n', data_coord)
end

% interpolate
if filetype(1)=='v' % volume
   
   % basic header check
   if matchNiiHeader(data_file, native_file)
      msg = bdp_linewrap(['data_file does not match specified coordinate. It can be because of '...
         'differences in nifti header. Please check these files: ' escape_filename(data_file) ' and ' ...
         escape_filename(native_file)]);
      error(msg);
   end
   
   temp_file = [tempname() '.nii.gz'];
   affine_transform_nii(data_file, M_mat, ref_loc, temp_file);
   interp3_nii(temp_file, target_file, output_file, method);
   delete(temp_file);
   
elseif filetype(1)=='s' % surface
   affine_transform_dfs(data_file, M_mat, ref_loc, output_file, native_file, target_file);
   
end

fprintf('\nFinished transformation. Saved file to disk: %s\n', escape_filename(output_file));
end


function out = matchNiiHeader(file1, file2)
% simple checks for consistency of nii volumes. Returns true when niis are
% consistent

out = 0;
nii1 = load_untouch_nii_gz(file1);
nii2 = load_untouch_nii_gz(file2);


sz1 = size(nii1.img); 
sz2 = size(nii2.img);
if ~isequal(sz1(1:3), sz2(1:3))
   out = 1; % image size don't match
   return;
end

if norm(nii1.hdr.dime.dim(2:4) - nii2.hdr.dime.dim(2:4)) > 1e-4
   out = 2; % Dimensions don't match
   return;
end

if norm(nii1.hdr.dime.pixdim(2:4)-nii2.hdr.dime.pixdim(2:4)) > 1e-3
   out = 3; % voxel size don't match
   return;
end

nx = norm(nii1.hdr.hist.srow_x - nii2.hdr.hist.srow_x) > 1e-2;
ny = norm(nii1.hdr.hist.srow_y - nii2.hdr.hist.srow_y) > 1e-2;
nz = norm(nii1.hdr.hist.srow_z - nii2.hdr.hist.srow_z) > 1e-2;

if nx || ny || nz
   out = 4; % sform does not match
   return;
end


end

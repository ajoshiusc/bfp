function [nii_in, reorient_matrix, sform_new, bMatrice_new] = reorient_nifti_sform(nifti_file, output_file, bmatrices)
% Reorients data matrix by applying all orthogonal rotation (& flipping) component of sform
% matrix (S). It tries to get modified sform matrix (Sn) as close as possible to identity
% matrix (before scaling). The output file should look same as in input file in a nifti viewer
% (except with minor differences due to nifti viewer's interpolation errors). It also adds
% appropriate q-form matrix after re-orientation. 
%
% This function does NOT perform any resampling or interpolation. 
%
%  nifti_file - Input nifti file with full path
%  output_file - Output nifti file with full path
%  bmatrices - (optional) bmat file or bmatrices for diffusion data, in voxel coordinate
%
%
%   Xw = S * Xi
%      = (Tr*Rf*Sc) * Xi
%      = (Tr*Rf*Sc) * (inv(Ro) * Ro) * Xi
%      = (Tr*Rf*Sc*inv(Ro)) * (Ro * Xi)
%      = Sn * Xn
%
%   Xw - world coordinate 
%   Xi - voxel coordinate (index starts from zero)
%   S  - Initial Sform matrix (4x4)
%   Xn - new voxel coordinate, after applying transforms; Xn = Ro*Xi
%   Sn - New Sform matrix; Sn = Tr*Rf*Sc*inv(Ro)
%   Tr - translation component of S
%   Rf - Rotation (& fliping) component of S
%   Sc - Scaling component of S
%   Ro - Re-orientation matrix which makes Sn as close to I as possible (before scaling)
%

output_file = [remove_extension(output_file) '.nii.gz'];
nii_in = load_untouch_nii_gz(nifti_file);

if nii_in.hdr.hist.sform_code ==0
   error('sform_code is unset! Try using add_sform first.')
end

sformT = zeros([4 4]);
sformT(1,:) = nii_in.hdr.hist.srow_x;
sformT(2,:) = nii_in.hdr.hist.srow_y;
sformT(3,:) = nii_in.hdr.hist.srow_z;
sformT(4,4) = 1;

% Seperate translation from transform
Mshift = eye(4);
Mshift(:,4) = sformT(:,4);

% normalize the matrix (to make abs determinant equal to 1)
Mrest = sformT(1:3, 1:3)*diag(1./abs(nii_in.hdr.dime.pixdim(2:4)));

% check det of Mrest
if is_canonical_mat(Mrest) % already aligned to canonical coordinates   
   [nii_in] = save_nii_with_qform(nii_in, output_file);
   
   reorient_matrix = eye(3);
   sform_new = sformT;
   
   if exist('bmatrices', 'var')
      bMatrice_new = bmatrices;
      write_bmat_file(bMatrice_new, [remove_extension(output_file) '.bmat']);
   end
   return;
end

% check if Mrest is rotation matrix 
if abs(abs(det(Mrest))-1)>1e-3
   sformT
   nii_in.hdr.dime.pixdim
   error('Absolute determinant of transformation matrix is not 1. Something is wrong!')
end


% find closest orthogonal coordinate
reorient_matrix = zeros(3,3);
B = sort(abs(Mrest(:)),'descend');
max_mask = abs(Mrest)>=B(3);
if sum(max_mask(:))>3
   sformT
   max_mask
   nii_in.hdr.dime.pixdim
   error('Something is wrong with sform matrix !')
end
temp = Mrest(max_mask);
reorient_matrix(max_mask) = temp./abs(temp);


% check orient_mat
if abs(abs(det(reorient_matrix))-1)>1e-3
   sformT
   reorient_matrix
   nii_in.hdr.dime.pixdim
   error('Absolute determinant of orientation matrix is not 1. Something is wrong!')
end

if is_canonical_mat(reorient_matrix) && abs(det(Mrest)-1)<1e-3  % already aligned to canonical coordinates
   [nii_in] = save_nii_with_qform(nii_in, output_file);
   sform_new = sformT;
   
   if exist('bmatrices', 'var')
      bMatrice_new = bmatrices;
      write_bmat_file(bMatrice_new, [remove_extension(output_file) '.bmat']);
   end
   return;
   
elseif is_canonical_mat(reorient_matrix) && abs(det(Mrest)-1)>1e-3
   sformT
   Mrest
   reorient_matrix
   nii_in.hdr.dime.pixdim
   error('Something went wrong with sform matrix!')
end


% Find residual matrix
Rot_residue = sformT(1:3, 1:3)*inv(reorient_matrix);
Mrest_rot_residue = Mrest*inv(reorient_matrix); % independent of resolution 
if abs(det(Mrest_rot_residue)-1)>1e-3
   sformT
   reorient_matrix
   Rot_residue
   Mrest_rot_residue
   nii_in.hdr.dime.pixdim
   error('Residue matrix is not rotation matrix. Something went wrong!')
end


% apply orientation transform to image matrix
%nii_out = nii_in;
%clear nii_in;

% First apply the permutation part of transform
temp = abs(reorient_matrix)'; % make row major
ind = find(temp);
permute_order = 1:ndims(nii_in.img);
permute_order(1:3) = ind(:)' - [0 3 6];
nii_in.img = permute(nii_in.img, permute_order);

temp = abs(nii_in.hdr.dime.pixdim(2:4));
nii_in.hdr.dime.pixdim(2:4) = temp(permute_order(1:3));

temp = abs(nii_in.hdr.dime.dim(2:4));
nii_in.hdr.dime.dim(2:4) = temp(permute_order(1:3));

% next apply the flipping for each dimension
orient_mat_diag = reorient_matrix*abs(reorient_matrix)';
if orient_mat_diag(1,1)<0
   nii_in.img = flipdim(nii_in.img, 1);
end

if orient_mat_diag(2,2)<0
   nii_in.img = flipdim(nii_in.img, 2);
end

if orient_mat_diag(3,3)<0
   nii_in.img = flipdim(nii_in.img, 3);
end

% get final sform matrix
size_data = size(nii_in.img);
trans_new = Rot_residue*((size_data(1:3)-1)'.*(diag(orient_mat_diag)<0));
Mshift(1:3,4) = Mshift(1:3,4) - trans_new;

Rot_residue(4,4) = 1;
sform_new = Mshift * Rot_residue;

nii_in.hdr.hist.srow_x = sform_new(1,:);
nii_in.hdr.hist.srow_y = sform_new(2,:);
nii_in.hdr.hist.srow_z = sform_new(3,:);

nii_in = save_nii_with_qform(nii_in, output_file);

% reorient bmat and save
if exist('bmatrices', 'var')
   if ischar(bmatrices)
      bmatrices = readBmat(bmatrices);
   end
   
   bMatrice_new =  zeros(size(bmatrices));
   for iDir = 1:size(bmatrices,3)
      bMatrice_new(:,:,iDir) = reorient_matrix*bmatrices(:,:,iDir)*(reorient_matrix');
   end   
   writeBmatFile(bMatrice_new, [remove_extension(output_file) '.bmat']);
end


end

% check for char/matrix and then save
function write_bmat_file(bMatrices, fname)
if ischar(bMatrices)
   if ~isequal(bMatrices, fname)
      copyfile(bMatrices, fname);
   end
else
   writeBmatFile(bMatrices, fname);
end
end


% check for canonical matrix
function b = is_canonical_mat(M)

if isequal(diag(abs(diag(M))),M)
   b = true;
else
   b = false;
end

end


function nii = save_nii_with_qform(nii, output_file)
% Adds qform structure to nii structure before saving. 
% * Overrides any qform parameters in header!
% * Assumes that sform is defined and is valid. 
% * Assumes that nii is in RAS after reorienting (probaly works without
%     this too). 

nii.hdr.hist.qform_code = 0;
nii = add_qform(nii);
save_untouch_nii_gz(nii, output_file);

end



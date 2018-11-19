function [vol_nii, Tnew] = affine_transform_nii(v1, F, origin_loc, outFile)
% Applies affine transform to nifti volume. As the location of voxels in nifti files are explicitly
% defined by an affine matrix (sform), any affine transform can be completely applied just by
% transforming the header affine matrix (without changing any voxel intensity). This function just
% changes the sform header of the volume appropriately. No reslicing/interpolation. Requires valid
% sform matrix/header. 
%
% v1 - nifti file name or nifti structure
% F - 4x4 affine matrix transformation
% origin_loc - 3x1 vector representing origin. See note below about origin.
% outfile - Optional, output file name
%
% Notes about Origin (& assumptions):
%   The sampling points in a nifti volume lies on a coordinate system defined by sform-matrix
%   in the nifti header. This coordinate system is often referred to as real-world coordinates.
%   See get_original_grid_data.m for details of getting the real-world coordinate.
%   It is assumed that the affine matrix 'F' corresponds to a coordinate system which we get
%   after shifting (no rotation) the real-world coordinate system by 'origin_loc'. Following
%   equations describe the transformation: 
%         Xw  = T*Xi
%         Xwn = Sinv*F*S*Xw
%             = Sinv*F*S*T*Xi
%             = Tnew*Xi
%
%    Xi - voxel coordinate 
%    T  - original (sform) affine matrix in nifti file
%    Xw - world coordinate of original nifti image
%    Xwn - world coordinate of the nifti image after affine transformation 
%    S  - translation matrix (4x4) which transfers origin to ORIGIN_LOC
%    F  - Affine transform to be applied
%    Sinv - inv(S); Puts back the origin from ORIGIN_LOC to [0; 0; 0]
%    Tnew - New sform matrix after affine transformation 
%

if ischar(v1)
   vol_nii = load_untouch_nii_gz(v1);
else
   if ~isfield(v1,'untouch') || v1.untouch ~= 1
      error('Please use ''load_untouch_nii.m'' to load the structure.');
   end
   vol_nii = v1;
   clear v1
end

if vol_nii.hdr.hist.sform_code<=0 
   error('sform_code is not set in the nifti header. Add sform matrix first!');
end

if length(origin_loc)~= 3
   error('ref_loc much be of length 3.');
else
   S = eye(4);
   S(1:3,end) = -1*origin_loc(:);
end

Tvol1 = eye(4);
Tvol1(1,:) = vol_nii.hdr.hist.srow_x(1:4);
Tvol1(2,:) = vol_nii.hdr.hist.srow_y(1:4);
Tvol1(3,:) = vol_nii.hdr.hist.srow_z(1:4);

Tnew = inv(S)*F*S*Tvol1;

vol_nii.hdr.hist.srow_x(1:4) = Tnew(1,:);
vol_nii.hdr.hist.srow_y(1:4) = Tnew(2,:);
vol_nii.hdr.hist.srow_z(1:4) = Tnew(3,:);

% unset qform 
vol_nii.hdr.hist.qform_code = 0;

if exist('outFile','var')
   save_untouch_nii_gz(vol_nii, outFile);
end

end


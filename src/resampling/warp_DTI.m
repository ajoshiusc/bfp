% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2019 The Regents of the University of California and the University of Southern California
% Created by Anand A. Joshi, Chitresh Bhushan, David W. Shattuck, Richard M. Leahy
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


function [L1out, L2out, L3out, W1, W2, W3] = warp_DTI(mapDTI, L1, L2, L3, V1, V2, V3, out_basename)
%
% Warps DTI data using voxel-to-voxel map.
%  - Assumes that tensor values are rotated in voxel coordinate (see
%      estimate_tensors.m)
%  - Uses the information from nifti header (sform) from all the files.
%  - Assumes that all eigen files have same nifti header (does not check)
%

if ischar(L1)
    L1 = load_untouch_nii_gz(L1);
    L2 = load_untouch_nii_gz(L2);
    L3 = load_untouch_nii_gz(L3);
    V1 = load_untouch_nii_gz(V1);
    V2 = load_untouch_nii_gz(V2);
    V3 = load_untouch_nii_gz(V3);
    mapDTI = load_untouch_nii_gz(mapDTI);
end


[len bre dep dim4] = size(mapDTI.img);

% sform transformation
Tmap = zeros([3 3]);
Tmap(1,:)= mapDTI.hdr.hist.srow_x(1:3);
Tmap(2,:)= mapDTI.hdr.hist.srow_y(1:3);
Tmap(3,:)= mapDTI.hdr.hist.srow_z(1:3);

Tdti = zeros([3 3]);
Tdti(1,:)= L1.hdr.hist.srow_x(1:3);
Tdti(2,:)= L1.hdr.hist.srow_y(1:3);
Tdti(3,:)= L1.hdr.hist.srow_z(1:3);

if ~isequal(Tmap, diag(diag(Tmap)))
    Tmap
    warning('Tmap is not diagonal matrix. Tensor reorientation may be wrong!\nContinuing anyway...')
end

if ~isequal(Tdti, diag(diag(Tdti)))
    Tdti
    warning('Tdti is not diagonal matrix. Tensor reorientation may be wrong!\nContinuing anyway...')
end

% Map in the local coordinate
% c(1,:) = reshape(mapDTI.img(:,:,:,1), 1, [ ]);
% c(2,:) = reshape(mapDTI.img(:,:,:,2), 1, [ ]);
% c(3,:) = reshape(mapDTI.img(:,:,:,3), 1, [ ]);
% c(4,:) = 1;
% c = (Tmap\Tdti)*c;
% mapDTI_local(:,:,:,1) = reshape(c(1,:), [len bre dep]);
% mapDTI_local(:,:,:,2) = reshape(c(2,:), [len bre dep]);
% mapDTI_local(:,:,:,3) = reshape(c(3,:), [len bre dep]);
% clear c

V1out = mapDTI;
V1out.img = [];
V2out = V1out;
V3out = V1out;

L1out = mapDTI;
L1out.img = [];
L1out.hdr.dime.dim(1) = 3;
L1out.hdr.dime.dim(5) = 1;
L2out = L1out;
L3out = L1out;

mapDTI.img = double(mapDTI.img);
% Interpolate
L1out.img = interp3(double(L1.img), mapDTI.img(:,:,:,2), mapDTI.img(:,:,:,1), mapDTI.img(:,:,:,3), 'nearest');
L2out.img = interp3(double(L2.img), mapDTI.img(:,:,:,2), mapDTI.img(:,:,:,1), mapDTI.img(:,:,:,3), 'nearest');
L3out.img = interp3(double(L3.img), mapDTI.img(:,:,:,2), mapDTI.img(:,:,:,1), mapDTI.img(:,:,:,3), 'nearest');

for m = 1:dim4
    V1out.img(:,:,:,m) = interp3(double(V1.img(:,:,:,m)), mapDTI.img(:,:,:,2), mapDTI.img(:,:,:,1), mapDTI.img(:,:,:,3), 'nearest');
    V2out.img(:,:,:,m) = interp3(double(V2.img(:,:,:,m)), mapDTI.img(:,:,:,2), mapDTI.img(:,:,:,1), mapDTI.img(:,:,:,3), 'nearest');
    V3out.img(:,:,:,m) = interp3(double(V3.img(:,:,:,m)), mapDTI.img(:,:,:,2), mapDTI.img(:,:,:,1), mapDTI.img(:,:,:,3), 'nearest');
end

L1out.img(~isfinite(L1out.img)) = 0;
L2out.img(~isfinite(L2out.img)) = 0;
L3out.img(~isfinite(L3out.img)) = 0;
V1out.img(~isfinite(V1out.img)) = 0;
V2out.img(~isfinite(V2out.img)) = 0;
V3out.img(~isfinite(V3out.img)) = 0;



% Compute Jacobian
J1 = zeros([len bre dep dim4]);
J2 = zeros([len bre dep dim4]);
J3 = zeros([len bre dep dim4]);
res = mapDTI.hdr.dime.pixdim(2:4);

for m = 1:dim4
    [temp1 temp2 temp3] = gradient(mapDTI.img(:,:,:,m), res(1), res(2), res(3));
    J2(:,:,:,m) = temp1;   % order of J1 & j2 are interchanged to accomodate the matlab's convention with ndgrid. check gradient & ndgrid.
    J1(:,:,:,m) = temp2;
    J3(:,:,:,m) = temp3;
end
clear temp1 temp2 temp3 mapDTI

% PPD to rotate eigenvectors
[W1 W2 W3] = PPD_FAST(V1out, V2out, V3out, J1, J2, J3);

if exist('out_basename', 'var')
    save_untouch_nii_gz(L1out, [out_basename '.L1.nii.gz']);
    save_untouch_nii_gz(L2out, [out_basename '.L2.nii.gz']);
    save_untouch_nii_gz(L3out, [out_basename '.L3.nii.gz']);
    save_untouch_nii_gz(W1, [out_basename '.V1.nii.gz']);
    save_untouch_nii_gz(W2, [out_basename '.V2.nii.gz']);
    save_untouch_nii_gz(W3, [out_basename '.V3.nii.gz']);
end


end


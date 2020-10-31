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

function [W1 W2 W3] = PPD_FAST(V1, V2, V3, G1, G2, G3, ~)
%PPD Rotates the eigenvectors of the diffusion matrix with Preservation of
%Principle component(PPD) algorithm & generates new eigenvectors. Saves the
%output in .nii.gz file.
%
%   V1 - Principal eigenvector (corresponding to eigenvalue L1)
%   V2 - Second eigenvector (corresponding to eigenvalue L2)
%   V3 - Third eigenvector (corresponding to eigenvalue L3)
%       such that L1>L2>L3
%
%   F1, F2, F3 - Rows of gradient matrix (of the map)
%
%   W1 - Rotated principal eigenvector
%   W2 - Rotated second eigenvector
%   W3 - Rotated third eigenvector
%

%fprintf('\n==== Applying PPD ====\n');

[len bre dep dim4] = size(V1.img);

% Invert the gradient matrix
% detG =   G1(:,:,:,1) .* G3(:,:,:,3) .* G2(:,:,:,2) ...
%        - G1(:,:,:,1) .* G3(:,:,:,2) .* G2(:,:,:,3) ...
%        - G2(:,:,:,1) .* G3(:,:,:,3) .* G1(:,:,:,2) ...
%        + G2(:,:,:,1) .* G3(:,:,:,2) .* G1(:,:,:,3) ...
%        + G3(:,:,:,1) .* G2(:,:,:,3) .* G1(:,:,:,2) ...
%        + G3(:,:,:,1) .* G2(:,:,:,2) .* G1(:,:,:,3);
%   
% F1(:,:,:,1) = (G3(:,:,:,3).*G2(:,:,:,2) - G3(:,:,:,2).*G2(:,:,:,3))./detG ;
% F1(:,:,:,2) = (G3(:,:,:,2).*G1(:,:,:,3) - G3(:,:,:,3).*G1(:,:,:,2))./detG ;
% F1(:,:,:,3) = (G2(:,:,:,3).*G1(:,:,:,2) - G2(:,:,:,2).*G1(:,:,:,3))./detG ;
% 
% F2(:,:,:,1) = (G3(:,:,:,1).*G2(:,:,:,3) - G3(:,:,:,3).*G2(:,:,:,1))./detG ;
% F2(:,:,:,2) = (G3(:,:,:,3).*G1(:,:,:,1) - G3(:,:,:,1).*G1(:,:,:,3))./detG ;
% F2(:,:,:,3) = (G2(:,:,:,1).*G1(:,:,:,3) - G2(:,:,:,3).*G1(:,:,:,1))./detG ;
% 
% F3(:,:,:,1) = (G3(:,:,:,2).*G2(:,:,:,1) - G3(:,:,:,1).*G2(:,:,:,2))./detG ;
% F3(:,:,:,2) = (G3(:,:,:,1).*G1(:,:,:,2) - G3(:,:,:,2).*G1(:,:,:,1))./detG ;
% F3(:,:,:,3) = (G2(:,:,:,2).*G1(:,:,:,1) - G2(:,:,:,1).*G1(:,:,:,2))./detG ;

F1=G1; F2=G2; F3=G3;
clear G1 G2 G3

% convert to double for higher precision
V1.img = double(V1.img);
V2.img = double(V2.img);
V3.img = double(V3.img);


% Normalizing V1, V2, V3
V1.img = norm_vec(V1.img);
V2.img = norm_vec(V2.img);
V3.img = norm_vec(V3.img);


% applying Jacobian to V1
n1 = zeros([len bre dep dim4]);
n1(:,:,:,1) = sum(F1.*V1.img, 4);
n1(:,:,:,2) = sum(F2.*V1.img, 4);
n1(:,:,:,3) = sum(F3.*V1.img, 4);
n1 = norm_vec(n1);

% applying Jacobian to V2
n2 = zeros([len bre dep dim4]);
n2(:,:,:,1) = sum(F1.*V2.img, 4);
n2(:,:,:,2) = sum(F2.*V2.img, 4);
n2(:,:,:,3) = sum(F3.*V2.img, 4);
n2 = norm_vec(n2);

cosTheta1 = sum(V1.img.*n1, 4);
axis1 = cross_vec(V1.img, n1);

% Projection of n2 perpendicular to n1
temp = zeros(size(n1));
n1_dot_n2 = sum(n1.*n2, 4);
temp(:,:,:,1) = n1(:,:,:,1).*n1_dot_n2;
temp(:,:,:,2) = n1(:,:,:,2).*n1_dot_n2;
temp(:,:,:,3) = n1(:,:,:,3).*n1_dot_n2;
proj_n2 = n2 - temp;
proj_n2 = norm_vec(proj_n2);
clear temp n2 n1_dot_n2

V2rot = rot_vec(V2.img, axis1, cosTheta1);

cosTheta2 = sum(V2rot.*proj_n2, 4);
axis2 = cross_vec(V2rot, proj_n2);


W1 = V1;
W1.img = n1;

W2 = V2;
W2.img = rot_vec(V2rot, axis2, cosTheta2);

W3 = V3;
V3rot = rot_vec(V3.img, axis1, cosTheta1);
W3.img = rot_vec(V3rot, axis2, cosTheta2);

end


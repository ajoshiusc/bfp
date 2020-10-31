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

function [ v_rot ] = rot_vec( v, axis, cosTheta )
% Rotates vector v about axis by angle theta. Function takes cosine of
% theta as argument. axis may not be a unit vector. Here vectors are stored
% in 4th dimesion of matrix a & b (usually the eigenvectors) 
% 
%  Uses Rodrigues' rotation formula

axis = norm_vec(axis);

term1 = zeros(size(v));
term1(:,:,:,1) = v(:,:,:,1).*cosTheta;
term1(:,:,:,2) = v(:,:,:,2).*cosTheta;
term1(:,:,:,3) = v(:,:,:,3).*cosTheta;


term2 = zeros(size(v));
sinTheta = sqrt(1 - (cosTheta.^2));
temp = cross_vec(axis, v);
term2(:,:,:,1) = temp(:,:,:,1).*sinTheta;
term2(:,:,:,2) = temp(:,:,:,2).*sinTheta;
term2(:,:,:,3) = temp(:,:,:,3).*sinTheta;
clear temp sinTheta


term3 = zeros(size(v));
temp = sum(axis.*v, 4).*(1-cosTheta);
term3(:,:,:,1) = axis(:,:,:,1).*temp;
term3(:,:,:,2) = axis(:,:,:,2).*temp;
term3(:,:,:,3) = axis(:,:,:,3).*temp;

v_rot = term1 + term2 + term3;

end


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

function [ v ] = cross_vec( a, b )
% Computes cross product axb, where vectors are stored in 4th
% dimesion of matrix a & b (usually the eigenvectors)
%
% a, b - 4D matrix (n x m x k x 3)
%
% Should be faster than matlab implementation in some cases.

v = zeros(size(a));

v(:,:,:,1) = (a(:,:,:,2).*b(:,:,:,3)) - (a(:,:,:,3).*b(:,:,:,2));
v(:,:,:,2) = (a(:,:,:,3).*b(:,:,:,1)) - (a(:,:,:,1).*b(:,:,:,3));
v(:,:,:,3) = (a(:,:,:,1).*b(:,:,:,2)) - (a(:,:,:,2).*b(:,:,:,1));

end


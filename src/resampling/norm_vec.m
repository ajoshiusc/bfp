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

function [ v_out ] = norm_vec( v )
%NORM_VEC Returns normalized vectors, where vectors are stored in 4th
%dimesion of a matrix v (usually the eigenvectors)

v_out = zeros(size(v));
v_norm = sqrt(sum(v.^2, 4));

v_out(:,:,:,1) = v(:,:,:,1) ./ v_norm;
v_out(:,:,:,2) = v(:,:,:,2) ./ v_norm;
v_out(:,:,:,3) = v(:,:,:,3) ./ v_norm;

v_out(isnan(v_out)) = 0;

end


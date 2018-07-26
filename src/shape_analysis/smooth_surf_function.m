% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2016 The Regents of the University of California and the University of Southern California
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


function f=smooth_surf_function(s,f0,a1,a2,aniso,normalize)
%function f=smooth_surf_function(s,f0,a1,a2,aniso)
%Performs smoothing of a function on a surface
%s : surface, f0 : input function defined at every vertex
% a1 a2: smoothing parameters for Laplacian and Gradient respectively
% recommended to be set equal to eachother
% aniso : anisotropy term for anisotropic diffusion

if ~exist('a1','var')
    a1=3.1;
end
if ~exist('a2','var')
    a2=3.1;
end
if ~exist('aniso','var')
    aniso=ones(size(s.vertices,1),1);
end
if ~exist('normalize','var')
    normalize=0;
end
%||f-f0||.^2 + a1*||Lf||^2 + a2*||grad f||^2

%L=loreta(s);

[S,Dx,Dy]=get_stiffness_matrix_tri_wt(s,aniso);

M=[a1*S;a2*[Dx;Dy]];
A=[speye(length(s.vertices));M]; b=[f0;zeros(size(M,1),1)];

f=mypcg(A'*A,A'*b,1e-100,3000,diag(A'*A));

if normalize>0
    f=f*norm(f0)/norm(f);
end

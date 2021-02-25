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


function O_trans = bspline_grid_generate_deformation(spacing, sizeI)
% Creates grid of control points for to parameterize deformation-field using uniform 3D
% b-spline. Output is different from bspline_grid_generate() which parameterized a map! To be
% used with bspline_repeated_endpoint_deformation_3d_double_only_x(). 
%
%  Spacing: vector of length 3, representing space between adjacent control points (in unit of
%           voxels) along three dimensions of the image.
%  sizeI: vector with the sizes of the image for which control point would be generated
%

if(length(spacing)==2)
   error('This function does not support 2D grid points.')
else
   % Determine grid spacing
   dx = spacing(1);
   dy = spacing(2);
   dz = spacing(3);
   
   if(mod(sizeI(1)-1,dx) + mod(sizeI(2)-1,dy) + mod(sizeI(3)-1,dz))~=0
      error('Size and spacing must be exact.');
   end
   
   sizeI = sizeI(:);
   spacing = spacing(:);
   
   O_size = ceil(sizeI./spacing);
   O_trans = zeros([O_size(:); 3]');
end

end


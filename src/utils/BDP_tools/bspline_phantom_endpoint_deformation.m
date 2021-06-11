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


function [t_x, t_y, t_z] = bspline_phantom_endpoint_deformation(O_grid, sizeI, Spacing)
% This function is interface to the mex function
% bspline_phantom_endpoint_deformation_3d_double.mex
% 
% This function adds the phantom points at the end, so that the end control
% points are interpolated. 
%
% Returns the deformation at each voxel in the image. 
%

if(mod(sizeI(1)-1,Spacing(1))+mod(sizeI(2)-1,Spacing(2))+mod(sizeI(3)-1,Spacing(3)))~=0
   error('SizeI and Spacing must be exact');
end

O_grid_ph_x = zeros(size(O_grid(:,:,:,1))+2);
O_grid_ph_y = zeros(size(O_grid(:,:,:,2))+2);
O_grid_ph_z = zeros(size(O_grid(:,:,:,3))+2);

O_grid_ph_x(2:end-1, 2:end-1, 2:end-1) = double(O_grid(:,:,:,1));
O_grid_ph_y(2:end-1, 2:end-1, 2:end-1) = double(O_grid(:,:,:,2));
O_grid_ph_z(2:end-1, 2:end-1, 2:end-1) = double(O_grid(:,:,:,3));

O_grid_ph_x = generate_phantom_grid_point(O_grid_ph_x);
O_grid_ph_y = generate_phantom_grid_point(O_grid_ph_y);
O_grid_ph_z = generate_phantom_grid_point(O_grid_ph_z);


global nthreads;
if isempty(nthreads), nthreads = 4; end

[t_x, t_y, t_z] = bspline_phantom_endpoint_deformation_3d_double(double(O_grid_ph_x), double(O_grid_ph_y), double(O_grid_ph_z), ...
                                    double(sizeI), double(Spacing), double(nthreads));

end



function O_grid_ph = generate_phantom_grid_point(O_grid_ph)
% Generate phantom points to interpolate the end points of actual grid

O_grid_ph(1,:,:) = 2*O_grid_ph(2,:,:)-O_grid_ph(3,:,:);
O_grid_ph(:,1,:) = 2*O_grid_ph(:,2,:)-O_grid_ph(:,3,:);
O_grid_ph(:,:,1) = 2*O_grid_ph(:,:,2)-O_grid_ph(:,:,3);

O_grid_ph(end,:,:) = 2*O_grid_ph(end-1,:,:)-O_grid_ph(end-2,:,:);
O_grid_ph(:,end,:) = 2*O_grid_ph(:,end-1,:)-O_grid_ph(:,end-2,:);
O_grid_ph(:,:,end) = 2*O_grid_ph(:,:,end-1)-O_grid_ph(:,:,end-2);

end

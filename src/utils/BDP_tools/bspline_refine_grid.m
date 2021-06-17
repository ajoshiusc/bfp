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


function [O_new, Spacing] = bspline_refine_grid(O_grid, Spacing, sizeI)
%
%  TEMP workaround - calls the refining code from toolbox - should be
%  checked and implemented correctly.
%
% Refine image transformation grid of 1D b-splines with use of spliting matrix
%       (Lane-Riesenfeld Algorithm - for uniform cubic B-spline)
% Msplit=(1/8)*[4 4 0 0;
%               1 6 1 0;
%               0 4 4 0;
%               0 1 6 1;
%               0 0 4 4];
%
%     [O_new,Spacing] = refine_grid(O_trans,Spacing,sizeI)
%

if(ndims(O_grid)==3)
   error('Not supported')
else
   
   O_grid_ph_x = zeros(size(O_grid(:,:,:,1))+4);
   O_grid_ph_y = zeros(size(O_grid(:,:,:,2))+4);
   O_grid_ph_z = zeros(size(O_grid(:,:,:,3))+4);
   
   O_grid_ph_x(2:end-3, 2:end-3, 2:end-3) = double(O_grid(:,:,:,1));
   O_grid_ph_y(2:end-3, 2:end-3, 2:end-3) = double(O_grid(:,:,:,2));
   O_grid_ph_z(2:end-3, 2:end-3, 2:end-3) = double(O_grid(:,:,:,3));
   
   O_grid_ph_x = generate_phantom_grid_point(O_grid_ph_x);
   O_grid_ph_y = generate_phantom_grid_point(O_grid_ph_y);
   O_grid_ph_z = generate_phantom_grid_point(O_grid_ph_z);
   
   O_trans = cat(4, O_grid_ph_x, O_grid_ph_y, O_grid_ph_z);
   
   [O_temp, Spacing] = refine_grid(O_trans, Spacing, sizeI);
   
   O_new = O_temp(2:end-2, 2:end-2, 2:end-2, :);
   
end

end

function O_grid_ph = generate_phantom_grid_point(O_grid_ph)
% Generate phantom points to interpolate the end points of actual grid

O_grid_ph(1,:,:) = 2*O_grid_ph(2,:,:)-O_grid_ph(3,:,:);
O_grid_ph(:,1,:) = 2*O_grid_ph(:,2,:)-O_grid_ph(:,3,:);
O_grid_ph(:,:,1) = 2*O_grid_ph(:,:,2)-O_grid_ph(:,:,3);

O_grid_ph(end-2,:,:) = 2*O_grid_ph(end-3,:,:)-O_grid_ph(end-4,:,:);
O_grid_ph(:,end-2,:) = 2*O_grid_ph(:,end-3,:)-O_grid_ph(:,end-4,:);
O_grid_ph(:,:,end-2) = 2*O_grid_ph(:,:,end-3)-O_grid_ph(:,:,end-4);

% dummy 2 more rows of points
O_grid_ph(end-1:end,:,:) = O_grid_ph([end-2 end-2],:,:);
O_grid_ph(:,end-1:end,:) = O_grid_ph(:,[end-2 end-2],:);
O_grid_ph(:,:,end-1:end) = O_grid_ph(:,:,[end-2 end-2]);

end



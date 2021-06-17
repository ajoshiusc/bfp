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


function [O_new, Spacing] = bspline_refine_grid_deformation(O_grid, Spacing, sizeI)
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
   sz = size(O_grid);
   
   xx = 0:(sz(1)+2);
   yy = 0:(sz(2)+2);
   zz = 0:(sz(3)+2);
   
   % repeat end points
   xx(xx<1) = 1;
   xx(xx>sz(1)) = sz(1);
   yy(yy<1) = 1;
   yy(yy>sz(2)) = sz(2);
   zz(zz<1) = 1;
   zz(zz>sz(3)) = sz(3);
   
   [xx, yy, zz] = ndgrid(xx, yy, zz);
   ind = sub2ind(sz, xx(:), yy(:), zz(:));
   
   temp = O_grid(:,:,:,1);
   O_trans_x = reshape(temp(ind), sz(1:3)+3);
   
   temp = O_grid(:,:,:,2);
   O_trans_y = reshape(temp(ind), sz(1:3)+3);
   
   temp = O_grid(:,:,:,3);
   O_trans_z = reshape(temp(ind), sz(1:3)+3);
   
   O_trans = cat(4, O_trans_x, O_trans_y, O_trans_z);   
   [O_temp, Spacing] = refine_grid(O_trans, Spacing, sizeI);
   
   O_new = O_temp(2:end-2, 2:end-2, 2:end-2, :);
   
end

end




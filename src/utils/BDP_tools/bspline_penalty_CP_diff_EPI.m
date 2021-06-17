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


function [O_penalty, O_grad] = bspline_penalty_CP_diff_EPI(O_grid1, spacing, res, penalty_wt)
% This function computes penalty based on difference of adjacent knot points. Similar to the 
% regularization proposed in (Chun and Fessler, IEEE STSP, 2009) but this function computes penalty
% for deformation which is only constrained in one dimensions (ie. Ogrid1 should be a 3D volume of
% x-coordinates of the control points). 
%
%  Ogrid1 : A 3D volume of x-coordinates of the b-spline control-point grid
%  res : 1x3 vector representing dimensions of voxels
%  Spacing : 1x3 vector of space between control points (in voxels)
%  penalty_wt : A vector in form [kx bigKx ky] or [kx bigKx]. Later form is equivalent to ky = kx.
%               Implicitly it is assumed that kz = ky in both cases. 
%               To encourage invertibility, all of following must be satisfied: 
%               (1) (kx + ky + kz)<1 
%               (2) kx, ky, kz must be in range [0, 0.5)
%               (3) bigKx >= -kx
%

if numel(penalty_wt)~=2 && numel(penalty_wt)~=3
   error('penalty_wt must be either be a vector of length 2 or 3.');
end

kx = penalty_wt(1);
bigkx = penalty_wt(2);
if numel(penalty_wt)~=3
   ky = penalty_wt(1);
else
   ky = penalty_wt(3);
end

if kx<0 || kx>=0.5 
   error('kx has to be in range [0, 0.5)');
elseif ky<0 || ky>=0.5 
   error('ky has to be in range [0, 0.5)');
elseif bigkx<-kx
   error('bigkx must be >= -kx')
end

if max(penalty_wt)>0
   spacing_mm = spacing .* res;
   c_knots = O_grid1 * res(1);
   clear O_grid1
   
   c_diff_x = diff(c_knots, 1, 1) - spacing_mm(1);
   c_diff_y = diff(c_knots, 1, 2);
   c_diff_z = diff(c_knots, 1, 3);
   clear c_knots
   
   if nargout==2
      [px, pxg] = calculate_penalty(c_diff_x, -spacing_mm(1)*kx, spacing_mm(1)*bigkx, 1);
      [py, pyg] = calculate_penalty(c_diff_y, -spacing_mm(1)*ky, spacing_mm(1)*ky, 2);
      [pz, pzg] = calculate_penalty(c_diff_z, -spacing_mm(1)*ky, spacing_mm(1)*ky, 3);
      O_penalty = mean(px(:)) + mean(py(:)) + mean(pz(:));
      O_grad = pxg./numel(px) + pyg./numel(py) + pzg./numel(pz);
   else
      px = calculate_penalty(c_diff_x, -spacing_mm(1)*kx, spacing_mm(1)*bigkx, 1);
      py = calculate_penalty(c_diff_y, -spacing_mm(1)*ky, spacing_mm(1)*ky, 2);
      pz = calculate_penalty(c_diff_z, -spacing_mm(1)*ky, spacing_mm(1)*ky, 3);
      O_penalty = mean(px(:)) + mean(py(:)) + mean(pz(:));
   end
   
else
   O_penalty = eps;
   O_grad = zeros(size(O_grid1));
end

end

function [PE, Pgrad] = calculate_penalty(c_diff, el, eu, dim)
m1 = c_diff<el;
m3 = c_diff>eu;

c_diff(m1) = c_diff(m1)-el;
c_diff(m3) = c_diff(m3)-eu;
c_diff(~m1 & ~m3) = 0;

PE = 0.5*(c_diff.^2);

if nargout==2
   padsz = [0 0 0];
   padsz(dim) = 1;
   Pgrad = -1 * diff(padarray(c_diff, padsz, 0, 'both'), 1, dim); % -1 b/c weights are opposite of what diff() computes
end

end


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


function [V, cp_ind, Q] = bspline_penalty_repeated_endpoint_metal_bending_operator(size_O, spacing, res)
% Computes operators for computing bending energy directly from control points for 3D volume.
% The operators applies reptition of control points on edges - ONLY works when control-points
% are parameterizing deformation (& not map). 
%
% This function uses the mathematical expressions given in 
%  Shackleford et al, "Analytic regularization of uniform cubic B-spline deformation fields",
%  MICCAI 2012, pp 122-129
%
% Inputs: 
%  size_O - size of bspline control points 
%  spacing - Spacing between control points
%  res - resolution along each direction 
%
% Bending energy can be computed from operators as following:
%
%   px = O(:,:,:,1); % get x-dim (or y,z dim) of bspline grid
%   p = px(cp_ind);
%   pV = p*V;
%   pV = pV(:);
%   S_penalty = transpose(pV)*p(:); % bending energy
%   G = reshape(Q*pV, size_O); % gradient wrt CP
%   clear pV p px
%
% See bspline_penalty_repeated_endpoint_metal_bending_EPI()
%

size_O = transpose(size_O(:));
spacing = transpose(spacing(:));
res = transpose(res(:));
spacing_mm = spacing.*res;

V_opr = sixMatOperators(spacing_mm);
V = sum(V_opr(:,:,1:3), 3) +  2*sum(V_opr(:,:,4:6), 3); % sum over all 6 different orders
clear V_opr

% Generate index of 64 control-points for each segment - also repeats end-points once

local_ind = -1:2;
[seg_x, seg_y, seg_z] = ndgrid(1:(size_O(1)-1), 1:(size_O(2)-1), 1:(size_O(3)-1));
seg_x = seg_x(:);
seg_y = seg_y(:);
seg_z = seg_z(:);

[lx, ly, lz] = ndgrid(local_ind);
lx = lx(:)';
ly = ly(:)';
lz = lz(:)';

n_CPneighbor = numel(lx); % should be 64
n_seg = numel(seg_x); % should be same as prod(size_O-1)

cp_x = seg_x(:, ones(1,n_CPneighbor)) + lx(ones(1,n_seg), :);
cp_y = seg_y(:, ones(1,n_CPneighbor)) + ly(ones(1,n_seg), :);
cp_z = seg_z(:, ones(1,n_CPneighbor)) + lz(ones(1,n_seg), :);
clear seg_x seg_y seg_z

% repeat end-points 
cp_x(cp_x<1) = 1;
cp_x(cp_x>size_O(1)) = size_O(1);
cp_y(cp_y<1) = 1;
cp_y(cp_y>size_O(2)) = size_O(2);
cp_z(cp_z<1) = 1;
cp_z(cp_z>size_O(3)) = size_O(3);

cp_ind = sub2ind(size_O, cp_x, cp_y, cp_z);
clear cp_x cp_y cp_z lx ly lz


% Operator for computing gradient wrt CP 
% - Sums up terms in p'*V corresponding to each CP 
% - Elements in p'*V is same as V*p as V is symmetric 
i = cp_ind(:);
j = 1:numel(cp_ind);
s = 2*ones(numel(cp_ind),1);
Q = sparse(i, j, s, prod(size_O), numel(cp_ind));
clear i j s 

end

function V = sixMatOperators(spacing)
% computes six operators corresponding to different partial derivative terms.
% spacing - vector of length 6

ord = [... % order of partial derivatives, column wise [x y z]
   2 0 0;...
   0 2 0;...
   0 0 2;...
   1 1 0;...
   1 0 1;...
   0 1 1;...
   ];

V = zeros(64, 64, 6);
for o = 1:size(ord, 1)
   SIG_x = bigSigmaMat(spacing(1), ord(o,1));
   SIG_y = bigSigmaMat(spacing(2), ord(o,2));
   SIG_z = bigSigmaMat(spacing(3), ord(o,3));
   V(:,:,o) = kron(SIG_x, kron(SIG_y, SIG_z));
end

end

function SIG = bigSigmaMat(sp, ord)
% Computes $\bar{\Sigma}$ for a particular order and spacing (along one dim)
% sp = scalar value representing spacing
% ord = derivative order, 0, 1 or 2

Q = QMat(sp, ord);
psi = [sp, (sp^2)/2, (sp^3)/3, (sp^4)/4, (sp^5)/5, (sp^6)/6, (sp^7)/7];

SIG = zeros(4,4);
for a = 1:4
   for b = 1:4
      SIG(a,b) = integralRegion(a, b);
   end
end

   function out = integralRegion(a, b)
      % dot product of small sigma and psi
      Xi = kron(vect(Q(a,:)),Q(b,:));
      sig = [...
         Xi(1); ...
         Xi(2)+Xi(5); ...
         Xi(3)+Xi(6)+Xi(9);...
         Xi(4)+Xi(7)+Xi(10)+Xi(13);...
         Xi(8)+Xi(11)+Xi(14);...
         Xi(12)+Xi(15);...
         Xi(16);...
         ];
      out = psi*sig;
   end
end

function Q = QMat(sp, ord)
% sp = scalar value representing spacing
% ord = derivative order, 0, 1 or 2

B = [1 -3 3 -1; 4 0 -6 3; 1 3 3 -3; 0 0 0 1]./6;
sp = 1/sp;
R = diag([1 sp sp^2 sp^3]);
Delta = cat(3, eye(4), [0 0 0 0; 1 0 0 0; 0 2 0 0; 0 0 3 0], ...
   [0 0 0 0; 0 0 0 0; 2 0 0 0; 0 6 0 0]);

Q = B*R*Delta(:,:,ord+1);
end

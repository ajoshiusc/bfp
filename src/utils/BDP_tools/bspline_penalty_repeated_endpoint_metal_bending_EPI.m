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


function [O_penalty, O_grad] = bspline_penalty_repeated_endpoint_metal_bending_EPI(O, V, cp_ind, Q)
% Computes bending energy of thin metal sheet from bspline control-points. It works only when
% control-points are paramterizing deformation (& NOT map). Check these functions: 
%    bspline_repeated_endpoint_deformation_3d_double_only_x()
%    bspline_grid_generate_deformation()
%
% See bspline_penalty_repeated_endpoint_metal_bending_operator() for details of V, cp_ind, Q
%

% use only x-dim for EPI
O = O(:,:,:,1); 

p = O(cp_ind);
pV = p*V;
pV = pV(:);
O_penalty = transpose(pV)*p(:); % bending energy

if nargout>1
   O_grad = reshape(Q*pV, size(O)); % gradient wrt CP
end

end

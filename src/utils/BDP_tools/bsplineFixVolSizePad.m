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


function [vol_out, st_vox, end_vox, sz_out] = bsplineFixVolSizePad(vol, spacing)
% fix image size by padding zeros - make image size to be a perfect fit for control points

input_size = size(vol);
vox_diff = (ceil(input_size./spacing).*spacing)+ 1 - input_size;
st_vox = floor(vox_diff/2)+1;
end_vox = st_vox + input_size - 1;

sz_out = input_size + vox_diff;
vol_out = zeros(sz_out);
vol_out(st_vox(1):end_vox(1), st_vox(2):end_vox(2), st_vox(3):end_vox(3)) = vol;

end

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


function overlay_volume_edge2png(vol, vol_overlay, clim, pngFileName, edge_param)

slices3 = 1:size(vol, 3);
slices1 = 1:size(vol, 1);
slices2 = 1:size(vol, 2);
nRows = NaN;
nCols = NaN;

if exist('edge_param', 'var')
   sigma = edge_param(1);
   threshold = edge_param(2);
else
   sigma = 2;
   threshold = 0.002;
end

vol_edge =  zeros(size(vol_overlay));
for k = 1:size(vol_overlay,3)
   vol_edge(:,:,k) = edge(vol_overlay(:,:,k), 'log', threshold, sigma);
end
A = overlay_volumes(vol, vol_edge, clim, 'edge');

E = permute(A, [2 1 4 3]);
E = flipdim(E, 1);
I1 = montage_image(E , 'size', [nRows nCols], 'Indices', slices3);
imwrite(I1, [pngFileName '.1.png'], 'png');
clear I1


vol_edge =  zeros(size(vol_overlay));
for k = 1:size(vol_overlay,2)
   vol_edge(:,k,:) = permute(edge(squeeze(vol_overlay(:,k,:)), 'log', threshold, sigma), [1 3 2]);
end
A = overlay_volumes(vol, vol_edge, clim, 'edge');

E = permute(A, [3 1 4 2]);
E = flipdim(E, 1);
I2 = montage_image(E , 'size', [nRows nCols], 'Indices', slices2);
imwrite(I2, [pngFileName '.2.png'], 'png');
clear I2


vol_edge =  zeros(size(vol_overlay));
for k = 1:size(vol_overlay,1)
   vol_edge(k,:,:) = permute(edge(squeeze(vol_overlay(k,:,:)), 'log', threshold, sigma), [3 1 2]);
end
A = overlay_volumes(vol, vol_edge, clim, 'edge');

E = permute(A, [3 2 4 1]);
E = flipdim(E, 1);
E = flipdim(E, 2);
I3 = montage_image(E , 'size', [nRows nCols], 'Indices', slices1);
imwrite(I3, [pngFileName '.3.png'], 'png');
clear I3 A E


end


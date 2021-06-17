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


function overlay_volumes2png(vol, vol_overlay, clim, pngFileName, type, type_clim)
% Generates colorFA in png format of the slices. 
%   CLIM - Vector of length two - lower limit & upper limit of intensity
%   PNGFILENAME - file name with full path
%   TYPE, TYPE_CLIM - Optional. See overlay_volumes.m
%

if ~exist('type', 'var')
   type = 'redgreen';
elseif strcmpi(type, 'edge')
   if exist('type_clim', 'var')
      overlay_volume_edge2png(vol, vol_overlay, clim, pngFileName, type_clim);
   else
      overlay_volume_edge2png(vol, vol_overlay, clim, pngFileName);
   end
   return;
end

if ~exist('type_clim', 'var')
   type_clim = [70 94];
end

A = overlay_volumes(vol, vol_overlay, clim, type, type_clim);
vol2png(A, [0 1], pngFileName, 'rgb');

end


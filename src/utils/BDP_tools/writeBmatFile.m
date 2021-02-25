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


function writeBmatFile(bMatrices, fname)
% bMatrices should be 3x3xnDir

fid = fopen(fname, 'w');
for iDir = 1:size(bMatrices,3)
   temp = squeeze(bMatrices(:,:,iDir))';
   temp = temp(:)';
   fprintf(fid, '%22.15f %22.15f %22.15f\n%22.15f %22.15f %22.15f\n%22.15f %22.15f %22.15f\n\n', temp);
end
fclose(fid);
end

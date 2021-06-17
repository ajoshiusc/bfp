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


function bMatrices = readBmat(fname)
% bMatrices are 3x3xnDir

fid = fopen(fname, 'r');

if fid<0
   error('BDP:FileDoesNotExist', ['Could not open the bmat file: ' escape_filename(fname)]);
end

data = fscanf(fid, '%f');
fclose(fid);

if mod(numel(data), 9)~=0
   err_msg = ['Number of elements in the bmat file seems to be invalid: ' escape_filename(fname) ...
      '\n Please make sure that the bmat has only 3x3 matrices. Total number of matrices should be '...
      'same as number of images in diffusion dataset.'];
   error('BDP:InvalidFile', bdp_linewrap(err_msg));
end

data = reshape(data, 3, 3, []);
bMatrices = permute(data, [2 1 3]);

end

% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2017 The Regents of the University of California and the University of Southern California
% Created by Anand A. Joshi, Chitresh Bhushan, David W. Shattuck, Richard M. Leahy 
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


function [fileName, ext] = remove_extension( fileName )
% Returns file name (with path) removing following extensions:
%

ext = '';

if length(fileName)>=3 && strcmpi(fileName(end-2:end), '.gz')
   ext = [fileName(end-2:end) ext];
   fileName = fileName(1:end-3);   
end

while length(fileName)>=4 && ( ...
      strcmpi(fileName(end-3:end), '.nii') || strcmpi(fileName(end-3:end), '.img') || ...
      strcmpi(fileName(end-3:end), '.hdr') || strcmpi(fileName(end-3:end), '.dfs') ||...
      strcmpi(fileName(end-3:end), '.dfc') || strcmpi(fileName(end-3:end), '.txt') ||...
      strcmpi(fileName(end-3:end), '.eig') || strcmpi(fileName(end-3:end), '.ext') ||...
      strcmpi(fileName(end-3:end), '.gii'));
   
   ext = [fileName(end-3:end) ext];
   fileName = fileName(1:end-4);
end

if length(fileName)>=5 && strcmpi(fileName(end-4:end), '.bmat')
   ext = [fileName(end-4:end) ext];
   fileName = fileName(1:end-5);
end

end

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


function str_out = escape_filename(filename)
% Replaces all occurance of \ by \\ in input string or cellstr.
% It is useful for escaping windows filenames and paths so that fprintf does not generate:
%      Warning: Escape sequence 'U' is not valid
%

if ischar(filename)
   str_out = strrep(filename, '\\', '\'); % to avoid repeating a previously escaped \
   str_out = strrep(str_out, '\', '\\');
   
elseif iscellstr(filename)
   cfuna = @(s) strrep(s, '\\', '\'); % to avoid repeating a previously escaped \
   str_out = cellfun(cfuna, filename, 'UniformOutput', false);
   cfun = @(s) strrep(s, '\', '\\');
   str_out = cellfun(cfun, str_out, 'UniformOutput', false);
else
   error('Input type must be string or cellstring.')
end

end

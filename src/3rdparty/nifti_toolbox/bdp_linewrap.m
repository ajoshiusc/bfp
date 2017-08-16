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


function strout = bdp_linewrap(str, maxchars, separator)
% Wrapper around linewrap and strjoin to simplify linewrap

if nargin == 2 && ischar(maxchars)
   separator = maxchars;
   clear maxchars
end

if ~exist('maxchars', 'var')
   maxchars = 80;
end

if ~exist('separator', 'var')
   separator = '\n';
end

if ischar(str)
   strout = strjoin(linewrap(str, maxchars), separator);
   
elseif iscellstr(str)
   strout = cellfun(@(x) bdp_linewrap(x, maxchars, separator), str, 'UniformOutput', false);
   strout = strjoin(strout, '');
   
else
   error('Input str must be string or cellstr.')
end

end

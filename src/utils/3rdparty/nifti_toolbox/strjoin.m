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


function output = strjoin(input, separator)
%STRJOIN Concatenate an array into a single string.

if nargin < 2, separator = ','; end
assert(ischar(separator), 'Invalid separator input: %s', class(separator));
separator = strrep(separator, '%', '%%');

output = '';
if ~isempty(input)
    if ischar(input)
        input = cellstr(input);
    end
    if isnumeric(input) || islogical(input)
        output = [repmat(sprintf(['%.15g', separator], input(1:end-1)), ...
            1, ~isscalar(input)), ...
            sprintf('%.15g', input(end))];
    elseif iscellstr(input)
        output = [repmat(sprintf(['%s', separator], input{1:end-1}), ...
            1, ~isscalar(input)), ...
            sprintf('%s', input{end})];
    elseif iscell(input)
        output = strjoin(cellfun(@(x)strjoin(x, separator), input, ...
            'UniformOutput', false), ...
            separator);
    else
        error('strjoin:invalidInput', 'Unsupported input: %s', class(input));
    end
end
end

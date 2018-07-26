% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2016 The Regents of the University of California and the University of Southern California
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


function disp1(ss,modn,flags)

if ~exist('flags','var')
    flags='';
end

if isempty(strfind(flags,'v'))
    verbosity=2;
else
    a=strfind(flags,'v');
    verbosity=flags(a(1)+1);
end
verbosity=str2double(verbosity);
if verbosity == 0
    return;
end

if ~exist('modn','var')
    modn='';
end

if isempty(strfind(flags,'m'))
    modn='';
end
v=fix(clock);
if isempty(strfind(flags,'m'))
    if isempty(strfind(flags,'t'))
        fprintf('%s\n',ss);
    else
        fprintf('%02d:%02d:%02d %s\n',v(4),v(5),v(6),ss);
    end
else
    if isempty(strfind(flags,'t'))
        fprintf('%s: %s\n',modn,ss);
    else
        fprintf('SVREG:%s: %02d:%02d:%02d %s\n',modn,v(4),v(5),v(6),ss);
    end
end

%v=load([subbasename,'.verbosity.txt']);

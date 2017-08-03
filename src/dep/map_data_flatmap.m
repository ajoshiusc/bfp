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


function tar_data=map_data_flatmap(sub,subdata,tar,method,cleanflag)

if ~exist('method','var');
    method='linear';
end
if ~exist('cleanflag','var');
    cleanflag=1;
end
if cleanflag
    sub.attributes=subdata;
    sub=clean_sqr_map(sub);
    subdata=sub.attributes;
end
tar_data=zeros(length(tar.u),size(subdata,2));
for jj=1:size(subdata,2)
    tar_data(:,jj)=mygriddata(sub.u',sub.v',subdata(:,jj),tar.u',tar.v',method);
end



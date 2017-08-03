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



function znew=mygriddata(x,y,z,xin,yin,interpm)
    warning off
    
    if ~exist('interpm','var')
        F = TriScatteredInterp(x,y,z);
    else
        F = TriScatteredInterp(x,y,z,interpm);
    end
    
    znew=xin;
    znew(:)=F(xin(:),yin(:));
    indx=find(isnan(znew));
    %znew(indx)=-100;
    %disp('mygriddata changed!!');
    F = TriScatteredInterp(x,y,z,'nearest');

    znew(indx)=F(xin(indx),yin(indx));
    
    warning on

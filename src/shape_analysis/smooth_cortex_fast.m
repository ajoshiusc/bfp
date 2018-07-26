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



function [surf_sm,A,VC]=smooth_cortex_fast(surf,a,nIterations,A,size_preserve)
% function [surf_sm,A]=smooth_cortex_fast(surf,a,nIterations,A)
%
% Smooths a surface
%
% Input:
%   surf: surface to smooth
%   a: scalar smooth weighting parameter (0-1 less-more smoothing)
%   nIterations: number of times to apply the smoothing

%calculate smoothing matrix

if ~exist('A','var') || isempty(A)
    [vertConn,VC] = vertices_connectivity_fast(surf);
    A=spones(VC);
    A=spdiags((a./sum(A))',0,size(A,1),size(A,2))*A;
    A=A+spdiags((1-a)*ones(size(A,1),1),0,size(A,1),size(A,2));
end

%smooth
surf_sm = surf;
for i = 1:nIterations
    surf_sm.vertices=A*surf_sm.vertices;
end


if exist('size_preserve','var')

maxxsm=max(surf_sm.vertices(:,1));minxsm=min(surf_sm.vertices(:,1));
maxysm=max(surf_sm.vertices(:,2));minysm=min(surf_sm.vertices(:,2));
maxzsm=max(surf_sm.vertices(:,3));minzsm=min(surf_sm.vertices(:,3));
maxx=max(surf.vertices(:,1));minx=min(surf.vertices(:,1));
maxy=max(surf.vertices(:,2));miny=min(surf.vertices(:,2));
maxz=max(surf.vertices(:,3));minz=min(surf.vertices(:,3));

surf_sm.vertices(:,1)=surf_sm.vertices(:,1)-minxsm;
surf_sm.vertices(:,1)=surf_sm.vertices(:,1)*((maxx-minx)/(maxxsm-minxsm))+minx;

surf_sm.vertices(:,2)=surf_sm.vertices(:,2)-minysm;
surf_sm.vertices(:,2)=surf_sm.vertices(:,2)*((maxy-miny)/(maxysm-minysm))+miny;

surf_sm.vertices(:,3)=surf_sm.vertices(:,3)-minzsm;
surf_sm.vertices(:,3)=surf_sm.vertices(:,3)*((maxz-minz)/(maxzsm-minzsm))+minz;
% 
%     
%     A=tri_area(surf_sm.faces,surf_sm.vertices);     A=sum(A);
%     Ao=tri_area(surf.faces,surf.vertices);          Ao=sum(Ao);
%     scalf=sqrt(Ao/A);
%     mean1=mean(surf_sm.vertices(:,1));mean2=mean(surf_sm.vertices(:,2));mean3=mean(surf_sm.vertices(:,3));
%     surf_sm.vertices(:,1)=scalf*(surf_sm.vertices(:,1)-mean1)+mean1;
%     surf_sm.vertices(:,2)=scalf*(surf_sm.vertices(:,2)-mean2)+mean2;
%     surf_sm.vertices(:,3)=scalf*(surf_sm.vertices(:,3)-mean3)+mean3;        
end


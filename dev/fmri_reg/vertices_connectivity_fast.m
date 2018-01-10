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



function [VertConn,C] = vertices_connectivity_fast(FV,VERBOSE);
%
% Authors: Dimitrios Pantazis, Anand Joshi, November 2007 

rowno=[FV.faces(:,1);FV.faces(:,1);FV.faces(:,2);FV.faces(:,2);FV.faces(:,3);FV.faces(:,3)];
colno=[FV.faces(:,2);FV.faces(:,3);FV.faces(:,1);FV.faces(:,3);FV.faces(:,1);FV.faces(:,2)];
data=ones(size(rowno));
C=sparse(rowno,colno,data);
C=spones(C);%(C>0)=1;
[rows,cols,vals] = find(C);
d = find(diff([cols' 0]));

%if there are vertices with no connections at all
nVertices = size(FV.vertices,1);
if length(d)~=nVertices
    VertConn = vertices_connectivity(FV);
    return
end

%fast calculation of vertices connectivity
VertConn = cell(nVertices,1); % one empty cell per vertex

if isempty(d)
    VertConn={};
    return
end
VertConn{1} = rows(1:d(1))';
for i = 2:nVertices
    VertConn{i} = rows((d(i-1)+1):d(i))';
end



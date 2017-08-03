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



function [TriConn] = faces_connectivity_fast(FV,VERBOSE);
% function [TriConn] = faces_connectivity_fast(FV,VERBOSE);
%
% Computes faces connectivity fast
% 
% Input:
%   FV: the standard matlab structure for faces and vertices,
%         where FV.faces is m x 3, and FV.vertices is n x 3.
%   VERBOSE (optional): redundant, used for compatibility with an earlier function
% Output:
%   TriConn:  vector of cells, one-to-one with each row of FV.vertices.
%       TriConn{i} returns a row vector of the face numbers (rows in FV.faces) that
%       are connected to the ith vertex row in FV.vertices.
%       Thus if we want the faces numbers for a given set of vertices, TriConn{i} will tell
%       us which faces to use in for instance a patch command.
%
% See also VERTICES_CONNECTIVITY
%
% Author: Dimitrios Pantazis, November 2007 

%initialize
nFaces = size(FV.faces,1);
nVertices = size(FV.vertices,1);
TriConn = cell(nVertices,1); % one empty cell per vertex

%get assignment of vertices to faces
[Y,I] = sort([FV.faces(:,1);FV.faces(:,2);FV.faces(:,3)]);
I = mod(I,nFaces);
I(I==0)=nFaces;
%optional: sort vertices
[dummy,ndx] = sort([Y*(nFaces+1)+I]);
I = I(ndx);

%find when vertex changes
d = find(diff([Y' Y(end)+1]));

%if there are vertices without any connections, do a slow loop
if length(d)~=nVertices
    for i = 1:length(Y)
        TriConn{Y(i)} = [TriConn{Y(i)} I(i)];
    end
    return
end

%fast calculation of faces connectivity (if no isolated vertices)
TriConn{1} = I(1:d(1))';
for i = 2:nVertices
    TriConn{i} = I((d(i-1)+1):d(i))';
end



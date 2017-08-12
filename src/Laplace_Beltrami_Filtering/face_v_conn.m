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


function [TriConn] = face_v_conn(FV,VERBOSE)
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

rows=[FV.faces(:,1);FV.faces(:,2);FV.faces(:,3)];
cols=[(1:nFaces)';(1:nFaces)';(1:nFaces)'];
data=ones(length(rows),1);

TriConn=sparse(rows,cols,data);

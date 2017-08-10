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



function TRIn = myclean_tri(TRI1)
A = sort(TRI1,2);
[B,I,J] = unique(A,'rows');
TRI = TRI1(I,:);
%disp(sprintf('%d duplicates!!',size(TRI1,1)-size(TRI,1)));
S=(abs(TRI(:,1)-TRI(:,2))>0)+(abs(TRI(:,2)-TRI(:,3))>0)+(abs(TRI(:,1)-TRI(:,3))>0);
I=find(S==3);
%disp(sprintf('%d duplicate vertices!!',size(TRI,1)-size(I,1)));
TRIn=TRI(I,:);

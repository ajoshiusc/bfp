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


function sub=clean_sqr_map(sub)
% This function cleans the square map. If there are extreme angles 
% or long triangles, they are deleted by this function and patch is cleaned
% it is useful in conjunction with mygriddata for better interpolation

EdTh=10;
AngTh1=5;AngTh2=170;
u1=sub.u(sub.faces(:,1));v1=sub.v(sub.faces(:,1));
u2=sub.u(sub.faces(:,2));v2=sub.v(sub.faces(:,2));
u3=sub.u(sub.faces(:,3));v3=sub.v(sub.faces(:,3));

ed1=sqrt((u1-u2).^2+(v1-v2).^2)';
ed2=sqrt((u1-u3).^2+(v1-v3).^2)';
ed3=sqrt((u2-u3).^2+(v3-v2).^2)';
ang1=angle([u2'-u1',v2'-v1'],[u3'-u1',v3'-v1']);
ang2=angle([u3'-u2',v3'-v2'],[u1'-u2',v1'-v2']);
ang3=angle([u1'-u3',v1'-v3'],[u2'-u3',v2'-v3']);
ind=(ed1./ed2>EdTh | ed2./ed1>EdTh | ed3./ed2>EdTh | ed2./ed3>EdTh | ed1./ed3>EdTh | ed3./ed1>EdTh | ang1<AngTh1 | ang2<AngTh1 | ang3<AngTh1 | ang1>AngTh2 | ang2>AngTh2 | ang3>AngTh2 );
sub.faces(ind,:)=[];
sub=myclean_patch_cc(sub);

function ang=angle(V1,V2)
ang=acos(dot(V1',V2')'./(mynorm2(V2).*mynorm2(V1)));
ang=ang*180/pi;


function nrm=mynorm2(W)
nrm=sqrt(W(:,1).^2 + W(:,2).^2);


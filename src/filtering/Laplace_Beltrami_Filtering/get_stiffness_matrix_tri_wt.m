
function [S,Dx,Dy]=get_stiffness_matrix_tri_wt(surf1,W)
%This finds stiffness matrix for triangular elements using techniques from
%book by Sadiku: Numerical Techniques Formulas 6.17,6.18

% A=tri_area(surf1.faces,surf1.vertices);
% 
% 
% f=surf1.faces;v=surf1.vertices;
% 
% C11e=(1./(4.*A)).*(v(f(:,2),2)- v(f(:,3),2)).^2
% 

X=surf1.vertices(:,1);
Y=surf1.vertices(:,2);
Z=surf1.vertices(:,3);
NumTri=size(surf1.faces,1);
NumVertx=size(X,1);
vertx_1=surf1.faces(:,1);vertx_2=surf1.faces(:,2);vertx_3=surf1.faces(:,3);
V1=[X(vertx_1),Y(vertx_1),Z(vertx_1)];
V2=[X(vertx_2),Y(vertx_2),Z(vertx_2)];
V3=[X(vertx_3),Y(vertx_3),Z(vertx_3)];

x1=zeros(NumTri,1); y1=zeros(NumTri,1);
v2_v1temp=V2-V1;
x2=sqrt(sum(v2_v1temp.^2,2));
y2=zeros(NumTri,1);
x3=dot((V3-V1),(v2_v1temp./[x2,x2,x2]),2);
mynorm = cross((v2_v1temp),V3-V1,2);
yunit = cross(mynorm,v2_v1temp,2);
y3 = dot(yunit,(V3-V1),2)./sqrt(sum(yunit.^2,2));
sqrt_DT= (abs((x1.*y2 - y1.*x2)+(x2.*y3 - y2.*x3)+(x3.*y1 - y3.*x1)));
Ar=.5*(abs((x1.*y2 - y1.*x2)+(x2.*y3 - y2.*x3)+(x3.*y1 - y3.*x1)));%Ar=ones(size(Ar));
y1 = y1./sqrt_DT; y2 = y2./sqrt_DT; y3 = y3./sqrt_DT;
x1 = x1./sqrt_DT; x2 = x2./sqrt_DT; x3 = x3./sqrt_DT;
tmp_A=[y2-y3;y3-y1;y1-y2]; tmp_B=[x3-x2;x1-x3;x2-x1];

rowno=[1:NumTri]';%rowno_1=rowno-1;
rowno_all=[rowno;rowno;rowno];
vertx_all=[vertx_1;vertx_2;vertx_3];
Dx=sparse(rowno_all,vertx_all,tmp_A,NumTri,size(X,1));
Dy=sparse(rowno_all,vertx_all,tmp_B,NumTri,size(X,1));


[TC] = face_v_conn(surf1);

Wt=(1/3)*(TC'*W);Wt=spdiags(Wt.*Ar,0,NumTri,NumTri);
S=Dx'*Wt*Dx +Dy'*Wt*Dy;


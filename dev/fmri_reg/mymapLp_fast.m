function [x_map,y_map,TRI]=mymapLp_fast(TRI, X, Y, Z, NumTri, pinnedVertx,xbdr,ybdr,Pp)

NumVertx=size(X,1);
NumPinnedVertx=size(pinnedVertx,1);   
NumfreeVertx=size(X,1)-NumPinnedVertx;
    vertx_1=TRI(:,1);vertx_2=TRI(:,2);vertx_3=TRI(:,3);
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
    sqrt_DT= (abs((x1.*y2 - y1.*x2)+(x2.*y3 - y2.*x3)+(x3.*y1 - y3.*x1))).^((Pp-1)/Pp);    
    y1 = y1./sqrt_DT; y2 = y2./sqrt_DT; y3 = y3./sqrt_DT;
    x1 = x1./sqrt_DT; x2 = x2./sqrt_DT; x3 = x3./sqrt_DT;    
    tmp_A=[y2-y3,y3-y1,y1-y2]; tmp_B=[x3-x2,x1-x3,x2-x1];
    
rowno=2*[1:NumTri]';rowno_1=rowno-1;
rowno_all=[rowno_1;rowno_1;rowno_1;rowno;rowno;rowno];
vertx_all=[vertx_1;vertx_2;vertx_3;vertx_1;vertx_2;vertx_3];
tmp_AB_all=[tmp_A(:,1);tmp_A(:,2);tmp_A(:,3);tmp_B(:,1);tmp_B(:,2);tmp_B(:,3)];
A=sparse(rowno_all,vertx_all,tmp_AB_all);
indxFV=[1:NumVertx]';indxFV(pinnedVertx)=[];
A_f=A(:,indxFV); b1=sparse(-A(:,pinnedVertx));
pinnedX = xbdr; pinnedY = ybdr;
b1_x=b1*pinnedX;b1_y = b1*pinnedY;
AtA=A_f'*A_f; M=diag(AtA);
clear *all X Y Z i A
x_map=mypcg(AtA,A_f'*b1_x,1e-10,50000,M);
y_map=mypcg(AtA,A_f'*b1_y,1e-10,50000,M);

if (Pp ~= 2)
    x_map=mypcgdlp(A_f,b1_x,x_map,1e-12,10000,Pp);%<----
    y_map=mypcgdlp(A_f,b1_y,y_map,1e-12,10000,Pp);%<----
end
x_map([indxFV;pinnedVertx])=[x_map;xbdr];
y_map([indxFV;pinnedVertx])=[y_map;ybdr];
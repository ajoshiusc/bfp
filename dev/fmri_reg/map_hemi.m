%||AUM||
%||Shree Ganeshaya Namaha||

function [xmap,ymap]=map_hemi(surf)

MINV_IHF=25;%Min No of vertices in interhemispheric fissure

surfvertConn=vertices_connectivity_fast(surf,1);

jj=1;bdr=[];
while(length(bdr)<MINV_IHF)
while (length(surfvertConn{jj})>3)
    jj=jj+1;
end
apst=jj;
bdr=trace_boundary(apst,surfvertConn,surf);
jj=jj+1;
end
view_patch(surf);



[temp,p1]=max(surf.vertices(bdr,2));
[temp,p2]=min(surf.vertices(bdr,2));
par=para_curve_circle_segmt(surf.vertices(bdr,:),[p1;p2],[0;pi]);
zb=surf.vertices(bdr,3);
if mean(zb(find(par<pi)))>mean(zb(find(par>pi)))
    par=2*pi-par;
end;

%par=2*pi-par;
Xb=surf.vertices(bdr,1);Yb=surf.vertices(bdr,2);Zb=surf.vertices(bdr,3);
view_patch(surf);fh=find(par<pi);sh=find(par>=pi);view(-90,0);
hold on;line(Xb(fh),Yb(fh),Zb(fh),'color','r');line(Xb(sh),Yb(sh),Zb(sh),'color','b');
%line(surf2.vertices(C2,1),surf2.vertices(C2,2),surf2.vertices(C2,3),'color','r');
title('Blue should be above and Red should be below!!')
xbdr=cos(par);ybdr=sin(par);plot(xbdr,ybdr);
%sqrmap
xbdr=sign(xbdr).*(xbdr.^2); ybdr=sign(ybdr).*(ybdr.^2);
tmpxybdr=[xbdr,ybdr]*[1,-1;1,1];
xbdr=tmpxybdr(:,1);ybdr=tmpxybdr(:,2);
[xmap,ymap]=mymapLp_fast(surf.faces,surf.vertices(:,1),surf.vertices(:,2),surf.vertices(:,3),size(surf.faces,1),bdr,xbdr,ybdr,2);


disp('mapping on unit square');



clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/3rdParty'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/MEX_Files'));
MINV_IHF=100;%Min No of vertices in interhemispheric fissure

load('/deneb_disk/studyforrest_bfp/sub-03/func/sub-03_ses-movie_task-movie_run-3_bold.32k.GOrd.mat')
load /home/ajoshi/coding_ground/bfp/supp_data/HCP_32k_Label.mat

sl=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kleft.dfs');
numVert=length(sl.vertices);


labs=brainstructure(1:numVert);
fmriL=dtseries(1:numVert,:);

view_patch(sl);


numnan=sum(isnan(labs(sl.faces)),2);
sl.faces(numnan==3,:)=[];
[sl,ind]=myclean_patch_cc(sl);

fmriL=fmriL(ind,:);
% Falt mapping of hemisphere
[xmap,ymap]=map_hemi(sl);

figure;
patch('faces',sl.faces,'vertices',[xmap,ymap],'facevertexcdata',sl.vertices(:,1),'edgecolor','k','facecolor','interp');
axis equal;axis off;camlight;material dull;

%%
[surf1vertConn,C1]=vertices_connectivity_fast(sl);
surf1facesVConn=faces_connectivity_fast(sl);
surf1facesConn=faces2faces_connectivity(sl,surf1facesVConn);
boundary1 = boundary_vertices(sl,surf1vertConn,surf1facesVConn,surf1facesConn);

view_patch(sl);
jj=1;
bdr1=[];

while(length(bdr1)<MINV_IHF) % There must be atleast MINV_IHF vertices in the interhemispherical fissure!!
   apst=boundary1(jj);
   bdr1=trace_boundary(apst,surf1vertConn,sl);
   jj=jj+1;
end
dim=2;
[~,p1]=max(sl.vertices(bdr1,dim));
[~,p2]=min(sl.vertices(bdr1,dim));
par=para_curve_circle_segmt(sl.vertices(bdr1,:),[p1;p2],[0;pi]);



xbdr=cos(par);ybdr=sin(par);


mysphere(sl.vertices(bdr1,:),1,'r',10);
mysphere(sl.vertices(bdr1(p1),:),2,'g',10);
mysphere(sl.vertices(bdr1(p2),:),2,'k',10);

xbdr=sign(xbdr).*(xbdr.^2); ybdr=sign(ybdr).*(ybdr.^2);

tmpxybdr=[xbdr,ybdr]*[1,-1;1,1];

xbdr1=tmpxybdr(:,1);ybdr1=tmpxybdr(:,2);
Verbose=1;
if(Verbose)
   figure; plot(xbdr1,ybdr1);title('This should be full circle!!');
end

%%

dtseries(isnan(labs),:)=300;
figure;
patch('faces',sl.faces,'vertices',sl.vertices,'facevertexcdata',dtseries(1:length(sl.vertices),6),'edgecolor','none','facecolor','interp');
axis equal;axis off;view(90,0),camlight;colormap hot;material dull;caxis([-300,300]);colormap;

sr=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kright.dfs');
funr=dtseries(1+length(sr.vertices):2*length(sr.vertices),6);
funr(isnan(labs))=300;
figure;
patch('faces',sr.faces,'vertices',sr.vertices,'facevertexcdata',funr,'edgecolor','none','facecolor','interp');
axis equal;axis off;view(90,0),camlight;colormap hot;material dull;caxis([-300,300]);colormap;




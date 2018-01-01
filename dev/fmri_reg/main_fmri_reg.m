clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/3rdParty'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/MEX_Files'));
MINV_IHF=100;%Min No of vertices in interhemispheric fissure
%add weights based on surface area in the cost function
load /home/ajoshi/coding_ground/bfp/supp_data/HCP_32k_Label.mat
NPTS=256;

sl=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kleft.dfs');
numVert=length(sl.vertices);


labs=brainstructure(1:numVert);

view_patch(sl);


numnan=sum(isnan(labs(sl.faces)),2);
sl.faces(numnan==3,:)=[];
[sl,ind]=myclean_patch_cc(sl);


% Flat mapping of hemisphere
[xmap,ymap]=map_hemi(sl);

figure;
patch('faces',sl.faces,'vertices',[xmap,ymap],'facevertexcdata',sl.vertices(:,1),'edgecolor','k','facecolor','interp');
axis equal;axis off;camlight;material dull;
%%
% Resample fmri data to square
load('/deneb_disk/studyforrest_bfp/sub-03/func/sub-03_ses-movie_task-movie_run-3_bold.32k.GOrd.mat');
% fMRI data
fmriL=dtseries(1:numVert,:);
fmriL=fmriL(ind,:);
ll=linspace(-1,1,NPTS);
[X,Y]=meshgrid(ll,ll);
fmriLSq=zeros(size(X,1),size(X,2),size(fmriL,2));
parfor jj=1:size(fmriL,2)
    fmriLSq(:,:,jj)=mygriddata(xmap,ymap,fmriL(:,jj),X,Y);
    jj
end
I1=fmriLSq;
%%

load('/deneb_disk/studyforrest_bfp/sub-02/func/sub-02_ses-movie_task-movie_run-3_bold.32k.GOrd.mat');
% fMRI data
fmriL=dtseries(1:numVert,:);
fmriL=fmriL(ind,:);
ll=linspace(-1,1,NPTS);
[X,Y]=meshgrid(ll,ll);
fmriLSq=zeros(size(X,1),size(X,2),size(fmriL,2));
parfor jj=1:size(fmriL,2)
    fmriLSq(:,:,jj)=mygriddata(xmap,ymap,fmriL(:,jj),X,Y);
    jj
end
I2=fmriLSq;
%% Compute first fundamental form
xx=mygriddata(xmap,ymap,sl.vertices(:,1),X,Y);
yy=mygriddata(xmap,ymap,sl.vertices(:,2),X,Y);
zz=mygriddata(xmap,ymap,sl.vertices(:,3),X,Y);

[xx_u,xx_v]=gradient(xx);
[yy_u,yy_v]=gradient(yy);
[zz_u,zz_v]=gradient(zz);

g11=xx_u.^2+yy_u.^2+zz_u.^2;
g22=xx_v.^2+yy_v.^2+zz_v.^2;
g12=xx_u.*xx_v+yy_u.*yy_v+zz_u.*zz_v;

g=g11.*g22-g12.^2;
figure;
imagesc(g);


%% Set static and moving image
S=I2; M=I1;
clear fmriLSq I2 dtseries fmriL g g11 g12 g22

% Alpha (noise) constant
alpha=2.5;

% Velocity field smoothing kernel
Hsmooth=fspecial('gaussian',[60 60],10);
%g=imfilter(g,Hsmooth);
%g(g<0)=0;g(isnan(g))=0;
%figure; imagesc(g);
%M=M.*sqrt(g);S=S.*sqrt(g);
% The transformation fields
Tx=zeros(size(M,1),size(M,2)); Ty=zeros(size(M,1),size(M,2));

[Sy,Sx] = gradient(S);
[X,Y]=meshgrid(1:NPTS);
NIT=10000;
costiter=zeros(NIT,1);
hh=tic;
for itt=1:NIT
    % Difference image between moving and static image
    Idiff=M-S;
    
    % Default demon force, (Thirion 1998)
    %Ux = -(Idiff.*Sx)./((Sx.^2+Sy.^2)+Idiff.^2);
    %Uy = -(Idiff.*Sy)./((Sx.^2+Sy.^2)+Idiff.^2);
    
    % Extended demon force. With forces from the gradients from both
    % moving as static image. (Cachier 1999, He Wang 2005)
    [My,Mx] = gradient(M);
    Ux = sum(-Idiff.*  ((Sx./(sum(Sx.^2+Sy.^2,3)+alpha^2*sum(Idiff.^2,3)))+(Mx./(sum(Mx.^2+My.^2,3)+alpha^2*sum(Idiff.^2,3)))),3);
    Uy = sum(-Idiff.*  ((Sy./(sum(Sx.^2+Sy.^2,3)+alpha^2*sum(Idiff.^2,3)))+(My./(sum(Mx.^2+My.^2,3)+alpha^2*sum(Idiff.^2,3)))),3);
    
    % When divided by zero
    Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;
    
    % Smooth the transformation field
    Uxs=imfilter(Ux,Hsmooth);
    Uys=imfilter(Uy,Hsmooth);
    
    % Add the new transformation field to the total transformation field.
    Tx=Tx+Uxs;
    Ty=Ty+Uys;
    %        M=movepixels(I1,Tx,Ty);
    for kk=1:size(M,3)
        M(:,:,kk)=interp2(I1(:,:,kk),max(min(X+Ty,size(X,1)),1),max(min(Y+Tx,size(Y,1)),1));
    end
    costiter(itt)=norm(Idiff(:));
    fprintf('iter = %d, diff=%g, def=%g\n',itt,costiter(itt),sqrt(mean((Tx(:)).^2+(Ty(:)).^2)));
end
t1=toc(hh)

%  [X,Y]=meshgrid(1:NPTS);%
YY1=min(max(1,Y+Tx),NPTS);XX1=min(max(1,X+Ty),NPTS);

WX1=((NPTS-1)/2)*(xmap+1)+1;
WY1=((NPTS-1)/2)*(ymap+1)+1;

WX1=min(max(WX1,1),NPTS);WY1=min(max(WY1,1),NPTS);

xmap2=interp2((XX1),WX1',WY1');
ymap2=interp2((YY1),WX1',WY1');

xmap2=(xmap2-1)*(2/(NPTS-1)) - 1;
ymap2=(ymap2-1)*(2/(NPTS-1)) - 1;

figure;
patch('faces',sl.faces,'vertices',[xmap,ymap],'facevertexcdata',sqrt((xmap2'-xmap).^2+(ymap2'-ymap).^2),'edgecolor','k','facecolor','interp');
axis equal;axis off;camlight;material dull;

figure;
patch('faces',sl.faces,'vertices',[xmap2',ymap2'],'facevertexcdata',sqrt((xmap2'-xmap).^2+(ymap2'-ymap).^2),'edgecolor','k','facecolor','interp');
axis equal;axis off;camlight;material dull;

figure;plot(costiter);
save temp1

%%
a=load('/home/ajoshi/coding_ground/bfp/supp_data/bci_grayordinates_surf_ind.mat');
atl=readdfs('/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.left.mid.cortex.dfs');
sl.labels=atl.labels(a.ind_left(ind));
writedfs('out1.dfs',sl);
recolor_by_label('out1.dfs','/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain');
sl=readdfs('out1.dfs');

slsm=smooth_cortex_fast(sl,.1,600);

figure;
patch('faces',sl.faces,'vertices',slsm.vertices,'facevertexcdata',sl.vcolor,'edgecolor','none','facecolor','interp');
axis equal;axis off;camlight;material dull;
slw=sl;
slw.labels=griddata(xmap,ymap,sl.labels,xmap2',ymap2','nearest');
writedfs('out1.dfs',slw);
recolor_by_label('out1.dfs','/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain');
slw=readdfs('out1.dfs');

figure;
patch('faces',slw.faces,'vertices',slsm.vertices,'facevertexcdata',slw.vcolor,'edgecolor','none','facecolor','interp');
axis equal;axis off;camlight;material dull;



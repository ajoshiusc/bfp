function [wsub,origmap,newmap,surfObj,costiter,ind] = fmri_demons(sub1,sub2,BFPPATH,hemi)

NPTS=256;
%% Alpha (noise) constant
alpha=2.5;
SMPARA=3;
%% Number of iterations
NIT=1000;
RDIM=5;
surfObj=readdfs(fullfile(BFPPATH,'supp_data',['bci32k',hemi,'.dfs']));
numVert=length(surfObj.vertices);
a=load(fullfile(BFPPATH,'supp_data',['HCP_32k_Label','.mat']));

labs=a.brainstructure(1:numVert);

numnan=sum(isnan(labs(surfObj.faces)),2);
surfObj.faces(numnan==3,:)=[];
[surfObj,ind]=myclean_patch_cc(surfObj);


%% Flat mapping of hemisphere
[xmap,ymap]=map_hemi(surfObj);


%a=load(sub1);
fMRI1=sub1.dtseries(1:numVert,:);

%a=load(sub2);
fMRI2=sub2.dtseries(1:numVert,:);

%% load inf for hemi
%fMRI1=fMRI1(ind,:);fMRI1(isnan(fMRI1))=0;
%fMRI2=fMRI2(ind,:);fMRI2(isnan(fMRI2))=0;
fMRI1 = normalizeData(fMRI1(ind,:)')';
fMRI2 = normalizeData(fMRI2(ind,:)')';
fMRI2 = brainSync(fMRI1',fMRI2')';
[W,Y]=pca([fMRI1;fMRI2]);

fMRI1=fMRI1/W';
fMRI2=fMRI2/W';
fMRI1=fMRI1(:,[1:RDIM]);
fMRI2=fMRI2(:,[1:RDIM]);

%% data 1
fMRIData=fMRI1;
lL=linspace(-1,1,NPTS);
[X,Y]=meshgrid(lL,lL);
fMRIDataSq=zeros(size(X,1),size(X,2),size(fMRIData,2));
parfor jj=1:size(fMRIData,2)
    fMRIDataSq(:,:,jj)=mygriddata(xmap,ymap,fMRIData(:,jj),X,Y);
    jj
end
I1=fMRIDataSq;

%% data 2
fMRIData=fMRI2;
%% Sim
if (1)
    Z=zeros(256,256);
    Z(128,128)=5000;
    H=fspecial('gaussian',round([6*10 6*10]),10);
    Zs=imfilter(Z,H)*2/256;
else
    Zs=0;
end

%%
lL=linspace(-1,1,NPTS);
[X,Y]=meshgrid(lL,lL);
fMRIDataSq=zeros(size(X,1),size(X,2),size(fMRIData,2));
parfor jj=1:size(fMRIData,2)
    fMRIDataSq(:,:,jj)=mygriddata(xmap,ymap,fMRIData(:,jj),X+Zs,Y);
    jj
end
I2=fMRIDataSq;

%% Compute first fundamental form
xx=mygriddata(xmap,ymap,surfObj.vertices(:,1),X,Y);
yy=mygriddata(xmap,ymap,surfObj.vertices(:,2),X,Y);
zz=mygriddata(xmap,ymap,surfObj.vertices(:,3),X,Y);

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
clear fMRIDataSq I2 dtseries g g11 g12 g22

costiter=zeros(NIT,1);
res1=256;%[64,128,256];
Tx=0;Ty=0;
Mo=M;So=S;
for r1=1:length(res1)

    NPTS=res1(r1);
    [X,Y]=meshgrid(1:NPTS);
    
    if r1>1
        % upsample the grid and the deformation field
        Tx=2*Tx;Ty=2*Ty;
        Tx=interp2(Tx,max(1,0.5*X),max(1,0.5*Y));
        Ty=interp2(Ty,max(1,0.5*X),max(1,0.5*Y));
    end
    
    M=zeros(NPTS,NPTS,size(Mo,3));S=M;I1=M;
    parfor kk=1:size(M,3)
        M(:,:,kk)=interp2(Mo(:,:,kk),max(min(1+(256/NPTS)*(X+Tx-1),256),1),max(min(1+(256/NPTS)*(Y+Ty-1),256),1));
        S(:,:,kk)=interp2(So(:,:,kk),max(min(1+(256/NPTS)*(X+Tx-1),256),1),max(min(1+(256/NPTS)*(Y+Ty-1),256),1));
        I1(:,:,kk)=interp2(Mo(:,:,kk),max(min(1+(256/NPTS)*(X-1),256),1),max(min(1+(256/NPTS)*(Y-1),256),1));        
    end
    [Sx,Sy] = gradient(S);
    ks=SMPARA*(NPTS/256);
    Hsmooth=fspecial('gaussian',round([6*ks 6*ks]),ks);

    for itt=1:NIT
        % Difference image between moving and static image
        Idiff=M-S;
            
        % Extended demon force. With forces from the gradients from both
        % moving as static image. (Cachier 1999, He Wang 2005)
        [Mx,My] = gradient(M);
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
        
        %  Warp the fmri data
        parfor kk=1:size(M,3)
            M(:,:,kk)=interp2(I1(:,:,kk),max(min(X+Tx,size(X,1)),1),max(min(Y+Ty,size(Y,1)),1));
        end
        costiter(itt+NIT*(r1-1))=norm(Idiff(:))/NPTS;
        fprintf('res=%d iter = %d, diff=%g, def=%g\n',res1(r1),itt,costiter(itt+NIT*(r1-1)),sqrt(mean((Tx(:)).^2+(Ty(:)).^2)));
    end
end

%% Warp the flat map
YY1=min(max(1,Y+Ty),NPTS);XX1=min(max(1,X+Tx),NPTS);

WX1=((NPTS-1)/2)*(xmap+1)+1;
WY1=((NPTS-1)/2)*(ymap+1)+1;

WX1=min(max(WX1,1),NPTS);WY1=min(max(WY1,1),NPTS);

xmap2=interp2((XX1),WX1',WY1');
ymap2=interp2((YY1),WX1',WY1');

xmap2=(xmap2-1)*(2/(NPTS-1)) - 1;
ymap2=(ymap2-1)*(2/(NPTS-1)) - 1;
origmap=[xmap,ymap];
newmap=[xmap2',ymap2'];

%% Warp fmri data
wsub=zeros(size(fMRIData));
parfor jj=1:size(fMRIData,2)
    wsub(:,jj)=mygriddata(xmap,ymap,fMRIData(:,jj),xmap2',ymap2');
    jj
end

function wsub2 = fmri_demons(sub1,sub2,hemi)

'load inf for hemi
fmriL1 = normalizeData(fmriL1(ind,:)')';
fmriL2 = normalizeData(fmriL2(ind,:)')';
fmriL2 = brainSync(fmriL1',fmriL2')';


load_flatmap32kmesh


fmriL=fmriL1;
ll=linspace(-1,1,NPTS);
[X,Y]=meshgrid(ll,ll);
fmriLSq=zeros(size(X,1),size(X,2),size(fmriL,2));
parfor jj=1:size(fmriL,2)
    fmriLSq(:,:,jj)=mygriddata(xmap,ymap,fmriL(:,jj),X,Y);
    jj
end
I1=fmriLSq;
%%

% load('/deneb_disk/studyforrest_bfp/sub-02/func/sub-02_ses-movie_task-movie_run-3_bold.32k.GOrd.mat');
% % fMRI data
% fmriL=dtseries(1:numVert,:);
% fmriL=fmriL(ind,:);
fmriL=fmriL2;

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

NIT=200;
costiter=zeros(NIT,1);
hh=tic;
res1=[64,128,256];
Tx=0;Ty=0;
Mo=M;So=S;
for r1=1:3
    %X,Y
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
        M(:,:,kk)=interp2(Mo(:,:,kk),max(min(1+(256/NPTS)*(X+Ty-1),256),1),max(min(1+(256/NPTS)*(Y+Tx-1),256),1));
        S(:,:,kk)=interp2(So(:,:,kk),max(min(1+(256/NPTS)*(X+Ty-1),256),1),max(min(1+(256/NPTS)*(Y+Tx-1),256),1));
        I1(:,:,kk)=interp2(Mo(:,:,kk),max(min(1+(256/NPTS)*(X-1),256),1),max(min(1+(256/NPTS)*(Y-1),256),1));        
    end
    [Sy,Sx] = gradient(S);
    ks=round(10*(NPTS/256));
    Hsmooth=fspecial('gaussian',[6*ks 6*ks],ks);

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
        parfor kk=1:size(M,3)
            M(:,:,kk)=interp2(I1(:,:,kk),max(min(X+Ty,size(X,1)),1),max(min(Y+Tx,size(Y,1)),1));
        end
        costiter(itt+NIT*(r1-1))=norm(Idiff(:))/NPTS;
        fprintf('res=%d iter = %d, diff=%g, def=%g\n',res1(r1),itt,costiter(itt+NIT*(r1-1)),sqrt(mean((Tx(:)).^2+(Ty(:)).^2)));
    end
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



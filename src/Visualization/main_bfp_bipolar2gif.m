clear all;close all;clc;
restoredefaultpath;
bfpdir='/home/ajoshi/coding_ground/bfp';
fmridatfile='/home/ajoshi/coding_ground/bfp/data/sub-01/func/sub-01_ses-movie_task-movie_run-1_bold.32k.GOrd.mat';
outgiffile='fmri_Movie.gif';
TR=2;

%% Do not edit below this line
addpath(genpath(fullfile(bfpdir,'src')));
cmap=bipolarcmapW(100,[-1.5,1.5],'linear','br');
dfs_refL = readdfs(fullfile(bfpdir,'supp_data/bci32kleft.dfs'));
nV=length(dfs_refL.vertices);

dfs_refR = readdfs(fullfile(bfpdir,'supp_data/bci32kright.dfs'));

%load('/deneb_disk/studyforrest/sub-01-run2/fmrit_reduce3_v2.mat');
load(fmridatfile)
% dataL and dataR are fMRI data on two hemispheres, N x T
dataL=dtseries(1:nV,:);dataR=dtseries((1+nV):(2*nV),:);

dataL=normalizeData(dataL')';dataL=dataL*sqrt(size(dataL,2));
dataR=normalizeData(dataR')';dataR=dataR*sqrt(size(dataR,2));
% interpolate data to 10 fps 
tMax = 30; % desired length of the video in seconds, here 1 mins
tItvOrg = TR; % TR, here for HCP 0.72
tAxisOrg = 0:tItvOrg:tMax;
fsItp = 10; % video fps
tItvItp = 1 / fsItp;
tAxisItp = 0:tItvItp:tMax;

numFOrg = length(tAxisOrg);

% interpolation
dataL = interp1(tAxisOrg(:), dataL(:, 1:numFOrg).', tAxisItp(:)).';
dataR = interp1(tAxisOrg(:), dataR(:, 1:numFOrg).', tAxisItp(:)).';

numF = length(tAxisItp);

% initiate figure and plot the first frame
hFig = figure;
%whitebg(1,'k');
% your patch command here, p1, p2 are the handlers
% you only need to plot the surface once
subaxis(1, 2, 1,'Margin',0.01,'Spacing', 0.01, 'Padding', 0);caxis([-1.5,1.5]);colormap(cmap);
p1 = patch('faces',dfs_refL.faces,'vertices',dfs_refL.vertices,'facevertexcdata', dataL(:,1),'edgecolor','none','facecolor','interp');axis equal;axis tight;axis off;view(-90,0);camlight;material dull;lighting phong;
subaxis(1, 2, 2,'Margin',0.01,'Spacing', 0.01, 'Padding', 0);caxis([-1.5,1.5]);colormap(cmap);
pr1 = patch('faces',dfs_refR.faces,'vertices',dfs_refR.vertices,'facevertexcdata', dataR(:,1),'edgecolor','none','facecolor','interp');axis equal;axis tight;axis off;view(90,0);camlight;material dull;lighting phong;
set(gcf,'color','w', 'Units', 'Inches', 'Position', [0, 0, 14, 4], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

for k = 1:numF
    % make sure plotting in the same figure
    figure(hFig);
    
    % just replace the vertex data
    p1.FaceVertexCData = dataL(:, k);
    pr1.FaceVertexCData = dataR(:, k);

    % draw the image and take the frame and put into the video object
    drawnow;
    frame = getframe(hFig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);    
    % Write to the GIF File 
    if k == 1 
        imwrite(imind,cm,outgiffile,'gif', 'DelayTime',0.1,'Loopcount',inf); 
    else 
        imwrite(imind,cm,outgiffile,'gif','DelayTime',0.1,'WriteMode','append'); 
    end     
    disp(k);
end

close(hFig);

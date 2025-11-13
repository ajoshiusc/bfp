clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/home/ajoshi/Projects/bfp/src'));
%    1050345 rest 2

s = dlmread('/home/ajoshi/Downloads/sub_0903_finger_tapping_time_series.csv');

activation = s(:,2);
plot(activation)

gord = '/home/ajoshi/Downloads/for_parsa/bfpoutput/sub_2_13_lf/func/sub_2_13_lf_task-fingertap_bold.32k.GOrd.mat';

g = load(gord);
c=zeros(size(g.dtseries,1),1);

for j = 1:size(g.dtseries,1)

    aa = corrcoef(g.dtseries(j,:),activation);
    c(j)=aa(1,2);
    j
end

%dtseries = g.dtseries;
bfpdir='/home/ajoshi/Projects/bfp';

dfs_refL = readdfs(fullfile(bfpdir,'supp_data/bci32kleft.dfs'));
dfs_refL=smooth_cortex_fast(dfs_refL,0.1,1000);

nV=length(dfs_refL.vertices);

dfs_refR = readdfs(fullfile(bfpdir,'supp_data/bci32kright.dfs'));
dfs_refR=smooth_cortex_fast(dfs_refR,0.1,1000);

lab=load(fullfile(bfpdir,'supp_data','HCP_32k_Label.mat'));
llab=lab.brainstructure(1:nV);
rlab=lab.brainstructure((1+nV):2*nV);
% dataL and dataR are fMRI data on two hemispheres, N x T
dataL=c(1:nV,:);dataR=c((1+nV):(2*nV),:);
dataL(isnan(llab),:)=0;dataR(isnan(rlab),:)=0;

figure;
patch('faces',dfs_refL.faces,'vertices',dfs_refL.vertices,'facevertexcdata',dataL,'facecolor','interp','edgecolor','none');
axis off;axis equal;camlight;material dull;


figure;
patch('faces',dfs_refR.faces,'vertices',dfs_refR.vertices,'facevertexcdata',dataR,'facecolor','interp','edgecolor','none');
axis off;axis equal;camlight;material dull;




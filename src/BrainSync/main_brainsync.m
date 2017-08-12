%||AUM||
%||Shree Ganeshaya Namaha||
clc;clear all;close all;
addpath(genpath('../'));
g1=load('/home/ajoshi/coding_ground/bfp/data/sub-01/func/sub-01_ses-movie_task-movie_run-1_bold.32k.GOrd.mat');
g2=load('/home/ajoshi/coding_ground/bfp/data/sub-02/func/sub-02_ses-movie_task-movie_run-1_bold.32k.GOrd.mat');
sl=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kleft.dfs');
msk=load('/big_disk/ajoshi/HCP100-fMRI-NLM/HCP100-fMRI-NLM/reference/100307.LR_mask.mat');

g1.dtseries=normalizeData(g1.dtseries')';
g2.dtseries=normalizeData(g2.dtseries')';
%g1.dtseries=g1.dtseries(1:32000,:);
%g2.dtseries=g2.dtseries(1:32000,:);

rho_before=sum(g1.dtseries.*g2.dtseries,2);
dtseries_sync=brainSync(g1.dtseries',g2.dtseries')';
rho_after=sum(g1.dtseries.*dtseries_sync,2);
fprintf('Correlation before=%g, after=%g\n',mean(rho_before), mean(rho_after));


cSZ=length(msk.LR_flag);
g1.dtseries=g1.dtseries(1:cSZ,:);
g2.dtseries=g2.dtseries(1:cSZ,:);

g1.dtseries=g1.dtseries(msk.LR_flag,:);
g2.dtseries=g2.dtseries(msk.LR_flag,:);
g1.dtseries=normalizeData(g1.dtseries')';
g2.dtseries=normalizeData(g2.dtseries')';

rho_before=sum(g1.dtseries.*g2.dtseries,2);
figure;
patch('faces',sl.faces,'vertices',sl.vertices,'facevertexcdata',rho_before,'facecolor','interp','edgecolor','none');
view(-90,0);axis equal;axis off;camlight;material dull;caxis([0,1]);

dtseries_sync=brainSync(g1.dtseries',g2.dtseries')';
rho_after=sum(g1.dtseries.*dtseries_sync,2);

figure;
patch('faces',sl.faces,'vertices',sl.vertices,'facevertexcdata',rho_after,'facecolor','interp','edgecolor','none');
view(-90,0);axis equal;axis off;camlight;material dull;caxis([0,1]);


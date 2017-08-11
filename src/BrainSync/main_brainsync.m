%||AUM||
%||Shree Ganeshaya Namaha||
clc;clear all;close all;
addpath(genpath('../'));
g1=load('/home/ajoshi/coding_ground/bfp/data/sub-01/func/sub-01_ses-movie_task-movie_run-1_bold.32k.GOrd.mat');
g2=load('/home/ajoshi/coding_ground/bfp/data/sub-02/func/sub-02_ses-movie_task-movie_run-1_bold.32k.GOrd.mat');

msk=load('/big_disk/ajoshi/HCP100-fMRI-NLM/HCP100-fMRI-NLM/reference/100307.LR_mask.mat');

g1.dtseries=normalizeData(g1.dtseries')';
g2.dtseries=normalizeData(g2.dtseries')';
%g1.dtseries=g1.dtseries(1:32000,:);
%g2.dtseries=g2.dtseries(1:32000,:);

rho_before=sum(g1.dtseries.*g2.dtseries,2);
g2.dtseries=brainSync(g1.dtseries',g2.dtseries')';
rho_after=sum(g1.dtseries.*g2.dtseries,2);
fprintf('Correlation before=%g, after=%g\n',mean(rho_before), mean(rho_after));
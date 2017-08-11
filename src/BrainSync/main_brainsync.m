%||AUM||
%||Shree Ganeshaya Namaha||
clc;clear all;close all;
addpath(genpath('../'));
g1=load_nii_BIG_Lab('/home/ajoshi/coding_ground/bfp/data/sub-02/func/sub-02_ses-movie_task-movie_run-1_bold.32k.GOrd.nii.gz');
%g2=load_nii_BIG_Lab('/home/ajoshi/coding_ground/bfp/data/sub-01/func/sub-01_ses-movie_task-movie_run-1_bold.32k.GOrd.filt.nii.gz');
msk=load('/big_disk/ajoshi/HCP100-fMRI-NLM/HCP100-fMRI-NLM/reference/100307.LR_mask.mat');

s=readdfs('')
clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
bfpdir='/home/ajoshi/coding_ground/bfp';
%fmridatfile='/home/ajoshi/coding_ground/fMRIReg/src/task_diff_dyn.mat';
fmridatfile='/home/ajoshi/coding_ground/brainsync/src/dynamic_fmri/diffafter_smooth.mat';
fmridatfile='/home/ajoshi/coding_ground/brainsync/src/dynamic_fmri/pval_sm_wt.mat';

%'/deneb_disk/studyforrest_bfp/sub-01/func/sub-01_ses-movie_task-movie_run-3_bold.32k.GOrd.filt.mat';
outgiffile='pval_sm_wt.gif';
TR=0.72;%2;

addpath(genpath(fullfile(bfpdir,'src')));

make_movie_fmri(bfpdir,fmridatfile,outgiffile,TR,0.1);

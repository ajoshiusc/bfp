clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
bfpdir='/home/ajoshi/coding_ground/bfp';
fmridatfile='/home/ajoshi/coding_ground/fMRIReg/src/task_diff_dyn.mat';
%'/deneb_disk/studyforrest_bfp/sub-01/func/sub-01_ses-movie_task-movie_run-3_bold.32k.GOrd.filt.mat';
outgiffile='task_dyn_diff.gif';
TR=0.72;%2;

addpath(genpath(fullfile(bfpdir,'src')));

make_movie_fmri(bfpdir,fmridatfile,outgiffile,TR);

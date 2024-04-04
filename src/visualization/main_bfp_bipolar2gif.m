clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/ajoshi/Projects/bfp/src'));

bfpdir='/home/ajoshi/Projects/bfp';
%fmridatfile='/home/ajoshi/coding_ground/fMRIReg/src/task_diff_dyn.mat';
fmridatfile='/home/ajoshi/Downloads/for_parsa/bfpoutput/sub1/func/sub1_task-restingstate_bold.32k.GOrd.mat';

%'/deneb_disk/studyforrest_bfp/sub-01/func/sub-01_ses-movie_task-movie_run-3_bold.32k.GOrd.filt.mat';
outgiffile='for_parsa.gif';
TR=1.3;

addpath(genpath(fullfile(bfpdir,'src')));

make_movie_fmri(bfpdir,fmridatfile,outgiffile,TR,6);

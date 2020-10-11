clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/ImagePTE1/ajoshi/code_farm/bfp/src'));
bfpdir='/ImagePTE1/ajoshi/code_farm/bfp';
%fmridatfile='/home/ajoshi/coding_ground/fMRIReg/src/task_diff_dyn.mat';
fmridatfile='/ImagePTE1/ajoshi/for_abhijit/bfp_out/sub-01/func/sub-01_task-rest_bold.32k.GOrd.filt.mat';

%'/deneb_disk/studyforrest_bfp/sub-01/func/sub-01_ses-movie_task-movie_run-3_bold.32k.GOrd.filt.mat';
outgiffile='for_abhijit1.gif';
TR=2;

addpath(genpath(fullfile(bfpdir,'src')));

make_movie_fmri(bfpdir,fmridatfile,outgiffile,TR,60);

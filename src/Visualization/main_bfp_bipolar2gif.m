clear all;close all;clc;
restoredefaultpath;
addpath(genpath(fullfile(bfpdir,'src')));

bfpdir='/home/ajoshi/coding_ground/bfp';
fmridatfile='/home/ajoshi/coding_ground/bfp/data/sub-01.old/func/sub-01_ses-movie_task-movie_run-1_bold.32k.GOrd.mat';
outgiffile='fmri_Movie.gif';
TR=2;

make_movie_fmri(bfpdir,fmridatfile,outgiffile,TR);

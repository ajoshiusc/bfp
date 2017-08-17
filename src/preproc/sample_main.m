clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/bfp/src'));
% Configuration file
configfile='/home/ajoshi/coding_ground/bfp/supp_data/config.ini';
%T1 and fMRI images in NIFTI-1 format (nii.gz)
 t1='/home/ajoshi/coding_ground/bfp/data/sub-01_T1w.nii.gz';
 fmri{1}='/home/ajoshi/coding_ground/bfp/data/sub-01_ses-movie_task-movie_run-1_bold.nii.gz';
% % Outputs will be saved in the studydir/subid for each subject
 studydir='/home/ajoshi/coding_ground/bfp/data/';
 subid='sub-01';
 % fMRI files will be saved in BIDS format with sessionid in fmri names
 sessionid{1}='ses-movie_task-movie_run-1';
 TR=2;
 
 bfp(configfile,t1,fmri,studydir,subid,sessionid,TR);
 
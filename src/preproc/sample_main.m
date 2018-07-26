clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/bfp/src'));

% Set the input arguments
configfile='/home/ajoshi/bfp_sample/bfpconfig.ini';
t1='/home/ajoshi/bfp_sample/input/sub08251/anat/mprage_anonymized.nii.gz';
fmri='/home/ajoshi/bfp_sample/input/sub08251/func/rest.nii.gz';
studydir='/home/ajoshi/bfp_sample/output';
subid='sub08251';
sessionid='rest';
TR='2';
 
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid, TR);


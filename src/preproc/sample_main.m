clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/bfp/src'));

% Set the input arguments
configfile='/home/ajoshi/coding_ground/bfp/supp_data/config.ini';
t1='/home/ajoshi/Desktop/for_SY/SCD062.nii.gz';
fmri='/home/ajoshi/Desktop/for_SY/SCD062.BOLD.resting.nii.gz';
studydir='/home/ajoshi/Desktop/for_SY/bfpout';
subid='SCD062';
sessionid='rest';
TR='2';
 
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid, TR);


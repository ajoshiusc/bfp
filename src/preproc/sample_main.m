clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/bfp/src'));

% Set the input arguments
configfile='/home/ajoshi/coding_ground/bfp/supp_data/config.ini';
t1='/big_disk/ajoshi/bfp_sample/input/sub08001/anat/mprage_anonymized.nii.gz';
fmri='/big_disk/ajoshi/bfp_sample/input/sub08001/func/rest.nii.gz';
studydir='/big_disk/ajoshi/bfp_sample/output';
subid='sub08001';
sessionid='rest';
TR='2';
 
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid, TR);


clc;clear all;close all;restoredefaultpath;
addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/src'));

% Set the input arguments
configfile='/big_disk/ajoshi/bfp_sample/bfpconfig.ini';
t1='/big_disk/ajoshi/bfp_sample/input/sub08001/anat/mprage_anonymized.nii.gz';
fmri='/big_disk/ajoshi/bfp_sample/input/sub08001/func/rest.nii.gz';
studydir='/big_disk/ajoshi/bfp_sample/output';
subid='sub08001';
sessionid='rest';
TR='2';
 
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid, TR);


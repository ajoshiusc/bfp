clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/home/ajoshi/projects/BrainSuite/bfp/src'));
%    1050345 rest 2

% Set the input arguments
configfile='/home/ajoshi/Downloads/bfp_sample/bfpconfig.ini';

t1='/home/ajoshi/Downloads/bfp_sample/input/sub08001/anat/mprage_anonymized.nii.gz';
fmri='/home/ajoshi/Downloads/bfp_sample/input/sub08001/func/rest.nii.gz';

%fmri='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005.BOLD.resting.nii.gz';
studydir='/home/ajoshi/Downloads/bfp_sample/output';
subid='sub08001';
sessionid='task-restingstate';
TR='2';
 
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid,TR);

%bfp.sh config.ini input/sub08001/anat/mprage_anonymized.nii.gz /big_disk/ajoshi/bfp_sample/input/sub08001/func/rest.nii.gz /big_disk/ajoshi/bfp_sample/output sub11 rest 2

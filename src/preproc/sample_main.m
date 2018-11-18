clc;clear all;close all;restoredefaultpath;
addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/src'));

% Set the input arguments
configfile='/big_disk/ajoshi/bfp_sample/config.ini';
t1='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005_T1w.nii.gz'
%t1='/big_disk/ajoshi/bfp_sample/input/sub08001/anat/mprage_anonymized.nii.gz';
fmri='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005.BOLD.resting.nii.gz';
studydir='/home/ajoshi/Downloads/BFP_issues';
subid='ACTL005';
sessionid='rest';
TR='2';
 
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid, TR);

%bfp.sh config.ini input/sub08001/anat/mprage_anonymized.nii.gz /big_disk/ajoshi/bfp_sample/input/sub08001/func/rest.nii.gz /big_disk/ajoshi/bfp_sample/output sub11 rest 2

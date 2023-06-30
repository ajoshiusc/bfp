clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/home/ajoshi/projects/BrainSuite/bfp/src'));
%    1050345 rest 2

% Set the input arguments
configfile='/home/ajoshi/sub-0009_ses-PIOP1/config_aaj.ini';

t1='/home/ajoshi/Downloads/for_david_debug_bfp/sub-0009/anat/sub-0009_T1w.nii.gz';
fmri='/home/ajoshi/Downloads/for_david_debug_bfp/sub-0009/func/sub-0009_task-restingstate_acq-mb3_bold.nii.gz';

%fmri='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005.BOLD.resting.nii.gz';
studydir='/home/ajoshi/Downloads/for_david_debug_bfp/mybfpout/';
subid='sub-0009';
sessionid='task-restingstate';
TR='';
 
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid,TR);

%bfp.sh config.ini input/sub08001/anat/mprage_anonymized.nii.gz /big_disk/ajoshi/bfp_sample/input/sub08001/func/rest.nii.gz /big_disk/ajoshi/bfp_sample/output sub11 rest 2

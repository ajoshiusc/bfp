clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/home/ajoshi/Projects/bfp/src'));
%    1050345 rest 2

% Set the input arguments
configfile='/home/ajoshi/Projects/bfp/supp_data/config_bfp_preproc_minimal.ini';

t1='/home/ajoshi/Downloads/for_parsa/low_field_fMRI_Finger_tapping_Scan_Jan29/low-field_fMRI_Finger-tapping_Scan_Jan29_T1_MPRAGE_SAG_20240129143535_2.nii.gz';
fmri='/home/ajoshi/Downloads/for_parsa/low_field_fMRI_Finger_tapping_Scan_Jan29/low-field_fMRI_Finger-tapping_Scan_Jan29_ep2d_bold_20240129143535_11.nii.gz';

%fmri='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005.BOLD.resting.nii.gz';
studydir='/home/ajoshi/Downloads/for_parsa/bfpoutput';
subid='sub1';
sessionid='task-restingstate';
TR='1.3';

setenv('BrainSuiteMCR','/home/ajoshi/Software/MATLAB/MATLAB_Runtime/R2023a');
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid,TR);

%bfp.sh config.ini input/sub08001/anat/mprage_anonymized.nii.gz /big_disk/ajoshi/bfp_sample/input/sub08001/func/rest.nii.gz /big_disk/ajoshi/bfp_sample/output sub11 rest 2

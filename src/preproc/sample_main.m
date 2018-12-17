clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/rcf-40/ajoshi/aaj/git_sandbox/bfp/src'));

% Set the input arguments
configfile='/home/rcf-40/ajoshi/aaj/git_sandbox/bfp/supp_data/hpcconfig.ini';

t1='/home/rcf-40/ajoshi/aaj/maryland_rao_v1/TBI_INVAP401RGR/T1mni.nii.gz'
%t1='/big_disk/ajoshi/bfp_sample/input/sub08001/anat/mprage_anonymized.nii.gz';
fmri='/home/rcf-40/ajoshi/aaj/maryland_rao_v1/TBI_INVAP401RGR/rest.nii'

%fmri='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005.BOLD.resting.nii.gz';
studydir='/home/rcf-40/ajoshi/aaj/maryland_rao_v1_bfp';
subid='TBI_INVAP401RGR';
sessionid='rest';
TR='2';
 
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid, TR);

%bfp.sh config.ini input/sub08001/anat/mprage_anonymized.nii.gz /big_disk/ajoshi/bfp_sample/input/sub08001/func/rest.nii.gz /big_disk/ajoshi/bfp_sample/output sub11 rest 2

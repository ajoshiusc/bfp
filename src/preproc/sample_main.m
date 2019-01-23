clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/src'));
%    1050345 rest 2

% Set the input arguments
configfile='/big_disk/ajoshi/coding_ground/bfp/supp_data/config.ini';

t1='/big_disk/ajoshi/bfp_sample/input/1050345/session_1/anat_1/mprage_noface.nii.gz';
%t1='/big_disk/ajoshi/bfp_sample/input/sub08001/anat/mprage_anonymized.nii.gz';
fmri='/big_disk/ajoshi/bfp_sample/input/1050345/session_1/rest_1/rest.nii.gz';

%fmri='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005.BOLD.resting.nii.gz';
studydir='/big_disk/ajoshi/bfp_sample/output';
subid='1050345';
sessionid='rest';
TR='2';
 
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid, TR);

%bfp.sh config.ini input/sub08001/anat/mprage_anonymized.nii.gz /big_disk/ajoshi/bfp_sample/input/sub08001/func/rest.nii.gz /big_disk/ajoshi/bfp_sample/output sub11 rest 2

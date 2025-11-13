clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/home/ajoshi/Projects/bfp/src'));
%    1050345 rest 2

% Set the input arguments
configfile='/home/ajoshi/Projects/bfp/supp_data/config_bfp_preproc_minimal.ini';

t1='/home/ajoshi/Downloads/vol0903-tak_feb12/vol0903-tak_feb12_SAG_MPRAGE_1.1_20240212153243_2.nii.gz';
fmri='/home/ajoshi/Downloads/vol0903-tak_feb12/vol0903-tak_feb12_ep2d_bold_4mm_TR3000_20240212153243_6.nii.gz';

%fmri='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005.BOLD.resting.nii.gz';
studydir='/home/ajoshi/Downloads/for_parsa/bfpoutput';
subid='sub_2_13_lf';
sessionid='task-fingertap';
TR='3';

setenv('BrainSuiteMCR','/home/ajoshi/Software/MATLAB/MATLAB_Runtime/R2023a');
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid,TR);

%bfp.sh config.ini input/sub08001/anat/mprage_anonymized.nii.gz /big_disk/ajoshi/bfp_sample/input/sub08001/func/rest.nii.gz /big_disk/ajoshi/bfp_sample/output sub11 rest 2

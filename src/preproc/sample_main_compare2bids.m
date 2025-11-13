clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/home/ajoshi/Projects/bfp/src'));
%    1050345 rest 2

% Set the input arguments
configfile='/home/ajoshi/Projects/BrainSuiteBIDSAppSampleData/config_mysoftware.ini';

t1='/home/ajoshi/Projects/BrainSuiteBIDSAppSampleData/AOMIC-PIOP2/sub-0015/anat/sub-0015_T1w.nii.gz';
fmri='/home/ajoshi/Projects/BrainSuiteBIDSAppSampleData/AOMIC-PIOP2/sub-0015/func/sub-0015_task-restingstate_acq-seq_bold.nii.gz';

%fmri='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005.BOLD.resting.nii.gz';
studydir='/home/ajoshi/Desktop/bfpoutput_local_funcproconly';
subid='sub-0015';
sessionid='task-restingstate';
TR='0';

setenv('BrainSuiteMCR','/home/ajoshi/Software/MATLAB/MATLAB_Runtime/R2023a');
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid,TR);

%bfp.sh config.ini input/sub08001/anat/mprage_anonymized.nii.gz /big_disk/ajoshi/bfp_sample/input/sub08001/func/rest.nii.gz /big_disk/ajoshi/bfp_sample/output sub11 rest 2

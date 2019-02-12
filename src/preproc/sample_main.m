clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/home/ajoshi/coding_ground/bfp/src'));
%    1050345 rest 2

% Set the input arguments
configfile='/home/ajoshi/coding_ground/bfp/supp_data/config.ini';

t1='/home/ajoshi/Desktop/002_S_4521/MPRAGE_SENSE2/2012-02-17_17_05_51.0/S141306/co20120217_165208MPRAGESENSE2SENSEs401a1004.nii.gz';
fmri='/home/ajoshi/Desktop/002_S_4521/Resting_State_fMRI/2012-02-17_16_52_08.0/S141311/20120217_165208RestingStatefMRIs501a1005.nii.gz';

%fmri='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005.BOLD.resting.nii.gz';
studydir='/deneb_disk';
subid='002_S_4521';
sessionid='rest';
TR='3';
 
% Call the bfp function
bfp(configfile, t1, fmri, studydir, subid, sessionid, TR);

%bfp.sh config.ini input/sub08001/anat/mprage_anonymized.nii.gz /big_disk/ajoshi/bfp_sample/input/sub08001/func/rest.nii.gz /big_disk/ajoshi/bfp_sample/output sub11 rest 2

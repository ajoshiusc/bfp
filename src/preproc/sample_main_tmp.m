clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/bfp/src'));
% Configuration file
configfile='/home/ajoshi/coding_ground/bfp/supp_data/config.ini';
%T1 and fMRI images in NIFTI-1 format (nii.gz)

for jj=1:21
    subid=sprintf('sub-%02d',jj);
    t1=sprintf('/deneb_disk/studyforrest_bfp/%s_T1w.nii.gz',subid);
    
    fmri{1}=sprintf('/deneb_disk/studyforrest_bfp/%s_ses-movie_task-movie_run-3_bold.nii.gz',subid);

    if ~exist(t1,'file') || ~exist(fmri{1},'file')
        continue;
    end

    % Outputs will be saved in the studydir/subid for each subject
    studydir='/deneb_disk/studyforrest_bfp';

    % fMRI files will be saved in BIDS format with sessionid in fmri names
    sessionid{1}='ses-movie_task-movie_run-3';
    TR=2;    
    bfp(configfile,t1,fmri,studydir,subid,sessionid,TR);
end

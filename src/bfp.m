%% BFP: BrainSuite fMRI Pipeline

%% User Inputs

t1 = '~/coding_ground/bfp/data/sub-01_ses-movie_task-movie_run-1_defacemask.nii.gz';
fmri{1} = '~/coding_ground/bfp/data/sub-01_ses-movie_task-movie_run-1_bold.nii.gz';
subbasename = '~/coding_ground/bfp/data/sub-01-run1/sub-01-run1';
BrainSuitePath='/home/ajoshi/BrainSuite17a';
RMFLAG=1;
fprintf('Processing T1=%s\n',t1)
for ind=1:length(fmri)
    fprintf('Processing fMRI=%s\n',fmri{ind})
end

%% Create Directory Structure
subdir=fileparts(subbasename);
if RMFLAG && exist(subdir,'dir')
    rmdir(subdir,'s')
end

mkdir(subbasename)
anatDir=fullfile(subbasename,'anat');
mkdir(anatDir);
copyfile(t1,anatDir)
for ind = 1:length(fmri)
    outdir=fullfile(subbasename,sprintf('func-%d',ind));
    mkdir(outdir);
    copyfile(fmri{ind},outdir);
end

%%
%% Skull Strip MRI
% 
%% Coregister fMRI to MNI Space
% 
%% Run BrainSuite
% 
%% Run Batch_Process Pipeline
% 
%% Transfer data to Surface and then to USCBrain atlas surface
% 
%% Produce Surface Grayordinates
% 
%% Apply fNIRT or Inverse map to map vol data to USCBrain atlas
% 
%% Generate Volumetric Grayordinates
% 
% 
%
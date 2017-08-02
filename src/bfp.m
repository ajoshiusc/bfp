%% BFP: BrainSuite fMRI Pipeline

%% User Inputs

t1 = '/home/ajoshi/coding_ground/bfp/data/sub-01_T1w.nii.gz';
fmri{1} = '/home/ajoshi/coding_ground/bfp/data/sub-01_ses-movie_task-movie_run-1_bold.nii.gz';
subbasename = '~/coding_ground/bfp/data/sub-01-run1/sub-01-run1';
BrainSuitePath='/home/ajoshi/BrainSuite17a';
RMFLAG=1;

fprintf('OS:%s\n',computer);
if ~strcmp(computer,'GLNXA64')
    error('OS is not supported!')
end

fprintf('Processing T1=%s\n',t1)
for ind=1:length(fmri)
    fprintf('Processing fMRI=%s\n',fmri{ind})
end

%% Create Directory Structure
fprintf('Creating Directory Structure\n');
subdir=fileparts(subbasename);
if RMFLAG && exist(subdir,'dir')
    rmdir(subdir,'s')
end

mkdir(subdir)
anatDir=fullfile(subdir,'anat');
fprintf('Creating Dir:%s\n',anatDir);
mkdir(anatDir);
t1new=fullfile(anatDir,'orig.nii.gz');
copyfile(t1,t1new);
for ind = 1:length(fmri)
    outdir=fullfile(subdir,sprintf('func-%d',ind));
    fprintf('Creating Dir:%s\n',outdir);
    mkdir(outdir);    
    copyfile(fmri{ind},outdir);
end
%% Skull Strip MRI
fprintf('Performing Skull Extraction\n');
bse=fullfile(BrainSuitePath,'bin','bse');
bseout=fullfile(anatDir,'orig.bse.nii.gz');
cmd=sprintf('%s --auto -i %s -o %s',bse,t1new,bseout);
unix(cmd);
 
%% Coregister t1 to MNI Space
bsenew=fullfile(anatDir,'mprage_skullstripped.nii.gz');

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
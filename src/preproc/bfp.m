%% BFP: BrainSuite fMRI Pipeline
% This pipeline takes fMRI and anatomical data and processes them using a series 
% of scripts from BrainSuite, FSL and AFNI. The functional processing script is 
% based on batch_process. sh script from fcon1000.
%% Inputs/Outputs
% * Inputs: Unprocessed T1/fMRI images in NIFTI format
% * Outputs: Preprocessed fMRI data in 32k Grayordinates
%% User Inputs [Edit this section]

% clc;clear all;close all;restoredefaultpath;
% addpath(genpath('/home/ajoshi/coding_ground/bfp/src'));
% % Configuration file
% configfile='/home/ajoshi/coding_ground/bfp/supp_data/config.ini';
% %T1 and fMRI images in NIFTI-1 format (nii.gz)
%  t1='~/coding_ground/bfp/data/sub-01_T1w.nii.gz';
%  fmri{1}='~/coding_ground/bfp/data/sub-01_ses-movie_task-movie_run-1_bold.nii.gz';
% % % Outputs will be saved in the studydir/subid for each subject
%  studydir='~/coding_ground/bfp/data/';
%  subid='sub-01';
%  % fMRI files will be saved in BIDS format with sessionid in fmri names
%  sessionid{1}='ses-movie_task-movie_run-1';
%  TR=2;
%% **DO NOT EDIT BELOW THIS****
%% Commonly used Variables
%%
function bfp(configfile,t1,fmri,studydir,subid,sessionid,TR)

if nargin ~=7
   disp('Incorrect number of input arguments (7 required)');
   disp('Usage:');
   disp('bfp configfile t1 fmri studydir subid sessionid TR');
   disp('t1: ' )
   
   
   error('exiting');
end
disp('hi');
p=inputParser;

addRequired(p,'configfile',@ischar);
addRequired(p,'t1',@ischar);
addRequired(p,'fmri',@iscellstr);
addRequired(p,'studydir',@ischar);
addRequired(p,'subid',@ischar);
addRequired(p,'sessionid',@ischar);
addRequired(p,'TR',@(x) isnumeric(x)||ischar(x));
parse(p,configfile,t1,fmri,studydir,subid,sessionid,TR);
%%
subdir=fullfile(studydir,subid);
anatDir=fullfile(subdir,'anat');
funcDir=fullfile(subdir,sprintf('func'));
subbasename=fullfile(anatDir,sprintf('%s_T1w',subid));

if ischar(class(fmri)) % This is for command line input when fmri data is cell string
    if contains(fmri,'{')
        eval(['fmri = ' fmri]);
    end
else
    fmri{1}=fmri;
end

if ischar(class(sessionid)) % This is for command line input when fmri data is cell string
    if contains(sessionid,'{')
        eval(['sessionid = ' fmri]);
    end
else
    sessionid{1}=fmri;
end

for i = 1:size(fmri,2)
    fprintf('fmri %d:%s\n',fmri{i});
end

for i = 1:size(sessionid,2)
    fprintf('sessionid %d:%s\n',sessionid{i});
end
%% Check if OS is supported
%%
fprintf('OS:%s\n',computer);
if ~strcmp(computer,'GLNXA64')
    error('OS %s is not supported, please use Linux 64 bit computer to run BFP!',computer)
end
%% Read configuration file and set environment variables
%%
fprintf("# Starting BFP Run");
if ~exist(configfile,'file')
    error('Config file: %s \n: File does not exist\n',configfile);
end

fprintf("## Reading config file\n");
config=ini2struct(configfile);
fprintf(" done\n");
%% Setting up the environment
%%
fprintf("## Setting up the environment\n");
setenv('PATH', [getenv('PATH'),':',config.FSLPATH,':',config.FSLPATH,'/bin']);
setenv('PATH', [getenv('PATH'),':',config.AFNIPATH,':',config.AFNIPATH,'/bin']);
setenv('FSLOUTPUTTYPE',config.FSLOUTPUTTYPE);
setenv('BrainSuiteDir',config.BrainSuitePath);
setenv('LD_LIBRARY_PATH', [config.LD_LIBRARY_PATH]);
BrainSuitePath=config.BrainSuitePath;
BFPPATH=config.BFPPATH;
bst_exe=fullfile(BFPPATH,'supp_data','cortical_extraction_nobse.sh');
svreg_exe=fullfile(BrainSuitePath,'svreg','bin','svreg.sh');
BCIbasename=fullfile(BrainSuitePath,'svreg','BCI-DNI_brain_atlas','BCI-DNI_brain');
BSA=fullfile(BrainSuitePath,'svreg','BrainSuiteAtlas1','mri.bfc.nii.gz');
GOrdSurfIndFile=fullfile(BFPPATH,'supp_data','bci_grayordinates_surf_ind.mat');
GOrdVolIndFile=fullfile(BFPPATH,'supp_data','bci_grayordinates_vol_ind.mat');
RMFLAG='0';%depricated
nuisance_template=fullfile(BFPPATH,'supp_data','nuisance.fsf');
func_prepro_script=fullfile(BFPPATH,'supp_data','func_preproc.sh');

fprintf(" done\n");
%% Create Directory Structure
% This directory structure is in BIDS format
%%
fprintf('## Creating Directory Structure\n');
if str2double(RMFLAG) && exist(subdir,'dir')
    %rmdir(subdir,'s');
    fprintf('RMFLAG is a depricated and is a placeholder.\nDeletion of subject directory by this flag is disabled for now.\n')
end
if exist(subdir,'dir')
    error('The subject directory (%s) exists!\n Please check that directory for previous runs and delete it if necessary\n',subdir);
end
    
mkdir(subdir)
fprintf('Creating Dir:%s\n',anatDir);
mkdir(anatDir);
t1hires=fullfile(anatDir,sprintf('%s.orig.hires.nii.gz',subid));
t1ds=fullfile(anatDir,sprintf('%s_T1w.ds.orig.nii.gz',subid)); 
copyfile(t1,t1hires);
fprintf('Creating Dir:%s\n',funcDir);
mkdir(funcDir);    
for ind = 1:length(fmri)
    copyfile(fmri{ind},fullfile(funcDir,sprintf('%s_%s_bold.nii.gz',subid,sessionid{ind})));
end
fprintf('done');
%% Generate 3mm BrainSuiteAtlas1 brain as a standard template
% This is used a template for fMRI data
%%
cmd=sprintf('flirt -ref %s -in %s -out %s -applyisoxfm 3',BSA,BSA,fullfile(funcDir,'standard.nii.gz'));
fprintf('Creating 3mm isotropic standard brain\n');
unix(cmd);
fprintf('done\n');
%% Resample T1w image to 1mm cubic resolution
% BrainSuite works best at this resolution
%%
fprintf('## Resample T1w image to 1mm cubic resolution \n')
cmd=sprintf('flirt -in %s -ref %s -out %s -applyisoxfm 1',t1hires,t1hires,t1ds);
unix(cmd);
fprintf('done\n');
%% Skull Strip MRI
%%
fprintf('## Performing Skull Extraction\n');
bse=fullfile(BrainSuitePath,'bin','bse');
bseout=fullfile(anatDir,sprintf('%s_T1w.ds.orig.bse.nii.gz',subid));
cmd=sprintf('%s --auto -i %s -o %s',bse,t1ds,bseout);
unix(cmd);
fprintf('done');
%% Coregister t1 to MNI Space
%%
fprintf('## Coregister t1 to MNI Space\n');
bsenew=fullfile(anatDir,sprintf('%s_T1w.nii.gz',subid));
cmd=sprintf('flirt -ref %s -in %s -out %s',BSA,bseout,bsenew);
unix(cmd);
bsemask=fullfile(anatDir,sprintf('%s_T1w.mask.nii.gz',subid));
cmd=sprintf('fslmaths %s -thr 0 -bin -mul 255 %s -odt char',bsenew,bsemask);
unix(cmd);
bsenew2=fullfile(anatDir,sprintf('%s_T1w.bse.nii.gz',subid));
copyfile(bsenew,bsenew2);
fprintf('done\n');
%% Run BrainSuite and SVReg
%%
fprintf('## Running BrainSuite CSE\n');
cmd=sprintf('%s %s',bst_exe,subbasename);
unix(cmd);
fprintf('Running SVReg');
cmd=sprintf('%s %s %s',svreg_exe,subbasename,BCIbasename);
unix(cmd);
fprintf('done');
%% Run Batch_Process Pipeline for fMRI
%%
fprintf('## Run fmri preprocessing script\n');
for ind=1:length(fmri)
    fmribasename=fullfile(funcDir,sprintf('%s_%s_bold',subid,sessionid{ind}));
    unix(sprintf('%s %s %s %s %s %s',func_prepro_script,subbasename,fmribasename,funcDir,num2str(TR),nuisance_template));
end
fprintf('done\n');
%% Grayordinate representation
% Transfer data to surface and then to USCBrain atlas surface and produce surface 
% grayordinates and then Transfer data to volumetric grayordinates.
% 
% The Grayordinate data is in the same format as HCP data on 32k surfaces. 
% The filename of grayordinate data is fmri_bold.32k.GOrd.nii.gz
%%
fprintf('## Transferring data from subject to atlas...\n');
for ind = 1:length(fmri)
    fmri2surfFile=fullfile(funcDir,sprintf('%s_%s_bold2surf.mat',subid,sessionid{ind}));    
    GOrdSurfFile=fullfile(funcDir,sprintf('%s_%s_bold2surf_GOrd.mat',subid,sessionid{ind}));
    fmri2standard=fullfile(funcDir,sprintf('%s_%s_bold_res2standard.nii.gz',subid,sessionid{ind})); 
    GOrdVolFile=fullfile(funcDir,sprintf('%s_%s_bold2Vol_GOrd.mat',subid,sessionid{ind}));
    GOrdFile=fullfile(funcDir,sprintf('%s_%s_bold.32k.GOrd.mat',subid,sessionid{ind}));
    resample2surf(subbasename,fmri2standard,fmri2surfFile);
    generateSurfGOrdfMRI(GOrdSurfIndFile,fmri2surfFile,GOrdSurfFile);
    generateVolGOrdfMRI(GOrdVolIndFile,subbasename,fmri2standard,GOrdVolFile);
    combineSurfVolGOrdfMRI(GOrdSurfFile,GOrdVolFile,GOrdFile);      
end
fprintf('The grayordinates file is: %s\n',GOrdFile);
%% tNLMPDF
% This part of the code takes grayordinate data generated by the previous pipeline 
% and performs tNLMPdf filtering. 
% 
% The output is stored in <fmri fmri>.32k.GOrd.filt.nii.gz
%%
fprintf('## tNLMPDF Filtering...\n');
for ind = 1:length(fmri)
    GOrdFile=fullfile(funcDir,sprintf('%s_%s_bold.32k.GOrd.mat',subid,sessionid{ind}));
    GOrdFiltFile=fullfile(funcDir,sprintf('%s_%s_bold.32k.GOrd.filt.mat',subid,sessionid{ind}));
    tNLMPDFGOrdfMRI(GOrdFile,GOrdFiltFile);
end
fprintf('The output filtered grayordinates file is: %s\n All done!! Good Night!\n',GOrdFiltFile);

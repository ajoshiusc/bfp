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
disp('Starting BFP Processing');
p=inputParser;

addRequired(p,'configfile',@ischar);
addRequired(p,'t1',@ischar);
addRequired(p,'fmri',@(x) ischar(x)||iscellstr(x));
addRequired(p,'studydir',@ischar);
addRequired(p,'subid',@ischar);
addRequired(p,'sessionid',@(x) ischar(x)||iscellstr(x));
addRequired(p,'TR',@(x) isnumeric(x)||ischar(x));
%addOptional(p,'contprevrun',@(x) isnumeric(x)||ischar(x));

parse(p,configfile,t1,fmri,studydir,subid,sessionid,TR);
%%
subdir=fullfile(studydir,subid);
anatDir=fullfile(subdir,'anat');
funcDir=fullfile(subdir,sprintf('func'));
subbasename=fullfile(anatDir,sprintf('%s_T1w',subid));

if ischar(fmri) % This is for command line input when fmri data is cell string
    if strfind(fmri,'{')
        eval(['fmri = ' fmri]);
    else
        tmp=fmri;clear fmri;
        fmri{1}=tmp;
    end
    
end

if ischar(sessionid) % This is for command line input when fmri data is cell string
    if strfind(sessionid,'{')
        eval(['sessionid = ' sessionid]);
    else
        tmp=sessionid;clear sessionid
        sessionid{1}=tmp;
    end
end

for i = 1:size(fmri,2)
    fprintf('fmri %d:%s\n',i,fmri{i});
end

for i = 1:size(sessionid,2)
    fprintf('sessionid %d:%s\n',i,sessionid{i});
end

%% Read configuration file and set environment variables
%%
fprintf('# Starting BFP Run\n');
if ~exist(configfile,'file')
    error('Config file: %s \n: File does not exist\n',configfile);
end

fprintf('## Reading config file\n');
config=ini2struct(configfile);
fprintf(' done\n');
%% Setting up the environment
%%
fprintf('## Setting up the environment\n');
setenv('PATH', [getenv('PATH'),':',config.FSLPATH,':',config.FSLPATH,'/bin']);
setenv('PATH', [getenv('PATH'),':',config.AFNIPATH,':',config.AFNIPATH,'/bin']);
setenv('FSLOUTPUTTYPE',config.FSLOUTPUTTYPE);
setenv('FSLDIR', config.FSLPATH);
setenv('BrainSuiteDir',config.BrainSuitePath);
setenv('LD_LIBRARY_PATH', [config.LD_LIBRARY_PATH]);



BrainSuitePath=config.BrainSuitePath;
BFPPATH=config.BFPPATH;
bst_exe=fullfile(BFPPATH,'supp_data','cortical_extraction_nobse.sh');
svreg_exe=fullfile(BrainSuitePath,'svreg','bin','svreg.sh');
thicknessPVC_exe=fullfile(BrainSuitePath,'svreg','bin','thicknessPVC.sh');

BCIbasename=fullfile(BrainSuitePath,'svreg','BCI-DNI_brain_atlas','BCI-DNI_brain');
ATLAS=fullfile(BrainSuitePath,'svreg','BCI-DNI_brain_atlas','BCI-DNI_brain.bfc.nii.gz');
GOrdSurfIndFile=fullfile(BFPPATH,'supp_data','bci_grayordinates_surf_ind.mat');
GOrdVolIndFile=fullfile(BFPPATH,'supp_data','bci_grayordinates_vol_ind.mat');
nuisance_template=fullfile(BFPPATH,'supp_data','nuisance.fsf');
func_prepro_script=fullfile(BFPPATH,'supp_data','func_preproc.sh');
fwhm=config.FWHM;
hp=config.HIGHPASS;
lp=config.LOWPASS;
continueRun=str2double(config.CONTINUERUN);
FSLRigid=1; %PUT THIS IN CONFIG FILE

if isfield(config, 'MultiThreading')
    config.MultiThreading=str2double(config.MultiThreading);
    if isnan(config.MultiThreading)
        config.MultiThreading = 0;
    end
else
    config.MultiThreading=1;
end

if isfield(config, 'T1SpaceProcessing')
    config.T1SpaceProcessing=str2double(config.T1SpaceProcessing);
    if isnan(config.T1SpaceProcessing)
        config.T1SpaceProcessing = 1;
    end
else
    config.T1SpaceProcessing = 1;
end

if isfield(config, 'EnabletNLMPdfFiltering')
    config.EnabletNLMPdfFiltering=str2double(config.EnabletNLMPdfFiltering);
    if isnan(config.EnabletNLMPdfFiltering)
        config.EnabletNLMPdfFiltering = 0;
    end
else
    config.EnabletNLMPdfFiltering=1;
end

if isfield(config, 'EnableShapeMeasures')
    config.EnableShapeMeasures=str2double(config.EnableShapeMeasures);
    if isnan(config.EnableShapeMeasures)
        config.EnableShapeMeasures = 0;
    end
else
    config.EnableShapeMeasures = 1;
end

if isfield(config, 'FSLRigid')
    config.FSLRigid=str2double(config.FSLRigid);
    if isnan(config.FSLRigid)
        config.FSLRigid = 0;
    end
else
    config.FSLRigid = 1;
end
fprintf(' done\n');

%% Binaries that are compiled in MATLAB
nii2int16_bin = fullfile(BFPPATH, 'nii2int16.sh');
resample2surf_bin = fullfile(BFPPATH, 'resample2surf.sh');
generateGOrdSCT_bin = fullfile(BFPPATH, 'generateGOrdSCT.sh');
generateSurfGOrdfMRI_bin = fullfile(BFPPATH, 'generateSurfGOrdfMRI.sh');
generateVolGOrdfMRI_bin = fullfile(BFPPATH, 'generateVolGOrdfMRI.sh');
combineSurfVolGOrdfMRI_bin = fullfile(BFPPATH, 'combineSurfVolGOrdfMRI.sh');
tNLMPDFGOrdfMRI_bin = fullfile(BFPPATH,'tNLMPDFGOrdfMRI.sh');

%% Create Directory Structure
% This directory structure is in BIDS format
%%
fprintf('## Creating Directory Structure\n');
if exist(subdir,'dir') && continueRun==0
    error('The subject directory (%s) exists!\n Please check that directory for previous runs and delete it if necessary\n',subdir);
end

if ~exist(subdir,'dir')
    mkdir(subdir)
end
fprintf('Creating Dir:%s\n',anatDir);

if ~exist(anatDir,'dir')
    mkdir(anatDir);
end

t1hires=fullfile(anatDir,sprintf('%s.orig.hires.nii.gz',subid));
t1ds=fullfile(anatDir,sprintf('%s_T1w.ds.orig.nii.gz',subid));

if config.T1SpaceProcessing
    ATLAS=t1hires;
end


if ~exist(t1hires,'file')
    copyfile(t1,t1hires);
end

fprintf('Creating Dir:%s\n',funcDir);
if ~exist(funcDir,'dir')
    mkdir(funcDir);
end

fprintf('Copying fMRI files\n');
for ind = 1:length(fmri)
    if ~exist(fullfile(funcDir,sprintf('%s_%s_bold.nii.gz',subid,sessionid{ind})),'file')
        copyfile(fmri{ind},fullfile(funcDir,sprintf('%s_%s_bold.nii.gz',subid,sessionid{ind})));
    else
        fprintf('subject=%s session=%s : Already done\n',subid,sessionid{ind});
    end
end
fprintf('done\n');


%% Generate 1mm BCI-DNI_brain brain as a standard template
% This is used a template for anatomical T1 data
%%
if ~config.T1SpaceProcessing
    ATLAS_DS=fullfile(anatDir,'standard1mm.nii.gz');
    cmd=sprintf('flirt -ref %s -in %s -out %s -applyisoxfm 1',ATLAS,ATLAS,ATLAS_DS);
    fprintf('Creating 1mm isotropic standard brain\n');
else
    ATLAS_DS=fullfile(anatDir,'standard.nii.gz');
    cmd=sprintf('cp %s %s', ATLAS, ATLAS_DS);
    fprintf('Creating a copy of the t1 image to be used as the standard brain\n');
end

%why is this created?
if ~exist(ATLAS_DS,'file')
    unix(cmd);
else
    fprintf('Already ');
end
fprintf('done\n');


%% Resample T1w image to 1mm cubic resolution
% BrainSuite works best at this resolution
%%
if ~(config.T1SpaceProcessing)
    fprintf('## Resample T1w image to 1mm cubic resolution \n')
    cmd=sprintf('flirt -in %s -ref %s -out %s -applyisoxfm 1',t1hires,t1hires,t1ds);
else
    t1ds = t1hires;
end

if ~exist(t1ds,'file')
    unix(cmd);
    cmd=sprintf('%s %s %s', nii2int16_bin, t1ds, t1ds);
    unix(cmd);
%    nii2int16(t1ds,t1ds);
else
    fprintf('Already ');
end
fprintf('done\n');

%% Skull Strip MRI
%%
fprintf('## Performing Skull Extraction\n');
bse=fullfile(BrainSuitePath,'bin','bse');

if config.T1SpaceProcessing
    bseout=fullfile(anatDir,sprintf('%s_T1w.orig.bse.nii.gz',subid));    
else
    bseout=fullfile(anatDir,sprintf('%s_T1w.ds.orig.bse.nii.gz',subid));
end

cmd=sprintf('%s --auto --trim -i %s -o %s',bse,t1ds,bseout);
if ~exist(bseout,'file')
    unix(cmd);
else
    fprintf('Already ');
end
fprintf('done\n');
%% Coregister t1 to BCI-DNI Space
%%
fprintf('## Coregister t1 to BCI-DNI Space\n');
bsenew=fullfile(anatDir,sprintf('%s_T1w.nii.gz',subid));

if config.T1SpaceProcessing
    cmd = sprintf('cp %s %s', t1, bsenew);
else
    cmd=sprintf('flirt -ref %s -in %s -out %s',ATLAS_DS,bseout,bsenew);
end

if ~exist(bsenew,'file')
    unix(cmd);
    
%    nii2int16(bsenew, bsenew, 0);
    cmd=sprintf('nii2int16.sh %s %s %s', bsenew, bsenew, '0');
    unix(cmd);
%    
end

bsenew2=fullfile(anatDir,sprintf('%s_T1w.bse.nii.gz',subid));

if ~exist(bsenew2,'file')
    if config.T1SpaceProcessing
        copyfile(bseout,bsenew2);

        %nii2int16(bsenew2, bsenew2,0);
        cmd=sprintf('nii2int16.sh %s %s %s', bsenew2, bsenew2, '0');
        unix(cmd);

    else
        copyfile(bsenew,bsenew2);
    end
else
    fprintf('Already ');
end
fprintf('done\n');

bsemask=fullfile(anatDir,sprintf('%s_T1w.mask.nii.gz',subid));
cmd=sprintf('fslmaths %s -thr 0 -bin -mul 255 %s -odt char',bsenew2,bsemask);
if ~exist(bsemask,'file')
    unix(cmd);
end

%% Run BrainSuite and SVReg
%%
fprintf('## Running BrainSuite CSE\n');
cmd=sprintf('%s %s',bst_exe,subbasename);
if ~exist([subbasename,'.right.pial.cortex.dfs'],'file')
    unix(cmd);
end
fprintf('Running SVReg');
if (config.MultiThreading == 1)
    cmd=sprintf('%s %s %s',svreg_exe,subbasename,BCIbasename);
else
    cmd=sprintf('%s %s %s -U',svreg_exe,subbasename,BCIbasename);
end

if ~exist([subbasename,'.svreg.label.nii.gz'],'file')
    unix(cmd);
else
    fprintf('Already ');
end
fprintf('done\n');

%% Generate 3mm BCI-DNI_brain brain as a standard template
% This is used a template for fMRI data
%%
if config.T1SpaceProcessing
    ATLAS=fullfile(anatDir,sprintf('%s_T1w.bfc.nii.gz',subid));
end

cmd=sprintf('flirt -ref %s -in %s -out %s -applyisoxfm 3',ATLAS,ATLAS,fullfile(funcDir,'standard.nii.gz'));
fprintf('Creating 3mm isotropic standard brain\n');
if ~exist(fullfile(funcDir,'standard.nii.gz'),'file')
    unix(cmd);
else
    fprintf('Already ');
end
fprintf('done\n');

%% Run Batch_Process Pipeline for fMRI
%%
fprintf('## Run fmri preprocessing script\n');
for ind=1:length(fmri)
    fmribasename=fullfile(funcDir,sprintf('%s_%s_bold',subid,sessionid{ind}));
    if ~exist([fmribasename,'_res2standard.nii.gz'],'file')
        func_preproc(subbasename,fmribasename,funcDir,num2str(TR),nuisance_template,fwhm,hp,lp,FSLRigid);

%        unix(sprintf('%s %s %s %s %s %s %s %s %s',func_prepro_script,subbasename,fmribasename,funcDir,num2str(TR),nuisance_template,fwhm,hp,lp));
    else
        fprintf('fMRI %s : Already done\n',fmribasename);
    end
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
    fprintf('Resampling fMRI to surface\n')
    if ~exist(fmri2surfFile,'file') && ~exist(GOrdSurfFile,'file') && ~exist(GOrdFile,'file')
%       resample2surf(subbasename,fmri2standard,fmri2surfFile,config.MultiThreading);
        cmd = sprintf('%s %s %s %d', resample2surf_bin, fmri2standard, fmri2surfFile, config.MultiThreading);
        unix(cmd);
    else
        fprintf('Already ');
    end
    
    fprintf('done\n');
    fprintf('Generating Surface Grayordinates\n');
    if ~exist(GOrdSurfFile,'file') && ~exist(GOrdFile,'file')
%       generateSurfGOrdfMRI(GOrdSurfIndFile,fmri2surfFile,GOrdSurfFile);

        cmd = sprintf('%s %s %s', generateSurfGOrdfMRI_bin, GOrdSurfIndFile, fmri2surfFile, GOrdSurfFile);
        unix(cmd);
        
        % The surf file is very large, deleting to save space
        delete(fmri2surfFile);
    else
        fprintf('Already ');
    end
    fprintf('done\n');
    fprintf('Generating Volume Grayordinates\n');
    if ~exist(GOrdVolFile,'file') && ~exist(GOrdFile,'file')       
%        generateVolGOrdfMRI(GOrdVolIndFile,subbasename,fmri2standard,GOrdVolFile);

        cmd = sprintf('%s %s %s %s %s', generateVolGOrdfMRI_bin, GOrdVolIndFile, subbasename, fmri2standard, GOrdVolFile);
        unix(cmd);

    else
        fprintf('Already ');
    end
    fprintf('done\n');
    fprintf('Combining Surface and Volume Grayordinates\n');
    if ~exist(GOrdFile,'file')
%        combineSurfVolGOrdfMRI(GOrdSurfFile,GOrdVolFile,GOrdFile);

        cmd = sprintf('%s %s %s %s', combineSurfVolGOrdfMRI_bin, GOrdSurfFile, GOrdVolFile, GOrdFile);
        unix(cmd);

        delete(GOrdSurfFile);
        delete(GOrdVolFile);
    else
        fprintf('Already ');
    end
    fprintf('done\n');
end
fprintf('The grayordinates file is: %s\n',GOrdFile);
%% tNLMPDF
% This part of the code takes grayordinate data generated by the previous pipeline
% and performs tNLMPdf filtering.
%
% The output is stored in <fmri fmri>.32k.GOrd.filt.nii.gz
%%
if config.EnabletNLMPdfFiltering>0
    
    fprintf('## tNLMPDF Filtering...\n');
    for ind = 1:length(fmri)
        GOrdFile=fullfile(funcDir,sprintf('%s_%s_bold.32k.GOrd.mat',subid,sessionid{ind}));
        GOrdFiltFile=fullfile(funcDir,sprintf('%s_%s_bold.32k.GOrd.filt.mat',subid,sessionid{ind}));
        fprintf('tNLMPdf filtering for subject = %s session = %s\n',subid,sessionid{ind});
        if ~exist(GOrdFiltFile,'file')
%          tNLMPDFGOrdfMRI(GOrdInFile,GOrdOutFile,config.fpr,config.memory,config.MultiThreading);
           cmd = sprintf('%s %s %s %s %s %s %d', tNLMPDFGOrdfMRI_bin, GOrdInFile, GOrdOutFile, config.fpr, config.memory, config.MultiThreading);
           unix(cmd);

        else
            fprintf('Already ');
        end
        fprintf('done\n');
    end
    fprintf('The output filtered grayordinates files are ready !! \n Good Night!\n');
end

if config.EnableShapeMeasures>0
    
    fprintf('Running thicknessPVC');
    cmd=sprintf('%s %s %s',thicknessPVC_exe,subbasename);
    subdir = fileparts(subbasename);
    if ~exist(fullfile(subdir, 'atlas.pvc-thickness_0-6mm.right.mid.cortex.dfs'),'file')
        unix(cmd);
    else
        fprintf('Already computed thicknessPVC');
    end
    
    if ~exist([subbasename,'.SCT.GOrd.mat'],'file')
        generateGOrdSCT(subbasename, GOrdSurfIndFile);
        generateGOrdSCT_bin
    else
        fprintf('Already done SCT');
    end
    
    fprintf('done\n');
end


fprintf('All done!\n');

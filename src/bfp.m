%% BFP: BrainSuite fMRI Pipeline
%% TBD make names BIDS compatible
% * I have to check this
%%


%% User Inputs
%%
t1='/home/ajoshi/coding_ground/bfp/data/sub-01_T1w.nii.gz';
fmri{1}='/home/ajoshi/coding_ground/bfp/data/sub-01_ses-movie_task-movie_run-1_bold.nii.gz';
sessionid{1}='ses-movie_task-movie_run-1';
studydir='~/coding_ground/bfp/data/';
subid='sub-01';
configfile='config.ini';
fprintf('OS:%s\n',computer);
if ~strcmp(computer,'GLNXA64')
    error('OS is not supported!')
end

%% Read configuration file 
%%
config=ini2struct(configfile);

setenv('LD_LIBRARY_PATH',config.ld_library_path);
BrainSuitePath=config.brainsuitepath;
bst_exe=config.bst_exe;
svreg_exe=config.svreg_exe;
BCIbasename=config.brainsuitepath;
BSA=config.bsa;
GOrdSurfIndFile=config.gordsurfindfile;
GOrdVolIndFile=config.gordvolindfile;
RMFLAG=config.rmflag;
fprintf('Processing T1=%s\n',t1)
for ind=1:length(fmri)
    fprintf('Processing fMRI=%s\n',fmri{ind})
end
%% Create Directory Structure
% This directory structure is in BIDS format
%%
fprintf('Creating Directory Structure\n');
subdir=outsubdir;
if RMFLAG && exist(subdir,'dir')
    rmdir(subdir,'s')
end

mkdir(subdir)
anatDir=fullfile(subdir,'anat');
fprintf('Creating Dir:%s\n',anatDir);
mkdir(anatDir);
t1hires=fullfile(anatDir,sprintf('%s.orig.hires.nii.gz',subid));
t1ds=fullfile(anatDir,sprintf('%s_T1w.ds.orig.nii.gz',subid)); 
copyfile(t1,t1hires);
%% Generate 3mm BrainSuiteAtlas1 brain as a standard template
% This is used a template for fMRI data
%%
cmd=sprintf('flirt -ref %s -in %s -out %s -applyisoxfm 3',BSA,BSA,fullfile(anatDir,'standard.nii.gz'));
fprintf('Creating 3mm isotropic standard brains\n');
unix(cmd);
%% Resample T1w image to 1mm cubic resolution. 
% BrainSuite works best at this resolution

cmd=sprintf('flirt -in %s -ref %s -out %s -applyisoxfm 1',t1hires,t1hires,t1ds);
unix(cmd);
outdir=fullfile(subdir,sprintf('func'));
fprintf('Creating Dir:%s\n',outdir);
mkdir(outdir);    

for ind = 1:length(fmri)
    copyfile(fmri{ind},fullfile(outdir,fprintf('%s_%s_bold.nii.gz',subid,sessionid{ind})));
end
%% Skull Strip MRI
%%
fprintf('Performing Skull Extraction\n');
bse=fullfile(BrainSuitePath,'bin','bse');
bseout=fullfile(anatDir,sprintf('%s_T1w.ds.orig.bse.nii.gz',subid));
cmd=sprintf('%s --auto -i %s -o %s',bse,t1ds,bseout);
unix(cmd);
%% Coregister t1 to MNI Space
%%
bsenew=fullfile(anatDir,sprintf('%s_T1w.nii.gz',subid));
cmd=sprintf('flirt -ref %s -in %s -out %s',BSA,bseout,bsenew);
unix(cmd);
bsemask=fullfile(anatDir,'t1.mask.nii.gz');
cmd=sprintf('fslmaths %s -thr 0 -bin -mul 255 %s -odt char',bsenew,bsemask);
unix(cmd);
bsenew2=fullfile(anatDir,sprintf('%s_T1w.bse.nii.gz',subid));
copyfile(bsenew,bsenew2);
%% Run BrainSuite and SVReg
%%
subbasename=fullfile(anatDir,sprintf('%s_T1w',subid));
cmd=sprintf('%s %s',bst_exe,subbasename);
unix(cmd)
cmd=sprintf('%s %s %s',svreg_exe,subbasename,BCIbasename);
unix(cmd);
%% Run Batch_Process Pipeline for fMRI
% 
%% Transfer data to surface and then to USCBrain atlas surface and produce surface grayordinates
%%
fprintf('Transferring data from subject to atlas\n');
for ind = 1:length(fmri)
    outdir=fullfile(subdir,'func');
    fmri2surfFile=fullfile(outdir,'fmri2surf.mat');
    GOrdSurfFile=fullfile(outdir,'fmri2surf_GOrd.mat');
    resample2surf(subbasename,fullfile(outdir,'fmri2standard.nii.gz'),fmri2surfFile);
    generateSurfGOrdfMRI(GOrdSurfIndFile,fmri2surfFile,GOrdSurfFile);
end
%% Transfer data to volumetric grayordinates
%%
fprintf('Generating Volumetric Grayordinates\n');
for ind = 1:length(fmri)
    outdir=fullfile(subdir,sprintf('func-%d',ind));
    GOrdVolFile=fullfile(outdir,'fmri2Vol_GOrd.mat');
    GOrdFile=fullfile(outdir,'fmri.32k.nii.gz');
    generateVolGOrdfMRI(GOrdVolIndFile,subbasename,fullfile(outdir,'fmri2standard.nii.gz'),GOrdVolFile);
    combineSurfVolGOrdfMRI(GOrdSurfFile,GOrdVolFile,GOrdFile);
end
%%
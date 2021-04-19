# BFP: BrainSuite fMRI Pipeline
The BrainSuite Functional Pipeline (BFP) processes 4D fMRI datasets using a combination of tools from AFNI, FSL, BrainSuite and additional in-house tools develeoped for BrainSuite. BFP additionally provides novel workflows for statistical analysis using [BrainSync](https://www.sciencedirect.com/science/article/pii/S1053811918300582), a tool that temporally aligns spatially registered fMRI datasets for direct timesries comparions between subjects.


BFP requires T1w-images that have been processed through BrainSuite in BIDS format where all subject prefixes are formatted as *subID_T1w* and all outputs are in an *anat* subdirectory. Output Directory structure should be as follows:
```
+ Main Directory
  + subID
    + anat
      + subID_T1w.nii.gz
      + subID_T1w.bfc.nii.gz
      + subID_T1w ...
```
If T1w-images have not been previously processed or are not in BIDSapp format, BFP will process them through BrainSuite's structural pipeline before processing fMRI data.  

BFP will then process fMRI images through a series modules. Final outputs are put in GrayOrdinate and BrainOrdinate space. Workflows are customized using a configuration file as described in the usage section.
![](https://github.com/ajoshiusc/bfp/blob/master/docs/BFP_workflow_ver4p01_t1distcorr.png)

## System Requirements
Current version of BFP supports Linux and MacOS. 
Linux: A Debian based linux distribution (ubuntu/mint/debian), 8GB RAM recommended. 
Windows: May run through a virtual machine, linux subsystem or docker. We have not done testing on such configurations, so we do not support it at the moment, contributions are welcome.

## Installation and Setup
 * Clone or Download BFP
	git clone https://github.com/ajoshiusc/bfp.git
	<required: git https://git-scm.com/>
 * Install [BrainSuite](http://forums.brainsuite.org/download/). We recommend version 18a or later.
 * Install AFNI (Ver. Jun 12 2017 or newer) and FSL ( Ver. 6.0 or newer). We recommend using [NeuroDebian](http://neuro.debian.net) as it makes the installation process easier. BFP has been tested with this configuration, but it should work for other versions of Linux and other softwares. However we recommend and support the above mentioned versions.

## Usage
### Matlab
bfp(configfile,t1,fmri,studydir,subid,sessionid,TR)
+ configfile: Configuration file indicates processing options. (<bfpdir>/src/preproc/config.ini) 
+ t1: T1 image in NIFTI-1 format e.g. 'subID_T1w.nii.gz'
+ fmri: string or cell array of fMRI NIFTI-1 formatted files e.g. 'subID_rest_bold.nii.gz'
+ studydir, subid: Outputs will be saved in the studydir/subid for each subject e.g. studydir='mystudy', subid='subID'
+ sessionid: string or cell array of session id for the input fmris e.g. 'rest' for resting state
 
 e.g.
bfp /home/ajoshi/bfp/supp_data/config.ini /home/ajoshi/sub-01_T1w.nii.gz /home/ajoshi/sub-01_rest.nii.gz /home/ajoshi/mystudy sub-01 rest 2


### Compiling
You need matlab compiler to create the binary from bfp matlab scripts.
In matlab, change directory to src/preproc and run
compile_bfp version number
where version number can be a string like ver2p23
compile_bfp ver2p23



### Binary
You can run BFP from the bash shell using the command mentioned below.

bfp.sh /home/ajoshi/bfp/supp_data/config.ini /home/ajoshi/sub-01_T1w.nii.gz /home/ajoshi/sub-01_rest.nii.gz /home/ajoshi/mystudy sub-01 rest 2


## Statistical analysis pipelines
The statistical analysis pipelines are available in src/stats directory as python scripts. All documentation is included within the sample config and python scripts.

required inputs:
* config file: input parameters and options for your statistical test. 
Sample config file is provided in (<bfpdir>/src/stats/sample_config_stats.ini) 
* demographics file: csv file with your study's demographics which will be used for identifying subjects and their respective demographic information. You can control for up to 2 variables through linear regression (such as age and sex). 
Sample csv file is provided in (<bfpdir>/src/stats/sample_demo_linear_regr.csv)

bfp_run_stat.py is the main statistical pipeline. All parameters are set using a config file as input found (<bfpdir>/src/stats/sample_config_stats.ini)

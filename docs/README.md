# BFP: BrainSuite fMRI Pipeline
 This pipeline takes fMRI and anatomical data and processes them using a series 
 of scripts from [AFNI](https://afni.nimh.nih.gov/) , [BrainSuite](http://brainsuite.org/) and [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki). The functional processing script is 
 based on `batch_process.sh` script from [fcon1000](http://fcon_1000.projects.nitrc.org/).
 
## System Requirements
Current version of BFP supports Linux and MacOS. 
Linux: A Debian based linux distribution (ubuntu/mint/debian), 8GB RAM recommended. 
Windows: May run through a virtual machine, linux subsystem or docker. We have not done testing on such configurations, so we do not support it at the moment, contributions are welcome.

## Installation and Setup
 * Clone or Download BFP
	git clone https://github.com/ajoshiusc/bfp.git
	<required: git https://git-scm.com/>
 * Install [BrainSuite18a](http://forums.brainsuite.org/download/) or later version.
 * Install AFNI (Ver. Jun 12 2017 or newer) and FSL ( Ver. 5.0 or newer). We recommend using [NeuroDebian](http://neuro.debian.net) as it makes the installation process easier. BFP has been tested with this configuration, but it should work for other versions of Linux and other softwares. However we recommend and support the above mentioned versions.

 * Open <bfpdir>/src/preproc/sample_main.m in Matlab.
 * Open <bfpdir>/src/preproc/bfp_preproc.ini and set paths for FSL, AFNI and BrainSuite.
 * For preprocessing fMRI, run src/preproc/sample_main.m from matlab. For binary, please follow instructions below.

## Usage
### Matlab
bfp(configfile,t1,fmri,studydir,subid,sessionid,TR)

 configfile: Configuration file that you edited during installation, (<bfpdir>/supp_data/config.ini) 
 
 t1: T1 image in NIFTI-1 format e.g. 'sub01_bT1w.nii.gz'

 fmri: string or cell array of fMRI NIFTI-1 formatted files e.g. 'sub_01_rest_bold.nii.gz'

 studydir, subid: Outputs will be saved in the studydir/subid for each subject e.g. studydir='mystudy', subid='sub-01'

 sessionid: string or cell array of session id for the input fmris e.g. 'rest' for resting state
 
 
 e.g.
bfp /home/ajoshi/bfp/supp_data/config.ini /home/ajoshi/sub-01_T1w.nii.gz /home/ajoshi/sub-01_rest.nii.gz /home/ajoshi/mystudy sub-01 rest 2


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

bfp_group_compare.py is the group difference pipeline.
bfp_linear_regr.py is the regression pipeline.
bfp_linear_regr_pairwise.py is the atlas-free regression pipeline.

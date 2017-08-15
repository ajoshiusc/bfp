# BFP: BrainSuite fMRI Pipeline
 This pipeline takes fMRI and anatomical data and processes them using a series 
 of scripts from [AFNI](https://afni.nimh.nih.gov/) , [BrainSuite](http://brainsuite.org/) and [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki). The functional processing script is 
 based on `batch_process.sh` script from [fcon1000](http://fcon_1000.projects.nitrc.org/).
 
## System Requirements
Current version of BFP is Linux only. A Debian based linux distribution (ubuntu/mint/debian), 8GB RAM recommended. 

## Installation and Setup
 * Install AFNI (Ver. Jun 12 2017) and FSL 5.0. We recommend using [NeuroDebian](http://neuro.debian.net) as it makes the installation process easier. 
 * Install [BrainSuite 17a](http://brainsuite.org).
 * Install Matlab 2017a if you want to run it from Matlab
 * Download bfp, unzip
 * Open bfp.mtx in Matlab Live Notebook and edit *User Inputs* section.
 * open <bfpdir>/supp_data/config.ini and set paths for FSL, AFNI and BrainSuite.
 * For preprocessing fMRI, run src/preproc/bfp.mtx

## Usage
### Matlab
bfp(configfile,t1,fmri,studydir,subid,sessionid,TR)

 configfile: Configuration file that you edited during installation, (<bfpdir>/supp_data/config.ini) 
 
 t1: T1 image in NIFTI-1 format e.g. 'sub01-bold.nii.gz'

 fmri: string or cell array of fMRI NIFTI-1 formatted files e.g. 'sub_01_rest_bold.nii.gz'

 studydir, subid: Outputs will be saved in the studydir/subid for each subject e.g. studydir='mystudy', subid='sub-01'

 sessionid: string or cell array of session id for the input fmris e.g. 'rest' for resting state
 
 
 e.g.
bfp /home/ajoshi/bfp/supp_data/config.ini /home/ajoshi/sub-01_T1w.nii.gz /home/ajoshi/sub-01_rest.nii.gz /home/ajoshi/mystudy sub-01 rest 2


### Binary
You can run BFP from the bash shell using the command mentioned below.

bfp.sh /home/ajoshi/bfp/supp_data/config.ini /home/ajoshi/sub-01_T1w.nii.gz /home/ajoshi/sub-01_rest.nii.gz /home/ajoshi/mystudy sub-01 rest 2

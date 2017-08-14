# BFP: BrainSuite fMRI Pipeline
 This pipeline takes fMRI and anatomical data and processes them using a series 
 of scripts from AFNI(https://afni.nimh.nih.gov/) , [BrainSuite](http://brainsuite.org/) and [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki). The functional processing script is 
 based on `batch_process.sh` script from [fcon1000](http://fcon_1000.projects.nitrc.org/).
 
## System Requirements
Current version of BFP is Linux only. A Debian based linux distribution (ubuntu/mint/debian), 8GB RAM recommended. 

## Installation and Setup
 * Install AFNI (Ver. Jun 12 2017) and FSL 5.0. We recommend using neurodebian (http://neuro.debian.net) as it makes the installation process easier. 
 * Install BrainSuite 17a (http://brainsuite.org).
 * Install Matlab 2017a
 * Download bfp, unzip
 * Open bfp.mtx in Matlab Live Notebook and edit *User Inputs* section.
 * open config.ini and set paths for FSL, AFNI and BrainSuite.
 * For preprocessing fMRI, run src/preproc/bfp.mtx


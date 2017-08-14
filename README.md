# BFP: BrainSuite fMRI Pipeline
 This pipeline takes fMRI and anatomical data and processes them using a series 
 of scripts from BrainSuite, FSL and AFNI. The functional processing script is 
 based on batch_process.sh script from fcon1000.
## System Requirements
 * Debian based linux distribution, 16GB RAM. 
## Installation and Setup
 * Install AFNI (Ver. Jun 12 2017), FSL 5.0 using neurodebian (http://neuro.debian.net). Install BrainSuite 17a (http://brainsuite.org).
 * Install Matlab 2017a
 * Download bfp, unzip
 * Open bfp.mtx in Matlab Live Notebook and edit *User Inputs* section
 * open config.ini and set paths make sure that the paths of the programs are 
 correct.
 * For preprocessing fMRI, run src/preproc/bfp.mtx


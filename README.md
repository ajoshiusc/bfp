# BFP: BrainSuite fMRI Pipeline
 This pipeline takes fMRI and anatomical data and processes them using a series 
 of scripts from BrainSuite, FSL and AFNI. The functional processing script is 
 based on batch_process. sh script from fcon1000.
 
## Issues/TBD
 * Global Signal Regression is on for now. We need to figure out what signals 
 to regress. The current pipeline regresses global, WM and GM signals.
 * Surface corresponding to 32k need to be finalized and distributed with the 
 BFP distribution.
 * Visualization of Grayordinate data: How to do this?
 * Current pipeline includes volumetric smoothing with 6mm FWHM. This corresponds 
 to 3mm isotropic voxels in fmri data. Is this required, considering that we 
 do tNLMPDF filtering at the end.
 * tNLMPDF filtering is done at the end. Is this appropriate place for this?
 * BrainSync is not included in this pipeline yet. We have to decide how to 
 include it since this pipeline is for single subject.
 * Most scans don't cover full cerebellum, so there are NaNs in the grayordinates. 
 We need to make sure that this is accounted in the further processing.
 * There are several outputs generated. Should we keep them all or delete them 
 once done.
 * Should we support partial runs or continuation of full runs of the pipeline?
 * Documentation needs to be done
## Installation and Setup
 * Install AFNI (Ver. Jun 12 2017), FSL 5.0 using neurodebian (<http://neuro.debian.net/ 
 http://neuro.debian.net/> ). Install BrainSuite 17a (<http://brainsuite.org 
 http://brainsuite.org>).
 * Install Matlab 2017a
 * Download bfp, unzip
 * Open bfp.mtx in Matlab Live Notebook and edit *User Inputs* section
 * open config.ini and set paths make sure that the paths of the programs are 
 correct.
 * Run bfp.mtx

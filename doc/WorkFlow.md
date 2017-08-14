# BFP Workflow

## Creating Directory Structure
BFP stores data in [BIDS](http://bids.neuroimaging.io/) format. The inputs are T1 and fMRI images. They are stored in a subject subdirectory inside the study directory. T1 images are stored in anat subdirectory and fmri images in func subdirectory. There can be only one T1 image but multiple fMRI images are supported. 

## Anatomical Processing
The analtomical processing pipeline in BFP is based on [BrainSuite](http://brainsuite.org). 
* First, we resample the T1 image to 1mm cubic resolution. This is because BrainSuite tools perform optimally at that resolution. 
* Next we perform skull extraction of the T1 image. This is done using [bse](http://brainsuite.org/processing/surfaceextraction/bse/) executable in BrainSuite.
* The extracted image is then coregistered using [FLIRT](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT) to the BrainSuiteAtlas1 atlas that is included with BrainSuite. This step is performed for both fMRI and T1 images so that both are in the common space. 
* The [CSE sequence](http://brainsuite.org/processing/surfaceextraction/) in BrainSuite is then executed on the coregistered image (skull stripping is skipped). This performs tissue classification and generate inner, mid and pial cortical surface representations.
* The brain registration and labeling is performed using [SVReg](http://brainsuite.org/processing/svreg/).

## Functional Processing
The functional data is processed using several tools from [AFNI](https://afni.nimh.nih.gov/) and [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) as detailed below. 
* We generate 3mm isotropic representation of BrainSuiteAtlas1 as a standard reference. fMRI data is coregistered to this template. This steps puts subjects T1 and fMRI data in a common space.
* The fMRI data is deobliqued. This step changes some of the information inside a fmri dataset's header and is reoriented in FSL friendly space. 
* The motion correction is performed to average of time series. This is followed by skull stripping.
* Spatial smoothing is performed with Gaussian kernel and FWHM of 6mm.
* This is followed by [grand mean scaling](http://dbic.dartmouth.edu/wiki/index.php/Global_Scaling) to account for signal differences across sessions.
* A temporal band pass filtering is then applied with bandwidth of 0.005-0.1 Hz.

## fMRI analysis tools



# BFP Workflow

## Anatomical Processing
The analtomical processing pipeline in BFP is based on [BrainSuite](http://brainsuite.org). 
* First, we perform skull extraction of the T1 image. This is done using [bse](http://brainsuite.org/processing/surfaceextraction/bse/) executable in BrainSuite.
* The extracted image is then coregistered using FLIRT(https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT) to the BrainSuiteAtlas1 atlas that is included with BrainSuite. This step is performed for both fMRI and T1 images so that both are in the common space. 

## Functional Processing

## fMRI analysis tools



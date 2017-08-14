# BFP Workflow

## Anatomical Processing
The analtomical processing pipeline in BFP is based on [BrainSuite](http://brainsuite.org). 
* First, we perform skull extraction of the T1 image. This is done using [bse](http://brainsuite.org/processing/surfaceextraction/bse/) executable in BrainSuite.
* The extracted image is then coregistered to the BrainSuiteAtlas1. This step is performed for both fMRI and T1 images so that both are in the common space. 
First the T1 image is coregistered to the BrainSuiteAtlas1 

## Functional Processing

## fMRI analysis tools



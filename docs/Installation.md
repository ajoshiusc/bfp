# BFP: Installation

 BFP pipeline takes fMRI and anatomical data and processes them using a series 
 of scripts from [AFNI](https://afni.nimh.nih.gov/) , [BrainSuite](http://brainsuite.org/) and [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki). The functional processing script is 
 based on `batch_process.sh` script from [fcon1000](http://fcon_1000.projects.nitrc.org/).
 
## System Requirements
Current version of BFP is Linux only. It may work on Mac, but is not supported on that platform. A Debian based linux distribution (ubuntu/mint/debian), 8GB RAM recommended. 

## Installation and Setup
 * Install AFNI (Ver. Jun 12 2017 or greater) and FSL 5.0 or greater. We recommend using [NeuroDebian](http://neuro.debian.net) as it makes the installation process easier. 
 * Install [BrainSuite](http://brainsuite.org).
 * Download bfp, unzip
 * open config.ini and set paths for FSL, AFNI and BrainSuite.


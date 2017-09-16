#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=11:59:59
#PBS -l mem=32gb
#PBS -l vmem=32gb
#PBS -m a
echo ${subid}

/usr/usc/matlab/default/bin/matlab -nodisplay -nosplash -r "addpath(genpath('/home/rcf-proj2/aaj/git_sandbox/bfp/src')); bfp /home/rcf-proj2/aaj/git_sandbox/bfp/supp_data/hpcconfig.ini /home/rcf-proj2/aaj/Beijing_Zang/${subid}/anat/mprage_anonymized.nii.gz {'/home/rcf-proj2/aaj/Beijing_Zang/${subid}/func/rest.nii.gz'} /home/rcf-proj2/aaj/git_sandbox/bfp/data/ ${subid} {'rest'} 2"

exit 0


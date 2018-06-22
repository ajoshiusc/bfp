#!/bin/bash
#SBATCH --mem=16000
#SBATCH --time=1:00:00
#SBATCH --ntasks=4
echo ${subid}

/usr/usc/matlab/default/bin/matlab -nodisplay -nosplash -r "addpath(genpath('/home/rcf-proj2/aaj/git_sandbox/bfp/src')); bfp /home/rcf-proj2/aaj/git_sandbox/bfp/supp_data/hpcconfig.ini /home/rcf-proj2/aaj/ADHD_Peking/data/Peking_all/${subid}/session_1/anat_1/mprage.nii.gz {'/home/rcf-proj2/aaj/ADHD_Peking/data/Peking_all/${subid}/session_1/rest_1/rest.nii.gz'} /home/rcf-proj2/aaj/ADHD_Peking_bfp/ ${subid} {'rest'} 2"

exit 0


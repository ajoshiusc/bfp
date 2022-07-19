clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/home/ajoshi/projects/bfp/src'));

%studydir='/ImagePTE1/ajoshi/ADHD_Peking_bfp';
studydir = '/home/ajoshi/tmp_study_dir';
sessionid='rest';

a=dir(studydir);

for j=1:length(a)
    subid = a(j).name;
    gen_brainordinates_alff('/home/ajoshi/BrainSuite19b', studydir, subid, sessionid);
end

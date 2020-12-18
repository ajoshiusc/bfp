clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/ImagePTE1/ajoshi/code_farm/bfp/src'));

%studydir='/ImagePTE1/ajoshi/ADHD_Peking_bfp';
studydir = '/ImagePTE1/ajoshi/maryland_rao_v1_bfp';
sessionid='rest';

a=dir(studydir);

for j=1:length(a)
    subid = a(j).name;
    gen_brainordinates('/home/ajoshi/BrainSuite19b', studydir, subid, sessionid);
end

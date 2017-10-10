function compile_bfp()

clc;clear all;close all;
restoredefaultpath;
addpath(genpath('../../src'));
mcc -m -v bfp.m
copyfile('../scripts/bfp_linux.sh','./bfp.sh');






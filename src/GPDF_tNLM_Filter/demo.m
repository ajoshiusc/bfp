clear; close all; clc;

T = 200;
data1 = randn(60000, T);
data2 = repmat(randn(1, T), 20000, 1);
data = [data2; data1; data2];
data = data + 3 * randn(100000, T);

option = tNLMGPDF();
option.isPlot = false;
option.isVerbose = true;
option.memoryLimit = 'auto';
% option.SCBFile = fullfile(pwd, 'SCB.mat');

%%
[dataSm, output] = tNLMGPDF(data, option);

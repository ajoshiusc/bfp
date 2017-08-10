clear; close all; clc;

T = 200;
data1 = randn(600, T);
data2 = repmat(randn(1, T), 200, 1);
data = [data2; data1; data2];
data = data + 3 * randn(1000, T);

option = tNLMPdf();
option.isPlot = true;
option.isVerbose = true;
option.memoryLimit = 'auto';
option.SCBFile = fullfile(pwd, 'SCB.mat');

%%
[dataSm, output] = tNLMPdf(data, option);

%
% [dataOut, meanVal, stdVal] = zeroMeanUnitVariance(dataIn, dim)
% 
% Description:
%     normalize data by removing mean and make it unit variance
% 
% Input:
%     dataIn - data input
%     dim - on which dimension to normalize data
% 
% Output:
%     dataOut - data output
%     meanVal - mean value that subtracted
%     stdVal - standard deviation that divided by
% 
% Copyright:
%     2013-2018 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     1.0.3
% Date:
%     2018/05/02
%

function [dataOut, meanVal, stdVal] = zeroMeanUnitVariance(dataIn, dim)
    [dataZM, meanVal] = zeroMean(dataIn, dim);
    [dataOut, stdVal] = unitVariance(dataZM, dim);
end
%
% [dataOut, meanVal] = zeroMean(dataIn, dim)
% 
% Description:
%     remove mean from data on a particular dimension
% 
% Input:
%     dataIn - data input
%     dim - on which dimension to remove mean
% 
% Output:
%     dataOut - data output
%     meanVal - mean value that subtracted
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

function [dataOut, meanVal] = zeroMean(dataIn, dim)

    meanVal = mean(dataIn, dim);
    
    if dim == 1
        dataOut = dataIn - meanVal(ones(size(dataIn, dim), 1), :, :, :);
    elseif dim == 2
        dataOut = dataIn - meanVal(:, ones(size(dataIn, dim), 1), :, :);
    elseif dim == 3
        dataOut = dataIn - meanVal(:, :, ones(size(dataIn, dim), 1), :);
    elseif dim == 4
        dataOut = dataIn - meanVal(:, :, :, ones(size(dataIn, dim), 1));
    end
    
end
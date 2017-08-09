%
% [dataOut, stdVal] = unitVariance(dataIn, dim)
% 
% Description:
%     normalize data to unit variance on a particular dimension
% 
% Input:
%     dataIn - data input
%     dim - on which dimension to unit variance
% 
% Output:
%     dataOut - data output
%     stdVal - the standard deviation that divided by
% 
% Copyright:
%     2013-2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     1.0.2
% Date:
%     2017/08/04
%

function [dataOut, stdVal] = unitVariance(dataIn, dim)

    stdVal = std(dataIn, 0, dim);
    
    if dim == 1
        dataOut = dataIn ./ stdVal(ones(size(dataIn, dim), 1), :, :, :);
    elseif dim == 2
        dataOut = dataIn ./ stdVal(:, ones(size(dataIn, dim), 1), :, :);
    elseif dim == 3
        dataOut = dataIn ./ stdVal(:, :, ones(size(dataIn, dim), 1), :);
    elseif dim == 4
        dataOut = dataIn ./ stdVal(:, :, :, ones(size(dataIn, dim), 1));
    end
    
end
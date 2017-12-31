%
% dataOut = unitNorm(dataIn, dim)
% 
% Description:
%     normalize the data to unit norm on a particular dimension
% 
% Input:
%     dataIn - data input
%     dim - on which dimension to unit norm, only support 1 or 2 for now
% 
% Output:
%     dataOut - data output
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

function dataOut = unitNorm(dataIn, dim)
    
    if dim == 1
        flag = false;
    elseif dim == 2
        flag = true;
    else
        error('dim can only be 1 or 2');
    end
    
    if flag
        X = dataIn';
    else
        X = dataIn;
    end
    
    for m = 1:size(X, 2)
        X(:, m) = X(:, m) ./ norm(X(:, m));
    end
    
    if flag
        dataOut = X';
    else
        dataOut = X;
    end
end
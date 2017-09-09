%
% numCores = findNumberOfCores()
% 
% Description:
%     find number of available cpu cores
% 
% Input:
% 
% Output:
%     numCores - number of cpu cores
% 
% Copyright:
%     2014-2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     1.0.2
% Date:
%     2017/08/16
%

function numCores = findNumberOfCores()
    numCores = feature('numCores');
end
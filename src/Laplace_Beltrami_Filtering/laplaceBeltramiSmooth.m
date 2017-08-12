%
% dataSm = laplaceBeltramiSmooth(LBO, data, sigma)
% 
% Description:
%     Laplace Beltrami filtering
% 
% Input:
%     LBO - Laplace Beltrami operator structure
%     data - data to be smoothed
%     sigma - Gaussian sigma
% 
% Output:
%     dataSm - filtered data
% 
% Copyright:
%     2016-2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Anand Joshi, Jian (Andrew) Li
% Revision:
%     1.0.3
% Date:
%     2017/06/01
%

function dataSm = laplaceBeltramiSmooth(LBO, data, sigma)
    
    idxNaN = any(isnan(data), 2);
    
    data2 = data;
    data2(idxNaN, :) = 0;
    
    E = LBO.E;
    L = LBO.L;
    kernel = exp(-L * sigma);

    dataSm = E * bsxfun(@times, kernel, E' * data2);
    dataSm(idxNaN, :) = nan;
    
end
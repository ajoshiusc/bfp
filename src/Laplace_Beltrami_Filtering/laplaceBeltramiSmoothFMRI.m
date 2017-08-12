%
% [dataSmL, dataSmR] = laplaceBeltramiSmoothFMRI(dataL, dataR, LBOL, LBOR, sigmaLB)
% 
% Description:
%     Laplace Beltrami filtering for fMRI data
% 
% Input:
%     dataL - data to be smoothed on left hemisphere
%     dataR - data to be smoothed on right hemisphere
%     LBOL - Laplace Beltrami operator structure on left hemisphere
%     LBOR - Laplace Beltrami operator structure on right hemisphere
%     sigmaLB - Gaussian sigma
% 
% Output:
%     dataSmL - filtered data on left hemisphere
%     dataSmR - filtered data on right hemisphere
% 
% Copyright:
%     2016-2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Anand Joshi, Jian (Andrew) Li
% Revision:
%     1.0.5
% Date:
%     2017/06/01
%

function [dataSmL, dataSmR] = laplaceBeltramiSmoothFMRI(dataL, dataR, LBOL, LBOR, sigmaLB)
    
    idxNaNL = any(~isfinite(dataL), 2);
    dataT1L = dataL;
    dataT1L(idxNaNL, :) = 0;
    dataT2L = laplaceBeltramiSmooth(LBOL, dataT1L, sigmaLB);
    dataT3L = zeroMean(dataT2L(~idxNaNL, :), 2);
    dataSmL = zeros(size(dataL));
    dataSmL(idxNaNL, :) = nan; 
    dataSmL(~idxNaNL, :) = dataT3L;
    
    idxNaNR = any(~isfinite(dataR), 2);
    dataT1R = dataR;
    dataT1R(idxNaNR, :) = 0;
    dataT2R = laplaceBeltramiSmooth(LBOR, dataT1R, sigmaLB);
    dataT3R = zeroMean(dataT2R(~idxNaNR, :), 2);
    dataSmR = zeros(size(dataR));
    dataSmR(idxNaNR, :) = nan; 
    dataSmR(~idxNaNR, :) = dataT3R;
    
end
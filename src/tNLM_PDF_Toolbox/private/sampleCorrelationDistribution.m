%
% P = sampleCorrelationDistribution(r, rho, N, isPlot)
% 
% Description:
%     calcualte the sample correlation distribution (Fisher 1915)
% 
% Input:
%     r - sample correlation
%     rho - population correlation
%     N - length of the time series
%     isPlot - flag to plot
% 
% Output:
%     P - the distribution evaluted at r
% 
% Copyright:
%     2016-2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     1.0.2
% Date:
%     2017/06/01
%

function P = sampleCorrelationDistribution(r, rho, N, isPlot)

    if ~exist('isPlot', 'var')
        isPlot = false;
    end
    
    a1 = (N-2) / sqrt(2*pi);
    a2 = gammaln(N-1);
    a3 = (N-1)/2 * log(1-rho^2);
    a4 = (N-4)/2 * log(1-r.^2);
    a5 = gammaln(N-0.5);
    a6 = (N-3/2) * log(1-rho.*r);
    
    a7 = hypergeom([0.5 0.5], (2*N-1)/2, (rho.*r+1)/2);

    P = a1 .* exp(a2 + a3 + a4 - a5 - a6) .* a7;
    P = P ./ sum(P);
    
    if isPlot
        figure, plot(r, P);
        grid on;
    end
    
end
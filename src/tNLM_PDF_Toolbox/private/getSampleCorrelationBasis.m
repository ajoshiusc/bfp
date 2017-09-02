%
% [basis, r, rho] = getSampleCorrelationBasis(T, fileName)
% 
% Description:
%     obtain sample correlation basis
% 
% Input:
%     T - number of time samples
%     fileName - file name for caching
% 
% Output:
%     basis - sample correlation basis
%     r - sample correlation sampling point
%     rho - population correlation sampling point
% 
% Copyright:
%     2016-2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     2.1.3
% Date:
%     2017/07/06
%

function [basis, r, rho] = getSampleCorrelationBasis(T, fileName)

    if ~exist('fileName', 'var') || isempty(fileName)
        fileName = fullfile(pwd, 'SCB.mat');
    end
    
    if exist(fileName, 'file')
        s = load(fileName);
        SCB = s.SCB;
        Ts = SCB.Ts;
        Bs = SCB.Basis;
        r = SCB.R;
        rho = SCB.Rho;
        
        idx = find(Ts == T);
        
        if isempty(idx)
            [basis, r, rho] = calculateBasis(T);
            num = length(Ts);
            Ts = [Ts, T];
            Bs{num+1} = basis;
            
            SCB.R = r;
            SCB.Rho = rho;
            SCB.Ts = Ts;
            SCB.Basis = Bs;
            
            save(fileName, 'SCB');
        else
            basis = Bs{idx};
        end
    else
        [basis, r, rho] = calculateBasis(T);
        SCB = struct();
        SCB.R = r;
        SCB.Rho = rho;
        SCB.Ts = T;
        SCB.Basis = {basis};
        save(fileName, 'SCB');
    end
end

function [basis, r, rhos] = calculateBasis(T)
    
    disp(['you do not have sample correlation basis for length ' num2str(T)]);
    disp(['calculating the basis and it needs to be done only once for length ' num2str(T)]);
    
    r = -1:0.001:1;
    rhos = -1:0.01:1;
    numRhos = length(rhos);

    basis = zeros(length(r), numRhos);
    
    openHere = false;
    if ~parpoolOperator('isopen')
        openHere = true;
        parpoolOperator('open');
    end
    
    parProgressTracker('s', numRhos-2);

    parfor m = 2:numRhos-1
        rho = rhos(m);
        P = sampleCorrelationDistribution(r, rho, T);
        basis(:, m) = P(:);
        parProgressTracker('p');
    end
    
    parProgressTracker('e');
    
    if openHere
        parpoolOperator('close');
    end
    
    % when rho is -1
    firstOne = zeros(length(r), 1);
    firstOne(1) = 1;
    basis(:, 1) = firstOne;
    
    % when rho = 1
    lastOne = zeros(length(r), 1);
    lastOne(end) = 1;
    basis(:, end) = lastOne;

end

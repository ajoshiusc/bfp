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
%     2016-2018 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     2.3.0
% Date:
%     2018/05/02
%

function [basis, r, rho] = getSampleCorrelationBasis(T, fileName, isVerbose)

    if ~exist('fileName', 'var') || isempty(fileName)
        fileName = fullfile(pwd, 'SCB.mat');
    end
    
    if ~exist('isVerbose', 'var') || isempty(isVerbose)
        isVerbose = true;
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
            [basis, r, rho] = calculateBasis(T, isVerbose);
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
        [basis, r, rho] = calculateBasis(T, isVerbose);
        SCB = struct();
        SCB.R = r;
        SCB.Rho = rho;
        SCB.Ts = T;
        SCB.Basis = {basis};
        save(fileName, 'SCB');
    end
end

function [basis, r, rhos] = calculateBasis(T, isVerbose)
    
    if isVerbose
        disp(['you do not have sample correlation basis for length ' num2str(T)]);
        disp(['calculating the basis and it needs to be done only once for length ' num2str(T)]);
    end
    
    r = -1:0.001:1;
    rhos = -1:0.01:1;
    numRhos = length(rhos);

    basis = zeros(length(r), numRhos);
    
    isParTB = checkMatlabToolbox('pct'); % if parallel processing toolbox installed
    
    if isParTB
        openHere = false;
        if ~parpoolOperator('isopen')
            openHere = true;
            parpoolOperator('open', 8);
        end
        
        if isVerbose, parProgressTracker('s', numRhos-2); end

        parfor m = 2:numRhos-1
            rho = rhos(m);
            P = sampleCorrelationDistribution(r, rho, T);
            basis(:, m) = P(:);
            if isVerbose, parProgressTracker('p'); end
        end

        if isVerbose, parProgressTracker('e'); end

        if openHere
            parpoolOperator('close');
        end
    else
        if isVerbose
            strLen = progressTracker(0, numRhos-2, 0, 50);
        end
        
        for m = 2:numRhos-1
            rho = rhos(m);
            P = sampleCorrelationDistribution(r, rho, T);
            basis(:, m) = P(:);
            
            if isVerbose
                strLen = progressTracker(m-1, numRhos-2, strLen, 50);
            end
        end
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

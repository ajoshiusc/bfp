%
% [dataSm, output] = tNLMPdf(data, option)
% 
% Description:
%     PDF-based temporal non-local means filtering
% 
% Input:
%     data - time series N x T
%     option - option structure
% 
% Output:
%     dataSm - filtered time series
%     output - output structure containing intermediate results
% 
% Copyright:
%     2016-2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     9.2.2
% Date:
%     2017/08/03
%

function [dataSm, output] = tNLMPdf(data, option)
    
    if nargin == 0 % return default option
        option = struct;
        option.selfWeightMode = 1;
        option.normalization = 1;
        option.FPR = 1e-3;
        option.memoryLimit = 'auto';
        option.isPlot = false;
        option.isVerbose = false;
        option.SCBFile = fullfile(pwd, 'SCB.mat');
        dataSm = option;
        return;
    end

    if ~exist('option', 'var') || isempty(option)
        option = tNLMPdf();
    end
    
    if option.isVerbose
        disp('estimating memory requirement');
    end
    
    [numV, numT] = size(data);
    
    % minimum and recommended number of samples for kernel estimation
    minNumSpKE = min(numV, 1e4);
    recomNumSpKE = min(numV, 2.3e4);
    
    numSpKE = recomNumSpKE;
    
    memReq = estimateVarSize('double', numV*numT*4 + minNumSpKE^2, 'GB');
    memRec = estimateVarSize('double', numV*numT*4 + recomNumSpKE^2, 'GB');
    memFullCorrMat = estimateVarSize('double', numV*numV, 'GB');
    
    % only work for linux now
    memInfo = getMemoryInfo();
    if ~isempty(memInfo)
        isDetected = true;
        memLim = memInfo.available;
    else
        isDetected = false;
    end
    
    if isscalar(option.memoryLimit) && (option.memoryLimit > 0)
        if isDetected && (memLim < option.memoryLimit)
            disp('you have less available memory than the specified limit, try to use as much as possible');
        else
            memLim = option.memoryLimit;
        end
    elseif ischar(option.memoryLimit) && strcmp(option.memoryLimit, 'auto')
        if ~isDetected
            disp('could not detect memory information automatically');
            disp('specify option.memoryLimit manually in GB and then try again');
            dataSm = []; output = [];
            return;
        end
    else
        error('set memory limit to a postive number in GB or auto');
    end
    
    if memLim < memReq
        disp('you do not have enough memory');
        disp(['the minimum required memory for this dataset is ' num2str(memReq) ' GB']);
        disp('try to free some cache or downsample your data, then run this function again');
        disp('process ceased');
        dataSm = []; output = [];
        return;
    elseif (memLim >= memReq) && (memLim < memRec)
        numSpKE = minNumSpKE;
        disp(['we recommend you have at least ' num2str(memRec) ' GB memory']);
        disp('this recommendation is not satisfied, will try to proceed in memory-saving mode');
        disp('however, in this mode, the estimation of the kernel may not be accurate and the filtering process may be very slow');
    end
    
    if option.isVerbose
        disp('pre-processing data');
    end
    
    dataN = zeroMeanUnitVariance(data, 2);
    idxNaN = any(~isfinite(dataN), 2);
    data2 = data(~idxNaN, :);
    dataN2 = dataN(~idxNaN, :);
    clear data dataN;
    
    if option.isVerbose
        disp('estimating the kernel');
        disp('downsample dataset and calculate the correlation matrix');
    end
    
    if numV > numSpKE
        dataN3 = dataN2(randsample(numV, numSpKE), :);
    else
        dataN3 = dataN2;
    end
    
    A = corr(dataN3'); clear dataN3;
    A2 = A(:); clear A;
    A2(A2 == 1) = [];
    
    if option.isVerbose
        disp('fit the conditional distribution');
    end
    
    % get basis <=> conditional distribution P(r|rho)
    [basis, r, rhos] = getSampleCorrelationBasis(numT, option.SCBFile);
    
    % boundary of H0
    thH0 = 0.011;
    idxH0 = (rhos >= -thH0) & (rhos <= thH0);
    
    binRes = 0.001;
    binEdge = (-1-binRes/2):binRes:(1+binRes/2);
    ct = histcounts(A2, binEdge); clear A2;
    ct2 = ct(:) ./ sum(ct); % empirical distribution
    
    AtA = basis' * basis;
    Atb = basis' * ct2;
    coef = fastNNLS(AtA, Atb); % prior distribution P(rho)
    
    PJoint = bsxfun(@times, basis, coef'); % joint dist P(r, rho)
    Pr = sum(PJoint, 2); % fitted P(r)
    PPost = bsxfun(@rdivide, PJoint, Pr); % posterior P(rho|r)
    
    res = abs(ct2 - Pr);
    err = res' * res;
    
    if option.isVerbose
        disp(['fitting error = ' num2str(err)]);
        disp('optimize the parameter based on the FPR');
    end
    
    % use 1/P(rho_0|r)-1 as the ratio
    PPostH0 = PPost(:, idxH0);
    PPostH0(isnan(PPostH0)) = 0;
    PPostH0 = sum(PPostH0, 2);
    t = 1./PPostH0 - 1;

    p0 = sum(PJoint(:, idxH0), 2);
    p0 = p0 ./ sum(p0);
    
    % binary search find optimal h
    idxNeg = (r < 0);
    hLeft = 0.01; hRight = 1e10;
    h = (hLeft + hRight) / 2;
    e = 1e10; th = 1e-3;
    while e > th
        w = 1 - exp(-t/(h^2));
        w(idxNeg) = 0; 
        fpr = sum(w .* p0);
        if fpr < option.FPR
            hRight = h;
        elseif fpr > option.FPR
            hLeft = h;
        elseif fpr == option.FPR
            break;
        end
        h = (hLeft + hRight) / 2;
        e = abs(fpr - option.FPR) / option.FPR;
    end
    
    if option.isVerbose
        disp(['best h = ' num2str(h)]);
    end
    
    if option.isPlot
        figure;
        ax1 = axes();
        maxVal = max(max(ct2), max(p0));
        plot(ax1, r, ct2/maxVal, 'k');
        hold on; grid on;
        plot(ax1, r, p0/maxVal);
        plot(ax1, r, w, 'r');
        ax1.YTick = 0:0.1:1;
        ax1.XTick = -1:0.2:1;
        ax1.YLim = [-0.1 1.1];
        xlabel('r');
        legend(ax1, 'Sample Correlation Histogram', ...
                    'Noise distribution', ...
                    'Filtering weights', ...
                    'Location', 'NorthWest');
    end
    
    clear PJoint PPost basis AtA;
    
    if option.isVerbose
        disp('tNLM-pdf filtering');
    end
    
    % subtract memory used for storing dataset
    memLim = memLim - estimateVarSize('double', numV*numT*4, 'GB');
    
    % use 95% avaiable memory to avoid freezing
    memLim = memLim * 0.95;
    
    % multiply by 2 because A and B simutaneously exist below for a short
    % time
    numBlk = ceil(memFullCorrMat / memLim * 2);
    if option.isVerbose
        disp(['based on the memory limitation, we will have ' num2str(numBlk) ' iterations']);
    end
    
    szBlk = ceil(numV / numBlk);
    
    data3 = zeros(size(data2));
    
    strLen = progressTracker(0, numBlk, 0, 50);
    
    for m = 1:numBlk
        idxS = (m-1)*szBlk+1;
        if m < numBlk
            idxE = m*szBlk;
        elseif m == numBlk
            idxE = size(dataN2, 1);
        end
        
        dataBlkN = dataN2(idxS:idxE, :);
        A = corr(dataBlkN', dataN2');
        B = interp1(r, w, A);
        
        szThisBlk = size(dataBlkN, 1);
        idxDiag = sub2ind(size(A), 1:szThisBlk, (1:szThisBlk)+(m-1)*szBlk);
        clear dataBlkN A;
        
        if option.selfWeightMode == 0
            % do nothing
        elseif option.selfWeightMode == 1
            B2 = B;
            B2(idxDiag) = 0;
            maxB2 = max(B2, [], 2);
            B(idxDiag) = maxB2; clear B2;
        elseif option.selfWeightMode == 2
            B(idxDiag) = 0;
        else
            error('self weight mode wrong');
        end
        
        d = sum(B, 2);
        B = spdiags(1./d, 0, size(B, 1), size(B, 1)) * B;
        data3(idxS:idxE, :) = B * data2;
        clear B d;
        
        strLen = progressTracker(m, numBlk, strLen, 50);
    end
    
    if option.isVerbose
        disp('normalize data as specified');
    end
    
    if option.normalization ==0
        % do nothing
    elseif option.normalization == 1
        data3 = zeroMean(data3, 2);
    elseif option.normalization == 2
        data3 = zeroMeanUnitVariance(data3, 2);
    else
        error('wrong normalization option');
    end
    
    if option.isVerbose
        disp('generate output');
    end
    
    dataSm = zeros(numV, numT);
    dataSm(~idxNaN, :) = data3;
    dataSm(idxNaN, :) = nan;
    
    output = struct();
    output.r = r;
    output.h = h;
    output.w = w;
    
    if option.isVerbose
        disp('done');
    end
end

%
% [dataSm, output] = tNLMGPDF(data, option)
% 
% Description:
%     Global PDF-based temporal non-local means filtering
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
%     2016-2019 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian Li (Andrew)
% Revision:
%     9.5.2
% Date:
%     2019/11/12
%

function [dataSm, output] = tNLMGPDF(data, option)
    
    if nargin == 0 % return default option
        option = struct;
        option.selfWeightMode = 1;
        option.normalization = 1;
        option.FPR = 1e-3;
        option.memoryLimit = 'auto';
        option.numCPU = 'auto';
        option.isPlot = false;
        option.isVerbose = false;
        option.SCBFile = '/home/andrew/Developer/Matlab/icMatBox/signal_processing/SCB.mat';
%         option.SCBFile = fullfile(pwd, 'SCB.mat');
        dataSm = option;
        return;
    end

    if ~exist('option', 'var') || isempty(option)
        option = tNLMGPDF();
    end
    
    if option.isVerbose
        disp('pre-processing data');
    end
    
    if ~(strcmp(option.numCPU, 'auto') || (option.numCPU >= 1))
        error('option numCPU needs to be either auto or a postive integer');
    end
    
    % make sure data is double type
    data = double(data);
    
    [numV, numT] = size(data);
    
    dataN = zeroMeanUnitVariance(data, 2);
    idxNaN = any(~isfinite(dataN), 2);
    data2 = data(~idxNaN, :);
    dataN2 = dataN(~idxNaN, :);
    clear data dataN;
    
    numV2 = size(data2, 1);
    
    if option.isVerbose
        disp('estimating memory requirement');
    end
    
    % minimum and recommended number of samples for kernel estimation
    minNumSpKE = min(numV2, 1e4);
    recomNumSpKE = min(numV2, 2.3e4);
    
    numSpKE = recomNumSpKE;
    
    memReq = estimateVarSize('double', minNumSpKE^2 * 3, 'GB');
    memRec = estimateVarSize('double', recomNumSpKE^2 * 3, 'GB');
    memFullCorrMat = estimateVarSize('double', numV2*numV2, 'GB');
    
    memInfo = getMemoryInfo();
    if ~isempty(memInfo.available)
        isDetected = true;
        memLim = memInfo.available;
    else
        isDetected = false;
    end
    
    if isscalar(option.memoryLimit) && (option.memoryLimit > 0)
        if isDetected && (memLim < option.memoryLimit)
            disp('you have less available memory than the specified limit');
        else
            memLim = option.memoryLimit;
        end
    elseif ischar(option.memoryLimit) && strcmp(option.memoryLimit, 'auto')
        if ~isDetected
            str = {'Can not detect memory information automatically.', ...
                   'You need to specify option.memoryLimit manually in GB and try again.', ...
                   sprintf('Or you may proceed with assumption that you have at least %s GB available.', num2str(memReq)), ...
                   'Note that this may be significantly slower than usual or even causing crash if the assumption is not satisfied', ...
                   'Would you like to proceed anyway? (Default is No)'};
            choice = questdlg(str, 'Memory Information Not Available', 'Yes', 'No', 'No');
            switch choice
                case 'Yes'
                    memLim = memReq * 1.0001;
                case 'No'
                    dataSm = []; output = [];
                    return;
            end
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
        disp('however, in this mode, the estimation of the kernel may not be accurate and the filtering process may be slower than usual');
    end
    
    if option.isVerbose
        disp('estimating the kernel');
        disp('downsample dataset and calculate the correlation matrix');
    end
    
    if numV2 > numSpKE
        dataN3 = dataN2(randsample(numV2, numSpKE), :);
    else
        dataN3 = dataN2;
    end
    
    A = corr(dataN3'); clear dataN3;
    A = A(:);
    A = A(A < 1 - 1e-6);
    
    if option.isVerbose
        disp('fit the conditional distribution');
    end
    
    % get basis <=> conditional distribution P(r|rho)
    [basis, r, rhos] = getSampleCorrelationBasis(numT, option.SCBFile, ...
                                        option.numCPU, option.isVerbose);
    
    % get boundary of H0 based on theoretical h0 distribution 50% prob
    p00 = basis(:, 101);
    t = cumsum(p00);
    thH0 = abs(r(find(t > 0.25, 1, 'first')));
    idxH0 = (rhos >= -thH0) & (rhos <= thH0);
    
    binRes = 0.001;
    binEdge = (-1-binRes/2):binRes:(1+binRes/2);
    ct = histcounts(A, binEdge); clear A;
    ct2 = ct(:) ./ sum(ct); % empirical distribution
    
    AtA = basis' * basis;
    Atb = basis' * ct2;
    coef = fastNNLS(AtA, Atb); % prior distribution P(rho)
    
    if sum(coef(idxH0)) < 1e-3
        disp('could not estimate noise distribution correctly, no filtering performed');
        dataSm = []; output = [];
        return;
    end
    
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
        legend(ax1, 'Sample Correlation Histogram', 'Noise distribution', ...
                    'Filtering weights', 'Location', 'NorthWest');
    end
    
    clear PJoint PPost basis AtA;
    
    if option.isVerbose
        disp('GPDF filtering');
    end
    
    % use 90% avaiable memory to avoid freezing
    memLim = memLim * 0.9;
    
    % multiply by 4 because A, B and other data simutaneously exist below
    % especially interp1 will double the memory usage of A or B
    % for a short time
    numBlk = ceil(memFullCorrMat / memLim * 4);
    if option.isVerbose
        str = ['based on the memory limitation, we will have ' num2str(numBlk)];
        if numBlk == 1
            str = [str ' iteration'];
        elseif numBlk > 1
            str = [str ' iterations'];
        end
        disp(str);
    end
    
    szBlk = ceil(numV2 / numBlk);
    
    data3 = zeros(size(data2));
    
    if option.isVerbose
        strLen = progressTracker(0, numBlk, 0, 50);
    end
    
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
        clear d;
        
        if option.isVerbose
            strLen = progressTracker(m, numBlk, strLen, 50);
        end
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
    output.B = B;
    output.PPrior = coef;
    
    if option.isVerbose
        disp('done');
    end
end

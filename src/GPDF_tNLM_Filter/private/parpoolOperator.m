%
% res = parpoolOperator(mode, numReqCores, isVerbose)
% 
% Description:
%     a helper function to operate parpool
% 
% Input:
%     mode - open, forceopen, reopendiff, close or isopen
%     numReqCores - number of requested cores
%     isVerbose - whether display verbose message
% 
% Output:
%     res - for isopen mode only, return binary result of whether there is an existing parpool
% 
% Copyright:
%     2015-2018 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     1.5.2
% Date:
%     2018/05/02
%

function res = parpoolOperator(mode, numReqCores, isVerbose)

    if ~exist('isVerbose', 'var') || isempty(isVerbose)
        isVerbose = false;
    end

    numSysCores = findNumberOfCores();

    if ~exist('numReqCores', 'var') || isempty(numReqCores)
        numReqCores = numSysCores - 1;
    else
        if numReqCores > numSysCores
            if isVerbose
                disp('You do not have that many cores as requested, set to the maximum number of cores in the system - 1');
            end
            numReqCores = numSysCores - 1;
        end
    end
    
    res = false;
    
    % get pool object
    poolobj = gcp('nocreate');
    
    % save one core for the system, 2 threads per core
    poolSize = numReqCores * 2;
    
    c = parcluster();
    prop = properties(c);
    if ismember('NumWorkers', prop)
        c.NumWorkers = numReqCores * 2;
    end
    
    if ismember('NumThreads', prop)
        c.NumThreads = 1;
    end
    
    if strcmpi(mode, 'open')
        if isempty(poolobj)
            parpool(c, poolSize);
        else
            if isVerbose
                disp('a parallel pool exists');
            end
        end
    elseif strcmpi(mode, 'forceopen')
        if ~isempty(poolobj)
            delete(poolobj);
            parpool(c, poolSize);
        end
    elseif strcmpi(mode, 'reopendiff')
        if ~isempty(poolobj)
           	if poolSize ~= poolobj.NumWorkers
                delete(poolobj);
                parpool(c, poolSize);
            else
                if isVerbose
                    disp('a parallel pool exists with same size');
                end
            end
        else
            parpool(c, poolSize);
        end
    elseif strcmpi(mode, 'close')
        if isempty(poolobj)
            if isVerbose
                disp('a parallel pool does not exist');
            end
        else
            delete(poolobj);
        end
    elseif strcmpi(mode, 'isopen')
        if isempty(poolobj)
            res = false;
        else
            res = true;
        end
    end
    
end
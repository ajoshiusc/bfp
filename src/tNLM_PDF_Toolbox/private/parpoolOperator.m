%
% res = parpoolOperator(mode)
% 
% Description:
%     a helper function to operate parpool
% 
% Input:
%     mode - open, forceopen, close or isopen
% 
% Output:
%     res - for isopen mode only, return binary result of whether there is an existing parpool
% 
% Copyright:
%     2015-2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     1.3.1
% Date:
%     2017/08/17
%

function res = parpoolOperator(mode)

    res = false;
    
    poolobj = gcp('nocreate');
    
    numCores = findNumberOfCores();
    poolSize = numCores * 2 - 2;
    
    c = parcluster();
    c.NumWorkers = numCores * 2;
    c.NumThreads = 1;

    if strcmpi(mode, 'open')
        
        if isempty(poolobj)
            parpool(c, poolSize);
        else
            warning('a parallel pool exists');
        end
        
    elseif strcmpi(mode, 'forceopen')
        
        if ~isempty(poolobj)
            delete(poolobj);
        end
        
        parpool(c, poolSize);
        
    elseif strcmpi(mode, 'close')
        
        if isempty(poolobj)
            warning('a parallel pool does not exist');
        else
            delete(poolobj);
        end
       
    elseif strcmpi(mode, 'isopen')
        res = ~isempty(poolobj);
    end
    
end
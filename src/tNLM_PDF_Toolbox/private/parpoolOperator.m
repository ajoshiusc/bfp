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
%     1.2.5
% Date:
%     2017/08/16
%

function res = parpoolOperator(mode)

    res = false;
    
    poolobj = gcp('nocreate');
    
    numCores = findNumberOfCores();
    poolSize = numCores * 2 - 2;
    
    if strcmpi(mode, 'open')
        
        if isempty(poolobj)
            parpool(poolSize);
        else
            warning('a parallel pool exists');
        end
        
    elseif strcmpi(mode, 'forceopen')
        
        if ~isempty(poolobj)
            delete(poolobj);
        end
        
        parpool(poolSize);
        
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
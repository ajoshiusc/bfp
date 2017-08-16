%
% memInfo = getMemoryInfo(unit)
% 
% Description:
%     get system memory information
% 
% Input:
%     unit - unit of Byte
% 
% Output:
%     memInfo - memory information structure
% 
% Copyright:
%     2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     1.1.3
% Date:
%     2017/08/16
%

function memInfo = getMemoryInfo(unit)
    
    if ~exist('unit', 'var') || isempty(unit)
        unit = 'GB';
    end
    
    units = getByteUnits();
    if ~ismember(unit, units)
        disp(units);
        error('unit has to be one of above:');
    end
    
    memInfo = struct;
    memInfo.total = [];
    memInfo.used = [];
    memInfo.free = [];
    memInfo.shared = [];
    memInfo.buffers = [];
    memInfo.cache = [];
    memInfo.available = [];

	[status, res] = unix('free -bw');
    
    if status ~= 0
        [status, res] = unix('free -b');
    end
    
    if status == 0
        mem = regexp(res, 'Mem:', 'split');
        ctg = strtrim(mem{1});
        ctg = strsplit(ctg);
        mem = strtrim(mem{2});
        mem = regexp(mem, 'Swap:', 'split');
        mem = strtrim(mem{1});
        mem = strsplit(mem);
        mem = cellfun(@str2double, mem, 'UniformOutput', false);
        
        for m = 1:length(ctg)
            if strncmpi(ctg{m}, 'total', 4)
                memInfo.total = convertByteUnit(mem{m}, 'B', unit);
            elseif strncmpi(ctg{m}, 'used', 4)
                memInfo.used = convertByteUnit(mem{m}, 'B', unit);
            elseif strncmpi(ctg{m}, 'free', 4)
                memInfo.free = convertByteUnit(mem{m}, 'B', unit);
            elseif strncmpi(ctg{m}, 'shared', 4)
                memInfo.shared = convertByteUnit(mem{m}, 'B', unit);
            elseif strncmpi(ctg{m}, 'buffers', 4)
                memInfo.buffers = convertByteUnit(mem{m}, 'B', unit);
            elseif strncmpi(ctg{m}, 'cache', 4)
                memInfo.cache = convertByteUnit(mem{m}, 'B', unit);
            elseif strncmpi(ctg{m}, 'available', 4)
                memInfo.available = convertByteUnit(mem{m}, 'B', unit);
            end
        end
        
        if isempty(memInfo.available)
            t = 0;
            if ~isempty(memInfo.cache)
                t = t + memInfo.cache;
            end
            
            if ~isempty(memInfo.free)
                t = t + memInfo.free;
            end
            
            if t > 0
                memInfo.available = t;
            end
        end
    end
end

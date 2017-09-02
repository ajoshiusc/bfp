%
% varargout = estimateVarSize(varType, numE, unit)
% 
% Description:
%     estimate variable size
% 
% Input:
%     varType - type of variable, e.g. double
%     numE - number of elements
%     unit - unit of Byte
% 
% Output:
%     varargout - print to screen if no output, otherwise output the estimated size and unit
% 
% Copyright:
%     2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     1.0.3
% Date:
%     2017/08/16
%

function varargout = estimateVarSize(varType, numE, unit)

    if ~exist('unit', 'var') || isempty(unit)
        unit = 'auto';
    end

    map = mapType2Size();
    
    if ~ismember(varType, map.keys)
        disp('primitive variable type has to be one of the following:');
        disp(map.keys);
        return;
    end
    
    if numE < 0
        disp('num of elememts must not be negative');
        return;
    end
    
    Bt = numE * map(varType);
    [estSize, unitOut] = convertByteUnit(Bt, 'B', unit);
    
    if nargout == 0
        fprintf('%g %s\n', estSize, unitOut);
        varargout = {};
    elseif nargout == 1
        varargout{1} = estSize;
    elseif nargout == 2
        varargout{1} = estSize;
        varargout{2} = unitOut;
    end
end
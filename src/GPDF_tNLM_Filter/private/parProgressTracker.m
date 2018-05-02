%
% parProgressTracker(mode, num)
% 
% Description:
%     Track progress for parallel loop
% 
% Input:
%     mode - 's' - start, 'p' - progress or 'e' - end
%     num - total number of iterations, only needed for mode 's'
% 
% Output:
% 
% Copyright:
%     2016-2018 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     2.1.8
% Date:
%     2018/05/02
%

function parProgressTracker(mode, num)

    if strcmpi(mode, 's')
        fprintf('Total number of iterations:\n');
        fprintf(['[' repmat('-', 1, num) ']\n']);
        fprintf('Progress:\n[ \n');
    elseif strcmpi(mode, 'p')
        fprintf('\b\b=>\n');
    elseif strcmpi(mode, 'e')
        fprintf('\b\b]\n');
        fprintf('Done\n');
    end
    
end
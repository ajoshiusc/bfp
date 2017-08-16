%
% strLen = progressTracker(counter, totalNum, prevStrLen, barLen)
% 
% Description:
%     tracking a progress of a for loop by printing a progress bar and the percentage
% 
% Input:
%     counter - current iteration counter
%     totalNum - total number of iterations
%     prevStrLen - previous output of PROGRESSTRACKER
%     barLen - the length of progress bar wanted to show on screen
% 
% Output:
%     strLen - the current length of string shown on screen used for the next input
% 
% Copyright:
%     2013-2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     1.1.2
% Date:
%     2017/08/16
%

function strLen = progressTracker(counter, totalNum, prevStrLen, barLen)

    % if the counter is greater than the total number of iteratoin,
    % terminate and pop error
    if counter > totalNum
        error('current counter can not be greater than the total iterations');
    end
    
    % generate string for progress bar
    pctg = counter/totalNum;
    numProc = floor(pctg * barLen);
    numLeft = barLen - numProc;
    strProgressBar = ['[', repmat('=', 1, numProc), repmat(' ', 1, numLeft), ']'];
    
    if numProc == barLen
        strProgressBar(numProc + 1) = '>';
    else
        strProgressBar(numProc + 2) = '>';
    end
    
    % generate string for percentages
    strNumExp = sprintf('%d/%d', counter, totalNum);
    strPctgExp = sprintf('%.1f%%', pctg*100);
    
    % change the arrow to = in last position of progress bar if all done
    if pctg == 1
        strProgressBar(barLen + 1) = '=';
    end
    
    % concatenate into output string
    output = [strProgressBar, '   ', strNumExp, '   ', strPctgExp];
    
    % prevent backspace in the very beginning
    if prevStrLen ~= 0
        prevStrLen = prevStrLen + 1;
    end
    
    % delete previous output
    fprintf(1, repmat('\b', 1, prevStrLen));
    
    % print current output
    disp(output);
    
    % calculate the length of current output which will be used for next
    % iteration
    strLen = numel(output);
    
    % if all done, print 'Done' and turn to a new line
    if pctg == 1
        fprintf('\b  Done.\n');
    end
    
end
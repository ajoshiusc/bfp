%RESTOREDEFAULTPATH Restores the MATLAB search path to installed products.
%   RESTOREDEFAULTPATH restores the MATLAB search path to its factory-
%   installed state. RESTOREDEFAULTPATH is only intended for situations
%   when the search path is corrupted and MATLAB is experiencing problems
%   during startup.
%
%   RESTOREDEFAULTPATH; MATLABRC sets the search path to include only
%   folders for installed products from MathWorks and corrects search path
%   problems encountered during startup.
%
%   MATLAB does not support issuing RESTOREDEFAULTPATH from a UNC path
%   name. Doing so might result in MATLAB being unable to find files on the
%   search path. If you do issue RESTOREDEFAULTPATH from a UNC path name,
%   restore the expected behavior by changing the current folder to an
%   absolute path, and then reissuing RESTOREDEFAULTPATH. 
%
%   See also ADDPATH, GENPATH, MATLABRC, RMPATH, SAVEPATH.

%   Copyright 2003-2011 The MathWorks, Inc.

% Get system path to Perl (MATLAB installs Perl on Windows)
if strncmp(computer,'PC',2)
    RESTOREDEFAULTPATH_perlPath = [matlabroot '\sys\perl\win32\bin\perl.exe'];
    RESTOREDEFAULTPATH_perlPathExists = exist(RESTOREDEFAULTPATH_perlPath,'file')==2;
else
    [RESTOREDEFAULTPATH_status, RESTOREDEFAULTPATH_perlPath] = unix('which perl');
    RESTOREDEFAULTPATH_perlPathExists = RESTOREDEFAULTPATH_status==0;
    RESTOREDEFAULTPATH_perlPath = (regexprep(RESTOREDEFAULTPATH_perlPath,{'^\s*','\s*$'},'')); % deblank lead and trail
end

% If Perl exists, execute "getphlpaths.pl"
if RESTOREDEFAULTPATH_perlPathExists
    RESTOREDEFAULTPATH_cmdString = sprintf('"%s" "%s" "%s"', ...
        RESTOREDEFAULTPATH_perlPath, which('getphlpaths.pl'), matlabroot);
    [RESTOREDEFAULTPATH_perlStat, RESTOREDEFAULTPATH_result] = system(RESTOREDEFAULTPATH_cmdString);
else
    error(message('MATLAB:restoredefaultpath:PerlNotFound'));
end

% Check for errors in shell command
if (RESTOREDEFAULTPATH_perlStat ~= 0)
    error(message('MATLAB:restoredefaultpath:PerlError',RESTOREDEFAULTPATH_result,RESTOREDEFAULTPATH_cmdString));
end

% Check that we aren't about to set the MATLAB path to an empty string
if isempty(RESTOREDEFAULTPATH_result)
    error(message('MATLAB:restoredefaultpath:EmptyPath'))
end

% Set the path, adding userpath if possible
if exist( 'userpath.m', 'file' ) == 2
    matlabpath([userpath ';' RESTOREDEFAULTPATH_result]);
else
    matlabpath(RESTOREDEFAULTPATH_result);
end

clear('RESTOREDEFAULTPATH_*');

% Create this variable so that if MATLABRC is run again, it won't try to
% use pathdef.m
RESTOREDEFAULTPATH_EXECUTED = true;

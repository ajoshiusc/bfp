%
% isInstalled = checkMatlabToolbox(tb)
% 
% Description:
%     check if a particular toolbox is installed
% 
% Input:
%     tb - toolbox short name, only support pct, smlt, ipt now
% 
% Output:
%     isInstalled - whether is installed or not
% 
% Copyright:
%     2013-2018 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     0.0.3
% Date:
%     2018/01/22
%

function isInstalled = checkMatlabToolbox(tb)

    switch tb
        case 'pct'
            v = ver('distcomp');
        case 'smlt'
            v = ver('stats');
        case 'ipt'
            v = ver('images');
        otherwise
            v = [];
    end
    
    if isempty(v)
        isInstalled = false;
    else
        isInstalled = true;
    end
    
end

% use this command to find the right word to check in ver
% ls(toolboxdir(''))

% bioinfo  dig  fixedpoint  images  optim  signal  slde  stateflow
% coder  distcomp  fixpoint  instrument  physmod  simulink  sldv  stats
% control  dsp  glee  local  realtime  simulinktest  sl_pir_cap  symbolic
% curvefit  eml  hdlcoder  matlab  rtw  sl3d  sltp  target
% iagram  finance  idelink  multisim  shared  slcontrol  slvnv
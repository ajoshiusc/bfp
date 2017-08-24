function varargout = gbfp(varargin)
% GBFP MATLAB code for gbfp.fig
%      GBFP, by itself, creates a new GBFP or raises the existing
%      singleton*.
%
%      H = GBFP returns the handle to a new GBFP or the handle to
%      the existing singleton*.
%
%      GBFP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GBFP.M with the given input arguments.
%
%      GBFP('Property','Value',...) creates a new GBFP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gbfp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
% 
% <<FILENAME.PNG>>
% 
%      stop.  All inputs are passed to gbfp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gbfp

% Last Modified by GUIDE v2.5 23-Aug-2017 17:18:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gbfp_OpeningFcn, ...
                   'gui_OutputFcn',  @gbfp_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gbfp is made visible.
function gbfp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gbfp (see VARARGIN)

% Choose default command line output for gbfp
handles.output = hObject;
fslPath = "/usr/share/fsl/5.0/bin"; 
set(handles.fslDir, 'string', fslPath);

afniPath = "/usr/lib/afni/bin"; 
set(handles.afniDir, 'string', afniPath);

BrainSuitePath = "~/BrainSuite/BrainSuite17a"; 
set(handles.BrainSuiteDir, 'string', BrainSuitePath);

current = pwd;
cd('../../')
bfpRoot = sprintf('%s',pwd);
cd(current);

set(handles.bfpDir,'string',bfpRoot);

scbdir = bfpRoot;
set(handles.scbDir,'string',scbdir);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gbfp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gbfp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fslDir_Callback(hObject, eventdata, handles)
% hObject    handle to fslDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fslDir as text
%        str2double(get(hObject,'String')) returns contents of fslDir as a double


% --- Executes during object creation, after setting all properties.
function fslDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fslDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function afniDir_Callback(hObject, eventdata, handles)
% hObject    handle to afniDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of afniDir as text
%        str2double(get(hObject,'String')) returns contents of afniDir as a double


% --- Executes during object creation, after setting all properties.
function afniDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to afniDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BrainSuiteDir_Callback(hObject, eventdata, handles)
% hObject    handle to BrainSuiteDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BrainSuiteDir as text
%        str2double(get(hObject,'String')) returns contents of BrainSuiteDir as a double


% --- Executes during object creation, after setting all properties.
function BrainSuiteDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BrainSuiteDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fslBrowse.
function fslBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to fslBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir(get(handles.fslDir,'String'));
textLabel = sprintf('%s', folder);
set(handles.fslDir, 'string', textLabel);




% --- Executes on button press in afniBrowse.
function afniBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to afniBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% User to select afni directory
folder = uigetdir(get(handles.afniDir,'String'));
% load folder directory in a string
textLabel = sprintf('%s', folder);
% display user selected directory in GUI
set(handles.afniDir, 'string', textLabel);




% --- Executes on button press in BrainSuiteBrowse.
function BrainSuiteBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to BrainSuiteBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir(get(handles.BrainSuiteDir,'String'));
textLabel = sprintf('%s', folder);
set(handles.BrainSuiteDir, 'string', textLabel); 


% --- Executes on button press in T1Browse.
function T1Browse_Callback(hObject, eventdata, handles)
% hObject    handle to T1Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path,~] = uigetfile({'*.nii.gz','Compressed NIFTI files (*.nii.gz)'}, ...
    'Choose Structural Image(s)','MultiSelect','on');
textLabel = sprintf('%s%s', path ,file);
set(handles.structural, 'string', textLabel);  



% --- Executes on button press in fmriBrowse.
function fmriBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to fmriBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path,~] = uigetfile({'*.nii.gz','Compressed NIFTI files (*.nii.gz)'}, ...
    'Choose Functional Image(s) ','MultiSelect','on');
textLabel = sprintf('%s%s', path, file);
set(handles.functional, 'string', textLabel); 

% --- Executes on button press in subjectBrowse.
function subjectBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to subjectBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function structural_Callback(hObject, eventdata, handles)
% hObject    handle to structural (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of structural as text
%        str2double(get(hObject,'String')) returns contents of structural as a double


% --- Executes during object creation, after setting all properties.
function structural_CreateFcn(hObject, eventdata, handles)
% hObject    handle to structural (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function functional_Callback(hObject, eventdata, handles)
% hObject    handle to functional (see GCBO)
% eventdata  reserved - to btextLabele defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of functional as text
%        str2double(get(hObject,'String')) returns contents of functional as a double


% --- Executes during object creation, after setting all properties.
function functional_CreateFcn(hObject, eventdata, handles)
% hObject    handle to functional (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SubjectID_Callback(hObject, eventdata, handles)
% hObject    handle to SubjectID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SubjectID as text
%        str2double(get(hObject,'String')) returns contents of SubjectID as a double


% --- Executes during object creation, after setting all properties.
function SubjectID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SubjectID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loadData(handles)


function loadData(handles) 
% Choose directory for configuration file
[config,configPath,~] = uigetfile('*.ini','Choose Configuration file');
addpath(genpath(uigetdir('~','Choose BFP Root Directory')));
% List all params and input them to the GUI
params = ini2struct([configPath config]);

set(handles.fslDir, 'string', params.FSLPATH)
set(handles.afniDir, 'string', params.AFNIPATH)
set(handles.BrainSuiteDir, 'string', params.BrainSuitePath)
set(handles.structural, 'string', params.T1)
set(handles.functional, 'string', params.fmri)
set(handles.SubjectID, 'string', params.subid)
set(handles.SessionID, 'string', params.sessionid)
set(handles.HighPass, 'string', params.HIGHPASS)
set(handles.LowPass, 'string', params.LOWPASS)
set(handles.fwhm, 'string', params.FWHM)
set(handles.TR, 'string', params.TR)
set(handles.Continue, 'value', str2double(params.CONTINUERUN))
set(handles.bfpDir, 'string',params.BFPPATH)
set(handles.studyDir,'string',params.studydir)
set(handles.config,'string',params.configName)
set(handles.fpr,'string',params.fpr)
set(handles.memory,'string',params.memory)
set(handles.scbDir,'string',params.scbPath)

% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Needs to execute the entie program. Include path to bfp suite
MyConfig(hObject,handles);


% --- Executes on button press in stage.
function stage_Callback(hObject, eventdata, handles)
% hObject    handle to stage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Continue.
function Continue_Callback(hObject, eventdata, handles)
% hObject    handle to Continue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Continue
status = get(hObject,'Value');
if status == 1 
    set(hObject,'Value',str2num('1'))
else
    set(hObject,'Value',str2num('0'))
end


% --- Executes on button press in FSLoutput.
function FSLoutput_Callback(hObject, eventdata, handles)
% hObject    handle to FSLoutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FSLoutput


% --- Executes on button press in LowPass.
function LowPass_Callback(hObject, eventdata, handles)
% hObject    handle to LowPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LowPass


% --- Executes on button press in HighPass.
function HighPass_Callback(hObject, eventdata, handles)
% hObject    handle to HighPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HighPass



function fwhm_Callback(hObject, eventdata, handles)
% hObject    handle to fwhm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fwhm as text
%        str2double(get(hObject,'String')) returns contents of fwhm as a double


% --- Executes during object creation, after setting all properties.
function fwhm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fwhm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
string


function SessionID_Callback(hObject, eventdata, handles)
% hObject    handle to SessionID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SessionID as text
%        str2double(get(hObject,'String')) returns contents of SessionID as a double


% --- Executes during object creation, after setting all properties.
function SessionID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SessionID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TR_Callback(hObject, eventdata, handles)
% hObject    handle to TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TR as text
%        str2double(get(hObject,'String')) returns contents of TR as a double


% --- Executes during object creation, after setting all properties.
function TR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MyConfig(hObject, handles)

% If your LD library path is in another directory change it here.
lib = 'LD_LIBRARY_PATH=/usr/lib/fsl/5.0';
fsloutput = 'FSLOUTPUTTYPE=NIFTI_GZ';
fslPath = sprintf('FSLPATH=%s', get(handles.fslDir,'string'));
afniPath = sprintf('AFNIPATH=%s', get(handles.afniDir,'string'));
brainSuitePath = sprintf('BrainSuitePath=%s', get(handles.BrainSuiteDir,'string'));
fwhm = sprintf('FWHM=%s', get(handles.fwhm,'string'));
high = sprintf('HIGHPASS=%s', get(handles.HighPass,'string'));
low  = sprintf('LOWPASS=%s', get(handles.LowPass,'string'));
proceed = sprintf('CONTINUERUN=%s', num2str(get(handles.Continue,'Value')));
bfpPath = sprintf('BFPPATH=%s', get(handles.bfpDir,'string'));
addpath(genpath(get(handles.bfpDir,'string')))
t1 = sprintf('T1=%s', get(handles.structural,'string'));
fmri{1} = sprintf('fmri=%s', get(handles.functional,'string'));
subid = sprintf('subid=%s', char(get(handles.SubjectID, 'string')));
TR = sprintf('TR=%s', char(get(handles.TR, 'string')));
sessionid{1} = sprintf('sessionid=%s',char(get(handles.SessionID,'string')));
memory = sprintf('memory=%s',char(get(handles.memory,'string')));
fpr = sprintf('fpr=%s',get(handles.fpr,'string'));
scbPath = sprintf('scbPath=%s',get(handles.scbDir,'string'));
% write to config.ini
configName = get(handles.config,'string');
configname = sprintf('configName=%s',configName);
studydir = get(handles.studyDir,'string');
studyDir = sprintf('studydir=%s', studydir);
keys = {'','','',char(bfpPath);'','','',char(afniPath);'','','',char(fslPath); ...
        '','','',char(brainSuitePath);'','','',char(studyDir);'','','',char(lib); ...
        '','','',char(fsloutput);'','','',char(fwhm);'','','',char(high);'','', '',char(low);...
        '','','',char(proceed);'','','',char(t1);'','','',char(fmri);'','','',char(subid); ...
        '','','',char(TR); '','','',char(sessionid);'','','',char(memory);'','','',char(fpr); ...
        '','','',char(configname);'','','',char(scbPath)};
PWD = pwd;
cd(studydir);
inifile(configName,'write',keys,'plain');
configfile = fullfile(studydir,configName);
cd(PWD);


bfp(configfile,char(get(handles.structural,'string')),cellstr(get(handles.functional,'string'))...
    ,studydir,char(get(handles.SubjectID, 'string')),cellstr(get(handles.SessionID,'string'))...
    ,char(get(handles.TR, 'string')))



function config_Callback(hObject, eventdata, handles)
% hObject    handle to config (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of config as text
%        str2double(get(hObject,'String')) returns contents of config as a double

% --- Executes during object creation, after setting all properties.
function config_CreateFcn(hObject, eventdata, handles)
% hObject    handle to config (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function studyDir_Callback(hObject, eventdata, handles)
% hObject    handle to studyDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of studyDir as text
%        str2double(get(hObject,'String')) returns contents of studyDir as a double


% --- Executes during object creation, after setting all properties.
function studyDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to studyDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in StudyBrowse.
function StudyBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to StudyBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir('~','Choose Output Study Directory');
textLabel = sprintf('%s', folder);
set(handles.studyDir, 'string', textLabel); 

% --- Executes on selection change in memory.
function memory_Callback(hObject, eventdata, handles)
% hObject    handle to memory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns memory contents as cell array
%        contents{get(hObject,'Value')} returns selected item from memory


% --- Executes during object creation, after setting all properties.
function memory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to memory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fpr_Callback(hObject, eventdata, handles)
% hObject    handle to fpr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fpr as text
%        str2double(get(hObject,'String')) returns contents of fpr as a double


% --- Executes during object creation, after setting all properties.
function fpr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fpr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bfpDir_Callback(hObject, eventdata, handles)
% hObject    handle to bfpDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bfpDir as text
%        str2double(get(hObject,'String')) returns contents of bfpDir as a double


% --- Executes during object creation, after setting all properties.
function bfpDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bfpDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bfp.
function bfp_Callback(hObject, eventdata, handles)
% hObject    handle to bfp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir('~','Choose BFP root Directory');
textLabel = sprintf('%s', folder);
set(handles.bfpDir, 'string', textLabel); 
addpath(genpath(textLabel));



function scbDir_Callback(hObject, eventdata, handles)
% hObject    handle to scbDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scbDir as text
%        str2double(get(hObject,'String')) returns contents of scbDir as a double


% --- Executes during object creation, after setting all properties.
function scbDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scbDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in scb.
function scb_Callback(hObject, eventdata, handles)
% hObject    handle to scb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir('~','Choose SCB root Directory');
textLabel = sprintf('%s', folder);
set(handles.scbDir, 'string', textLabel); 

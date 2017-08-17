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
%      stop.  All inputs are passed to gbfp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gbfp

% Last Modified by GUIDE v2.5 16-Aug-2017 10:46:15

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


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir(get(handles.fslDir,'String'));
textLabel = sprintf('%s', folder);
set(handles.fslDir, 'string', textLabel);
% create string to append to config.ini
fslpath = sprintf('fslpath=%s',textLabel);
% write to config.ini
fileID = fopen('supp_data/config.ini','a');
fprintf(fileID,'%s\n',fslpath);
fclose(fileID);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% User to select afni directory
folder = uigetdir(get(handles.afniDir,'String'));
% save folder directory in a string
textLabel = sprintf('%s', folder);
% display user selected directory in GUI
set(handles.afniDir, 'string', textLabel);
% create string to append to config.ini
afnipath = sprintf('afnipath=%s',textLabel);
% write to config.ini
fileID = fopen('supp_data/config.ini','a');
fprintf(fileID,'%s\n',afnipath);
fclose(fileID);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir(get(handles.BrainSuiteDir,'String'));
textLabel = sprintf('%s', folder);
set(handles.BrainSuiteDir, 'string', textLabel); 
% create string to append to config.ini
brainsuitepath = sprintf('brainsuitepath=%s',textLabel);
% write to config.ini
fileID = fopen('supp_data/config.ini','a');
fprintf(fileID,'%s\n',brainsuitepath);
fclose(fileID);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path,~] = uigetfile({'*.nii.gz','Compressed NIFTI files (*.nii.gz)'}, ...
    'Choose Structural Image(s)','MultiSelect','on');
textLabel = sprintf('%s%s', path ,file);
set(handles.structural, 'string', textLabel);  
% T1 = sprintf('t1=%s',textLabel);
% % write to config.ini
% fileID = fopen('src/preproc/bfp.mlx','a');
% fprintf(fileID,'%s\n',T1);
% fclose(fileID);
%% This section is supposed to save the T1 and fmri files to the bfp.mlx file.. need to discuss with Anand.

% FileName = 'src/preproc/bfp.mlx';
% S = fileread(FileName);
% S = [textLabel, '\n', S];
% FID = fopen(FileName, 'w');
% if FID == -1, error('Cannot open file %s', FileName); end
% fwrite(FID, S, 'char');
% fclose(FID);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path,~] = uigetfile({'*.nii.gz','Compressed NIFTI files (*.nii.gz)'}, ...
    'Choose Functional Image(s) ','MultiSelect','on');
textLabel = sprintf('%s%s', path, file);
set(handles.functional, 'string', textLabel); 

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
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
% eventdata  reserved - to be defined in a future version of MATLAB
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



function ID_Callback(hObject, eventdata, handles)
% hObject    handle to ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ID as text
%        str2double(get(hObject,'String')) returns contents of ID as a double


% --- Executes during object creation, after setting all properties.
function ID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Needs to get data from updated fields. 6 total. 


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Needs to execute the entie program. Include path to bfp suite
run('src/preproc/bfp.mlx')




% --- Executes on button press in stage.
function stage_Callback(hObject, eventdata, handles)
% hObject    handle to stage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

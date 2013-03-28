function varargout = calim_gui(varargin)
% CALIM_GUI MATLAB code for calim_gui.fig
%      CALIM_GUI, by itself, creates a new CALIM_GUI or raises the existing
%      singleton*.
%
%      H = CALIM_GUI returns the handle to a new CALIM_GUI or the handle to
%      the existing singleton*.
%
%      CALIM_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIM_GUI.M with the given input arguments.
%
%      CALIM_GUI('Property','Value',...) creates a new CALIM_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calim_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calim_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calim_gui

% Last Modified by GUIDE v2.5 17-Oct-2012 19:13:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calim_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @calim_gui_OutputFcn, ...
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


% --- Executes just before calim_gui is made visible.
function calim_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calim_gui (see VARARGIN)

% Choose default command line output for calim_gui
handles.output = hObject;


% ---- CODE ADDED BY PEEYUSH ----
% Initialize handles structure with parameters from GUI. 
addpath ../ % Add calibration code to the path.
handles.path = '~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr.bin';
set (handles.binfile, 'String', handles.path);

handles.do_calib = 1;
set(handles.calib_on, 'String', 'Calib. on');
set(handles.calib_on,'Value',1); % 1 indicates pressed.

handles.do_precalib = 1;
set (handles.precalib_on, 'String', 'PreCal. on');
set (handles.precalib_on, 'Value', 1);

handles.do_image = 1;
set (handles.image_on, 'String', 'Imaging. on');
set (handles.image_on, 'Value', 1);

handles.do_flag = 1;
set (handles.flag_on, 'String', 'Flag. on');
set (handles.flag_on, 'Value', 1);

handles.do_wsfpos = 1;
handles.model_sky_src = 5;
handles.taper = 'uniform';
handles.windowsize = 10;

handles.acc = zeros (288, 288);
handles.t_obs = 0;
handles.freq = 0;
handles.fid = fopen (handles.path, 'rb');
handles.start = 0;

% Update handles structure
guidata(hObject, handles);
master_ctrl(hObject);

% UIWAIT makes calim_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Function to image a visibility set as per parameters in handles.
function showimage (hObject)

function master_ctrl (hObject)
while 1
    pause (1);
    handles = guidata (hObject);
    if handles.start == 1
        [handles.acc, handles.t_obs, handles.freq] = readms2float (handles.fid, -1, -1);
        disp (['Processing timeslice: ' num2str(handles.t_obs)]);
        
    end
end




% --- Outputs from this function are returned to the command line.
function varargout = calim_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in catalog.
function catalog_Callback(hObject, eventdata, handles)
% hObject    handle to catalog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns catalog contents as cell array
%        contents{get(hObject,'Value')} returns selected item from catalog


% --- Executes during object creation, after setting all properties.
function catalog_CreateFcn(hObject, eventdata, handles)
% hObject    handle to catalog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calib_on.
function calib_on_Callback(hObject, eventdata, handles)
% hObject    handle to calib_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calib_on


% --- Executes on key press with focus on calib_on and none of its controls.
function calib_on_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to calib_on (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in precalib_on.
function precalib_on_Callback(hObject, eventdata, handles)
% hObject    handle to precalib_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of precalib_on


% --------------------------------------------------------------------
function uipushtool2_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.start = 1;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in flag_on.
function flag_on_Callback(hObject, eventdata, handles)
% hObject    handle to flag_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flag_on


% --- Executes on button press in image_on.
function image_on_Callback(hObject, eventdata, handles)
% hObject    handle to image_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of image_on


% --- Executes during object creation, after setting all properties.

% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes5


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


% --- Executes during object creation, after setting all properties.
function axes6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes6


% --- Executes during object creation, after setting all properties.
function axes11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes11


% --- Executes during object creation, after setting all properties.
function axes13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes13


% --- Executes during object creation, after setting all properties.
function axes7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes7


% --- Executes during object creation, after setting all properties.
function axes8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes8

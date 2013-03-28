function varargout = imaging(varargin)
% IMAGING MATLAB code for imaging.fig
%      IMAGING, by itself, creates a new IMAGING or raises the existing
%      singleton*.
%
%      H = IMAGING returns the handle to a new IMAGING or the handle to
%      the existing singleton*.
%
%      IMAGING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGING.M with the given input arguments.
%
%      IMAGING('Property','Value',...) creates a new IMAGING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imaging_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imaging_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imaging

% Last Modified by GUIDE v2.5 24-Oct-2012 18:20:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imaging_OpeningFcn, ...
                   'gui_OutputFcn',  @imaging_OutputFcn, ...
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

% --- Executes just before imaging is made visible.
function imaging_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imaging (see VARARGIN)

% Choose default command line output for imaging
handles.output = hObject;

% Code added by Peeyush.
addpath ('~/WORK/AARTFAAC/Afaac_matlab_calib');

% Initialize various parameters, 
handles.path = '~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr_cal.bin';
set (handles.binfile, 'String', handles.path);

handles.fid = fopen (handles.path, 'rb');
% Initial read to determine array sizes.
[handles.acc, handles.tobs, handles.freq] = readms2float (handles.fid, -1, -1);
set(handles.image_on, 'String', 'Start Imaging');
set(handles.image_on,'Value',0); % 0 indicates not pressed.

handles.flagant = [6, 103];
handles.nelem = 288;

% Local horizon based coordinates of array in ITRF
load ('poslocal.mat', 'posITRF', 'poslocal');
% Generate uv coordinates in local horizon coord. system, needed for imaging
handles.normal = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
handles.uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
handles.vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
handles.wloc = meshgrid (poslocal(:,3)) - meshgrid (poslocal (:,3)).';
uvw = [handles.uloc(:), handles.vloc(:), handles.wloc(:)];
handles.uvdist = sqrt(sum(uvw.^2, 2) - (uvw * handles.normal).^2);
clear poslocal, posITRF;
load ('../srclist3CR.mat');
handles.catalog = srclist3CR;
handles.abovejy = 25;
clear srclist3CR;

% Initialize imaging related parameters
handles.duv = 2.5;
handles.Nuv = 500;
handles.uvpad = 512;
handles.calmap = zeros (handles.uvpad);
handles.prevmap = zeros (handles.uvpad);
handles.integrecs = 1;
set(handles.edit1,'String', [num2str(handles.duv) ',' num2str(handles.Nuv) ...
    ',' num2str(handles.uvpad) ',' num2str(handles.integrecs)]);
handles.radec = 0;
handles.overplot3cr = 0;
handles.statcell = 20; % Pixels, to form a region 2*statcell for generating stats
handles.sigthresh = 3; % Threshold for sigma clipping, to generate stats.
handles.mask = eye (handles.nelem);
handles.taper_mask = ones (size (handles.mask)); % Uniform by default
handles.fwhm = max (handles.uvdist(:))/5;
handles.imageinfo = sprintf ('Res: %2.4f deg., DR: %5.3f, rms: %8.5f.', ...
             (2/handles.uvpad)*180, 0, 0);
set (handles.text3, 'String', handles.imageinfo);

handles.play = 0;
% handles.timer = timer('ExecutionMode', 'fixedRate', ... % Run timer repeatedly 
%     'Period', 1, ... % Initial period is 1 sec.
%    'TimerFcn', {@update_display,hObject}); % Specify callback

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using imaging.
if strcmp(get(hObject,'Visible'),'off')

end

% UIWAIT makes imaging wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = imaging_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% while handles.play == 1
         % read in an ACM.
         handles.tobs = 0; handles.acc = zeros (size (handles.acc));
         for ind = 1:handles.integrecs           
           [acc, tobs, freq] = readms2float (handles.fid, -1, -1);
           handles.acc = handles.acc + acc;
           handles.tobs = handles.tobs + tobs;
         end
         handles.acc = handles.acc ./ handles.integrecs;
         handles.tobs = handles.tobs ./ handles.integrecs;
         
         
         [handles.radecmap, handles.calmap, handles.gridvis, handles.l, handles.m] = ... 
         fft_imager_sjw_radec (handles.acc (:).* handles.taper_mask(:), handles.uloc(:), handles.vloc(:), ... 
                              handles.duv, handles.Nuv, handles.uvpad, handles.tobs, handles.freq, handles.radec);
         
         
         axes (handles.axes1);
         imagesc (handles.l, handles.m, abs(handles.calmap));
        
         title (num2str(handles.tobs));
         axis equal;
         axis tight;
         xlabel('South <- m -> North');
         ylabel('East <- l -> West');
         set(colorbar, 'FontSize', 14);
         % Overplotting
         if handles.overplot3cr == 1
             [sel, sel_l, sel_m] = overplot3cr (handles.tobs, handles.catalog, handles.abovejy);
         end
         
         [dr, noise] = getimagedr(handles.calmap, handles.statcell, handles.sigthresh);
         handles.imageinfo = sprintf ('Res: %2.2f deg., DR: %5.3f, rms: %8.5f.', ...
             (2/handles.uvpad)*180, dr, noise);
         set (handles.text3, 'String', handles.imageinfo);
         
         axes (handles.axes3);
         axis tight;
         plot (handles.l, sum(handles.calmap, 1));
         
         axes (handles.axes4);
         axis tight;
         % axis (,  -1, 1);
         plot (sum(handles.calmap, 2), handles.m);         
                 
         
         % Obtain current choice of plot on axes2
         popup_sel_index = get(handles.popupmenu2, 'Value');
         axes (handles.axes2);
         axis equal;
         axis tight;
         switch popup_sel_index
            case 1 % Plot uvdistance                
                t =  handles.acc - diag(diag(handles.acc));
                plot (handles.uvdist(:), abs(t(:) .* handles.taper_mask(:)), '.');
                xlabel('UV Distance (m)');
                ylabel('Visibility amp');
                
            case 2 % Plot real Vs. Imag. visibilities               
                plot (real(handles.acc(:)), imag(handles.acc(:)),'.');
                xlabel('Real');
                ylabel('Imaginary');
                                
            case 3 % Plot visibility taper function
                plot (handles.uvdist(:), handles.taper_mask(:), '.');
                xlabel('UV Distance (m)');
                ylabel('Normalized taper');
                
            case 4 % Plot the PSF due to current taper etc.
               [t,psf,gvis, l, m] = ... 
                  fft_imager_sjw_radec (handles.taper_mask(:), handles.uloc(:), handles.vloc(:), ... 
                              handles.duv, handles.Nuv, handles.uvpad, handles.tobs, handles.freq, handles.radec); 
               imagesc (handles.l, handles.m, 20*log10(abs(psf)/max(psf(:))));
               xlabel('South <- m -> North');
               ylabel('East <- l -> West');
               set(colorbar, 'FontSize', 14);
               
           case 5 % Plot the difference image
               imagesc (handles.l, handles.m, abs(handles.calmap-handles.prevmap));
               xlabel('South <- m -> North');
               ylabel('East <- l -> West');
               set(colorbar, 'FontSize', 14);
                
                
         end
         
         handles.prevmap = handles.calmap;
% end
        


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.path = uigetfile('*.bin');
if ~isequal(handles.path, 0)
    handles.fid = open(handles.path, 'rb');
    set (handles.binfile, 'String', handles.path);
    % Update handles structure
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
% Obtain current taper value.
    popup_sel_index = get(handles.popupmenu1, 'Value');
    switch popup_sel_index
        case 1
            handles.taper_mask = ones (size (handles.mask)); % Uniform
            
        case 2 % Gaussian taper specified by a FWHM
            % fwhm = max(handles.uvdist(:))/10;           
            c = handles.fwhm/(8*log(2));
            handles.taper_mask = exp(-handles.uvdist(:).^2/(2*c.^2));        
            
        case 3 % Triangular taper specified by an end point
            duvdist = -1/max(handles.uvdist(:));
            handles.taper_mask = handles.uvdist(:).*duvdist;
    end
    % Update handles structure
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

% set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on button press in image_on.
function image_on_Callback(hObject, eventdata, handles)
% hObject    handle to image_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of image_on
if (get(hObject, 'Value') == 1)
    % start (handles.timer);
    set(handles.image_on, 'String', 'Stop Imaging');
    set(handles.image_on,'Value',1); % 1 indicates not pressed.
    disp ('---------- >START IMAGING!');
    handles.play = 1;
elseif (get(hObject, 'Value') == 0)
    % stop (handles.timer);
    set(handles.image_on, 'String', 'Start Imaging');
    set(handles.image_on,'Value',0); % 1 indicates not pressed.
    disp ('STOP IMAGING!');
    handles.play = 0;
end
% Update handles structure
guidata(hObject, handles);


function update_display (hObject, eventdata, handles)
        [handles.acc, handles.tobs, handles.freq] = readms2float (handles.fid, -1, -1);
        % [handles.radecmap, handles.calmap, handles.gridvis] = ... 
        % fft_imager_sjw_radec (handles.acc (:), handles.uloc(:), handles.vloc(:), ... 
        %                      handles.duv, handles.Nuv, handles.uvpad, handles.t_obs, handles.freq, handles.radec);        
        axes (handles.axes1);
        imagesc (abs(handles.acc));
        title (num2str(handles.tobs));



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
user_string = get(hObject,'String');
tmp = zeros (1,3);
for i=1:4 % && length(user_string) > 0 
    [tok, user_string] = strtok (user_string, ',');
    tmp(i) = str2double(tok);    
    % disp(['User entered: ' user_string, 'tok:' tok 'rem:' user_string]);
end
handles.duv = tmp(1); handles.Nuv = tmp(2); 
handles.uvpad = tmp(3); handles.integrecs = tmp(4);
% Update handles structure
guidata(hObject, handles);

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


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
if (get(hObject, 'Value') == 1)
    handles.overplot3cr = 1;    
else 
    handles.overplot3cr = 0;
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname] = uigetfile('*.*','Enter data file');
handles.path = [pname fname];
if handles.fid ~= -1
    fclose (handles.fid);
end
handles.fid = fopen (handles.path, 'rb');
set (handles.binfile, 'String', handles.path);
guidata(hObject, handles);

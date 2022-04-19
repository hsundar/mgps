function varargout = sview(varargin)
% SVIEW M-file for sview.fig
%      SVIEW, by itself, creates a new SVIEW or raises the existing
%      singleton*.
%
%      H = SVIEW returns the handle to a new SVIEW or the handle to
%      the existing singleton*.
%
%      SVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SVIEW.M with the given input arguments.
%
%      SVIEW('Property','Value',...) creates a new SVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sview_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sview

% Last Modified by GUIDE v2.5 27-Dec-2004 15:05:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sview_OpeningFcn, ...
                   'gui_OutputFcn',  @sview_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before sview is made visible.
function sview_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sview (see VARARGIN)

% Choose default command line output for sview
handles.output = hObject;

if (length(varargin) == 1 )
    vol = varargin{1};
    minval = min(min(min(vol)));
    maxval = max(max(max(vol)));
    U = (vol - minval)/maxval;
    handles.volume = U;
    sz = size(handles.volume);
    handles.numSlices = sz(3);
    handles.curSlice = 1;
    % Update handles structure
    guidata(hObject, handles);
    % Show the Image ...
    i = floor(handles.curSlice);
    colormap(jet);
    colorbar
    imagesc(handles.volume(:,:,i), [0 1]);
else
    errordlg('Please provide volume file as the only argument')
end

% UIWAIT makes sview wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sview_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

slider_value = get(hObject,'Value');
handles.curSlice = slider_value*handles.numSlices;
guidata(hObject, handles)
set(handles.SliceNumber, 'String', sprintf('%d', floor(handles.curSlice)));
i = floor(handles.curSlice);
if (i<1)
    i = 1;
end
if (i>handles.numSlices)
    i = handles.numSlices;
end
imagesc(handles.volume(:,:,i), [0 1]);
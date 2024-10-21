function varargout = Approximation(varargin)
% APPROXIMATION MATLAB code for Approximation.fig
%      APPROXIMATION, by itself, creates a new APPROXIMATION or raises the existing
%      singleton*.
%
%      H = APPROXIMATION returns the handle to a new APPROXIMATION or the handle to
%      the existing singleton*.
%
%      APPROXIMATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APPROXIMATION.M with the given input arguments.
%
%      APPROXIMATION('Property','Value',...) creates a new APPROXIMATION or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Approximation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Approximation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Approximation

% Last Modified by GUIDE v2.5 16-Feb-2017 14:10:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Approximation_OpeningFcn, ...
                   'gui_OutputFcn',  @Approximation_OutputFcn, ...
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

% --- Executes just before Approximation is made visible.
function Approximation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Approximation (see VARARGIN)

% Choose default command line output for Approximation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global do_app_even_finite; do_app_even_finite = 0;


% --- Outputs from this function are returned to the command line.
function varargout = Approximation_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function No_order_Callback(hObject, eventdata, handles)
global order_No; order_No = str2double(get(hObject,'String'));
function No_order_CreateFcn(hObject, eventdata, handles)
global order_No; order_No = 7;
set(hObject,'String',num2str(order_No));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function No_rel_order_Callback(hObject, eventdata, handles)
global rel_order_No; rel_order_No = str2double(get(hObject,'String'));
function No_rel_order_CreateFcn(hObject, eventdata, handles)
global rel_order_No; rel_order_No = [];
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Hd_order_Callback(hObject, eventdata, handles)
global order_Hd; order_Hd = str2double(get(hObject,'String'));
function Hd_order_CreateFcn(hObject, eventdata, handles)
global order_Hd; order_Hd = 7;
set(hObject,'String',num2str(order_Hd));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Hd_rel_order_Callback(hObject, eventdata, handles)
global rel_order_Hd; rel_order_Hd = str2double(get(hObject,'String'));
% --- Executes during object creation, after setting all properties.
function Hd_rel_order_CreateFcn(hObject, eventdata, handles)
global rel_order_Hd; rel_order_Hd = [];
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function do_app_Callback(hObject, eventdata, handles)
global do_app_even_finite;  do_app_even_finite = get(hObject,'Value');

function continue_to_optimal_Callback(hObject, eventdata, handles)
global order_Hd rel_order_Hd order_No rel_order_No do_app_even_finite;
DESIRED_ORDERS.Hd_order = order_Hd;
DESIRED_ORDERS.Hd_rel_order = rel_order_Hd;
DESIRED_ORDERS.No_order = order_No;
DESIRED_ORDERS.No_rel_order = rel_order_No;
DESIRED_ORDERS.DO_APP_EVEN_FINITE = do_app_even_finite;
assignin('base','DESIRED_ORDERS',DESIRED_ORDERS);
main_gui_object = findobj('Name','Finds Hinf optimal controller for neutral/retarded systems');
if isempty(main_gui_object)
    main_gui_object = findobj('Name','hinfconOnly');
    main_gui_data = guidata(main_gui_object); closereq; drawnow;
    hinfconOnly('controllerCalculator',main_gui_data);
else
    main_gui_data = guidata(main_gui_object); closereq; drawnow;
    Yalta_GUI('controllerCalculator',main_gui_data);
end
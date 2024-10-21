function varargout = hinfconOnly(varargin)
% HINFCONONLY MATLAB code for hinfconOnly.fig
%      HINFCONONLY, by itself, creates a new HINFCONONLY or raises the existing
%      singleton*.
%
%      H = HINFCONONLY returns the handle to a new HINFCONONLY or the handle to
%      the existing singleton*.
%
%      HINFCONONLY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HINFCONONLY.M with the given input arguments.
%
%      HINFCONONLY('Property','Value',...) creates a new HINFCONONLY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before hinfconOnly_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to hinfconOnly_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help hinfconOnly

% Last Modified by GUIDE v2.5 13-Oct-2017 21:30:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @hinfconOnly_OpeningFcn, ...
                   'gui_OutputFcn',  @hinfconOnly_OutputFcn, ...
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


% --- Executes just before hinfconOnly is made visible.
function hinfconOnly_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hinfconOnly (see VARARGIN)
% Choose default command line output for hinfconOnly
global W1affected;
global W2affected;
% global b1;
syms s
W1affected=0;
W2affected=0;
% b1 = 0.0017446*(s+1.992e05)*(s^2 + 0.02931*s + 0.00667)/((s+0.3455)*(s^2 + 0.0002462*s + 0.01279));
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes hinfconOnly wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = hinfconOnly_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function W1_NUM_Callback(hObject, eventdata, handles)
global W1Num;
W1Num_str = get(hObject,'String');
W1Num_sym = str2sym(W1Num_str);
W1Num = sym2poly(W1Num_sym);


% --- Executes during object creation, after setting all properties.
function W1_NUM_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Enable', 'Inactive');



function W1_DEN_Callback(hObject, eventdata, handles)
global W1Den;
W1Den_str = get(hObject,'String');
W1Den_sym = str2sym(W1Den_str);
W1Den = sym2poly(W1Den_sym);


% --- Executes during object creation, after setting all properties.
function W1_DEN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Enable', 'Inactive');



function W2_NUM_Callback(hObject, eventdata, handles)
global W2Num;
qNum_str = get(hObject,'String');
qNum_sym = str2sym(qNum_str);
W2Num = sym2poly(qNum_sym);

% --- Executes during object creation, after setting all properties.
function W2_NUM_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Enable', 'Inactive');



function W2_DEN_Callback(hObject, eventdata, handles)
global W2Den;
W2Den_str = get(hObject,'String');
W2Den_sym = str2sym(W2Den_str);
W2Den = sym2poly(W2Den_sym);

function W2_DEN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Enable', 'Inactive');

function plantButton_Callback(hObject, eventdata, handles)
global Mn_Num;
global Md_Num;
global No_Num;
global Mn_Den;
global Md_Den;
global No_Den;
global Mn;
global Md;
global No;
syms s
try
    Mn = Mn_Num/Mn_Den;
    Md = Md_Num/Md_Den;
    No = No_Num/No_Den;
    Plant = Mn*No/Md;
    Plant = simplify(Plant);
    set(handles.plant_TEXT,'FontName','Courier New');
    set(handles.plant_TEXT,'FontSize',7);
    set(handles.plant_TEXT,'String',evalc('pretty(vpa(Plant,3))'));
    set(handles.gamma_opt_BUTTON,'Enable','On');
    set(handles.plant_PANEL,'Visible','On');
    assignin('base','Mn',Mn);
    assignin('base','Md',Md);
    assignin('base','No',No);
catch
    msgbox('Check numerators and denominators entered for Mn, Md & No')
end

function slider1_Callback(hObject, eventdata, handles)
global wmin;
wmin=round(get(hObject,'Value'));
set(handles.lw1_EDIT,'String',num2str(wmin));
slmin = wmin;
slmax = wmin+10;
set(handles.slider3,'Min',slmin+1);
set(handles.slider3,'Max',slmax);
set(handles.slider3,'SliderStep',[1 1]./(slmax-slmin));
set(handles.slider3,'Value',slmin+1);
set(handles.lw2_EDIT,'String',num2str(slmin+1));

function slider1_CreateFcn(hObject, eventdata, handles)
slmin = -5;
slmax = 3;
set(hObject,'Min',-5);
set(hObject,'Max',3);
set(hObject,'SliderStep',[1 1]./(slmax-slmin));
set(hObject,'Value',-5);
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider3_Callback(hObject, eventdata, handles)
global wmax;
global W1Num;
global W2Num;
global W1Den;
global W2Den;
global W1affected;
global W2affected;
syms s;
wmax=round(get(hObject,'Value'));
wmax=round(get(hObject,'Value'));
if(length(W1Num)~=length(W1Den))||W1affected==1
    W1affected=1;
    if(length(W1Num)~=length(W1Den))    % assumed W1 is strictly proper
        W1Num=[10^(-wmax) W1Num];
    else
        W1Num(1)=10^(-wmax);
    end
    W1Den(end)=10^(-wmax);
end
if(length(W2Num)~=length(W2Den))||W2affected==1
    W2affected=1;
    if length(W2Num)~=length(W2Den)     % assumed W2 is improper
        W2Den=[10^(-wmax) W2Den];
    else
        W2Den(1)=10^(-wmax);
    end
    W2Num(end)=W2Num(1)*10^(-wmax);
end
set(handles.lw2_EDIT,'String',num2str(wmax));
set(handles.W1_NUM,'String',num2str(char(vpa(poly2sym(W1Num,s),2))));
set(handles.W1_DEN,'String',num2str(char(vpa(poly2sym(W1Den,s),2))));
set(handles.W2_NUM,'String',num2str(char(vpa(poly2sym(W2Num,s),2))));
set(handles.W2_DEN,'String',num2str(char(vpa(poly2sym(W2Den,s),2))));

function slider3_CreateFcn(hObject, eventdata, handles)
set(hObject,'Min',0);
set(hObject,'Max',10);
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function lw1_EDIT_Callback(hObject, eventdata, handles)
global wmin;
wmin = str2double(get(hObject,'String'));
set(handles.slider1,'Value',wmin);

function lw1_EDIT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lw2_EDIT_Callback(hObject, eventdata, handles)
global wmax;
wmax = str2double(get(hObject,'String'));
set(handles.slider3,'Value',wmax);

function lw2_EDIT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nuOfPoints_Callback(hObject, eventdata, handles)
global fPoints;
fPoints=str2double(get(hObject,'String'));

function nuOfPoints_CreateFcn(hObject, eventdata, handles)
global fPoints;
fPoints=1000;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minGamma_Callback(hObject, eventdata, handles)
global minGamma;
minGamma = str2double(get(hObject,'String'));

function minGamma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxGamma_Callback(hObject, eventdata, handles)
global maxGamma;
maxGamma = str2double(get(hObject,'String'));

function maxGamma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nuOfGamma_Callback(hObject, eventdata, handles)
global gPoints;
gPoints=str2double(get(hObject,'String'));

function nuOfGamma_CreateFcn(hObject, eventdata, handles)
global gPoints;
gPoints=1000;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function threshold_Callback(hObject, eventdata, handles)
global threshold;
threshold = str2double(get(hObject,'String'));

function threshold_CreateFcn(hObject, eventdata, handles)
global threshold;
threshold = 1e-6;
set(hObject,'String','1e-6');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function W2_NUM_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w2_num=get(hObject,'String');
if strcmp(w2_num,'W2 Numerator')
   set(hObject,'String',''); 
end
uicontrol(hObject);

function W1_NUM_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w1_num=get(hObject,'String');
if strcmp(w1_num,'W1 Numerator')
   set(hObject,'String',''); 
end
uicontrol(hObject);

function W1_DEN_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w1_den=get(hObject,'String');
if strcmp(w1_den,'W1 Denominator')
   set(hObject,'String',''); 
end
uicontrol(hObject);

function W2_DEN_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w2_den=get(hObject,'String');
if strcmp(w2_den,'W2 Denominator')
   set(hObject,'String',''); 
end
uicontrol(hObject);

function Mn_NUM_Callback(hObject, eventdata, handles)
global Mn_Num;
Mn_Num_str = get(hObject,'String');
Mn_Num = str2sym(Mn_Num_str);
Mn_Num = simplify(Mn_Num);

function Mn_NUM_CreateFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'Inactive');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Mn_DEN_Callback(hObject, eventdata, handles)
global Mn_Den;
Mn_Den_str = get(hObject,'String');
Mn_Den = str2sym(Mn_Den_str);
Mn_Den = simplify(Mn_Den);

function Mn_DEN_CreateFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'Inactive');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Md_NUM_Callback(hObject, eventdata, handles)
global Md_Num;
global alpha;
Md_Num_str = get(hObject,'String');
Md_Num = str2sym(Md_Num_str);
Md_Num = simplify(Md_Num);
alpha = roots(sym2poly(Md_Num));

function Md_NUM_CreateFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'Inactive');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Md_DEN_Callback(hObject, eventdata, handles)
global Md_Den;
Md_Den_str = get(hObject,'String');
Md_Den = str2sym(Md_Den_str);
Md_Den = simplify(Md_Den);

function Md_DEN_CreateFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'Inactive');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function No_NUM_Callback(hObject, eventdata, handles)
global No_Num;
No_Num_str = get(hObject,'String');
No_Num = str2sym(No_Num_str);
No_Num = simplify(No_Num);

function No_NUM_CreateFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'Inactive');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function No_DEN_Callback(hObject, eventdata, handles)
global No_Den;
No_Den_str = get(hObject,'String');
No_Den = str2sym(No_Den_str);
No_Den = simplify(No_Den);

function No_DEN_CreateFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'Inactive');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gamma_opt_BUTTON_Callback(hObject, eventdata, handles)
global W1Num;
global W1Den;
global W2Num;
global W2Den;
global alpha;
global gPoints;
global minGamma;
global maxGamma;
global threshold;
global Mn;
global gamma_opt;
set(hinfconOnly, 'HandleVisibility', 'off');
close all;
set(hinfconOnly, 'HandleVisibility', 'on');
    assignin('base','nW1',W1Num);
    assignin('base','dW1',W1Den);
    assignin('base','nW2',W2Num);
    assignin('base','dW2',W2Den);
[gamma_opt,xvalue,yvalue,newFunc,betas] = gammaOpt(W1Num,W1Den,W2Num,W2Den,alpha,gPoints,minGamma,maxGamma,threshold,Mn);
plot(xvalue,yvalue); grid on;
xlabel('\gamma'); ylabel('smin(M)'); title('Optimal case');
set(handles.uipanel5,'Visible','On');
set(handles.gamma_opt_TEXT,'String',['Optimal Gamma     =     ',num2str(gamma_opt,'%1.9f\n')]);
set(handles.gamma_opt_TEXT,'FontSize',14);
set(handles.gamma_opt_TEXT,'FontWeight','bold');
figure(3)
plot(xvalue,yvalue);
xlabel('\gamma');
ylabel('smin(R)');
title('Optimal case');
figure(4)
plot(fliplr(xvalue),newFunc);
xlabel('\gamma');
ylabel('-smin(M)');
title('Negative smin(M) versus \gamma');

function Mn_NUM_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w1_den=get(hObject,'String');
if strcmp(w1_den,'Mn Numerator')
   set(hObject,'String',''); 
end
uicontrol(hObject);

function Mn_DEN_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w1_den=get(hObject,'String');
if strcmp(w1_den,'Mn Denominator')
   set(hObject,'String',''); 
end
uicontrol(hObject);

function Md_NUM_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w1_den=get(hObject,'String');
if strcmp(w1_den,'Md Numerator')
   set(hObject,'String',''); 
end
uicontrol(hObject);

function Md_DEN_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w1_den=get(hObject,'String');
if strcmp(w1_den,'Md Denominator')
   set(hObject,'String',''); 
end
uicontrol(hObject);

function No_NUM_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w1_den=get(hObject,'String');
if strcmp(w1_den,'No Numerator')
   set(hObject,'String',''); 
end
uicontrol(hObject);

function No_DEN_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w1_den=get(hObject,'String');
if strcmp(w1_den,'No Denominator')
   set(hObject,'String',''); 
end
uicontrol(hObject);

function optC_BUTTON_Callback(hObject, eventdata, handles)
run('Approximation.m');

function controllerCalculator(handles)
set(handles.uipanel7,'Visible','off');
global W1Num;
global W1Den;
global W2Num;
global W2Den;
global alpha;
global fPoints;
global wmin;
global wmax;
global threshold;
global Mn;
global Md;
global No;
global gamma_opt;
order_inputs = evalin('base','DESIRED_ORDERS');
[warning, E, G, F, L, Ro, Hn, appHd, Hd, hat_G_gamma, ~, Ca, appNo, Copt]=optController(gamma_opt,W1Num,W1Den,W2Num,W2Den,alpha,threshold,wmin,wmax,fPoints,Mn,Md,No,order_inputs);
if strcmp(warning(1),'C')
    str = evalc('sym2zpk(E)');
    set(handles.E_GAMMA,'FontName','Courier New')
    set(handles.E_GAMMA,'FontSize',8)
    set(handles.E_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(G)');
    set(handles.G_GAMMA,'FontName','Courier New')
    set(handles.G_GAMMA,'FontSize',8)
    set(handles.G_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(F)');
    set(handles.F_GAMMA,'FontName','Courier New')
    set(handles.F_GAMMA,'FontSize',8)
    set(handles.F_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(L)');
    set(handles.L_GAMMA,'FontName','Courier New')
    set(handles.L_GAMMA,'FontSize',8)
    set(handles.L_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(hat_G_gamma)');
    set(handles.G_hat_GAMMA,'FontName','Courier New')
    set(handles.G_hat_GAMMA,'FontSize',8)
    set(handles.G_hat_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(Ro)');
    set(handles.Ro_GAMMA,'FontName','Courier New')
    set(handles.Ro_GAMMA,'FontSize',8)
    set(handles.Ro_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(Hn)');
    set(handles.Hn_GAMMA,'FontName','Courier New')
    set(handles.Hn_GAMMA,'FontSize',8)
    set(handles.Hn_GAMMA,'String',str(10:length(str)-42))
    Co = 1/(1+Hn);
    str = evalc('sym2zpk(Co)');
    set(handles.Co_GAMMA,'FontName','Courier New')
    set(handles.Co_GAMMA,'FontSize',8)
    set(handles.Co_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(appHd)');
    set(handles.Hd_GAMMA,'FontName','Courier New')
    set(handles.Hd_GAMMA,'FontSize',8)
    set(handles.Hd_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(appNo)');
    set(handles.No_GAMMA,'FontName','Courier New')
    set(handles.No_GAMMA,'FontSize',7)
    set(handles.No_GAMMA,'String',str(10:length(str)-42))
    set(handles.uipanel7,'Visible','on');
    set(handles.uipanel12,'Visible','on');
    app_Copt = sym2tf(Ca);
    nume = real(app_Copt.num{1});
    deno = real(app_Copt.den{1});
    app_Copt = tf(nume,deno);
    app_Copt = minreal(app_Copt,1e-3);
    Ca = tf2sym(app_Copt);
    app_Copt = sym2zpk(Ca)
    format long
    
    assignin ('base','E',E)
    assignin ('base','G',G)
    assignin ('base','F',F)
    assignin ('base','L',L)
    assignin ('base','K_opt',1/L)
    assignin ('base','Ro',Ro)
    assignin ('base','Hn',Hn)
    assignin ('base','appHd',appHd)
    assignin ('base','hat_G_gamma',hat_G_gamma)
    assignin ('base','Ca',Ca)
    assignin ('base','Hd',Hd)
    assignin ('base','appNo',appNo)
    assignin ('base','Copt',Copt)
    s= 1e+45;
    dInf = eval(Ro)/eval(L)/gamma_opt;
    assignin ('base','dInf',dInf)
else
    errordlg(warning,'Error')
end

function import_all_variables(handles)
global wmin wmax fpoints minGamma maxGamma gpoints threshold;
global save_qNum alphaNum hNum save_qDen alphaDen hDen ho W1Num W1Den W2Num W2Den;
% try
    inputs = evalin('base','INITIAL_VARIABLES');
    save_qNum = take_struct_field(inputs,'q_numerator',save_qNum,handles.q_NUM,2); 
    save_qDen = take_struct_field(inputs,'q_denominator',save_qDen,handles.q_DEN,2); 
    ho = take_struct_field(inputs,'plant_delay',ho,handles.h_PLANT,0); 
    alphaNum = take_struct_field(inputs,'plant_num_fractional_power',alphaNum,handles.fPower_NUM,0); 
    alphaDen = take_struct_field(inputs,'plant_den_fractional_power',alphaDen,handles.fPower_DEN,0); 
    hNum = take_struct_field(inputs,'plant_num_tau_YALTA',hNum,handles.h_NUM,0); 
    hDen = take_struct_field(inputs,'plant_den_tau_YALTA',hDen,handles.h_DEN,0); 
    W1Num = take_struct_field(inputs,'W1_numerator',W1Num,handles.W1_NUM,1); 
    W1Den = take_struct_field(inputs,'W1_denominator',W1Den,handles.W1_DEN,1); 
    W2Num = take_struct_field(inputs,'W2_numerator',W2Num,handles.W2_NUM,1); 
    W2Den = take_struct_field(inputs,'W2_denominator',W2Den,handles.W2_DEN,1); 
    wmin = take_struct_field(inputs,'log_w_min',wmin,handles.lw1_EDIT,0); 
    wmax = take_struct_field(inputs,'log_w_max',wmax,handles.lw2_EDIT,0); 
    fpoints = take_struct_field(inputs,'Bode_plot_sample_no',fpoints,handles.nuOfPoints,0); 
    minGamma = take_struct_field(inputs,'gamma_min',minGamma,handles.minGamma,0); 
    maxGamma = take_struct_field(inputs,'gamma_max',maxGamma,handles.maxGamma,0); 
    gpoints = take_struct_field(inputs,'Gamma_sample_no',gpoints,handles.nuOfGamma,0); 
    threshold = take_struct_field(inputs,'threshold',threshold,handles.threshold,0); 
% catch err
%     err
%     disp('You need to define a structure called INITIAL_VARIABLES with specific fields defined in the user guide.')
% end
    
    
function stru1 = export_all_variables(myBool)
global wmin wmax fpoints minGamma maxGamma gpoints threshold;
global save_qNum alphaNum hNum save_qDen alphaDen hDen ho W1Num W1Den W2Num W2Den Mn Md No;
global alpha gamma_opt betas;
global E F G L Ro Hn hat_G_gamma Ca appHd Hd appNo Copt dInf K_opt;
stru1 = struct;
stru1 = create_struct_field(stru1,'q_numerator',save_qNum); 
stru1 = create_struct_field(stru1,'q_denominator',save_qDen); 
stru1 = create_struct_field(stru1,'plant_delay',ho); 
stru1 = create_struct_field(stru1,'plant_num_fractional_power',alphaNum); 
stru1 = create_struct_field(stru1,'plant_den_fractional_power',alphaDen); 
stru1 = create_struct_field(stru1,'plant_num_tau_YALTA',hNum); 
stru1 = create_struct_field(stru1,'plant_den_tau_YALTA',hDen); 
stru1 = create_struct_field(stru1,'W1_numerator',W1Num); 
stru1 = create_struct_field(stru1,'W1_denominator',W1Den); 
stru1 = create_struct_field(stru1,'W2_numerator',W2Num); 
stru1 = create_struct_field(stru1,'W2_denominator',W2Den); 
stru1 = create_struct_field(stru1,'log_w_min',wmin); 
stru1 = create_struct_field(stru1,'log_w_max',wmax); 
stru1 = create_struct_field(stru1,'Bode_plot_sample_no',fpoints); 
stru1 = create_struct_field(stru1,'gamma_min',minGamma); 
stru1 = create_struct_field(stru1,'gamma_max',maxGamma); 
stru1 = create_struct_field(stru1,'Gamma_sample_no',gpoints); 
stru1 = create_struct_field(stru1,'threshold',threshold); 
if myBool==0
    assignin('base','INITIAL_VARIABLES',stru1);
else
    stru1 = create_struct_field(stru1,'Mn',Mn); 
    stru1 = create_struct_field(stru1,'Md',Md); 
    stru1 = create_struct_field(stru1,'No',No); 
    stru1 = create_struct_field(stru1,'Noa',appNo);
    stru1 = create_struct_field(stru1,'unstable_poles',alpha);
    stru1 = create_struct_field(stru1,'gamma_opt',gamma_opt);
    stru1 = create_struct_field(stru1,'zeros_of_E',betas);
    stru1 = create_struct_field(stru1,'E',E);
    stru1 = create_struct_field(stru1,'G',G);
    stru1 = create_struct_field(stru1,'F',F);
    stru1 = create_struct_field(stru1,'hat_G',hat_G_gamma);
    stru1 = create_struct_field(stru1,'L',L);
    stru1 = create_struct_field(stru1,'Kopt',K_opt);
    stru1 = create_struct_field(stru1,'Ro',Ro);
    stru1 = create_struct_field(stru1,'d_infty',dInf);
    stru1 = create_struct_field(stru1,'Hn',Hn);
    stru1 = create_struct_field(stru1,'Hd',Hd);
    stru1 = create_struct_field(stru1,'Hda',appHd);
    stru1 = create_struct_field(stru1,'Copt',Copt);
    stru1 = create_struct_field(stru1,'Ca',Ca);
    assignin('base','ALL_VARIABLES',stru1);
end

function structure = create_struct_field(structure,fieldName,variable)
if ~isempty(variable)
    structure = setfield(structure, fieldName, variable);
end

function global_var = take_struct_field(structure,field,global_var,gui_obj,type)
syms s;
if isfield(structure, field)
    global_var = getfield(structure,field);
    if(type==0) % double
        if strcmp(field,'plant_delay')
            set(gui_obj,'String',['-',num2str(global_var)]);
        else
            set(gui_obj,'String',num2str(global_var));
        end
    elseif(type==1) % polynomial
        set(gui_obj,'String',num2str(char(vpa(poly2sym(global_var,s),2))));
    elseif(type==2) % matrix
        result = '';
        for k=size(global_var,1):-1:1
            result = strcat(result,char(vpa(poly2sym(global_var(k,:),s),2)));
            if k~=1
                result = strcat(result,', ');
            end
        end
        set(gui_obj,'String',result);
    end
end


% --------------------------------------------------------------------
function program_options_Callback(hObject, eventdata, handles)
% %


% --------------------------------------------------------------------
function hinfcon_only_Callback(hObject, eventdata, handles)
% hObject    handle to hinfcon_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function import_from_workspace_Callback(hObject, eventdata, handles)
import_all_variables(handles);


% --------------------------------------------------------------------
function export_to_workspace_Callback(hObject, eventdata, handles)
% % 


% --------------------------------------------------------------------
function export_initials_Callback(hObject, eventdata, handles)
export_all_variables(0);


% --------------------------------------------------------------------
function export_all_Callback(hObject, eventdata, handles)
export_all_variables(1);


% --------------------------------------------------------------------
function user_guide_ClickedCallback(hObject, eventdata, handles)
winopen('UserGuide_GUI.pdf')


% --------------------------------------------------------------------
function link_to_src_ClickedCallback(hObject, eventdata, handles)
web('http://www.sciencedirect.com/science/article/pii/S2405896316307376/pdf?md5=14453f144077ac3aa385e435cd797479&amp;pid=1-s2.0-S2405896316307376-main.pdf');


% --------------------------------------------------------------------
function reset_program_ClickedCallback(hObject, eventdata, handles)
clear all;
close all;
run('hinfconOnly.m');


% --------------------------------------------------------------------
function save_all_ClickedCallback(hObject, eventdata, handles)
folder_name = uigetdir;
h = findall(0,'Type','figure');
[fig_titles,file_names] = original_figure_titles();
for k=1:size(h,1)
    if ~isempty(h(k).Number)
        figure(h(k).Number); h1=get(gca,'title'); fig_title=get(h1,'string');
        if isempty(fig_title)
            fig_axes = gca;
            if(strcmp(fig_axes(1).YLabel.String,'Number of unstable poles') && strcmp(fig_axes(1).XLabel.String,'Delay (s)'))
                savefig(h(k),[folder_name,'\',file_names{2}]);
            end
        end
        for n=1:length(fig_titles)
            if strcmp(fig_title,fig_titles{n})
                savefig(h(k),[folder_name,'\',file_names{n}]);
            elseif strcmp(h(k).Name,fig_titles{n})
                savefig(h(k),[folder_name,'\',file_names{n}]);
            end
        end
    else
        if(strcmp(h(k).Name,fig_titles{1}))
            savefig(h(k),[folder_name,'\',file_names{1}]);
        end
    end
end
ALL_VARIABLES = export_all_variables(1);
save([folder_name,'\','ALL_VARIABLES.mat'],'ALL_VARIABLES');


function [fig_titles,file_names] = original_figure_titles()
fig_titles{1} = 'Finds Hinf optimal controller for neutral/retarded systems'; file_names{1} = 'program_window';
fig_titles{2} = 'Stability Window';  file_names{2} = 'stability_window';
fig_titles{3} = 'Root Loci'; file_names{3} = 'root_loci';
fig_titles{4} = 'Optimal case'; file_names{4} = 'gamma_plot';
fig_titles{5} = 'Bode Plots of Controllers'; file_names{5} = 'bode_controller';
fig_titles{6} = 'Nyquist plot with Copt(jw)'; file_names{6}='nyquist_Copt';
fig_titles{7} = 'Nyquist plot with Ca(jw)'; file_names{7}='nyquist_Ca';
fig_titles{8} = 'Bode Plot of Hd(s)'; file_names{8}='bode_approximation_Hd';
fig_titles{9} = 'Bode Plot of No(s)'; file_names{9}='bode_approximation_No';
fig_titles{10} = 'Performance plot'; file_names{10}='performance_plot';

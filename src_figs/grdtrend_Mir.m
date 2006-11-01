function varargout = grdtrend_Mir(varargin)
% M-File changed by desGUIDE 
%   The output is a structure with the following fields (or empty):
%   out.opt_what -> It wil be '-T','-D' or 'W' and optionaly it may have a 'r' apended
%   out.opt_N    -> Number of model parameters (allways)

%	Copyright (c) 2004-2006 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

hObject = figure('Tag','figure1','Visible','off');
handles = guihandles(hObject);
guidata(hObject, handles);
grdtrend_Mir_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
 
movegui(hObject,'center')

% Give a Pro look (3D) to the frame boxes 
bgcolor = get(0,'DefaultUicontrolBackgroundColor');
framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
set(hObject,'Units','pixels')    % Pixels are easier to reason with
h_f = findobj(hObject,'Style','Frame');
for i=1:length(h_f)
    frame_size = get(h_f(i),'Position');
    f_bgc = get(h_f(i),'BackgroundColor');
    usr_d = get(h_f(i),'UserData');
    if abs(f_bgc(1)-bgcolor(1)) > 0.01           % When the frame's background color is not the default's
        frame3D(hObject,frame_size,framecolor,f_bgc,usr_d)
    else
        frame3D(hObject,frame_size,framecolor,'',usr_d)
        delete(h_f(i))
    end
end

% Recopy the text fields on top of previously created frames (uistack is to slow)
h_t = findobj(hObject,'Style','Text');
for i=1:length(h_t)
    usr_d = get(h_t(i),'UserData');
    t_size = get(h_t(i),'Position');   t_str = get(h_t(i),'String');    fw = get(h_t(i),'FontWeight');
    bgc = get (h_t(i),'BackgroundColor');   fgc = get (h_t(i),'ForegroundColor');
    uicontrol('Parent',hObject, 'Style','text', 'Position',t_size,'String',t_str, ...
        'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw,'UserData',usr_d);
end
delete(h_t)

% Choose default command line output for grdtrend_Mir_export
handles.output = hObject;
guidata(hObject, handles);

set(hObject,'Visible','on');
% UIWAIT makes grdtrend_Mir_export wait for user response (see UIRESUME)
uiwait(handles.figure1);

handles = guidata(hObject);
out = grdtrend_Mir_OutputFcn(hObject, [], handles);
varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = grdtrend_Mir_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;
% The figure can be deleted now
delete(handles.figure1);

% -------------------------------------------------------------------------------------
function radiobutton_trend_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_residuals,'Value',0)
    set(handles.radiobutton_weights,'Value',0)
else
    set(hObject,'Value',1)
end

% -------------------------------------------------------------------------------------
function radiobutton_residuals_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_trend,'Value',0)
    set(handles.radiobutton_weights,'Value',0)
else
    set(hObject,'Value',1)
end

% -------------------------------------------------------------------------------------
function radiobutton_weights_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_trend,'Value',0)
    set(handles.radiobutton_residuals,'Value',0)
    set(handles.checkbox_RobustFit,'Value',1)
else
    set(hObject,'Value',1)
end

% -------------------------------------------------------------------------------------
function popup_Nmodel_Callback(hObject, eventdata, handles)
% Nothing to do, the OK button will do the rest

% -------------------------------------------------------------------------------------
function checkbox_RobustFit_Callback(hObject, eventdata, handles)
% Nothing to do, the OK button will do the rest

% -------------------------------------------------------------------------------------
function pushbutton_Help_Nmodel_Callback(hObject, eventdata, handles)
message = {'The trend surface is defined by:'
    ' '
    'm1  +  m2*x + m3*y + m4*x*y + m5*x*x + m6*y*y + m7*x*x*x + m8*x*x*y + m9*x*y*y + m10*y*y*y'
    ''
    'The user must specify "Number of mode parameters", the number of model'
    'parameters to use; thus, 4 fits a bilinear trend, 6 a quadratic surface,'
    'and so on. Optionally, select "Robust fit" to perform a robust fit.'};
helpdlg(message,'Help on model parameters');

% -------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
% See what to compute
if (get(handles.radiobutton_trend,'Value'))
    out.opt_what = '-T';
elseif (get(handles.radiobutton_residuals,'Value'))
    out.opt_what = '-D';
elseif (get(handles.radiobutton_weights,'Value'))
    out.opt_what = '-W';
else
    errordlg('Nothing selected in "What to compute"','Error');  return
end

% Find the model, but first check if "Robust" was choosen
if (get(handles.checkbox_RobustFit,'Value'))
    out.opt_N = '-Nr';
else
    out.opt_N = '-N';
end

val = get(handles.popup_Nmodel,'Value');
str = get(handles.popup_Nmodel, 'String');
switch str{val};
    case '1',       out.opt_N = [out.opt_N '1'];
    case '2',       out.opt_N = [out.opt_N '2'];
    case '3',       out.opt_N = [out.opt_N '3'];
    case '4',       out.opt_N = [out.opt_N '4'];
    case '5',       out.opt_N = [out.opt_N '5'];
    case '6',       out.opt_N = [out.opt_N '6'];
    case '7',       out.opt_N = [out.opt_N '7'];
    case '8',       out.opt_N = [out.opt_N '8'];
    case '9',       out.opt_N = [out.opt_N '9'];
    case '10',      out.opt_N = [out.opt_N '10'];
end

handles.output = out;
guidata(hObject, handles);    uiresume(handles.figure1);

% -------------------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
handles.output = '';        % User gave up, return nothing
guidata(hObject, handles);  uiresume(handles.figure1);

% -------------------------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    handles.output = [];        % User gave up, return nothing
    guidata(hObject, handles);    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    handles.output = [];        % User gave up, return nothing
    guidata(hObject, handles);    delete(handles.figure1);
end

% -------------------------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% Check for "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = '';    % User said no by hitting escape
    guidata(hObject, handles);    uiresume(handles.figure1);
end

% --- Creates and returns a handle to the GUI figure. 
function grdtrend_Mir_LayoutFcn(h1,handles);
set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','grdtrend',...
'NumberTitle','off',...
'Position',[520 686 290 114],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

h2 = uicontrol('Parent',h1,...
'Position',[10 70 271 35],...
'String',{  '' },...
'Style','frame',...
'Tag','frame1');

h3 = uicontrol('Parent',h1,...
'Callback',{@grdtrend_Mir_uicallback,h1,'radiobutton_trend_Callback'},...
'Position',[24 79 51 15],...
'String','Trend',...
'Style','radiobutton',...
'TooltipString','Compute the trend surface resulting from the choosen model',...
'Tag','radiobutton_trend');

h4 = uicontrol('Parent',h1,...
'Callback',{@grdtrend_Mir_uicallback,h1,'radiobutton_residuals_Callback'},...
'Position',[112 79 71 15],...
'String','Residuals',...
'Style','radiobutton',...
'TooltipString','Compute the residuals. That is, The difference (input data - trend)',...
'Tag','radiobutton_residuals');

h5 = uicontrol('Parent',h1,...
'Callback',{@grdtrend_Mir_uicallback,h1,'radiobutton_weights_Callback'},...
'Position',[212 79 61 15],...
'String','Weights',...
'Style','radiobutton',...
'TooltipString','Compute the weights of the fit between the data and the model',...
'Tag','radiobutton_weights');

h6 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdtrend_Mir_uicallback,h1,'popup_Nmodel_Callback'},...
'Position',[10 18 53 22],...
'String',{  '1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_Nmodel');

h7 = uicontrol('Parent',h1,...
'Callback',{@grdtrend_Mir_uicallback,h1,'checkbox_RobustFit_Callback'},...
'Position',[67 23 66 15],...
'String','Robust Fit',...
'Style','checkbox',...
'TooltipString','This is a matematical thing. Either you know what it is or not.',...
'Tag','checkbox_RobustFit');

h8 = uicontrol('Parent',h1,...
'Callback',{@grdtrend_Mir_uicallback,h1,'pushbutton_Help_Nmodel_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[151 18 22 22],...
'String','?',...
'Tag','pushbutton_Help_Nmodel');

h9 = uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Enable','inactive',...
'Position',[11 43 138 15],...
'String','Number of model parameters',...
'Style','text',...
'Tag','text1');

h10 = uicontrol('Parent',h1,...
'Callback',{@grdtrend_Mir_uicallback,h1,'pushbutton_cancel_Callback'},...
'Position',[214 40 66 23],...
'String','Cancel',...
'Tag','pushbutton_cancel');

h11 = uicontrol('Parent',h1,...
'Callback',{@grdtrend_Mir_uicallback,h1,'pushbutton_OK_Callback'},...
'Position',[214 10 66 23],...
'String','OK',...
'Tag','pushbutton_OK');

h12 = uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Enable','inactive',...
'Position',[20 95 91 17],...
'String','What to compute',...
'Style','text',...
'Tag','text2');

function grdtrend_Mir_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

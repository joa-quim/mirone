function varargout = about_box(varargin)
% M-File changed by desGUIDE  
% varargin   command line arguments to about_box (see VARARGIN)

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
about_box_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
% Reposition the window on screen
movegui(hObject,'northeast')

global home_dir

if isempty(home_dir)        % Case when this function was called directly
    f_path = [pwd filesep 'data' filesep];
else
    f_path = [home_dir filesep 'data' filesep];
end

% Load background image
axes(handles.axes1)
logo = imread([f_path 'logo.png']);
image(logo)
set(handles.axes1,'Visible', 'off');

axes(handles.axes2)
fundo = imread([f_path 'cafe_ico.png']);
image(fundo)
set(handles.axes2,'Visible', 'off');

if (~isempty(varargin))
    rem = varargin{1};    text = '';
    while (length(rem) ~= 0)
        [tok,rem] = strtok(rem,'_');
        text = [text ' ' tok];
    end
    set(handles.text_prog,'String',text)
end
if (length(varargin) == 2) 
    str_prog = 'Mirone. The ultimate indescrete grid viewer';
    str_analpha = 'Mirone, Version 1.0.1';
    str_url = 'w3.ualg.pt/~jluis/mirone';
    h1 = findobj(hObject,'Style','text','Tag','text_ProgName');
    h2 = findobj(hObject,'Style','text','Tag','text_AnalphaBeta');
    h3 = findobj(hObject,'Style','text','Tag','text_url');
    set(h1,'String',str_prog)
    set(h2,'String',str_analpha)
    set(h3,'String',str_url)
end

%h = findobj('Tag','text_url');
%makeurl(h,'http://w3.ualg.pt/~jluis/m_gmt');

% Choose default command line output for about_box_export
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes about_box_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = about_box_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = about_box_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
delete(handles.figure1);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% Hint: delete(hObject) closes the figure
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    handles.output = '';        % User gave up, return nothing
    guidata(hObject, handles);
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    handles.output = '';        % User gave up, return nothing
    guidata(hObject, handles);
    delete(handles.figure1);
end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% Check for "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = '';    % User said no by hitting escape
    guidata(hObject, handles);
%    uiresume(handles.figure1);
    delete(handles.figure1);
end

% --- Creates and returns a handle to the GUI figure. 
function about_box_LayoutFcn(h1,handles);

set(h1,...
'PaperUnits','centimeters',...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','about_box',...
'NumberTitle','off',...
'Position',[303 228 269 380],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1');

h2 = uicontrol('Parent',h1,...
'Max',10,...
'Position',[13 99 241 30],...
'String',{  'M-GMT is a MATLAB GUI to the Generic Mapping'; 'Tools (GMT) software'; '' },...
'Style','text',...
'Tag','text_ProgName');

h3 = axes('Parent',h1,'Units','pixels',...
'Color',get(0,'defaultaxesColor'),...
'Position',[15 139 240 240],...
'Tag','axes1');

h8 = uicontrol('Parent',h1,...
'FontSize',11,...
'HorizontalAlignment','left',...
'Position',[11 69 293 20],...
'String','M-GMT, Version Analpha-Beta 0.5 (Built 4)',...
'Style','text',...
'Tag','text_AnalphaBeta');

h9 = axes('Parent',h1,'Units','pixels',...
'Color',get(0,'defaultaxesColor'),...
'Position',[10 49 75 18],'Tag','axes2');

h14 = uicontrol('Parent',h1,...
'FontSize',12,...
'HorizontalAlignment','left',...
'Position',[86 45 45 21],...
'String','J Luis',...
'Style','text',...
'Tag','text3');

h15 = uicontrol('Parent',h1,...
'ForegroundColor',[0 0.501960784313725 0.250980392156863],...
'Position',[138 47 113 17],...
'String','w3.ualg.pt/~jluis/m_gmt',...
'Style','text',...
'Tag','text_url');

h16 = uicontrol('Parent',h1,'Position',[10 91 251 1],...
'String',{''},'Style','frame','Tag','frame3');

h17 = uicontrol('Parent',h1,'Position',[10 41 251 1],...
'String',{''},'Style','frame','Tag','frame4');

h18 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[10 16 201 15],...
'Style','text','Tag','text_prog');

h19 = uicontrol('Parent',h1,...
'Callback',{@about_box_uicallback,h1,'pushbutton_OK_Callback'},...
'Position',[216 11 46 23],...
'String','OK',...
'Tag','pushbutton_OK');

function about_box_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

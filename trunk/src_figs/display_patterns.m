function varargout = display_patterns(varargin)
% M-File changed by desGUIDE 

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
display_patterns_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
 
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to display_patterns_export (see VARARGIN)

global home_dir

if isempty(home_dir)        % Case when this function was called directly
    f_path = ['data' filesep];
else
    f_path = [home_dir filesep 'data' filesep];
end

% Load background image
fundo = imread([f_path 'GMT_patterns.jpg']);
position = get(hObject, 'Position');
image(fundo)
set(handles.axes1,'Visible', 'off');

% Reposition the window on screen
movegui(hObject,'northwest')

% Choose default command line output for display_patterns_export
handles.output = hObject;
guidata(hObject, handles);
set(hObject,'Visible','on');

% UIWAIT makes display_patterns_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = display_patterns_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = display_patterns_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% Check for "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = '';    % User said no by hitting escape
    guidata(hObject, handles);
    delete(hObject);
end

% --- Creates and returns a handle to the GUI figure. 
function display_patterns_LayoutFcn(h1,handles);
set(h1,...
'PaperUnits','centimeters',...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','display_patterns',...
'NumberTitle','off',...
'Position',[114.47204765826 76.6480317721734 497 601],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

setappdata(h1, 'GUIDEOptions',struct(...
'active_h', 107.001831054688, ...
'taginfo', struct(...
'figure', 2, ...
'axes', 2), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'lastSavedFile', 'D:\m_gmt2\display_patterns.m'));

h2 = axes(...
'Parent',h1,...
'Units','pixels',...
'Position',[0 6 495 595],...
'Tag','axes1');

function display_patterns_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

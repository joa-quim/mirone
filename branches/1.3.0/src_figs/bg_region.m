function varargout = bg_region(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to bg_region (see VARARGIN)
% It accepts inputs in dd:mm or dd:mm:ss format

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
bg_region_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
 
global nError
handles.command = cell(15,1);
movegui(hObject,'center');                  % Reposition the window on screen

if ~isempty(varargin) && strcmp(varargin{1},'empty')
    handles.x_min = [];                     handles.x_max = [];
    handles.y_min = [];                     handles.y_max = [];
    handles.command{3} = [];                handles.command{5} = [];
    handles.command{7} = [];                handles.command{9} = [];
    set(handles.edit_Xmin,'String','');     set(handles.edit_Xmax,'String','');
    set(handles.edit_Ymin,'String','');     set(handles.edit_Ymax,'String','');
    set(handles.figure1,'Name','Limits')
elseif ~isempty(varargin) && strcmp(varargin{1},'with_limits')
    tmp = varargin{2};
    handles.x_min = tmp(1);                 handles.x_max = tmp(2);
    handles.y_min = tmp(3);                 handles.y_max = tmp(4);
    handles.command{3} = num2str(tmp(1));   handles.command{5} = num2str(tmp(2));
    handles.command{7} = num2str(tmp(3));   handles.command{9} = num2str(tmp(4));
    set(handles.edit_Xmin,'String',num2str(tmp(1)));   set(handles.edit_Xmax,'String',num2str(tmp(2)));
    set(handles.edit_Ymin,'String',num2str(tmp(3)));   set(handles.edit_Ymax,'String',num2str(tmp(4)));
    set(handles.figure1,'Name','Limits')
    if ((tmp(2) - tmp(1)) > 360 || (tmp(4) - tmp(3)) > 180)  % See if limits rule out "geog"
        set(handles.checkbox_IsGeog,'Value',0)
    end
else
    handles.x_min = -180;                   handles.x_max = 180;
    handles.y_min = -90;                    handles.y_max = 90;
    handles.command{3} = '-180';            handles.command{5} = '180';
    handles.command{7} = '-90';             handles.command{9} = '90';
end
% Choose default command line output for bg_region_export
handles.output = hObject;
guidata(hObject, handles);

set(hObject,'Visible','on');
% UIWAIT makes bg_region_export wait for user response (see UIRESUME)
uiwait(handles.figure1);

handles = guidata(hObject);
out = bg_region_OutputFcn(hObject, [], handles);
varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = bg_region_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure

% Get default command line output from handles structure
varargout{1} = handles.output;
% The figure can be deleted now
delete(handles.figure1);

% --------------------------------------------------------------------------------------------------
function edit_Xmin_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)            % when dd:mm or dd:mm:ss was given
    x_min = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val),  x_min = x_min + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val),  x_min = x_min - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.x_min = x_min;
    if ~isempty(handles.x_max) && x_min >= handles.x_max
        errordlg('West Longitude >= East Longitude ','Error in Longitude limits')
        handles.command{3} = '';
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        handles.command{3} = xx;    % Save entered string
        guidata(hObject, handles);
    end
    % Guess if we are probably dealing with geog coordinates
    if ~isempty(get(handles.edit_Xmax,'String')) && (handles.x_max - x_min) <= 360
        set(handles.checkbox_IsGeog,'Value',1)
    else
        set(handles.checkbox_IsGeog,'Value',0)
    end
else                % box is empty
    handles.command{3} = '';
    set(hObject,'String','');   guidata(hObject, handles);
end

% --------------------------------------------------------------------------------------------------
function edit_Xmax_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)            % when dd:mm or dd:mm:ss was given
    x_max = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val),  x_max = x_max + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val),  x_max = x_max - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.x_max = x_max;
    if ~isempty(handles.x_min) && x_max <= handles.x_min
        errordlg('East Longitude <= West Longitude','Error in Longitude limits')
        handles.command{5} = '';
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        handles.command{5} = xx;    % Save entered string
        guidata(hObject, handles);
    end
    % Guess if we are probably dealing with geog coordinates
    if ~isempty(get(handles.edit_Xmin,'String')) && (x_max - handles.x_min) <= 360
        set(handles.checkbox_IsGeog,'Value',1)
    else
        set(handles.checkbox_IsGeog,'Value',0)
    end
else                % box is empty
    handles.command{5} = '';
    set(hObject,'String','');   guidata(hObject, handles);
end

% --------------------------------------------------------------------------------------------------
function edit_Ymin_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)            % when dd:mm or dd:mm:ss was given
    y_min = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val),  y_min = y_min + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val),  y_min = y_min - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.y_min = y_min;
    if ~isempty(handles.y_max) && y_min >= handles.y_max
        errordlg('South Latitude >= North Latitude','Error in Latitude limits')
        handles.command{7} = '';
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        handles.command{7} = xx;    % Save entered string
        guidata(hObject, handles);
    end
    % Guess if we are probably dealing with geog coordinates
    if ~isempty(get(handles.edit_Ymax,'String')) && (handles.y_max - y_min) <= 180
        set(handles.checkbox_IsGeog,'Value',1)
    else
        set(handles.checkbox_IsGeog,'Value',0)
    end
else                % box is empty
    handles.command{7} = '';
    set(hObject,'String','');   guidata(hObject, handles);
end

% --------------------------------------------------------------------------------------------------
function edit_Ymax_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)            % when dd:mm or dd:mm:ss was given
    y_max = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val),   y_max = y_max + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val),   y_max = y_max - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.y_max = y_max;
    if ~isempty(handles.y_min) && y_max <= handles.y_min
        errordlg('North Latitude <= South Latitude','Error in Latitude limits')
        handles.command{9} = '';
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        handles.command{9} = xx;    % Save entered string
        guidata(hObject, handles);
    end
    % Guess if we are probably dealing with geog coordinates
    if ~isempty(get(handles.edit_Ymin,'String')) && (y_max - handles.y_min) <= 180
        set(handles.checkbox_IsGeog,'Value',1)
    else
        set(handles.checkbox_IsGeog,'Value',0)
    end
else                % box is empty
    handles.command{9} = '';
    set(hObject,'String','');   guidata(hObject, handles);
end

% --------------------------------------------------------------------------------------------------
function VerifyCommand(handles)
global nError
% ERROR TESTING
error = 0;
if isempty(handles.command{3})
    errordlg('Lon Min box is empty','Error');
    error = error + 1;    
end
if isempty(handles.command{5})
    errordlg('Lon Max box is empty','Error');
    error = error + 1;    
end
if isempty(handles.command{7})
    errordlg('Lat Min box is empty','Error');
    error = error + 1;    
end
if isempty(handles.command{9})
    errordlg('Lat Max box is empty','Error');
    error = error + 1;    
end
nError = error;

% --------------------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
VerifyCommand(handles)
global nError
if nError == 0
    is_geog = get(handles.checkbox_IsGeog,'Value');
    handles.output = [handles.x_min handles.x_max handles.y_min handles.y_max is_geog];
    guidata(hObject,handles);
    uiresume(handles.figure1);
else
    nError = 0;
end

% --------------------------------------------------------------------------------------------------
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
handles.output = [];        % User gave up, return nothing
guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit8_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit9_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit10_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit11_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------------------------------------
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

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = [];    % User said no by hitting escape
    guidata(hObject, handles);    uiresume(handles.figure1);
end

% --------------------------------------------------------------------------------------------------
function checkbox_IsGeog_Callback(hObject, eventdata, handles)
% Don't need any code because the output will be checked by pushbutton_OK and
% the only thing that matters is if this is checked or not.


% --- Creates and returns a handle to the GUI figure. 
function bg_region_LayoutFcn(h1,handles)
set(h1,...
'PaperUnits','centimeters',...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','bg_region',...
'NumberTitle','off',...
'Position',[520 659 258 141],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bg_region_uicallback,h1,'edit_Ymax_Callback'},...
'FontSize',10,...
'Position',[123 113 79 26],...
'String','90',...
'Style','edit',...
'TooltipString','Enter nothern map limit',...
'Tag','edit_Ymax');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bg_region_uicallback,h1,'edit_Xmin_Callback'},...
'CData',[],...
'FontSize',10,...
'Position',[79 78 79 26],...
'String','-180',...
'Style','edit',...
'TooltipString','Enter western map limit',...
'Tag','edit_Xmin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bg_region_uicallback,h1,'edit_Xmax_Callback'},...
'FontSize',10,...
'Position',[168 78 79 26],...
'String','180',...
'Style','edit',...
'TooltipString','Enter eastern map limit',...
'Tag','edit_Xmax');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bg_region_uicallback,h1,'edit_Ymin_Callback'},...
'FontSize',10,...
'Position',[124 43 79 26],...
'String','-90',...
'Style','edit',...
'TooltipString','Enter southern map limit',...
'Tag','edit_Ymin');

uicontrol('Parent',h1,...
'Callback',{@bg_region_uicallback,h1,'pushbutton_OK_Callback'},...
'Position',[185 8 66 23],...
'String','OK',...
'Tag','pushbutton_OK');

uicontrol('Parent',h1,...
'Callback',{@bg_region_uicallback,h1,'pushbutton_Cancel_Callback'},...
'Position',[110 8 66 23],...
'String','Cancel',...
'Tag','pushbutton_Cancel');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bg_region_uicallback,h1,'edit8_Callback'},...
'Enable','off',...
'FontSize',10,...
'Position',[25 104 22 22],...
'String','N',...
'Style','edit',...
'Tag','edit8');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bg_region_uicallback,h1,'edit9_Callback'},...
'Enable','off',...
'FontSize',10,...
'Position',[9 80 24 22],...
'String','W',...
'Style','edit',...
'Tag','edit9');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bg_region_uicallback,h1,'edit10_Callback'},...
'Enable','off',...
'FontSize',10,...
'Position',[39 80 24 22],...
'String','E',...
'Style','edit',...
'Tag','edit10');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bg_region_uicallback,h1,'edit11_Callback'},...
'Enable','off',...
'FontSize',10,...
'Position',[25 57 22 22],...
'String','S',...
'Style','edit',...
'Tag','edit11');

uicontrol('Parent',h1,...
'Position',[10 32 91 15],...
'String','Is Geographic?',...
'Style','checkbox',...
'Value',1,...
'TooltipString','Check if limits are in geographic coordinates',...
'Tag','checkbox_IsGeog');

function bg_region_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

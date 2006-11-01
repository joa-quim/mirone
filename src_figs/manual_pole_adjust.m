function varargout = manual_pole_adjust(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to manual_pole_adjust (see VARARGIN)

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
manual_pole_adjust_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

movegui(hObject,'center');

if (~isempty(varargin))
    handles.h_calling_fig = varargin{1};
    handles.mirone_axes = get(varargin{1},'CurrentAxes');
else
    errordlg('EULER_STUFF: wrong number of arguments.','Error')
end

handles.h_line = [];
handles.line_x = [];
handles.line_y = [];
handles.have_pole = 0;
handles.p_lon = 0;      % Default to this non-sense
handles.p_lat = 0;
handles.p_omega = 0;
handles.h_active_line_str = findobj(handles.figure1,'Tag','text_activeLine');      % Get this handle
handles.path_continent = [pwd filesep 'continents' filesep];

set(handles.edit_lon, 'String', num2str(handles.p_lon))
set(handles.edit_lat, 'String', num2str(handles.p_lat))
set(handles.edit_omega, 'String', num2str(handles.p_omega))

set(handles.slider_lon, 'Value', handles.p_lon)
set(handles.slider_lat, 'Value', handles.p_lat)
set(handles.slider_omega, 'Value', handles.p_omega)

% Choose default command line output for manual_pole_adjust_export
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes earthquakes wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = manual_pole_adjust_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --------------------------------------------------------------------------
function varargout = manual_pole_adjust_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------------
function edit_lon_Callback(hObject, eventdata, handles)
x = str2double(get(hObject,'String'));
set(handles.slider_lon,'Value',x)
handles.p_lon = x;
if (isempty(handles.p_lat) | isempty(handles.p_omega))
    handles.have_pole = 0;
else                            % OK, we have all the pole parameters
    handles.have_pole = 1;
end
guidata(hObject, handles);
if (~handles.have_pole)         return;     end     % We don't still have a pole
if (isempty(handles.h_line))    return;     end     % There is nothing to do yet
[rlon,rlat] = rot_euler(handles.line_x,handles.line_y,handles.p_lon,handles.p_lat,handles.p_omega);
try,    set(handles.h_line,'XData',rlon,'YData',rlat);      end     % Use a 'try' because line may have been killed

% --------------------------------------------------------------------------
function slider_lon_Callback(hObject, eventdata, handles)
if (isempty(handles.h_line))
    set(hObject,'Value',0);         return;
end
if (~handles.have_pole)     return;     end     % We don't still have a pole
val = get(hObject,'Value');
set(handles.edit_lon,'String',num2str(val))
handles.p_lon = val;
[rlon,rlat] = rot_euler(handles.line_x,handles.line_y,handles.p_lon,handles.p_lat,handles.p_omega);
try,    set(handles.h_line,'XData',rlon,'YData',rlat);      end     % Use a 'try' because line may have been killed
guidata(hObject, handles);

% --------------------------------------------------------------------------
function edit_lat_Callback(hObject, eventdata, handles)
x = str2double(get(hObject,'String'));
set(handles.slider_lat,'Value',x)
handles.p_lat = x;
if (isempty(handles.p_lon) | isempty(handles.p_omega))
    handles.have_pole = 0;
else                            % OK, we have all the pole parameters
    handles.have_pole = 1;
end
guidata(hObject, handles);
if (~handles.have_pole)         return;     end     % We don't still have a pole
if (isempty(handles.h_line))    return;     end     % There is nothing to do yet
[rlon,rlat] = rot_euler(handles.line_x,handles.line_y,handles.p_lon,handles.p_lat,handles.p_omega);
try,    set(handles.h_line,'XData',rlon,'YData',rlat);      end     % Use a 'try' because line may have been killed

% --------------------------------------------------------------------------
function slider_lat_Callback(hObject, eventdata, handles)
if (isempty(handles.h_line))
    set(hObject,'Value',0);         return;
end
if (~handles.have_pole)     return;     end     % We don't still have a pole
val = get(hObject,'Value');
set(handles.edit_lat,'String',num2str(val))
handles.p_lat = val;
[rlon,rlat] = rot_euler(handles.line_x,handles.line_y,handles.p_lon,handles.p_lat,handles.p_omega);
try,    set(handles.h_line,'XData',rlon,'YData',rlat);      end     % Use a 'try' because line may have been killed
guidata(hObject, handles);

% --------------------------------------------------------------------------
function edit_omega_Callback(hObject, eventdata, handles)
x = str2double(get(hObject,'String'));
set(handles.slider_omega,'Value',x)
handles.p_omega = x;
if (isempty(handles.p_lon) | isempty(handles.p_lat))
    handles.have_pole = 0;
else                            % OK, we have all the pole parameters
    handles.have_pole = 1;
end
guidata(hObject, handles);
if (~handles.have_pole)         return;     end     % We don't still have a pole
if (isempty(handles.h_line))    return;     end     % There is nothing to do yet
[rlon,rlat] = rot_euler(handles.line_x,handles.line_y,handles.p_lon,handles.p_lat,handles.p_omega);
try,    set(handles.h_line,'XData',rlon,'YData',rlat);      end     % Use a 'try' because line may have been killed

% --------------------------------------------------------------------------
function slider_omega_Callback(hObject, eventdata, handles)
if (isempty(handles.h_line))
    set(hObject,'Value',0);         return;
end
if (~handles.have_pole)     return;     end     % We don't still have a pole
val = get(hObject,'Value');
set(handles.edit_omega,'String',num2str(val))
handles.p_omega = val;
[rlon,rlat] = rot_euler(handles.line_x,handles.line_y,handles.p_lon,handles.p_lat,handles.p_omega);
try,    set(handles.h_line,'XData',rlon,'YData',rlat);      end     % Use a 'try' because line may have been killed
guidata(hObject, handles);

% --------------------------------------------------------------------------
function pushbutton_polesList_Callback(hObject, eventdata, handles)
fid = fopen([handles.path_continent 'lista_polos.dat'],'rt');
c = fread(fid,'*char').';
fclose(fid);
s = strread(c,'%s','delimiter','\n');

[s,v] = choosebox('Name','One Euler list',...
                    'PromptString','List of poles:',...
                    'SelectString','Selected poles:',...
                    'ListSize',[380 300],...
                    'ListString',s);

if (v == 1)         % Finite pole
    handles.p_lon = s(1);
    handles.p_lat = s(2);
    handles.p_omega = s(3);
    set(handles.edit_lon, 'String', num2str(s(1)))
    set(handles.edit_lat, 'String', num2str(s(2)))
    set(handles.edit_omega, 'String', num2str(s(3)))
    handles.have_pole = 1;
    guidata(hObject,handles)
else                % Stage poles or cancel
    handles.have_pole = 0;
    return;
    %set(handles.edit_polesFile,'String',s)
end

set(handles.slider_lon, 'Value', s(1))
set(handles.slider_lat, 'Value', s(2))
set(handles.slider_omega, 'Value', s(3))

% --------------------------------------------------------------------------
function togglebutton_pickLine_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    % Test if we have potential target lines and their type
    h_mir_lines = findobj(handles.h_calling_fig,'Type','line');     % Fish all objects of type line in Mirone figure
    if (isempty(h_mir_lines))                                       % We don't have any lines
        str = ['If you hited this button on purpose, than you deserve the following insult.',...
                'You #!|"*!%!?~^)--$&.',... 
                'THERE ARE NO LINES IN THAT FIGURE.'];
        errordlg(str,'Chico Clever');   set(hObject,'Value',0);     return;
    end
    % The above test is not enough. For exemple, coastlines are not eligible neither,
    % but is very cumbersome to test all the possibilities of pure non-eligible lines.
    set(handles.h_calling_fig,'pointer','crosshair')
    h_line = get_polygon(handles.h_calling_fig);          % Get the line handle
    if (~isempty(h_line))
        x = get(h_line,'XData');        y = get(h_line,'YData');
        handles.line_x = x(:);          handles.line_y = y(:);
        % Create a empty line handle that will hold the rotated line
        handles.h_line = line('parent',get(handles.h_calling_fig,'CurrentAxes'),'XData',[],'YData',[], ...
            'LineStyle','-.','LineWidth',2);
    else
        handles.line_x = [];                handles.line_y = [];
        set(hObject,'Value',0)
    end
    set(handles.h_calling_fig,'pointer','arrow')
    set(hObject,'Value',0)
    set(handles.h_active_line_str,'String','GOT A LINE TO WORK WITH')
    figure(handles.figure1)                 % Bring this figure to front again
else        % What should I do?
    %handles.do_graphic = 0;
end
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')
    delete(handles.figure1);
end


% --- Creates and returns a handle to the GUI figure. 
function manual_pole_adjust_LayoutFcn(h1,handles);

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','Manual pole adjust',...
'NumberTitle','off',...
'Position',[520 660 650 140],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

h2 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@manual_pole_adjust_uicallback,h1,'edit_lon_Callback'},...
'Position',[558 73 81 21],...
'Style','edit',...
'Tag','edit_lon');

h3 = uicontrol('Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',{@manual_pole_adjust_uicallback,h1,'slider_lon_Callback'},...
'Max',360,...
'Min',-180,...
'Position',[78 75 471 18],...
'String',{  '' },...
'Style','slider',...
'SliderStep',[0.001 0.01],...
'Tag','slider_lon');

h4 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@manual_pole_adjust_uicallback,h1,'edit_lat_Callback'},...
'Position',[558 42 81 21],...
'Style','edit',...
'Tag','edit_lat');

h5 = uicontrol('Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',{@manual_pole_adjust_uicallback,h1,'slider_lat_Callback'},...
'Max',90,...
'Min',-90,...
'Position',[78 44 471 18],...
'String',{  '' },...
'Style','slider',...
'SliderStep',[0.001 0.01],...
'Tag','slider_lat');

h6 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@manual_pole_adjust_uicallback,h1,'edit_omega_Callback'},...
'Position',[558 13 81 21],...
'Style','edit',...
'Tag','edit_omega');

h7 = uicontrol('Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',{@manual_pole_adjust_uicallback,h1,'slider_omega_Callback'},...
'Max',90,...
'Min',-90,...
'Position',[78 15 471 18],...
'String',{  '' },...
'Style','slider',...
'SliderStep',[0.0005 0.005],...
'Tag','slider_omega');

h8 = uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[10 75 55 16],...
'String','Longitude',...
'Style','text',...
'Tag','text1');

h9 = uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[10 45 55 16],...
'String','Latitude',...
'Style','text',...
'Tag','text2');

h10 = uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[11 17 55 16],...
'String','Angle',...
'Style','text',...
'Tag','text3');

h11 = uicontrol('Parent',h1,...
'Callback',{@manual_pole_adjust_uicallback,h1,'pushbutton_polesList_Callback'},...
'Position',[520 107 121 23],...
'String','Poles selector',...
'TooltipString','Select a pole from the default list',...
'Tag','pushbutton_polesList');

h12 = uicontrol('Parent',h1,...
'Callback',{@manual_pole_adjust_uicallback,h1,'togglebutton_pickLine_Callback'},...
'Position',[10 107 121 23],...
'String','Pick line from Figure',...
'TooltipString','Allows you to mouse select one line from a Mirone figure',...
'Tag','togglebutton_pickLine');

h13 = uicontrol('Parent',h1,...
'FontSize',10,'Position',[240 109 231 17],...
'String','NO ACTIVE LINE',...
'Style','text',...
'Tag','text_activeLine');

function manual_pole_adjust_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

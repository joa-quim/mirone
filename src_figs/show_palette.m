function varargout = show_palette(varargin)
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
show_palette_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
 
global home_dir
if (isempty(home_dir))        % Case when this function was called directly
    handles.home_dir = pwd;
else
    handles.home_dir = home_dir;
end

h_m_gmt = findobj('Name','m_gmt');
if ishandle(h_m_gmt)
    handles.work_dir = getappdata(h_m_gmt,'current_dir');
    handles.last_dir = handles.work_dir;
else
    handles.work_dir = pwd;
    handles.last_dir = pwd;
end

movegui(hObject,'northeast');     % Reposition the window on screen
cmaplim = varargin{1};      datalim = varargin{2};
datainc = varargin{3};      cmap = varargin{4};

set(hObject,'Name','');     set(handles.figure1,'Colormap',cmap)
% Determine tick labels
tickvalues = datainc * [ceil(datalim(1)/datainc) : floor(datalim(2)/datainc)];
for i = 1 : length(tickvalues)
  ticklabel{i} = num2str(tickvalues(i));
end

% Determine tick positions and image array
if cmaplim(1) < cmaplim(2)
    tick = cmaplim(1) + (tickvalues - datalim(1)) * diff(cmaplim) / diff(datalim);
    I = cmaplim(1):cmaplim(2);    ilim = cmaplim;
else
    tick = cmaplim(2) - (tickvalues - datalim(1)) * diff(cmaplim) / diff(datalim);
    I = cmaplim(1):-1:cmaplim(2);    ilim = [cmaplim(2) cmaplim(1)];
end

% Create the color bar axes and image object
h = axes('Units','Pixels','Visible','off');
image(I','CDataMapping','direct','YData',ilim);
set(h,'YTick',tick,'XTick',[],'YTickLabel',ticklabel);

pos = get(h,'Position');
% Set axes additional properties
set(h,'YDir','normal'); % Override the default YDir set by IMAGE
set(h,'Position', pos+[55 -30 -45 +40]); % Separate in case we have {...,'YDir','reverse',...}
set(h,'Visible','On');  % Always make it visible

set(hObject,'Visible','on');

% Choose default command line output for show_palette_export
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes shading_params wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = show_palette_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = show_palette_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')      % Check for "escape"
    delete(hObject);
end


% --- Creates and returns a handle to the GUI figure. 
function show_palette_LayoutFcn(h1,handles);

set(h1,...
'PaperUnits','centimeters',...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','show_palette',...
'NumberTitle','off',...
'Position',[100 100 106 420],...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

setappdata(h1, 'GUIDEOptions',struct(...
'active_h', 102.003295898438, ...
'taginfo', struct(...
'figure', 2), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'lastSavedFile', 'D:\m_gmt\show_palette.m'));

function show_palette_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

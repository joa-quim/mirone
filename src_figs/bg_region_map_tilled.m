function varargout = bg_region_map_tilled(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to bg_region_map_tilled (see VARARGIN)

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
bg_region_map_tilled_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
 
global home_dir

% Case when this function was called directly
if isempty(home_dir),   f_path = ['data' filesep];
else                    f_path = [home_dir filesep 'data' filesep];   end
handles.f_path = f_path;
handles.first_WT = 1;       % Flag to signal the first time World Map is used
handels.h_ll     = [];      % Will contail the handles to the Lower Left pushbutton

set(hObject,'Name','World Topo Tiles')
load([f_path 'WorldMapTilled.mat'])
set_tiles(handles)              % Put the logo tiles in the pusbuttons

movegui(hObject,'center');      % Reposition the window on screen

% Choose default command line output for bg_region_map_tilled_export
handles.output = hObject;
guidata(hObject, handles);

set(hObject,'Visible','on');
% UIWAIT makes bg_region_map_tilled_export wait for user response (see UIRESUME)
uiwait(handles.figure1);


handles = guidata(hObject);
out = bg_region_map_tilled_OutputFcn(hObject, [], handles);
varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = bg_region_map_tilled_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in pushbutton_90N180W.
function pushbutton_90N180W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '90N180W.jpg']),1 );
out.X = [-180 -135];  out.Y = [45 90];    out.Iname = [handles.f_path '90N180W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_90N135W.
function pushbutton_90N135W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '90N135W.jpg']),1 );
out.X = [-135 -90];  out.Y = [45 90];    out.Iname = [handles.f_path '90N135W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_90N090W.
function pushbutton_90N090W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '90N090W.jpg']),1 );
out.X = [-90 -45];  out.Y = [45 90];    out.Iname = [handles.f_path '90N090W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_90N045W.
function pushbutton_90N045W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '90N045W.jpg']),1 );
out.X = [-45 0];  out.Y = [45 90];    out.Iname = [handles.f_path '90N045W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_90N000E.
function pushbutton_90N000E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '90N000E.jpg']),1 );
out.X = [0 45];  out.Y = [45 90];    out.Iname = [handles.f_path '90N000E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_90N045E.
function pushbutton_90N045E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '90N045E.jpg']),1 );
out.X = [45 90];  out.Y = [45 90];    out.Iname = [handles.f_path '90N045E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_90N090E.
function pushbutton_90N090E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '90N090E.jpg']),1 );
out.X = [90 135];  out.Y = [45 90];    out.Iname = [handles.f_path '90N090E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_90N135E.
function pushbutton_90N135E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '90N135E.jpg']),1 );
out.X = [135 180];  out.Y = [45 90];    out.Iname = [handles.f_path '90N135E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45N180W.
function pushbutton_45N180W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45N180W.jpg']),1 );
out.X = [-180 -135];  out.Y = [0 45];    out.Iname = [handles.f_path '45N180W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45N135W.
function pushbutton_45N135W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45N135W.jpg']),1 );
out.X = [-135 -90];  out.Y = [0 45];    out.Iname = [handles.f_path '45N135W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45N090W.
function pushbutton_45N090W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45N090W.jpg']),1 );
out.X = [-90 -45];  out.Y = [0 45];    out.Iname = [handles.f_path '45N090W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45N045W.
function pushbutton_45N045W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45N045W.jpg']),1 );
out.X = [-45 0];  out.Y = [0 45];    out.Iname = [handles.f_path '45N045W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45N000E.
function pushbutton_45N000E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45N000E.jpg']),1 );
out.X = [0 45];  out.Y = [0 45];    out.Iname = [handles.f_path '45N000E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45N045E.
function pushbutton_45N045E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45N045E.jpg']),1 );
out.X = [45 90];  out.Y = [0 45];    out.Iname = [handles.f_path '45N045E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45N090E.
function pushbutton_45N090E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45N090E.jpg']),1 );
out.X = [90 135];  out.Y = [0 45];    out.Iname = [handles.f_path '45N090E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45N135E.
function pushbutton_45N135E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45N135E.jpg']),1 );
out.X = [135 180];  out.Y = [0 45];    out.Iname = [handles.f_path '45N135E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_00N180W.
function pushbutton_00N180W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '00N180W.jpg']),1 );
out.X = [-180 -135];  out.Y = [-45 0];    out.Iname = [handles.f_path '00N180W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_00N135W.
function pushbutton_00N135W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '00N135W.jpg']),1 );
out.X = [-135 -90];  out.Y = [-45 0];    out.Iname = [handles.f_path '00N135W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_00N090W.
function pushbutton_00N090W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '00N090W.jpg']),1 );
out.X = [-90 -45];  out.Y = [-45 0];    out.Iname = [handles.f_path '00N090W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_00N045W.
function pushbutton_00N045W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '00N045W.jpg']),1 );
out.X = [-45 0];  out.Y = [-45 0];    out.Iname = [handles.f_path '00N045W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_00N000E.
function pushbutton_00N000E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '00N000E.jpg']),1 );
out.X = [0 45];  out.Y = [-45 0];    out.Iname = [handles.f_path '00N000E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_00N045E.
function pushbutton_00N045E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '00N045E.jpg']),1 );
out.X = [45 90];  out.Y = [-45 0];    out.Iname = [handles.f_path '00N045E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_00N090E.
function pushbutton_00N090E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '00N090E.jpg']),1 );
out.X = [90 135];  out.Y = [-45 0];    out.Iname = [handles.f_path '00N090E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_00N135E.
function pushbutton_00N135E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '00N135E.jpg']),1 );
out.X = [135 180];  out.Y = [-45 0];    out.Iname = [handles.f_path '00N135E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45S180W.
function pushbutton_45S180W_Callback(hObject, eventdata, handles)
if (get(handles.radiobutton_WorldMap,'Value'))
    out.img = flipdim( imread([handles.f_path 'etopo2.jpg']),1 );
    out.X = [-180 179.967];  out.Y = [-90 90];    out.Iname = [handles.f_path 'etopo2.jpg'];
else
    out.img = flipdim( imread([handles.f_path '45S180W.jpg']),1 );
    out.X = [-180 -135];  out.Y = [-90 -45];    out.Iname = [handles.f_path '45S180W.jpg'];
end
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45S135W.
function pushbutton_45S135W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45S135W.jpg']),1 );
out.X = [-135 -90];  out.Y = [-90 -45];    out.Iname = [handles.f_path '45S135W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45S090W.
function pushbutton_45S090W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45S090W.jpg']),1 );
out.X = [-90 -45];  out.Y = [-90 -45];    out.Iname = [handles.f_path '45S090W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45S045W.
function pushbutton_45S045W_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45S045W.jpg']),1 );
out.X = [-45 0];  out.Y = [-90 -45];    out.Iname = [handles.f_path '45S045W.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45S000E.
function pushbutton_45S000E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45S000E.jpg']),1 );
out.X = [0 45];  out.Y = [-90 -45];    out.Iname = [handles.f_path '45S000E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45S045E.
function pushbutton_45S045E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45S045E.jpg']),1 );
out.X = [45 90];  out.Y = [-90 -45];    out.Iname = [handles.f_path '45S045E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45S090E.
function pushbutton_45S090E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45S090E.jpg']),1 );
out.X = [90 135];  out.Y = [-90 -45];    out.Iname = [handles.f_path '45S090E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);

% --- Executes on button press in pushbutton_45S135E.
function pushbutton_45S135E_Callback(hObject, eventdata, handles)
out.img = flipdim( imread([handles.f_path '45S135E.jpg']),1 );
out.X = [135 180];  out.Y = [-90 -45];    out.Iname = [handles.f_path '45S135E.jpg'];
handles.output = out;   guidata(hObject, handles);  uiresume(handles.figure1);


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

% --- Executes on button press in radiobutton_MapTiles.
function radiobutton_MapTiles_Callback(hObject, eventdata, handles)
if ~(get(hObject,'Value'))
    set(handles.radiobutton_MapTiles, 'Value', 1);
    set(handles.radiobutton_WorldMap, 'Value', 0);
    return
end
set(handles.radiobutton_WorldMap, 'Value', 0);
set(gcf,'Name','World Topo Tiles')
set(handles.h_ll,'Position',[1 1 64 64],'TooltipString','45S180W')
set_tiles(handles)
set(handles.h_all,'Visible','on')
%guidata(hObject, handles);

% --- Executes on button press in radiobutton_WorldMap.
function radiobutton_WorldMap_Callback(hObject, eventdata, handles)
set(handles.radiobutton_WorldMap, 'Value', 1);
set(handles.radiobutton_MapTiles, 'Value', 0);
set(gcf,'Name','World Topo')
if (handles.first_WT)           % Don't need to do this all times (if the user plays)
    h_all = findall(gcf,'Style','pushbutton');
    zz=get(h_all,'Position');   xx = [];
    for i=1:32
        xx = [xx; zz{i}];
    end
    x_min = min(xx(:,1));   x_max = max(xx(:,1));
    y_min = min(xx(:,2));   y_max = max(xx(:,2));
    x_max = x_max + 65;     y_max = y_max + 65;
    handles.h_ll = findobj(h_all,'Tag','pushbutton_45S180W');
    handles.first_WT = 0;   handles.h_all = h_all;
    handles.pushbutonWT_pos = [x_min y_min x_max y_max];
    handles.logo = imread([handles.f_path 'logo_etopo2.jpg']);
    guidata(hObject, handles);
end
set(handles.h_all,'Visible','off')
set(handles.h_ll,'Position',handles.pushbutonWT_pos,'CData',handles.logo,'Visible','on','TooltipString','')

% --------------------------------------------------------------------
function set_tiles(handles)
load([handles.f_path 'WorldMapTilled.mat'])
set(handles.pushbutton_90N180W,'CData',i90N180W);   set(handles.pushbutton_90N135W,'CData',i90N135W)
set(handles.pushbutton_90N090W,'CData',i90N090W);   set(handles.pushbutton_90N045W,'CData',i90N045W)
set(handles.pushbutton_90N000E,'CData',i90N000E);   set(handles.pushbutton_90N045E,'CData',i90N045E)
set(handles.pushbutton_90N090E,'CData',i90N090E);   set(handles.pushbutton_90N135E,'CData',i90N135E)

set(handles.pushbutton_45N180W,'CData',i45N180W);   set(handles.pushbutton_45N135W,'CData',i45N135W)
set(handles.pushbutton_45N090W,'CData',i45N090W);   set(handles.pushbutton_45N045W,'CData',i45N045W)
set(handles.pushbutton_45N000E,'CData',i45N000E);   set(handles.pushbutton_45N045E,'CData',i45N045E)
set(handles.pushbutton_45N090E,'CData',i45N090E);   set(handles.pushbutton_45N135E,'CData',i45N135E)

set(handles.pushbutton_00N180W,'CData',i00N180W);   set(handles.pushbutton_00N135W,'CData',i00N135W)
set(handles.pushbutton_00N090W,'CData',i00N090W);   set(handles.pushbutton_00N045W,'CData',i00N045W)
set(handles.pushbutton_00N000E,'CData',i00N000E);   set(handles.pushbutton_00N045E,'CData',i00N045E)
set(handles.pushbutton_00N090E,'CData',i00N090E);   set(handles.pushbutton_00N135E,'CData',i00N135E)

set(handles.pushbutton_45S180W,'CData',i45S180W);   set(handles.pushbutton_45S135W,'CData',i45S135W)
set(handles.pushbutton_45S090W,'CData',i45S090W);   set(handles.pushbutton_45S045W,'CData',i45S045W)
set(handles.pushbutton_45S000E,'CData',i45S000E);   set(handles.pushbutton_45S045E,'CData',i45S045E)
set(handles.pushbutton_45S090E,'CData',i45S090E);   set(handles.pushbutton_45S135E,'CData',i45S135E)


% --- Creates and returns a handle to the GUI figure. 
function bg_region_map_tilled_LayoutFcn(h1,handles);

set(h1,'PaperUnits','centimeters',...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','bg_region_map_tilled',...
'NumberTitle','off',...
'Position',[520 525 505 280],...
'Resize','off',...
'Tag','figure1');

h2 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_90N180W_Callback'},...
'Position',[1 190 64 64],...
'TooltipString','90N180W',...
'Tag','pushbutton_90N180W');

h3 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_90N135W_Callback'},...
'Position',[64 190 64 64],...
'TooltipString','90N135W',...
'Tag','pushbutton_90N135W');

h4 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_90N090W_Callback'},...
'Position',[127 190 64 64],...
'TooltipString','90N090W',...
'Tag','pushbutton_90N090W');

h5 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_90N045W_Callback'},...
'Position',[190 190 64 64],...
'TooltipString','90N045W',...
'Tag','pushbutton_90N045W');

h6 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_90N000E_Callback'},...
'Position',[253 190 64 64],...
'TooltipString','90N000E',...
'Tag','pushbutton_90N000E');

h7 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_90N045E_Callback'},...
'Position',[316 190 64 64],...
'TooltipString','90N045E',...
'Tag','pushbutton_90N045E');

h8 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_90N090E_Callback'},...
'Position',[379 190 64 64],...
'TooltipString','90N090E',...
'Tag','pushbutton_90N090E');

h9 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_90N135E_Callback'},...
'Position',[442 190 64 64],...
'TooltipString','90N135E',...
'Tag','pushbutton_90N135E');

h10 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45N180W_Callback'},...
'Position',[1 127 64 64],...
'TooltipString','45N180W',...
'Tag','pushbutton_45N180W');

h11 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45N135W_Callback'},...
'Position',[64 127 64 64],...
'TooltipString','45N135W',...
'Tag','pushbutton_45N135W');

h12 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45N090W_Callback'},...
'Position',[127 127 64 64],...
'TooltipString','45N090W',...
'Tag','pushbutton_45N090W');

h13 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45N045W_Callback'},...
'Position',[190 127 64 64],...
'TooltipString','45N045W',...
'Tag','pushbutton_45N045W');

h14 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45N000E_Callback'},...
'Position',[253 127 64 64],...
'TooltipString','45N000E',...
'Tag','pushbutton_45N000E');

h15 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45N045E_Callback'},...
'Position',[316 127 64 64],...
'TooltipString','45N045E',...
'Tag','pushbutton_45N045E');

h16 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45N090E_Callback'},...
'Position',[379 127 64 64],...
'TooltipString','45N090E',...
'Tag','pushbutton_45N090E');

h17 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45N135E_Callback'},...
'Position',[442 127 64 64],...
'TooltipString','45N135E',...
'Tag','pushbutton_45N135E');

h18 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_00N180W_Callback'},...
'Position',[1 64 64 64],...
'TooltipString','00N180W',...
'Tag','pushbutton_00N180W');

h19 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_00N135W_Callback'},...
'Position',[64 64 64 64],...
'TooltipString','00N135W',...
'Tag','pushbutton_00N135W');

h20 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_00N090W_Callback'},...
'Position',[127 64 64 64],...
'TooltipString','00N090W',...
'Tag','pushbutton_00N090W');

h21 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_00N045W_Callback'},...
'Position',[190 64 64 64],...
'TooltipString','00N045W',...
'Tag','pushbutton_00N045W');

h22 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_00N000E_Callback'},...
'Position',[253 64 64 64],...
'TooltipString','00N000E',...
'Tag','pushbutton_00N000E');

h23 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_00N045E_Callback'},...
'Position',[316 64 64 64],...
'TooltipString','00N045E',...
'Tag','pushbutton_00N045E');

h24 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_00N090E_Callback'},...
'Position',[379 64 64 64],...
'TooltipString','00N090E',...
'Tag','pushbutton_00N090E');

h25 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_00N135E_Callback'},...
'Position',[442 64 64 64],...
'TooltipString','00N135E',...
'Tag','pushbutton_00N135E');

h26 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45S180W_Callback'},...
'Position',[1 1 64 64],...
'TooltipString','45S180W',...
'Tag','pushbutton_45S180W');

h27 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45S135W_Callback'},...
'Position',[64 1 64 64],...
'TooltipString','45S135W',...
'Tag','pushbutton_45S135W');

h28 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45S090W_Callback'},...
'Position',[127 1 64 64],...
'TooltipString','45S090W',...
'Tag','pushbutton_45S090W');

h29 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45S045W_Callback'},...
'Position',[190 1 64 64],...
'TooltipString','45S045W',...
'Tag','pushbutton_45S045W');

h30 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45S000E_Callback'},...
'Position',[253 1 64 64],...
'TooltipString','45S000E',...
'Tag','pushbutton_45S000E');

h31 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45S045E_Callback'},...
'Position',[316 1 64 64],...
'TooltipString','45S045E',...
'Tag','pushbutton_45S045E');

h32 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45S090E_Callback'},...
'Position',[379 1 64 64],...
'TooltipString','45S090E',...
'Tag','pushbutton_45S090E');

h33 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'pushbutton_45S135E_Callback'},...
'Position',[442 1 64 64],...
'TooltipString','45S135E',...
'Tag','pushbutton_45S135E');

h34 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'radiobutton_MapTiles_Callback'},...
'Position',[10 259 94 15],...
'String','World Map Tiles',...
'Style','radiobutton',...
'Value',1,...
'Tag','radiobutton_MapTiles');

h35 = uicontrol('Parent',h1,...
'Callback',{@bg_region_map_tilled_uicallback,h1,'radiobutton_WorldMap_Callback'},...
'Position',[128 259 79 15],...
'String','World Map',...
'Style','radiobutton',...
'Tag','radiobutton_WorldMap');

function bg_region_map_tilled_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

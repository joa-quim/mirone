function varargout = mirone_pref(varargin)
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
mirone_pref_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
 
movegui(hObject,'northwest');

global home_dir
if isempty(home_dir),   handles.d_path = [pwd filesep 'data' filesep];    % Case when this function was called directly
else                    handles.d_path = [home_dir filesep 'data' filesep];   end

directory_list = [];
gridNCformat = 1;           % Default to new netCDF format
load([handles.d_path 'mirone_pref.mat']);
handles.geog = geog;        % Just to not be empty.
handles.ForceInsitu = 0;    % Just to not be empty.
handles.gridNCformat = gridNCformat;

j = logical(zeros(1,length(directory_list)));           % vector for eventual cleaning non-existing dirs
if iscell(directory_list)                               % When exists a dir list in mirone_pref
    for (i = 1:length(directory_list))
        try,        cd(directory_list{i});              % NOTE. I don't use 'exist' anymore because
        catch,      j(i) = 1;       cd(home_dir);       % the stupid compiler allways return something > 0
        end
    end
    cd(home_dir);                           % Need to come back home because it was probably somewere out there
    directory_list(j) = [];                             % clean eventual non-existing directories
    if (~isempty(directory_list))                       % If there is one left
        set(handles.popup_directory_list,'String',directory_list)
        handles.last_directories = directory_list;
    else
        handles.last_directories = {[home_dir filesep 'tmp']; home_dir};    % Let it have something existent
        set(handles.popup_directory_list,'String',handles.last_directories)
    end
else                                                    % mirone_pref had no dir list
    handles.last_directories = {[home_dir filesep 'tmp']; home_dir};    % Let it have something existent
    set(handles.popup_directory_list,'String',handles.last_directories)
end

if (nargin >= 3)                    % That's the case when called from mirone
    if varargin{1} == 1             % Signals a geographic grid/image
        set(handles.radiobutton_geog,'Value',1)
        set(handles.radiobutton_cart,'Value',0)
        handles.geog = 1;
    elseif varargin{1} == 0
        set(handles.radiobutton_geog,'Value',0)
        set(handles.radiobutton_cart,'Value',1)
        handles.geog = 0;
    else
        handles.geog = 1;
    end
    set(handles.edit_GridMaxSize,'String',num2str(varargin{2}))
    set(handles.edit_swathRatio,'String',num2str(varargin{3}))
    set(handles.checkbox_ForceInsitu,'Value',varargin{4})
    handles.ForceInsitu = varargin{4};
end    

% Well this is split from the above because it was written later and I don't want to mess
% with what is working. Wrap in a try-catch because the first time the variables are not
% yet in mirone_pref.mat
try     % Goes here all other times
    set(handles.popupmenu_DefLineThickness,'String',DefLineThick)
    set(handles.popupmenu_DefLineColor,'String',DefLineColor)
    set(handles.popup_MeasureUnites,'String',DefineMeasureUnit)
    set(handles.popup_ellipsoide,'String',DefineEllipsoide)
    set(handles.checkbox_NewGrids,'Value',out_in_NewWindow)
    set(handles.checkbox_SaveAsInt16,'Value',saveAsInt16)
    set(handles.checkbox_newNCformat,'Value',gridNCformat)
catch       % Comes here in first call before variables are stored in mirone_pref.mat
    DefLineThick = {'2 pt'; '1 pt'; '3 pt'; '4 pt'};
    DefLineColor = {'White'; 'Black'; 'Red'; 'Green'; 'Blue'; 'Cyan'; 'Yellow'; 'Magenta'};
    DefineMeasureUnit = {'nautic miles'; 'kilometers'; 'meters'; 'user'};
    set(handles.popupmenu_DefLineThickness,'String',DefLineThick)
    set(handles.popupmenu_DefLineColor,'String',DefLineColor)
    set(handles.popup_MeasureUnites,'String',DefineMeasureUnit)
    set(handles.checkbox_NewGrids,'Value',1)
end

% This is the default ellipsoide order. It will be changed (and saved as so in mirone_pref) by the user
% Note, however, that only the ellipsoide names are stored in mirone_pref. The parameters of the selected
% ellipsoid are found by a "case" loop runned by the OK pushbutton.
handles.ellipsoide = {'WGS-84 - 1984', 6378137.0, 0.0, 1.0/298.2572235630;
		'OSU91A - 1991', 6378136.3, 0.0, 1.0/298.25722;
		'OSU86F - 1986', 6378136.2, 0.0, 1.0/298.25722;
		'Engelis - 1985', 6378136.05, 0.0, 1.0/298.2566;
		'SGS-85 - 1985', 6378136.0, 0.0, 1.0/298.257;
		'MERIT-83 - 1983', 6378137.0, 0.0, 1.0/298.257;
		'GRS-80 - 1980', 6378137.0, 0.0, 1.0/298.257222101;
		'Lerch - 1979', 6378139.0, 0.0, 1.0/298.257;
		'ATS77 - 1977', 6378135.0, 0.0, 1.0/298.257;
		'IAG-75 - 1975', 6378140.0, 0.0, 1.0/298.257222;
		'Indonesian - 1974', 6378160.0, 0.0, 1.0/298.247;
		'WGS-72 - 1972', 6378135.0, 0.0, 1.0/298.26;
		'NWL-10D - 1972', 6378135.0, 0.0, 1.0/298.26;
		'South-American - 1969', 6378160.0, 0.0, 1.0/298.25;
		'Fischer-1968', 6378150.0, 0.0, 1.0/298.3;
		'Modified-Mercury-1968', 6378150.0, 0.0, 1.0/298.3;
		'GRS-67 - 1967', 6378160.0, 0.0, 1.0/298.247167427;
		'International-1967', 6378157.5, 0.0, 1.0/298.25;
		'WGS-66 - 1966', 6378145.0, 0.0, 1.0/298.25;
		'NWL-9D - 1966', 6378145.0, 0.0, 1.0/298.25;
		'Australian - 1965', 6378160.0, 0.0, 1.0/298.25;
		'APL4.9 - 1965', 6378137.0, 0.0, 1.0/298.25;
		'Kaula - 1961', 6378163.0, 0.0, 1.0/298.24;
		'Hough - 1960', 6378270.0, 0.0, 1.0/297.0;
		'WGS-60 - 1960', 6378165.0, 0.0, 1.0/298.3;
		'Fischer-1960', 6378166.0, 0.0, 1.0/298.3;
		'Mercury-1960', 6378166.0, 0.0, 1.0/298.3;
		'Modified-Fischer-1960', 6378155.0, 0.0, 1.0/298.3;
		'Fischer-1960-SouthAsia', 6378155.0, 0.0, 1.0/298.3;
		'Krassovsky - 1940', 6378245.0, 0.0, 1.0/298.3;
		'War-Office - 1926', 6378300.583, 0.0, 1.0/296.0;
		'International-1924', 6378388.0, 0.0, 1.0/297.0;
		'Hayford-1909', 6378388.0, 0.0, 1.0/297.0;
		'Helmert-1906', 6378200.0, 0.0, 1.0/298.3;
		'Clarke-1880', 6378249.145, 0.0, 1.0/293.465;
		'Clarke-1880-Arc1950 - 1880', 6378249.145326, 0.0, 1.0/293.4663076;
		'Clarke-1880-IGN', 6378249.2, 0.0, 1.0/293.4660213;
		'Clarke-1880-Jamaica', 6378249.136, 0.0, 1.0/293.46631;
		'Clarke-1880-Merchich', 6378249.2, 0.0, 1.0/293.46598;
		'Clarke-1880-Palestine', 6378300.79, 0.0, 1.0/293.46623;
		'Andrae - 1876', 6377104.43, 0.0, 1.0/300.0;
		'Clarke-1866', 6378206.4, 0.0, 1.0/294.9786982;
		'Clarke-1866-Michigan', 6378450.047484481, 0.0, 1.0/294.9786982;
		'Struve - 1860', 6378297.0, 0.0, 1.0/294.73;
		'Clarke-1858', 6378293.639, 0.0, 1.0/294.26068;
		'Airy - 1830', 6377563.396, 0.0, 1.0/299.3249646;
		'Airy-Ireland - 1830', 6377340.189, 0.0, 1.0/299.3249646;
		'Modified-Airy - 1830', 6377340.189, 0.0, 1.0/299.3249646;
		'Bessel - 1841', 6377397.155, 0.0, 1.0/299.1528128;
		'Bessel-Schwazeck - 1841', 6377483.865, 0.0, 1.0/299.1528128;
		'Bessel-Namibia - 1841', 6377483.865, 0.0, 1.0/299.1528128;
		'Bessel-NGO1948 - 1841', 6377492.0176, 0.0, 1.0/299.15281;
		'Everest-1830', 6377276.345, 0.0, 1.0/300.8017;
		'Everest-1830-Kalianpur', 6377301.243, 0.0, 1.0/300.80174;
		'Everest-1830-Kertau', 6377304.063, 0.0, 1.0/300.8017;
		'Everest-1830-Timbalai', 6377298.556, 0.0, 1.0/300.8017;
		'Everest-1830-Pakistan', 6377309.613, 0.0, 1.0/300.8017;
		'Walbeck - 1819', 6376896.0, 0.0, 1.0/302.78;
		'Plessis - 1817', 6376523.0, 0.0, 1.0/308.64;
		'Delambre - 1810', 6376428.0, 0.0, 1.0/311.5;
		'CPM - 1799', 6375738.7, 0.0, 1.0/334.29;
		'Maupertius - 1738', 6397300.0, 0.0, 1.0/191.0;
		'Sphere - 1980', 6371008.7714, 0.0, 0.0};

try
	if (handles.geog == 0)      % For cartesian coords the following is no applyable
        set(handles.popup_ellipsoide,'Enable','off')
        set(handles.popup_MeasureUnites,'Enable','off')
	end
catch   % In case of error, set the default list.
    set(handles.popup_ellipsoide,'String',handles.ellipsoide(:,1))
    set(handles.popup_MeasureUnites,'String',{'nautic miles'; 'kilometers'; 'meters'; 'user'})
end

% Create the "ForceInsitu" TooltipString
h = findobj(hObject,'Tag','checkbox_ForceInsitu','Style','checkbox');
str = sprintf(['Importing grids implies a conversion that uses\n'...
    'matrix transposition. This operations is fast if\n'...
    'we make a copy of the importing grid. However,\n'...
    'this requires twice the grid size on memory.\n'...
    'If you don''t have enough memory to import a\n'...
    'large grid, use this option that do the transposition\n'...
    '"insitu". That is, it uses only one time the grid\n'...
    'size in memory. The price you will pay, however,\n'...
    'is in speed because it runs about 10 times slower.']);
set(h,'TooltipString',str)

% Choose default command line output for mirone_pref_export
handles.output = hObject;
guidata(hObject, handles);
set(hObject,'Visible','on');
% UIWAIT makes mirone_pref_export wait for user response (see UIRESUME)
uiwait(handles.figure1);

handles = guidata(hObject);
out = mirone_pref_OutputFcn(hObject, [], handles);
varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = mirone_pref_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure

% Get default command line output from handles structure
varargout{1} = handles.output;
% The figure can be deleted now
delete(handles.figure1);

% ------------------------------------------------------------------------------------
function checkbox_newNCformat_Callback(hObject, eventdata, handles)
if get(hObject,'Value')     handles.gridNCformat = 1;
else                        handles.gridNCformat = 0;   end
guidata(hObject,handles)
    
% ------------------------------------------------------------------------------------
function radiobutton_geog_Callback(hObject, eventdata, handles)
if get(hObject,'Value')     handles.geog = 1;   set(handles.radiobutton_cart,'Value',0)
else                        handles.geog = 0;   set(handles.radiobutton_cart,'Value',1); end
guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function radiobutton_cart_Callback(hObject, eventdata, handles)
if get(hObject,'Value')     handles.geog = 0;   set(handles.radiobutton_geog,'Value',0)
else                        handles.geog = 1;   set(handles.radiobutton_geog,'Value',1); end
guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function edit_GridMaxSize_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if isnan(str2double(xx)) | isempty(xx)    % Just a stupid user error
    set(hObject, 'String', '20');    return
else
    if (str2double(xx) >= 10)   % Numbers bigger than 10 Mb don't need decimal places
        set(hObject,'String',num2str(fix(str2double(xx))))
    end
end

% ------------------------------------------------------------------------------------
function edit_swathRatio_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if isnan(str2double(xx)) | isempty(xx)    % Just a stupid user error
    set(hObject, 'String', '3');    return
end

% ------------------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
handles.output = [];        % User gave up, return nothing
guidata(hObject, handles);  uiresume(handles.figure1);

% ------------------------------------------------------------------------------------
function popup_directory_list_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
% Put the selected field on top of the String list. This is necessary because the "OK" button will
% read this list and save it in mirone_pref_export, so next time the selected field will show up first.
tmp = str(val);         str(val) = [];
new_str = [tmp; str];   set(hObject,'String',new_str); 
set(hObject,'Value',1)

% ------------------------------------------------------------------------------------
function pushbutton_change_dir_Callback(hObject, eventdata, handles)
contents = get(handles.popup_directory_list,'String');
if (ispc)
    work_dir = uigetfolder_standalone('Select a directory',contents{get(handles.popup_directory_list,'Value')});
else            % This guy says it cannot be compiled
    work_dir = uigetdir;
end

if (isempty(work_dir) | isequal(work_dir,0))    return;     end

handles.last_directories = [cellstr(work_dir); handles.last_directories];
if length(handles.last_directories) > 15            % Keep only 15 adresses
    handles.last_directories(16:end) = [];
end
set(handles.popup_directory_list,'String',handles.last_directories)
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function checkbox_NewGrids_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    handles.out_in_NewWindow = 1;
else
    handles.out_in_NewWindow = 0;
end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function checkbox_ForceInsitu_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    handles.ForceInsitu = 1;
else
    handles.ForceInsitu = 0;
end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function checkbox_SaveAsInt16_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    handles.saveAsInt16 = 1;
else
    handles.saveAsInt16 = 0;
end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function popupmenu_DefLineThickness_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
% Put the selected field on top of the String list. This is necessary because the "OK" button will
% read this list and save it in mirone_pref, so next time the selected field will show up first.
tmp = str(val);         str(val) = [];
new_str = [tmp; str];   set(hObject,'String',new_str); 
set(hObject,'Value',1)

% ------------------------------------------------------------------------------------
function popupmenu_DefLineColor_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
% Put the selected field on top of the String list. This is necessary because the "OK" button will
% read this list and save it in mirone_pref, so next time the selected field will show up first.
tmp = str(val);         str(val) = [];
new_str = [tmp; str];   set(hObject,'String',new_str); 
set(hObject,'Value',1)

% ------------------------------------------------------------------------------------
function popup_MeasureUnites_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
% Put the selected field on top of the String list. This is necessary because the "OK" button will
% read this list and save it in mirone_pref, so next time the selected field will show up first.
tmp = str(val);         str(val) = [];
new_str = [tmp; str];   set(hObject,'String',new_str); 
set(hObject,'Value',1)

% ------------------------------------------------------------------------------------
function popup_ellipsoide_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
% Put the selected field on top of the String list. This is necessary because the "OK" button will
% read this list and save it in mirone_pref, so next time the selected field will show up first.
tmp = str(val);         str(val) = [];
new_str = [tmp; str];   set(hObject,'String',new_str); 
set(hObject,'Value',1)

% ------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
Out.geog = handles.geog;
Out.grdMaxSize = str2double(get(handles.edit_GridMaxSize,'String'));
Out.swathRatio = str2double(get(handles.edit_swathRatio,'String'));
directory_list = get(handles.popup_directory_list, 'String');
Out.last_dir   = directory_list{1};
DefLineThick = get(handles.popupmenu_DefLineThickness, 'String');
DefLineColor = get(handles.popupmenu_DefLineColor, 'String');
Out.out_in_NewWindow = get(handles.checkbox_NewGrids,'Value');
Out.saveAsInt16 = get(handles.checkbox_SaveAsInt16,'Value');
Out.ForceInsitu = handles.ForceInsitu;
Out.gridNCformat = handles.gridNCformat;
% Decode the line thickness string into a number
Out.DefLineThick = str2num(DefLineThick{1}(1));
% Decode the line color string into the corresponding char (e.g. k,w, etc...)
switch DefLineColor{1}
    case 'Black'
        Out.DefLineColor = 'k';
    case 'White'
        Out.DefLineColor = 'w';
    case 'Red'
        Out.DefLineColor = 'r';
    case 'Green'
        Out.DefLineColor = 'g';
    case 'Blue'
        Out.DefLineColor = 'b';
    case 'Cyan'
        Out.DefLineColor = 'c';
    case 'Yellow'
        Out.DefLineColor = 'y';
    case 'Magenta'
        Out.DefLineColor = 'm';
end

% Decode the Measure units into a char code (e.g n, k, m, u)
DefineMeasureUnit = get(handles.popup_MeasureUnites, 'String');
switch DefineMeasureUnit{1}
    case 'nautic miles'
        Out.DefineMeasureUnit = 'n';
    case 'kilometers'
        Out.DefineMeasureUnit = 'k';
    case 'meters'
        Out.DefineMeasureUnit = 'm';
    case 'user'
        Out.DefineMeasureUnit = 'u';
end

% Decode the Ellipsoide into a var containg a,b,f
DefineEllipsoide = get(handles.popup_ellipsoide, 'String');
if (handles.geog == 1)
    for i=1:length(handles.ellipsoide)
        switch DefineEllipsoide{1}
            case handles.ellipsoide(i)
                Out.DefineEllipsoide(1) = handles.ellipsoide{i,2};
                Out.DefineEllipsoide(2) = handles.ellipsoide{i,3};
                Out.DefineEllipsoide(3) = handles.ellipsoide{i,4};
        end
    end
else        % For the time beeing default to WGS-84
    Out.DefineEllipsoide(1) = handles.ellipsoide{1,2};
    Out.DefineEllipsoide(2) = handles.ellipsoide{1,3};
    Out.DefineEllipsoide(3) = handles.ellipsoide{1,4};
end

fname = [handles.d_path 'mirone_pref.mat'];
% Save the preferences to a mat file under the data directory
% Note: for the ellipsoide we save it's parameters (a,b,f) instead of the name
DefineEllipsoide_params = Out.DefineEllipsoide;    % For saving purposes
geog = Out.geog;      grdMaxSize = Out.grdMaxSize;    swathRatio = Out.swathRatio;
out_in_NewWindow = Out.out_in_NewWindow;
saveAsInt16 = Out.saveAsInt16;
gridNCformat = handles.gridNCformat;
%ForceInsitu = handles.ForceInsitu;     % We don't save it because the user must choose it every time
save(fname,'geog','grdMaxSize','swathRatio','directory_list','DefLineThick','DefLineColor',...
    'DefineMeasureUnit','DefineEllipsoide','DefineEllipsoide_params', 'out_in_NewWindow',...
    'saveAsInt16', 'gridNCformat', '-append')
handles.output = Out;           guidata(hObject,handles)
uiresume(handles.figure1);

% ------------------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    handles.output = [];        % User gave up, return nothing
    guidata(hObject, handles);    uiresume(handles.figure1);
else    % The GUI is no longer waiting, just close it
    handles.output = [];        % User gave up, return nothing
    guidata(hObject, handles);    delete(handles.figure1);
end

% ------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = [];        % User said no by hitting escape
    guidata(hObject, handles);    uiresume(handles.figure1);
end

% --- Creates and returns a handle to the GUI figure. 
function mirone_pref_LayoutFcn(h1,handles);

set(h1, 'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','mirone_pref',...
'NumberTitle','off',...
'Position',[520 437 272 363],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[10 335 261 2],'String',{''},...
'Style','frame','Tag','frame2');

h2 = uicontrol('Parent',h1,'Position',[10 269 111 50],'String',{''},...
'Style','frame','Tag','frame1');

h3 = uicontrol('Parent',h1,...
'Callback',{@mirone_pref_uicallback,h1,'radiobutton_geog_Callback'},...
'Position',[20 293 79 15],...
'String','Geographic',...
'Style','radiobutton',...
'TooltipString','GMT grid is in geographical coordinates',...
'Value',1,...
'Tag','radiobutton_geog');

h4 = uicontrol('Parent',h1,...
'Callback',{@mirone_pref_uicallback,h1,'radiobutton_cart_Callback'},...
'Position',[20 274 79 15],...
'String','Cartesian',...
'Style','radiobutton',...
'TooltipString','GMT grid is in cartesian coordinates',...
'Tag','radiobutton_cart');

h5 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'edit_GridMaxSize_Callback'},...
'HorizontalAlignment','left',...
'Position',[134 297 36 20],...
'String','20',...
'Style','edit',...
'TooltipString','Grid max size that will be stored in memory',...
'Tag','edit_GridMaxSize');

h6 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'edit_swathRatio_Callback'},...
'HorizontalAlignment','left',...
'Position',[134 272 36 20],...
'String','3',...
'Style','edit',...
'TooltipString','Swath Ratio for multibeam planing',...
'Tag','edit_swathRatio');

h7 = uicontrol('Parent',h1,'Position',[21 311 90 15],...
'String','Grid coordinates','Style','text','Tag','text1');

h8 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[174 300 91 15],...
'String','Grid max size (Mb)','Style','text','Tag','text2');

h9 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[174 275 71 15],...
'String','Swath ratio','Style','text','Tag','text3');

h10 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'popup_MeasureUnites_Callback'},...
'Position',[10 217 101 22],...
'String',{  'nautic miles'; 'kilometers'; 'meters'; 'user' },...
'Style','popupmenu',...
'TooltipString','Select the default measure units',...
'Value',1,...
'Tag','popup_MeasureUnites');

h11 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'popup_ellipsoide_Callback'},...
'Position',[120 217 145 22],...
'String','WGS84',...
'Style','popupmenu',...
'TooltipString','Select the default ellipsoide',...
'Value',1,...
'Tag','popup_ellipsoide');

h12 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'popup_directory_list_Callback'},...
'Position',[10 166 236 22],...
'Style','popupmenu',...
'TooltipString','Select the default initial directory from list',...
'Value',1,...
'Tag','popup_directory_list');

h13 = uicontrol('Parent',h1,...
'Callback',{@mirone_pref_uicallback,h1,'pushbutton_change_dir_Callback'},...
'FontSize',10,...
'FontWeight','bold',...
'Position',[246 167 18 21],...
'String','...',...
'TooltipString','Select a different directory',...
'Tag','pushbutton_change_dir');

h14 = uicontrol('Parent',h1,...
'Callback',{@mirone_pref_uicallback,h1,'checkbox_NewGrids_Callback'},...
'Position',[10 140 147 15],...
'String','New grids in new window',...
'Style','checkbox',...
'TooltipString','If not checked new grids are written on disk',...
'Value',1,...
'Tag','checkbox_NewGrids');

h15 = uicontrol('Parent',h1,...
'Callback',{@mirone_pref_uicallback,h1,'checkbox_ForceInsitu_Callback'},...
'Position',[10 115 147 15],...
'String','Force "Insitu" transposition',...
'Style','checkbox',...
'Tag','checkbox_ForceInsitu');

h16 = uicontrol('Parent',h1,...
'Callback',{@mirone_pref_uicallback,h1,'checkbox_SaveAsInt16_Callback'},...
'Position',[10 90 160 15],...
'String','When possible, save as Int16',...
'Style','checkbox',...
'TooltipString','Grids that were originaly on short int format (2 bytes) will be saved as so',...
'Tag','checkbox_SaveAsInt16');

h17 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'popupmenu_DefLineThickness_Callback'},...
'Position',[10 41 100 22],...
'String',{'1 pt'; '2 pt'; '3 pt' },...
'Style','popupmenu',...
'TooltipString','All drawn lines will have this thickness',...
'Value',1,...
'Tag','popupmenu_DefLineThickness');

h18 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'popupmenu_DefLineColor_Callback'},...
'Position',[144 41 100 22],...
'String',{  'Black'; 'White'; 'Red'; 'Green'; 'Blue'; 'Cyan'; 'Yellow'; 'Magenta' },...
'Style','popupmenu',...
'TooltipString','All drawn lines will have this color',...
'Value',1,...
'Tag','popupmenu_DefLineColor');

h19 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[10 188 80 15],...
'String','Default directory','Style','text','Tag','text4');

h20 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[10 63 102 15],...
'String','Default line thickness','Style','text','Tag','text5');

h21 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[145 63 86 15],...
'String','Default line color','Style','text','Tag','text6');

h22 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[10 240 80 15],...
'String','Measure units','Style','text','Tag','text7');

h23 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[121 240 90 15],...
'String','Default ellipsoid','Style','text','Tag','text8');

h24 = uicontrol('Parent',h1,...
'Callback',{@mirone_pref_uicallback,h1,'pushbutton_OK_Callback'},...
'Position',[120 10 66 23],'String','OK','Tag','pushbutton_OK');

h25 = uicontrol('Parent',h1,...
'Callback',{@mirone_pref_uicallback,h1,'pushbutton_cancel_Callback'},...
'Position',[198 10 66 23],'String','Cancel','Tag','pushbutton_cancel');

h26 = uicontrol('Parent',h1,...
'Callback',{@mirone_pref_uicallback,h1,'checkbox_newNCformat_Callback'},...
'Position',[10 342 191 15],...
'String','Use new GMT netCDF grid format',...
'Style','checkbox',...
'TooltipString','Use the new GMT 4.1 grid format as default for writting grid files',...
'Tag','checkbox_newNCformat');

function mirone_pref_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

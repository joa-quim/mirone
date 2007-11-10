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
	mirone_pref_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'northwest');

    handMir = varargin{1};
    home_dir = handMir.home_dir;
    handles.d_path = handMir.path_data;
    handles.handMir = handMir;

	directory_list = [];
	load([handMir.path_data 'mirone_pref.mat']);
	handles.geog = geog;		% Just to not be empty.
	handles.ForceInsitu = 0;	% 		"
	handles.flederPlanar = 1;	%		"
	handles.flederBurn = 1;
	handles.whichFleder = 1;

	% The next are new (20-1-07) and therefore we need to wrap it in try because old prefs do not have it yet
	try
        handles.flederPlanar = flederPlanar;
        set(handles.radio_planar,'Val',handles.flederPlanar)
        set(handles.radio_spherical,'Val',~handles.flederPlanar)
        handles.flederBurn = flederBurn;
        if (handles.flederBurn == 0),       set(handles.radio_noBurnAtAll,'Val',1); set(handles.radio_coastsOnly,'Val',0)
        elseif (handles.flederBurn == 1),   set(handles.radio_coastsOnly,'Val',1);
        else                                set(handles.radio_burnAll,'Val',1);     set(handles.radio_coastsOnly,'Val',0)
        end
		handles.whichFleder = whichFleder;	% whichFleder = 1 for the free iview3d or 0 for the true thing (fledermaus)
	end

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

    if (handMir.geog)             % Signals a geographic grid/image
        set(handles.radiobutton_geog,'Value',1)
        set(handles.radiobutton_cart,'Value',0)
        handles.geog = 1;
    elseif (handMir.geog == 0)
        set(handles.radiobutton_geog,'Value',0)
        set(handles.radiobutton_cart,'Value',1)
        handles.geog = 0;
    else
        handles.geog = 1;
    end
    set(handles.edit_GridMaxSize,'String',sprintf('%d',fix(handMir.grdMaxSize / (2^20))))
    set(handles.edit_swathRatio,'String',sprintf('%g',handMir.swathRatio))
    set(handles.checkbox_ForceInsitu,'Value',handMir.ForceInsitu)
    handles.ForceInsitu = handMir.ForceInsitu;

	% Well this is split from the above because it was written later and I don't want to mess
	% with what is working. Wrap in a try-catch because the first time the variables are not
	% yet in mirone_pref.mat
	try     % Goes here all other times
        set(handles.popupmenu_DefLineThickness,'String',DefLineThick)
        set(handles.popupmenu_DefLineColor,'String',DefLineColor)
        set(handles.popup_MeasureUnites,'String',DefineMeasureUnit)
        set(handles.popup_ellipsoide,'String',DefineEllipsoide)
        set(handles.checkbox_meanLat,'Value',scale2meanLat)
	catch       % Comes here in first call before variables are stored in mirone_pref.mat
        DefLineThick = {'2 pt'; '1 pt'; '3 pt'; '4 pt'};
        DefLineColor = {'White'; 'Black'; 'Red'; 'Green'; 'Blue'; 'Cyan'; 'Yellow'; 'Magenta'};
        DefineMeasureUnit = {'nautic miles'; 'kilometers'; 'meters'; 'user'};
        set(handles.popupmenu_DefLineThickness,'String',DefLineThick)
        set(handles.popupmenu_DefLineColor,'String',DefLineColor)
        set(handles.popup_MeasureUnites,'String',DefineMeasureUnit)
        set(handles.checkbox_meanLat,'Value',1)
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
str = sprintf(['Importing grids implies a conversion that uses\n'...
    'matrix transposition. This operations is fast if\n'...
    'we make a copy of the importing grid. However,\n'...
    'this requires twice the grid size on memory.\n'...
    'If you don''t have enough memory to import a\n'...
    'large grid, use this option that do the transposition\n'...
    '"insitu". That is, it uses only one time the grid\n'...
    'size in memory. The price you will pay, however,\n'...
    'is in speed because it runs about 10 times slower.']);
set(handles.checkbox_ForceInsitu,'TooltipString',str)

str = sprintf(['Since 1 degree of longitude and latitude do not cover the same\n'...
        'arc length at Earth surface, isometric plotting of geographical\n'...
        'grids squeezes the image vertically. Scaling the image to the\n'...
        'cosinus of the mean lat minimizes this effect.']);
set(handles.checkbox_meanLat,'TooltipString',str)

% ------------------ TABPANEL SECTION ----------------------------------------
% This is the tag that all tab push buttons share.  If you have multiple
% sets of tab push buttons, each group should have unique tag.
group_name = 'tab_group';

% This is a list of the UserData values used to link tab push buttons and
% the components on their linked panels.  To add a new tab panel to the group
%  Add the button using GUIDE
%  Assign the Tag based on the group name - in this case tab_group
%  Give the UserData a unique name - e.g. another_tab_panel
%  Add components to GUIDE for the new panel
%  Give the new components the same UserData as the tab button
%  Add the new UserData name to the below cell array
panel_names = {'general','fleder'};

% tabpanelfcn('makegroups',...) adds new fields to the handles structure,
% one for each panel name and another called 'group_name_all'.  These fields
% are used by the tabpanefcn when tab_group_handler is called.
handles = tabpanelfcn('make_groups',group_name, panel_names, handles, 1);
% ------------------------------------------------------------------------------

% Choose default command line output for mirone_pref
guidata(hObject, handles);
set(hObject,'Visible','on');

if (nargout)
    varargout{1} = hObject;
end

% -------------------------------------------------------------------------------------
function tab_group_ButtonDownFcn(hObject, eventdata, handles)
% Call the tab_group_handler.  This updates visiblity of components as needed to
% hide the components from the previous tab and show components on this tab.
% This also updates the last_tab field in the handles structure to keep track
% of which panel was hidden.
    handles = tabpanelfcn('tab_group_handler',hObject, handles, get(hObject, 'Tag'));
    guidata(hObject, handles);
    
% ------------------------------------------------------------------------------------
function radiobutton_geog_Callback(hObject, eventdata, handles)
	if get(hObject,'Value')     handles.geog = 1;   set(handles.radiobutton_cart,'Value',0)
	else                        handles.geog = 0;   set(handles.radiobutton_cart,'Value',1);
	end
	guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function radiobutton_cart_Callback(hObject, eventdata, handles)
	if get(hObject,'Value')     handles.geog = 0;   set(handles.radiobutton_geog,'Value',0)
	else                        handles.geog = 1;   set(handles.radiobutton_geog,'Value',1);
	end
	guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function edit_GridMaxSize_Callback(hObject, eventdata, handles)
	xx = get(hObject,'String');
	if isnan(str2double(xx)) || isempty(xx)    % Just a stupid user error
        set(hObject, 'String', '50')
	else
        if (str2double(xx) >= 10)   % Numbers bigger than 10 Mb don't need decimal places
            set(hObject,'String',sprintf('%d',fix(str2double(xx))))
        end
	end

% ------------------------------------------------------------------------------------
function edit_swathRatio_Callback(hObject, eventdata, handles)
	xx = get(hObject,'String');
	if isnan(str2double(xx)) || isempty(xx)    % Just a stupid user error
        set(hObject, 'String', '3')
	end

% ------------------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
    delete(handles.figure1)

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

if (isempty(work_dir) || isequal(work_dir,0))    return;     end

handles.last_directories = [cellstr(work_dir); handles.last_directories];
if length(handles.last_directories) > 15            % Keep only 15 adresses
    handles.last_directories(16:end) = [];
end
set(handles.popup_directory_list,'String',handles.last_directories)
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function checkbox_meanLat_Callback(hObject, eventdata, handles)
	if (get(hObject,'Value')),      handles.scale2meanLat = 1;
	else                            handles.scale2meanLat = 0;
	end
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function checkbox_ForceInsitu_Callback(hObject, eventdata, handles)
	if (get(hObject,'Value')),      handles.ForceInsitu = 1;
	else                            handles.ForceInsitu = 0;
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
	handles.handMir.geog = handles.geog;
	handles.handMir.grdMaxSize = str2double(get(handles.edit_GridMaxSize,'String')) * 2^20;
	handles.handMir.swathRatio = str2double(get(handles.edit_swathRatio,'String'));
	directory_list = get(handles.popup_directory_list, 'String');
	handles.handMir.last_dir   = directory_list{1};
	handles.handMir.work_dir   = handles.handMir.last_dir;
	DefLineThick = get(handles.popupmenu_DefLineThickness, 'String');
	DefLineColor = get(handles.popupmenu_DefLineColor, 'String');
	handles.handMir.scale2meanLat = get(handles.checkbox_meanLat,'Value');
	handles.handMir.ForceInsitu = handles.ForceInsitu;
	handles.handMir.flederPlanar = handles.flederPlanar;		flederPlanar = handles.flederPlanar;
	handles.handMir.flederBurn = handles.flederBurn;			flederBurn = handles.flederBurn;
	handles.handMir.whichFleder = handles.whichFleder;			whichFleder = handles.whichFleder;
	% Decode the line thickness string into a number
	handles.handMir.DefLineThick = str2num(DefLineThick{1}(1));
	% Decode the line color string into the corresponding char (e.g. k,w, etc...)
	switch DefLineColor{1}
		case 'Black',       handles.handMir.DefLineColor = 'k';
		case 'White',       handles.handMir.DefLineColor = 'w';
		case 'Red',         handles.handMir.DefLineColor = 'r';
		case 'Green',       handles.handMir.DefLineColor = 'g';
		case 'Blue',        handles.handMir.DefLineColor = 'b';
		case 'Cyan',        handles.handMir.DefLineColor = 'c';
		case 'Yellow',      handles.handMir.DefLineColor = 'y';
		case 'Magenta',     handles.handMir.DefLineColor = 'm';
	end

	% Decode the Measure units into a char code (e.g n, k, m, u)
	DefineMeasureUnit = get(handles.popup_MeasureUnites, 'String');
	switch DefineMeasureUnit{1}
		case 'nautic miles',    handles.handMir.DefineMeasureUnit = 'n';
		case 'kilometers',      handles.handMir.DefineMeasureUnit = 'k';
		case 'meters',          handles.handMir.DefineMeasureUnit = 'm';
		case 'user',            handles.handMir.DefineMeasureUnit = 'u';
	end

	% Decode the Ellipsoide into a var containg a,b,f
	DefineEllipsoide = get(handles.popup_ellipsoide, 'String');
	if (handles.geog == 1)
		for i=1:length(handles.ellipsoide)
			switch DefineEllipsoide{1}
				case handles.ellipsoide(i)
					handles.handMir.DefineEllipsoide(1) = handles.ellipsoide{i,2};
					handles.handMir.DefineEllipsoide(2) = handles.ellipsoide{i,3};
					handles.handMir.DefineEllipsoide(3) = handles.ellipsoide{i,4};
			end
		end
	else        % For the time beeing default to WGS-84
		handles.handMir.DefineEllipsoide(1) = handles.ellipsoide{1,2};
		handles.handMir.DefineEllipsoide(2) = handles.ellipsoide{1,3};
		handles.handMir.DefineEllipsoide(3) = handles.ellipsoide{1,4};
	end

	fname = [handles.d_path 'mirone_pref.mat'];
	% Save the preferences to a mat file under the data directory
	% Note: for the ellipsoide we save it's parameters (a,b,f) instead of the name
	DefineEllipsoide_params = handles.handMir.DefineEllipsoide;    % For saving purposes
	geog = handles.handMir.geog;      grdMaxSize = handles.handMir.grdMaxSize;
	swathRatio    = handles.handMir.swathRatio;
	scale2meanLat = handles.handMir.scale2meanLat;
	%ForceInsitu = handles.ForceInsitu;     % We don't save it because the user must choose it every time

	% Detect which matlab version is beeing used. For the moment I'm only interested to know if R13 or >= R14
	version7 = version;
	if (str2double(version7(1)) > 6),   version7 = 1;
	else                                version7 = 0;
	end

	if (~version7)                  % R<=13
		save(fname,'geog','grdMaxSize','swathRatio','directory_list','DefLineThick','DefLineColor',...
			'DefineMeasureUnit','DefineEllipsoide','DefineEllipsoide_params', 'scale2meanLat',...
			'flederPlanar', 'flederBurn', 'whichFleder', '-append')
	else
		save(fname,'geog','grdMaxSize','swathRatio','directory_list','DefLineThick','DefLineColor',...
			'DefineMeasureUnit','DefineEllipsoide','DefineEllipsoide_params', 'scale2meanLat',...
			'flederPlanar', 'flederBurn', 'whichFleder', '-append', '-v6')
	end

	% Save the Mirone handles, on the Mirone fig obviously
	guidata(handles.handMir.figure1, handles.handMir)
	setappdata(handles.handMir.figure1,'swathRatio',swathRatio);    % We need this in getline_mb
	set(handles.handMir.ToolsMeasureDist,'Label',['Distance in ' handles.handMir.DefineMeasureUnit])

	delete(handles.figure1)

% ------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
        delete(hObject);
	end

% ------------------------------------------------------------------------------------
function tab_group_Callback(hObject, eventdata, handles)

% ------------------------------------------------------------------------------------
function radio_iview_CB(hObject, eventdata, handles)
    if (get(hObject,'Val'))
        set(handles.radio_fleder,'Val',0)
        handles.whichFleder = 1;
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
function radio_fleder_CB(hObject, eventdata, handles)
    if (get(hObject,'Val'))
        set(handles.radio_iview,'Val',0)
        handles.whichFleder = 0;
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
function radio_planar_Callback(hObject, eventdata, handles)
    if (get(hObject,'Val'))
        handles.flederPlanar = 1;
        set(handles.radio_spherical,'Val',0)
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end
    
% ------------------------------------------------------------------------------------
function radio_spherical_CB(hObject, eventdata, handles)
    if (get(hObject,'Val'))
        handles.flederPlanar = 0;
        set(handles.radio_planar,'Val',0)
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
function radio_coastsOnly_CB(hObject, eventdata, handles)
    if (get(hObject,'Val'))
        set([handles.radio_noBurnAtAll handles.radio_burnAll],'Val',0)
        handles.flederBurn = 1;
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
function radio_noBurnAtAll_CB(hObject, eventdata, handles)
    if (get(hObject,'Val'))
        set([handles.radio_coastsOnly handles.radio_burnAll],'Val',0)
        handles.flederBurn = 0;
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
function radio_burnAll_CB(hObject, eventdata, handles)
    if (get(hObject,'Val'))
        set([handles.radio_coastsOnly handles.radio_noBurnAtAll],'Val',0)
        handles.flederBurn = 2;
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
% ----------------------- Creates and returns a handle to the GUI figure. 
function mirone_pref_LayoutFcn(h1);

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Preferences',...
'NumberTitle','off',...
'Position',[520 437 275 363],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[4 334 61 23],...
'Callback',{@mirone_pref_uicallback,h1,'tab_group_Callback'},...
'Enable','inactive',...
'String','General',...
'ButtonDownFcn',{@mirone_pref_uicallback,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','general');

uicontrol('Parent',h1,'Position',[65 334 70 23],...
'Callback',{@mirone_pref_uicallback,h1,'tab_group_Callback'},...
'Enable','inactive',...
'String','Fledermaus',...
'ButtonDownFcn',{@mirone_pref_uicallback,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[4 5 266 331],'Enable','off','Tag','pushbutton_tab_bg');
uicontrol('Parent',h1,'Position',[10 269 111 50],'Style','frame','Tag','frame1','UserData','general');

uicontrol('Parent',h1,'Position',[20 293 90 15],...
'Callback',{@mirone_pref_uicallback,h1,'radiobutton_geog_Callback'},...
'String','Geographic',...
'Style','radiobutton',...
'TooltipString','GMT grid is in geographical coordinates',...
'Value',1,...
'Tag','radiobutton_geog',...
'UserData','general');

uicontrol('Parent',h1,'Position',[20 274 90 15],...
'Callback',{@mirone_pref_uicallback,h1,'radiobutton_cart_Callback'},...
'String','Cartesian',...
'Style','radiobutton',...
'TooltipString','GMT grid is in cartesian coordinates',...
'Tag','radiobutton_cart',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'edit_GridMaxSize_Callback'},...
'HorizontalAlignment','left',...
'Position',[134 297 36 20],...
'String','20',...
'Style','edit',...
'TooltipString','Grid max size that will be stored in memory',...
'Tag','edit_GridMaxSize',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'edit_swathRatio_Callback'},...
'HorizontalAlignment','left',...
'Position',[134 272 36 20],...
'String','3',...
'Style','edit',...
'TooltipString','Swath Ratio for multibeam planing',...
'Tag','edit_swathRatio',...
'UserData','general');

uicontrol('Parent',h1,'Position',[21 311 90 15],...
'String','Grid coordinates',...
'Style','text',...
'Tag','textGridCoords',...
'UserData','general');

uicontrol('Parent',h1,'Position',[174 300 90 15],...
'HorizontalAlignment','left',...
'String','Grid max size (Mb)',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[174 275 90 15],...
'HorizontalAlignment','left',...
'String','Swath ratio',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'popup_MeasureUnites_Callback'},...
'Position',[10 217 101 22],...
'String',{  'nautic miles'; 'kilometers'; 'meters'; 'user' },...
'Style','popupmenu',...
'TooltipString','Select the default measure units',...
'Value',1,...
'Tag','popup_MeasureUnites',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'popup_ellipsoide_Callback'},...
'Position',[120 217 145 22],...
'String','WGS84',...
'Style','popupmenu',...
'TooltipString','Select the default ellipsoide',...
'Value',1,...
'Tag','popup_ellipsoide',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'popup_directory_list_Callback'},...
'Position',[10 166 236 22],...
'Style','popupmenu',...
'TooltipString','Select the default initial directory from list',...
'Value',1,...
'Tag','popup_directory_list',...
'UserData','general');

uicontrol('Parent',h1,'Position',[246 167 18 21],...
'Callback',{@mirone_pref_uicallback,h1,'pushbutton_change_dir_Callback'},...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'TooltipString','Select a different directory',...
'Tag','pushbutton_change_dir',...
'UserData','general');

uicontrol('Parent',h1,'Position',[10 140 230 15],...
'Callback',{@mirone_pref_uicallback,h1,'checkbox_meanLat_Callback'},...
'String','Scale geog images at mean lat',...
'Style','checkbox',...
'Value',1,...
'Tag','checkbox_meanLat',...
'UserData','general');

uicontrol('Parent',h1,...
'Callback',{@mirone_pref_uicallback,h1,'checkbox_ForceInsitu_Callback'},...
'Position',[10 115 230 15],...
'String','Force "Insitu" transposition',...
'Style','checkbox',...
'Tag','checkbox_ForceInsitu',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'popupmenu_DefLineThickness_Callback'},...
'Position',[10 41 100 22],...
'String',{'1 pt'; '2 pt'; '3 pt'},...
'Style','popupmenu',...
'TooltipString','All drawn lines will have this thickness',...
'Value',1,...
'Tag','popupmenu_DefLineThickness',...
'UserData','general');

h21 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@mirone_pref_uicallback,h1,'popupmenu_DefLineColor_Callback'},...
'Position',[144 41 100 22],...
'String',{  'Black'; 'White'; 'Red'; 'Green'; 'Blue'; 'Cyan'; 'Yellow'; 'Magenta' },...
'Style','popupmenu',...
'TooltipString','All drawn lines will have this color',...
'Value',1,...
'Tag','popupmenu_DefLineColor',...
'UserData','general');

uicontrol('Parent',h1,'Position',[10 188 180 15],...
'HorizontalAlignment','left',...
'String','Default directory',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[10 63 102 15],...
'HorizontalAlignment','left',...
'String','Default line thickness',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[145 63 120 15],...
'HorizontalAlignment','left',...
'String','Default line color',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[10 240 100 15],...
'HorizontalAlignment','left',...
'String','Measure unites',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[121 240 120 15],...
'HorizontalAlignment','left',...
'String','Default ellipsoide',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[120 10 66 23],...
'Callback',{@mirone_pref_uicallback,h1,'pushbutton_OK_Callback'},...
'String','OK',...
'Tag','pushbutton_OK');

uicontrol('Parent',h1,'Position',[198 10 66 23],...
'Callback',{@mirone_pref_uicallback,h1,'pushbutton_cancel_Callback'},...
'String','Cancel',...
'Tag','pushbutton_cancel');

% -------------------- 	FLEDERMAUS TAB ---------------------------------------- 
uicontrol('Parent',h1,'Position',[27 280 220 30],...
'FontAngle','italic','FontSize',12,'FontWeight','demi',...
'HorizontalAlignment','left',...
'String','Which one?',...
'Style','text',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 265 120 15],...
'Callback',{@mirone_pref_uicallback,h1,'radio_iview_CB'},...
'String','The free viewer',...
'Style','radiobutton',...
'TooltipString','View the fleder files using the iview3d free viewer',...
'Value',1,...
'Tag','radio_iview',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[146 265 120 15],...
'Callback',{@mirone_pref_uicallback,h1,'radio_fleder_CB'},...
'String','The real thing',...
'Style','radiobutton',...
'TooltipString','Open the fleder files using the (not free) fledermaus',...
'Tag','radio_fleder',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[27 205 220 30],...
'FontAngle','italic','FontSize',12,'FontWeight','demi',...
'HorizontalAlignment','left',...
'String','How to deal with Images',...
'Style','text',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 190 110 15],...
'Callback',{@mirone_pref_uicallback,h1,'radio_planar_Callback'},...
'String','Planar Image',...
'Style','radiobutton',...
'TooltipString','Create planar (2D) images',...
'Value',1,...
'Tag','radio_planar',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[146 190 110 15],...
'Callback',{@mirone_pref_uicallback,h1,'radio_spherical_CB'},...
'String','Spherical Image',...
'Style','radiobutton',...
'TooltipString','Create spherical images (wraped on the sphere)',...
'Tag','radio_spherical',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 130 220 30],...
'FontAngle','italic','FontSize',12,'FontWeight','demi',...
'HorizontalAlignment','left',...
'String','How to deal with Vectors',...
'Style','text',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 116 180 15],...
'Callback',{@mirone_pref_uicallback,h1,'radio_coastsOnly_CB'},...
'String','Burn coastlines only',...
'Style','radiobutton',...
'TooltipString','Burn coast lines into the image before creating the fleder file',...
'Value',1,...
'Tag','radio_coastsOnly',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 86 180 15],...
'Callback',{@mirone_pref_uicallback,h1,'radio_burnAll_CB'},...
'String','Burn them all',...
'Style','radiobutton',...
'TooltipString','Burn all lines into the image before creating the fleder file',...
'Tag','radio_burnAll',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 56 180 15],...
'Callback',{@mirone_pref_uicallback,h1,'radio_noBurnAtAll_CB'},...
'String','No burning at all',...
'Style','radiobutton',...
'TooltipString','All lines are converted into vectors when creating the fleder file',...
'Tag','radio_noBurnAtAll',...
'UserData','fleder');

function mirone_pref_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

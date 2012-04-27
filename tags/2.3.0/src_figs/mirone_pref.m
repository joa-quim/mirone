function varargout = mirone_pref(varargin)
% Helper window to select some importat defaults

%	Copyright (c) 2004-2012 by J. Luis
%
% 	This program is part of Mirone and is free software; you can redistribute
% 	it and/or modify it under the terms of the GNU Lesser General Public
% 	License as published by the Free Software Foundation; either
% 	version 2.1 of the License, or any later version.
% 
% 	This program is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% 	Lesser General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------
 
	hObject = figure('Vis','off');
	mirone_pref_LayoutFcn(hObject);
	handles = guihandles(hObject);

    handMir = varargin{1};
    home_dir = handMir.home_dir;
    handles.handMir = handMir;
	move2side(handMir.figure1, hObject, 'left');

	directory_list = [];
	load([handMir.path_data 'mirone_pref.mat']);
	handles.geog = geog;		% Just to not be empty.
	handles.ForceInsitu = 0;	% 		"
	handles.moveDoubleClick = 1;% 		"
	handles.flederPlanar = 1;	%		"
	handles.flederBurn = 1;
	handles.whichFleder = 1;
	handles.bg_color = [1 1 1];		% Default is white, but should be update by mirone_pref contents
	handles.proxyAddress = [];
	handles.proxyPort = [];
	handles.proxyAddressPort = [];

	handles.flederPlanar = flederPlanar;
	set(handles.radio_planar,'Val',handles.flederPlanar)
	set(handles.radio_spherical,'Val',~handles.flederPlanar)
	handles.flederBurn = flederBurn;
	if (handles.flederBurn == 0)
		set(handles.radio_noBurnAtAll,'Val',1);		set(handles.radio_coastsOnly,'Val',0)
	elseif (handles.flederBurn == 1)
		set(handles.radio_coastsOnly,'Val',1);
	else
		set(handles.radio_burnAll,'Val',1);			set(handles.radio_coastsOnly,'Val',0)
	end
	handles.whichFleder = whichFleder;	% whichFleder = 1 for the free iview3d or 0 for the true thing (fledermaus)
	if (~whichFleder)
		set(handles.radio_fleder,'Val',1),		set(handles.radio_iview,'Val',0)
	end
	handles.moveDoubleClick = moveDoubleClick;

	try
		handles.bg_color = nanColor;		% Wrap it into a try while in probation period
	end

	if iscell(directory_list)							% When exists a dir list in mirone_pref		
		j = false(1,numel(directory_list));				% vector for eventual cleaning non-existing dirs
		for (i = 1:numel(directory_list))				% Check that all dirs in last_directories exist
			j(i) = (exist(directory_list{i},'dir') ~= 7);
		end
		directory_list(j) = [];							% clean non-existing directories

		if (~isempty(directory_list))					% If there is one left
			set(handles.popup_directory_list,'String',directory_list)
			handles.last_directories = directory_list;
		else
			handles.last_directories = {[home_dir filesep 'tmp']; home_dir};    % Let it have something existent
			set(handles.popup_directory_list,'String',handles.last_directories)
		end
	else												% mirone_pref had no dir list
		handles.last_directories = {[home_dir filesep 'tmp']; home_dir};    % Let it have something existent
		set(handles.popup_directory_list,'String',handles.last_directories)
	end

    if (handMir.geog)             % Signals a geographic grid/image
		set(handles.radio_geog,'Value',1)
		set(handles.radio_cart,'Value',0)
		handles.geog = 1;
    elseif (handMir.geog == 0)
		set(handles.radio_geog,'Value',0)
		set(handles.radio_cart,'Value',1)
		handles.geog = 0;
    else
		handles.geog = 1;
    end
    set(handles.edit_GridMaxSize,'String',sprintf('%d',fix(handMir.grdMaxSize / (2^20))))
    set(handles.edit_swathRatio,'String',sprintf('%g',handMir.swathRatio))
    set(handles.checkbox_ForceInsitu,'Value',handMir.ForceInsitu)
	set(handles.check_movePolyg,'Val',handles.moveDoubleClick)
    handles.ForceInsitu = handMir.ForceInsitu;

	% Well this is split from the above because it was written later and I don't want to mess
	% with what is working. Wrap in a try-catch because the first time the variables are not
	% yet in mirone_pref.mat
	try     % Goes here all other times
        set(handles.popup_DefLineThickness,'String',DefLineThick)
        set(handles.popup_DefLineColor,'String',DefLineColor)
        set(handles.popup_MeasureUnites,'String',DefineMeasureUnit)
        set(handles.popup_ellipsoide,'String',DefineEllipsoide)
        set(handles.checkbox_meanLat,'Value',scale2meanLat)
	catch       % Comes here in first call before variables are stored in mirone_pref.mat
        DefLineThick = {'2 pt'; '1 pt'; '3 pt'; '4 pt'};
        DefLineColor = {'White'; 'Black'; 'Red'; 'Green'; 'Blue'; 'Cyan'; 'Yellow'; 'Magenta'};
        DefineMeasureUnit = {'nautic miles'; 'kilometers'; 'meters'; 'user'};
        set(handles.popup_DefLineThickness,'String',DefLineThick)
        set(handles.popup_DefLineColor,'String',DefLineColor)
        set(handles.popup_MeasureUnites,'String',DefineMeasureUnit)
        set(handles.checkbox_meanLat,'Value',1)
	end

	% This is the default ellipsoide order. It will be changed (and saved as so in mirone_pref) by the user
	% The parameters of the selected ellipsoid are found by a "case" loop runned by the OK pushbutton.
	was_ellips_from_file = false;
	if (exist([handMir.path_data 'ellipsoids.txt'], 'file') == 2)		% If we have the ellipsoids provied in a file
		fid = fopen([handMir.path_data 'ellipsoids.txt'],'rt');
		if (fid > 0)				% Otherwise, non-existent file, revert to just WGS-84
			ellips = strread(fread(fid,'*char').','%s','delimiter','\n');		fclose(fid);
			% Remove comment and also eventual empty lines
			m = numel(ellips);		c = false(m,1);
			for (k = 1:m)
				if ( isempty(ellips{k}) || ellips{k}(1) == '#' ),		c(k) = true;	end
			end
			ellips(c) = [];			n_ellips = numel(ellips);
			handles.ellipsoide = cell(n_ellips,4);
			for (k = 1:n_ellips)					% Loop over number of ellipsoides
				ind = strfind(ellips{k},',');
				handles.ellipsoide{k,1} = ellips{k}(1:ind(1)-1);		% Ellipsoid name
				handles.ellipsoide{k,2} = str2double(ellips{k}(ind(1)+1:ind(2)-1));		% Major axis
				handles.ellipsoide{k,3} = 0;
				tmp = ellips{k}(ind(2)+1:end);
				id = strfind(tmp,'/');
				if (~isempty(id))
					handles.ellipsoide{k,4} = 1 / str2double(tmp(id+1:end));		% Flattening given as 1/...
				else
					handles.ellipsoide{k,4} = str2double(tmp);
				end
			end
			was_ellips_from_file = true;
			set(handles.popup_ellipsoide,'String',handles.ellipsoide(:,1))
		else
			handles.ellipsoide = ellips_list;		% Use builtin list
		end

	else
		handles.ellipsoide = ellips_list;			% Use builtin list
	end

	% This test should resolve the case when we used to have a ellipsoids.txt file but it
	% was removed so now we gat back to the builtin ellipsoids. However, it will fail if the
	% ellipsoids file used to have the same number of ellipsoids as the builtins (30)
	if (~was_ellips_from_file && numel(DefineEllipsoide) ~= size(handles.ellipsoide, 1) )
		set(handles.popup_ellipsoide,'String',handles.ellipsoide(:,1))
	end

	try
		if (handles.geog == 0)      % For cartesian coords the following is no applyable
			set(handles.popup_ellipsoide,'Enable','off')
			set(handles.popup_MeasureUnites,'Enable','off')
		end
	catch   % In case of error, set the default list.
		set(handles.popup_ellipsoide,'String',handles.ellipsoide(:,1))
		set(handles.popup_MeasureUnites,'String',{'nautic miles'; 'kilometers'; 'meters'; 'user'})
	end

	% Create the "ForceInsitu" Tooltip
	str = sprintf(['Importing grids implies a conversion that uses\n'...
		'matrix transposition. This operations is fast if\n'...
		'we make a copy of the importing grid. However,\n'...
		'this requires twice the grid size on memory.\n'...
		'If you don''t have enough memory to import a\n'...
		'large grid, use this option that do the transposition\n'...
		'"insitu". That is, it uses only one time the grid\n'...
		'size in memory. The price you will pay, however,\n'...
		'is in speed because it runs about 10 times slower.']);
	set(handles.checkbox_ForceInsitu,'Tooltip',str)

	str = sprintf(['Since 1 degree of longitude and latitude do not cover the same\n'...
		'arc length at Earth surface, isometric plotting of geographical\n'...
		'grids squeezes the image vertically. Scaling the image to the\n'...
		'cosinus of the mean lat minimizes this effect.']);
	set(handles.checkbox_meanLat,'Tooltip',str)

	str = sprintf(['Controls what mouse selection is used to move polylines/patches\n'...
			'If checked lines are moved with a left click drag-n-drop. Though\n',...
			'easier to operate this has often anoying side effects (move the polyline\n',...
			'when in fact we wanted to edit one of its vertex).\n\n',...
			'Uncheck if you want to use a Shift-click left mouse button or click both\n',...
			'left and right mouse buttons to move the line. A bit more cumbersome, but safer']);
	set(handles.check_movePolyg,'Tooltip',str)

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, handles.txt_GC)
	%------------- END Pro look (3D) ------------------------------

	set(hObject,'Vis','on');

	% ------------ Paint the NaN color pushbutton -----------------
	siz = get(handles.push_NaNcolor, 'Pos');
	img = uint8(zeros(siz(4)-3,siz(3)-3,3));
	for (k = 1:3),	img(:,:,k) = round(handles.bg_color(k)*255);		end
	set(handles.push_NaNcolor,'CData', img)
	% -------------------------------------------------------------

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
	panel_names = {'general','fleder','more'};

	% tabpanelfcn('makegroups',...) adds new fields to the handles structure,
	% one for each panel name and another called 'group_name_all'.  These fields
	% are used by the tabpanefcn when tab_group_handler is called.
	handles = tabpanelfcn('make_groups',group_name, panel_names, handles, 1);
	% ------------------------------------------------------------------------------

	% Choose default command line output for mirone_pref
	guidata(hObject, handles);

	if (nargout),	varargout{1} = hObject;		end

% -------------------------------------------------------------------------------------
function tab_group_ButtonDownFcn(hObject, handles)
% Call the tab_group_handler.  This updates visiblity of components as needed to
% hide the components from the previous tab and show components on this tab.
% This also updates the last_tab field in the handles structure to keep track
% of which panel was hidden.
    handles = tabpanelfcn('tab_group_handler',hObject, handles, get(hObject, 'Tag'));
    guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function ellipsoids = ellips_list()
	ellipsoids = {'WGS-84 - 1984', 6378137.0, 0.0, 1.0/298.2572235630;
			'OSU91A - 1991', 6378136.3, 0.0, 1.0/298.25722;
			'Engelis - 1985', 6378136.05, 0.0, 1.0/298.2566;
			'SGS-85 - 1985', 6378136.0, 0.0, 1.0/298.257;
			'MERIT-83 - 1983', 6378137.0, 0.0, 1.0/298.257;
			'GRS-80 - 1980', 6378137.0, 0.0, 1.0/298.257222101;
			'IAG-75 - 1975', 6378140.0, 0.0, 1.0/298.257222;
			'Indonesian - 1974', 6378160.0, 0.0, 1.0/298.247;
			'WGS-72 - 1972', 6378135.0, 0.0, 1.0/298.26;
			'WGS-66 - 1966', 6378145.0, 0.0, 1.0/298.25;
			'WGS-60 - 1960', 6378165.0, 0.0, 1.0/298.3;
			'South-American - 1969', 6378160.0, 0.0, 1.0/298.25;
			'Fischer-1968', 6378150.0, 0.0, 1.0/298.3;
			'GRS-67 - 1967', 6378160.0, 0.0, 1.0/298.247167427;
			'International-1967', 6378157.5, 0.0, 1.0/298.25;
			'Australian - 1965', 6378160.0, 0.0, 1.0/298.25;
			'Hough - 1960', 6378270.0, 0.0, 1.0/297.0;
			'Krassovsky - 1940', 6378245.0, 0.0, 1.0/298.3;
			'International-1924', 6378388.0, 0.0, 1.0/297.0;
			'Hayford-1909', 6378388.0, 0.0, 1.0/297.0;
			'Helmert-1906', 6378200.0, 0.0, 1.0/298.3;
			'Clarke-1880', 6378249.145, 0.0, 1.0/293.465;
			'Andrae - 1876', 6377104.43, 0.0, 1.0/300.0;
			'Airy - 1830', 6377563.396, 0.0, 1.0/299.3249646;
			'Modified-Airy - 1830', 6377340.189, 0.0, 1.0/299.3249646;
			'Bessel - 1841', 6377397.155, 0.0, 1.0/299.1528128;
			'Bessel-Namibia - 1841', 6377483.865, 0.0, 1.0/299.1528128;
			'Everest-1830', 6377276.345, 0.0, 1.0/300.8017;
			'Sphere - 1980', 6371008.7714, 0.0, 0.0;
			'Mars Sphere', 3396000.0, 0.0, 0.0};
	
% ------------------------------------------------------------------------------------
function radio_geog_CB(hObject, handles)
	if get(hObject,'Value')     handles.geog = 1;   set(handles.radio_cart,'Value',0)
	else                        handles.geog = 0;   set(handles.radio_cart,'Value',1);
	end
	guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function radio_cart_CB(hObject, handles)
	if get(hObject,'Value')     handles.geog = 0;   set(handles.radio_geog,'Value',0)
	else                        handles.geog = 1;   set(handles.radio_geog,'Value',1);
	end
	guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function edit_GridMaxSize_CB(hObject, handles)
	xx = get(hObject,'String');
	if isnan(str2double(xx)) || isempty(xx)    % Just a stupid user error
        set(hObject, 'String', '50')
	else
        if (str2double(xx) >= 10)   % Numbers bigger than 10 Mb don't need decimal places
            set(hObject,'String',sprintf('%d',fix(str2double(xx))))
        end
	end

% ------------------------------------------------------------------------------------
function edit_swathRatio_CB(hObject, handles)
	xx = get(hObject,'String');
	if isnan(str2double(xx)) || isempty(xx)    % Just a stupid user error
        set(hObject, 'String', '3')
	end

% ------------------------------------------------------------------------------------
function popup_directory_list_CB(hObject, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	% Put the selected field on top of the String list. This is necessary because the "OK" button will
	% read this list and save it in mirone_pref_export, so next time the selected field will show up first.
	tmp = str(val);         str(val) = [];
	new_str = [tmp; str];   set(hObject,'String',new_str); 
	set(hObject,'Value',1)

% ------------------------------------------------------------------------------------
function push_change_dir_CB(hObject, handles)
	contents = get(handles.popup_directory_list,'String');
	if (strcmp(computer, 'PCWIN'))
		work_dir = uigetfolder_win32('Select a directory', contents{get(handles.popup_directory_list,'Value')});
	else            % This guy says it cannot be compiled
		work_dir = uigetdir(contents{get(handles.popup_directory_list,'Value')}, 'Select a directory');
	end
	
	if (isempty(work_dir) || isequal(work_dir,0))    return;     end
	
	handles.last_directories = [cellstr(work_dir); handles.last_directories];
	if length(handles.last_directories) > 15            % Keep only 15 adresses
		handles.last_directories(16:end) = [];
	end
	set(handles.popup_directory_list,'String',handles.last_directories)
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function checkbox_meanLat_CB(hObject, handles)
	handles.scale2meanLat = get(hObject,'Value');
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function checkbox_ForceInsitu_CB(hObject, handles)
	handles.ForceInsitu = get(hObject,'Value');
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function check_movePolyg_CB(hObject, handles)
	handles.moveDoubleClick = get(hObject,'Value');
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function popup_DefLineThickness_CB(hObject, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	% Put the selected field on top of the String list. This is necessary because the "OK" button will
	% read this list and save it in mirone_pref, so next time the selected field will show up first.
	tmp = str(val);         str(val) = [];
	new_str = [tmp; str];   set(hObject,'String',new_str); 
	set(hObject,'Value',1)

% ------------------------------------------------------------------------------------
function popup_DefLineColor_CB(hObject, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	% Put the selected field on top of the String list. This is necessary because the "OK" button will
	% read this list and save it in mirone_pref, so next time the selected field will show up first.
	tmp = str(val);         str(val) = [];
	new_str = [tmp; str];   set(hObject,'String',new_str); 
	set(hObject,'Value',1)

% ------------------------------------------------------------------------------------
function popup_MeasureUnites_CB(hObject, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	% Put the selected field on top of the String list. This is necessary because the "OK" button will
	% read this list and save it in mirone_pref, so next time the selected field will show up first.
	tmp = str(val);         str(val) = [];
	new_str = [tmp; str];   set(hObject,'String',new_str); 
	set(hObject,'Value',1)

% ------------------------------------------------------------------------------------
function popup_ellipsoide_CB(hObject, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	% Put the selected field on top of the String list. This is necessary because the "OK" button will
	% read this list and save it in mirone_pref, so next time the selected field will show up first.
	tmp = str(val);         str(val) = [];
	new_str = [tmp; str];   set(hObject,'String',new_str); 
	set(hObject,'Value',1)

% ------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
	handles.handMir.geog = handles.geog;
	handles.handMir.grdMaxSize = str2double(get(handles.edit_GridMaxSize,'String')) * 2^20;
	handles.handMir.swathRatio = str2double(get(handles.edit_swathRatio,'String'));
	directory_list = get(handles.popup_directory_list, 'String');
	handles.handMir.last_dir   = directory_list{1};
	handles.handMir.work_dir   = handles.handMir.last_dir;
	DefLineThick = get(handles.popup_DefLineThickness, 'String');
	DefLineColor = get(handles.popup_DefLineColor, 'String');
	handles.handMir.scale2meanLat = get(handles.checkbox_meanLat,'Value');
	handles.handMir.ForceInsitu = handles.ForceInsitu;
	handles.handMir.flederPlanar = handles.flederPlanar;		flederPlanar = handles.flederPlanar;
	handles.handMir.flederBurn = handles.flederBurn;			flederBurn = handles.flederBurn;
	handles.handMir.whichFleder = handles.whichFleder;			whichFleder = handles.whichFleder;
	% Decode the line thickness string into a number
	handles.handMir.DefLineThick = str2double(DefLineThick{1}(1));
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
	for i = 1:length(handles.ellipsoide)
		switch DefineEllipsoide{1}
			case handles.ellipsoide(i)
				handles.handMir.DefineEllipsoide(1) = handles.ellipsoide{i,2};
				handles.handMir.DefineEllipsoide(2) = handles.ellipsoide{i,3};
				handles.handMir.DefineEllipsoide(3) = handles.ellipsoide{i,4};
		end
	end

	delete(handles.figure1)		% Killing it here will make look that the tool runs faster

	fname = [handles.handMir.path_data 'mirone_pref.mat'];
	% Save the preferences to a mat file under the data directory
	% Note: for the ellipsoid we save it's parameters (a,b,f) and name on separate vars
	DefineEllipsoide_params = handles.handMir.DefineEllipsoide;    % For saving purposes
	geog = handles.handMir.geog;      grdMaxSize = handles.handMir.grdMaxSize;
	swathRatio    = handles.handMir.swathRatio;
	scale2meanLat = handles.handMir.scale2meanLat;
	%ForceInsitu = handles.ForceInsitu;     % We don't save it because the user must choose it every time
	moveDoubleClick = handles.moveDoubleClick;
	nanColor = handles.bg_color;
	handles.handMir.bg_color = handles.bg_color;

	% Detect which matlab version is beeing used. For the moment I'm only interested to know if R13 or >= R14
	version7 = version;
	V7 = (sscanf(version7(1),'%f') > 6);

	if (~V7)                  % R <= 13
		save(fname,'geog','grdMaxSize','swathRatio','directory_list','DefLineThick','DefLineColor',...
			'DefineMeasureUnit','DefineEllipsoide','DefineEllipsoide_params', 'scale2meanLat',...
			'flederPlanar', 'flederBurn', 'whichFleder', 'moveDoubleClick', 'nanColor', '-append')
	else
		save(fname,'geog','grdMaxSize','swathRatio','directory_list','DefLineThick','DefLineColor',...
			'DefineMeasureUnit','DefineEllipsoide','DefineEllipsoide_params', 'scale2meanLat',...
			'flederPlanar', 'flederBurn', 'whichFleder', 'moveDoubleClick', 'nanColor', '-append', '-v6')
	end

	% Save the Mirone handles, on the Mirone fig obviously
	guidata(handles.handMir.figure1, handles.handMir)
	setappdata(handles.handMir.figure1,'swathRatio',swathRatio);    % We need this in getline_mb
	set(handles.handMir.ToolsMeasureDist,'Label',['Distance in ' handles.handMir.DefineMeasureUnit])
	if (moveDoubleClick)			% this info is used by UI_EDIT_POLYGON()
		setappdata(handles.handMir.axes1,'MovPolyg',[])				% Move lines with a click drag-n-drop
	else
		setappdata(handles.handMir.axes1,'MovPolyg','extend')		% Move lines with a Shift-click drag-n-drop
	end

% ------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
        delete(hObject);
	end

% ------------------------------------------------------------------------------------
function tab_group_CB(hObject, handles)

% ------------------------------------------------------------------------------------
function radio_iview_CB(hObject, handles)
    if (get(hObject,'Val'))
        set(handles.radio_fleder,'Val',0)
        handles.whichFleder = 1;
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
function radio_fleder_CB(hObject, handles)
    if (get(hObject,'Val'))
        set(handles.radio_iview,'Val',0)
        handles.whichFleder = 0;
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
function radio_planar_CB(hObject, handles)
    if (get(hObject,'Val'))
        handles.flederPlanar = 1;
        set(handles.radio_spherical,'Val',0)
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end
    
% ------------------------------------------------------------------------------------
function radio_spherical_CB(hObject, handles)
    if (get(hObject,'Val'))
        handles.flederPlanar = 0;
        set(handles.radio_planar,'Val',0)
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
function radio_coastsOnly_CB(hObject, handles)
    if (get(hObject,'Val'))
        set([handles.radio_noBurnAtAll handles.radio_burnAll],'Val',0)
        handles.flederBurn = 1;
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
function radio_noBurnAtAll_CB(hObject, handles)
    if (get(hObject,'Val'))
        set([handles.radio_coastsOnly handles.radio_burnAll],'Val',0)
        handles.flederBurn = 0;
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
function radio_burnAll_CB(hObject, handles)
    if (get(hObject,'Val'))
        set([handles.radio_coastsOnly handles.radio_noBurnAtAll],'Val',0)
        handles.flederBurn = 2;
        guidata(handles.figure1, handles);
    else
        set(hObject,'Val',1)
    end

% ------------------------------------------------------------------------------------
function push_NaNcolor_CB(hObject, handles)
% ...
	c = uisetcolor;
	if (numel(c) > 1)
		siz = get(hObject, 'Pos');
		img = uint8(zeros(siz(4)-3,siz(3)-3,3));
		for (k = 1:3),	img(:,:,k) = round(c(k)*255);		end
		set(hObject,'CData', img)
		handles.bg_color = c;
		guidata(handles.figure1, handles)
	end

% -------------------------------------------------------------------------------------
function check_proxy_CB(hObject, handles)
	if ( get(hObject, 'Val') )
		set([handles.edit_proxyAddress handles.edit_proxyPort handles.text_proxyPort handles.text_proxyAddress], 'Enable', 'on')
		ind = strfind(handles.proxyAddressPort, ':');
		if (~isempty(ind))					% We already have a full proxy adress
			set_gmt(['http_proxy=' handles.proxyAddressPort]);
		end
	else
		set([handles.edit_proxyAddress handles.edit_proxyPort handles.text_proxyPort handles.text_proxyAddress], 'Enable', 'off')
		set_gmt('http_proxy=');
	end
	guidata(handles.figure1, handles)		% as a request to change the proxy settings

% -------------------------------------------------------------------------------------
function edit_proxyAddress_CB(hObject, handles)
	handles.proxyAddress = get(hObject, 'String');
	handles.proxyAddressPort = [handles.proxyAddress ':' handles.proxyPort];
	if (~isempty(handles.proxyAddressPort))
		set_gmt(['http_proxy=' handles.proxyAddressPort]);
	end
	guidata(handles.figure1, handles)
	
% -------------------------------------------------------------------------------------
function edit_proxyPort_CB(hObject, handles)
	handles.proxyPort = get(hObject, 'String');
	handles.proxyAddressPort = [handles.proxyAddress ':' handles.proxyPort];
	if (~isempty(handles.proxyPort))
		set_gmt(['http_proxy=' handles.proxyAddressPort]);
	end
	guidata(handles.figure1, handles)

% ------------------------------------------------------------------------------------
% ----------------------- Creates and returns a handle to the GUI figure. 
function mirone_pref_LayoutFcn(h1)

fUiBgColor = get(0,'factoryUicontrolBackgroundColor');
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

uicontrol('Parent',h1,'Position',[4 333 61 22],...
'Call',{@mirone_pref_uiCB,h1,'tab_group_CB'},...
'Enable','inactive',...
'String','General',...
'ButtonDownFcn',{@mirone_pref_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','general');

uicontrol('Parent',h1,'Position',[65 333 80 22],...
'Call',{@mirone_pref_uiCB,h1,'tab_group_CB'},...
'Enable','inactive',...
'String','Fledermaus',...
'ButtonDownFcn',{@mirone_pref_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[145 333 80 22],...
'Call',{@mirone_pref_uiCB,h1,'tab_group_CB'},...
'Enable','inactive',...
'String','More',...
'ButtonDownFcn',{@mirone_pref_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','more');

uicontrol('Parent',h1,'Position',[4 5 266 331],'Enable','off','BackgroundColor',fUiBgColor);
uicontrol('Parent',h1,'Position',[10 269 111 50],'Style','frame','UserData','general');

uicontrol('Parent',h1,'Position',[20 293 90 15],...
'Call',{@mirone_pref_uiCB,h1,'radio_geog_CB'},...
'String','Geographic',...
'Style','radiobutton',...
'Tooltip','Grid is in geographical coordinates',...
'Value',1,...
'Tag','radio_geog',...
'UserData','general');

uicontrol('Parent',h1,'Position',[20 274 90 15],...
'Call',{@mirone_pref_uiCB,h1,'radio_cart_CB'},...
'String','Cartesian',...
'Style','radiobutton',...
'Tooltip','Grid is in cartesian coordinates',...
'Tag','radio_cart',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@mirone_pref_uiCB,h1,'edit_GridMaxSize_CB'},...
'HorizontalAlignment','center',...
'Position',[134 297 36 21],...
'String','20',...
'Style','edit',...
'Tooltip','Grid max size that will be stored in memory',...
'Tag','edit_GridMaxSize',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@mirone_pref_uiCB,h1,'edit_swathRatio_CB'},...
'HorizontalAlignment','center',...
'Position',[134 272 36 21],...
'String','3',...
'Style','edit',...
'Tooltip','Swath Ratio for multibeam planing',...
'Tag','edit_swathRatio',...
'UserData','general');

%uicontrol('Parent',h1,'Units','characters','Position',[4 23.846 18.0 1.154],...
uicontrol('Parent',h1,'Position',[21 311 90 15],...
'BackgroundColor',fUiBgColor,...
'String','Grid coordinates',...
'Style','text',...
'Tag','txt_GC',...
'UserData','general');

uicontrol('Parent',h1,'Position',[172 300 84 15],...
'HorizontalAlignment','left',...
'String','Grid max size (Mb)',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[172 275 70 15],...
'HorizontalAlignment','left',...
'String','Swath ratio',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@mirone_pref_uiCB,h1,'popup_MeasureUnites_CB'},...
'Position',[10 217 101 22],...
'String',{'nautic miles'; 'kilometers'; 'meters'; 'user'},...
'Style','popupmenu',...
'Tooltip','Select the default measure units',...
'Value',1,...
'Tag','popup_MeasureUnites',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@mirone_pref_uiCB,h1,'popup_ellipsoide_CB'},...
'Position',[120 217 145 22],...
'String','WGS84',...
'Style','popupmenu',...
'Tooltip','Select the default ellipsoide',...
'Value',1,...
'Tag','popup_ellipsoide',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@mirone_pref_uiCB,h1,'popup_directory_list_CB'},...
'Position',[10 166 236 22],...
'Style','popupmenu','String',' ', ...
'Tooltip','Select the default initial directory from list',...
'Value',1,...
'Tag','popup_directory_list',...
'UserData','general');

uicontrol('Parent',h1,'Position',[246 167 18 21],...
'Call',{@mirone_pref_uiCB,h1,'push_change_dir_CB'},...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'Tooltip','Select a different directory',...
'Tag','push_change_dir',...
'UserData','general');

uicontrol('Parent',h1,'Position',[10 140 230 15],...
'Call',{@mirone_pref_uiCB,h1,'checkbox_meanLat_CB'},...
'String','Scale geog images at mean lat',...
'Style','checkbox',...
'Value',1,...
'Tag','checkbox_meanLat',...
'UserData','general');

uicontrol('Parent',h1,...
'Call',{@mirone_pref_uiCB,h1,'checkbox_ForceInsitu_CB'},...
'Position',[10 115 230 15],...
'String','Force "Insitu" transposition',...
'Style','checkbox',...
'Tag','checkbox_ForceInsitu',...
'UserData','general');

uicontrol('Parent',h1,...
'Call',{@mirone_pref_uiCB,h1,'check_movePolyg_CB'},...
'Position',[10 90 230 15],...
'String','Move lines with a left-click',...
'Style','checkbox',...
'Value',1,...
'Tag','check_movePolyg',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@mirone_pref_uiCB,h1,'popup_DefLineThickness_CB'},...
'Position',[10 41 100 22],...
'String',{'1 pt'; '2 pt'; '3 pt'},...
'Style','popupmenu',...
'Tooltip','All drawn lines will have this thickness',...
'Value',1,...
'Tag','popup_DefLineThickness',...
'UserData','general');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@mirone_pref_uiCB,h1,'popup_DefLineColor_CB'},...
'Position',[144 41 100 22],...
'String',{  'Black'; 'White'; 'Red'; 'Green'; 'Blue'; 'Cyan'; 'Yellow'; 'Magenta' },...
'Style','popupmenu',...
'Tooltip','All drawn lines will have this color',...
'Value',1,...
'Tag','popup_DefLineColor',...
'UserData','general');

uicontrol('Parent',h1,'Position',[15 188 100 15],...
'HorizontalAlignment','left',...
'String','Default directory',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[10 63 125 15],...
'HorizontalAlignment','left',...
'String','Default line thickness',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[148 63 100 15],...
'HorizontalAlignment','left',...
'String','Default line color',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[13 240 95 15],...
'HorizontalAlignment','left',...
'String','Measure unites',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[123 240 105 15],...
'HorizontalAlignment','left',...
'String','Default ellipsoide',...
'Style','text',...
'UserData','general');

uicontrol('Parent',h1,'Position',[198 10 66 21],...
'Call',{@mirone_pref_uiCB,h1,'push_OK_CB'},...
'String','OK',...
'Tag','push_OK');

% -------------------- 	FLEDERMAUS TAB ---------------------------------------- 
uicontrol('Parent',h1,'Position',[27 280 220 30],...
'FontAngle','italic','FontSize',12,'FontWeight','demi',...
'HorizontalAlignment','left',...
'String','Which one?',...
'Style','text',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 265 120 15],...
'Call',{@mirone_pref_uiCB,h1,'radio_iview_CB'},...
'String','The free viewer',...
'Style','radiobutton',...
'Tooltip','View the fleder files using the iview3d free viewer',...
'Value',1,...
'Tag','radio_iview',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[146 265 120 15],...
'Call',{@mirone_pref_uiCB,h1,'radio_fleder_CB'},...
'String','The real thing',...
'Style','radiobutton',...
'Tooltip','Open the fleder files using the (not free) fledermaus',...
'Tag','radio_fleder',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[27 205 220 30],...
'FontAngle','italic','FontSize',12,'FontWeight','demi',...
'HorizontalAlignment','left',...
'String','How to deal with Images',...
'Style','text',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 190 110 15],...
'Call',{@mirone_pref_uiCB,h1,'radio_planar_CB'},...
'String','Planar Image',...
'Style','radiobutton',...
'Tooltip','Create planar (2D) images',...
'Value',1,...
'Tag','radio_planar',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[146 190 110 15],...
'Call',{@mirone_pref_uiCB,h1,'radio_spherical_CB'},...
'String','Spherical Image',...
'Style','radiobutton',...
'Tooltip','Create spherical images (wraped on the sphere)',...
'Tag','radio_spherical',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 130 220 30],...
'FontAngle','italic','FontSize',12,'FontWeight','demi',...
'HorizontalAlignment','left',...
'String','How to deal with Vectors',...
'Style','text',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 116 180 15],...
'Call',{@mirone_pref_uiCB,h1,'radio_coastsOnly_CB'},...
'String','Burn coastlines only',...
'Style','radiobutton',...
'Tooltip','Burn coast lines into the image before creating the fleder file',...
'Value',1,...
'Tag','radio_coastsOnly',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 86 180 15],...
'Call',{@mirone_pref_uiCB,h1,'radio_burnAll_CB'},...
'String','Burn them all',...
'Style','radiobutton',...
'Tooltip','Burn all lines into the image before creating the fleder file',...
'Tag','radio_burnAll',...
'UserData','fleder');

uicontrol('Parent',h1,'Position',[28 56 180 15],...
'Call',{@mirone_pref_uiCB,h1,'radio_noBurnAtAll_CB'},...
'String','No burning at all',...
'Style','radiobutton',...
'Tooltip','All lines are converted into vectors when creating the fleder file',...
'Tag','radio_noBurnAtAll',...
'UserData','fleder');

% -------------------- The 'MORE' TAB ---------------------------------------- 
uicontrol('Parent',h1,'Position',[27 280 120 23],...
'FontSize',11,...
'HorizontalAlignment','left',...
'String','NaN Color',...
'Tooltip','Background color used to paint NaNs',...
'Style','text',...
'UserData','more');

uicontrol('Parent',h1,'Position',[200 280 21 21],...
'Call',{@mirone_pref_uiCB,h1,'push_NaNcolor_CB'},...
'UserData','more',...
'Tag','push_NaNcolor');

uicontrol('Parent',h1, 'Position',[25 223 55 15],...
'Call',{@mirone_pref_uiCB,h1,'check_proxy_CB'},...
'String','proxy?',...
'Style','checkbox',...
'UserData','more',...
'Tag','check_proxy');

uicontrol('Parent',h1, 'Position',[85 220 96 20],...
'Call',{@mirone_pref_uiCB,h1,'edit_proxyAddress_CB'},...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','proxy adress here',...
'UserData','more',...
'Tag','edit_proxyAddress');

uicontrol('Parent',h1, 'Position',[190 220 51 20],...
'Call',{@mirone_pref_uiCB,h1,'edit_proxyPort_CB'},...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','port here',...
'UserData','more',...
'Tag','edit_proxyPort');

uicontrol('Parent',h1,'Position',[105 240 60 15],...
'Enable','off',...
'String','Address',...
'Style','text',...
'UserData','more',...
'Tag','text_proxyAddress');

uicontrol('Parent',h1,'Position',[200 240 28 15],...
'Enable','off',...
'String','Port',...
'Style','text',...
'UserData','more',...
'Tag','text_proxyPort');


function mirone_pref_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

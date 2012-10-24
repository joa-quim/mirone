function varargout = aquamoto(varargin)
% A general purpose 3D grid/image file viewer, specialy taylored for tsunami (ANUGA) files
%
%	To read a file and call aquaPlugin and direct it to use a "control script"
%		aquamoto file.nc 'file_name_of_control_script'
%	To read a file and tell aquaPlugin to search the control script name in the OPTcontrol.txt file:
%		aquamoto('file.nc', 0)

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

% $Id$

% For compiling one need to include the aqua_suppfuns.m aquaPlugin.m files.

	hObject = figure('Tag','figure1','Visible','off');
	aquamoto_LayoutFcn(hObject);
	handles = guihandles(hObject);
 
	got_a_file_to_start = [];		run_aquaPlugin = false;		handles.hMirFig = [];
	if ( numel(varargin) > 0 && ~ischar(varargin{1}) )	% Expects Mirone handles as first arg
		handMir = varargin{1};
		if (numel(varargin) == 2)
			if ( exist(varargin{2}, 'file') == 2 )		% Optional file with an aquaPlugin control script
				got_a_file_to_start = varargin{2};
				handles.hMirFig = handMir.hMirFig;		% This fig already has first layer
			end
		end
		handles.home_dir = handMir.home_dir;
		handles.last_dir = handMir.last_dir;
		handles.work_dir = handMir.work_dir;
		handles.DefineEllipsoide = handMir.DefineEllipsoide;	% Potentially need in aquaPlugin
		handles.DefineMeasureUnit = handMir.DefineMeasureUnit;	%				"
		handles.IamCompiled = handMir.IamCompiled;
        d_path = handMir.path_data;
	else
		if (numel(varargin) >= 1)		% File name in input
			if ( exist(varargin{1}, 'file') == 2 )
				got_a_file_to_start = varargin{1};
			end
			if (numel(varargin) == 2)		% Run aquaPlugin after loading input file
				if (ischar(varargin{2}))
					run_aquaPlugin = varargin{2};	% Name of the control script
				else
					run_aquaPlugin = true;			% Search the control script in OPTcontrol.txt
				end
			end
		end
		handles.home_dir = cd;
		handles.last_dir = handles.home_dir;
		handles.work_dir = handles.home_dir;
		handles.DefineEllipsoide = [6378137, 0, 1/298.2572235630];	% Defaults to WGS-84
 		handles.DefineMeasureUnit = 'u';							% Defaults to 'user' units
		d_path = [handles.home_dir filesep 'data' filesep];
		% Need to know if "IamCompiled". Since that info is in Mirone handles, we need to find it out here
		try			which('mirone');			handles.IamCompiled = false;
		catch,		handles.IamCompiled = true;
		end
	end

	% -------------- Import/set icons --------------------------------------------
	load([d_path 'mirone_icons.mat'],'Mfopen_ico','Marrow_ico','um_ico','dois_ico','color_ico');
	set(handles.push_swwName,'CData',Mfopen_ico)
	ind = isnan(Marrow_ico(:,:,:));
	Marrow_ico(~ind) = 0.6;
	set(handles.push_vel,'CData',Marrow_ico)
	set(handles.toggle_1,'CData',um_ico)
	set(handles.toggle_2,'CData',dois_ico)
	set(handles.push_landPhoto,'CData',Mfopen_ico)
	set(handles.push_batGrid,'CData',Mfopen_ico)
	set(handles.push_singleWater,'CData',Mfopen_ico)
	set(handles.push_namesList,'CData',Mfopen_ico)
	set(handles.push_movieName,'CData',Mfopen_ico)
	set(handles.push_profile,'CData',Mfopen_ico)
	set(handles.push_maregs,'CData',Mfopen_ico)
	set(handles.push_palette,'CData',color_ico)
	clear Mfopen_ico Marrow_ico um_ico dois_ico color_ico;

	handles.handMir = [];		handles.fname = [];
	handles.sliceNumber = 0;	handles.one_or_zero = 1;
	handles.first = true;
	handles.volumes = [];		handles.illumComm = [];
	handles.dms_xinc = 0;		handles.dms_yinc = 0;
	handles.hQuiver = [];
	handles.geoPhoto = [];		handles.indMaxWater = [];
	handles.runinPoly = [];		handles.xyData = [];
	handles.head_bat = [];
	handles.head = [];			% Header from a sww file. POOR NAMER -> CHANGE IT
	handles.imgBat = [];		% Bathymetry only (that is, dry) image
	handles.cmapBat = [];		% To hold the color map with a discontinuity at the shoreline
	handles.geog = 0;			% For the time being ANUGA no geoga
	handles.is_sww = false;		% It will be true if a SWW file is loaded
	handles.is_coards = false;	% Updated, if the case, in aqua_suppfuns
	handles.is_otherMultiband = false;	%			"--"Show slice
	handles.minWater = [];		% Used when "scale to global min/max". If empty on call, it will be computed than
	handles.plotVector = false;	% Use to know if plot arrows
	handles.vecScale = 1;		% To scale eventual vector plots
	handles.vecSpacing = 10;	% If vector plot, plot every this interval
	handles.spacingChanged = false;		% To signal quiver when old arrows should be deleted
	handles.useLandPhoto = false;
	handles.firstLandPhoto = true;

	set(handles.popup_derivedVar, 'String', ...
		{'Absolute Velocity (V)'; ...
		'Absolute Momentum (VxD)'; ...
		'Water Depth'; ...
		'Elevation'; ...
		'Max Water'; ...
		'Max Depth'; ...
		'Froude Number'; ...
		'Velocity Head (V^2 / (2g))'; ...
		'Total Energy (Stage + V^2 / (2g))'; ...
		'Specific Energy (D + V^2 / (2g))'; ...
		'Taylor-V (D*(1+V+V**2))'; ...
		'Hazard-RVD (D(1+V^2))'; ...
		'Hazard-RVD (DV^2)'})

	% In push_swwName_CB() individual ranges MUST be assigned to handles.ranges
	% following exactly the variables order used in the popup_derivedVar.
	handles.ranges = cell(16,1);		% 16 = 3 (direct vars) + 13 (derived vars)
	handles.elevRange = [];

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
 	frames = [handles.frame1 handles.frame2 handles.frame3 handles.frame4 handles.frame5 handles.frame6 handles.frame7 handles.frame8];
	h_t = [handles.text_Pq handles.text_GLG handles.text_MovType handles.text_MovSize handles.text_testFile handles.text_globalMM handles.text_montage];
	new_frame3D(hObject, h_t, frames)
	%------------- END Pro look (3D) -----------------------------------------------------

	%------------- Resize the Figure to its final dimensions -------------------------------
	figPos = get(hObject, 'Pos');
	figPos(3:4) = [391 430];
	set(hObject, 'Pos', figPos)
	move2side(hObject,'right');

	% ------------------ TABPANEL SECTION ----------------------------------------
	% For a total unknown reason the folowing line (less the 'Parent') in tabpanelfcn would cause
	% the creation of a new figure. So I do here, where without the 'Parent', would have the same effect ???? 
	handles.blank_patch = uicontrol('Parent',hObject,'style','text', 'units','pixels','Position',[0,0,1,1]);
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
	panel_names = {'anuga','t_grids','misc','shade','cinema','plug'};
	
	% tabpanelfcn('makegroups',...) adds new fields to the handles structure,
	% one for each panel name and another called 'group_name_all'.  These fields
	% are used by the tabpanefcn when tab_group_handler is called.
	handles = tabpanelfcn('make_groups',group_name, panel_names, handles, 1);
	% ------------------------------------------------------------------------------
	
	set(hObject,'Visible','on');	drawnow

	%------ Ok, we still have a couple of things to do, but we can do them with the figure already visible
	% -- Move the uicontrols that are still 'somewhere' to their correct positions
	hTab = findobj(handles.tab_group,'UserData','anuga'); 		% Find the handle of the "ANUGA" tab push button
	handles.hTabAnuga = hTab;				% save this in case we will need it later
	
	hhs = findobj(hObject, 'UserData', 't_grids');
	hTab = findobj(handles.tab_group,'UserData','t_grids'); 	% remove the "Time Grids" tab push button from the hhs list
    hhs = setdiff(hhs, hTab);
	for (k = 1:numel(hhs))
		hhsPos = get(hhs(k), 'Pos');		hhsPos = hhsPos - [660 110 0 0];
		set(hhs(k), 'Pos', hhsPos)
	end

	hhs = findobj(hObject, 'UserData', 'misc');
	hTab = findobj(handles.tab_group,'UserData','misc'); 	% remove the "Misc" tab push button from the hhs list
	handles.hTabMisc = hTab;				% save this for we will need it later on
    hhs = setdiff(hhs, hTab);
	for (k = 1:numel(hhs))
		hhsPos = get(hhs(k), 'Pos');		hhsPos = hhsPos - [660 -160 0 0];
		set(hhs(k), 'Pos', hhsPos)
	end

	hhs = findobj(hObject, 'UserData', 'shade');
	hTab = findobj(handles.tab_group,'UserData','shade'); 	% remove the "Shading OR Image" tab push button from the hhs list
	handles.hTabShade = hTab;				% save this for we will need it later on
    hhs = setdiff(hhs, hTab);
	for (k = 1:numel(hhs))
		hhsPos = get(hhs(k), 'Pos');		hhsPos = hhsPos - [350 210 0 0];
		set(hhs(k), 'Pos', hhsPos)
	end

	hhs = findobj(hObject, 'UserData', 'cinema');
	hTab = findobj(handles.tab_group,'UserData','cinema'); 	% remove the "Cinema" tab push button from the hhs list
    hhs = setdiff(hhs, hTab);
	for (k = 1:numel(hhs))
		hhsPos = get(hhs(k), 'Pos');		hhsPos = hhsPos - [320 -80 0 0];
		set(hhs(k), 'Pos', hhsPos)
	end
	
	% Move the Plug buttons
	hhsPos = get(handles.push_plugFun, 'Pos');		hhsPos = hhsPos - [0 260 0 0];
	set(handles.push_plugFun, 'Pos', hhsPos)
	hhsPos = get(handles.check_plugFun, 'Pos');		hhsPos = hhsPos - [0 140 0 0];
	set(handles.check_plugFun, 'Pos', hhsPos)

	% ---------------- Transited from tsunamovie --------------------------------------------
	handles.dither = 'nodither';% Default
	handles.checkedMM = 0;      % To signal if need or not to get the ensemble water Min/Max
	handles.usrMM = 0;          % To signal if user has changed the ensemble Min|Max
	handles.lambCteComm = '/0.55/0.6/0.4/10';   % These Lambertian params are here const
	handles.waterIllumComm = '-E0/30/0.55/0.6/0.4/10';      % Starting values
	handles.landIllumComm = '-A0';
	handles.landCurrIllumType  = 'grdgradient_A';
	handles.waterCurrIllumType = 'lambertian';
	handles.fps = 5;            % Frames per second
	handles.dt = 0.2;           % 1/fps
	handles.scaleFactor = 1;    % To shrink or increase the movie dimensions
	handles.Z_bat   = [];		handles.Z_water = [];
	handles.nameList = [];		handles.testTime = [];
	handles.multiLayerInc = 1;
	handles.reinterpolated_bat = false; % To when we need to reinterpolate bat to fit with water
	handles.flederize = false;	% Temporary way of controling if flederize

	% The rest is donne in push_swwName_CB() because only than we have all the necessary info
	S = load([d_path 'gmt_other_palettes.mat'],'DEM_screen');
	handles.cmapLand = S.DEM_screen;
% 	S = load([d_path 'gmt_other_palettes.mat'],'Terre_Mer');
% 	handles.terraMar = S.Terre_Mer;

	% By default use a blue only colormap for water
	handles.cmapWater = [0 0 1; 0 0 1];
% 	S = load([d_path 'gmt_other_palettes.mat'],'polar');
% 	handles.cmapWater = S.polar;
	handles.cmapWater_bak = handles.cmapWater;      % Make a copy for cmaps reseting
	handles.cmapLand_bak = handles.cmapLand;
	
	% --------------- Import background image for axes1 ------------------------------------
	astrolabio = imread([d_path 'astrolabio.jpg']);
	image(astrolabio,'parent',handles.axes1);

	pos = get(handles.axes1,'Position');
	set(handles.axes1,'Visible','off','XTick',[],'YTick',[])
	
	% Draw everything that may be needed for all options. Later, depending on the
	% option selected, only the allowed features will be let visible
	x0 = pos(3)/2;      y0 = pos(4)/2;      radius = pos(3)/2 - 3;
	h_line(1) = line('parent',handles.axes1,'XData',[x0 x0],'YData',[y0 0],'Color','r','LineWidth',3,'Userdata',radius,'HitTest','off');
	% Now draw, on axes2, a quarter of circle and a line
	t = 0:0.02:pi/2;    x = [0 cos(t) 0];     y = [0 sin(t) 0];
	line('parent',handles.axes2,'XData',x,'YData',y,'HitTest','off','Color','k','LineWidth',1);
	xdataElev = [0 cos(30*pi/180)];     ydataElev = [0 sin(30*pi/180)];
	h_line(2) = line('parent',handles.axes2,'XData',xdataElev,'YData',ydataElev,'Color','k','LineWidth',3,'Visible','off');
	set(h_line(2),'Tag','Elev','Userdata',1)        % save radius of circumscribed circle

	% --------------- Backup the graphical info
	handles.landLineAzBack = [x0 x0 y0 0];
	handles.waterLineAzBack = [x0 x0 y0 0];
	handles.landLineElevBack = [xdataElev ydataElev];
	handles.waterLineElevBack = [xdataElev ydataElev];
	handles.landAzStrBack  = '0';			handles.waterAzStrBack = '0';
	handles.landElevStrBack  = '30';		handles.waterElevStrBack = '30';

	handles.ciclePar = [x0 y0 radius];
	handles.h_line = h_line;
	guidata(hObject, handles);
	show_needed(handles,'grdgradient_A')
	set(handles.edit_azim,'Visible','off')	% We need to this to invisible again because is too soon to show it up
	set(handles.toggle_1,'Value',1)         % Start it in a pressed state
	set(hObject,'WindowButtonDownFcn',{@ButtonDown,h_line,handles});

	% If we got a file in input
	if (~isempty(got_a_file_to_start))
		push_swwName_CB(handles.push_swwName, [], handles, got_a_file_to_start)
		handles = guidata(handles.figure1);		% Get updated handles
		% And now, if asked for, run the aquaPlugin
		if (run_aquaPlugin),	aquaPlugin(handles, run_aquaPlugin),	end
	else
		% This will cause a silent error but it also load the mex file in memory so it will be fast on "first" use
		try		mexnc('open', 'lixoxo', 0 );	end
	end

	guidata(hObject, handles);
	if (nargout),   varargout{1} = hObject;     end

% -------------------------------------------------------------------------------------
function tab_group_ButtonDownFcn(hObject, eventdata, handles)
% Call the tab_group_handler.  This updates visiblity of components as needed to
% hide the components from the previous tab and show components on this tab.
% This also updates the last_tab field in the handles structure to keep track
% of which panel was hidden.
	handles = tabpanelfcn('tab_group_handler',hObject, handles, get(hObject, 'Tag'));
	if ( hObject == handles.hTabShade )
		set(handles.push_bg, 'Vis', 'off')		% It has axes, which would be otherwise hidden
		set(handles.axes2,'Visible','off')		% tabpanelfcn had made them visible
	elseif ( hObject == handles.hTabMisc )
		set([handles.edit_globalWaterMin handles.edit_globalWaterMax handles.textMin handles.textMax ...
				handles.text_globalMM handles.push_bg], 'Vis', 'on');
	else
		set(handles.push_bg, 'Vis', 'on')
	end
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function tab_group_CB(hObject, eventdata, handles)
%     handles = tabpanelfcn('tab_group_handler',hObject, handles, get(hObject, 'Tag'));
% 	if ( hObject == handles.hTabShade )
% 		set(handles.push_bg, 'Vis', 'off')		% It has axes, which would be otherwise hidden
% 		set(handles.axes2,'Visible','off')		% tabpanelfcn had made them visible
% 	else
% 		set(handles.push_bg, 'Vis', 'on')
% 	end
%     guidata(hObject, handles);

% -----------------------------------------------------------------------------------------
function slider_layer_CB(hObject, eventdata, handles)
	handles.sliceNumber = round(get(handles.slider_layer,'Value')) - 1;
	set(handles.edit_sliceNumber,'String', handles.sliceNumber+1)		% Update slice n? box
	set(handles.figure1,'pointer','watch')
	if (handles.is_sww)
		push_showSlice_CB([], [], handles)		% and update image (also saves handles)
	elseif (handles.is_coards)
		aqua_suppfuns('coards_slice', handles)
	elseif (handles.is_otherMultiband)
		aqua_suppfuns('forGDAL_slice', handles)
	end
	set(handles.figure1,'pointer','arrow')

% -----------------------------------------------------------------------------------------
function edit_swwName_CB(hObject, eventdata, handles)
    fname = get(hObject,'String');
	if ( ~isempty(fname) ),    push_swwName_CB(handles.push_swwName, [], handles, fname),	end

% -----------------------------------------------------------------------------------------
function push_swwName_CB(hObject, eventdata, handles, opt)
% This function does quite some work. It reads and extract relevant info from the netCDF file

	if (nargin == 3)		% Direct call
		[FileName, PathName, handles] = put_or_get_file(handles, ...
			{'*.sww;*.SWW;*.nc;*.NC', 'Data files (*.sww,*.SWW,*.nc,*.NC)';'*.*', 'All Files (*.*)'},'sww file','get');
		if isequal(FileName,0),		return,		end
		
	else					% File name on input
		[PathName,FNAME,EXT] = fileparts(opt);
		if (~isempty(PathName)),	PathName = [PathName filesep];	end		% To be coherent with the 'if' branch
		FileName = [FNAME EXT];
	end
	pause(0.01);	handles.fname = [PathName FileName];

	if (exist(handles.fname, 'file') ~= 2)
		errordlg(['File: ' handles.fname ' does not exist.'],'Error')
		handles.fname = [];
		return
	end
	set(handles.edit_swwName,'String',handles.fname)

	% ----------------- this make sense when one are reloading a(nother) file ----------
	handles.sliceNumber = 0;
    set(handles.slider_layer,'Value',1)
    set(handles.edit_sliceNumber,'String','1')
	
	% ---- Maybe the dimensions should be fished out of the "s" structurem as well -----
	% But for now I'll just test that 'number_of_volumes' exists, otherwise ... street

	% Check if it's a netCDF file before decide what to do
	fid = fopen(handles.fname, 'r');
	ID = fread(fid,3,'*char');      ID = ID';      fclose(fid);
	if (strcmpi(ID,'CDF'))
		s = nc_funs('info',handles.fname);
	else			% Some other format. Get it opened with gdalread in aqua_suppfuns
		set(handles.check_splitDryWet,'Val',0)
		set(handles.push_showMesh,'Enable','off')
		aqua_suppfuns('forGDAL_hdr',handles)
		return
	end
	if (~isempty(s.Attribute))
		attribNames = {s.Attribute.Name};
	else
		attribNames = [];
	end
	ind = strcmp({s.Dimension.Name},'number_of_volumes');
	if (~any(ind))
		[X,Y,Z,head,misc] = nc_io(handles.fname,'R');
		if (numel(head) == 9 && isfield(misc,'z_id'))			% INTERCEPT POINT FOR PLAIN COARDS NETCDF FILES
			if (numel(misc.z_dim) <= 2)
				errordlg('This netCDF file is not 3D. Use Mirone directly to read it.','Error')
				return
			end
			handles.nc_info = s;		% Save the nc file info
			set(handles.check_splitDryWet,'Val',0)
			aqua_suppfuns('coards_hdr',handles,X,Y,head,misc)
			set(handles.push_showMesh,'Enable','off')
			return
		end
		errordlg('ERROR: This .sww file is not of recognizable type. For example: "number_of_volumes" was not found.','Error')
		set(handles.figure1,'pointer','arrow')
		return
	end

	ind = strcmp(attribNames,'xllcorner');		xllcorner = 0;
	if (any(ind)),	xllcorner = s.Attribute(ind).Value;		end
	ind = strcmp(attribNames,'yllcorner');		yllcorner = 0;
	if (any(ind)),	yllcorner = s.Attribute(ind).Value;		end
	
	% Fill the grid size boxes which imply some parameter guessings
	st = nc_funs('getdiminfo', handles.fname,'number_of_volumes');
	handles.number_of_volumes = st.Length;
	st = nc_funs('getdiminfo', handles.fname,'number_of_points');
	handles.number_of_points = st.Length;
	st = nc_funs('getdiminfo', handles.fname,'number_of_timesteps');
	handles.number_of_timesteps = st.Length;

	% ------------------ OK, Get numerics now -----------------------------------
	handles.x = double(nc_funs('varget', handles.fname, 'x')) + xllcorner;
	handles.y = double(nc_funs('varget', handles.fname, 'y')) + yllcorner;
	handles.time = nc_funs('varget', handles.fname, 'time');
	handles.volumes = nc_funs('varget', handles.fname, 'volumes');
	if (~isa(handles.volumes, 'int32')),	handles.volumes = int32(handles.volumes);	end
	set(handles.figure1,'pointer','arrow')

	head = [min(handles.x) max(handles.x) min(handles.y) max(handles.y) 0 1 0];
	set( handles.edit_x_min,'String',sprintf('%.8g',head(1)) )
	set( handles.edit_x_max,'String',sprintf('%.8g',head(2)) )
	set( handles.edit_y_min,'String',sprintf('%.8g',head(3)) )
	set( handles.edit_y_max,'String',sprintf('%.8g',head(4)) )
	handles.x_min = head(1);			handles.x_max = head(2);
	handles.y_min = head(3);			handles.y_max = head(4);
	handles.x_min_or = head(1);			handles.x_max_or = head(2);
	handles.y_min_or = head(3);			handles.y_max_or = head(4);

	% ------------------- Get/compute the global min/max --------------------------------
	varNames = {s.Dataset.Name};						% WHAT IF THIS FAILS??? 
	ind = strcmp(varNames,'stage_range');
	if (any(ind))
		handles.ranges{1} = double(nc_funs('varget', handles.fname, 'stage_range'));
		handles.minWater = handles.ranges{1}(1);
		handles.maxWater = handles.ranges{1}(2);
		set(handles.edit_globalWaterMin, 'String', handles.ranges{1}(1))
		set(handles.edit_globalWaterMax, 'String', handles.ranges{1}(2))
	end
	ind = strcmp(varNames,'xmomentum_range');
	if (any(ind))
		handles.ranges{2} = double(nc_funs('varget', handles.fname, 'xmomentum_range'));
	end
	ind = strcmp(varNames,'ymomentum_range');
	if (any(ind))
		handles.ranges{3} = double(nc_funs('varget', handles.fname, 'ymomentum_range'));
	end
	ind = strcmp(varNames,'elevation_range');
	if (any(ind))				% Not unlikely that some of the following are wrong. F... devided by (nearly) zero
		handles.elevRange = double(nc_funs('varget', handles.fname, 'elevation_range'));
		% ok, here we are going to compute the derived vars min/max as well
		D = [eps; max(handles.ranges{1} - handles.elevRange)];		% Why the hell handles.ranges{1}(2) == handles.elevRange(2) !!!!
		handles.ranges{5}  = [0; max(sqrt(handles.ranges{2} .^ 2 + handles.ranges{3} .^ 2))];	% momentumRange
		handles.ranges{6}  = D;											% waterDepthRange
		handles.ranges{4}  = handles.ranges{5} ./ D;					% vRange
		handles.ranges{10} = handles.ranges{4} ./ sqrt(D * 9.8);		% V / sqrt(gD) - froudeRange
		handles.ranges{11} = handles.ranges{4} .^ 2 / (2*9.8);			% vHeadRange
		handles.ranges{12} = [min(handles.ranges{1}); max(handles.ranges{1})] + handles.ranges{4} .^ 2 / (2*9.8);		% total energyRange
		handles.ranges{13} = D + handles.ranges{4} .^ 2 / (2*9.8);		% specific energyRange
		handles.ranges{14} = D .* (1 + handles.ranges{4} + handles.ranges{4} .^ 2);		% D(1+V+V^2) - taylorVRange
		handles.ranges{15} = D .* (1 + handles.ranges{4} .^ 2);			% D(1+V^2) -- hazard_rvd_Range
		handles.ranges{16} = D .* (handles.ranges{4} .^ 2);				% DV^2 -- hazard2_rvd_Range
	end

	% --------------- Estimate a "reasonable" proposition for grid dimensions------------------
	n = round( sqrt(double(handles.number_of_volumes)) );
	inc = ( diff(head(1:2)) + diff(head(3:4)) ) / (2*(n-1));	% A mean dx dy
	set( handles.edit_x_inc,'String',sprintf('%.8g',inc) )
	set( handles.edit_y_inc,'String',sprintf('%.8g',inc) )
	set( handles.edit_Ncols,'String',sprintf('%d',n) )
	set( handles.edit_Nrows,'String',sprintf('%d',n) )
	% Call dim_funs to compute & update the correct size
	dim_funs('xInc', handles.edit_x_inc, handles)
	dim_funs('yInc', handles.edit_y_inc, handles)

	% ----------------- Remainings ... --------------------------------------------------------
	handles.head = head;		% INCOMPLETE HEAD. The rest is computed in push_showSlice_CB()

	if ( ~isempty(handles.elevRange) )
		head(5:6) = handles.elevRange;
		if ( head(5) == head(6) ),		head(6) = head(6) + 0.1;	end		% not screw const depth cases
		handles.cmapBat = makeCmapBat(handles, head, handles.cmapLand, 1);		% Put the cmap discontinuity at the zero of bat
	else
		warndlg('Could not find elevation (bathymetry) min/max. End of the world is NEEEAAAAR','WARNING')	% SCREEEEEMMMMMMMMM
	end

	% ----------------- Finish slider configurations ------------------------------------------
	st = [1 10] / (handles.number_of_timesteps - 1);
	set(handles.slider_layer,'Min',1,'Max',handles.number_of_timesteps,'Val',1,'SliderStep',st) 	
	set(handles.slider_layer,'Enable','on')

	set(handles.edit_sliceNumber,'Enable','on')
	set(handles.text_Info,'String',sprintf('Triangles = %d & Time steps = %d',handles.number_of_volumes,handles.number_of_timesteps))

	handles.illumComm = [];					% New file. Reset illum state.
	handles.imgBat = [];
	handles.is_sww = true;
	set(handles.radio_multiLayer, 'Val', 1)
	set(handles.edit_multiLayerInc, 'Enable', 'on')
	set(handles.radio_timeGridsList,'Val',0)
	set([handles.textResize handles.popup_resize], 'Enable', 'off')

	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function push_showSlice_CB(hObject, eventdata, handles)
	if ( ~get(handles.radio_multiLayer, 'Val') && ~isempty(handles.nameList) )
		errordlg('This button is for exclusive use of ANUGA files. Not for "Time Grids" list.','Error'),	return
	end
	if ( isempty(handles.fname) )
		errordlg('Hey Lou. What about a walk on the Wild Side? Maybe you''ll find a file there that you can use here!','Chico clever')
		return
	end

	nx = str2double(get(handles.edit_Ncols,'String'));
	ny = str2double(get(handles.edit_Nrows,'String'));

	if (handles.is_coards)			% We are dealing with a coards netCDF file
		aqua_suppfuns('coards_slice', handles),		return
	elseif (handles.is_otherMultiband)
		aqua_suppfuns('forGDAL_slice', handles),	return
	elseif (handles.is_sww)
		[theVar, U, V, indVar, indWater] = get_swwVar(handles);
	end
	if (~isa(theVar, 'double')),	theVar = double(theVar);	end		% While we don't f... these doubles as well

	% ------------ Create a grid in x and y. But before see if any of the limits has changed --------------
	x = str2double(get(handles.edit_x_min,'String'));
	if (abs( (x - handles.head(1)) ) > 1e-5),		handles.head(1) = x;	end
	x = str2double(get(handles.edit_x_max,'String'));
	if (abs( (x - handles.head(2)) ) > 1e-5),		handles.head(2) = x;	end
	y = str2double(get(handles.edit_y_min,'String'));
	if (abs( (y - handles.head(3)) ) > 1e-5),		handles.head(3) = y;	end
	y = str2double(get(handles.edit_y_max,'String'));
	if (abs( (y - handles.head(4)) ) > 1e-5),		handles.head(4) = y;	end

	x = linspace(handles.head(1),handles.head(2),nx);
	y = linspace(handles.head(3),handles.head(4),ny);
	Z = mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) theVar(:)],x,y);
	% -----------------------------------------------------------------------------------------------------

	handles.head(5:6) = [min(Z(:)) max(Z(:))];
	handles.head(8) = str2double(get(handles.edit_x_inc,'String'));
	handles.head(9) = str2double(get(handles.edit_y_inc,'String'));

	splitDryWet = get(handles.check_splitDryWet,'Val');		% See if we need to build a wet and dry images, or only one

	if ( splitDryWet )
		indLand = get_landInd(handles, x, y, indWater);		% Compute Dry/Wet indexes - have to recompute, ... since water ... moves
		if ( indVar == 8 )										% Max Water option. Calculate once the Water delimiting polygon
			handles.indMaxWater = ~indLand;
		end

		if ( isempty(handles.imgBat) || size(handles.imgBat,1) ~= numel(y) || size(handles.imgBat,2) ~= numel(x) )	% First time, compute it (not shaded)
			handles.imgBat = do_imgBat(handles, indVar, x, y);	% IMG is always RGB
		end

		if (handles.useLandPhoto)
			alfa = round((1 - get(handles.slider_transparency, 'Val')) * 255);	% (1 - val) since it will be applyied to Land
			if ( ~isempty(handles.hMirFig) && ishandle(handles.handMir.hImg) )	% A Mirone figure with image already exists
				alphaMask = get(handles.handMir.hImg,'AlphaData');
				if ( numel(alphaMask) == 1 )			% No AlphaMask yet, but we'll need one further down
					alphaMask = alloc_mex(size(indLand),'uint8');	% Create an image mask of Dry/Wets
					alphaMask(~indLand) = alfa;
				else							% AlphaMask exists, but we need to update it to reflect this slice water level
					alphaMask(indLand) = 0;		alphaMask(~indLand) = alfa;
				end
			else								% The Mirone figure will be created later. Compute the AlphaMask
				alphaMask = alloc_mex(size(indLand),'uint8');	% Create an image mask of Dry/Wets
				alphaMask(~indLand) = alfa;
			end
		end
	end

	% ------------------ Open or update a Mirone window with the slice display -----------------------------
	if (isempty(handles.hMirFig) || ~ishandle(handles.hMirFig))			% First run or killed Mirone window
		tmp.X = x;		tmp.Y = y;		tmp.head = handles.head;	tmp.cmap = handles.cmapLand;
		tmp.name = sprintf('SWW time = %g', handles.time(handles.sliceNumber+1));
		handles.hMirFig = mirone(Z,tmp);
		move2side(handles.figure1,handles.hMirFig,'left')
		handles.handMir = guidata(handles.hMirFig);			% Get the handles of the now existing Mirone fig
		handles.firstLandPhoto = true;
		if ( splitDryWet && handles.useLandPhoto )
			h = image('XData',handles.geoPhotoX,'YData',handles.geoPhotoY, 'CData',handles.geoPhoto, 'Parent',handles.handMir.axes1);
			uistack(h,'bottom')
			handles.firstLandPhoto = false;
			set(handles.handMir.hImg,'AlphaData',alphaMask)	% 'alphaMask' was updated above
		end

	else									% We already have a Mirone image. Update it with this new slice
		handles.handMir = guidata(handles.hMirFig);			% Get updated handles to see if illum has changed
		setappdata(handles.handMir.figure1,'dem_z',Z);		% Update grid so that coursor display correct values
															% Have to do it here because minmax arg to scalet8 CHANGES Z

		if ( splitDryWet && handles.useLandPhoto )			% Land/Water spliting with external Land image
			if (handles.firstLandPhoto)						% First time, create the background image
				h = image('XData',handles.geoPhotoX,'YData',handles.geoPhotoY, 'CData',handles.geoPhoto, 'Parent',handles.handMir.axes1);
				uistack(h,'bottom')
				handles.firstLandPhoto = false;
			end
			set(handles.handMir.hImg,'AlphaData',alphaMask)	% 'alphaMask' was updated above

		elseif ( splitDryWet && get(handles.radio_shade, 'Val') && ...		% Land/Water spliting with Illumination
				~isequal(handles.illumComm, handles.landIllumComm) )	% Recompute 'imgBat' if illumination has changed
			handles.imgBat = do_imgBat(handles, indVar, x, y);
			handles.illumComm = handles.landIllumComm;		% save illum command for future comparison

		elseif ( ~splitDryWet )								% No Land/Water spliting
			if ( ~get(handles.check_globalMinMax, 'Val') ),		minmax = [];		% Use Slice's min/max
			else				minmax = [handles.minWater handles.maxWater];
			end
			
			if ( ~isempty(minmax) ),		img = scaleto8(Z, 8, minmax);
			else							img = scaleto8(Z);
			end

			if ( get(handles.radio_shade, 'Val') )
				img = ind2rgb8(img, handles.cmapLand);		% img is now RGB
				head = handles.head;
				if ( ~isempty(handles.ranges{indVar}) ),	head(5:6) = handles.ranges{indVar};		end
				R = illumByType(handles, Z, head, handles.landIllumComm);
				img = shading_mat(img,R,'no_scale');		% and now it is illuminated
			end
			set(handles.handMir.hImg, 'CData', img)
		end

		setappdata(handles.handMir.figure1,'dem_x',x);		% Don't get bad surprises if space increments have changed
		setappdata(handles.handMir.figure1,'dem_y',y);
		set(handles.handMir.figure1, 'Name', sprintf('SWW time = %g',handles.time(handles.sliceNumber+1)))
	end

	if ( splitDryWet )				% Get the final mixed land/water image 
		if ( handles.useLandPhoto && ~handles.firstLandPhoto )
			img = do_imgWater(handles, indVar, Z, [], indLand);		% [] means "no need to mix water/land sub-images"
		else
			img = do_imgWater(handles, indVar, Z, handles.imgBat, indLand);
		end
		set(handles.handMir.hImg, 'CData', img)				% IMG is always RGB
	end

	if ( ~isempty(U) )		% Plot vectors
		x = linspace(handles.head(1),handles.head(2), numel(1:handles.vecSpacing:nx) );
		y = linspace(handles.head(3),handles.head(4), numel(1:handles.vecSpacing:ny) );
		U = double( mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) U(:)],x,y) );
		V = double( mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) V(:)],x,y) );
		s.hAx = handles.handMir.axes1;		s.hQuiver = handles.hQuiver;	s.spacingChanged = handles.spacingChanged;
		handles.hQuiver = loc_quiver(s, x, y, U, V, handles.vecScale);
		set(handles.hQuiver(1), 'UserData', {x,y,U,V})		% Store for eventual future file saving
		cmenuHand = uicontextmenu('Parent',handles.handMir.figure1);
		set(handles.hQuiver(1), 'UIContextMenu', cmenuHand);
		set(handles.hQuiver(2), 'UIContextMenu', cmenuHand);
		uimenu(cmenuHand, 'Label', 'Save arrows', 'Call', {@save_arrows, handles.hQuiver(1)});
		if (handles.spacingChanged),	handles.spacingChanged = false;		end		% Worked once, reset it
	elseif (~isempty(handles.hQuiver))			% If we had plotted vectors, delete them
		delete(handles.hQuiver),			handles.hQuiver = [];
	end

	guidata(handles.figure1,handles)

	% Save also the updated header in Mirone handles
	handles.handMir.head = handles.head;
	guidata(handles.handMir.figure1,handles.handMir)

% -----------------------------------------------------------------------------------------
function push_runIn_CB(hObject, eventdata, handles)
% Find and plot the polygon which delimits the inundation zone

	if ( isempty(handles.hMirFig) || ~ishandle(handles.hMirFig) || isempty(handles.indMaxWater) )		% Do a lot of tricks 
		% We dont have a valid Mirone figure with data displayed. Try to create one from here.
		% But since we need "Max Water", set the popup to that before calling "push_showSlice"
		if ( isempty(handles.indMaxWater) )				% If it has not yet been computed
			set(handles.check_derivedVar, 'Val', 1)		% Simulate the all process
			check_derivedVar_CB(handles.check_derivedVar, eventdata, handles)
			set(handles.popup_derivedVar, 'Val', 5)
		else
			set(handles.check_derivedVar, 'Val', 0)		% Uncheck it to avoid a second maxwater calculation
			check_derivedVar_CB(handles.check_derivedVar, eventdata, handles)
		end
		push_showSlice_CB(handles.push_showSlice, eventdata, handles)
		handles = guidata(handles.figure1);
		% And teeeestttt again. Who knows
		if (isempty(handles.hMirFig) || ~ishandle(handles.hMirFig))
			warndlg('It is the second time (without your knowledge) that I test for a valid Mirone figure. Giving up','Warning'),	return
		end
	end
	if ( isempty(handles.indMaxWater) || ~ishandle(handles.hMirFig) )
		errordlg('Something screewed up before. It should never come here. Please inform base','Error'),	return
	end

	if ( ~isempty(handles.runinPoly) )		% We already know the Run in polygon
		h = findobj(handles.handMir.axes1, 'type', 'line', 'Tag','inundPoly');
		if ( ~isempty(h) ),		return,		end		% We already have it. Don't draw it again
		h = line('XData',handles.runinPoly(:,1), 'YData',handles.runinPoly(:,2), 'Parent',handles.handMir.axes1, ...
			'Linewidth',handles.handMir.DefLineThick, 'Color',handles.handMir.DefLineColor, 'Tag','inundPoly');
		draw_funs(h,'line_uicontext')        % Set lines's uicontextmenu
		return
	end

	% Calculate the still water mask
	nx = str2double(get(handles.edit_Ncols,'String'));
	ny = str2double(get(handles.edit_Nrows,'String'));
	x = linspace(handles.head(1),handles.head(2),nx);
	y = linspace(handles.head(3),handles.head(4),ny);
	Z = nc_funs('varget', handles.fname, 'elevation')';
	if (~isa(Z, 'double')),		Z = double(Z);		end
	Z = mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) Z(:)],x,y);
	indLand = (Z >= 0);
	indInund = (handles.indMaxWater & indLand);				% Now we have the inundation mask
	B = img_fun('bwboundaries',indInund,'noholes');
	if ( isempty(B) )
		warndlg('Could not find any inundation zone','Warning'),	return
	end
	
	y = (B{1}(:,1)-1) * handles.handMir.head(9) + handles.handMir.head(3);
	x = (B{1}(:,2)-1) * handles.handMir.head(8) + handles.handMir.head(1);
	h = line('XData',x, 'YData',y, 'Parent',handles.handMir.axes1, 'Tag','inundPoly', ...
		'Linewidth',handles.handMir.DefLineThick, 'Color',handles.handMir.DefLineColor);
	draw_funs(h,'line_uicontext')        % Set lines's uicontextmenu
	handles.runinPoly = [x y];
    guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function [theVar, U, V, indVar, indWater, qual] = get_derivedVar(handles)
	% Compute a derived quantity from a combination of the primary quantities
	% INDVAR is the index of the variable in the popup. It starts at 4 because
	% the first 3 are 'stage', 'xmomentum' and 'ymomentum'.
	% INDWATER ~= [] when THEVAR is maximum ... (e.g. Max Water)
	% QUAL is the name of the derived quantity as it appears in the popup
	%
	% 4  'Absolute Velocity (V)'; ...
	% 5  'Absolute Momentum (VxD)'; ...
	% 6  'Water Depth'; ...
	% 7  'Elevation'; ...
	% 8  'Max Water'; ...
	% 9  'Max Depth'; ...
	% 10 'Froude Number'; ...
	% 11 'Velocity Head (V^2 / (2g))'; ...
	% 12 'Total Energy (Stage + V^2 / (2g))'; ...
	% 13 'Specific Energy (D + V^2 / (2g))'; ...
	% 14 'Taylor-V (D*(1+V+V**2))'
	% 15 'Hazard-RVD (D(1+V^2))'; ...
	% 16 'Hazard-RVD (DV^2)'; ...

	U = [];		V = [];		indWater = [];
	contents = get(handles.popup_derivedVar, 'String');
	qual = contents{get(handles.popup_derivedVar,'Value')};
	switch qual(1:min(numel(qual),11))
		case {'Absolute Ve' 'Froude Numb'}			% Absolute Velocity (V) || Froude Number
			x = nc_funs('varget', handles.fname, 'xmomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			y = nc_funs('varget', handles.fname, 'ymomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			if (~isa(x, 'double')),		x = double(x);		y = double(y);		end
			if ( qual(1) == 'A' && handles.plotVector )
				U = x;		V = y;					% momentums copy
			end
			theVar = sqrt(x.^2 + y.^2);				% |M|
			x = nc_funs('varget', handles.fname, 'stage', [handles.sliceNumber 0], [1 handles.number_of_points]);
			y = nc_funs('varget', handles.fname, 'elevation')';
			if (~isa(x, 'double')),		x = double(x);		y = double(y);		end
			D = x - y + 1e-10;		clear x y;
			ind_0 = (D < 1e-8);			% To get arround a devide-by-nearly-zero and Anuga bug in velocity problem 
			theVar = theVar ./ D;
			theVar(ind_0) = 0;
			if ( qual(1) == 'A' && handles.plotVector )
				U = U ./ D;			V = V ./ D;
			end
			if ( ~isempty(U) )				% Now U,V are vx,vy
				U = U ./ D;		V = V ./ D;
				U(ind_0) = 0;	V(ind_0) = 0;
			end
			indVar = 4;
			if ( strcmp(qual(1:6), 'Froude') )		% V / sqrt(gD)
				theVar = theVar ./ sqrt( 9.8 * D );
				indVar = 10;
			end

		case 'Absolute Mo'			% Absolute Momentum (VxD)
			x = nc_funs('varget', handles.fname, 'xmomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			y = nc_funs('varget', handles.fname, 'ymomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			if (~isa(x, 'double')),		x = double(x);		y = double(y);		end
			theVar = sqrt(x.^2 + y.^2);
			if ( qual(1) == 'A' && handles.plotVector )
				U = x;		V = y;
			end
			indVar = 5;

		case 'Elevation'			% Elevation
			theVar = nc_funs('varget', handles.fname, 'elevation')';
			indVar = 7;

		case 'Water Depth'			% Water Depth
			stage = nc_funs('varget', handles.fname, 'stage', [handles.sliceNumber 0], [1 handles.number_of_points]);
			elevation = nc_funs('varget', handles.fname, 'elevation')';
			if (~isa(stage, 'double')),		stage = double(stage);		elevation = double(elevation);		end
			theVar = (stage - elevation);
			indVar = 6;

		case 'Max Water'			% Maximum water height (we need to compute it)
			aguentabar(0, 'title','Computing max water height ...', 'CreateCancelBtn');
			theVar = zeros(1, handles.number_of_points);
			stage = nc_funs('varget', handles.fname, 'elevation')';
			indWater = (stage < 0);	% Still water indices 
			for (k = 0:handles.number_of_timesteps - 1)
				stage = nc_funs('varget', handles.fname, 'stage', [k 0], [1 handles.number_of_points]);
				ind = (stage > theVar);			% ones at new max heights
				theVar(ind) = stage(ind);
				% Compute maximum inundation (runin). Do not use k = 0 because than all IND == 1, since theVar was still == 0
				if (k > 0),		indWater(ind) = 1;		end
				h = aguentabar((k+1) / handles.number_of_timesteps);
				if ( isnan(h) ),	break,	end
			end
			indVar = 8;

		case 'Max Depth'			% Maximum water depth (we need to compute it)
			aguentabar(0, 'title','Computing max water depth ...', 'CreateCancelBtn');
			theVar = zeros(1, handles.number_of_points);
			elevation = nc_funs('varget', handles.fname, 'elevation')';
			indWater = (elevation < 0);			% Still water indices 
			for (k = 0:handles.number_of_timesteps - 1)
				stage = nc_funs('varget', handles.fname, 'stage', [k 0], [1 handles.number_of_points]);
				cvlib_mex('sub', stage, elevation);
				ind = (stage > theVar);			% ones at new max depths
				theVar(ind) = stage(ind);
				% Compute maximum inundation (runin). Do not use k = 0 because than all IND == 1, since theVar was still == 0
				if (k > 0),		indWater(ind) = 1;		end
				h = aguentabar((k+1) / handles.number_of_timesteps);
				if ( isnan(h) ),	break,	end
			end
			indVar = 9;

		case {'Velocity He' 'Specific En' 'Total Energ' 'Taylor-V (D' 'Hazard-RVD '}
			% V^2 / (2g) || Stage + V^2 / (2g) || D + V^2 / (2g) || D*(1+V+V^2) || D * (1 + V^2) || D * V^2
			x = nc_funs('varget', handles.fname, 'xmomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			y = nc_funs('varget', handles.fname, 'ymomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			if (~isa(x, 'double')),		x = double(x);		y = double(y);		end
			theVar = (x.^2 + y.^2);		% = D^2 * V^2 = (DV)^2
			x = nc_funs('varget', handles.fname, 'stage', [handles.sliceNumber 0], [1 handles.number_of_points]);
			y = nc_funs('varget', handles.fname, 'elevation')';
			if (~isa(x, 'double')),		x = double(x);		y = double(y);		end
			D = (x - y + 1e-10);
			clear y;
			if (~strncmp(qual, 'To', 2)),	clear x,	end
			ind_0 = (D < 1e-8);			% To get arround a devide-by-nearly-zero and Anuga bug in velocity problem 
			if (qual(1) == 'V')
				D = D .* D;
				theVar = theVar ./ (19.6 * D);		% 19.6 = 2 * 9.8
				indVar = 11;
			elseif (qual(1) == 'T')		% Total Energy -- Stage + V^2 / (2g)
				x(x <= 0) = 0;			% We don't want negative energies
				theVar = x + theVar ./ (19.6 * (D .* D));
				indVar = 12;
			elseif (qual(1) == 'E')		% D + V^2 / (2g)
				theVar = D + theVar ./ (19.6 * (D .* D));
				indVar = 13;
			elseif (strncmp(qual, 'Ta', 2))
				% theVar = D .* (1 + sqrt(theVar ./ (D .* D)) + theVar ./ (D .* D) );	% D(1+V+V^2) 
				% theVar = D + theVar ./ D + sqrt((D .* D) .* theVar ./ (D .* D));		% D + DV^2 + sqrt(D^2 * (DV)^2 / D^2)
				theVar = D + theVar ./ D + sqrt(theVar);								% D + DV^2 + DV = D(1+V+V^2)
				indVar = 14;
			elseif (qual(1) == 'H')
				theVar = D + theVar ./ D;		% = D + (DV)^2/D = D + DV^2 = D * (1 + V^2)
				indVar = 15;
			else
				theVar = theVar ./ D;			% = (DV)^2/D = DV^2
				indVar = 16;
			end
			theVar(ind_0) = 0;

		case 'Shear Stress'						% ....
			
	end

% -----------------------------------------------------------------------------------------
function img = do_imgWater(handles, indVar, Z, imgBat, indLand)
% Compute Water (wet) image and illuminate it if that's the case
%	INDVAR -> indice of the currently processing variable (either a direct or derived var)
%
%	Output IMG is always RGB

% 	minmax = handles.ranges{indVar};		minmax = minmax(:)';
	minmax = [handles.minWater handles.maxWater];		% Despite the name it respects the "indVar"
	if ( ~get(handles.check_globalMinMax, 'Val') ),		minmax = [];	end		% Use Slice's min/max
	if ( ~isempty(minmax) ),		imgWater = scaleto8(Z, 8, minmax);
	else							imgWater = scaleto8(Z);
	end

	handles.handMir = guidata(handles.hMirFig);			% Get updated handles to see if illum has changed
	if ( handles.handMir.Illumin_type >= 1 && handles.handMir.Illumin_type <= 4 )
		illumComm = handles.waterIllumComm;
		pal = get(handles.handMir.figure1,'Colormap');
		if ( get(handles.check_splitDryWet, 'Val') ),	pal = handles.cmapWater;	end
		imgWater = ind2rgb8(imgWater, pal);						% image is now RGB
		R = illumByType(handles, Z, handles.head, illumComm);
		imgWater = shading_mat(imgWater,R,'no_scale');			% and now it is illuminated
		clear R;
	end
	if ( ndims(imgWater) == 2 ),	imgWater = ind2rgb8(imgWater,handles.cmapWater);	end		% Like promissed above
	alfa = get(handles.slider_transparency, 'Val');
	if ( ~isempty(imgBat) )
		img = mixe_images(handles, imgBat, imgWater, indLand, alfa);
	else							% When we are using a Land external image the spliting is done via an AlphaData
		img = imgWater;
	end

% ----------------------------------------------------------------------------------------------	
function imgBat = do_imgBat(handles, indVar, x, y)	
% Compute Land (dry) image and illuminate it if that's the case
%	INDVAR -> indice of the currently processing variable (either a direct or derived var)
%	X,Y regular coordinates to interpolate the triangles into regular grid
%
%	Output IMGBAT is always RGB

	bat = nc_funs('varget', handles.fname, 'elevation')';
	if (~isa(bat,'double') ),		bat = double(bat);	end
	bat = mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) bat(:)], x, y);

	imgBat = scaleto8(bat);
    imgBat = ind2rgb8(imgBat, handles.cmapBat);
	if ( get(handles.radio_shade, 'Val') )
		illumComm = handles.landIllumComm;
		head = handles.head;
		if ( ~isempty(handles.ranges{indVar}) ),	head(5:6) = handles.ranges{indVar};		end
		R = illumByType(handles, bat, head, illumComm);
		imgBat = shading_mat(imgBat, R, 'no_scale');
	end

% -----------------------------------------------------------------------------------------
function R = illumByType(handles, Z, head, illumComm)
% Compute the illuminance matrix in function of the illumination type

	if ( get(handles.toggle_1, 'Val') )
		if (handles.geog),  R = grdgradient_m(Z,head,'-M',illumComm,'-Nt');
		else                R = grdgradient_m(Z,head,illumComm,'-Nt');
		end
	else
		R = grdgradient_m(Z,head,illumComm);
	end

% -----------------------------------------------------------------------------------------
function ind = get_landInd(handles, x, y, indWater)
% Compute indices such that 1 -> Dry; 0 -> Wet
% Note that "Land" (stage - elevation) is gridded
% INDWATER, if not empty, has ones at the indices of a max quantity areal extent
% The minimalist form (first one tried) is inherited from 'tsunamovie'

	if (nargin == 2)
		zBat = handles;		zWater = x;				% Just for human understanding
		dife = cvlib_mex('absDiff',zBat,zWater);
		ind = (dife < 1e-3);						% The 1e-3 SHOULD be parameterized
		return
	end

	if (nargin == 3),	indWater = [];		end
	if ( isempty(indWater) )		% Compute Dry/Wet indices for a specific time slice
		stage = nc_funs('varget', handles.fname, 'stage', [handles.sliceNumber 0], [1 handles.number_of_points]);
		elevation = nc_funs('varget', handles.fname, 'elevation')';
		dife = cvlib_mex('absDiff',stage,elevation);
		if (~isa(dife,'double') ),		dife = double(dife);	end
		dife = mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) dife(:)],x,y);
		ind = (dife < 1e-7);
	else							% A Max quantity. INDWATER has ones at the indices of that quantity areal extent
		indWater = mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) double(indWater(:))],x,y);
		ind = (indWater < 0.5);		% Test with < 0.5 since we dont't know if mxgridtrimesh has changed the ones
	end

% -----------------------------------------------------------------------------------------
function imgWater = mixe_images(handles, imgBat, imgWater, ind, alfa)
% Mixes land and water images simulating transparency.
% It also resizes the image if the scale factor is ~= 1
    
	try
		if (alfa > 0.01)    % Only if transparency is greater than 1%
			cvlib_mex('addweighted',imgWater,(1 - alfa),imgBat,alfa);     % In-place
		end

		tmpW = imgWater(:,:,1);     tmpB = imgBat(:,:,1);           % R
		tmpW(ind) = tmpB(ind);      imgWater(:,:,1) = tmpW;

		tmpW = imgWater(:,:,2);     tmpB = imgBat(:,:,2);           % G
		tmpW(ind) = tmpB(ind);      imgWater(:,:,2) = tmpW;

		tmpW = imgWater(:,:,3);     tmpB = imgBat(:,:,3);           % B
		tmpW(ind) = tmpB(ind);      imgWater(:,:,3) = tmpW;

		if (handles.scaleFactor ~= 1)
			imgWater = cvlib_mex('resize',imgWater,handles.scaleFactor);
		end
	catch
		errordlg(['Aquamoto:mixe_images ' lasterr],'Error')
	end

% -----------------------------------------------------------------------------------------
function new_cmap = makeCmapBat(handles, head, cmap, orig)
% Put the cmap discontinuity at the zero of bathymetry (coastline)

	if (~orig),		new_cmap = cmap;	return,		end		% Untill I know better what to do

	% Isto assume que a bat tem partes neg e pos (testar)
	% cmap = handles.terraMar;
	% ind_old = 147;			% Discontinuity in the handles.terraMar cmap
	cmap(end-4:end,:) = [];		% Remove last five colors (nearly white)
	cmap = [repmat([196 156 104]/255,5,1); cmap];		% Add a yelowish color
	ind_old = 5;				% New discontinuity in cmap

	nc = length(cmap);
	z_inc = (head(6) - head(5)) / (nc - 1);
	ind_c = round(abs(0 - head(5)) / z_inc + 1);

	nl = ind_old;		nu = ind_c;
	new_cmap_l = interp1(linspace(0,1,nl), cmap(1:nl,:), linspace(0,1,nu));
	new_cmap_u = interp1(linspace(0,1,nc-ind_old), cmap(ind_old+1:nc,:), linspace(0,1,nc-ind_c));
	new_cmap = [new_cmap_l; new_cmap_u];

% -------------------------------------------------------------------------------------
function edit_x_min_CB(hObject, eventdata, handles)
	dim_funs('xMin', hObject, handles)

% -------------------------------------------------------------------------------------
function edit_x_max_CB(hObject, eventdata, handles)
	dim_funs('xMax', hObject, handles)

% --------------------------------------------------------------------
function edit_y_min_CB(hObject, eventdata, handles)
	dim_funs('yMin', hObject, handles)

% --------------------------------------------------------------------
function edit_y_max_CB(hObject, eventdata, handles)
	dim_funs('yMax', hObject, handles)

% --------------------------------------------------------------------
function edit_x_inc_CB(hObject, eventdata, handles)
	dim_funs('xInc', hObject, handles)

% --------------------------------------------------------------------
function edit_Ncols_CB(hObject, eventdata, handles)
	dim_funs('nCols', hObject, handles)

% --------------------------------------------------------------------
function edit_y_inc_CB(hObject, eventdata, handles)
	dim_funs('yInc', hObject, handles)

% --------------------------------------------------------------------
function edit_Nrows_CB(hObject, eventdata, handles)
	dim_funs('nRows', hObject, handles)

% -----------------------------------------------------------------------------------------
function push_Help_CB(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------------
function push_showMesh_CB(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------------
function edit_sliceNumber_CB(hObject, eventdata, handles)
	xx = fix(str2double(get(hObject,'String')));		% Make sure its an int
	if (isnan(xx) || xx < 1 || xx > handles.number_of_timesteps)
		handles.sliceNumber = 0;		set(hObject,'String','1')
	else
		set(hObject,'String',xx);		set(handles.slider_layer,'Val',xx)		% Update slider
		handles.sliceNumber = xx - 1;
	end
	set(handles.figure1,'pointer','watch')
	push_showSlice_CB([], [], handles)		% Update image (also saves handles)
	set(handles.figure1,'pointer','arrow')

% -----------------------------------------------------------------------------------------
function radio_stage_CB(hObject, eventdata, handles)
	if (get(hObject,'Value'))
		set([handles.radio_xmoment handles.radio_ymoment], 'Value', 0)
		handles.minWater = handles.ranges{1}(1);		handles.maxWater = handles.ranges{1}(2);
		set(handles.edit_globalWaterMin, 'String', handles.minWater)
		set(handles.edit_globalWaterMax, 'String', handles.maxWater)
		guidata(handles.figure1, handles)
	else
		set(hObject,'Value', 1)
	end

% -----------------------------------------------------------------------------------------
function radio_xmoment_CB(hObject, eventdata, handles)
	if (get(hObject,'Value'))
		set([handles.radio_stage handles.radio_ymoment], 'Value', 0)
		handles.minWater = handles.ranges{2}(1);		handles.maxWater = handles.ranges{2}(2);
		set(handles.edit_globalWaterMin, 'String', handles.minWater)
		set(handles.edit_globalWaterMax, 'String', handles.maxWater)
		guidata(handles.figure1, handles)
	else
		set(hObject,'Value', 1)
	end

% -----------------------------------------------------------------------------------------
function radio_ymoment_CB(hObject, eventdata, handles)
	if (get(hObject,'Value'))
		set([handles.radio_stage handles.radio_xmoment], 'Value', 0)
		handles.minWater = handles.ranges{3}(1);		handles.maxWater = handles.ranges{3}(2);
		set(handles.edit_globalWaterMin, 'String', handles.minWater)
		set(handles.edit_globalWaterMax, 'String', handles.maxWater)
		guidata(handles.figure1, handles)
	else
		set(hObject,'Value', 1)
	end

% -----------------------------------------------------------------------------------------
function popup_derivedVar_CB(hObject, eventdata, handles)
	contents = get(handles.popup_derivedVar, 'String');
	val = get(handles.popup_derivedVar,'Value');
	qual = contents{val};
	ico = get(handles.push_vel,'CData');
	ind = ~isnan(ico(:,:,:));
	switch qual(1:min(numel(qual),11))
		case {'Absolute Ve' 'Absolute Mo'}			% Absolute Velocity || Momentum
			set(handles.push_vel,'Enable', 'on'),	ico(ind) = 0.0;
		otherwise
			set(handles.push_vel,'Enable', 'off'),	ico(ind) = 0.6;
	end
	set(handles.push_vel,'CData',ico)
	if ( ~isempty(handles.ranges{val}) )
		handles.minWater = handles.ranges{val}(1);		handles.maxWater = handles.ranges{val}(2);
		set(handles.edit_globalWaterMin, 'String', handles.minWater)
		set(handles.edit_globalWaterMax, 'String', handles.maxWater)
	else
		set([handles.edit_globalWaterMin handles.edit_globalWaterMax], 'String', '')
		handles.minWater = [];		handles.maxWater = [];
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function check_derivedVar_CB(hObject, eventdata, handles)
	ico = get(handles.push_vel,'CData');
	if (get(hObject,'Value'))
		set([handles.radio_stage handles.radio_xmoment handles.radio_ymoment], 'Enable', 'off')
		set(handles.popup_derivedVar, 'Enable', 'on')
		val = get(handles.popup_derivedVar,'Val');
		if (val == 1 || val == 2)		% Activate the vector controls
			ind = ~isnan(ico(:,:,:));		ico(ind) = 0.0;
			set(handles.push_vel,'Enable', 'on', 'CData',ico)
		end
	else
		set([handles.radio_stage handles.radio_xmoment handles.radio_ymoment], 'Enable', 'on')
		set(handles.popup_derivedVar, 'Enable', 'off')
		ind = ~isnan(ico(:,:,:));		ico(ind) = 0.6;
		set(handles.push_vel,'Enable', 'off', 'CData',ico)
	end

% -----------------------------------------------------------------------------------------
function push_vel_CB(hObject, eventdata, handles)
	vector_plot(handles.figure1,handles.vecScale, handles.vecSpacing, handles.plotVector);

% -----------------------------------------------------------------------------------------
function toggle_1_CB(hObject, eventdata, handles)
	show_needed(handles,'grdgradient_A')
	if (get(handles.radio_land,'Value')),   handles.landCurrIllumType = 'grdgradient_A';
	else                                    handles.waterCurrIllumType = 'grdgradient_A';
	end
	ButtonUp([],[],handles.h_line,handles)		% This call updates the handles. water|land IllumComm

% -----------------------------------------------------------------------------------------
function toggle_2_CB(hObject, eventdata, handles)
	show_needed(handles,'lambertian')
	if (get(handles.radio_land,'Value')),   handles.landCurrIllumType = 'lambertian';
	else                                    handles.waterCurrIllumType = 'lambertian';
	end
	ButtonUp([],[],handles.h_line,handles)		% This call updates the handles. water|land IllumComm

% -----------------------------------------------------------------------------------------
function radio_land_CB(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end

	set(handles.radio_water,'Value',0)
	set(handles.h_line(1),'XData',handles.landLineAzBack(1:2),'YData',handles.landLineAzBack(3:4))
	set(handles.h_line(2),'XData',handles.landLineElevBack(1:2),'YData',handles.landLineElevBack(3:4))
	show_needed(handles,handles.landCurrIllumType)
	set(handles.edit_azim,'String', handles.landAzStrBack)
	set(handles.edit_elev,'String', handles.landElevStrBack)
	str = sprintf(['Set parametrs for Land illumination\nCurrent selection (in GMT parlance) is:\n',...
		handles.landIllumComm]);
	set(hObject,'TooltipString',str)
	set(handles.push_palette,'TooltipString','Choose another Land palette')

% -----------------------------------------------------------------------------------------
function radio_water_CB(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end

	set(handles.radio_land,'Value',0)
	set(handles.h_line(1),'XData',handles.waterLineAzBack(1:2),'YData',handles.waterLineAzBack(3:4))
	set(handles.h_line(2),'XData',handles.waterLineElevBack(1:2),'YData',handles.waterLineElevBack(3:4))
	show_needed(handles,handles.waterCurrIllumType)
	set(handles.edit_azim,'String', handles.waterAzStrBack)
	set(handles.edit_elev,'String', handles.waterElevStrBack)
	str = sprintf(['Set parametrs for Water illumination\nCurrent selection (in GMT parlance) is:\n',...
		handles.waterIllumComm]);
	set(hObject,'TooltipString',str)
	set(handles.push_palette,'TooltipString','Choose another Water palette')

% -----------------------------------------------------------------------------------------
function push_palette_CB(hObject, eventdata, handles)
% Get the new color map and assign it to either Land or Water cmaps
	if ( ~get(handles.check_splitDryWet, 'Val') && get(handles.radio_water,'Val') )
		warndlg('Without the "Split Dry/Wet" option selected, what you asked has no effect.','Warning'),	return
	end
	handles.no_file = 1;		% We don't want color_palettes trying to update a ghost image (DON'T KILL THIS LINE)

	cmap = color_palettes(handles);
	if (~isempty(cmap))
		if (get(handles.radio_land,'Val'))
			handles.cmapLand = cmap;        % This copy will be used if user loads another bat grid
			if ( ~isempty(handles.head_bat) )		head = handles.head_bat;	% ---> DESENRASQUE DA TRETA
			else									head = handles.head;
			end
			handles.cmapBat = makeCmapBat(handles, head, cmap, 1);
			handles.imgBat = [];			% Force recomputing on next call
		else
			handles.cmapWater = cmap;
		end
		set(handles.check_resetCmaps,'Enable','on')
		guidata(handles.figure1,handles)
	end

% -----------------------------------------------------------------------------------------
function check_resetCmaps_CB(hObject, eventdata, handles)
% Reset to the default cmaps and make this ui invisible
	if (get(hObject,'Value'))
		if ( ~isempty(handles.head_bat) )		head = handles.head_bat;	% ---> DESENRASQUE DA TRETA
		else									head = handles.head;
		end
		handles.cmapBat = makeCmapBat(handles, head, handles.cmapLand_bak, 1);    % Put the discontinuity at the zero of bat
		handles.cmapWater = handles.cmapWater_bak;
		handles.imgBat = [];			% Force recomputing on next call
		set(hObject, 'Val',0, 'Enable','off')
		guidata(handles.figure1,handles)
	end

% -----------------------------------------------------------------------------------------
function radio_noShade_CB(hObject, eventdata, handles)
	if ( ~get(hObject,'Val'))
		set(hObject,'Val', 1)
		return
	end
	set(handles.radio_shade, 'Val', 0)
	set([handles.edit_elev handles.edit_azim handles.toggle_1 handles.toggle_2], 'Enable', 'off')
	set(handles.h_line, 'HitTest', 'off')

% -----------------------------------------------------------------------------------------
function radio_shade_CB(hObject, eventdata, handles)
	if ( ~get(hObject,'Val'))
		set(hObject,'Val', 1)
		return
	end
	set(handles.radio_noShade, 'Val', 0)
	set([handles.edit_azim handles.toggle_1 handles.toggle_2], 'Enable', 'on')
	set(handles.h_line, 'HitTest', 'on')
	if ( ~get(handles.toggle_1, 'Val') )
		set(handles.edit_elev, 'Enable', 'on')
	end

% -----------------------------------------------------------------------------------------
function push_apply_CB(hObject, eventdata, handles)
	% The job to be done is the same as in "Show slice"
	if ( ishandle(handles.hMirFig) )
		push_showSlice_CB(handles.push_showSlice, eventdata, handles)
	end

% -----------------------------------------------------------------------------------------
function radio_gif_CB(hObject, eventdata, handles)
	if (get(hObject,'Value')),      set([handles.radio_avi handles.radio_mpg],'Value',0)
	else                            set(hObject,'Value',1)
	end
	mname = get(handles.edit_movieName,'String');
	if (~isempty(mname))
		mname = [handles.moviePato handles.movieName '.gif'];
		set(handles.edit_movieName,'String',mname)
	end

% -----------------------------------------------------------------------------------------
function radio_avi_CB(hObject, eventdata, handles)
	if (get(hObject,'Value')),      set([handles.radio_gif handles.radio_mpg],'Value',0)
	else                            set(hObject,'Value',1)
	end
	mname = get(handles.edit_movieName,'String');
	if (~isempty(mname))
		mname = [handles.moviePato handles.movieName '.avi'];
		set(handles.edit_movieName,'String',mname)
	end
    
% -----------------------------------------------------------------------------------------
function radio_mpg_CB(hObject, eventdata, handles)
	if (get(hObject,'Value')),		set([handles.radio_avi handles.radio_gif],'Value',0)
	else							set(hObject,'Value',1)
	end
	mname = get(handles.edit_movieName,'String');
	if (~isempty(mname))
		mname = [handles.moviePato handles.movieName '.mpg'];
		set(handles.edit_movieName,'String',mname)
	end

% -----------------------------------------------------------------------------------------
function check_dither_CB(hObject, eventdata, handles)
	if (get(hObject,'Value')),		handles.dither = 'dither';
	else							handles.dither = 'nodither';
	end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_fps_CB(hObject, eventdata, handles)
% Frames per second
	fps = round(str2double(get(hObject,'String')));
	if (isnan(fps))
		set(hObject,'String',num2str(handles.fps)),		return
	end
	set(hObject,'String',num2str(fps))      % In case there were decimals
	handles.fps = fps;
	handles.dt = 1/fps;
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function popup_resize_CB(hObject, eventdata, handles)
	contents = get(hObject,'String');
	handles.scaleFactor = str2double(contents{get(hObject,'Value')});
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_movieName_CB(hObject, eventdata, handles)
	fname = get(hObject,'String');
	push_movieName_CB([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_movieName_CB(hObject, eventdata, handles, opt)
	if (nargin == 3)        % Direct call
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.gif;*.avi;*.mpg', 'Grid files (*.gif,*.avi,*.mpg)'},'Select Movie name','put');
		if isequal(FileName,0),		return,		end
		[dumb,FNAME,EXT]= fileparts(FileName);
		fname = [PathName FileName];
	else        % File name on input
		[PathName,FNAME,EXT] = fileparts(opt);
		fname = opt;
		if (~isempty(PathName)),	PathName = [PathName filesep];	end		% To be coherent with the 'if' branch
	end
	if ( ~any(strcmpi(EXT,{'.gif' '.avi' '.mpg' '.mpeg'})) )
		errordlg('Ghrrrrrrrr! Don''t be smart. Only ''.gif'', ''.avi'', ''.mpg'' or ''mpeg'' extensions are acepted.', ...
            'Chico Clever');
		return
	end
	set(handles.edit_movieName, 'String', fname)
	
	handles.moviePato = PathName;			handles.movieName = FNAME;
	if (strcmpi(EXT,'.gif'))
		set(handles.radio_gif,'Value',1)
		radio_gif_CB(handles.radio_gif, [], handles)
	elseif (strcmpi(EXT,'.avi'))
		set(handles.radio_avi,'Value',1)
		radio_avi_CB(handles.radio_avi, [], handles)
	else
		set(handles.radio_mpg,'Value',1)
		radio_mpg_CB(handles.radio_mpg, [], handles)
	end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_multiLayerInc_CB(hObject, eventdata, handles)
% Get the increment for movie of a multi-layer grid.
	str = get(hObject, 'String');
	xx = round( str2double(str) );
	if ( isnan(xx) )		% As will be the case if we have a start:inc:stop
		ind = strfind(str,':');			% Search for a start:inc:end form
		if (numel(ind) == 2)
			start = round( str2double( str(1:(ind(1)-1)) ) );
			inc   = round( str2double( str((ind(1)+1):(ind(2)-1)) ) );
			fim   = round( str2double( str((ind(2)+1):end) ) );
			if ( isnan(fim) && strcmp( str((ind(2)+1):end), 'end') ),	fim = handles.number_of_timesteps;	end
		elseif (numel(ind) == 1)		% Only start:stop
			start = round( str2double( str(1:(ind(1)-1)) ) );
			fim   = round( str2double( str((ind(1)+1):end) ) );
			if ( isnan(fim) && strcmp( str((ind(1)+1):end), 'end') ),	fim = handles.number_of_timesteps;	end
			inc = 1;
		end
		try
			handles.multiLayerInc = start:inc:fim;
		catch
			set(hObject, 'String', 1),			handles.multiLayerInc = 1;
		end
	elseif ( xx < 1)
		set(hObject, 'String', 1)
		handles.multiLayerInc = 1:handles.number_of_timesteps;
	else
		handles.multiLayerInc = 1:xx:handles.number_of_timesteps;
	end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function radio_multiLayer_CB(hObject, eventdata, handles)
% Movie selection button. Make movie from a netCDF file(s)
	if ( ~get(hObject, 'Val') ),		set(hObject, 'Val', 1),		return,		end
	set(handles.radio_timeGridsList, 'Val', 0)
	set(handles.edit_multiLayerInc, 'Enable', 'on')

% -----------------------------------------------------------------------------------------
function radio_timeGridsList_CB(hObject, eventdata, handles)
% Movie selection button. Make movie from a list of files
	if ( ~get(hObject, 'Val') ),		set(hObject, 'Val', 1),		return,		end
	set(handles.radio_multiLayer, 'Val', 0)
	set(handles.edit_multiLayerInc, 'Enable', 'off')

% -----------------------------------------------------------------------------------------
function push_OK_CB(hObject, eventdata, handles, opt)
% Do whatever it has to do (and that is a lot) to make a movie

	is_gif = get(handles.radio_gif,'Value');
	is_avi = get(handles.radio_avi,'Value');
	is_mpg = get(handles.radio_mpg,'Value');
	is_montage = false;
	do_logo = false;				% Insert a logo image but it needs to be manually positioned in fig
	if (nargin == 4)
		is_montage = true;		is_gif = false;		is_avi = false;		is_mpg = false;
	end

	% ------------------------- Do a movie (if so) from a "packed" netcdf file ------------------------
	if ( get(handles.radio_multiLayer, 'Val') )
		if ( isempty(handles.fname) )
			errordlg('Sorry, I can''t make movies out of a ... NOTHIIIINGGGG. Percebes???','ERROR');
			return
		end
		if ( numel(handles.multiLayerInc) == 1 )		% A const increment. Build the slice vector
			slice_vec = 1:handles.multiLayerInc:handles.number_of_timesteps;
		else
			slice_vec = handles.multiLayerInc;			% Increments were already in a vector form
		end

		aguentabar(0,'title','Relax and wait ... I''m flooding','createcancelbtn');
		n_slices = numel(slice_vec);
		M.cdata = [];		M.colormap = [];
		for (i = 1:n_slices )
			handles.sliceNumber = slice_vec(i) - 1;
			if (handles.is_sww)
				push_showSlice_CB([], [], handles)	% and update image (also saves handles)
			elseif (handles.is_coards)
				aqua_suppfuns('coards_slice', handles)
			elseif (handles.is_otherMultiband)
				aqua_suppfuns('forGDAL_slice', handles)
			else
				errordlg('What kind of file is this?','ERROR')
				return
			end
			handles = guidata(handles.figure1);			% We need to update for what changed in push_showSlice_CB()
			hAguenta = aguentabar(i / n_slices);		% Show visualy the processing advance
			if (isnan(hAguenta)),	break,		end		% User hit cancel
			if (ishandle(hAguenta)),	figure(hAguenta),	end		% Make sure the aguentabar is always on top
			if ( ~handles.useLandPhoto )
				img = get(handles.handMir.hImg, 'CData');
			else
				%img = imcapture(handles.handMir.axes1,'img',150);		% Do a screen capture
				F = getframe(handles.handMir.axes1);		img = F.cdata;
				if ( strcmp(get(handles.handMir.axes1,'Ydir'),'normal') ),	img = flipdim(img,1);	end
			end
			if (ndims(img) == 3),	map = [];
			else					map = get(handles.handMir.figure1,'Colormap');
			end
			[M, map] = aux_movie(handles, is_gif, is_mpg, is_avi, img, i, M, map);
		end
		if (is_avi)
			mname = [handles.moviePato handles.movieName '.avi'];
			movie2avi_j(M,mname,'compression','none','fps',handles.fps)
		elseif (is_mpg)
			mname = [handles.moviePato handles.movieName '.mpg'];
			opt = [1, 0, 1, 0, 10, 5, 5, 5];
			mpgwrite(M,map,mname,opt)
		elseif (is_montage)
			montage(M)
		end
		return
	end
	% --------------------------- END MOVIE FROM netCDF PACK FILE (e.g SWW) section ------------------

	% --------------------------- Input came from a list of grids ------------------------------------
	if (isempty(handles.Z_bat))
		errordlg('Noooo! Where is the bathymetry file? Do you think I''m bruxo?','ERROR');  return
	end

	% 'surface elevation' and 'water depth' grids are treated diferently
	is_surfElev = (get(handles.popup_surfType,'Val') == 1);

	geog = guessGeog(handles.head_bat(1:4));	% Guess if grids are geogs

	if (isempty(handles.Z_water))				% We don't have a testing grid
		if (isempty(handles.nameList))			% Neither water list nor single water grid
			errordlg('Where is the water to make the WaterWorld movie?','ERROR'),	return
		end
		[handles,X,Y,Z_water,handles.head_water] = read_gmt_type_grids(handles,handles.nameList{1});
	end

	% See if we need to (and can) reinterpolate the bat to go along with the water grids
	if (~handles.reinterpolated_bat && (any(handles.head_water(1:4) - handles.head_bat(1:4)) || ...
			( ~isequal( size(handles.Z_bat), size(handles.Z_water)) && ~isequal( size(handles.Z_bat), size(Z_water))) ) )

		h = warndlg('Ai, Ai, Bathymetry and Water grids are not compatible. Trying to fix that ...','Warning');
		opt_R = sprintf( '-R%.12f/%.12f/%.12f/%.12f', handles.head_water(1:4) );
		opt_I = sprintf( '-I%.12f/%.12f',handles.head_water(8:9) );
		handles.Z_bat = grdsample_m(handles.Z_bat,handles.head_bat,opt_R,opt_I);
		handles.head_bat = handles.head_water;
		handles.head_bat(5) = double(min(handles.Z_bat(:)));
		handles.head_bat(6) = double(max(handles.Z_bat(:)));
		% We reach here when reinterpolation was sucesseful
		h_txt = findobj(h,'Type','text');
		txt = get(h_txt, 'String');
		if (~iscell(txt)),  txt = {txt};    end
		txt{2} = '...';
		txt{3} = '... Luckyly for you that I am so good. We can proceed.';
		set(h_txt, 'String', txt)
		handles.reinterpolated_bat = true;
		guidata(handles.figure1, handles)
	end

	alfa = get(handles.slider_transparency,'Value');
	head_struct.head = [handles.head_bat(1:4) 0 255 handles.head_bat(7:9)];
	head_struct.X = handles.X_bat;		head_struct.Y = handles.Y_bat;
	head_struct.geog = geog;			head_struct.name = 'Tsu frame';

	% Do ground illum
	imgBat = ind2rgb8(scaleto8(handles.Z_bat), handles.cmapBat);
	opt_M = ' ';
	if (geog),		opt_M = '-M';	end
	R = grdgradient_m(handles.Z_bat,handles.head_bat,opt_M,handles.landIllumComm,'-Nt');
	imgBat = shading_mat(imgBat,R,'no_scale');

    % ------------------- If we have a testing grid, work on it and return -------------------------- 
	if ( ~isempty(handles.Z_water) )
		imgWater = ind2rgb8(scaleto8(handles.Z_water), handles.cmapWater);
		R = grdgradient_m(handles.Z_water,handles.head_water,handles.waterIllumComm);
		imgWater = shading_mat(imgWater,R,'no_scale');    	clear R;

		% ------ Compute indices of Land
		if (is_surfElev),	indLand = (handles.Z_bat >= 0);		% Surface height (not water depth)
		else				indLand = (handles.Z_water == 0);
		end
		imgWater = mixe_images(handles, imgBat, imgWater, indLand, alfa);
		if (~isempty(handles.strTimes))
			imgWater = cvlib_mex('text',imgWater,handles.testTime,[10 30]);
		end
		hFig = mirone(imgWater, head_struct);
		handMir = guidata(hFig);
		if ( handles.useLandPhoto )
			h = image('XData',handles.geoPhotoX,'YData',handles.geoPhotoY, 'CData',handles.geoPhoto, 'Parent',handMir.axes1);
			uistack(h,'bottom')			% Send to bottom because we want the alphamask applyied to data derived image
			alphaMask = alloc_mex(size(indLand),'uint8');	% Create an image mask of Dry/Wets
			alfa = round((1 - get(handles.slider_transparency, 'Val')) * 255);	% (1 - val) since it will be applyied to Land
			alphaMask(~indLand) = alfa;
			set(handMir.hImg,'AlphaData',alphaMask)
		end
		return
	end
	% ---------------------------------------------------------------------------------------------

	% ---------------------------- If we have a list of names -------------------------------------
	if ( ~isempty(handles.nameList) )
		if (isempty(handles.movieName))
			errordlg('Hei! what shoult it be the movie name?','ERROR'),		return
		end

		nGrids = numel(handles.nameList);
		if (~handles.checkedMM)				% We don't know yet the water ensemble Min/Max
			[minWater, maxWater,heads] = push_checkGlobalMM_CB([], [], handles);
		else								% Get what we know
			minWater = handles.minWater;
			maxWater = handles.maxWater;
			heads = handles.gridHeaders;
		end

 		if (do_logo),	logo = flipdim( imread('C:\SVN\mironeWC\data\mirone.tif'), 1);	end

		if ( ~handles.useLandPhoto )		% Otherwise the bar is either hiden by the Mirone window, or screws the transparency  
			aguentabar(0,'title','Please wait ... I''m flooding','createcancelbtn');
		end
		
		% ------------ Compute indices of Land
		if (is_surfElev),		indLand = (handles.Z_bat >= 0);		end		% Surface height (not water depth)

		M = struct('cdata',[], 'colormap',[]);
		for (i = 1:nGrids)
			if ( ~handles.useLandPhoto )
				hAguenta = aguentabar(i / nGrids);			% Show visualy the processing advance
				if (isnan(hAguenta)),	break,		end		% User hit cancel
			end

			[handles,X,Y,Z] = read_gmt_type_grids(handles,handles.nameList{i});
			imgWater = ind2rgb8(scaleto8(Z,8,[minWater maxWater]), handles.cmapWater);	% Would change Z if min/max were not global
			R = grdgradient_m(Z,heads{i},handles.waterIllumComm);
			imgWater = shading_mat(imgWater,R,'no_scale');

			% ------------ Compute indices of Land
			if (~is_surfElev),		indLand = (Z == 0);		end		% Water depth
			imgWater = mixe_images(handles, imgBat, imgWater, indLand, alfa);
			if (~isempty(handles.strTimes))
				cvlib_mex('text',imgWater,handles.strTimes{i},[10 30]);
			end

			if (do_logo)
				imgWater(30:29+size(logo,1),670:669+size(logo,2),:) = logo;
			end
			if (handles.flederize)
				flederize(handles.nameList{i}, i, Z, imgWater, indLand, [handles.head_water(1:4) minWater maxWater])
			end

			if ( handles.useLandPhoto )
				if (i == 1)			% First round, open a new Mirone window
					hFig = mirone(imgWater, head_struct);
					handMir = guidata(hFig);
					alphaMask = alloc_mex(size(indLand),'uint8');	% Create an image mask of Dry/Wets
					h = image('XData',handles.geoPhotoX,'YData',handles.geoPhotoY, 'CData',handles.geoPhoto, 'Parent',handMir.axes1);
					uistack(h,'bottom')		% Send to bottom because we want the alphamask applyied to data derived image
					alfaMasca = round((1 - get(handles.slider_transparency, 'Val')) * 255);	% (1 - val) ... see above
				end
				alphaMask(~indLand) = alfaMasca;
 				set(handMir.hImg,'CData',imgWater,'AlphaData',alphaMask)
				F = getframe(handMir.axes1);		% Image capture. Unfortunately, without any control
				imgWater = F.cdata;
				if ( strcmp(get(handMir.axes1,'Ydir'),'normal') ),	imgWater = flipdim(imgWater,1);	end
			end

			[M, map] = aux_movie(handles, is_gif, is_mpg, is_avi, imgWater, i, M);
		end			% for (i=1:nGrids)
	else
		error('aquamoto:computing movie. Should not have passed here')
	end

	% This section is common to both Multi-layer grid or List of grids
	if (is_avi)
		mname = [handles.moviePato handles.movieName '.avi'];
		movie2avi_j(M,mname,'compression','none','fps',handles.fps)
	elseif (is_mpg)
		mname = [handles.moviePato handles.movieName '.mpg'];
		opt = [1, 0, 1, 0, 10, 5, 5, 5];
		mpgwrite(M,map,mname,opt)
	elseif (is_montage)
		montage(M)
	end

% -----------------------------------------------------------------------------------------
function [M, map] = aux_movie(handles, is_gif, is_mpg, is_avi, img, i, M, map)
% Auxiliary function to write movie frames
	if (nargin == 7),	map = [];	end
	if ( (ndims(img) == 3) && (is_gif || is_mpg) )
		[img,map] = img_fun('rgb2ind',img,256,handles.dither);
	end
	img = flipdim(img,1);			% The stupid UL origin

	if (is_gif)
		mname = [handles.moviePato handles.movieName '.gif'];
		if (i == 1)
			writegif(img,map,mname,'loopcount',Inf)
		else
			writegif(img,map,mname,'WriteMode','append','DelayTime',handles.dt)
		end
	elseif (is_avi)			% AVI
		M(i) = im2frame(img);
	else					% MPEG
		M(i) = im2frame(img,map);
	end

% -----------------------------------------------------------------------------------------
function slider_transparency_CB(hObject, eventdata, handles)
	val = get(hObject,'Value');
	set(handles.text_WaterTrans,'String',['Water transparency  ' sprintf('%.0f',val*100) '%'])

% -----------------------------------------------------------------------------------------
function [minWater, maxWater, heads] = push_checkGlobalMM_CB(hObject, eventdata, handles)
	if (numel(handles.nameList) == 0)		% netCDF files don't use this mechanism
		set(hObject,'Visible','off'),	return
	end
	[minWater, maxWater, heads] = get_globalMinMax(handles);
	set(handles.edit_globalWaterMin,'String',sprintf('%.2f',minWater))
	set(handles.edit_globalWaterMax,'String',sprintf('%.2f',maxWater))
	set(hObject,'Visible','off')
	handles.checkedMM = 1;          % Signal that ensemble Min/Max is now known
	handles.minWater = minWater;    % Make a copy that can be changed by a user
	handles.maxWater = maxWater;    % direct intervenction on the edit boxes
	handles.gridHeaders = heads;    % a cell array
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function [minWater, maxWater, heads] = get_globalMinMax(handles)
% Run trough all water level grids and find the ensemble Min/Max
	minWater = 1000;        maxWater = -1000;
	nGrids = numel(handles.nameList);
	heads = cell(nGrids,1);
	for (i=1:nGrids)
		[handles,heads{i}] = read_gmt_type_grids(handles,handles.nameList{i},'hdr');
		minWater = min(minWater,heads{i}(5));
		maxWater = max(maxWater,heads{i}(6));
	end

% -----------------------------------------------------------------------------------------
function edit_globalWaterMax_CB(hObject, eventdata, handles)
% User can change the ensemble limits, but we must keep track of that
	xx = str2double(get(hObject,'String'));
	if (~isnan(xx))
		handles.maxWater = xx;
		handles.minWater = str2double(get(handles.edit_globalWaterMin,'String'));
		handles.usrMM = 1;			% Flag that user has changed the ensemble Min|Max
		guidata(handles.figure1,handles)
	end

% -----------------------------------------------------------------------------------------
function edit_globalWaterMin_CB(hObject, eventdata, handles)
% User can change the ensemble limits, but we must keep track of that
	xx = str2double(get(hObject,'String'));
	if (~isnan(xx))
		handles.minWater = xx;
		handles.maxWater = str2double(get(handles.edit_globalWaterMax,'String'));
		handles.usrMM = 1;			% Flag that user has changed the ensemble Min|Max
		guidata(handles.figure1,handles)
	end

% -----------------------------------------------------------------------------------------
function edit_landPhoto_CB(hObject, eventdata, handles)
	fname = get(hObject,'String');
	if ( ~isempty(fname) )
		push_landPhoto_CB([], [], handles, fname)
	else					% Reset
		handles.geoPhoto = [];
		handles.useLandPhoto = false;
		guidata(handles.figure1,handles)
	end

% -----------------------------------------------------------------------------------------
function push_landPhoto_CB(hObject, eventdata, handles, opt)
	if ( isempty(handles.head) && isempty(handles.head_bat) )
		errordlg('You need to load the grid data first. Only than I''l know what to do with this image.','Error')
		return
	end

	str1 = {'*.tif;*.TIF;*.tiff;*.TIFF', 'GeoTiff or TIF *.tfw (*.tif,*.tiff,*.TIF,*.TIFF)'; ...
		'*.sid;*.SID', 'Mr Sid (*.sid,*.SID)'; ...
		'*.ecw;*.ECW', 'ECW (*.ecw,*.ECW)'; ...
		'*.jp2;*.JP2', 'Jpeg 2000 (*.jp2,*.JP2)'; ...
		'*.jpg;*.JPG;*.jpeg;*.JPEG', 'JPEG with *.jgw (*.jpg,*.jpeg)'; ...
		'*.png;*.PNG', 'PNG with *.pgw (*.png,*PNG)'; ...
		'*.gif;*.GIF', 'GIF with *.gfw (*.gif,*.GIF)'; ...
		'*.kap;*.KAP;*.nos;*NOS', 'BSB Nautical Chart (*.kap,*.KAP,*.nos,*.NOS)'; ...
		'*.*', 'All Files (*.*)'};
	if (nargin == 3)
		[FileName,PathName,handles] = put_or_get_file(handles,str1,'Select image format','get');
		if isequal(FileName,0),		return,		end
	else                % Filename was transmited in input
		PathName = [];      FileName = opt;
	end
	fileName = [PathName FileName];

	check_fw = true;       head = [];
	[PATH,FNAME,EXT] = fileparts(fileName);
	
	if ( any(strcmpi(EXT,{'.tif' '.tiff' '.jp2' '.ecw' '.sid' '.gif' '.kap'})) )		% GDAL territory
		try		att = gdalread(fileName,'-M','-C');
		catch,	errordlg(lasterr,'Error'),	return
		end
		if ( att.RasterCount == 0 || att.RasterCount > 3 || ~strcmp(att.Band(1).DataType,'Byte') )
			errordlg('Hmm, I don''t think this is an apropriate file to use here. Bye Bye','ERROR'),	return
		end
		opt_U = '-U';		if (isempty(att.GeoTransform)),		opt_U = ' ';	end
		I = gdalread(fileName, opt_U);
		if (~isempty(att.GeoTransform))     % Georeferenced image
			head = att.GMT_hdr;
			head(8) = diff(head(1:2)) / (size(I,2) - 1);      head(9) = diff(head(3:4)) / (size(I,1) - 1);
			check_fw = false;				% Don't need to search for an auxiliary *fw file
		else
			errordlg('I''m sorry but your file doesn''t seam to be georerefrenced. At least one that GDAL understands.','Error')
			return
		end
		if (strcmpi(att.ColorInterp,'palette') && ~isempty(att.Band(1).ColorMap) )
			try		pal = att.Band(1).ColorMap.CMap;     pal = pal(:,1:3);       % GDAL creates a Mx4 colormap
			catch,	errordlg('Figure ColorMap had troubles. Quiting.','ERROR'),		return
			end
			I = ind2rgb8(I, pal);
		elseif ( strcmpi(att.ColorInterp,'gray') )
			pal = repmat( (att.GMT_hdr(5):att.GMT_hdr(6))' / 255, 1, 3);
			I = ind2rgb8(I, pal);
		end
		if ( ndims(I) == 2 )
			errordlg('Unknow error in file. Something to do with an uncorrect color palette.','ERROR'),		return
		end

	else
		try		[I,pal] = imread(fileName);
		catch,	errordlg(lasterr,'Error'),	return	% It realy may happen
		end
		I = flipdim(I,1);
		if ( ndims(I) == 2 && ~isempty(pal) )
			I = ind2rgb8(I, pal);
		elseif ( ndims(I) == 2 && isempty(pal) )
			errordlg('Unknow error in file. The image is indexed but the color palette is empty.','ERROR'),		return
		end
	end
	
	if ( check_fw )
		[head,err_msg] = tfw_funs('inquire',[size(I,1) size(I,2)],PATH,FNAME,EXT);  % See if we have .*fw file
		if (~isempty(err_msg))
			warndlg(['A registering world file was found but the following error occured: ' err_msg],'Warning')
		end
	end
	
	if ( isempty(head) )
		errordlg('Most probably a BUG in calculating the file''s header. Please inform base.','ERROR'),		return
	end

	handles.headGeoImage = head;

	if ( ~isempty(handles.head) )		% SWW file present
		rect_crop = [handles.head(1) handles.head(3) diff(handles.head(1:2)) diff(handles.head(3:4))];
	else								% List of grids
		rect_crop = [handles.head_bat(1) handles.head_bat(3) diff(handles.head_bat(1:2)) diff(handles.head_bat(3:4))];
	end
	[I,r_c] = cropimg(head(1:2),head(3:4), I, rect_crop, 'out_grid');
	[m,n] = size(I);

	if (m < 11 || n < 11)
		warndlg('The area common to both the image and the modelling region is ridiculous small. Ignoring it','Warning')
		return
	end

	set(handles.edit_landPhoto, 'String', fileName)

	handles.geoPhoto = I;
	handles.geoPhotoX = [head(1) + (r_c(3)-1)*head(8) head(1) + (r_c(4)-1)*head(8)];
	handles.geoPhotoY = [head(3) + (r_c(1)-1)*head(9) head(3) + (r_c(2)-1)*head(9)];
	handles.useLandPhoto = true;

	if ( ishandle(handles.hMirFig) )
		push_showSlice_CB([], [], handles)	% and update image (also saves handles)
	end
	set([handles.textResize handles.popup_resize], 'Enable', 'off')
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function show_needed(handles,opt)
	h_all = handles.h_line;
	if (strncmp(opt,'grdgradient',11))
		set(handles.edit_elev,'Enable','off');			set(handles.edit_azim,'Visible','on')
		set(handles.text_elev,'Enable','on');			set(handles.edit_azim,'Enable','on');
		set(handles.text_azim,'Enable','on');
		set(h_all(1),'Visible','on');					set(h_all(2),'Visible','off')
		set(handles.toggle_1,'Value',1);				set(handles.toggle_2,'Value',0)
	else                            % Lambertian
		set(handles.edit_elev,'Enable','on');			set(handles.edit_azim,'Visible','on')
		set(handles.text_elev,'Enable','on');			set(handles.edit_azim,'Enable','on');
		set(handles.text_azim,'Enable','on');			set(h_all,'Visible','on');
		set(handles.toggle_1,'Value',0);				set(handles.toggle_2,'Value',1)
	end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function ButtonDown(obj,eventdata,h_all,handles)
	% It could be cleverer.
	pt = get(gca, 'CurrentPoint');
	x_lim = get(gca,'xlim');      y_lim = get(gca,'ylim');
	% check if x,y is inside of axis
	if ~((pt(1,1)>=x_lim(1)) && (pt(1,1)<=x_lim(2)) && (pt(1,2)>=y_lim(1)) && (pt(1,2)<=y_lim(2)))    % outside axis limits
		return
	end
	if any(h_all == gco)
		h = h_all(h_all == gco);    % When more than one line handle exists, find only the selected one
		set(handles.figure1,'WindowButtonMotionFcn',{@ButtonMotion,h,handles},'WindowButtonUpFcn',{@ButtonUp,h_all,handles},...
			'Pointer', 'crosshair');
	end

% -----------------------------------------------------------------------------------------
function ButtonMotion(obj,eventdata,h,handles)
	selectionType = get(handles.figure1, 'SelectionType');
	pt = get(gca, 'CurrentPoint');
	if strcmp(selectionType, 'normal')      % right-cick
		xx = get(h,'XData');		yy = get(h,'YData');
		theta = cart2pol(pt(1,1)-xx(1),pt(1,2)-yy(1));
		radius = get(h,'Userdata');
		x2 = xx(1) + radius * cos(theta);      y2 = yy(1) + radius * sin(theta);
		if strcmp(get(h,'Tag'),'Elev') && (theta >= 0 && theta <= pi/2)   % Elevation line
			set(h,'XData',[xx(1) x2],'YData',[yy(1) y2]);
			set(handles.edit_elev,'String',sprintf('%.0f',theta *180/pi) )
		elseif ~strcmp(get(h,'Tag'),'Elev')     % Azimuth line(s)
			set(h,'XData',[xx(1) x2],'YData',[yy(1) y2]);
			% NOTE to if I ever want to reuse this code. Normally ang_2pi should be = pi/2 - (pi*.....)
			% for the normal y origin at bottm left corner. However, due to the stupid habit of using y=0
			% at top left corner when dealing with images, to get an azimuth angle we have to do like following. 

			% truncate angles into [-pi pi] range
			ang_2pi = pi/2 + ( pi*((abs(theta)/pi) - 2*ceil(((abs(theta)/pi)-1)/2)) * sign(theta) );
			epsilon = -1e-7;        %  Allow points near zero to remain there
			indx = find(ang_2pi < epsilon);
			%  Shift the points in the [-pi 0] range to [pi 2pi] range
			if ~isempty(indx);  ang_2pi(indx) = ang_2pi(indx) + 2*pi;  end;
			set(handles.edit_azim,'String',sprintf('%.0f',ang_2pi *180/pi) )
		end
	end

% -----------------------------------------------------------------------------------------
function ButtonUp(obj,eventdata,h,handles)
	handles = guidata(handles.figure1);     % We need an updated version
	set(handles.figure1,'WindowButtonMotionFcn','','WindowButtonDownFcn', ...
		{@ButtonDown,h,handles},'WindowButtonUpFcn','');
	set(handles.figure1,'Pointer', 'arrow')
	azim = get(handles.edit_azim,'String');
	xdataAz = get(h(1),'XData');        ydataAz = get(h(1),'YData');
	xdataElev = get(h(2),'XData');      ydataElev = get(h(2),'YData');
	if (get(handles.radio_land,'Value'))    % We are setting the LAND illum params
		if (get(handles.toggle_1,'Value'))  % grdgradient classic
			handles.landIllumComm = ['-A' azim];
		else                                % Lambertian
			elev = get(handles.edit_elev,'String');
			handles.landIllumComm = ['-E' azim '/' elev handles.lambCteComm];
		end
		handles.landLineAzBack = [xdataAz ydataAz];				% Backup the graphical info
		handles.landLineElevBack = [xdataElev ydataElev];
		handles.landAzStrBack = get(handles.edit_azim,'String'); % And string-numeric
		handles.landElevStrBack = get(handles.edit_elev,'String');
	else                                    % WATER
		if (get(handles.toggle_1,'Value'))  % grdgradient classic
			handles.waterIllumComm = ['-A' azim];
		else                                % Lambertian
			elev = get(handles.edit_elev,'String');
			handles.waterIllumComm = ['-E' azim '/' elev handles.lambCteComm];
		end
		handles.waterLineAzBack = [xdataAz ydataAz];			% Backup the graphical info
		handles.waterLineElevBack = [xdataElev ydataElev];
		handles.waterAzStrBack = get(handles.edit_azim,'String'); % And string-numeric
		handles.waterElevStrBack = get(handles.edit_elev,'String');
	end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_batGrid_CB(hObject, eventdata, handles)
    fname = get(hObject,'String');
    push_batGrid_CB([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_batGrid_CB(hObject, eventdata, handles, opt)
    if (nargin == 3)        % Direct call
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*','All Files (*.*)'},'Select GMT grid','get');
		if isequal(FileName,0),		return,		end
	else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
	fname = [PathName FileName];
	
	[handles,handles.X_bat,handles.Y_bat,handles.Z_bat,handles.head_bat] = read_gmt_type_grids(handles,fname);
	if (isempty(handles.X_bat)),	return,		end

    handles.cmapBat = makeCmapBat(handles, handles.head_bat, handles.cmapLand, 1);	% Put the cmap discontinuity at the zero of bat

	handles.reinterpolated_bat = false;		% In case we need to reinterpolate bat to be compatible with water grid
	set(handles.edit_batGrid,'String',fname)
	[dump,FNAME] = fileparts(FileName);
	if ( get(handles.radio_gif,'Value') ),		EXT = '.gif';
	elseif ( get(handles.radio_mpg,'Value') )   EXT = '.mpg';
	else										EXT = '.avi';
	end
	handles.moviePato = PathName;
	handles.movieName = FNAME;
	set(handles.edit_movieName,'String',[PathName handles.movieName EXT])
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_namesList_CB(hObject, eventdata, handles)
    fname = get(hObject,'String');
    push_namesList_CB([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_namesList_CB(hObject, eventdata, handles, opt)
    if (nargin == 3)        % Direct call
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.DAT,*.txt,*.TXT)';'*.*', 'All Files (*.*)'},'File with grids list','get');
		if isequal(FileName,0),		return,		end		
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        if (~isempty(PathName)),	PathName = [PathName filesep];		end      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
	fname = [PathName FileName];

    [bin,n_column] = guess_file(fname);
    % If error in reading file
    if isempty(bin)
		errordlg(['Error reading file ' fname],'Error');    return
    end

	fid = fopen(fname);
	c = fread(fid,'*char')';      fclose(fid);
	names = strread(c,'%s','delimiter','\n');   clear c fid;
	m = length(names);
	
	handles.strTimes = [];          % To hold time steps as strings
	if (n_column > 1)
		handles.strTimes = cell(m,1);
		c = false(m,1);
		for (k=1:m)
			[t,r] = strtok(names{k});
			if (t(1) == '#'),  c(k) = true;  continue;   end
			names{k} = t;
			handles.strTimes{k} = r;
		end
		% Remove eventual commented lines
		if (any(c))
			names(c) = [];			handles.strTimes(c) = [];
			m = length(names);      % Count remaining ones
		end
	end
    
	handles.shortNameList = cell(m,1);      % To hold grid names with path striped
	c = false(m,1);
	for (k=1:m)
		if ( isempty(names{k}) ),	continue,		end		% Jump empty lines
		if (n_column == 1 && names{k}(1) == '#')    % If n_column > 1, this test was already done above
			c(k) = true;    continue;
		end
		[PATH,FNAME,EXT] = fileparts(names{k});
		if (isempty(PATH))
			handles.shortNameList{k} = names{k};
			names{k} = [PathName names{k}];
		else
			handles.shortNameList{k} = [FNAME EXT];
		end
		if (any(c))
			names(c) = [];		handles.shortNameList(c) = [];
		end
	end
    
	% Check that at least the files in provided list do exist
	c = false(m,1);
	for (k=1:m)
		c(k) = (exist(names{k},'file') ~= 2);
	end
	names(c) = [];      handles.shortNameList(c) = [];

	handles.nameList = names;
	handles.checkedMM = 0;
	set(handles.listbox1,'String',handles.shortNameList)
	set(handles.edit_namesList,'String',[PathName FileName])
	set(handles.push_checkGlobalMM,'Visible','on')
	set(handles.radio_timeGridsList,'Val',1)
	set(handles.radio_multiLayer,'Val',0)
	set(handles.edit_multiLayerInc, 'Enable', 'off')
	set([handles.textResize handles.popup_resize], 'Enable', 'on')
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function listbox1_CB(hObject, eventdata, handles)
    % if this is a doubleclick,
    if ( strcmp(get(gcbf,'SelectionType'),'open') && ~isempty(handles.nameList) )
        val = get(hObject,'Value');
        if (~isempty(handles.strTimes))
            handles.testTime = handles.strTimes{val};
        end
        push_singleWater_CB([], [], handles, handles.nameList{val})
    end

% -----------------------------------------------------------------------------------------
function push_clearTestBat_CB(hObject, eventdata, handles)
	set(handles.edit_singleWater,'String','')
	set(hObject,'Visible','off')
	handles.Z_water = [];
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_singleWater_CB(hObject, eventdata, handles)
	fname = get(hObject,'String');
	push_singleWater_CB([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_singleWater_CB(hObject, eventdata, handles, opt)
    if (nargin == 3)        % Direct call
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'},'Select GMT grid','get');
		if isequal(FileName,0),		return,		end		
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
    fname = [PathName FileName];

	[handles,X,Y,handles.Z_water,handles.head_water] = read_gmt_type_grids(handles,fname);
	if (isempty(X)),    return;     end
	set(handles.push_clearTestBat,'Visible','on')
	
	set(handles.edit_singleWater,'String',fname)
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_maregs_CB(hObject, eventdata, handles)
    fname = get(hObject,'String');
	if ( ~isempty(fname) ),    push_maregs_CB([], [], handles, fname),	end

% -----------------------------------------------------------------------------------------
function push_maregs_CB(hObject, eventdata, handles, opt)

    if (nargin == 3)        % Direct call
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.txt)';'*.*', 'All Files (*.*)'},'(x,y) file','get');
		if isequal(FileName,0),		return,		end
	else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        if (~isempty(PathName)),	PathName = [PathName filesep];	end		% To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
	fname = [PathName FileName];
	
	if (exist(fname, 'file') ~= 2)
		errordlg(['File: ' fname ' does not exist.'],'Error'),		return
	end
    set(handles.edit_maregs,'String',fname)
	
	% Set a default name for interpolated file
 	[PATH,FNAME] = fileparts(fname);
	set(handles.edit_profile,'String',[PATH filesep FNAME '_intp.dat'])
   
	[bin,n_column,multi_seg,n_headers] = guess_file(fname);
	% If error in reading file
	if ( isempty(bin) ),	errordlg(['Error reading file ' fname],'Error'),	return,		end
	if (isa(bin,'struct') || bin ~= 0)
        errordlg('Sorry, reading binary files is not programed','Error'),	return
	end
	if (n_column < 2)
        errordlg('This file doesn''t have at least 2 columns','Error'),	return
	end
	if (isempty(n_headers)),    n_headers = NaN;    end
	if (multi_seg),		numeric_data = text_read(fname,NaN,n_headers,'>');
	else				numeric_data = text_read(fname,NaN,n_headers);
	end

	% NEED TO FORESEE THE CASE OF MULTISEGS
	handles.xyData = numeric_data;
	set([handles.push_interpolate handles.edit_miscLayerInc], 'Enable', 'on') 
	guidata(handles.figure1, handles)
	
% -----------------------------------------------------------------------------------------
function push_interpolate_CB(hObject, eventdata, handles)
	if (isempty(handles.fname))
		errordlg('Yea I do. But what? Things work better if you load a file first','Error'),	return
	end

	layerInc = edit_miscLayerInc_CB(handles.edit_miscLayerInc, eventdata, handles);
	n_layers = numel(layerInc);
	ncols = size(handles.xyData,1);		% Each point will be a column in output file, so number of cols is = npoints

	Z = single( zeros(n_layers, ncols) );
	if (handles.is_coards)
		z_id = handles.netcdf_z_id;
		s = handles.nc_info;			% Retrieve the .nc info struct
		head = handles.head;
	end

	aguentabar(0,'title','Interpolating multi-layer file','CreateCancelBtn')
	theVarName = 'none';		% Default
	for (k = layerInc)
		handles.sliceNumber = k - 1;
		if (handles.is_sww)
			[theVar, U, V, indVar, indWater, theVarName] = get_swwVar(handles);		% Does not save handles
			if (~isa(theVar, 'double')),	theVar = double(theVar);	end
			z = mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) theVar(:)],handles.xyData(:,1),handles.xyData(:,2), 'p');
		elseif (handles.is_coards)
			z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [handles.sliceNumber 0 0], [1 s.Dataset(z_id).Size(end-1:end)]);
			z = grdtrack_m(z,head,handles.xyData(:,1:2),'-Z')';
		end
		Z(k,:) = z';
		h = aguentabar(k/n_layers);
		if (isnan(h)),	break,	end
	end
	if (isnan(h)),	return,		end

	fname = get(handles.edit_profile,'String');
	if (isempty(fname))			% If it's empty ask for it
		[FileName,PathName] = put_or_get_file(handles,{'*.dat;*.DAT','ASCII file'; '*.*', 'All Files (*.*)'},'Output file','put');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	end
	
	%Open and write to ASCII file
	if (ispc),		fid = fopen(fname,'wt');
	elseif (isunix),fid = fopen(fname,'w');
	else			error('aquamoto: Unknown platform.');
	end
 	if ( get(handles.check_miscWriteHeader, 'Val') )		% Write an header
		fprintf(fid,'%s\n', ['# Variable: ' theVarName ]);
		fprintf(fid, ['#  \t', repmat('%g(X)\t', [1,ncols]) '\n'], handles.xyData(:,1));
		fprintf(fid, ['# T\t', repmat('%g(Y)\t', [1,ncols]) '\n'], handles.xyData(:,2));
	end
	t = handles.time(layerInc);		% Layers's times
	fprintf(fid,['%.2f\t' repmat('%f\t',[1,ncols]) '\n'], [t(:) double(Z)]');
	fclose(fid);

% -----------------------------------------------------------------------------------------
function save_arrows(hObject, evt, hQuiver)
% Save the arrow field on file
	handles = guidata(hObject);
	str1 = {'*.dat;*.DAT', 'Symbol file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select Symbol File name','put','.dat');
	if isequal(FileName,0),		return,		end
	f_name = [PathName FileName];
	
	setas = get(hQuiver, 'UserData');
	x = repmat(setas{1},numel(setas{2}),1);
	y = repmat(setas{2}(:),numel(setas{1}),1);
	double2ascii(f_name, [x(:) y(:) setas{3}(:) setas{4}(:)], '%.6f\t%.6f\t%.6f\t%.6f');

% -----------------------------------------------------------------------------------------
function push_plugFun_CB(hObject, eventdata, handles)
% THIS IS A SPECIAL CALLBACK THAT CALLS A FUNCTION NAMED 'aquaPlugin' THAT MAY RESIDE
% ANYWHERE IN THE PATH WORLD. IT'S UP TO THE USER TO DIFFINE ITS CONTENTS.
	aquaPlugin(handles)		% That's all it should be needed

% -----------------------------------------------------------------------------------------
function layerInc = edit_miscLayerInc_CB(hObject, eventdata, handles)
% This function's code reuses the code of edit_multiLayerInc_CB() because they
% are supposed to do exactly the same thing, which is - creating a vector map of layers to use
	bak = handles.multiLayerInc;		% Backup copy of this field
	edit_multiLayerInc_CB(hObject, eventdata, handles)
	handles = guidata(handles.figure1);	% Get updated handles
	layerInc = handles.multiLayerInc;
	handles.multiLayerInc = bak;		% Restore original value
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function [theVar, U, V, indVar, indWater, theVarName] = get_swwVar(handles)
% Get a primary or derived variable from a sww file. The slice is controled by handles.sliceNumber
% Since there is no saving of handles, the calling function can safely apply a teporary change to that handles member

	if ( ~get(handles.check_derivedVar,'Val') )			% Get one of the primary quantities
		% Get the ploting variable
		if (get(handles.radio_stage, 'Val')),			theVarName = 'stage';		indVar = 1;
		elseif (get(handles.radio_xmoment, 'Val')),		theVarName = 'xmomentum';	indVar = 2;
		elseif (get(handles.radio_ymoment, 'Val')),		theVarName = 'ymomentum';	indVar = 3;
		end
		theVar = nc_funs('varget', handles.fname, theVarName, [handles.sliceNumber 0], [1 handles.number_of_points]);
		U = [];		V = [];
		indWater = [];					% Not empty when a Max quantity (obvioulsly, not the case here)
	else
		[theVar, U, V, indVar, indWater, theVarName] = get_derivedVar(handles);
	end

% -----------------------------------------------------------------------------------------
function geog = guessGeog(lims)
   % Make a good guess if LIMS are geographic
    geog = double( ( (lims(1) >= -180 && lims(2) <= 180) || (lims(1) >= 0 && lims(2) <= 360) )...
        && (lims(3) >= -90 && lims(4) <= 90) );
	
% -----------------------------------------------------------------------------------------
function flederize(fname,n, Z, imgWater, indLand, limits)
% Write a .sd fleder file with z_max smashed to (?) times min water height
	pato = fileparts(fname);
	%fname = [pato filesep name '.sd'];
	fname = [pato filesep sprintf('z_%.2d.sd',n)];
	
	s = 4;				% Smash land to ? times max water height
	maxWater = 10;
	minWater = -17;		% We could use limits(5), but it's not sure min is not on land
	
	fact = abs( (s * maxWater) / limits(6) );		% smashing factor
	Z_smashed = Z;
	Z_smashed(indLand) = single(double(Z(indLand)) * fact);
	write_flederFiles('main_SD', fname, 'Planar', Z_smashed, imgWater, [limits(1:4) minWater maxWater*s])

% --------------------------------------------------------------------
function hh = loc_quiver(struc_in,varargin)
%QUIVER Quiver plot.
%   QUIVER(X,Y,U,V) plots velocity vectors as arrows with components (u,v)
%   at the points (x,y).  The matrices X,Y,U,V must all be the same size
%   and contain corresponding position and velocity components (X and Y
%   can also be vectors to specify a uniform grid).  QUIVER automatically
%   scales the arrows to fit within the grid.
%
%   QUIVER(X,Y,U,V,S) automatically scales the arrows to fit within the grid and
%   then stretches them by S. Use  S=0 to plot the arrows without the automatic scaling.
%
%   H = QUIVER(...) returns a vector of line handles.

	% Arrow head parameters
	alpha = 0.33;		% Size of arrow head relative to the length of the vector
	beta = 0.33;		% Width of the base of the arrow head relative to the length
	autoscale = 1;		% Autoscale if ~= 0 then scale by this.
	subsample = 1;		% Plot one every other grid node vector

	nin = nargin - 1;

	% Check numeric input arguments
	if (nin < 4)					% quiver(u,v) or quiver(u,v,s)
		[msg,x,y,u,v] = xyzchk(varargin{1:2});
	else
		[msg,x,y,u,v] = xyzchk(varargin{1:4});
	end
	if ~isempty(msg), error(msg); end

	if (nin == 5)		% quiver(x,y,u,v,s)
		autoscale = varargin{nin};
	elseif  (nin == 6)
		autoscale = varargin{nin-1};
		subsample = abs(round(varargin{nin}));
	end

	% Scalar expand u,v
	if (numel(u) == 1),     u = u(ones(size(x))); end
	if (numel(v) == 1),     v = v(ones(size(u))); end

	if (subsample > 1)
		x = x(1:subsample:end,1:subsample:end);		y = y(1:subsample:end,1:subsample:end);
		u = u(1:subsample:end,1:subsample:end);		v = v(1:subsample:end,1:subsample:end);
	end

	if (autoscale)
		% Base autoscale value on average spacing in the x and y directions.
		% Estimate number of points in each direction as either the size of the
		% input arrays or the effective square spacing if x and y are vectors.
		if min(size(x))==1, n=sqrt(numel(x)); m=n; else [m,n]=size(x); end
		delx = diff([min(x(:)) max(x(:))])/n;
		dely = diff([min(y(:)) max(y(:))])/m;
		del = delx.^2 + dely.^2;
		if (del > 0)
			len = (u.^2 + v.^2)/del;
			maxlen = sqrt(max(len(:)));
		else
			maxlen = 0;
		end
		
		if maxlen > 0
			autoscale = autoscale*0.9 / maxlen;
		else
			autoscale = autoscale*0.9;
		end
		u = u*autoscale; v = v*autoscale;
	end

	% Make velocity vectors
	x = x(:).';		y = y(:).';
	u = u(:).';		v = v(:).';
	uu = [x; x+u; ones(size(u))*NaN];
	vv = [y; y+v; ones(size(u))*NaN];
	% Make arrow heads
	hu = [x+u-alpha*(u+beta*(v+eps)); x+u; x+u-alpha*(u-beta*(v+eps)); ones(size(u))*NaN];
	hv = [y+v-alpha*(v-beta*(u+eps)); y+v; y+v-alpha*(v+beta*(u+eps)); ones(size(v))*NaN];

	if (struc_in.spacingChanged)
		try		delete(struc_in.hQuiver),	struc_in.hQuiver = [];	end		% Remove previous arrow field
	end

	if ( isempty(struc_in.hQuiver) || ~ishandle(struc_in.hQuiver(1)) )		% No arrows yet.
		h1 = line('XData',uu(:), 'YData',vv(:), 'Parent',struc_in.hAx, 'Color','k');
		h2 = line('XData',hu(:), 'YData',hv(:), 'Parent',struc_in.hAx, 'Color','k');
		if (nargout > 0),	hh = [h1;h2];	end
	else
		% We have the arrows and only want to change them
		set(struc_in.hQuiver(1),'XData',uu(:), 'YData',vv(:))
		set(struc_in.hQuiver(2),'XData',hu(:), 'YData',hv(:))
		if (nargout > 0),	hh = [struc_in.hQuiver(1); struc_in.hQuiver(2)];	end
	end

% ------------------------------------------------------------------------------
function changecolor(obj,evt,hPatch)
% Change the color of the progress bar patch
	colorlim = 2.8;				% Must be <= 3.0 - This keeps the color from being too light
	thiscolor = rand(1,3);
	while (sum(thiscolor) > colorlim)
		thiscolor = rand(1,3);
	end
	set(hPatch,'FaceColor',thiscolor);

% ------------------------------------------------------------------------------
function timestr = sec2timestr(sec)
% Convert a time measurement from seconds into a human readable string.

	% Convert seconds to other units
	d = floor(sec/86400);		sec = sec - d*86400;		% Days and remaing seconds
	h = floor(sec/3600);		sec = sec - h*3600;			% Hours and remaing seconds
	m = floor(sec/60);			sec = sec - m*60;			% Minutes and remaing seconds
	s = floor(sec);				% Seconds

	% Create time string
	if (d > 0)
		timestr = sprintf('%d day, %.1f hr', d, (h+m/60));
	elseif (h > 0)
		timestr = sprintf('%d hr, %d min',h, m);
	elseif (m > 0)
		if (m > 9),		timestr = sprintf('%d min',m);
		else			timestr = sprintf('%d min, %d sec',m,s);
		end
	else
		timestr = sprintf('%d sec',s);
	end

% ----- Insert those other guys here --------------

%--------------------------------------------------

% --- Creates and returns a handle to the GUI figure. 
function aquamoto_LayoutFcn(h1)

set(h1,'Position',[520 347 1030 601],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'MenuBar','none',...
'Name','Aquamoto',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'tab_group_CB'},...
'Enable','inactive',...
'Position',[132 399 50 21],...
'String','Misc',...
'ButtonDownFcn',{@aquamoto_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','misc');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'tab_group_CB'},...
'Enable','inactive',...
'Position',[10 399 50 21],...
'String','ANUGA',...
'ButtonDownFcn',{@aquamoto_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'tab_group_CB'},...
'Enable','inactive',...
'Position',[183 399 100 21],...
'String','Shading OR Image',...
'ButtonDownFcn',{@aquamoto_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','shade');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'tab_group_CB'},...
'Enable','inactive',...
'Position',[284 399 50 21],...
'String','Cinema',...
'ButtonDownFcn',{@aquamoto_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','cinema');

uicontrol('Parent',h1, 'Position',[335 399 46 21],...
'Call',{@aquamoto_uiCB,h1,'tab_group_CB'},...
'Enable','inactive',...
'String','Plug',...
'ButtonDownFcn',{@aquamoto_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','plug');

uicontrol('Parent',h1, 'Position',[61 399 70 21],...
'Call',{@aquamoto_uiCB,h1,'tab_group_CB'},...
'Enable','inactive',...
'String','Time Grids',...
'ButtonDownFcn',{@aquamoto_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','t_grids');

uicontrol('Parent',h1, 'Position',[10 10 371 391],...
'Enable','inactive',...
'Tag','push_bg');

uicontrol('Parent',h1, 'Position',[420 173 181 41],'Style','frame','Tag','frame5','UserData','cinema');
uicontrol('Parent',h1, 'Position',[690 190 221 41],'Style','frame','Tag','frame6','UserData','t_grids');
uicontrol('Parent',h1, 'Position',[420 238 181 58],'Style','frame','Tag','frame4','UserData','cinema');
uicontrol('Parent',h1, 'Position',[610 238 71 58], 'Style','frame','Tag','frame8','UserData','cinema');

uicontrol('Parent',h1, 'Position',[690 455 311 21],...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_batGrid_CB'},...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','Name of bathymetry grid',...
'Tag','edit_batGrid',...
'UserData','t_grids');

uicontrol('Parent',h1, 'Position',[1001 455 21 21],...
'Call',{@aquamoto_uiCB,h1,'push_batGrid_CB'},...
'TooltipString','Browse for a bathymetry grid',...
'Tag','push_batGrid',...
'UserData','t_grids');

uicontrol('Parent',h1, 'Position',[699 195 180 21],...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_singleWater_CB'},...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','Name of one water level grid for test illumination purposes',...
'Tag','edit_singleWater',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'push_singleWater_CB'},...
'Position',[879 195 21 21],...
'Tag','push_singleWater',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_namesList_CB'},...
'HorizontalAlignment','left',...
'Position',[690 405 311 21],...
'Style','edit',...
'TooltipString','Name of a file with the water level grids list',...
'Tag','edit_namesList',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'push_namesList_CB'},...
'Position',[1001 405 21 21],...
'TooltipString','Browse for a water level grids list file',...
'Tag','push_namesList',...
'UserData','t_grids');

uicontrol('Parent',h1, 'Position',[690 246 221 158],...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'listbox1_CB'},...
'Style','listbox',...
'Value',1,...
'Tag','listbox1',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'toggle_1_CB'},...
'Position',[470 537 20 20],...
'Enable', 'off',...
'Style','togglebutton',...
'TooltipString','GMT grdgradient classic Illumination',...
'Tag','toggle_1',...
'UserData','shade');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'toggle_2_CB'},...
'Position',[491 537 20 20],...
'Enable', 'off',...
'Style','togglebutton',...
'TooltipString','Lambertian with lighting Illumination',...
'Tag','toggle_2',...
'UserData','shade');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[616 419 30 21],...
'String','30',...
'Style','edit',...
'TooltipString','Elevation light direction',...
'Tag','edit_elev',...
'UserData','shade');

uicontrol('Parent',h1, 'Position',[620 255 51 21],...
'Call',{@aquamoto_uiCB2,h1,[],'push_OK_CB'},...
'FontName','Helvetica',...
'FontSize',10,...
'String','Do it',...
'TooltipString','Display multiple images as a montage of subplots',...
'Tag','push_montage',...
'UserData','cinema');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'push_OK_CB'},...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[596 10 66 21],...
'String','OK',...
'Tag','push_OK',...
'UserData','cinema');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[599 500 50 14],...
'String','Elev',...
'Style','text',...
'Tag','text_elev',...
'UserData','shade');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[513 419 34 21],...
'String','0',...
'Style','edit',...
'TooltipString','Azimuth direction',...
'Tag','edit_azim',...
'UserData','shade');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[469 420 42 16],...
'String','Azimuth',...
'Style','text',...
'Tag','text_azim',...
'UserData','shade');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'radio_avi_CB'},...
'FontName','Helvetica',...
'Position',[431 246 50 15],...
'String','AVI',...
'Style','radiobutton',...
'TooltipString','Write movie file in RGB AVI format',...
'Tag','radio_avi',...
'UserData','cinema');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'radio_land_CB'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[519 537 50 17],...
'String','Land',...
'Style','radiobutton',...
'TooltipString','Set parametrs for Land illumination',...
'Value',1,...
'Tag','radio_land',...
'UserData','shade');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'radio_water_CB'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[618 535 55 17],...
'String','Water',...
'Style','radiobutton',...
'TooltipString','Set parametrs for Water illumination',...
'Tag','radio_water',...
'UserData','shade');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'radio_gif_CB'},...
'FontName','Helvetica',...
'Position',[431 266 50 15],...
'String','GIF',...
'Style','radiobutton',...
'TooltipString','Write movie file in animated GIF format',...
'Value',1,...
'Tag','radio_gif',...
'UserData','cinema');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'check_dither_CB'},...
'FontName','Helvetica',...
'Position',[550 269 55 15],...
'String','Dither',...
'Style','checkbox',...
'TooltipString','If you don''t know what this is, ask google',...
'Tag','check_dither',...
'UserData','cinema');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[432 288 70 15],...
'String','Movie type',...
'Style','text',...
'Tag','text_MovType',...
'UserData','cinema');

uicontrol('Parent',h1, 'Position',[617 288 55 15],...
'FontName','Helvetica', 'FontSize',9,...
'String','Montage',...
'Style','text',...
'Tag','text_montage',...
'UserData','cinema');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_fps_CB'},...
'Position',[550 246 30 18],...
'String','5',...
'Style','edit',...
'TooltipString','Frames per second',...
'Tag','edit_fps',...
'UserData','cinema');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Position',[493 248 55 15],...
'String','Frames p/s',...
'Style','text',...
'Tag','text15',...
'UserData','cinema');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[690 477 82 17],...
'String','Bathymetry file',...
'Style','text',...
'Tag','text16',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[690 426 82 17],...
'String','Water files list',...
'Style','text',...
'Tag','text17',...
'UserData','t_grids');

axes('Parent',h1,...
'Units','pixels',...
'Position',[470 440 91 91],...
'XTick',[],...
'YTick',[],...
'Tag','axes1',...
'UserData','shade',...
'Visible','off');

axes('Parent',h1,...
'Units','pixels',...
'Position',[600 440 61 60],...
'Tag','axes2',...
'UserData','shade',...
'Visible','off');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[699 221 100 17],...
'String','Test with this file',...
'Style','text',...
'Tag','text_testFile',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'push_clearTestBat_CB'},...
'FontName','Helvetica',...
'FontSize',7,...
'Position',[860 214 40 16],...
'String','Clear',...
'TooltipString','Remove the test grid so that you can compute a movie',...
'Tag','push_clearTestBat',...
'UserData','t_grids',...
'Visible','off');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'popup_resize_CB'},...
'Enable','off',...
'Position',[607 175 50 22],...
'String',{  '0.25'; '0.33'; '0.5'; '1.0'; '2.0' },...
'Style','popupmenu',...
'TooltipString','Resize output images by this value',...
'Value',4,...
'Tag','popup_resize',...
'UserData','cinema');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_movieName_CB'},...
'HorizontalAlignment','left',...
'Position',[350 45 311 21],...
'Style','edit',...
'TooltipString','Name of movie file',...
'Tag','edit_movieName',...
'UserData','cinema');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'push_movieName_CB'},...
'Position',[661 44 21 21],...
'TooltipString','Browse for a movie file name (extention is ignored)',...
'Tag','push_movieName',...
'UserData','cinema');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[355 66 115 17],...
'String','Output movie name',...
'Style','text',...
'Tag','text23',...
'UserData','cinema');

uicontrol('Parent',h1,...
'Enable','off',...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[609 197 50 17],...
'String','Scale it?',...
'Style','text',...
'Tag','textResize',...
'UserData','cinema');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[920 370 101 21],...
'String',{  'Surface height'; 'Water depth' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_surfType',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'radio_mpg_CB'},...
'FontName','Helvetica',...
'Position',[492 266 58 16],...
'String','MPEG',...
'Style','radiobutton',...
'TooltipString','Write movie file in mpeg format',...
'Tag','radio_mpg',...
'UserData','cinema');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'push_palette_CB'},...
'Position',[599 514 18 18],...
'TooltipString','Choose another Land palette',...
'Tag','push_palette',...
'UserData','shade');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'check_resetCmaps_CB'},...
'FontName','Helvetica',...
'Position',[620 516 60 15],...
'Enable', 'off', ...
'String','Reset',...
'Style','checkbox',...
'TooltipString','Reset to default color palettes',...
'Tag','check_resetCmaps',...
'UserData','shade');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_landPhoto_CB'},...
'HorizontalAlignment','left',...
'Position',[370 330 311 21],...
'Style','edit',...
'TooltipString','Terrain image name (of a georeferenced image file)',...
'Tag','edit_landPhoto',...
'UserData','shade');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'push_landPhoto_CB'},...
'Position',[680 329 21 21],...
'TooltipString','Browse for a Terrain image name (of a georeferenced image)',...
'Tag','push_landPhoto',...
'UserData','shade');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'Position',[445 355 170 34],...
'String',{  'OR (fantastic)'; 'Satelitte/Aerial Photograph' },...
'Style','text',...
'Tag','text28',...
'UserData','shade');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[463 180 40 20],...
'String','300',...
'Style','edit',...
'TooltipString','Image height in pixels',...
'Tag','edit_imgHeight',...
'UserData','cinema');

uicontrol('Parent',h1, 'Position',[551 182 40 20],...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'String','600',...
'Style','edit',...
'TooltipString','Image width in pixels',...
'Tag','edit_imgWidth',...
'UserData','cinema');

uicontrol('Parent',h1, 'Position',[460 204 100 17],...
'Enable','off',...
'FontName','Helvetica',...
'FontSize',9,...
'String','Movie frame size',...
'Style','text',...
'Tag','text_MovSize',...
'UserData','cinema');

uicontrol('Parent',h1, 'Position',[431 181 31 15],...
'Enable','off',...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Height',...
'Style','text',...
'Tag','text38',...
'UserData','cinema');

uicontrol('Parent',h1, 'Position',[518 184 30 15],...
'Enable','off',...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Width',...
'Style','text',...
'Tag','text39',...
'UserData','cinema');

uicontrol('Parent',h1, 'Position',[370 401 331 2], 'ForegroundColor',[0 0.501961 0],...
'Style','frame', 'Tag','frame_h1', 'UserData','shade');

uicontrol('Parent',h1, 'Position',[30 271 331 38], 'Style','frame', 'Tag','frame2', 'UserData','anuga');
uicontrol('Parent',h1, 'Position',[30 309 331 35], 'Style','frame', 'Tag','frame1', 'UserData','anuga');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[319 253 40 15],...
'String','Slice n?',...
'Style','text',...
'Tag','text1',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[0.99 0.99 0.99],...
'Call',{@aquamoto_uiCB,h1,'slider_layer_CB'},...
'Enable','inactive',...
'Position',[30 231 291 17],...
'Style','slider',...
'Tag','slider_layer',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_swwName_CB'},...
'HorizontalAlignment','left',...
'Position',[30 355 310 21],...
'Style','edit',...
'TooltipString','Enter a ANUGA .sww file name',...
'Tag','edit_swwName',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'push_swwName_CB'},...
'Position',[340 356 21 21],...
'TooltipString','Browse for a sww file name',...
'Tag','push_swwName',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',get(0,'factoryUicontrolBackgroundColor'),...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[27 377 135 17],...
'String','Input a SWW or NC file',...
'Style','text',...
'Tag','text2',...
'UserData','anuga');

uicontrol('Parent',h1, 'Position',[30 108 331 97], 'Style','frame', 'Tag','frame3', 'UserData','anuga');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[50 198 121 15],...
'String','Griding Line Geometry',...
'Style','text',...
'Tag','text_GLG',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_x_max_CB'},...
'HorizontalAlignment','left',...
'Position',[146 166 75 21],...
'Style','edit',...
'Tag','edit_x_max',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_y_max_CB'},...
'HorizontalAlignment','left',...
'Position',[146 140 75 21],...
'Style','edit',...
'Tag','edit_y_max',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[38 170 30 15],...
'String','X Dir',...
'Style','text',...
'Tag','text4',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[37 144 30 15],...
'String','Y Dir',...
'Style','text',...
'Tag','text5',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[166 187 41 13],...
'String','Max',...
'Style','text',...
'Tag','text6',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[87 187 41 13],...
'String','Min',...
'Style','text',...
'Tag','text7',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_x_inc_CB'},...
'HorizontalAlignment','left',...
'Position',[226 166 71 21],...
'Style','edit',...
'TooltipString','DX grid spacing',...
'Tag','edit_x_inc',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_y_inc_CB'},...
'HorizontalAlignment','left',...
'Position',[226 140 71 21],...
'Style','edit',...
'TooltipString','DY grid spacing',...
'Tag','edit_y_inc',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_Ncols_CB'},...
'HorizontalAlignment','left',...
'Position',[302 166 50 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Ncols',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_Nrows_CB'},...
'HorizontalAlignment','left',...
'Position',[302 140 50 21],...
'Style','edit',...
'TooltipString','Number of rows in the grid',...
'Tag','edit_Nrows',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[243 189 41 13],...
'String','Spacing',...
'Style','text',...
'Tag','text8',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[305 189 51 13],...
'String','# of lines',...
'Style','text',...
'Tag','text9',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Call',{@aquamoto_uiCB,h1,'push_Help_CB'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[295 115 61 18],...
'String','?',...
'Tag','push_Help',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_x_min_CB'},...
'HorizontalAlignment','left',...
'Position',[65 166 75 21],...
'Style','edit',...
'Tag','edit_x_min',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_y_min_CB'},...
'HorizontalAlignment','left',...
'Position',[65 141 75 21],...
'Style','edit',...
'Tag','edit_y_min',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'push_showMesh_CB'},...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','bold',...
'Position',[146 72 100 21],...
'String','Show mesh',...
'TooltipString','Show the mesh triangulation used in simulation',...
'Tag','push_showMesh',...
'UserData','anuga');

uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[175 378 200 16],...
'String','Info',...
'Style','text',...
'Tag','text_Info',...
'UserData','anuga');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_sliceNumber_CB'},...
'Enable','inactive',...
'Position',[320 231 40 20],...
'String','1',...
'Style','edit',...
'TooltipString','Slice number (to go to a specific one, enter a valid slice number here)',...
'Tag','edit_sliceNumber',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'radio_stage_CB'},...
'Position',[37 320 60 15],...
'String','Stage',...
'Style','radiobutton',...
'TooltipString','Plot the "stage" variable',...
'Value',1,...
'Tag','radio_stage',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'radio_xmoment_CB'},...
'Position',[136 320 95 15],...
'String','Xmomentum',...
'Style','radiobutton',...
'TooltipString','Plot the "xmomentum" variable',...
'Tag','radio_xmoment',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'radio_ymoment_CB'},...
'Position',[263 320 95 15],...
'String','Ymomentum',...
'Style','radiobutton',...
'TooltipString','Plot the "ymomentum" variable',...
'Tag','radio_ymoment',...
'UserData','anuga');

uicontrol('Parent',h1, 'Position',[131 279 190 22],...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'popup_derivedVar_CB'},...
'Enable','off',...
'String',{  'Absolute Velocity (V)'; 'Absolute Momentum (VxD)'; 'Water Depth'; 'Topography'; 'Froude Number'; 'Velocity Head (V^2 / (2g))'; 'Shear Stress'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_derivedVar',...
'UserData','anuga');

uicontrol('Parent',h1, 'Position',[38 283 92 15],...
'Call',{@aquamoto_uiCB,h1,'check_derivedVar_CB'},...
'String','Derived var',...
'Style','checkbox',...
'TooltipString','Select a derived quantity from the side popup menu',...
'Tag','check_derivedVar',...
'UserData','anuga');

uicontrol('Parent',h1, 'Position',[60 336 100 15],...
'Enable','inactive',...
'String','Primary quantities',...
'Style','text',...
'Tag','text_Pq',...
'UserData','anuga');

uicontrol('Parent',h1, 'Position',[328 279 25 23],...
'Call',{@aquamoto_uiCB,h1,'push_vel_CB'},...
'Enable','off',...
'Style','pushbutton',...
'TooltipString','Open a new window with plot arrow field controls',...
'Tag','push_vel',...
'UserData','anuga');

uicontrol('Parent',h1,'Position',[212 50 110 15],...
'String','Split Dry/wet',...
'Style','checkbox',...
'TooltipString','If checked, water and land parts of the image are built separately - Nice with shadings.',...
'Value',1,...
'Tag','check_splitDryWet');

uicontrol('Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Call',{@aquamoto_uiCB,h1,'slider_transparency_CB'},...
'Position',[20 28 180 16],...
'Style','slider',...
'Value',0.0,...
'Tag','slider_transparency');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[20 45 150 17],...
'String','Water transparency  0%',...
'Style','text',...
'Tag','text_WaterTrans');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'push_showSlice_CB'},...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','bold',...
'Position',[260 72 100 21],...
'String','Show slice',...
'TooltipString','Extract the slice selected in "Slice n?" and shot it in a Mirone window',...
'Tag','push_showSlice',...
'UserData','anuga');

uicontrol('Parent',h1,...
'Position',[212 29 165 15],...
'String','Scale color to global min/max',...
'Style','checkbox',...
'TooltipString','If checked, water color (if other than plain blue) is calculated using global min/max',...
'Tag','check_globalMinMax');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'Position',[449 572 134 17],...
'String','Shading Illumination',...
'Style','text',...
'Tag','text40',...
'UserData','shade');

uicontrol('Parent',h1,...
'Position',[920 190 101 61],...
'Style','frame',...
'Tag','frame7',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_globalWaterMax_CB'},...
'Position',[960 218 50 19],...
'Style','edit',...
'TooltipString','Global maximum water level',...
'Tag','edit_globalWaterMax',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[932 219 25 16],...
'String','Max',...
'Style','text',...
'Tag','textMax',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[932 199 25 16],...
'String','Min',...
'Style','text',...
'Tag','textMin',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_globalWaterMin_CB'},...
'Position',[960 197 50 19],...
'Style','edit',...
'TooltipString','Global minimum water level',...
'Tag','edit_globalWaterMin',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[925 243 90 17],...
'String','Global min/max',...
'Style','text',...
'Tag','text_globalMM',...
'UserData','t_grids');

uicontrol('Parent',h1,...
'Call',{@aquamoto_uiCB,h1,'push_checkGlobalMM_CB'},...
'FontName','Helvetica',...
'Position',[920 260 100 18],...
'String','Get global Min/Max',...
'TooltipString','Run trough all grids to find ensemble Min/Max',...
'Tag','push_checkGlobalMM',...
'UserData','t_grids',...
'Visible','off');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[690 50 311 21],...
'Style','edit',...
'TooltipString','Name of a file where the interpolation results are stored. First column is time.',...
'Tag','edit_profile',...
'UserData','misc');

uicontrol('Parent',h1, 'Position',[1000 50 21 21],...
'TooltipString','Browse for a file name',...
'Tag','push_profile',...
'UserData','misc');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[690 72 100 17],...
'String','Interpolated file',...
'Style','text',...
'Tag','text44',...
'UserData','misc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_maregs_CB'},...
'HorizontalAlignment','left',...
'Position',[690 140 311 21],...
'Style','edit',...
'TooltipString','Name of a file with x,y positions where to extract time series',...
'Tag','edit_maregs',...
'UserData','misc');

uicontrol('Parent',h1, 'Position',[1000 140 21 21],...
'Call',{@aquamoto_uiCB,h1,'push_maregs_CB'},...
'TooltipString','Browse for a file with x,y positions where to extract time series',...
'Tag','push_maregs',...
'UserData','misc');

uicontrol('Parent',h1, 'Position',[690 162 120 19],...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Positions (x,y) file',...
'Style','text',...
'Tag','text45',...
'UserData','misc');

uicontrol('Parent',h1, 'Position',[30 72 100 21],...
'Call',{@aquamoto_uiCB,h1,'push_runIn_CB'},...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','bold',...
'String','Plot Run In',...
'TooltipString','Plot a line delimiting the inundation zone',...
'Tag','push_runIn',...
'UserData','anuga');

uicontrol('Parent',h1, 'Position',[370 503 85 17],...
'Call',{@aquamoto_uiCB,h1,'radio_noShade_CB'},...
'FontName','Helvetica',...
'FontSize',9,...
'String','No shading',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_noShade',...
'UserData','shade');

uicontrol('Parent',h1, 'Position',[371 473 70 17],...
'Call',{@aquamoto_uiCB,h1,'radio_shade_CB'},...
'FontName','Helvetica',...
'FontSize',9,...
'String','Shading',...
'Style','radiobutton',...
'Tag','radio_shade',...
'UserData','shade');

uicontrol('Parent',h1, 'Position',[370 420 81 21],...
'Call',{@aquamoto_uiCB,h1,'push_apply_CB'},...
'String','Apply',...
'TooltipString','Apply changes made to the illumination settings',...
'Tag','push_apply',...
'UserData','shade');

uicontrol('Parent',h1, 'Position',[385 91 70 20],...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_multiLayerInc_CB'},...
'Enable','off',...
'String','1',...
'Style','edit',...
'TooltipString','Slice increment. ''1'' means all layers, ''2'', every other layer. The ML start:inc:end form is acepted as well',...
'Tag','edit_multiLayerInc',...
'UserData','cinema');

uicontrol('Parent',h1, 'Position',[353 94 30 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Slices',...
'Style','text',...
'Tag','text47',...
'UserData','cinema');

uicontrol('Parent',h1, 'Position',[355 115 125 16],...
'Call',{@aquamoto_uiCB,h1,'radio_multiLayer_CB'},...
'FontName','Helvetica',...
'String','Multi-layer source',...
'Style','radiobutton',...
'TooltipString','Build movie from a multilayer (netCDF) source file',...
'Tag','radio_multiLayer',...
'UserData','cinema');

uicontrol('Parent',h1, 'Position',[490 115 104 16],...
'Call',{@aquamoto_uiCB,h1,'radio_timeGridsList_CB'},...
'FontName','Helvetica',...
'String','Time grids list',...
'Style','radiobutton',...
'TooltipString','Build movie from a list of time grids (as loaded in "Time Grids" tab)',...
'Tag','radio_timeGridsList',...
'UserData','cinema');

uicontrol('Parent',h1, 'Position',[690 110 110 21],...
'Call',{@aquamoto_uiCB,h1,'push_interpolate_CB'},...
'Enable','off',...
'String','Interpolate',...
'TooltipString','Push this button to do the interpolation',...
'Tag','push_interpolate',...
'UserData','misc');

uicontrol('Parent',h1, 'Position',[841 111 50 21],...
'BackgroundColor',[1 1 1],...
'Call',{@aquamoto_uiCB,h1,'edit_miscLayerInc_CB'},...
'Enable','off',...
'String','1',...
'Style','edit',...
'TooltipString','Slice increment. ''1'' means all layers, ''2'', every other layer. The ML start:inc:end form is acepted as well',...
'Tag','edit_miscLayerInc',...
'UserData','misc');

uicontrol('Parent',h1, 'Position',[810 114 30 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Slices',...
'Style','text',...
'Tag','text48',...
'UserData','misc');

uicontrol('Parent',h1, 'Position',[910 112 105 19],...
'FontName','Helvetica',...
'String','Write header?',...
'Style','checkbox',...
'TooltipString','Write an header with: Variable name and coordinates of each point',...
'Value',1,...
'Tag','check_miscWriteHeader',...
'UserData','misc');

uicontrol('Parent',h1, 'Position',[140 450 121 81],...
'Call',{@aquamoto_uiCB,h1,'push_plugFun_CB'},...
'String','Run Plugin function',...
'TooltipString','Execute the external ''aquaPlugin'' function ',...
'Tag','push_plugFun',...
'UserData','plug');

uicontrol('Parent',h1, 'Position',[140 450 131 21],...
'String','Seek OPTcontrol.txt',...
'Style','checkbox',...
'TooltipString','Scan the OPTcontrol.txt file for a "Control script" file',...
'Tag','check_plugFun',...
'UserData','plug');

function aquamoto_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,[],guidata(h1));

function aquamoto_uiCB2(hObject, eventdata, h1, opt, callback_name)
	feval(callback_name,hObject,[], guidata(h1), opt);

% ----------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------
function varargout = vector_plot(varargin)
% 
	hObject = figure('Tag','figure1','Visible','off');
	vector_plot_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(handles.figure1,'center')

	handles.hAquaFig = varargin{1};			% Aquamoto fig handle
	handles.vecScale = varargin{2};
	handles.vecSpacing = varargin{3};
	plotVector = varargin{4};

	if (~plotVector)
		set(handles.radio_plotVec,'Val',0),		set(handles.radio_noPlot,'Val',1)
	else
		set(handles.radio_plotVec,'Val',1),		set(handles.radio_noPlot,'Val',0)
	end
	set(handles.edit_vecScale, 'String', handles.vecScale) 
	set(handles.edit_gridSpacing, 'String', sprintf('1:%d',handles.vecSpacing) )
	set(handles.slider_scale, 'Val', handles.vecScale) 
	set(handles.slider_spacing, 'Val', handles.vecSpacing)

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),		varargout{1} = hObject;		end

% ---------------------------------------------------------------------
function radio_plotVec_CB(hObject, eventdata, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_noPlot, 'Val', 0)
	handAqua = guidata(handles.hAquaFig);
	handAqua.plotVector = true;
	guidata(handles.hAquaFig, handAqua)

% ---------------------------------------------------------------------
function radio_noPlot_CB(hObject, eventdata, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_plotVec, 'Val', 0)
	handAqua = guidata(handles.hAquaFig);
	handAqua.plotVector = false;
	if (~isempty(handAqua.hQuiver))			% If we had plotted vectors, delete them
		delete(handAqua.hQuiver)
		handAqua.hQuiver = [];
	end
	guidata(handles.hAquaFig, handAqua)

% ---------------------------------------------------------------------
function slider_scale_CB(hObject, eventdata, handles)
	val = round(get(hObject,'Val'));
	set(handles.edit_vecScale, 'String', val)

% ---------------------------------------------------------------------
function slider_spacing_CB(hObject, eventdata, handles)
	val = round(get(hObject,'Val'));
	set(handles.edit_gridSpacing, 'String', sprintf('1:%d',val))
	handAqua = guidata(handles.hAquaFig);
	handAqua.spacingChanged = true;			% Signal quiver that old arrows are to delete
	guidata(handles.hAquaFig,handAqua)			% Save result in the Aquamoto handles

% ---------------------------------------------------------------------
function push_vp_OK_CB(hObject, eventdata, handles)
	if ( get(handles.radio_plotVec, 'Val') )
		vscale = round(get(handles.slider_scale,'Val'));
		vspacing = round(get(handles.slider_spacing,'Val'));
		handAqua = guidata(handles.hAquaFig);
		handAqua.vecScale = vscale;
		handAqua.vecSpacing = vspacing;
		handAqua.plotVector = logical(get(handles.radio_plotVec,'Val'));
		guidata(handles.hAquaFig,handAqua)			% Save result in the Aquamoto handles
	end
	delete(handles.figure1)

% ---------------------------------------------------------------------
function push_vp_cancel_CB(hObject, eventdata, handles)
	delete(handles.figure1)


% ------------ Creates and returns a handle to the GUI figure. --------
function vector_plot_LayoutFcn(h1)

set(h1,'Position',[520 620 161 180],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none', 'NumberTitle','off',...
'Resize','off', 'Name','Vector Plot',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[90 99 61 20],...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'String','1',...
'Style','edit',...
'TooltipString','Scale arrow sizes by this amount',...
'Tag','edit_vecScale');

uicontrol('Parent',h1,'Position',[90 45 61 20],...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'String','1:10',...
'Style','edit',...
'TooltipString','Plot 1 vector for every n grid nodes',...
'Tag','edit_gridSpacing');

uicontrol('Parent',h1,'Position',[11 101 65 15],...
'HorizontalAlignment','left',...
'String','Arrow Scale', 'Style','text');

uicontrol('Parent',h1,'Position',[10 46 70 15],...
'HorizontalAlignment','left',...
'String','Grid Spacing', 'Style','text');

uicontrol('Parent',h1,'Position',[10 9 60 23],...
'Call',{@vector_plot_uiCB,h1,'push_vp_OK_CB'},...
'FontName','Helvetica',...
'FontSize',9,...
'String','OK',...
'Tag','push_vp_OK');

uicontrol('Parent',h1,'Position',[91 9 60 23],...
'Call',{@vector_plot_uiCB,h1,'push_vp_cancel_CB'},...
'FontName','Helvetica',...
'FontSize',9,...
'String','Cancel',...
'Tag','push_vp_cancel');

uicontrol('Parent',h1,'Position',[10 156 79 15],...
'Call',{@vector_plot_uiCB,h1,'radio_plotVec_CB'},...
'String','Plot vectors',...
'Style','radiobutton',...
'Tag','radio_plotVec');

uicontrol('Parent',h1,'Position',[100 156 55 15],...
'Call',{@vector_plot_uiCB,h1,'radio_noPlot_CB'},...
'String','No plot',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_noPlot');

uicontrol('Parent',h1,'Position',[10 120 141 15],...
'BackgroundColor',[0.96 0.96 0.96],...
'Call',{@vector_plot_uiCB,h1,'slider_scale_CB'},...
'Style','slider',...
'Max',20, 'Min',1,...
'SliderStep',[0.05 0.1], 'Value',1,...
'Tag','slider_scale');

uicontrol('Parent',h1,'Position',[10 66 141 15],...
'BackgroundColor',[0.96 0.96 0.96],...
'Call',{@vector_plot_uiCB,h1,'slider_spacing_CB'},...
'Style','slider',...
'Max',100, 'Min',1,...
'SliderStep',[0.01 0.1], 'Value',10,...
'Tag','slider_spacing');

function vector_plot_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,[],guidata(h1));

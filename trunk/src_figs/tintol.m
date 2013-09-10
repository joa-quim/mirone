function tintol(handles,axis_t,X,Y,I)
% Transform a Mirone window into a Ground Control Points selection and
% Image registration tool. If the master image has coordinates this will work
% as Image-to-Map rectification, otherwise it works in the Image-to-Image mode.

%	Copyright (c) 2004-2013 by J. Luis
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

% 	delete([handles.Image handles.Tools handles.Plates handles.MagGrav ...
% 		handles.Seismology handles.Tsunamis handles.GMT handles.GridTools])

% 	delete([handles.NewFigure handles.ImportKnownTypes handles.LoadBGMap handles.SaveGMTgrid handles.Preferences ...
% 		handles.Print handles.DrawText handles.DrawGeogCirc handles.DrawLine handles.DrawRect ...
% 		handles.DrawPolyg handles.DrawArrow handles.ColorPal handles.Shading handles.Anaglyph ...
% 		handles.toGE handles.MBplaning handles.FlederPlanar handles.ImageInfo handles.Refresh])

	% ------------- Cleverer deletion of unwanted uicontrols
% 	h1 = get(handles.File,'Children');
% 	h2 = setxor(h1,handles.OpenGI);
% 	delete(h2)

	handles.origFig = [];		% We don't need the image copy anymore

	if (~handles.validGrid)
		handles.head = [-20 0 25 45 0 255 0 20/511 20/511];
		set(handles.figure1, 'Vis', 'off')
		mirone('FileNewBgFrame_CB',handles,[-20 0 25 45 1 0], [512 512], 'TINTOL');
		resizetrue(handles,[512 512], 'xy', [350 550]);

		handles.head_bat = [];
		handles.head_src = [];
		handles.maregraph_xy = [];
		handles.BothGridsInMemory = 0;  % Flag to signal that bat & deform arrays are already in memory
		handles.BatGridInMemory = 0;    % Flag to signal that bat array is already in memory
		handles.MaregraphInMemory = 1;  % Flag to signal that the maregraphs locations are already in memory
		handles.n_jump = 0;             % If > 0, the maregraphs array will start at n_jump + 1. OR netCDF files start
		handles.dt = [];
		handles.grn = 10;
		handles.cumint = 1;

		h = tintol_buttons_LayoutFcn(handles.figure1);
		handTintButt = local_guihandles(h);				% THIS WILL BE SEEN AS 'HANDLES' inside _CBs
		handTintButt.figure1 = handles.figure1;			% Make a copy so that we can fish the main handles in _CBs
		handTintButt.axes1 = handles.axes1;				% Simplify access to this
		handTintButt.home_dir = handles.home_dir;		% Start in values
		handTintButt.work_dir = handles.work_dir;
		handTintButt.last_dir = handles.last_dir;
		handTintButt.outPato = handles.last_dir;

% 		[X,Y,Z,head] = load_grd(handles);				% Copy Z/head into the nested_level cell array
% 		handTintButt.nested_level{1,1} = Z;
% 		handTintButt.nested_level{1,2} = head;
		handTintButt.last_nested_level = 1;

		if (handles.last_dir(end) ~= filesep),		handTintButt.outPato = [handles.last_dir filesep];	end

		local_guidata(handles.figure1, handTintButt)	% Store the local handles
		set(handles.figure1, 'Vis', 'on')
	end

%--------------------------------------------------------------------------------
function data = local_guidata(h, data_in)
% A local version of guidata that uses a different key for appdata
	fig = getParentFigure(h);
	if (nargin == 1)
		data = getappdata(fig, 'L_UsedByGUIData');
	else
		setappdata(fig, 'L_UsedByGUIData', data_in);
	end

%--------------------------------------------------------------------------------
function fig = getParentFigure(fig)
	while ~isempty(fig) && ~strcmp('figure', get(fig,'Type'))
	  fig = get(fig,'Parent');
	end

%--------------------------------------------------------------------------------
function handles = local_guihandles(hAll)
% Local version of guihandles that builds a handles structure from a
% a vector list of uicontrols handles.
	handles = [];

	for i = (1:numel(hAll))
		if isprop(hAll(i), 'Tag')
			tag = get(hAll(i), 'Tag');
			if (~isempty(tag))
				if isfield(handles, tag)
					prev_h = handles.(tag);
				else
					prev_h = [];
				end
				handles.(tag) = [prev_h hAll(i)];	% In case a previous handle with same name already exists
			end
		end
	end

%--------------------------------------------------------------------------------
function push_nswing_CB(hObject, handles)


%--------------------------------------------------------------------------------
function push_anuga_CB(hObject, handles)

%--------------------------------------------------------------------------------
function push_changeR_CB(hObject, handles)
% Remember, this handles is NOT the Mirone handles a local one so we need
% to fetch the old Mirone one to transmit to its function.
	handMirMain = guidata(handles.figure1);
	mirone('FileNewBgFrame_CB', handMirMain)

%--------------------------------------------------------------------------------
function edit_SourceGrid_CB(hObject, handles)
	fname = get(hObject,'String');
	if (isempty(fname)),	return,		end
	push_SourceGrid_CB(handles.push_SourceGrid, handles, fname)

%--------------------------------------------------------------------------------
function push_SourceGrid_CB(hObject, handles, opt)
	if (nargin == 2)
		str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
		[FileName,PathName,handles] = put_or_get_file(handles,str1,'Select grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
		set(handles.edit_SourceGrid,'Str',fname)
	else
		fname = opt;
	end

	[handles,X,Y,handles.Z_src,handles.head_src] = read_gmt_type_grids(handles,fname);
	if (isempty(X)),    return,		end
	%geog = guessGeog(handles.head_src(1:4));
	%set(handles.radio_cartesian,'Value', ~double(geog))
	%set(handles.radio_geog,'Value', double(geog))
	local_guidata(handles.figure1,handles)

%--------------------------------------------------------------------------------
function edit_NestGrids_CB(hObject, handles)
	fname = get(hObject,'String');
	if (isempty(fname)),	return,		end
	push_NestGrids_CB(handles.push_NestGrids, handles, fname)

%--------------------------------------------------------------------------------
function push_NestGrids_CB(hObject, handles, opt)
% Read either base level (level 0) or descendeny nesting grids
	if (nargin == 2)
		str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
		[FileName,PathName,handles] = put_or_get_file(handles,str1,'Select grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
		set(handles.edit_NestGrids,'Str',fname)
	else
		fname = opt;
	end

	str = get(handles.popup_nestings,'Str');
	val = get(handles.popup_nestings,'Val');
	if (val > handles.last_nested_level + 1)		% User Screw up, probably a joke/esticanço
		errordlg('You cannot jump steps. Be modest and select one nesting level at a time.','Chico Clever')
		return
	elseif (val == handles.last_nested_level)
		val = val + 1;
		handles.last_nested_level = val;
	elseif (val == 1 && (numel(str{1}) == 1 && str{1} == '0'))
		handMain = guidata(handles.figure1);
		[X,Y,Z,head_bat] = load_grd(handMain,'silent');
		if (isempty(X))
			errordlg('You can not read a nesting grid before the base level. Use "File -> Open" to do that.','Error')
			return
		else
			% Since is complicated to do this when "File->Open" we do it here (saving base level)
			handles.nested_level{1,1} = Z;
			handles.nested_level{1,2} = head;	
			str{1} = '0 -- level ready to use';
			val = val + 1;
			handles.last_nested_level = val;
		end
	end		% ELSE just replace one previous level (DUMB and DANGEROUS thing to do -- should not be allowed)

	[handles,X,Y,Z,head] = read_gmt_type_grids(handles, fname);
	if (isempty(X)),    return,		end

	handles.nested_level{val,1} = Z;
	handles.nested_level{val,2} = head;

	str{val} = sprintf('%d -- level ready to use', val-1);
	set(handles.popup_nestings, 'Str', str, 'Val', val);		% Make clear what nesting level we are on

	% Draw the BoundingBox of this nested grid
	x = [head(1) head(1) head(2) head(2) head(1)];
	y = [head(3) head(4) head(4) head(3) head(3)];
	line('parent',handles.axes1, 'XData',x, 'YData',y);

	local_guidata(handles.figure1,handles)

%--------------------------------------------------------------------------------
function popup_nestings_CB(hObject, handles)
% Show what has been loaded so far or prepare for the -- eventual -- next level.
	%str = get(hObject,'Str');

%--------------------------------------------------------------------------------
function radio_outGrids_CB(hObject, handles)
	if (~get(hObject,'Value'))
		set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
				handles.check_velocity handles.check_momentum], 'Enable','off')
		set(handles.edit_gridNameStem, 'String','', 'Tooltip','')
		return
	end
	set([handles.radio_anuga handles.radio_most], 'Value',0)
	set(handles.edit_gridNameStem,'String',[handles.outPato 'tsu_time_'])
	set(handles.edit_gridNameStem,'Tooltip','Grids are numbered using this name stem')
	set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
			handles.check_velocity handles.check_momentum], 'Enable','on')
	set(handles.radio_surfLevel, 'Val', 1)

%--------------------------------------------------------------------------------
function radio_anuga_CB(hObject, handles)
	if (~get(hObject,'Value'))
		set(handles.edit_gridNameStem, 'String','', 'Tooltip','')
		return
	end
	set([handles.radio_outGrids handles.radio_most], 'Value',0)
	set(handles.edit_gridNameStem,'String',[handles.outPato 'swan_tsu.sww'])
	set(handles.edit_gridNameStem,'Tooltip','Output file name')
	set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
			handles.check_velocity handles.check_momentum], 'Val',0)
	set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
			handles.check_velocity handles.check_momentum], 'Enable','off')

%--------------------------------------------------------------------------------
function radio_most_CB(hObject, handles)
	if (~get(hObject,'Value'))
		set(handles.edit_gridNameStem, 'String','', 'Tooltip','')
		return
	end
	set([handles.radio_outGrids handles.radio_anuga], 'Value',0)
	set(handles.edit_gridNameStem,'String',[handles.outPato 'swan_tsu'])
	set(handles.edit_gridNameStem,'Tooltip','Basename for MOST files (do not give extension)')
	set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
			handles.check_velocity handles.check_momentum], 'Val',0)
	set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
			handles.check_velocity handles.check_momentum], 'Enable','off')

%--------------------------------------------------------------------------------
function edit_gridNameStem_CB(hObject, handles)


%--------------------------------------------------------------------------------
function radio_surfLevel_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_totalWater handles.radio_maxWater],'Val',0)

%--------------------------------------------------------------------------------
function radio_totalWater_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_surfLevel handles.radio_maxWater],'Val',0)

%--------------------------------------------------------------------------------
function radio_maxWater_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_surfLevel handles.radio_totalWater],'Val',0)

%--------------------------------------------------------------------------------
function check_velocity_CB(hObject, handles)


%--------------------------------------------------------------------------------
function check_momentum_CB(hObject, handles)


%--------------------------------------------------------------------------------
function check_wantMaregs_CB(hObject, handles)
	if (get(hObject,'Value'))
		set([handles.edit_MaregraphPosFile handles.push_MaregraphPosFile],'Enable','on')
		set([handles.edit_MaregraphDataFile handles.push_MaregraphDataFile],'Enable','on')
	else
		set(handles.edit_MaregraphPosFile,'String','','Enable','off')
		set(handles.edit_MaregraphDataFile,'String','','Enable','off')
		set(handles.push_MaregraphPosFile,'Enable','off')
		set(handles.push_MaregraphDataFile,'Enable','off')
	end

%--------------------------------------------------------------------------------
function edit_cumint_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)
		set(hObject,'String',handles.cumint),	return
	end
	handles.cumint = xx;
	if (handles.got_params),	handles.params(9) = xx;		end
	local_guidata(handles.figure1,handles)

%--------------------------------------------------------------------------------
function edit_MaregraphPosFile_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),  return,		end
	push_MaregraphPosFile_CB(handles.push_MaregraphPosFile, handles, fname)

%--------------------------------------------------------------------------------
function push_MaregraphPosFile_CB(hObject, handles)
	if (nargin == 3)
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.dat;*.DAT;*.xy', 'Maregraph location (*.dat,*.DAT,*.xy)';'*.*', 'All Files (*.*)'},'Select Maregraphs position','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	else
		fname = opt;
	end

	[handles, msg] = getMaregsPos(handles, fname);
	if (~isempty(msg))
		errordlg(msg,'Error')
		return
	end
	
	set(handles.edit_Number_of_cycles,'String',size(handles.maregraph_xy,1))
	set(handles.edit_MaregraphPosFile,'String',fname)
	local_guidata(handles.figure1,handles)

%--------------------------------------------------------------------------------
function edit_MaregraphDataFile_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),  return,		end
	push_MaregraphDataFile_CB(handles.push_MaregraphDataFile, handles, fname)

%--------------------------------------------------------------------------------
function push_MaregraphDataFile_CB(hObject, handles)
	if (nargin == 3)
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.dat;*.DAT;*.xy', 'Maregraph data file (*.dat,*.DAT,*.xy)';'*.*', 'All Files (*.*)'},'Select Maregraph','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
		
		if (~handles.is_tsun2)		% We are done with this case
			set(handles.edit_MaregraphDataFile,'String',fname)
			return
		end
	else
		fname = opt;
	end

	% If we get here we are in the tsun2 oputput maregraphs option. We must read that file
	bak = handles.maregraph_xy;		% To reuse the code from getMaregsPos() we must backup this
	[handles, msg] = getMaregsPos(handles, fname);
	if (~isempty(msg))
		errordlg(msg,'Error')
		handles.maregraph_xy = bak;			% Get back the original data
		local_guidata(handles.figure1,handles)
		return
	end
	% OK, we are still playing
	handles.maregraph_xy = bak;			% Get back the original data
	set(handles.edit_MaregraphDataFile,'String',fname)
	local_guidata(handles.figure1,handles)

%--------------------------------------------------------------------------------
function edit_Number_of_cycles_CB(hObject, handles)
% The OK button will get the number of cycles from here. So ensure it gives a valid number
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 2)
		set(hObject,'String','1010')
	else
		set(hObject,'String',sprintf('%.0f',xx)) % Ensure that it is an integer
	end

%--------------------------------------------------------------------------------
function edit_jumpInitial_CB(hObject, handles)
% Remove this number of seconds from the beguining of the maregraphs (if they were loaded)
	jmp = str2double(get(hObject,'String'));
	if (isnan(jmp) || jmp < 0)
		set(hObject,'String','0');		handles.n_jump = 0;
		local_guidata(handles.figure1,handles)
		return
	end

	if ( ~isempty(handles.maregraph_xy) )
		m = size(handles.maregraph_xy,1);
     	dt1 = diff(handles.maregraph_xy(:,1));    dt2 = diff(dt1);
		if any(dt2 ~= 0)				% First column doesn't have the time
		else							% First column has the time.
			dt = dt1(1);				% This is the time increment
			n_jmp = fix(jmp / dt);		% Number of records to jump
			if ((m - n_jmp) < 10)		% Stupid choice for jump time
				set(hObject,'String','0')
				n_jmp = 0;
			end
			set(handles.edit_Number_of_cycles,'String',sprintf('%d',m-n_jmp))    % Uppdate max pssible
			handles.n_jump = n_jmp;
		end
	else
		handles.n_jump = jmp;			% To be used in ANUGA or MOST outputs. It's a jump time, not cycles
	end
	local_guidata(handles.figure1,handles)

%--------------------------------------------------------------------------------
function edit_dt_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)
		set(hObject,'String',handles.dt),	return
	end
	% Check that xx is a valid (CFL condition) dt
	dx = handles.head_bat(9);
	if (handles.geog),		dx = dx * 111000;		end
	dt = dx / sqrt(abs(handles.head_bat(5)) * 9.8);
	if ( xx > dt )
		warndlg(['The value entered here doesn''t conform to the CFL condition (' ...
				sprintf('%.3f',dt) ') and will probably make the simulation diverge.'],'Warning')
	end
	handles.dt = xx;
	if (handles.got_params),	handles.params(1) = xx;		end
	local_guidata(handles.figure1,handles)

%--------------------------------------------------------------------------------
function edit_grn_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)
		set(hObject,'String',handles.grn),	return
	end
	handles.grn = xx;
	if (handles.got_params),	handles.params(2) = xx;		end
	local_guidata(handles.figure1,handles)

%--------------------------------------------------------------------------------
function push_RUN_CB(hObject, handles)
% ...
	[X,Y,Z,head_bat] = load_grd(handles,'silent');
	if (~isempty(X))		% So we have a bat grid. Check the remaining conditions for this case.
		err_str = check_errors(handles);
		if (~isempty(err_str))
			errordlg(err_str, 'Error'),		return
		end
	else
		% Print something
		return
	end

	
%---------------------------------------------------------------------------------------------------
function err_str = check_errors(handles)
% ...
	err_str = '';

	[X,Y,Z,head_bat] = load_grd(handles,'silent');

	if (isempty(handles.head_src))
		err_str = 'No tsunami source grid provided. And I NEED ONE!';	return
	end

	% Cheeck that at lest one operation was selected
	if ( ~(get(handles.check_wantMaregs,'Val') || get(handles.radio_outGrids,'Val') || ...
			get(handles.radio_anuga,'Val') || get(handles.radio_most,'Val')) )
		err_str = 'Comput what? You need to select at least one thing to do';
		return
	end

	if (get(handles.check_wantMaregs,'Value'))       % That is, requested maregraphs computation
		if (isempty(handles.maregraph_xy) && handles.MaregraphInMemory == 0)
			err_str = 'Where are your maregraphs? On the Moon?';
			return
		end
		if (isempty(get(handles.edit_MaregraphDataFile,'String')))
			err_str = 'You need to tell me the file name where I''ll write the maregraphs water hight';
			return
		end
	end

	% If we reach here, we can now check if grids are compatible.
	small = 1e-6;   % Used in postion comparation. It places the accuracy at sub-meter level
	difes = handles.head_src(1:4) - head_bat(1:4);
	if (any(abs(difes) > small))
		err_str{1} = 'Bathymetry & Source grids do not cover the same region';
		err_str{2} = sprintf('x_min diff = %.10g', difes(1));	err_str{3} = sprintf('x_max diff = %.10g', difes(2));
		err_str{4} = sprintf('y_min diff = %.10g', difes(3));	err_str{5} = sprintf('y_max diff = %.10g', difes(4));
		return
	end
	if (abs(handles.head_src(8) - head_bat(8)) > small || ...
		abs(handles.head_src(9) - head_bat(9)) > small )
		err_str = 'Bathymetry & Source grids have different nx and/or ny';
	end

% ----------------------------------------
function h = tintol_buttons_LayoutFcn(h1)

h = zeros(36);		% Change this if uis are added or removed
n = 1;
h(n) = uicontrol('Parent',h1, 'Position',[10 512 66 23],...
'Call',@tintolButtons_uiCB,...
'String','NSWING',...
'UserData','main',...
'Tag','push_nswing');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[80 512 80 23],...
'Call',@tintolButtons_uiCB,...
'String','ANUGA',...
'HitTest','off',...
'UserData','params',...
'Tag','push_anuga');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[0 140 348 375],...
'Enable','inactive',...
'Tag','push_bg');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[49 490 110 22],...
'Call',@tintolButtons_uiCB,...
'String','Change Region',...
'Tag','push_changeR');	n = n + 1;


% ------------------------- FRAME 1 -------------------------
h(n) = uicontrol('Parent',h1, 'Position',[10 400 331 100],...
'Style','frame',...
'UserData','main',...
'Tag','frame1');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[11 464 38 15],...
'HorizontalAlignment','right',...
'String','Source',...
'Style','text',...
'UserData','main',...
'Tag','text_Source');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[51 463 268 21],...
'Call',@tintolButtons_uiCB,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Style','edit',...
'UserData','main',...
'Tag','edit_SourceGrid');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[316 462 23 23],...
'Call',@tintolButtons_uiCB,...
'UserData','main',...
'Tag','push_SourceGrid');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[11 439 37 15],...
'HorizontalAlignment','right',...
'String','Nest ',...
'Style','text',...
'UserData','main',...
'Tag','text_SNests');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[51 437 268 21],...
'Call',@tintolButtons_uiCB,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Style','edit',...
'UserData','main',...
'Tag','edit_NestGrids');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[316 436 23 23],...
'Call',@tintolButtons_uiCB,...
'UserData','main',...
'Tag','push_NestGrids');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[50 410 282 20],...
'Call',@tintolButtons_uiCB,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Style','popupmenu',...
'String',{'0','1','2','3','4','5'},...
'Tooltip','Select nesting grid level',...
'UserData','main',...
'Tag','popup_nestings');	n = n + 1;
% -----------------------------------------------------------

% ------------------------- FRAME 2 -------------------------
h(n) = uicontrol('Parent',h1, 'Position',[10 350 331 40],...
'Style','frame',...
'UserData','main',...
'Tag','frame2');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[49 360 200 22],...
'Call',@tintolButtons_uiCB,...
'String','Bordering',...
'Tag','push_bordering');	n = n + 1;
% -----------------------------------------------------------

% ------------------------- FRAME 3 -------------------------
h(n) = uicontrol('Parent',h1, 'Position',[10 280 331 63],...
'Style','frame',...
'UserData','main',...
'Tag','frame3');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[21 318 89 15],...
'Call',@tintolButtons_uiCB,...
'String','Output grids',...
'Style','radiobutton',...
'Tooltip','Save output  as a series of Surfer format grids',...
'UserData','main',...
'Tag','radio_outGrids');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[131 318 99 15],...
'Call',@tintolButtons_uiCB,...
'String','ANUGA .sww',...
'Style','radiobutton',...
'Tooltip','Save output as a single netCDF file using the ANUGA''s .sww format',...
'UserData','main',...
'Tag','radio_anuga');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[251 318 76 15],...
'Call',@tintolButtons_uiCB,...
'String','MOST .nc',...
'Style','radiobutton',...
'Tooltip','Save output as 3 netCDF files as produced by MOST',...
'UserData','main',...
'Tag','radio_most');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[20 292 30 16],...
'HorizontalAlignment','left',...
'String','Name',...
'Style','text',...
'UserData','main',...
'Tag','text_Name');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[50 290 271 21],...
'Call',@tintolButtons_uiCB,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Grids are created using this name stem',...
'UserData','main',...
'Tag','edit_gridNameStem');	n = n + 1;
% -----------------------------------------------------------

h(n) = uicontrol('Parent',h1, 'Position',[9 254 93 15],...
'Call',@tintolButtons_uiCB,...
'Enable','off',...
'String','Surface level',...
'Style','radiobutton',...
'UserData','main',...
'Tag','radio_surfLevel');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[248 254 82 15],...
'Call',@tintolButtons_uiCB,...
'Enable','off',...
'String','Total water',...
'Style','radiobutton',...
'UserData','main',...
'Tag','radio_totalWater');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[130 254 82 15],...
'Call',@tintolButtons_uiCB,...
'Enable','off',...
'String','Max water',...
'Style','radiobutton',...
'UserData','main',...
'Tag','radio_maxWater');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[130 234 82 15],...
'Call',@tintolButtons_uiCB,...
'Enable','off',...
'String','Velocity',...
'Style','checkbox',...
'Tooltip','Write velocity grids (u and v with sufixes _U, _V)',...
'UserData','main',...
'Tag','check_velocity');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[248 234 82 15],...
'Call',@tintolButtons_uiCB,...
'Enable','off',...
'String','Momentum',...
'Style','checkbox',...
'Tooltip','Write momentum grids (with sufixes _Uh, _Vh)',...
'UserData','main',...
'Tag','check_momentum');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[10 202 87 15],...
'Call',@tintolButtons_uiCB,...
'String','Maregraphs',...
'Style','checkbox',...
'Tooltip','Check this if you want to compute water height at maregraphs',...
'UserData','main',...
'Tag','check_wantMaregs');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[109 202 180 16],...
'Enable','off',...
'HorizontalAlignment','right',...
'String','Saving step (time = Time step * this)',...
'Style','text',...
'Tag','text_SveStepTime');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[290 200 31 21],...
'Call',@tintolButtons_uiCB,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'String','1',...
'Style','edit',...
'Tooltip','cumint * dt = time at which the virtual tide gauges are writen',...
'UserData','params',...
'Tag','edit_cumint');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[50 175 271 21],...
'Call',@tintolButtons_uiCB,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Name of the file with maregraph locations',...
'UserData','main',...
'Tag','edit_MaregraphPosFile');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[321 174 23 23],...
'Call',@tintolButtons_uiCB,...
'UserData','main',...
'Tag','push_MaregraphPosFile');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[50 149 271 21],...
'Call',@tintolButtons_uiCB,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Name of the file that will contain the maregraphs water height',...
'UserData','main',...
'Tag','edit_MaregraphDataFile');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[321 148 23 23],...
'Call',@tintolButtons_uiCB,...
'UserData','main',...
'Tag','push_MaregraphDataFile');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[12 179 36 16],...
'HorizontalAlignment','right',...
'String','In file',...
'Style','text',...
'UserData','main',...
'Tag','text_InFile');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[10 152 38 16],...
'HorizontalAlignment','right',...
'String','Out file',...
'Style','text',...
'UserData','main',...
'Tag','text_OutFile');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[80 109 51 21],...
'Call',@tintolButtons_uiCB,...
'BackgroundColor',[1 1 1],...
'String','1010',...
'Style','edit',...
'Tooltip','Use this number of cycles from the Maregraph file',...
'Tag','edit_Number_of_cycles' );	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[280 109 41 21],...
'Call',@tintolButtons_uiCB,...
'BackgroundColor',[1 1 1],...
'String','0',...
'Style','edit',...
'Tooltip','Output is produced only after this modeling time (seconds) has elapsed',...
'Tag','edit_jumpInitial');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[90 79 41 21],...
'Call',@tintolButtons_uiCB,...
'BackgroundColor',[1 1 1],...
'Style','edit',...
'Tooltip','Time step for the modeling run',...
'UserData','params',...
'Tag','edit_dt');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[280 79 41 21],...
'Call',@tintolButtons_uiCB,...
'BackgroundColor',[1 1 1],...
'String','10',...
'Style','edit',...
'Tooltip','save grids at this cycle interval (it depends on dt)',...
'UserData','params',...
'Tag','edit_grn');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[16 113 63 15],...
'HorizontalAlignment','right',...
'String','Nº of cycles',...
'Style','text',...
'Tag','text_Ncycles');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[221 113 57 15],...
'HorizontalAlignment','right',...
'String','Jump initial ',...
'Style','text',...
'Tag','text_JumpInit');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[7 81 81 16],...
'HorizontalAlignment','right',...
'String','Time step (sec)',...
'Style','text',...
'UserData','params',...
'Tag','text_TimeStep');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[151 81 126 16],...
'HorizontalAlignment','right',...
'String','Saving step (cycle units)',...
'Style','text',...
'Tag','text_SaveCycle');	n = n + 1;

h(n) = uicontrol('Parent',h1, 'Position',[260 30 81 21],...
'Call',@tintolButtons_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'String','RUN',...
'Tag','push_RUN' );

function tintolButtons_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, local_guidata(hObject));

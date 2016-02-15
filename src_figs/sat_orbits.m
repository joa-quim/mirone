function varargout = sat_orbits(varargin)
% Helper window to call the SGP4 orbit propagator (calculate Satellite orbits)

%	Copyright (c) 2004-2016 by J. Luis
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

% $Id: sat_orbits.m 7799 2016-02-15 18:24:26Z j $

% For compiling one need to include the orbits.m file.

	hObject = figure('Vis','off');
	sat_orbits_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right');

	if (nargin == 0),	varargin{1} = [];	end		% So por agora
	handles.handMir = varargin{1};
	handles.rect = [];
	if (nargin == 2)		% Second arg must be an handle to a rectangle, but not tested.
		x = get(varargin{2}, 'XData');
		y = get(varargin{2}, 'YData');
		handles.rect = [min(x) max(x) min(y) max(y)];
	end
	
	handles.log_pass = [];		% To store login credentials valid for one session.

	% By default select current day
	t = date;
	set(handles.edit_dateStart, 'Str', [t ' 00:00:00'])
	set(handles.edit_dateStop , 'Str', [t ' 23:59:59'])

	% See if we have traces in preferences of the last used TLE file
	s = load([handles.handMir.path_data 'mirone_pref.mat'], 'lastTLE');	% Warning: This is NEVER empty ... even when it is
	if (isfield(s,'lastTLE') && exist(s.lastTLE,'file'))		% Yes, we found one
		set(handles.edit_TLE, 'Str', s.lastTLE)
	end

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) ------------------------------
	
	% ------- Add this figure handle to the carraças list ---------
	plugedWin = getappdata(handles.handMir.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.handMir.figure1,'dependentFigs',plugedWin);

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% ------------------------------------------------------------------------
function edit_dateStart_CB(hObject, handles)
% If manually set check likelyhood.
	d  = datenum(get(hObject, 'Str'));
	if (abs(now - d) > 10000)	% ~30 yrs!
		warndlg('I won''t censure this date but it''s highly probable that you are inventing.','Fiu Fiu')
	end

% ------------------------------------------------------------------------
function push_calendarStart_CB(hObject, handles)
% Select starting orbit time
	new_date = uisetdate;
	set(handles.edit_dateStart, 'Str', new_date)

% ------------------------------------------------------------------------
function edit_dateStop_CB(hObject, handles)
% If manually set check likelyhood.
	d  = datenum(get(hObject, 'Str'));
	if (abs(now - d) > 10000)	% ~30 yrs!
		warndlg('I won''t censure this date but it''s highly probable that you are inventing.','Fiu Fiu')
	end

% ------------------------------------------------------------------------
function push_calendarStop_CB(hObject, handles)
	new_date = uisetdate;
	set(handles.edit_dateStop, 'Str', new_date)

% ------------------------------------------------------------------------
function edit_TLE_CB(hObject, handles, fname)
% Manually entered file or via push_getTLE_CB(). Also saves the file name in preferences

	lastTLE = get(handles.edit_TLE, 'Str');
	if (check_TLE(lastTLE))			% Check TLE looks good/exists
		set(hObject, 'Str', '')
		errordlg('Non existing or badly formated TLE file. Ignoring it.','Error')
		return
	end

	% Ok, since we got here we will accept the TLE as good and save its name in preferences
	version7 = version;
	V7 = (sscanf(version7(1),'%f') > 6);
	if (~V7)		% R <= 13. That's all I need to know for now.
		save([handles.handMir.path_data 'mirone_pref.mat'], 'lastTLE', '-append')
	else
		save([handles.handMir.path_data 'mirone_pref.mat'], 'lastTLE', '-append', '-v6')
	end

% ------------------------------------------------------------------------
function push_getTLE_CB(hObject, handles)
% ...
	str1 = {'*.tle;*.TLE;*.txt;*.TXT;*.dat;*.DAT;', 'TLE files (*.tle,*.txt,*.dat)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select TLE File','get');
	if isequal(FileName,0),		return,		end
	set(handles.edit_TLE, 'Str', [PathName FileName])
	edit_TLE_CB(handles.edit_TLE, handles, [PathName FileName])

% ------------------------------------------------------------------------
function push_callSpaceTrack_CB(hObject, handles)
% Send the browser to the SpaceTrack site
	url = 'https://www.space-track.org/#/tle';
	web(url, '-browser');

% ------------------------------------------------------------------------
function push_gotoSpaceTrack_CB(hObject, handles)
% Select from a known (to me) satellite and get (try) the TLE via SpaceTrack API

	val = get(handles.popup_satellite, 'Val');
	if (val == 1)
		errordlg('Cm''on select a satellite. It''s not hard, is it?', 'Error')
		return
	end

	if (isempty(handles.log_pass))
		resp = inputdlg({'login: (you must have an account at www.space-track.org)' 'password'}, 'Space-Track login data',[1 60]);
		log = resp{1};		pass = resp{2};
		if (isempty(log) || isempty(pass))
			msgbox('OK, bye, bye', 'No I don''t have one.'),		return
		end
	end

	str = get(handles.popup_satellite, 'Str');
	if (strncmp(str{val}, 'AQUA', 4))
		ID = '27424';
		fname = [handles.handMir.path_data 'example_data/AQUA.tle'];
	elseif (strncmp(str{val}, 'TERRA', 5))
		ID = '25994';
		fname = [handles.handMir.path_data 'example_data/TERRA.tle'];
	end
	if (ispc)
		cmd = ['wget --post-data "identity=' log '&password=' pass '&query=' ...
				'https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/NORAD_CAT_ID/' ...
				ID ...
				'/orderby/TLE_LINE1 ASC/format/3le" --keep-session-cookies --no-check-certificate --save-cookies=cookies.txt' ...
				' https://www.space-track.org/ajaxauth/login -O ' fname];
		dos(cmd)
	else
		cmd = ['wget --post-data ''identity=' log '&password=' pass '&query=' ...
				'https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/NORAD_CAT_ID/' ...
				ID ...
				'/orderby/TLE_LINE1 ASC/format/3le'' --keep-session-cookies --no-check-certificate --save-cookies=cookies.txt' ...
				' ''https://www.space-track.org/ajaxauth/login'' -O ' fname];
		unix(cmd)
	end

	good = check_TLE(fname);
	if (~good)
		warndlg('Sorry, seams that the direct TLE download has failed. Either try again or do it manually.','Warning')
	else
		set(handles.edit_TLE, 'Str', fname)
		edit_TLE_CB(handles.edit_TLE, handles, fname)
		handles.log_pass = {log pass};
		guidata(handles.figure1, handles)	% Save credentials that will be valid for this session
	end

% ------------------------------------------------------------------------
function popup_satellite_CB(hObject, handles)
% ...
	val = get(hObject, 'Val');
	if (val > 1 ),		set(handles.push_tracksTiles, 'Vis', 'on')
	else				set(handles.push_tracksTiles, 'Vis', 'off')
	end

% ------------------------------------------------------------------------
function push_seeGE_CB(hObject, handles)
	warndlg('Ah,ah! not yet','')

% ------------------------------------------------------------------------
function push_tracksTiles_CB(hObject, handles)
% Plot the scenes as patches

	[t_start, t_stop, fname, msg] = check_input(handles);
	if (~isempty(msg)),		errordlg(msg, 'Error'),		return,		end

	if (isempty(handles.rect))	% If not called via a rectangle clicked-callback
		BB = [get(handles.handMir.axes1, 'XLim') get(handles.handMir.axes1, 'YLim')] + [-10 10 -10 10];	% Plus a pad
	else
		BB = handles.rect;
	end
	if (BB(2) - BB(1) > 359),	BB = [];	end		% If global image, no BoundingBox
	tracks = orbits(fname, t_start, t_stop, 0.5, BB);
	if (isempty(tracks))
		warndlg('No track data inside this region.', 'Warning'),	return
	end

	% --------------- Retain only the DAY or NIGHT parts of the orbits -----------------------
	str = get(handles.popup_satellite, 'Str');		val = get(handles.popup_satellite, 'Val');
	if (strfind(str{val}, '(day)'))
		ind = (tracks.date(:,4) < 7 | tracks.date(:,4) > 19);
	else
		ind = (tracks.date(:,4) > 7 | tracks.date(:,4) < 19);
	end
	Sat = str{val}(1);				% First char in sat name, used in uicontext to build file name
	tracks.date(ind,:) = [];		tracks.xyz(ind,:) = [];
	if (isnan(tracks.xyz(1,1)))		% Shit, bad luck with the cut above. We don't want to start with a NaN
		tracks.xyz(1, :)  = [];		tracks.date(1, :)  = [];
	elseif (isnan(tracks.xyz(2,1)))	% Eve more shity bad luck. We can't have it either in 2 row (because azims)
		tracks.xyz(1:2, :)  = [];	tracks.date(1:2, :)  = [];
	end
	% ----------------------------------------------------------------------------------------

	% Now suffer to compute the tiles
	rng = 2326958 / 2;			% This is if for AQUA (and got by measuring over a L2 grid

	ind_NaN = [0; find(isnan(tracks.xyz(:,1))); size(tracks.xyz,1)+1];		% First and last pts are for algorithmic reasons
	n_tracks = numel(ind_NaN) - 1;	% True number of tracks (number of NaNs + 1)

	for (nt = 1:n_tracks)
		new_track = true;
		cor = rand(1,3);
		x = tracks.xyz(ind_NaN(nt)+1:ind_NaN(nt+1)-1,1);
		y = tracks.xyz(ind_NaN(nt)+1:ind_NaN(nt+1)-1,2);
		data = tracks.date(ind_NaN(nt)+1:ind_NaN(nt+1)-1, :);		% Means 'date' in Tuguese
		data = datevec(datenum(data));	% Round trip to get rid of times like 29 min 60 sec which is probably an orbits() bug
		azims = azimuth_geo(y(1:end-1), x(1:end-1), y(2:end), x(2:end));
		t = data(:,5) + data(:,6) / 60;
		t = str2num(sprintf('%.1f\n', t)) * 10;		% Trick to round to one decimal only. Times 10 to use in rem()
		ind = find(rem(t, 50) == 0);	% Find all the times multiples of 5 min, which mark the start of a new scene
		if (ind(end) ~= numel(x)),	ind(end+1) = numel(x);	end		% To plot also the chunk of last (partially outside) patch
		for (k = 1:numel(ind)-1)
			if (new_track)				% Deal with the cases where a patch is partially in but starts outside.
				x1 = x(1);				y1 = y(1);
				x2 = x(ind(k));			y2 = y(ind(k));
				d = data(1,:);
				r = rem(d(5), 5);		% First patch only by chance has the correct minute we need for building file's name.
				if (r ~= 0),	d(5) = d(5) - r;	end		% This insures we have the 'right minute' for file's name.
				plot_patches(handles, x1, y1, x2, y2, azims(1), azims(ind(k)), rng, cor, Sat, d)
				new_track = false;
			end

			dmin = t(ind(k+1)) - t(ind(k)); 
			if ((dmin == 50 || dmin == -550) && (data(ind(k+1),4) - data(ind(k),4)) <= 1)
				x1 = x(ind(k));			y1 = y(ind(k));
				x2 = x(ind(k+1));		y2 = y(ind(k+1));
				az1 = azims(ind(k));			az2 = azims(ind(k+1));
			else				% Patches that End outside the displayed region, but we want those too.
				x1 = x(ind(k));			y1 = y(ind(k));
				x2 = x(end);			y2 = y(end);
				az1 = azims(ind(k));	az2 = azims(end);	% az2 is out of view so don't care to be very accurate.
				new_track = true;
			end
			plot_patches(handles, x1, y1, x2, y2, az1, az2, rng, cor, Sat, data(ind(k),:))
			if (new_track)
				cor = rand(1,3);		% New track's new color
			end
		end
	end

% ------------------------------------------------------------------------
function plot_patches(handles, x1, y1, x2, y2, az1, az2, rng, cor, Sat, data)
% AZ1 and AZ2 are the track's azimuth at location X1,Y1 and X2,Y2

	[lat1,lon1] = vreckon(y1, x1, rng, az1+90, 1);
	[lat4,lon4] = vreckon(y1, x1, rng, az1-90, 1);

	[lat2,lon2] = vreckon(y2, x2, rng, az2+90, 1);
	[lat3,lon3] = vreckon(y2, x2, rng, az2-90, 1);
	h = patch('XData',[lon1 lon2 lon3 lon4 lon1], 'YData',[lat1 lat2 lat3 lat4 lat1], 'parent', handles.handMir.axes1, ...
		'FaceColor',cor, 'EdgeColor',handles.handMir.DefLineColor,'LineWidth',handles.handMir.DefLineThick, ...
		'FaceAlpha', 0.5, 'Tag', 'L2_Scene');
	set_uictx(h, Sat, data)

% ------------------------------------------------------------------------
function push_tracks_CB(hObject, handles)
% Plot the grand tracks

	[t_start, t_stop, fname, msg] = check_input(handles);
	if (~isempty(msg)),		errordlg(msg, 'Error'),		return,		end

	if (isempty(handles.rect))	% If not called via a rectangle clicked-callback
		BB = [get(handles.handMir.axes1, 'XLim') get(handles.handMir.axes1, 'YLim')] + [-3 3 -6 6];	% Plus a pad
	else
		BB = handles.rect;
	end
	if (BB(2) - BB(1) > 359),	BB = [];	end		% If global image, bo BoundingBox
	tracks = orbits(fname, t_start, t_stop, 0.5, BB);
	if (isempty(tracks))
		warndlg('No track data inside this region.', 'Warning'),	return
	end

	% --------If so requested, retain only the DAY or NIGHT parts of the orbits --------------
	str = get(handles.popup_satellite, 'Str');		val = get(handles.popup_satellite, 'Val');
	if (val > 1)
		if (strfind(str{val}, '(day)'))
			ind = (tracks.date(:,4) < 7 | tracks.date(:,4) > 19);
		else
			ind = (tracks.date(:,4) > 7 | tracks.date(:,4) < 19);
		end
		tracks.date(ind,:) = [];	tracks.xyz(ind,:) = [];
	end
	% ----------------------------------------------------------------------------------------

	h = line('XData',tracks.xyz(:,1), 'YData',tracks.xyz(:,2), ...
		'parent',handles.handMir.axes1, 'Color',handles.handMir.DefLineColor,'LineWidth',handles.handMir.DefLineThick);
	setappdata(h, 'ZData', tracks.xyz(:,3))		% We can't just put 'ZData' in the plot because it f... the patches. Damn TMW BUGS
	draw_funs(h,'line_uicontext')

% ------------------------------------------------------------------------
function good = check_TLE(fname)
% Check that the given TLE file exists and looks good. For the second point
% we only that it has 2 or 3 non empty lines starting by [0] 1 2 but don't
% care about its real contents. That's a job for orbits.m

	fid = fopen(fname,'r');
	if (fid == -1),		good = false;	return,		end		% It doesn't even exists

	todos = fread(fid,'*char');		fclose(fid);
	todos = strread(todos','%s','delimiter','\n');		% Cell array
	n_lines = numel(todos);

	if ((n_lines < 2) || (n_lines > 4)),	good = false;	return,		end		% Empty or shity TLE file
	if ((n_lines == 2) && (todos{1}(1) == '1') && (todos{2}(1) == '2'))
		good = true;
	elseif ((n_lines >= 3) && (todos{1}(1) == '#' || todos{1}(1) == '0') && (todos{1}(2) == '1') && (todos{3}(1) == '2'))
		good = true;
	else
		good = false;
	end

% ------------------------------------------------------------------------
function [t_start, t_stop, fname, msg] = check_input(handles)
% Check that entered parameters are reasonable. Notice that TLE
% validity was already checked by check_TLE()
	msg = '';
	t_start = get(handles.edit_dateStart, 'Str');
	t_stop  = get(handles.edit_dateStop, 'Str');
	fname   = get(handles.edit_TLE, 'Str');
	
	if (isempty(fname))
		msg = 'Error: TLE file name cannot be empty';	return
	end
	if (datenum(t_stop) <= datenum(t_start))
		msg = 'Error: Starting date is later or equal to end date.';
	end

% ------------------------------------------------------------------------
function set_uictx(h, Sat, data)
% Set the UIcontext of the scene patches
% DATA = [YYYY MM DD hh mm ss]
	handles = guidata(h);
	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand)
	uimenu(cmenuHand, 'Label', 'Delete this patch', 'Call', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Delete all patches', 'Call', {@del_patches, h});
	dy = doy(data(1), data(2), data(3));	
	L2_name_SST = sprintf('%s%d%.3d%.2d%.2d00.L2_LAC_SST.nc', Sat, data(1), dy, data(4), data(5));
	L2_name_OC  = sprintf('%s%d%.3d%.2d%.2d00.L2_LAC_OC.nc', Sat, data(1), dy, data(4), data(5));
 	url = ['http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/' L2_name_SST];
	uimenu(cmenuHand, 'Label', L2_name_SST, 'Call', {@display_url,url}, 'Sep','on');
 	url = ['http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/' L2_name_OC];
	uimenu(cmenuHand, 'Label', L2_name_OC, 'Call', {@display_url,url}, 'Sep','on');
% 	dest_fiche = [handles.path_tmp L2_name];
% 	if (ispc)
% 		dos(['wget "' url '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche ' &']);
% 	else
% 		unix(['wget "' url '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche ' &']);
% 	end

% ------------------------------------------------------------------------
function del_patches(obj, evt, h)
	handles = guidata(h);
	delete(findobj(handles.axes1, 'Type', 'patch', 'Tag', 'L2_Scene'))

% ------------------------------------------------------------------------
function display_url(obj, evt, url)
% Need to have output from inputdlg otherwise the compiler screwes.
	resp = inputdlg({'Copy-paste this address into your browser to dowload it'}, 'File web address',[1 120], {url});

% ------- TO COMPILE INCLUDE ORBITS.M HERE -------------------------------

% ------- END OF INCLUDED CODE -------------------------------------------

% --- Executes on key press over figure1 with no controls selected.%
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(hObject);
	end

% ------------------------------------------------------------------------
function h1 = sat_orbits_LayoutFcn(h1)

set(h1, 'Position',[520 490 320 300],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Satellite orbits',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[12 49 297 141], 'Style','frame');

uicontrol('Parent',h1, 'Position',[20 273 221 15],...
'HorizontalAlignment','left',...
'String','Date Start UTC (e.g. 01-Apr-2016 00:30:00)',...
'Style','text',...
'Tag','text_dateStart');

uicontrol('Parent',h1, 'Position',[20 250 221 21],...
'BackgroundColor',[1 1 1],...
'Call',@sat_orbits_uiCB,...
'String','',...
'Style','edit',...
'Tag','edit_dateStart');

uicontrol('Parent',h1, 'Position',[240 249 61 23],...
'Call',@sat_orbits_uiCB,...
'String','Calendar',...
'TooltipString','Pick a date using the calendar tool',...
'Tag','push_calendarStart');

uicontrol('Parent',h1, 'Position',[20 223 81 15],...
'HorizontalAlignment','left',...
'String','Date Stop UTC',...
'Style','text',...
'Tag','text_dateStop');

uicontrol('Parent',h1, 'Position',[20 200 221 21],...
'BackgroundColor',[1 1 1],...
'Call',@sat_orbits_uiCB,...
'String','',...
'Style','edit',...
'Tag','edit_dateStop');

uicontrol('Parent',h1, 'Position',[240 199 61 23],...
'Call',@sat_orbits_uiCB,...
'String','Calendar',...
'TooltipString','Pick a date using the calendar tool',...
'Tag','push_calendarStop');

uicontrol('Parent',h1, 'Position',[20 110 221 21],...
'BackgroundColor',[1 1 1],...
'Call',@sat_orbits_uiCB,...
'String','',...
'Style','edit',...
'TooltipString','Full name of a TLE file',...
'Tag','edit_TLE');

uicontrol('Parent',h1, 'Position',[240 109 61 23],...
'Call',@sat_orbits_uiCB,...
'String','Select TLE',...
'TooltipString','Browse for the wished TLE file',...
'Tag','push_getTLE');

uicontrol('Parent',h1, 'Position',[140 92 31 15],...
'FontSize',10,...
'FontWeight','bold',...
'String','AND',...
'Style','text',...
'Tag','text1');

uicontrol('Parent',h1, 'Position',[20 159 281 23],...
'Call',@sat_orbits_uiCB,...
'String','Open Space-Track.org for TLE',...
'TooltipString','Tell your browser to open https://www.space-track.org/#/tle',...
'Tag','push_callSpaceTrack');

uicontrol('Parent',h1, 'Position',[140 134 31 15],...
'FontSize',10,...
'FontWeight','bold',...
'String','OR',...
'Style','text',...
'Tag','text2');

uicontrol('Parent',h1, 'Position',[20 60 111 22],...
'BackgroundColor',[1 1 1],...
'Call',@sat_orbits_uiCB,...
'String',{'Select satellite'; 'AQUA (day)'; 'AQUA (night)'; 'TERRA (day)'; 'TERRA (night)'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_satellite');

uicontrol('Parent',h1, 'Position',[130 60 171 23],...
'Call',@sat_orbits_uiCB,...
'String','Get it from Space-Track.org',...
'TooltipString','Get TLE from www.space-track.org (SLOW and may fail)' ,...
'Tag','push_gotoSpaceTrack');

uicontrol('Parent',h1, 'Position',[250 10 61 21],...
'Call',@sat_orbits_uiCB,...
'String','Plot tracks',...
'TooltipString','Plot the satellite tracks on the Mirone fig.',...
'Tag','push_tracks');

uicontrol('Parent',h1, 'Position',[90 10 143 21],...
'Call',@sat_orbits_uiCB,...
'String','Plot tracks and scene tiles',...
'TooltipString','Plot the the tracks plus the scene tiles',...
'Vis', 'off', ...
'Tag','push_tracksTiles');

uicontrol('Parent',h1, 'Position',[10 10 61 21],...
'Call',@sat_orbits_uiCB,...
'String','See in GE',...
'TooltipString','See the tracks in GoogleEarth',...
'Vis','off', ...
'Tag','push_seeGE');

function sat_orbits_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'], hObject, guidata(hObject));

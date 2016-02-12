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

% $Id: sat_orbits.m 7792 2016-02-12 19:02:19Z j $

	hObject = figure('Vis','off');
	sat_orbits_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right');

	if (nargin == 0),	varargin{1} = [];	end		% So por agora
	handles.handMir = varargin{1};
	handles.hRect = [];
	if (nargin == 2)
		handles.hRect = varargin{2};
	end

	% By default select current day
	t = date;
	set(handles.edit_dateStart, 'Str', [t ' 00:00:00'])
	set(handles.edit_dateStop , 'Str', [t ' 23:59:59'])

	% See if we have traces in preferences of the last used TLE file
	s = load([handles.handMir.path_data 'mirone_pref.mat'], 'lastTLE');	% Warning: This is NEVER empty ... even when it is
	if (~isempty(fields(s)) && exist(s.lastTLE,'file'))		% Yes, we found one
		set(handles.edit_TLE, 'Str', s.lastTLE)
	end

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) ------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% ------------------------------------------------------------------------
function edit_dateStart_CB(hObject, handles)


% ------------------------------------------------------------------------
function push_calendarStart_CB(hObject, handles)
% ...
	new_date = uisetdate;
	set(handles.edit_dateStart, 'Str', new_date)

% ------------------------------------------------------------------------
function edit_dateStop_CB(hObject, handles)


% ------------------------------------------------------------------------
function push_calendarStop_CB(hObject, handles)
	new_date = uisetdate;
	set(handles.edit_dateStop, 'Str', new_date)

% ------------------------------------------------------------------------
function edit_TLE_CB(hObject, handles)
% Manually entered file or via push_getTLE_CB(). Also saves the file name in preferences

	lastTLE = get(hObject, 'Str');
	if (check_TLE(lastTLE))			% Check TLE looks good/exists
		set(hObject, 'Str', '')
		errordlg('Non existing or badly formated TLE file. Ignoring it.','Error')
		return
	end

	% Ok, since we got here we will accept the TLE as good and save its name in preferences
	version7 = version;
	V7 = (sscanf(version7(1),'%f') > 6);
	if (~V7)		% R <= 13. That's all I need to know for now.
		save([handles.handMir.path_data 'mirone_pref.mat'], 'lastTLE', '-append', '-v6')
	else
		save([handles.handMir.path_data 'mirone_pref.mat'], 'lastTLE', '-append')
	end

% ------------------------------------------------------------------------
function push_getTLE_CB(hObject, handles)
% ...
	str1 = {'*.tle;*.TLE;*.txt;*.TXT;*.dat;*.DAT;', 'TLE files (*.tle,*.txt,*.dat)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select TLE File','get');
	if isequal(FileName,0),		return,		end
	set(handles.edit_TLE, 'Str', [PathName FileName])

% ------------------------------------------------------------------------
function push_callSpaceTrack_CB(hObject, handles)


% ------------------------------------------------------------------------
function popup_satellite_CB(hObject, handles)


% ------------------------------------------------------------------------
function push_gotoSpaceTrack_CB(hObject, handles)


% ------------------------------------------------------------------------
function push_tracks_CB(hObject, handles)
% ...
	[t_start, t_stop, fname, msg] = check_input(handles);
	if (~isempty(msg))
		errordlg(msg, 'Error'),		return
	end

	%profile on
	tic
	tracks = orbits(fname, t_start, t_stop);
	toc
	%profile viewer
	h = line('XData',tracks.xyz(:,1), 'YData',tracks.xyz(:,2), 'ZData',tracks.xyz(:,3), ...
		'parent',handles.handMir.axes1);
	draw_funs(h,'line_uicontext')

% ------------------------------------------------------------------------
function push_tracksTiles_CB(hObject, handles)


% ------------------------------------------------------------------------
function push_seeGE_CB(hObject, handles)


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
'Call',@sat_orbits_uiCB,...
'BackgroundColor',[1 1 1],...
'String',{'Select satellite'; 'AQUA'; 'TERRA'},...
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
'Tag','push_tracksTiles');

uicontrol('Parent',h1, 'Position',[10 10 61 21],...
'Call',@sat_orbits_uiCB,...
'String','See in GE',...
'TooltipString','See the tracks in GoogleEarth',...
'Tag','push_seeGE');

function sat_orbits_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'], hObject, guidata(hObject));

function varargout = earth_tides(varargin)
% Helper figure to plot Earth tides using the 'earthtides' MEX derived from the 'solid.f' program

%	Copyright (c) 2004-2018 by J. Luis
%
%             DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%                     Version 2, December 2004
% 
%  Everyone is permitted to copy and distribute verbatim or modified
%  copies of this license document, and changing it is allowed as long
%  as the name is changed.
% 
%             DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%    TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
% 
%   0. You just DO WHAT THE FUCK YOU WANT TO.
%
% BUT NOTE: This license does not apply to the included argofiles() and argodata()
%			functions that have their own author and license (BSD)
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id: earth_tides.m 11312 2018-05-31 01:13:15Z j $

	if (nargin > 1 && ischar(varargin{1}))
		gui_CB = str2func(varargin{1});
		[varargout{1:nargout}] = feval(gui_CB,varargin{2:end});
	else
		h = earth_tides_OF(varargin{:});
		if (nargout),	varargout{1} = h;   end
	end

% ---------------------------------------------------------------------------------
function hObject = earth_tides_OF(varargin)
	hObject = earth_tides_LayoutFcn;
	handles = guihandles(hObject);
	move2side(hObject,'right')

	d1 = datestr(now);		d1(end-1:end) = '00';	% put seconds to 0
	d2 = datestr(now+10);	d2(end-1:end) = '00';	% put seconds to 0
	set(handles.edit_date_start, 'Str', d1)
	set(handles.edit_date_end, 'Str', d2)
	
	if (nargin)
		handles.hMirFig = varargin{1};
	else
		handles.hMirFig = [];
		set(handles.push_clickPT, 'Visible','off')
	end

	set(hObject,'Vis','on');
	guidata(hObject, handles);
	
	if (nargin > 1),	external_drive(handles, 'earth_tides', varargin{2:end}),	end

% -------------------------------------------------------------------------
function edit_date_start_CB(hObject, handles)
	if (datenum(get(hObject, 'Str') < 723181))
		warndlg('Dates older than 1-Jan-1980 have loss precision','Warning')
	end

% -------------------------------------------------------------------------
function edit_date_end_CB(hObject, handles)
	d1 = datenum(get(handles.edit_date_start, 'Str'));
	d2 = datenum(get(hObject, 'Str'));
	if (d2 <= d1)
		errordlg('End date older than start date', 'Error')
	end

% -------------------------------------------------------------------------
function edit_lon_CB(hObject, handles)
	x = str2double(get(hObject, 'Str'));
	if (isnan(x) || x < -180 || x > 360),	set(hObject, 'Str', -8),	end

% -------------------------------------------------------------------------
function edit_lat_CB(hObject, handles)
	x = str2double(get(hObject, 'Str'));
	if (isnan(x) || x < -90 || x > 90),	set(hObject, 'Str', 37),		end

% -------------------------------------------------------------------------
function radio_time_series_CB(hObject, handles)
	if (~get(hObject,'Val')),		set(hObject,'Val',1),	return,		end
	set(handles.radio_grids,'Val',0)
	set(handles.push_clickPT, 'Enable', 'on')

% -------------------------------------------------------------------------
function radio_grids_CB(hObject, handles)
	if (~get(hObject,'Val')),		set(hObject,'Val',1),	return,		end
	set(handles.radio_time_series,'Val',0)
	set(handles.push_clickPT, 'Enable', 'off')

% -------------------------------------------------------------------------
function push_clickPT_CB(hObject, handles)
% Plot a symbol in the Mirone fig where compute the tides
	if (~ishandle(handles.hMirFig))
		errordlg('You killed it (the Mir window). Bye Bye','Error'),	return
	end
	handMir = guidata(handles.hMirFig);
	if (handMir.no_file || ~handMir.geog)
		errordlg('No file loaded or it''s not a Geographical coordinates one. Bye Bye','Error'),	return
	end
	figure(handles.hMirFig)
	h = mirone('Draw_CB', handMir, 'Symbol', 'o');
	set(h, 'Tag', 'ETide')

% -------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...
	do_up = false;		do_east = false;	do_north = false;
	date1 = datevec(get(handles.edit_date_start, 'Str'));
	hdr = [-180 180 -90 90 0.5 0.5];

	if (get(handles.check_up,    'Val')),	do_up = true;		end
	if (get(handles.check_east,  'Val')),	do_east = true;		end
	if (get(handles.check_north, 'Val')),	do_north = true;	end
	ind_comp = [do_north do_east do_up];
	if (sum(ind_comp) == 0)					% If nothing selected, default to vertical
		do_up = true;
		ind_comp(3) = true;
	end

	if (get(handles.radio_time_series, 'Val'))
		opt_TG = '-T';
		date2 = datevec(get(handles.edit_date_end, 'Str'));
		if (datenum(date2) <= datenum(date1))
			errordlg('End date older than start date', 'Error'),	return
		end
		handMir = guidata(handles.hMirFig);
		hPt = findobj(handMir.axes1, 'Tag', 'ETide');
		if (~isempty(hPt))
			lon = get(hPt(1), 'XData');			lat = get(hPt(1), 'YData');
		else
			lon = str2double(get(handles.edit_lon, 'Str'));
			lat = str2double(get(handles.edit_lat, 'Str'));
		end
		hdr(1:2) = [lon lat];		% For time series only this two elements are read in MEX
	else
		opt_TG = '-G';
		if (do_up),		opt_TG = [opt_TG 'u'];	end
		if (do_east),	opt_TG = [opt_TG 'e'];	end
		if (do_north),	opt_TG = [opt_TG 'n'];	end
		tmp.X = -180:0.5:180;		tmp.Y = -90:0.5:90;		tmp.geog = 1;
		head = [tmp.X(1) tmp.X(end) -90 90 0 0 0 0.5 0.5];		% Header vector les the z_min|max
		tmp.head = head;
		names = {'Earth tide North component' 'Earth tide East component' 'Earth tide Vertical component'};
		date2 = [];			% Need to exist as a var
	end

	O = earthtide(date1, date2, hdr, opt_TG);

	if (isa(O, 'struct'))
		if (~isempty(O.north))
			zzz = grdutils(O.north, '-L');	tmp.head(5:6) = double(zzz(1:2));	tmp.name = names{1};
			mirone(O.north, tmp)
		end
		if (~isempty(O.east))
			zzz = grdutils(O.east, '-L');	tmp.head(5:6) = double(zzz(1:2));	tmp.name = names{2};
			mirone(O.east, tmp)
		end
		if (~isempty(O.up))
			zzz = grdutils(O.up, '-L');		tmp.head(5:6) = double(zzz(1:2));	tmp.name = names{1};
			mirone(O.up, tmp)
		end
	else
		if (opt_TG(2) == 'G')
			zzz = grdutils(O, '-L');		tmp.head(5:6) = double(zzz(1:2));	tmp.name = names{ind_comp};
			hFig = mirone(O, tmp);
			datasets_funs('CoastLines',guidata(hFig),'l')
		else
			ecran(O(:,1), O(:,(find(ind_comp)+1)), 'Earth tide')	% Plot the 1 to 3 components in a single Fig.
		end
	end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, evt)
% Check for "escape"
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles = guidata(hObject);
		delete(handles.figure1)
	end

% --- Creates and returns a handle to the GUI figure. 
function h1 = earth_tides_LayoutFcn()
	h1 = figure('Units','pixels', 'Position',[748 979 340 161],...
	'Color',get(0,'factoryUicontrolBackgroundColor'),...
	'KeyPressFcn',@figure1_KeyPressFcn,...
	'MenuBar','none',...
	'Name','earth_tides',...
	'NumberTitle','off',...
	'Resize','off',...
	'HandleVisibility','callback',...
	'Vis','off',...
	'Tag','figure1');

	uicontrol('Parent',h1, 'Position',[60 131 135 21],...
	'Style','edit',...
	'BackgroundColor',[1 1 1],...
	'Callback',@earthtides_uiCB,...
	'TooltipString','Start date',...
	'Tag','edit_date_start');

	uicontrol('Parent',h1, 'Position',[60 101 135 21],...
	'Style','edit',...
	'BackgroundColor',[1 1 1],...
	'Callback',@earthtides_uiCB,...
	'TooltipString','End date',...
	'Tag','edit_date_end');

	uicontrol('Parent',h1, 'Position',[231 131 100 21],...
	'String','-8',...
	'Style','edit',...
	'BackgroundColor',[1 1 1],...
	'Callback',@earthtides_uiCB,...
	'TooltipString','Longitude in decimal degrees',...
	'Tag','edit_lon');

	uicontrol('Parent',h1, 'Position',[231 101 100 21],...
	'String','37',...
	'Style','edit',...
	'BackgroundColor',[1 1 1],...
	'Callback',@earthtides_uiCB,...
	'TooltipString','Latitude in decimal degrees',...
	'Tag','edit_lat');

	uicontrol('Parent',h1, 'Position',[10 64.5 75.5 17.5],...
	'String','Time series',...
	'Style','radiobutton',...
	'Value',1,...
	'Callback',@earthtides_uiCB,...
	'Tag','radio_time_series');

	uicontrol('Parent',h1, 'Position',[10 39.5 55.5 17.5],...
	'String','Grid(s)',...
	'Style','radiobutton',...
	'Callback',@earthtides_uiCB,...
	'Tag','radio_grids');

	uicontrol('Parent',h1, 'Position',[110 61.5 57 17.5],...
	'String','Vertical',...
	'Style','checkbox',...
	'Value',1,...
	'TooltipString','Compute vertical component',...
	'Tag','check_up');

	uicontrol('Parent',h1, 'Position',[110 37.5 47 17.5],...
	'String','East',...
	'Style','checkbox',...
	'TooltipString','Compute East component',...
	'Tag','check_east');

	uicontrol('Parent',h1, 'Position',[110 12.5 47 17.5],...
	'String','North',...
	'Style','checkbox',...
	'TooltipString','Compute North component',...
	'Tag','check_north');

	uicontrol('Parent',h1, 'Position',[3.5 135 55.5 14],...
	'HorizontalAlignment','right',...
	'String','Start date',...
	'Style','text');

	uicontrol('Parent',h1, 'Position',[3.5 104 55.5 14],...
	'HorizontalAlignment','right',...
	'String','End date',...
	'Style','text');

	uicontrol('Parent',h1, 'Position',[203.5 134 25.5 14],...
	'HorizontalAlignment','right',...
	'String','Lon',...
	'Style','text');

	uicontrol('Parent',h1, 'Position',[203.5 104 25.5 14],...
	'HorizontalAlignment','right',...
	'String','Lat',...
	'Style','text');

	uicontrol('Parent',h1, 'Position',[230.5 61 101.5 21],...
	'String','Click point on map',...
	'Style',get(0,'defaultuicontrolStyle'),...
	'Callback',@earthtides_uiCB,...
	'TooltipString','Get coordinates where to compute tides by clicking on the Mirone Map.',...
	'Tag','push_clickPT');

	uicontrol('Parent',h1, 'Position',[261 6.5 71.5 21],...
	'Callback',@earthtides_uiCB,...
	'String','OK',...
	'Style',get(0,'defaultuicontrolStyle'),...
	'Tag','push_OK',...
	'FontSize',9,...
	'FontWeight','bold');

function earthtides_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

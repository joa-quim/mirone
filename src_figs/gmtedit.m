function out = gmtedit(varargin)
% Revival of the ancient Sunview gmtedit.
% gmtedit lets the user edit a gmtfile by mousing. The use should be
% fairly obvious from the menues and buttons' labels. Clicking on a point
% changes its status from good (green) to bad (red), a second click turns
% it back into a good one again.
%
%   Usage:
%       GMTEDIT with no arguments opens a empty GUI from where files can be loaded.
%       GMTEDIT(varargin)
%           All varargin(s) are optional
%           varargin{i} = 'FILE' opens the .gmt file 'FILE' (FILE must include the .gmt extension)
%           varargin{i} = '-L<width>' sets the <width> windowwidth in km (default is 300 km)
%           varargin{i} = '-G' force geodetical longitudes (0-360) output [Default is -180/+180];
%                           I'm sorry, this is against GMT defaults but I hate the [0-360] range.
%       OUT, if given will contain this figure handle.
%

%	Copyright (c) 2004-2009 by J. Luis
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

f_name = '';
got_inFile = 0;
def_width_km = 200;         % Default width in km (approximatly 1/2 of a day at 10 knots)
opt_G = '-G';               % Default output to [-180/+180] range
begin = 1;                  % This means that the display starts at the begining of track
center_win = [];
is_gmt = false;
is_mgd77 = false;
vars = [];					% Container to eventual user selected fields to plot (MGD77+ only)

	if (nargin > 0)
		f_name = varargin{1};
		[PATH,FNAME,EXT] = fileparts(f_name);
		if (isempty(EXT))       % Remember that here we need the extension
			[fid, msg] = fopen([f_name '.nc'], 'r');
			if (fid < 0)		% Ok, since it is not a .nc we'll assume it is a .gmt file
				EXT = '.gmt';	f_name = [f_name EXT];
			else				% Confirm that it's a netCDF file
				ID = fread(fid,3,'*char');		ID = ID';		fclose(fid);
				if (strcmp(ID,'CDF')),			EXT = '.nc';	f_name = [varargin{1} EXT];	end
			end
		end
		if (exist(f_name,'file') == 2)
			got_inFile = 1;
			varargin(1) = [];
			if (strcmpi(EXT,'.nc')),		is_mgd77 = true;		is_gmt = false;
			else							is_gmt = true;			is_mgd77 = false;
			end
		end
		for (k = 1:numel(varargin))
            switch (varargin{k}(1:2))
				case '-L'
					def_width_km = str2double(varargin{k}(3:end));
					if (isnan(def_width_km))        % Nonsense use of -L option
						def_width_km = 200;
					end
				case '-G'
					opt_G = ' ';
				case '-P'
					str = varargin{k}(3:end);
					[tok,rem] = strtok(str,'/');
					lon = str2double(tok);
					lat = str2double(rem(2:end));
					begin = 0;
					center_win = [lon lat];
				case '-V'
					% For new MGD77+ files only. It contains a user selection of the fields to be plot (instead of GMT)
					str = varargin{k}(3:end);
					ind = strfind(str,',');
					vars{1} = str(1:ind(1)-1);
					vars{2} = str(ind(1)+1:ind(2)-1);
					vars{3} = str(ind(2)+1:end);
			end
		end
	end

	sc_size = get(0,'ScreenSize');
	fig_height = fix(sc_size(4)*.9);
	fp = [1 sc_size(4)-fig_height sc_size(3) fig_height];
	marg_l = 60;                            % Left margin
	marg_r = 15;                            % Right margin
	marg_tb = 10;                           % Top & Bottom margins
	marg_ax = 15;                           % Margin between axes
	ax_height = fix(fig_height * .295);
	ax_width = sc_size(3) - marg_l - marg_r;

	hf = figure('name','gmtedit','resize','off','numbertitle','off', 'visible','off', ...
        'position',fp, 'DoubleBuffer','on', 'Tag','figure1', 'closerequestfcn','delete(gcbf)');

	% Use system color scheme for figure:
	set(hf,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Apply this trick to get the icons. Let's hope that this is not version/OS dependent 
	hA = findall(hf);
	hh = findobj(hA,'TooltipString','Open File');
	openFile_img = get(hh(1),'CData');
	hh = findobj(hA,'TooltipString','Save Figure');
	saveFile_img = get(hh(1),'CData');
	hh = findobj(hA,'TooltipString','Zoom In');
	zoomIn_img = get(hh(1),'CData');
	hh = findobj(hA,'TooltipString','Zoom Out');
	zoomOut_img = get(hh(1),'CData');
	set(hf,'menubar','none')            % Set the menubar to none

	pos = [marg_l fig_height-ax_height-marg_tb ax_width ax_height];
	h_a1 = axes('Parent',hf, 'Units','pixels', 'Position',pos, 'XLim',[0 def_width_km], 'Tag','axes1');

	pos = [marg_l fig_height-2*(ax_height+marg_tb)-marg_ax ax_width ax_height];
	h_a2 = axes('Parent',hf, 'Units','pixels', 'Position',pos, 'XLim',[0 def_width_km], 'Tag','axes2');

	pos = [marg_l fig_height-3*(ax_height+marg_tb)-2*marg_ax ax_width ax_height];
	h_a3 = axes('Parent',hf,'Units','pixels', 'Position',pos, 'XLim',[0 def_width_km], 'Tag','axes3');

	handles = guihandles(hf);
	handles.opt_G = opt_G;
	handles.info = [];
	handles.h_broken = [];
	handles.begin = begin;
	handles.center_win = center_win;
	handles.is_gmt = is_gmt;
	handles.is_mgd77 = is_mgd77;
	handles.vars = vars;

	MIRONE_DIRS = getappdata(0,'MIRONE_DIRS');
	if (~isempty(MIRONE_DIRS))							% Should not be empty, but ...
		handles.home_dir = MIRONE_DIRS.home_dir;		% False info if not called from Mir root dir
	else
		handles.home_dir = cd;
	end

	if (got_inFile)
		handles.last_dir = PATH;
	else
		handles.last_dir = handles.home_dir;			% This means, last_dir is not saved between sessions
	end
	handles.work_dir = handles.last_dir;

	% Load some icons from mirone_icons.mat
	load([handles.home_dir filesep 'data' filesep 'mirone_icons.mat'],'rectang_ico','info_ico','trincha_ico');

	%SetAxesNumericType(handles.axes1);        % Set axes uicontextmenus
	if (isempty(vars))		% Default GMT axes names
		set(get(h_a1,'YLabel'),'String','Gravity anomaly (mGal)')
		set(get(h_a2,'YLabel'),'String','Magnetic anomaly (nT)')
		set(get(h_a3,'YLabel'),'String','Bathymetry (m)')
	else
		set(get(h_a1,'YLabel'),'String',vars{1})
		set(get(h_a2,'YLabel'),'String',vars{2})
		set(get(h_a3,'YLabel'),'String',vars{3})
	end

	handles.def_width_km = def_width_km;
	handles.max_x_data = 10000;
	scroll_plots(def_width_km,[0 10000])     % The [0 10000] will be reset inside scroll_plots to a more appropriate val
	movegui(hf,'north')

h_toolbar = uitoolbar('parent',hf,'Clipping', 'on', 'BusyAction','queue','HandleVisibility','on',...
   'Interruptible','on','Tag','FigureToolBar','Visible','on');
uipushtool('parent',h_toolbar,'Click',{@import_clickedCB,f_name},'Tag','import',...
   'cdata',openFile_img,'Tooltip','Open gmt file');
uipushtool('parent',h_toolbar,'Click',@save_clickedCB,'Tag','save', 'cdata',saveFile_img,'Tooltip','Save gmt file');
% uitoggletool('parent',h_toolbar,'Click',@zoom_clickedcallback,'Tag','zoom',...
%    'cdata',zoom_img,'TooltipString','Zoom');
uipushtool('parent',h_toolbar,'Click',@info_clickedCB,'Tag','info','cdata',info_ico, 'Tooltip','Cruise Info');
uipushtool('parent',h_toolbar,'Click',@rectang_clickedCB,'Tag','rectang','cdata',rectang_ico,...
   'Tooltip','Rectangular region','Sep','on');
uipushtool('parent',h_toolbar,'Click',@rectangMove_clickedCB,'cdata',rectang_ico,'Tooltip','Select for moving');
uipushtool('parent',h_toolbar,'Click',{@changeScale_clickedCB,'inc'}, 'cdata',zoomIn_img,'Tooltip','Increase scale','Sep','on');
uipushtool('parent',h_toolbar,'Click',{@changeScale_clickedCB,'dec'}, 'cdata',zoomOut_img,'Tooltip','Decrease scale');
uipushtool('parent',h_toolbar,'Click',@outliers_clickedCB, 'cdata',trincha_ico,'Tooltip','Outliers detector','Sep','on');

	% Create empty lines just for the purpose of having their handles
	handles.h_gl = line('XData',[],'YData',[],'Color','k','Parent',h_a1);
	handles.h_gm = line('XData',[],'YData',[],'LineStyle','none','Marker','s', ...
		'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',4,'Parent',h_a1);
	handles.h_ml = line('XData',[],'YData',[],'Color','k','Parent',h_a2);
	handles.h_mm = line('XData',[],'YData',[],'LineStyle','none','Marker','s', ...
		'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',4,'Parent',h_a2);
	handles.h_tl = line('XData',[],'YData',[],'Color','k','Parent',h_a3);
	handles.h_tm = line('XData',[],'YData',[],'LineStyle','none','Marker','s', ...
		'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',4,'Parent',h_a3);

	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set([h_a1 h_a2 h_a3], 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Overly interpolation', 'Call',@interp_on_grid);

	guidata(hf, handles);

	% Now that we have a figure and its handles, we can open the file if it was requested on input
	if (got_inFile),	import_clickedCB(hf,[],f_name),		end

	% Add or remove red Markers
	set(hf,'WindowButtonDownFcn',@add_MarkColor,'Visible','on')

	% Choose default command line output for gmtedit
	if (nargout == 1),	out = hf;   end

% --------------------------------------------------------------------
function import_clickedCB(hObject, eventdata, opt)
	handles = guidata(hObject);     % get handles
	if (isempty(opt))
		[FileName,PathName] = put_or_get_file(handles,{'*.gmt;*.GMT;*.nc', 'gmt files (*.gmt,*.GMT,*.nc)'},'Select gmt File','get');
		if isequal(FileName,0),		return,		end
		f_name = [PathName FileName];
	else
		f_name = opt;
	end

	[PATH,FNAME,EXT] = fileparts(f_name);
	if (strcmpi(EXT,'.nc'))
		handles.f_name = [PATH filesep FNAME '.nc'];		handles.is_mgd77 = true;		handles.is_gmt = false;
		[handles, track] = read_mgd77_plus(handles, f_name);
	else
		handles.f_name = [FNAME '.gmt'];	handles.is_gmt = true;			handles.is_mgd77 = false;
		track = gmtlist_m([PATH filesep FNAME], '-Fsxygmtd', handles.opt_G);
	end
	set(handles.figure1,'Name',['gmtedit  ' FNAME EXT])	

	% Save those for use when saving into a new file
	if (handles.is_gmt),	handles.time = track.time;		end
	handles.lon = track.longitude;
	handles.lat = track.latitude;
	handles.year = track.year;
	handles.info = track.info;
	if (length(track.agency) ~= 10)         % Ensures that agency is exactly 10 chars
		agency = '          ';              % 10 blanks
		len = min(length(track.agency),10);
		agency(1:len) = track.agency(1:len);
		handles.agency = agency;
	else
		handles.agency = track.agency;
	end

	% Search for de-activated (red marked) points handles - They may exist if another file was already loaded
	h_gn = findobj(handles.figure1,'Type','Line','tag','GravNull');
	h_mn = findobj(handles.figure1,'Type','Line','tag','MagNull');
	h_tn = findobj(handles.figure1,'Type','Line','tag','TopNull');
	if (~isempty(h_gn)),	set(h_gn,'XData',[],'YData',[]);	end
	if (~isempty(h_mn)),	set(h_mn,'XData',[],'YData',[]);	end
	if (~isempty(h_tn)),	set(h_tn,'XData',[],'YData',[]);	end

	% See if any "broken line" was left behind - They may exist if another file was already loaded and processed
	if (~isempty(handles.h_broken))         % If we have broken lines, delete them
		for (k=1:length(handles.h_broken))
			set(handles.h_broken(k),'Xdata',[],'YData',[])
		end
		rmfield(handles,'h_broken');
		handles.h_broken = [];
	end

	if ( numel(track.gravity) > 1 && ~all(isnan(track.gravity)) )
		set(handles.h_gm,'XData',track.distance,'YData',track.gravity, 'Tag','orig_grav')
	else
		set(handles.h_gm,'XData',[],'YData',[])    
	end
	if ( numel(track.magnetics) > 1 && ~all(isnan(track.magnetics)) )
		set(handles.h_mm,'XData',track.distance,'YData',track.magnetics, 'Tag','orig_mag')
	else
		set(handles.h_mm,'XData',[],'YData',[])    
	end
	if ( numel(track.topography) > 1 && ~all(isnan(track.topography)) )
		set(handles.h_tm,'XData',track.distance,'YData',track.topography, 'Tag','orig_topo')
	else
		set(handles.h_tm,'XData',[],'YData',[])    
	end

	% Update the slider Max propertie
	hs = findobj(handles.figure1,'style','slider');
	handles.max_x_data = track.distance(end);
	max_s = handles.max_x_data - handles.def_width_km;
	if (max_s < 0),		max_s = handles.max_x_data;		end			% I already had one case like this. A very short track
	val = track.distance(1);
	%set(hs,'Max',handles.max_x_data-handles.def_width_km,'value',track.distance(1))

	if (~handles.begin)         % Start the display at a user selected coordinate
		x = handles.lon - handles.center_win(1);    [zz,id1] = min(abs(x));     clear x;
		y = handles.lat - handles.center_win(2);    [zz,id2] = min(abs(y));     clear y;
		% id1 and id2 are not forcedly equal. Find out the "best"
		r1 = sqrt((handles.lon(id1) - handles.center_win(1))^2 + (handles.lat(id1) - handles.center_win(2))^2);
		r2 = sqrt((handles.lon(id2) - handles.center_win(1))^2 + (handles.lat(id2) - handles.center_win(2))^2);
		id = id1;
		if (r2 < r1),	id = id2;	end
		x_lim = track.distance(id) + [-handles.def_width_km/2 handles.def_width_km/2];
		set(findall(handles.figure1,'Type','axes'),'xlim',x_lim)
		val0 = track.distance(id)-handles.def_width_km;
		if (val0 > track.distance(1)),		val = val0;		end
	end

	% Update the slider properties
	r = 0.90 * handles.def_width_km / handles.max_x_data;		% Give a 10% overlap between display area when click on slider arrow
	set(hs,'Max',max_s,'value',val,'SliderStep',[r r*10])
	guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function interp_on_grid(obj, event)
	handles = guidata(obj);
	hAx = get(handles.figure1, 'CurrentAxes');		% Find which axes are we working on
	gmtedit_track(handles.lon, handles.lat, handles.year, hAx, handles.h_tm)

% --------------------------------------------------------------------
function [handles, track] = read_mgd77_plus(handles, fname)
	s = nc_funs('info',fname);

	% ------------------ OK, Get numerics now -----------------------------------
	tempo = [];
	track.longitude = double(nc_funs('varget', fname, 'lon'));
	track.latitude  = double(nc_funs('varget', fname, 'lat'));
	if (isempty(handles.vars))
		track.magnetics = nc_funs('varget', fname, 'mtf1');
		track.gravity = nc_funs('varget', fname, 'faa');
		track.topography = nc_funs('varget', fname, 'depth');
	else
		try
			ind_v = strmatch('vel', handles.vars);		% Was velocity asked for?
			if ( any(ind_v) )		% A bit boring below because 'vel' is not a MGD77+ Variable. So we need to avoid error
				tempo = nc_funs('varget', fname, 'time');
				if (ind_v(1) == 1)
					track.magnetics  = nc_funs('varget', fname, handles.vars{2});
					track.topography = nc_funs('varget', fname, handles.vars{3});
				elseif (ind_v(1) == 2)
					track.gravity  = nc_funs('varget', fname, handles.vars{1});
					track.topography = nc_funs('varget', fname, handles.vars{3});
				else
					track.gravity  = nc_funs('varget', fname, handles.vars{1});
					track.magnetics = nc_funs('varget', fname, handles.vars{2});
				end
			else
				track.gravity = nc_funs('varget', fname, handles.vars{1});
				track.magnetics  = nc_funs('varget', fname, handles.vars{2});
				track.topography = nc_funs('varget', fname, handles.vars{3});
			end
		catch
			errordlg('GMTEDIT: At least one custom selected field name does ot exist in dataset','Error')
			error('GMTEDIT:read_mgd77_plus', lasterr)
		end
	end
	
	if ( (any(strmatch('lon', handles.vars)) | any(strmatch('lat', handles.vars))) & isempty(tempo) )
		track.distance = 1:numel(track.gravity);		% We plot lon/lat against the record number
	else												% Need to compute the acumulated distance along profile
		D2R = pi / 180;		KMPRDEG = 111.1949;
		co = cos(track.latitude * D2R);
		dx = [0; diff(track.longitude)];		dy = [0; diff(track.latitude)] .* co;
		track.distance = cumsum(sqrt(dx .^2 + dy .^2) * KMPRDEG);
		if (~isempty(tempo))			% Compute velocity and find where to store it
			vel = [4; diff(track.distance) ./ diff(tempo) * 1e3];		% Speed in m/s
			if (ind_v(1) == 1),			track.gravity = vel;		% Remember that grav, mag, topo are default field names 
			elseif (ind_v(1) == 2),		track.magnetics = vel;
			else						track.topography = vel;
			end
			track.distance = 1:numel(track.gravity);	% We plot against the record number
		end
	end

	% -------------- Get the nodata-values ---------------------------------------
	ind = strcmp({s.Dataset.Name}, 'mtf1');
	id = strcmp({s.Dataset(ind).Attribute.Name}, 'missing_value');
	handles.magNoValue = s.Dataset(ind).Attribute(id).Value;
	id = strcmp({s.Dataset(ind).Attribute.Name}, 'scale_factor');
	handles.magScaleF = s.Dataset(ind).Attribute(id).Value;

	ind = strcmp({s.Dataset.Name}, 'faa');
	id = strcmp({s.Dataset(ind).Attribute.Name}, 'missing_value');
	handles.gravNoValue = s.Dataset(ind).Attribute(id).Value;
	id = strcmp({s.Dataset(ind).Attribute.Name}, 'scale_factor');
	handles.gravScaleF = s.Dataset(ind).Attribute(id).Value;
	
	ind = strcmp({s.Dataset.Name}, 'depth');
	id = strcmp({s.Dataset(ind).Attribute.Name}, 'missing_value');
	handles.topoNoValue = s.Dataset(ind).Attribute(id).Value;
	id = strcmp({s.Dataset(ind).Attribute.Name}, 'scale_factor');
	handles.topoScaleF = s.Dataset(ind).Attribute(id).Value;
	% -----------------------------------------------------------------------------

	ind = strcmp({s.Attribute.Name}, 'Survey_Departure_Year');
	if (~any(ind))
		ind = strcmp({s.Attribute.Name}, 'Survey_Departure_Year_REVISED');
		track.year = str2double(s.Attribute(ind).Value);
		ind = strcmp({s.Attribute.Name}, 'Survey_Departure_Month_REVISED');
		track.month = str2double(s.Attribute(ind).Value);
		ind = strcmp({s.Attribute.Name}, 'Source_Institution_REVISED');
		track.agency = s.Attribute(ind).Value;
	else
		track.year = str2double(s.Attribute(ind).Value);
		ind = strcmp({s.Attribute.Name}, 'Survey_Departure_Month');
		track.month = str2double(s.Attribute(ind).Value);
		ind = strcmp({s.Attribute.Name}, 'Source_Institution');
		track.agency = s.Attribute(ind).Value;
	end
	track.info = 'Bla Bla';				% PRECISAMOS ACABAR AQUI

% --------------------------------------------------------------------
function save_clickedCB(hObject, eventdata)
	handles = guidata(hObject);     % get handles
	if (handles.is_gmt)
		NODATA(1:3) = -32000;
		[FileName,PathName] = put_or_get_file(handles, handles.f_name,'Select gmt File', 'put','.gmt');
		if isequal(FileName,0),		return,		end
		f_name = [PathName FileName];
	else
		NODATA(1) = handles.gravNoValue;	NODATA(2) = handles.magNoValue;		NODATA(3) = handles.topoNoValue;
		f_name = handles.f_name;	% Since we only update, there is no choice here.
	end

	% Search for de-activated points handles (the reds)
	h_gn = findobj(handles.figure1,'Type','Line','tag','GravNull');
	h_mn = findobj(handles.figure1,'Type','Line','tag','MagNull');
	h_tn = findobj(handles.figure1,'Type','Line','tag','TopNull');

	% And the corresponding red marker values values
	x_gn = get(h_gn,'XData');			y_gn = get(h_gn,'YData');
	x_mn = get(h_mn,'XData');			y_mn = get(h_mn,'YData');
	x_tn = get(h_tn,'XData');			y_tn = get(h_tn,'YData');

	% Get G,M,T values
	x_g = get(handles.h_gm,'XData');	y_g = get(handles.h_gm,'YData');
	x_m = get(handles.h_mm,'XData');	y_m = get(handles.h_mm,'YData');
	x_t = get(handles.h_tm,'XData');	y_t = get(handles.h_tm,'YData');

	if (~isempty(handles.h_broken))         % If we have broken lines we must join them
		x_broken = [];    y_broken = [];
		for (k=1:length(handles.h_broken))
			x_broken = [x_broken get(handles.h_broken(k),'XData')];
			y_broken = [y_broken get(handles.h_broken(k),'YData')];
		end
		[x_m,id] = sort([x_m x_broken]);
		y_m = [y_m y_broken];
		y_m = y_m(id);                      % Otherwise the y's would be out of order
	end

	set(handles.figure1,'Pointer','watch')
	x_lim = get(get(handles.figure1,'CurrentAxes'),'XLim');
	if (~isempty(x_gn))
		id_x = zeros(numel(x_gn),1);
		for (k=1:numel(x_gn))
			tmp = find((x_g - x_gn(k)) == 0);   % Find the gravity points that were marked
			id_x(k) = tmp(1);				% Old files often have repeated coords
		end
		y_g(id_x) = NODATA(1);				% Remove them
	end
	if (~isempty(x_mn))
		for (k=1:numel(x_mn))
			tmp = find((x_m - x_mn(k)) == 0);   % Find the magnetic points that were marked
			y_m(tmp) = NODATA(2);				% Remove them
		end
	end
	if (~isempty(x_tn))
		id_x = zeros(numel(x_tn),1);
		for (k=1:numel(x_tn))
			tmp = find((x_t - x_tn(k)) == 0);   % Find the topo points that were marked
			id_x(k) = tmp(1);				% Old files often have repeated coords
		end
		y_t(id_x) = NODATA(3);					% Remove them
	end

	n_rec = numel(handles.lon);
	if (handles.is_gmt && (isempty(y_gn) || isempty(y_mn) || isempty(y_tn)) )
		if (isempty(y_g)),		y_g = repmat(int16(NODATA(1)),1,n_rec);		end		% Original had not this type of data
		if (isempty(y_m)),		y_m = repmat(int16(NODATA(1)),1,n_rec);		end
		if (isempty(y_t)),		y_t = repmat(int16(NODATA(1)),1,n_rec);		end
	end

	if (~isempty(y_g))
		if (handles.is_gmt),	y_g(isnan(y_g*10)) = NODATA(1);		y_g = int16(y_g);
		elseif (~isempty(x_gn))	y_g = y_g / handles.gravScaleF;		y_g(isnan(y_g)) = NODATA(1);	y_g = int16(y_g);
		end
	end
	if (~isempty(y_m))
		if (handles.is_gmt),	y_m(isnan(y_m)) = NODATA(2);		y_m = int16(y_m);
		elseif (~isempty(x_mn))	y_m = y_m / handles.magScaleF;		y_m(isnan(y_m)) = NODATA(2);	y_m = int32(y_m);
		end
	end
	if (~isempty(y_t))
		if (handles.is_gmt),	y_t(isnan(y_t)) = NODATA(3);		y_t = int16(y_t);
		elseif (~isempty(x_tn))	y_t = y_t / handles.topoScaleF;		y_t(isnan(y_t)) = NODATA(3);	y_t = int32(y_t);
		end
	end

	if (handles.is_gmt)			% Old style .gmt files
		tempo = int32(handles.time);
		lat  = int32(handles.lat * 1e6);    % And convert back to millidegrees
		lon  = int32(handles.lon * 1e6);

		fid = fopen(f_name,'wb');
		fwrite(fid,[int32(handles.year) n_rec],'int32');
		fwrite(fid,handles.agency,'schar');
		% This is STUPIDLY slow but I didn't find any other way to do it.
		for (k=1:n_rec)
			fwrite(fid,[tempo(k) lat(k) lon(k)],'int32');
			fwrite(fid,[y_g(k) y_m(k) y_t(k)],'int16');
		end
		fclose(fid);
	else						% New mgf77+ netCDF style files
		if (~isempty(y_gn)),		nc_funs('varput', f_name, 'faa', y_g, 0);		end		% The 0 flags nc_funs for not try to scale
		if (~isempty(y_mn)),		nc_funs('varput', f_name, 'mtf1', y_m, 0);		end
		if (~isempty(y_tn)),		nc_funs('varput', f_name, 'depth', y_t, 0);		end
	end

	set(handles.figure1,'Pointer','arrow')

% --------------------------------------------------------------------
function add_MarkColor(hObject, eventdata)
% Add a red Marker over the closest (well, near closest) clicked point.
	handles = guidata(hObject);     % get handles

	button = get(handles.figure1, 'SelectionType');
	if (~strcmp(button,'normal')),		return,		end		% Accept only left-clicks

	in_grav = 0;    in_mag = 0;     in_topo = 0;
	ax = get(handles.figure1,'CurrentAxes');
	pt = get(ax, 'CurrentPoint');
	if (strcmp(get(ax,'Tag'),'axes1'))
		in_grav = 1;        opt = 'GravNull';
	elseif (strcmp(get(ax,'Tag'),'axes2'))
		in_mag = 1;         opt = 'MagNull';
	elseif (strcmp(get(ax,'Tag'),'axes3'))
		in_topo = 1;        opt = 'TopNull';
	end
	if ((in_grav + in_mag + in_topo) == 0),		return,		end		% Click was outside axes

	if (in_grav)
		hM = findobj(handles.figure1,'Type','Line','tag','GravNull');
		x = get(handles.h_gm,'XData');      y = get(handles.h_gm,'YData');
	elseif (in_mag)
		hM = findobj(handles.figure1,'Type','Line','tag','MagNull');
		x = get(handles.h_mm,'XData');      y = get(handles.h_mm,'YData');
	else
		hM = findobj(handles.figure1,'Type','Line','tag','TopNull');
		x = get(handles.h_tm,'XData');      y = get(handles.h_tm,'YData');
	end

	x_lim = get(ax,'XLim');					y_lim = get(ax,'YLim');
	dx = diff(x_lim) / 20;					% Search only betweem +/- 1/10 of x_lim
	dx = max(dx, 10);
	id = (x < (pt(1,1)-dx) | x > (pt(1,1)+dx));
	x(id) = [];					y(id) = [];	% Clear outside-2*dx points to speed up the search code
	XScale = diff(x_lim);		YScale = diff(y_lim)*6;		% The 6 factor compensates the ~6:1 horizontal/vertical axes dimension

	r = sqrt(((pt(1,1)-x)./XScale).^2+((pt(1,2)-y)./YScale).^2);
	[temp,i] = min(r);
	pt_x = x(i);				pt_y = y(i);

	xr = get(hM,'XData');		yr = get(hM,'YData');		% Red markers
	id = find(xr == pt_x);
	if (isempty(id))			% New Marker
		if (isempty(hM))		% First red Marker on this axes
			line(pt_x,pt_y,'Marker','s','MarkerFaceColor','r','MarkerSize',4,'LineStyle','none','Tag',opt);
		else
			xr = [xr pt_x];		yr = [yr pt_y];
			set(hM,'XData',xr, 'YData', yr)
		end
	else						% Marker already exists. Kill it
		xr(id) = [];			yr(id) = [];
		set(hM,'XData',xr, 'YData', yr)
	end

% --------------------------------------------------------------------------------------------------
function zoom_clickedcallback(obj,eventdata)
	if (strcmp(get(obj,'State'),'on'))
		zoom_j xon;
	else
		zoom_j off;
	end

% --------------------------------------------------------------------------------------------------
function info_clickedCB(obj,eventdata)
	handles = guidata(obj);     % get handles
	if (isempty(handles.info)),		return,		end
	data = handles.info;
	str{1} = sprintf('N_recs = %d, N_grav = %d, N_mag = %d, N_topo = %d', data(1:4));
	str{2} = ['E: = ' num2str(data(5)) '  W: = ' num2str(data(6))];
	str{3} = ['S: = ' num2str(data(7)) '  N: = ' num2str(data(8))];
	str{4} = ['Start day,month,year: = ' num2str(data(9)) '  ' num2str(data(10)) '  ' num2str(data(11))];
	str{5} = ['End   day,month,year: = ' num2str(data(12)) '  ' num2str(data(13)) '  ' num2str(data(14))];
	msgbox(str,'Cruise Info')

% --------------------------------------------------------------------------------------------------
function rectang_clickedCB(obj,eventdata)
handles = guidata(obj);     % get handles
try
	[p1,p2] = rubberbandbox;
catch		% Don't know why but uisuspend sometimes breaks
	set(handles.figure1,'Pointer','arrow');
	return
end

in_grav = 0;	in_mag = 0;		in_topo = 0;
ax = get(handles.figure1,'CurrentAxes');
if (strcmp(get(ax,'Tag'),'axes1')),			in_grav = 1;	opt = 'GravNull';
elseif (strcmp(get(ax,'Tag'),'axes2')),		in_mag = 1;		opt = 'MagNull';
elseif (strcmp(get(ax,'Tag'),'axes3')),		in_topo = 1;	opt = 'TopNull';
end

if (in_grav)
    hM = findobj(handles.figure1,'Type','Line','tag','GravNull');
    x = get(handles.h_gm,'XData');      y = get(handles.h_gm,'YData');
elseif (in_mag)
    hM = findobj(handles.figure1,'Type','Line','tag','MagNull');
    x = get(handles.h_mm,'XData');      y = get(handles.h_mm,'YData');
else
    hM = findobj(handles.figure1,'Type','Line','tag','TopNull');
    x = get(handles.h_tm,'XData');      y = get(handles.h_tm,'YData');
end

% Search for "regular" points inside the rectangle
id = find(x >= p1(1,1) & x <= p2(1,1) & y >= p1(1,2) & y <= p2(1,2));
if (isempty(id))    return;     end     % Nothing inside rect
golo_x = x(id);             golo_y = y(id);

% Now search also for eventual pre-existing red markers inside the rectangle
xr = get(hM,'XData');       yr = get(hM,'YData');
if (~isempty(xr))           % We have red markers. See if any is inside rect
    id = find(xr >= p1(1,1) & xr <= p2(1,1) & yr >= p1(1,2) & yr <= p2(1,2));
    if (~isempty(id))       % Yes we have, so remove them for not creating duplicates
        xr(id) = [];        yr(id) = [];
    end
end

% Finaly plot the points cought inside the rect
if (isempty(hM))        % First red Markers on this axes
    line(golo_x,golo_y,'Marker','s','MarkerFaceColor','r','MarkerSize',4,'LineStyle','none','Tag',opt);
else
    golo_x = [golo_x xr];    golo_y = [golo_y yr];
    set(hM,'XData',golo_x, 'YData', golo_y)
end

% --------------------------------------------------------------------------------------------------
function rectangMove_clickedCB(obj,eventdata)
handles = guidata(obj);     % get handles
try
    [p1,p2] = rubberbandbox;
catch
    return
end

if (~strcmp(get(get(handles.figure1,'CurrentAxes'),'Tag'),'axes2'))		% This option is meant only to magnetic data
    return
end

x = get(handles.h_mm,'XData');      y = get(handles.h_mm,'YData');

% Search for "regular" points inside the rectangle
id = find(x >= p1(1,1) & x <= p2(1,1) & y >= p1(1,2) & y <= p2(1,2));
if (isempty(id)),	return,		end		% Nothing inside rect

% create a new line and set it the line_uicontext in order that they may be moved
golo_x = x(id);             golo_y = y(id);

x(id) = [];         y(id) = [];         % Remove these from the main line
set(handles.h_mm,'XData',x,'YData',y)

n = length(handles.h_broken);
handles.h_broken(n+1) = line(golo_x,golo_y,'Marker','s','MarkerFaceColor','b','MarkerSize',4,'LineStyle','-');
ui_edit_polygon(handles.h_broken(n+1))      % Set edition functions

guidata(handles.figure1, handles)

% --------------------------------------------------------------------------------------------------
function changeScale_clickedCB(obj,eventdata,opt)
handles = guidata(obj);     % get handles

if (strcmp(opt,'inc'))
	handles.def_width_km = handles.def_width_km - 50;
	if (handles.def_width_km < 50),		return,		end
else            % Decrease scale
	handles.def_width_km = handles.def_width_km + 50;
end

x_lim = get(get(handles.figure1,'CurrentAxes'),'Xlim');
set([handles.axes1 handles.axes2 handles.axes3],'Xlim',x_lim(1)+[0 handles.def_width_km])

h_slider = findobj(handles.figure1,'style','slider');
cb = get(h_slider,'callback');
% Here we make use of the knowledge that the "cb" is a string of the form:
% set(findall(gcf,'Type','axes'),'xlim',get(gcbo,'value')+[0 ???])
% where '???' is the width of the currently displyed axes. And that's what we need to change
new_cb = [cb(1:59) num2str(handles.def_width_km) '])'];
set(h_slider,'callback',new_cb)

% Now update the slider 'Max' propertie
new_max = handles.max_x_data - handles.def_width_km;
set(h_slider,'Max',new_max)

val = get(h_slider,'Value');
if (val > new_max),		set(h_slider,'Value',new_max),		end

guidata(handles.figure1, handles)

% --------------------------------------------------------------------------------------------------
function outliers_clickedCB(obj,eventdata,opt)
% Detect outliers using a spline smooth technique.
	handles = guidata(obj);					% get handles
	outliersdetect(handles.figure1, handles.axes1, handles.axes2, handles.axes3, ...
		[handles.h_gm handles.h_mm handles.h_tm]);

% --------------------------------------------------------------------
function scroll_plots(width,x);
%   This function uses the idea of the scrollplotdemo from Steven Lord (http://www.mathworks.com/matlabcentral)
%   scroll_plots(width,x);
%   width: window width in x units
%   x: absicssae vector

% Generate constants for use in uicontrol initialization
pos = get(gca,'position');

% This will create a slider which is just underneath the axis
% but still leaves room for the axis labels above the slider
Newpos=[pos(1) 5 pos(3) 15];

S=['set(findall(gcf,''Type'',''axes''),''xlim'',get(gcbo,''value'')+[0 ' num2str(width) '])'];

% Creating Uicontrol with initial value of the minimum of x
uicontrol('style','slider','units','pixels','position',Newpos, ...
    'callback',S,'min',x(1),'max',x(end)-width,'value',x(1));

% --------------------------------------------------------------------------
function varargout = outliersdetect(varargin)
% Do automatic outliers detection by comparison with spline sooth version data
 
	hObject = figure('Tag','figure1','Visible','off');
	outliersdetect_LayoutFcn(hObject);
	handles = guihandles(hObject);

	handles.hCallingFig = varargin{1};
	handles.hCallingAx1 = varargin{2};
	handles.hCallingAx2 = varargin{3};
	handles.hCallingAx3 = varargin{4};
	handles.hChannel = varargin{5};
	
	handles.thresh(1) = 0.4;
	handles.thresh(2) = 4;
	handles.thresh(3) = 10;
	handles.smooth(1) = 1;
	handles.smooth(2) = 1;
	handles.smooth(3) = 1;
	warning off SPLINES:CHCKXYWP:NaNs

	set(hObject,'Visible','on');	drawnow
	
	x = get(handles.hChannel(1),'XData');		y = get(handles.hChannel(1),'YData');
	have_g = any(~isnan(y));
	if (have_g)
		[pp,p] = spl_fun('csaps',x,y);		% This is just to get csaps's p estimate
		handles.smooth(1) = p;
	else
		set(handles.radio_G, 'Enable','off')
	end
	x = get(handles.hChannel(2),'XData');		y = get(handles.hChannel(2),'YData');
	have_m = any(~isnan(y));
	if (have_m)
		[pp,p] = spl_fun('csaps',x,y);
		handles.smooth(2) = p;
	else
		set(handles.radio_M, 'Enable','off')
	end
	x = get(handles.hChannel(3),'XData');		y = get(handles.hChannel(3),'YData');
	have_t = any(~isnan(y));
	if (have_t)
		[pp,p] = spl_fun('csaps',x,y);
		handles.smooth(3) = p;
	else
		set(handles.radio_T, 'Enable','off')
	end

	% Fill the edit boxes with appropriate values (priority is Mag, than Grav and last is Topo)
	if (have_g),		handles.id_gmt = 1;		set(handles.radio_G,'Val',1)
	elseif (have_m),	handles.id_gmt = 2;		set(handles.radio_M,'Val',1)
	else				handles.id_gmt = 3;		set(handles.radio_T,'Val',1)
	end
	set(handles.edit_thresh,'String',handles.thresh(handles.id_gmt));
	set(handles.edit_SmoothParam,'String',num2str(handles.smooth(handles.id_gmt)));

	str = sprintf(['Residues greater or equal than this are outliers.\n' ...
		'Notice that we use small numbers because the spline\n' ...
		'smoothing will do only a mild smoothing, so the residues\n', ...
		'are naturally small. Unless you decrease the p parameter']);
	set(handles.edit_thresh,'Tooltip', str)

	guidata(hObject, handles);	
	if (nargout),	varargout{1} = hObject;		end

% ----------------------------------------------------------------------------
function edit_SmoothParam_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (xx < 0 || xx > 1 || isnan(xx))
		xx = 1;		set(hObject,'String',xx)
	end
	handles.smooth(handles.id_gmt) = xx;
	guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------
function edit_thresh_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (xx < 0 || isnan(xx))
		xx = 0;		set(hObject,'String',xx)
	end
	handles.thresh(handles.id_gmt) = xx;
	guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------
function radio_G_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_M handles.radio_T],'Val', 0)
	handles.id_gmt = 1;		guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function radio_M_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_G handles.radio_T],'Val', 0)
	handles.id_gmt = 2;		guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function radio_T_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_G handles.radio_M],'Val', 0)
	handles.id_gmt = 3;		guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function push_Apply_Callback(hObject, eventdata, handles)
% Detect outliers using a spline smooth technique. Note that the threshold is normaly low
% because the smoothing is very mild and therefore the residues are small.

	% ------------- Get the working channel ---------------------------
	id_gmt = handles.id_gmt;
	hChannel = handles.hChannel(id_gmt);
	markers_tag = {'GravNull' 'MagNull' 'TopNull'};
	gmtedit_axes = [handles.hCallingAx1 handles.hCallingAx2 handles.hCallingAx3];
	% -----------------------------------------------------------------

	set(handles.figure1,'pointer','watch')
	x = get(hChannel,'XData');				y = get(hChannel,'YData');
	yy = spl_fun('csaps',x,y,handles.smooth(id_gmt),x);
	difa = abs(y - yy);
	ind = (difa >= handles.thresh(id_gmt));
	clear difa
	xx = x(ind);		yy = y(ind);

	hM = findobj(handles.hCallingFig,'Type','Line','tag',markers_tag{id_gmt});
	if (isempty(hM))
		line(xx, yy,'Parent',gmtedit_axes(id_gmt),'Marker','s','MarkerFaceColor','r','MarkerSize',4,'LineStyle','none','Tag',markers_tag{id_gmt});
	else
		set(hM, 'XData',xx, 'YData',yy)
	end
	set(handles.figure1,'pointer','arrow')

% ----------------------------------------------------------------------------
function push_applyNreturn_Callback(hObject, eventdata, handles)
	push_Apply_Callback(handles.push_Apply, eventdata, handles)
	delete(handles.figure1)

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
function push_clear_Callback(hObject, eventdata, handles)
% Remove all eventually detected outlaws from current channel
	if (handles.id_gmt == 1)
		hM = findobj(handles.hCallingFig,'Type','Line','tag','GravNull');
	elseif (handles.id_gmt == 2)
		hM = findobj(handles.hCallingFig,'Type','Line','tag','MagNull');
	else
		hM = findobj(handles.hCallingFig,'Type','Line','tag','TopNull');
	end
	set(hM, 'XData',[], 'YData',[])
	

% --- Creates and returns a handle to the GUI figure. 
function outliersdetect_LayoutFcn(h1)

set(h1, 'Position',[520 420 311 65],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Detect outliers',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 27 121 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@outliersdetect_CB,h1,'edit_SmoothParam_Callback'},...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','Enter a Smoothing Parameter between [0 1]',...
'Tag','edit_SmoothParam');

uicontrol('Parent',h1, 'Position',[9 48 125 15],...
'FontName','Helvetica',...
'String','Smoothing parameter (p)',...
'Style','text');

uicontrol('Parent',h1, 'Position',[212 31 90 21],...
'Callback',{@outliersdetect_CB,h1,'push_Apply_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'String','Apply',...
'TooltipString','Use this for testing',...
'Tag','push_Apply');

uicontrol('Parent',h1, 'Position',[146 27 51 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@outliersdetect_CB,h1,'edit_thresh_Callback'},...
'String','4',...
'Style','edit',...
'Tag','edit_thresh');

uicontrol('Parent',h1, 'Position',[143 48 56 15],...
'FontName','Helvetica',...
'String','Threshold',...
'Style','text');

uicontrol('Parent',h1, 'Position',[11 6 30 15],...
'Callback',{@outliersdetect_CB,h1,'radio_G_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'String','G',...
'Style','radiobutton',...
'TooltipString','Select Gravity channel',...
'Tag','radio_G');

uicontrol('Parent',h1, 'Position',[60 6 30 15],...
'Callback',{@outliersdetect_CB,h1,'radio_M_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'String','M',...
'Style','radiobutton',...
'TooltipString','Select Magnetic channel',...
'Tag','radio_M');

uicontrol('Parent',h1, 'Position',[108 6 30 15],...
'Callback',{@outliersdetect_CB,h1,'radio_T_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'String','T',...
'Style','radiobutton',...
'TooltipString','Select Topography channel',...
'Tag','radio_T');

uicontrol('Parent',h1, 'Position',[212 7 90 21],...
'Callback',{@outliersdetect_CB,h1,'push_applyNreturn_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'String','Apply n return',...
'TooltipString','Do it and go away',...
'Tag','push_applyNreturn');

uicontrol('Parent',h1, 'Position',[147 0 50 21],...
'Callback',{@outliersdetect_CB,h1,'push_clear_Callback'},...
'FontName','Helvetica',...
'String','Clear',...
'TooltipString','Clear detections from current selected channel',...
'Tag','push_clear');

function outliersdetect_CB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

% -------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------
function gmtedit_track(varargin)
% M-File changed by desGUIDE 
 
	hObject = figure('Tag','figure1','Visible','off');
	gmtedit_track_LayoutFcn(hObject);
	handles = guihandles(hObject);
 
	handles.lon = varargin{1};
	handles.lat = varargin{2};
	handles.year = varargin{3};
	handles.hCallingAxes = varargin{4};
	handles.h_tm = varargin{5};		% Handle of the Topography line

	handles.fileName = [];

	guidata(hObject, handles);
	set(hObject,'Visible','on');

% ----------------------------------------------------------------------------
function edit_gridToInterp_Callback(hObject, eventdata, handles)
	fname = get(hObject,'String');
	if isempty(fname),		return,		end
	% Let the push_gridToInterp_Callback do all the work
	push_gridToInterp_Callback(handles.push_gridToInterp, [], handles, fname)

% ----------------------------------------------------------------------------
function push_gridToInterp_Callback(hObject, eventdata, handles, opt)
	if (nargin == 3),		opt = [];		end
	if (nargin == 4),		fname = opt;	end

	if (isempty(opt))       % Otherwise 'opt' already transmited the file name.
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.grd;*.GRD', 'Grid files (*.grd,*.GRD,*.nc)';'*.*', 'All Files (*.*)'},'Select GMT grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
		set(handles.fname,'String',fname)
	end

	handles.fname = fname;
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function zz = push_OK_Callback(hObject, eventdata, handles)
% Read the grid and interpolate this track
	if (isempty(handles.fname)),	return,		end
	if (~exist(handles.fname,'file'))
		errordlg('The grid name provided does not exist.','Error'),		return
	end

	[handles, X, Y, Z, head] = read_gmt_type_grids(handles,handles.fname);
	if (isempty(X)),		zz = [];	return,		end
	zz = grdtrack_m(Z,head,[handles.lon handles.lat],'-Z');

	if (get(handles.check_addIGRF, 'Val'))
		out = igrf_m(handles.lon, handles.lat, 0, handles.year + 0.5, '-Ft');	% Date doesn't need to be very accurate here
		zz = zz + out;
	end

	if (strcmp(get(handles.hCallingAxes,'Tag'), 'axes1'))
		h = findobj('Type','line','Tag','orig_grav');
	elseif (strcmp(get(handles.hCallingAxes,'Tag'), 'axes2'))
		h = findobj('Type','line','Tag','orig_mag');
	else
		h = findobj('Type','line','Tag','orig_topo');
	end
	x = get(h,'XData');
	line('XData',x,'YData',zz,'Color','b','Parent',handles.hCallingAxes, 'LineWidth', 2, 'HitTest', 'off');
	delete(handles.figure1)

% ----------------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function gmtedit_track_LayoutFcn(h1)

set(h1, 'Position',[520 491 351 80],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Track from grid',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[9 39 311 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@gmtedit_track_uicallback,h1,'edit_gridToInterp_Callback'},...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','Name (full name) of grid to sample along track coords',...
'Tag','edit_gridToInterp');

uicontrol('Parent',h1, 'Position',[319 38 23 23],...
'Callback',{@gmtedit_track_uicallback,h1,'push_gridToInterp_Callback'},...
'FontSize',12,...
'FontWeight','bold',...
'String','...',...
'Tag','push_gridToInterp');

uicontrol('Parent',h1, 'Position',[10 15 70 15],...
'FontName','Helvetica',...
'String','Add IGRF',...
'Style','checkbox',...
'TooltipString','For magnetic anomalies only, add an IGRF to interpolation',...
'Tag','check_addIGRF');

uicontrol('Parent',h1, 'Position',[276 9 66 23],...
'Callback',{@gmtedit_track_uicallback,h1,'push_OK_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1, 'Position',[10 62 180 16],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Grid to sample along track coords',...
'Style','text');

function gmtedit_track_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

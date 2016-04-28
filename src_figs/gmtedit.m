function out = gmtedit(varargin)
% Revival of the ancient Sunview gmtedit.
%
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

%	Copyright (c) 2004-2015 by J. Luis
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

% $Id: gmtedit.m 4760 2015-07-28 23:51:57Z j $

	f_name = '';
	got_inFile = 0;
	def_width_km = 200;         % Default width in km (approximatly 1/2 of a day at 10 knots)
	opt_G = '-G';               % Default output to [-180/+180] range
	begin = 1;                  % This means that the display starts at the begining of track
	center_win = [];
	vars = [];					% Container to eventual user selected fields to plot (MGD77+ only)
	multi_plot = [];			% Container to eventual user selected extra fields to overlay in plot (MGD77+ only)
	xISdist = true;				% Default absicssae to distance in km
	is_gmt = false;
	is_mgd77 = false;
	hMirAxes = [];				% To hold the Mirone axes handle (case exists)
	MIRONE_DIRS = getappdata(0,'MIRONE_DIRS');

	if (nargin > 0)
		f_name = varargin{1};
		[PATH,FNAME,EXT] = fileparts(f_name);
		if (isempty(EXT))       % Remember that here we need the extension
			fid = fopen([f_name '.nc'], 'r');
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

		if (is_mgd77)
			[vars, multi_plot, xISdist] = parse_optV(MIRONE_DIRS.home_dir);
		end

		if (~isempty(varargin) && ishandle(varargin{end}) && strcmp(get(varargin{end},'type'),'axes'))
			hMirAxes = varargin{end};
			varargin(end) = [];				% Remove it to not complicate switch case below
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
					[tok,r] = strtok(str,'/');
					lon = str2double(tok);
					lat = str2double(r(2:end));
					begin = 0;
					center_win = [lon lat];
			end
		end
	end

	sc_size = get(0,'ScreenSize');
	fig_height = fix(sc_size(4)*.9);
	marg_l = 60;                            % Left margin
	marg_r = 15;                            % Right margin
	marg_tb = 10;                           % Top & Bottom margins
	marg_ax = 15;                           % Margin between axes
	ax_height = fix(fig_height * .295);
	ax_width = sc_size(3) - marg_l - marg_r;
	fp = [0 0.035 1 0.965];

	hf = figure('name','gmtedit','numbertitle','off', 'visible','off', 'units','normalized',...
        'outerposition',fp, 'DoubleBuffer','on', 'Tag','figure1', 'closerequestfcn',@fig_CloseRequestFcn, ...
		'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Apply this trick to get the icons. Let's hope that this is not version/OS dependent 
	hA = findall(hf);
	hh = findobj(hA,'Tooltip','Open File');			openFile_img = get(hh(1),'CData');
	hh = findobj(hA,'Tooltip','Save Figure');		saveFile_img = get(hh(1),'CData');
	hh = findobj(hA,'Tooltip','Zoom In');			zoomIn_img   = get(hh(1),'CData');
	hh = findobj(hA,'Tooltip','Zoom Out');			zoomOut_img  = get(hh(1),'CData');
	set(hf,'menubar','none','units','pixel')		% Set the menubar to none

	handles = guihandles(hf);

	if (~isempty(MIRONE_DIRS))							% Should not be empty, but ...
		handles.home_dir = MIRONE_DIRS.home_dir;		% False info if not called from Mir root dir
		handles.last_dir = MIRONE_DIRS.last_dir;
		handles.work_dir = MIRONE_DIRS.work_dir;
	else
		handles.home_dir = cd;
		handles.last_dir = cd;
		handles.work_dir = handles.last_dir;
	end

	if (got_inFile),		handles.last_dir = PATH;	end

	% Load some icons from mirone_icons.mat
	s = load([handles.home_dir filesep 'data' filesep 'mirone_icons.mat'],'rectang_ico','info_ico','trincha_ico');
	link_ico = make_link_ico;

	h_toolbar = uitoolbar('parent',hf,'Clipping', 'on', 'BusyAction','queue','HandleVisibility','on',...
		'Interruptible','on','Tag','FigureToolBar','Visible','on');
	uipushtool('parent',h_toolbar,'Click',{@import_clickedCB,f_name},'Tag','import',...
		'cdata',openFile_img,'Tooltip','Open gmt file');
	uipushtool('parent',h_toolbar,'Click',@save_clickedCB,'Tag','save', 'cdata',saveFile_img,'Tooltip','Save gmt file');
	uipushtool('parent',h_toolbar,'Click',@info_clickedCB,'Tag','info','cdata',s.info_ico, 'Tooltip','Cruise Info');
	uipushtool('parent',h_toolbar,'Click',@rectang_clickedCB,'Tag','rectang', 'cdata',s.rectang_ico, ...
		'Tooltip','Rectangular region','Sep','on');
	s.rectang_ico(7:10,8,:) = 0;
	uipushtool('parent',h_toolbar,'Click',@rectangMove_clickedCB,'cdata',s.rectang_ico,'Tooltip','Select for moving');
	uipushtool('parent',h_toolbar,'Click',{@changeScale_clickedCB,'inc'}, 'cdata',zoomIn_img,'Tooltip','Increase scale','Sep','on');
	uipushtool('parent',h_toolbar,'Click',{@changeScale_clickedCB,'dec'}, 'cdata',zoomOut_img,'Tooltip','Decrease scale');
	uipushtool('parent',h_toolbar,'Click',@outliers_clickedCB, 'cdata',s.trincha_ico,'Tooltip','Outliers detector','Sep','on');
	if (~is_gmt)				% Not an old .gmt file
		uipushtool('parent',h_toolbar,'Click',{@NavFilters_ClickedCB,f_name}, ...
			'cdata',flipdim(s.trincha_ico,1),'Tooltip','Find Nav/Grad troubles','Sep','on');
	end
	if (~isempty(hMirAxes))
		uitoggletool('parent',h_toolbar,'Click',@ptcoords_clickedCB, 'cdata',link_ico, ...
			'Tooltip','Pick data point in curve and plot it the mirone figure','Sep','on');
	end

	pos = [marg_l fig_height-ax_height-marg_tb ax_width ax_height];
	handles.axes1 = axes('Parent',hf, 'Units','pixels', 'Position',pos, 'XLim',[0 def_width_km], 'Tag','axes1');

	pos(2) = fig_height-2*(ax_height+marg_tb)-marg_ax;
	handles.axes2 = axes('Parent',hf, 'Units','pixels', 'Position',pos, 'XLim',[0 def_width_km], 'Tag','axes2');

	pos(2) = fig_height-3*(ax_height+marg_tb)-2*marg_ax;
	handles.axes3 = axes('Parent',hf,'Units','pixels', 'Position',pos, 'XLim',[0 def_width_km], 'Tag','axes3');

	if (isempty(vars))		% Default GMT axes names
		set(get(handles.axes1,'YLabel'),'String','Gravity anomaly (mGal)')
		set(get(handles.axes2,'YLabel'),'String','Magnetic field (nT)')
		set(get(handles.axes3,'YLabel'),'String','Bathymetry (m)')
	else
		set(get(handles.axes1,'YLabel'),'String',vars{1})
		set(get(handles.axes2,'YLabel'),'String',vars{2})
		set(get(handles.axes3,'YLabel'),'String',vars{3})
	end

	handles.def_width_km = def_width_km;
	handles.max_x_data = 10000;
	scroll_plots(def_width_km,[0 10000])     % The [0 10000] will be reset inside scroll_plots to a more appropriate val

	% Create empty lines just for the purpose of having their handles
	handles.h_gm = line('XData',[],'YData',[],'LineStyle','none','Marker','s', ...
		'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',4,'Parent',handles.axes1);
	handles.h_mm = line('XData',[],'YData',[],'LineStyle','none','Marker','s', ...
		'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',4,'Parent',handles.axes2);
	handles.h_tm = line('XData',[],'YData',[],'LineStyle','none','Marker','s', ...
		'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',4,'Parent',handles.axes3);

	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set([handles.axes1 handles.axes2 handles.axes3], 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Overlay interpolation', 'Call',@interp_on_grid);

	handles.opt_G = opt_G;
	handles.info = [];
	handles.h_broken = [];
	handles.begin = begin;
	handles.center_win = center_win;
	handles.is_gmt = is_gmt;
	handles.is_mgd77 = is_mgd77;
	handles.vars = vars;
	handles.xISdist = xISdist;
	handles.multi_plot = multi_plot;
	handles.force_gmt = false;		% If TRUE save in old .gmt format
	handles.got_vel = 0;			% To tell if we have a Vel channel and where.
	handles.adjustedVel = false;	% Will be set to true if any velocity (coords) is recalculated
	handles.hMirAxes = hMirAxes;	% Non-empty when called from draw_funs

	guidata(hf, handles);

	% Now that we have a figure and its handles, we can open the file if it was requested on input
	if (got_inFile),	import_clickedCB(hf,[],f_name),		end

	% Add or remove red Markers
	set(hf,'WindowButtonDownFcn',@add_MarkColor,'Visible','on')

	% Choose default command line output for gmtedit
	if (nargout == 1),	out = hf;   end

% --------------------------------------------------------------------------
function [vars, multi_plot, xISdist] = parse_optV(home_dir)

	opt_V = false;
	vars = [];					% Container to eventual user selected fields to plot (MGD77+ only)
	multi_plot = [];			% Container to eventual user selected extra fields to overlay in plot (MGD77+ only)
	xISdist = true;				% Default absicssae to distance in km

	opt_file = [home_dir filesep 'data' filesep 'OPTcontrol.txt'];
	if ( exist(opt_file, 'file') == 2 )
		fid = fopen(opt_file, 'r');
		c = (fread(fid,'*char'))';      fclose(fid);
		lines = strread(c,'%s','delimiter','\n');   clear c fid;
		m = numel(lines);
		for (k = 1:m)
			if (~strncmp(lines{k},'MIR_GMTEDIT',11)),	continue,	end
			if (numel(lines{k}) <= 14),	continue,	end		% The minimum it takes to have a -? switch
			t = strtok(lines{k}(13:end));
			if ( strcmp(t(1:2), '-V') )		% Here we only check for a -V... and do not check for errors
				opt_V = t(3:end);
				break
			end
		end
	end

	if (opt_V)
		% For new MGD77+ files only. It contains a user selection of the fields to be plot (instead of GMT)
		% USAGE (example) 
		%	-Vvar1,var2,var3
		%	-Vvar1,var2,var3[:varI+slotI[/varI+slotI[/varI+slotI]]][|]
		% Where var1,var2,var3 is the name of any of the variables in the .nc file
		% The varI+slotI form will plot (add) the varI (again, one of the .nc variables) to the slot I (1:3) axes
		% Append the optional end '|' char to indicate that absicssae is record number (fids) instead of distances in km
		if (opt_V(end) == '|')	% Absicssae will be the "fiducial" numbers
			xISdist = false;	opt_V(end) = '';
		end
		ind2 = strfind(opt_V,':');
		if (isempty(ind2)),		last = numel(opt_V);
		else					last = ind2 - 1;
		end

		vars = {'faa' 'mtf1' 'depth'};	% Used when (ex) '-V:anom+2' or -Vvar1,var2,var3[:...] forms
		if (~isempty(ind2) && ind2(1) == 1)		% The form '-V:anom+2' is also aceptable
		else
			ind = strfind(opt_V,',');
			if (numel(ind) == 2)
				vars{1} = opt_V(1:ind(1)-1);
				vars{2} = opt_V(ind(1)+1:ind(2)-1);
				vars{3} = opt_V(ind(2)+1:last);
			end
		end

		if (~isempty(ind2))				% Have a multiple plots request
			str = opt_V(ind2+1:end);	% Retain the multi-plot info string only
			ind = strfind(str,'/');		% How many extra curves?
			multi_plot = cell(numel(ind)+1,2);		% Pre-allocate
			ind3 = [1 ind numel(str)+1];
			for (n = 1:numel(ind3)-1)
				str2 = str(ind3(n):ind3(n+1)-1);
				ind = strfind(str2,'+');% To find in which slot will go this plot
				multi_plot{n,1} = str2(1:ind-1);
				multi_plot{n,2} = str2(ind+1:end);
			end
		end
	end

% -----------------------------------------------------------------------------
function ico = make_link_ico()
% construct the link icon (not yet in mirone_icons.mat)
	ico	= uint8([	226  239  238  238  238  238  238  238  238  238  238  238  238  238  240  221; ...
					240  254  253  253  253  253  253  253  253  255  255  255  253  253  254  236; ...
					236  251  250  250  250  250  250  250  253  244  156  183  255  250  251  233; ...
					235  247  247  247  247  247  247  250  245  132  145  124  167  254  248  232; ...
					234  244  244  244  244  244  245  240  127  147  250  224  119  176  251  231; ...
					232  240  239  239  239  239  248  148  145  189  243  255  141  144  248  230; ...
					232  238  238  238  239  244  247  198  151  107  189  149  113  223  241  228; ...
					229  232  232  233  236  137  187  154  107  139  137  107  217  236  232  227; ...
					225  221  223  228  110  124  139  100  131  188  129  207  226  221  221  223; ...
					226  223  227  106  114  199   84  129  174  233  231  224  222  222  222  224; ...
					227  233  135  105  230  238  201  116  121  233  225  225  225  225  225  225; ...
					228  233  154   81  206  244  116   88  219  227  225  225  225  225  225  226; ...
					227  227  246  126   77  103   85  223  231  228  228  228  228  228  228  225; ...
					228  230  229  247  149  119  221  233  230  230  230  230  230  230  230  226; ...
					230  232  232  232  239  241  234  232  232  232  232  232  232  232  232  228; ...
					218  229  228  228  228  228  228  228  228  228  228  228  228  228  229  216]);
	  ico = repmat(ico,[1 1 3]);

% -----------------------------------------------------------------------------
function fig_CloseRequestFcn(hObj, event)
	handles = guidata(hObj);
	try		h = getappdata(handles.figure1,'Filhas');
	catch,	delete(handles.figure1),	return
	end
	delete(handles.figure1);		delete(h(ishandle(h)))      % Delete also any eventual 'carra?as'

% --------------------------------------------------------------------
function import_clickedCB(hObject, evt, opt)
	handles = guidata(hObject);     % get handles
	if (isempty(opt))
		[FileName,PathName] = put_or_get_file(handles,{'*.gmt;*.GMT;*.nc', 'gmt files (*.gmt,*.GMT,*.nc)'},'Select gmt File','get');
		if isequal(FileName,0),		return,		end
		f_name = [PathName FileName];
	else
		f_name = opt;
	end

	[PATH,FNAME,EXT] = fileparts(f_name);
	if (isempty(PATH)),		PATH = cd;		end		% Play safe
	if (strcmpi(EXT,'.nc'))
		handles.f_name = [PATH filesep FNAME '.nc'];		handles.is_mgd77 = true;		handles.is_gmt = false;
		[handles, track] = read_mgd77_plus(handles, f_name);
		MIRONE_DIRS = getappdata(0,'MIRONE_DIRS');
		[handles.vars, handles.multi_plot, handles.xISdist] = parse_optV(MIRONE_DIRS.home_dir);
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
	if (length(track.agency) ~= 10)			% Ensures that agency is exactly 10 chars
		agency = '          ';				% 10 blanks
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
		handles = rmfield(handles,'h_broken');
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

	if (~handles.begin)         % Start the display at a user selected coordinate
		r1 = 1;		r2 = 1;		n_iter = 0;
		x = handles.lon - handles.center_win(1);		y = handles.lat - handles.center_win(2);
		while ( (min(r1, r2) > 0.02) && (n_iter < 10) )	% We need this loop for cases when both id1 & id2 are bad
			[zz,id1] = min(abs(x));			[zz,id2] = min(abs(y));
			% id1 and id2 are not forcedly equal. Find out the "best"
			r1 = sqrt((handles.lon(id1) - handles.center_win(1))^2 + (handles.lat(id1) - handles.center_win(2))^2);
			r2 = sqrt((handles.lon(id2) - handles.center_win(1))^2 + (handles.lat(id2) - handles.center_win(2))^2);
			id = id1;
			if (r2 < r1),	id = id2;	end
			x(id1) = 1e2;	y(id2) = 1e2;		% If they are no good, don't reuse them, otherwise it doesn't matter
			n_iter = n_iter + 1;
		end
		clear x y;
		x_lim = track.distance(id) + [-handles.def_width_km/2 handles.def_width_km/2];
		set([handles.axes1 handles.axes2 handles.axes3],'xlim',x_lim)
		val0 = track.distance(id)-handles.def_width_km;
		if (val0 > track.distance(1)),		val = val0;		end
		line('XData',[track.distance(id) track.distance(id)], 'YData',get(handles.axes1,'ylim'), 'Parent', handles.axes1)
		line('XData',[track.distance(id) track.distance(id)], 'YData',get(handles.axes2,'ylim'), 'Parent', handles.axes2)
		line('XData',[track.distance(id) track.distance(id)], 'YData',get(handles.axes3,'ylim'), 'Parent', handles.axes3)
	end

	% If we have multi-plots requests
	if (~isempty(handles.multi_plot))
		for (k = 1:size(handles.multi_plot,1))
			switch handles.multi_plot{k,2}
				case '1'
					if (numel(track.distance) == numel(track.multi{k}))
						line('XData',track.distance, 'YData',track.multi{k}, 'Parent', handles.axes1)
					end
				case '2'
					if (numel(track.distance) == numel(track.multi{k}))
						line('XData',track.distance, 'YData',track.multi{k}, 'Parent', handles.axes2)
					end
				case '3'
					if (numel(track.distance) == numel(track.multi{k}))
						line('XData',track.distance, 'YData',track.multi{k}, 'Parent', handles.axes3)
					end
			end
		end
	end

	% Update the slider properties
	r = 0.90 * handles.def_width_km / handles.max_x_data;		% Give a 10% overlap between display area when click on slider arrow
	if (r > 1),		set(hs,'Vis','off')			% Would otherwise generate an error
	else			set(hs,'Max',max_s,'value',val,'SliderStep',[r r*10])
	end
	guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function interp_on_grid(obj, event)
	handles = guidata(obj);
	hAx = get(handles.figure1, 'CurrentAxes');		% Find which axes are we working on
	h = gmtedit_track(handles.lon, handles.lat, handles.year, hAx, handles.h_tm, handles.last_dir, handles.work_dir, handles.hMirAxes);
	filhas = getappdata(handles.figure1,'Filhas');
	setappdata(handles.figure1,'Filhas',[filhas(:); h])

% --------------------------------------------------------------------
function NavFilters_ClickedCB(obj, event, fname)
% Rise a Window with controls for filtering Nav/Grads
	handles = guidata(obj);
	if ( isempty(get(handles.h_mm,'YData')) )
		errordlg('Warning. This option operates only on magnetic data ... and you don''t have any here.','Warnerr')
		return
	end
	h = gmtedit_NavFilters(handles.axes2, handles.h_mm, fname);
	filhas = getappdata(handles.figure1,'Filhas');
	setappdata(handles.figure1,'Filhas',[filhas(:); h])

% --------------------------------------------------------------------
function [handles, track] = read_mgd77_plus(handles, fname)
	s = nc_funs('info',fname);

	% ------------------ OK, Get numerics now -----------------------------------
	tempo = [];
	track.longitude = double(nc_funs('varget', fname, 'lon'));
	track.latitude  = double(nc_funs('varget', fname, 'lat'));
	if (isempty(handles.vars))
		track.magnetics  = nc_funs('varget', fname, 'mtf1');
		track.gravity    = nc_funs('varget', fname, 'faa');
		track.topography = nc_funs('varget', fname, 'depth');
	else
		try
			ind_v = strncmp('vel', handles.vars, 3);		% Was velocity asked for?
			if ( any(ind_v) )		% A bit boring below because 'vel' is not a MGD77+ Variable. So we need to avoid error
				ind_v = find(ind_v);
				tempo = nc_funs('varget', fname, 'time');
				if (ind_v == 1)
					track.magnetics  = nc_funs('varget', fname, handles.vars{2});
					track.topography = nc_funs('varget', fname, handles.vars{3});
				elseif (ind_v == 2)
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

			% If there is a multi plot request, try to fetch the requested extra fields
			if (~isempty(handles.multi_plot))
				try
					for (k = 1:size(handles.multi_plot,1))
						track.multi{k} = nc_funs('varget', fname, handles.multi_plot{k,1});
					end
				catch
					errordlg('GMTEDIT: At least one of custom selected EXTRA field name does not exist in dataset','Error')
					if (isfield(track,'multi')),	track = rmfield(track,'multi');		end
				end
			end
		catch
			errordlg('GMTEDIT: At least one of custom selected field name does not exist in dataset','Error')
			error('GMTEDIT:read_mgd77_plus', lasterr)
		end
	end

	D2R = pi / 180;		KMPRDEG = 111.1949;
	co = cos(track.latitude * D2R);
	dx = [0; diff(track.longitude)] .* co;		dy = [0; diff(track.latitude)];
	dist = cumsum(sqrt(dx .^2 + dy .^2) * KMPRDEG);

	if (~isempty(tempo))				% Compute velocity and find where to store it
		vel = [diff(dist(1:2)) ./ diff(tempo(1:2)); diff(dist) ./ diff(tempo)] * 1e3 / 1852 * 3600;		% Speed in knots
		if (ind_v == 1),			track.gravity = vel;		% Remember that grav, mag, topo are default field names 
		elseif (ind_v == 2),		track.magnetics = vel;
		else						track.topography = vel;
		end
		handles.got_vel = ind_v;		% Store in which axes we have the velocity
	end
	if (handles.xISdist)				% Plot against accumulated distance
		track.distance = dist;
	else								% Plot against the record number (the 'fids')
		track.distance = 1:numel(track.gravity);
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
	track.info = 'Bla Bla';

% --------------------------------------------------------------------
function save_clickedCB(hObject, evt)
	handles = guidata(hObject);     % get handles
	if (handles.is_gmt || handles.force_gmt)
		NODATA(1:3) = -32000;
		f_name = handles.f_name;
		if (handles.force_gmt && ~handles.is_gmt)		% Conversion from the mgd77+ format
			[PATH,FNAME] = fileparts(handles.f_name);
			f_name = [PATH filesep FNAME '.gmt'];
		end
		[FileName,PathName] = put_or_get_file(handles, f_name,'Select gmt File', 'put','.gmt');
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

	% Find out which variables we need to save (and NOT save by mistake)
	saveGRAV = false;		saveMAG = false;		saveTOPO = false;
	if (isempty(handles.vars))
		saveGRAV = ~isempty(y_gn);		saveMAG = ~isempty(y_mn);		saveTOPO = ~isempty(y_tn);
	else
		for (k = 1:numel(handles.vars))		% Would be more elegant with a strmatch but this is lighter
			if ( strcmp(handles.vars{k}, 'faa')   && ~isempty(y_gn) && k == 1 ),		saveGRAV = true;	end
			if ( strcmp(handles.vars{k}, 'mtf1') && k == 2 && (~isempty(y_mn) || ~isempty(handles.h_broken)) )
				saveMAG = true;
			end
			if ( strcmp(handles.vars{k}, 'depth') && ~isempty(y_tn) && k == 3 ),		saveTOPO = true;	end
		end
	end

	if (~isempty(handles.h_broken))			% If we have broken lines we must join them
		x_broken = [];    y_broken = [];
		for (k = 1:numel(handles.h_broken))
			x_broken = [x_broken get(handles.h_broken(k),'XData')];
			y_broken = [y_broken get(handles.h_broken(k),'YData')];
		end
		[x_m,id] = sort([x_m x_broken]);
		y_m = [y_m y_broken];
		y_m = y_m(id);						% Otherwise the y's would be out of order
		saveMAG = true;
		if (isempty(y_mn)),		y_mn = 0;	end	% It can't be empty on the test to update the .nc file
	end

	set(handles.figure1,'Pointer','watch')
	if (~isempty(x_gn))						% Find and remove marked GRAV points
		for (k=1:numel(x_gn)),		y_g((x_g - x_gn(k)) == 0) = NaN;		end
	end
	if (~isempty(x_mn))						% Find and remove marked MAG points
		for (k=1:numel(x_mn)),		y_m((x_m - x_mn(k)) == 0) = NaN;		end
	end
	if (~isempty(x_tn))						% Find and remove marked TOPO points
		for (k=1:numel(x_tn)),		y_t((x_t - x_tn(k)) == 0) = NaN;		end
	end
	
	if (~isempty(y_g))
		if (handles.is_gmt),	y_g(isnan(y_g*10)) = NODATA(1);		y_g = int16(y_g);
		elseif (~isempty(x_gn)),y_g = y_g / handles.gravScaleF;		y_g(isnan(y_g)) = NODATA(1);	y_g = int16(y_g);
		end
	end
	if (~isempty(y_m))
		hD = findobj(handles.figure1,'Type','Line','tag','despikado');	% Do we have "Despiked" points?
		if (~isempty(hD))			% Yes, we have them so update the y_m vector with their values
			yD = get(hD, 'YData');			ud = get(hD,'UserData');
			y_m(ud) = yD;
			saveMAG = true;
		end
		is_gmt = handles.is_gmt;	% Local copy to handle also the force_gmt case 
		if (handles.force_gmt && ~is_gmt),		y_m = y_m - 40000;	is_gmt = true;		end		% Conversion from the mgd77+ format
		if (is_gmt)
			y_m(isnan(y_m)) = NODATA(2);		y_m = int16(y_m);
		elseif ( ~isempty(x_mn) || ~isempty(handles.h_broken) || ~isempty(hD) )
			y_m = y_m / handles.magScaleF;
			y_m(isnan(y_m)) = NODATA(2);
			y_m = int32(y_m);
		end
	end
	if (~isempty(y_t))
		if (handles.is_gmt),	y_t(isnan(y_t)) = NODATA(3);		y_t = int16(y_t);
		elseif (~isempty(x_tn)),y_t = y_t / handles.topoScaleF;		y_t(isnan(y_t)) = NODATA(3);	y_t = int32(y_t);
		end
	end

	n_rec = numel(handles.lon);
	if ( (handles.is_gmt || handles.force_gmt) && (isempty(y_gn) || isempty(y_mn) || isempty(y_tn)) )
		if (isempty(y_g)),		y_g = repmat(int16(NODATA(1)),1,n_rec);		end		% Original had not this type of data
		if (isempty(y_m)),		y_m = repmat(int16(NODATA(1)),1,n_rec);		end
		if (isempty(y_t)),		y_t = repmat(int16(NODATA(1)),1,n_rec);		end
	end

	if (handles.is_gmt || handles.force_gmt)			% Old style .gmt files
		if (handles.force_gmt && ~handles.is_gmt)		% Conversion from the mgd77+ format
			tempo = nc_funs('varget', handles.f_name, 'time');
			handles.time = tempo - (date2jd(handles.year) - date2jd(1970)) * 86400;	% We need time in seconds since begining of year
			[PATH,FNAME] = fileparts(f_name);
			f_name = [PATH filesep FNAME '.gmt'];
		end
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
	else						% New mgd77+ netCDF style files
		if (saveGRAV),		nc_funs('varput', f_name, 'faa',   y_g, 0);		end		% The 0 flags nc_funs for not try to scale
		if (saveMAG),		nc_funs('varput', f_name, 'mtf1',  y_m, 0);		end
		if (saveTOPO),		nc_funs('varput', f_name, 'depth', y_t, 0);		end
		if (handles.adjustedVel)	% We had Nav changes
			nc_funs('varput', f_name, 'lon',  handles.lon);
			nc_funs('varput', f_name, 'lat',  handles.lat);
		end
	end

	set(handles.figure1,'Pointer','arrow')

% --------------------------------------------------------------------
function out = add_MarkColor(hObject, evt)
% Add a red Marker over the closest (well, near closest) clicked point.
% When used with an output arg, interpret as "Get index of clicked point only"
% and that index is returned in OUT

	handles = guidata(hObject);     % get handles

	do_despika   = false;			% For spike detection/killing
	only_find_pt = false;			% To just find the index of closes clicked pt

	if (nargout)
		only_find_pt = true;	out = [];
	end
	button = get(handles.figure1, 'SelectionType');
	if (strcmp(button,'extend'))		% Shift-click to despike a point
		do_despika = true;
	elseif (~strcmp(button,'normal'))	% Than accept only left-clicks
		return
	end

	in_grav = 0;    in_mag = 0;     in_topo = 0;
	hAx = get(handles.figure1,'CurrentAxes');
	pt = get(hAx, 'CurrentPoint');
	if (strcmp(get(hAx,'Tag'),'axes1'))
		in_grav = 1;		opt = 'GravNull';		currHline = handles.h_gm;
	elseif (strcmp(get(hAx,'Tag'),'axes2'))
		in_mag = 1;			opt = 'MagNull';		currHline = handles.h_mm;
	elseif (strcmp(get(hAx,'Tag'),'axes3'))
		in_topo = 1;		opt = 'TopNull';		currHline = handles.h_tm;
	end
	if ((in_grav + in_mag + in_topo) == 0),		return,		end		% Click was outside axes

	hM = findobj(hAx,'Type','Line','tag',opt);
	x = get(currHline,'XData');		y = get(currHline,'YData');

	x_lim = get(hAx,'XLim');		y_lim = get(hAx,'YLim');
	dx = diff(x_lim) / 20;			% Search only betweem +/- 1/10 of x_lim
	dx = max(dx, 10);
	id = (x < (pt(1,1)-dx) | x > (pt(1,1)+dx));
	if (do_despika || only_find_pt)
		fids = find(~id);		% W'll need this info to reconstruct the true fid of the clicked point 
	end
	x(id) = [];					y(id) = [];	% Clear outside-2*dx points to speed up the search code
	XScale = diff(x_lim);		YScale = diff(y_lim)*6;		% The 6 factor compensates the ~6:1 horizontal/vertical axes dimension

	r = sqrt(((pt(1,1)-x)./XScale).^2+((pt(1,2)-y)./YScale).^2);
	[mini,i] = min(r);
	if (mini > 0.04),			return,		end				% Do not accept click if its too far away of the closest point
	ptClicked_x = x(i);			ptClicked_y = y(i);

	if (do_despika)
		if ( handles.got_vel && strcmp( get(hAx,'Tag'), sprintf('axes%d',handles.got_vel) ) )	% Velocity recomp
			despikeNav(hAx, currHline, (fids(1) + i - 1))
		else
			despika(hAx, currHline, (fids(1) + i - 1))
		end
		return
	elseif (only_find_pt)		% We are done here.
		out = fids(1) + i - 1;
		return
	end

	xr = get(hM,'XData');		yr = get(hM,'YData');		% Red markers
	if ( ~isempty(xr) && ~isempty(ptClicked_x) )
		id = find(xr == ptClicked_x);
		if (~isempty(id))		% Marker already exists. Kill it and return
			xr(id) = [];		yr(id) = [];
			set(hM,'XData',xr, 'YData', yr)
			return
		end
	end

	% New Marker
	if (isempty(hM))			% First red Marker on this axes
		line(ptClicked_x,ptClicked_y,'Marker','s','MarkerFaceColor','r','MarkerSize',4,'LineStyle','none','Tag',opt);
	else
		xr = [xr ptClicked_x];		yr = [yr ptClicked_y];
		set(hM,'XData',xr, 'YData', yr)
	end

% --------------------------------------------------------------------------------------------------
function despikeNav(hAx, hLine, n)
% Try to remove a spike in Vel by interpolation of its faulty coords
	handles = guidata(hAx);
	if (~handles.xISdist)
		warndlg('Sorry. Cannot recompute velocity when X axes has record numbers instead of ditances. Bye.','Warning')
		return
	end
	if (n == 1 || n == numel(handles.lon))
		warndlg('First or last point cannot be recalculated. Bye.','Warning'),		return
	end
	ind1 = max(n-3,1);			ind2 = min(n+3,numel(handles.lon));
	lon = handles.lon(ind1:ind2);			lat = handles.lat(ind1:ind2);
	meanDLon = abs(lon(end) - lon(1)) / (numel(lon) - 1);
	meanDLat = abs(lat(end) - lat(1)) / (numel(lat) - 1);
	lonInt = handles.lon(n+1) - handles.lon(n-1);
	latInt = handles.lat(n+1) - handles.lat(n-1);
	if (lonInt > 4*meanDLon || latInt > 4*meanDLat)		% Ad-hoc test, but we must test for something resonable.
		warndlg('Sorry. Coordinate data arround clicked point does not allow a reasonable estimate of new Lon/Lat. Likely close to a data gap','Warning')
		return
	end
	newLon = handles.lon(n-1) + lonInt/2;
	newLat = handles.lat(n-1)  + latInt/2;
	handles.lon(n) = newLon;
	handles.lat(n)  = newLat;
	handles.adjustedVel = true;
	guidata(hAx, handles)

	% Recompute the velocity for this point but we must use a trick to recover the time which we don't store anywhere
	co = cos(newLat *  pi / 180);
	dx = (newLon - handles.lon(n+1)) .* co;		dy = newLat - handles.lat(n-1);
	dist = sqrt(dx .^2 + dy .^2) * 111.1949;
	x = get(hLine,'XData');			y = get(hLine,'YData');
	dt = (x(n) - x(n-1)) / y(n);			% Recover the time from the previous velocity value
	newVel = dist / dt;						% Speed in knots. Notice that dt above was obtained by km / (NM / h)
	x(n) = x(n-1) + dist;			y(n) = newVel;
	set(hLine, 'XData', x, 'YData', y)

% --------------------------------------------------------------------------------------------------
function despika(hAx, hLine, n)
% Move one isolated point (spike) onto its most likely position as computed from neighbors
	if (~strcmp(get(hAx,'Tag'), 'axes2'))
		warndlg('Sorry, relocation by interpolation only works when the magnetic channel is plotted in the middle axes','Warning')
		return
	end

	x = get(hLine,'XData');		y = get(hLine,'YData');
	ind1 = max(n-4,1);			ind2 = min(n+4,numel(x));
	Xs = x(ind1:ind2);			Ys = y(ind1:ind2);
	ind = find(Xs == x(n));
	Xs(ind) = [];				Ys(ind) = [];		% Remove the clicked point. Easier this way
	ind = isnan(Ys);
	if (~isempty(ind))			% We cannot have NaNs in the vectors sent to interp1
		Xs(ind) = [];			Ys(ind) = [];
		if (numel(Xs) < 4),		return,		end
	end
	
	newY = interp1(Xs,Ys,x(n), 'spline');
	
	hD = findobj(hAx,'Type','Line','tag','despikado');
	if (isempty(hD))				% First despiked Marker on this axes
		hD = line(x(n),newY,'Parent',hAx,'Marker','o','MarkerFaceColor','r','MarkerSize',5,'LineStyle','none','Tag','despikado');
		set(hD,'UserData',false(1,numel(x)))
	else
		xd = get(hD,'XData');		yd = get(hD,'YData');		% Circle markers
		[xd,id] = sort([xd x(n)]);	yd = [yd newY];
		set(hD,'XData',xd, 'YData', yd(id))
	end
	
	ud = get(hD,'UserData');		% Store the FID of the relocated point in UserData
	ud(n) = true;
	set(hD,'UserData',ud)

	y(n) = NaN;
	set(hLine, 'YData', y)
	

% --------------------------------------------------------------------------------------------------
function info_clickedCB(obj,evt)
	handles = guidata(obj);     % get handles
	if (isempty(handles.info)),		return,		end
	if (handles.is_gmt)
		data = handles.info;
	else
		[data, agency] = aux_funs('mgd77info',handles.f_name);
	end

	str{1} = sprintf('N_recs = %d, N_grav = %d, N_mag = %d, N_topo = %d', data(1:4));
	str{2} = ['E: = ' num2str(data(5)) '  W: = ' num2str(data(6))];
	str{3} = ['S: = ' num2str(data(7)) '  N: = ' num2str(data(8))];
	str{4} = ['Start day,month,year: = ' num2str(data(9)) '  ' num2str(data(10)) '  ' num2str(data(11))];
	str{5} = ['End   day,month,year: = ' num2str(data(12)) '  ' num2str(data(13)) '  ' num2str(data(14))];
	if (~handles.is_gmt),	str{6} = '';	str{7} = agency;	end
	msgbox(str,'Cruise Info')

% --------------------------------------------------------------------------------------------------
function rectang_clickedCB(obj,evt)
	handles = guidata(obj);     % get handles
	try
		[p1,p2] = rubberbandbox;
	catch		% Don't know why but uisuspend sometimes breaks
		set(handles.figure1,'Pointer','arrow')
		return
	end

	in_grav = 0;	in_mag = 0;
	hAx = get(handles.figure1,'CurrentAxes');
	if (strcmp(get(hAx,'Tag'),'axes1')),			in_grav = 1;	opt = 'GravNull';
	elseif (strcmp(get(hAx,'Tag'),'axes2')),		in_mag = 1;		opt = 'MagNull';
	else	opt = 'TopNull';
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
	if (isempty(id)),	return,		end		% Nothing inside rect
	golo_x = x(id);				golo_y = y(id);

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
function rectangMove_clickedCB(obj,evt)
	handles = guidata(obj);		% get handles
	
	try			[p1,p2] = rubberbandbox(handles.axes2);
	catch,		return
	end
	
	if (~strcmp(get(get(handles.figure1,'CurrentAxes'),'Tag'),'axes2'))		% This option is meant only to magnetic data
		return
	end

	x = get(handles.h_mm,'XData');		y = get(handles.h_mm,'YData');

	% Search for "regular" points inside the rectangle
	id = find(x >= p1(1,1) & x <= p2(1,1) & y >= p1(1,2) & y <= p2(1,2));
	if (isempty(id)),	return,		end		% Nothing inside rect

	% create a new line and set it with the ui_edit_polygon so that it can be moved
	golo_x = x(id);				golo_y = y(id);

	x(id) = [];         y(id) = [];			% Remove these from the main line
	set(handles.h_mm,'XData',x,'YData',y)

	n = numel(handles.h_broken);
	handles.h_broken(n+1) = line(golo_x,golo_y,'Marker','s','MarkerFaceColor','b','MarkerSize',4,'LineStyle','-');
	ui_edit_polygon(handles.h_broken(n+1), 'y')		% Set edition functions

	guidata(handles.figure1, handles)

% --------------------------------------------------------------------------------------------------
function changeScale_clickedCB(obj,evt,opt)
	handles = guidata(obj);		% get handles

	if (strcmp(opt,'inc'))
		handles.def_width_km = handles.def_width_km - 50;
		handles.def_width_km = max(handles.def_width_km, 25);
		if (handles.def_width_km < 25),		return,		end
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

	r = 0.90 * handles.def_width_km / handles.max_x_data;
	set(h_slider,'SliderStep', [r r*10])
	
	guidata(handles.figure1, handles)

% --------------------------------------------------------------------------------------------------
function outliers_clickedCB(obj,evt,opt)
% Detect outliers using a spline smooth technique.
	handles = guidata(obj);					% get handles
	h = gmtedit_outliersdetect(handles.figure1, handles.axes1, handles.axes2, handles.axes3, ...
		[handles.h_gm handles.h_mm handles.h_tm]);
	filhas = getappdata(handles.figure1,'Filhas');
	setappdata(handles.figure1,'Filhas',[filhas(:); h])

% --------------------------------------------------------------------------------------------------
function ptcoords_clickedCB(obj,evt,opt)
% Get geog coordinates of nearest clicked point.
	handles = guidata(obj);					% get handles
	h = findobj(handles.figure1, 'type','line');
	if ( strcmp(get(obj,'State'),'on') )	% 
		set(handles.figure1,'WindowButtonDownFcn','', 'pointer','crosshair')
		set(h,'ButtonDownFcn', @bdn_ptcoords)
	else
		set(handles.figure1,'WindowButtonDownFcn',@add_MarkColor, 'pointer','arrow')
		set(h,'ButtonDownFcn', '')
	end

function bdn_ptcoords(obj, evt)
	handles = guidata(obj);					% get handles
	if (~ishandle(handles.hMirAxes))
		% User f up by killing the Mirone fig. Limit the damages
		set(handles.figure1,'WindowButtonDownFcn',@add_MarkColor,  'pointer','arrow')
		set(findobj(handles.figure1, 'type','line'),'ButtonDownFcn', '')
		% falta resetar o togglebutao
		return
	end
	ind = add_MarkColor(obj);				% Use this function to get the index of clicked pt
	if (isempty(ind)),	return,		end		% A bad click
	h = line(handles.lon(ind), handles.lat(ind),'Parent',handles.hMirAxes,'Marker','o', ...
		'MarkerFaceColor','k', 'MarkerSize',6,'LineStyle','none', 'Tag','LinkedSymb');
	setappdata(h,'box',handles.f_name)		% Used by draw_funs to append this file name to x,y when saving pts
	draw_funs(h,'DrawSymbol')				% Set uicontexts
% --------------------------------------------------------------------------------------------------
	
% --------------------------------------------------------------------
function scroll_plots(width, x)
%   This function uses the idea of the scrollplotdemo from Steven Lord (http://www.mathworks.com/matlabcentral)
%   scroll_plots(width,x);
%   width: window width in x units
%   x: absicssae vector

	pos = get(gca,'position');
	Newpos = [pos(1) 2 pos(3) 13];
	S = ['set(findall(gcf,''Type'',''axes''),''xlim'',get(gcbo,''value'')+[0 ' num2str(width) '])'];

	% Creating Uicontrol with initial value of the minimum of x
	uicontrol('style','slider','units','pixels','Pos',Newpos, 'Call',S,'min',x(1),'max',x(end)-width,'value',x(1));

% --------------------------------------------------------------------------
function varargout = gmtedit_outliersdetect(varargin)
% Do automatic outliers detection by comparison with spline sooth version data
 
	hObject = figure('Vis','off');
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

	set(hObject,'Vis','on');	drawnow
	
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
function edit_SmoothParam_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (xx < 0 || xx > 1 || isnan(xx))
		xx = 1;		set(hObject,'String',xx)
	end
	handles.smooth(handles.id_gmt) = xx;
	guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------
function edit_thresh_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (xx < 0 || isnan(xx))
		xx = 0;		set(hObject,'String',xx)
	end
	handles.thresh(handles.id_gmt) = xx;
	guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------
function radio_G_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_M handles.radio_T],'Val', 0)
	handles.id_gmt = 1;		guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function radio_M_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_G handles.radio_T],'Val', 0)
	handles.id_gmt = 2;		guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function radio_T_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_G handles.radio_M],'Val', 0)
	handles.id_gmt = 3;		guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function push_Apply_CB(hObject, handles)
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
function push_applyNreturn_CB(hObject, handles)
	push_Apply_CB(handles.push_Apply, handles)
	delete(handles.figure1)

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
function push_clear_CB(hObject, handles)
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
'Callback',@outliersdetect_CB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Enter a Smoothing Parameter between [0 1]',...
'Tag','edit_SmoothParam');

uicontrol('Parent',h1, 'Position',[9 48 125 15],...
'FontName','Helvetica',...
'String','Smoothing parameter (p)',...
'Style','text');

uicontrol('Parent',h1, 'Position',[212 31 90 21],...
'Callback',@outliersdetect_CB,...
'FontName','Helvetica',...
'FontSize',9,...
'String','Apply',...
'Tooltip','Use this for testing',...
'Tag','push_Apply');

uicontrol('Parent',h1, 'Position',[146 27 51 21],...
'BackgroundColor',[1 1 1],...
'Callback',@outliersdetect_CB,...
'String','4',...
'Style','edit',...
'Tag','edit_thresh');

uicontrol('Parent',h1, 'Position',[143 48 56 15],...
'FontName','Helvetica',...
'String','Threshold',...
'Style','text');

uicontrol('Parent',h1, 'Position',[11 6 30 15],...
'Callback',@outliersdetect_CB,...
'FontName','Helvetica',...
'FontSize',9,...
'String','G',...
'Style','radiobutton',...
'Tooltip','Select Gravity channel',...
'Tag','radio_G');

uicontrol('Parent',h1, 'Position',[60 6 30 15],...
'Callback',@outliersdetect_CB,...
'FontName','Helvetica',...
'FontSize',9,...
'String','M',...
'Style','radiobutton',...
'Tooltip','Select Magnetic channel',...
'Tag','radio_M');

uicontrol('Parent',h1, 'Position',[108 6 30 15],...
'Callback',@outliersdetect_CB,...
'FontName','Helvetica',...
'FontSize',9,...
'String','T',...
'Style','radiobutton',...
'Tooltip','Select Topography channel',...
'Tag','radio_T');

uicontrol('Parent',h1, 'Position',[212 7 90 21],...
'Callback',@outliersdetect_CB,...
'FontName','Helvetica',...
'FontSize',9,...
'String','Apply n return',...
'Tooltip','Do it and go away',...
'Tag','push_applyNreturn');

uicontrol('Parent',h1, 'Position',[147 0 50 21],...
'Callback',@outliersdetect_CB,...
'FontName','Helvetica',...
'String','Clear',...
'Tooltip','Clear detections from current selected channel',...
'Tag','push_clear');

function outliersdetect_CB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

% -------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------
function h = gmtedit_track(varargin)
% Little figure that asks for the name of grid to be interpolated and profile added to main fig
	h = figure('Vis','off');
	gmtedit_track_LayoutFcn(h);
	handles = guihandles(h);
 
	handles.lon = varargin{1};
	handles.lat = varargin{2};
	handles.year = varargin{3};
	handles.hCallingAxes = varargin{4};
	handles.h_tm = varargin{5};			% Handle of the Topography line
	handles.last_dir = varargin{6};
	handles.work_dir = varargin{7};
	handles.hMirAxes = varargin{8};
	handles.fname = [];

	if (ishandle(handles.hMirAxes))		% See if we have a previously used grid. If yes repeat its offer
		fname = getappdata(handles.hMirAxes, 'lastGrid');
		if (~isempty(fname))
			if (exist(fname,'file'))
				set(handles.edit_gridToInterp, 'Str', fname)
				handles.fname = fname;
			else
				rmappdata(handles.hMirAxes, 'lastGrid')
			end
		else
			set(handles.edit_gridToInterp, 'Str', handles.last_dir)
		end
	end

	guidata(h, handles);
	set(h,'Visible','on');

% ----------------------------------------------------------------------------
function edit_gridToInterp_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),		return,		end
	% Let the push_gridToInterp_CB do all the work
	push_gridToInterp_CB(handles.push_gridToInterp, handles, fname)

% ----------------------------------------------------------------------------
function push_gridToInterp_CB(hObject, handles, opt)
	if (nargin == 2),		opt = [];		end
	if (nargin == 3),		fname = opt;	end
	handles = guidata(hObject);		% Get udated handles

	if (isempty(opt))       % Otherwise 'opt' already transmited the file name.
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.grd;*.GRD', 'Grid files (*.grd,*.GRD,*.nc)';'*.*', 'All Files (*.*)'},'Select GMT grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
		set(handles.edit_gridToInterp,'String',fname)
	end

	handles.fname = fname;
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function zz = push_OK_CB(hObject, handles)
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
		h = findobj('Parent',handles.hCallingAxes, 'Type','line','Tag','orig_grav');
	elseif (strcmp(get(handles.hCallingAxes,'Tag'), 'axes2'))
		h = findobj('Parent',handles.hCallingAxes, 'Type','line','Tag','orig_mag');
	else
		h = findobj('Parent',handles.hCallingAxes, 'Type','line','Tag','orig_topo');
	end
	x = get(h,'XData');
	line('XData',x,'YData',zz,'Color','b','Parent',handles.hCallingAxes, 'LineWidth', 2, 'HitTest', 'off');
	setappdata(handles.hMirAxes, 'lastGrid', handles.fname)		% Store it so we can restart by it
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
'Callback',@gmtedit_track_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Name (full name) of grid to sample along track coords',...
'Tag','edit_gridToInterp');

uicontrol('Parent',h1, 'Position',[319 38 23 23],...
'Callback',@gmtedit_track_uiCB,...
'FontSize',12,...
'FontWeight','bold',...
'String','...',...
'Tag','push_gridToInterp');

uicontrol('Parent',h1, 'Position',[10 15 70 15],...
'FontName','Helvetica',...
'String','Add IGRF',...
'Style','checkbox',...
'Tooltip','For magnetic anomalies only, add an IGRF to interpolation',...
'Tag','check_addIGRF');

uicontrol('Parent',h1, 'Position',[276 9 66 23],...
'Callback',@gmtedit_track_uiCB,...
'FontName','Helvetica',...
'FontSize',9,...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1, 'Position',[10 62 180 16],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Grid to sample along track coords',...
'Style','text');

function gmtedit_track_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

% -------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------
function h = gmtedit_NavFilters(varargin)
 % Window with params for Nav/grads filtering
 
	h = figure('Vis','off');
	NavFilters_LayoutFcn(h);
	handles = guihandles(h);
 
	handles.hAxes = varargin{1};
	handles.hMag  = varargin{2};
	handles.fname = varargin{3};
	handles.minSpeed = 1;
	handles.maxSpeed = 15;
	handles.maxSlope = 250;

	guidata(h, handles);
	set(h,'Visible','on');


% ----------------------------------------------------------------------------
function edit_MaxSpeed_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (xx < 0 || isnan(xx))
		xx = 15;		set(hObject,'String',xx)
	end
	handles.maxSpeed = xx;
	guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------
function edit_MinSpeed_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (xx < 0 || isnan(xx))
		xx = 0;		set(hObject,'String',xx)
	end
	handles.minSpeed = xx;
	guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------
function edit_maxSlope_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (xx < 0 || isnan(xx))
		xx = 250;		set(hObject,'String',xx)
	end
	handles.maxSlope = xx;
	guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------
function push_navFiltClean_CB(hObject, handles)
% Remove all eventually detected bad Nav/grads from Mag channel
	hM = findobj(handles.hAxes,'Type','Line','tag','MagNull');
	set(hM, 'XData',[], 'YData',[])
	oldTit = get(handles.figure1, 'Name');
	set(handles.figure1, 'Name', oldTit(1:22))

% ----------------------------------------------------------------------------
function push_navFiltApply_CB(hObject, handles)
% Detect bad Nav and excessive gradients in data

	x = get(handles.hMag,'XData');			y = get(handles.hMag,'YData');
	
	try			tempo = nc_funs('varget', handles.fname, 'time');
	catch,		errordlg('No Time in this file (???)','Error'),		return
	end
	
	vel = diff(x) ./ diff(tempo') * (1000 / 1852 * 3600);		% x -> km; tempo -> sec. Now vel -> knots
	ind = (vel < handles.minSpeed | vel > handles.maxSpeed);
	indNaN = ~isnan(y(2:end));
	ind = ind & indNaN;			% We don't want the numbers reflect bad nav on ... non existing data
	ind = ind | ( abs(diff(y) ./ diff(x)) > handles.maxSlope);
	
	ind = find(ind);
	howmany = numel(ind);
	oldTit = get(handles.figure1, 'Name');
	set(handles.figure1, 'Name', sprintf('%s (found %d)', oldTit(1:22), howmany))
	if (~howmany),		return,		end			% Nothing

	ind = ind + 1;
	xx = x(ind);		yy = y(ind);

	hM = findobj(handles.hAxes,'Type','Line','tag','MagNull');
	if (isempty(hM))
		line(xx, yy,'Parent',handles.hAxes,'Marker','s','MarkerFaceColor','r','MarkerSize',4,'LineStyle','none','Tag','MagNull');
	else
		set(hM, 'XData',xx, 'YData',yy)
	end

% ----------------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function NavFilters_LayoutFcn(h1)

set(h1, 'Position',[520 500 356 61],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Speed and Slope filter', 'NumberTitle','off',...
'Resize','off', 'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[5 32 38 22], 'Style','edit',...
'BackgroundColor',[1 1 1],...
'Call',@navFilters_uiCB,...
'String','15',...
'Tooltip','Flag speeds higher than this',...
'Tag','edit_MaxSpeed');

uicontrol('Parent',h1, 'HorizontalAlignment','left', 'Position',[44 31 95 20],...
'String','Max speed (knots)','Style','text')

uicontrol('Parent',h1, 'Position',[151 31 38 22], 'Style','edit',...
'BackgroundColor',[1 1 1],...
'Call',@navFilters_uiCB,...
'String','1',...
'Tag','edit_MinSpeed');

uicontrol('Parent',h1, 'HorizontalAlignment','left', 'Position',[191 31 60 19],'String','Min speed','Style','text')

uicontrol('Parent',h1, 'Position',[5 4 38 22], 'Style','edit',...
'BackgroundColor',[1 1 1],...
'Call',@navFilters_uiCB,...
'String','250',...
'Tag','edit_maxSlope');

uicontrol('Parent',h1, 'HorizontalAlignment','left', 'Position',[44 3 109 20],...
'String','Max Slope (nT/km)','Style','text')

uicontrol('Parent',h1, 'Position',[260 32 90 21],...
'Call',@navFilters_uiCB,...
'String','Clean', 'Tag','push_navFiltClean')

uicontrol('Parent',h1, 'Position',[260 5 90 21],...
'Call',@navFilters_uiCB,...
'String','Apply', 'Tag','push_navFiltApply')

function navFilters_uiCB(hObject, event)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

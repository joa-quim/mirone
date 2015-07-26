function varargout = mag_synthetic(varargin)
% Extract a chunk of a mgd77+ file dictated from a 'guidingline' drawn interactively
% WARNING:	When the user selects to compute a RTP as well, this is be the curve
%			that will be used when doing correlation Fit in 'ecran'
%
% This function has the magmodel repeated from 'ecran', which is bad practice

%	Copyright (c) 2004-2012 by J. Luis
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

% $Id$

	if (nargin ~= 4)
		errordlg('MAG_SYNTHETIC: Called with wrong number of args (4)'),	return
	end
	
	hObject = figure('Vis','off');
	mag_synthetic_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'left');

	handles.hCallingFig = varargin{1};
	handles.hTrack = varargin{2};
	handles.hGuideLine = varargin{3};
	handles.chunkAzim = varargin{4};
	f_name = getappdata(varargin{2},'FullName');
	handles.f_name = f_name;
	handles.handlesMir = guidata(handles.hCallingFig);
	handles.last_dir = handles.handlesMir.last_dir;
	handles.work_dir = handles.handlesMir.work_dir;
	handles.home_dir = handles.handlesMir.home_dir;
	handles.d_path   = handles.handlesMir.path_data;
	handles.path_tmp = handles.handlesMir.path_tmp;

	% Set an uicontext for profile extraction of the mgd77++ file
	cmenuHand = uicontextmenu('Parent',handles.hCallingFig);
	set(handles.hGuideLine, 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Extract mgd77+ profile', 'Call', {@getProfile, hObject});
	uimenu(cmenuHand, 'Label', 'Delete me', 'Call', {@del_line, hObject}, 'Sep','on');
	ui_edit_polygon(handles.hGuideLine)		% Set edition functions

	handles.hZeroAge = [];
	handles.batProfile = [];
	handles.hEcran = [];
	handles.hSynthet = [];
	handles.spreadDir = '';
	handles.decl = '';
	handles.inc = '';
	handles.ageInit = '';
	handles.depthZero = 3;
	handles.contamin = 1;
	handles.speed = 2.5;
	handles.spreadObliqDir = 0;

	% See if we have a init file from a previous run to make our life a bit easier
	fid = fopen([handles.path_tmp 'mag_synt.ini'],'r');
	if (fid > 0)
		try
			fgetl(fid);		% Jump first line (it's a comment line)
			t = fread(fid,'*char');
			[handles.speed handles.spreadDir handles.spreadObliqDir handles.decl handles.inc handles.contamin] = ...
				strread(t,'%f %f %f %f %f %f');
			fclose(fid);
			set(handles.edit_speed, 'Str',handles.speed),	set(handles.edit_spreadDir, 'Str',handles.spreadDir)
			set(handles.edit_decl, 'Str',handles.decl),		set(handles.edit_spreadObliqDir, 'Str',handles.spreadObliqDir)
			set(handles.edit_inc, 'Str',handles.inc),		set(handles.edit_contamin, 'Str',handles.contamin)
		end
	end
	
	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, handles.text_synthet)
	%------------- END Pro look (3D) ------------------------------

	plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
	plugedWin = [plugedWin hObject];		% Add this figure handle to the carra?as list
	setappdata(handles.hCallingFig,'dependentFigs',plugedWin);

	set(hObject, 'Vis', 'on')

	% Choose default command line output for show_mag_mgd_export
	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end

% -----------------------------------------------------------------------------------------
function edit_speed_CB(hObject, handles)
	xx = sscanf(get(hObject,'Str'), '%f');
	if (isnan(xx)),		set(hObject,'Str',sprintf('%.2f',handles.speed)),	return,		end
	handles.speed = abs(xx);
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function edit_spreadDir_CB(hObject, handles)
	xx = sscanf(get(hObject,'Str'), '%f');
	if (isnan(xx)),		set(hObject,'Str',sprintf('%.1f',handles.spreadDir)),	return,		end
	handles.spreadDir = xx;
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function edit_spreadObliqDir_CB(hObject, handles)
% When the ridge opening is oblique
	xx = sscanf(get(hObject,'Str'), '%f');
	if (isnan(xx)),		set(hObject,'Str',sprintf('%.1f',handles.spreadObliqDir)),	return,		end
	handles.spreadObliqDir = xx;
	guidata(handles.figure1, handles)
	
% ----------------------------------------------------------------------------
function edit_ageInit_CB(hObject, handles)
% Age at the begining of the profile
	xx = sscanf(get(hObject,'Str'), '%f');
	if (isnan(xx)),		set(hObject,'Str',sprintf('%.1f',handles.ageInit)),	return,		end
	handles.ageInit = xx;
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function edit_depthZero_CB(hObject, handles)
	xx = sscanf(get(hObject,'Str'), '%f');
	if (isnan(xx)),		set(hObject,'Str',sprintf('%.1f',handles.depthZero)),	return,		end
	handles.depthZero = abs(xx);
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function edit_contamin_CB(hObject, handles)
	xx = sscanf(get(hObject,'Str'), '%f');
	if (isnan(xx) || xx > 1 || xx <= 0)
		set(hObject,'Str',sprintf('%.1f',handles.contamin)),	return
	end
	handles.contamin = abs(xx);
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function edit_decl_CB(hObject, handles)
	xx = sscanf(get(hObject,'Str'), '%f');
	if (isnan(xx)),		set(hObject,'Str',sprintf('%.1f',handles.decl)),	return,		end
	handles.decl = xx;
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function edit_inc_CB(hObject, handles)
	xx = sscanf(get(hObject,'Str'), '%f');
	if (isnan(xx)),		set(hObject,'Str',sprintf('%.1f',handles.inc)),		return,		end
	handles.inc = xx;
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function check_batGrid_CB(hObject, handles)
	onoff = {'off' 'on'};
	set([handles.edit_batGrid handles.push_batGrid],'Enable', onoff{get(hObject,'Val')+1})

% ----------------------------------------------------------------------------
function edit_batGrid_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),		return,		end
	% Let the push_batGrid_CB do all the work
	push_batGrid_CB(hObject,guidata(gcbo),fname)

% ----------------------------------------------------------------------------
function push_batGrid_CB(hObject, handles, opt)
% Get age grid name
	if (nargin == 3),	fname = opt;
	else				opt = [];
	end

	if (isempty(opt))           % Otherwise we already know fname from the 3th input argument
		str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
		[FileName,PathName] = put_or_get_file(handles,str1,'Select bathym grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	end
	set(handles.edit_batGrid,'String',fname)

% ----------------------------------------------------------------------------
function push_run_CB(hObject, handles)
% ...
	if (isempty(handles.decl))
		errordlg('Insult-declina','Error'),	return
	elseif (isempty(handles.inc))
		errordlg('Insult-inclina','Error'),	return
	elseif (isempty(handles.spreadDir))
		errordlg('Insult-spreadDir','Error'),	return
	elseif (isempty(handles.hZeroAge) && isempty(handles.ageInit))
		errordlg('Insult-zeroAge','Error'),	return
	elseif (~isempty(handles.hZeroAge) && ~ishandle(handles.hZeroAge) && isempty(handles.ageInit))
		errordlg('Who told you to kill the Zero Age marker? Bye-Bye.','Error'),	return
	end

	if ( isempty(handles.hEcran) || ~ishandle(handles.hEcran) )
		handles = getProfile([], [], handles.figure1);
	end

	if (ishandle(handles.hZeroAge))				% Get the Zero age location
		pt_x = get(handles.hZeroAge,'XData');		pt_y = get(handles.hZeroAge,'YData');
	else										% Nope, no Zero age location. Use begining of guiding line instead
		pt_x = get(handles.hGuideLine,'XData');		pt_y = get(handles.hGuideLine,'YData');
		pt_x = pt_x(1);		pt_y = pt_y(1);
	end
	ind_0 = find_nearpoint(handles.trk_chunk.lon, handles.trk_chunk.lat, pt_x, pt_y);

	dist = handles.trk_chunk.distance - handles.trk_chunk.distance(ind_0);
	if (get(handles.check_batGrid,'Val') && isempty(handles.batProfile))
		fname = get(handles.edit_batGrid,'str');
		if (isempty(fname))
			errordlg('Insult-bat','Error'),	return
		elseif (exist(fname,'file') ~= 2)
			errordlg('Bathymetry grid does not exist.','Error'),	return
		end
		[handles, X, Y, Z, head] = read_gmt_type_grids(handles,fname);
		if (isempty(X)),	return,		end
		x = get(handles.hGuideLine,'XData');		y = get(handles.hGuideLine,'YData');
		handles.batProfile = abs(double(grdtrack_m(Z,head,[x(:) y(:)],'-Z')) / 1000 );
		handles.trk_chunk.topography = handles.batProfile(:);
	elseif (isempty(handles.trk_chunk.topography))
		handles.trk_chunk.topography = ones(numel(dist),1)*2.5;
	end

	magRTP = [];
	if ( get(handles.check_RTP,'Val') )			% See if we need to compute the RTP profile as well
		[theta,ampfac] = nskew(0, handles.spreadDir-90, handles.decl, handles.inc, pt_x);
		magRTP = -rtp2d(handles.trk_chunk.mag, theta, ampfac);
	end

	reversalsFile = [handles.d_path 'Cande_Kent_95.dat'];
	handEcran = guidata(handles.hEcran);		% Need to access the handle of second axes of the 'Ecran' figure
	dxyz = [dist(:) handles.trk_chunk.lon(:) handles.trk_chunk.lat(:) handles.trk_chunk.topography(:)];

	if (isempty(handles.hSynthet) || ~ishandle(handles.hSynthet) )	% If model line does not exist yet
		a1 = abs(dist(1)) / (handles.speed / 2 * 10);				% Remember, half speed
		a2 = abs(dist(end)) / (handles.speed / 2 * 10);
		handEcran.ageEnd = a2;
		if (a1 > a2),		handEcran.ageEnd = a1;		end			% Pick the abs longer of both half lengths
		if (isempty(handles.hZeroAge) && ~isempty(handles.ageInit))	% No ZeroAge marker, but start age from edit box
			handEcran.ageStart = handles.ageInit;
		else
			handEcran.ageStart = 0;
		end

		% Get obliquity that takes into account the spreading obliquity for eventual Mag Bars stretch/shrink
		[dir_spread, obliq] = obliquity_care(handles.spreadDir, handles.chunkAzim, handles.spreadObliqDir);
		if (obliq > 0),			handEcran.stretchMagBar = 1 / cos(obliq * pi / 180);
		elseif (obliq < 0)		handEcran.stretchMagBar = cos(obliq * pi / 180);
		end
		guidata(handEcran.figure1, handEcran)						% Update Ecran handles so it knows hot to plot the mag bars
		cb = get(handEcran.push_magBar,'Call');						% Clever trick to fish the callback
		feval(cb, handEcran.push_magBar, [])						% Call the ecran_uiCB function thats plots the mag bars
		set(handEcran.push_syntheticRTP, 'Vis','off')				% 'Ecran' button that we don't want to be visible here
		set(handEcran.popup_selectSave, 'Vis','off')
		set(handEcran.push_ageFit, 'Vis','on')

		% -------------- Have to go over this again. 'Ecran' needs to know about these fields -----------
		handEcran = guidata(handles.hEcran);
		syntPar = iquire_OPTcontrol([handles.d_path 'OPTcontrol.txt']);
		handEcran.syntPar.agePad = 1.5;				% Probably too short for wide isochrons
		if (isempty(syntPar))
			set(handEcran.edit_ageFit, 'Vis','on')
		else
			handEcran.syntPar.ageMarkers = syntPar.ageMarkers;
			handEcran.syntPar.agePad = syntPar.agePad;
			handEcran.syntPar.ageStretch = syntPar.ageStretch;
			set( handEcran.popup_ageFit,'Vis','on','Str',handEcran.syntPar.ageMarkers )
		end
		[anoma, handEcran.age_line] = magmodel(handEcran.axes2, reversalsFile, dxyz, handles.decl, handles.inc, ...
			handles.speed, handles.spreadDir, handles.chunkAzim, handles.spreadObliqDir, handles.contamin);
		handles.hSynthet = line('XData',dist-dist(1), 'YData',anoma, 'Parent',handEcran.axes1);
		if (~isempty(magRTP))
			handles.hRTP = line('XData',dist-dist(1), 'YData',magRTP, 'Color','b', 'LineStyle','--','Parent',handEcran.axes1);
			y_lim = get(handEcran.axes1, 'YLim');	% Update YLim because magRTP can go out of old limits.
			y_lim(1) = min(y_lim(1), min(magRTP));		y_lim(2) = max(y_lim(2), max(magRTP));
			set(handEcran.axes1, 'YLim', y_lim);
			handEcran.hLine = handles.hRTP;			% From now one this the current line in Ecran ... which may give troubles
		end
		handEcran.hSynthetic = handles.hSynthet;
		if (obliq ~= 0),	handEcran.age_line = handEcran.age_line / handEcran.stretchMagBar;	end		% Yes, divide.
		guidata(handEcran.figure1, handEcran)		% handEcran.age_line is needed in Ecran's edit_ageFit_CB()
		% -------------------------------------------------------------------------------------------------

		y_lim = get(handEcran.axes1,'ylim');
		set(handEcran.axes1,'ylim',[min(y_lim(1),min(anoma)) max(y_lim(2),max(anoma))])

		uistack_j(handEcran.axes2,'up')				% Bring it up relatively to axes1 so that we can see the ticks
		ages = get(handEcran.axes2, 'UserData');
		if (isempty(ages))
			errordlg('This time you exagerated with the experiments. Something screw the Mag Bars. Bye.','Error'),	return
		end
		hT = findobj(handEcran.axes1,'Tag','chroneName');
		set(hT,'Vis','off')
		if (ages(1) == 0),		ages(1) = [];	hT(end) = [];	end	% We don't want the Zero tick
		ages = [-ages(end:-1:1) ages];								% Will be new ticks
		set(handEcran.axes2,'XTick', ages(:))
		str = get(hT,'Str');
		set(handEcran.axes2,'XTickLabel',[str; str(end:-1:1)]);		% And finally, the new tickmarks

	else				% Line exists. Update it
		anoma = magmodel(handEcran.axes2, reversalsFile, dxyz, handles.decl, handles.inc, handles.speed, ...
			handles.spreadDir, handles.chunkAzim, handles.spreadObliqDir, handles.contamin);
		set(handles.hSynthet, 'XData',dist-dist(1), 'YData',anoma)
	end

	guidata(handles.figure1, handles);
% ---------------------------------------------------------------------
function syntPar = iquire_OPTcontrol(opt_file)
% See if the opt_file (OPTcontrol.txt) file has relevant info for this run
% If yes, the parametrs will be stored in the SYNTPAR struct. Otherwise the struct is empty

	syntPar = [];		% The default value
	if ( exist(opt_file, 'file') == 2 )
		fid = fopen(opt_file, 'r');
		c = (fread(fid,'*char'))';      fclose(fid);
		lines = strread(c,'%s','delimiter','\n');   clear c fid;
		m = numel(lines);		kk = 2;		% kk starts at 2 because first line in popup is empty
		for (k = 1:m)
			if (~strncmp(lines{k},'MIR_MAGPROF',7)),	continue,	end
			if (numel(lines{k}) <= 14),	continue,	end		% The minimum it takes to have a -? switch
			[t, r] = strtok(lines{k}(13:end));
			switch t
				case 'ISOC'			% The list of isochrons (by its age) that will be tentatively fit
					ind = strfind(r, '_');
					if (~isempty(ind))
						if (numel(ind) == 2)
							syntPar.ageStretch(kk) = fix(str2double(r(ind(2)+1:end)));	% Used to expand/shrink in corr
							r(ind(2):end) = [];		% Wipe it out so the "...agePad(kk)" assignement works for any 'ind'
						else
							syntPar.ageStretch(kk) = 0;
						end
						syntPar.ageMarkers{kk} = r(1:ind(1)-1);
						syntPar.agePad(kk) = str2double(r(ind(1)+1:end));
					else
						syntPar.ageMarkers{kk} = r;
						syntPar.agePad(kk) = 1.5;		% Default (and possibly too short for old isocs) value
					end
					kk = kk + 1;
			end
		end
	end

% ----------------------------------------------------------------------------
function handles = getProfile(obj, evt, hFig)
% Get a chunk of the mgd77++ profile roughly deffined by the guideLine
	handles = guidata(hFig);
	trk_chunk = get_track_chunk(handles);
	if (~isempty(trk_chunk))
		bak = handles.handlesMir.DefineMeasureUnit;		% Save this for restoring
		handles.handlesMir.DefineMeasureUnit = 'k';
		handles.hEcran = ecran(handles.handlesMir,trk_chunk.lon,trk_chunk.lat,trk_chunk.mag,['Track from ' handles.f_name]);
		handles.trk_chunk = trk_chunk;
		handles.handlesMir.DefineMeasureUnit = bak;		% Reset orig value
		set(handles.text_gotDepth, 'Vis', 'on')
	end

% ----------------------------------------------------------------------------
function del_line(obj, evt, hFig)
% Delete the guiding line and this figure too. OBJ is a child of Mirone fig
	handles = guidata(hFig);
	delete([handles.hGuideLine handles.figure1])

% ----------------------------------------------------------------------------
function push_zeroAge_CB(hObject, handles)
% Plot a star on the Mirone fig indicating the position of the Zero age
	figure(handles.hCallingFig)
	handles.hZeroAge = mirone('Draw_CB',guidata(handles.hCallingFig),'Symbol','p');
	set(handles.hZeroAge,'MarkerFaceColor','b','Tag','ZeroAgePt');
	% Use the time of first point to compute IGRF dec and dip
	lon = get(handles.hZeroAge,'XData');	lat = get(handles.hZeroAge,'YData');
	if ( isempty(get(handles.edit_decl, 'Str')) && isempty(get(handles.edit_inc, 'Str')) )
		% Only fill with IGRF values if User didn't go there first
		out = igrf_m(lon, lat, 0, 1990, '-Fdi');
		handles.decl = out(1);		handles.inc = out(2);
		set(handles.edit_decl, 'Str', out(1)),		set(handles.edit_inc, 'Str', out(2))
	end
	guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function trk_chunk = get_track_chunk(handles)
% Working function called by getProfile

	[PATH,FNAME,EXT] = fileparts(handles.f_name);
	if (isempty(PATH)),		PATH = cd;		end		% Play safe
	if (strcmpi(EXT,'.nc'))
		handles.f_name = [PATH filesep FNAME '.nc'];		handles.is_mgd77 = true;		handles.is_gmt = false;
		[handles, track] = read_mgd77_plus(handles, handles.f_name);
	else
		handles.f_name = [FNAME '.gmt'];	handles.is_gmt = true;			handles.is_mgd77 = false;
		track = c_gmtlist([PATH filesep FNAME], '-Fsxygmtd', handles.opt_G);
	end

	if ~( numel(track.magnetics) > 1 && ~all(isnan(track.magnetics)) )
		trk_chunk = [];
		return
	end

	% Start the display at a user selected coordinate
	pt_hGuideLine = [get(handles.hGuideLine,'XData')' get(handles.hGuideLine,'YData')'];
	if (size(pt_hGuideLine,1) > 2)		% We don't want more than two vertex
		pt_hGuideLine = [pt_hGuideLine(1,:); pt_hGuideLine(end,:)];
	end
	id = zeros(1,2);
	for (k = 1:2)
		id(k) = find_nearpoint(track.lon, track.lat, pt_hGuideLine(k,1), pt_hGuideLine(k,2));
	end

	step = 1;
	if (id(1) > id(2)),		step = -1;		end	% 50% likely cases
	ind = id(1):step:id(2);

	trk_chunk.lon = track.lon(ind);
	trk_chunk.lat = track.lat(ind);
	trk_chunk.mag = double(track.magnetics(ind));
	if (step == 1)				% Picked track runs in the same sense as data in file
		trk_chunk.distance = track.distance(ind) - track.distance(ind(1));		% distances start at 0
	else						% Track and data in file run in oposit sense. Adjust acumulated dist accordingly.
		t = track.distance(ind);
		trk_chunk.distance = [0; cumsum( diff(t(end:-1:1)) )];
	end

	guideLineLength = vdist(pt_hGuideLine(1,2), pt_hGuideLine(1,1), pt_hGuideLine(2,2), pt_hGuideLine(2,1)) / 1000;
	if ( (trk_chunk.distance(end) - trk_chunk.distance(1)) > guideLineLength * 1.2)
		warndlg('Extracted chunk suspiciously longer than guide line. Picked on a track cross-over?','Warning')
	end

	indNan = isnan(track.topography);
	if ( numel(track.topography) > 1 && ~all(indNan) )
		trk_chunk.topography = double(track.topography(ind)) / 1000;
		indNan = isnan(trk_chunk.topography);
		if (any(indNan))		% Shit. Let's hope they are not too many
			xx = trk_chunk.distance;
			trk_chunk.topography(indNan) = interp1(xx(~indNan), trk_chunk.topography(~indNan), xx(indNan), 'linear', 'extrap');
		end
	else
		trk_chunk.topography = [];
	end

	% Reinterpolate at a constant step because that may be needed for the FFT (RTP)
	x = linspace(0, trk_chunk.distance(end), numel(trk_chunk.distance))';
	trk_chunk.mag = interp1(trk_chunk.distance, trk_chunk.mag, x);
	trk_chunk.distance = x;
	% We should do the topography too, but I think it doesn't worth it because the apprximation is good enough

% --------------------------------------------------------------------
function id = find_nearpoint(vec_x, vec_y, pt_x, pt_y)
% Find the nearest point on vector [vec_x vec_y] to given point [pt_x pt_y] and return its index
	x = vec_x - pt_x;		y = vec_y - pt_y;
	[zzz, id] = min(x.*x + y.*y);

% --------------------------------------------------------------------
function [handles, track] = read_mgd77_plus(handles, fname)
	s = nc_funs('info',fname);

	% ------------------ OK, Get numerics now -----------------------------------
	track.lon  = double(nc_funs('varget', fname, 'lon'));
	track.lat  = double(nc_funs('varget', fname, 'lat'));
	track.magnetics  = nc_funs('varget', fname, 'anom_cm4');
	track.topography = nc_funs('varget', fname, 'depth');
	tempo = nc_funs('varget', fname, 'time');
	track.time = tempo(1);		% We will use this to compute IGRF

	D2R = pi / 180;		KMPRDEG = 111.1949;
	co = cos(track.lat * D2R);
	dx = [0; diff(track.lon)] .* co;		dy = [0; diff(track.lat)];
	dist = cumsum(sqrt(dx .^2 + dy .^2) * KMPRDEG);

	track.distance = dist;

	% -------------- Get the nodata-values ---------------------------------------
	ind = strcmp({s.Dataset.Name}, 'mtf1');
	id = strcmp({s.Dataset(ind).Attribute.Name}, 'missing_value');
	handles.magNoValue = s.Dataset(ind).Attribute(id).Value;
	id = strcmp({s.Dataset(ind).Attribute.Name}, 'scale_factor');
	handles.magScaleF = s.Dataset(ind).Attribute(id).Value;
	
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

% ----------------------------------------------------------------------------
function fout = rtp2d(fin, theta, ampfac) 
% RTP2D - Reduce magnetic anomaly profile to the pole
% given the phase angle (theta) of the observed profile
% see nskew for calculating phase angle
%
%  fin   : input profile
%  theta : skew angle (degrees) of observed profile
% ampfac : amplitude factor
%   fout : deskewed output profile
%
% Maurice A. Tivey Feb 8 1993 MATLAB

	nn = numel(fin);
	nn2 = nn/2;
	ind = isnan(fin);
	if ( any(ind) )			% Sometimes we have NaNs in mag. If yes, fill them
		x = 1:nn;
		fin(ind) = interp1(x(~ind),fin(~ind),x(ind),'linear','extrap');
	end
	meanfld = mean(fin);
	M = fft(fin-meanfld);
	%------------------ PHASE filter  --------------
	phase = exp(1i * theta * pi/180 .* sign(nn2-(1:nn)));
	%------------------ INVERSE fft ----------------
	fout = ifft(M ./ (ampfac .* phase)');
	fout = real(fout) + meanfld;		% Add back the mean value to result
	if ( any(ind) )						% If we had NaNs, put them back
		fout(ind) = NaN;
	end

% ----------------------------------------------------------------------------
function [theta,ampfac] = nskew(zobs, slin, fdec, fdip, rlat, sdec, sdip)
% NSKEW - Compute skewness parameter and amplitude factor following Schouten (1971)
%  Computes GEOCENTRIC DIPOLE unless given declination and dip of magnetization      
% Usage: 
%    [theta,ampfac]=nskew(zobs, slin, fdec, fdip)
% or
%    [theta,ampfac]=nskew(zobs, slin, fdec, fdip, rlat, sdec, sdip)
%
%  Input variables:
%   rlat : regional latitude in decimal degrees
%   zobs : level of observation in km above sealevel
%   slin : strike of lineations normal to profile (+cw degrees from north)
%   fdec,fdip : field declination,inclination 
%   sdec,sdip : magnetization declination,inclination 
% Output variables
%   theta : phase angle
%   ampfac : amplitude factor
%
% Maurice A. Tivey February 3, 1993
%    checked April 1996
%
% MIRONEFIED 24-Jan-2011	J. Luis
%---------------------------------------------------------

	D2R = pi/180;
	if (nargin == 5)
		%  NOTE FOR GEOCENTRIC DIPOLE TAN(INC)=2*TAN(LAT)
		sdip = atan2( 2.*sin(rlat*D2R),cos(rlat*D2R) ) / D2R;
		sdec = 0;
	end
	% compute phase and amplitude factors
	ra1 = fdip * D2R;
	rb1 = (fdec-slin) * D2R;
	ra2 = sdip * D2R;
	rb2 = (sdec-slin) * D2R;
	% compute phase and amplitude factors
	inclm  = atan2(tan(ra2),sin(rb2));
	inclf  = atan2(tan(ra1),sin(rb1));
	ampfac = ((sin(ra2))*(sin(ra1))) / ((sin(inclm))*(sin(inclf)));
	theta  = (inclm/D2R) + (inclf/D2R) - 180;
	if (theta <= -360),		theta = theta + 360;
	elseif (theta >= 360),	theta = theta - 360;
	end

% ------------------------------------------------------------------------------
function [anoma, age_line, obliquity] = magmodel(hAxesMagBar, reversalsFile, dxyp, fDec, fInc, speed, ...
	dir_spread, dir_profil, spreadObliqDir, contam)
% hAxesMagBar	Must contain the handle of axes where to plot the Mag Bars (axes2 of the 'Ecran' fig)
% spreadObliqDir	Is the spreaging obliquity (often, but not always, is zero)
% contam		The Tissot & Patriat's contamination factor
%
% NOTE: dir_spread can be set via OPTcontrol.txt, if it was not than it defaults to equal dir_profil
%
% This function was adapted from the MAGMOD program
% "MODMAG, a MATLAB program to model marine magnetic anomalies."
% Véronique Mendel, Marc Munschy and Daniel Sauter. Computers & Geosciences 31 (2005) 589–597
% Intention is to clean it more or eventually re-write it (for example it doesn't account
% for possibly different Remanent & Inducted mags).

	if (nargin < 10),	contam = 1;		end
	if (size(dxyp,2) == 3)		% If depth is not provided, default to const = 2.5 km
		dxyp = [dxyp(:,1:3) ones(size(dxyp,1),1)*2.5];
	end

	zObs = 0;		magAtAxe = 18;		magFlatOfEachBlock = 4;
	psubsi = 0.35;			% Coefficient of thermal subsidence equation
	thickness = 0.5;
	speed = speed * 10;		% From full rate in cm/year to rate in km/Ma 

	fid = fopen(reversalsFile,'r');
	todos = fread(fid,'*char');
	[chron age_start age_end s.age_txt] = strread(todos,'%s %f %f %s');
	fclose(fid);    clear todos
	BA = zeros(numel(age_start)*2,2);
	BA(1:2:end-1,1) = age_start;	BA(2:2:end,1) = age_end;
	BA(1:2:end-1,2) = age_end;
	BA(2:2:end-2,2) = age_start(2:end);
	BA(end,:) = [];

	distAlongProfile = dxyp(:,1);
	profileDepth = dxyp(:,4);
	ind = isnan(profileDepth);
	if ( any(ind) )			% Many files have gaps in bathymetry. Reinvent it
		profileDepth(ind) = interp1(distAlongProfile(~ind), profileDepth(~ind), distAlongProfile(ind));
	end
	
	nBlocks = size(BA,1);
	nPts = numel(distAlongProfile);		% Number of points
	stat_z = zeros(nPts,1) + zObs;		% Depth of the points where the magnetic anomaly will be compute

	% Determination of the age limits of normal and inverse polarity bodies and
	% of their respecting magnetization

	twiceNBlocks = (nBlocks*2) - 1;

	blockAge = zeros(twiceNBlocks,2);
	blockMag = zeros(twiceNBlocks,1);
	blockAge(1:nBlocks-1,1) = -BA(nBlocks:-1:1+1,2);
	blockAge(1:nBlocks-1,2) = -BA(nBlocks:-1:1+1,1);
	blockAge(nBlocks,:) = [-BA(1,2) BA(1,2)];
	blockAge(nBlocks+1:twiceNBlocks,:) = BA(1+1:nBlocks,:);
	blockMag(nBlocks) = magAtAxe;
	blockMag(nBlocks-2:-2:1) = magFlatOfEachBlock;
	blockMag(nBlocks+2:2:twiceNBlocks) = magFlatOfEachBlock;
	blockMag(nBlocks-1:-2:1) = -magFlatOfEachBlock;
	blockMag(nBlocks+1:2:twiceNBlocks) = -magFlatOfEachBlock;

	% Calculation of magnetized bodies position in km
	polygon_x = zeros(4,twiceNBlocks);
	polygon_x(1,:) = blockAge(:,1)' * speed / 2;
	polygon_x(2,:) = blockAge(:,2)' * speed / 2;
	polygon_x(3,:) = polygon_x(2,:);
	polygon_x(4,:) = polygon_x(1,:);

	% Before calculation of the magnetic anomaly, the spreading direction has
	% to be less than 90° away from the profile direction in order to be plot
	% in the right sense of distance. To do this we first calculate the obliquity
	% between the direction of the profile and the spreading direction.
	[dir_spread, obliquity] = obliquity_care(dir_spread, dir_profil, spreadObliqDir);
	if (obliquity ~= 0)
		if (obliquity > 0)
			polygon_x = polygon_x / cos (obliquity * pi / 180);
		else
			% Rarer cases where ridge spreading is oblique and the profile lies inside the obliquity cone
			polygon_x = polygon_x * cos (obliquity * pi / 180);
		end
	end

	hMagBar = findobj(hAxesMagBar, 'type', 'patch');
	xx = get(hMagBar,'XData');
	ind = find((xx(1) - blockAge(:,1)) > 0);	% Find index of closest block start of displyed bricks and those from file
	ind = ind(end);								% The last one is closest from the left (starting) side
	f = (xx(1) - blockAge(ind,1)) / (blockAge(ind,2) - blockAge(ind,1));
	distAlongProfile = distAlongProfile + polygon_x(1,ind) + f * (polygon_x(2,ind) - polygon_x(1,ind));

	% Calculation of the magnetized bodies depth
	polygon_z = zeros(4,twiceNBlocks);

	ind = find(polygon_x(1,:) < distAlongProfile(1));	% Find start block limits to the left of profile start
	if ~isempty(ind)
		maxlignm11 = max(ind);
		if (polygon_x(1,maxlignm11) < 0)
			polygon_z(1,ind) = profileDepth(1) + psubsi*sqrt(blockAge(maxlignm11,1)-blockAge(ind,1)');
		else
			polygon_z(1,ind) = profileDepth(1);
		end
		polygon_z(4,ind) = polygon_z(1,ind) + thickness;
	end

	ind = find(polygon_x(2,:) < distAlongProfile(1));	% Find end block limits to the left of profile start
	if ~isempty(ind)
		maxlignm12 = max(ind);
		if polygon_x(2,maxlignm12) < 0
			polygon_z(2,ind) = profileDepth(1) + psubsi*sqrt(blockAge(maxlignm12,2) - blockAge(ind,2)');
		else
			polygon_z(2,ind) = profileDepth(1);
		end
		polygon_z(3,ind) = polygon_z(2,ind) + thickness;
	else
		maxlignm12 = 0;
	end

	ind = find(polygon_x(1,:) >= distAlongProfile(end));	% Find start block limits to the right of profile end
	if ~isempty(ind)
		minlignm21 = min(ind);
		if (polygon_x(1,minlignm21) > 0)
			polygon_z(1,ind) = profileDepth(end) + psubsi*sqrt(blockAge(ind,1)' - blockAge(minlignm21,1));	
		else
			polygon_z(1,ind) = profileDepth(end);
		end
		polygon_z(4,ind) = polygon_z(1,ind) + thickness;
	else
		minlignm21 = twiceNBlocks + 1;
	end

	ind = find(polygon_x(2,:) >= distAlongProfile(end));	% Find end block limits to the right of profile end
	if ~isempty(ind)
		minlignm22 = min(ind);
		if polygon_x(2,minlignm22) > 0
			polygon_z(2,ind) = profileDepth(end) + psubsi*sqrt(blockAge(ind(:),2)' - blockAge(minlignm22,2));	
		else
			polygon_z(2,ind) = profileDepth(end);
		end
		polygon_z(3,ind) = polygon_z(2,ind) + thickness;
	end

	% Now find the block limits that are contained inside the profile
	ind = find(polygon_x(1,:) >= distAlongProfile(1) & polygon_x(1,:) < distAlongProfile(end));
	if ~isempty(ind)
		polygon_z(1,ind) = interp1(distAlongProfile, profileDepth, polygon_x(1,ind));
		polygon_z(4,ind) = polygon_z(1,ind) + thickness;
	end
	ind = find(polygon_x(2,:) >= distAlongProfile(1) & polygon_x(2,:) < distAlongProfile(end));
	if ~isempty(ind)
		polygon_z(2,ind) = interp1(distAlongProfile, profileDepth, polygon_x(2,ind));
		polygon_z(3,ind) = polygon_z(2,ind) + thickness;
	end

	PolXX = cell(1,twiceNBlocks);
    PolZZ = cell(1,twiceNBlocks);

	for i = 1:twiceNBlocks
		PolXX{i} = [polygon_x(:,i); polygon_x(1,i)];
		PolZZ{i} = [polygon_z(:,i); polygon_z(1,i)];
	end

	for i = maxlignm12+1:minlignm21-1
		ptdist = find(distAlongProfile > polygon_x(1,i)+0.00001 & distAlongProfile < polygon_x(2,i)-0.00001);
		if ~isempty(ptdist)
			tempox1 = distAlongProfile(ptdist);
			tempoz1 = profileDepth(ptdist);
			tempox2 = flipud(tempox1);
			tempoz2 = flipud(tempoz1+thickness);
			tempox3 = [polygon_x(1,i);tempox1;polygon_x(2,i);polygon_x(3,i);tempox2;polygon_x(4,i)];
			tempoz3 = [polygon_z(1,i);tempoz1;polygon_z(2,i);polygon_z(3,i);tempoz2;polygon_z(4,i)];
			PolXX{i}= [tempox3; polygon_x(1,i)];
			PolZZ{i}= [tempoz3; polygon_z(1,i)];
		end
	end

	% Re-positionning of the magnetized bodies if the contamination coefficient is different from 1
	if (contam ~= 1)
		for (i = 1:twiceNBlocks)
			PolXX{i} = PolXX{i} * contam;
		end
		distAlongProfile = distAlongProfile * contam;
	end

	% Calculation of the magnetic anomaly created by the magnetized bodies
	anoma = calcmag(twiceNBlocks,fInc,fDec,dir_spread,blockMag,distAlongProfile,stat_z,nPts,PolXX,PolZZ);

	if (nargout >= 2)						% Calculate the ages. 
		age_line = distAlongProfile / speed * 2;
	end

% ---------------------------------------------------------------------------------------
function [dir_spread, obliquity] = obliquity_care(dir_spread, dir_profil, ridgeObliquity)
% Compute the obliquity between the direction of the profile and the opening direction
% but taking into account that the Ridge my open aslant as well.
% OBLIQUITY Is positive when DIR_PROFIL is outside the obliquity cone and negative otherwise.
% The 'obliquity cone' exists only for oblique spreading (ridgeObliquity ~= 0) and is deffined
% (by Me) as having an angle twice that of the spreading direction and the Ridge normal. That
% is, the 'ridgeObliquity', and whose axial line is the ridge normal.

	if (nargin < 3),	ridgeObliquity = 0;		end

	% First bring the profile and spreading directions to the same side of the 180 deg barrier.
	obliquity = dir_profil - dir_spread;
	if (abs(obliquity) > 90)
		temp1 = mod(obliquity,90);
		temp2 = mod(obliquity,-90);
		if (abs(temp1) <= abs(temp2)),		obliquity = temp1;
		else								obliquity = temp2;
		end
		dir_spread = dir_profil - obliquity;
	end
	obliquity = abs(obliquity);				% Remember, signal counts. Positive ==> 'outer' obliquity

	% End of the story?
	if (~ridgeObliquity),		return,		end		% We are done, bye

	if (dir_spread > 180)		% Let's call it the left side branch
		if ( (dir_profil <= dir_spread) || (dir_profil >= dir_spread + 2*ridgeObliquity) )
			obliquity = abs(dir_profil - dir_spread);	% Simpler case.
		else
			obliquity = -abs(ridgeObliquity - mod(dir_profil - dir_spread, ridgeObliquity));
		end
	else						% The right side branch
		if ( (dir_profil >= dir_spread + 2*ridgeObliquity) || (dir_profil <= dir_spread) )
			obliquity = abs(dir_profil - dir_spread);	% Simpler case.
		else
			obliquity = -abs(ridgeObliquity - mod(dir_profil - dir_spread, ridgeObliquity));
		end
	end

% ---------------------------------------------------------------------
function anoma = calcmag(nb_struct,fInc,fDec,dir_spread,blockMag,stat_x,stat_z,nPts,PolXX,PolZZ)
% Calculation of the magnetic anomaly created by magnetized bodies

	D2R = pi / 180;
	magInc(1:nb_struct,1) = fInc;
	magDec(1:nb_struct,1) = fDec;
    dir_anomag = dir_spread + 90;

	c1 = sin(fInc*D2R);    
	c2 = cos(fInc*D2R) * cos(dir_spread*D2R - fDec*D2R);

	d1 = sin(magInc*D2R);
	d2 = cos(magInc*D2R).*cos((dir_anomag-90)*D2R - magDec*D2R);
	d3 = 200 * blockMag;

	anomax = 0;		anomaz = 0;

	for (k = 1:nb_struct)
		n = numel(PolXX{k});
		[amx,amz] = fcalcmagpt(n,stat_x,stat_z,nPts,PolXX{k},PolZZ{k},d1(k),d2(k),d3(k));
		anomax = anomax + amx;
		anomaz = anomaz + amz;
	end

	anoma  = c2*anomax + c1*anomaz;

% --------------------------------------------------------------------
function [amxx,amzz] = fcalcmagpt(nbps,stax,staz,nbsta,polxxm,polzzm,dd1,dd2,dd3)
% Function to calculate the magnetized anomaly created
% by each point of a magnetized body

	matstax = repmat(stax,1,nbps-1);
	matstaz = repmat(staz,1,nbps-1);
	matpolxx = repmat(polxxm',nbsta,1);
	matpolzz = repmat(polzzm',nbsta,1);

	x1 = matpolxx(:,1:(nbps-1)) - matstax(:,1:(nbps-1));
	z1 = matpolzz(:,1:(nbps-1)) - matstaz(:,1:(nbps-1));
	x2 = matpolxx(:,2:nbps) - matstax(:,1:(nbps-1));
	z2 = matpolzz(:,2:nbps) - matstaz(:,1:(nbps-1));

	indx1 = find(x1==0);
	if (~isempty(indx1)),		x1(indx1) = 1e-11;		end
	indz1 = find(z1==0);
	if (~isempty(indz1)),		z1(indz1) = 1e-11;		end
	indx2 = find(x2==0);
	if (~isempty(indx2)),		x2(indx2) = 1e-11;		end
	indz2 = find(z2==0);
	if (~isempty(indz2)),		z2(indz2) = 1e-11;		end

	th1 = atan2(z1,x1);
	th2 = atan2(z2,x2);
	t12 = th1-th2;
	z21=z2-z1;
	x21=x2-x1;
	xz12=x1.*z2-x2.*z1;
	r21s=x21.*x21+z21.*z21;
	r1s = x1.^2+z1.^2;
	r2s = x2.^2+z2.^2;
	rln = 0.5*log(r2s./r1s);

	p=(xz12./r21s).*((x1.*x21-z1.*z21)./r1s-(x2.*x21-z2.*z21)./r2s);
	q=(xz12./r21s).*((x1.*z21+z1.*x21)./r1s-(x2.*z21+z2.*x21)./r2s);

	f1 = (t12.*z21-rln.*x21)./r21s;
	f2 = (t12.*x21 + rln.*z21)./r21s;

	dxx = p + z21.*f1;
	dxz = q - x21.*f1;
	dzz = -p + x21.*f2;
	dzx = q - z21.*f2;

	amxx = dd3*(dd1*sum(dxz,2)+dd2*sum(dxx,2));
	amzz = dd3*(dd1*sum(dzz,2)+dd2*sum(dzx,2));
% ------------------------------------------------------------------------------

% ----------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	try		delete([handles.hGuideLine handles.hZeroAge]),	end
	fid = fopen([handles.path_tmp 'mag_synt.ini'],'w');
	if (fid > 0)
		fprintf(fid,'# Config file holding the parameters of last usage of mag_synthetic.m\n');
		fprintf(fid,'%g %g %g %g %g %g\n', handles.speed, handles.spreadDir, handles.spreadObliqDir, ...
			handles.decl, handles.inc, handles.contamin);
		fclose(fid);
	end
	delete(handles.figure1);

% ----------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	handles = guidata(hObject);
	if isequal(get(hObject,'CurrentKey'),'escape')
		try		delete([handles.hGuideLine handles.hZeroAge]),	end
		delete(handles.figure1);
	end


% --- Creates and returns a handle to the GUI figure. 
function mag_synthetic_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'MenuBar','none',...
'Name','Mag Anom profiles',...
'NumberTitle','off',...
'Position',[520 298 319 401],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[10 10 301 381],'Style','frame');
uicontrol('Parent',h1,'Position',[20 270 281 103],'Style','frame');
uicontrol('Parent',h1,'Position',[20 179 281 83],'Style','frame');
uicontrol('Parent',h1,'Position',[20 69 281 101],'Style','frame');

uicontrol('Parent',h1,'Position',[35 327 61 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','2.5',...
'Style','edit',...
'Tooltip','Full spreading in cm/year',...
'Tag','edit_speed');

uicontrol('Parent',h1,'Position',[120 327 71 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','',...
'Style','edit',...
'Tooltip','Azimuth of spreading.',...
'Tag','edit_spreadDir');

uicontrol('Parent',h1,'Position',[210 327 71 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','0',...
'Style','edit',...
'Tooltip','Spreading obliquity in degrees (not always zero).',...
'Tag','edit_spreadObliqDir');

uicontrol('Parent',h1,'Position',[30 277 61 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','',...
'Style','edit',...
'Tooltip','Geomagnetic field declination',...
'Tag','edit_decl');

uicontrol('Parent',h1,'Position',[120 277 61 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','',...
'Style','edit',...
'Tooltip','Geomagnetic field inclination',...
'Tag','edit_inc');

uicontrol('Parent',h1,'Position',[210 277 71 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','1',...
'Style','edit',...
'Tag','edit_contamin');

uicontrol('Parent',h1,'Position',[75 234 170 21],...
'Callback',@mag_synthetic_uiCB,...
'String','Insert marker at Zero Age',...
'Tooltip','Insert a marker at the ridge intersection with the chunk line.',...
'Tag','push_zeroAge');

uicontrol('Parent',h1,'Position',[203 189 71 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','',...
'Style','edit',...
'Tooltip','Age of crust at the start of profile',...
'Tag','edit_ageInit');

uicontrol('Parent',h1,'Position',[30 128 85 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','3',...
'Style','edit',...
'Tooltip','Depth in km of sea-floor at ridge axis.',...
'Tag','edit_depthZero');

uicontrol('Parent',h1,'Position',[30 97 130 21],...
'Callback',@mag_synthetic_uiCB,...
'String','Get depth from grid',...
'Style','checkbox',...
'Tooltip','Get bathymetry profile by interpolation of a bat grid.',...
'Tag','check_batGrid');

uicontrol('Parent',h1,'Position',[30 77 241 21],...
'Callback',@mag_synthetic_uiCB,...
'HorizontalAlignment','left',...
'BackgroundColor',[1 1 1],...
'String','',...
'Style','edit',...
'Tooltip','Full name of bathymetry grid where to extract the depth profile',...
'Tag','edit_batGrid');

uicontrol('Parent',h1,'Position',[270 76 23 23],...
'Callback',@mag_synthetic_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'String','...',...
'Tag','push_batGrid');

uicontrol('Parent',h1,'Pos',[23  352 84 15], 'String','Spreading rate', 'Style','text');
uicontrol('Parent',h1,'Pos',[100 352 110 15],'String','Spreading dir', 'Style','text');
uicontrol('Parent',h1,'Pos',[193 352 110 15],'String','Spreading obliquity', 'Style','text');
uicontrol('Parent',h1,'Pos',[28  298 65 15], 'String','Declination', 'Style','text');
uicontrol('Parent',h1,'Pos',[120 298 60 15], 'String','Inclination', 'Style','text');
uicontrol('Parent',h1,'Pos',[204 298 80 15], 'String','Contamination', 'Style','text',...
'Tooltip','Contamination factor (< 1) in the sense of Tissot and Patriat');
uicontrol('Parent',h1,'Pos',[66 193 135 15], 'String','Age at begining of profile', 'Style','text');
uicontrol('Parent',h1,'Pos',[30 149 85 16],  'String','Depth at age 0', 'Style','text');

uicontrol('Parent',h1,'Pos',[152 115 25 18],...
'FontSize',10,...
'String','OR',...
'Style','text');

 uicontrol('Parent',h1,'Pos',[152 211 25 18],...
'FontSize',10,...
'String','OR',...
'Style','text');

uicontrol('Parent',h1,'Position',[135 134 151 16],...
'FontSize',10,...
'ForegroundColor',[0.7 0 0],...
'String','PROFILE HAS DEPTH',...
'Style','text',...
'Vis', 'off',...
'Tag','text_gotDepth');

uicontrol('Parent',h1,'Position',[20 379 270 20],...
'FontSize',9,...
'String','Generate a companion synthetic profile (opt)',...
'Style','text',...
'Tag','text_synthet');

uicontrol('Parent',h1,'Position',[20 15 137 23],...
'String','Compute RTP too',...
'Style','checkbox',...
'Tooltip','Compute also an Reduced To the Pole anomaly too',...
'Tag','check_RTP');

uicontrol('Parent',h1,'Position',[220 18 80 21],...
'Callback',@mag_synthetic_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'String','Run',...
'Tag','push_run');

function mag_synthetic_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

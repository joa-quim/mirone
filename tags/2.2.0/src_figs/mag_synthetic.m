function varargout = mag_synthetic(varargin)
% Extract a chunk of a mgd77++ file dictated from 'guideline' drawn interactively
%
% This function has the magmodel repeated from 'ecran'. This is bad practice

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

	% Set an uicontext for profile extraction of the mgd77++ file
	cmenuHand = uicontextmenu('Parent',handles.hCallingFig);
	set(handles.hGuideLine, 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Extract mgd77++ profile', 'Call', {@getProfile,hObject});
	uimenu(cmenuHand, 'Label', 'Delete me', 'Call', {@del_line,hObject}, 'Sep','on');
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
function edit_spreadingDir_CB(hObject, handles)
	xx = sscanf(get(hObject,'Str'), '%f');
	if (isnan(xx)),		set(hObject,'Str',sprintf('%.1f',handles.spreadDir)),	return,		end
	handles.spreadDir = xx;
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function edit_ageInit_CB(hObject, handles)
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
		errordlg('insulto-declina','Error'),	return
	elseif (isempty(handles.inc))
		errordlg('insulto-inclina','Error'),	return
	elseif (isempty(handles.spreadDir))
		errordlg('insulto-spreadDir','Error'),	return
	elseif (isempty(handles.hZeroAge))
		errordlg('insulto-zeroAge','Error'),	return
	elseif (~ishandle(handles.hZeroAge))
		errordlg('Who told you to kill the Zero Age marker? Bye-Bye.','Error'),	return
	end

	pt_x = get(handles.hZeroAge,'XData');		pt_y = get(handles.hZeroAge,'YData');	% Get the Zero age location
	ind_0 = find_nearpoint(handles.trk_chunk.lon, handles.trk_chunk.lat, pt_x, pt_y);

	if ( ~isempty(handles.trk_chunk) && ~isempty(handles.trk_chunk.topography) )
		distance = handles.trk_chunk.distance - handles.trk_chunk.distance(ind_0);
	elseif (get(handles.check_batGrid,'Val') && isempty(handles.batProfile))
		if (isempty(handles.hGuideLine))
			errordlg('You need to draw a profile line first. The "Select line chunk" button','Error')
			return
		end
		fname = get(handles.edit_batGrid,'str');
		if (isempty(fname))
			errordlg('insulto-bat','Error'),	return
		elseif (exist(fname,'file') ~= 2)
			errordlg('Bathymetry grid does not exist.','Error'),	return
		end
		[handles, X, Y, Z, head] = read_gmt_type_grids(handles,fname);
		if (isempty(X)),	return,		end
		x = get(handles.hGuideLine,'XData');		y = get(handles.hGuideLine,'YData');
		handles.batProfile = grdtrack_m(Z,head,[x(:) y(:)],'-Z');
	end

	magRTP = [];
	if ( get(handles.check_RTP,'Val') )			% See if we need to compute the RTP profile as well
		[theta,ampfac] = nskew(0, handles.spreadDir-90, handles.decl, handles.inc, pt_x);
		magRTP = -rtp2d(handles.trk_chunk.mag, theta, ampfac);
	end

	dxyza = [distance(:) handles.trk_chunk.lon(:) handles.trk_chunk.lat(:) handles.trk_chunk.topography(:)/1000 handles.trk_chunk.mag(:)];
	anoma = magmodel(dxyza, handles.decl, handles.inc, handles.speed, handles.spreadDir, handles.chunkAzim, handles.contamin);
	if ( ~isempty(handles.hEcran) && ishandle(handles.hEcran))			% If we are going to append to a ecran Fig
		if (isempty(handles.hSynthet) || ~ishandle(handles.hSynthet) )	% If model line does not exist yet
			handEcran = guidata(handles.hEcran);
			handles.hSynthet = line('XData',distance-distance(1), 'YData',anoma, 'Parent',handEcran.axes1);
			a1 = abs(distance(1)) / (handles.speed / 2 * 10);			% Remember, half speed
			a2 = abs(distance(end)) / (handles.speed / 2 * 10);
			handEcran.ageEnd = a2;
			if (a1 > a2),		handEcran.ageEnd = a1;		end			% Pick the abs longer of both half lengths
			handEcran.ageStart = 0;

			if (~isempty(magRTP))
				handles.hRTP = line('XData',distance-distance(1), 'YData',magRTP, 'Color','b', 'LineStyle','--','Parent',handEcran.axes1);
			end

			guidata(handEcran.figure1, handEcran)
			cb = get(handEcran.push_magBar,'Call');						% Clever trick to fish the callback
			feval(cb{1}, handEcran.push_magBar, [], cb{2}, cb{3})		% Call the ecran_uiCB function
			set(handEcran.axes2,'xlim',[-a1 a2])
			y_lim = get(handEcran.axes1,'ylim');
			set(handEcran.axes1,'ylim',[min(y_lim(1),min(anoma)) max(y_lim(2),max(anoma))])
			hP = findobj(handEcran.axes2,'Type','patch');				% Fish prev patch and plot it mirrored.
			vert = get(hP,'Vertices');	faces = get(hP,'Faces');	cor = get(hP,'FaceVertexCData');
			patch('Parent',handEcran.axes2,'Faces',faces,'Vertices',[-vert(:,1) vert(:,2)],'FaceVertexCData',cor,'FaceColor','flat');

			uistack_j(handEcran.axes2,'up')				% Bring it up relatively to axes1 so that we can see the ticks
			ages = get(handEcran.axes2, 'UserData');
			hT = findobj(handEcran.axes1,'Tag','chroneName');
			set(hT,'Vis','off')
			if (ages(1) == 0),		ages(1) = [];	hT(end) = [];	end	% We don't want the Zero tick
			ages = [-ages(end:-1:1) ages];								% Will be new ticks
			set(handEcran.axes2,'XTick', ages(:))
			str = get(hT,'Str');
			set(handEcran.axes2,'XTickLabel',[str; str(end:-1:1)]);		% And finally, the new tickmarks

		else				% Line exists. Update it
			set(handles.hSynthet, 'XData',distance-distance(1), 'YData',anoma)
		end
	else
		figure;		plot(distance-distance(1), anoma)
	end
	guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------
function getProfile(obj, evt, hFig)
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
		guidata(handles.figure1, handles);
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
	out = igrf_m(lon, lat, 0, 1990, '-Fdi');
	handles.decl = out(1);		handles.inc = out(2);
	set(handles.edit_decl, 'Str', out(1)),		set(handles.edit_inc, 'Str', out(2))
	guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------
function check_RTP_CB(hObject, handles)
% ...

% --------------------------------------------------------------------
function trk_chunk = get_track_chunk(handles)
% ...

	[PATH,FNAME,EXT] = fileparts(handles.f_name);
	if (isempty(PATH)),		PATH = cd;		end		% Play safe
	if (strcmpi(EXT,'.nc'))
		handles.f_name = [PATH filesep FNAME '.nc'];		handles.is_mgd77 = true;		handles.is_gmt = false;
		[handles, track] = read_mgd77_plus(handles, handles.f_name);
	else
		handles.f_name = [FNAME '.gmt'];	handles.is_gmt = true;			handles.is_mgd77 = false;
		track = gmtlist_m([PATH filesep FNAME], '-Fsxygmtd', handles.opt_G);
	end

	if (length(track.agency) ~= 10)			% Ensures that agency is exactly 10 chars
		agency = '          ';				% 10 blanks
		len = min(length(track.agency),10);
		agency(1:len) = track.agency(1:len);
		trk_chunk.agency = agency;
	else
		trk_chunk.agency = track.agency;
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
	trk_chunk.mag = track.magnetics(ind);
	if (step == 1)				% Picked track runs in the same sense as data in file
		trk_chunk.distance = track.distance(ind);
	else						% Track and data in file run in oposit sense. Adjust acumulated dist accordingly.
		t = track.distance(ind);
		trk_chunk.distance = [0; cumsum( diff(t(end:-1:1)) )];
	end

	if ( numel(track.topography) > 1 && ~all(isnan(track.topography)) )
		trk_chunk.topography = track.topography(ind);
	else
		trk_chunk.topography = [];
	end

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
function [anoma, age_line] = magmodel(dxypa, chtdec, chtinc, speed, dir_spread, dir_profil, contam)
% This function is in large part taken from the MAGMOD program
% Intention is to clean it more or eventually re-write it (for example it doesn't account
% for possibly different Remanent & Inducted maags).

	if (nargin == 6),	contam = 1;		end
	if (size(dxypa,2) == 4)		% If depth not provided, default to const = 2.5 km
		dxypa = [dxypa(:,1:3) ones(size(dxypa,1),1)*2.5 dxypa(:,4)];
	end

	levelobs = 0;	aimaaxe = 18;		aimanta = 4;
	psubsi = 0.35;		% Coefficient of thermal subsidence equation
	subsi = 1;			% PORQUÊ????
	profond = 0;		% PORQUÊ????
	thickness = 0.5;
	stat_x = dxypa(:,1);
	txouv = speed * 10;		% From full rate in cm/year to rate in km/Ma 

	mir_dirs = getappdata(0,'MIRONE_DIRS');
	if (~isempty(mir_dirs))
		chron_file = [mir_dirs.home_dir filesep 'data' filesep 'Cande_Kent_95.dat'];
	else
		chron_file = [cd filesep 'data' filesep 'Cande_Kent_95.dat'];
	end

	fid = fopen(chron_file,'r');
	todos = fread(fid,'*char');
	[chron age_start age_end age_txt] = strread(todos,'%s %f %f %s');
	fclose(fid);    clear todos
	BA = zeros(numel(age_start)*2,2);
	BA(1:2:end-1,1) = age_start;	BA(2:2:end,1) = age_end;
	BA(1:2:end-1,2) = age_end;
	BA(2:2:end-2,2) = age_start(2:end);
	BA(end,:) = [];
	
	debut_bloc = 1;
	fin_bloc = size(BA,1);

	if (nargin == 0)
		Ficp = fopen ('C:\j\artigos\review\cg_modmag\publicado\swir.dxypa');
		MX = (fscanf(Ficp,'%f',[5,inf]))';
		fclose(Ficp);
		MXS = sortrows(MX,1);
	else
		MXS = dxypa;
	end
	distfic = MXS(:,1);
	proffic = MXS(:,4);
	ind = isnan(proffic);
	if ( any(ind) )			% Many files have gaps in bathymetry. Reinvent it
		proffic(ind) = interp1(distfic(~ind),proffic(~ind),distfic(ind));
	end
	Nficprof = numel(distfic);
	
	% Number of points
	nb_stat_mod = numel(stat_x);

	% Depth of the points where the magnetic anomaly will be compute
	stat_z = zeros(nb_stat_mod,1)+levelobs;

	% Calculation of the position of the magnetized bodies

	% Determination of the age limits of normal and inverse polarity bodies and
	% of their respecting magnetization

	nb_struct = (fin_bloc*2)-1;

	blocag = zeros(nb_struct,2);
	aimbloc = zeros(nb_struct,1);
	blocag(debut_bloc:fin_bloc-1,1) = -BA(fin_bloc:-1:debut_bloc+1,2);
	blocag(debut_bloc:fin_bloc-1,2) = -BA(fin_bloc:-1:debut_bloc+1,1);
	blocag(fin_bloc,:) = [-BA(1,2) BA(1,2)];
	blocag(fin_bloc+1:nb_struct,:) = BA(debut_bloc+1:fin_bloc,:);
	aimbloc(fin_bloc,1) = aimaaxe;
	aimbloc(fin_bloc-2:-2:1,1) = aimanta;
	aimbloc(fin_bloc+2:2:nb_struct,1) = aimanta;
	aimbloc(fin_bloc-1:-2:1,1) = aimanta*(-1.);
	aimbloc(fin_bloc+1:2:nb_struct,1) = aimanta*(-1.);

	% Initialisation of magnetized bodies position in km
	polygon_x = zeros(4,nb_struct);

	% Calculation of magnetized bodies position in km
	polygon_x(1,:) = blocag(:,1)' * txouv(1) / 2;
	polygon_x(2,:) = blocag(:,2)' * txouv(1) / 2;
	polygon_x(3,:) = polygon_x(2,:);
	polygon_x(4,:) = polygon_x(1,:);

	% Before calculation of the magnetic anomaly, the spreading direction has
	% to be less than 90° away from the profile direction in order to be plot
	% in the right sense of distance. To do this we first calculate the obliquity
	% between the direction of the profile and the spreading direction.

	obliquity = dir_profil - dir_spread;
	if (abs(obliquity) > 90)
		temp1 = mod(obliquity,90);
		temp2 = mod(obliquity,-90);
		if (abs(temp1) <= abs(temp2)),		obliquity = temp1;
		else								obliquity = temp2;
		end
		dir_spread = dir_profil - obliquity;
	end

	% Calculation of the magnetized bodies depth
	if (subsi == 1),	cosubsi = psubsi;
	else				cosubsi = 0.0;
	end
	polygon_z = zeros(4,nb_struct);
	if profond ~= 0.0
		polygon_z(1,:) = profond + cosubsi*sqrt(abs(blocag(:,1)'));
		polygon_z(2,:) = profond + cosubsi*sqrt(abs(blocag(:,2)'));
		polygon_z(3,:) = polygon_z(2,:) + thickness;
		polygon_z(4,:) = polygon_z(1,:) + thickness;
		if cosubsi ~= 0.0
			indbmil = polygon_x(1,:) < 0 & polygon_x(2,:) > 0;
			shift = profond - polygon_z(1,indbmil);
			polygon_z(:,:) = polygon_z(:,:) + shift;
		end
	else
		PolX = cell(1,nb_struct);
		PolZ = cell(1,nb_struct);
		if  (obliquity ~= 0)
			polygon_x = polygon_x / cos (obliquity * pi / 180);
		end
		lignm11 = find(polygon_x(1,:) < distfic(1));
		if ~isempty(lignm11)
			maxlignm11 = max(lignm11);
			if polygon_x(1,maxlignm11) < 0
				polygon_z(1,lignm11(:)) = proffic(1) + cosubsi*sqrt(blocag(maxlignm11,1)-blocag(lignm11(:),1)');
			else
				polygon_z(1,lignm11(:)) = proffic(1);
			end
			polygon_z(4,lignm11(:)) = polygon_z(1,lignm11(:)) + thickness;
		end
		lignm12 = find(polygon_x(2,:) < distfic(1));
		if ~isempty(lignm12)
			maxlignm12 = max(lignm12);
			if polygon_x(2,maxlignm12) < 0
				polygon_z(2,lignm12(:)) = proffic(1) + cosubsi*sqrt(blocag(maxlignm12,2)-blocag(lignm12(:),2)');
			else
				polygon_z(2,lignm12(:)) = proffic(1);
			end
			polygon_z(3,lignm12(:)) = polygon_z(2,lignm12(:)) + thickness;
		else
			maxlignm12=0;
		end

		lignm21 = find(polygon_x(1,:) >= distfic(Nficprof));
		if ~isempty(lignm21)
			minlignm21 = min(lignm21);
			if (polygon_x(1,minlignm21) > 0)
				polygon_z(1,lignm21(:)) = proffic(Nficprof) + cosubsi*sqrt(blocag(lignm21(:),1)' - blocag(minlignm21,1));	
			else
				polygon_z(1,lignm21(:)) = proffic(Nficprof);
			end
			polygon_z(4,lignm21(:)) = polygon_z(1,lignm21(:)) + thickness;
		else
			minlignm21=nb_struct+1;
		end
		lignm22 = find(polygon_x(2,:) >= distfic(Nficprof));
		if ~isempty(lignm22)
			minlignm22 = min(lignm22);
			if polygon_x(2,minlignm22) > 0
				polygon_z(2,lignm22(:)) = proffic(Nficprof) + cosubsi*sqrt(blocag(lignm22(:),2)' - blocag(minlignm22,2));	
			else
				polygon_z(2,lignm22(:)) = proffic(Nficprof);
			end
			polygon_z(3,lignm22(:)) = polygon_z(2,lignm22(:)) + thickness;
		end

		lignm31  = find(polygon_x(1,:) >= distfic(1) & polygon_x(1,:) < distfic(Nficprof));
		if ~isempty(lignm31)
			polygon_z(1,lignm31(:)) = interp1(distfic,proffic,polygon_x(1,lignm31(:)));
			polygon_z(4,lignm31(:)) = polygon_z(1,lignm31(:)) + thickness;
		end
		lignm32  = find(polygon_x(2,:) >= distfic(1) & polygon_x(2,:) < distfic(Nficprof));
		if ~isempty(lignm32)
			polygon_z(2,lignm32(:)) = interp1(distfic,proffic,polygon_x(2,lignm32(:)));
			polygon_z(3,lignm32(:)) = polygon_z(2,lignm32(:)) + thickness;
		end
	end

	PolXX = cell(1,nb_struct);
    PolZZ = cell(1,nb_struct);

	for i=1:nb_struct
		PolX(1,i)  = {polygon_x(:,i)};
		PolXX(1,i) = {[polygon_x(:,i);polygon_x(1,i)]};
		PolZ(1,i)  = {polygon_z(:,i)};
		PolZZ(1,i) = {[polygon_z(:,i);polygon_z(1,i)]};
	end
	if profond==0
		for i=maxlignm12+1:minlignm21-1
			ptdist=find(distfic > polygon_x(1,i)+0.00001 & distfic < polygon_x(2,i)-0.00001);
			if ~isempty(ptdist)
				tempox1 = distfic(ptdist(:),1);
				tempoz1 = proffic(ptdist(:),1);
				tempox2 = flipud(tempox1);
				tempoz2 = flipud(tempoz1+thickness);
				tempox3 = [polygon_x(1,i);tempox1;polygon_x(2,i);polygon_x(3,i);tempox2;polygon_x(4,i)];
				tempoz3 = [polygon_z(1,i);tempoz1;polygon_z(2,i);polygon_z(3,i);tempoz2;polygon_z(4,i)];
				PolX(1,i)={tempox3};
				PolXX(1,i)={[tempox3;polygon_x(1,i)]};
				PolZ(1,i)={tempoz3};
				PolZZ(1,i)={[tempoz3;polygon_z(1,i)]};
			end
		end
	end

	% Re-positionning of the magnetized bodies if the contamination coefficient is different from 1

	if contam ~= 1
		for i=1:nb_struct
			PolXX{1,i}=PolXX{1,i} * contam;
		end
		stat_x(:,1) = stat_x(:,1) * contam;
	end

	% Calculation of the magnetic anomaly created by the magnetized bodies

	anoma = calcmag(nb_struct,chtinc,chtdec,dir_spread,aimbloc,stat_x,stat_z,nb_stat_mod,PolXX,PolZZ);

	if (nargout == 2)						% Calculate the ages. Note that this will be RELATIVE
		age_line = distfic / txouv * 2;		% to the age of ofirst point in ptofile (r = 0)
	end

% ---------------------------------------------------------------------
function anoma = calcmag(nb_struct,chtinc,chtdec,dir_spread,aimbloc,stat_x,stat_z,nb_stat_mod,PolXX,PolZZ)
% Calculation of the magnetic anomaly created by magnetized bodies

	D2R = pi / 180;
	aiminc(1:nb_struct,1) = chtinc;
	aimdir(1:nb_struct,1) = chtdec;
    dir_anomag = dir_spread+90;

	c1 = sin(chtinc*D2R);    
	c2 = cos(chtinc*D2R)*cos(dir_spread*D2R - chtdec*D2R);

	d1 = sin(aiminc*D2R);
	d2 = cos(aiminc*D2R).*cos((dir_anomag-90)*D2R - aimdir*D2R);
	d3 = 200.*aimbloc;

	anomax = 0;		anomaz = 0;

	for ik = 1:nb_struct
		nbpoints = length(PolXX{1,ik});
		[amx,amz] = fcalcmagpt(nbpoints,stat_x,stat_z,nb_stat_mod,PolXX{1,ik},PolZZ{1,ik},d1(ik),d2(ik),d3(ik));
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

uicontrol('Parent',h1,'Position',[30 327 72 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','2.5',...
'Style','edit',...
'TooltipString','Full spreading in cm/year',...
'Tag','edit_speed');

uicontrol('Parent',h1,'Position',[190 327 81 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','',...
'Style','edit',...
'TooltipString','Azimuth of spreading.',...
'Tag','edit_spreadingDir');

uicontrol('Parent',h1,'Position',[30 277 61 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','',...
'Style','edit',...
'TooltipString','Geomagnetic field declination',...
'Tag','edit_decl');

uicontrol('Parent',h1,'Position',[113 277 61 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','',...
'Style','edit',...
'TooltipString','Geomagnetic field inclination',...
'Tag','edit_inc');

uicontrol('Parent',h1,'Position',[199 277 71 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','1',...
'Style','edit',...
'Tag','edit_contamin');

uicontrol('Parent',h1,'Position',[75 234 170 21],...
'Callback',{@mag_synthetic_uiCB,h1,'push_zeroAge_CB'},...
'String','Insert marker at Zero Age',...
'Tooltip','Insert a marker at the ridge intersection with the chunk line.',...
'Tag','push_zeroAge');

uicontrol('Parent',h1,'Position',[203 189 71 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','',...
'Style','edit',...
'TooltipString','Age of crust at the start of profile',...
'Tag','edit_ageInit');

uicontrol('Parent',h1,'Position',[30 128 85 21],...
'Callback',@mag_synthetic_uiCB,...
'BackgroundColor',[1 1 1],...
'String','3',...
'Style','edit',...
'TooltipString','Depth in km of sea-floor at ridge axis.',...
'Tag','edit_depthZero');

uicontrol('Parent',h1,'Position',[30 97 130 21],...
'Callback',@mag_synthetic_uiCB,...
'String','Get depth from grid',...
'Style','checkbox',...
'TooltipString','Get bathymetry profile by interpolation of a bat grid.',...
'Tag','check_batGrid');

uicontrol('Parent',h1,'Position',[30 77 241 21],...
'Callback',@mag_synthetic_uiCB,...
'HorizontalAlignment','left',...
'BackgroundColor',[1 1 1],...
'String','',...
'Style','edit',...
'TooltipString','Full name of bathymetry grid where to extract the depth profile',...
'Tag','edit_batGrid');

uicontrol('Parent',h1,'Position',[270 76 23 23],...
'Callback',@mag_synthetic_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'String','...',...
'Tag','push_batGrid');

uicontrol('Parent',h1,'Position',[174 352 110 15],...
'String','Spreading direction',...
'Style','text');

uicontrol('Parent',h1,'Position',[23 352 84 15],...
'String','Spreading rate',...
'Style','text');

uicontrol('Parent',h1,'Position',[66 193 135 15],...
'String','Age at begining of profile',...
'Style','text');

uicontrol('Parent',h1,'Position',[30 149 85 16],...
'String','Depth at age 0',...
'Style','text');

uicontrol('Parent',h1,'Position',[152 115 25 18],...
'FontSize',10,...
'String','OR',...
'Style','text');

 uicontrol('Parent',h1,'Position',[152 211 25 18],...
'FontSize',10,...
'String','OR',...
'Style','text');

uicontrol('Parent',h1,'Position',[28 298 65 15],...
'String','Declination',...
'Style','text');

uicontrol('Parent',h1,'Position',[113 298 60 15],...
'String','Inclination',...
'Style','text');

uicontrol('Parent',h1,'Position',[194 298 80 15],...
'String','Contamination',...
'Style','text',...
'TooltipString','Contamination factor (< 1) in the sense of Tissot and Patriat');

uicontrol('Parent',h1,'Position',[138 134 151 16],...
'FontSize',10,...
'ForegroundColor',[0 1 0],...
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
'Callback',@mag_synthetic_uiCB,...
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

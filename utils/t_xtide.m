function pred = t_xtide(varargin)
% T_XTIDE Tidal prediction
% YOUT=T_XTIDE(LONG,LAT) makes a tidal prediction for the current day using the
% harmonics file from XTIDE. LONG,LAT is used to find the closest station to that position
%
% The times of predicted tides are given by the next numerical argument
% (if any), e.g. [...]=T_XTIDE(LONG,LAT,TIM). 
% TIM can be: a vector of matlab-format decimal days (from DATENUM).
%           : a scalar <1000, taken as the number of days from present
%           : a scalar >1000, taken as the starting time in matlab-format 
%             for a 2 day time series. 
%           : not given, in which case the current time is used as a start time 
%
% Times are usually taken to be in standard time for the given location (no daylight savings adjustments); 
% if in doubt use the 'info' or 'full' options where offset from UTC is given.
%
% 
% Other optional arguments can be specified following this using property/value pairs: 
%
%     'format'     'raw' (default)
%                    YOUT is just a time series to match the time in TIM
%
%                  'times'
%                    YOUT is a structure of high/low tide information
%                    between times min(TIM) and max(TIM).
%
%                  'info'
%                    YOUT is a structure giving station information
%                    (location, time zone, units, harmonic constituents)
%
%                  'full'
%                    Combination of 'raw' and 'info' in a structure YOUT.
%
%     'units'     {'meters' | 'feet' | 'm/s' | 'knots' | 'original' }
%                    Units of result (default is original units)
%
% If no output argument is specified data is plotted and/or displayed.
%
%  Requires the xtide harmonics file  - get this from http://bel-marduk.unh.edu/xtide/files.html

% R. Pawlowicz 1/Dec/2001
% Version 1.0
%          16/May/02 - added lat/long options (thanks to Richard Dewey).
%
%   As usual, Mironified by J. Luis

	% Star by finding if we are being reused or it is the first time this function is called
	hFig = findobj('Type','figure','Tag','TidalPred');
	if (isempty(hFig))     % First time use. Need to load harmonics file
		load(['data' filesep 't_xtide.mat']);
	else                    % We already have the harmonics somwhere. Get them
		xharm = getappdata(hFig,'XHARM');    xtide = getappdata(hFig,'XTIDE');
		if (isempty(xharm)),	load(['data' filesep 't_xtide.mat']);   end     % Something wrong occured before
	end

	dist = t_gcdist(xharm.latitude,xharm.longitude,varargin{2},varargin{1});
	[mind,ista] = min(dist);
	curr_pt = [varargin{1:2}];
	varargin(1:2)=[];

	% Time vector (if available) otherwise take current time.
	if (~isempty(varargin) && ~isa(varargin{1},'char'))
		tim = varargin{1};
		varargin(1) = [];
		if (length(tim) == 1)
			if (tim < 1000)
				dat = clock;
				tim = datenum(dat(1),dat(2),dat(3))-0.5+(0:1/48:tim);
			else
				tim = tim+(0:1/48:2)-0.5;       % 2 days worth.
			end	 	
		end
	else 
		dat = clock;    tim = datenum(dat(1),dat(2),dat(3))-0.5+(0:.25:48)/24;  clear dat;
	end

	format = 'raw';     unt = 'original';

	while (~isempty(varargin))
		switch lower(varargin{1}(1:3)),
			case 'for',     format = lower(varargin{2});
			case 'uni',     unt = lower(varargin{2}); 
			otherwise,      error(['Can''t understand property:' varargin{1}]);
		end
		varargin([1 2])=[]; 
	end

	% If we want a time series
	pred = [];
	[units,convf] = convert_units(unt,xharm.units(ista,:));       % Convert units if requested.
	if strcmp(format(1:2),'ra') || strcmp(format(1:2),'fu') || strcmp(format(1:2),'ti')
		if strcmp(format(1:2),'ti')             % Data every minute for hi/lo forecasting.
			tim = tim(1):(1/1440):tim(end); 
		end

		% Convert into time since the beginning of year
		mid = datevec(mean(tim));		iyr = mid(1)-xtide.startyear+1;
		lt  = length(tim);
		xtim = (tim-datenum(mid(1),1,1))*24;    % Hours since beginning of year

		%----- Sum up everything for the prediction! ------------
		pred = xharm.datum(ista)+sum(repmat(xtide.nodefactor(:,iyr).*xharm.A(ista,:)',1,lt).* ...
			cos( ( xtide.speed*xtim + repmat(xtide.equilibarg(:,iyr)-xharm.kappa(ista,:)',1,lt) )*(pi/180) ),1);
		pred = pred * convf;
		%-----------------------------------------------------

		% Compute times of hi/lo from every-minute data
		if strcmp(format(1:2),'ti')
			% Check if this is a current station
			%if ~isempty(strfind('Current',xharm.station(ista,:))), currents=1; else currents=0; end;
			ddpred = diff(diff(pred) > 0);
			flat = find(ddpred ~= 0) + 1;

			hi.mtime = tim(flat);           hi.value = pred(flat);
			hi.type = zeros(size(flat));
			hi.type(ddpred(flat-1) < 0) = 1;  % 0=lo, 1=hi
			hi.units = deblank(units);      pred = hi;
		end
	end

	% Create information structure
	if strcmp(format(1:2),'in') || strcmp(format(1:2),'fu'),
		if (~isempty(pred))
			pred.yout=pred;     pred.mtime=tim; 
		else
			kk = find(xharm.A(ista,:) ~= 0);
			pred.freq = xtide.name(kk,:);
			pred.A = full(xharm.A(ista,kk)')*convf;
			pred.kappa = full(xharm.kappa(ista,kk)'); 
		end
		pred.longitude = xharm.longitude(ista);     pred.latitude = xharm.latitude(ista);
		pred.timezone = xharm.timezone(ista);       pred.units = deblank(units);
		pred.datum = xharm.datum(ista)*convf;       pred.station = deblank(xharm.station(ista,:));
	end

% 	% ------------
% 	if (strcmp(format(1:2),'ti') & nargout > 0)
% 		lh_times = datestr(hi.mtime);
% 		lh_heights = num2str(hi.value,'%.1f\n')';
% 		data_ini = lh_times(1,1:11);
% 		c = calendar(data_ini);
% 		time_tides = lh_times(:,13:17);
% 		str = sprintf('Sunday\tMonday\tTuesday\tWednesday\tThursday\tFriday\tSaturday\n\n');
% 		ind = find(c > 0);
% 		tmp = [time_tides ' ' lh_heights];
% 		tmp = repmat(tmp,6,7);
% 	end
% 	% ------------

	% If no output parameters then we plot things
	if (nargout == 0)
		hFig = findobj('Type','figure','Tag','TidalPred');
		if (isempty(hFig))
			[hFig,hAxes] = t_fig;     % Create a new figure
			setappdata(hFig,'XHARM',xharm);    setappdata(hFig,'XTIDE',xtide);    % First time, save those
			reuse_fig = 0;
		else
			hAxes = get(hFig,'CurrentAxes');  reuse_fig = 1;
		end
		low_high = t_xtide(curr_pt(1),curr_pt(2),'format','ti');    % Get the next hi/low tides
		lh_times = datestr(low_high.mtime(1:2),16);
		lh_heights = low_high.value(1:2);
		%line(tim,pred)
		setappdata(hFig,'CurrentPoint',curr_pt);   % Save current point (the original click on the map window)

		start_time = datevec(tim(1));
		rel_hr = (tim - datenum(start_time(1),start_time(2),start_time(3)))*24; % Hours since begining of display date
		setappdata(hFig,'TideData',[rel_hr; pred]);    % Save the tide data for eventual exporting

		ddpred = diff(diff(pred) > 0);
		flat = find(ddpred ~= 0) + 1;
		X = cell(length(flat)+1,1);    Y = X;
		X{1} = tim(1:flat(1));    Y{1} = pred(1:flat(1));
		for (k=2:length(flat))
			X{k} = tim(flat(k-1):flat(k));
			Y{k} = pred(flat(k-1):flat(k));
		end
		X{end} = tim(flat(end):end);    Y{end} = pred(flat(end):end);
		y_pos = get(hAxes,'Ylim');

		% When a different date was selected (with the calendar) it's easier to start over again
		if (reuse_fig),		delete(get(hAxes,'Children'));     end

		hold on
		for (k=1:length(X))
			if (Y{k}(end) > Y{k}(1)),	cor = [0 0.6275 1];
			else						cor = [0 1 0.4];
			end
			patch('Parent',hAxes,'XData',[X{k} X{k}(end) X{k}(1)],'YData',[Y{k} y_pos(1) y_pos(1)], ...
			      'FaceColor',cor,'EdgeColor','none')
		end
		hold off

		dat = clock;    exact_now = datenum(dat(1),dat(2),dat(3),dat(4),dat(5),dat(6));     clear dat;
		if (exact_now >= tim(1) && exact_now <= tim(end))    % If current time fits inside the TIM vector
			y_now = interp1(tim,pred,exact_now);
			line(exact_now,y_now,'Marker','+','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',2)
			txt_now = ['Time now ' datestr(exact_now,16) ' UTC ' num2str(y_now,'%.1f') 'm'];
			if (lh_heights(1) > lh_heights(2))      % Next is a high tide
				txt_h = ['Next High Tide ' lh_times(1,:) ' UTC ' num2str(lh_heights(1),'%.1f') 'm'];
				txt_l = ['Next Low Tide ' lh_times(2,:) ' UTC ' num2str(lh_heights(2),'%.1f') 'm'];
			else                                    % Next is a low tide
				txt_h = ['Next High Tide ' lh_times(2,:) ' UTC ' num2str(lh_heights(2),'%.1f') 'm'];
				txt_l = ['Next Low Tide ' lh_times(1,:) ' UTC ' num2str(lh_heights(1),'%.1f') 'm'];
			end
			x_pos = tim(fix(end/2));        y_pos = get(hAxes,'Ylim');
			ht1 = text(x_pos,y_pos(2),txt_h,'Parent',hAxes,'HorizontalAlignment','center','VerticalAlignment','top');
			ht2 = text(x_pos,y_pos(1)+.1,txt_l,'Parent',hAxes,'HorizontalAlignment','center','VerticalAlignment','bot');
			ht3 = text(x_pos,y_pos(2)-.15,txt_now,'Parent',hAxes,'HorizontalAlignment','center','VerticalAlignment','top');
			set(ht1, 'ButtonDownFcn', 'move_obj(1)')	% Set the dragging option
			set(ht2, 'ButtonDownFcn', 'move_obj(1)')
			set(ht3, 'ButtonDownFcn', 'move_obj(1)')
		end

		tit = ['Tidal prediction for ',deblank(xharm.station(ista,:)) ' beginning ' datestr(tim(1))];
		setappdata(hFig,'StartDate',tit);       % Save the starting time for eventual exporting
		datetick(hAxes);
		set(hFig,'Name',tit)
		ylabel(deblank(xharm.units(ista,:)));
		setappdata(hAxes, 'LabelFormatType', 'Date');		% Tell pixval_stsbar to display XX coords as Date-Time
		set(hFig,'Visible','on')
	end
  
% --------------------------------------------------------------------------
function [units,convf] = convert_units(unt,origunits)
% Conversion factors from origianl units if requested and possible
% (no conversions from knots to feet).

	if strcmp(unt(1:3),origunits(1:3)) || strcmp(unt(1:3),'ori'),
		units = origunits;
		convf = 1;
	else
		switch unt(1:3),
			case 'fee',
				if strcmp(origunits(1:3), 'met'),
					units = 'feet';     convf = 3.2808399;
				else
					units = origunits;  convf = 1;
				end
			case 'met',
				if strcmp(origunits(1:3), 'fee'),
					units = 'meters';   convf = 0.3048;
				else
					units = origunits;  convf = 1;
				end
			case 'm/s',
				if strcmp(origunits(1:3), 'kno'),
					units = 'meters/sec';   convf = 0.51444444;
				else
					units = origunits;      convf = 1;
				end
			case 'kno',
				if strcmp(origunits(1:3), 'm/s'),
					units = 'knots';        convf = 1.9438445;
				else
					units = origunits;      convf = 1;
				end
			otherwise
				error('Unknown units')
			end
	end

% ---------------------------------------------------------------------------
function [d,hdg] = t_gcdist(lat1,lon1,lat2,lon2)
% Function to calculate distance in kilometers and heading between two
% positions in latitude and longitude.
% Assumes -90 > lat > 90  and  -180 > long > 180
%    north and east are positive
% Uses law of cosines in spherical coordinates to calculate distance
% calculate conversion constants
%
%  Code from Richard Dewey.

	raddeg = 180/pi;
	degrad = 1/raddeg;
	% convert latitude and longitude to radians
	lat1 = lat1 * degrad;   lat2 = lat2 * degrad;
	in1 = find(lon1>180);   lon1(in1) = lon1(in1)-360;
	in2 = find(lon2>180);   lon2(in2) = lon2(in2)-360;
	lon1 = -lon1.*degrad;   lon2 = -lon2.*degrad;
	% calculate some basic functions
	coslat1=cos(lat1);      sinlat1=sin(lat1);
	coslat2=cos(lat2);      sinlat2=sin(lat2);
	%calculate distance on unit sphere
	dtmp=cos(lon1-lon2);
	dtmp=sinlat1.*sinlat2 + coslat1.*coslat2.*dtmp;

	% check for invalid values due to roundoff errors
	in1= dtmp > 1;     dtmp(in1)=1.0;
	in2= dtmp < -1;    dtmp(in2)=-1.0;

	% convert to meters for earth distance
	ad = acos(dtmp);
	d=(111.112) .* raddeg .* ad;

	% now find heading (If it was required)
	if (nargout == 2)
		hdgcos = (sinlat2-sinlat1.*cos(ad))./(sin(ad).*coslat1);

		% check value to be legal range
		in1 = hdgcos > 1.0;   hdgcos(in1) = 1.0;
		in2 = hdgcos < -1.0;  hdgcos(in2) = -1.0;
		hdg = acos(hdgcos).*raddeg;

		% if longitude is decreasing then heading is between 180 and 360
		test = sin(lon2-lon1);
		in1 = find(test > 0);
		hdg(in1) = 360-hdg(in1);
	end

% ---------------------------------------------------------------------------
function [hFig,hAxes] = t_fig()
% Create a new figure
	hFig = ecran;
	handles = guidata(hFig);
	hAxes = handles.axes1;
	set(hFig, 'Tag','TidalPred')
% 	hFig = figure('NumberTitle','off','Visible','off','Color',get(0,'factoryUicontrolBackgroundColor'),...
% 		'Units','normalized','Position',[0.3 0.3 0.57 0.45],'Tag','TidalPred');
% 	hAxes = axes('position',[0.06, 0.07, 0.913, 0.909]);
	delete(findobj(hFig, 'Label','Analysis'))		% Better delete this menu where most options do not apply here
	delete(findobj(hFig, 'Tag', 'isocs_but'))		% And neither show this.
	options = uimenu(hFig,'Label','Options');
	uimenu(options,'Label','Calendar','callback',{@calendario,hFig});
	uimenu(options,'Label','Export tides','callback',{@export_mare,hFig});
	uimenu(options,'Label','Harmonics to Mat file','callback',@harmonics2mat);

% -------------------------------------------------------------------------------
function calendario(obj,eventdata,hFig)
	new_date = uisetdate;
	if (~isempty(new_date))
		yr = str2double(new_date(8:end));
		if (yr < 1970 || yr > 2037)
			errordlg('Tide prediction is not possible for this date.','Error');     return
		end
		tim = datenum(new_date);
		pt = getappdata(hFig,'CurrentPoint');
		t_xtide(pt(1),pt(2),tim);
	end

% -------------------------------------------------------------------------------
function export_mare(obj,eventdata,hFig)
	% Export the displayed tidal heights
	start_date = getappdata(hFig,'StartDate');
	xy = getappdata(hFig,'TideData');
	str1 = {'*.dat', 'Data file (*.dat)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = uiputfile(str1,'Tidal data file');
	if isequal(FileName,0);     return;     end
	% Open and write to ASCII file
	if ispc;        fid = fopen([PathName FileName],'wt');
	elseif isunix;  fid = fopen([PathName FileName],'w');
	else    errordlg('Unknown platform.','Error');
	end
	fprintf(fid,'# %s\n', start_date);
	fprintf(fid,'%9.5f\t%6.3f\n', xy);
	fclose(fid);

% -------------------------------------------------------------------------------
function harmonics2mat(obj,eventdata)
% Attempt to generate one mat-file from an xtide harmonics file....
% Latest version available from http://bel-marduk.unh.edu/xtide/files.html

	str1 = {'*.txt;*.dat', 'Harmonics file (*.txt,*.dat)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = uigetfile(str1,'Select Harmonics file');
	pause(0.01)
	if isequal(FileName,0);     return;     end
	fid = fopen([PathName,FileName],'r');

	fprintf('Reading harmonics file (this will take a while)\n');
	[xtide,xharm] = read_xtidefile(fid);

	[FileName,PathName] = uiputfile({'*.mat;*.MAT', 'Data files (*.mat,*.MAT)'},'Select Harmonics Mat-file');
	pause(0.01)
	if isequal(FileName,0);     return;     end
	h = msgbox('Saving harmonic information to t_xtide.mat');
	[PATH,FNAME,EXT] = fileparts([PathName FileName]);
	if isempty(EXT),	fname = [PathName FNAME '.mat'];
	else				fname = [PathName FNAME EXT];
	end
	save(fname,'xtide', 'xharm')
	delete(h)
  
% --------------------------------------------------------------------------
function [xtide,xharm] = read_xtidefile(fid)
% Reads the xtide harmonics file and creates a data structure
% with all that info for faster access

	l = fgetl_nocom(fid);       ncon = sscanf(l,'%d');

	xtide=struct('name',repmat(' ',ncon,8),'speed',zeros(ncon,1),...
			 'startyear',0,'equilibarg',zeros(ncon,68),'nodefactor',zeros(ncon,68));

	for (k=1:ncon)
		l = fgetl_nocom(fid);
		xtide.name(k,:) = l(1:8);
		xtide.speed(k) = sscanf(l(9:end),'%f');
	end

	xtide.startyear = sscanf(fgetl_nocom(fid),'%d');
	nyear = sscanf(fgetl_nocom(fid),'%d');

	for (k=1:ncon)
		fgetl(fid);
		xtide.equilibarg(k,:) = fscanf(fid,'%f',nyear);
		fgetl(fid);
	end
	fgetl(fid); % Skip *END*

	nyear=sscanf(fgetl_nocom(fid),'%d');

	for (k=1:ncon)
		fgetl(fid);
		xtide.nodefactor(k,:) = fscanf(fid,'%f',nyear);
		fgetl(fid);
	end
	l = fgetl(fid);		% Skip *END*

	% Now read in all harmonic data
	%nsta=1754; % This is number of stations in harmonics (1998-07-18)
	%nsta=3351; % This is number of stations in v1.42 or harmonics file
	nsta=3316; % This is number in v1.51

	xharm = struct('station',repmat(' ',nsta,79),'units',repmat(' ',nsta,8),...
			 'longitude',zeros(nsta,1),'latitude',zeros(nsta,1),...
			 'timezone',zeros(nsta,1),'datum',zeros(nsta,1),...
			 'A',zeros(nsta,ncon),'kappa',zeros(nsta,ncon));

	nh=0;
	while (~isempty(l) && l(1) ~= -1)
		l=[l '   '];
		nh=nh+1;
		while ~strcmp(l(1:3),'# !'),
			l=[fgetl(fid) '   '];
		end
		while strcmp(l(1:3),'# !'),
			switch l(4:7),
				case 'unit',
					tmp=deblank(l(strfind(l,':')+2:end));
					xharm.units(nh,1:length(tmp))=tmp;
				case 'long',
					xharm.longitude(nh)=sscanf(l(strfind(l,':')+1:end),'%f');
				case 'lati'  
					xharm.latitude(nh)=sscanf(l(strfind(l,':')+1:end),'%f');
			end
			l=fgetl(fid);
		end
		tmp = deblank(l);
		if (tmp(1) ~= '#')    % Not commented out
			xharm.station(nh,1:length(tmp)) = tmp;
			tmp = fgetl(fid);
			k = min(strfind(tmp,':'));
			tim = sscanf(tmp(1:k-1),'%d')+sscanf(tmp(k+(1:2)),'%d')/60;
			xharm.timezone(nh) = tim;
			xharm.datum(nh) = sscanf(fgetl(fid),'%f');

			for (k=1:ncon)
				l = fgetl(fid);
				if (l(1) ~= 'x')
					ll = min([strfind(' ',l) find(abs(l)==9)]); % space or tab
					tmp = sscanf(l(ll+1:end),'%f',2);
					xharm.A(nh,k) = tmp(1);
					xharm.kappa(nh,k) = tmp(2);
				end
			end
			l = fgetl(fid);
		else
			nh = nh - 1;  
		end
		if (rem(nh,50) == 0),	fprintf('.');   end
	end
	fprintf('\n');

	% Convert internally to sparse matrix storage (much smaller).
	xharm.A = sparse(xharm.A);
	xharm.kappa = sparse(xharm.kappa);
  
% --------------------------------------------------------------------------
function l = fgetl_nocom(fid)
% Gets a line that isn't a comment line
	l = fgetl(fid);
	while (~isempty(l) && l(1) == '#'),		l = fgetl(fid);  end

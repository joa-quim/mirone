function pole2neighbor(obj, evt, hLine)
% Compute finite poles from isochron N to isochron N+1 of the data/isochrons.dat file (must be loaded)
%
	if (nargin == 2),		hLine = gco;	end
	nNext = hLine;
	while (~isempty(nNext))
		nNext = compute_pole2neighbor(nNext);
	end

% -----------------------------------------------------------------------------------------------------------------
function hNext = compute_pole2neighbor(hLine)
% It inserts/updates the header line of current object (gco) or of HLINE
% HNEXT is the handle of nearest older isochron line

% $Id: $

	if (nargin == 0),		hLine = gco;	end % blabla
	D2R = pi/180;
	this_LineInfo = getappdata(hLine,'LineInfo');
	[this_isoc, r] = strtok(this_LineInfo);
	[this_plate_pair, r] = strtok(r);
	this_plate_pair_first_name = '';
	if (strncmp(this_plate_pair,'NOR',3) || strncmp(this_plate_pair,'SOU',3))		% NORTH or SOUTH is no good
		this_plate_pair_first_name = this_plate_pair;		% Either 'NORTH' or 'SOUTH'
		this_plate_pair = strtok(r);
	end
	this_isoc = get_isoc_numeric(this_isoc);

	p = parse_finite_pole(this_LineInfo);	% Get the Euler pole parameters of this isoc
	if (isempty(p))
		errordlg(sprintf('This line has no pole info. See:\n\n%s',this_LineInfo),'Error')
		hNext = [];		return
	end

	hAllIsocs = findobj('Tag', get(gco,'Tag'));
	closest_so_far = 5000;		% Just a pretend-to-be oldest anomaly (must be > 1000 + oldest M)
	ind_closest = [];
	for (i = 1:numel(hAllIsocs))
		LineInfo = getappdata(hAllIsocs(i),'LineInfo');
		if (isempty(LineInfo)),		continue,	end
		[isoc, r] = strtok(LineInfo);
		[plates, r] = strtok(r);
		if (strncmp(plates,'NOR',3) || strncmp(plates,'SOU',3))		% NORTH or SOUTH is no good, but wee need a further test
			if (~strcmp(this_plate_pair_first_name, plates)),	continue,	end		% A 'NORTH' and a 'SOUTH'
			plates = strtok(r);
		end
		if (~strcmp(this_plate_pair, plates)),	continue,	end		% This isoc is not in the same plate
		if (isoc(1) == 'A'),	continue,	end		% The ones named Axx are of no interest for the moment
		isoc = get_isoc_numeric(isoc);
		if (isoc <= this_isoc),	continue,	end		% We are only interested on older-than-this-isoc
		if (isoc < closest_so_far)
			closest_so_far = isoc;
			ind_closest = i;
		end
	end

	if (isempty(ind_closest))		% If nothing older is found, we are done
		hNext = [];
		return
	end

	% OK, so here we are supposed to have found the next older isoc
	closest_LineInfo = getappdata(hAllIsocs(ind_closest),'LineInfo');
	disp(closest_LineInfo)
	x = get(hLine,'Xdata');		y = get(hLine,'Ydata');
	x2 = get(hAllIsocs(ind_closest),'XData');
	y2 = get(hAllIsocs(ind_closest),'YData');
	n_pts = numel(x);
	inds = 1:fix(n_pts / 10):n_pts;		% Indices of one every other 10% of number of total points in working isoc

	xcFirst = [];		k = 0;
	% Find the first point in Working isoc hose circle about his Euler pole crosses its nearest (older) neighbor
	while (isempty(xcFirst))
		k = k + 1;
		if (k > numel(inds))
			errordlg(sprintf('Something screw.\n\n%s\n\nand\n\n%s\n\n are clearly not neighbors', ...
				this_LineInfo,closest_LineInfo),'Error')
			hNext = [];		return
		end
		indF = inds(k);
		c = sin(p.lat*D2R)*sin(y(indF)*D2R) + cos(p.lat*D2R)*cos(y(indF)*D2R)*cos( (p.lon-x(indF))*D2R );
		rad = acos(c) / D2R;
		[latc,lonc] = circ_geo(p.lat, p.lon, rad, [], 360);
		[lonc,I] = sort(lonc);	latc = latc(I);		% To avoid a false cross detection (due to circularity not taking into acount)
		[xcFirst,ycFirst] = intersections(lonc,latc,x2,y2);
	end

	% Now do the same but find last (but rounded to 10% N_pts) point in Working isoc that fills the same condition
	idr = numel(inds):-1:1;
	inds = inds(idr);			% Reversing oder makes our life easier
	xcLast = [];	k = 0;
	while (isempty(xcLast))
		k = k + 1;
		indL = inds(k);
		c = sin(p.lat*D2R)*sin(y(indL)*D2R) + cos(p.lat*D2R)*cos(y(indL)*D2R)*cos( (p.lon-x(indL))*D2R );
		rad = acos(c) / D2R;
		[latc,lonc] = circ_geo(p.lat, p.lon, rad, [], 360);
		[lonc,I] = sort(lonc);	latc = latc(I);
		[xcLast,ycLast] = intersections(lonc,latc,x2,y2);
	end

	% So at this point we should be able to get a first good estimate of the seeked rotation pole using Bonin Method
	[plon,plat,omega] = calc_bonin_euler_pole([x(indF); x(indL)],[y(indF); y(indL)],[xcFirst; xcLast],[ycFirst; ycLast]);
	%sprintf('Lon = %.1f  Lat = %.1f   Ang = %.3f', plon, plat, omega)

	% Get the neighbor pole parameters, though we will only use its age
	p_closest = parse_finite_pole(closest_LineInfo);	% Get the pole of neighbor because we need its age
	if (isempty(p_closest))
		errordlg(sprintf('This line has no pole info. See:\n\n%s',closest_LineInfo),'Error')
		hNext = [];		return
	end

	indS = strfind(this_LineInfo, 'STG1"');		% See if we aready have one INSITU STG
	if (isempty(indS))							% No, create one
		LineInfo = sprintf('%s STG1"%.1f %.1f %.1f %.2f %.3f"', this_LineInfo, plon, plat, p_closest.age, p.age, omega);
	else									% Yes, update it
		ind2 = strfind(this_LineInfo(indS+5:end),'"') + indS + 5 - 1;	% So that indS refers to string start too
		str = sprintf('%.1f %.1f %.1f %.2f %.3f', plon, plat, p_closest.age, p.age, omega);
		LineInfo = [this_LineInfo(1:indS+4) str this_LineInfo(ind2(1):end)];
	end
	setappdata(hLine,'LineInfo',LineInfo)	% Update this isoc with just computed STG pole

	if (nargout),	hNext = hAllIsocs(ind_closest);		end

% -----------------------------------------------------------------------------------------------------------------
function p = parse_finite_pole(LineInfo)
% Find the Euler pole parameters from the header line of an isochron
	indF = strfind(LineInfo, 'FIN"');
	if (isempty(indF)),		p = [];		return,		end
	ind2 = strfind(LineInfo(indF+4:end),'"') + indF + 4 - 1;	% So that ind refers to the begining string too
	[t, r] = strtok(LineInfo(indF(1)+4:ind2(1)-1));
	p.lon = sscanf(t, '%f');
	[t, r] = strtok(r);			p.lat = sscanf(t, '%f');
	[t, r] = strtok(r);			p.ang = sscanf(t, '%f');
	t      = strtok(r);			p.age = sscanf(t, '%f');

% -----------------------------------------------------------------------------------------------------------------
function isoc = get_isoc_numeric(isoc_str)
% Get the isoc name in numeric form. We need a fun because isoc names sometimes have letters
	if (isoc_str(end) == 'a' || isoc_str(end) == 'r')		% We have 4a, 33r, etc...
		isoc_str = [isoc_str(1:end-1) '.1'];	% Just to be older than the corresponding 'no a' version
	end
	if (isoc_str(1) == 'M'),	isoc_str = ['10' isoc_str(2:end)];	end		% The M's are older. Pretend it's 100, etc.
	isoc = str2double(isoc_str);


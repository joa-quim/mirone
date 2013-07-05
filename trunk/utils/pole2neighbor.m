function pole2neighbor(obj, evt, hLine, mode, opt)
% Compute finite poles from isochron N to isochron N+1 of the data/isochrons.dat file (must be loaded)

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

% $Id$

	if (isempty(hLine)),	hLine = gco;	end
	hNext = hLine;
	if (strcmpi(mode, 'Bonin'))
		while (~isempty(hNext))
			hNext = compute_pole2neighbor_Bonin(hNext);
		end
	else								% Best-fit (Brute-Force)
		if (nargin == 4)				% Make one full round on all poles in same plate of the seed
			nNewPoles = 0;
			while (~isempty(hNext))
				[hNext, is_pole_new] = compute_pole2neighbor_BruteForce(hNext);
				nNewPoles = nNewPoles + is_pole_new;
			end
			fprintf('Computed %d new poles', nNewPoles)
		else							% Compute the pole of the seed isoc only, iterating OPT times or not improving
			n = 1;		is_pole_new = true;
			while (n <= opt && is_pole_new)
				[hNext, is_pole_new] = compute_pole2neighbor_BruteForce(hLine);
				n = n + 1;
			end
		end
	end

% -----------------------------------------------------------------------------------------------------------------
function [hNext, is_pole_new] = compute_pole2neighbor_BruteForce(hLine)
% ...
	[hNext, pole, LineInfo, p, p_closest, CA, CB] = compute_pole2neighbor_Bonin(hLine);
	if (isempty(hNext))		% Suposedly the oldest isochron in its town
		is_pole_new = false;
		return
	end

	% OK, brute force now it is
	opt_D = '-D10/10/0.1';
	opt_I = '-I101/101/11';
	xA = get(hLine, 'XData');	yA = get(hLine, 'YData');
	xB = get(hNext, 'XData');	yB = get(hNext, 'YData');
	new_pole = compute_euler([xA(:) yA(:)], [xB(:) yB(:)], pole(1), pole(2), pole(3), opt_D, opt_I);
	% new_pole is [pLon, pLat, pAng, resid]
	[LineInfo, is_pole_new] = set_stg_info('bestfit', hLine, LineInfo, new_pole(1), new_pole(2), new_pole(3), ...
                                           p, p_closest, CA, CB, new_pole(4));

% -----------------------------------------------------------------------------------------------------------------
function [hNext, pole, LineInfo, p, p_closest, CA, CB] = compute_pole2neighbor_Bonin(hLine)
% Inserts/update the header line of current object (gco) or of HLINE
% HNEXT is the handle of nearest older isochron line
%
% However, if nargout > 1 return also the Bonin pole and the names of current and its
% closest old neighbor as strings Ca & CB.
	
	hNext = [];		pole = [];	LineInfo = '';		p_closest = [];		CB = '';

	if (nargin == 0),		hLine = gco;	end
	D2R = pi/180;
	this_LineInfo = getappdata(hLine,'LineInfo');
	[this_isoc, r] = strtok(this_LineInfo);
	[this_plate_pair, r] = strtok(r);
	this_plate_pair_first_name = '';
	if (strncmp(this_plate_pair,'NOR',3) || strncmp(this_plate_pair,'SOU',3))		% NORTH or SOUTH is no good
		this_plate_pair_first_name = this_plate_pair;		% Either 'NORTH' or 'SOUTH'
		this_plate_pair = strtok(r);
	end
	CA = this_isoc;								% A copy as a string for eventual use in Best-Fit info data
	this_isoc = get_isoc_numeric(this_isoc);	% Convert the isoc name into a numeric form

	p = parse_finite_pole(this_LineInfo);		% Get the Euler pole parameters of this isoc
	if (isempty(p))
		errordlg(sprintf('This line has no pole info. See:\n\n%s',this_LineInfo),'Error')
		return
	end

	hAllIsocs = findobj('Tag', get(gco,'Tag'));
	closest_so_far = 5000;			% Just a pretend-to-be oldest anomaly (must be > 1000 + oldest M)
	ind_closest = [];
	for (i = 1:numel(hAllIsocs))	% Loop to find the closest and older to currently active isoc
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
		isocN = get_isoc_numeric(isoc);				% Convert the isoc name into a numeric form
		if (isocN <= this_isoc),	continue,	end	% We are only interested on older-than-this-isoc
		if (isocN < closest_so_far)
			closest_so_far = isocN;					% Note: for example C3a is given an age = 3.1 to tell from C3
			CB = isoc;								% A copy as a string for eventual use in Best-Fit info data
			ind_closest = i;
		end
	end

	if (isempty(ind_closest))		% If nothing older is found, we are done
		return
	end

	% ------------- OK, so here we are supposed to have found the next older isoc -------------------
	hNext = hAllIsocs(ind_closest);
	
	if (nargout > 1)
		[plon,plat,omega] = get_last_pole(this_LineInfo);
		if (~isempty(plon))				% Remember, this is the case where we are going to do Brute-Force
			pole = [plon,plat,omega];	% So if a previous pole exists it's all we want, OTHERWISE we must
			LineInfo = this_LineInfo;	% proceed till the end of this fucntion to compute a Bonin pole.
			return
		end
	end

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
			return
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
		return
	end

	new_LineInfo = set_stg_info('bonin', hLine, this_LineInfo, plon, plat, omega, p, p_closest);

	if (nargout > 1)
		[plon,plat,omega] = get_last_pole(new_LineInfo);
		pole = [plon,plat,omega];
		LineInfo = new_LineInfo;
	end

% -----------------------------------------------------------------------------------------------------------------
function [new_LineInfo, is_pole_new] = set_stg_info(mode, hLine, LineInfo, plon, plat, omega, p, p_closest, CA, CB, residue)
% Update (or insert) the stage pole info in the isochron header data.
% MODE = 'bonin' insert/update the bonin type pole
% MODE = 'whatever' insert/update the Best-Fit type pole
% Actually the pole is a finite pole but due to how it was compute ... it is fact a stage pole (not convinced?)
%
% P & P_CLOSEST are structures with the finite poles of isochron CA and CB. Only P.AGE is used here
% CA, CB, RESIDUE are used only with the second generation Brute-Force (best-fit) poles and denote
% the isochron names the best fit residue
%
% On output NEW_LINEINFO holds the eventualy updated pole info (if a better best-fit pole was achieved)
%           IS_NEW_POLE -> true if the new best-fit pole has lower residue and hence used to update NEW_LINEINFO

	is_pole_new = false;
	new_LineInfo = LineInfo;

	if (strcmpi(mode,'bonin'))
		indS = strfind(LineInfo, 'STG1"');			% See if we aready have one INSITU STG
		if (isempty(indS))							% No, create one
			new_LineInfo = sprintf('%s STG1"%.1f %.1f %.1f %.2f %.3f"', LineInfo, plon, plat, p_closest.age, p.age, omega);
		else										% Yes, update it
			ind2 = strfind(LineInfo(indS+5:end),'"') + indS + 5 - 1;	% So that indS refers to string start too
			str = sprintf('%.1f %.1f %.1f %.2f %.3f', plon, plat, p_closest.age, p.age, omega);
			new_LineInfo = [LineInfo(1:indS+4) str LineInfo(ind2(1):end)];
		end
	else											% Set/update the Best-Fit pole info
		indS = strfind(LineInfo, 'STG');
		if (isempty(indS))
			errordlg(['This Isoc ', LineInfo, ' doesn''t have a Bonin STG, so it should not pass here.'],'WarnError')
			return
		elseif (numel(indS) > 1 && LineInfo(indS(1)+3) == ' ') % The annoying old STG not quoted case. Strip it.
			indS = indS(2:end);
		end

		% We don't always know the age of the older chron sent as argument (p_closest) so get it directly from LineInfo
		ind2 = strfind(LineInfo(indS(end):end), '"') + indS(end) - 1;
		[t, r] = strtok(LineInfo(ind2(1)+1:end));
		[t, r] = strtok(r);
		[t, r] = strtok(r);		ageB = t;
		t = strtok(r);			ageA = t;

		if (numel(indS) == 1)							% --- Only a Bonin type STG, insert the new Brute-Force one
			new_LineInfo = sprintf('%s STG2_%s-%s"%.1f %.1f %s %s %.3f %.3f"', ...
				LineInfo, CA, CB, plon, plat, ageB, ageA, omega, residue);
		else											% --- Already have a Brute-Force type STG, update it
			% the following gymnastic is because the STG string has the form STGXX_YY-ZZ"
			indS = indS(end);
			ind2 = ind2(2) - 1;		ind1 = ind2;
			while (LineInfo(ind1) ~= ' '),	ind1 = ind1 - 1;	end
			old_res = str2double(LineInfo(ind1:ind2));	% Get old residue
			if ( (old_res - residue) < 1e-3 )			% Old one is better, don't change it
				return
			end
			fprintf('Anom --> %s  Res antes = %.6f  Res depois = %.6f\n',CA, old_res, residue)
			%ind2 = strfind(LineInfo(indS:end), '_') + indS - 1;		% So it refers to start of string
			%old_ver = sscanf(LineInfo(indS+3:ind2-1),'%d');	% Old 'release' number (may grow as far as it gets better)
			str = sprintf('"%.1f %.1f %s %s %.3f %.3f"', plon, plat, ageB, ageA, omega, residue);
			new_LineInfo = [LineInfo(1:indS+2) sprintf('2_%s-%s', CA, CB) str];
		end	
		is_pole_new = true;
	end
	setappdata(hLine,'LineInfo',new_LineInfo)			% Update this isoc with just computed STG pole
		
% -----------------------------------------------------------------------------------------------------------------
function [plon,plat,ang] = get_last_pole(LineInfo)
% Get the finite pole parameters from the pseudo stage pole stored in STGXX
% By getting the last STG in LineInfo we are getting the most updated.

	plon = [];		plat = [];		ang = [];

	indS = strfind(LineInfo, 'STG');
	if (isempty(indS) || ((numel(indS) == 1) && (LineInfo(indS+3) == ' ')) )
		% An old STG not quoted case (to be deleted some time later)
		return
	end

	str = LineInfo(indS(end):end);		% This now has only the STGXX_??"....." string
	ind2 = strfind(str,'"') + 1;		% Add 1 to start after the quote
	[plon, r] = strtok(str(ind2(1):end));
	[plat, r] = strtok(r);
	[t, r] = strtok(r);
	[t, r] = strtok(r);
	ang = strtok(r);
	if (ang(end) == '"'),	ang(end) = [];	end		% STG1 doesn't have the last element as residue, hence the '"'
	plon = str2double(plon);
	plat = str2double(plat);
	ang = str2double(ang);

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


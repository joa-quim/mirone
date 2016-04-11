function [sun_params, lon, lat] = solar_params(lon_pt, lat_pt, TZ, termin, day, month, year, UTChour, UTCmin, UTCsec)
% Compute the day-night terminator and the civil, nautical and astronomical twilights
%
% [sun_params, lon, lat] = solar_params(lon_pt, lat_pt, TZ, termin, day, month, year, UTChour, UTCmin, UTCsec)
%	Inputs:
%		LON_PT,LAT_PT -> Location where data is reported
%		TZ            -> Time Zone
%		TERMIN        -> One of 'daylight', 'civil', 'nautical' or 'astronomical' for selecting which terminator
%		remainings    -> Optional date and time. If not provided compute them from 'now'
%
%	Outputs:
%		SUN_PARAMS    -> A structure with several parametrs os Solar data (Sunset, Sunrise, Noon, Equation of time, etc...)
%		LON,LAT       -> Terminator coordinates.
%
% Particular cases
% [...] = solar_params(lonp, latp)
%       Run computations for location LONP, LATP, time 'now' and Time Zone = 0 
% [...] = solar_params(lonp, latp, termin)
%       Run computations for location LONP, LATP, time 'now', Time Zone = 0 and terminator time TERMIN 
% [...] = solar_params(lonp, latp, TZ)
%       Run computations for location LONP, LATP, for the time 'now' and Time Zone = TZ
% [...] = solar_params(lonp, latp, TZ, termin)
%       Run computations for location LONP, LATP, time 'now', Time Zone = TZ and terminator time TERMIN 
%
% This function was partially inspired in 'plotdaynightterminator' found in FEX-54875
% 'Geostationary' package and, according to it, written by Mattia Rossi <mrossi@swin.edu.au>
% but the link to the original no longer exists. That code however, produced wrong results,
% so I used the Excell spreadshits in http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
% to verify the equations and added all computations in solar_calcs() as well as computations
% of the various daylight terminators, which is done with circ_geo();
%
% This function uses equinox information found at http://www.timeanddate.com/calendar/seasons.html
% and is set for the year 2016 in function get_hemisphere(). The consequence of not updating it is
% that terminators calculated in either the March or September equinox may be closed by the wrong pole.

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

% $Id$

	if (nargin == 0)		% Undocumented option. Defaults to Faro, Portugal.
		lon_pt = -7.92;		lat_pt = 37.073;
		TZ = 0;		termin = '';
	elseif (nargin == 2)
		TZ = 0;		termin = '';
	elseif (nargin == 3)
		if (isa(termin,'char'))
			TZ = 0;
		else
			termin = '';
		end
	end
	if (nargin <= 4)
		UTC_time = now - TZ / 24;
		year  = str2double(datestr(UTC_time,'yyyy'));
		month = str2double(datestr(UTC_time,'mm'));
		day   = str2double(datestr(UTC_time,'dd'));

		UTChour = str2double(datestr(UTC_time,'HH'));
		UTCmin  = str2double(datestr(UTC_time,'MM'));
		UTCsec  = str2double(datestr(UTC_time,'SS'));
	else
		error('Wrong number of arguments')
	end
	%day=26;UTChour = 17;	UTCmin = 49;	UTCsec = 30;

	D2R = pi / 180;
	UT = UTChour + UTCmin / 60 + UTCsec / 3600;
	
	% From http://scienceworld.wolfram.com/astronomy/JulianDate.html
	JD = 367 * year - floor(7 * ( year + floor((month + 9) / 12) ) / 4) - ...
		floor(3 * ( floor((year + (month - 9) / 7) / 100) + 1 ) / 4) + ...
		floor((275 * month) / 9) + day + 1721028.5 + (UT / 24);

	JC = (JD - 2451545) / 36525;		% number of Julian centuries since Jan 1, 2000, 12 UT
	L = 280.46645 + 36000.76983 * JC + 0.0003032 .* JC .* JC;		% Sun mean longitude, degree
	L = mod(L,360);		
	if (L < 0),		L = L + 360;	end
	M = 357.5291 + 35999.0503 * JC - 0.0001559 .* JC .* JC - 0.00000048 .* JC .* JC .* JC;		% mean anomaly, degrees
	M = mod(M,360);
	if (M < 0),		M = M + 360;	end

	C = (1.914602 - 0.004817*JC - 0.000014 .* JC .* JC) * sin(M * D2R);
	C = C + (0.019993 - 0.000101 * JC) * sin(2 * M * D2R) + 0.000289 * sin(3 * M * D2R);		% Sun ecliptic longitude
	theta = L + C;												% Sun true longitude, degree
	meanObliqEclipt =  23 + 26/60 + 21.448/3600 - (46.815*JC + 0.00059 .* JC .* JC - 0.001813 .* JC .* JC .* JC) / 3600;
	obliqCorr = meanObliqEclipt + 0.00256 * cos((125.04 - 1934.136 * JC) * D2R);	% Oblique correction
	lambda = theta - 0.00569 - 0.00478 * sin((125.04 - 1934.136*JC) * D2R);			% Sun apparent longitude
	SUNdec = asin(sin(obliqCorr * D2R) * sin(lambda * D2R)) / D2R;					% Sun declination

	sun_params = solar_calcs(JC, L, M, UT, SUNdec, TZ, obliqCorr, lon_pt, lat_pt, termin);

	if (nargout > 1)
		[lat,lon] = circ_geo(SUNdec, -(sun_params.HourAngle - lon_pt), sun_params.radius, [], 181);
		[mi, ind] = min(lon);
		lon = [lon(ind:end) lon(1:ind-1)];			% Force that longitudes always start at -pi or min lon.
		lat = [lat(ind:end) lat(1:ind-1)];
		lon = [lon 180 180 -180 -180 lon(1)];		% The extra points are to close the polygon cleanly around [-180 180]
		% Close the polygon (either from North or South, depending on the value of baseline)
		N_or_S = get_hemisphere(year, month, day, UTChour, UTCmin);
		lat = [lat lat(end) N_or_S N_or_S lat(1) lat(1)];
	end

% 	Sunrise_H = fix(sun_params.Sunrise * 24);
% 	Sunrise_M = (sun_params.Sunrise * 24 - Sunrise_H) * 60;
% 	Sunset_H  = fix(sun_params.Sunset * 24);
% 	Sunset_M  = (sun_params.Sunset * 24 - Sunset_H) * 60;

% -----------------------------------------------------------------------------------------------------------
function out = solar_calcs(JC, L, M, UT, SUNdec, TZ, obliqCorr, lon_pt, lat_pt, termin)

	% Type of terminator
	if (isempty(termin) || strncmp(termin, 'daylight', 3))
		radius = 90.833;
	elseif (strncmp(termin, 'civil', 3))
		radius = 90+6;
	elseif (strncmp(termin, 'nautical', 3))
		radius = 90+12;
	elseif (strncmp(termin, 'astronomical', 3))
		radius = 90+18;
	end

	D2R = pi/180;
	lat_pt = lat_pt * D2R;		% lon_pt doesn't need conversion
	SUNdec = SUNdec * D2R;
	M = M * D2R;
	L = L * D2R;

	var_y = tan(D2R*obliqCorr/2) * tan(D2R*obliqCorr/2);
	EEO   = 0.016708634 - JC .* (0.000042037 + 0.0000001267 * JC);		% Earth Eccentric Orbit
	out.EQ_time = 4*(var_y*sin(2*L)-2*EEO*sin(M)+4*EEO*var_y*sin(M)*cos(2*L)-0.5*var_y*var_y*sin(4*L)-1.25*EEO*EEO*sin(2*M)) / D2R;	% minutes
	HA_Sunrise  = (acos(cos(radius*D2R)/(cos(lat_pt) * cos(SUNdec)) - tan(lat_pt) * tan(SUNdec))) / D2R;
	out.SolarNoon = (720 - 4 * lon_pt - out.EQ_time + TZ * 60) / 1440;
	out.Sunrise = out.SolarNoon - HA_Sunrise * 4 / 1440;
	out.Sunset  = out.SolarNoon + HA_Sunrise * 4 / 1440;
	out.Sunlight_duration = 8 * HA_Sunrise;
	TrueSolarTime = mod(UT/24 * 1440 + out.EQ_time + 4 * lon_pt - 60 * TZ, 1440);
	if (TrueSolarTime < 0)
		out.HourAngle = TrueSolarTime / 4 + 180;
	else
		out.HourAngle = TrueSolarTime / 4 - 180;
	end
	SolarZenit = (acos(sin(lat_pt) * sin(SUNdec) + cos(lat_pt) * cos(SUNdec) * cos(out.HourAngle*D2R))) / D2R;
	SolarElevation = 90 - SolarZenit;

	if (SolarElevation > 85)
		r = 0;
	else
		if (SolarElevation > 5)
			r = 58.1 / tan(SolarElevation * D2R) -0.07 / (tan(SolarElevation*D2R)^3) + 0.000086 / (tan(SolarElevation*D2R)^5);
		else
			if (SolarElevation > -0.575)	% If > -34.5 arcmin (refraction index effect?)
				r = 1735 + SolarElevation * (-518.2 + SolarElevation * (103.4 + SolarElevation * (-12.79 + SolarElevation * 0.711)));
			else
				r = -20.772 / tan(SolarElevation*D2R);
			end
		end
	end

	refraction_index = r / 3600;
	out.SolarElevationCorrected = SolarElevation + refraction_index;

	sz = SolarZenit * D2R;		% Just to make below lines a bit simpler to read
	if (out.HourAngle > 0)
		out.SolarAzim = mod(acos(((sin(lat_pt)*cos(sz)) - sin(SUNdec))/(cos(lat_pt)*sin(sz))) / D2R + 180, 360);
	else
		out.SolarAzim = mod(540 - acos(((sin(lat_pt)*cos(sz)) - sin(SUNdec))/(cos(lat_pt)*sin(sz))) / D2R, 360);
	end
	out.SolarElevation = SolarElevation;
	out.radius = radius;

% ----------------------------------------------------------------------------------
function hemi = get_hemisphere(year, month, day, UTChour, UTCmin)
% Hemisphere calculations for YYYY equinox times. 
% What we need to find is if the terminator polygon will be closed from North or South pole
% 2016 equinox times: (20/03/2016 04:30 and 22/09/2016 15:21)
% 2017 equinox times: (20/03/2017 10:29 and 22/09/2017 21:02)
% 2018 equinox times: (20/03/2018 16:15 and 23/03/2018 02:54)
% 2019 equinox times: (20/03/2019 21:59 and 23/03/2019 08:50)
% 2020 equinox times: (20/03/2019 03:50 and 22/03/2019 14:31)
% TODO select the equinox times directly from the year variable (currently unused)

	Eq1_day = 20;		Eq1_hour =  4;		Eq1_min = 30;
	Eq2_day = 22;		Eq2_hour = 15;		Eq2_min = 21;

	if (month < 3 || month > 9 || day < Eq1_day || day > Eq2_day)
		hemi = 90;
		return
	elseif (month > 3 && month < 9 || day > Eq1_day || day < Eq2_day)
		hemi = -90;
		return
	end

	if ((month == 3) && (day == Eq1_day))		% March equinox
		if (UTChour < Eq1_hour)
			hemi = 90;
		elseif (UTChour > Eq1_hour)
			hemi = -90;
		else
			if (UTCmin <= Eq1_min),		hemi = 90;
			else						hemi = -90;
			end
		end
	elseif ((month == 9) && (day == Eq2_day))	% September equinox
		if (UTChour > Eq2_hour)
			hemi = 90;
		elseif (UTChour < Eq2_hour)
			hemi=-90;
		else		% UTChour == 15
			if (UTCmin >= Eq2_min),		hemi = 90;
			else						hemi = -90;
			end
		end
	end

function [rlon,rlat] = rot_euler(lon,lat,p_lon,p_lat,omega,units, ecc)
% ROT_EULER  Rotate points (or lines) about an Euler pole
%
%  [rlon,rlat] = ROT_EULER(lat,lon,p_lon,p_lat,omega)   # Two column vectors output
%   OR
%  rlon_rlat = ROT_EULER(lat,lon,p_lon,p_lat,omega)     # one Mx2 matrix output
%  P_LON, P_LAT & OMEGA are the obvious parameters that define the Euler pole
%  UNITS -> by default angles are expected to be in degrees. Use UNITS = 'radians'
%  if input AND output angles are in radians
%
%	GEOCENTRIC LATITUDES.
%	Convert to geocentric latitudes -> rotate -> convert back to geodetic latitudes
%	(..., units) UNITS	-> Numeric value. The numeric value contain the first exccentricity
%	(...,units, ecc)	-> Units selects either 'radians' or 'degrees' as above
%	If ECC = -1, use the default value for the WGS84 datum
%	If ECC = 0,  Do not convert to geocentrics
%
%   Based on code extracted from GMT

%	Copyright (c) 2004-2012 by J. Luis
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

	is_radians = false;		do_geocentric = false;
	if (nargin < 5)
		error('ROT_EULER: Incorrect number of arguments')
	elseif (nargin == 6)
		if (ischar(units))
			if (strncmp(units,'rad',3))
				is_radians = true;
			end
		else
			ecc = units;	% Excentricity
			do_geocentric = true;
			if (ecc < 0)
				ecc = 0.0818191908426215;	% WGS84
			elseif (ecc == 0)
				do_geocentric = false;		% No need to go for a useless conversion 
			end
		end
	elseif (nargin == 7)
		if (strncmp(units,'rad',3))
			is_radians = true;
		end
		do_geocentric = true;
		if (ecc < 0)
			ecc = 0.0818191908426215;		% WGS84
		elseif (ecc == 0)
			do_geocentric = false;			% No need to go for a useless conversion 
		end
	end

	D2R = pi/180;			R2D = 180 / pi;
	if (~is_radians)		% Angles are in degrees. We must convert them
		lon = lon(:) * D2R;			lat = lat(:) * D2R;
		p_lon = p_lon * D2R;		p_lat = p_lat * D2R;
		omega = omega * D2R;
	end

	if (do_geocentric)		% Convert to geocentric latitudes
		lat = atan2( (1-ecc^2)*sin(lat), cos(lat) );
		%p_lat = atan2( (1-ecc^2)*sin(p_lat), cos(p_lat) );
	end

	p_sin_lat = sin(p_lat);			p_cos_lat = cos(p_lat);
	s_lat = sin(lat);				c_lat = cos(lat);
	s_lon = sin(lon - p_lon);		c_lon = cos(lon - p_lon);
	cc = c_lat .* c_lon;

	tlon = atan2(c_lat .* s_lon, p_sin_lat * cc - p_cos_lat * s_lat);
	s_lat = p_sin_lat * s_lat + p_cos_lat * cc;
	c_lat = sqrt(1 - s_lat .* s_lat);

	s_lon = sin(tlon + omega);		c_lon = cos(tlon + omega);
	cc = c_lat .* c_lon;			clear tlon c_lon;

	rlat = asin(p_sin_lat * s_lat - p_cos_lat * cc);
	rlon = p_lon + atan2(c_lat .* s_lon, p_sin_lat * cc + p_cos_lat * s_lat);

	if (do_geocentric)		% Convert back to geodetic latitudes
		rlat = atan2( sin(rlat), (1-ecc^2)*cos(rlat) );
	end

	ind = (rlon > pi);
	if(any(ind)),	rlon(ind) = rlon(ind) - 2*pi;	end
	if (~is_radians)				% User wants angles in degrees
		rlon = rlon * R2D;			rlat = rlat * R2D;
	end

	if (nargout == 1),		rlon = [rlon rlat];		end

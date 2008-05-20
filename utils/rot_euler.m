function [rlon,rlat] = rot_euler(lon,lat,p_lon,p_lat,omega,units)
% ROT_EULER  Rotate points (or lines) about an Euler pole
%
%  [rlon,rlat] = ROT_EULER(lat,lon,p_lon,p_lat,omega)   # Two column vectors output
%   OR
%  rlon_rlat = ROT_EULER(lat,lon,p_lon,p_lat,omega)     # one Mx2 matrix output
%  P_LON, P_LAT & OMEGA are the obvious parameters that define the Euler pole
%  UNITS -> by default angles are expected to be in degrees. Use UNITS = 'radians'
%  if input AND output angles are in radians
%
%   Based on code extracted from GMT

%	Copyright (c) 2004-2006 by J. Luis
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

if (nargin < 5)
    error('ROT_EULER: Incorrect number of arguments')
elseif (nargin == 6 && strcmp(units,'radians'))
    is_radians = 1;
else
    is_radians = 0;
end

D2R = pi/180;			R2D = 180 / pi;
if (~is_radians)		% Angles are in degrees. We must convert them
	lon = lon(:) * D2R;			lat = lat(:) * D2R;
	p_lon = p_lon * D2R;		p_lat = p_lat * D2R;
	omega = omega * D2R;
end

p_sin_lat = sin(p_lat);			p_cos_lat = cos(p_lat);
s_lat = sin(lat);				c_lat = cos(lat);
s_lon = sin(lon - p_lon);		c_lon = cos(lon - p_lon);
clear lon lat;
cc = c_lat .* c_lon;

tlon = atan2(c_lat .* s_lon, p_sin_lat * cc - p_cos_lat * s_lat);
s_lat = p_sin_lat * s_lat + p_cos_lat * cc;
c_lat = sqrt(1 - s_lat .* s_lat);

s_lon = sin(tlon + omega);		c_lon = cos(tlon + omega);
cc = c_lat .* c_lon;			clear tlon c_lon;

rlat = asin(p_sin_lat * s_lat - p_cos_lat * cc);
rlon = p_lon + atan2(c_lat .* s_lon, p_sin_lat * cc + p_cos_lat * s_lat);

ind = (rlon > pi);				rlon(ind) = rlon(ind) - 2*pi;
if (~is_radians)				% User wants angles in degrees
	rlon = rlon * R2D;			rlat = rlat * R2D;
end

if (nargout == 1),		rlon = [rlon rlat];		end

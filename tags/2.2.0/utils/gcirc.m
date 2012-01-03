function [lon,lat] = gcirc(lon1,lat1,lon2,lat2,npts)
% Compute greate circle between (lon1,lat1) and (lon2,lat2)
% IF NPTS is not given, use a default of 100

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

	if (nargin < 4 || nargin > 5)
		error('Incorrect number of arguments')
	elseif (nargin == 4)
		npts = 100;
	end

	D2R = pi / 180;
	lat1 = lat1 * D2R;    lon1 = lon1 * D2R;
	lat2 = lat2 * D2R;    lon2 = lon2 * D2R;
	small = 1e-6*D2R;

	%  If a track starts at a pole, make sure it traverses the starting meridian
	indx = find(abs(lat1) >= pi/2-small);
	lon1(indx) = lon2(indx);

	%  Compute azimuth and distance
	azim  = azimuth_geo(lat1,lon1,lat2,lon2,'radians');
	tmp = sin(lat1).*sin(lat2) + cos(lat1).*cos(lat2).*cos(lon2-lon1);
	tmp(tmp >  1) =  1;     tmp(tmp < -1) = -1;     % Just in case
	teta = acos(tmp);

	teta = linspace(0,teta,npts);
	lat = repmat(lat1,1,npts);
	lon = repmat(lon1,1,npts);
	azim = repmat(azim,1,npts);

	% Vectors are nomally short sized, so we can aford a litle memory wasting in favor of spead.
	sinLat = sin(lat);      cosLat = cos(lat);
	sinTeta = sin(teta);    cosTeta = cos(teta);
	sinAzim = sin(azim);    cosAzim = cos(azim);
	lat = asin(sinLat.*cosTeta + cosLat.*sinTeta.*cosAzim) / D2R;
	lon = (lon + atan2(sinTeta.*sinAzim,cosLat.*cosTeta-sinLat.*sinTeta.*cosAzim)) / D2R;

function [lon_s,lat_s,ang_s] = add_poles(lon1, lat1, ang1, lon2, lat2, ang2, ecc)
% Add two finite rotation Euler poles.
%
% All angles are assumed to be in degrees
%	ECC is the excentricity for the computations using geocentric latitudes.
%	That is, if lats are in geodetic coords ECC is used to convert to geocentric Latitude.
%	If not provided, use the provided lats as they come in.
%	ECC = -1	 ->		use the default value for the WGS84 datum
%	ECC = 0		 ->		Ignore the to geocentric conversion request

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

	do_geocentric = false;
	if (nargin == 7)
		do_geocentric = true;
		if (ecc < 0)
			ecc = 0.0818191908426215;		% WGS84
		elseif (ecc == 0)
			do_geocentric = false;			% No need to go for a useless conversion 
		end
	end

	D2R = pi / 180;			R2D = 180 / pi;
	lon1 = lon1 * D2R;		lat1 = lat1 * D2R;		ang1 = ang1 * D2R;
	lon2 = lon2 * D2R;		lat2 = lat2 * D2R;		ang2 = ang2 * D2R;

	if (do_geocentric)		% Convert to geocentric latitudes
		lat1 = atan2( (1-ecc^2)*sin(lat1), cos(lat1) );
		lat2 = atan2( (1-ecc^2)*sin(lat2), cos(lat2) );
	end

	% Build the first rotation matrix
	[px,py,pz] = sph2cart(lon1,lat1,1);				% Convert first pole to cartesian
	num = 1.0 - cos(ang1);
	r1(1,1) = px * px * num + cos(ang1);
	r1(1,2) = px * py * num - pz * sin(ang1);
	r1(1,3) = px * pz * num + py * sin(ang1);

	r1(2,1) = py * px * num + pz * sin(ang1);
	r1(2,2) = py * py * num + cos(ang1);
	r1(2,3) = py * pz * num - px * sin(ang1);

	r1(3,1) = pz * px * num - py * sin(ang1);
	r1(3,2) = pz * py * num + px * sin(ang1);
	r1(3,3) = pz * pz * num + cos(ang1);


	% Build the second rotation matrix
	[px,py,pz] = sph2cart(lon2,lat2,1);				% Convert the second pole to cartesian
	num = 1.0 - cos(ang2);
	r2(1,1) = px * px * num + cos(ang2);
	r2(1,2) = px * py * num - pz * sin(ang2);
	r2(1,3) = px * pz * num + py * sin(ang2);

	r2(2,1) = py * px * num + pz * sin(ang2);
	r2(2,2) = py * py * num + cos(ang2);
	r2(2,3) = py * pz * num - px * sin(ang2);

	r2(3,1) = pz * px * num - py * sin(ang2);
	r2(3,2) = pz * py * num + px * sin(ang2);
	r2(3,3) = pz * pz * num + cos(ang2);

	% Compute the total rotation matrix T = R2 * R1
	T = r2 * r1;

	% Compute the Euler pole whose rotation matrix is T
	lon_s = atan2( (T(1,3) - T(3,1)), (T(3,2) - T(2,3)) ) * R2D;		% Lon
	tmp = sqrt( (T(3,2)-T(2,3))^2 + (T(1,3)-T(3,1))^2 + (T(2,1)-T(1,2))^2);
	lat_s = asin( (T(2,1) - T(1,2)) / tmp);								% Lat

	if (do_geocentric)		% Convert back to geodetic latitudes
		lat_s = atan2( sin(lat_s), (1-ecc^2)*cos(lat_s) ) * R2D;
	else
		lat_s = lat_s * R2D;
	end
	
	ang_s = atan(tmp / (T(1,1)+T(2,2)+T(3,3) -1 )) * R2D;

	% Make sure that 0 <= ang_s <= 180
	if (ang_s < 0),     ang_s = ang_s + 180;    end

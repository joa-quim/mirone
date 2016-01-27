function [latc,lonc] = circ_geo(lat, lon, rng, azim, np, noforce_pipi)
% Small circle defined by its center, range and azimuth
% circ_geo(lat,lon,rng,azim,np)
% azim and np are optional arguments. "azim" is a one or two-column vector. 
% For single column, returns the arc between 0 and azim. For two columns, returns
% the arc between azim(1) and azim(2).
% RNG is the radius of the circle. If it si a vector than compute as many circles as elements of RNG. 
% "np" specifies the number of output points [default = 180].
% NOFORCE_PIPI		If TRUE angles are not truncated into the [-pi pi] interval (default)
% All angles are in degrees.

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

% $Id: circ_geo.m 7759 2016-01-27 23:23:24Z j $

	D2R = pi/180;   npts  = 180;       az    = [];		n_args = nargin;
	if (n_args < 3)
		errordlg('Error calling circ_geo. Must give at least 3 arguments','Error')
	elseif (n_args == 3)
		noforce_pipi = false;
	elseif (n_args == 4)
		az = azim;					noforce_pipi = false;
	elseif (n_args == 5)
		az = azim;    npts = np;	noforce_pipi = false;
	else
		az = azim;    npts = np;
	end
	force_pipi = ~noforce_pipi;		% I'm fed up of tests on negations of negations to know if true

	%  Allow for multiple circles starting from the same point
	if (numel(lat) == 1 && numel(lon) == 1 && numel(rng) > 1)
		lat = lat(ones(size(rng)));		lon = lon(ones(size(rng)));
	end

	if (isempty(az))
		az = [0 360];			az = az(ones([size(lat,1) 1]), :);
	elseif (isequal(size(az), [1 2]))
		az = az(ones([size(lat,1) 1]), :);
	end
	if (numel(rng) == 1),		rng = ones(size(az,1), 1) * rng;	end

	% convert to radians
	lat = lat * D2R;    lon = lon * D2R;    az = az * D2R;    rng = rng * D2R;

	%  Expand the azimuth inputs
	if size(az,2) == 1              % Single column azimuth inputs
		az_neg = zeros(size(az));    az_pos = az;
	else                            % Two column azimuth inputs
		az_neg = az(:,1);            az_pos = az(:,2);
	end

	az = zeros([size(az_neg,1) npts]);
	for i = 1:size(az,1)
		% Handle case where limits give more than half of the circle. Go clockwise from start to end.
		if az_neg(i) > az_pos(i);     az_pos(i) = az_pos(i)+2*pi;	end
		az(i,:) = linspace(az_neg(i),az_pos(i),real(npts));	
	end

	% Each circle occupies a row of the output matrices.
	lat = lat(:,ones([1,size(az,2)]) );
	lon = lon(:,ones([1,size(az,2)]) );
	rng = rng(:,ones([1,size(az,2)]) );

	%  Ensure correct azimuths at either pole.
	epsilon = 1e-07;     % Set tolerance
	ind = (lat >= pi/2-epsilon);		az(ind) = pi;    % starting at north pole
	ind = (lat <= epsilon-pi/2);		az(ind) = 0;     % starting at south pole

	temp1  = sin(lat).*cos(rng);		temp2  = cos(lat).*sin(rng).*cos(az);
	latc   = asin(temp1+temp2);
	temp1  = sin(rng).*sin(az);			temp2  = cos(lat).*cos(rng);
	temp3  = sin(lat).*sin(rng).*cos(az);
	lonc   = lon + atan2(temp1,temp2-temp3);

	if (force_pipi)		% Truncate angles into the [-pi pi] range
		%lonc = pi*((abs(lonc)/pi) - 2*ceil(((abs(lonc)/pi)-1)/2)) .* sign(lonc);
		ind = (lonc > pi);
		if (any(ind)),		lonc(ind) = lonc(ind) - 2*pi;	end
		ind = (lonc < -pi);
		if (any(ind)),		lonc(ind) = lonc(ind) + 2*pi;	end
		[mi, ind] = min(lonc);
	else				% Truncate angles into the [0 2pi] range
		ind = (lonc > 2*pi);
		if (any(ind)),		lonc(ind) = lonc(ind) - 2*pi;	end
		ind = (lonc < 0);
		if (any(ind)),		lonc(ind) = lonc(ind) + 2*pi;	end
	end
	lonc = [lonc(ind:end) lonc(1:ind-1)];
	latc = [latc(ind:end) latc(1:ind-1)];

	%  Convert back to degrees
	latc = latc / D2R;      lonc = lonc / D2R;

function [latb,lonb] = buffer_j(lat, lon, dist, direction, npts, geog, outputformat)
%BUFFERM2 Computes buffer zone around a polygon
%
% [latb,lonb] = buffer_j(lat,lon,dist,direction)
% [latb,lonb] = buffer_j(lat,lon,dist,direction,npts)
% [latb,lonb] = buffer_j(lat,lon,dist,direction,npts,geog)
% [latb,lonb] = buffer_j(lat,lon,dist,direction,npts,geog,outputformat)
%
% This function is a replacement for the Mapping Toolbox function bufferm, 
% which calculates a buffer zone around a polygon. It uses the same concept as
% the original, constructing a buffer by placing rectangles around each edge
% and circles around each vertex of the input polygon.  However, it relies
% on the PolygonClip MEX library to calculate polygon unions/differences.
% The output is identical in format to that for bufferm but it has more inputs.
% This function is a modification of Kelly Kearney's BUFFERM2 which still
% depended on the Mapping Toolbox. This one doesn't.
%
% Input variables:
%
%   lon, lat        coordinates values defining the polygon to be buffered.
%                   This can be either a NaN-delimited vector, or a cell
%                   array containing individual polygonal contours.
%
%   dist:           Width of buffer, in degrees of arc along the surface
%
%   direction:      'in' or 'out' or 'both'. Default is 'both'
%
%   npts:           Number of points used to contruct the circles around
%                   each polygon vertex.  If omitted, default is 13.
%	geog 			= [a b], or [a f] uses geodetic computation, where a & b|f are ellipsoid params
%					= [] uses WGS-84
%					= 1  uses spherical aproximation
%					= 0  do cartesian calculation
%					If omitted, guesses if is geog from the X,Y coords (but never ellipsoidal)
%
%   outputformat:   'vector' (NaN-delimited vectors), 'cutvector'
%                   (NaN-clipped vectors with cuts connecting holes to the
%                   exterior of the polygon), or 'cell' (cell arrays in
%                   which each element of the cell array is a separate
%                   polygon), defining format of output.  If omitted,
%                   default is 'vector'.
%
% Output variables:
%
%   latb:           Latitude values for buffer polygon
%   lonb:           Longitude values for buffer polygon
%
% Example:
%
%   load conus
%   tol = 0.1; % Tolerance for simplifying polygon outlines
%   [reducedlat, reducedlon] = reducem(gtlakelat, gtlakelon, tol);
%   dist = 1;  % Buffer distance in degrees
%   [latb, lonb] = buffer_j(reducedlat, reducedlon, dist, 'out');
%   figure('Renderer','painters')
%   usamap({'MN','NY'})
%   geoshow(latb, lonb, 'DisplayType', 'polygon', 'FaceColor', 'yellow')
%   geoshow(gtlakelat, gtlakelon,...
%                       'DisplayType', 'polygon', 'FaceColor', 'blue')
%   geoshow(uslat, uslon)
%   geoshow(statelat, statelon)
%
% NOTE: Not all of the above was tested and there still misses adding a non-geog operability

%	Coffeeright (c) 2004-2012 by J. Luis
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

	% Check input
	nin = nargin;
	if (nin < 3 || nin > 7),	error('Wrong number of input arguments'),	end

	% Set defaults if not provided as input
	if (nin < 4)
		direction = 'both';		geog = guess_geog(lon, lat);
		npts = 13;				outputformat = 'vector';
	elseif (nin < 5)
		npts = 13;				geog = guess_geog(lon, lat);
		outputformat = 'vector';
	elseif (nin < 6)
		geog = guess_geog(lon, lat);	outputformat = 'vector';
	else
		outputformat = 'vector';
	end

	% Check format and dimensions of input
	if (direction(1) ~= 'i' && direction(1) ~= 'o' && direction(1) ~= 'b'),
		error('Direction must be either ''in'' or ''out'' or ''both'' '),
	end
	if (~isnumeric(dist) || numel(dist) > 1),	error('Distance must be a scalar.'),	end
	if (~isnumeric(npts) || numel(npts) > 1),	error('Number of points must be a scalar.'),	end

	if ( ~(strcmp(outputformat(1:2), 've') || strcmp(outputformat(1:2), 'cu') || strcmp(outputformat(1:2), 'ce')) )   
		error('Unrecognized output format flag.');
	end

	% Split (and test) NaN-polygon(s) into separate elements of cell arrays 
	if ~(isa(lat,'cell'))
		if (~isequal(size(lat), size(lon))),	error('Lat and lon must have identical dimensions'),	end
		[latcells, loncells] = map_funs('polysplit',lat(:), lon(:));
	else
		[lat, lon] = map_funs('polyjoin',lat, lon);		% In case multiple faces in one cell.
		lat = lat(:);	lon = lon(:);
		[latcells, loncells] = map_funs('polysplit',lat, lon);
		latsizes = cellfun(@size, lat, 'UniformOutput', false);
		lonsizes = cellfun(@size, lon, 'UniformOutput', false);
		if (~isequal(latsizes, lonsizes)),		error('Lat and lon must have identical dimensions'),	end
	end

	% Create buffer shapes on each face
	xout = cell(numel(latcells), 1);	yout = cell(numel(latcells), 1);
	for (ipoly = 1:numel(latcells))
		[xb, yb] = buffer(loncells{ipoly}, latcells{ipoly}, dist, npts, direction, geog);
		xout{ipoly} = xb;		yout{ipoly} = yb;
	end

	% Reformat output
	switch outputformat
        case 'vector'
			if (isa(yout,'cell'))
            	[latb, lonb] = map_funs('polyjoin',yout, xout);
			end
        case 'cutvector'
            [latb, lonb] = map_funs('polycut',yout, xout);
        case 'cell'
			latb = yout;	lonb = xout;
	end

% ------------------------------------------------------------------------------------------
function [xb, yb] = buffer(x, y, dist, npts, direction, geog)
% BUFFER handles one non-interrupted polygon at a time
% GEOG = [a b], or [a f] uses geodetic computation, where a & b are ellipsoid params
% GEOG = [] uses WGS-84
% GEOG = 1  uses spherical aproximation
% GEOG = 0  do cartesian calculation

	if ( numel(geog) == 2 || numel(geog) == 0 )
		fhandle_circ = @vreckon;
		fhandle_azim = @vdist;
		dist = dist * 111194.927;
	elseif (geog)
		fhandle_circ = @circ_geo;
		fhandle_azim = @azimuth_geo;
	else
		fhandle_circ = @circ_cart;
		fhandle_azim = @azim_cart;	
	end

	xb = [];	yb = [];
	perim = [0 360];
	range = dist;
 	[y_circ, x_circ] = feval(fhandle_circ, y, x, range, perim, npts, geog == 2);
	
	% Rectangles around each edge
	[dumb, az] = feval(fhandle_azim, y(1:end-1), x(1:end-1), y(2:end), x(2:end));
	[latbl1,lonbl1] = feval( fhandle_circ, y(1:end-1), x(1:end-1), range, az-90, 1, geog == 2);
	[latbr1,lonbr1] = feval( fhandle_circ, y(1:end-1), x(1:end-1), range, az+90, 1, geog == 2);
	[latbl2,lonbl2] = feval( fhandle_circ, y(2:end),   x(2:end),   range, az-90, 1, geog == 2);
	[latbr2,lonbr2] = feval( fhandle_circ, y(2:end),   x(2:end),   range, az+90, 1, geog == 2);

	y_rect = [latbl1 latbl2 latbr2 latbr1 latbl1]';
	x_rect = [lonbl1 lonbl2 lonbr2 lonbr1 lonbl1]';

	P1.x = x_circ(1,:);			P1.y = y_circ(1,:);	P1.hole = 0;
	P2.x = x_rect(:,1);			P2.y = y_rect(:,1);	P2.hole = 0;
	Pprev = PolygonClip(P1, P2, 3);				% Union of first circle and rectangle

	n_vertex = numel(y);
	for (i_vertex = 2:n_vertex-1)
		P1.x = x_circ(i_vertex,:);		P1.y = y_circ(i_vertex,:);
		P2.x = x_rect(:,i_vertex);		P2.y = y_rect(:,i_vertex);
		P3 = PolygonClip(P1, P2, 3);			% Union of circle and rectangle
		Pprev = PolygonClip(Pprev, P3, 3);		% Union of previous unions and new circle+rectangle union
	end
	P1.x = x_circ(n_vertex,:);		P1.y = y_circ(n_vertex,:);
	P3 = PolygonClip(Pprev, P1, 3);				% Union of previous unions and last circle

	% Now deal with any particular 'direction' request
	n_polygs = numel(P3);
	if (n_polygs == 1)
		if (direction(1) == 'b' || direction(1) == 'o')		% If only one, it will always be an outside
			xb = P3.x;		yb = P3.y;
		end
	else
		if (direction(1) == 'b')			% Want both 'in' and 'out'
			for (k = 1:n_polygs-1)
				xb = [xb; P3(k).x; NaN];		yb = [yb; P3(k).y; NaN];
			end
			xb = [xb; P3(k+1).x];				yb = [yb; P3(k+1).y];		% We didn't want the last NaN
		else								% Either 'out' or 'in'
			n_in_polygs = zeros(n_polygs,1);
			for (k = 1:n_polygs),	n_in_polygs(k) = numel(P3(k).x);	end	% Count np in each polygon
			[dumb, m] = max(n_in_polygs);
			if (direction(1) == 'o')		% 'out'
				xb = P3(m).x;		yb = P3(m).y;
			else							% 'in'
				u = 1:n_polygs;
				u(m) = [];					% Remove the larger polygon count (supposedly the outer polygon)
				for (k = u)
					xb = [xb; P3(k).x; NaN];		yb = [yb; P3(k).y; NaN];
				end
				xb(end) = [];						yb(end) = [];		% We dont't want the last NaN
			end
		end
	end

% ------------------------------------------------------------------------------------------
function [x_circ, y_circ] = circ_cart(x, y, range, perim, npts, dumb)
% DUMB is a variable never used except on the calling side (to make the nargin equal to the geog case)
	perim = perim * pi / 180;
	t = linspace(perim(1),perim(end),npts);
	x_circ = zeros(numel(x), npts);		y_circ = zeros(numel(x), npts);
	xc = range * cos(t);				yc = range * sin(t);
	for (k = 1:numel(x))
		x_circ(k,:) = x(k) + xc;		y_circ(k,:) = y(k) + yc;
	end

% ------------------------------------------------------------------------------------------
function [dumb, az] = azim_cart(x1, y1, x2, y2)
	dumb = [];
	dx = x2 - x1;   dy = y2 - y1;
	angs = atan2(dy,dx) * 180/pi;		% and convert to degrees
	az = (90 - angs);					% convert to azim (cw from north)
	ind = find(az < 0);
	az(ind) = 360 + az(ind);

% --------------------------------------------------------------------
function geog = guess_geog(x, y)
% Make a good guess if LIMS are geographic
	lims = [min(x) max(x) min(y) max(y)];
	geog = ( (lims(1) >= -180 && lims(2) <= 180) || (lims(1) >= 0 && lims(2) <= 360) )...
        && (lims(3) >= -90 && lims(4) <= 90);

function [latb,lonb] = buffer_j(lat,lon,dist,direction,npts,outputformat)
%BUFFERM2 Computes buffer zone around a polygon
%
% [latb,lonb] = buffer_j(lat,lon,dist,direction)
% [latb,lonb] = buffer_j(lat,lon,dist,direction,npts)
% [latb,lonb] = buffer_j(lat,lon,dist,direction,npts,outputformat)
%
% This function is a replacement for the Mapping Toolbox function bufferm, 
% which calculates a buffer zone around a polygon. It uses the same concept as
% the original, constructing a buffer by placing rectangles around each edge
% and circles around each vertex of the input polygon.  However, it relies
% on the PolygonClip MEX library to calculate polygon unions/differences.
% The input and output is identical in format to that for bufferm.
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
%   direction:      'in' or 'out'
%
%   npts:           Number of points used to contruct the circles around
%                   each polygon vertex.  If omitted, default is 13. 
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

% Coffeeright (c) 2002-2008 by J. Luis
% Joaquim Luis	27-Aug-2008

	% Check input
	nin = nargin;
	if (nin < 3 || nin > 6),	error('Wrong number of input arguments'),	end

	% Set defaults if not provided as input
	if (nin < 4),	direction = 'out';		end
	if (nin < 5)	npts = 13;				end
	if (nin < 6)	outputformat = 'vector';end

	% Check format and dimensions of input
	if (direction(1) ~= 'i' && direction(1) ~= 'o'),	error('Direction must be either ''in'' or ''out''.'),	end
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
		[xb, yb] = buffer(loncells{ipoly}, latcells{ipoly}, dist, npts, direction);
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
function [xb, yb] = buffer(x, y, dist, npts, direction)
% BUFFER handles one non-interrupted polygon at a time

	n_vertex = numel(y);
	range = zeros(n_vertex, 1) + dist;
 	[y_circ, x_circ] = circ_geo(y, x, range, [], npts);
	
	% Rectangles around each edge
	range(end) = [];
	az = azimuth_geo(y(1:end-1), x(1:end-1), y(2:end), x(2:end));
	[latbl1,lonbl1] = circ_geo(y(1:end-1), x(1:end-1), range, az-90, 1);
	[latbr1,lonbr1] = circ_geo(y(1:end-1), x(1:end-1), range, az+90, 1);
	[latbl2,lonbl2] = circ_geo(y(2:end),   x(2:end),   range, az-90, 1);
	[latbr2,lonbr2] = circ_geo(y(2:end),   x(2:end),   range, az+90, 1);

	y_rect = [latbl1 latbl2 latbr2 latbr1 latbl1]';
	x_rect = [lonbl1 lonbl2 lonbr2 lonbr1 lonbl1]';

	P1.x = x_circ(1,:);			P1.y = y_circ(1,:);	P1.hole = 0;
	P2.x = x_rect(:,1);			P2.y = y_rect(:,1);	P2.hole = 0;
	Pprev = PolygonClip(P1, P2, 3);				% Union of first circle and rectangle

	for (i_vertex = 2:n_vertex-1)
		P1.x = x_circ(i_vertex,:);		P1.y = y_circ(i_vertex,:);
		P2.x = x_rect(:,i_vertex);		P2.y = y_rect(:,i_vertex);
		P3 = PolygonClip(P1, P2, 3);			% Union of circle and rectangle
		Pprev = PolygonClip(Pprev, P3, 3);		% Union of previous unions and new circle+rectangle union
	end
	P1.x = x_circ(n_vertex,:);		P1.y = y_circ(n_vertex,:);
	P3 = PolygonClip(Pprev, P1, 3);				% Union of previous unions and last circle

	if (numel(P3) == 1)
		xb = P3.x;		yb = P3.y;
	else
		xb = [];		yb = [];
		for (k = 1:numel(P3)-1)
			xb = [xb; P3(k).x; NaN];		yb = [yb; P3(k).y; NaN];
		end
		xb = [xb; P3(k+1).x];				yb = [yb; P3(k+1).y];		% We didn't want the last NaN
	end

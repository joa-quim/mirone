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
%   lat:            Latitude values defining the polygon to be buffered.
%                   This can be either a NaN-delimited vector, or a cell
%                   array containing individual polygonal contours.
%                   External contours should be listed in a clockwise
%                   direction, and internal contours (holes) in a
%                   counterclockwise direction.
%
%   lon:            Longitude values defining the polygon to be buffered.
%                   Same format as lat. 
%
%   dist:           Width of buffer, in degrees of arc along the surface
%
%   direction:      'in' or 'out'
%
%   npts:           Number of points used to contruct the circles around
%                   each polygon vertex.  If omitted, default is 25. 
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

%---------------------------
% Check input
%---------------------------

nin = nargin;
if (nin < 3 || nin > 6)
    error('Wrong number of input arguments');
end

% Set defaults if not provided as input
if (nin < 4),	direction = 'out';		end
if (nin < 5)	npts = 25;				end
if (nin < 6)	outputformat = 'vector';end

% Check format and dimensions of input
if (direction(1) ~= 'i' && direction(1) ~= 'o'),	error('Direction must be either ''in'' or ''out''.'),	end
if (~isnumeric(dist) || numel(dist) > 1),	error('Distance must be a scalar.'),	end
if (~isnumeric(npts) || numel(npts) > 1),	error('Number of points must be a scalar.'),	end

if ( ~(strcmp(outputformat(1:2), 've') || strcmp(outputformat(1:2), 'cu') || strcmp(outputformat(1:2), 'ce')) )   
	error('Unrecognized output format flag.');
end

if iscell(lat)
	latsizes = cellfun(@size, lat, 'UniformOutput', false);
	lonsizes = cellfun(@size, lon, 'UniformOutput', false);
	if ~isequal(latsizes, lonsizes)
		error('Lat and lon must have identical dimensions');
	end
else
	if ~isequal(size(lat), size(lon))
		error('Lat and lon must have identical dimensions');
	end
end

%---------------------------
% Split polygon(s) into separate faces 
%---------------------------

if iscell(lat)
	[lat, lon] = map_funs('polyjoin',lat, lon);		% In case multiple faces in one cell.
end

lat = lat(:);
lon = lon(:);
[latcells, loncells] = map_funs('polysplit',lat, lon);

%---------------------------
% Create buffer shapes
%---------------------------

latcrall = cell(0);		loncrall = cell(0);

for (ipoly = 1:length(latcells))		% Circles around each vertex
    
	range = repmat(dist, size(latcells{ipoly}));
% 	[lattemp, lontemp] = scircle1(latcells{ipoly}, loncells{ipoly}, range, [], [], [], npts);
 	[lattemp, lontemp] = circ_geo(latcells{ipoly}, loncells{ipoly}, range, [], npts);
%	lattemp = lattemp';		lontemp = lontemp';

% 	latc = lattemp(:,1);    lonc = lontemp(:,1);
	latc = lattemp(1,:);    lonc = lontemp(1,:);
%	for (icirc = 2:size(lattemp,2))
	for (icirc = 2:size(lattemp,1))
		P1.x = lonc;	P1.y = latc;	P1.hole = 0;
% 		P2.x = lontemp(:,icirc);		P2.y = lattemp(:,icirc);	P2.hole = 0;
		P2.x = lontemp(icirc,:);		P2.y = lattemp(icirc,:);	P2.hole = 0;
		P3 = PolygonClip(P1, P2, 3);
		lonc = [P3(1).x; P3(1).x(1)];	latc = [P3(1).y;  P3(1).y(1)];		% Add first pt to close (bug in PolygonClip?)
		if (numel(P3) > 1)			% The two circles did not intercept
			lonc = [lonc; P3(2).x; P3(2).x(1)];
			latc = [latc; P3(2).y; P3(2).y(1)];
		end
	end 
    
	% Rectangles around each edge
	range(end) = [];
	az = azimuth_geo(latcells{ipoly}(1:end-1), loncells{ipoly}(1:end-1), latcells{ipoly}(2:end), loncells{ipoly}(2:end));
% 	az = azimuth(latcells{ipoly}(1:end-1), loncells{ipoly}(1:end-1), latcells{ipoly}(2:end), loncells{ipoly}(2:end));
% 	[latbl1,lonbl1] = reckon(latcells{ipoly}(1:end-1), loncells{ipoly}(1:end-1), range, az-90);
% 	[latbr1,lonbr1] = reckon(latcells{ipoly}(1:end-1), loncells{ipoly}(1:end-1), range, az+90);
% 	[latbl2,lonbl2] = reckon(latcells{ipoly}(2:end),   loncells{ipoly}(2:end),   range, az-90);
% 	[latbr2,lonbr2] = reckon(latcells{ipoly}(2:end),   loncells{ipoly}(2:end),   range, az+90);

	[latbl1,lonbl1] = circ_geo(latcells{ipoly}(1:end-1), loncells{ipoly}(1:end-1), range, az-90, 1);
	[latbr1,lonbr1] = circ_geo(latcells{ipoly}(1:end-1), loncells{ipoly}(1:end-1), range, az+90, 1);
	[latbl2,lonbl2] = circ_geo(latcells{ipoly}(2:end),   loncells{ipoly}(2:end),   range, az-90, 1);
	[latbr2,lonbr2] = circ_geo(latcells{ipoly}(2:end),   loncells{ipoly}(2:end),   range, az+90, 1);

    lattemp = [latbl1 latbl2 latbr2 latbr1 latbl1]';
    lontemp = [lonbl1 lonbl2 lonbr2 lonbr1 lonbl1]';
    
	latr = lattemp(:,1);    lonr = lontemp(:,1);
	c = 0;
    for (irect = 2:size(lattemp,2))
		P1.x = lonr;	P1.y = latr;	P1.hole = 0;
		P2.x = lontemp(:,irect);		P2.y = lattemp(:,irect);	P2.hole = 0;
		P3 = PolygonClip(P1, P2, 3);
		lonr = P3(1).x;		latr = P3(1).y;			% WHAT IS IN THE P3(I>1) ?????
		N = numel(P3);
		if (N > 1)
			c = c + 1;
		end
    end

	lonb = lonc;	latb = latc;
    % Union of circles and rectangles
    
	P1.x = lonr;	P1.y = latr;	P1.hole = 0;
	P2.x = lonc;	P2.y = latc;	P2.hole = 0;
	P3 = PolygonClip(P1, P2, 3);
	loncr = P3(1).x;		latcr = P3(1).y;
	if (ipoly == 1)
		P_1.x = loncr;		P_1.y = latcr;	P_1.hole = 0;
	else
		P2.x = loncr;	P2.y = latcr;	P2.hole = 0;
		P3 = PolygonClip(P_1, P2, 3);
		P_1 = P3;
	end
    
end

	N = numel(P3);
	if (N == 1)
		loncrall = P3.x;		latcrall = P3.y;
		if ( latcrall(1) ~= latcrall(end) || loncrall(1) ~= loncrall(end) )		% open polygon -> BUG somewhere
			loncrall = [loncrall; loncrall(1)];
			latcrall = [latcrall; latcrall(1)];
		end
	else
		loncrall = [];		latcrall = [];
		for (k = 1:N)
			loncrall = [loncrall; P3(k).x; NaN];
			latcrall = [latcrall; P3(k).y; NaN];
		end
		loncrall(end) = [];
		latcrall(end) = [];
	end

%---------------------------
% Calculate union/difference
%---------------------------

% switch direction
%     case 'out'
% 		[lonb, latb] = polybool('+', loncells, latcells, loncrall, latcrall);
%     case 'in'
% 		[lonb, latb] = polybool('-', loncells, latcells, loncrall, latcrall);
% end
lonb = loncrall;
latb = latcrall;

%---------------------------
% Reformat output
%---------------------------

switch outputformat
    case 'vector'
		if (isa(latb,'cell'))
        	[latb, lonb] = map_funs('polyjoin',latb, lonb);
		end
    case 'cutvector'
        [latb, lonb] = map_funs('polycut',latb, lonb);
    case 'cell'
end

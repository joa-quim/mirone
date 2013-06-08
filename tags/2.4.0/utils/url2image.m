function  varargout = url2image(opt, varargin)
% Mosaic image tiles from the Microsoft Virtual Earth web map servers
%
%	[img, lon_mm, lat_mm, x, y] = url2image('tile2img', lon, lat, zoom, PN/PV)
%			LON, LAT , ZOOM are mandatory. PN, PV are property name/value pairs
%			If LON = [x_min x_max] & LAT = [y_min y_max] build a 2D MOSAIC (see below)
%
%		url2image(...,'mosaic',mosaic).
%			If MOSAIC == [] process only the tile that contain the lon, lat coords
%			IF MOSAIC is a scalar, process a mosaic of size ones(mosaic,mosaic) 
%			IF MOSAIC is a 2D matrix, process a mosaic of size ones(size(mosaic)) 
%
%		url2image(...,'source', server).
%			The Microsoft image servers are used by default. However, there are others, namely the Google's
%			The problem with Google images is that they are not free to use (see their terms of use)
%			Since the Mars images are free, we'll use that as an example
%			SERVER is cell array with two fields
%				SERVER{1} = url-to-the-map-server (for example for Mars: mw1.google.com/mw-planetary/mars/elevation/t )
%				SERVSR{2} = the key characters used to code the tiles names ( example MS '0132' or Google's 'qrst' )
%
%		url2image(...,'what', whatKind).
%			whatKind is any of 'aerial', 'road' or 'hybrid'
%			Note: this option works only with Virtual Earth servers
%
%		url2image(...,'lonlat', 'yes')
%			If OPT == 'callmir' || OPT == 'url2image' returns coordinates in geogs instead of Mercator.
%			Ignored otherwise. NOTICE that image is still mercator so the lat coords are correct but they
%			don't grow linearly from lat_min to lat_max
%
%		url2image(...,'verbose', 'yes')
%			print out info while downloading the image files (silent when geting file from cache)
%			If value is 'all' instead of 'yes', prints info either when retrieving from cache or downloading
%
%		url2image(...,'cache', cache).
%			Search in the CACHE directory for tile files before trying to download them
%
%		VERY IMPORTANT INFO about CACHE.
%			CACHE dir can be provided into two diferent forms:
%				1) NASA World Wind type cache. There is a plugin for WW that allows downloading of Virtual Earth tiles.
%					In such cases, tiles are stored in a directory tree typicaly like:
%					C:\Program Files\NASA\World Wind 1.4\Cache\Virtual Earth
%					So, if you have WW installed with the VE plugin, give the above adress to be used as a search
%					path to VE files before trying to dowload them.
%				2) C:\whatever\you\want\to\call\this\path\cache
%					In this case 3 further subdirectories are appended to the CACHE dir (see also below) to acomodate
%					the 'what' propertie.
%						IF 'what' == 'aerial' cache = [cache 'kh'];
%						IF 'what' == 'road'   cache = [cache 'mt'];
%						IF 'what' == 'hybrid' cache = [cache 'tt'];
%					Furthermore, a subsequent directory is still appended based on the ZOOM level required.
%					Example of an absolute dir of a 12 zoom level aerial request: 
%						C:\lixo\Cache\kh\12
%
%				CACHE summary. CACHE is a base name directory of which subdirectories are assumed to exist (or
%				created if they do not) to hold proguessive refinement zoom level images.
%				Please follow EXACTLY one of the two possible forms as explained above. ... Otherwise cache is ignored. 
%
%	[img, lon_mm, lat_mm, x, y] = url2image('tile2img',...)
%			IMG is the mosaiced image
%			LON_MM, LAT_MM, -> min/max of the region long, lat. 
%			X, Y are the spherical mercator tile/ensemble of tiles coordinates
%
%	[url, lon_mm, lat_mm, x, y] = url2image('tile2url', lon, lat, zoom, PN/PV)
%			This form works pretty much like 'tile2img' except that it returns tile(s) url(s)
%			The only meaning options for PN/PV are (...,'mosaic',mosaic).
%			Or, to a less extent (see below), (...,'quadonly',whatever)
%			If MOSAIC implyies more than one tile (see above), URL is a cell array, otherwise it's a char string
%
%	[quad, lon_mm, lat_mm] = url2image('tile2url', lon, lat, zoom,'quadonly', 'whatever')
%			QUAD contains only the quadtree coding name, that is without the prepending server adress (e.g, 'a121')
%			This is usefull in conjunction with a posterior call like [lims, tiles_bb]  = url2image('quadcoord', quad);
%			that will compute the limits of a set of quadtrees names.
%
%	[lims, tiles_bb] = url2image('quadcoord', quad)
%			Compute the coordinates of a QUAD quadtree string, url, or cell array of one the previous
%			LIMS contains the BoundingBox coords in geogs
%			TILES_BB contains a Mx4 matrix with rectangle BB coords of each tile [lonMin lonMax latMin latMax]
%			This matrix is written by columns. That is like coords(quad{:})
%
%	url2image('mappcache', quad, cache, ext)
% 		Writes a mapping file of the files listed by QUAD in the CACHE directory.
%		If QUAD == [] the CACHE dir is scanned for *.EXT files. If EXT is not provided, it defaults to .jpg
%		The mapping file has a fix name of tilesMapping.mat and contains three variables
% 			region -> all points bounding box
%			zoomL ->  zoom level
% 			tiles_midpt -> a Mx2 matrix with tiles mid point coordinates
%
%
%	WARNING: If you access the internet via a proxy you might need to declare this environmental variable
%			set http_proxy=proxy_adress:port
%
%	EXAMPLES:
%	---------
%			- Get the whole world (well, not higher than +- 85 degrees N/S)
%		[img, lon, lat] = url2image('tile2img',[-180 180],[-80 80], 4);
%
%			- Get a 3x3 mosaic of level 19 over a part of l'ile de France in Paris
%		[img, lon, lat, x, y] = url2image('tile2img',2+20.869/60,48+51.29/60,19,'mosaic',3,'lonlat', 'yes');
%
%			- Sow the image in its original Mercator projection
%		image(x,y,img), axis xy, axis image
%
%			- Or for Mirone usage
%		url2image('callmir',-8,37.014,14,'mosaic',3);
%		url2image('callmir',[-180 180],[-80 80], 4);
%
%
%   AUTHOR
%       Joaquim Luis (jluis@ualg.pt)   01-Feb-2008

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

	quadkey = {'0' '1'; '2' '3'};				% Default to Virtual Earth
	prefix = 'http://a0.ortho.tiles.virtualearth.net/tiles/a';
	
	if (~nargin),		return,		end
	% ---------------------------- Test args ----------------------------
	varargs = [];		cache = 'D:\lixo\GoogleMV2\GoogleMV\Cache';		%cache = [];
	if ( strcmp(opt(1:3), 'til') || strcmp(opt(1:3), 'cal') )
		if (nargin < 4)
			error('url2image: minimum number of argins is 4')
		end
		
		lon = varargin{1};		lat = varargin{2};		zoomL = varargin{3};
		if ( ~isnumeric(lon) || ~isnumeric(lat) || ~isnumeric(zoomL) )
			error('url2image: First 3 arguments must be numeric (lon, lat, zoomlevel)')
		end
	
		if ( numel(varargin) > 3 ),		varargs = varargin(4:end);		end
	elseif ( strcmp(opt(1:3), 'qua') )		% 'quadcoord'
		if (numel(varargin) < 1)
			error('url2image: minimum number of argins is 1')
		else
			varargs = varargin(2:end);		% To use if needed to change 'quadkey'
		end
	end
	if ( ~isempty(varargs) && rem(numel(varargs), 2) )
		error('url2image:tile2url', 'Property names - property values must came in pairs')
	end

	% -------------------------------------------------------------------------

	if (~isempty(varargs))
		for (k = 1:2:numel(varargs))
			if ( strncmp(varargs{k}, 'cache', 2) )
				cache = varargs{k+1};
				if (~isempty(cache) && cache(end) == filesep),		cache(end) = [];	end		% No '/' at the end
			elseif ( strncmp(varargs{k}, 'source', 3) )
				dumb = varargs{k+1};		msg = [];
				if (~isa(dumb, 'cell'))
					msg = '"source" property value must be a cell array';
				elseif (numel(dumb) ~= 2)
					msg = '"source" property value must have two elements';
				end
				if (~isempty(msg)),		error('url2image:tile2url', msg),	end
				prefix  = varargs{k+1}{1};
				if (~strcmp(prefix(1:7),'http://')),	prefix = ['http://' prefix];	end
				quadkey = varargs{k+1}{2};
				if ( isa(quadkey, 'char') )			% Need to conver to the cell array format
					if (numel(quadkey) ~= 4)
						error('url2image:tile2url', '"quadkee" keeyword must contain 4 characters')
					end
					dumb = quadkey;
					quadkey = {dumb(1) dumb(2); dumb(4) dumb(3)};
				elseif ( isa(quadkey, 'cell') )
					if (numel(quadkey) ~= 4)
						error('url2image:tile2url', '"quadkee" keeyword must contain 4 characters')
					end
				else
					error('url2image:tile2url', 'Wrong type for the "quadkee" keeyword')
				end
			elseif ( strncmp(varargs{k}, 'proxy', 2) )
				if (strncmp(computer,'PC',2))
					dos(['set http_proxy=' varargs{k+1}])
				else
					unix(['set http_proxy=' varargs{k+1}])
				end
			end
		end
	else
		varargs = {''};			% So that varargs{:} below can work peacefully
	end
	% -------------------------------------------------------------------------
	
	%geoid = [6378137, 1/298.2572235630];	% WGS84
	%geoid = [6371008.7714, 0];				% Spherical
	geoid = [6378137, 0];					% Spherical - I don't understand why the WGS82 major axis is used to specify a spherical Earth
											%				but it seams that everybody is using this (e.g. EPSG:3785)
	switch lower(opt)
		case {'tile2url' 'tile2img'}
			[varargout{1:nargout}] = tile2url(opt, geoid, quadkey, prefix, lon, lat, zoomL, cache, varargs{:});
		case 'callmir'					% After composing the image, feed it into Mirone
			tile2url(opt, geoid, quadkey, prefix, lon, lat, zoomL, cache, varargs{:});
		case 'quadcoord'
			switch nargout
				case 1
					varargout{1} = quadcoord(geoid(2), quadkey, varargin{1});
				case 2
					[varargout{1:nargout}] = quadcoord(geoid(2), quadkey, varargin{1});
				otherwise
					error('url2image:quadcoord', 'Wrong number (min 1, max 2) of output args')
			end
		case 'mappcache'				% Call quadcoord() and save a mappimg of the existing tiles in the 'cache' dir
			mappcache(geoid(2), quadkey, varargin{:})
		otherwise
			disp('Valid keewords are: ''tile2url'' ''callmir'' ''quadcoord'' ''mappcache'' ')
			error('url2image: Unknown operating option')
	end

% --------------------------------------------------------------------------------------	
function [url, lon_mm, lat_mm, x, y] = tile2url(opt, geoid, quadkey, prefix, lon, lat, zoomL, cache, varargin)
% Compute the URL of the tile with the given coordinates.

	eq_radius = geoid(1);		flatness = geoid(2);

	% ---------------------------- Test args ----------------------------
	n_varargin = numel(varargin);
	mosaic = [];		verbose = false;	whatKind = 'a';		reportMerc = true;		quadonly = false;

	for (k = 1:2:n_varargin)
		if ( strncmp(varargin{k}, 'mosaic', 1) )
			mosaic = varargin{k+1};
		elseif ( strncmp(varargin{k}, 'what', 1) )
			whatKind = varargin{k+1}(1);
			if (quadkey{1} == '0')			% Virtual Earth
				%if ( isempty(strfind('r',whatKind)) && isempty(strfind('h',whatKind)) ),	whatKind = 'a';	continue,	end		% Was error
				if ( (whatKind ~= 'r') && (whatKind ~= 'h') ),		whatKind = 'a';		continue,	end		% Was error
				prefix(8) = whatKind;			prefix(end) = whatKind;
			end
		elseif ( strncmp(varargin{k}, 'lonlat', 1) )
			reportMerc = false;
		elseif ( strncmp(varargin{k}, 'verbose', 1) )
			verbose = true;
			if (varargin{k+1}(1) == 'a'),		verbose = 2;	end		% Speak up also when reading file from cache
		elseif ( strncmp(varargin{k}, 'quadonly', 1) )		% NOT TO BE USED WITH IMGs
			quadonly = true;
		end
	end

	if ( ~isequal(size(lon), size(lat)) )
		error('lon & lat must be of the same size')
	end

	% ----------------------- Test for well behaved coords values.
	if (any(lat > 85)),			lat(lat > 85)  = 85;		end 
	if (any(lat < -85)),		lat(lat < -85) = -85;		end
	
	ind = (lon < -180);
	if (any(ind))
		while (lon(ind) < -180),	lon(ind) = lon(ind) + 360;	end
	end
	ind = (lon > 180);
	if (any(ind))
		while (lon(ind) > 180),		lon(ind) = lon(ind) - 360;	end
	end
		
	if (numel(mosaic) == 1),	mosaic = ones(mosaic);			% "mosaic" is a scalar
	else						mosaic = ones(size(mosaic));	% "mosaic" is a matrix
	end

	lat = geod2isometric(lat, flatness);
	[x, y, xmm, ymm] = getPixel(lon, lat, zoomL);	% x,y are the fractional number of the 256 bins counting from origin
	x = floor(x);		y = floor(y);
	[lon_mm, latiso_mm] = getLonLat(xmm, ymm, zoomL+8);
	latiso_mm = flipud(latiso_mm(:,end:-1:1));		% Revert because origin was at top -> down.	flipud to have south at 1 row and growing north with rows
	lat_mm = isometric2geod(latiso_mm, flatness);

	% ---------------------- Case when rectangle BB was given as coords
	if (numel(lon) == 2)
		% See if have other tiles in between the ones deffined by lon_min and lon_max 
		Dtile = xmm(2,1) - xmm(1,2);			% > 0, have tiles in the midle; == 0, two contiguous tiles; < 0, same tile
		if (Dtile > 0),			nInTilesX = Dtile / 256;
		elseif (Dtile == 0),	nInTilesX = 0;
		else					nInTilesX = -1;
		end
		Dtile = ymm(2,1) - ymm(1,2);			% Idem for lat
		if (Dtile > 0),			nInTilesY = Dtile / 256;
		elseif (Dtile == 0),	nInTilesY = 0;
		else					nInTilesY = -1;
		end
		mosaic = ones(nInTilesY+2, nInTilesX+2);		% Create the mosaic matrix

		% lon_mm, lat_mm & latiso_mm have > 2 elements. Reduce them to the [min max] form
		lon_mm = [min(lon_mm(:)) max(lon_mm(:))];
		lat_mm = [min(lat_mm(:)) max(lat_mm(:))];
		latiso_mm = [min(latiso_mm(:)) max(latiso_mm(:))];

		% Calculate center so that the below MxN mosaic code case can be reused
		x = fix((x(1) + x(2)) / 2);			y = fix((y(1) + y(2)) / 2);
	end

	% ---------------------- CORE THING ---- Calculate the quadtree string
	quad = repmat(' ', 1, zoomL-1);
	for (i = 1:(zoomL-1))
		[x, rx] = rema(x, 2);		rx = rx+1;
		[y, ry] = rema(y, 2);		ry = ry+1;
		quad(i) = quadkey{ry,rx};
	end
	quad = fliplr(quad);
	
	% ----------------------- Find the tile decimal adress. x counts from 0 at -180 to 2^(zoomL - 1) - 1 at +180 
	decimal_adress = getQuadLims([prefix(end) quad], quadkey, 1);		% prefix(end) == 't' || 'a'
	
	isMapForFree = false;
	pref_bak = prefix;
	if ( quadkey{1} ~= '0' && ~isempty(strfind('rh',whatKind)) )		% hiden cat with outside tail
		prefix = [prefix sprintf('%d&y=%d&z=%d',decimal_adress,zoomL-1)];
	elseif ( quadkey{1} == 'm' && quadkey{2} == 'e' && quadkey{3} == 'f' && quadkey{4} == 'r')
		isMapForFree = true;
		prefix = [prefix sprintf( '%d/row%d/%d_%d-%d.jpg',zoomL-1,decimal_adress(2),zoomL-1,decimal_adress(1),decimal_adress(2) )];
	end

	if (quadonly)		% NOT TO BE USED WITH IMGs. The character used is irrelevant
		prefix = prefix(end);
		quadkey{1} = 'z';		% Dumb attribution but it forces the go trough the NOT == '0' branches bellow
	end
	
	if (isempty(mosaic))				% One single tile
		url = prefix;
		quad_ = quad;					% A copy
		if (quadkey{1} == '0'),			url = [url quad '?g=244'];		% VE
		elseif (whatKind(1) == 'a'),	url = [url quad];
		end

	else
		% Mosaic is a MxN matrix
		[mMo, nMo] = size(mosaic);
		url = cell(mMo,nMo);	quad_ = cell(mMo,nMo);		% quad_ is a copy with the quadtree string only
		mc = ceil(mMo/2);		nc = ceil(nMo/2);			% Find central point
		mm = [1 mMo] - mc;		nn = [1 nMo] - nc;			% Shift vector about central point
		for (i=mm(1):mm(end))
			for (j=nn(1):nn(end))
 				quad_{i+mc,j+nc} = getNext(quad, quadkey, i, j);
				if (quadkey{1} == '0')						% VE
					url{i+mc,j+nc} = [prefix quad_{i+mc,j+nc} '?g=244'];
				else
					decAdr = getQuadLims([quadkey{2} quad_{i+mc,j+nc}], quadkey, 1);		% Get decimal adress
					if (whatKind(1) == 'a')					% ~VE and 'aerial'
						if (isMapForFree)					% We already know the name (e.g. maps-for-free) - no QuadTree in it
							url{i+mc,j+nc} = [pref_bak sprintf('%d/row%d/%d_%d-%d.jpg',zoomL-1,decAdr(2),zoomL-1,decAdr(1),decAdr(2))];
						else
							url{i+mc,j+nc} = [prefix quad_{i+mc,j+nc}];
							if ( numel(prefix) > 9 && strcmp(prefix(8:10), 'mw1') )			% Good for Mars but fails for Moon (different algo)
								url{i+mc,j+nc} = [url{i+mc,j+nc} '.jpg'];
							elseif ( numel(prefix) > 9)
								url{i+mc,j+nc} = [pref_bak sprintf('%d&y=%d&z=%d',decAdr, zoomL-1)];
							end
						end
					else									% ~VE and road|hybrid
						url{i+mc,j+nc} = [pref_bak sprintf('%d&y=%d&z=%d',decAdr, zoomL-1)];
					end
				end
			end
		end
		if (numel(lon) == 1)		% Otherwise (rectangle limits on input) we already know them
			lon_mm 	= [(lon_mm(1) + nn(1)*diff(lon_mm)) (lon_mm(2) + nn(2)*diff(lon_mm))];
			latiso_mm = [(latiso_mm(1) - mm(2)*diff(latiso_mm)) (latiso_mm(2) - mm(1)*diff(latiso_mm))];
			lat_mm = isometric2geod(latiso_mm, flatness);
		end
	end
	
	% ----------------------- Compute mercator meters (even if don't need they are cheap) ------------ 
	D2R = pi / 180;
	if (flatness ~= 0)		% Copyed from GMT
		meridionalRadius = meridionalRad(eq_radius,flatness);
		eq_radius = eq_radius * (eq_radius / meridionalRadius);
	end
	x = lon_mm * D2R * eq_radius;
	if (flatness ~= 0)
		y = log (tan (pi/4 + 0.5 * D2R * latiso_mm)) * eq_radius;
	else
		y = log (tan (pi/4 + 0.5 * D2R * lat_mm)) * eq_radius;
	end
	% -------------------------------------------------------------------------------------------------

	if ( strcmp(opt, 'callmir') || strcmp(opt, 'tile2img') )

		nTilesX = size(mosaic,2);			nTilesY = size(mosaic,1);
		if (reportMerc)			% Output coords in Mercator
			dx = diff(x)/(256*nTilesX);			dy = diff(y)/(256*nTilesY);		% Here we are still in the Pixel reg world
			x = x + [1 -1] * dx / 2;			y = y + [1 -1] * dy / 2;		% But now we are on Grid reg
			tmp.head = [x y 0 255 0 dx dy];
			tmp.X = x;			tmp.Y = y;
		else					% Output coords in Geogs (image is unchanged, that is it remain mercatorized)
			dx = diff(lon_mm)/(256*nTilesX);	dy = diff(lat_mm)/(256*nTilesY);		% Here we are still in the Pixel reg world
			lon_mm = lon_mm + [1 -1] * dx / 2;	lat_mm = lat_mm + [1 -1] * dy / 2;		% But now we are on Grid reg
			tmp.head = [lon_mm lat_mm 0 255 0 dx dy];
			tmp.X = lon_mm;		tmp.Y = lat_mm;
		end
		tmp.geog = 0;		tmp.name = 'EuGooglo';
		ext = 'jpg';		% Default image extension
		if ( (whatKind == 'r') || (whatKind == 'h') ),	ext = 'png';	end		% Maps use png files

		if (~iscell(url))
			[cache, cache_supp, ext, isWW] = completeCacheName(cache, ext, zoomL, decimal_adress, whatKind);
			[img, cmap] = getImgTile(quadkey, quad_, url, cache, cache_supp, ext, decimal_adress, isWW, verbose);
		else
			img = repmat( repmat(uint8(200), [256 256 3]), mMo, nMo );		% Pre-allocate the final image size
			for (i=1:mMo)
				[cache, cache_supp, ext, isWW] = completeCacheName(cache, ext, zoomL, decimal_adress, whatKind);	% WW changes cache dir per rows
				for (j=1:nMo)
					if ( isWW && quadkey{1} == '0' )		% VE and WW cache - We need to recompute the dec adress for each quad
						decimal_adress = getQuadLims(['a' quad_{mMo-i+1,j}], quadkey, 1);	% Get decimal adress
					end
					[img((i-1)*256+1:i*256, (j-1)*256+1:j*256, :), cmap] = ...
						getImgTile(quadkey, quad_{mMo-i+1,j}, url{mMo-i+1,j}, cache, cache_supp, ext, decimal_adress, isWW, verbose);
				end
			end
		end
		if ( strcmp(opt, 'callmir') )
			if (~isempty(cmap)),	tmp.cmap = cmap;	end
			if (reportMerc),		tmp.srsWKT = ogrproj(['+proj=merc +R=' num2str(geoid(1))]);	end
			mirone(img, tmp)
		else
			url = img;
			lon_mm = tmp;
		end
	end

% --------------------------------------------------------------------------------------	
function [cache, cache_supp, ext, isWW] = completeCacheName(cache, ext, zoomL, decimal_adress, whatKind)
% See if we have a "World Wind" with "Virtual Earth" cache. If yes complete the dir name for
% the tile in question, which is expressed in DECIMAL_ADRESS.
% For the type 2) cache refered at the main help, if cache dir doesn't end with 'cache', append it as well
	isWW = false;				% Will be true if CACHE is a WW cache
	cache_supp = [];
	if ( ~isempty(cache) )
		[pato, nome] = fileparts(cache);
		fsep = filesep;
		if ( strcmp(nome, 'Virtual Earth') )
			%cache = [cache fsep sprintf('%d%c%c%c%.4d',zoomL-1,fsep,whatKind(1), fsep, decimal_adress(2) )];
			if (whatKind(1) == 'a'),	ext = 'jpeg';		% Ghrr, could they not have used .jpg ?
			else						ext = 'png';
			end
			isWW = true;
			cache_supp = [fsep sprintf('%d%c%c%c%.4d',zoomL-1,fsep,whatKind(1), fsep, decimal_adress(2) )];
		else
			if ( ~strcmpi(nome, 'cache') )				% If the cache dir doesn't end by 'cache', add it
				cache = [cache fsep 'cache'];
			end
			plusZLnumber = sprintf('%s%.2d',fsep, zoomL);		% Append the zoomlevel to the cache dir tree
			if (whatKind(1) == 'a')
				cache_supp = [fsep 'kh' plusZLnumber];			ext = 'jpg';
			elseif (whatKind(1) == 'r')
				cache_supp = [fsep 'mt' plusZLnumber];			ext = 'png';
			else
				cache_supp = [fsep 'tt' plusZLnumber];			ext = 'png';
			end
		end
	end

% --------------------------------------------------------------------------------------	
function mappcache(flatness, quadkey, quad, cache, ext)
% Calculate & save a mappimg of the existing tiles in the 'cache' dir
% The mapping file has a fix name of tilesMapping.mat and will contain two variables
% region -> all points bounding box
% tiles_midpt -> a Mx2 matrix with tiles mid point coordinates

	if (nargin < 3)
		error('url2image:mappcache', 'number of argins must be >= 3')
	end
	if (~isempty(quad) && ~(ischar(quad) || iscell(quad)) )
		error('url2image:mappcache','"QUAD" must be a string or a cell array of strings')
	end
	if (~exist(cache, 'dir'))
		error('url2image:mappcache','The ''cache'' directory name does not exist. Bye Bye.')
	end
	if (nargin == 4),		ext = 'jpg';	end

	if (isempty(quad))		% Scan 'cache' dir for existing tiles and use them to build the mapping file
		files = dir([cache filesep '*.' ext]);
		if (isempty(files))
			warning('url2image:mappcache', 'There is nothing useful in the provided directory')
			return
		end
		quad = {files.name};
	end
	
	[lims, tiles_bb, zoomL] = quadcoord(quad, quadkey, flatness);
	
	% Calculate mid point of each tile
	tiles_midpt = [(tiles_bb(:,1) + tiles_bb(:,2)) / 2 (tiles_bb(:,3) + tiles_bb(:,4)) / 2];
	
	region = [min(tiles_midpt(:,1)) max(tiles_midpt(:,1)) min(tiles_midpt(:,2)) max(tiles_midpt(:,2))];
	% Make bb_limits pix reg in longitude and apply the same margin to lat (too complicated to do correctly)
	region = region + [-0.5 0.5 -0.5 0.5] * diff(tiles_bb(1:2));
	if ( region(3) < 85.0841 ),		region(3) = -85.0841;	end		% Isometric 180
	if ( region(4) > 85.0841 ),		region(4) =  85.0841;	end

	save([cache filesep 'tilesMapping'], 'region', 'zoomL', 'tiles_midpt')

% --------------------------------------------------------------------------------------
function [lims, tiles_bb, zoomL] = quadcoord(flatness, quadkey, quad)
	% Compute the coordinates of a QUAD quadtree file or string
	% LIMS contains the BoundingBox coords in geogs
	% TILES_BB contains the rectangle BB coords of each tile [lonMin lonMax latMin latMax]
	% This matrix is written by columns. That is like coords(quad{:})

	if ( ~(isa(quad,'char') || isa(quad,'cell')) )
		error('url2image:quadcoord','First arg must be a string or a cell array of strings')
	end
	
	n_argout = nargout;
	if (isa(quad,'char'))		% A single quadtree
		[lims, zoomL] = getQuadLims(quad, quadkey);
		tiles_bb = lims([1 2 4 3]);			% In case idiot choice of 2 argouts
	else						% Several in a cell array
		if (n_argout == 2),		tiles_bb = zeros(numel(quad),4);	end		% Pre-allocation
		lims = [361 -361 -361 361];
		for (k = 1:numel(quad))
			[lim, zoomL] = getQuadLims(quad{k}, quadkey);		% Remember that lim(3) = y_max
			lims =  [min(lims(1), lim(1)) max(lims(2), lim(2)) max(lims(3), lim(3)) min(lims(4), lim(4))];
			if (n_argout >= 2)
				tiles_bb(k,:) = lim;
			end
		end
	end
	lims(3:4) = isometric2geod(lims(4:-1:3), flatness);
	if (n_argout >= 2)
		tiles_bb(:,3:4) = isometric2geod(tiles_bb(:,4:-1:3), flatness);
	end

% ----------------------------------------------------------------------------	
function new_quad = getNext(quad, quadkey, v, h)
% To navigate up/down/left/right, we can translate the string of q-t this way: q=0 r=1 t=2 s=3
% We can treat each letter as a 2-bit number, where bit 0 means left/right and bit 1 means up/down.
% Therefore, a string gets translated like: 
%    tqstqrrqstsrtsqsqr =
%  = 203201103231230301
%  H 001001101011010101 (horizontal)
%  V 101100001110110100 (vertical)
%
% Now, we have two 18-bit binary numbers that encode position. To navigate east, for example,
% we simply increment the horizontal number, then reencode: 
%  H 001001101011010110
%  V 101100001110110100
%  = 203201103231230310
%
%	From: web.media.mit.edu/~nvawter/projects/googlemaps/index.html

	if (isequal([h v], [0 0]))		% Trivial case. Nothing to do
		new_quad = quad;		return
	end

	N = numel(quad);			quad_num = zeros(1, N);
	ind = (quad == quadkey{1,2});	quad_num(ind) = 1;		% 1|r 		Convert to numeric
	ind = (quad == quadkey{2,2});	quad_num(ind) = 3;		% 3|s
	ind = (quad == quadkey{2,1});	quad_num(ind) = 2;		% 2|t
	
	H = repmat(' ', 1, N);		new_quad = H;	% pre-allocations
	tmp = dec2bin(quad_num,2);
	V = tmp(:,1);				H = tmp(:,2);
	
	if (h),		H = dec2bin(bin2dec(H')+h, N);		end		% East-West
	if (v),		V = dec2bin(bin2dec(V')+v, N);		end		% Up-Down
	% Test if we are getting out of +180. If yes, wrap the exccess it to -180
	if (numel(H) > N),		H = dec2bin(h-1, N);	end
	% Test if we are getting out of -180. If yes, wrap the exccess it to +180
	if ( H(1) == '/' ),		H = dec2bin(bin2dec(tmp(:,2)')+(2^N - 1), N);	end
	
	new_tile_bin = [V(:) H(:)];
	quad_num = bin2dec(new_tile_bin);

	% Now reencode the new numeric quadtree
	ind = (quad_num == 0);	new_quad(ind) = quadkey{1,1};		ind = (quad_num == 1);	new_quad(ind) = quadkey{1,2};
	ind = (quad_num == 2);	new_quad(ind) = quadkey{2,1};		ind = (quad_num == 3);	new_quad(ind) = quadkey{2,2};

% ----------------------------------------------------------------------------	
function [lims, zoomL] = getQuadLims(quad, quadkey, opt)
% Compute the limits of a quadtree file or string.
% OPT == [] => lims = [lon1 lon2 lat1 lat2];	ELSE,	lims = [pixelX pixelY];

	if (nargin == 2),	opt = [];	end
	[pato, quad] = fileparts(quad);		% If quad is a filename, retain only what matters
	if (strncmpi(pato,'http',4))				% QUAD is embeded in a URL
		if (quadkey{1} ~= '0')					% filename starts after a '='
			ind = strfind(quad,'=');
			quad = quad(ind(end)+1:end);
		end
	end

	zoomL = numel(quad);		quad_x = zeros(1, zoomL-1);		quad_y = quad_x;

	quad = (quad(2:end));
	% MAPPING: X -> 0=0 1=1 3=1 2=0;		Y -> 0=0 1=0 3=1 2=1
	% x component
	ind = (quad == quadkey{1,2});	quad_x(ind) = 1;		% Convert to numeric
	ind = (quad == quadkey{2,1});	quad_x(ind) = 0;
	ind = (quad == quadkey{2,2});	quad_x(ind) = 1;
	% And y component
	ind = (quad == quadkey{1,2});	quad_y(ind) = 0;
	ind = (quad == quadkey{2,1});	quad_y(ind) = 1;
	ind = (quad == quadkey{2,2});	quad_y(ind) = 1;
	
 	pixelX = 0;		pixelY = 0;
	for k=1:numel(quad)
		pixelX = pixelX*2 + quad_x(k);
		pixelY = pixelY*2 + quad_y(k);
	end

	if (isempty(opt))
		[lon1, lat1] = getLonLat(pixelX, pixelY, zoomL);
		[lon2, lat2] = getLonLat(pixelX+1, pixelY+1, zoomL);
		lims = [lon1 lon2 lat1 lat2];
	else
		lims = [pixelX pixelY];
	end

% ----------------------------------------------------------------------------	
function [x, y, xmm, ymm] = getPixel(lon, lat, zoomL)
% Compute the pixel coordinate for a given longitude and latitude.
% In fact x,y are the fractional number of the 256 bins counting from origin
	pixPerDeg = 2^(zoomL - 1) / 360;
	x = (lon + 180) * pixPerDeg;
 	y = (180 - lat) * pixPerDeg;
	y = y(end:-1:1);
	if (numel(lon) == 1)
		xmm = [fix(x) fix(x)+1] * 256;		% [x_min x_max] pixel coords of the rectangle
		ymm = [fix(y) fix(y)+1] * 256;
	elseif (numel(lon) == 2)				% lon, lat contain a rectangle limits
		xmm = [fix(x(1)) fix(x(1))+1; fix(x(2)) fix(x(2))+1] * 256;
		ymm = [fix(y) fix(y)+1] * 256;
		ymm = reshape(ymm, 2, 2);
	else
		error('url2image:getPixel','Lon must be a scalar or a 2 elements vector')
	end
	
% ----------------------------------------------------------------------------	
function [lon, lat] = getLonLat(pixelX, pixelY, zoomL)
% Compute the longitude and isometric latitude for a given pixel coordinate.
	pixPerDeg = 2^(zoomL - 1) / 360;
	lon = pixelX / pixPerDeg - 180;
	lat = 180 - pixelY / pixPerDeg;

% ----------------------------------------------------------------------------	
function [img, cmap] = getImgTile(quadkey, quad, url, cache, cache_supp, ext, decAdr, isWW, verbose)
% Get the image either from a local cache or by url 
% cache_supp contains (for example) the /../13 subdirs
	try
		if (quadkey{1} == '0'),		quad = [url(8) quad];			% (Virtual Earth -> url(8) = 'a'|'r'|'h')
		else						quad = [quadkey{2} quad];
		end

		cmap = [];		att = [];	fname = [];		img = [];
		if (~isempty(cache))		% We have a cache dir to look for

			if ( isWW )		% We guessed it before as beeing a WW cache
				fname = [cache cache_supp filesep sprintf('%.4d_%.4d', decAdr(2), decAdr(1)) '.' ext];
			else
				fname = [cache cache_supp filesep quad '.' ext];
			end
			if (exist(fname, 'file'))
				if (verbose > 1),	disp(['Retrieving file from cache: ' fname]),	end
				[img, cmap] = imread(fname);
			else										% File is not in cache. Need to download
				[img, att] = netFetchTile(url, [cache cache_supp], quad, ext, verbose, fname);
			end
		else						% No cache, download it
			[img, att] = netFetchTile(url, [], quad, ext, verbose);
		end
		img = flipdim(img, 1);
		
		if ( isempty(img) )
			img = repmat(uint8(200), [256 256 3]);
		elseif ( ndims(img) == 2 && ext(1) == 'p')
			if ( ~isempty(cmap) )
				img = ind2rgb8(img, cmap);
			elseif ( isempty(cmap) && ~isempty(att) )
				cmap = att.Band.ColorMap.CMap(:,1:3);
				img = ind2rgb8(img, cmap);
			else
				disp('url2image:getImgTile  Unknown error in retrieving image')
				img = repmat(uint8(200), [256 256 3]);
			end
		elseif ( ndims(img) == 2 && ext(1) == 'j')		% VE returns a 256x256 when the tile has no image
			img = repmat(uint8(200), [256 256 3]);
		end
	catch
		disp(lasterr)
		if ( isempty(img) && ~isempty(fname) )			% An empty screwed tile exists in the cache. Delete it
			builtin('delete',fname);
		end
		img = repmat(uint8(200), [256 256 3]);
	end

% --------------------------------------------------------------------------------------	
function [img, att] = netFetchTile(url, cache, quad, ext, verbose, fname)
% Fetch a file from the web either using gdal or wget (when gdal is not able to)
	if (verbose),	disp(['Downloading file ' url]),	end
	if (true || strcmp(url(8:9),'kh') || strcmp(url(8:9),'mt'))		% --> FORCE USE OF WGET 
		if ( ~isempty(cache) )
			if ( exist(cache,'dir') ~= 7 ),		make_dir(cache),	end		% Need to create a new dir
			dest_fiche = [cache filesep quad '.' ext];		% Save the file directly in susitu (cache)
		else
			dest_fiche = 'lixogrr';
		end

		if (ispc),		dos(['wget "' url '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
		else			unix(['wget ''' url ''' -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
		end
		finfo = dir(dest_fiche);
		if (finfo.bytes < 100)					% Delete the file anyway because it exists but is empty
			disp(['Failed to download file: ' url])
			builtin('delete',dest_fiche);
			img = [];							% At least this way it won't generate an error
			return
		end
		[img, att] = gdalread(dest_fiche);
		if (size(img,3) == 4),		img(:,:,4) = [];	end		% We don't deal yet with alpha-maps in images
	else
		[img, att] = gdalread(url);				% Don't read flipped ('-U') because of saving
		if ( ~isempty(cache) )
			saveInCache(cache, fname, att, img, ext)	% Save tile in the cache directory
		end
	end

% --------------------------------------------------------------------------------------	
function saveInCache(cache, fname, att, img, ext)
% This function gets called when a file was downloaded. If a cache dir info is available, save it there. 
% Still, if cache (main)dir exists but not the required cache (sub)dir, than create it (or raise error exception)
% CACHE is here the full path to the cache directory name. As just said, if it doesn't exist creat it
% FNAME is the full file name, which is = CACHE/image_name (image_name has some degrees of freedom)
% ATT is the gdalread attribute structure
% IMG is ... the image (2D or 3D)
% EXT is either 'jpg', 'jpeg' (3D) or 'png' (2D)

	if ( isempty(cache) ),		return,		end			% No CACHE, no money
	
	if ( exist(cache,'dir') ~= 7 )						% Ghrr, make a new dir
		make_dir(cache)
	end

	if ( ndims(img) == 3 && (strcmp(ext,'jpg') || strcmp(ext,'jpeg')) )
		imwrite(img, fname, 'Quality', 100);
	elseif ( ndims(img) == 3 && strcmp(ext,'png') )
		imwrite(img, fname);
	elseif ( ndims(img) == 2 && strcmp(ext,'png') )
		cmap = att.Band.ColorMap.CMap(:,1:3);
		imwrite(img, cmap, fname);
	elseif ( ndims(img) == 2 && strcmp(ext,'jpg') )		% Do nothing. VE returns a 256x256 when the tile has no image
	else
		error('url2image:saveInCache', 'Shit. Unknown error while trying to save image in cache')
	end

% --------------------------------------------------------------------------------------	
function make_dir(cache)
% Create the directory CACHE
% At least on R13 MKDIR is a brain-dead fun. WHY CAN'T the F. WE CREAT AN ABSOLUTE PATH 'DIR' ???????
	try
		if (strncmp(computer,'PC',2))
			if (strfind(cache, ' ')),	cache = ['"' cache '"'];	end
			disp(['Creating DIR: ' cache])
			Status = dos(['mkdir ' cache]);
		else								
			Status = unix(['mkdir -p ' cache]);
		end
		if (Status),	error('url2image:make_dir', 'Error while creating new cache (sub)directory'),	end
	catch
		error('url2image:make_dir', 'Error while creating new cache directory')
	end

% ----------------------------------------------------------------------------	
function [a, b]	= rema(x, y)
	a = fix(x / y);
	b = x - a*y;

% -------------------------------------------------------------------------------
function lat = isometric2geod(latin, flat)
% Computes the geodetic latitude given the isometric latitude

	if (nargin == 1)
		flat = 1/298.257223563;			% WGS84
	end

	latcnf = 2 * atan(exp(latin * pi / 180)) - pi/2;       %  Compute the conformal lat
	lat = conf2geod(latcnf,flat);		%  Transform to geodetic
	lat = lat * 180 / pi; 

% -------------------------------------------------------------------------------
function lat = conf2geod(latin,f)
% Converts from conformal latitude to geodetic latitude

	e2 = f * (2.0 - f);		e4 = e2 * e2;
	e6 = e4 * e2;			e8 = e4 * e4;

	%  Compute the series expansion terms
	c0 = e2 /2 + 5*e4 /24 + e6 /12 + 13*e8 /360;
	c1 = 7*e4 /48 + 29*e6 /240 + 811*e8 /11520;
	c2 = 7*e6 / 120 + 81*e8 /1120;
	c3 = 4279*e8 / 161280;

	rl2 = 2.0 * latin;
	sin2phi = sin(rl2);		cos2phi = cos(rl2);
	lat = latin + sin2phi .* (c0 + cos2phi .* (c1 + cos2phi .* (c2 + cos2phi * c3)));

% -------------------------------------------------------------------------------
function lat = geod2isometric(latin, flat)
% Computes the isometric latitude given the geodetic latitude

	if (nargin == 1)
		flat = 1/298.257223563;			% WGS84
	end
	e = sqrt(2*flat - flat.^2);

	latin = latin * pi / 180;

	latcnf = geod2cnf(latin,e);   %  Compute the conformal lat
	lat = log ( tan(pi/4 + latcnf/2) );      %  Transform to isometric
	lat = lat * 180 / pi; 

% -------------------------------------------------------------------------------
function latout = geod2cnf(latin,e)
% Converts from geodetic latitude to conformal latitude

	f1  = 1 - e*sin(latin);     f2 = 1 + e*sin(latin);
	f3  = 1 - sin(latin);       f4 = 1 + sin(latin);
	latout = 2 * atan(sqrt((f4./f3) .* ((f1./f2).^e)) ) - pi/2;

% -------------------------------------------------------------------------------
function meridionalRadius = meridionalRad(a,f)
% Compute Meridional radius as deffined in GMT_lar_swap_init() of gmt_map.c
%	A Equatorial radius, F flatness

	e2 = f * (2.0 - f);		e4 = e2 * e2;
	e6 = e4 * e2;			e8 = e4 * e4;
	
	xx0 = 1 / 4;
	xx1 = xx0 * 3 / 16;
	xx2 = xx1 * 3 * 5 / 36;
	xx3 = xx2 * 5 * 7 / 64;
	x = xx0 * e2 + ( xx1 * e4 + ( xx2 * e6 + xx3 * e8));
	meridionalRadius = a * (1 - x);	


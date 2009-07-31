function shapenc(fname, data, varargin)
% Create netCDF shape type files
%
% FNAME -> Name of the netCDF output file
% DATA  -> Mx2, or Mx3, or filename of a Mx2, Mx3 table, or a shapefile, or a cell of Mx2, Mx3 data
% shapenc(..., 'outer', polyg)	where POLYG is a polygon embracing the DATA data
% shapenc(..., 'inner', polyg)	where POLYG is a polygon representing a hole in the DATA
%								Several polygons/holes are allowed but if > 1 they must be in a cell array
% shapenc(..., 'srs','proj4string')  DATA is in projected coords described by PROJ4STRING 
% shapenc(..., 'geog',1|0)  If keyword GEOG is used  PROJ4STRING -> +proj=latlong  ---- NOT CONFIRMED
% shapenc(..., 'desc','description')  Create a global attribute with contents DESCRIPTION
%
% possibilidades
% data				<== Point swarm 
% data, data		<== Point swarm + external polygon
% data, data, data	<== Point swarm + 1 external polygon + 1 internal polygon
% data, data, cell	<== Um poly externo e varios internos
% cell				<== Cada cellula tem um molho
% cell, cell		<== Cada cellula tem um molho e cada molho tem 0 ou 1 poly externo

	% ------------------ Parsing and defaults settings------------------ --------------
	if (nargin < 2),	error('SHAPENC:error', 'Please, minimum 2 arguments -- FNAME & DATA'),	end

	is_Pts_cell = false;	n_swarms = 1;	is_file_new = true;		ultimo = 0;		is_geog = true;
	outer_polygs = [];		inner_polygs = [];	nOuterPolys = 0;	nInnerPolys = 0;
	desc = [];		spatial_ref = '+proj=longlat';
	is_point_2D = false;	is_point_3D = false;	is_polygon_2D = false;	is_polygon_3D = false;
	is_polyline_2D = false;	is_polyline_3D = false;
	nMaxChunks = 50;		% Apply only to shapefiles and multisegment ascii files
	multiSegPos = [];		% Positions of NaN multiseg separators in a packed multisegment data array (mostly for shapefiles)
	fver = [];				% File version

	for (k = 1:2:numel(varargin))
		if (~ischar(varargin{k})),		error('SHAPENC:error','property name must be a character string'),	end
		if ( strncmpi(varargin{k}, 'out', 3) )
			outer_polygs = varargin{k+1};
		elseif ( strncmpi(varargin{k}, 'inner', 2) )
			inner_polygs = varargin{k+1};
		elseif ( strncmpi(varargin{k}, 'geog', 4) )
			if (~varargin{k+1}),	is_geog = false;	end
			spatial_ref = '+proj=xy';
		elseif ( strcmpi(varargin{k}, 'srs') )
			spatial_ref = varargin{k+1};
		elseif ( strncmpi(varargin{k}, 'desc', 4) )
			desc = varargin{k+1};
		elseif ( strcmpi(varargin{k}, 'point2d') )
			is_point_2D = true;
		elseif ( strcmpi(varargin{k}, 'polygon2D') )
			is_polygon_2D = true;
		elseif ( strcmpi(varargin{k}, 'polygon3D') )
			is_polygon_3D = true;
		elseif ( strcmpi(varargin{k}, 'polyline2D') )
			is_polyline_2D = true;
		elseif ( strcmpi(varargin{k}, 'polyline3D') )
			is_polyline_3D = true;
		elseif ( strcmpi(varargin{k}, 'maxpoly') )
			nMaxChunks = varargin{k+1};
		elseif ( strncmpi(varargin{k}, 'version', 3) )
			fver = varargin{k+1};
		end
	end
	if ( ~is_geog && strcmp(spatial_ref(7:end), 'longlat') ),	spatial_ref = '+proj=xy';	end

	if (ischar(data))		% Try to load data from file (including shapefiles)
		try
			[data, multiSegPos, is_polygon_2D, is_polygon_3D, is_polyline_2D, is_polyline_3D, is_point_2D] = ...
				readFile(data, nMaxChunks, is_polygon_2D, is_polygon_3D, is_polyline_2D, is_polyline_3D, is_point_2D);
		catch
			error(['SHAPENC:reading DATA ' lasterr])
		end
	end
	if (ischar(outer_polygs))		% Try to load data from file
		try			outer_polygs = readFile(outer_polygs, nMaxChunks);
		catch		error(['SHAPENC:reading OUTER_POLYGS ' lasterr])
		end
	end
	if (ischar(inner_polygs))		% Try to load data from file
		try			inner_polygs = readFile(inner_polygs, nMaxChunks);
		catch		error(['SHAPENC:reading INNER_POLYGS ' lasterr])
		end
	elseif (isa(inner_polygs,'cell'))	% Hmm, must see if it is a cell of fnames to load
		if (ischar(inner_polygs{1}))	% Yep, so they seam
			try
				for (k = 1:numel(inner_polygs))
					inner_polygs{k} = readFile(inner_polygs{k}, nMaxChunks);
				end
			catch
				error( sprintf('SHAPENC:screw_reading %d element of cellular INNER_POLYGS \n%s', k, lasterr) )
			end
		end
	end

	if (isa(data,'cell')),		is_Pts_cell = true;		n_swarms = numel(data);		end

	if ( is_Pts_cell && ~isempty(outer_polygs) && ~isa(outer_polygs,'cell') )
		error('SHAPENC:error','When DATA is a cell array so must be the POLYGON')
	end

	if (is_geog)
		long_name = {'Longitude' 'Latitude'};	units = {'degrees_east' 'degrees_north'};
	else
		long_name = {'unknown' 'unknown'};		units = {'xunits' 'yunits'};
	end
	% -------------------------- END PARSING & DEFAULT SETTINGS -------------------------

	% ------------------ Guessing based on number of columns ----------------------------
	% Point_2D and Point_2D are guessed (if not provided earlier) here
	% All other cases must be explicitly declared
	if (~(is_polygon_2D || is_polygon_3D || is_polyline_2D || is_polyline_3D) )
		if (is_Pts_cell)
			if ( size(data{1},2) == 2 ),		is_point_3D = false;	is_point_2D = true;
			elseif ( size(data{1},2) == 3 ),	is_point_3D = true;		is_point_2D = false;
			else 	error('SHAPENC:error','DATA must be a Mx2 OR a Mx3 array')
			end
		else
			if ( size(data,2) == 2 ),			is_point_3D = false;	is_point_2D = true;
			elseif ( size(data,2) == 3 ),		is_point_3D = true;		is_point_2D = false;
			else 	error('SHAPENC:error','DATA must be a Mx2 OR a Mx3 array')
			end
		end
	end
	% -----------------------------------------------------------------------------------

	% ------------------ Count number of points in each (if > 1) ensemble --------------
	if (is_Pts_cell)
		n_pts = zeros(numel(data));
		for (k = 1:n_swarms),		n_pts(k) = size(data{k},1);			end
		if (isa(outer_polygs,'cell'))
			nOuterPolys = zeros(numel(outer_polygs));
			for (k = 1:numel(outer_polygs)),	nOuterPolys(k) = size(outer_polygs{k},1);		end
		end
		if (isa(inner_polygs,'cell'))
			nInnerPolys = zeros(numel(inner_polygs));
			for (k = 1:numel(inner_polygs)),	nInnerPolys(k) = size(inner_polygs{k},1);		end
		end
	else
		n_pts = size(data,1);
		if (~isempty(outer_polygs)),		nOuterPolys = size(outer_polygs,1);		end
		if (~isempty(inner_polygs))
			if (~isa(inner_polygs, 'cell'))
				nInnerPolys = size(inner_polygs,1);
			else
				nInnerPolys = zeros(1, numel(inner_polygs));
				for (k = 1:numel(inner_polygs)),	nInnerPolys(k) = size(inner_polygs{k},1);	end
			end
		end
	end
	% -----------------------------------------------------------------------------------

	% -----------------------------------------------------------------------------------
	if (exist(fname,'file') ~= 2)			% If file does not exist, create it.
		nc_funs('create_empty', fname)
	else
		is_file_new = false;
		s = nc_funs('info',fname);
		if (~isempty(s.Attribute)),		attribNames = {s.Attribute.Name};
		else							attribNames = [];
		end
		ind = strcmp(attribNames,'SHAPENC_type');
		if ( any(ind) )
			dimNames = {s.Dimension.Name};
			ind = strncmp(dimNames,'dimpts_',7);		% Find the dimpts_??? variables
			for (k = 1:numel(ind))
				if (~ind(k)),	continue,	end			% Jump other dimension names
				ultimo = max(ultimo, str2double( dimNames{k}(8:end) ));		% Get the last used name
			end
		end
	end
	% -----------------------------------------------------------------------------------

	% ---------------------------- Write the dimensions --------------------------------
	for (k = 1:n_swarms)		% Start with the Points dims
		nc_funs('add_dimension', fname, sprintf('dimpts_%d',k+ultimo), n_pts(k) )
	end
	if (~isempty(outer_polygs))		% Write dim of external polygons
		for (k = 1:n_swarms)
			if (~nOuterPolys(k)),	continue,	end		% This swarm has no polygon. Jump it
			nc_funs('add_dimension', fname, sprintf('dimpolyOUT_%d',k+ultimo), nOuterPolys(k) )			
		end
	end
	if (~isempty(inner_polygs))		% Write dim of internal polygons
		for (k = 1:numel(nInnerPolys))
			if (~nInnerPolys(k)),	continue,	end		% This inner polygon is empty. Jump it
			nc_funs('add_dimension', fname, sprintf('dimpolyIN_%d_%d', ultimo+1, k+ultimo), nInnerPolys(k) )			
		end
	end
	% -----------------------------------------------------------------------------------

	% ---------------------------- Write Globals attributes -----------------------------
	global_BB = [Inf -Inf Inf -Inf];
	if (is_file_new)
		nc_funs('attput', fname, -1, 'title', 'SHAPENC: A netCDF extended storage version of some types of shapefiles');
		if (~isempty(desc)),	nc_funs('attput', fname, -1, 'Description', desc);		end
		nc_funs('attput', fname, -1, 'version', '1.0');
		if (is_point_3D)
			nc_funs('attput', fname, -1, 'SHAPENC_type', 'PointZ');
		elseif (is_point_2D)
			nc_funs('attput', fname, -1, 'SHAPENC_type', 'Point');
		elseif (is_polygon_3D)
			nc_funs('attput', fname, -1, 'SHAPENC_type', 'PolygonZ');
		elseif (is_polygon_2D)
			nc_funs('attput', fname, -1, 'SHAPENC_type', 'Polygon');
		elseif (is_polyline_3D)
			nc_funs('attput', fname, -1, 'SHAPENC_type', 'PolyLineZ');
		elseif (is_polyline_2D)
			nc_funs('attput', fname, -1, 'SHAPENC_type', 'PolyLine');
		end
		nc_funs('attput', fname, -1, 'spatial_ref', spatial_ref);
		nc_funs('attput', fname, -1, 'BoundingBox', global_BB);
	else
		% Get old global BB
		ind = strcmp(attribNames,'BoundingBox');
		if ( any(ind) ),				global_BB = s.Attribute(ind).Value;		end
		if (numel(global_BB) ~= 4),		global_BB = [Inf -Inf Inf -Inf];		end
	end
	nc_funs('attput', fname, -1, 'Number_of_main_ensembles', ultimo+n_swarms );
	if (~isempty(fver))
		nc_funs('attput', fname, -1, 'File Version', fver );
	end
	% -----------------------------------------------------------------------------------

	% ------------------------------ Write the variables --------------------------------
	if (is_point_2D || is_point_3D)
		write_pt_type(fname, data, outer_polygs, nOuterPolys, inner_polygs, nInnerPolys, is_point_3D, is_geog, is_Pts_cell, ...
						long_name, units, ultimo, global_BB, n_swarms)
	else
		write_poly_type(fname, data, is_polygon_2D, is_polygon_3D, is_polyline_2D, is_polyline_3D, is_geog, is_Pts_cell, ...
						long_name, units, ultimo, global_BB, n_swarms, multiSegPos)
	end


% ----------------------------------------------------------------------------------------------
function write_pt_type(fname, data, outer_polygs, nOuterPolys, inner_polygs, nInnerPolys, is_point_3D, is_geog, is_Pts_cell, ...
						long_name, units, ultimo, global_BB, n_swarms)
% Write variables for the 2D & 3D point swarm cases

	if (is_point_3D),	prefix = 'PointZ_';
	else				prefix = 'Point_';
	end

	if (is_geog),		Xname = 'lon';		Yname = 'lat';
	else				Xname = 'X';		Yname = 'Y';
	end

	for (k = 1:n_swarms)
		kk = k + ultimo;
		if (is_Pts_cell),	x = data{k}(:,1);		y = data{k}(:,2);
		else				x = data(:,1);			y = data(:,2);
		end
		if ( ~isa(x,'double') && any(isnan(x)) )	% BLOODY ML BUGS (R13)
			BB = [min(double(x)) max(double(x)) min(double(y)) max(double(y))];
		else
			BB = double([min(x) max(x) min(y) max(y)]);
		end
		if (is_point_3D)		% Split case between 3d and 2D Point ensembles
			if (is_Pts_cell),	z = data{k}(:,3);
			else				z = data(:,3);
			end
			BB = [BB double([min(z) max(z)])];
		end
		if (k == 1),	Nctype = getDataType(x);		end		% Get data type

		write_var(fname, sprintf('%s%s%d', Xname, prefix, kk), Nctype, sprintf('dimpts_%d',kk), long_name{1}, units{1})
		nc_funs('varput', fname, sprintf('%s%s%d', Xname, prefix, kk), x);
		write_var(fname, sprintf('%s%s%d', Yname, prefix, kk), Nctype, sprintf('dimpts_%d',kk), long_name{2}, units{2})
		nc_funs('varput', fname, sprintf('%s%s%d', Yname, prefix ,kk), y);
		if (is_point_3D)
			write_var(fname, sprintf('z_%d',kk), Nctype, sprintf('dimpts_%d',kk), 'z', 'meters')
			nc_funs('varput', fname, sprintf('z_%d',kk), z);
		end

		% Create a container variable to hold the BB info
		write_var(fname, sprintf('%s%d', prefix, kk), 2)
		nc_funs('attput', fname, sprintf('%s%d', prefix, kk), 'BoundingBox', BB);
		
		% Update the global BB
		global_BB = [min(BB(1), global_BB(1)) max(BB(2), global_BB(2)) min(BB(3), global_BB(3)) max(BB(4), global_BB(4))];
	end
	nc_funs('attput', fname, -1, 'BoundingBox', global_BB);

	if (~isempty(outer_polygs))
		for (k = 1:n_swarms)
			if (~nOuterPolys(k)),	continue,	end		% This outer polygon is empty. Jump it
			kk = k + ultimo;
			if (~isa(outer_polygs, 'cell')),	x = outer_polygs(:,1);		y = outer_polygs(:,2);
			else								x = outer_polygs{k}(:,1);	y = outer_polygs{k}(:,2);
			end
			write_var(fname, sprintf('%sPolyOUT_%d', Xname, kk), Nctype, sprintf('dimpolyOUT_%d',kk), long_name{1}, units{1})
			nc_funs('varput', fname, sprintf('%sPolyOUT_%d', Xname, kk), x);
			write_var(fname, sprintf('%sPolyOUT_%d', Yname, kk), Nctype, sprintf('dimpolyOUT_%d',kk), long_name{2}, units{2})
			nc_funs('varput', fname, sprintf('%sPolyOUT_%d', Yname, kk), y);
		end
	end

	if (~isempty(inner_polygs))
		for (k = 1:numel(nInnerPolys))
			if (~nInnerPolys(k)),	continue,	end		% This inner polygon is empty. Jump it
			kk = k + ultimo;
			if (~isa(inner_polygs, 'cell')),	x = inner_polygs(:,1);		y = inner_polygs(:,2);
			else								x = inner_polygs{k}(:,1);	y = inner_polygs{k}(:,2);
			end
			write_var(fname, sprintf('%sPolyIN_%d_%d', Xname, ultimo+1, kk), Nctype, sprintf('dimpolyIN_%d_%d',ultimo+1, kk), long_name{1}, units{1})
			nc_funs('varput', fname, sprintf('%sPolyIN_%d_%d', Xname, ultimo+1, kk), x);
			write_var(fname, sprintf('%sPolyIN_%d_%d', Yname, ultimo+1, kk), Nctype, sprintf('dimpolyIN_%d_%d',ultimo+1, kk), long_name{2}, units{2})
			nc_funs('varput', fname, sprintf('%sPolyIN_%d_%d', Yname, ultimo+1, kk), y);
		end
	end

% ----------------------------------------------------------------------------------------------
function write_poly_type(fname, data, is_polygon_2D, is_polygon_3D, is_polyline_2D, is_polyline_3D, is_geog, is_Pts_cell, ...)
						long_name, units, ultimo, global_BB, n_ensembles, multiSegPos)
% Write variables for the 2D & 3D polygon|polyline cases

	is_3D = false;
	if (is_polygon_2D || is_polygon_3D)
		prefix = 'Polygon_';
		if (is_polygon_3D),		prefix = [prefix(1:end-1) 'Z_'];	is_3D = true;	end
	elseif (is_polyline_2D || is_polyline_3D)
		prefix = 'PolyLine_';
		if (is_polyline_3D),	prefix = [prefix(1:end-1) 'Z_'];	is_3D = true;	end
	end

	if (is_geog),		Xname = 'lon';		Yname = 'lat';
	else				Xname = 'X';		Yname = 'Y';
	end

	for (k = 1:n_ensembles)
		kk = k + ultimo;
		if (is_Pts_cell),	x = data{k}(:,1);		y = data{k}(:,2);
		else				x = data(:,1);			y = data(:,2);
		end
		if ( ~isa(x,'double') && any(isnan(x)) )	% BLOODY ML BUGS (R13)
			BB = [min(double(x)) max(double(x)) min(double(y)) max(double(y))];
		else
			BB = double([min(x) max(x) min(y) max(y)]);
		end
		if (is_3D)		% Split cases between 3d and 2D
			if (is_Pts_cell),	z = data{k}(:,3);
			else				z = data(:,3);
			end
			BB = [BB double([min(z) max(z)])];
		end
		if (k == 1),	Nctype = getDataType(x);		end		% Get data type

		write_var(fname, sprintf('%s%s%d', Xname, prefix, kk), Nctype, sprintf('dimpts_%d',kk), long_name{1}, units{1})
		nc_funs('varput', fname, sprintf('lon%s%d', prefix, kk), x);
		write_var(fname, sprintf('%s%s%d', Yname, prefix, kk), Nctype, sprintf('dimpts_%d',kk), long_name{2}, units{2})
		nc_funs('varput', fname, sprintf('lat%s%d', prefix ,kk), y);
		if (is_3D)
			write_var(fname, sprintf('Z%s%d', prefix, kk), Nctype, sprintf('dimpts_%d',kk), 'z', 'meters')
			nc_funs('varput', fname, sprintf('Z%s%d', prefix, kk), z);
		end

		% Create a container variable to hold the BB info
		write_var(fname, sprintf('%s%d', prefix, kk), 2)
		nc_funs('attput', fname, sprintf('%s%d', prefix, kk), 'BoundingBox', BB);
		if (k == 1 && ~isempty(multiSegPos))
			nc_funs('attput', fname, sprintf('%s%d', prefix, kk), 'SegmentBoundaries', multiSegPos);
		end
		
		% Update the global BB
		global_BB = [min(BB(1), global_BB(1)) max(BB(2), global_BB(2)) min(BB(3), global_BB(3)) max(BB(4), global_BB(4))];
	end
	nc_funs('attput', fname, -1, 'BoundingBox', global_BB);

% ----------------------------------------------------------------------------------------------
function write_var(fname, name, tipo, dim, long_name, units, actual_range, comment, fillValue, missing_value, scale_factor)	
% NARGIN == 3 & NARGIN == 6 are special handled cases
	if (nargin == 3)
		dim = [];	long_name = [];	units = [];	actual_range = [];	comment = [];	fillValue = [];	missing_value = [];	scale_factor = [];
	elseif (nargin == 6)
		actual_range = [];	comment = [];	fillValue = [];	missing_value = [];	scale_factor = [];
	end
	varstruct.Name = name;
	varstruct.Nctype = tipo;
	if (~isempty(dim)),		varstruct.Dimension = {dim};	end
	nc_funs('addvar', fname, varstruct)
	if (~isempty(long_name)),		nc_funs('attput', fname, varstruct.Name, 'long_name', long_name ),	end
	if (~isempty(units)),			nc_funs('attput', fname, varstruct.Name, 'units', units),			end
	if (~isempty(comment)),			nc_funs('attput', fname, varstruct.Name, 'comment', comment),		end
	if (~isempty(actual_range)),	nc_funs('attput', fname, varstruct.Name, 'actual_range', actual_range),	end
	if (~isempty(fillValue)),		nc_funs('attput', fname, varstruct.Name, '_FillValue', fillValue),		end
	if (~isempty(missing_value)),	nc_funs('attput', fname, varstruct.Name, 'missing_value', missing_value),	end
	if (~isempty(scale_factor)),	nc_funs('attput', fname, varstruct.Name, 'scale_factor', scale_factor),	end

% ----------------------------------------------------------------------------------------------
function Nctype = getDataType(data)	
	
	switch ( class(data) )
		case 'single',		Nctype = 5;			% NC_FLOAT
		case 'int16',		Nctype = 3;			% NC_SHORT
		case 'int32',		Nctype = 4;			% NC_INT
		case 'int8',		Nctype = 1;			% NC_BYTE
		case 'uint8',		Nctype = 1;			% NC_CHAR
		case 'double',		Nctype = 6;			% NC_DOUBLE
		otherwise
			error('SHAPENC:getDataType', ['Unsuported data type: ' class(data)])
	end


% ----------------------------------------------------------------------------------------------
function [data, pos, is_polygon_2D, is_polygon_3D, is_polyline_2D, is_polyline_3D, is_point_2D] = ...
		readFile(fname, nMaxChunks, is_polygon_2D, is_polygon_3D, is_polyline_2D, is_polyline_3D, is_point_2D)
% Try to read from a file (if ASCII or from a shapefile is binary)
% The POS output variable (when not empty) contains the positions of the multiseg separators (NaNs)
% when DATA was packed into a single multisegment array.
%
% The several IS_XXX variables are currently only used/checked when FNAME is a shapefile file name

	pos = [];
	[bin,n_column,multi_seg,n_headers] = guess_file(fname);

	if (bin ~= 0)   % NOT ASCII
		% Try to see if it is a shapefile
		try		[s,t] = mex_shape(fname);
		catch
			error('SHAPENC:readFile',[lasterr ' - NOT a shapefile and reading other binary files is not available'])
		end

		% So, it was a shapefile. Minimalist support
		is_3D = false;
		switch lower(t)
			case 'polygon'
				is_polygon_2D = true;
			case 'polygonz'
				is_polygon_3D = true;		is_3D = true;
			case 'polyline'
				is_polyline_2D = true;
			case 'polylinez'
				is_polyline_3D = true;		is_3D = true;
			case 'point'
				is_point_2D = true;
		end

		nChunks = numel(s);
		if (nChunks > nMaxChunks)
			pos = zeros(1,nChunks+1);
			pos(1) = 1;
			n = 0;
			for (k = 1:nChunks)						% Count total number of data points
				n = n + numel(s(k).X);
				pos(k+1) = n + k+1;					% Positions that will hold the NaN multiseg separators
			end

			nCols = 2;
			if (is_3D),		nCols = 3;	end
			data = single(zeros(n+nChunks-1,1)*NaN);		% Pre-allocate
			data = repmat(data,1,nCols);

			for (k = 1:nChunks)
				ini = pos(k)+1;		fim = pos(k+1)-1;
				data(ini:fim,1) = single(s(k).X);
				data(ini:fim,2) = single(s(k).Y);
				if (is_3D),		data(ini:fim,3) = single(s(k).Z);	end
			end
			pos(end) = [];		% Last point was in excess (for algo convenience only)
			pos = int32(pos);

		else
			% pack the individual chunks into a cell array
			data = cell(nChunks);
			for (k = 1:nChunks)
				x = single(s(k).X);		y = single(s(k).Y);
				data{k}(:,1) = x(:);	data{k}(:,2) = y(:);
				if (is_3D)
					y = single(s(k).Z);
					data{k}(:,3) = y(:);
				end
			end
		end

		return
	end

	if isempty(bin) && isempty(n_column) && isempty(multi_seg) && isempty(n_headers)
		error('SHAPENC:readFile',['Error reading file ' fname])
	end
	if (n_column < 2)
		error('SHAPENC:readFile','File error. The file doesn''t have at least 2 columns')
	end
	if (isempty(n_headers)),    n_headers = NaN;    end
	if (multi_seg)
		data = text_read(fname,NaN,n_headers,'>');
		for (k = 1:numel(data))
			data{k} = single(data{k});
		end
	else
		%fid = fopen('lisboa_igoe_geo.bin');		data=fread(fid,'*float');	fclose(fid);	data=reshape(data,3,numel(data)/3)';
		data = single(text_read(fname,NaN,n_headers));
	end

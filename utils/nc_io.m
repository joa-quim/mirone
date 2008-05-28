function varargout = nc_io(fname, mode, handles, data, misc)
% HANDLES, DATA & MISC only matters on writing
% On reading -> MODE == 'r' & varargout = [X,Y,Z,head,misc] = nc_io(fname, mode);
%	where:	X, Y, Z are the coords row vectors and the data 2D array
%			HEAD = [x_min x_max y_min y_max z_min z_max node_offset x_inc y_inc]
%			MISC is a struct with 'desc', 'title', 'history', 'srsWKT', 'strPROJ4' fields
%
%	Special cases on reading
% 		MODE == 'R' ==> varargout = [X,Y,[],head,misc] = nc_io(fname, mode);
% 		MODE == 'h' ==> varargout = [head,misc] = nc_io(fname, mode);
%
% On writing -> nc_io(fname, mode, handles, data, misc)
%			MODE == 'w' and DATA the 2D array to be saved
%			MISC is optional and if provided it must contain the fields as declared below.
%
%	Special case of multi-Layer writing:
%			MODE == 'wN/levelsName' initialize the netCDF file for writing a 3D file.
%					Where 'N' is the number of levels and 'levelsName' the name of thirth
%					dimension variable.
%					If N == 0 or negative 3rth dim is UNLIMITED
%					If N < 0 it is taken as the first 3rth dim numeric value
%					Example 'w10/time' or 'w-4.5/time.
%			MODE == 'wK' or 'wK\thisLevel' use in subsequent call to save the K'th level.
%					The second form transmits the 3rth dim numeric value. If this is not
%					provided thisLevel = K is used.
%					ATTENTION: K is zero based
%
%				In case of a LIMITED var the levels vector may be transmitted via the handles structure.
%				If absent, a default one with (1:N) is created. This information may be transmited 
%				only on the last call	(for the cases of only than all values of that vector are known)

%	AUTHOR
%		Joaquim Luis  - 15-October-2007
%		jluis@ualg.pt - Universidade do Algarve
	
	if ( nargin == 4 && ~isempty(data) )
		misc = struct('x_units',[],'y_units',[],'z_units',[],'z_name',[],'desc',[], ...
			'title',[],'history',[],'srsWKT',[], 'strPROJ4',[]);
	end

	if (mode(1) == 'w')			% Write file
		if (numel(mode) == 1)	% Called in the normal, write once, mode
			write_nc(fname, handles, data, misc)
		else					% Called in the 'append' mode.
			ind = strfind(mode,'/');
			if (~isempty(ind))					% Initialize a new 3D file
				page.nLevels = str2double(mode(2:ind(1)-1));	% If it is 0 or a negative float -> 3D unlimited
				page.levelName = mode(ind(1)+1:end);
			else								% Append a new page to an existing 3D file

				indLev = strfind(mode,'\');		% Search for a level (3rth dim) value
				if (~isempty(indLev))
					page{2} = str2double(mode(indLev(1)+1:end));
					mode(indLev(1):end) = [];		% Rip the level info from the mode string
				end
				tmp = abs( round(str2double(mode(2:end))) );
				if (isnan(tmp)),	error('NC_IO:write_nc','Nonsense in PAGE number.'),		end
				page{1} = tmp;
			end

			write_nc(fname, handles, data, misc, page)
		end
	elseif (mode(1) == 'R')		% Get all but Z
		[varargout{1} varargout{2} varargout{3} varargout{4} varargout{5}] = read_nc(fname, 1);
	elseif (mode(1) == 'h')		% Get header[misc]
		[varargout{1} varargout{2}] = read_nc(fname);
	else
		[varargout{1} varargout{2} varargout{3} varargout{4} varargout{5}] = read_nc(fname);
	end

% _________________________________________________________________________________________________	
% -*-*-*-*-*-*-$-$-$-$-$-$-#-#-#-#-#-#-%-%-%-%-%-%-@-@-@-@-@-@-(-)-(-)-(-)-&-&-&-&-&-&-{-}-{-}-{-}-
function write_nc(fname, handles, data, misc, page)

	if (nargin > 4),	
		persistent levelName z_name is_unlimited
	end
	if (nargin == 4)							% Normal, write once, mode
		is3D = false;
	elseif (nargin == 5 && isa(page,'struct'))	% Initialize a new 3D file
		is_unlimited = false;
		levelName = page.levelName;
		nLevels = page.nLevels;
		if (nLevels <= 0)						% UNLIMITED
			is_unlimited = true;
			if (nLevels < 0)
				levelVec = -nLevels;			% First value of the unlimited var
				nLevels = 0;
			else		% == 0
				levelVec = 1;					% At least it won't be 9...10^36
			end
		elseif (isfield(handles,'levelVec') && ~isempty(handles.levelVec))
			levelVec = handles.levelVec;		% Only for not UNLIMITED
		else
			levelVec = 1:nLevels;				%			"
		end
		is3D = true;
	elseif (nargin == 5 && isa(page,'cell'))	% 'append' layer mode. Work and return
		if (is_unlimited)
			if (numel(page) == 2)
				level = page{2};
			else
				level = page{1};
			end
			nc_funs('varput', fname, levelName, level, page{1}, 1 );	% The UNLIMITED var value
		end

		if (ndims(data) == 2)
			[ny, nx] = size(data);
			nc_funs('varput', fname, z_name, reshape(data,[1 ny nx]), [page{1} 0 0], [1 ny nx] );
		elseif (ndims(data) == 3 && size(data,1) == 1)
			nc_funs('varput', fname, z_name, data, [page{1} 0 0], [1 ny nx] );
		else
			error('NC_IO:write_nc','input array on append mode must be MxN or 1xMxN')
		end

		% Update global min/max
		if ( isa(data, 'double') )
			mima = [min(data(:)) max(data(:))];
		else		% min/max are bugged when NaNs in singles
			zz = grdutils(data,'-L');  mima = [zz(1) zz(2)];
		end
		z_actual_range = nc_funs('attget', fname, z_name, 'actual_range');
		if (~isempty(z_actual_range))
			mima = [min(z_actual_range(1),mima(1)) max(z_actual_range(2),mima(2))];
		end
		nc_funs('attput', fname, z_name, 'actual_range', mima);
		if (isfield(handles,'levelVec') && ~isempty(handles.levelVec))	% We may have updated info on this
			nc_funs('varput', fname, levelName, handles.levelVec );
		end
		return
	end

	nc_funs('create_empty', fname)

	if (handles.geog)
		x_var = 'longitude';		y_var = 'latitude';
		misc.x_units = 'degrees_east';		misc.y_units = 'degrees_north';		% Last word is here
	else
		x_var = 'x';				y_var = 'y';
		if (isempty(misc.x_units))
			misc.x_units = 'unknown';		misc.y_units = 'unknown';
		end
	end

	% ---------------------- Make sure we always have something assigned to these vars --------------
	x_units = misc.x_units;			y_units = misc.y_units;
	if (isempty(misc.z_name)),		z_name = 'z';
	else							z_name = misc.z_name;
	end
	if (isempty(misc.z_units)),		z_units = 'unknown';
	else							z_units = misc.z_units;
	end
	if (isempty(misc.desc)),		desc = 'File written from Matlab';
	else							desc = misc.desc;
	end
	if (isempty(misc.title)),		tit = 'Grid created by Mirone';
	else							tit = misc.title;
	end
	% ----------------------------------------------------------------
	
	% See if we have them in appdata
	if ( isempty(misc.srsWKT) ),		misc.srsWKT = getappdata(handles.figure1,'ProjWKT');	end
	if ( isempty(misc.strPROJ4) ),		misc.strPROJ4 = getappdata(handles.figure1,'Proj4');	end

	% Create the coordinates vectors
	nx = round(diff(handles.head(1:2)) / handles.head(8) + ~handles.head(7));
	ny = round(diff(handles.head(3:4)) / handles.head(9) + ~handles.head(7));
	X = linspace(handles.head(1), handles.head(2), nx);
	Y = linspace(handles.head(3), handles.head(4), ny);

	% ---------------------------- Write the dimensions --------------------------------
	nc_funs('add_dimension', fname, x_var, nx )
	nc_funs('add_dimension', fname, y_var, ny )

	if (is3D)		% Initialize a 3D file
 		nc_funs('add_dimension', fname, levelName, nLevels)
	end

	x_varstruct.Name = x_var;		x_varstruct.Dimension = {x_var};
	y_varstruct.Name = y_var;		y_varstruct.Dimension = {y_var};
	nc_funs('addvar', fname, x_varstruct)
	nc_funs('addvar', fname, y_varstruct)
	varstruct.Dimension = {y_var, x_var};

	if (is3D)		% Initialize a 3D file
		t_varstruct.Name = levelName;		t_varstruct.Dimension = {levelName};
		nc_funs('addvar', fname, t_varstruct)
		varstruct.Dimension = {levelName, y_var, x_var};
	end

	varstruct.Name = z_name;
	add_off = [];
	switch ( class(data) )
		case 'single'			% NC_FLOAT
			if (handles.was_int16 && ~handles.computed_grid)
				if (handles.have_nans),		data(isnan(data)) = -32768;		end
				data = int16(data);
				varstruct.Nctype = 3;		no_val = -32768;
			else
				varstruct.Nctype = 5;		no_val = single(nan);
			end
		case 'int16'			% NC_SHORT
			varstruct.Nctype = 3;		no_val = -32768;
		case 'int32'			% NC_INT
			varstruct.Nctype = 4;		no_val = -2147483648;
		case 'int8'				% NC_BYTE
			varstruct.Nctype = 1;		no_val = [];	add_off = 128;
		case 'uint8'			% NC_CHAR
			data = grdutils(data,'-C');
			varstruct.Nctype = 1;		no_val = [];	add_off = 128;
		case 'double'			% NC_DOUBLE
			varstruct.Nctype = 6;		no_val = nan;
		otherwise
			error('NC_IO:write_nc', ['Unsuported data type: ' class(data)])
	end
	nc_funs('addvar', fname, varstruct )

	% -------------------------- Globals --------------------------
	nc_global = -1;
	nc_funs('attput', fname, nc_global, 'Conventions', 'COARDS' );
	nc_funs('attput', fname, nc_global, 'title', tit );
	if (~isempty(misc.history)),	nc_funs('attput', fname, nc_global, 'history', misc.history);	end
	nc_funs('attput', fname, nc_global, 'description', desc );
	nc_funs('attput', fname, nc_global, 'node_offset', handles.head(7) );
	nc_funs('attput', fname, nc_global, 'Source_Software', 'Mirone' );

	% ------------------------------ Write the variables ------------------------------------
	% ------- Put the coords vectors ---
	nc_funs('varput', fname, x_var, X );
	nc_funs('varput', fname, y_var, Y );
	if (is3D && ~is_unlimited),	nc_funs('varput', fname, levelName, levelVec );		end		% not UNLIMITED

	nc_funs('attput', fname, x_var, 'long_name', x_var );
	nc_funs('attput', fname, x_var, 'units', x_units);
	nc_funs('attput', fname, x_var, 'actual_range', [X(1) X(end)] );

	nc_funs('attput', fname, y_var, 'long_name', y_var );
	nc_funs('attput', fname, y_var, 'units', y_units);
	nc_funs('attput', fname, y_var, 'actual_range', [Y(1) Y(end)] );

	nc_funs('attput', fname, z_name, 'long_name', z_name );
	if (~isempty(no_val)),		nc_funs('attput', fname, z_name, '_FillValue', no_val );	end
	if (~isempty(add_off)),		nc_funs('attput', fname, z_name, 'add_offset', add_off);	end
	nc_funs('attput', fname, z_name, 'actual_range', handles.head(5:6) );
	nc_funs('attput', fname, z_name, 'units', z_units);
	
	if ( ~isempty(misc.srsWKT) || ~isempty(misc.strPROJ4) )
		% Create a container variable named "grid_mapping" to hold the projection info
		nc_funs('attput', fname, z_name, 'grid_mapping', 'grid_mapping');
		nc_funs('addvar', fname, struct('Name','grid_mapping', 'Nctype',2))		% 2 -> char
		
		if (~isempty(misc.srsWKT))
			nc_funs('attput', fname, 'grid_mapping', 'spatial_ref', misc.srsWKT);
		end
		if (~isempty(misc.strPROJ4))
			nc_funs('attput', fname, 'grid_mapping', 'spatial_ref', misc.strPROJ4);
		end
	end

	if (is3D && is_unlimited),	nc_funs('varput', fname, levelName, levelVec(1), 0, 1 );		end		% The UNLIMITED var value
	if (ndims(data) == 2)
		nc_funs('varput', fname, z_name, data, [0 0], [ny nx] );
	elseif (ndims(data) == 3)
		nc_funs('varput', fname, z_name, data, [0 0 0], [1 ny nx] );
	else
		error('NC_IO: array must be 2 or 3D only')
	end

% _________________________________________________________________________________________________	
% -*-*-*-*-*-*-$-$-$-$-$-$-#-#-#-#-#-#-%-%-%-%-%-%-@-@-@-@-@-@-(-)-(-)-(-)-&-&-&-&-&-&-{-}-{-}-{-}-
function [X,Y,Z,head,misc] = read_nc(fname, opt)
	% Read ...
	% If OPT, do not load Z (return Z = []) -- Used by Aquamoto to read multi-levels grids

	if (nargin == 1),		get_Z = true;
	else					get_Z = false;		Z = [];
	end

	s = nc_funs('info',fname);
	dimNames = {s.Dimension.Name};
	if ( any(strcmp(dimNames, 'xysize')) )			% Old style netCDF gmt grid format
		[X,Y,Z,head,misc] = read_old_cdf(fname, s);
		return
	end

	nvars = numel(s.Dataset);		z_id = 0;
	for (k = 1:nvars)			% Find first 2-dimensional (z) variable
		if (numel(s.Dataset(k).Size) >= 2 )
			z_id = k;
			break
		end
	end

	if (~z_id),		error('NC_IO:read_nc', 'Didn''t find any (at least) 2D variable.'),		end

	ndims = numel(s.Dataset(z_id).Size);	% z number of dimensions
	dims = s.Dataset(z_id).Dimension;		% Variable names of dimensions z variable - ORDER IS CRUTIAL

	% --------------------- Fish in the attribs of the Z var --------------------
	if (~isempty(s.Dataset(z_id).Attribute))
		attribNames = {s.Dataset(z_id).Attribute.Name};
	else
		attribNames = [];
	end
	dataNames = {s.Dataset.Name};			% example: {'x' 'y' 'z'}; Order may change
	ind = strcmp(attribNames,'actual_range');		z_actual_range = [];
	if (any(ind)),	z_actual_range = s.Dataset(z_id).Attribute(ind).Value;	end
	ind = strcmp(attribNames,'scale_factor');		scale_factor = 1;
	if (any(ind)),	scale_factor = s.Dataset(z_id).Attribute(ind).Value;	end
	ind = strcmp(attribNames,'add_offset');			add_offset = 0;
	if (any(ind)),	add_offset = s.Dataset(z_id).Attribute(ind).Value;		end
	% GDAL and ESRI put the eventual georeferencing on the layer
	grid_mapping = [];		srsWKT = [];		strPROJ4 = [];
	ind = strcmp(attribNames,'grid_mapping');
	if (any(ind))
		grid_mapping = s.Dataset(z_id).Attribute(ind).Value;
		ind = strcmp(dataNames, grid_mapping);
		containor_id = find(ind);
		containerNames = {s.Dataset(containor_id).Attribute.Name};
		ind = strcmp(containerNames,'spatial_ref');
		if (any(ind)),	srsWKT = s.Dataset(containor_id).Attribute(ind).Value;		end
		
		ind = strcmp(containerNames,'strPROJ4');
		if (any(ind)),	strPROJ4 = s.Dataset(containor_id).Attribute(ind).Value;		end
	end
	if (isempty(srsWKT))
		ind = strcmp(attribNames,'esri_pe_string');			% ESRI uses this
		if (any(ind)),	srsWKT = s.Dataset(z_id).Attribute(ind).Value;		end
	end
	% ----------------------------------------------------------------------------

	% ----------- Get the ids of the x and y (and depth and time) variables ------
	ind = strcmp(dataNames, dims{end});
	x_id = find(ind);
	ind = strcmp(dataNames, dims{end-1});
	y_id = find(ind);
	
	x_actual_range = [];		y_actual_range = [];
	if ( ~isempty(x_id) )
		% ------------------ Get the X & Y ranges ------------------------------------
		attribNames = {s.Dataset(x_id).Attribute.Name};
		ind = strcmp(attribNames,'actual_range');
		if (any(ind)),	x_actual_range = s.Dataset(x_id).Attribute(ind).Value;	end
		%ind = strcmp(attribNames,'units');		units = [];
		%if (any(ind))
			%units = s.Dataset(x_id).Attribute(ind).Value;
			%data.geog = 0;
			%if (strncmp(units,'degree',6)),		data.geog = 1;	end
		%end
	end
	if ( ~isempty(y_id) )
		attribNames = {s.Dataset(y_id).Attribute.Name};
		ind = strcmp(attribNames,'actual_range');
		if (any(ind)),	y_actual_range = s.Dataset(y_id).Attribute(ind).Value;	end
	end
	% ----------------------------------------------------------------------------
	
	% --------------------- Fish the Global attributes ---------------------------
	if (~isempty(s.Attribute))
		attribNames = {s.Attribute.Name};
	else
		attribNames = [];
	end
	ind = strcmp(attribNames,'node_offset');
	if (any(ind)),	node_offset = s.Attribute(ind).Value;
	else			node_offset = 0;				% Hmmm, ...
	end
	ind = strcmp(attribNames,'title');					titulo = [];
	if (any(ind)),	titulo = s.Attribute(ind).Value;	end
	ind = strcmp(attribNames,'history');				history = [];
	if (any(ind)),	history = s.Attribute(ind).Value;	end
	ind = strcmp(attribNames,'description');			description = [];
	if (any(ind)),	description = s.Attribute(ind).Value;	end
	
	misc = struct('desc',description, 'title',titulo, 'history',history, 'srsWKT',srsWKT, 'strPROJ4',[]);
	if (~isempty(strPROJ4)),	misc.strPROJ4 = strPROJ4;	end
	% ----------------------------------------------------------------------------
	
	% --------------------- Finally, get the Data --------------------------------
	% Now z may be > 2D. If it is, get the first layer
	z_dim = s.Dataset(z_id).Size;
	if (get_Z)
		nD = numel(z_dim);
		if (nD == 2 )
			Z = nc_funs('varget', fname, s.Dataset(z_id).Name);
		else
			if ( ~all([s.Dimension.Length]) )
				warndlg('One ore more dimensions has ZERO length. Expect #&%%&&$.','Warning')
			end
			Z = nc_funs('varget', fname, s.Dataset(z_id).Name, zeros(1,nD), [ones(1,nD-2) s.Dataset(z_id).Size(end-1:end)]);
		end
		%[ny, nx] = size(Z);
	else
		misc.z_id = z_id;		% Return the z_id so that we don't have to repeat the fishing process
		misc.z_dim = z_dim;
	end
	nx = z_dim(end);		ny = z_dim(end-1);

	if (~isempty(x_id)),	X = double(nc_funs('varget', fname, s.Dataset(x_id).Name));
	else					X = 1:nx;
	end
	if (~isempty(y_id)),	Y = double(nc_funs('varget', fname, s.Dataset(y_id).Name));
	else					Y = 1:ny;
	end

	if (~isempty(x_actual_range)),		head(1:2) = x_actual_range;
	else								head(1:2) = [X(1) X(end)];
	end
	if (~isempty(y_actual_range)),		head(3:4) = y_actual_range;
	else								head(3:4) = [Y(1) Y(end)];
	end

	if (get_Z && (scale_factor ~= 1 || add_offset ~= 0) )		% If we have scale or offset convert to single the "lower"
		if ( ~(isa(Z, 'double') || isa(Z, 'single')) )			% types. Maybe not always apropriate but not general rule
			Z = single(Z);
		end
		cvlib_mex('CvtScale',Z,scale_factor,add_offset)			% Do inplace
	end

	if (get_Z && isempty(z_actual_range))
		if ( isa(Z, 'double') )
			z_actual_range = [min(Z(:)) max(Z(:))];
		else		% min/max are bugged when NaNs in singles
			zz = grdutils(Z,'-L');  z_actual_range = [zz(1) zz(2)];
		end
	elseif (~get_Z)
		z_actual_range = [0 1];		% Dumb, but it is not meant to be used anywhere
	end
	head(5:7) = double([z_actual_range node_offset]);
	head(8) = diff(head(1:2)) / (nx - ~node_offset);
	head(9) = diff(head(3:4)) / (ny - ~node_offset);
	
	if (nargout <= 2)
		X = head;	Y = misc;
	end

% _________________________________________________________________________________________________	
% -*-*-*-*-*-*-$-$-$-$-$-$-#-#-#-#-#-#-%-%-%-%-%-%-@-@-@-@-@-@-(-)-(-)-(-)-&-&-&-&-&-&-{-}-{-}-{-}-
function [X,Y,Z,head,misc] = read_old_cdf(fname, s)
	% Read old GMT netCDF format

	% Fish the IDs of	x_range, y_range, z_range, spacing, dimension & z
	dataNames = {s.Dataset.Name};
	x_range_id = find(strcmp(dataNames, 'x_range'));
	y_range_id = find(strcmp(dataNames, 'y_range'));
	z_range_id = find(strcmp(dataNames, 'z_range'));
	spacing_id = find(strcmp(dataNames, 'spacing'));
	dimension_id = find(strcmp(dataNames, 'dimension'));
	z_id = find(strcmp(dataNames, 'z'));

	% --------------------- Fish in the attribs of the Z var ----------------------
	attribNames = {s.Dataset(z_id).Attribute.Name};
	ind = strcmp(attribNames,'node_offset');		node_offset = 0;
	if (any(ind)),	node_offset = s.Dataset(z_id).Attribute(ind).Value;		end
	% ------------------------------------------------------------------------------

	% --------------------- Fish the Global attributes -----------------------------
	attribNames = {s.Attribute.Name};
	ind = strcmp(attribNames,'title');					titulo = [];
	if (any(ind)),	titulo = s.Attribute(ind).Value;	end
	ind = strcmp(attribNames,'source');					source = [];
	if (any(ind)),	source = s.Attribute(ind).Value;	end
	
	% --------------------- Finally, get the Data --------------------------------
	x_range = double(nc_funs('varget', fname, s.Dataset(x_range_id).Name));
	y_range = double(nc_funs('varget', fname, s.Dataset(y_range_id).Name));
	z_range = double(nc_funs('varget', fname, s.Dataset(z_range_id).Name));
	spacing = double(nc_funs('varget', fname, s.Dataset(spacing_id).Name));
	dimension = double(nc_funs('varget', fname, s.Dataset(dimension_id).Name));
	Z = nc_funs('varget', fname, s.Dataset(z_id).Name);
	
	nx = round(diff(x_range) / spacing(1) + ~node_offset);
	ny = round(diff(y_range) / spacing(2) + ~node_offset);
	Z = reshape(Z,nx,ny);
	Z = Z.';
	Z = flipud(Z);

	if (node_offset)		% Pixel registration
		X = linspace(x_range(1)+spacing(1)/2, x_range(2)-spacing(1)/2, dimension(1));
		Y = linspace(y_range(1)+spacing(2)/2, y_range(2)-spacing(2)/2, dimension(2));
	else
		X = linspace(x_range(1), x_range(2), dimension(1));
		Y = linspace(y_range(1), y_range(2), dimension(2));
	end
	head = double([x_range(:)' y_range(:)' z_range(:)' node_offset spacing(:)']);
	misc = struct('desc',[], 'title',titulo, 'history',source, 'srsWKT',[], 'strPROJ4',[]);


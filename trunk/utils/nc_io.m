function varargout = nc_io(fname, mode, handles, data, misc)
% HANDLES, DATA & MISC only matters on writing
% On reading -> MODE == 'r' & varargout = [X,Y,Z,head,misc] = nc_io(fname, mode);
%	where:	X, Y, Z are the coords row vectors and the data 2D array
%			HEAD = [x_min x_max y_min y_max z_min z_max node_offset x_inc y_inc]
%			MISC is a struct with 'desc', 'title', 'history', 'srsWKT', 'strPROJ4' fields
%
% On writing -> nc_io(fname, mode, handles, data, misc)
%			MODE == 'w' and DATA the 2D array to be saved
%			MISC is optional and if provided it must contain the fields as declared below.

%	AUTHOR
%		Joaquim Luis  - 15-October-2007
%		jluis@ualg.pt - Universidade do Algarve
	
	if ( nargin == 4 && ~isempty(data) )
		misc = struct('x_units',[],'y_units',[],'z_units',[],'z_name',[],'desc',[], ...
			'title',[],'history',[],'srsWKT',[], 'strPROJ4',[]);
	end

	if (mode(1) == 'w')
		write_nc(fname, handles, data, misc)
	elseif (mode(1) == 'R')		% Get all but Z
		[varargout{1} varargout{2} varargout{3} varargout{4} varargout{5}] = read_nc(fname, 1);
	else
		[varargout{1} varargout{2} varargout{3} varargout{4} varargout{5}] = read_nc(fname);
	end

% _________________________________________________________________________________________________	
% -*-*-*-*-*-*-$-$-$-$-$-$-#-#-#-#-#-#-%-%-%-%-%-%-@-@-@-@-@-@-(-)-(-)-(-)-&-&-&-&-&-&-{-}-{-}-{-}-
function write_nc(fname, handles, data, misc)

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
	nc_funs('add_dimension', fname, 'x', nx )
	nc_funs('add_dimension', fname, 'y', ny )

	x_varstruct.Name = 'x';		x_varstruct.Dimension = {'x'};
	y_varstruct.Name = 'y';		y_varstruct.Dimension = {'y'};
	nc_funs('addvar', fname, x_varstruct)
	nc_funs('addvar', fname, y_varstruct)

	varstruct.Dimension = {'y', 'x'};
	varstruct.Name = 'z';
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
	nc_funs('varput', fname, 'x', X );
	nc_funs('varput', fname, 'y', Y );

	nc_funs('varput', fname, 'z', data, [0 0], [ny nx] );

	nc_funs('attput', fname, 'x', 'long_name', x_var );
	nc_funs('attput', fname, 'x', 'units', x_units);
	nc_funs('attput', fname, 'x', 'actual_range', [X(1) X(end)] );

	nc_funs('attput', fname, 'y', 'long_name', y_var );
	nc_funs('attput', fname, 'y', 'units', y_units);
	nc_funs('attput', fname, 'y', 'actual_range', [Y(1) Y(end)] );

	nc_funs('attput', fname, 'z', 'long_name', z_name );
	if (~isempty(no_val)),		nc_funs('attput', fname, 'z', '_FillValue', no_val );	end
	if (~isempty(add_off)),		nc_funs('attput', fname, 'z', 'add_offset', add_off);	end
	nc_funs('attput', fname, 'z', 'actual_range', handles.head(5:6) );
	nc_funs('attput', fname, 'z', 'units', z_units);
	
	if ( ~isempty(misc.srsWKT) || ~isempty(misc.strPROJ4) )
		% Create a container variable named "grid_mapping" to hold the projection info
		nc_funs('attput', fname, 'z', 'grid_mapping', 'grid_mapping');
		nc_funs('addvar', fname, struct('Name','grid_mapping', 'Nctype',2))		% 2 -> char
		
		if (~isempty(misc.srsWKT))
			nc_funs('attput', fname, 'grid_mapping', 'spatial_ref', misc.srsWKT);
		end
		if (~isempty(misc.strPROJ4))
			nc_funs('attput', fname, 'grid_mapping', 'spatial_ref', misc.strPROJ4);
		end
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
% 	ind = strcmp(attribNames,'scale_factor');		scale_factor = 1;
% 	if (any(ind)),	scale_factor = s.Dataset(z_id).Attribute(ind).Value;	end
% 	ind = strcmp(attribNames,'add_offset');			add_offset = 0;
% 	if (any(ind)),	add_offset = s.Dataset(z_id).Attribute(ind).Value;		end
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
	attribNames = {s.Attribute.Name};
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

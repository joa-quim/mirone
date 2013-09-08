function varargout = nc_io(fname, mode, handles, data, misc)
% Wraper function to netCDF I/O
%
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
%	Special case of multi-Layer writing:	PAY ATTENTION TO THE SLASHES
%
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
%				This is done by setting the the calling function
%					handles.levelVec = levelVec;
%				If absent, a default one with (1:N) is created. This information may be transmitted 
%				only on the last call (for the cases of only than all values of that vector are known)
%
%				The above is extremely important for the standalone version. For an highly mysterious
%				reason the mexnc call in nc_funs/write_the_data() errors when writing UNLIMITED variables
%				so the only way out is to send in the levelVec (vector of times, most of times)

%	Copyright (c) 2004-2013 by J. Luis
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

% $Id$

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

	if (nargin > 4)
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
			if (numel(page) == 2),		level = page{2};
			else						level = page{1};
			end
			nc_funs('varput', fname, levelName, level, page{1}, 1 );	% The UNLIMITED var value
		end

		ny = size(data, 1);		nx = size(data, 2);
		if (ndims(data) == 2)
			nc_funs('varput', fname, z_name, reshape(data,[1 ny nx]), [page{1} 0 0], [1 ny nx] );
		elseif (ndims(data) == 3 && size(data,1) == 1)
			nc_funs('varput', fname, z_name, data, [page{1} 0 0], [1 ny nx] );
		else
			error('NC_IO:write_nc','input array on append mode must be MxN or 1xMxN')
		end

		% Update global min/max
		if ( isa(data, 'single') )			% min/max is bugged when NaNs in singles
			zz = grdutils(data,'-L');		mima = [zz(1) zz(2)];
		else
			mima = [double(min(data(:))) double(max(data(:)))];
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

	% ---------------- Decide if netCDF-4 or classic based on deflation level -----------------------
	deflation_level = get_deflation(handles);
	if (deflation_level)
		nc_funs('create_empty', fname, 'netcdf4-classic')	% Despite 'classic' it will be compressed
	else
		nc_funs('create_empty', fname)
	end
	% -----------------------------------------------------------------------------------------------

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
	try			% Wrap with a try because we want to be able to use this fun even when no figure1
		if ( isempty(misc.srsWKT) ),		misc.srsWKT = getappdata(handles.figure1,'ProjWKT');	end
		if ( isempty(misc.strPROJ4) ),		misc.strPROJ4 = getappdata(handles.figure1,'Proj4');	end
	end

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
	add_off = [];	scale_factor = [];
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
			varstruct.Nctype = 1;		no_val = [];
			if ( min(data(:)) < 0 ),	add_off = 128;	end		% Apply the offset trick only when needed
		case 'uint8'			% NC_CHAR
			if ( (ndims(data) == 3) && (size(data,1) == 1) )
				thisSize = size(data);
				data = grdutils(squeeze(data),'-C');
				data = reshape(data, thisSize);
			else
				data = grdutils(data,'-C');
			end
			varstruct.Nctype = 1;		no_val = [];	add_off = 128;
			handles.head(5:6) = double([min(data(:)) max(data(:))]);		% We dont save handles
		case 'double'			% NC_DOUBL
			varstruct.Nctype = 6;		no_val = nan;
		otherwise
			error('NC_IO:write_nc', ['Unsuported data type: ' class(data)])
	end

	% -------------------- The Attributes of the 'z' variable ----------------------------
	% We now do this here so that the attributes are written at the same time that the 'z'
	% variable is created and thus avoid the BAD netCDF 187 bug that ends up in a crash.
	n = 0;
	if (~isempty(no_val))
		if (isnan(no_val) && varstruct.Nctype == 5)		% The following chunk it's because TMW is so shameless
			try											% that for the ML compiler (the TRUE compiler)
				IamCompiled = handles.IamCompiled;		% class(single(NaN)) = double
			catch
				try			s.s = which('mirone');		IamCompiled = false;
				catch,		IamCompiled = true;
				end
			end
			if (IamCompiled)	% OK here we have it all fucked because it's a double NaN, which is unaceptable
				no_val = alloc_mex(1,1,'single',nan);
			end
		end
		varstruct.Attribute(1).Name = '_FillValue';
		varstruct.Attribute(1).Value = no_val;
		n = n + 1;
	end
	varstruct.Attribute(n+1).Name = 'long_name';		varstruct.Attribute(n+1).Value = z_name;
	varstruct.Attribute(n+2).Name = 'actual_range';		varstruct.Attribute(n+2).Value = handles.head(5:6);
	varstruct.Attribute(n+3).Name = 'units';			varstruct.Attribute(n+3).Value = z_units;
	if (~isempty(add_off))
		varstruct.Attribute(end+1).Name = 'add_offset';
		varstruct.Attribute(n+4).Value = add_off;
		if (isempty(scale_factor))						% We don't use the scale_factor yet in Mirone but CF
			varstruct.Attribute(end+1).Name = 'scale_factor';	% recommends that one of scale_factor or add_offset
			varstruct.Attribute(n+5).Value = 1;				% is present the other is set as well (with neutral value)
		end
	end
	% ------------------------------------------------------------------------------------

	if (deflation_level)
		varstruct.Deflate = deflation_level;
	end

	nc_funs('addvar', fname, varstruct )

	% ------------------------------ Globals ---------------------------------
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
	if (is3D && ~is_unlimited)				% For LIMITED files, save the vector of times
		nc_funs('varput', fname, levelName, levelVec );
	end

	nc_funs('attput', fname, x_var, 'long_name', x_var );
	nc_funs('attput', fname, x_var, 'units', x_units);
	nc_funs('attput', fname, x_var, 'actual_range', [X(1) X(end)] );

	nc_funs('attput', fname, y_var, 'long_name', y_var );
	nc_funs('attput', fname, y_var, 'units', y_units);
	nc_funs('attput', fname, y_var, 'actual_range', [Y(1) Y(end)] );
	
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

	if (is3D && is_unlimited),	nc_funs('varput', fname, levelName, levelVec(1), 0, 1 );	end		% The UNLIMITED var value
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

	dims = s.Dataset(z_id).Dimension;		% Variable names of dimensions z variable - ORDER IS CRUCIAL

	% --------------------- Fish in the attribs of the Z var --------------------
	if (~isempty(s.Dataset(z_id).Attribute))
		attribNames = {s.Dataset(z_id).Attribute.Name};
	else
		attribNames = [];
	end
	dataNames = {s.Dataset.Name};			% example: {'x' 'y' 'z'}; Order may change
	ind = strcmp(attribNames,'actual_range');		z_actual_range = [];
	if (any(ind)),	z_actual_range = s.Dataset(z_id).Attribute(ind).Value;	end
	% GDAL and ESRI put the eventual georeferencing on the layer
	srsWKT = [];		strPROJ4 = [];
	ind = strcmp(attribNames,'grid_mapping');
	if (any(ind))
		grid_mapping = s.Dataset(z_id).Attribute(ind).Value;
		ind = strcmp(dataNames, grid_mapping);
		containor_id = find(ind);
		containerNames = {s.Dataset(containor_id).Attribute.Name};
		ind = strcmp(containerNames,'spatial_ref');
		if (any(ind)),	srsWKT = s.Dataset(containor_id).Attribute(ind).Value;		end
		
		ind = strcmp(containerNames,'strPROJ4');
		if (any(ind)),	strPROJ4 = s.Dataset(containor_id).Attribute(ind).Value;	end
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
	if ( ~isempty(x_id) && ~isempty(s.Attribute) )
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
	if ( ~isempty(y_id) && ~isempty(s.Attribute) )
		attribNames = {s.Dataset(y_id).Attribute.Name};
		ind = strcmp(attribNames,'actual_range');
		if (any(ind)),	y_actual_range = s.Dataset(y_id).Attribute(ind).Value;	end
	end
	% ----------------------------------------------------------------------------
	
	% --------------------- Fish the Global attributes ---------------------------
	if (~isempty(s.Attribute)),		attribNames = {s.Attribute.Name};
	else							attribNames = [];
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
			try
				Z = nc_funs('varget', fname, s.Dataset(z_id).Name);
			catch
				if ( ~isempty(strfind(lasterr, 'A memory allocation request failed')) )
					h = warndlg('We got an Out of memory error. Trying again with a MUCH SLOWER method.','WARNING');
					try
						Z = nc_funs('varget_t', fname, s.Dataset(z_id).Name);
					catch
						if (ishandle(h)),	delete(h),	end
						error(lasterr)		% Error out so it can be catch by other try-catch
					end
					if (ishandle(h)),	delete(h),	end
				end
			end
		else
			if ( ~all([s.Dimension.Length]) )
				warndlg('One ore more dimensions has ZERO length. Expect #&%%&&$.','Warning')
			end
			Z = nc_funs('varget', fname, s.Dataset(z_id).Name, zeros(1,nD), [ones(1,nD-2) s.Dataset(z_id).Size(end-1:end)]);
		end
		if (nD == 3 && z_dim(1) == 1)
			z_dim = z_dim(2:3);	% Remove singleton
		end
	else
		misc.z_id = z_id;		% Return the z_id so that we don't have to repeat the fishing process
	end
	
	misc.z_dim = z_dim;			% Use this if calling code needs to know number of layers
	nx = z_dim(end);		ny = z_dim(end-1);

	if (~isempty(x_id)),	X = double(nc_funs('varget', fname, s.Dataset(x_id).Name));
	else					X = 1:nx;
	end
	if (~isempty(y_id)),	Y = double(nc_funs('varget', fname, s.Dataset(y_id).Name));
	else					Y = 1:ny;
	end

	if (Y(2) < Y(1))
		Y = Y(end:-1:1);		Z = flipud(Z);
		if (~isempty(y_actual_range)),		y_actual_range = y_actual_range(2:-1:1);	end
	end

	if (~isempty(x_actual_range)),		head(1:2) = x_actual_range;
	else								head(1:2) = [X(1) X(end)];
	end
	if (~isempty(y_actual_range)),		head(3:4) = y_actual_range;
	else								head(3:4) = [Y(1) Y(end)];
	end

	[X, Y, Z, head] = deal_exceptions(Z, X, Y, head, s, attribNames);	% Currently, deal with Ifremer hosted SST stupidities (no coords and Kelvins)

	if (get_Z && isempty(z_actual_range))
		if ( isa(Z, 'single') )			% min/max are bugged when NaNs in singles
			zz = grdutils(Z,'-L');  z_actual_range = [zz(1) zz(2)];
		else
			z_actual_range = double([min(Z(:)) max(Z(:))]);
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
	x_range_id = strcmp(dataNames, 'x_range');
	y_range_id = strcmp(dataNames, 'y_range');
	z_range_id = strcmp(dataNames, 'z_range');
	spacing_id = strcmp(dataNames, 'spacing');
	dimension_id = strcmp(dataNames, 'dimension');
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


% _________________________________________________________________________________________________	
% -*-*-*-*-*-*-$-$-$-$-$-$-#-#-#-#-#-#-%-%-%-%-%-%-@-@-@-@-@-@-(-)-(-)-(-)-&-&-&-&-&-&-{-}-{-}-{-}-
function [X, Y, Z, head] = deal_exceptions(Z, X, Y, head, s, attribNames)
% For one reason or another some files need special treatment.
% This function deals with the known (to me) cases

	if (isempty(attribNames)),		return,		end

	% SST files hosted at ftp.ifremer.fr/pub/ifremer/cersat/SAFOSI/Products/NARSST/
	% The guys claim that they are COARDS compliant, which is a shameless lie.
	ind = strcmp(attribNames,'producer_agency');		agency = [];
	if (any(ind)),	agency = s.Attribute(ind).Value;	end
	if (strncmp(agency, 'O&SI SAF', 8))
		cvlib_mex('addS', Z, -273.15)		% I want temperatures in Celsius
		
		dirs = getappdata(0,'MIRONE_DIRS');	% Not called through a Mirone session
		if (isempty(dirs)),		return,		end
		fid = fopen([dirs.home_dir filesep 'data' filesep 'sst_by_ifremer.txt'],'r');
		if (fid < 0),			return,		end		% Parameter file does not exist

	    todos = fread(fid,'*char');		fclose(fid);
		txt = strread(todos,'%s','delimiter','\n');
		fname_coords = [];		nGrids = 0;
		for (k = 1:numel(txt))
			if (isempty(txt{k}) || txt{k}(1) == '#'),	continue,	end		% Jump comment lines
			if (strncmpi(txt{k}, 'sst_by_ifremer_coords', 21))
				[t, r] = strtok(txt{k});
				nGrids = nGrids + 1;
				fname_coords{nGrids} = ddewhite(r);
			end
		end

		if (isempty(fname_coords)),		return,		end		% Parameter file has not the necessary info

		lat = nc_funs('varget', fname_coords{1}, 'latitude');			% Load the LAT array
		lon = nc_funs('varget', fname_coords{1}, 'longitude');			% Load the LON array
		ind = strcmp(attribNames,'north_west_latitude');		NW_lat_check = double(s.Attribute(ind).Value);
		ind = strcmp(attribNames,'south_east_longitude');		SE_lon_check = double(s.Attribute(ind).Value);

		NW_lat = double(lat(end,1));	NE_lat = double(lat(end,end));		SW_lat = double(lat(1));		SE_lat = double(lat(1,end));
		NW_lon = double(lon(end,1));	NE_lon = double(lon(end,end));		SW_lon = double(lon(1));		SE_lon = double(lon(1,end));
		
		% Confirm that the given LAT/LON arrays are likely to be the right one for this file
		k = 2;
		while ( abs(NW_lat - NW_lat_check) > 0.1 || abs(SE_lon - SE_lon_check) > 0.1 && (k <= nGrids) )	% If TRUE try next file on list
			lat = nc_funs('varget', fname_coords{k}, 'latitude');			% Try with this coords file
			lon = nc_funs('varget', fname_coords{k}, 'longitude');
			NW_lat = double(lat(end,1));	NE_lat = double(lat(end,end));		SW_lat = double(lat(1));		SE_lat = double(lat(1,end));
			NW_lon = double(lon(end,1));	NE_lon = double(lon(end,end));		SW_lon = double(lon(1));		SE_lon = double(lon(1,end));
			k = k + 1;
		end
		if (k > nGrids + 1)
			warndlg('WARNING: Very likely NONE of the provided coordinate files are correct for this file.','Warning')
		end
	
		% Find out if we are using Matlab or compiled code
		try
			which('mirone');
			opt_e = '';
		catch
			opt_e = '-e';
		end

		opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', min(NW_lon,SW_lon), max(NE_lon,SE_lon), min(SW_lat,SE_lat), max(NW_lat,NE_lat));
		lon = single(lon);		lat = (single(lat));
		[Z, head] = nearneighbor_m(lon(:), lat(:), Z(:), opt_R, opt_e, '-N2', '-I0.02', '-S0.06');
		if ( any( isnan(head(5:6)) ) )
			zz = grdutils(Z,'-L');  head(5:6) = zz(1:2);
		end
		X = linspace(head(1), head(2), size(Z,2));
		Y = linspace(head(3), head(4), size(Z,1));
	end

% -----------------------------------------------------------------------------------------------------
function deflation_level = get_deflation(handles)
% Try to fish the current deflation level from several sources

	deflation_level = 0;		% We cannot trust that handles.deflation_level always exists
	if (isfield(handles, 'deflation_level'))
		deflation_level = handles.deflation_level;
	elseif (isfield(handles, 'path_data'))
		try
			prf = load([handles.path_data 'mirone_pref.mat']);
			if (isfield(prf, 'deflation_level')),	deflation_level = prf.deflation_level;	end
		catch
			mir_dirs = getappdata(0,'MIRONE_DIRS');
			if (~isempty(mir_dirs))
				home_dir = mir_dirs.home_dir;
				path_data = [home_dir filesep 'data' filesep];
				try
					prf = load([path_data 'mirone_pref.mat']);
					if (isfield(prf, 'deflation_level')),	deflation_level = prf.deflation_level;	end
				end
			end
		end
	end

function  varargout = nc_funs(opt,varargin)
% This contains some of the MEXNC (John Evans) functions packed in a single file
% It is not the whole package. Just the functions we use in Mirone
% Data is returned on its native type. That is, it is not compliant with the doubles tirany
% Except in the case of - ints and scale_factor ~= 1 OR have _FillValue and array data have _FillValue
% In that case we return the array variable as a single
% Its is also striped of the Java shit
% Joaquim Luis

	switch opt
		case 'getdiminfo'
			varargout{1} = nc_getdiminfo(varargin{:});
		case 'varget'
			varargout{1} = nc_varget(varargin{:});
		case 'varget_t'
			varargout{1} = nc_varget_t(varargin{:});
		case 'attget'
			varargout{1} = nc_attget(varargin{:});
		case 'info'
			varargout{1} = nc_info(varargin{:});
		case 'add_dimension'
			nc_add_dimension(varargin{:})
		case 'addvar'
			nc_addvar(varargin{:})
		case 'attput'
			nc_attput(varargin{:})
		case 'addhist'
			nc_addhist(varargin{:})
		case 'varput'
			nc_varput(varargin{:})
		case 'create_empty'
			nc_create_empty(varargin{:})
		case 'cat'
			nc_cat(varargin{:})
		case 'dump'
			if (nargout),		varargout{1} = nc_dump(varargin{:});
			else				nc_dump(varargin{:});
			end
	end

% --------------------------------------------------------------------
function fileinfo = nc_info(ncfile)
% This function is the core function of the snctools/nc_info

fileinfo.Filename = ncfile;

[ncid, status] = mexnc('open', ncfile, nc_nowrite_mode);
if status ~= 0
	if (status == -101)			% Damn dirty trick to avoid a situation where a file (from Ocean Color only?)
		fileinfo.Dataset = [];	% was previously open by gdalread and something in it didn't copletely close.
		return					% As a consequence accessing it with mecnc fails. But that's what we want anyway.
	end							% so return right now and aux_funs(findFileType,...) will send it to gdal again.
    snc_error ( 'NC_INFO:MEXNC:OPEN', mexnc('strerror', status) );
end

[ndims, nvars, ngatts, record_dimension, status] = mexnc('INQ', ncid);
if status ~= 0
    mexnc('close',ncid);
    snc_error('NC_FUNS:NC_INFO:MEXNC:INQ', mexnc('strerror', status));
end

% Get the dimensions
if ndims == 0
	Dimension = struct ( [] );
else
	if ndims > 0
		Dimension(1) = nc_getdiminfo( ncid, 0 );
	end
	Dimension = repmat( Dimension, ndims,1 );
	for dimid = 1:ndims-1
		Dimension(dimid+1) = nc_getdiminfo( ncid, dimid );
	end
end

% Get the global attributes.
if ngatts == 0
	fileinfo.Attribute = struct([]);
else
	if ngatts > 0
		Attribute(1) = nc_get_attribute_struct ( ncid, nc_global, 0 );
	end
	Attribute = repmat ( Attribute, ngatts, 1 );
	for attnum = 1:ngatts-1
		Attribute(attnum+1) = nc_get_attribute_struct ( ncid, nc_global, attnum );
	end
	fileinfo.Attribute = Attribute;
end

% Get the variable information.
if nvars == 0
	Dataset = struct([]);
else
	if ( nvars > 0 )
		Dataset(1) = nc_getvarinfo( ncid, 0 );
	end
	Dataset = repmat ( Dataset, nvars, 1 );
	for varid=1:nvars-1
		Dataset(varid+1) = nc_getvarinfo( ncid, varid );
	end
end

fileinfo.Dimension = Dimension;
fileinfo.Dataset = Dataset;
mexnc('close',ncid);

% --------------------------------------------------------------------
function dinfo = nc_getdiminfo ( arg1, arg2 )
% NC_GETDIMINFO:  returns metadata about a specific NetCDF dimension
%
% DINFO = NC_GETDIMINFO(NCFILE,DIMNAME) returns information about the
% dimension DIMNAME in the netCDF file NCFILE.
%
% DINFO = NC_GETDIMINFO(NCID,DIMID) returns information about the
% dimension with numeric Id DIMID in the already-opened netCDF file
% with file Id NCID.  This form is not recommended for use from the
% command line.
%
% Upon output, DINFO will have the following fields.
%
%    Name:  
%        a string containing the name of the dimension.
%    Length:  
%        a scalar equal to the length of the dimension
%    Unlimited:  
%        A flag, either 1 if the dimension is an unlimited dimension
%        or 0 if not.
%
% In case of an error, an exception is thrown.
%

snc_nargchk(2,2,nargin);
snc_nargoutchk(1,1,nargout);

% If we are here, then we must have been given something local.
if ischar(arg1) && ischar(arg2)
    dinfo = handle_char_nc_getdiminfo(arg1,arg2);
elseif isnumeric( arg1 ) && isnumeric( arg2 )
	dinfo = handle_numeric_nc_getdiminfo(arg1,arg2);
else
	snc_error( 'NC_FUNS:NC_GETDIMINFO:badInput', ...
	            'Must supply either two character or two numeric arguments.' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dinfo = handle_char_nc_getdiminfo ( ncfile, dimname )

[ncid,status ]=mexnc('open', ncfile, nc_nowrite_mode );
if status ~= 0
	snc_error ( 'NC_FUNS:NC_GETDIMINFO:handle_char_nc_getdiminfo:openFailed', mexnc('strerror', status) );
end

[dimid, status] = mexnc('INQ_DIMID', ncid, dimname);
if ( status ~= 0 )
	mexnc('close',ncid);
	snc_error ( 'NC_FUNS:NC_GETDIMINFO:handle_char_nc_getdiminfo:inq_dimidFailed', mexnc('strerror', status) );
end

dinfo = handle_numeric_nc_getdiminfo ( ncid,  dimid );
mexnc('close',ncid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dinfo = handle_numeric_nc_getdiminfo ( ncid, dimid )

[unlimdim, status] = mexnc ( 'inq_unlimdim', ncid );
if status ~= 0
	mexnc('close',ncid);
	snc_error ( 'NC_FUNS:NC_GETDIMINFO:MEXNC:inq_ulimdimFailed', mexnc ( 'strerror', status ) );
end

[dimname, dimlength, status] = mexnc('INQ_DIM', ncid, dimid);
if status ~= 0
	mexnc('close',ncid);
	snc_error ( 'NC_FUNS:NC_GETDIMINFO:MEXNC:inq_dimFailed', mexnc ( 'strerror', status ) );
end

dinfo.Name = dimname;
dinfo.Length = dimlength;

if (dimid == unlimdim),		dinfo.Unlimited = true;
else						dinfo.Unlimited = false;
end

% --------------------------------------------------------------------
function Dataset = nc_getvarinfo ( arg1, arg2 )
% NC_GETVARINFO:  returns metadata about a specific NetCDF variable
%
% VINFO = NC_GETVARINFO(NCFILE,VARNAME) returns a metadata structure VINFO about
% the variable VARNAME in the netCDF file NCFILE.
%
% VINFO = NC_GETVARINFO(NCID,VARID) returns a metadata structure VINFO about
% the variable whose netCDF variable-id is VARID, and whose parent file-id is 
% NCID.  The netCDF file is assumed to be open, and in this case the file will
% not be closed upon completion.
%
% VINFO will have the following fields:
%
%    Name:  
%       a string containing the name of the variable.
%    Nctype:  
%       a string specifying the NetCDF datatype of this variable.
%    Unlimited:  
%       Flag, either 1 if the variable has an unlimited dimension or 0 if not.
%    Dimensions:  
%       a cell array with the names of the dimensions upon which this variable 
%       depends.
%    Attribute:  
%       An array of structures corresponding to the attributes defined for the 
%       specified variable.
%                         
%    Each "Attribute" element contains the following fields.
%
%       Name:  
%           a string containing the name of the attribute.
%       Nctype:  
%           a string specifying the NetCDF datatype of this attribute.
%       Attnum:  
%           a scalar specifying the attribute id
%       Value: 
%           either a string or a double precision value corresponding to the 
%           value of the attribute
%
% In case of an error, an exception is thrown.
%

% Show usage if too few arguments.
snc_nargchk(2,2,nargin);
snc_nargoutchk(1,1,nargout);

% If we are here, then we must have been given something local.
if ischar(arg1) && ischar(arg2)

	ncfile = arg1;
	varname = arg2;

	[ncid,status ]=mexnc('open',ncfile,nc_nowrite_mode);
	if status ~= 0
	    snc_error( 'NC_FUNS:NC_VARGET:MEXNC:OPEN', mexnc('strerror', status) );
	end

	[varid, status] = mexnc('INQ_VARID', ncid, varname);
	if ( status ~= 0 )
	    snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:INQ_VARID', mexnc('strerror', status) );
	end
	
	Dataset = get_varinfo ( ncid,  varid );

	% close whether or not we were successful.
	mexnc('close',ncid);

elseif isnumeric ( arg1 ) && isnumeric ( arg2 )
	ncid = arg1;
	varid = arg2;
	Dataset = get_varinfo ( ncid,  varid );

else
	snc_error ( 'NC_FUNS:NC_GETVARINFO:badTypes', 'Must have either both character inputs, or both numeric.' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dataset = get_varinfo ( ncid, varid )

[record_dimension, status] = mexnc ( 'INQ_UNLIMDIM', ncid );
if status ~= 0
    mexnc('close',ncid);
    snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:INQ_UNLIMDIM', mexnc('strerror', status) );
end

[varname, datatype, ndims, dims, natts, status] = mexnc('INQ_VAR', ncid, varid);
if status ~= 0 
    mexnc('close',ncid);
    snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:INQ_VAR', mexnc('strerror', status) );
end

Dataset.Name = varname;
Dataset.Nctype = datatype;

% Assume the current variable does not have an unlimited dimension until we know that it does.
Dataset.Unlimited = false;

if ndims == 0
	Dataset.Dimension = {};
	Dataset.Size = 1;
else

	for j = 1:ndims
		[dimname, dimlength, status] = mexnc('INQ_DIM', ncid, dims(j));
		if ( status ~= 0 )
		    mexnc('close',ncid);
		    snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:INQ_DIM', mexnc('strerror', status) );
		end
	
		Dataset.Dimension{j} = dimname; 
		Dataset.Size(j) = dimlength;
		if (dims(j) == record_dimension),		Dataset.Unlimited = true;		end
	end
end

% get all the attributes
if natts == 0
	Dataset.Attribute = struct([]);
else
	for attnum = 0:natts-1
		Dataset.Attribute(attnum+1) = nc_get_attribute_struct ( ncid, varid, attnum );
	end
end

% --------------------------------------------------------------------------------------------
function tf = nc_isunlimitedvar(ncfile,varname)
%NC_ISUNLIMITEDVAR determine if variable has unlimited dimension.
%
%   TF = NC_ISUNLIMITEDVAR(NCFILE,VARNAME) returns true if the netCDF
%   variable VARNAME in the netCDF file NCFILE has an unlimited dimension,
%   and false otherwise.
%
%   Example:
%       nc_dump('example.nc');
%       tf = nc_isunlimitedvar('example.nc','time_series')
%
%   See also NC_ISCOORDVAR, NC_DUMP.

	try
		info = nc_getvarinfo(ncfile, varname);
	catch 
		e = lasterror;
		switch ( e.identifier )
			case 'NC_FUNS:NC_VARGET:MEXNC:INQ_VAR'
				tf = false;
				return
			otherwise
				error(e);
		end
	end

	tf = info.Unlimited;

% --------------------------------------------------------------------------------------------
function values = nc_varget_t(ncfile, varname, varargin)
% Same as nc_varget() but reads already transposed (10 times SLOWER)

	parse_and_validate_args_varget(ncfile,varname,varargin{:});

	[ncid, status] = mexnc('open',ncfile,'NOWRITE');
	if (status ~= 0)
		snc_error('NC_FUNS:NC_VARGET_T:MEXNC:OPEN', mexnc('strerror', status) );
	end

	[varid, status] = mexnc('inq_varid', ncid, varname);
	if (status ~= 0)
		mexnc('close',ncid);
		snc_error ( 'NC_FUNS:NC_VARGET_T:MEXNC:INQ_VARID', mexnc('strerror', status) );
	end

	[dud, var_type, nvdims, dimids, dud, status] = mexnc('inq_var', ncid, varid);
	if (status ~= 0)
		mexnc('close',ncid);
		snc_error ( 'NC_FUNS:NC_VARGET_T:MEXNC:INQ_VAR', mexnc('strerror',status) );
	end

	the_var_size = determine_varsize_mex( ncid, dimids, nvdims );
    [start, count, stride] = snc_get_indexing(nvdims, the_var_size, varargin{:});

	% Check that START index is ok.
	if ~isempty(start) && any(start > the_var_size)
		snc_error('NC_FUNS:NC_VARGET:badStartIndex', 'The START index argument exceeds the size of the variable.');
	end

	% If the user had set non-positive numbers in "count", then we replace them
	% with what we need to get the rest of the variable.
	negs = find(count < 0);
	count(negs) = the_var_size(negs) - start(negs);

	% Check that the start, count, stride parameters have appropriate lengths.
	% Otherwise we get confusing error messages later on.
	check_index_vectors(start,count,stride,nvdims,ncid,varname);

	% What mexnc operation will we use?
	if (numel(count) > 1)
		prefix = 'get_varm';
		switch (var_type)
			case nc_int,		funcstr = [prefix '_int'];
			case nc_float,		funcstr = [prefix '_float'];
			case nc_double,		funcstr = [prefix '_double'];
			case nc_short,		funcstr = [prefix '_short'];
			case nc_char,		funcstr = [prefix '_text'];
			case nc_byte,		funcstr = [prefix '_schar'];
			otherwise
				snc_error ( 'NC_FUNS:NC_VARGET_T:badDatatype', sprintf('Unhandled datatype %d.', var_type) );
		end
		funcstr = funcstr_as_float ( ncid, varid, var_type, funcstr);	% If we'll need to convert to float, better read as such

		imap_coord = [1 count(1)];
		[values, status] = mexnc(funcstr, ncid, varid, start, count, stride, imap_coord);

	else
		[funcstr_family, funcstr] = determine_funcstr(var_type, nvdims, start, count, stride);
		funcstr = funcstr_as_float ( ncid, varid, var_type, funcstr);	% If we'll need to convert to float, better read as such
		switch funcstr_family
			case 'get_var',		[values, status] = mexnc( funcstr, ncid, varid );
			case 'get_var1',	[values, status] = mexnc( funcstr, ncid, varid, 0 );
			case 'get_vara',	[values, status] = mexnc( funcstr, ncid, varid, start, count );
			case 'get_vars',	[values, status] = mexnc( funcstr, ncid, varid, start, count, stride );
			otherwise
				snc_error('NC_FUNS:NC_VARGET:unhandledType', sprintf ('Unhandled function string type ''%s''\n', funcstr_family) );
		end
	end

	if ( status ~= 0 )
		mexnc('close',ncid);
		snc_error ( 'NC_FUNS:NC_VARGET_T', mexnc('strerror', status) );
	end

	% Test for situations like the Ifremer SST files where when we get here we have
	% the_var_size = [1 1024 1024]		and values = [1024 x 1024]
	% That is, the variable on the netCDF file is singleton on the leading dimension
	if ( (the_var_size(1) == 1) && (numel(the_var_size) > 2) && (ndims(values) == (numel(the_var_size) - 1)) )
		the_var_size = the_var_size(2:end);
	end

	% If it's a 1D vector, make it a column vector.  Otherwise permute the data
	% to make up for the row-major-order-vs-column-major-order issue.
	if (numel(the_var_size) == 1),		values = values(:);		end

	[values, status] = handle_fill_value_mex ( ncid, varid, var_type, values );	% Do both '_FillValue' & 'missing_value'
	if (status)		% '_FillValue' was not found. Try the 'missing_value'
		values = handle_mex_missing_value ( ncid, varid, var_type, values );
	end

	values = handle_scaling_mex(ncid, varid, values);
	values = squeeze( values );		% remove any singleton dimensions.
	mexnc('close',ncid);

% --------------------------------------------------------------------------------------------
function values = nc_varget(ncfile, varname, varargin )
% NC_VARGET:  Retrieve data from a netCDF variable.
%
% DATA = NC_VARGET(NCFILE,VARNAME) retrieves all the data from the 
% variable VARNAME in the netCDF file NCFILE.
%
% DATA = NC_VARGET(NCFILE,VARNAME,START,COUNT) retrieves the contiguous
% portion of the variable specified by the index vectors START and 
% COUNT.  Remember that NC_FUNS indexing is zero-based, not 
% one-based.  Specifying a -1 in COUNT means to retrieve everything 
% along that dimension from the START coordinate.
%
% DATA = NC_VARGET(NCFILE,VARNAME,START,COUNT,STRIDE) retrieves 
% a non-contiguous portion of the dataset.  The amount of
% skipping along each dimension is given through the STRIDE vector.
%
% NCFILE can also be an OPeNDAP URL if the proper NC_FUNS backend is
% installed.  See the README for details.
% 
% NC_VARGET tries to be intelligent about retrieving the data.
% Since most general matlab operations are done in double precision,
% retrieved numeric data will be cast to double precision, while 
% character data remains just character data.  
%
% Singleton dimensions are removed from the output data.  
%
% A '_FillValue' attribute is honored by flagging those datums as NaN.
% A 'missing_value' attribute is honored by flagging those datums as 
% NaN.  The exception to this is for NC_CHAR variables, as mixing 
% character data and NaN doesn't really seem to work in matlab.
%
% If the named NetCDF variable has valid scale_factor and add_offset 
% attributes, then the data is scaled accordingly.  
%
% EXAMPLE:
% #1.  In this case, the variable in question has rank 2, and has size 
%      500x700.  We want to retrieve starting at row 300, column 250.
%      We want 100 contiguous rows, 200 contiguous columns.
% 
%      vardata = nc_varget ( file, variable_name, [300 250], [100 200] );

% 	values = nc_varget_t(ncfile, varname, varargin{:} );
% 	return

	snc_nargchk(2,5,nargin);
	snc_nargchk(1,1,nargout);

	[start, count, stride] = parse_and_validate_args_varget(ncfile,varname,varargin{:});

	[ncid,status] = mexnc('open',ncfile,'NOWRITE');
	if status ~= 0
		snc_error('NC_FUNS:NC_VARGET:MEXNC:OPEN', mexnc('strerror', status) );
	end

	[varid, status] = mexnc('inq_varid', ncid, varname);
	if (status ~= 0)
		mexnc('close',ncid);
		snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:INQ_VARID', mexnc('strerror', status) );
	end

	[dud,var_type,nvdims,dimids,dud,status] = mexnc('inq_var', ncid, varid);
	if (status ~= 0)
		mexnc('close',ncid);
		snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:INQ_VAR', mexnc('strerror',status) );
	end

	the_var_size = determine_varsize_mex( ncid, dimids, nvdims );
    %[start, count, stride] = snc_get_indexing(nvdims, the_var_size, varargin{:});	% Only would worth if fvd == true

	% Check that the start, count, stride parameters have appropriate lengths.
	% Otherwise we get confusing error messages later on.
	check_index_vectors(start,count,stride,nvdims,ncid,varname);

	% What mexnc operation will we use?
	[funcstr_family, funcstr] = determine_funcstr( var_type, nvdims, start, count, stride );

	% Check that START index is ok.
	if ~isempty(start) && any(start > the_var_size)
		snc_error('NC_FUNS:NC_VARGET:badStartIndex', 'The START index argument exceeds the size of the variable.');
	end
	
	funcstr = funcstr_as_float ( ncid, varid, var_type, funcstr);	% If we'll need to convert to float, better read as such

	% If the user had set non-positive numbers in "count", then we replace them
	% with what we need to get the rest of the variable.
	negs = find(count < 0);
	count(negs) = the_var_size(negs) - start(negs);

	% At long last, retrieve the data.
	switch funcstr_family
		case 'get_var',		[values, status] = mexnc( funcstr, ncid, varid );
		case 'get_var1',	[values, status] = mexnc( funcstr, ncid, varid, 0 );
		case 'get_vara',	[values, status] = mexnc( funcstr, ncid, varid, start, count );
		case 'get_vars',	[values, status] = mexnc( funcstr, ncid, varid, start, count, stride );
		otherwise
			snc_error('NC_FUNS:NC_VARGET:unhandledType', sprintf ('Unhandled function string type ''%s''\n', funcstr_family) );
	end

	if ( status ~= 0 )
		mexnc('close',ncid);
		snc_error ( 'NC_FUNS:NC_VARGET:%s', mexnc('strerror', status) );
	end

	% Test for situations like the Ifremer SST files where when we get here we have
	% the_var_size = [1 1024 1024]		and values = [1024 x 1024]
	% That is, the variable on the netCDF file is singleton on the leading dimension
	if ( (the_var_size(1) == 1) && (numel(the_var_size) > 2) && (ndims(values) == (numel(the_var_size) - 1)) )
		the_var_size = the_var_size(2:end);
	end

	% If it's a 1D vector, make it a column vector.  Otherwise permute the data
	% to make up for the row-major-order-vs-column-major-order issue.
	if (length(the_var_size) == 1)
		values = values(:);
	else
        % Ok it's not a 1D vector.  If we are not preserving the fastest
        % varying dimension, we should permute the data.
		pv = numel(the_var_size):-1:1;
		values = permute(values, pv);
	end                                                                                   

	[values, status] = handle_fill_value_mex ( ncid, varid, var_type, values );	% Do both '_FillValue' & 'missing_value'
	if (status)		% '_FillValue' was not found. Try the 'missing_value'
		values = handle_mex_missing_value ( ncid, varid, var_type, values );
	end
	values = handle_scaling_mex ( ncid, varid, values );

	% remove any singleton dimensions.
	values = squeeze( values );

	mexnc('close',ncid);

% --------------------------------------------------------------------------------
function [start, count, stride] = parse_and_validate_args_varget(ncfile,varname,varargin)
%
% Set up default outputs.
start = [];		count = [];		stride = [];

	switch nargin
		case 4
			start = varargin{1};
			count = varargin{2};
		case 5
			start = varargin{1};
			count = varargin{2};
			stride = varargin{3};
	end

	% Error checking on the inputs.
	if ~ischar(ncfile)
		snc_error ( 'NC_FUNS:NC_VARGET:badInput', 'the filename must be character.' );
	end
	if ~ischar(varname)
		snc_error ( 'NC_FUNS:NC_VARGET:badInput', 'the variable name must be character.' );
	end

	if ~isnumeric ( start )
		snc_error ( 'NC_FUNS:NC_VARGET:badInput', 'the ''start'' argument must be numeric.' );
	end
	if ~isnumeric ( count )
		snc_error ( 'NC_FUNS:NC_VARGET:badInput', 'the ''count'' argument must be numeric.' );
	end
	if ~isnumeric ( stride )
		snc_error ( 'NC_FUNS:NC_VARGET:badInput', 'the ''stride'' argument must be numeric.' );
	end
	
% --------------------------------------------------------------------------------
function [start, count, stride] = snc_get_indexing(nvdims, var_size, varargin)
% Common private function for setting up indexing for NC_VARGET.

	if (nvdims == 0)	% This will happen in the case of a singleton.  No need to go further.
		start = 0;		count = 1;		stride = 1;
		return
	end

	switch(numel(varargin))
		case 0			% retrieve everything.
			start  = zeros(1,nvdims);
			count  = var_size;
			stride = ones(1,nvdims);

		case 1			% if only start was provided, then the count is implied to be one.
			start  = varargin{1};
			count  = ones(1,nvdims);
			stride = ones(1,nvdims);

		case 2			% just a contiguous hyperslab.
			start  = varargin{1};
			count  = varargin{2};
			stride = ones(1,nvdims);

		case 3
			start  = varargin{1};
			count  = varargin{2};
			stride = varargin{3};
	end

	start = double(start(:)');
	count = double(count(:)');
	stride = double(stride(:)');

	% If the user had set non-positive numbers in "count", then we replace them
	% with what we need to get the rest of the variable.
	negs = find((count<0) | isinf(count));
	count(negs) = (var_size(negs) - start(negs)) ./ stride(negs);

	% Ok, now do some final validation.
	if (any(start < 0))
		snc_error('NC_FUNS:INDEXING:badStartIndex', 'The START argument should be nonnegative.');
	end

	if (any(count <= 0))
		snc_error('NC_FUNS:INDEXING:badStartIndex', 'The COUNT argument should be positive.');
	end

	if (any(stride <= 0))
		snc_error('NC_FUNS:INDEXING:badStartIndex', 'The STRIDE argument should be positive.');
	end

	if (~isnumeric(start) || ~isnumeric(count) || ~isnumeric(stride))
		snc_error('NC_FUNS:INDEXING:badIndexType', 'Any index arguments should be numeric');
	end

	if (numel(start) ~= numel(count)) || (numel(count) ~= numel(stride)) || (numel(stride) ~= nvdims)
		snc_error('NC_FUNS:INDEXING:badIndexLength', sprintf('The lengths of the index arguments should be %d.', nvdims));
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prefix,funcstr] = determine_funcstr ( var_type, nvdims, start, count, stride )
% DETERMINE_FUNCSTR
%     Determines if we are to use, say, 'get_var1_text', or 'get_vars_double', or whatever.

	% Determine if we are retriving a single value, the whole variable, a 
	% contiguous portion, or a strided portion.
	if (nvdims == 0)													% It is a singleton variable.
		prefix = 'get_var1';
	elseif isempty(start)												% retrieving the entire variable.
		prefix = 'get_var';
	elseif (~isempty(start) && ~isempty(count) && isempty(stride))		% retrieving a contiguous portion
		prefix = 'get_vara';
	elseif (~isempty(start) && ~isempty(count) && ~isempty(stride))		% retrieving a contiguous portion
		prefix = 'get_vars';
	else
		snc_error ( 'NC_FUNS:NC_VARGET:FUNCSTR', 'Could not determine funcstr prefix.' );
	end

	switch (var_type)
		case nc_int,		funcstr = [prefix '_int'];
		case nc_float,		funcstr = [prefix '_float'];
		case nc_double,		funcstr = [prefix '_double'];
		case nc_short,		funcstr = [prefix '_short'];
		case nc_char,		funcstr = [prefix '_text'];
		case nc_byte,		funcstr = [prefix '_schar'];
		case nc_ubyte,		funcstr = [prefix '_uchar'];
		otherwise
			snc_error ('NC_FUNS:NC_VARGET:badDatatype', sprintf('Unhandled datatype %d.', var_type));
	end
	
% ------------------------------------------------------------------------------
function funcstr = funcstr_as_float ( ncid, varid, var_type, funcstr)
% Check for the conditions upon which a short int variable (the grid array) will
% later be converted to single. If any of those conditions are met, change the
% reading function string FUNCSTR to its 'float' type which will force the reading
% into a single variable in first place.
%
% Currently tested conditions are that a SCALE_FACTOR ~= 1 or that the 
% '_FillValue' or 'missing_value' attributes exist

	if (var_type ~= nc_short && var_type ~= nc_int)
		return
	end

	have_scale = false;
	[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'scale_factor' );
	if (status == 0),    have_scale = true;			end

	% Return early if we don't have it.
	if (~have_scale),    return,	end

	[scale_factor, status] = mexnc ( 'get_att_double', ncid, varid, 'scale_factor' );
	if ( status ~= 0 )
		mexnc('close',ncid);
		snc_error('NC_FUNS:NC_VARGET:MEXNC:GET_ATT_DOUBLE', mexnc('strerror', status) );
	end

	do_replace = false;

	if ( scale_factor ~= 1)
		do_replace = true;
	else
		[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, '_FillValue' );
		[dud, dud, status2] = mexnc('INQ_ATT', ncid, varid, 'missing_value' );
		if (status == 0 || status2 == 0)
			do_replace = true;
		end
	end

	if (do_replace)
		ind = strfind(funcstr, '_');
		funcstr = [funcstr(1:ind(end)) 'float'];
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [values, status] = handle_fill_value_mex ( ncid, varid, var_type, values )
% HANDLE_MEX_FILL_VALUE: If there is a fill value, then replace such values with NaN.
% And, since the Job is exactly the same as for the missing_value case we now deal
% both cases here (J. Luis 11-06-2009).
% However, a situation may arise where there is no '_FillValue'. In that case STATUS
% contains that ifo and the calling routine should call handle_mex_missing_value()

miss_value = [];

[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, '_FillValue' );
if ( status == 0 )
	
	[dud, dud, status2] = mexnc('INQ_ATT', ncid, varid, 'missing_value' );

    switch ( var_type )
        case {nc_char, nc_byte}
			% For now, do nothing.  Does a fill value even make sense with char data?
			% If it does, please tell me so.
        case { nc_int, nc_short}
			[fill_value, status] = mexnc( 'get_att_double', ncid, varid, '_FillValue' );
			if ( ~status2 )			% Get also the 'missing_value' and do the job for both
				miss_value = mexnc( 'get_att_double', ncid, varid, 'missing_value' );
			end
			if (~isnan(fill_value) || ~isnan(miss_value))
				ind = (values == fill_value);
				if ( ~isempty(miss_value) && (fill_value ~= miss_value) )	% If 'missing_value' is different from '_FillValue'
					ind2 = (values == miss_value);
					ind = (ind | ind2);				% Ensemble of fill_ and miss_ values
				end
				if (any(ind(:)))
					values = single(values);		% Here I give up and convert to singles
					values(ind) = NaN;
				end
			else
				% An idiotic case. A _FillValue = NAN in a array of integers.
			end
        case { nc_double, nc_float }
			[fill_value, status] = mexnc( 'get_att_double', ncid, varid, '_FillValue' );
			if ( ~status2 )			% Get also the 'missing_value' and do the job for both
				miss_value = mexnc( 'get_att_double', ncid, varid, 'missing_value' );
			end
			if (~isnan(fill_value))
				values(values == fill_value) = NaN;
			end
			if (  ~isempty(miss_value) && (~isnan(miss_value) && (miss_value ~= fill_value)) )
				values(values == miss_value) = NaN;
			end
        otherwise
			mexnc('close',ncid);
			snc_error ( 'NC_FUNS:fill_value_mex',  sprintf('unhandled datatype %d\n', var_type) );
    end

    if ( status ~= 0 )
        mexnc('close',ncid);
        snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:GET_ATT', mexnc ( 'strerror', status ) );
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function values = handle_mex_missing_value ( ncid, varid, var_type, values )
% HANDLE_MEX_MISSING_VALUE
%     If there is a missing value, then replace such values with NaN.

[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'missing_value' );
if ( status == 0 )

    switch ( var_type )
        case {nc_char, nc_byte}
			% For now, do nothing.  Does a fill value even make sense with char data?
			% If it does, please tell me so.
        case { nc_int, nc_short }
			[fill_value, status] = mexnc('get_att_double', ncid, varid, 'missing_value');
			if (~isnan(fill_value))
				ind = (values == fill_value);
				if (any(ind(:)))
					values = single(values);		% Here I give up and convert to singles
					values(ind) = NaN;
				end
			else
				% An idiotic case. A _FillValue = NAN in a array of integers.
			end
        case { nc_double, nc_float }
			[fill_value, status] = mexnc('get_att_double', ncid, varid, 'missing_value');
			if (status ~= 0)		% It can be a string !!!!!!!!!!!!!!!!!!!!
				[fill_value, status] = mexnc('get_att_text', ncid, varid, 'missing_value');
				if (status == 0 && strncmpi(fill_value, 'nan', 3))
					fill_value = NaN;
				else
					fill_value = str2double(fill_value);
				end
			end
			if (~isnan(fill_value))
				values(values == fill_value) = NaN;
			end
        otherwise
			mexnc('close',ncid);
			snc_error ( 'NC_FUNS:mex_missing_value', sprintf('unhandled datatype %d\n', mfilename, var_type) );
    end

    if ( status ~= 0 )
        mexnc('close',ncid);
        snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:GET_ATT', mexnc('strerror', status) );
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function values = handle_scaling_mex ( ncid, varid, values )
% HANDLE_MEX_SCALING
%     If there is a scale factor and/or  add_offset attribute, convert the data
%     to double precision and apply the scaling.

have_scale = false;
have_addoffset = false;
[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'scale_factor' );
if (status == 0),    have_scale = true;			end
[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'add_offset' );
if (status == 0),    have_addoffset = true;		end

% Return early if we don't have either one.
if (~(have_scale || have_addoffset)),    return,	end

if (have_scale)
    [scale_factor, status] = mexnc ( 'get_att_double', ncid, varid, 'scale_factor' );
    if ( status ~= 0 )
        mexnc('close',ncid);
        snc_error('NC_FUNS:NC_VARGET:MEXNC:GET_ATT_DOUBLE', mexnc('strerror', status) );
    end
	% If data is other than single or double, change it to single
	if ( scale_factor ~= 1 && ~( isa(values,'single') || isa(values,'double')) )
		values = single(values);
	end
	if ( scale_factor == 1 ),	have_scale = false;		end
end

if (have_addoffset)
    [add_offset, status] = mexnc ( 'get_att_double', ncid, varid, 'add_offset' );
    if ( status ~= 0 )
        mexnc('close',ncid);
        snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:GET_ATT_DOUBLE', mexnc('strerror', status) );
    end
	if ( add_offset == 0 ),		have_addoffset = false;		end
end

if (ndims(values) > 2)			% OpenCV is bugged with 3D grids and sometime we get here a singleton 3D array
	values = squeeze(values);
end

if (have_scale && have_addoffset)
	cvlib_mex('CvtScale',values, scale_factor, add_offset)
elseif (have_scale)
	cvlib_mex('CvtScale',values, scale_factor)
elseif (have_addoffset)
	if (~isa(values,'int8'))		% OpenCV doesn't add to int8 arrays
		cvlib_mex('addS',values, add_offset)
	else
		if (add_offset > 0)			% We take it to mean do [-128 127] => [0 255] conversion
			values = uint8(cvlib_mex('addS', int16(values), add_offset));
		else						% Whatever
			values = int8(cvlib_mex('addS', int16(values), add_offset));
		end
	end
end
	

% -----------------------------------------------------------------------------
function dump = nc_dump(file_name, varargin )
% NC_DUMP:  a Matlab counterpart to the NetCDF utility 'ncdump'.
%     NC_DUMP(NCFILE) prints metadata about the netCDF file NCFILE.  
%     NC_DUMP(NCFILE,VARNAME) prints metadata about just the one netCDF variable
%     named VARNAME.

snc_nargchk(1,2,nargin);
%snc_nargoutchk(0,0,nargout);

if nargin == 2
	do_restricted_variable = true;
	restricted_variable = varargin{1};
else
	do_restricted_variable = false;	
	restricted_variable = [];
end

metadata = nc_info ( file_name );

% print out name of file
if (nargout)
	dump = sprintf('netcdf %s { \n', metadata.Filename );
	dump = dump_dimension_metadata ( metadata, dump );
	dump = dump_variable_metadata ( metadata, restricted_variable, dump );
else
	fprintf(1, 'netcdf %s { \n', metadata.Filename );
	dump_dimension_metadata ( metadata );
	dump_variable_metadata ( metadata, restricted_variable );
end

	
if ( do_restricted_variable == false )
	if (nargout),	dump = dump_global_attributes ( metadata, dump );
	else			dump_global_attributes ( metadata );
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dump = dump_dimension_metadata ( metadata, dump )

if isfield( metadata, 'Dimension' )
	num_dims = length(metadata.Dimension);
else
	num_dims = 0;
end

n_argin = nargin;
if (n_argin == 2),	dump = [dump sprintf('dimensions:\n')];
else				fprintf ( 1, 'dimensions:\n' );
end
for j = 1:num_dims
	if metadata.Dimension(j).Unlimited
		str = sprintf('\t%s = UNLIMITED ; (%i currently)\n', deblank(metadata.Dimension(j).Name), metadata.Dimension(j).Length );
		if (n_argin == 2),	dump = [dump str];
		else				fprintf( 1, str );
		end
	else
		str = sprintf('\t%s = %i ;\n', metadata.Dimension(j).Name, metadata.Dimension(j).Length);
		if (n_argin == 2),	dump = [dump str];
		else				fprintf( 1, str );
		end
	end
end
if (n_argin == 2),	dump = [dump sprintf('\n')];
else				fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dump = dump_variable_metadata ( metadata, restricted_variable, dump )

if isfield ( metadata, 'Dataset' )
	num_vars = length(metadata.Dataset);
else
	num_vars = 0;
end

n_argin = nargin;
if (n_argin == 3),	dump = [dump sprintf('variables:\n')];
else				fprintf ( 1, 'variables:\n' );
end
for j = 1:num_vars
	if ~isempty(restricted_variable)
		if ~strcmp ( restricted_variable, metadata.Dataset(j).Name )
			continue
		end
	end
	if (n_argin == 3),		dump = dump_single_variable( metadata.Dataset(j), dump );
	else					dump_single_variable( metadata.Dataset(j) );
	end
end
if (n_argin == 3),	dump = [dump sprintf('\n')];
else				fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dump = dump_single_variable ( var_metadata, dump )

switch ( var_metadata.Nctype )
	case 1,		str = sprintf('\tbyte ' );
	case 2,		str = sprintf('\tchar ' );
	case 3,		str = sprintf('\tshort ' );
	case 4,		str = sprintf('\tlong ' );
	case 5,		str = sprintf('\tfloat ' );
	case 6,		str = sprintf('\tdouble ' );
end

n_argin = nargin;
str = [str sprintf('%s', var_metadata.Name)];
if (n_argin == 1),		fprintf ( 1, str );
else					dump = [dump str];
end

if isempty(var_metadata.Dimension)
	str = sprintf ('([]),');
else
	str = sprintf('(%s', var_metadata.Dimension{1} );
	for (k = 2:length(var_metadata.Size))
		str = [str sprintf(',%s', var_metadata.Dimension{k} )];
	end
	str = [str sprintf('), ')];
end

if isempty(var_metadata.Dimension)
	str = [str sprintf('shape = [1]\n')];
else
	str = [str sprintf('shape = [%d', var_metadata.Size(1) )];
	for k = 2:length(var_metadata.Size)
		str = [str sprintf(' %d', var_metadata.Size(k) )];
	end
	str = [str sprintf(']\n')];
end

if (n_argin == 2),	dump = [dump str];
else				fprintf(1, str);
end

% Now do all attributes for each variable.
num_atts = length(var_metadata.Attribute);
for k = 1:num_atts
	if (n_argin == 1)
		dump_single_attribute ( var_metadata.Attribute(k), var_metadata.Name );
	else
		dump = dump_single_attribute ( var_metadata.Attribute(k), var_metadata.Name, dump );
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dump = dump_single_attribute ( attribute, varname, dump )

switch ( attribute.Nctype )
	case 0
		att_val = '';		att_type = 'NC_NAT';
	case 1
		att_type = 'x';		att_val = sprintf ('%d ', fix(attribute.Value) );
	case 2
		att_type = '';		att_val = sprintf ('"%s" ', attribute.Value );
	case 3
		att_type = 's';		att_val = sprintf ('%i ', attribute.Value );
	case 4
		att_type = 'd';		att_val = sprintf ('%i ', attribute.Value );
	case 5
		att_type = 'f';		att_val = sprintf ('%f ', attribute.Value );
	case 6
		att_type = '';		att_val = sprintf ('%g ', attribute.Value );
end

n_argin = nargin;
if (nargout)
	n_argin = n_argin - 1;
	if (n_argin == 1),	dump = varname;		end
end
if (n_argin == 1)
	str = sprintf('\t\t:%s = %s%s\n', attribute.Name, att_val, att_type);
else
	str = sprintf('\t\t%s:%s = %s%s\n', varname, attribute.Name, att_val, att_type);
end

if (nargout),	dump = [dump str];
else			fprintf(1, str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dump = dump_global_attributes ( metadata, dump )

if isfield ( metadata, 'Attribute' )
	num_atts = length(metadata.Attribute);
else
	num_atts = 0;
end

if (num_atts > 0)
	str = sprintf('//global attributes:\n' );
	if (nargin == 2),	dump = [dump str];
	else				fprintf(1, str);
	end
end

for (k = 1:num_atts)
	if (nargin == 1)
		dump_single_attribute ( metadata.Attribute(k) );
	else
		dump = dump_single_attribute ( metadata.Attribute(k), dump );
	end
end

if (nargout),	dump = [dump sprintf('}\n' )];
else			fprintf ( 1, '}\n' );
end
		
% --------------------------------------------------------------------------------		
function attribute = nc_get_attribute_struct ( cdfid, varid, attnum )
% NC_GET_ATTRIBUTE_STRUCT:  Returns a NetCDF attribute as a structure
%
% You don't want to be calling this routine directly.  Just don't use 
% it.  Use nc_attget instead.  Go away.  Nothing to see here, folks.  
% Move along, move along.
%
% USAGE:  attstruct = nc_get_attribute_struct ( cdfid, varid, attnum );
%
% PARAMETERS:
% Input:
%     cdfid:  NetCDF file id
%     varid:  NetCDF variable id
%     attnum:  number of attribute
% Output:
%     attstruct:  structure with "Name", "Nctype", "Attnum", and "Value" fields
%
% In case of an error, an exception is thrown.
%
% USED BY:  nc_getinfo.m, nc_getvarinfo.m
%

% Fill the attribute struct with default values
attribute.Name = '';
attribute.Nctype = NaN;
attribute.Attnum = attnum;   % we know this at this point
attribute.Value = NaN;       % In case the routine fails?

[attname, status] = mexnc('INQ_ATTNAME', cdfid, varid, attnum);
if status < 0 
	snc_error ( 'NC_FUNS:nc_get_attribute_struct:INQ_ATTNAME', sprintf('%s:  mexnc:inq_attname failed on varid %d.\n', mfilename, varid) );
end
attribute.Name = attname;

[att_datatype, status] = mexnc('INQ_ATTTYPE', cdfid, varid, attname);
if status < 0 
	snc_error ( 'NC_FUNS:nc_get_attribute_struct:INQ_ATTTYPE', sprintf('%s:  mexnc:inq_att failed on varid %d, attribute %s.\n', mfilename, varid, attname) );
end
attribute.Nctype = att_datatype;

switch att_datatype
	case 0
		attval = NaN;
	case nc_char
		[attval, status] = mexnc('get_att_text',cdfid,varid,attname);
	case { nc_double, nc_float, nc_int, nc_short, nc_byte }
		[attval, status] = mexnc('get_att_double',cdfid,varid,attname);
	otherwise
		snc_error ( 'NC_FUNS:nc_get_attribute_struct', sprintf('att_datatype is %d.\n', att_datatype) );
end
if status < 0 
	snc_error ( 'NC_FUNS:nc_get_attribute_struct', sprintf('%s:  mexnc:attget failed on varid %d, attribute %s.\n', mfilename, varid, attname) );
end

% this puts the attribute into the variable structure
attribute.Value = attval;

% ----------------------------------------------------------------------------------------
function values = nc_attget(ncfile, varname, attribute_name )
% NC_ATTGET: Get the values of a NetCDF attribute.
%
% USAGE:  att_value = nc_attget(ncfile, varname, attribute_name);
%
% PARAMETERS:
% Input:
%   ncfile:  
%       name of netcdf file in question
%   varname:  
%       name of variable in question.  Specify nc_global to retrieve a 
%       global attribute.  Do NOT use 'global'.
%   attribute_name:  
%       name of attribute in question
% Output:    
%   values:  
%       value of attribute asked for.  Returns the empty matrix 
%       in case of an error.  There is an ambiguity in the case of 
%       NC_BYTE data, so it is always retrieved as an int8 datatype.
%       If you wanted uint8, then you must cast it yourself.
%
% You can specify that java be used instead of the mex-file by setting
% the appropriate preference, i.e.
%     >> setpref('NC_FUNS','USE_JAVA',true);
%
% Example:
%    values = nc_attget('foo.nc', 'x', 'scale_factor')
%
% SEE ALSO:  NC_GLOBAL

snc_nargchk(3,3,nargin);
snc_nargoutchk(1,1,nargout);

[ncid, status] = mexnc('open', ncfile, nc_nowrite_mode );
if ( status ~= 0 )
	snc_error ( 'NC_FUNS:NC_ATTGET:MEXNC:OPEN', mexnc('STRERROR', status) );
end

switch class(varname)
	case {'double'},	varid = varname;
	case 'char',		varid = figure_out_varid ( ncid, varname );
	otherwise
		snc_error ( 'NC_FUNS:NC_ATTGET:badType', 'Must specify either a variable name or NC_GLOBAL' );
end

funcstr = determine_funcstr_attget(ncid,varid,attribute_name);

% And finally, retrieve the attribute.
[values, status] = mexnc(funcstr,ncid,varid,attribute_name);
if ( status ~= 0 )
	snc_error(['NC_FUNS:NC_ATTGET:MEXNC:' funcstr ], mexnc('STRERROR', status) );
end

status = mexnc('close',ncid);
if ( status ~= 0 )
	snc_error ( 'NC_FUNS:NC_ATTGET:MEXNC:CLOSE', mexnc('STRERROR', status) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function funcstr = determine_funcstr_attget(ncid,varid,attribute_name)
% This function is for the mex-file backend.  Determine which netCDF function
% string we invoke to retrieve the attribute value.

[dt, status] = mexnc('inq_atttype',ncid,varid,attribute_name);
if ( status ~= 0 )
	mexnc('close',ncid);
	snc_error ( 'NC_FUNS:NC_ATTGET:MEXNC:INQ_ATTTYPE', mexnc('STRERROR', status) );
end

switch ( dt )
	case nc_double,		funcstr = 'GET_ATT_DOUBLE';
	case nc_float,		funcstr = 'GET_ATT_FLOAT';
	case nc_int,		funcstr = 'GET_ATT_INT';
	case nc_short,		funcstr = 'GET_ATT_SHORT';
	case nc_byte,		funcstr = 'GET_ATT_SCHAR';
	case nc_char,		funcstr = 'GET_ATT_TEXT';
	otherwise
		mexnc('close',ncid);
		snc_error ( 'NC_FUNS:NC_ATTGET:badDatatype', sprintf ( 'unhandled datatype ID %d', dt ) );
end

%===============================================================================
% Did the user do something really stupid like say 'global' when they meant NC_GLOBAL?
function varid = figure_out_varid ( ncid, varname )

if isempty(varname)
	varid = nc_global;
	return
end

if ( strcmpi(varname,'global') )
    [varid, status] = mexnc ( 'inq_varid', ncid, varname );
	if status 
		% Ok, the user meant NC_GLOBAL
		warning ( 'NC_FUNS:nc_attget:doNotUseGlobalString', ...
		          'Please consider using the m-file NC_GLOBAL.M instead of the string ''%s''.', varname );
		varid = nc_global;
		return
	end
end

[varid, status] = mexnc ( 'inq_varid', ncid, varname );
if ( status ~= 0 )
	mexnc('close',ncid);
	snc_error ( 'NC_FUNS:NC_ATTGET:MEXNC:INQ_VARID', mexnc('STRERROR', status) );
end

%----------------------------------------------------------------------------
function nc_varput( ncfile, varname, data, varargin )
% NC_VARPUT:  Writes data into a netCDF file.
%
% NC_VARPUT(NCFILE,VARNAME,DATA) writes the matlab variable DATA to
% the variable VARNAME in the netCDF file NCFILE.  The main requirement
% here is that DATA have the same dimensions as the netCDF variable.
%
% NC_VARPUT(NCFILE,VARNAME,DATA,START,COUNT) writes DATA contiguously, 
% starting at the zero-based index START and with extents given by COUNT.
%
% NC_VARPUT(NCFILE,VARNAME,DATA,START,COUNT,STRIDE) writes DATA  
% starting at the zero-based index START with extents given by
% COUNT, but this time with strides given by STRIDE.  If STRIDE is not
% given, then it is assumes that all data is contiguous.
%
% EXAMPLES:
%    Suppose you have a netcdf variable called 'x' of size 6x4.  If you 
%    have an array of data called 'mydata' that is 6x4, then you can 
%    write to the entire variable with 
% 
%        >> nc_varput ( 'foo.nc', 'x', mydata );
%
%    If you wish to only write to the first 2 rows and three columns,
%    you could do the following
%
%        >> subdata = mydata(1:2,1:3);
%        >> nc_varput ( 'foo.nc', 'x', subdata, [0 0], [2 3] );
%

snc_nargchk(3,6,nargin);
snc_nargoutchk(0,0,nargout);

n_in = nargin;		try_to_scale = true;
if (n_in == 4)
	n_in = 3;
	varargin = {[]};
	try_to_scale = false;
end

[start, count, stride] = parse_and_validate_args_varput(ncfile,varname,varargin{:});

[ncid, status] = mexnc('open', ncfile, nc_write_mode);
if (status ~= 0)
    snc_error ( 'NC_FUNS:NC_VARPUT:MEXNC:OPEN', mexnc('STRERROR', status) );
end

% check to see if the variable already exists.  
[varid, status] = mexnc('INQ_VARID', ncid, varname );
if ( status ~= 0 )
    mexnc ( 'close', ncid );
    snc_error ( 'NC_FUNS:NC_VARPUT:MEXNC:INQ_VARID', mexnc('STRERROR', status) );
end

[dud,var_type,nvdims,var_dim,dud, status] = mexnc('INQ_VAR',ncid,varid);
if status ~= 0 
    mexnc ( 'close', ncid );
    snc_error ( 'NC_FUNS:NC_VARPUT:MEXNC:INQ_VAR', mexnc('STRERROR', status) );
end

nc_count = determine_varsize_mex ( ncid, var_dim, nvdims );
validate_input_data_rank ( ncid, varname, data, nvdims, nc_count );
[start, count] = validate_indexing (ncid,nvdims,data,start,count,stride);

% check that the length of the start argument matches the rank of the variable.
if length(start) ~= length(nc_count)
    mexnc ( 'close', ncid );
    fmt = 'Length of START index (%d) does not make sense with a variable rank of %d.\n';
    msg = sprintf ( fmt, length(start), length(nc_count) );
    snc_error ( 'NC_FUNS:NC_VARPUT:badIndexing', msg );
end


% Figure out which write routine we will use.  If the target variable is a singleton, then we must use
% VARPUT1.  If a stride was given, we must use VARPUTG.  Otherwise just use VARPUT.
if nvdims == 0
    write_op = 'put_var1';
elseif n_in == 3
    write_op = 'put_var';
elseif n_in == 5
    write_op = 'put_vara';
elseif n_in == 6
    write_op = 'put_vars';
else
    snc_error ( 'NC_FUNS:NC_VARPUT', 'unhandled write op.  How did we come to this??\n' );
end

validate_input_size_vs_netcdf_size(ncid,data,nc_count,count,write_op);

if ( try_to_scale && ~(isa(data,'int8') || isa(data,'uint8')) )
	[data, did_scale] = handle_scaling(ncid,varid,data);	% Use cvlib_mex
	data = handle_fill_value ( ncid, varid, data );			% WARNING: Operates only in singles or doubles
	if ( did_scale )					% Dangerous case. (J. LUIS)
		var_type = mexnc('INQ_VARTYPE',ncid,varid);
		switch var_type
			case nc_byte,	    if (~isa(data,'int8'))		data = int8(data);		end
			case nc_char,	    if (~isa(data,'uint8'))		data = uint8(data);		end
			case nc_short,	    if (~isa(data,'int16'))		data = int16(data);		end
			case nc_int,	    if (~isa(data,'int32'))		data = int32(data);		end
			case nc_float,	    if (~isa(data,'single'))	data = single(data);	end
		end
	end
end

write_the_data(ncid,varid,start,count,stride,write_op,data);

status = mexnc ( 'close', ncid );
if ( status ~= 0 )
    snc_error ( 'NC_FUNS:NC_VARPUT:CLOSE', mexnc('STRERROR',status));
end

% ----------------------------------------------------------------------------------
function [start, count, stride] = parse_and_validate_args_varput(ncfile,varname,varargin)
% Set up default outputs.
start = [];
count = [];
stride = [];

switch length(varargin)
	case 2
		start = varargin{1};
		count = varargin{2};
	case 3
		start = varargin{1};
		count = varargin{2};
		stride = varargin{3};
end

% Error checking on the inputs.
if ~ischar(ncfile)
    snc_error ( 'NC_FUNS:NC_VARPUT:badInput', 'the filename must be character.' );
end
if ~ischar(varname)
    snc_error ( 'NC_FUNS:NC_VARPUT:badInput', 'the variable name must be character.' );
end

if ~isnumeric ( start )
    snc_error ( 'NC_FUNS:NC_VARPUT:badInput', 'the ''start'' argument must be numeric.' );
end
if ~isnumeric ( count )
    snc_error ( 'NC_FUNS:NC_VARPUT:badInput', 'the ''count'' argument must be numeric.' );
end
if ~isnumeric ( stride )
    snc_error ( 'NC_FUNS:NC_VARPUT:badInput', 'the ''stride'' argument must be numeric.' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function validate_input_data_rank ( ncid, varname, data, nvdims, nc_count )
% VALIDATE_INPUT_DATA_RANK:
%     Make sure that the rank of the input data matches up with the netCDF
%     variable.  There are a few different cases to consider.

if (nvdims == 0)
    % singleton case.
    if numel(data) ~= 1
        mexnc ( 'close', ncid );
        snc_error ( 'NC_FUNS:NC_VARPUT:badRank', 'Input data size for this variable should be a scalar.' );
    end

elseif (nvdims == 1)
    % This is a special case.  Matlab reports 1D vectors as having at least rank 2.
    if ndims(data) > 2
        mexnc ( 'close', ncid );
        snc_error ( 'NC_FUNS:NC_VARPUT:badRank', 'Input data size for this variable should be a 1D vector.' );
    end

    % At least one of the input data sizes must be 1
    sz = size(data);
    if ~any(find(sz==1))
        mexnc ( 'close', ncid );
        snc_error ( 'NC_FUNS:NC_VARPUT:badRank', 'Input data size for this variable should be a 1D vector.' );
    end

else
	% Trim any trailing singleton dimension.s
	n = length(nc_count);
	for k = n:-1:ceil(n/2)
		if ( nc_count(k) ~= 1 )
			break;
		end
		nc_count(k) = 0;
	end

	%for k = 1:-1:floor(n/2)
	%	if ( nc_count(k) ~= 1 )
	%		break;
	%	end
	%	nc_count(k) = 0;
	%end
	effective_nc_rank = numel(find(nc_count));
	effective_nc_rank = max(2,effective_nc_rank);

    % 2D and higher case.  These cases are easier.
    if ( ndims(data) ~= effective_nc_rank)
% 		mexnc ( 'close', ncid );
% 		efmt = 'Rank of input data (%d) does not work for netCDF variable %s.\n';
% 		snc_error ( 'NC_FUNS:NC_VARPUT:badRank', sprintf(efmt, ndims(data), varname) );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [start, count] = validate_indexing(ncid,nvdims,data,start,count,stride)
% Singletons are a special case.  We need to set the start and count carefully.
if nvdims == 0

    if isempty(start) && isempty(count) && isempty(stride)
        % This is the case of "nc_varput ( file, var, data );"
        start = 0;
        count = 1;

    elseif ~isempty(start) && ~isempty(count) && isempty(stride)
        % This is the case of "nc_varput ( file, var, data, start, count );"
        % So the user gave us "start" and "count".  Did they do so correctly?
        if ( start ~= 0 ) || (count ~= 1 )
            mexnc ( 'close', ncid );
            err_msg = 'In case of singleton variable,  ''start'' must be [0] and ''count'' must be [1].';
            snc_error ( 'NC_FUNS:NC_VARPUT:MEXNC:badIndexing', err_msg );
        end

    elseif ~isempty(start) && ~isempty(count) && ~isempty(stride)
        mexnc ( 'close', ncid );
        err_msg = 'Strides make no sense for a singleton variable.';
        snc_error ( 'NC_FUNS:NC_VARPUT:MEXNC:badIndexing', err_msg );
    end

end

% If START and COUNT not given, and if not a singleton variable, then START is [0,..] and COUNT is 
% the size of the data.  
if isempty(start) && isempty(count) && ( nvdims > 0 )
    start = zeros(1,nvdims);
    count = zeros(1,nvdims);
	for j = 1:nvdims
		count(j) = size(data,j);
	end
end

% Check that the start, count, and stride arguments have the same length.
if ( numel(start) ~= numel(count) )
    mexnc ( 'close', ncid );
    err_msg = 'START and COUNT arguments must have the same length.';
    snc_error ( 'NC_FUNS:NC_VARPUT:MEXNC:badIndexing', err_msg );
end
if ( ~isempty(stride) && (length(start) ~= length(stride)) )
    mexnc ( 'close', ncid );
    err_msg = 'START, COUNT, and STRIDE arguments must have the same length.';
    snc_error ( 'NC_FUNS:NC_VARPUT:MEXNC:badIndexing', err_msg );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, did_scale] = handle_scaling(ncid,varid,data)
% HANDLE_MEX_SCALING
%	If there is a scale factor and/or  add_offset attribute, convert the data
%	to double precision and apply the scaling. The DID_SCALE informs the caller
%	that a scalling and/or offset was performed

did_scale = false;
[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'scale_factor' );
if ( status == 0 ),		have_scale_factor = 1;
else					have_scale_factor = 0;
end

[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'add_offset' );

if ( status == 0 ),		have_add_offset = 1;
else					have_add_offset = 0;
end

% Return early if we don't have either one.
if ~(have_scale_factor || have_add_offset),		return,		end

scale_factor = 1.0;
add_offset = 0.0;

if have_scale_factor
	[scale_factor, status] = mexnc ( 'get_att_double', ncid, varid, 'scale_factor' );
	if ( status ~= 0 )
	    mexnc ( 'close', ncid );
	    snc_error ( 'NC_FUNS:NC_VARPUT:MEXNC:GET_ATT_DOUBLE', mexnc('STRERROR', status) );
	end
end

if have_add_offset
	[add_offset, status] = mexnc ( 'get_att_double', ncid, varid, 'add_offset' );
	if ( status ~= 0 )
	    mexnc ( 'close', ncid );
	    snc_error ( 'NC_FUNS:NC_VARPUT:MEXNC:GET_ATT_DOUBLE', mexnc('STRERROR', status) );
	end
end

if (add_offset == 0 && scale_factor == 1),		return,		end

[var_type,status] = mexnc('INQ_VARTYPE',ncid,varid);
if status ~= 0 
    mexnc ( 'close', ncid );
    snc_error ( 'NC_FUNS:NC_VARPUT:MEXNC:INQ_VARTYPE', mexnc('STRERROR', status) );
end

if (have_scale_factor && have_add_offset)
	data = cvlib_mex('CvtScale',data, 1 / scale_factor, add_offset);
	did_scale = true;
elseif (have_scale_factor)
	data = cvlib_mex('CvtScale',data, 1 / scale_factor);
	did_scale = true;
else
	data = cvlib_mex('addS',data, add_offset);
end

% data = (double(data) - add_offset) / scale_factor;
% did_scale = true;
% 
% % When scaling to an integer, we should add 0.5 to the data.  Otherwise
% % there is a tiny loss in precision, e.g. 82.7 should round to 83, not .
% switch var_type
% 	case { nc_int, nc_short, nc_byte, nc_char }
% 		data = round(data);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = handle_fill_value(ncid,varid,data)
% Handle the fill value.  We do this by changing any NaNs into
% the _FillValue.  That way the netcdf library will recognize it.
[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, '_FillValue' );

if ( status == 0 )
    switch ( class(data) )
        case 'double',		funcstr = 'get_att_double';
        case 'single',		funcstr = 'get_att_float';
        case 'int32',		funcstr = 'get_att_int';
        case 'int16',		funcstr = 'get_att_short';
        case 'int8',		funcstr = 'get_att_schar';
        case 'uint8',		funcstr = 'get_att_uchar';
        case 'char',		funcstr = 'get_att_text';
        otherwise
			mexnc ( 'close', ncid );
			snc_error ('NC_FUNS:NC_VARPUT:unhandledDatatype', sprintf('Unhandled datatype for fill value, ''%s''.', class(data)) );
    end

    [fill_value, status] = mexnc(funcstr,ncid,varid,'_FillValue' );
    if ( status ~= 0 )
	    mexnc( 'close', ncid );
	    snc_error( ['NC_FUNS:NC_VARPUT:MEXNC:' funcstr], mexnc('STRERROR', status) );
    end

	if ( ~isnan(fill_value) && (isa(data,'single') || isa(data,'double')) )
    	data(isnan(data)) = fill_value;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function validate_input_size_vs_netcdf_size(ncid,data,nc_count,count,write_op)
% VALIDATE_INPUT_SIZE_VS_NETCDF_SIZE
%
% We are now in a position to do a check on the input var size vs. the known 
% size in the netcdf file.  If we don't do this, it would be possible to send
% a larger then expected chunk of data to the netcdf file, have parts of it
% get lopped off in order to fit, and  never be the wiser.

switch ( write_op )

	case 'put_var1'
		% Just check that the length of the input data was 1.
		if numel(data) ~= 1
			mexnc ( 'close', ncid );
			snc_error ( 'NC_FUNS:NC_VARPUT:badInput', 'Length of input data must be 1 for singleton variable.'  );
		end

	case 'put_var'
        % Since 'put_var' writes all the data, check that the extents match up exactly.  
        if ( numel(data) ~= prod(nc_count) )
			% Added the following (stupid) test TO LET GO WITH THE UNLIMITED VAR -- J. LUIS
			rec_dim = mexnc( 'INQ_UNLIMDIM', ncid );
			if ( ~(rec_dim == 2 && numel(data) == 1) )
				mexnc ( 'close', ncid );
				fmt = 'Total number of input datums was %d, but the netcdf variable size is %d elements.';
 				snc_error ( 'NC_FUNS:NC_VARPUT:badInput', sprintf ( fmt, numel(data), prod(nc_count) ) );
			end
        end

	case { 'put_vara', 'put_vars' }
        % Just check that the chunk of data the user gave us is the same
        % size as the given count.  This works for put_vars as well.
        if ( numel(data) ~= prod(count) )
			mexnc ( 'close', ncid );
			fmt = 'Total number of input datums was %d, but the count parameter indicated %d elements.\n';
			snc_error ( 'NC_FUNS:NC_VARPUT:badInput', sprintf( fmt, numel(data), prod(count)) );
        end

	otherwise 
		mexnc ( 'close', ncid );
		snc_error ( 'NC_FUNS:NC_VARPUT:unhandledCase', sprintf('Unhandled write operation family, ''%s''.', write_op) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_the_data(ncid,varid,start,count,stride,write_op,data)

% write the data
pdata = permute(data, fliplr( 1:ndims(data) ));
switch (write_op)
    case 'put_var1'
        switch (class(data)),
            case 'double',		funcstr = 'put_var1_double';
            case 'single',		funcstr = 'put_var1_float';
            case 'int32',		funcstr = 'put_var1_int';
            case 'int16',		funcstr = 'put_var1_short';
            case 'int8',		funcstr = 'put_var1_schar';
            case 'uint8',		funcstr = 'put_var1_uchar';
            case 'char',		funcstr = 'put_var1_text';
            otherwise
                mexnc('close',ncid);
                snc_error('NC_FUNS:NC_VARPUT:unhandledMatlabType', sprintf('unhandled data class %s\n', class(pdata)) );
        end
        status = mexnc (funcstr, ncid, varid, start, pdata );

    case 'put_var'
        switch (class(data))
            case 'double',		funcstr = 'put_var_double';
            case 'single',		funcstr = 'put_var_float';
            case 'int32',		funcstr = 'put_var_int';
            case 'int16',		funcstr = 'put_var_short';
            case 'int8',		funcstr = 'put_var_schar';
            case 'uint8',		funcstr = 'put_var_uchar';
            case 'char',		funcstr = 'put_var_text';
            otherwise
                mexnc('close',ncid);
                snc_error ( 'NC_FUNS:NC_VARPUT:unhandledMatlabType', sprintf('unhandled data class %s\n', class(pdata)) );
        end
        status = mexnc (funcstr, ncid, varid, pdata );
    
    case 'put_vara'
        switch (class(data))
            case 'double',		funcstr = 'put_vara_double';
            case 'single',		funcstr = 'put_vara_float';
            case 'int32',		funcstr = 'put_vara_int';
            case 'int16',		funcstr = 'put_vara_short';
            case 'int8',		funcstr = 'put_vara_schar';
            case 'uint8',		funcstr = 'put_vara_uchar';
            case 'char',		funcstr = 'put_vara_text';
            otherwise
                mexnc('close',ncid);
                snc_error ('NC_FUNS:NC_VARPUT:unhandledMatlabType', sprintf('unhandled data class %s\n', class(pdata)) );
        end
        status = mexnc(funcstr, ncid, varid, start, count, pdata );

    case 'put_vars'
        switch ( class(data) )
			case 'double',		funcstr = 'put_vars_double';
			case 'single',		funcstr = 'put_vars_float';
			case 'int32',		funcstr = 'put_vars_int';
			case 'int16',		funcstr = 'put_vars_short';
			case 'int8',		funcstr = 'put_vars_schar';
			case 'uint8',		funcstr = 'put_vars_uchar';
			case 'char',		funcstr = 'put_vars_text';
			otherwise
				mexnc('close',ncid);
				snc_error ('NC_FUNS:NC_VARPUT:unhandledMatlabType', sprintf('unhandled data class %s\n', class(pdata)) );
        end
        status = mexnc(funcstr, ncid, varid, start, count, stride, pdata );

    otherwise 
        mexnc('close', ncid);
        snc_error('NC_FUNS:NC_VARPUT:unhandledWriteOp', sprintf('unknown write operation''%s''.\n', write_op) );

end

if (status ~= 0)
	mexnc('close', ncid);
	snc_error('NC_FUNS:NC_VARPUT:unhandledMatlabType', sprintf('write operation ''%s'' failed with error ''%s''.\n', write_op, mexnc('STRERROR', status) ) );
end

% ------------------------------------------------------------------------------------------
function nc_addvar ( ncfile, varstruct )
% NC_ADDVAR:  adds a variable to a NetCDF file
%
% USAGE:  nc_addvar ( ncfile, varstruct );
%
% PARAMETERS:
% Input
%    ncfile:
%    varstruct:
%        This is a structure with four fields:
%
%        Name
%        Nctype
%        Dimension
%        Attribute
%
%      "Name" is just that, the name of the variable to be defined.
%
%      "Nctype" should be 
%          'double', 'float', 'int', 'short', or 'byte', or 'char'
%          'NC_DOUBLE', 'NC_FLOAT', 'NC_INT', 'NC_SHORT', 'NC_BYTE', 'NC_CHAR'
%
%      "Dimension" is a cell array of dimension names.
%
%      "Attribute" is also a structure array.  Each element has two
%      fields, "Name", and "Value".
%
% Output: 
%     None.  In case of an error, an exception is thrown.
%
% AUTHOR:
%    john.g.evans.ne@gmail.com

snc_nargchk(2,2,nargin);

if  ~ischar(ncfile) 
	snc_error ('NC_FUNS:NC_ADDVAR:badInput', 'file argument must be character');
end

if ( ~isstruct(varstruct) )
	snc_error ('NC_FUNS:NC_ADDVAR:badInput', '2nd argument must be a structure');
end

varstruct = validate_varstruct (varstruct);

[ncid, status] = mexnc ('open', ncfile, nc_write_mode);
if ( status ~= 0 )
	msg = sprintf ('OPEN failed on %s, ''%s''', ncfile, mexnc('STRERROR', status));
	snc_error ('NC_FUNS:NC_ADDVAR:MEXNC:OPEN', msg);
end

% determine the dimids of the named dimensions
num_dims = length(varstruct.Dimension);
dimids = zeros(num_dims,1);
for j = 1:num_dims
	[dimids(j), status] = mexnc ( 'inq_dimid', ncid, varstruct.Dimension{j} );
	if ( status ~= 0 )
		mexnc ( 'close', ncid );
		snc_error ( 'NC_FUNS:NC_ADDVAR:MEXNC:INQ_DIMID', mexnc('STRERROR', status) );
	end
end

% go into define mode
status = mexnc ( 'redef', ncid );
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	snc_error ( 'NC_FUNS:NC_ADDVAR:MEXNC:REDEF', mexnc('strerror', status) );
end

% We prefer to use 'Datatype' instead of 'Nctype', but we'll try to be backwards compatible.
if isfield(varstruct,'Datatype')
	[varid, status] = mexnc('DEF_VAR', ncid, varstruct.Name, varstruct.Datatype, num_dims, dimids );
else
	[varid, status] = mexnc('DEF_VAR', ncid, varstruct.Name, varstruct.Nctype, num_dims, dimids );
end

%[varid, status] = mexnc('DEF_VAR', ncid, varstruct.Name, varstruct.Nctype, num_dims, dimids );

if ( status ~= 0 )
	mexnc('enddef', ncid);
	mexnc('close', ncid);
	snc_error('NC_FUNS:NC_ADDVAR:MEXNC:DEF_VAR', mexnc('STRERROR', status));
end

if (varstruct.Shuffle || varstruct.Deflate)
    status = mexnc('DEF_VAR_DEFLATE',ncid,varid, varstruct.Shuffle,varstruct.Deflate,varstruct.Deflate);
    if ( status ~= 0 )
        ncerr = mexnc('strerror', status);
        mexnc('enddef', ncid);
        mexnc('close', ncid);
        snc_error('NC_FUNS:NC_ADDVAR:MEXNC:DEF_VAR_DEFLATE', ncerr );
    end
end

% if (~isempty(varstruct.Attribute) && strcmp(varstruct.Attribute(1).Name, '_FillValue'))
% 	attval = varstruct.Attribute(1).Value;
% 	nc_attput_while_open ( ncid, varstruct.Name, '_FillValue', attval );
% 	varstruct.Attribute(1) = [];		% Remove this and let the eventual others follow the normal nc_attput path
% end

% Now just use nc_attput to put in the attributes
for j = 1:length(varstruct.Attribute)
	attname = varstruct.Attribute(j).Name;
	attval = varstruct.Attribute(j).Value;
	nc_attput_while_open ( ncid, varstruct.Name, attname, attval );
end

status = mexnc ( 'enddef', ncid );
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	snc_error ( 'NC_FUNS:NC_ADDVAR:MEXNC:ENDDEF', mexnc('STRERROR', status) );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 )
	snc_error ( 'NC_FUNS:NC_ADDVAR:MEXNC:CLOSE', mexnc('STRERROR', status) );
end

% % Now just use nc_attput to put in the attributes
% for j = 1:length(varstruct.Attribute)
% 	attname = varstruct.Attribute(j).Name;
% 	attval = varstruct.Attribute(j).Value;
% 	nc_attput ( ncfile, varstruct.Name, attname, attval );
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varstruct = validate_varstruct(varstruct)
% Check that required fields are there.
% Must at least have a name.
	if ~isfield (varstruct, 'Name')
		snc_error('NC_FUNS:NC_ADDVAR:badInput', 'structure argument must have at least the ''Name'' field.');
	end

	% Check that required fields are there. Default Nctype is double.
	if (~isfield (varstruct, 'Nctype')),	varstruct.Nctype = 'double';	end

	if (~isfield(varstruct,'Datatype'))
		if ~isfield (varstruct, 'Nctype'),	varstruct.Datatype = 'double';
		else								varstruct.Datatype = varstruct.Nctype;
		end
	end

	% Are there any unrecognized fields?
	fnames = fieldnames(varstruct);
	for (j = 1:length(fnames))
		fname = fnames{j};
		switch (fname)
			case {'Datatype', 'Nctype', 'Name', 'Dimension', 'Attribute', ...
					'FillValue', 'Storage', 'Chunking', 'Shuffle', 'Deflate', 'DeflateLevel'}
				% These are used to create the variable.  They are ok.

			case { 'Unlimited', 'Size', 'Rank' }
				% These come from the output of nc_getvarinfo.  We don't use
				% them, but let's not give the user a warning about	them either.

			otherwise
				fprintf (2, '%s:  unrecognized field name ''%s''.  Ignoring it...\n', mfilename, fname);
		end
	end

	% If the datatype is not a string.
	% Change suggested by Brian Powell
	if (isa(varstruct.Nctype, 'double') && varstruct.Nctype < 8 )
		types={'byte' 'char' 'short' 'int' 'float' 'double' 'ubyte'};
		varstruct.Nctype = char(types(varstruct.Nctype));
		varstruct.Datatype = char(types(varstruct.Datatype));
	end

	% Check that the datatype is known.
	switch (varstruct.Nctype)
		case {'NC_DOUBLE', 'double', ...
			'NC_FLOAT', 'float', ...
			'NC_INT', 'int', ...
			'NC_SHORT', 'short', ...
			'NC_BYTE', 'byte', ...
			'NC_UBYTE', 'ubyte', ...
			'NC_CHAR', 'char'}
			% Do nothing
		case 'single'
			varstruct.Datatype = 'float';
		case 'int32'
			varstruct.Datatype = 'int';
		case 'int16'
			varstruct.Datatype = 'short';
		case 'int8'
			varstruct.Datatype = 'byte';
		case 'uint8'
			varstruct.Datatype = 'ubyte';
		case {'uint16', 'uint32', 'int64', 'uint64'}
			snc_error('NC_FUNS:NC_ADDVAR:notClassicDatatype', ...
				'Datatype ''%s'' is not a classic model datatype.', varstruct.Datatype); 
	otherwise
		snc_error( 'NC_FUNS:NC_ADDVAR:unknownDatatype', sprintf('unknown type ''%s''\n', mfilename, varstruct.Nctype) );
	end

	% Check that required fields are there. Default Dimension is none.  Singleton scalar.
	if (~isfield ( varstruct, 'Dimension' )),	varstruct.Dimension = [];	end

	% Check that required fields are there. Default Attributes are none
	if (~isfield ( varstruct, 'Attribute' )),	varstruct.Attribute = [];	end

	if (~isfield(varstruct,'Storage')),		varstruct.Storage = 'contiguous';	end
	if (~isfield(varstruct,'Chunking')),	varstruct.Chunking = [];	end
	if (~isfield(varstruct,'Shuffle')),		varstruct.Shuffle = 0;		end
	if (~isfield(varstruct,'Deflate')),		varstruct.Deflate = 0;		end
	if (~isfield(varstruct,'DeflateLevel')),varstruct.DeflateLevel = 0;	end

% ---------------------------------------------------------------------------------------------------
function nc_attput_while_open ( ncid, varname, attribute_name, attval )
% Same as nc_attput but works on a ncid of a file that is already open
% We need this to workaround the netCDF bug 187 that would otherwise result in a crash
% when adding teh _FillValue attribute

	if isnumeric(varname)
		varid = varname;
	else
		[varid, status] = mexnc ( 'inq_varid', ncid, varname );
		if ( status ~= 0 )
			mexnc ( 'close', ncid );
			snc_error ( 'NC_FUNS:NC_ATTPUT_WHILE_OPEN:MEXNC:INQ_VARID', mexnc('STRERROR', status) );
		end
	end

	% Figure out which mexnc operation to perform.
	switch class(attval)

		case 'double'
			funcstr = 'put_att_double';
			atttype = nc_double;
		case 'single'
			funcstr = 'put_att_float';
			atttype = nc_float;
		case 'int32'
			funcstr = 'put_att_int';
			atttype = nc_int;
		case 'int16'
			funcstr = 'put_att_short';
			atttype = nc_short;
		case 'int8'
			funcstr = 'put_att_schar';
			atttype = nc_byte;
		case 'uint8'
			funcstr = 'put_att_uchar';
			atttype = nc_byte;
		case 'char'
			funcstr = 'put_att_text';
			atttype = nc_char;
		otherwise
			msg = sprintf ('attribute class %s is not handled by %s', class(attval), mfilename );
			snc_error ( 'NC_FUNS:NC_ATTPUT_WHILE_OPEN:unhandleDatatype', msg );
	end

	status = mexnc ( funcstr, ncid, varid, attribute_name, atttype, length(attval), attval);
	if ( status ~= 0 )
		mexnc ( 'close', ncid );
		snc_error ( ['NC_FUNS:NC_ATTPUT_WHILE_OPEN:MEXNC:' upper(funcstr)], mexnc('STRERROR', status) );
	end

% ---------------------------------------------------------------------------------------------------
function nc_attput ( ncfile, varname, attribute_name, attval )
% NC_ATTPUT:  writes an attribute into a netCDF file
%     NC_ATTPUT(NCFILE,VARNAME,ATTNAME,ATTVAL) writes the data in ATTVAL to
%     the attribute ATTNAME of the variable VARNAME of the netCDF file NCFILE.
%     VARNAME should be the name of a netCDF VARIABLE, but one can also use the
%     mnemonic nc_global to specify a global attribute.  Do not use 'global'.
%
% The attribute datatype will match that of the class of ATTVAL.  So if
% if you want to have a 16-bit short integer attribute, make class of
% ATTVAL to be INT16.

snc_nargchk(4,4,nargin);
snc_nargoutchk(0,0,nargout);

[ncid, status] = mexnc( 'open', ncfile, nc_write_mode );
if  status ~= 0 
	snc_error ( 'NC_FUNS:NC_ATTPUT:MEXNC:badFile', mexnc('STRERROR', status) );
end

% Put into define mode.
status = mexnc ( 'redef', ncid );
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	snc_error ( 'NC_FUNS:NC_ATTPUT:MEXNC:REDEF', mexnc('STRERROR', status) );
end

if isnumeric(varname)
	varid = varname;
else
	[varid, status] = mexnc ( 'inq_varid', ncid, varname );
	if ( status ~= 0 )
		mexnc ( 'close', ncid );
		snc_error ( 'NC_FUNS:NC_ATTPUT:MEXNC:INQ_VARID', mexnc('STRERROR', status) );
	end
end

% Figure out which mexnc operation to perform.
switch class(attval)

	case 'double'
		funcstr = 'put_att_double';
		atttype = nc_double;
	case 'single'
		funcstr = 'put_att_float';
		atttype = nc_float;
	case 'int32'
		funcstr = 'put_att_int';
		atttype = nc_int;
	case 'int16'
		funcstr = 'put_att_short';
		atttype = nc_short;
	case 'int8'
		funcstr = 'put_att_schar';
		atttype = nc_byte;
	case 'uint8'
		funcstr = 'put_att_uchar';
		atttype = nc_byte;
	case 'char'
		funcstr = 'put_att_text';
		atttype = nc_char;
	otherwise
		msg = sprintf ('attribute class %s is not handled by %s', class(attval), mfilename );
		snc_error ( 'NC_FUNS:NC_ATTPUT:unhandleDatatype', msg );
end

status = mexnc ( funcstr, ncid, varid, attribute_name, atttype, length(attval), attval);
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	snc_error ( ['NC_FUNS:NC_ATTPUT:MEXNC:' upper(funcstr)], mexnc('STRERROR', status) );
end

% End define mode.
status = mexnc ( 'enddef', ncid );
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	snc_error ( 'NC_FUNS:NC_ATTPUT:MEXNC:ENDDEF', mexnc('STRERROR', status) );
end

status = mexnc('close',ncid);
if ( status ~= 0 )
	snc_error ( 'NC_FUNS:NC_ATTPUT:MEXNC:CLOSE', mexnc('STRERROR', status) );
end

%--------------------------------------------------------------------------------------
function nc_addhist ( ncfile, attval )
% NC_ADDHIST:  adds text to a global history attribute
%
% NC_ADDHIST(NCFILE,TEXT) adds the TEXT string to the standard convention
% "history" global attribute of the netCDF file NCFILE.  The string is 
% prepended, rather than appended.

snc_nargchk(2,2,nargin);

if ~exist(ncfile,'file')
	snc_error('NC_FUNS:NC_ADDHIST:badFilename', sprintf('%s does not exist', ncfile) );
end
if ~ischar(attval)
	snc_error('NC_FUNS:NC_ADDHIST:badDatatype', 'history attribute addition must be character.' );
end

try
	old_hist = nc_attget ( ncfile, nc_global, 'history' );
catch
	% The history attribute must not have existed.  That's ok.
	old_hist = '';
end

if isempty(old_hist)
	new_history = sprintf('%s:  %s', datestr(now), attval );
else
	new_history = sprintf('%s:  %s\n%s', datestr(now), attval, old_hist );
end
nc_attput ( ncfile, nc_global, 'history', new_history );

% ---------------------------------------------------------------------------------------
function nc_add_dimension ( ncfile, dimension_name, dimension_length )
% NC_ADD_DIMENSION:  adds a dimension to an existing netcdf file
%
% USAGE:  nc_add_dimension ( ncfile, dimension_name, dimension_size );
%
% PARAMETERS:
% Input:
%     ncfile:  path to netcdf file
%     dimension_name:  name of dimension to be added
%     dimension_size:  length of new dimension.  If zero, it will be an
%         unlimited dimension.
% Output:
%     none
%
% In case of an error, an exception is thrown.

snc_nargchk(3,3,nargin);

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if status
	snc_error ( 'NC_FUNS:NC_ADD_DIMENSION:openFailed', mexnc('STRERROR', status) );
end

status = mexnc ( 'redef', ncid );
if status
	mexnc ( 'close', ncid );
	snc_error ( 'NC_FUNS:NC_ADD_DIMENSION:redefFailed', mexnc('STRERROR', status) );
end

[dimid, status] = mexnc ( 'def_dim', ncid, dimension_name, dimension_length );
if status
	mexnc ( 'close', ncid );
	snc_error ( 'NC_FUNS:NC_ADD_DIMENSION:defdimFailed', mexnc('STRERROR', status) );
end

status = mexnc ( 'enddef', ncid );
if status
	mexnc('close', ncid);
	snc_error('NC_FUNS:NC_ADD_DIMENSION:enddefFailed', mexnc('STRERROR', status));
end

status = mexnc('close', ncid);
if status 
	snc_error('NC_FUNS:NC_ADD_DIMENSION:closeFailed', mexnc('STRERROR', status));
end

% ---------------------------------------------------------------------------------------
function nc_create_empty (ncfile, mode)
% NC_CREATE_EMPTY:  creates an empty netCDF file 
%     NC_CREATE_EMPTY(NCFILE,MODE) creates the empty netCDF file NCFILE
%     with the given MODE.  MODE is optional, defaulting to nc_clobber_mode.
%     MODE can be one of the following strings:
%   
%       'clobber'         - deletes existing file, creates netcdf-3 file
%       'noclobber'       - creates netcdf-3 file if it does not already
%                           exist.
%       '64bit_offset'    - creates a netcdf-3 file with 64-bit offset
%       'netcdf4-classic' - creates a netcdf-3 file with 64-bit offset
%
%   MODE can also be a numeric value that corresponds either to one of the
%   named netcdf modes or a numeric bitwise-or of them.
%
%   EXAMPLE:  Create an empty classic netCDF file.
%       nc_create_empty('myfile.nc');
%
%   EXAMPLE:  Create an empty netCDF file with the 64-bit offset mode, but
%   do not destroy any existing file with the same name.
%       mode = bitor(nc_noclobber_mode,nc_64bit_offset_mode);
%       nc_create_empty('myfile.nc',mode);
%
%   EXAMPLE:  Create a netCDF-4 file.  This assumes that you have a 
%   netcdf-4 enabled mex-file.
%       nc_create_empty('myfile.nc','netcdf4-classic');  
%
%   SEE ALSO:  nc_adddim, nc_addvar.

	snc_nargchk(1,2,nargin)

	% Set the default mode if necessary.
	if (nargin == 1),	mode = nc_clobber_mode;		end

	switch(mode)
		case {'clobber', 'noclobber', '64bit_offset'}		% Do nothing, this is ok
		case 'netcdf4-classic'
			mode = nc_netcdf4_classic;
	end

	[ncid, status] = mexnc('CREATE', ncfile, mode);
	if (status ~= 0)
		snc_error('NC_FUNS:NC_CREATE_EMPTY:MEXNC:CREATE', mexnc('STRERROR', status));
	end

	status = mexnc('CLOSE', ncid);
	if (status ~= 0)
		snc_error('NC_FUNS:NC_CREATE:MEXNC:CLOSE', mexnc('STRERROR', status));
	end

% --------------------------------------------------------------------
function nc_cat ( input_ncfiles, output_ncfile, abscissa_var )
% NC_CAT_A:  concatentates a set of netcdf files into ascending order
%
% The concatenation is done only along unlimited variable, which by
% definition have an unlimited dimension.  Variables which do NOT have
% an unlimited dimension are copied over from the first of the input
% netcdf input files.
%
% This m-file is not meant as a replacement for ncrcat or any of Charles
% Zender's terrific NCO tools.  If you need NCO functionality, you should
% get NCO tools from http://nco.sourceforge.net
% 
% USAGE:  nc_cat_a ( input_ncfiles, output_ncfile, abscissa_var )
% 
% PARAMETERS:
%   Input:
%       input_ncfiles:
%           This can be either a cell array of netcdf files, or a text
%           file with one netcdf file per line
%       output_ncfile:
%           This file will be generated from scratch.
%       abscissa_var:
%           Name of an unlimited variable.  Supposing we are dealing
%           with time series, then a good candidate for this would
%           be a variable called, oh, I don't know, maybe "time".  
%   Output:
%       None.  An exception is thrown in case of an error.
%
% The best way to explain this is with simple examples.  Suppose that
% the abscissa_var is "time" and that the other netcdf variable is "tsq".
% Suppose that the first netcdf file has files for "time" and "tsq" of
%
%      time: 0 2  4
%      tsq:  0 4 16
%
% Suppose the 2nd netcdf file has values of
%
%      time:  4  6  8
%      tsq:  18 36 64
%
% Note that the 2nd time series has a different value of "tsq" for the 
% abscissa value of 4.
%
% Running nc_cat_asc will produce a single time series of
% 
%      time:  0   2   4   6   8
%      tsq:   0   4  18  36  64
%
% In other words, the 2nd netcdf file's abscissa/ordinate values take
% precedence.  So the order of your netcdf files matter, and the output
% netcdf file will have unique abscissa values.

snc_nargchk(3,3,nargin);
snc_nargoutchk(0,0,nargout);

% If the first input is of type char and is a file, then read it in.
% At the end of this process, the list of netcdf files to be 
% concatenated is in a cell array.
if ischar(input_ncfiles) && exist(input_ncfiles,'file')
	snc_error('NC_FUNS:NC_CAT','Input files by file was not implemented')
% 	afid = fopen ( input_ncfiles, 'r' );
% 	input_ncfiles = textscan ( afid, '%s' );
elseif iscell ( input_ncfiles )
	% Do nothing
else
	snc_error ('NC_FUNS:NC_CAT', 'first input must be either a text file or a cell array\n' );
end

num_input_files = length(input_ncfiles);

% This is how close the abscissa variable values have to be before they
% are considered to be the same value.
tol = 10*eps;

% Now construct the empty output netcdf file.
ncm = nc_info ( input_ncfiles{1} );
mode = nc_clobber_mode;
[ncid, status] = mexnc ( 'CREATE', output_ncfile, mode );
if status ~= 0
	snc_error ( 'NC_FUNS:NC_CAT:CREATE', mexnc ( 'STRERROR', status ) );
end

status = mexnc ( 'CLOSE', ncid );
if status ~= 0
	snc_error ( 'NC_FUNS:NC_CAT:CLOSE', mexnc ( 'STRERROR', status ) );
end

% Add the dimensions.
% ncm.Dimension(3).Unlimited=1;
for d = 1:numel(ncm.Dimension)
	if ( ncm.Dimension(d).Unlimited )
		nc_add_dimension ( output_ncfile, ncm.Dimension(d).Name, 0 );
	else
		nc_add_dimension ( output_ncfile, ncm.Dimension(d).Name, ncm.Dimension(d).Length );
	end
end

% Add the variables
% ncm.Dataset(3).Unlimited=1;
% ncm.Dataset(4).Unlimited=1;
for v = 1:numel(ncm.Dataset)
	nc_addvar ( output_ncfile, ncm.Dataset(v) );

	% If the variable is NOT unlimited, then we can copy over its data now
	if ~ncm.Dataset(v).Unlimited
		vardata = nc_varget ( input_ncfiles{1}, ncm.Dataset(v).Name );
		nc_varput ( output_ncfile, ncm.Dataset(v).Name, vardata );
	end
end

% Add the attributes
for (j = 1:numel(ncm.Attribute))
	nc_attput(output_ncfile, nc_global, ncm.Attribute(j).Name, ncm.Attribute(j).Value );
end

% Go thru and figure out how much data we are looking at, then pre-allocate for speed.
total_length = 0;
for (j = 1:num_input_files)
	sz = nc_varsize ( input_ncfiles{j}, abscissa_var );
	total_length = total_length + sz;
end

abscissa_vardata = NaN*ones(total_length,1);
file_index = NaN*ones(total_length,1);
infile_abscissa_varindex = NaN*ones(total_length,1);

% Now read in the abscissa variable for each file.
start_index = 1;
use_fake_time = false;			% To signal when "time" var repeats between two contiguous files
last_first_v = nan;		end_index = nan;		% Just to shut up the compiler
for j = 1:num_input_files
	v = double(nc_varget ( input_ncfiles{j}, abscissa_var ));
	if (j > 1 && v(1) == last_first_v)		% Poor patch for when the two contiguous files have the the same "time" origin
		v = v + abscissa_vardata(end_index) + diff(abscissa_vardata(1:2));
		use_fake_time = true;
	end
	nv = length(v);

	end_index = start_index + nv - 1;
	inds = start_index:end_index;
	if (use_fake_time),		jump_index = inds(1);	end		% If not two input files result is unknown

	abscissa_vardata(inds) = v;
	file_index(inds) = j*ones(nv,1);
	infile_abscissa_varindex(inds) = (0:nv-1)';
	start_index = start_index + nv;
	last_first_v = v(1);
end

% Sort the ascissa_vardata into ascending order.  
[abscissa_vardata,I] = sort ( abscissa_vardata );
file_index = file_index(I);
infile_abscissa_varindex = infile_abscissa_varindex(I);

% Are there any duplicates?
ind = find ( diff(abscissa_vardata) < tol );
if ~isempty(ind)
	abscissa_vardata(ind) = [];
	file_index(ind) = [];
	infile_abscissa_varindex(ind) = [];
end

% So now go thru each record and append it to the output file and we are done.
for j = 1:length(abscissa_vardata)
	ncfile = input_ncfiles{file_index(j)};
	start = infile_abscissa_varindex(j);
	input_record = nc_getbuffer ( ncfile, start, 1 );
	if (use_fake_time && j >= jump_index)			% PATCH for the case that the "time" variable repeats in the two files
		input_record.(abscissa_var) = abscissa_vardata(j);
	end
	nc_addnewrecs ( output_ncfile, input_record, abscissa_var );
end

% --------------------------------------------------------------------
function new_data = nc_addnewrecs ( ncfile, input_buffer, record_variable )
% NC_ADDNEWRECS:  Tacks on new data from simple matlab structure to an unlimited-dimension netcdf file
% 
% The difference between this m-file and nc_add_recs is that this 
% routine assumes that the unlimited dimension has a monotonically
% increasing coordinate variable, e.g. time series.  This routine
% actually calls nc_add_recs with suitable arguments.
%
% If the length of the record variable data that is to be appended is
% just one, then a check is made for the rest of the incoming data to
% make sure that they also have the proper rank.  This addresses the
% issue of squeezed-out leading singleton dimensions.
%
% From this point foreward, assume we are talking about time series.
% It doesn't have to be that way (the record variable could be 
% monotonically increasing spatially instead ), but talking about it
% in terms of time series is just easier.  If a field is present in 
% the structure, but not in the netcdf file, then that field is 
% ignored.  Only data that is more recent than the last record 
% currently in the NetCDF file is written.   Older data is discarded.
%
% USAGE:  new_data = nc_addnewrecs ( ncfile, input_buffer, record_variable )
% 
% PARAMETERS:
%   Input:
%      ncfile:  
%          netcdf file that we write information to
%      input_buffer:  
%          structure of time series data.  
%      record_variable:
%          Coordinate variable that is monotonically increasing.  
%          In ROMS, it is "ocean_time".  For purposes of backwards
%          compatibility, if this is not provided, it is assumed
%          to be "time".
%   Output:
%      new_data:  
%          Matlab structure of data corresponding in structure to "input_buffer", but
%          consisting only of those records which were actually written to file.
%  
%   The dimensions of the data should match that of the target netcdf file.  For example, 
%   suppose an ncdump of the
%   NetCDF file looks something like
%
%       netcdf a_netcdf_file {
%       dimensions:
%       	lat = 1 ;
%       	lon = 2 ;
%       	depth = 2 ; 
%       	time = UNLIMITED ; // (500 currently)
%       variables:
%       	double time(time) ;
%       	float var1(time, depth) ;
%       	float var2(time, depth, lat, lon) ;
%       	float var3(time, depth, lat) ;
%       
%       // global attributes:
%       }
% 
%   The "input_input_buffer" should look something like the following:
%
%       >> input_input_buffer
%
%       input_input_buffer =
%
%           time: [3x1 double]
%           var1: [3x2 double]
%           var2: [4-D double]
%           var3: [3x2 double]
%
% The reason for the possible size discrepency here is that matlab will
% ignore trailing singleton dimensions (but not interior ones, such as
% that in var2.
%
% If a netcdf variable has no corresponding field in the input input_buffer,
% then the corresponding NetCDF variable will populate with the appropriate
% _FillValue for each new time step.
%          
% In case of an error, an exception is thrown.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: nc_funs.m 7888 2016-05-07 21:55:25Z j $
% $LastChangedDate: 2007-04-23 09:05:21 -0400 (Mon, 23 Apr 2007) $
% $LastChangedRevision: 2178 $
% $LastChangedBy: johnevans007 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_data = [];
snc_nargchk(2,3,nargin);
snc_nargoutchk(0,1,nargout);

if nargin == 2
	record_variable = 'time';
end

if isempty ( input_buffer ),	return,		end

% Check that the record variable is present in the input buffer.
if ~isfield ( input_buffer, record_variable )
	snc_error ( 'NC_FUNS:NC_ADDNEWRECS:missingRecordVariable', ...
	        sprintf('input structure is missing the record variable ''%s''.\n', record_variable) );
end

% check to see that all fields are actually there.
nc = nc_info ( ncfile );
num_nc_vars = length(nc.Dataset);


fnames = fieldnames ( input_buffer );
num_fields = length(fnames);
for j = 1:num_fields
	not_present = 1;
	for k = 1:num_nc_vars
		if strcmp(fnames{j}, nc.Dataset(k).Name)
			not_present = 0;
		end
	end
	if not_present
		fprintf ( 1, '  %s not present in file %s.  Ignoring it...\n', fnames{j}, ncfile );
		input_buffer = rmfield ( input_buffer, fnames{j} );
	end
end

% If the length of the record variable data to be added is just one,
% then we may have a special corner case.  The leading dimension might
% have been squeezed out of the other variables.  MEXNC wants the rank
% of the incoming data to match that of the infile variable.  We address 
% this by forcing the leading dimension in these cases to be 1.
input_buffer = force_leading_dimension ( ncfile, input_buffer, record_variable );

% Retrieve the dimension id of the unlimited dimension upon which
% all depends.  It must be the first dimension listed.
varinfo = nc_getvarinfo ( ncfile, record_variable );
unlimited_dimension_name = varinfo.Dimension{1};

% Get the last time value.   If the record variable is empty, then
% only take datums that are more recent than the latest old datum
input_buffer_time_values = input_buffer.(record_variable);
if varinfo.Size > 0
	last_time = nc_getlast ( ncfile, record_variable, 1 );
    recent_inds = find( input_buffer_time_values > last_time );
else
    recent_inds = 1:length(input_buffer_time_values);
end

% if no data is new enough, just return.  There's nothing to do.
if isempty(recent_inds),	return,		end

% Go thru each variable.  Restrict to what's new.
varnames = fieldnames ( input_buffer );
for j = 1:numel(varnames)
	data = input_buffer.(varnames{j});
	current_varsize = size(data);
	new_output_size = [length(recent_inds) current_varsize(2:end)];
	restricted_data = data(recent_inds,:);
	restricted_data = reshape ( restricted_data, new_output_size );
	input_buffer.(varnames{j}) = restricted_data;
end

% Write the records out to file.
nc_add_recs ( ncfile, input_buffer, unlimited_dimension_name );

new_data = input_buffer;

%==============================================================================
function input_buffer = force_leading_dimension ( ncfile, input_buffer, record_variable )

varnames = fieldnames ( input_buffer );
num_vars = length(varnames);
if length(input_buffer.(record_variable)) == 1 
	for j = 1:num_vars
		% Skip the record variable, it's irrelevant at this stage.
		if strcmp ( varnames{j}, record_variable ),		continue,	end

		infile_vsize = nc_varsize(ncfile, varnames{j} );

		% Disregard any trailing singleton dimensions.
		effective_nc_rank = calculate_effective_nc_rank(infile_vsize);

		% The input data has to have at least 2 as a rank.
		mlrank = calculate_mlrank ( input_buffer, varnames, j );

		if (effective_nc_rank == mlrank)
			% Do nothing.  The data is fine in this case.
			% This would be an example of a 1D variable, where we
			% are tacking on a single value. Or it is the case of a row vector.

		elseif ( effective_nc_rank == mlrank+1 )
			% In this case, the infile definition is, e.g.,
			% 10x5x5, or in other words, has rank 3.
			% But the incoming data is just a single timestep,
			% e.g. 5x5.  So we change it into a 1x5x5.
			tmp(1,:,:) = input_buffer.(varnames{j});
			input_buffer.(varnames{j}) = tmp;

		end
	end
end

%==============================================================================
function effective_nc_rank = calculate_effective_nc_rank(infile_vsize)

	n = length(infile_vsize);

	% Trim any trailing singleton dimensions. Do this by zeroing them out.
	for k = n:-1:ceil(n/2)
		if (infile_vsize(k) ~= 1),		break,		end
		infile_vsize(k) = 0;
	end

	% Don't get fooled if there is no data in the file.
	if ( infile_vsize(1) == 0 )
		infile_vsize(1) = -1;
	end
	effective_nc_rank = numel(find(infile_vsize));

%==============================================================================
function mlrank = calculate_mlrank ( input_buffer, varnames, j )
% If the rank of the file variable and the data is different,  then we assume two
% conditions have to hold before we augment the data with a leading singleton dimension.
% The extent of the incoming data must not be one, and the length of the size of the
% incoming data must not match up with the length of the size of the file variable.
	if ndims(input_buffer.(varnames{j})) == 2
		sz = size(input_buffer.(varnames{j}));
		if (sz(1) == 1) && (sz(2) == 1)
			mlrank = 1;
		else
			mlrank = 2;
		end
	else
		mlrank = length(size(input_buffer.(varnames{j})));
	end

% --------------------------------------------------------------------
function nc_add_recs ( ncfile, new_data, varargin )
% NC_ADD_RECS:  add records onto the end of a netcdf file
%
% USAGE:  nc_add_recs ( ncfile, new_data, unlimited_dimension );
% 
% INPUT:
%   ncfile:  netcdf file
%   new_data:  Matlab structure.  Each field is a data array
%      to be written to the netcdf file.  Each array had
%      better be the same length.  All arrays are written
%      in the same fashion.
%   unlimited_dimension:
%      Optional.  Name of the unlimited dimension along which the data 
%      is written.  If not provided, we query for the first unlimited 
%      dimension (looking ahead to HDF5/NetCDF4).
%     
% OUTPUT:
%   None.  In case of an error, an exception is thrown.
%
% AUTHOR: 
%   johnevans@acm.org

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: nc_funs.m 7888 2016-05-07 21:55:25Z j $
% $LastChangedDate: 2007-08-31 16:30:56 -0400 (Fri, 31 Aug 2007) $
% $LastChangedRevision: 2309 $
% $LastChangedBy: johnevans007 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snc_nargchk(2,3,nargin);

% Check that we were given good inputs.
if ~isstruct ( new_data )
	snc_error ( 'NC_FUNS:NC_ADD_RECS:badStruct', '2nd input argument must be a structure .\n' );
end

% Check that each field of the structure has the same length.
varnames = fieldnames ( new_data );
num_fields = length(varnames);
if ( num_fields <= 0 )
	snc_error ( 'NC_FUNS:NC_ADD_RECS:badRecord', 'data record cannot be empty' );
end
field_length = zeros(num_fields,1);
for j = 1:num_fields
	command = sprintf ( 'field_length(j) = size(new_data.%s,1);', varnames{j} );
	eval ( command );
end
if any(diff(field_length))
	snc_error ( 'NC_FUNS:NC_ADD_RECS:badFieldLengths', 'Some of the fields do not have the same length.\n' );
end

% So we have this many records to write.
record_count(1) = field_length(1);


[unlim_dimname, unlim_dimlen, unlim_dimid] = get_unlimdim_info ( ncfile, varargin{:} );

varsize = get_all_varsizes ( ncfile, new_data, unlim_dimid );

% So we start writing here.
record_corner(1) = unlim_dimlen;

% write out each data field, as well as the minimum and maximum
input_variable = fieldnames ( new_data );
num_vars = length(input_variable);
for i = 1:num_vars

	current_var = input_variable{i};
	%fprintf ( 1, '%s:  processing %s...\n', mfilename, current_var );

	current_var_data = new_data.(current_var);
	var_buffer_size = size(current_var_data);

	netcdf_var_size = varsize.(current_var);

	corner = zeros( 1, length(netcdf_var_size) );
	count = ones( 1, length(netcdf_var_size) );

	corner(1) = record_corner(1);
	count(1) = record_count(1);
	
	for j = 2:numel(var_buffer_size)
		if ( var_buffer_size(j) > 1 )
			count(j) = var_buffer_size(j);
		end
	end

	% Ok, we are finally ready to write some data.
	nc_varput ( ncfile, current_var, current_var_data, corner, count );
end

% ===============================================================================
function varsize = get_all_varsizes ( ncfile, new_data,unlimited_dimension_dimid )

[ncid,status ]=mexnc( 'open', ncfile, nc_nowrite_mode );
if status ~= 0
	ncerr = mexnc ( 'strerror', status );
	error_id = 'NC_FUNS:NC_ADD_RECS:openFailed';
	snc_error ( error_id, ncerr );
end

% For each field of "new_data" buffer, inquire as to the dimensions in the
% NetCDF file.  We need this data to properly tell nc_varput how to write the data
input_variable = fieldnames ( new_data );
num_vars = length(input_variable);
varsize = [];
for j = 1:num_vars
	[varid, status] = mexnc('INQ_VARID', ncid, input_variable{j} );
	if ( status ~= 0 )
		mexnc('close',ncid);
		snc_error ( 'NC_FUNS:NC_ADD_RECS:inq_varidFailed', mexnc ( 'strerror', status ) );
	end

	[dimids, status] = mexnc('INQ_VARDIMID', ncid, varid);
	if ( status ~= 0 )
		mexnc('close',ncid);
		snc_error ( 'NC_FUNS:NC_ADD_RECS:inq_vardimidFailed', mexnc ( 'strerror', status ) );
	end
	ndims = length(dimids);
	dimsize = zeros(ndims,1);

	% make sure that this variable is defined along the unlimited dimension.
	if ~any(find(dimids==unlimited_dimension_dimid))
		mexnc('close',ncid);
		format = 'variable %s must be defined along unlimited dimension.\n';
		snc_error ( 'NC_FUNS:NC_ADD_RECS:missingUnlimitedDimension', format, input_variable{j} );
	end

	for k = 1:ndims
		[dim_length, status] = mexnc('INQ_DIMLEN', ncid, dimids(k) );
		if ( status ~= 0 )
			mexnc('close',ncid);
			snc_error ( 'NC_FUNS:NC_ADD_RECS:inq_dimlenFailed', mexnc ( 'strerror', status ) );
		end
		dimsize(k) = dim_length;
	end
	varsize.(input_variable{j}) = dimsize;
end

status = mexnc('close',ncid);
if status ~= 0 
	snc_error ( 'NC_FUNS:NC_ADD_RECS:closeFailed', mexnc ( 'strerror', status ) );
end

% =======================================================================
function [dimname, dimlen, dimid] = get_unlimdim_info ( ncfile, varargin )

[ncid,status ]=mexnc( 'open', ncfile, nc_nowrite_mode );
if status ~= 0
	mexnc('close',ncid);
	ncerr = mexnc ( 'strerror', status );
	snc_error ( 'NC_FUNS:NC_ADD_RECS:openFailed', ncerr );
end

% If we were not given the name of an unlimited dimension, get it now
if nargin < 2
	[dimid, status] = mexnc ( 'inq_unlimdim', ncid );
	if status ~= 0
		mexnc('close',ncid);
		ncerr = mexnc ( 'strerror', status );
		snc_error ( 'NC_FUNS:NC_ADD_RECS:inq_unlimdimFailed', ncerr );
	end

	[dimname, status] = mexnc ( 'INQ_DIMNAME', ncid, dimid );
	if status ~= 0
		mexnc('close',ncid);
		ncerr = mexnc ( 'strerror', status );
		snc_error ( 'NC_FUNS:NC_ADD_RECS:inq_dimnameFailed', ncerr );
	end

	if dimid == -1
		error_id = 'NC_FUNS:NC_ADD_RECS:noUnlimitedDimension';
		snc_error ( error_id, sprintf('%s is missing an unlimited dimension, %s requires it', ncfile, mfilename) );
	end

else
	dimname = varargin{1};
	[dimid, status] = mexnc ( 'inq_dimid', ncid, dimname );
	if status ~= 0
		mexnc('close',ncid);
		ncerr = mexnc ( 'strerror', status );
		snc_error ( 'NC_FUNS:NC_ADD_RECS:OPEN', ncerr );
	end
	
end
	
[dimlen, status] = mexnc ( 'INQ_DIMLEN', ncid, dimid );
if status ~= 0
	mexnc('close',ncid);
	ncerr = mexnc ( 'strerror', status );
	snc_error ( 'NC_FUNS:NC_ADD_RECS:inq_dimlenFailed', ncerr );
end

status = mexnc('close',ncid);
if status ~= 0 
	ncerr = mexnc ( 'strerror', status );
	snc_error ( 'NC_FUNS:NC_ADD_RECS:closeFailed', ncerr );
end

% --------------------------------------------------------------------
function theBuffer = nc_getbuffer ( ncfile, varargin )
% NC_GETBUFFER:  read the unlimited variables of a netcdf file into a structure
%
% USAGE:  theBuffer = nc_getbuffer ( ncfile );
% USAGE:  theBuffer = nc_getbuffer ( ncfile, varlist );
% USAGE:  theBuffer = nc_getbuffer ( ncfile, start, count );
% USAGE:  theBuffer = nc_getbuffer ( ncfile, varlist, start, count );
%
% PARAMETERS:
% INPUT:
%     ncfile:  
%        Input netcdf file name.
%     varlist:
%        cell array of named variables.  Only data for these variables 
%        will be retrieved from the file
%     start, count:
%        starting index and number of records to retrieve.  This is 
%        optional.  If not provided, all of the record variables will
%        be retrieved.
%        
%        If start is negative, then the last few records (total of
%        "count") are retrieved.
%
%        If count is negative, then everything beginning at "start" 
%        and going all the way to the end is retrieved.
% 
% OUTPUT:
%     theBuffer:
%        Structure with fields corresponding to each netcdf record variable.  
%        Each such field contains the data for that variable.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: nc_funs.m 7888 2016-05-07 21:55:25Z j $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% assume failure until success is known
theBuffer = [];

snc_nargchk(1,4,nargin);
snc_nargoutchk(1,1,nargout);

% check that the first argument is a char
if ~ischar ( ncfile )
   	snc_error (  'NC_FUNS:NC_GETBUFFER:badInput', 'filename argument must be character.' );
end

[varlist,start,count] = parse_inputs_getbuf(varargin{:});
metadata = nc_info ( ncfile );
% metadata.Dimension(3).Unlimited=1;
% metadata.Dataset(3).Unlimited=1;
% metadata.Dataset(4).Unlimited=1;
num_datasets = length(metadata.Dataset);
skip_this_variable = construct_skip_list(varlist,metadata);

%
% Find the unlimited dimension and it's length
record_length = -1;
num_dims = length(metadata.Dimension);
for j = 1:num_dims
	if metadata.Dimension(j).Unlimited
		record_length = metadata.Dimension(j).Length;
	end
end
if record_length < 0
   	snc_error (  'NC_FUNS:NC_GETBUFFER:noUnlimitedDimension', 'An unlimited dimension is required.');
end

% figure out what the start and count really are.
if ~isempty(start) && ~isempty(count)
	if start < 0
		start = record_length - count;
	end
	if count < 0
		count = record_length - start;
	end
	if (start < 0) && (count < 0)
   		snc_error (  'NC_FUNS:NC_GETBUFFER:badIndexing', 'both start and count cannot be less than zero.');
	end
end

for (j = 1:num_datasets)
	% Did we restrict retrieval to a few variables?
	if (~isempty(varlist) && skip_this_variable(j)),	continue,		end

	% If it is not an unlimited variable, we don't want it.
	if (~metadata.Dataset(j).Unlimited),	continue,	end

	if ~isempty(start) && ~isempty(count) 
		varstart = zeros(size(metadata.Dataset(j).Size));
		varstart(1) = start;
		varcount = metadata.Dataset(j).Size;
		varcount(1) = count;
		vardata = nc_varget ( ncfile, metadata.Dataset(j).Name, varstart, varcount );
	else
		vardata = nc_varget ( ncfile, metadata.Dataset(j).Name );
	end

	theBuffer.(metadata.Dataset(j).Name) = vardata;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varlist, start, count] = parse_inputs_getbuf( varargin )

varlist = {};	start = [];		count = [];

% figure out what the inputs actually were
switch nargin
case 1
	if iscell(varargin{1})
		varlist = varargin{1};
	else
		snc_error ( 'NC_FUNS:NC_GETBUFFER:badInput', '2nd of two input arguments must be a cell array.' );
	end
case 2
	if isnumeric(varargin{1}) && isnumeric(varargin{2})
		start = varargin{1};
		count = varargin{2};
	else
		snc_error ( 'NC_FUNS:NC_GETBUFFER:badInput', '2nd and 3rd of three input arguments must be numeric.' );
	end
case 3
	if iscell(varargin{1})
		varlist = varargin{1};
	else
		snc_error ( 'NC_FUNS:NC_GETBUFFER:badInput', '2nd of four input arguments must be a cell array.' );
	end
	if isnumeric(varargin{2}) && isnumeric(varargin{3})
		start = varargin{2};
		count = varargin{3};
	else
		snc_error ( 'NC_FUNS:NC_GETBUFFER:badInput', '3rd and 4th of four input arguments must be numeric.' );
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function skip_it = construct_skip_list(varlist,metadata)

num_datasets = length(metadata.Dataset);
if ~isempty(varlist) 
	skip_it = ones(num_datasets,1);
	% Go thru and quickly set up a flag for each Dataset
	for j = 1:num_datasets
		for k = 1:length(varlist)
			if strcmp(varlist{k}, metadata.Dataset(j).Name)
				skip_it(j) = 0;
			end
		end
	end
else
	skip_it = zeros(num_datasets,1);
end

retrievable_datasets = find(1 - skip_it);
if ~any(retrievable_datasets)
	snc_error ( 'NC_FUNS:skip_list',  'No datasets found.\n' );
end

% -------------------------------------------------------------------
function varsize = nc_varsize(ncfile, varname)
% NC_VARSIZE:  return the size of the requested netncfile variable
%
% VARSIZE = NC_VARSIZE(NCFILE,NCVAR) returns the size of the netCDF variable 
% NCVAR in the netCDF file NCFILE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: nc_funs.m 7888 2016-05-07 21:55:25Z j $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snc_nargchk(2,2,nargin);
snc_nargoutchk(1,1,nargout);

if ~ischar(ncfile)
	snc_error ( 'NC_FUNS:NC_VARSIZE:badInputType', 'The input filename must be a string.' );
end
if ~ischar(varname)
	snc_error ( 'NC_FUNS:NC_VARSIZE:badInputType', 'The input variable name must be a string.' );
end

v = nc_getvarinfo ( ncfile, varname );
varsize = v.Size;

% -------------------------------------------------------------------
function values = nc_getlast(ncfile, var, num_datums)
% NC_GETLAST:  Retrieves records at the end of an unlimited netCDF file
%
% DATA = NC_GETLAST(NCFILE,VARNAME,NUM_DATUMS) retrieves NUM_DATUMS 
% datums from the netCDF variable VARNAME in the netCDF file NCFILE.
% If NUM_DATUMS is not supplied, the default value is 1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: nc_funs.m 7888 2016-05-07 21:55:25Z j $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snc_nargchk(2,3,nargin);
snc_nargoutchk(1,1,nargout);

if ~ischar(ncfile) 
	snc_error ( 'NC_FUNS:NC_GETLAST:badInput', 'The netCDF file argument must be char.' );
end

if ~ischar(var) 
	snc_error ( 'NC_FUNS:NC_GETLAST:badInput', 'The netCDF variable argument must be char.' );
end

if ( nargin == 2 )
	num_datums = 1;
else
	if ~isnumeric(num_datums) 
	    snc_error ( 'NC_FUNS:NC_GETLAST:badInput', 'The num_datums argument must be numeric.' );
	end
	if num_datums <= 0
	    snc_error ( 'NC_FUNS:NC_GETLAST:badInput', 'The num_datums argument must be positive.' );
	end
end

nb = nc_getbuffer (ncfile, {var}, -1, num_datums);
values = nb.(var);


% --------------------------------------------------------------------
function snc_error (error_id, error_msg)
	disp ([error_id ' ' error_msg]);			% On compiled version that's the only message we'll get
	error (error_id, error_msg);

% --------------------------------------------------------------------
function snc_nargchk(low,high,N)
% SNC_NARGCHK:  wrapper for NARGCHK, which changed functionality at R???
	error (nargchk(low,high,N));

% --------------------------------------------------------------------
function snc_nargoutchk(low,high,N)
% SNC_NARGOUTCHK:  wrapper for NARGOUTCHK, which changed functionality at R???
	error (nargoutchk(low,high,N));

% --------------------------------------------------------------------
function check_index_vectors(start,count,stride,nvdims,ncid,varname)
% CHECK_INDEX_VECTORS
%    We need to check the lengths of the index vectors before calling the
%    netCDF library.  A bad length can confuse the mex-file, and it's really 
%    not the mex-file's responsibility to do this.  We can't do this when 
%    parsing the input arguments since we need the number of dimensions 
%    corresponding to the input variable.

	if ~isempty(start) && (length(start) ~= nvdims)
		mexnc('close',ncid);
		fmt = 'length of the start index vector (%d) does not equal the number of dimensions (%d) for %s';
		snc_error ( 'NC_FUNS:NC_VARGET:badStartIndex', sprintf(fmt, length(start), nvdims, varname) );
	end
	if ~isempty(count) && (length(count) ~= nvdims)
		mexnc('close',ncid);
		fmt = 'length of the count index vector (%d) does not equal the number of dimensions (%d) for %s';
		snc_error ( 'NC_FUNS:NC_VARGET:badCountIndex', sprintf ( fmt, length(count), nvdims, varname ) );
	end
	if ~isempty(stride) && (length(stride) ~= nvdims)
		mexnc('close',ncid);
		fmt = 'length of the stride index vector (%d) does not equal the number of dimensions (%d) for %s';
		snc_error ( 'NC_FUNS:NC_VARGET:badStrideIndex', sprintf ( fmt, length(count), nvdims, varname ) );
	end

% --------------------------------------------------------------------
function the_var_size = determine_varsize_mex ( ncid, dimids, nvdims )
% DETERMINE_VARSIZE_MEX: Need to figure out just how big the variable is.
%
% VAR_SIZE = DETERMINE_VARSIZE_MEX(NCID,DIMIDS,NVDIMS);

% If not a singleton, we need to figure out how big the variable is.
if nvdims == 0
    the_var_size = 1;
else
    the_var_size = zeros(1,nvdims);
    for (j = 1:nvdims)
        dimid = dimids(j);
        [dim_size,status] = mexnc('inq_dimlen', ncid, dimid);
        if ( status ~= 0 )
			mexnc('close',ncid);
            snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:INQ_DIM_LEN', mexnc('strerror', status) );
        end
        the_var_size(j)=dim_size;
    end
end

% --------------------------------------------------------------------
function nc_datatype = nc_nat()
% NC_NAT:  returns constant corresponding to NC_NAT enumerated constant in netcdf.h
	nc_datatype = 0;

% --------------------------------------------------------------------
function nc_datatype = nc_byte()
% NC_BYTE:  returns constant corresponding to NC_BYTE enumerated constant in netcdf.h
	nc_datatype = 1;

% --------------------------------------------------------------------
function nc_datatype = nc_char()
% NC_CHAR:  returns constant corresponding to NC_CHAR enumerated constant in netcdf.h
	nc_datatype = 2;

% --------------------------------------------------------------------
function nc_datatype = nc_short()
% NC_SHORT:  returns constant corresponding to NC_SHORT enumerated constant in netcdf.h
	nc_datatype = 3;

% --------------------------------------------------------------------
function nc_datatype = nc_int()
% NC_INT:  returns constant corresponding to NC_INT enumerated constant in netcdf.h
	nc_datatype = 4;

% --------------------------------------------------------------------
function nc_datatype = nc_float()
% NC_FLOAT:  returns constant corresponding to NC_FLOAT enumerated constant in netcdf.h
	nc_datatype = 5;

% --------------------------------------------------------------------
function nc_datatype = nc_double()
% NC_DOUBLE:  returns constant corresponding to NC_DOUBLE enumerated constant in netcdf.h
	nc_datatype = 6;

% --------------------------------------------------------------------
function nc_datatype = nc_ubyte()
% NC_UBYTE:  returns constant corresponding to NC_UBYTE enumerated constant in netcdf.h
	nc_datatype = 7;

% --------------------------------------------------------------------
function nc_datatype = nc_ushort()
% NC_USHORT:  returns constant corresponding to NC_USHORT enumerated constant in netcdf.h
	nc_datatype = 8;

% --------------------------------------------------------------------
function nc_datatype = nc_uint()
% NC_UINT:  returns constant corresponding to NC_UINT enumerated constant in netcdf.h
	nc_datatype = 9;

% --------------------------------------------------------------------
function nc_datatype = nc_int64()
% NC_INT64:  returns constant corresponding to NC_INT64 enumerated constant in netcdf.h
	nc_datatype = 10;

% --------------------------------------------------------------------
function nc_datatype = nc_uint64()
% NC_UINT64:  returns constant corresponding to NC_UINT64 enumerated constant in netcdf.h
	nc_datatype = 11;

% --------------------------------------------------------------------
function nc_datatype = nc_string()
% NC_STRING:  returns constant corresponding to NC_STRING enumerated constant in netcdf.h
	nc_datatype = 12;

% --------------------------------------------------------------------
function the_value = nc_global ()
	the_value = -1;	% returns enumerated constant NC_GLOBAL in netcdf.h

% --------------------------------------------------------------------
function mode = nc_nowrite_mode()
	mode = 0;	% returns integer mnemonic for NC_NOWRITE

% --------------------------------------------------------------------
function mode = nc_clobber_mode()
	mode = 0;	% returns integer mnemonic for NC_CLOBBER

% --------------------------------------------------------------------
function mode = nc_noclobber_mode()
	mode = 4;	% returns integer mnemonic for NC_NOCLOBBER
	
% --------------------------------------------------------------------
function mode = nc_write_mode()
	mode = 1;	% returns integer mnemonic for NC_WRITE

% --------------------------------------------------------------------
function mode = nc_netcdf4_classic()
% NC_NETCDF4_CLASSIC:  returns integer mnemonic for BITOR(NC_NETCDF4,NC_CLASSIC_MODEL)
% Use this when you wish to create a netcdf-4 file with the classic model.
	mode = 4352;

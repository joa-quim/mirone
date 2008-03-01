function  varargout = nc_funs(opt,varargin)
% This contains some of the MEXNC (John Evans) functions packed in a single file
% It is not the whole package. Just the functions we use in Mirone
% Data is returned on its native type. That is, it is not complient with the doubles tirany
% Except in the case of - ints and scale_factor ~= 1 OR have _FillValue and array data have _FillValue
% In that case we return the array variable as a single
% Its is also striped of the Java shit
% Joaquim Luis

	switch opt
		case 'getdiminfo'
			varargout{1} = nc_getdiminfo(varargin{:});
		case 'varget'
			varargout{1} = nc_varget(varargin{:});
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
		case 'dump'
			if (nargout),		varargout{1} = nc_dump(varargin{:});
			else				nc_dump(varargin{:});
			end
	end
		

% --------------------------------------------------------------------
function fileinfo = nc_info( ncfile )
% This function is the core function of the snctools/nc_info

fileinfo.Filename = ncfile;

[ncid, status] = mexnc('open', ncfile, nc_nowrite_mode );
if status ~= 0
    snc_error ( 'NC_INFO:MEXNC:OPEN', mexnc('strerror', status) );
end

[ndims, nvars, ngatts, record_dimension, status] = mexnc('INQ', ncid);
if status ~= 0
    mexnc('close',ncid);
    snc_error( 'NC_FUNS:NC_INFO:MEXNC:INQ', mexnc('strerror', status) );
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
	snc_error ( 'NC_FUNS:NC_GETVARINFO:badTypes', ...
	        'Must have either both character inputs, or both numeric.' );
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

% Assume the current variable does not have an unlimited dimension until
% we know that it does.
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

snc_nargchk(2,5,nargin);
snc_nargchk(1,1,nargout);

[start, count, stride] = parse_and_validate_args_varget(ncfile,varname,varargin{:});

[ncid,status] = mexnc('open',ncfile,'NOWRITE');
if status ~= 0
    snc_error('NC_FUNS:NC_VARGET:MEXNC:OPEN', mexnc('strerror', status) );
end

[varid, status]=mexnc('inq_varid', ncid, varname);
if status ~= 0
    mexnc('close',ncid);
    snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:INQ_VARID', mexnc('strerror', status) );
end

[dud,var_type,nvdims,dimids,dud,status] = mexnc('inq_var', ncid, varid);
if status ~= 0
    mexnc('close',ncid);
    snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:INQ_VAR', mexnc('strerror',status) );
end

% Check that the start, count, stride parameters have appropriate lengths.
% Otherwise we get confusing error messages later on.
check_index_vectors(start,count,stride,nvdims,ncid,varname);

% What mexnc operation will we use?
[funcstr_family, funcstr] = determine_funcstr( var_type, nvdims, start, count, stride );

the_var_size = determine_varsize_mex( ncid, dimids, nvdims );

% Check that START index is ok.
if ~isempty(start) && any(start > the_var_size)
    snc_error('NC_FUNS:NC_VARGET:badStartIndex', 'The START index argument exceeds the size of the variable.');
end

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
        error('NC_FUNS:NC_VARGET:unhandledType', sprintf ('Unhandled function string type ''%s''\n', funcstr_family) );
end

if ( status ~= 0 )
	mexnc('close',ncid);
	snc_error ( 'NC_FUNS:NC_VARGET:%s', mexnc('strerror', status) );
end

% If it's a 1D vector, make it a column vector.  Otherwise permute the data
% to make up for the row-major-order-vs-column-major-order issue.
if (length(the_var_size) == 1)
    values = values(:);
else
    pv = numel(the_var_size):-1:1;
    values = permute(values,pv);
end                                                                                   

values = handle_fill_value_mex ( ncid, varid, var_type, values );
values = handle_mex_missing_value ( ncid, varid, var_type, values );
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

switch ( var_type )
	case nc_int,		funcstr = [prefix '_int'];
	case nc_float,		funcstr = [prefix '_float'];
	case nc_double,		funcstr = [prefix '_double'];
	case nc_short,		funcstr = [prefix '_short'];
	case nc_char,		funcstr = [prefix '_text'];
	case nc_byte,		funcstr = [prefix '_schar'];
	otherwise
		error ( 'NC_FUNS:NC_VARGET:badDatatype', sprintf('Unhandled datatype %d.', var_type) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function values = handle_fill_value_mex ( ncid, varid, var_type, values )
% HANDLE_MEX_FILL_VALUE: If there is a fill value, then replace such values with NaN.

[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, '_FillValue' );
if ( status == 0 )

    switch ( var_type )
        case {nc_char, nc_byte}
			% For now, do nothing.  Does a fill value even make sense with char data?
			% If it does, please tell me so.
        case { nc_int, nc_short}
			[fill_value, status] = mexnc( 'get_att_double', ncid, varid, '_FillValue' );
			if (~isnan(fill_value))
				ind = (values == fill_value);
				if (any(ind))
					values = single(values);		% Here I give up and convert to singles
					values(ind) = NaN;
				end
			else
				% An idiotic case. A _FillValue = NAN in a array of integers.
			end
        case { nc_double, nc_float }
			[fill_value, status] = mexnc( 'get_att_double', ncid, varid, '_FillValue' );
			if (~isnan(fill_value))
				values(values == fill_value) = NaN;
			end
        otherwise
			mexnc('close',ncid);
			error ( sprintf('unhandled datatype %d\n', var_type) );
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
			[fill_value, status] = mexnc( 'get_att_double', ncid, varid, 'missing_value' );
			if (~isnan(fill_value))
				ind = (values == fill_value);
				if (any(ind))
					values = single(values);		% Here I give up and convert to singles
					values(ind) = NaN;
				end
			else
				% An idiotic case. A _FillValue = NAN in a array of integers.
			end
        case { nc_double, nc_float }
			[fill_value, status] = mexnc( 'get_att_double', ncid, varid, 'missing_value' );
			if (~isnan(fill_value))
				values(values == fill_value) = NaN;
			end
        otherwise
			mexnc('close',ncid);
			error ( sprintf('unhandled datatype %d\n', mfilename, var_type) );
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
	if ( scale_factor ~= 1 && (~isa(values,'single') || ~isa(values,'double')) )
		values = single(values);
	end
end

if (have_addoffset)
    [add_offset, status] = mexnc ( 'get_att_double', ncid, varid, 'add_offset' );
    if ( status ~= 0 )
        mexnc('close',ncid);
        snc_error ( 'NC_FUNS:NC_VARGET:MEXNC:GET_ATT_DOUBLE', mexnc('strerror', status) );
    end
end

if (have_scale && have_addoffset)
	cvlib_mex('CvtScale',values, scale_factor, add_offset)
elseif (have_scale)
	cvlib_mex('CvtScale',values, scale_factor)
elseif (have_addoffset)
	cvlib_mex('addS',values, add_offset)
end
	
%values = values * scale_factor + add_offset;

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
	else				fprintf(1, str)
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
	error ( sprintf('%s:  mexnc:inq_attname failed on varid %d.\n', mfilename, varid) );
end
attribute.Name = attname;

[att_datatype, status] = mexnc('INQ_ATTTYPE', cdfid, varid, attname);
if status < 0 
	error ( sprintf('%s:  mexnc:inq_att failed on varid %d, attribute %s.\n', mfilename, varid, attname) );
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
		error ( sprintf('att_datatype is %d.\n', att_datatype) );
end
if status < 0 
	error ( sprintf('%s:  mexnc:attget failed on varid %d, attribute %s.\n', mfilename, varid, attname) );
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

[dud,var_type,nvdims,var_dim,dud, status]=mexnc('INQ_VAR',ncid,varid);
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
elseif nargin == 3
    write_op = 'put_var';
elseif nargin == 5
    write_op = 'put_vara';
elseif nargin == 6
    write_op = 'put_vars';
else
    error ( 'unhandled write op.  How did we come to this??\n' );
end

validate_input_size_vs_netcdf_size(ncid,data,nc_count,count,write_op);

data = handle_fill_value ( ncid, varid, data );
data = handle_scaling(ncid,varid,data);

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
        mexnc ( 'close', ncid );
        efmt = 'Rank of input data (%d) does not work for netCDF variable %s.\n';
        snc_error ( 'NC_FUNS:NC_VARPUT:badRank', sprintf(efmt, ndims(data), varname) );
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
function data = handle_scaling(ncid,varid,data)
% HANDLE_MEX_SCALING
%     If there is a scale factor and/or  add_offset attribute, convert the data
%     to double precision and apply the scaling.
%

[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'scale_factor' );
if ( status == 0 )
    have_scale_factor = 1;
else
    have_scale_factor = 0;
end
[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'add_offset' );
if ( status == 0 )
    have_add_offset = 1;
else
    have_add_offset = 0;
end

% Return early if we don't have either one.
if ~(have_scale_factor || have_add_offset)
    return;
end

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

[var_type,status]=mexnc('INQ_VARTYPE',ncid,varid);
if status ~= 0 
    mexnc ( 'close', ncid );
    snc_error ( 'NC_FUNS:NC_VARPUT:MEXNC:INQ_VARTYPE', mexnc('STRERROR', status) );
end

data = (double(data) - add_offset) / scale_factor;

% When scaling to an integer, we should add 0.5 to the data.  Otherwise
% there is a tiny loss in precision, e.g. 82.7 should round to 83, not .
switch var_type
	case { nc_int, nc_short, nc_byte, nc_char }
	    data = round(data);
end

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
			error ('NC_FUNS:NC_VARPUT:unhandledDatatype', sprintf('Unhandled datatype for fill value, ''%s''.', class(data)) );
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
% We are now in a position to do a check on the input var size vs. 
% the known size in the netcdf file.  If we don't do this, it would 
% be possible to send a larger then expected chunk of data to the 
% netcdf file, have parts of it get lopped off in order to fit, and 
% never be the wiser.
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
	    mexnc ( 'close', ncid );
        fmt = 'Total number of input datums was %d, but the netcdf variable size is %d elements.';
	    msg = sprintf ( fmt, numel(data), prod(nc_count) );
	    snc_error ( 'NC_FUNS:NC_VARPUT:badInput', msg );
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
switch ( write_op )

    case 'put_var1'
        switch ( class(data) ),
        case 'double',		funcstr = 'put_var1_double';
        case 'single',		funcstr = 'put_var1_float';
        case 'int32',		funcstr = 'put_var1_int';
        case 'int16',		funcstr = 'put_var1_short';
        case 'int8',		funcstr = 'put_var1_schar';
        case 'uint8',		funcstr = 'put_var1_uchar';
        case 'char',		funcstr = 'put_var1_text';
        otherwise
            mexnc('close',ncid);
            snc_error ( 'NC_FUNS:NC_VARPUT:unhandledMatlabType', sprintf('unhandled data class %s\n', class(pdata)) );
        end
        status = mexnc (funcstr, ncid, varid, start, pdata );

    case 'put_var'
        switch ( class(data) )
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
        switch ( class(data) )
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
        status = mexnc (funcstr, ncid, varid, start, count, stride, pdata );

    otherwise 
        mexnc ( 'close', ncid );
        snc_error('NC_FUNS:NC_VARPUT:unhandledWriteOp', sprintf('unknown write operation''%s''.\n', write_op) );

end

if ( status ~= 0 )
    mexnc ( 'close', ncid );
    error ( sprintf('write operation ''%s'' failed with error ''%s''.\n', write_op, mexnc('STRERROR', status) ) );
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
	snc_error ( 'NC_FUNS:NC_ADDVAR:badInput', 'file argument must be character' );
end

if ( ~isstruct(varstruct) )
	snc_error ( 'NC_FUNS:NC_ADDVAR:badInput', '2nd argument must be a structure' );
end

varstruct = validate_varstruct ( varstruct );

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 )
	msg = sprintf ( 'OPEN failed on %s, ''%s''', ncfile, mexnc('STRERROR', status) );
	snc_error ( 'NC_FUNS:NC_ADDVAR:MEXNC:OPEN', msg );
end

% determine the dimids of the named dimensions
num_dims = length(varstruct.Dimension);
dimids = zeros(num_dims,1);
for j = 1:num_dims
	[dimids(j), status] = mexnc ( 'dimid', ncid, varstruct.Dimension{j} );
	if ( status ~= 0 )
		mexnc ( 'close', ncid );
		snc_error ( 'NC_FUNS:NC_ADDVAR:MEXNC:DIMID', mexnc('STRERROR', status) );
	end
end

% go into define mode
status = mexnc ( 'redef', ncid );
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	snc_error ( 'NC_FUNS:NC_ADDVAR:MEXNC:REDEF', mexnc('strerror', status) );
end

[varid, status] = mexnc('DEF_VAR', ncid, varstruct.Name, varstruct.Nctype, num_dims, dimids );
if ( status ~= 0 )
	mexnc ( 'endef', ncid );
	mexnc ( 'close', ncid );
	snc_error ( 'NC_FUNS:NC_ADDVAR:MEXNC:DEF_VAR', mexnc('STRERROR', status) );
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

% Now just use nc_attput to put in the attributes
for j = 1:length(varstruct.Attribute)
	attname = varstruct.Attribute(j).Name;
	attval = varstruct.Attribute(j).Value;
	nc_attput ( ncfile, varstruct.Name, attname, attval );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varstruct = validate_varstruct ( varstruct )
% Check that required fields are there.
% Must at least have a name.
if ~isfield ( varstruct, 'Name' )
	snc_error( 'NC_FUNS:NC_ADDVAR:badInput', 'structure argument must have at least the ''Name'' field.' );
end

% Check that required fields are there. Default Nctype is double.
if (~isfield ( varstruct, 'Nctype' )),	varstruct.Nctype = 'double';	end

% Are there any unrecognized fields?
fnames = fieldnames ( varstruct );
for j = 1:length(fnames)
	fname = fnames{j};
	switch ( fname )

	case { 'Nctype', 'Name', 'Dimension', 'Attribute' }
		% These are used to create the variable.  They are ok.
		
	case { 'Unlimited', 'Size', 'Rank' }
		% These come from the output of nc_getvarinfo.  We don't 
		% use them, but let's not give the user a warning about
		% them either.

	otherwise
		fprintf ( 2, '%s:  unrecognized field name ''%s''.  Ignoring it...\n', mfilename, fname );
	end
end

% If the datatype is not a string.
% Change suggested by Brian Powell
if ( isa(varstruct.Nctype, 'double') && varstruct.Nctype < 7 )
	types={ 'byte' 'char' 'short' 'int' 'float' 'double'};
	varstruct.Nctype = char(types(varstruct.Nctype));
end

% Check that the datatype is known.
switch ( varstruct.Nctype )
case { 'NC_DOUBLE', 'double', ...
	'NC_FLOAT', 'float', ...
	'NC_INT', 'int', ...
	'NC_SHORT', 'short', ...
	'NC_BYTE', 'byte', ...
	'NC_CHAR', 'char'  }
	% Do nothing
otherwise
	snc_error( 'NC_FUNS:NC_ADDVAR:unknownDatatype', sprintf('unknown type ''%s''\n', mfilename, varstruct.Nctype) );
end

% Check that required fields are there. Default Dimension is none.  Singleton scalar.
if (~isfield ( varstruct, 'Dimension' )),	varstruct.Dimension = [];	end

% Check that required fields are there. Default Attributes are none
if (~isfield ( varstruct, 'Attribute' )),	varstruct.Attribute = [];	end

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

[ncid, status] =mexnc( 'open', ncfile, nc_write_mode );
if  status ~= 0 
	snc_error ( 'NC_FUNS:NC_ATTGET:MEXNC:badFile', mexnc('STRERROR', status) );
end

% Put into define mode.
status = mexnc ( 'redef', ncid );
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	snc_error ( 'NC_FUNS:NC_ATTGET:MEXNC:REDEF', mexnc('STRERROR', status) );
end

if isnumeric(varname)
	varid = varname;
else
	[varid, status] = mexnc ( 'inq_varid', ncid, varname );
	if ( status ~= 0 )
		mexnc ( 'close', ncid );
		snc_error ( 'NC_FUNS:NC_ATTGET:MEXNC:INQ_VARID', mexnc('STRERROR', status) );
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
		snc_error ( 'NC_FUNS:NC_ATTGET:unhandleDatatype', msg );
end

status = mexnc ( funcstr, ncid, varid, attribute_name, atttype, length(attval), attval);
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	snc_error ( ['NC_FUNS:NC_ATTGET:MEXNC:' upper(funcstr)], mexnc('STRERROR', status) );
end

% End define mode.
status = mexnc ( 'enddef', ncid );
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	snc_error ( 'NC_FUNS:NC_ATTGET:MEXNC:ENDDEF', mexnc('STRERROR', status) );
end

status = mexnc('close',ncid);
if ( status ~= 0 )
	snc_error ( 'NC_FUNS:NC_ATTGET:MEXNC:CLOSE', mexnc('STRERROR', status) );
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
	mexnc ( 'close', ncid );
	snc_error ( 'NC_FUNS:NC_ADD_DIMENSION:enddefFailed', mexnc('STRERROR', status) );
end

status = mexnc ( 'close', ncid );
if status 
	snc_error ( 'NC_FUNS:NC_ADD_DIMENSION:closeFailed', mexnc('STRERROR', status) );
end

% ---------------------------------------------------------------------------------------
function nc_create_empty ( ncfile, mode )
% NC_CREATE_EMPTY:  creates an empty netCDF file 
%     NC_CREATE_EMPTY(NCFILE,MODE) creates the empty netCDF file NCFILE
%     with the given MODE.  MODE is optional, defaulting to nc_clobber_mode.
%
% EXAMPLE:
%     Suppose you wish to create a netCDF file with large file support.  
%     This would do it.
%
%     >> my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
%     >> nc_create_empty ( 'test.nc', my_mode );
%
% SEE ALSO:  
%     nc_noclobber_mode, nc_clobber_mode, nc_64bit_offset_mode

snc_nargchk(1,2,nargin);

% Set the default mode if necessary.
if nargin == 1
	mode = nc_clobber_mode;
end

[ncid, status] = mexnc ( 'CREATE', ncfile, mode );
if ( status ~= 0 )
	snc_error ( 'NC_FUNS:NC_CREATE_EMPTY:MEXNC:CREATE', mexnc('STRERROR', status) );
end

status = mexnc ( 'CLOSE', ncid );
if ( status ~= 0 )
	snc_error ( 'NC_FUNS:NC_CREATE:MEXNC:CLOSE', mexnc('STRERROR', status) );
end








% --------------------------------------------------------------------
function snc_error ( error_id, error_msg )
	error ( error_id, error_msg );

% --------------------------------------------------------------------
function snc_nargchk(low,high,N)
% SNC_NARGCHK:  wrapper for NARGCHK, which changed functionality at R???
	error ( nargchk(low,high,N ) );

% --------------------------------------------------------------------
function snc_nargoutchk(low,high,N)
% SNC_NARGOUTCHK:  wrapper for NARGOUTCHK, which changed functionality at R???
	error ( nargoutchk(low,high,N ) );

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
    error ( 'NC_FUNS:NC_VARGET:badCountIndex', sprintf ( fmt, length(count), nvdims, varname ) );
end
if ~isempty(stride) && (length(stride) ~= nvdims)
    mexnc('close',ncid);
    fmt = 'length of the stride index vector (%d) does not equal the number of dimensions (%d) for %s';
    error ( 'NC_FUNS:NC_VARGET:badStrideIndex', sprintf ( fmt, length(count), nvdims, varname ) );
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
    for j=1:nvdims,
        dimid = dimids(j);
        [dim_size,status]=mexnc('inq_dimlen', ncid, dimid);
        if ( status ~= 0 )
			mexnc('close',ncid);
            error ( 'NC_FUNS:NC_VARGET:MEXNC:INQ_DIM_LEN', mexnc('strerror', status) );
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
	
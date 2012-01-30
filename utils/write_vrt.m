function [name_vrt, comp_type] = write_vrt(full_name, opt, names, varargin)
% Write an GDAL .vrt header companion of the raw file FULL_NAME
%
% If FULL_NAME is a two elements cell array, first elements contains the VRT full file name
% and the second one the full filename (possibly a URL) of the pointed to file name.
%
% OPT contains the type of file we will deal with
%   Predefined values are 'SRTM30' or 'SRTM1'
%   Optionaly it may have [x_min y_max n_cols n_rows x_inc y_inc]
%
% NAMES	When writing a composit VRT holds the names if the individual VRTs. Otherwise empty.
%
% VARARGIN - See the params defaults below (to be completed)
%
% It returns the header full name (e.g with absolute path)
% It also tests if file is zip or gzip compressed. In case it is, the
% compression type is returned in COMP_TYPE. Otherwise it is empty.

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

	params.raw = false;
	params.band = 1;
	params.BYTEORDER = 'LSB';
	params.PixelOffset = 2;
	params.relative2VRT = 0;
	params.type = 'Int16';
	params.source = 'simple';
	params.nodata = false;
	params.interleave = 'BSQ';
	params.ImageOffset = 0;
	params = parse_pv_pairs(params, varargin);

	only_one = true;
	bare = true;
	comp_type = [];

	if (~isempty(names)),	only_one = false;	end

	if (ischar(opt))
		opt = upper(opt);
	else
		x_min = opt(1);		y_max = opt(2);
		n_cols = opt(3);    n_rows = opt(4);
		x_inc = opt(5);     y_inc = opt(6);
	end

	% Check if VTR and pointed file live in different places
	if (isa(full_name,'cell'))
		pointed_file = full_name{2};
		if (numel(full_name) == 3)
			simple = strncmpi(full_name{3},'simple',3);
			bare = false;
		end
		full_name = full_name{1};
	else
		pointed_file = full_name;	% The live together
	end

	% Test if we have a compressed file
	[PATH,fname,EXT] = fileparts(full_name);
	if (~isempty(PATH)),	PATH = [PATH filesep];	end		% This way we won't have a leading '/' on non-path names
	name_copy = fname;
	fname = lower(name_copy);       % Turn it to lower so that we can search for both W & w, N & n, etc ...
	if (strcmpi(EXT,'.zip')),		comp_type = 'zip';
	elseif (strcmpi(EXT,'.gz'))		comp_type = 'gzip';
	end

	if (~isempty(comp_type))		% File is compressed. Need to remove the true extension
		[pato,name_copy] = fileparts(name_copy);
		[pato,fname] = fileparts(fname);
	end

	if (ischar(opt))				% SRTM1|3|30
		[x_min y_max n_cols n_rows x_inc y_inc params.nodata] = known_cases(fname, opt);
		params.BYTEORDER = 'MSB';
	end
	y_inc = -abs(y_inc);			% GDAL wants it like that
	
	params.nXSize = n_cols;			params.nYSize = n_rows;

	name_vrt = [PATH name_copy '.vrt'];		% The extension has already been removed above
	fid = fopen(name_vrt,'wt');

	% Calculate the <dims> that go into first line
	if (only_one)
		n_totalRows = n_rows;		n_totalCols = n_cols;
	else
		n_totalRows = (max([names{:,2}]) + 1) * n_rows;		% Total size of the master virtual file
		n_totalCols = (max([names{:,3}]) + 1) * n_cols;
	end
	
	fprintf(fid,'<VRTDataset rasterXSize="%d" rasterYSize="%d">;\n', n_totalCols, n_totalRows);
	if (~params.raw)
		fprintf(fid, '%s\n', SRS_block);
		fprintf(fid,'  <GeoTransform>%.18g, %.19g,  0.0,  %.18g,  0.0, %.19g</GeoTransform>\n', x_min,x_inc,y_max, y_inc);

		if (only_one)			% A VRT for a single file
			if (bare)
				RasterBand_bare(fid, pointed_file, params)
			elseif (simple)
				RasterBand_simple(fid, pointed_file, params)
			end
		else					% A master VRT with individual VRTS as 'childrens'
			RasterBand_complex(fid, names, params)
		end
	else
		RasterBand_raw(fid, pointed_file, params)
	end

	fclose(fid);

% ----------------------------------------------------------------------------------
function txt = SRS_block(opt)
% Create a SRS block. OPT will be used to decide for different SRS when need comes.
	txt = ['  <SRS>GEOGCS[&quot;WGS 84&quot;,DATUM[&quot;WGS_1984&quot;,SPHEROID[&quot;WGS 84&quot;,' ...
		'6378137,298.257223563,AUTHORITY[&quot;EPSG&quot;,&quot;7030&quot;]],AUTHORITY[&quot;EPSG&quot;,' ...
		'&quot;6326&quot;]],PRIMEM[&quot;Greenwich&quot;,0,AUTHORITY[&quot;EPSG&quot;,&quot;8901&quot;]],'...
		'UNIT[&quot;degree&quot;,0.0174532925199433,AUTHORITY[&quot;EPSG&quot;,&quot;9122&quot;]],'...
		'AUTHORITY[&quot;EPSG&quot;,&quot;4326&quot;]]</SRS>'];

% ----------------------------------------------------------------------------------
function RasterBand_bare(fid, name, params)
% Create a RasterBand block. Last 3 args are optional under the condition on next lines
	fprintf(fid,'  <VRTRasterBand dataType="%s" band="%d" subClass="VRTRawRasterBand">\n', params.type, params.band(1));
	fprintf(fid,'    <SourceFilename relativetoVRT="%d">%s</SourceFilename>\n', params.relative2VRT, name);
	fprintf(fid,'    <ImageOffset>%d</ImageOffset>\n', params.ImageOffset);
	fprintf(fid,'    <PixelOffset>%d</PixelOffset>\n', params.PixelOffset);
	fprintf(fid,'    <LineOffset>%d</LineOffset>\n', params.PixelOffset * params.nXSize);
	fprintf(fid,'    <ByteOrder>%s</ByteOrder>\n', params.BYTEORDER);
	fprintf(fid,'    <NoDataValue>%f</NoDataValue>\n', params.nodata);
	fprintf(fid,'  </VRTRasterBand>\n');
	fprintf(fid,'</VRTDataset>\n');

% ----------------------------------------------------------------------------------
function RasterBand_simple(fid, name, params)
% Create a RasterBand block. Last 3 args are optional under the condition on next lines
	fprintf(fid,'  <VRTRasterBand dataType="%s" band="%d">\n', params.type, params.band(1));
	fprintf(fid,'    <NoDataValue>%f</NoDataValue>\n', params.nodata);
	fprintf(fid,'    <SimpleSource>\n');
	fprintf(fid,'      <SourceFilename relativetoVRT="%d">%s</SourceFilename>\n', params.relative2VRT, name);
	fprintf(fid,'      <SourceBand>%d</SourceBand>\n', params.band(1));
	fprintf(fid,'      <SourceProperties RasterXSize="%d" RasterYSize="%d" DataType="%s" BlockXSize="%d" BlockYSize="1"/>\n', ...
		params.nXSize, params.nXSize, params.type, params.nXSize);
	fprintf(fid,'      <SrcRect xOff="0" yOff="0" xSize="%d" ySize="%d"/>\n', params.nXSize, params.nYSize);
	fprintf(fid,'      <DstRect xOff="0" yOff="0" xSize="%d" ySize="%d"/>\n', params.nXSize, params.nYSize);
	fprintf(fid,'    </SimpleSource>\n');
	fprintf(fid,'  </VRTRasterBand>\n');
	fprintf(fid,'</VRTDataset>\n');

% ----------------------------------------------------------------------------------
function RasterBand_complex(fid, names, params)
% Create a RasterBand block.

	fprintf(fid,'  <VRTRasterBand dataType="%s" band="%d">\n', params.type, params.band(1));
	fprintf(fid,'    <NoDataValue>%f</NoDataValue>\n', params.nodata);
	for (k = 1:size(names,1))
		fprintf(fid,'    <ComplexSource>\n');
		fprintf(fid,'      <SourceFilename relativetoVRT="%d">%s</SourceFilename>\n', params.relative2VRT, names{k,1});
		fprintf(fid,'      <SourceBand>%d</SourceBand>\n', params.band(1));
		fprintf(fid,'      <SourceProperties RasterXSize="%d" RasterYSize="%d" DataType="%s" BlockXSize="1800" BlockYSize="1"/>\n', ...
			params.nXSize, params.nYSize, params.type);
		fprintf(fid,'      <SrcRect xOff="0" yOff="0" xSize="%d" ySize="%d"/>\n', params.nXSize, params.nYSize);
		fprintf(fid,'      <DstRect xOff="%d" yOff="%d" xSize="%d" ySize="%d"/>\n', ...
			names{k,3}*params.nXSize, names{k,2}*params.nYSize, params.nXSize, params.nYSize);
		fprintf(fid,'      <NODATA>%f</NODATA>\n', params.nodata);
		fprintf(fid,'    </ComplexSource>\n');
	end
	fprintf(fid,'  </VRTRasterBand>\n');
	fprintf(fid,'</VRTDataset>\n');

% ----------------------------------------------------------------------------------
function RasterBand_raw(fid, name, params)
% ... params.nTotalBands
	for (k = 1:numel(params.band))
		fprintf(fid,'  <VRTRasterBand dataType="%s" band="%d" subClass="VRTRawRasterBand">\n', bands.type, params.band(k));
		fprintf(fid,'    <SourceFilename relativetoVRT="%d">%s</SourceFilename>\n', params.relative2VRT, name);
		if (strcmp(params.interleave,'BSQ'))
			fprintf(fid,'    <ImageOffset>%d</ImageOffset>\n', params.ImageOffset + (k-1) * params.nXSize * params.nYSize);
			fprintf(fid,'    <PixelOffset>%d</PixelOffset>\n', params.PixelOffset);
			fprintf(fid,'    <LineOffset>%d</LineOffset>\n',   params.nXSize);
		elseif (strcmp(params.interleave,'BIL'))
			fprintf(fid,'    <ImageOffset>%d</ImageOffset>\n', params.ImageOffset + (k-1) * params.nXSize);
			fprintf(fid,'    <PixelOffset>%d</PixelOffset>\n', params.PixelOffset);
			fprintf(fid,'    <LineOffset>%d</LineOffset>\n',   params.nXSize * params.nTotalBands);
		else			% BIP and RGB or RGBA only
			fprintf(fid,'    <ImageOffset>%d</ImageOffset>\n', params.ImageOffset + k - 1);
			fprintf(fid,'    <PixelOffset>%d</PixelOffset>\n', numel(params.band));
			fprintf(fid,'    <LineOffset>%d</LineOffset>\n',   params.nXSize * numel(params.band));
		end
		fprintf(fid,'    <ByteOrder>%s</ByteOrder>\n', params.BYTEORDER);
		fprintf(fid,'  </VRTRasterBand>\n');
	end
	fprintf(fid,'</VRTDataset>\n');

% ---------------------------------------------
function [x_min, y_max, n_cols, n_rows, x_inc, y_inc, nodata] = known_cases(fname, tipo)
% ...
	lon_sng = 0;    lat_sng = 0;
	% Find the tile coordinates from the file name
	if (strcmp(tipo,'SRTM30'))
		x_w = strfind(fname(1:7),'w');
		x_e = strfind(fname(1:7),'e');     % fname(1:7), we only want to shearch here.
		y_s = strfind(fname(1:7),'s');
		y_n = strfind(fname(1:7),'n');
		if ~isempty(x_w),		ind_x = x_w(1);		lon_sng = -1;
		elseif  ~isempty(x_e),	ind_x = x_e(1);		lon_sng = 1;
		end

		if ~isempty(y_n),		ind_y = y_n(1);		lat_sng = 1;
		elseif ~isempty(y_s),	ind_y = y_s(1);		lat_sng = -1;
		end

		x_inc = 30 / 3600;      y_inc = x_inc;
		x_min = str2double(fname(ind_x+1:ind_x+3)) * lon_sng;
		y_max = str2double(fname(ind_y+1:ind_y+2)) * lat_sng;
		nodata = -32768;

		if (y_max > -60)      % Lower row tiles have different size
			n_rows = 6000;      n_cols = 4800;
		else
			n_rows = 3600;      n_cols = 7200;
		end
		x_min = x_min + x_inc / 2;      % Remember that those are pixel registered grids.
		y_max = y_max - y_inc / 2;
	elseif (strcmp(tipo,'SRTM1') || strcmp(tipo,'SRTM3'))
		x_w = strfind(fname(1:7),'w');
		x_e = strfind(fname(1:7),'e');
		y_s = strfind(fname(1:7),'s');
		y_n = strfind(fname(1:7),'n');
		if ~isempty(x_w),       ind_x = x_w(1);     lon_sng = -1;
		elseif  ~isempty(x_e),  ind_x = x_e(1);     lon_sng = 1;
		end

		if ~isempty(y_n),        lat_sng = 1;
		elseif ~isempty(y_s),   lat_sng = -1;    end
		x_min = str2double(fname(ind_x+1:ind_x+3)) * lon_sng;
		y_max = str2double(fname(2:ind_x-1)) * lat_sng + 1;
		nodata = -32768;
		if (strcmp(tipo,'SRTM1'))
			n_rows = 3601;      n_cols = 3601;      x_inc = 1 / 3600;
		else
			n_rows = 1201;      n_cols = 1201;      x_inc = 3 / 3600;
		end
		y_inc = x_inc;
	end

% ----------------------------------------------------------------------------
function params = parse_pv_pairs(params, pv_pairs)
% parse_pv_pairs: parses sets of property value pairs, allows defaults
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.
%
% Example usage:
% First, set default values for the parameters. Assume we have four
% parameters that we wish to use optionally in the function examplefun.
%
%  - 'viscosity', which will have a default value of 1
%  - 'volume', which will default to 1
%  - 'pie' - which will have default value 3.141592653589793
%  - 'description' - a text field, left empty by default
%
% The first argument to examplefun is one which will always be supplied.
%
%   function examplefun(dummyarg1,varargin)
%   params.Viscosity = 1;
%   params.Volume = 1;
%   params.Pie = 3.141592653589793
%
%   params.Description = '';
%   params=parse_pv_pairs(params,varargin);
%   params
%
% Use examplefun, overriding the defaults for 'pie', 'viscosity'
% and 'description'. The 'volume' parameter is left at its default.
%
%   examplefun(rand(10),'vis',10,'pie',3,'Description','Hello world')
%
% params = 
%     Viscosity: 10
%        Volume: 1
%           Pie: 3
%   Description: 'Hello world'
%
% Note that capitalization was ignored, and the property 'viscosity' was truncated
% as supplied. Also note that the order the pairs were supplied was arbitrary.

n = length(pv_pairs) / 2;

if (n == 0),	return,		end     % just return the defaults

if n ~= floor(n)
    error 'Property/value pairs must come in PAIRS.'
end

if ~isstruct(params)
    error 'No structure for defaults was supplied'
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for (i = 1:n)
	p_i = lower(pv_pairs{2*i-1});
	v_i = pv_pairs{2*i};
	
	ind = strcmp(p_i, lpropnames);
    if (~any(ind))
	    ind = find(strncmp(p_i, lpropnames, length(p_i)));
        if isempty(ind)
            error(['No matching property found for: ', pv_pairs{2*i-1}])
	    elseif (length(ind) > 1)
            error(['Ambiguous property name: ', pv_pairs{2*i-1}])
        end
    end
    p_i = propnames{ind};
	%params = setfield(params,p_i,v_i);      % override the corresponding default in params
	params.(p_i) = v_i;
end

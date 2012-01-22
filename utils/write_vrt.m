function [name_vrt, comp_type] = write_vrt(full_name, opt, rel, names)
% Write an GDAL .vrt header companion of the raw file FULL_NAME
%
% If FULL_NAME is a two elements cell array, first elements contains the VRT full file name
% and the second one the full filename (possibly a URL) of the pointed to file name.
%
% OPT contains the type of file we will deal with
%   Predefined values are 'SRTM30' or 'SRTM1'
%   Optionaly it may have [lon_min lat_max n_cols n_rows x_inc y_inc nodata_value]
%   OR [lon_min lat_max n_cols n_rows x_inc y_inc nodata_value byte_order n_bytes]
%       Where BYTE_ORDER is = 1 to indicate 'Intel' or = 0 'Motorola' byte order
%       N_BYTES is what it says (per grid value)
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

	BYTEORDER = 'MSB';      % Make these defaults to work with SRTM cases
	pixOff = 2;
	only_one = true;
	bare = true;
	relative = 0;
	comp_type = [];
	
	if (nargin == 4),	only_one = false;	end
	if (nargin >= 3),	relative = rel;		end

	if (ischar(opt))
		opt = upper(opt);
	else
		lon_min = opt(1);   lat_max = opt(2);
		n_cols = opt(3);    n_rows = opt(4);
		x_inc = opt(5);     y_inc = opt(6);
		nodata = opt(7);
		if (numel(opt) == 9)
			if (opt(8)),    BYTEORDER = 'LSB';
			else            BYTEORDER = 'MSB';
			end
			pixOff = opt(9);
		end
	end

	% Check if VTR and pointed file live in different places
	if (isa(full_name,'cell'))
		pointed_file = full_name{2};
		if (numel(full_name) == 3)
			simple = strncmpi(full_name{3},'simple',3);
			complx = strncmpi(full_name{3},'complex',3);
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
		[lon_min lat_max n_cols n_rows x_inc y_inc nodata] = known_cases(fname, opt);
	end
	y_inc = -abs(y_inc);			% GDAL wants it like that

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
	fprintf(fid, '%s\n', SRS_block);
	fprintf(fid,'  <GeoTransform>%.18g, %.19g,  0.0,  %.18g,  0.0, %.19g</GeoTransform>\n', lon_min,x_inc,lat_max, y_inc);

	if (only_one)			% A VRT for a single file
		if (bare)
			RasterBand_bare(fid, pointed_file, 0, pixOff, n_cols, BYTEORDER, nodata, 'Int16', 1, relative)
		elseif (simple)
			RasterBand_simple(fid, pointed_file, n_rows, n_cols, BYTEORDER, nodata, 'Int16', 1, relative)
		end
	else					% A master VRT with individual VRTS as 'childrens'
		RasterBand_complex(fid, names, n_rows, n_cols, BYTEORDER, nodata, 'Int16', 1, 1)
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

% --------------------------1-----2------3-------4-------5------6--------7------8----9-------10---
function RasterBand_bare(fid, name, imgOff, pixOff, n_cols, endian, nodata, tipo, band, relative)
% Create a RasterBand block. Last 3 args are optional under the condition on next lines
	if (nargin <= 7)
		tipo = 'Int16';		% Default to this data type
		band = 1;			% Default band number
		relative = 0;		% relativetoVRT number
	elseif (nargin <= 8)
		band = 1;		relative = 0;
	elseif (nargin <= 9)
		relative = 0;
	end
	fprintf(fid,'  <VRTRasterBand dataType="%s" band="%d" subClass="VRTRawRasterBand">\n', tipo, band);
	fprintf(fid,'    <SourceFilename relativetoVRT="%d">%s</SourceFilename>\n', relative, name);
	fprintf(fid,'    <ImageOffset>%d</ImageOffset>\n', imgOff);
	fprintf(fid,'    <PixelOffset>%d</PixelOffset>\n', pixOff);
	fprintf(fid,'    <LineOffset>%d</LineOffset>\n', pixOff * n_cols);
	fprintf(fid,'    <ByteOrder>%s</ByteOrder>\n', endian);
	fprintf(fid,'    <NoDataValue>%f</NoDataValue>\n', nodata);
	fprintf(fid,'  </VRTRasterBand>\n');
	fprintf(fid,'</VRTDataset>\n');

% --------------------------1-----2------3-------4-------5------6--------7------8----9----
function RasterBand_simple(fid, name, n_rows, n_cols, endian, nodata, tipo, band, relative)
% Create a RasterBand block. Last 3 args are optional under the condition on next lines
	fprintf(fid,'  <VRTRasterBand dataType="%s" band="%d">\n', tipo, band);
	fprintf(fid,'    <NoDataValue>%f</NoDataValue>\n', nodata);
	fprintf(fid,'    <SimpleSource>\n');
	fprintf(fid,'      <SourceFilename relativetoVRT="%d">%s</SourceFilename>\n', relative, name);
	fprintf(fid,'      <SourceBand>%d</SourceBand>\n', band);
	fprintf(fid,'      <SourceProperties RasterXSize="%d" RasterYSize="%d" DataType="%s" BlockXSize="%d" BlockYSize="1"/>\n',n_cols,n_rows,tipo,n_cols);
	fprintf(fid,'      <SrcRect xOff="0" yOff="0" xSize="%d" ySize="%d"/>\n', n_cols, n_rows);
	fprintf(fid,'      <DstRect xOff="0" yOff="0" xSize="%d" ySize="%d"/>\n', n_cols, n_rows);
	fprintf(fid,'    </SimpleSource>\n');
	fprintf(fid,'  </VRTRasterBand>\n');
	fprintf(fid,'</VRTDataset>\n');

% ----------------------------------------------------------------------------------
function RasterBand_complex(fid, names, n_rows, n_cols, endian, nodata, tipo, band, relative)
% Create a RasterBand block. Last 3 args are optional under the condition on next lines

	fprintf(fid,'  <VRTRasterBand dataType="%s" band="%d">\n', tipo, band);
	fprintf(fid,'    <NoDataValue>%f</NoDataValue>\n', nodata);
	for (k = 1:size(names,1))
		fprintf(fid,'    <ComplexSource>\n');
		fprintf(fid,'      <SourceFilename relativetoVRT="%d">%s</SourceFilename>\n', relative, names{k,1});
		fprintf(fid,'      <SourceBand>%d</SourceBand>\n', band);
		fprintf(fid,'      <SourceProperties RasterXSize="%d" RasterYSize="%d" DataType="%s" BlockXSize="1800" BlockYSize="1"/>\n',n_cols,n_rows,tipo);
		fprintf(fid,'      <SrcRect xOff="0" yOff="0" xSize="%d" ySize="%d"/>\n', n_cols, n_rows);
		fprintf(fid,'      <DstRect xOff="%d" yOff="%d" xSize="%d" ySize="%d"/>\n', names{k,2}*n_cols, names{k,3}*n_rows, n_cols, n_rows);
		fprintf(fid,'      <NODATA>%f</NODATA>\n', nodata);
		fprintf(fid,'    </ComplexSource>\n');
	end
	fprintf(fid,'  </VRTRasterBand>\n');
	fprintf(fid,'</VRTDataset>\n');

% ---------------------------------------------
function [lon_min, lat_max, n_cols, n_rows, x_inc, y_inc, nodata] = known_cases(fname, tipo)
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
		lon_min = str2double(fname(ind_x+1:ind_x+3)) * lon_sng;
		lat_max = str2double(fname(ind_y+1:ind_y+2)) * lat_sng;
		nodata = -32768;

		if (lat_max > -60)      % Lower row tiles have different size
			n_rows = 6000;      n_cols = 4800;
		else
			n_rows = 3600;      n_cols = 7200;
		end
		lon_min = lon_min + x_inc / 2;      % Remember that those are pixel registered grids.
		lat_max = lat_max - y_inc / 2;
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
		lon_min = str2double(fname(ind_x+1:ind_x+3)) * lon_sng;
		lat_max = str2double(fname(2:ind_x-1)) * lat_sng + 1;
		nodata = -32768;
		if (strcmp(tipo,'SRTM1'))
			n_rows = 3601;      n_cols = 3601;      x_inc = 1 / 3600;
		else
			n_rows = 1201;      n_cols = 1201;      x_inc = 3 / 3600;
		end
		y_inc = x_inc;
	end

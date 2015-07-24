function [name_hdr,comp_type] = write_esri_hdr(full_name, opt)
% Write an ESRI .hdr header companion of the raw file FULL_NAME
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

BYTEORDER = 'BYTEORDER     M';      % Make these defaults to work with SRTM cases
to_bytes = 2;

if (ischar(opt))
    opt = upper(opt);
else
    lon_min = opt(1);   lat_max = opt(2);
    n_cols = opt(3);    n_rows = opt(4);
    x_inc = opt(5);     y_inc = opt(6);
    nodata = opt(7);
    if (length(opt) == 9)
        if (opt(8)),    BYTEORDER = 'BYTEORDER     I';
        else            BYTEORDER = 'BYTEORDER     M';      end
        to_bytes = opt(9);
    end
end

comp_type = [];

% Test if we have a compressed file
[PATH,fname,EXT] = fileparts(full_name);
name_copy = fname;
fname = lower(name_copy);       % Turn it to lower so that we can search for both W & w, N & n, etc ...
if (strcmpi(EXT,'.zip'))
    comp_type = 'zip';
elseif (strcmpi(EXT,'.gz'))
    comp_type = 'gzip';
end

if (~isempty(comp_type))     % File is compressed. Need to remove the true extension
    [pato,name_copy] = fileparts(name_copy);
    [pato,fname] = fileparts(fname);
end

lon_sng = 0;    lat_sng = 0;
% Find the tile coordinates from the file name
if (strcmp(opt,'SRTM30'))
	x_w = strfind(fname(1:7),'w');   x_e = strfind(fname(1:7),'e');     % fname(1:7), we only want to shearch here.
	y_s = strfind(fname(1:7),'s');   y_n = strfind(fname(1:7),'n');
    if ~isempty(x_w),       ind_x = x_w(1);     lon_sng = -1;
    elseif  ~isempty(x_e),  ind_x = x_e(1);     lon_sng = 1;
    end
	
    if ~isempty(y_n),       ind_y = y_n(1);     lat_sng = 1;
    elseif ~isempty(y_s),   ind_y = y_s(1);     lat_sng = -1;
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
elseif (strcmp(opt,'SRTM1') || strcmp(opt,'SRTM3'))
	x_w = strfind(fname(1:7),'w');   x_e = strfind(fname(1:7),'e');
	y_s = strfind(fname(1:7),'s');   y_n = strfind(fname(1:7),'n');
    if ~isempty(x_w),       ind_x = x_w(1);     lon_sng = -1;
    elseif  ~isempty(x_e),  ind_x = x_e(1);     lon_sng = 1;
    end
	
	if ~isempty(y_n),        lat_sng = 1;
	elseif ~isempty(y_s),   lat_sng = -1;    end
    lon_min = str2double(fname(ind_x+1:ind_x+3)) * lon_sng;
    lat_max = str2double(fname(2:ind_x-1)) * lat_sng + 1;
    nodata = -32768;
    if (strcmp(opt,'SRTM1'))
        n_rows = 3601;      n_cols = 3601;      x_inc = 1 / 3600;
    else
        n_rows = 1201;      n_cols = 1201;      x_inc = 3 / 3600;
    end
    y_inc = x_inc;
end

% Write the ESRI header
tmp = cell(15,1);
tmp{1} = BYTEORDER;
tmp{2} = 'LAYOUT        BIL';
tmp{3} = ['NROWS       ' num2str(n_rows)];
tmp{4} = ['NCOLS       ' num2str(n_cols)];
tmp{5} = 'NBANDS        1';
tmp{6} = ['NBITS         ' num2str(to_bytes*8)];
tmp{7} = 'PIXELTYPE  SIGNEDINT';
tmp{8} = ['BANDROWBYTES  ' num2str(n_cols*to_bytes)];
tmp{9} = ['TOTALROWBYTES ' num2str(n_cols*to_bytes)];
tmp{10} = 'BANDGAPBYTES  0';
tmp{11} = ['NODATA        ' num2str(nodata)];
tmp{12} = ['ULXMAP        ' num2str(lon_min,15)];
tmp{13} = ['ULYMAP        ' num2str(lat_max,15)];
tmp{14} = ['XDIM         ' num2str(x_inc,'%.15f')];
tmp{15} = ['YDIM         ' num2str(y_inc,'%.15f')];

name_hdr = [PATH filesep name_copy '.hdr'];      % The extension has already been removed above
fid = fopen(name_hdr,'wt');
for (i = 1:15),   fprintf(fid,'%s\n',tmp{i});   end
fclose(fid);

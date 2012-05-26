function [Z, X, Y, srsWKT, handles, att] = read_grid(handles, fullname, tipo, opt)
% Loads grid files that may contain DEMs or other grid (not images (byte)) types
%
% HANDLES	-> Normally, the Mirone's handles structure but can actually be any structure with these fields:
%			-	grdMaxSize, ForceInsitu
%			-	path_tmp	(only need if TIPO == GMT and file is GMT grid not in netCDF)
%			-  Can be empty ([]) when TIPO ~= GMT.
% FULLNAME	-> Like it says. The file's full name
% TIPO		-> 'GMT' read grids using the read_gmt_type_grids() function
%			   'MOLA_lbl' To read Mars MOLA .img with a .lbl header file (cannot be compressed)
%			   'MOLA' To read Mars MOLA .img with a .lbl header file *OR* a V3 PDS .img
%			   'OVR' Used only by the overview tool together with -P option
%			   'IN' Used to send in an internally computed grid, transmitted in OPT
%			-> 'whatever'	Let GDAL guess what to do (it means, any string)
% OPT		-> -R<...> or -P<...> options of gdalread OR the "att" attributes structure
%				got from att = gdalread(fname,'-M',...).
%			-> It can also hold a structure with fields 'X','Y','Z','head' & 'name' (optional).
%				We use this construct to take advantage to all testing/setting machinery of this function.
%				With all of that still note that OPT is optional

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

% $Id$

	if (nargin == 3)	opt = ' ';	end
	opt_I = ' ';	srsWKT = [];	att = [];	attVRT = [];	Z = [];		X = [];		Y = [];
	if (isa(fullname, 'cell') && numel(fullname) == 2 )
		fname = [fullname{1} fullname{2}];
	else
		fname = fullname;
	end
	if (isempty(handles))
		if (strcmp(tipo, 'GMT')),	error('read_grid: handles cannot be empty when reading a GMT type'),	end
		handles.ForceInsitu = false;
		handles.grdMaxSize = 1e15;
	end
	try
		if (handles.ForceInsitu),	opt_I = '-I';		end	% Use only in desperate cases.
	end
	handles.was_int16 = 0;			% To make sure that it wasnt left = 1 from a previous use.

	if (strncmp(tipo,'GMT',3))		% GMT_relatives - Reading is done by the read_gmt_type_grids function
		[handles, X, Y, Z, head, misc] = read_gmt_type_grids(handles, fname);
		if (isempty(X))		return,		end
		if (isfield(misc,'z_dim') && numel(misc.z_dim) == 3),	handles.nLayers = misc.z_dim(3);	end
		if (~isempty(misc) && ~isempty(misc.srsWKT) ),			srsWKT = misc.srsWKT;	end
	elseif (strncmpi(tipo,'IN',2))
		Z = opt.Z;
		X = opt.X;
		Y = opt.Y;
		head = opt.head;
		if (isfield(opt,'name'))
			fname = opt.name;
		end
	else
		grdMaxSize = handles.grdMaxSize;
		if (nargin >= 4 && isa(opt,'struct'))		% From FileOpenGeoTIFF_CB
			att = opt;
			opt = ' ';								% So that we can use it harmlesslly below
		else
			if (strncmp(tipo,'MOLA',4))			% This means that a .img MOLA file with a .lbl cannot be cmpressed
				fname_t = [fname(1:end-3) 'lbl'];
				if (strcmp(tipo,'MOLA_lbl') || exist(fname_t, 'file'))
					fname = fname_t;
					[limits, n_cols, n_rows, n_bytes, A_rad, B_rad, inc] = parse_MOLA_lbl(fname);
					fnameVRT = write_vrt(fname, [limits(1) limits(4) n_cols n_rows inc inc], [], ...
								'source', 'simple','BYTEORDER','MSB','PixelOffset',2);
					attVRT = gdalread(fnameVRT, '-M','-C');
					builtin('delete',fnameVRT);
				end
			else
				[PATH,FNAME,EXT] = fileparts(fname);
				if (strcmpi(EXT,'.zip'))
					fname = ['/vsizip/' fname filesep FNAME];
				elseif (strcmpi(EXT,'.gz'))
					fname = ['/vsigzip/' fname];
				end
			end
			try
				att = gdalread(fname,'-M','-C', opt);
				if (~isempty(attVRT))
					att.GeoTransform = attVRT.GeoTransform;		att.GMT_hdr = attVRT.GMT_hdr;
					att.GEOGCorners = attVRT.GEOGCorners;		att.Corners = attVRT.Corners;
					att.ProjectionRef = ogrproj(sprintf('+proj=longlat +a=%.3f +b=%.3f +no_defs',A_rad,B_rad));
				end
			catch
				errordlg(['GDALREAD: Unable to open file ' fname],'Error'),		return
			end
			if (strcmp(tipo,'OVR') && strncmp(opt,'-P',2))	% A call from the overview tool
				grdMaxSize = 1e15;							% Preview mode does not care of grdMaxSize
			end
		end
		if ((att.RasterXSize * att.RasterYSize * 4) > grdMaxSize)
			if ( strcmp(yes_or_no('title','Warning'),'Yes')),	return,		end		% Advise accepted
		end
		if (strcmp(att.Band(1).DataType,'Int16')),		handles.was_int16 = 1;	end

		if (~strncmp(att.DriverShortName, 'HDF4', 4))
			Z = gdalread(att.Name, '-U', opt_I, opt);
		elseif (strcmp(att.Band(1).DataType,'L3Bin'))
			[Z, handles.have_nans, att] = empilhador('getZ', fname, att, false, false, false, 1, 0, []);
		else								% HDF files need a special care. Search for an offset and scale factor, etc...
			[head, slope, intercept, base, is_modis, is_linear, is_log, att] = empilhador('getFromMETA', att);
			[Z, handles.have_nans, att] = empilhador('getZ', fname, att, is_modis, is_linear, is_log, slope, intercept, base);
			fname = att.fname;				% We'll need better but for now this ensures that no subdataset name is taken as the whole.
			if (~isa(Z,'int16')),		handles.was_int16 = 0;		end
		end
		[Z, did_scale] = handle_scaling(Z, att);	% See if we need to apply a scale/offset
		if (~did_scale),	handles.Nodata_int16 = att.Band(1).NoDataValue;		end
		handles.image_type = 4;
		srsWKT = att.ProjectionRef;
	end
	handles.fileName = fname;

	if (~isa(Z,'single')),		Z = single(Z);		end

	if ( ~strncmp(tipo,'GMT',3) && ~strncmpi(tipo,'IN',2) )
		if ( ~isempty(att.Band(1).NoDataValue) && ~isnan(att.Band(1).NoDataValue) && att.Band(1).NoDataValue ~= 0 )
			% Do this because some formats (e.g MOLA PDS v3) are so dumb that they declare a NoDataValue
			% and than don't use it !!!!!!
			if (att.Band(1).NoDataValue < 0)
				ind = (Z <= single(att.Band(1).NoDataValue));
			else
				ind = (Z >= single(att.Band(1).NoDataValue));
			end
			if (any(ind(:)))
				Z(ind) = NaN;		handles.have_nans = 1;
			end
			clear ind;
		end
		% GDAL wrongly reports the corners as [0 nx] [0 ny] when no SRS
		if ( isequal((att.Corners.LR - att.Corners.UL),[att.RasterXSize att.RasterYSize]) && ~all(att.Corners.UL) )
			att.GMT_hdr(1:4) = [1 att.RasterXSize 1 att.RasterYSize];
		end
		head = att.GMT_hdr;
		if (isequal(head(5:6),[0 0]) || any(isnan(head(5:6))))		% It can happen with GeoTiff_DEM
			zz = grdutils(Z,'-L');			head(5:6) = zz(1:2);		att.GMT_hdr(5:6) = head(5:6);
		end
		X = linspace(head(1),head(2),size(Z,2));	Y = linspace(head(3),head(4),size(Z,1));
	end
	
	handles.head = head;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, did_scale] = handle_scaling(data, att)
%	If there is a scale factor and/or  add_offset attribute, convert the data
%	to single precision and apply the scaling. 

	have_scale_factor = false;			have_add_offset = false;	did_scale = false;
	if ( att.Band(1).ScaleOffset(1) ~= 1 )
		have_scale_factor = true;
		scale_factor = att.Band(1).ScaleOffset(1);
	end
	if ( att.Band(1).ScaleOffset(2) ~= 0 )
		have_add_offset = true;
		add_offset = att.Band(1).ScaleOffset(2);
	end

	% Return early if we don't have either one.
	if ~(have_scale_factor || have_add_offset),		return,		end

	if (~isa(data,'single')),		data = single(data);	end

	did_scale = true;
	if (have_scale_factor && have_add_offset)
		data = cvlib_mex('CvtScale',data, scale_factor, add_offset);
	elseif (have_scale_factor)
		data = cvlib_mex('CvtScale',data, scale_factor);
	else
		data = cvlib_mex('addS',data, add_offset);
	end

% --------------------------------------------------------------------------------------
function [limits, n_cols, n_rows, n_bytes, A_rad, B_rad, inc] = parse_MOLA_lbl(fname)
% This a resseruction of my old way because I found GDAL is bugged in decoding
% those Martian PSD files and I'm afraid a fix will take forever to come.
% LIMITS	[x_min x_max y_min y_max] suposedly in grid-registration
	fid = fopen(fname,'rt');
	if (fid < 0)
		msg = ['ERROR: Could not open format descriptive file: ' fname];
		errordlg(msg,'Error');		error(msg);
	end
	s = strread(fread(fid,'*char').','%s','delimiter','\n');	fclose(fid);

	res = findcell('LINES', s);					ind = strfind(s{res.cn},'=');
	n_rows = str2double(s{res.cn}(ind(1)+1:end));

	s = s(res.cn+1:end);		% Strip the lines before this. They would only waste memory.

	res = findcell('LINE_SAMPLES', s);			ind = strfind(s{res.cn},'=');
	n_cols = str2double(s{res.cn}(ind(1)+1:end));

	res = findcell('A_AXIS_RADIUS', s);			ind = strfind(s{res.cn},'=');
	A_rad = str2double( strtok(s{res.cn}(ind(1)+1:end)) ) * 1000;

	res = findcell('B_AXIS_RADIUS', s);			ind = strfind(s{res.cn},'=');
	B_rad = str2double( strtok(s{res.cn}(ind(1)+1:end)) ) * 1000;

	res = findcell('SAMPLE_BITS', s);			ind = strfind(s{res.cn},'=');
	n_bytes = str2double(s{res.cn}(ind(1)+1:end)) / 8;

	res = findcell('CENTER_LATITUDE', s);		ind = strfind(s{res.cn},'=');
	lat0 = str2double(strtok(s{res.cn}(ind(1)+1:end)));

	res = findcell('CENTER_LONGITUDE', s);		ind = strfind(s{res.cn},'=');
	lon0 = str2double(strtok(s{res.cn}(ind(1)+1:end)));

	res = findcell('LINE_FIRST_PIXEL', s);		ind = strfind(s{res.cn},'=');
	line_first_pix = str2double(s{res.cn}(ind(1)+1:end));

	res = findcell('LINE_LAST_PIXEL', s);		ind = strfind(s{res.cn},'=');
	line_last_pix = str2double(s{res.cn}(ind(1)+1:end));

	res = findcell('SAMPLE_FIRST_PIXEL', s);	ind = strfind(s{res.cn},'=');
	sample_first_pix = str2double(s{res.cn}(ind(1)+1:end));

	res = findcell('SAMPLE_LAST_PIXEL', s);		ind = strfind(s{res.cn},'=');
	sample_last_pix = str2double(s{res.cn}(ind(1)+1:end));

	res = findcell('MAP_RESOLUTION', s);		ind = strfind(s{res.cn},'=');
	inc = 1 / str2double(strtok(s{res.cn}(ind(1)+1:end)));

	res = findcell('WESTERNMOST_LONGITUDE', s);	ind = strfind(s{res.cn},'=');
	W0 = str2double(strtok(s{res.cn}(ind(1)+1:end)));

	res = findcell('LINE_PROJECTION_OFFSET', s);	ind = strfind(s{res.cn},'=');
	line_off = str2double(s{res.cn}(ind(1)+1:end));

	res = findcell('SAMPLE_PROJECTION_OFFSET', s);	ind = strfind(s{res.cn},'=');
	sample_off = str2double(s{res.cn}(ind(1)+1:end));

	limits = [(sample_first_pix - sample_off)*inc+lon0 (sample_last_pix - sample_off)*inc+lon0 ...
		(line_first_pix - line_off)*inc+lat0 (line_last_pix - line_off)*inc+lat0];

	% Because I don't know if the PDS (lbl) files are always pix-reg as the cases I tested
	% I'm using this crazy schema to determine that. Couldn't they just report that on file?
	% What the hell are those LINE_PROJECTION_OFFSET meant for????
	if ( (abs(W0 - limits(1)) - inc) < 1e-5 )
		limits = limits + [-inc inc -inc inc] / 2;
	end

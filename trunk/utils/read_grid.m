function [Z, X, Y, srsWKT, handles, att] = read_grid(handles, fullname, tipo, opt)
% Loads grid files that may contain DEMs or other grid (not images (byte)) types
%
% HANDLES	-> Normally, the Mirone's handles structure but can actually be any structure with these fields:
%			-	grdMaxSize
%			-	path_tmp	(only need if TIPO == GMT and file is GMT grid not in netCDF)
% FULLNAME	-> Like it says. The file's full name
% TIPO		-> 'GMT' read grids using the read_gmt_type_grids() function
%			   'MOLA_lbl' To read Mars MOLA .img with a .lbl header file (cannot be compressed)
%			   'MOLA' To read Mars MOLA .img with a .lbl header file *OR* a V3 PDS .img
%			   'OVR' Used only by the overview tool together with -P option
%			-> 'whatever'	Let GDAL guess what to do (it means, any string)
% OPT		-> -R<...> or -P<...> options of gdalread OR the "att" attributes structure
%				got from att = gdalread(fname,'-M',...). This argument is optional

%	Copyright (c) 2004-2010 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

	if (nargin == 3)	opt = ' ';	end
	opt_I = ' ';	srsWKT = [];	att = [];	Z = [];		X = [];		Y = [];
	if (isa(fullname, 'cell') && numel(fullname) == 2 )
		fname = [fullname{1} fullname{2}];
	else
		fname = fullname;
	end
	if (handles.ForceInsitu),	opt_I = '-I';		end	% Use only in desperate cases.
	handles.was_int16 = 0;			% To make sure that it wasnt left = 1 from a previous use.

	if (strncmp(tipo,'GMT',3))		% GMT_relatives - Reading is done by the read_gmt_type_grids function
		[handles, X, Y, Z, head, misc] = read_gmt_type_grids(handles, fname);
		if (isempty(X))		return,		end
		if (~isempty(misc) && ~isempty(misc.srsWKT) ),		srsWKT = misc.srsWKT;	end
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
			catch
				errordlg(['GDALREAD: Unable to open file ' fname],'Error'),		return
			end
			if (strcmp(tipo,'OVR') && strncmp(opt,'-P',2))	% A call from the overview tool
				grdMaxSize = 1e15;						% Preview mode does not care of grdMaxSize
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

	if (~strncmp(tipo,'GMT',3))
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

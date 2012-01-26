function varargout = empilhador(varargin)
% Stacks a bunch of grids into a single 3D file
%
% WARNING: FOR COMPILING THIS WE NEED TO INCLUDE THE HDF_FUNS.M SRC
%
% NOTE: The gotFromMETA and getZ functions are called directly by mirone

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

	if (nargin > 1 && ischar(varargin{1}))
		gui_CB = str2func(varargin{1});
		[varargout{1:nargout}] = feval(gui_CB,varargin{2:end});
	else
		h = empilhador_OpeningFcn(varargin{:});
		if (nargout)    varargout{1} = h;   end
	end

% ---------------------------------------------------------------------------------
function hObject = empilhador_OpeningFcn(varargin)
	hObject = figure('Vis','off');
	empilhador_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	if (numel(varargin) > 0)
		handMir = varargin{1};
		handles.home_dir = handMir.home_dir;
		handles.last_dir = handMir.last_dir;
		handles.work_dir = handMir.work_dir;
		handles.IamCompiled = handMir.IamCompiled;		% Need to know due to crazy issue of nc_funs
        handles.path_tmp = handMir.path_tmp;
        d_path = handMir.path_data;
	else
		handles.home_dir = cd;
		handles.last_dir = handles.home_dir;
		handles.work_dir = handles.home_dir;
		handles.IamCompiled = false;
        handles.path_tmp = [pwd filesep 'tmp' filesep];
        d_path = [pwd filesep 'data' filesep];
	end
	handles.nameList = [];
	handles.OneByOneNameList = [];	% For when files are loaded one by one (risky)
	handles.OneByOneFirst = true;	% Safety valve to deal with the load one by one case
	handles.testedDS = false;		% To test if a Sub-Dataset request is idiot
	handles.Interactive = false;	% Used when need to reinterpolate L2 files. FALSE means do not
									% call the helper window that asks questions (use default ans)

	% -------------- Import/set icons --------------------------------------------
	load([d_path 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_namesList, 'CData',Mfopen_ico)

	%------------ Give a Pro look (3D) to the frame box -----------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) ------------------------------

	set(hObject,'Visible','on');
	guidata(hObject, handles);

% -----------------------------------------------------------------------------------------
function edit_namesList_CB(hObject, handles)
    fname = get(hObject,'String');
    push_namesList_CB([], handles, fname)

% -----------------------------------------------------------------------------------------
function push_namesList_CB(hObject, handles, opt)
    if (nargin == 2)        % Direct call
    	str1 = {'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.DAT,*.txt,*.TXT)';'*.*', 'All Files (*.*)'};
        [FileName,PathName,handles] = put_or_get_file(handles, str1,'File with grids list','get');
	    if isequal(FileName,0),		return,		end
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
	fname = [PathName FileName];

    [bin,n_column] = guess_file(fname);

    if isempty(bin)					% If error in reading file
		errordlg(['Error reading file ' fname],'Error'),	return
	elseif (bin)					% Binary file. Assume it's a target file, not a name list
		[PATH,FNAME,EXT] = fileparts(fname);
		handles.OneByOneNameList{end+1} = fname;		% Save the full name
		str = get(handles.listbox_list, 'Str');
		str{end+1} = [FNAME EXT];
		set(handles.listbox_list, 'Str', str);
		handles.OneByOneFirst = true;					% Repetitive but ensures that things are always updated
		guidata(handles.figure1, handles)
		return
    end

	fid = fopen(fname);
	c = fread(fid,'*char')';      fclose(fid);
	names = strread(c,'%s','delimiter','\n');   clear c fid;
	m = length(names);

	handles.strTimes = cell(m,1);		% To hold time steps as strings
	SDSinfo = cell(m,1);				% To hold Sub Datasets info
	handles.SDSinfo = [];				% If above exists, it will be copied here
	c = false(m,1);
	caracol = false(m,1);				% Case name list has '@' to pause for a CD change

	n_msg = 1;							% Will hold the "change DV messages" counter
	if (n_column > 1)					% When 2nd column holds the 3D numbering
		for (k = 1:m)
			[t,r] = strtok(names{k});
			if (t(1) == '#' || numel(t) < 2),	c(k) = true;	continue,	end		% Jump empty and comment lines
			if ( t(1) == '@')
				caracol(k) = true;
				if (~isempty(t(2:end))),	handles.changeCD_msg{n_msg} = t(2:end);		% The '@' was glued with the message
				else						handles.changeCD_msg{n_msg} = r;
				end
				n_msg = n_msg + 1;
				continue
			end

			names{k} = ddewhite(t);
			if (n_column == 2)			% Names & numeric label format
				r = ddewhite(r);
				handles.strTimes{k} = r;
			else						% Names, numeric label & SDS info format
				[t,r] = strtok(r);
				t = ddewhite(t);
				if (t(1) == '?')		% Means get the numeric label as time extracted from file name (OceanColor products)
					% Example names: A2012024021000.L2_LAC_SST4 S1998001130607.L2_MLAC_OC.x.hdf
					[PATH,FNAME,EXT] = fileparts(names{k});
					indDot = strfind(FNAME,'.');
					if (~isempty(indDot) && strcmpi(FNAME(16:17), 'L2'))	% Second case type name
						FNAME(indDot(1):end) = [];
					elseif (~isempty(EXT) && strcmpi(EXT(2:3), 'L2'))		% First case type name (nothing to do)
					else
						errordlg(sprintf('This "%s" is not a MODIS type name',names{k}),'ERROR'),	return
					end
					% Compose name as YYYY.xxxxx where 'xxxxx' is the decimal day of year truncated to hour precision
					%t = sprintf('%s.%f',FNAME(2:5),sscanf(FNAME(6:8),'%f') + sscanf(FNAME(9:10),'%f')/24); 
					t = sprintf('%f',sscanf(FNAME(6:8),'%f') + sscanf(FNAME(9:10),'%f')/24); 
				end
				handles.strTimes{k} = t;
				SDSinfo{k} = ddewhite(r);
			end
		end
	else								% Only one column with fnames
		for (k = 1:m)
			if ( isempty(names{k}) ),	continue,			end		% Jump empty lines
			if ( names{k}(1) == '#'),	c(k) = true;		continue,	end
			if ( names{k}(1) == '@')
				caracol(k) = true;
				handles.changeCD_msg{n_msg} = names{k}(2:end);
				n_msg = n_msg + 1;
				continue
			end
			handles.strTimes{k} = sprintf('%d',k);
		end
	end

	if (any(c))					% Remove eventual comment lines
		names(c) = [];			handles.strTimes(c) = [];		caracol(c) = [];
	end
	m = numel(names);			% Count remaining ones
	
	if (m == 0)
		errordlg('Beautiful service. Your list file had nothing but a bunch of trash.','Error'),	return
	end

	% -------------- Check if we have a Sub-Datasets request ------------------
	if (n_column == 3 && strncmpi(SDSinfo{1}, 'sds', 3))
		if (any(c)),	SDSinfo(c) = [];	end
		handles.SDSinfo = SDSinfo;
	elseif (n_column == 2 && strncmpi(handles.strTimes{1}, 'sds', 3))		% Two cols with SDS info in the second
		handles.SDSinfo = handles.strTimes;
		for (k = 1:m),		handles.strTimes{k} = sprintf('%d',k);		end
	end
	if (~isempty(handles.SDSinfo))
		handles.SDSthis = str2double(handles.SDSinfo{1}(4:end));
	end
	% --------------------------------------------------------------------------

	handles.shortNameList = cell(m,1);      % To hold grid names with path striped
	for (k = 1:m)
		[PATH,FNAME,EXT] = fileparts(names{k});
		if (isempty(PATH))
			handles.shortNameList{k} = names{k};
			names{k} = [PathName names{k}];
		else
			handles.shortNameList{k} = [FNAME EXT];
		end
	end

	% Check that the files provided in list do exist
	if (~any(caracol))			% It makes no sense to test existance of files in other DVDs
		c = false(m,1);
		for (k = 1:m)
			c(k) = (exist(names{k},'file') ~= 2);		% Flag to kill all non-existant files
		end
		names(c) = [];		handles.shortNameList(c) = [];
	end

	handles.nameList = names;
	handles.caracol = caracol;
	set(handles.edit_namesList, 'String', fname)
	set(handles.listbox_list,'String',handles.shortNameList)
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function radio_conv2netcdf_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_multiBand handles.radio_conv2vtk],'Val',0)

% -----------------------------------------------------------------------------------------
function radio_conv2vtk_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_multiBand handles.radio_conv2netcdf],'Val',0)

% -----------------------------------------------------------------------------------------
function radio_multiBand_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_conv2netcdf handles.radio_conv2vtk],'Val',0)

% -----------------------------------------------------------------------------------------
function check_region_CB(hObject, handles)
	if (get(hObject,'Val'))
		set([handles.edit_north handles.edit_south handles.edit_west handles.edit_east],'Enable','on')
	else
		set([handles.edit_north handles.edit_south handles.edit_west handles.edit_east],'Enable','off')
	end

% -----------------------------------------------------------------------------------------
function edit_north_CB(hObject, handles)
	x1 = str2double(get(hObject,'String'));
	if (isnan(x1)),		set(hObject,'String',''),	return,		end
	x2 = get(handles.edit_south,'String');
	if (~isempty(x2))
		x2 = str2double(x2);
		if (x2 >= x1) 
			errordlg('North Latitude <= South Latitude','Error in Latitude limits')
			set(hObject,'String','')
		end
	end

% -----------------------------------------------------------------------------------------
function edit_south_CB(hObject, handles)
	x1 = str2double(get(hObject,'String'));
	if (isnan(x1)),		set(hObject,'String',''),	return,		end
	x2 = get(handles.edit_north,'String');
	if (~isempty(x2))
		x2 = str2double(x2);
		if (x2 <= x1) 
			errordlg('South Latitude >= North Latitude','Error in Latitude limits')
			set(hObject,'String','')
		end
	end

% -----------------------------------------------------------------------------------------
function edit_west_CB(hObject, handles)
	x1 = str2double(get(hObject,'String'));
	if (isnan(x1)),		set(hObject,'String',''),	return,		end
	x2 = get(handles.edit_east,'String');
	if (~isempty(x2))
		x2 = str2double(x2);
		if (x2 <= x1) 
			errordlg('East Longitude <= West Longitude','Error in Longitude limits')
			set(hObject,'String','')
		end
	end

% -----------------------------------------------------------------------------------------
function edit_east_CB(hObject, handles)
	x1 = str2double(get(hObject,'String'));
	if (isnan(x1)),		set(hObject,'String',''),	return,		end
	x2 = get(handles.edit_west,'String');
	if (~isempty(x2))
		x2 = str2double(x2);
		if (x2 >= x1) 
			errordlg('East Longitude <= West Longitude','Error in Longitude limits')
			set(hObject,'String','')
		end
	end

% -----------------------------------------------------------------------------------------
function edit_stripeWidth_CB(hObject, handles)
	x1 = str2double(get(hObject,'String'));
	if (isnan(x1)),		set(hObject,'String','0.5'),	end

% -----------------------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
% ...

	if (~isempty(handles.OneByOneNameList) && handles.OneByOneFirst)	% Files we entered one by one. Must trick to reuse code
		lixoName = [handles.path_tmp 'listName_lixo.txt'];
		fid = fopen(lixoName, 'w');
		for (k = 1:numel(handles.OneByOneNameList))
			fprintf(fid, '%s\n', handles.OneByOneNameList{k});	% Bloody thing doesn't let write all at once
		end
		fclose(fid);
		% Now call push_namesList_CB as if we had a file with the names list. Clever me, no?
		push_namesList_CB([], handles, lixoName);
		builtin('delete', lixoName);
		handles = guidata(handles.figure1);		% Get the updated version
		handles.OneByOneFirst = false;			% To get out of the otherwise dead-end logic
		guidata(handles.figure1, handles)
	end

	if (isempty(handles.nameList))
		errordlg('No files to work on. You either didn''t give them or all names are wrong.','Error')
		return
	end

	got_R = false;	west = [];			east = [];		south = [];		north = [];
	if (get(handles.check_region,'Val'))
		north = str2double(get(handles.edit_north,'String')); 
		west = str2double(get(handles.edit_west,'String')); 
		east = str2double(get(handles.edit_east,'String')); 
		south = str2double(get(handles.edit_south,'String')); 
		if ( any(isempty([west east south north])) )
			errordlg('One of the region limits was not provided','Error'),	return
		end
		got_R = true;
	end

	cut2cdf(handles, got_R, west, east, south, north)

% -----------------------------------------------------------------------------------------
function cut2tif(handles, got_R, west, east, south, north, FileName)
% Save into a multi-band GeoTIFF file

	[pato, fname, EXT] = fileparts(FileName);
	if (isempty(EXT)),		FileName = [FileName '.tiff'];	end
	fname = FileName;

	att = gdalread(handles.nameList{1}, '-M');

	opt_R = ' ';		head = att.GMT_hdr;
	% If user wants a sub-region
	if (got_R)			% We must give the region in pixels since the image is trully not georeferenced (comment for nasa HDF)
		cp = round(([west east] - head(1)) / head(8));
		rp = round(([south north] - head(3)) / head(9));
		if (cp(1) < 0 || cp(2) > att.RasterXSize)		% Almost sure it should be >=
			msg = 'Sub-region West/Est is outside that grid''s limits';
			errordlg(msg, 'ERROR'),		error(msg)
		end
		if (rp(1) < 0 || rp(2) > att.RasterYSize)		% Almost sure it should be >=
			msg = 'Sub-region South/North is outside that grid''s limits';
			errordlg(msg, 'ERROR'),		error(msg)
		end
		head(1) = head(1) + cp(1)*head(8);		head(2) = head(1) + cp(2)*head(8);
		head(3) = head(3) + rp(1)*head(9);		head(4) = head(3) + rp(2)*head(9);
		rows = att.RasterYSize;
		rp = rows - rp -1;		rp = [rp(2) rp(1)];
		opt_R = sprintf('-r%d/%d/%d/%d',cp(1:2),rp(1:2));
	end

    nSlices = numel(handles.nameList);
    img = gdalread(handles.nameList{1}, opt_R);
	n_row = size(img, 1);
	n_col = size(img, 2);
	for (k = 2:nSlices)
		set(handles.listbox_list,'Val',k),		pause(0.01)			% Show advance
    	Z = gdalread(handles.nameList{k}, opt_R);
		ny = size(Z, 1);	nx = size(Z, 2);
		if ((nx ~= n_col) || (ny ~= n_row))
			errordlg('This image has not the same size as precedentes.','ERROR'),	return
		end
		n_col = nx;		n_row = ny;
		img = cat(3, img, Z);
	end
	clear Z
	
	hdr.name = fname;		hdr.driver = 'GTiff';
	hdr.projWKT = att.ProjectionRef;
	hdr.Xinc = head(8);		hdr.Yinc = head(9);
	hdr.ULx = head(1);		hdr.ULy = head(4);
	if (~isempty(att.GCPvalues)),	hdr.gcp = att.GCPvalues;	end
	
	gdalwrite(img,hdr)
	set(handles.listbox_list,'Val',1)

% -----------------------------------------------------------------------------------------
function cut2cdf(handles, got_R, west, east, south, north)
% Save into a multi-layer netCDF file

	if (get(handles.radio_conv2netcdf,'Val'))		% netCDF format
		this_ext = '.nc';							txt0 = '*.nc;*.grd';
		txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
	elseif (get(handles.radio_conv2vtk,'Val'))		% VTK
		this_ext = '.vtk';							txt0 = '*.vtk';
		txt1 = 'VTK format (*.vtk)';				txt2 = 'Select output VRT file';
	else											% Multi-band
		txt0 = '*.tiff;*.tif';
		txt1 = '(Geo)Tiff format (*.tiff)';			txt2 = 'Select output Tiff file';
	end
	[FileName,PathName] = put_or_get_file(handles,{txt0,txt1; '*.*', 'All Files (*.*)'},txt2,'put');
	if isequal(FileName,0),		return,		end
	[pato, fname, EXT] = fileparts(FileName);
	if (isempty(EXT)),		FileName = [fname this_ext];	end
	grd_out = [PathName FileName];

	if (get(handles.radio_multiBand,'Val'))		% Multi-band. We now pass the hand to its own function
		cut2tif(handles, got_R, west, east, south, north, grd_out)
		return
	end

	% Read relevant metadata. Attention, if we have a subdataset request ATT holds the attribs of the SDS
	[head, opt_R, slope, intercept, base, is_modis, is_linear, is_log, att, do_SDS] = ...
		get_headerInfo(handles, handles.nameList{1}, got_R, west, east, south, north);

	handles.geog = 1;			handles.head = head;
	handles.was_int16 = 0;		handles.computed_grid = 0;

	if (get(handles.radio_conv2vtk,'Val')),		fid = write_vtk(handles, grd_out, 'hdr');	end

	nSlices = numel(handles.nameList);
	n_cd = 1;
	for (k = 1:nSlices)
		set(handles.listbox_list,'Val',k),		pause(0.01)			% Show advance

		if (handles.caracol(k))				% Ai, we need to change CD
			msg = handles.changeCD_msg{n_cd};
			resp = yes_or_no('string',['OK, this one is over. ' msg '  ... and Click "Yes" to continue']);
			n_cd = n_cd + 1;
			if (strcmp(resp, 'No')),	return
			else						continue
			end
		end

		if (do_SDS && k > 1)		% If we have an SDS request, get the attribs of that SDS (needed in getZ)
			[att, do_SDS] = get_headerInfo(handles, handles.nameList{k}, got_R, west, east, south, north);
		end

		% In the following, if any of slope, intercept or base changes from file to file ... f
		NoDataValue = att.Band(1).NoDataValue;		% Backup it because it might be changed for other (now unknow) reasons.
		[Z, handles.have_nans, att] = ...
			getZ(handles.nameList{k}, att, is_modis, is_linear, is_log, slope, intercept, base, opt_R, handles);
		att.Band(1).NoDataValue = NoDataValue;

		% Check if all grids have the same size
		if (k == 1)
			n_rows = size(Z,1);		n_cols = size(Z,2);
		elseif (~isequal([n_rows n_cols],[size(Z,1) size(Z,2)]))
			warndlg(['The grid ' handles.nameList{k} ' has different size than precedents. Jumping it.'],'Warning')
			continue
		end

		if (isfield(att, 'hdrModisL2') && ~isempty(att.hdrModisL2))	% Grid was very likely reinterpolated. Update header
			handles.head = att.GMT_hdr;
		end

		if (get(handles.radio_conv2vtk,'Val'))				% Write this layer of the VTK file and continue
			write_vtk(fid, grd_out, Z);
			continue
		end
		
		if ( isa(Z,'int8') && (min(Z(:)) >= 0) )
			grdutils(Z,'-c');								% Shift by -128 so it goes well with the uint8 add_off elsewere
		end

		if (isa(Z,'single'))
			zz = grdutils(Z,'-L');		handles.head(5:6) = [zz(1) zz(2)];
		else			% min/max is bugged when NaNs in singles
			handles.head(5:6) = [double(min(Z(:))) double(max(Z(:)))];
		end

		% Must treate compiled version differently since, biggest of misteries, nc_funs
		% than hangs when writting unlimited variables
		if (~handles.IamCompiled)
			if (~isempty(handles.strTimes)),		t_val = handles.strTimes{k};
			else									t_val = sprintf('%d',k - 1);
			end
			if (k == 1)
				nc_io(grd_out, ['w-' t_val '/time'], handles, reshape(Z,[1 size(Z)]))
			else
				kk = k - n_cd;			% = k - 1 - (n_cd - 1)	We need this when we had "@ change CD" messages
				nc_io(grd_out, sprintf('w%d\\%s', kk, t_val), handles, Z)
			end
		else
			if (k == 1)
				handles.levelVec = str2double(handles.strTimes);
     			nc_io(grd_out,sprintf('w%d/time',nSlices), handles, reshape(Z,[1 size(Z)]))
			else
				kk = k - n_cd;			% = k - 1 - (n_cd - 1)	We need this when we had "@ change CD" messages
				nc_io(grd_out, sprintf('w%d', kk), handles, Z)
			end
		end
	end

	if (get(handles.radio_conv2vtk,'Val')),		fclose(fid);	end
	set(handles.listbox_list,'Val',1)

% -----------------------------------------------------------------------------------------
function [head, opt_R, slope, intercept, base, is_modis, is_linear, is_log, att, do_SDS] = ...
			get_headerInfo(handles, name, got_R, west, east, south, north)
% Get several direct and inderect (computed) informations about the file NAME or one of its subdatasets.
% The [att, do_SDS] = get_headerInfo(...) form is also supported and used when processing L2 files.

	[att, do_SDS] = get_att(handles, name);

	% GDAL wrongly reports the corners as [0 nx] [0 ny] when no SRS
	if ( isequal((att.Corners.LR - att.Corners.UL), [att.RasterXSize att.RasterYSize]) && ~all(att.Corners.UL) )
		att.GMT_hdr(1:4) = [1 att.RasterXSize 1 att.RasterYSize];
	end
	
	att.fname = name;			% This case needs it
	[head , slope, intercept, base, is_modis, is_linear, is_log, att, opt_R] = ...
		getFromMETA(att, got_R, handles, west, east, south, north);
	
	if (nargout <= 2)			% Short form
		head = att;
		if (nargout == 2),		opt_R = do_SDS;		end
	end

% -----------------------------------------------------------------------------------------
function [att, indSDS] = get_att(handles, name)
% Get the attributes of the root file or, in case we have one, of the requested subdataset

	indSDS = 0;
	att = get_baseNameAttribs(name);

	if ( att.RasterCount == 0 && ~isempty(att.Subdatasets) )	
		indSDS = 1;
		if (~isempty(handles.SDSinfo))
			indSDS = handles.SDSthis * 2 - 1;
			ind = strfind(att.Subdatasets{indSDS}, '=');
		elseif (strncmp(att.DriverShortName, 'HDF4', 4))	% Some MODIS files
			ind = strfind(att.Subdatasets{1}, '=');
		else
			errordlg('File has Sub-Datasets but you told me nothing about it.','ERROR')
			error('File has Sub-Datasets but you told me nothing about it.')
		end
		AllSubdatasets = att.Subdatasets;				% Copy this for keeping it as a subdataset field too
		FileName = att.Subdatasets{indSDS}(ind+1:end);	% First "ind" chars are of the form SUBDATASET_1_NAME=
		att = gdalread(FileName,'-M','-C');				% Try again (it will probably fail on ziped files)
		att.AllSubdatasets = AllSubdatasets;			% A non-standard that is also in some cases set in Mirone
	end

% -----------------------------------------------------------------------------------------
function [head , slope, intercept, base, is_modis, is_linear, is_log, att, opt_R] = ...
	getFromMETA(att, got_R, handles, west, east, south, north)
% Get complementary data from an att struct. This is mostly a helper function to get_headerInfo()
% but it is detached from it because in this way it can be called by exterir code. Namelly by Mirone.
% The cases addressed here are some ones raised by HDF files.
% WARNING: the ATT struct must have and extra field att.fname (to be eventualy used by hdfread)
% NOTE: When called from outside (e.g Mirone) use only the [...] = getFromMETA(att) form
% NOTE: For HDF files the ATT struct will be added the field 'hdrInfo' (returned when hdrfinfo)

	if (nargin == 1),	got_R = false;		end

	opt_R = ' ';	is_modis = false;		is_linear = false;		is_log = false;
	slope = 1;		intercept = 0;			base = 1;
	modis_or_seawifs = false;				is_HDFEOS = false;		is_ESA = false;
	att.hdrInfo = [];
	head = att.GMT_hdr;

	if ( ~isempty(att.Metadata) && ~isempty(search_scaleOffset(att.Metadata, 'HDFEOSVersion')) )
		is_HDFEOS = true;
		if ( isnan(search_scaleOffset(att.Metadata, 'ENVISAT')) )	% Poor trick to find ESA (well, ENVISAT) products
			is_ESA = true;
		end
	end

	if ( ~is_HDFEOS && ~isempty(att.Metadata) && ...
			~isempty(search_scaleOffset(att.Metadata, 'MODIS')) || ~isempty(search_scaleOffset(att.Metadata, 'SeaWiFS')) )
		modis_or_seawifs = true;
	end

	if ( modis_or_seawifs && strncmp(att.DriverShortName, 'HDF4', 4) && ~isempty(search_scaleOffset(att.Metadata, 'Level-2')) )
		out = search_scaleOffset(att.Metadata, 'slope');
		if (~isempty(out))		% Otherwise, no need to search for a 'intercept'
			slope = out;
			if (slope ~= 1)
				out = search_scaleOffset(att.Metadata, 'intercept');
				if (~isempty(out)),		intercept = out;	end
				is_linear = true;
			end
		end
		head(1:4) = [1 att.RasterXSize 1 att.RasterYSize];
		att.Band(1).NoDataValue = -32767;	% Shity format doesn't declare this null part (good for SST)
		if (isfield(att, 'subDsName'))		% Known values (found by file inspection)
			NoDataValue = guess_nodataval(att.subDsName);
			if (~isempty(NoDataValue)),		att.Band(1).NoDataValue = NoDataValue;		end
		end
		is_modis = true;					% We'll use this knowledge to 'avoid' Land pixels = -32767  
		att.hdrModisL2 = hdf_funs('hdfinfo', att.fname);
		if (got_R)			% For L2 files we cannot find -R here (before sensor to geog coords conversion)
							% So we store the croping info in 'att' to use later after reinterpolation
			att.crop_info.opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g',west,east,south,north);
			att.crop_info.limits = [west east south north];
			got_R = false;	% So that last block in this function won't try to execute.
		end
	elseif ( modis_or_seawifs && strncmp(att.DriverShortName, 'HDF4', 4) )
		x_max = search_scaleOffset(att.Metadata, 'Easternmost');	% Easternmost Latitude=180
		x_min = search_scaleOffset(att.Metadata, 'Westernmost');	% Westernmost Latitude=-180
		y_max = search_scaleOffset(att.Metadata, 'Northernmost');	% Northernmost Latitude=90
		y_min = search_scaleOffset(att.Metadata, 'Southernmost');	% Southernmost Latitude=-90
		dx = (x_max - x_min) / att.RasterXSize;
		dy = dx;
		x_min = x_min + dx/2;		x_max = x_max - dx/2;	% Orig data was pixel registered
		y_min = y_min + dx/2;		y_max = y_max - dx/2;
		head(1:4) = [x_min x_max y_min y_max];
		head(8:9) = dx;
		att.Corners.UL = [x_min y_max];			
		att.Corners.LR = [x_max y_min];			
		if (got_R)			% We must give the region in pixels since the image is trully not georeferenced
			rows = att.RasterYSize;
		end
		head(7) = 0;		% Make sure that grid reg is used

		% Get the the scaling equation and its parameters
		if ( ~isempty(search_scaleOffset(att.Metadata, 'linear')) )
			slope = search_scaleOffset(att.Metadata, 'Slope');	% att.Metadata{47} -> Slope=0.000717185
			intercept = search_scaleOffset(att.Metadata, 'Intercept');
			is_linear = true;
		elseif ( ~isempty(search_scaleOffset(att.Metadata, 'logarithmic')) )
			base = search_scaleOffset(att.Metadata, 'Base');	% att.Metadata{41} -> Base=10
			slope = search_scaleOffset(att.Metadata, 'Slope');	% att.Metadata{41} -> Slope=0.000717185
			intercept = search_scaleOffset(att.Metadata, 'Intercept');
			is_log = true;
		end

		att.Band(1).NoDataValue = 65535;						% Shity format doesn't declare this.
		nv = search_scaleOffset(att.Metadata, 'Fill');			% But sometimes (some SeaWifs) it exists
		if (~isempty(nv)),	att.Band(1).NoDataValue = nv;	end
		if (head(5) == att.Band(1).NoDataValue),	head(5) = NaN;	end		% Force later recomputing of array min/max
		att.GMT_hdr = head;			% We need this updated
		is_modis = true;			% We'll use this knowledge to 'avoid' Land pixels = 65535

	elseif ( ~is_HDFEOS && ~modis_or_seawifs && strncmp(att.DriverShortName, 'HDF4', 4)  )		% TEMP -> SST PATHFINDER
		finfo = hdf_funs('hdfinfo', att.fname);
		if (strcmpi(finfo.SDS.Attributes(11).Name, 'slope'))
			slope = double(finfo.SDS.Attributes(11).Value);		% = 0.075;
			intercept = double(finfo.SDS.Attributes(12).Value);	% = -3.0;
		else
			out = search_scaleOffset(finfo.SDS.Attributes, 'slope');
			if (~isempty(out))		% Otherwise, no need to search for a 'intercept'
				slope = out;
				out = search_scaleOffset(finfo.SDS.Attributes, 'intercept');
				if (~isempty(out)),		intercept = out;	end
			end
			if (slope == 1)			% We may have a netCDF style naming. Check
				out = search_scaleOffset(finfo.SDS.Attributes, 'scale_factor');
				if (~isempty(out))
					slope = out;
					out = search_scaleOffset(finfo.SDS.Attributes, 'add_off', 7);	% Sometimes is ADD_OFF, others ADD_OFFSET and no-one is killed for that
					if (~isempty(out)),		intercept = out;	end
				end
			end
		end
		lat = finfo.SDS.Dims(1).Scale;		% Get the latitudes
		lon = finfo.SDS.Dims(2).Scale;		% Get the longitudes
		if (isnumeric(lat))
			x_min = lon(1);			x_max = lon(end);
			y_min = min(lat(1),lat(end));
			y_max = max(lat(1),lat(end));
			att.GMT_hdr(1:4) = [x_min x_max y_min y_max];	% We need this updated
			att.Corners.UL = [x_min y_max];			
			att.Corners.LR = [x_max y_min];			
			dx = lon(3) - lon(2);	dy = abs(lat(3) - lat(2));
			att.GMT_hdr(7:9) = [0 dx dy];
		else				% If not, use array size as coordinates
			x_min = 1;		x_max = finfo.SDS.Dims(2).Size;
			y_min = 1;		y_max = finfo.SDS.Dims(1).Size;
			dx = 1;			dy = 1;
		end
		head(1:4) = [x_min x_max y_min y_max];
		head(8:9) = [dx dy];
		if (got_R)			% We must give the region in pixels since the image is trully not georeferenced
			rows = finfo.SDS.Dims(1).Size;					% Number of rows
		end
		head(7) = 0;		% Make sure that grid reg is used
		is_linear = true;
		att.hdrInfo = finfo;
	elseif ( is_HDFEOS && ~is_ESA)		% This case might not be complete as yet.
		x_min = search_scaleOffset(att.Metadata, 'WESTBOUNDINGCOORDINATE');
		x_max = search_scaleOffset(att.Metadata, 'EASTBOUNDINGCOORDINATE');
		y_min = search_scaleOffset(att.Metadata, 'SOUTHBOUNDINGCOORDINATE');
		y_max = search_scaleOffset(att.Metadata, 'NORTHBOUNDINGCOORDINATE');
		rows = search_scaleOffset(att.Metadata, 'DATAROWS');
		cols = search_scaleOffset(att.Metadata, 'DATACOLUMNS');
		dx = (x_max - x_min) / cols;
		dy = (y_max - y_min) / rows;
		att.GMT_hdr(1:4) = [x_min x_max y_min y_max];	% We need this updated
		att.GMT_hdr(7:9) = [0 dx dy];
		head = att.GMT_hdr;
		att.Corners.UL = [x_min y_max];			
		att.Corners.LR = [x_max y_min];			
		is_linear = true;
	elseif ( is_HDFEOS && is_ESA)		% One more incredible messy HDF product. Nothing (spatialy) reliable inside.
		ind = strfind(att.fname, '_');
		tmp = att.fname(ind(end)+1:end);
		indH = strfind(tmp, 'H');		indV = strfind(tmp, 'V');		indDot = strfind(tmp, '.');
		if (isempty(indH) || isempty(indV))
			errordlg('Sorry but this is not a 5 degrees tile of the super non-documented ESA product. Don''t know how to proceed.','Error')
			return
		end
		% The following is crazy. ESA actually uses a grid registration schema but calls it pixel reg.
		% Hence the left side of each tile is aligned with a multiple of 5 but lacks the last col/row.
		x_min = sscanf(tmp(indH+1:indV-1), '%d') * 5 - 180;
		y_max = 90 - sscanf(tmp(indV+1:indDot-1), '%d') * 5;
		rows = att.RasterYSize;
		dx = 5 / att.RasterXSize;	dy = 5 / att.RasterYSize;
		x_max = x_min + 5 - dx;		y_min = y_max - 5 + dx;
		att.GMT_hdr(1:4) = [x_min x_max y_min y_max];	% We need this updated
		att.GMT_hdr(7:9) = [0 dx dy];
		head = att.GMT_hdr;
		att.Corners.UL = [x_min y_max];		att.Corners.LL = [x_min y_min];
		att.Corners.LR = [x_max y_min];		att.Corners.UR = [x_max y_max];
		att.DriverShortName = sprintf('ATENTION: Spatial info displayed here may NOT be\n reliable due to ESA awfull lack of info in file\n\n%s\n',att.DriverShortName);
		att.ProjectionRef = ogrproj('+proj=longlat +ellps=wgs84 +nodefs');
	else					% Other types
		if (got_R),		rows = att.RasterYSize;		end
		x_min = head(1);	y_min = head(3);
		dx = head(8);		dy = head(9);
	end

	% If user wants a sub-region
	if (got_R)			% We must give the region in pixels since the image is trully not georeferenced (comment for nasa HDF)
		cp = round(([west east] - x_min) / dx);
		rp = round(([south north] - y_min) / dx);
		if (cp(1) < 0 || cp(2) > att.RasterXSize)		% Almost sure it should be >=
			msg = 'Sub-region West/Est is outside that grid''s limits';
			errordlg(msg, 'ERROR'),		error(msg)
		end
		if (rp(1) < 0 || rp(2) > att.RasterYSize)		% Almost sure it should be >=
			msg = 'Sub-region South/North is outside that grid''s limits';
			errordlg(msg, 'ERROR'),		error(msg)
		end
		head(1) = x_min + cp(1)*dx;		head(2) = x_min + cp(2)*dx;
		head(3) = y_min + rp(1)*dy;		head(4) = y_min + rp(2)*dy;
		rp = rows - rp -1;		rp = [rp(2) rp(1)];
		opt_R = sprintf('-r%d/%d/%d/%d',cp(1:2),rp(1:2));
	end

% ----------------------------------------------------------------------------------------
function out = search_scaleOffset(attributes, what, N)
% Search for the WHAT attribute in ATTRIBUTES. If find return its VALUE.
% Used to search for slope/intercept or scale_factor/add_offset in HDF files
	out = [];
	if (isa(attributes, 'struct'))
		if (nargin == 2)						% Exact search for WHAT
			for (k = numel(attributes):-1:1)					% Start from the bottom because they are likely close to it 
				if ( strcmpi(attributes(k).Name, what) )
					out = double(attributes(k).Value);
					break
				end
			end
		else									% Search with a strncmp. Motivated by the uterly stupid play with ADD_OFF & ADD_OFFSET
			for (k = numel(attributes):-1:1)				% Start from the bottom because they are likely close to it 
				if ( strncmpi(attributes(k).Name, what, N) )
					out = double(attributes(k).Value);
					break
				end
			end
		end
	else					% Not tested but it must be a cell array (the att.Metadata)
		for (k = 1:numel(attributes))
			id = strfind(attributes{k}, what);
			if ( ~isempty(id) )
				id_eq = strfind(attributes{k},'=');			% Find the '=' sign
				if (numel(id_eq) > 1),	continue,	end		% Sometimes comments have also the '=' char
				out = str2double(attributes{k}(id_eq+1:end));
				break
			end
		end
	end

% -----------------------------------------------------------------------------------------
function NoDataValue = guess_nodataval(DsName)
% Make an educated guess of the no-data value of HDF4 (MODIS) files
	NoDataValue = [];
	if (strncmp(DsName, 'sst', 3)),			NoDataValue = -32767;
	elseif (strcmp(DsName, 'l3m_data')),	NoDataValue = 65535;		% L3 SST
	elseif (strcmp(DsName, 'chlor_a')),		NoDataValue = -1;
	elseif (strncmp(DsName, 'nLw_',4)),		NoDataValue = 0;			% ???
	elseif (strcmp(DsName, 'K_490')),		NoDataValue = -5000;
	elseif (strcmp(DsName, 'tau_869')),		NoDataValue = 0;
	elseif (strcmp(DsName, 'eps_78')),		NoDataValue = 0;
	elseif (strcmp(DsName, 'angstrom_531')),NoDataValue = -32767;
	end

% -----------------------------------------------------------------------------------------
function [Z, have_nans, att] = getZ(fname, att, is_modis, is_linear, is_log, slope, intercept, base, opt_R, handles)
% ATT may be still unknown (empty). In that case it will be returned by read_gdal()
% HANDLES, is transmitted only within internal calls of this function (that is, not from outside calls)

	str_d = [];			IamInteractive = true;
	if (nargin < 9),	opt_R = ' ';	end
	if (nargin == 10 && ~handles.Interactive),	IamInteractive = false;		end		% External calls may be interactive

	if (nargin == 10 && ~isempty(handles.SDSinfo))		% We have a Sub-Dataset request
		clear_att = false;
		if (~isempty(att) && ~isempty(att.Subdatasets))
			ind = strfind(att.Subdatasets{handles.SDSthis * 2 - 1}, '=');
			fname = att.Subdatasets{handles.SDSthis * 2 - 1}(ind+1:end);	% First "ind" chars are of the form SUBDATASET_?_NAME=
		elseif (~isempty(att) && isempty(att.Subdatasets))
			ind = strfind(att.AllSubdatasets{handles.SDSthis * 2 - 1}, '=');
			fname = att.AllSubdatasets{handles.SDSthis * 2 - 1}(ind+1:end);
		else		% isempty(att) = true
			[fname, str_d] = deal_with_compressed(fname);		% MUST GET RID OF THIS (read compressed directly)
			att = gdalread(fname, '-M');	clear_att = true;
			if (~isempty(str_d)),	str_d = fname;		end
		end
		if (clear_att),		att = [];	end
		IamCompiled = handles.IamCompiled;
	else
		% Need to know if "IamCompiled". Since that info is in handles, we need to find it out here
		try			dumb = which('mirone');			IamCompiled = false;
		catch,		IamCompiled = true;
		end
	end

	saveNoData = false;
	if (is_modis && ~isempty(att))
		NoDataValue = att.Band(1).NoDataValue;		% Save this value to reset it after read_gdal()
		saveNoData = true;
	end

	GMT_hdr = [];
	if (~isempty(att))					% Make copies of these
		GMT_hdr = att.GMT_hdr;
		Corners.LR = att.Corners.LR;
		Corners.UL = att.Corners.UL;
	end

	[Z, att, known_coords, have_nans] = read_gdal(fname, att, IamCompiled, IamInteractive, '-C', opt_R, '-U');

	% See if we knew the image coordinates but that knowedge is lost in new att
	if ( ~known_coords && ~isempty(GMT_hdr) )	% For the moment we only know for sure for L2 georeferenced products
		if ( (isequal(att.GMT_hdr(8:9), [1 1]) && ~isequal(GMT_hdr(8:9), [1 1])) || (diff(GMT_hdr(1:2)) > 359.9) )
			att.GMT_hdr = GMT_hdr;		% Recover the header info
			att.Corners.LR = Corners.LR;
			att.Corners.UL = Corners.UL;
		end
	end

	if (~isempty(str_d)),	delete(str_d);		end		% Delete uncompressed file.

	if ( is_modis && ~isempty(att.Band(1).NoDataValue) && saveNoData)	% att.Band... is isempty when all work has been done before
		att.Band(1).NoDataValue = NoDataValue;		% Shity format doesn't inform on the no-data.
	end

	if (isempty(have_nans))			% It's only non-empty when processing L2 products with reinterpolation
		[Z, have_nans, att] = sanitizeZ(Z, att, is_modis, is_linear, is_log, slope, intercept, base);
	end

% -----------------------------------------------------------------------------------------
function [Z, have_nans, att] = sanitizeZ(Z, att, is_modis, is_linear, is_log, slope, intercept, base)
% Take care of possible scaling/offset transformations.
% Have it separate in a function because when processing L2 products we need to apply this right
% before the interpolation (in read_gdal()). This is too early with regard to the normal work-flow
% where it is applied at the end of function getZ

	if (nargin == 2)	% Go again to getFromMETA but without changing again att
		[head, slope, intercept, base, is_modis, is_linear, is_log] = getFromMETA(att);
	end

	ind = [];
	if ( ~isempty(att.Band(1).NoDataValue) && (att.Band(1).NoDataValue == -9999) )		% TEMP -> PATHFINDER
		if ( ~isempty(att.Metadata) && ~isempty(search_scaleOffset(att.Metadata, 'dsp_SubImageName=QUAL')) )
			% Quality flags files cannot be NaNified. 
			% However, they should NOT have a scaling equation either. If quality is [0 7] why scalling?
			is_linear = false;
			% Furthermore, we also need to recast the flag array into int8 (it was uint8) because netCDF doesn't know UINT8
			Z = int8(Z);
		else
			ind = (Z == 0);
		end
	elseif ( ~isempty(att.Band(1).NoDataValue) && ~isnan(att.Band(1).NoDataValue) )
		ind = (Z == (att.Band(1).NoDataValue));
	elseif (isnan(att.Band(1).NoDataValue))		% The nodata is NaN, replace NaNs in Z by zero
		ind = isnan(Z);
	end

	if ( is_modis && (isa(Z, 'int8') || isa(Z, 'uint8')) )
		is_linear = false;		ind = [];
		if (isa(Z, 'uint8')),	Z = int8(Z);	end		% netCDF doesn't know UINT8
	end

	% See if we must apply a scaling equation
	if (is_linear && (slope ~= 1 || intercept ~= 0))
		if (~isa(Z,'single')),		Z = single(Z);		end
		cvlib_mex('CvtScale',Z, slope, intercept)
	elseif (is_log)
		Z = single(base .^ (double(Z) * slope + intercept));
	end
	have_nans = 0;
	if (~isempty(ind))
		if (~isa(Z,'single') && ~isa(Z,'double'))		% Otherwise NaNs would be converted to 0
			Z = single(Z);
		end
		Z(ind) = NaN;		have_nans = 1;
	elseif (isempty(att.Band(1).NoDataValue) && isa(Z,'single'))		% It comes from the interpolation
		att.Band(1).NoDataValue = NaN;
		have_nans = grdutils(Z, '-N');
	end

% -----------------------------------------------------------------------------------------
function [Z, att, known_coords, have_nans] = read_gdal(full_name, att, IamCompiled, IamInteractive, varargin)
% Help function to gdalread that deals with cases when file is compressed.
% ATT is the GDALREAD returned attributes. If empty, we'll get it here
% VARARGIN will normally contain one or more of '-C', opt_R, '-U'
% WARNING: If exist(att.hdfInfo) than att.fname should exist as well (both non standard)
%
% KNOWN_COORDS	Is a logical that whn true the informs the caller that we already know the coordinates for sure
%				and no attempt should be made to fish them from the matadata info.

	have_nans = [];			% Will only become ~[] when input is a L2 product to be interpolated and referenced
	NoDataValue = [];		% Some defaults
	known_coords = false;	% If we know for sure the coords (as for georefed L2 products) tell that to caller
	[full_name, str_d, uncomp_name] = deal_with_compressed(full_name);

	opt_e = '';
	if (IamCompiled),	opt_e = '-e';	end		% Use aguentabar.dll

	if (isempty(att))
		att = get_baseNameAttribs(full_name);
	end

	if (att.RasterCount == 0 && ~isempty(att.Subdatasets) && strncmp(att.DriverShortName, 'HDF4', 4))		% Some MODIS files
		sds_num = 1;
		handles = guidata(gcf);
		if (~isempty(handles.SDSinfo))		sds_num = handles.SDSthis * 2 - 1;		end		% We have a SubDataset request
		ind = strfind(att.Subdatasets{sds_num}, '=');
		full_name = att.Subdatasets{sds_num}(ind+1:end);				% First "ind" chars are of the form SUBDATASET_1_NAME=
		ind = strfind(att.Subdatasets{sds_num+1}, ' ');					% Need to guess the no-data value
		subDsName = att.Subdatasets{sds_num+1}(ind(1)+1:ind(2)-1);		% Get dataset name
		NoDataValue = guess_nodataval(subDsName);
	end

	% att.hdrInfo and att.hdrModisL2 are not default fields of the ATT struct
	try		fname = att.fname;
	catch,	fname = [];
	end
	if (isfield(att, 'hdrInfo') && ~isempty(att.hdrInfo) && (strcmp(att.hdrInfo.SDS.Name,'sst')) )
		% Only particular case dealt now
		Z = hdf_funs('hdfread', att.fname, att.hdrInfo.SDS.Name, 'index', {[1 1],[1 1], [att.RasterYSize att.RasterXSize]});
		Z = flipud(Z);
	else
		opt_L = ' ';	GCPvalues = [];
		if (isfield(att, 'hdrModisL2') && ~isempty(att.hdrModisL2))
			% First check if lon, lat are of the same size of Z. If yes, than there is no need to reinterpolate
			% the location grids ... to their same positions. This has the further advantage that all readings are
			% done with GDAL and so will work also with the stand-alone version.
			lonID = find_in_subdatasets(att.AllSubdatasets, 'longitude');
			ind = strfind(att.AllSubdatasets{lonID},'=');			% Still must rip the 'SUBDATASET_XX_NAME='
			lon_full = gdalread(att.AllSubdatasets{lonID}(ind+1:end), '-L');
			if ( isequal(size(lon_full), [att.RasterYSize att.RasterXSize]) )
				latID = find_in_subdatasets(att.AllSubdatasets, 'latitude');
				ind = strfind(att.AllSubdatasets{latID},'=');
				lat_full = gdalread(att.AllSubdatasets{latID}(ind+1:end), '-L');

			else			% Bad luck. We need to go through the hdfread way.
				% These boys think they are very funy. Oh it's so cute to write the file from right-to-left !!!
				Vg_index = numel(att.hdrModisL2.Vgroup);	% The uncomprehensible MESS never ends. Assume last has things of interst
				lon = fliplr( hdf_funs('hdfread', att.fname, att.hdrModisL2.Vgroup(Vg_index).SDS(1).Name, 'index', {[],[], []}) );
				lat = fliplr( hdf_funs('hdfread', att.fname, att.hdrModisL2.Vgroup(Vg_index).SDS(2).Name, 'index', {[],[], []}) );
				cntl_pt_cols = hdf_funs('hdfread', att.fname, att.hdrModisL2.Vgroup(Vg_index).SDS(3).Name, 'index', {[],[], []});

				cntl_pt_cols = double(cntl_pt_cols);
				lat = double(lat);		lon = double(lon);
				lat_full = zeros(size(lon,1), cntl_pt_cols(end));
				lon_full = zeros(size(lon,1), cntl_pt_cols(end));
				cols_vec = 1:cntl_pt_cols(end);
				for (k = 1:size(lon,1))
					lon_full(k,:) = akimaspline(cntl_pt_cols, lon(k,:), cols_vec);
					lat_full(k,:) = akimaspline(cntl_pt_cols, lat(k,:), cols_vec);
				end
				clear lon lat cols_vec
			end
			x_min = double(min([lon_full(1) lon_full(1,end) lon_full(end,1) lon_full(end)]));	% double because R6.5
			x_max = double(max([lon_full(1) lon_full(1,end) lon_full(end,1) lon_full(end)]));
			y_min = double(min([lat_full(1) lat_full(1,end) lat_full(end,1) lat_full(end)]));
			y_max = double(max([lat_full(1) lat_full(1,end) lat_full(end,1) lat_full(end)]));
			if ( any([x_min x_max y_min y_max] == -999) )		% CAN WE BELIVE THIS???? BUT HAPPENS!!!!!!!!!!!
				% YES, L2 MODIS files can have -999 as coordinates. This is unbelievable but happens.
				% The remedy is to recompute limits and forget the f.. coordinates that will be trimmed by -R below
				indNotFckCoords = (lon_full ~= -999);
				x_min = double( min(min(lon_full(indNotFckCoords))) );
				x_max = double( max(max(lon_full(indNotFckCoords))) );
				y_min = double( min(min(lat_full(indNotFckCoords))) );
				y_max = double( max(max(lat_full(indNotFckCoords))) );
				clear indNotFckCoords
			end
			opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', x_min, x_max, y_min, y_max);

			AllSubdatasets = att.AllSubdatasets;		% Copy of all subdatsets names in this sub-dataset att
			full_name = att.Name;
			NoDataValue = att.Band(1).NoDataValue;		% Save this that has been ESTIMATED before
			if (strcmp(varargin{end}, '-U')),	varargin(end) = [];		end		% We also don't want to UpDown
			opt_L = '-L';
		end

		if (nargout == 2)
			[Z, att] = gdalread(full_name, varargin{:}, opt_L);	% This ATT may be of a subdataset
		else
			Z = gdalread(full_name, varargin{:}, opt_L);
		end
		if (isempty(att.GCPvalues) && ~isempty(GCPvalues)),		att.GCPvalues = GCPvalues;		end
		if ( (att.Band(1).NoDataValue == -1) && (min(Z(:)) == -32767 ) )	% DIRTY PATCH to avoid previous bad nodata guessing
			att.Band(1).NoDataValue = -32767;
			if (~isempty(NoDataValue)),		NoDataValue = -32767;	end
		end
		if (~isempty(NoDataValue)),		att.Band(1).NoDataValue = NoDataValue;	end		% Recover ESTIMATED value

		% For MODIS L2 products, we may still have many things to do
		if (strcmp(opt_L,'-L'))		% We are using this as an indication that the file is MODIS (need clever solution)

			% Here we must sanitize Z in case we must apply a scale/offset transform and/or NaNifying NodataValues
			[Z, have_nans, att] = sanitizeZ(Z, att);
			if (have_nans),		NoDataValue = NaN;		end

			if (IamInteractive)
				what = l2_choices(AllSubdatasets);		% Call secondary GUI to select what to do next
			else
				what = struct('georeference',1,'nearneighbor',0,'mask',0,'coastRes',0,'quality','');	% sensor coords
				ID = find_in_subdatasets(AllSubdatasets, 'qual_sst', 8);	% Check if we have a quality flags array
				if (ID)
					ind = strfind(AllSubdatasets{ID}, '=');			% Yes we have. Use it if not overruled by info in OPTcontrol
					what.qualSDS = AllSubdatasets{ID}(ind+1:end);
					what.quality = 0;
				end
				if (isfield(att, 'crop_info'))			% We have a crop request
					opt_R = att.crop_info.opt_R;
				end
			end

			% Go check if -R or quality flags request exists in OPTcontrol.txt file
			[opt_R, opt_I, opt_C, bitflags, flagsID] = sniff_in_OPTcontrol(opt_R, att);		% Output opt_R gets preference

			if (isempty(what))							% User killed the window, but it's too late to stop so pretend ...
				what =  struct('georeference',1,'nearneighbor',1,'mask',0,'coastRes',0,'quality','');	% sensor coords
			end
			if ( isempty(bitflags) && ~isempty(what.quality) && what.quality < 2 )		% We have a GUI quality request
				qual = gdalread(what.qualSDS, opt_L);
				Z(qual > what.quality) = NoDataValue;
			end

			if (~isempty(what) && what.georeference)	% OK, let's interpolate it into a regular geog grid
				if (~isempty(bitflags))
					ind = strfind(att.AllSubdatasets{flagsID},'=');	% Still must rip the 'SUBDATASET_XX_NAME='
					Zf = gdalread(att.AllSubdatasets{flagsID}(ind+1:end), varargin{:}, opt_L);
					c = false(size(Zf));
					Zf = uint32(Zf);
					for (k = 1:numel(Zf))
						if ( any(bitget(Zf(k),bitflags)) ),	c(k) = true;	end
					end
					clear Zf
					Z(c) = [];		lon_full(c) = [];		lat_full(c) = [];
					clear c
				end
				ind = (Z == (att.Band(1).NoDataValue));
				Z(ind) = [];		lon_full(ind) = [];		lat_full(ind) = [];
				if (isempty(Z))
					errordlg('As a result of applying (probably wrongly) quality flags (in OPTcontrol.txt), no data was left. Aborting','Error')
					error('As a result of applying (probably wrongly) quality flags, no data was left. Aborting')
				end

				if (isempty(opt_I)),	opt_I = '-I0.01';	end
				if (isempty(opt_C)),	opt_C = '-C3';		end		% For gmtmbgrid only
				if (what.nearneighbor)
					lon_full = single(lon_full);			lat_full = single(lat_full);	Z = single(Z);
					[Z, head] = nearneighbor_m(lon_full(:), lat_full(:), Z(:), opt_R, opt_e, '-N2', opt_I, '-S0.04');
				else
					if (~isa(lon_full,'double')),	lon_full = double(lon_full);	lat_full = double(lat_full);	end
					[Z, head] = gmtmbgrid_m(lon_full(:), lat_full(:), double(Z(:)), opt_I, opt_R, '-Mz', opt_C);
					Z = single(Z);
				end
				if (what.mask)
					opt_D = {'-Dc' '-Dl' '-Di' '-Dh' '-Df'};
					opt_D = opt_D{what.coastRes};
					mask = grdlandmask_m(opt_R, opt_D, opt_e, '-I0.01', '-V');
					Z(mask) = NaN;
				end
				att.GMT_hdr = head;
				known_coords = true;				% Signal that coordinates are known and should not be guessed again
				att.Band(1).NoDataValue = [];		% Don't waist time later trying to NaNify again
				x_min = head(1) - head(8)/2;		x_max = head(2) + head(8)/2;		% Goto pixel registration
				y_min = head(3) - head(9)/2;		y_max = head(4) + head(9)/2;		% But not for att.GMT_hdr(7)
			end
			att.RasterXSize = size(Z,2);		att.RasterYSize = size(Z,1);
			att.Band.XSize = size(Z,2);			att.Band.YSize = size(Z,1);
			att.Corners.LL = [x_min y_min];		att.Corners.UL = [x_min y_max];		% CONFIRMAR
			att.Corners.UR = [x_max y_max];		att.Corners.LR = [x_max y_min];
		end
	end

	if (~isempty(fname)),	att.fname = fname;		end		% Needed in some HDF cases
	if (~isempty(str_d)),	delete(uncomp_name);	end		% Delete uncompressed file

% -----------------------------------------------------------------------------------------
function att = get_baseNameAttribs(full_name)
% Get the file's metadata and also tests if an SDS was required but does not exist
	att = gdalread(full_name, '-M');
	if (att.RasterCount > 0)
		handles = guidata(gcf);
		if (~isempty(handles.SDSinfo) && handles.SDSthis > 1 && ~handles.testedDS)
			handles.testedDS = true;
			guidata(handles.figure1, handles)
			errordlg('Input File has no Sub-Datasets so the silly sub-dataset request forced me to abort.','Error')
			error('empilhador: File has no Sub-Datasets so the silly sub-dataset request forced me to abort')
		end
	end

% -----------------------------------------------------------------------------------------
function [opt_R, opt_I, opt_C, bitflags, flagsID] = sniff_in_OPTcontrol(old_R, att)
% Check the OPTcontrol file for particular requests in terms of -R, -I or quality flags
% OPT_R is what the OPTcontrol has in
% BITFLAGS is a vector with the bit number corresponding to the flgs keys in OPTcontrol.
%			Returns [] when no bitflags keywords are found.
% FLAGSID is the subdatset number adress of the flags array (l2_flags)
%			Return 0 when no l2_flags array is found.

	got_flags = false;		bitflags = [];		flagsID = 0;
	opt_I = [];				opt_C = [];
	opt_R = old_R;			% In case we return without finding a new -R
	mir_dirs = getappdata(0,'MIRONE_DIRS');
	if (~isempty(mir_dirs))
		opt_file = [mir_dirs.home_dir '/data/OPTcontrol.txt'];
	else
		return				% Since we cannot find the OPTcontrol file
	end
	if (~(exist(opt_file, 'file') == 2)),		return,		end		% Nickles
	
	fid = fopen(opt_file, 'r');
	c = (fread(fid,'*char'))';      fclose(fid);
	lines = strread(c,'%s','delimiter','\n');   clear c fid;
	m = numel(lines);
	for (k = 1:m)
		if (~strncmp(lines{k},'MIR_EMPILHADOR',8)),	continue,	end
		opt = ddewhite(lines{k}(15:end));
		got_one = false;
		[t,r] = strtok(opt);
		while (t)
			if (strncmp(t,'-R',2)),		opt_R = t;		got_one = true;
			elseif (strncmp(t,'-I',2)),	opt_I = t;		got_one = true;
			elseif (strncmp(t,'-C',2)),	opt_C = t;		got_one = true;
			end
			[t,r] = strtok(ddewhite(r));
		end
		if (got_one),	continue,	end		% Done with this line
		got_flags = true;			% If it comes here means that we have a flags request
		flaglist = opt;
	end

	if (got_flags)
		% Before anything else find the 'fl_flags' array ID. If not found go away right away
		flagsID = find_in_subdatasets(att.AllSubdatasets, 'l2_flags');
		if (~flagsID)
			warndlg('You requested for a FLAGS masking (via OPTcontrol.txt) but this file does not have one "l2_flags" array','Warning')
			return
		end

		fmap = {'ATMFAIL' 1;
				'LAND' 2;
				'HIGLINT' 4;
				'HILT' 5;
				'HISATZEN' 6;
				'STRAYLIGHT' 9;
				'CLDICE' 10;
				'COCCOLITH' 11;
				'HISOLZEN' 13;
				'LOWLW' 15;
				'CHLFAIL' 16;
				'NAVWARN' 17;
				'MAXAERITER' 20;
				'CHLWARN' 22;
				'ATMWARN' 23;
				'NAVFAIL' 26;
				'FILTER' 27;
				'SSTWARN' 28;
				'SSTFAIL' 29};
		
		loc = [0 strfind(flaglist, ',') numel(flaglist)+1];		% Find the ',' separator. Add 2 to easy algo
		c = false(1, size(fmap,1));								% Vector with as many elements as input flags
		fmap_names = fmap(:,1);
		for (k = 1:numel(loc)-1)		% Loop over number of input flag keys
			ind = strcmp(flaglist(loc(k)+1 : loc(k+1)-1), fmap_names);
			n = find(ind);
			if (~isempty(n)),	c(n) = true;	end
		end
		bitflags = [fmap{c,2}];
	end

% -----------------------------------------------------------------------------------------
function ID = find_in_subdatasets(AllSubdatasets, name, ncmp)
% Find the position in the subdatasets array containing the array called 'name'
% NCMP	is an optional argument containing the N if the strNcmp function.
%		If not provided, defaults to numel(name)

	if (nargin == 2),	ncmp = numel(name);		end
	got_it = false;		ID = 0;
	for (k = 2:2:numel(AllSubdatasets))
		ind = strfind(AllSubdatasets{k}, ' ');
		if ( strncmp(AllSubdatasets{k}(ind(1)+1:end), name, ncmp) )
			got_it = true;		break
		end
	end
	if (got_it),	ID = k - 1;		end		% -1 because the subdataset name is one position before

% -----------------------------------------------------------------------------------------
function [full_name, str_d, uncompressed_name] = deal_with_compressed(full_name)
% Check if FULL_NAME is a compressed file. If it is, uncompress it and return the uncompressed
% name as FULL_NAME.
% STR_D informs if decompressing was done. If it is not empty, it means file was decompressed.
% Use this info to eventualy remove the FULL_NAME (note that in this case this not the original file name)

	str_d = [];		do_warn = true;		cext = [];
	[PATH,fname,EXT] = fileparts(full_name);
	uncompressed_name = [PATH filesep fname];		% Only used if file is compressed
	if (strcmpi(EXT,'.bz2'))
		str_d = ['bzip2 -d -q -f -c ' full_name ' > ' uncompressed_name];		cext = EXT;
	elseif (strcmpi(EXT,'.zip') || strcmpi(EXT,'.gz'))
		str_d = ['gunzip -q -N -f -c ' full_name ' > ' uncompressed_name];		cext = EXT;
	end

	if (~isempty(str_d))     % File is compressed.
		[pato,fname,EXT] = fileparts(fname);	% Need to remove the true extension
	
		if (do_warn),	aguentabar(0.5,'title',['Uncompressing ' fname EXT cext]);	end
		if (isunix),	s = unix(str_d);
		elseif ispc,	s = dos(str_d);
		else			errordlg('Unknown platform.','Error');	error('Unknown platform.')
		end
		if ~(isequal(s,0))                  % An error as occured
			errordlg(['Error decompressing file ' full_name],'Error');
			if (do_warn),   aguentabar(1,'title','By'),		end
			error(['Error decompressing file ' full_name])
		end
		if (do_warn),	aguentabar(1,'title','Donne'),		end
		full_name = uncompressed_name;				% The uncompressed file name
	end

% -----------------------------------------------------------------------------------------
function fid = write_vtk(handles, grd_out, arg3)
% Help function to write in the VTK format
	if (isa(arg3,'char') && strcmp(arg3,'hdr'))
		nx = (handles.head(2) - handles.head(1)) / handles.head(8) + ~handles.head(7);
		ny = (handles.head(4) - handles.head(3)) / handles.head(9) + ~handles.head(7);
		nz = numel(handles.nameList);
% 		fid = fopen(grd_out, 'wa');
		fid = fopen(grd_out, 'wb','b');
		fprintf(fid, '# vtk DataFile Version 2.0\n');
		fprintf(fid, 'converted from A B\n');
% 		fprintf(fid, 'ASCII\n');
		fprintf(fid, 'BINARY\n');
		fprintf(fid, 'DATASET RECTILINEAR_GRID\n');
		fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
		fprintf(fid, 'X_COORDINATES %d float\n', nx);
		X = linspace(handles.head(1), handles.head(2), nx);
% 		fprintf(fid, '%.6f\n', X);
		fwrite(fid, X, 'real*4');
		fprintf(fid, 'Y_COORDINATES %d float\n', ny);
		X = linspace(handles.head(3), handles.head(4), ny);
% 		fprintf(fid, '%.6f\n', X);
		fwrite(fid, X, 'real*4');
		fprintf(fid, 'Z_COORDINATES %d float\n', nz);
		X = str2double(handles.strTimes);
% 		fprintf(fid, '%.6f\n', X);
		fwrite(fid, X, 'real*4');
		fprintf(fid, 'POINT_DATA %d\n', nx * ny * nz);
		fprintf(fid, 'SCALARS dono float 1\n');
		fprintf(fid, 'LOOKUP_TABLE default\n');
	else
		fid = handles;			% First arg actually contains the file id
		Z = double((arg3)');
% 		ind = isnan(Z);
% 		if (any(ind(:)))
% 			Z(ind) = 0;
% 		end
% 		fprintf(fid, '%.6f\n', Z(:));		
		fwrite(fid, Z(:), 'real*4');		
	end

% ---------------------------------------------------------------------
function empilhador_LayoutFcn(h1)

set(h1,'Position',[520 532 440 280],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Toolbar','none',...
'Name','Empilhador',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[240 39 191 117],'Style','frame');
uicontrol('Parent',h1,'Position',[358 40 20 15],'String','S','Style','text');
uicontrol('Parent',h1,'Position',[250 98 20 15],'String','W','Style','text');
uicontrol('Parent',h1,'Position',[402 99 20 15],'String','E','Style','text');
uicontrol('Parent',h1,'Position',[356 122 20 15],'String','N','Style','text');

uicontrol('Parent',h1,'Position',[6 254 401 21],...
'BackgroundColor',[1 1 1],...
'Call','empilhador(''edit_namesList_CB'',gcbo,guidata(gcbo))',...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Name of an ascii file with the grids list. One grid name per row',...
'Tag','edit_namesList');

uicontrol('Parent',h1,'Position',[407 252 23 23],...
'Call','empilhador(''push_namesList_CB'',gcbo,guidata(gcbo))',...
'Tooltip','Browse for a grids list file',...
'Tag','push_namesList');

uicontrol('Parent',h1, 'Position',[244 233 165 15],...
'Call','empilhador(''radio_conv2netcdf_CB'',gcbo,guidata(gcbo))',...
'String','Convert to 3D netCDF file',...
'Style','radiobutton',...
'Tooltip','Take a list of files and create a single 3D netCDF file',...
'Value',1,...
'Tag','radio_conv2netcdf');

uicontrol('Parent',h1, 'Position',[244 214 150 15],...
'Call','empilhador(''radio_conv2vtk_CB'',gcbo,guidata(gcbo))',...
'String','Convert to 3D VTK file',...
'Style','radiobutton',...
'Tooltip','Take a list of files and create a single 3D VTK file',...
'Value',0,...
'Tag','radio_conv2vtk');

uicontrol('Parent',h1, 'Position',[244 195 140 15],...
'Call','empilhador(''radio_multiBand_CB'',gcbo,guidata(gcbo))',...
'String','Make multi-band image',...
'Style','radiobutton',...
'Tooltip','Take a list of files and compute a zonal average file',...
'Tag','radio_multiBand');

uicontrol('Parent',h1,'Position',[250 133 115 15],...
'Call','empilhador(''check_region_CB'',gcbo,guidata(gcbo))',...
'String','Use sub-region?',...
'Style','checkbox',...
'Tooltip','Perform computations inside a data sub-region',...
'Tag','check_region');

uicontrol('Parent',h1,'Position',[300 104 71 21],...
'BackgroundColor',[1 1 1],...
'Call','empilhador(''edit_north_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_north');

uicontrol('Parent',h1,'Position',[250 79 71 21],...
'BackgroundColor',[1 1 1],...
'Call','empilhador(''edit_west_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_west');

uicontrol('Parent',h1,'Position',[350 79 71 21],...
'BackgroundColor',[1 1 1],...
'Call','empilhador(''edit_east_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_east');

uicontrol('Parent',h1,'Position',[300 54 71 21],...
'BackgroundColor',[1 1 1],...
'Call','empilhador(''edit_south_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_south');

uicontrol('Parent',h1,'Position',[6 9 225 236],...
'BackgroundColor',[1 1 1],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_list');

uicontrol('Parent',h1,'Position',[340 5 90 21],...
'Call','empilhador(''push_compute_CB'',gcbo,guidata(gcbo))',...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'String','Compute',...
'Tag','push_compute');

% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
function varargout = l2_choices(varargin)
% ... 
	hObject = figure('Vis','off');
	l2_choices_LayoutFcn(hObject);
	handles = guihandles(hObject);

	handles.out.georeference = 0;
	handles.out.nearneighbor = 1;
	handles.out.mask = 0;
	handles.out.coastRes = 0;
	handles.out.quality = '';

	AllSubdatasets = varargin{1};
	ID = find_in_subdatasets(AllSubdatasets, 'qual_sst', 8);
	if (ID)
		set([handles.popup_quality handles.text_quality],'Enable','on')
		ind = strfind(AllSubdatasets{ID}, '=');
		handles.out.qualSDS = AllSubdatasets{ID}(ind+1:end);
	end

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	% UIWAIT makes l2_choices wait for user response
	uiwait(handles.figure1);
	handles = guidata(handles.figure1);
	varargout{1} = handles.out;
	delete(handles.figure1),	drawnow

% -----------------------------------------------------------------------
function radio_sensor_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_interpMin handles.radio_interpNear handles.check_landMask],'Enable','off')
	set(handles.radio_georef,'Val',0)

% -----------------------------------------------------------------------
function radio_georef_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_interpMin handles.radio_interpNear handles.check_landMask],'Enable','on')
	set(handles.radio_sensor,'Val',0)

% -----------------------------------------------------------------------
function radio_interpNear_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set(handles.radio_interpMin,'Val',0)

% -----------------------------------------------------------------------
function radio_interpMin_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set(handles.radio_interpNear,'Val',0)

% -----------------------------------------------------------------------
function check_landMask_CB(hObject, handles)
	if ( get(hObject,'Val') )
		info = getappdata(0,'gmt_version');		% It should never be empty
		if (info.crude ~= 'y')
			errordlg('You do not have GMT installed. So, no coastlines no masking.','Error')
			set(hObject,'Val',0),	return
		elseif (info.full == 'y'),			handles.out.coastRes = 5;
		elseif (info.high == 'y'),			handles.out.coastRes = 4;
		elseif (info.intermediate == 'y'),	handles.out.coastRes = 3;
		elseif (info.low == 'y'),			handles.out.coastRes = 2;
		else								handles.out.coastRes = 1;
		end
		guidata(handles.figure1, handles)
	end

% -----------------------------------------------------------------------
function push_OK_CB(hObject, handles)	
	handles.out.georeference = get(handles.radio_georef, 'Val');
	handles.out.nearneighbor = get(handles.radio_interpNear, 'Val');
	handles.out.mask = get(handles.check_landMask, 'Val');
	if ( strcmp(get(handles.popup_quality,'Enable'),'on') )		% Only if we are using it
		ind = get(handles.popup_quality, 'Val');
		if (ind ~= 1)
			handles.out.quality = 3 - ind;
		end
	end
	guidata(handles.figure1, handles)
	uiresume(handles.figure1);

% -----------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
% The GUI is still in UIWAIT, do UIRESUME
	handles = guidata(hObject);
	handles.out = '';		% User gave up, return nothing
	guidata(handles.figure1, handles);
	uiresume(handles.figure1);

% -----------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
% Check for "escape"
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles = guidata(hObject);
		handles.out = '';	% User said no by hitting escape
		guidata(handles.figure1, handles);
		uiresume(handles.figure1);
	end

% --- Creates and returns a handle to the GUI figure. 
function l2_choices_LayoutFcn(h1)

set(h1, 'Position',[520 400 291 148],...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','L2 product choices',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[2 124 265 22],...
'FontAngle','oblique',...
'FontName','Helvetica',...
'FontSize',11,...
'FontWeight','demi',...
'ForegroundColor',[1 0.50196 0],...
'String','L2 Level product - Needs decisions',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 101 161 16],...
'Call',@l2_choices_uiCB,...
'FontName','Helvetica',...
'String','Plot in Sensor coordinates',...
'Style','radiobutton',...
'Tooltip','Plot data as it is in file (faster, but deformed)',...
'Value',1,...
'Tag','radio_sensor');

uicontrol('Parent',h1, 'Position',[10 78 231 16],...
'Call',@l2_choices_uiCB,...
'FontName','Helvetica',...
'String','Compute georeferenced grid (takes time)',...
'Style','radiobutton',...
'Tooltip','Reinterpolate data to get a georeferenced grid. (It may take 1 minute)',...
'Tag','radio_georef');

uicontrol('Parent',h1, 'Position',[31 55 131 16],...
'Call',@l2_choices_uiCB,...
'Enable','off',...
'FontName','Helvetica',...
'String','Minimum curvature',...
'Style','radiobutton',...
'Value',1,...
'Tooltip','Interpolation method',...
'Tag','radio_interpMin');

uicontrol('Parent',h1, 'Position',[31 36 110 16],...
'Call',@l2_choices_uiCB,...
'Enable','off',...
'FontName','Helvetica',...
'String','Nearneighbor',...
'Style','radiobutton',...
'Tooltip','Interpolation method',...
'Value',0,...
'Tag','radio_interpNear');

uicontrol('Parent',h1, 'Position',[10 6 51 22],...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'ListboxTop',0,...
'String',{'2'; '1'; '0'},...
'Style','popupmenu',...
'Tooltip','Select the least quality level. 2 - worst - means all values. 0 - only the best',...
'Value',3,...
'Tag','popup_quality');

uicontrol('Parent',h1, 'Position',[64 10 75 15],...
'Enable','off',...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Quality factor',...
'Style','text',...
'Tag','text_quality');

uicontrol('Parent',h1, 'Position',[180 48 110 15],...
'Call',@l2_choices_uiCB,...
'Enable','off',...
'String','Apply Land Mask',...
'Style','checkbox',...
'Tooltip','Mask reinterpolated Land pixels (need high definition coastlines installed)',...
'Tag','check_landMask');

uicontrol('Parent',h1, 'Position',[215 7 66 23],...
'Call',@l2_choices_uiCB,...
'FontName','Helvetica',...
'FontSize',9,...
'String','OK',...
'Tag','push_OK');

function l2_choices_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

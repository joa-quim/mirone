function varargout = empilhador(varargin)
% Stacks a bunch of grids into a single 3D file
% (It used to be called zonal_integrator)
%
% WARNING: FOR COMPILING THIS WE NEED TO INCLUDE THE HDF_FUNS.M SRC
%
% NOTE: The gotFromMETA and getZ functions are called directly by mirone

%	Copyright (c) 2004-2009 by J. Luis
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

	if (nargin > 1 && ischar(varargin{1}))
		gui_Callback = str2func(varargin{1});
		[varargout{1:nargout}] = feval(gui_Callback,varargin{2:end});
	else
		h = empilhador_OpeningFcn(varargin{:});
		if (nargout)    varargout{1} = h;   end
	end

% ---------------------------------------------------------------------------------
function hObject = empilhador_OpeningFcn(varargin)
	hObject = figure('Tag','figure1','Visible','off');
	empilhador_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'east')

	if (numel(varargin) > 0)
		handMir = varargin{1};
		handles.home_dir = handMir.home_dir;
		handles.last_dir = handMir.last_dir;
		handles.work_dir = handMir.work_dir;
		handles.IamCompiled = handMir.IamCompiled;		% Need to know due to crazy issue of nc_funs
        d_path = handMir.path_data;
	else
		handles.home_dir = cd;
		handles.last_dir = handles.home_dir;
		handles.work_dir = handles.home_dir;
		handles.IamCompiled = false;
        d_path = [pwd filesep 'data' filesep];
	end
	handles.nameList = [];
	handles.testedDS = false;		% To test if a Sub-Dataset request is idiot

	set([handles.edit_stripeWidth handles.radio_lon handles.radio_lat],'Enable','off')

	% -------------- Import/set icons --------------------------------------------
	load([d_path 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_namesList, 'CData',Mfopen_ico)

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	bgcolor = get(0,'DefaultUicontrolBackgroundColor');
	framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
	frame_size = get(handles.frame1,'Position');
	frame3D(hObject,frame_size,framecolor,'',get(handles.frame1,'UserData'))
	delete(handles.frame1)
	%------------- END Pro look (3D) -------------------------------------------------------

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
        [FileName,PathName] = put_or_get_file(handles, str1,'File with grids list','get');
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
    end

	fid = fopen(fname);
	c = fread(fid,'*char')';      fclose(fid);
	names = strread(c,'%s','delimiter','\n');   clear c fid;
	m = length(names);

	handles.strTimes = cell(m,1);	% To hold time steps as strings
	SDSinfo = cell(m,1);			% To hold Sub Datasets info
	handles.SDSinfo = [];			% If above exists, it will be copied here
	c = false(m,1);
	caracol = false(m,1);

	n_msg = 1;							% Will hold the "change DV messages" counter
	if (n_column > 1)					% When 2nd column holds the 3D numbering
		for (k = 1:m)
			[t,r] = strtok(names{k});
			if (t(1) == '#' || isempty(t)),  	c(k) = true;	continue,	end		% Jump empty and comment lines
			if ( t(1) == '@')
				caracol(k) = true;
				if (~isempty(t(2:end))),		handles.changeCD_msg{n_msg} = t(2:end);		% The '@' was glued with the message
				else							handles.changeCD_msg{n_msg} = r;
				end
				n_msg = n_msg + 1;
				continue
			end

			names{k} = t;
			if (n_column == 2)			% Names & numeric label format
				r = ddewhite(r);
				handles.strTimes{k} = r;
			else						% Names, numeric label & SDS info format
				[t,r] = strtok(r);
				t = ddewhite(t);
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
		handles.SDSthis = str2num(handles.SDSinfo{1}(4:end));
	end
	% --------------------------------------------------------------------------

	handles.shortNameList = cell(m,1);      % To hold grid names with path striped
	for (k=1:m)
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
	set([handles.edit_stripeWidth handles.radio_lon handles.radio_lat],'Enable','off')
	set([handles.radio_zonalInteg handles.radio_conv2vtk],'Val',0)

% -----------------------------------------------------------------------------------------
function radio_conv2vtk_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.edit_stripeWidth handles.radio_lon handles.radio_lat],'Enable','off')
	set([handles.radio_zonalInteg handles.radio_conv2netcdf],'Val',0)

% -----------------------------------------------------------------------------------------
function radio_zonalInteg_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.edit_stripeWidth handles.radio_lon handles.radio_lat],'Enable','on')
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
function radio_lat_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_lon,'Val',0)

% -----------------------------------------------------------------------------------------
function radio_lon_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_lat,'Val',0)

% -----------------------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
% ...

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

	% to netCDF or VTK conversion?
	if ( get(handles.radio_conv2netcdf,'Val') || get(handles.radio_conv2vtk,'Val') )
		cut2cdf(handles, got_R, west, east, south, north)
		return
	end

	% MORE OR LESS DEPRECATED CODE (IT WILL PROBABLY FAIL) - FUNCTIONALITY MOVED TO AQUAMOTO SUPP FUNS
	[head, opt_R, slope, intercept, base, is_modis, is_linear, is_log, N_spatialSize, integDim] = ...
			get_headerInfo(handles, got_R, west, east, south, north);

	att = read_gdal(handles.nameList{1}, [], '-M','-C');

	if ( strcmp(att.DriverShortName, 'HDF4') && att.RasterCount == 0 && ~isempty(att.Subdatasets) )
		ind = strfind(att.Subdatasets{1}, '=');
		FileName = att.Subdatasets{1}(ind+1:end);		% First "ind" chars are of the form SUBDATASET_1_NAME=
		att = read_gdal(FileName,[],'-M','-C');			% Try again
		handles.nameList{1} = FileName;					% This way the first time in loop for below will access the subdataset
	end

	aguentabar(0,'title','Computing zonal means','CreateCancelBtn')

	% Build the vectors to deal with the zonal integration
	dlat = str2double(get(handles.edit_stripeWidth,'String'));
	if (get(handles.radio_lon, 'Val'))
		vecD = floor(head(3)):dlat:ceil(head(4));
		Y = linspace(head(3),head(4), N_spatialSize);
	else
		vecD = floor(head(1)):dlat:ceil(head(2));
		Y = linspace(head(1),head(2), N_spatialSize);
	end
	nStripes = numel(vecD) - 1;
	indStripe = ones(numel(vecD),1);
	for (k = 2:nStripes)
		ind = find(Y >= vecD(k));
		indStripe(k) = ind(1);
	end
	indStripe(end) = N_spatialSize;

	nSeries = numel(handles.nameList);
	allSeries = zeros(nStripes, nSeries);
	for (k = 1:nSeries)
		att = gdalread(handles.nameList{k},'-M','-C');				% First, get only the attribs
		if ( strcmp(att.DriverShortName, 'HDF4') && att.RasterCount == 0 && ~isempty(att.Subdatasets) )
			ind = strfind(att.Subdatasets{1}, '=');
			FileName = att.Subdatasets{1}(ind+1:end);		% First "ind" chars are of the form SUBDATASET_1_NAME=
			[Z,att] =  read_gdal(FileName, [], '-C', opt_R, '-U');
		else
			Z =  read_gdal(handles.nameList{k}, att, '-C', opt_R, '-U');
		end
		this_has_nans = false;
		if (is_modis)
			ind = (Z == 65535);
			if (any(ind(:))),		Z(ind) = 0;		this_has_nans = true;		end
		elseif ( ~isempty(att.Band(1).NoDataValue) && ~isnan(att.Band(1).NoDataValue) )
			ind = (Z == single(att.Band(1).NoDataValue));
			if (any(ind(:))),		Z(ind) = 0;		this_has_nans = true;		end
		elseif (isnan(att.Band(1).NoDataValue))		% The nodata is NaN, replace NaNs in Z by zero
			ind = isnan(Z);
			if (any(ind(:))),		Z(ind) = 0;		this_has_nans = true;		end
		end

		tmp = sum(Z,integDim);			% Add along integration dim
		if (get(handles.radio_lon, 'Val')),		N = size(Z,2);
		else									N = size(Z,1);
		end
		if (this_has_nans)
			tmp2 = sum(ind,integDim);
			tmp = tmp ./ (N - tmp2);
		else
			tmp = tmp / N;
		end
		% Now add all inside each stripe
		for (m = 1:nStripes)
			tmp2 = tmp( indStripe(m):indStripe(m+1) );
			allSeries(m,k) = sum(tmp2) / numel(tmp2);
		end

		h = aguentabar(k/nSeries);
		if (isnan(h)),	break,	end
	end
	if (isnan(h)),	return,		end

	% See if we must apply a scaling equation
	if (is_modis && is_linear)
		allSeries = allSeries * slope + intercept;
	elseif (is_modis && is_log)
		allSeries = base .^ (allSeries * slope + intercept);
	end

	allSeries = single(allSeries);
	zz = grdutils(allSeries,'-L');
	head = [1 nSeries vecD(1) vecD(end) zz(1) zz(2) 0 1 dlat];
	tmp.X = 1:nSeries;		tmp.Y = linspace( (vecD(1)+dlat/2), (vecD(end)-dlat/2), nStripes );
	tmp.head = [head(1:2) tmp.Y(1) tmp.Y(end) head(5:end)];
	tmp.geo = 0;			tmp.name = 'Zonal integration';
	mirone(allSeries, tmp)

% -----------------------------------------------------------------------------------------
function cut2cdf(handles, got_R, west, east, south, north)
% Save into a multi-layer netCDF file

% 	grd_out = 'teste2.nc';
	if (get(handles.radio_conv2netcdf,'Val'))		% netCDF format
		this_ext = '.nc';							txt0 = '*.nc;*.grd';
		txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
	else											% VTK
		this_ext = '.vtk';							txt0 = '*.vtk';
		txt1 = 'VTK format (*.vtk)';				txt2 = 'Select output VRT file';
	end
	[FileName,PathName] = put_or_get_file(handles,{txt0,txt1; '*.*', 'All Files (*.*)'},txt2,'put');
	if isequal(FileName,0),		return,		end
	[pato, fname, EXT] = fileparts(FileName);
	if (isempty(EXT)),		FileName = [fname this_ext];	end
	grd_out = [PathName FileName];

	% Read relevant metadata
	[head, opt_R, slope, intercept, base, is_modis, is_linear, is_log, N_spatialSize, integDim, att] = ...
		get_headerInfo(handles, got_R, west, east, south, north);

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

		% In the following, if any of slope, intercept or base changes from file to file ... f
		[Z, handles.have_nans, att] = getZ(handles.nameList{k}, [], is_modis, is_linear, is_log, slope, intercept, base, opt_R, handles);

		if (get(handles.radio_conv2vtk,'Val'))				% Write this layer of the VTK file and continue
			write_vtk(fid, grd_out, Z);
			continue
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
function [head, opt_R, slope, intercept, base, is_modis, is_linear, is_log, N_spatialSize, integDim, att] = ...
			get_headerInfo(handles, got_R, west, east, south, north)
% ...

	att = read_gdal(handles.nameList{1}, [], '-M','-C');	% Do this because it deals also with ziped files

	if ( att.RasterCount == 0 && ~isempty(att.Subdatasets) )
		if (~isempty(handles.SDSinfo))
			ind = strfind(att.Subdatasets{handles.SDSthis}, '=');
		elseif (strncmp(att.DriverShortName, 'HDF4', 4))	% Some MODIS files
			ind = strfind(att.Subdatasets{1}, '=');
		else
			errordlg('File has Sub-Datasets but you told me nothing about it.','ERROR'),	return
		end
		FileName = att.Subdatasets{1}(ind+1:end);		% First "ind" chars are of the form SUBDATASET_1_NAME=
		att = gdalread(FileName,'-M','-C');				% Try again (it will probably fail on ziped files)
	end

	% GDAL wrongly reports the corners as [0 nx] [0 ny] when no SRS
	if ( isequal([att.Corners.LR - att.Corners.UL],[att.RasterXSize att.RasterYSize]) && ~all(att.Corners.UL) )
		att.GMT_hdr(1:4) = [1 att.RasterXSize 1 att.RasterYSize];
	end

	N_spatialSize = att.RasterYSize;		integDim = 2;
	if (get(handles.radio_lon, 'Val'))
		N_spatialSize = att.RasterXSize;		% Number of points in the spatial dim
		integDim = 1;							% Dimension along which we are going to integrate
	end

	att.fname = handles.nameList{1};			% This case needs it
	[head , slope, intercept, base, is_modis, is_linear, is_log, att, opt_R] = ...
		getFromMETA(att, got_R, handles, west, east, south, north);

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
	att.hdrInfo = [];
	head = att.GMT_hdr;

	if ( ~isempty(att.Metadata) && ~isempty(strfind(att.Metadata{2}, 'MODIS')))
		modis_or_seawifs = true;
	elseif ( ~isempty(att.Metadata) && ~isempty(strfind(att.Metadata{2}, 'SeaWiFS')))
		modis_or_seawifs = true;
	else
		modis_or_seawifs = false;
	end

	if ( modis_or_seawifs && strncmp(att.DriverShortName, 'HDF4', 4) && ~isempty(strfind(att.Metadata{1}, 'Level-2')) )
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
		is_modis = true;						% We'll use this knowledge to 'avoid' Land pixels = -32767  
		att.hdrModisL2 = hdf_funs('hdfinfo', att.fname);
	elseif ( modis_or_seawifs && strncmp(att.DriverShortName, 'HDF4', 4) )
		y_max = str2double(att.Metadata{39}(23:end));		% att.Metadata{39} -> Northernmost Latitude=90
		y_min = str2double(att.Metadata{40}(23:end));		% att.Metadata{40} -> Southernmost Latitude=90
		x_min = str2double(att.Metadata{41}(23:end));		% att.Metadata{41} -> Westernmost Longitude=-180
		x_max = str2double(att.Metadata{42}(23:end));		% att.Metadata{41} -> Easternmost Longitude=-180
		dx = (x_max - x_min) / att.RasterXSize;
		dy = dx;
		x_min = x_min + dx/2;		x_max = x_max - dx/2;	% Orig data was pixel registered
		y_min = y_min + dx/2;		y_max = y_max - dx/2;
		head(1:4) = [x_min x_max y_min y_max];
		head(8:9) = dx;
		att.GMT_hdr(1:4) = head(1:4);	% We need this updated
		att.Corners.UL = [x_min y_max];			
		att.Corners.LR = [x_max y_min];			
		if (got_R)			% We must give the region in pixels since the image is trully not georeferenced
			rows = att.RasterYSize;
		end
		head(7) = 0;		% Make sure that grid reg is used

		% Get the the scaling equation and its parameters
		if ( strcmp(att.Metadata{53}(9:end), 'linear') )
			slope = str2double(att.Metadata{55}(7:end));		% att.Metadata{41} -> Slope=0.000717185
			intercept = str2double(att.Metadata{56}(11:end));	% att.Metadata{41} -> Intercept=-2
			is_linear = true;
		elseif ( strcmp(att.Metadata{53}(9:end), 'logarithmic') )
			base = str2double(att.Metadata{55}(6:end));		 	% att.Metadata{41} -> Base=10
			slope = str2double(att.Metadata{56}(7:end));		% att.Metadata{41} -> Slope=0.000717185
			intercept = str2double(att.Metadata{57}(11:end));	% att.Metadata{41} -> Intercept=-2
			is_log = true;
		end
		att.Band(1).NoDataValue = 65535;		% Shity format doesn't declare this null part
		is_modis = true;			% We'll use this knowledge to 'avoid' Land pixels = 65535
	elseif ( strncmp(att.DriverShortName, 'HDF4', 4) && ~modis_or_seawifs )		% TEMP -> SST PATHFINDER
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
					out = search_scaleOffset(finfo.SDS.Attributes, 'add_offset');
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
		if (get(handles.radio_lon, 'Val')),		N_spatialSize = round(diff(rp) + 1);
		else									N_spatialSize = round(diff(cp) + 1);
		end
	end

% ----------------------------------------------------------------------------------------
function out = search_scaleOffset(attributes, what)
% Search for the WHAT attribute in ATTRIBUTES. If find return its VALUE.
% Used to search for slope/intercept or scale_factor/add_offset in HDF files
	out = [];
	if (isa(attributes, 'struct'))
		for (k = numel(attributes):-1:1)					% Start from the bottom because they are likely close to it 
			if ( strcmpi(attributes(k).Name, what) )
				out = double(attributes(k).Value);			break
			end
		end
	else					% Not tested but it must be a cell array (the att.Metadata)
		for (k = numel(attributes):-1:1)					% Start from the bottom because they are likely close to it
			id = strfind(attributes{k}, what);
			if ( ~isempty(id) )
				id_eq = findstr(attributes{k},'=');			% Find the '=' sign
				out = str2double(attributes{k}(id_eq+1:end));			break
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

	if (nargin < 9),	opt_R = ' ';	end
	str_d = [];
	if (nargin == 10 && ~isempty(handles.SDSinfo))	% We have a Sub-Dataset request
		clear_att = false;
		if (isempty(att))
			[fname, str_d] = deal_with_compressed(fname);
			att = gdalread(fname, '-M');	clear_att = true;
			str_d = fname;
		end
		ind = strfind(att.Subdatasets{handles.SDSthis * 2 - 1}, '=');
		fname = att.Subdatasets{handles.SDSthis * 2 - 1}(ind+1:end);	% First "ind" chars are of the form SUBDATASET_?_NAME=
		if (clear_att),		att = [];	end
	end
	saveNoData = false;
	if (is_modis && ~isempty(att))
		NoDataValue = att.Band(1).NoDataValue;		% Save this value to reset it after read_gdal()
		saveNoData = true;
	end
	[Z,att] = read_gdal(fname, att, '-C', opt_R, '-U');
	if (~isempty(str_d)),	delete(str_d);		end		% Delete uncompressed file.

	if ( is_modis && ~isempty(att.Band(1).NoDataValue) && saveNoData)	% att.Band... is isempty when all work has been done before
		att.Band(1).NoDataValue = NoDataValue;		% Shity format doesn't inform on the no-data.
	end

	ind = [];
	if ( ~isempty(att.Band(1).NoDataValue) && (att.Band(1).NoDataValue == -9999) )		% TEMP -> PATHFINDER
		if ( ~isempty(att.Metadata) && strcmp(att.Metadata{2}, 'dsp_SubImageName=QUAL') )
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
		if (~isa(Z,'single') || ~isa(Z,'double'))		% Otherwise NaNs would be converted to 0
			Z = single(Z);
		end
		Z(ind) = NaN;		have_nans = 1;
	elseif (isempty(att.Band(1).NoDataValue) && isa(Z,'single'))		% It comes from the interpolation
		att.Band(1).NoDataValue = NaN;
		have_nans = grdutils(Z, '-N');
	end

% -----------------------------------------------------------------------------------------
function [Z, att] = read_gdal(full_name, att, varargin)
% Help function to gdalread that deals with cases when file is compressed.
% ATT is the GDALREAD returned attributes. If empty, we'll get it here
% VARARGIN will normally contain one or more of '-C', '-M', opt_R, '-U'
% WARNING: If exist(att.hdfInfo) than att.fname should exist as well (both non standard)

	NoDataValue = [];	% Some defaults
	[full_name, str_d, uncomp_name] = deal_with_compressed(full_name);

	if (isempty(att))
		att = gdalread(full_name, '-M');	% This first call is used in the next test
		if ( att.RasterCount > 0 && numel(varargin) <= 2)		% CONVOLUTED TEST, BUT THIS IS A CONFUSED PROGRAM
			handles = guidata(gcf);
			if ( ~isempty(handles.SDSinfo) && handles.SDSthis > 1 && ~handles.testedDS)
				errordlg('Input File has no Sub-Datasets so the silly sub-dataset request forced me to abort.','Error')
				error('empilhador: File has no Sub-Datasets so the silly sub-dataset request forced me to abort')
				handles.testedDS = true;
				guidata(handles.figure1, handles)
			end
		end
	end

	if ( att.RasterCount == 0 && ~isempty(att.Subdatasets) && strncmp(att.DriverShortName, 'HDF4', 4) )		% Some MODIS files
		sds_num = 1;
		handles = guidata(gcf);
		if (~isempty(handles.SDSinfo))		sds_num = handles.SDSthis * 2 - 1;		end		% We have a SubDataset request
		ind = strfind(att.Subdatasets{sds_num}, '=');
		full_name = att.Subdatasets{sds_num}(ind+1:end);				% First "ind" chars are of the form SUBDATASET_1_NAME=
		ind = strfind(att.Subdatasets{sds_num+1}, ' ');					% Need to guess the no-data value
		subDsName = att.Subdatasets{sds_num+1}(ind(1)+1:ind(2)-1);		% Get dataset name
		NoDataValue = guess_nodataval(subDsName);
	end

	% att.hdrInfo and att.hdrModisL2 are not a default fields of the ATT struct
	fname = [];
	if (isfield(att, 'hdrInfo') && ~isempty(att.hdrInfo) && (strcmp(att.hdrInfo.SDS.Name,'sst')) )
		% Only particular case dealt now
		Z = hdf_funs('hdfread', att.fname, att.hdrInfo.SDS.Name, 'index', {[1 1],[1 1], [att.RasterYSize att.RasterXSize]});
		Z = flipud(Z);
		fname = att.fname;
	else
		opt_L = ' ';	GCPvalues = [];
		if (isfield(att, 'hdrModisL2') && ~isempty(att.hdrModisL2) )
			% These boys think they are very funy. Oh it's so cute to write the file from right-to-left !!!
			lon = fliplr( hdf_funs('hdfread', att.fname, att.hdrModisL2.Vgroup(4).SDS(1).Name, 'index', {[],[], []}) );
			lat = fliplr( hdf_funs('hdfread', att.fname, att.hdrModisL2.Vgroup(4).SDS(2).Name, 'index', {[],[], []}) );
			cntl_pt_cols = hdf_funs('hdfread', att.fname, att.hdrModisL2.Vgroup(4).SDS(3).Name, 'index', {[],[], []});

			cntl_pt_cols = double(cntl_pt_cols);
			lat = double(lat);		lon = double(lon);
			lat_full = zeros(size(lon,1), cntl_pt_cols(end));
			lon_full = zeros(size(lon,1), cntl_pt_cols(end));
			cols_vec = 1:cntl_pt_cols(end);
			for (k = 1:size(lon,1))
				lon_full(k,:) = akimaspline(cntl_pt_cols, lon(k,:), cols_vec);
				lat_full(k,:) = akimaspline(cntl_pt_cols, lat(k,:), cols_vec);
			end
			x_min = lon(end,1);		x_max = lon(1,end);
			y_min = lat(1);			y_max = lat(end,end);
			clear lon lat cols_vec
			opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', x_min, x_max, y_min, y_max);

			fname = att.fname;
			AllSubdatasets = att.AllSubdatasets;
			NoDataValue = att.Band(1).NoDataValue;		% Save this that has been ESTIMATED before
			if (strcmp(varargin{end}, '-U')),	varargin(end) = [];		end		% We also don't want to UpDown
			opt_L = '-L';
		end

		if (nargout == 2)
			[Z, att] = gdalread(full_name, varargin{:}, opt_L);	% This ATT might be of a subdataset
		else
			Z = gdalread(full_name, varargin{:}, opt_L);
		end
		if (isempty(att.GCPvalues) && ~isempty(GCPvalues)),		att.GCPvalues = GCPvalues;		end
		if (~isempty(NoDataValue)),		att.Band(1).NoDataValue = NoDataValue;	end		% Recover ESTIMATED value

		% For MODIS L2 products, we may still have many things to do
		if (strcmp(opt_L,'-L'))		% We are using this as an indication that the file is MODIS (need clever solution)
			what = l2_choices(AllSubdatasets);			% Call secondary GUI to select what to do next
			if ( ~isempty(what.quality) && what.quality < 2 )		% We have a quality request
				qual = gdalread(what.qualSDS, opt_L);
				Z(qual > what.quality) = NoDataValue;
			end
			if ( ~isempty(what) && what.georeference)
				ind = (Z == (att.Band(1).NoDataValue));
				Z(ind) = [];		lon_full(ind) = [];		lat_full(ind) = [];
				if (what.nearneighbor)
					lon_full = single(lon_full);			lat_full = single(lat_full);	Z = single(Z);
					[Z, head] = nearneighbor_m(lon_full(:), lat_full(:), Z(:), opt_R, '-N2', '-I0.01', '-S0.04');
				else
					[Z, head] = gmtmbgrid_m(lon_full(:), lat_full(:), double(Z(:)), '-I0.01', opt_R, '-Mz', '-C5');
					Z = single(Z);
				end
				if (what.mask)
					opt_D = {'-Dc' '-Dl' '-Di' '-Dh' '-Df'};
					opt_D = opt_D{what.coastRes};
					mask = grdlandmask_m(opt_R, opt_D, '-I0.01', '-V');
					Z(mask) = NaN;
				end
				att.GMT_hdr = head;
				att.Band(1).NoDataValue = [];		% Don't wast time later trying to NaNify again
				x_min = head(1) - head(8)/2;		x_max = head(2) + head(8)/2;		% Goto pixel registration
				y_min = head(3) - head(9)/2;		y_max = head(4) + head(9)/2;
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
'Name','Empilhador',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[241 170 30 15],'String','Delta','Style','text');
uicontrol('Parent',h1,'Position',[240 39 191 117],'Style','frame','Tag','frame1');
uicontrol('Parent',h1,'Position',[358 40 20 15],'String','S','Style','text');
uicontrol('Parent',h1,'Position',[250 98 20 15],'String','W','Style','text');
uicontrol('Parent',h1,'Position',[402 99 20 15],'String','E','Style','text');
uicontrol('Parent',h1,'Position',[356 122 20 15],'String','N','Style','text');

uicontrol('Parent',h1,'Position',[6 254 401 21],...
'BackgroundColor',[1 1 1],...
'Callback','empilhador(''edit_namesList_CB'',gcbo,guidata(gcbo))',...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','Name of an ascii file with the grids list. One grid name per row',...
'Tag','edit_namesList');

uicontrol('Parent',h1,'Position',[407 252 23 23],...
'Callback','empilhador(''push_namesList_CB'',gcbo,guidata(gcbo))',...
'TooltipString','Browse for a grids list file',...
'Tag','push_namesList');

uicontrol('Parent',h1, 'Position',[244 233 150 15],...
'Callback','empilhador(''radio_conv2netcdf_CB'',gcbo,guidata(gcbo))',...
'String','Convert to 3D netCDF file',...
'Style','radiobutton',...
'TooltipString','Take a list of files and create a single 3D netCDF file',...
'Value',1,...
'Tag','radio_conv2netcdf');

uicontrol('Parent',h1, 'Position',[244 214 150 15],...
'Callback','empilhador(''radio_conv2vtk_CB'',gcbo,guidata(gcbo))',...
'String','Convert to 3D VTK file',...
'Style','radiobutton',...
'TooltipString','Take a list of files and create a single 3D VTK file',...
'Value',0,...
'Tag','radio_conv2vtk');

uicontrol('Parent',h1, 'Position',[244 195 130 15],...
'Callback','empilhador(''radio_zonalInteg_CB'',gcbo,guidata(gcbo))',...
'String','Do zonal integration',...
'Style','radiobutton',...
'TooltipString','Take a list of files and compute a zonal average file',...
'Tag','radio_zonalInteg');

uicontrol('Parent',h1,'Position',[271 166 51 21],...
'BackgroundColor',[1 1 1],...
'Callback','empilhador(''edit_stripeWidth_CB'',gcbo,guidata(gcbo))',...
'String','0.5',...
'Style','edit',...
'TooltipString','Width of the stripe over which integration is carried on',...
'Tag','edit_stripeWidth');

uicontrol('Parent',h1,'Position',[335 169 50 15],...
'Callback','empilhador(''radio_lon_CB'',gcbo,guidata(gcbo))',...
'String','Long',...
'Style','radiobutton',...
'TooltipString','Integrate in Longitude',...
'Value',1,...
'Tag','radio_lon');

uicontrol('Parent',h1,'Position',[390 169 40 15],...
'Callback','empilhador(''radio_lat_CB'',gcbo,guidata(gcbo))',...
'String','Lat',...
'Style','radiobutton',...
'TooltipString','Integrate in Latitude',...
'Tag','radio_lat');

uicontrol('Parent',h1,'Position',[250 133 110 15],...
'Callback','empilhador(''check_region_CB'',gcbo,guidata(gcbo))',...
'String','Use sub-region?',...
'Style','checkbox',...
'TooltipString','Perform computations inside a data sub-region',...
'Tag','check_region');

uicontrol('Parent',h1,'Position',[300 104 71 21],...
'BackgroundColor',[1 1 1],...
'Callback','empilhador(''edit_north_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_north');

uicontrol('Parent',h1,'Position',[250 79 71 21],...
'BackgroundColor',[1 1 1],...
'Callback','empilhador(''edit_west_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_west');

uicontrol('Parent',h1,'Position',[350 79 71 21],...
'BackgroundColor',[1 1 1],...
'Callback','empilhador(''edit_east_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_east');

uicontrol('Parent',h1,'Position',[300 54 71 21],...
'BackgroundColor',[1 1 1],...
'Callback','empilhador(''edit_south_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_south');

uicontrol('Parent',h1,'Position',[6 9 225 236],...
'BackgroundColor',[1 1 1],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_list');

uicontrol('Parent',h1,'Position',[340 5 90 23],...
'Callback','empilhador(''push_compute_CB'',gcbo,guidata(gcbo))',...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'String','Compute',...
'Tag','push_compute');

% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
function varargout = l2_choices(varargin)
% M-File changed by desGUIDE  
 
	hObject = figure('Tag','figure1','Visible','off');
	l2_choices_LayoutFcn(hObject);
	handles = guihandles(hObject);

	handles.out.georeference = 0;
	handles.out.nearneighbor = 1;
	handles.out.mask = 0;
	handles.out.coastRes = 0;
	handles.out.quality = '';

	AllSubdatasets = varargin{1};
	for (k = numel(AllSubdatasets):-2:2)
		ind = strfind(AllSubdatasets{k}, ' ');
		subDsName = AllSubdatasets{k}(ind(1)+1:ind(2)-1);		% Get dataset name
		if (strcmp(subDsName,'qual_sst'))
			set([handles.popup_quality handles.text_quality],'Enable','on')
			ind = strfind(AllSubdatasets{k-1}, '=');
			handles.out.qualSDS = AllSubdatasets{k-1}(ind+1:end);
			break
		end
	end

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	% UIWAIT makes l2_choices wait for user response
	uiwait(handles.figure1);
	handles = guidata(handles.figure1);
	varargout{1} = handles.out;
	delete(handles.figure1),	drawnow

% -----------------------------------------------------------------------
function radio_sensor_Callback(hObject, eventdata, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_interpMin handles.radio_interpNear handles.check_landMask],'Enable','off')
	set(handles.radio_georef,'Val',0)

% -----------------------------------------------------------------------
function radio_georef_Callback(hObject, eventdata, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_interpMin handles.radio_interpNear handles.check_landMask],'Enable','on')
	set(handles.radio_sensor,'Val',0)

% -----------------------------------------------------------------------
function radio_interpNear_Callback(hObject, eventdata, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set(handles.radio_interpMin,'Val',0)

% -----------------------------------------------------------------------
function radio_interpMin_Callback(hObject, eventdata, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set(handles.radio_interpNear,'Val',0)

% -----------------------------------------------------------------------
function check_landMask_Callback(hObject, eventdata, handles)
	if ( get(hObject,'Val') )
		info = getappdata(0,'gmt_version');		% It should never be empty
		if (info.crude ~= 'y')
			errordlg('You do not have GMT installed. So, no coastlines no masking.','Error')
			set(hObject,'Val',0),	return
		elseif (info.full == 'y'),			handles.out.coastRes = 5;
		elseif (info.high == 'y'),			handles.out.coastRes = 4;
		elseif (info.intermediate == 'y'),	handles.out.coastRes = 3;
		elseif (info.low == 'y'),			handles.out.coastRes = 2;
		else,								handles.out.coastRes = 1;
		end
		guidata(handles.figure1, handles)
	end

% -----------------------------------------------------------------------
function push_OK_Callback(hObject, eventdata, handles)	
	handles.out.georeference = get(handles.radio_georef, 'Val');
	handles.out.nearneighbor = get(handles.radio_interpNear, 'Val');
	handles.out.mask = get(handles.check_landMask, 'Val');
	ind = get(handles.popup_quality, 'Val');
	if (ind ~= 1)
		handles.out.quality = 3 - ind;
	end
	guidata(handles.figure1, handles)
	uiresume(handles.figure1);

% -----------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
% The GUI is still in UIWAIT, us UIRESUME
	handles = guidata(hObject);
	handles.out = '';        % User gave up, return nothing
	guidata(handles.figure1, handles);
	uiresume(handles.figure1);

% -----------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
% Check for "escape"
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles = guidata(hObject);
		handles.out = '';    % User said no by hitting escape
		guidata(handles.figure1, handles);
		uiresume(handles.figure1);
	end

% --- Creates and returns a handle to the GUI figure. 
function l2_choices_LayoutFcn(h1)

set(h1, 'Position',[520 400 271 148],...
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
'Callback',{@l2_choices_uicallback,h1,'radio_sensor_Callback'},...
'FontName','Helvetica',...
'String','Plot in Sensor coordinates',...
'Style','radiobutton',...
'TooltipString','Plot data as it is in file (faster, but deformed)',...
'Value',1,...
'Tag','radio_sensor');

uicontrol('Parent',h1, 'Position',[10 78 231 16],...
'Callback',{@l2_choices_uicallback,h1,'radio_georef_Callback'},...
'FontName','Helvetica',...
'String','Compute georeferenced grid (takes time)',...
'Style','radiobutton',...
'TooltipString','Reinterpolate data to get a georeferenced grid. (It may take 1 minute)',...
'Tag','radio_georef');

uicontrol('Parent',h1, 'Position',[31 55 101 16],...
'Callback',{@l2_choices_uicallback,h1,'radio_interpNear_Callback'},...
'Enable','off',...
'FontName','Helvetica',...
'String','Nearneighbor',...
'Style','radiobutton',...
'TooltipString','Interpolation method',...
'Value',1,...
'Tag','radio_interpNear');

uicontrol('Parent',h1, 'Position',[31 36 121 16],...
'Callback',{@l2_choices_uicallback,h1,'radio_interpMin_Callback'},...
'Enable','off',...
'FontName','Helvetica',...
'String','Minimum curvature',...
'Style','radiobutton',...
'TooltipString','Interpolation method',...
'Tag','radio_interpMin');

uicontrol('Parent',h1, 'Position',[10 6 51 22],...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'ListboxTop',0,...
'String',{'2'; '1'; '0'},...
'Style','popupmenu',...
'TooltipString','Select the least quality level. 2 - worst - means all values. 0 - only the best',...
'Value',1,...
'Tag','popup_quality');

uicontrol('Parent',h1, 'Position',[64 10 75 15],...
'Enable','off',...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Quality factor',...
'Style','text',...
'Tag','text_quality');

uicontrol('Parent',h1, 'Position',[160 48 110 15],...
'Callback',{@l2_choices_uicallback,h1,'check_landMask_Callback'},...
'Enable','off',...
'String','Apply Land Mask',...
'Style','checkbox',...
'TooltipString','Mask reinterpolated Land pixels (need high definition coastlines installed)',...
'Tag','check_landMask');

uicontrol('Parent',h1, 'Position',[195 7 66 23],...
'Callback',{@l2_choices_uicallback,h1,'push_OK_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'String','OK',...
'Tag','push_OK');

function l2_choices_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

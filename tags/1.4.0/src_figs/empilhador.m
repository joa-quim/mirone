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
	if (nargout),	varargout{1} = hObject;     end

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
	c = false(m,1);
	caracol = false(m,1);

	n_msg = 1;							% Will hold the "change DV messages" counter
	if (n_column > 1)					% When 2nd column holds the 3D numbering
		for (k=1:m)
			if ( isempty(names{k}) ),	continue,		end		% Jump empty lines
			[t,r] = strtok(names{k});
			if (t(1) == '#'),  c(k) = true;			continue,	end
			if ( t(1) == '@')
				caracol(k) = true;
				if (~isempty(t(2:end))),		handles.changeCD_msg{n_msg} = t(2:end);		% The '@' was glued with the message
				else							handles.changeCD_msg{n_msg} = r;
				end
				n_msg = n_msg + 1;
				continue
			end

			names{k} = t;
			r = ddewhite(r);
			handles.strTimes{k} = r;
		end
	else								% Only one column with fnames
		for (k=1:m)
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
	m = length(names);			% Count remaining ones

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
		for (k=1:m)
			c(k) = (exist(names{k},'file') ~= 2);		% Flag to kill all non-existant files
		end
		names(c) = [];		handles.shortNameList(c) = [];
	end

	handles.nameList = names;
	handles.caracol = caracol;
	set(handles.edit_namesList, 'String', fname)
	set(handles.listbox_list,'String',handles.shortNameList)
	guidata(handles.figure1,handles)
	set(handles.figure1,'pointer','arrow')
	
	if (isempty(names))
		warndlg('As you may have already realized, the goodness of the name list provied is ... fiu, fiu, fiu!','Warning')
	end

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
function listbox_list_CB(hObject, handles)

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
	if (isempty(handles.nameList))
		errordlg('Yes, Com-Pute and Sem-Pute -- as you like. Empty filename list.','ERROR'),		return
	end

	got_R = false;		is_modis = false;		is_linear = false;		is_log = false;
	west = [];			east = [];		south = [];		north = [];
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
			[Z,att] =  read_gdal(FileName, [], '-U', '-C', opt_R);
		else
			Z =  read_gdal(handles.nameList{k}, att, '-U', '-C', opt_R);
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

	set(handles.figure1,'pointer','watch')
	% Read relevant metadata
	[head, opt_R, slope, intercept, base, is_modis, is_linear, is_log] = get_headerInfo(handles, got_R, west, east, south, north);

	handles.geog = 1;			handles.head = head;
	handles.was_int16 = 0;		handles.computed_grid = 0;

	if (get(handles.radio_conv2vtk,'Val'))
		fid = write_vtk(handles, grd_out, 'hdr');
	end

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
		[Z, handles.have_nans, att] = getZ(handles.nameList{k}, [], is_modis, is_linear, is_log, slope, intercept, base, opt_R);

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
	set(handles.figure1,'pointer','arrow')

% -----------------------------------------------------------------------------------------
function [head, opt_R, slope, intercept, base, is_modis, is_linear, is_log, N_spatialSize, integDim, att] = ...
			get_headerInfo(handles, got_R, west, east, south, north)
% ...

	att = read_gdal(handles.nameList{1}, [], '-M','-C');	% Do this because it deals also with ziped files

	if ( att.RasterCount == 0 && ~isempty(att.Subdatasets) && strncmp(att.DriverShortName, 'HDF4', 4) )		% Some MODIS files
		ind = strfind(att.Subdatasets{1}, '=');
		FileName = att.Subdatasets{1}(ind+1:end);		% First "ind" chars are of the form SUBDATASET_1_NAME=
		att = gdalread(FileName,'-M','-C');				% Try again (it will probabçy fail on ziped files)
	end

	% GDAL wrongly reports the corners as [0 nx] [0 ny] when no SRS
	if ( isequal([att.Corners.LR - att.Corners.UL],[att.RasterXSize att.RasterYSize]) && ~all(att.Corners.UL) )
		att.GMT_hdr(1:4) = [1 att.RasterXSize 1 att.RasterYSize];
	end

	if (get(handles.radio_lon, 'Val'))
		N_spatialSize = att.RasterYSize;		% Number of points in the spatial dim
		integDim = 2;							% Dimension along which we are going to integrate
	else
		N_spatialSize = att.RasterXSize;
		integDim = 1;
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

	if ( ~isempty(att.Metadata) && ~isempty(strfind(att.Metadata{2}, 'MODIS')) && strfind(att.Metadata{2}, 'MODIS'))
		modis_or_seawifs = true;
	elseif ( ~isempty(att.Metadata) && ~isempty(strfind(att.Metadata{2}, 'SeaWiFS')) && strfind(att.Metadata{2}, 'SeaWiFS'))
		modis_or_seawifs = true;
	else
		modis_or_seawifs = false;
	end
		
	if ( strncmp(att.DriverShortName, 'HDF4', 4) && ~isempty(att.Metadata) && modis_or_seawifs )
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
	for (k = numel(attributes):-1:1)					% Start from the bottom because they are likely close to it 
		if ( strcmpi(attributes(k).Name, what) )
			out = double(attributes(k).Value);			break
		end
	end

% -----------------------------------------------------------------------------------------
function [Z, have_nans, att] = getZ(fname, att, is_modis, is_linear, is_log, slope, intercept, base, opt_R)
% ATT may be still unknown (empty). In that case it will be returned by read_gdal()

	if (nargin < 9),	opt_R = ' ';	end
	[Z,att] = read_gdal(fname, att, '-U', '-C', opt_R);

	ind = [];
	if (is_modis)
		ind = (Z == 65535);
	elseif ( ~isempty(att.Band(1).NoDataValue) && (att.Band(1).NoDataValue == -9999) )		% TEMP -> PATHFINDER
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
	end

% -----------------------------------------------------------------------------------------
function [Z, att] = read_gdal(full_name, att, varargin)
% Help function to gdalread that deals with cases when file is compressed.
% ATT is the GDALREAD returned attributes. If empty, we'll get it here
% VARARGIN will normally contain one or more of '-U', '-C', '-M', opt_R
% WARNING: If exist(att.hdfInfo) than att.fname should exist as well (both non standard)

	str_d = [];		do_warn = 'true';		cext = [];		% Some defaults
	
	[PATH,fname,EXT] = fileparts(full_name);
	out_name = [PATH filesep fname];		% Only used if file is compressed
	if (strcmpi(EXT,'.bz2'))
		str_d = ['bzip2 -d -q -f -c ' full_name ' > ' out_name];		cext = EXT;
	elseif (strcmpi(EXT,'.zip') || strcmpi(EXT,'.gz'))
		str_d = ['gunzip -q -N -f -c ' full_name ' > ' out_name];		cext = EXT;
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
		full_name = out_name;				% The uncompressed file name
	end

	if (isempty(att))
		att = gdalread(full_name, '-M');		% This first call is used in the next test
	end

	if ( att.RasterCount == 0 && ~isempty(att.Subdatasets) && strncmp(att.DriverShortName, 'HDF4', 4) )		% Some MODIS files
		ind = strfind(att.Subdatasets{1}, '=');
		full_name = att.Subdatasets{1}(ind+1:end);		% First "ind" chars are of the form SUBDATASET_1_NAME=
	end

	% att.hdrInfo is not a default field of the ATT struct
	if (isfield(att, 'hdrInfo') && ~isempty(att.hdrInfo) && (strcmp(att.hdrInfo.SDS.Name,'sst')) )
		% Only particular case dealt now
		Z = hdf_funs('hdfread', att.fname, att.hdrInfo.SDS.Name, 'index', {[1 1],[1 1], [att.RasterYSize att.RasterXSize]});
		Z = flipud(Z);
	else
		if (nargout == 2)
			[Z, att] = gdalread(full_name, varargin{:});	% This ATT might be of a subdataset
		else
			Z = gdalread(full_name, varargin{:});
		end
	end

	if (~isempty(str_d)),	delete(out_name);	end		% Delete uncompressed file.

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
function empilhador_LayoutFcn(h1);

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
'Callback','empilhador(''listbox_list_CB'',gcbo,guidata(gcbo))',...
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

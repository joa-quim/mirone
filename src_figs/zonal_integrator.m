function varargout = zonal_integrator(varargin)
% M-File changed by desGUIDE 

	hObject = figure('Tag','figure1','Visible','off');
	zonal_integrator_LayoutFcn(hObject);
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
	if (nargout),   varargout{1} = hObject;     end

% -----------------------------------------------------------------------------------------
function edit_namesList_Callback(hObject, eventdata, handles)
    fname = get(hObject,'String');
    push_namesList_Callback([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_namesList_Callback(hObject, eventdata, handles, opt)
    if (nargin == 3)        % Direct call
        cd(handles.last_dir)
    	str1 = {'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.DAT,*.txt,*.TXT)';'*.*', 'All Files (*.*)'};
        [FileName,PathName] = uigetfile(str1,'File with grids list');
        cd(handles.home_dir);
	    if isequal(FileName,0);     return;     end
        if (PathName ~= 0),         handles.last_dir = PathName;    end
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
	fname = [PathName FileName];

    [bin,n_column] = guess_file(fname);
    % If error in reading file
    if isempty(bin)
		errordlg(['Error reading file ' fname],'Error'),	return
    end

	fid = fopen(fname);
	c = fread(fid,'*char')';      fclose(fid);
	names = strread(c,'%s','delimiter','\n');   clear c fid;
	m = length(names);
	
	handles.strTimes = [];          % To hold time steps as strings
	if (n_column > 1)
		handles.strTimes = cell(m,1);
		c = false(m,1);
		for (k=1:m)
			[t,r] = strtok(names{k});
			if (t(1) == '#'),  c(k) = true;		continue,	end
			names{k} = t;
			r = ddewhite(r);
			handles.strTimes{k} = r;
		end
		% Remove eventual commented lines
		if (any(c))
			names(c) = [];			handles.strTimes(c) = [];
			m = length(names);      % Count remaining ones
		end
	end

	handles.shortNameList = cell(m,1);      % To hold grid names with path striped
	c = false(m,1);
	for (k=1:m)
		if ( isempty(names{k}) ),	continue,		end		% Jump empty lines
		if (n_column == 1 && names{k}(1) == '#')    % If n_column > 1, this test was already done above
			c(k) = true;    continue;
		end
		[PATH,FNAME,EXT] = fileparts(names{k});
		if (isempty(PATH))
			handles.shortNameList{k} = names{k};
			names{k} = [PathName names{k}];
		else
			handles.shortNameList{k} = [FNAME EXT];
		end
		if (any(c))
			names(c) = [];		handles.shortNameList(c) = [];
		end
	end

	% Check that at least the files in provided list do exist
	c = false(m,1);
	for (k=1:m)
		c(k) = (exist(names{k},'file') ~= 2);
	end
	names(c) = [];      handles.shortNameList(c) = [];

	handles.nameList = names;
	set(handles.edit_namesList, 'String', fname)
	set(handles.listbox_list,'String',handles.shortNameList)
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function radio_conv2netcdf_Callback(hObject, eventdata, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.edit_stripeWidth handles.radio_lon handles.radio_lat],'Enable','off')
	set(handles.radio_zonalInteg,'Val',0)

% -----------------------------------------------------------------------------------------
function radio_zonalInteg_Callback(hObject, eventdata, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.edit_stripeWidth handles.radio_lon handles.radio_lat],'Enable','on')
	set(handles.radio_conv2netcdf,'Val',0)

% -----------------------------------------------------------------------------------------
function check_region_Callback(hObject, eventdata, handles)
	if (get(hObject,'Val'))
		set([handles.edit_north handles.edit_south handles.edit_west handles.edit_east],'Enable','on')
	else
		set([handles.edit_north handles.edit_south handles.edit_west handles.edit_east],'Enable','off')
	end

% -----------------------------------------------------------------------------------------
function edit_north_Callback(hObject, eventdata, handles)
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
function edit_south_Callback(hObject, eventdata, handles)
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
function edit_west_Callback(hObject, eventdata, handles)
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
function edit_east_Callback(hObject, eventdata, handles)
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
function listbox_list_Callback(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------------
function edit_stripeWidth_Callback(hObject, eventdata, handles)
	x1 = str2double(get(hObject,'String'));
	if (isnan(x1)),		set(hObject,'String','0.5'),	end

% -----------------------------------------------------------------------------------------
function radio_lat_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_lon,'Val',0)

% -----------------------------------------------------------------------------------------
function radio_lon_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_lat,'Val',0)

% -----------------------------------------------------------------------------------------
function push_compute_Callback(hObject, eventdata, handles)
% ...
	if (isempty(handles.nameList))
		errordlg('Yes, ComPute and SemPute. Empty filename list.','ERROR'),		return
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
	
	% HDF to netCDF conversion?
	if (get(handles.radio_conv2netcdf,'Val'))
		cut2cdf(handles, got_R, west, east, south, north)
		return
	end

	% MORE OR LESS DEPRECATED CODE (IT WILL PROBABLY FAIL) - FUNCTIONALITY MOVED TO AQUAMOTO SUPP FUNS
	[head, opt_R, slope, intercept, base, is_modis, is_linear, is_log, N_spatialSize, integDim] = ...
			get_headerInfo(handles, got_R, west, east, south, north);

	att = read_gdal(handles.nameList{1},'-M','-C');

	if ( strcmp(att.DriverShortName, 'HDF4') && att.RasterCount == 0 && ~isempty(att.Subdatasets) )
		ind = strfind(att.Subdatasets{1}, '=');
		FileName = att.Subdatasets{1}(ind+1:end);		% First "ind" chars are of the form SUBDATASET_1_NAME=
		att = read_gdal(FileName,'-M','-C');			% Try again
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
			[Z,att] =  read_gdal(FileName, '-U', '-C', opt_R);
		else
			Z =  read_gdal(handles.nameList{k}, '-U', '-C', opt_R);
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
	txt1 = 'netCDF grid format (*.nc,*.grd)';
	txt2 = 'Select output netCDF grid';
	[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put');
	if isequal(FileName,0),		return,		end
	[pato, fname, EXT] = fileparts(FileName);
	if (isempty(EXT)),		FileName = [fname '.nc'];	end
	grd_out = [PathName FileName];

	set(handles.figure1,'pointer','watch')
	% Read relevant metadata
	[head, opt_R, slope, intercept, base, is_modis, is_linear, is_log] = get_headerInfo(handles, got_R, west, east, south, north);

	handles.geog = 1;
	handles.head = head;
	handles.was_int16 = 0;
	handles.computed_grid = 0;
	
	nSlices = numel(handles.nameList);
	for (k = 1:nSlices)
		set(handles.listbox_list,'Val',k),		pause(0.01)			% Show advance

		[Z,att] = read_gdal(handles.nameList{k}, '-U', '-C', opt_R);

		ind = [];
		if (is_modis)
			ind = (Z == 65535);
		elseif ( ~isempty(att.Band(1).NoDataValue) && (att.Band(1).NoDataValue == -9999) )		% TEMP -> PATHFINDER
			ind = (Z == 0);
		elseif ( ~isempty(att.Band(1).NoDataValue) && ~isnan(att.Band(1).NoDataValue) )
			ind = (Z == (att.Band(1).NoDataValue));
		elseif (isnan(att.Band(1).NoDataValue))		% The nodata is NaN, replace NaNs in Z by zero
			ind = isnan(Z);
		end

		% See if we must apply a scaling equation
		if (is_linear)
			Z = Z * slope + intercept;
		elseif (is_log)
			Z = base .^ (Z * slope + intercept);
		end
		Z = single(Z);
		handles.have_nans = 0;
		if (~isempty(ind)),		Z(ind) = NaN;	handles.have_nans = 1;		end
		%Z(Z > 5) = 5;			% <==== CLIPING
		zz = grdutils(Z,'-L');  handles.head(5:6) = [zz(1) zz(2)];

		% Must treate compiled version differently since, biggest of misteries, nc_funs
		% than hangs when writting unlimited variables
		if (~handles.IamCompiled)
			if (~isempty(handles.strTimes)),		t_val = handles.strTimes{k};
			else									t_val = sprintf('%d',k - 1);
			end
			if (k == 1)
				nc_io(grd_out, ['w-' t_val '/time'], handles, reshape(Z,[1 size(Z)]))
			else
				nc_io(grd_out, sprintf('w%d\\%s', k-1, t_val), handles, Z)
			end
		else
			if (k == 1)
				handles.levelVec = str2double(handles.strTimes);
     			nc_io(grd_out,sprintf('w%d/time',nSlices), handles, reshape(Z,[1 size(Z)]))
			else
				nc_io(grd_out, sprintf('w%d', k-1), handles, Z)
			end
		end
	end
	set(handles.listbox_list,'Val',1)
	set(handles.figure1,'pointer','arrow')
	
% -----------------------------------------------------------------------------------------
function [head, opt_R, slope, intercept, base, is_modis, is_linear, is_log, N_spatialSize, integDim] = ...
			get_headerInfo(handles, got_R, west, east, south, north)
% ...
	opt_R = ' ';	is_modis = false;		is_linear = false;		is_log = false;		base = 0;

	att = read_gdal(handles.nameList{1},'-M','-C');

	if ( att.RasterCount == 0 && ~isempty(att.Subdatasets) && strcmp(att.DriverShortName, 'HDF4') )		% Some MODIS files
		ind = strfind(att.Subdatasets{1}, '=');
		FileName = att.Subdatasets{1}(ind+1:end);		% First "ind" chars are of the form SUBDATASET_1_NAME=
		att = gdalread(FileName,'-M','-C');				% Try again
	end

	% GDAL wrongly reports the corners as [0 nx] [0 ny] when no SRS
	if ( isequal([att.Corners.LR - att.Corners.UL],[att.RasterXSize att.RasterYSize]) && ~all(att.Corners.UL) )
		att.GMT_hdr(1:4) = [1 att.RasterXSize 1 att.RasterYSize];
	end
	head = att.GMT_hdr;

	if (get(handles.radio_lon, 'Val'))
		N_spatialSize = att.RasterYSize;		% Number of points in the spatial dim
		integDim = 2;							% Dimension along which we are going to integrate
	else
		N_spatialSize = att.RasterXSize;
		integDim = 1;
	end

	if (~isempty(strfind(att.Metadata{2}, 'MODIS')) && strfind(att.Metadata{2}, 'MODIS'))
		modis_or_seawifs = true;
	elseif (~isempty(strfind(att.Metadata{2}, 'SeaWiFS')) && strfind(att.Metadata{2}, 'SeaWiFS'))
		modis_or_seawifs = true;
	else
		modis_or_seawifs = false;
	end
		
	if ( strncmp(att.DriverShortName, 'HDF4', 4) && ~isempty(att.Metadata) && modis_or_seawifs )
		y_max = str2double(att.Metadata{39}(23:end));		% att.Metadata{39} -> Northernmost Latitude=90
		y_min = str2double(att.Metadata{40}(23:end));		% att.Metadata{40} -> Southernmost Latitude=90
		x_min = str2double(att.Metadata{41}(23:end));		% att.Metadata{41} -> Westernmost Longitude=-180
		x_max = str2double(att.Metadata{42}(23:end));		% att.Metadata{41} -> Easternmost Longitude=-180
		columns = str2double(att.Metadata{49}(19:end));		% att.Metadata{49} -> Number of Columns=4320
		dx = (x_max - x_min) / columns;
		x_min = x_min + dx/2;		x_max = x_max - dx/2;	% Orig data was pixel registered
		y_min = y_min + dx/2;		y_max = y_max - dx/2;
		head(1:4) = [x_min x_max y_min y_max];
		head(8:9) = dx;
		if (got_R)			% We must give the region in pixels since the image is trully not georeferenced
			rows = str2double(att.Metadata{48}(17:end));	% att.Metadata{48} -> Number of Lines=2160
			cp = round(([west east] - x_min) / dx);
			rp = round(([south north] - y_min) / dx);
			rp(rp < 0) = 0;			cp(cp < 0) = 0;
			head(1) = x_min + cp(1)*dx;		head(2) = x_min + cp(2)*dx;
			head(3) = y_min + rp(1)*dx;		head(4) = y_min + rp(2)*dx;
			rp = rows - rp -1;		rp = [rp(2) rp(1)];
			opt_R = sprintf('-r%d/%d/%d/%d',cp(1:2),rp(1:2));
			if (get(handles.radio_lon, 'Val')),		N_spatialSize = round(diff(rp) + 1);
			else									N_spatialSize = round(diff(cp) + 1);
			end
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
	else							% TEMP -> SST PATHFINDER
		finfo = hdfinfo(handles.nameList{1});
		slope = double(finfo.SDS.Attributes(11).Value);		% = 0.075;
		intercept= double(finfo.SDS.Attributes(12).Value);	% = -3.0;
		lat = finfo.SDS.Dims(1).Scale;		% Get the latitudes
		lon = finfo.SDS.Dims(2).Scale;		% Get the longitudes
		x_min = lon(1);			x_max = lon(end);
		y_min = lat(end);		y_max = lat(1);
		head(1:4) = [x_min x_max y_min y_max];
		dx = lon(3) - lon(2);
		head(8:9) = dx;
		if (got_R)			% We must give the region in pixels since the image is trully not georeferenced
			rows = finfo.SDS.Dims(1).Size;					% Number of rows
			cp = round(([west east] - x_min) / dx);
			rp = round(([south north] - y_min) / dx);
			rp(rp < 0) = 0;			cp(cp < 0) = 0;
			head(1) = x_min + cp(1)*dx;		head(2) = x_min + cp(2)*dx;
			head(3) = y_min + rp(1)*dx;		head(4) = y_min + rp(2)*dx;
			rp = rows - rp -1;		rp = [rp(2) rp(1)];
			opt_R = sprintf('-r%d/%d/%d/%d',cp(1:2),rp(1:2));
			if (get(handles.radio_lon, 'Val')),		N_spatialSize = round(diff(rp) + 1);
			else									N_spatialSize = round(diff(cp) + 1);
			end
		end
		head(7) = 0;		% Make sure that grid reg is used
		is_linear = true;
	end

% -----------------------------------------------------------------------------------------
function [Z, att] = read_gdal(full_name, varargin)
% Help function to gdalread that deals with cases when file is compressed.
% VARARGIN will normally contain one or more of '-U', '-C', '-M', opt_R

	str_d = [];		do_warn = 'true';		cext = [];
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
			return
		end
		if (do_warn),	aguentabar(1,'title','Donne'),		end
		full_name = out_name;				% The uncompressed file name
	end

	att = gdalread(full_name, '-M');		% This first call serves to be used in the next test
	if ( att.RasterCount == 0 && ~isempty(att.Subdatasets) && strcmp(att.DriverShortName, 'HDF4') )		% Some MODIS files
		ind = strfind(att.Subdatasets{1}, '=');
		full_name = att.Subdatasets{1}(ind+1:end);		% First "ind" chars are of the form SUBDATASET_1_NAME=
	end

	if (nargout == 2)
		[Z, att] = gdalread(full_name, varargin{:});
		Z = double(Z);
	else
		Z = gdalread(full_name, varargin{:});
		if ( ~isa(Z,'struct') )			% Get Z, otherwise get only the metadata
			Z = double(Z);
		end
	end

	if (~isempty(str_d)),	delete(out_name);	end		% Delete uncompressed file.
	
% ---------------------------------------------------------------------
function zonal_integrator_LayoutFcn(h1);

set(h1,'Position',[520 532 440 280],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Zonal Integrator',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[241 175 30 15],...
'String','Delta',...
'Style','text');

uicontrol('Parent',h1,'Position',[240 40 191 117],...
'Style','frame',...
'Tag','frame1');

uicontrol('Parent',h1,'Position',[358 40 20 15],...
'String','S','Style','text');

uicontrol('Parent',h1,'Position',[250 98 20 15],...
'String','W','Style','text');

uicontrol('Parent',h1,'Position',[402 99 20 15],...
'String','E','Style','text');

uicontrol('Parent',h1,'Position',[356 122 20 15],...
'String','N','Style','text');

uicontrol('Parent',h1,'Position',[6 254 401 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@zonal_integrator_uicallback,h1,'edit_namesList_Callback'},...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','Name of an ascii file with the grids list. One grid name per row',...
'Tag','edit_namesList');

uicontrol('Parent',h1,'Position',[407 252 23 23],...
'Callback',{@zonal_integrator_uicallback,h1,'push_namesList_Callback'},...
'TooltipString','Browse for a grids list file',...
'Tag','push_namesList');

uicontrol('Parent',h1,'Position',[271 171 51 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@zonal_integrator_uicallback,h1,'edit_stripeWidth_Callback'},...
'String','0.5',...
'Style','edit',...
'TooltipString','Width of the stripe over which integration is carried on',...
'Tag','edit_stripeWidth');

uicontrol('Parent',h1,'Position',[390 174 40 15],...
'Callback',{@zonal_integrator_uicallback,h1,'radio_lat_Callback'},...
'String','Lat',...
'Style','radiobutton',...
'TooltipString','Integrate in Latitude',...
'Tag','radio_lat');

uicontrol('Parent',h1,'Position',[335 174 50 15],...
'Callback',{@zonal_integrator_uicallback,h1,'radio_lon_Callback'},...
'String','Long',...
'Style','radiobutton',...
'TooltipString','Integrate in Longitude',...
'Value',1,...
'Tag','radio_lon');

uicontrol('Parent',h1,'Position',[250 133 110 15],...
'Callback',{@zonal_integrator_uicallback,h1,'check_region_Callback'},...
'String','Use sub-region?',...
'Style','checkbox',...
'TooltipString','Perform computations inside a data sub-region',...
'Tag','check_region');

uicontrol('Parent',h1,'Position',[300 104 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@zonal_integrator_uicallback,h1,'edit_north_Callback'},...
'Enable','off',...
'Style','edit',...
'Tag','edit_north');

uicontrol('Parent',h1,'Position',[250 79 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@zonal_integrator_uicallback,h1,'edit_west_Callback'},...
'Enable','off',...
'Style','edit',...
'Tag','edit_west');

uicontrol('Parent',h1,'Position',[350 79 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@zonal_integrator_uicallback,h1,'edit_east_Callback'},...
'Enable','off',...
'Style','edit',...
'Tag','edit_east');

uicontrol('Parent',h1,'Position',[300 54 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@zonal_integrator_uicallback,h1,'edit_south_Callback'},...
'Enable','off',...
'Style','edit',...
'Tag','edit_south');

uicontrol('Parent',h1,'Position',[6 9 225 236],...
'BackgroundColor',[1 1 1],...
'Callback',{@zonal_integrator_uicallback,h1,'listbox_list_Callback'},...
'Style','listbox',...
'Value',1,...
'Tag','listbox_list');

uicontrol('Parent',h1,'Position',[340 5 90 23],...
'Callback',{@zonal_integrator_uicallback,h1,'push_compute_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'String','Compute',...
'Tag','push_compute');

uicontrol('Parent',h1, 'Position',[244 210 150 15],...
'Callback',{@zonal_integrator_uicallback,h1,'radio_conv2netcdf_Callback'},...
'String','Convert to 3D netCDF file',...
'Style','radiobutton',...
'TooltipString','Take a list of .hdf files and create a single 3D netCDF file',...
'Value',1,...
'Tag','radio_conv2netcdf');

uicontrol('Parent',h1, 'Position',[244 229 130 15],...
'Callback',{@zonal_integrator_uicallback,h1,'radio_zonalInteg_Callback'},...
'String','Do zonal integration',...
'Style','radiobutton',...
'TooltipString','Take a list of .hdf files and compute a zonal average file',...
'Tag','radio_zonalInteg');

function zonal_integrator_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

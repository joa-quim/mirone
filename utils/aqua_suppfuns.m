function varargout = aqua_suppfuns(opt, varargin)
% Supplement functions to allow using Aquamoto with plain netCDF coards grids

%	Copyright (c) 2004-2013 by J. Luis
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

	switch opt
		case 'coards_hdr',		[varargout{1:nargout}] = init_header_params(varargin{:});
		case 'coards_slice',	coards_sliceShow(varargin{:})
		case 'forGDAL_hdr',		[varargout{1:nargout}] = init_header_gdal(varargin{:});
		case 'forGDAL_slice',	gdal_sliceShow(varargin{:})
	end

% --------------------------------------------------------------------------
function out = init_header_params(handles,X,Y,head,misc,getAllMinMax)
% Use the OUT option when using this function to get several usefull info 
% about the netCDF file but NOT using a GUI.
%
% The 'getAllMinMax' when set to TRUE (or not provided) will scan the entire file to
% compute all individual layers 'actual_range'. However, this may slow quite a bit the
% loading of big files and apparently is not very used/useful because most of the time
% is spent loading the layer and the min/max computation is extremely fast.

	if (nargin < 6),		getAllMinMax = true;		end
	handles.x = X;		handles.y = Y;
	handles.time = [];
	handles.number_of_timesteps = misc.z_dim(1);		% ... NEEDS THINKING
	
	handles.x_min = head(1);			handles.x_max = head(2);
	handles.y_min = head(3);			handles.y_max = head(4);

	% ------------- Finish slider configurations -------------
	s = handles.nc_info;
	if (handles.number_of_timesteps > 1)
		st = [1 10] / (handles.number_of_timesteps - 1);
		id = find(strcmpi('time',{s.Dataset.Name}));			% ONLY WHEN 3RTH DIM IS CALLED time
		if (~isempty(id))
			handles.time = double(nc_funs('varget', handles.fname, s.Dataset(id).Name));
		else
			handles.time = 1:handles.number_of_timesteps;
		end
		slMax = handles.number_of_timesteps;
	else
		slMax = 1+eps;	st = [1 1];		handles.time = 1;		% Defaults for no crashing
	end

	% ------ Compute individual and global min/maxs ----------------------------------
	handles.zMinMaxs = zeros(handles.number_of_timesteps,2);
	if (handles.IamTSU)			% This one needs special doing since we want min/max over ocean only
		aguentabar(0,'title','Computing water surface min/max')
		zBat = nc_funs('varget', handles.fname, 'bathymetry');
		ind = (zBat < 0);	clear zBat
		for (k = 1:handles.number_of_timesteps)
			Z = nc_funs('varget', handles.fname, s.Dataset(misc.z_id).Name, [(k-1) 0 0], [1 s.Dataset(misc.z_id).Size(end-1:end)]);
			Z = Z(ind);
			indNaN = isnan(Z);
			if (any(indNaN)),	Z(indNaN) = [];	end
			if (~isempty(Z))
				handles.zMinMaxs(k,:) = [min(Z) max(Z)];
			end
			aguentabar(k/handles.number_of_timesteps);
		end
		handles.zMinMaxsGlobal = [min(handles.zMinMaxs(:,1)) max(handles.zMinMaxs(:,2))];
		k = 1;
		while ((k <= handles.number_of_timesteps) && isequal(handles.zMinMaxs(k,:), [0 0]))
			k = k + 1;
		end
		if (k == handles.number_of_timesteps + 1)
			k = 1;		% To not error below
			warndlg('WARNING: All layers are zero. You better go home.', 'Warning')
		end
		head(5:6) = handles.zMinMaxs(k,:);			% Take the first non zero slice min/max
	elseif (getAllMinMax)
		aguentabar(0,'title','Computing global min/max')
		for (k = 1:handles.number_of_timesteps)
			Z = nc_funs('varget', handles.fname, s.Dataset(misc.z_id).Name, [(k-1) 0 0], [1 s.Dataset(misc.z_id).Size(end-1:end)]);
			if ( isa(Z, 'single') )			% min/max are bugged when NaNs in singles
				zz = grdutils(Z,'-L');
				handles.zMinMaxs(k,:) = [zz(1) zz(2)];
			else
				handles.zMinMaxs(k,:) = [double(min(Z(:))) double(max(Z(:)))];
			end
			aguentabar(k/handles.number_of_timesteps);
		end
		handles.zMinMaxsGlobal = [min(handles.zMinMaxs(:,1)) max(handles.zMinMaxs(:,2))];
		while (isequal(handles.zMinMaxs(k,:), [0 0]) && (k <= handles.number_of_timesteps))
			k = k + 1;
		end
		if (k > handles.number_of_timesteps)
			warndlg('WARNING: All layers are zero. You better go home.', 'Warning')
		end
		head(5:6) = handles.zMinMaxs(k,:);			% Take the first non zero slice min/max
	else
		ind = strcmp({s.Dataset(misc.z_id).Attribute.Name},'actual_range');
		if (any(ind))
			handles.zMinMaxsGlobal = s.Dataset(misc.z_id).Attribute(ind).Value;
		else
			warndlg('Non CF compliant file. Missing ''actual_range'' attribute. Expect color screws.','Warnerror')
			handles.zMinMaxsGlobal = [0 1];
		end
		head(5:6) = handles.zMinMaxsGlobal;
	end
	handles.minWater = handles.zMinMaxsGlobal(1);
	handles.maxWater = handles.zMinMaxsGlobal(2);
	handles.geog = aux_funs('guessGeog',head(1:4));
	% ---------------------------------------------------------------------------------

	if (~handles.IamTSU)
		handles.cmapLand = jet(256);		% Reset the default colormap (default's Aquamoto is a specific one)
	end
	handles.head = head;
	handles.illumComm = [];					% New file. Reset illum state.
	handles.imgBat = [];
	handles.netcdf_z_id = misc.z_id;
	handles.is_coards = true;

	% -------------------- See if we have a projection ----------------------------------
	if (~isempty(misc.strPROJ4)),	handles.strPROJ4 = misc.strPROJ4;
	else							handles.strPROJ4 = [];
	end
	if (~isempty(misc.srsWKT)),		handles.srsWKT = misc.srsWKT;
	else							handles.srsWKT = [];
	end

	if (nargout)
		out = handles;
	else
		set( handles.edit_Ncols,'String',sprintf('%d',misc.z_dim(end)) )
		set( handles.edit_Nrows,'String',sprintf('%d',misc.z_dim(end-1)) )
		set(handles.slider_layer,'Min',1,'Max',slMax,'Val',1,'SliderStep',st) 	
		set(handles.edit_globalWaterMin,'String',handles.zMinMaxsGlobal(1))
		set(handles.edit_globalWaterMax,'String',handles.zMinMaxsGlobal(2))
		set(handles.hTabAnuga,'String','netCDF')
		set_common(handles, handles.head)
		guidata(handles.figure1,handles)
		% Store the nc_info and z_id in Mirone handles so we can access it from there as well
		if (~isempty(handles.hMirFig))	% This is another patch for the mess of the hMirFig existing or not
			hFig = handles.hMirFig;		% and how the code flow in places relies on that.
		else
			hFig = handles.hMirFig_safe;
		end
		handMir = guidata(hFig);
		handMir.netcdf_z_id = misc.z_id;
		handMir.nc_info = handles.nc_info;
		handMir.time_z = handles.time;
		guidata(hFig, handMir)
	end

% --------------------------------------------------------------------------
function out = init_header_gdal(handles)
% Read a multiband file with gdal and fill the header parameters

	handles.illumComm = [];					% New file. Reset illum state.
	handles.imgBat = [];
	handles.is_otherMultiband = true;
	att = gdalread(handles.fname,'-M','-C');
	X = linspace(att.GMT_hdr(1), att.GMT_hdr(2), att.RasterXSize);
	Y = linspace(att.GMT_hdr(3), att.GMT_hdr(4), att.RasterYSize);
	handles.number_of_timesteps = att.RasterCount;
	handles.time = 1:att.RasterCount;
	handles.head = att.GMT_hdr;
	handles.x_min = X(1);		handles.x_max = X(2);
	handles.y_min = Y(1);		handles.y_max = Y(2);
	handles.x = X;				handles.y = Y;
	handles.flip_on_read = true;
	if (isempty(att.GeoTransform)),		handles.flip_on_read = false;	end

	st = [1 10] / (att.RasterCount - 1);

	% ------ Compute individual and global min/maxs ------------------------
	handles.zMinMaxs = zeros(att.RasterCount, 2);
	for (k = 1:att.RasterCount)
		if ( isnan(att.Band(1).MinMax(1)) )		% Shit, nothing usable here
			if (~nargout),	set(handles.check_globalMinMax,'Enable', 'off'),	end
			break
		end
		handles.zMinMaxs(k,:) = att.Band(k).MinMax(:);
	end
	handles.zMinMaxsGlobal = [min(handles.zMinMaxs(:,1)) max(handles.zMinMaxs(:,2))];
	handles.minWater = handles.zMinMaxsGlobal(1);
	handles.maxWater = handles.zMinMaxsGlobal(2);
	handles.geog = aux_funs('guessGeog',att.GMT_hdr(1:4));
	% ---------------------------------------------------------------------------------

	handles.cmapLand = jet(256);			% Reset the default colormap (default's Aquamoto is a specific one)

	% -------------------- See if we have a projection ----------------------------------
	if (~isempty(att.ProjectionRef)),	handles.srsWKT = att.ProjectionRef;
	else								handles.srsWKT = [];
	end
	handles.strPROJ4 = [];

	if (nargout)
		out = handles;
	else
		set(handles.slider_layer,'Min',1,'Max',att.RasterCount,'Val',1,'SliderStep',st)
		set( handles.edit_Ncols,'String',sprintf('%d',att.RasterXSize) )
		set( handles.edit_Nrows,'String',sprintf('%d',att.RasterYSize) )
		set( handles.check_splitDryWet,'Enable', 'off' )
		set( handles.push_runIn,'Enable', 'off' )
		set( handles.slider_transparency,'Enable', 'off' )
		set(handles.edit_globalWaterMin,'String',handles.zMinMaxsGlobal(1))
		set(handles.edit_globalWaterMax,'String',handles.zMinMaxsGlobal(2))
		set(handles.hTabAnuga,'String','GDALish')
		set_common(handles, handles.head)
		guidata(handles.figure1,handles)
	end

% --------------------------------------------------------------------------------------------
function gdal_sliceShow(handles, att)
% Read the slice with gdalread and send in the array to coards_sliceShow() to do the rest

	opt_U = ' ';
	if (handles.flip_on_read),		opt_U = '-U';	end
	opt_B = sprintf('-B%d', handles.sliceNumber + 1);
	Z = gdalread(handles.fname, opt_B, opt_U, '-C');
	if (isa(Z, 'uint16') || isa(Z, 'int16') || isa(Z, 'double'))
		Z = single(Z);
	end
	coards_sliceShow(handles, Z)

% --------------------------------------------------------------------------------------------
function coards_sliceShow(handles, Z)
% ...

	if (nargin == 1)		% Otherwise we suposedly already know Z (from gdalread)
		if (isempty(handles.fname))
			errordlg('Hey Lou. What about a walk on the Wild Side? Maybe you''ll find a little file there that you can use here!','Chico clever')
			return
		end

		z_id = handles.netcdf_z_id;
		s = handles.nc_info;					% Retrieve the .nc info struct 
		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [handles.sliceNumber 0 0], [1 s.Dataset(z_id).Size(end-1:end)]);
		dims = s.Dataset(z_id).Dimension;		% Variable names of dimensions z variable
		y_id = find(strcmp({s.Dataset.Name}, dims{end-1}));
		if (~isempty(y_id))						% Check if we need to flipud(Z)
			Y = double(nc_funs('varget', handles.fname, s.Dataset(y_id).Name));
			if (Y(2) < Y(1)),	Z = flipud(Z);	end
		end
		try		handles.head(5:6) =  handles.zMinMaxs(handles.sliceNumber,:);	end		%  handles.head was from 1st layer
	end

	have_nans = 0;
	if (isa(Z,'single') && ~handles.IamTSU)		% Actually I don't know if this is used in any circumnstances.
		have_nans = grdutils(Z,'-N');			% No worry, very fast
	end
	if (have_nans && handles.useLandPhoto)
		alphaMask = alloc_mex(size(Z),'uint8');	% Create an image mask of Dry/Wets
		alphaMask(~isnan(Z)) = 255;				% nan pixeis will be transparent
	end

	try			% MUST BE TURN INTO A PROPER TEST (AQUAMOTO OR NOT-AQUAMOTO)
		splitDryWet = get(handles.check_splitDryWet,'Val');		% See if we need to build wet and dry images, or only one
	catch
		splitDryWet = false;
		handles.cmapLand = jet(256);			% Was not yet deffined if a TSU file opened in SLICES
	end

	if (splitDryWet && handles.IamTSU)
		zBat = nc_funs('varget', handles.fname, 'bathymetry');
		dife = cvlib_mex('absDiff', zBat, Z);
		indLand = (dife < 1e-2);				% The 1e-2 SHOULD be parameterized
		Z(indLand) = 0;							% Smash these to to steal the color dynamics

		if (handles.useLandPhoto)
			alfa = 255;		% Means land will be completely opac and water 100% transparent
			if (~isempty(handles.hMirFig) && ishandle(handles.handMir.hImg))	% A Mirone figure with image already exists
				alphaMask = get(handles.handMir.hImg,'AlphaData');
				if (numel(alphaMask) == 1)		% No AlphaMask yet, but we'll need one further down
					alphaMask = alloc_mex(size(indLand),'uint8');	% Create an image mask of Dry/Wets
					alphaMask(~indLand) = alfa;
				else							% AlphaMask exists, but we need to update it to reflect this slice water level
					alphaMask(indLand) = 0;		alphaMask(~indLand) = alfa;
				end
			else								% The Mirone figure will be created later. Compute the AlphaMask
				alphaMask = alloc_mex(size(indLand),'uint8');	% Create an image mask of Dry/Wets
				alphaMask(~indLand) = alfa;
			end
		end
	end

	% ----- Open or update a Mirone window with the slice display ----
	if (isempty(handles.hMirFig) || ~ishandle(handles.hMirFig))			% First run or killed Mirone window
		tmp.X = handles.x;		tmp.Y = handles.y;		tmp.head = handles.head;	tmp.cmap = handles.cmapLand;
		tmp.name = sprintf('Layer = %g',handles.time(handles.sliceNumber+1));
		if (~isempty(handles.srsWKT)),		tmp.srsWKT = handles.srsWKT;	end
		hFigs = findobj(0,'type','figure');
		if (numel(hFigs) == 2)	% Often we have an empty Mir fig but that is very difficult to use here. So blow it
			inds = [isempty(getappdata(hFigs(1), 'IAmAMirone')) isempty(getappdata(hFigs(2), 'IAmAMirone'))];
			hFigs = hFigs(~inds);							% Only one of them is a Mirone fig
			handThis = guidata(hFigs);
			if (~handThis.validGrid),		delete(hFigs),	clear handThis,	end
		end
		
		reset = false;
		if (~ishandle(handles.hMirFig))						% Save this before handles.hMirFig is reset
			nLayers = handles.handMir.nLayers;
			reset = true;
		end
		
		handles.hMirFig = mirone(Z, tmp);
		move2side(handles.figure1,handles.hMirFig,'left')
		handles.handMir = guidata(handles.hMirFig);			% Get the handles of the now existing Mirone fig

		if (reset)											% Figure was killed but the folowing info exists and is needed
			handles.handMir.nLayers = nLayers;
			handles.handMir.netcdf_z_id = handles.netcdf_z_id;
			handles.handMir.nc_info = handles.nc_info;
			handles.handMir.time_z = handles.time;
			handles.handMir.grdname = handles.fname;		% At least grid_profiler() needs this for the 3D interpolations
		end

		handles.firstLandPhoto = true;
		if (handles.useLandPhoto)
			h = image('XData',handles.geoPhotoX,'YData',handles.geoPhotoY, 'CData',handles.geoPhoto, 'Parent',handles.handMir.axes1);
			uistack(h,'bottom')
			handles.firstLandPhoto = false;
			set(handles.handMir.hImg,'AlphaData',alphaMask)	% 'alphaMask' was updated ... maybe somewhere
		end
		handles.imgBat = [];				% Make sure this one is always reset when user kills fig

	else									% We already have a Mirone image. Update it with this new slice
		handles.handMir = guidata(handles.hMirFig);			% Get updated handles to see if illum has changed
		if (~(isa(Z,'uint8') || isa(Z,'int8')))
			setappdata(handles.handMir.figure1,'dem_z',Z);	% Update grid so that coursor display correct values
		end													% Have to do it here because minmax arg to scalet8 CHANGES Z

		if (handles.useLandPhoto)							% External Land image
			if (handles.firstLandPhoto)						% First time, create the background image
				h = image('XData',handles.geoPhotoX,'YData',handles.geoPhotoY, 'CData',handles.geoPhoto, 'Parent',handles.handMir.axes1);
				uistack(h,'bottom')
				handles.firstLandPhoto = false;
			end
			set(handles.handMir.hImg,'AlphaData',alphaMask)	% 'alphaMask' was updated ... somewhere
		end

		if (~get(handles.check_globalMinMax, 'Val')),		minmax = [];		% Use Slice's min/max
		else							minmax = handles.zMinMaxsGlobal;
		end
		if (isa(Z,'int8') && ~isempty(minmax)),	minmax = [];	end	% We don't want to scale a 1 byte array

		if (~isa(Z,'uint8'))
			if (~isempty(minmax)),		img = scaleto8(Z, 8, minmax);
			else						img = scaleto8(Z);
			end
		else
			img = Z;
		end

		if (handles.IamTSU && splitDryWet)
			recomp = false;
			if (~isempty(handles.landIllumComm_bak) && ~isequal(handles.landIllumComm_bak, handles.landIllumComm))
				recomp = true;
			end
			if (recomp || isempty(handles.imgBat))		% First time, compute it (not necessarily shaded)
				% Put the cmap discontinuity at the zero of bat
				zz = grdutils(zBat,'-L');	% I'm lazy to fish this inside the s info struct
				head = handles.head;	head(5:6) = [zz(1) zz(2)];
				handles.cmapBat = aquamoto('makeCmapBat', handles, head, handles.cmapLand, 1);
				handles.imgBat  = ind2rgb8(scaleto8(zBat), handles.cmapBat);
				if (get(handles.radio_shade, 'Val'))
					R = aquamoto('illumByType', handles, zBat, head, handles.landIllumComm);
					handles.imgBat = shading_mat(handles.imgBat, R, 'no_scale');
				end
			end
			img = aquamoto('do_imgWater', handles, 1, Z, handles.imgBat, indLand);

		elseif (get(handles.radio_shade, 'Val'))
			indVar = 1;									% FAR FROM SURE THAT THIS IS CORRECT
			img = ind2rgb8(img, handles.cmapLand);		% img is now RGB
			head = handles.head;
			if (~isempty(handles.ranges{indVar})),		head(5:6) = handles.ranges{indVar};		end
			if (head(5) == 0 && head(6) == 0)
				warndlg('Something screwed up. Don''t know this grid min/max and must comput it now.','Warning')
				zz = grdutils(Z, '-L');		head(5:6) = [zz(1) zz(2)];
			end
			R = aquamoto('illumByType', handles, Z, head, handles.landIllumComm);
			img = shading_mat(img,R,'no_scale');		% and now it is illuminated
		end

		set(handles.handMir.hImg, 'CData', img)
		set(handles.handMir.figure1, 'Name', sprintf('Level = %.10g',handles.time(handles.sliceNumber+1)))
		setappdata(handles.handMir.figure1,'dem_x',handles.x);		% Don't get bad surprises (like load another file)
		setappdata(handles.handMir.figure1,'dem_y',handles.y);
	end

    guidata(handles.figure1,handles)

	% Save also the updated header in Mirone handles
	handles.handMir.head = handles.head;
    guidata(handles.handMir.figure1,handles.handMir)

% --------------------------------------------------------------------------
function set_common(handles, head)
% Common settingd to both 'init_header' functions
	set( handles.edit_x_min,'String',sprintf('%.8g',head(1)) )
	set( handles.edit_x_max,'String',sprintf('%.8g',head(2)) )
	set( handles.edit_y_min,'String',sprintf('%.8g',head(3)) )
	set( handles.edit_y_max,'String',sprintf('%.8g',head(4)) )
	set( handles.edit_x_inc,'String',sprintf('%.8g',head(8)) )
	set( handles.edit_y_inc,'String',sprintf('%.8g',head(9)) )

	set(handles.slider_layer,'Enable','on')
	set(handles.edit_sliceNumber,'Enable','on')
	set(handles.text_Info,'String',sprintf('Time steps = %d',handles.number_of_timesteps))

	set(handles.radio_multiLayer, 'Val', 1)
	set(handles.edit_multiLayerInc, 'Enable', 'on')
	set(handles.radio_timeGridsList,'Val',0)
	set([handles.textResize handles.popup_resize], 'Enable', 'off')
	set([handles.radio_stage handles.radio_xmoment handles.radio_ymoment handles.check_derivedVar], 'Enable', 'off')
	set([handles.edit_x_min handles.edit_x_max handles.edit_y_min handles.edit_y_max ...
		handles.edit_x_inc handles.edit_y_inc handles.edit_Ncols handles.edit_Nrows], 'Enable', 'inactive')

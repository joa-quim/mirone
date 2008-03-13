function aqua_suppfuns(opt, varargin)
% Supplement functions to allow using Aquamoto with plain netCDF coards grids

	switch opt
		case 'coards_hdr',		init_header_params(varargin{:})
		case 'coards_slice',	coards_sliceShow(varargin{:})
	end

% --------------------------------------------------------------------------
function init_header_params(handles,X,Y,head,misc)
% 
	handles.x = X;			handles.y = Y;
	handles.time = [];
	handles.number_of_timesteps = misc.z_dim(1);		% ... NEEDS THINKING
	
	set( handles.edit_x_min,'String',sprintf('%.8g',head(1)) )
	set( handles.edit_x_max,'String',sprintf('%.8g',head(2)) )
	set( handles.edit_y_min,'String',sprintf('%.8g',head(3)) )
	set( handles.edit_y_max,'String',sprintf('%.8g',head(4)) )
	handles.x_min = head(1);			handles.x_max = head(2);
	handles.y_min = head(3);			handles.y_max = head(4);
	
	set( handles.edit_x_inc,'String',sprintf('%.8g',head(8)) )
	set( handles.edit_y_inc,'String',sprintf('%.8g',head(9)) )
	set( handles.edit_Ncols,'String',sprintf('%d',misc.z_dim(end)) )
	set( handles.edit_Nrows,'String',sprintf('%d',misc.z_dim(end-1)) )

	% ------------- Finish slider configurations -------------
	s = handles.nc_info;
	if (handles.number_of_timesteps > 1)
		st = [1 10] / (handles.number_of_timesteps - 1);
		id = strmatch('time',{s.Dataset.Name});				% ONLY WHEN 3RTH DIM IS CALLED time
		if (~isempty(id))
			handles.time = double(nc_funs('varget', handles.fname, s.Dataset(id).Name));
		else
			handles.time = [1:handles.number_of_timesteps];
		end
		slMax = handles.number_of_timesteps;
	else
		slMax = 1+eps;	st = [1 1];		handles.time = 1;		% Defaults for no crashing
	end
	set(handles.slider_layer,'Min',1,'Max',slMax,'Val',1,'SliderStep',st) 	
	set(handles.slider_layer,'Enable','on')
	set(handles.edit_sliceNumber,'Enable','on')
	set(handles.text_Info,'String',sprintf('Time steps = %d',handles.number_of_timesteps))

	% ------ Compute individual and global min/maxs ----------------------------------
	handles.zMinMaxs = zeros(handles.number_of_timesteps,2);
	aguentabar(0,'title','Computing global min/max')
	for (k = 1:handles.number_of_timesteps)
		Z = nc_funs('varget', handles.fname, s.Dataset(misc.z_id).Name, [(k-1) 0 0], [1 s.Dataset(misc.z_id).Size(end-1:end)]);
		if ( isa(Z, 'double') )
			handles.zMinMaxs(k,:) = [min(Z(:)) max(Z(:))];
		else					% min/max are bugged when NaNs in singles
			zz = grdutils(Z,'-L');
			handles.zMinMaxs(k,:) = [zz(1) zz(2)];
		end
		aguentabar(k/handles.number_of_timesteps);
	end
	handles.zMinMaxsGlobal = [min(handles.zMinMaxs(:,1)) max(handles.zMinMaxs(:,2))];
	set(handles.edit_globalWaterMin,'String',handles.zMinMaxsGlobal(1))
	set(handles.edit_globalWaterMax,'String',handles.zMinMaxsGlobal(2))
	handles.minWater = handles.zMinMaxsGlobal(1);
	handles.maxWater = handles.zMinMaxsGlobal(2);
	head(5:6) = handles.zMinMaxs(1,:);				% Take the first slice min/max
	% ---------------------------------------------------------------------------------

	handles.cmapLand = jet(256);			% Reset the default colormap (default's Aquamoto is a specific one)

	handles.head = head;
	handles.illumComm = [];					% New file. Reset illum state.
	handles.imgBat = [];
	handles.netcdf_z_id = misc.z_id;
	handles.is_coards = true;
	set(handles.radio_multiLayer, 'Val', 1)
	set(handles.edit_multiLayerInc, 'Enable', 'on')
	set(handles.radio_timeGridsList,'Val',0)
	set([handles.textResize handles.popup_resize], 'Enable', 'off')
	set([handles.radio_stage handles.radio_xmoment handles.radio_ymoment handles.check_derivedVar], 'Enable', 'off')
	set([handles.edit_x_min handles.edit_x_max handles.edit_y_min handles.edit_y_max ...
		handles.edit_x_inc handles.edit_y_inc handles.edit_Ncols handles.edit_Nrows], 'Enable', 'inactive')
	set(handles.hTabAnuga,'String','netCDF')

	guidata(handles.figure1,handles)

% --------------------------------------------------------------------------------------------
function coards_sliceShow(handles)

	if ( isempty(handles.fname) )
		errordlg('Hey Lou. What about a walk on the Wild Side? Maybe you''ll find a little file there that you can use here!','Chico clever')
		return
	end

	z_id = handles.netcdf_z_id;
	s = handles.nc_info;			% Retrieve the .nc info struct 
	Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [handles.sliceNumber 0 0], [1 s.Dataset(z_id).Size(end-1:end)]);
	have_nans = grdutils(Z,'-N');	% No worry, very fast
	if ( have_nans && handles.useLandPhoto )
		alphaMask = alloc_mex(size(Z),'uint8');	% Create an image mask of Dry/Wets
		alphaMask(~isnan(Z)) = 255;				% nan pixeis will be transparent
	end

	% ----- Open or update a Mirone window with the slice display ----
	if (isempty(handles.hMirFig) || ~ishandle(handles.hMirFig))			% First run or killed Mirone window
		tmp.X = handles.x;		tmp.Y = handles.y;		tmp.head = handles.head;	tmp.cmap = handles.cmapLand;
		tmp.name = sprintf('Layer = %g',handles.time(handles.sliceNumber+1));
		handles.hMirFig = mirone(Z, tmp);
		move2side(handles.figure1,handles.hMirFig,'left')
		handles.handMir = guidata(handles.hMirFig);			% Get the handles of the now existing Mirone fig
		handles.firstLandPhoto = true;
		if ( handles.useLandPhoto )
			h = image('XData',handles.geoPhotoX,'YData',handles.geoPhotoY, 'CData',handles.geoPhoto, 'Parent',handles.handMir.axes1);
			uistack(h,'bottom')
			handles.firstLandPhoto = false;
			set(handles.handMir.hImg,'AlphaData',alphaMask)	% 'alphaMask' was updated ... maybe somewhere
		end

	else									% We already have a Mirone image. Update it with this new slice
		handles.handMir = guidata(handles.hMirFig);			% Get updated handles to see if illum has changed
		setappdata(handles.handMir.figure1,'dem_z',Z);		% Update grid so that coursor display correct values
															% Have to do it here because minmax arg to scalet8 CHANGES Z

		if ( handles.useLandPhoto )							% External Land image
			if (handles.firstLandPhoto)						% First time, create the background image
				h = image('XData',handles.geoPhotoX,'YData',handles.geoPhotoY, 'CData',handles.geoPhoto, 'Parent',handles.handMir.axes1);
				uistack(h,'bottom')
				handles.firstLandPhoto = false;
			end
			set(handles.handMir.hImg,'AlphaData',alphaMask)	% 'alphaMask' was updated ... somewhere
		end

		if ( ~get(handles.check_globalMinMax, 'Val') ),		minmax = [];		% Use Slice's min/max
		else							minmax = handles.zMinMaxsGlobal;
		end

		if ( ~isempty(minmax) ),		img = scaleto8(Z, 8, minmax);
		else							img = scaleto8(Z);
		end

		if ( get(handles.radio_shade, 'Val') )
			indVar = 1;								% FAR FROM SURE THAT THIS IS CORRECT
			img = ind2rgb8(img, handles.cmapLand);		% img is now RGB
			head = handles.head;
			if ( ~isempty(handles.ranges{indVar}) ),	head(5:6) = handles.ranges{indVar};		end
			R = illumByType(handles, Z, head, handles.landIllumComm);
			img = shading_mat(img,R,'no_scale');		% and now it is illuminated
		end

		set(handles.handMir.hImg, 'CData', img)
		set(handles.handMir.figure1, 'Name', sprintf('Level = %g',handles.time(handles.sliceNumber+1)))
		setappdata(handles.handMir.figure1,'dem_x',handles.x);		% Don't get bad surprises (like loaded another file)
		setappdata(handles.handMir.figure1,'dem_y',handles.y);
	end
	
    guidata(handles.figure1,handles)

	% Save also the updated header in Mirone handles
	handles.handMir.head = handles.head;
    guidata(handles.handMir.figure1,handles.handMir)
	

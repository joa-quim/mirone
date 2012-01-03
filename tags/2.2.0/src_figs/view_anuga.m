function varargout = view_anuga(varargin)
% M-File changed by desGUIDE

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

	hObject = figure('Tag','figure1','Visible','off');
	view_anuga_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	if (numel(varargin) > 0)
		handMir = varargin{1};
		handles.home_dir = handMir.home_dir;
		handles.last_dir = handMir.last_dir;
		handles.work_dir = handMir.work_dir;
        d_path = handMir.path_data;
	else
		handles.home_dir = cd;
		handles.last_dir = handles.home_dir;
		handles.work_dir = handles.home_dir;
        d_path = [pwd filesep 'data' filesep];
	end

	% Import icon
	load([d_path 'mirone_icons.mat'],'Mfopen_ico','Marrow_ico');
	set(handles.push_swwName,'CData',Mfopen_ico)
	set(handles.toggle_vel,'CData',Marrow_ico)
	clear Mfopen_ico Marrow_ico;

	handles.handMir = [];
	handles.fname = [];
	handles.sliceNumber = 0;
	handles.one_or_zero = 1;
	handles.first = true;
	handles.volumes = [];		handles.illumComm = [];
	handles.dms_xinc = 0;		handles.dms_yinc = 0;
	handles.hQuiver = [];		handles.hMirFig = [];
	handles.imgBat = [];		% Bathymetry only (that is, dry) image
	handles.cmapBat = [];		% To hold the color map with a discontinuity at the shoreline
	handles.geog = 0;			% For the time being ANUGA no geoga

	set(handles.popup_derivedVar, 'String', ...
		{'Absolute Velocity (V)'; ...
		'Absolute Momentum (VxD)'; ...
		'Water Depth'; ...
		'Elevation'; ...
		'Max Water'; ...
		'Froude Number'; ...
		'Velocity Head (V^2 / (2g))'; ...
		'Hazard-RVD (D(1+V^2))'; ...
		'Taylor-V (D*(1+V+V**2))' })

	% In push_swwName_CB() individual ranges MUST be assigned to handles.ranges
	% following exactly the variables order used in the popup_derivedVar.
	handles.ranges = cell(12,1);		% 12 = 3 (direct vars) + 9 (derived vars)
	handles.elevRange = [];

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, handles.text_GLG, handles.frame1)
	%------------- END Pro look (3D) -----------------------------------------------------

	set(hObject,'Visible','on');

	% The rest is donne in push_swwName_CB() because only than we have all the necessary info
	S = load([d_path 'gmt_other_palettes.mat'],'DEM_screen');
	handles.cmapLand = S.DEM_screen;
% 	S = load([f_path 'gmt_other_palettes.mat'],'Terre_Mer');
% 	handles.terraMar = S.Terre_Mer;

	% By default use a blue only colormap for water
	handles.cmapWater = [0 0 1; 0 0 1];
% 	S = load([d_path 'gmt_other_palettes.mat'],'polar');
% 	handles.cmapWater = S.polar;

	% This will cause a silent error but it also load the mex file in memory so it will be fast on "first" use
	try     lili = nc_funs('info','lixoxo');	end

	guidata(hObject, handles);
	if (nargout),   varargout{1} = hObject;     end

% --------------------------------------------------------------------
function radio_stage_CB(hObject, handles)
	if (get(hObject,'Value'))
		set([handles.radio_xmoment handles.radio_ymoment], 'Value', 0)
	else
		set(hObject,'Value', 1)
	end

% --------------------------------------------------------------------
function radio_xmoment_CB(hObject, handles)
	if (get(hObject,'Value'))
		set([handles.radio_stage handles.radio_ymoment], 'Value', 0)
	else
		set(hObject,'Value', 1)
	end

% --------------------------------------------------------------------
function radio_ymoment_CB(hObject, handles)
	if (get(hObject,'Value'))
		set([handles.radio_stage handles.radio_xmoment], 'Value', 0)
	else
		set(hObject,'Value', 1)
	end

% --------------------------------------------------------------------
function check_derivedVar_CB(hObject, handles)
	if (get(hObject,'Value'))
		set([handles.radio_stage handles.radio_xmoment handles.radio_ymoment], 'Enable', 'off')
		set(handles.popup_derivedVar, 'Enable', 'on')
	else
		set([handles.radio_stage handles.radio_xmoment handles.radio_ymoment], 'Enable', 'on')
		set(handles.popup_derivedVar, 'Enable', 'off')
	end

% --------------------------------------------------------------------
function popup_derivedVar_CB(hObject, handles)
	contents = get(handles.popup_derivedVar, 'String');
	qual = contents{get(handles.popup_derivedVar,'Value')};
	switch qual(1:min(numel(qual),11))
		case {'Absolute Ve' 'Absolute Mo'}			% Absolute Velocity || Momentum
			set(handles.toggle_vel,'Vis', 'on')
		otherwise
			set(handles.toggle_vel,'Vis', 'off', 'Value', 0)
	end

% -----------------------------------------------------------------------------------------
function toggle_vel_CB(hObject, handles)
	if ( ~get(hObject, 'Val') && ~isempty(handles.hQuiver) )
		try     delete(handles.hQuiver),	end
		handles.hQuiver = [];
		guidata(handles.figure1, handles);
	end

% -----------------------------------------------------------------------------------------
function checkbox_splitDryWet_CB(hObject, handles)
	if ( get(hObject, 'Val') )
		set(handles.checkbox_globalMinMax, 'Val', 0)		% Spliting impliyes RGB final image, so no use scaling
	end

% -----------------------------------------------------------------------------------------
function slider1_CB(hObject, handles)
	handles.sliceNumber = round(get(handles.slider1,'Value')) - 1;
	set(handles.edit_sliceNumber,'String', handles.sliceNumber+1)		% Update slice nº box
	set(handles.figure1,'pointer','watch')
	push_showSlice_CB([], [], handles)		% and update image (also saves handles)
	set(handles.figure1,'pointer','arrow')

% -----------------------------------------------------------------------------------------
function edit_sliceNumber_CB(hObject, handles)
	xx = fix(str2double(get(hObject,'String')));		% Make sure its an int
	if (isnan(xx) || xx < 1 || xx > handles.number_of_timesteps)
		handles.sliceNumber = 0;		set(hObject,'String','1')
	else
		set(hObject,'String',xx);		set(handles.slider1,'Val',xx)		% Update slider
		handles.sliceNumber = xx - 1;
	end
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_swwName_CB(hObject, handles)
    fname = get(hObject,'String');
    push_swwName_CB([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_swwName_CB(hObject, handles, opt)
% This function does quite some work. It reads and extract relevant info from the netCDF file

    if (nargin == 3)        % Direct call
        cd(handles.last_dir)
    	str1 = {'*.sww;*.SWW;', 'Data files (*.sww,*.SWW)';'*.*', 'All Files (*.*)'};
        [FileName,PathName] = uigetfile(str1,'sww file');
        cd(handles.home_dir);
	    if isequal(FileName,0),		return,		end
        if (PathName ~= 0),         handles.last_dir = PathName;    end
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
	pause(0.01);	handles.fname = [PathName FileName];
	
	if (exist(handles.fname, 'file') ~= 2)
		errordlg(['File: ' handles.fname ' does not exist.'],'Error')
		handles.fname = [];
		return
	end
    set(handles.edit_swwName,'String',handles.fname)

	% ----------------- this make sense when one are reloading a(nother) file ----------
	handles.sliceNumber = 0;
    set(handles.slider1,'Value',1)
    set(handles.edit_sliceNumber,'String','1')
	
	% ---- Maybe the dimensions should be fished out of the "s" structurem as well -----
	% But for now I'll just test that 'number_of_volumes' exists, otherwise ... street
	set(handles.figure1,'pointer','watch')
	s = nc_funs('info',handles.fname);
	attribNames = {s.Attribute.Name};
	ind = strcmp({s.Dimension.Name},'number_of_volumes');
	if (~any(ind))
		errordlg('ERROR: This .sww file is not of recognizable type. For example: "number_of_volumes" was not found.','Error')
		set(handles.figure1,'pointer','arrow')
		return
	end

	ind = strcmp(attribNames,'xllcorner');		xllcorner = 0;
	if (any(ind)),	xllcorner = s.Attribute(ind).Value;		end
	ind = strcmp(attribNames,'yllcorner');		yllcorner = 0;
	if (any(ind)),	yllcorner = s.Attribute(ind).Value;		end
	
	% Fill the grid size boxes which imply some parameter guessings
	st = nc_funs('getdiminfo', handles.fname,'number_of_volumes');
	handles.number_of_volumes = st.Length;
	st = nc_funs('getdiminfo', handles.fname,'number_of_points');
	handles.number_of_points = st.Length;
	st = nc_funs('getdiminfo', handles.fname,'number_of_timesteps');
	handles.number_of_timesteps = st.Length;

	% ------------------ OK, Get numerics now -----------------------------------
	handles.x = double(nc_funs('varget', handles.fname, 'x')) + xllcorner;
	handles.y = double(nc_funs('varget', handles.fname, 'y')) + yllcorner;
	handles.time = nc_funs('varget', handles.fname, 'time');
	handles.volumes = nc_funs('varget', handles.fname, 'volumes');
	if (~isa(handles.volumes, 'int32')),	handles.volumes = int32(handles.volumes);	end
	set(handles.figure1,'pointer','arrow')

	head = [min(handles.x) max(handles.x) min(handles.y) max(handles.y) 0 1 0];
	set( handles.edit_x_min,'String',sprintf('%.8g',head(1)) )
	set( handles.edit_x_max,'String',sprintf('%.8g',head(2)) )
	set( handles.edit_y_min,'String',sprintf('%.8g',head(3)) )
	set( handles.edit_y_max,'String',sprintf('%.8g',head(4)) )
	handles.x_min = head(1);			handles.x_max = head(2);
	handles.y_min = head(3);			handles.y_max = head(4);
	handles.x_min_or = head(1);			handles.x_max_or = head(2);
	handles.y_min_or = head(3);			handles.y_max_or = head(4);

	% ------------------- Get/compute the global min/max --------------------------------
	varNames = {s.Dataset.Name};						% WHAT IF THIS FAILS??? 
	ind = strcmp(varNames,'stage_range');
	if (any(ind))
		handles.ranges{1} = double(nc_funs('varget', handles.fname, 'stage_range'));
	end
	ind = strcmp(varNames,'xmomentum_range');
	if (any(ind))
		handles.ranges{2} = double(nc_funs('varget', handles.fname, 'xmomentum_range'));
	end
	ind = strcmp(varNames,'ymomentum_range');
	if (any(ind))
		handles.ranges{3} = double(nc_funs('varget', handles.fname, 'ymomentum_range'));
	end
	ind = strcmp(varNames,'elevation_range');
	if (any(ind))				% Not unlikely that some of the following are wrong. F... devided by (nearly) zero
		handles.elevRange = double(nc_funs('varget', handles.fname, 'elevation_range'));
		% ok, here we are going to compute the derived vars min/max as well
		D = [eps; max(handles.ranges{1} - handles.elevRange)];		% Why the hell handles.ranges{1}(2) == handles.elevRange(2) !!!!
		handles.ranges{5}  = [0; max(sqrt(handles.ranges{2} .^ 2 + handles.ranges{3} .^ 2))];	% momentumRange
		handles.ranges{6}  = D;											% waterDepthRange
		handles.ranges{4}  = handles.ranges{5} ./ D;					% vRange
		handles.ranges{9}  = handles.ranges{4} ./ sqrt(D * 9.8);		% V / sqrt(gD) - froudeRange
		handles.ranges{10} = handles.ranges{4} .^ 2 / (2*9.8);			% vHeadRange
		handles.ranges{11} = D .* (1 + handles.ranges{4} .^ 2);			% D(1+V^2) - hazard_rvd_Range
		handles.ranges{12} = D .* (1 + handles.ranges{4} + handles.ranges{4} .^ 2);		% D(1+V+V^2) - taylorVRange
	end

	% --------------- Estimate a "reasonable" proposition for grid size ---------------------
	n = round( sqrt(double(handles.number_of_volumes)) );
	inc = ( diff(head(1:2)) + diff(head(3:4)) ) / (2*(n-1));	% A mean dx dy
	set( handles.edit_x_inc,'String',sprintf('%.8g',inc) )
	set( handles.edit_y_inc,'String',sprintf('%.8g',inc) )
	set( handles.edit_Ncols,'String',sprintf('%d',n) )
	set( handles.edit_Nrows,'String',sprintf('%d',n) )
	% Call dim_funs to compute & update the correct size
	dim_funs('xInc', handles.edit_x_inc, handles)
	dim_funs('yInc', handles.edit_y_inc, handles)

	% ----------------- Remainings ... --------------------------------------------------------
	handles.head = head;		% INCOMPLETE HEAD. The rest is computed in push_showSlice_CB()

	if ( ~isempty(handles.elevRange) )
		head(5:6) = handles.elevRange;
		handles.cmapBat = makeCmapBat(handles, head, handles.cmapLand, 1);		% Put the cmap discontinuity at the zero of bat
	else
		warndlg('Could not find elevation (bathymetry) min/max. End of the world is NEEEAAAAR','WARNING')	% SCREEEEEMMMMMMMMM
	end

	% ----------------- Finish slider configurations ------------------------------------------
	st = [1 10] / (handles.number_of_timesteps - 1);
	set(handles.slider1,'Min',1,'Max',handles.number_of_timesteps,'Val',1,'SliderStep',st) 	
	set(handles.slider1,'Enable','on')
	
	set(handles.edit_sliceNumber,'Enable','on')
	set(handles.text_Info,'String',sprintf('Triangles = %d & Time steps = %d',handles.number_of_volumes,handles.number_of_timesteps))

	handles.illumComm = [];					% New file. Reset illum state.
	handles.imgBat = [];
 
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function push_showSlice_CB(hObject, handles)
	if (isempty(handles.fname))
		errordlg('Hey Lou. What about a walk on the Wild Side? Maybe you''ll find a file there that you can use here!','Chico clever')
		return
	end

	nx = str2double(get(handles.edit_Ncols,'String'));
	ny = str2double(get(handles.edit_Nrows,'String'));

	if ( ~get(handles.check_derivedVar,'Val') )			% Get one of the primary quantities
		% Get the ploting variable
		if (get(handles.radio_stage, 'Val')),			theVarName = 'stage';		indVar = 1;
		elseif (get(handles.radio_xmoment, 'Val')),		theVarName = 'xmomentum';	indVar = 2;
		elseif (get(handles.radio_ymoment, 'Val')),		theVarName = 'ymomentum';	indVar = 3;
		end
		theVar = nc_funs('varget', handles.fname, theVarName, [handles.sliceNumber 0], [1 handles.number_of_points]);
		U = [];		V = [];
		indWater = [];					% Not empty when a Max quantity (obvioulsly, not the case here)
	else
		[theVar, U, V, indVar, indWater] = get_derivedVar(handles);
	end
	if (~isa(theVar, 'double')),	theVar = double(theVar);	end		% While we don't f... these doubles as well

	% create a grid in x and y
	x = linspace(handles.head(1),handles.head(2),nx);
	y = linspace(handles.head(3),handles.head(4),ny);
	Z = mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) theVar(:)],x,y);
	handles.head(5:6) = [min(Z(:)) max(Z(:))];

	splitDryWet = get(handles.checkbox_splitDryWet, 'Val');		% See if we need to build a wet and dry images, or only one

	if ( splitDryWet )
		indLand = get_landInd(handles, x, y, indWater);			% Compute Dry/Wet indexes - have to recompute, ... since water ... moves
		if ( isempty(handles.imgBat) )							% First time, compute it (not shaded)
			handles.imgBat = do_imgBat(handles, indVar, x, y);	% IMG is always RGB
		end
	end

	handles.head(8) = str2double(get(handles.edit_x_inc,'String'));
	handles.head(9) = str2double(get(handles.edit_y_inc,'String'));

	if (isempty(handles.hMirFig) || ~ishandle(handles.hMirFig))			% First run or killed Mirone window
		tmp.X = x;		tmp.Y = y;		tmp.head = handles.head;
		tmp.name = sprintf('SWW time = %g',handles.time(handles.sliceNumber+1));
		handles.hMirFig = mirone(Z,tmp);
		handles.handMir = guidata(handles.hMirFig);			% Get the handles of the now existing Mirone fig

	else									% We already have a Mirone image. Update it with this new slice
		handles.handMir = guidata(handles.hMirFig);			% Get updated handles to see if illum has changed
		illumComm = getappdata(handles.handMir.figure1,'illumComm');
		setappdata(handles.handMir.figure1,'dem_z',Z);		% Update grid so that cursor display correct values
															% Have to do it here because minmax arg to scalet8 CHANGES Z
		if ( splitDryWet && handles.handMir.Illumin_type )	% Land/Water spliting with Illumination
			if ( ~isequal(handles.illumComm, illumComm) )	% Recomputed 'imgBat' if illumination has changed
				handles.imgBat = do_imgBat(handles, indVar, x, y);
				handles.illumComm = illumComm;				% save illum command for future comparison
			end
		elseif ( ~splitDryWet )								% No Land/Water spliting
			minmax = handles.ranges{indVar};		minmax = minmax(:)';
			if ( ~get(handles.checkbox_globalMinMax, 'Val') ),		minmax = [];	end		% Use Slice's min/max
			if ( ~isempty(minmax) ),		img = scaleto8(Z, 8, minmax);
			else							img = scaleto8(Z);
			end

			if ( handles.handMir.Illumin_type >= 1 && handles.handMir.Illumin_type <= 4 )%&& ...	%  with Illumination
					%~isequal(handles.illumComm, illumComm) )
				img = ind2rgb8(img,get(handles.handMir.figure1,'Colormap'));	% img is now RGB	
				head = handles.head;
				if ( ~isempty(handles.ranges{indVar}) ),	head(5:6) = handles.ranges{indVar};		end
				R = illumByType(handles, Z, head, illumComm);
				img = shading_mat(img,R,'no_scale');		% and now it is illuminated
				handles.illumComm = illumComm;				% save illum command for future comparison
			end
			set(handles.handMir.hImg, 'CData', img)
		end

		setappdata(handles.handMir.figure1,'dem_x',x);		% Don't get bad surprises if space increments have changed
		setappdata(handles.handMir.figure1,'dem_y',y);
		set(handles.handMir.figure1, 'Name', sprintf('SWW time = %g',handles.time(handles.sliceNumber+1)))
	end
	
	if ( splitDryWet )
		img = do_imgWater(handles, indVar, Z, handles.imgBat, indLand);		% IMG is always RGB
		set(handles.handMir.hImg, 'CData', img)
	end

	if ( ~isempty(U) )		% Plot vectors
		if ( ~isempty(handles.hQuiver) ),	try     delete(handles.hQuiver),	end,	end
		x = linspace(handles.head(1),handles.head(2), min(fix(nx/10), 15) );
		y = linspace(handles.head(3),handles.head(4), min(fix(ny/10), 15) );
		U = double( mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) U(:)],x,y) );
		V = double( mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) V(:)],x,y) );
		handles.hQuiver = loc_quiver(handles.handMir.axes1, x, y, U, V);
	end

    guidata(handles.figure1,handles)

	% Save also the updated header in Mirone handles
	handles.handMir.head = handles.head;
    guidata(handles.handMir.figure1,handles.handMir)

% --------------------------------------------------------------------
function [theVar, U, V, indVar, indWater] = get_derivedVar(handles)
	% Compute a derived quantity from a combination of the primary quantities
	% INDVAR is the index of the variable in the popup. It starts at 4 because
	% the first 3 are 'stage', 'xmomentum' and 'ymomentum'.
	% INDWATER ~= [] when THEVAR is maximum ... (e.g. Max Water)
	%
	% 4  'Absolute Velocity (V)'; ...
	% 5  'Absolute Momentum (VxD)'; ...
	% 6  'Water Depth'; ...
	% 7  'Elevation'; ...
	% 8  'Max Water'; ...
	% 9  'Froude Number'; ...
	% 10 'Velocity Head (V^2 / (2g))'; ...
	% 11 'Hazard-RVD (D(1+V^2))'; ...
	% 12 'Taylor-V (D*(1+V+V**2))'

	U = [];		V = [];		indWater = [];
	contents = get(handles.popup_derivedVar, 'String');
	qual = contents{get(handles.popup_derivedVar,'Value')};
	switch qual(1:min(numel(qual),11))
		case {'Absolute Ve' 'Froude Numb'}			% Absolute Velocity (V) || Froude Number
			x = nc_funs('varget', handles.fname, 'xmomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			y = nc_funs('varget', handles.fname, 'ymomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			if (~isa(x, 'double')),		x = double(x);		y = double(y);		end
			if ( qual(1) == 'A' && get(handles.toggle_vel, 'Value') )
				U = x;		V = y;					% momentums copy
			end
			theVar = sqrt(x.^2 + y.^2);				% |M|
			x = nc_funs('varget', handles.fname, 'stage', [handles.sliceNumber 0], [1 handles.number_of_points]);
			y = nc_funs('varget', handles.fname, 'elevation')';
			if (~isa(x, 'double')),		x = double(x);		y = double(y);		end
			D = x - y + 1e-10;		clear x y;
			ind_0 = (D < 1e-8);			% To get arround a devide-by-nearly-zero and Anuga bug in velocity problem 
			theVar = theVar ./ D;
			theVar(ind_0) = 0; 
			if ( ~isempty(U) )				% Now U,V are vx,vy
				U = U ./ D;		V = V ./ D;
				U(ind_0) = 0;	V(ind_0) = 0;
			end
			indVar = 4;
			if ( strcmp(qual(1:6), 'Froude') )		% V / sqrt(gD)
				theVar = theVar ./ sqrt( 9.8 * D );
				indVar = 9;
			end
		case 'Absolute Mo'			% Absolute Momentum (VxD)
			x = nc_funs('varget', handles.fname, 'xmomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			y = nc_funs('varget', handles.fname, 'ymomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			if (~isa(x, 'double')),		x = double(x);		y = double(y);		end
			theVar = sqrt(x.^2 + y.^2);
			if ( qual(1) == 'A' && get(handles.toggle_vel, 'Value') )
				U = x;		V = y;
			end
			indVar = 5;
		case 'Elevation'			% Elevation
			theVar = nc_funs('varget', handles.fname, 'elevation')';
			indVar = 7;
		case 'Water Depth'			% Water Depth
			stage = nc_funs('varget', handles.fname, 'stage', [handles.sliceNumber 0], [1 handles.number_of_points]);
			elevation = nc_funs('varget', handles.fname, 'elevation')';
			if (~isa(stage, 'double')),		stage = double(stage);		elevation = double(elevation);		end
			theVar = (stage - elevation);
			indVar = 6;
		case 'Max Water'			% Maximum water height (we need to compute it)
			h = waitbar(0, 'Computing max water height ...');
			theVar = zeros(1, handles.number_of_points);
			stage = nc_funs('varget', handles.fname, 'elevation')';
			indWater = (stage < 0);	% Still water indices 
			for (k = 0:handles.number_of_timesteps - 1)
				stage = nc_funs('varget', handles.fname, 'stage', [k 0], [1 handles.number_of_points]);
				ind = (stage > theVar);
				theVar(ind) = stage(ind);
				% Compute maximum inundation (runin). Do not use k = 0 because than all IND == 1, since theVar was still == 0
				if (k > 0),		indWater(ind) = 1;		end
				waitbar((k+1) / handles.number_of_timesteps)
			end
			if (ishandle(h)),	close(h);	end
			indVar = 8;
% 		case 'Max Depth'			% Maximum water Depth (we need to compute it)
% 			h = waitbar(0, 'Computing max depth ...');
% 			theVar = zeros(1, handles.number_of_points);
% 			elevation = nc_funs('varget', handles.fname, 'elevation')';
% 			for (k = 0:handles.number_of_timesteps - 1)
% 				stage = nc_funs('varget', handles.fname, 'stage', [handles.sliceNumber 0], [1 handles.number_of_points]);
% 				theVar = cvlib_mex('sub',stage,elevation);		% Water depth
% 				ind = (theVar == 0);
%				if ( isequal(ind) )
% 				waitbar((k+1) / handles.number_of_timesteps)
% 			end
% 			if (ishandle(h)),	close(h);	end
% 
% 			stage = nc_funs('varget', handles.fname, 'stage', [handles.sliceNumber 0], [1 handles.number_of_points]);
% 			theVar = cvlib_mex('sub',stage,elevation);
% 			indVar = 8;
			
		case {'Velocity He' 'Hazard-RVD ' 'Taylor-V (D'}	% Velocity Head (V^2 / (2g)) || D * (1 + V^2) || D*(1+V+V^2)
			x = nc_funs('varget', handles.fname, 'xmomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			y = nc_funs('varget', handles.fname, 'ymomentum', [handles.sliceNumber 0], [1 handles.number_of_points]);
			if (~isa(x, 'double')),		x = double(x);		y = double(y);		end
			theVar = (x.^2 + y.^2);		% = D^2 * V^2 = (DV)^2
			x = nc_funs('varget', handles.fname, 'stage', [handles.sliceNumber 0], [1 handles.number_of_points]);
			y = nc_funs('varget', handles.fname, 'elevation')';
			if (~isa(x, 'double')),		x = double(x);		y = double(y);		end
			D = (x - y + 1e-10);
			clear x y;
			ind_0 = (D < 1e-8);			% To get arround a devide-by-nearly-zero and Anuga bug in velocity problem 
			if (qual(1) == 'V')
				D = D .* D;
				theVar = theVar ./ (2 * 9.8 * D);
				indVar = 10;
			elseif (qual(1) == 'H')
				theVar = D + theVar ./ D;		% = D + (DV)^2/D = D + DV^2 = D * (1 + V^2)
				indVar = 11;
			else
				% theVar = D .* (1 + sqrt(theVar ./ (D .* D)) + theVar ./ (D .* D) );	% D(1+V+V^2) 
				% theVar = D + theVar ./ D + sqrt((D .* D) .* theVar ./ (D .* D));		% D + DV^2 + sqrt(D^2 * (DV)^2 / D^2)
				theVar = D + theVar ./ D + sqrt(theVar);								% D + DV^2 + DV = D(1+V+V^2)
				indVar = 12;
			end
			theVar(ind_0) = 0; 
	end

% -----------------------------------------------------------------------------------------
function img = do_imgWater(handles, indVar, Z, imgBat, indLand)
% Compute Water (wet) image and illuminate it if that's the case
%	INDVAR -> indice of the currently processing variable (either a direct or derived var)
%	X,Y regular coordinates to interpolate the triangles into regular grid
%
%	Output IMG is always RGB

	minmax = handles.ranges{indVar};		minmax = minmax(:)';
% 	Z(indLand) = 0;
	if ( ~get(handles.checkbox_globalMinMax, 'Val') ),		minmax = [];	end		% Use Slice's min/max
	if ( ~isempty(minmax) ),		imgWater = scaleto8(Z, 8, minmax);
	else							imgWater = scaleto8(Z);
	end

	handles.handMir = guidata(handles.hMirFig);			% Get updated handles to see if illum has changed
	if ( handles.handMir.Illumin_type >= 1 && handles.handMir.Illumin_type <= 4 )
		illumComm = getappdata(handles.handMir.figure1,'illumComm');
		pal = get(handles.handMir.figure1,'Colormap');
		if ( get(handles.checkbox_splitDryWet, 'Val') ),	pal = handles.cmapWater;	end
		imgWater = ind2rgb8(imgWater, pal);						% image is now RGB
		R = illumByType(handles, Z, handles.head, illumComm);
		imgWater = shading_mat(imgWater,R,'no_scale');			% and now it is illuminated
		clear R;
	end
	if ( ndims(imgWater) == 2 ),	imgWater = ind2rgb8(imgWater,handles.cmapWater);	end		% Like promissed above
	alfa = 0.0;
	img = mixe_images(handles, imgBat, imgWater, indLand, alfa);

% ----------------------------------------------------------------------------------------------	
function imgBat = do_imgBat(handles, indVar, x, y)	
% Compute Land (dry) image and illuminate it if that's the case
%	INDVAR -> indice of the currently processing variable (either a direct or derived var)
%	X,Y regular coordinates to interpolate the triangles into regular grid
%
%	Output IMGBAT is always RGB

	bat = nc_funs('varget', handles.fname, 'elevation')';
	if (~isa(bat,'double') ),		bat = double(bat);	end
	bat = mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) bat(:)], x, y);

	imgBat = scaleto8(bat);
    imgBat = ind2rgb8(imgBat, handles.cmapBat);
	if ( ~isempty(handles.handMir) && handles.handMir.Illumin_type >= 1 && handles.handMir.Illumin_type <= 4 )
		illumComm = getappdata(handles.handMir.figure1,'illumComm');
		head = handles.head;
		if ( ~isempty(handles.ranges{indVar}) ),	head(5:6) = handles.ranges{indVar};		end
		R = illumByType(handles, bat, head, illumComm);
		imgBat = shading_mat(imgBat, R, 'no_scale');
	end

% -----------------------------------------------------------------------------------------
function R = illumByType(handles, Z, head, illumComm)
% Compute the illuminance matrix in function of the illumination type

	if (handles.handMir.Illumin_type == 1)
		if (handles.geog),  R = grdgradient_m(Z,head,'-M',illumComm,'-Nt');
		else                R = grdgradient_m(Z,head,illumComm,'-Nt');
		end
	else
		R = grdgradient_m(Z,head,illumComm);
	end

% -----------------------------------------------------------------------------------------
function ind = get_landInd(handles, x, y, indWater)
    % Compute indices such that 1 -> Dry; 0 -> Wet
	% Note that "Land" (stage - elevation) is gridded
	% INDWATER, if not empty, has ones at the indices of a max quantity areal extent

	if ( isempty(indWater) )		% Compute Dry/Wet indices for a specific time slice
		stage = nc_funs('varget', handles.fname, 'stage', [handles.sliceNumber 0], [1 handles.number_of_points]);
		elevation = nc_funs('varget', handles.fname, 'elevation')';
		dife = cvlib_mex('absDiff',stage,elevation);
		if (~isa(dife,'double') ),		dife = double(dife);	end
		dife = mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) dife(:)],x,y);
		ind = (dife < 1e-7);
	else							% A Max quantity. INDWATER has ones at the indices of that quantity areal extent
		indWater = mxgridtrimesh(handles.volumes, [handles.x(:) handles.y(:) double(indWater(:))],x,y);
		ind = (indWater < 0.5);		% Test with < 0.5 since we dont't know if mxgridtrimesh has changed the ones
	end

% -----------------------------------------------------------------------------------------
function imgWater = mixe_images(handles, imgBat, imgWater, ind, alfa)
% Mixes land and water images simulating transparency.
% It also resizes the image if the scale factor is ~= 1
    
    try
        if (alfa > 0.01)    % Only if transparency is greater than 1%
            cvlib_mex('addweighted',imgWater,(1 - alfa),imgBat,alfa);     % In-place
        end
        
        tmpW = imgWater(:,:,1);     tmpB = imgBat(:,:,1);           % R
        tmpW(ind) = tmpB(ind);      imgWater(:,:,1) = tmpW;
        
        tmpW = imgWater(:,:,2);     tmpB = imgBat(:,:,2);           % G
        tmpW(ind) = tmpB(ind);      imgWater(:,:,2) = tmpW;
        
        tmpW = imgWater(:,:,3);     tmpB = imgBat(:,:,3);           % B
        tmpW(ind) = tmpB(ind);      imgWater(:,:,3) = tmpW;

    catch
        errordlg(['View Anuga:mixe_images ' lasterr],'Error')
    end

% -----------------------------------------------------------------------------------------
function new_cmap = makeCmapBat(handles, head, cmap, orig)
% Put the cmap discontinuity at the zero of bathymetry (coastline)

	if (~orig),		new_cmap = cmap;	return,		end		% Untill I know better what to do

	% Isto assume que a bat tem partes neg e pos (testar)
	% cmap = handles.terraMar;
	% ind_old = 147;			% Discontinuity in the handles.terraMar cmap
	cmap(end-4:end,:) = [];		% Remove last five colors (nearly white)
	cmap = [repmat([196 156 104]/255,5,1); cmap];		% Add a yelowish color
	ind_old = 5;				% New discontinuity in cmap

	nc = length(cmap);
	z_inc = (head(6) - head(5)) / (nc - 1);
	ind_c = round(abs(0 - head(5)) / z_inc + 1);

	nl = ind_old;		nu = ind_c;
	new_cmap_l = interp1(linspace(0,1,nl), cmap(1:nl,:), linspace(0,1,nu));
	new_cmap_u = interp1(linspace(0,1,nc-ind_old), cmap(ind_old+1:nc,:), linspace(0,1,nc-ind_c));
	new_cmap = [new_cmap_l; new_cmap_u];

% --------------------------------------------------------------------
function push_showMesh_CB(hObject, handles)
	
	if (isempty(handles.fname)),	errordlg('Go to ... #!?-#$%}&*"*_yu','#!?-#$%}'),	return,		end

	% Get the 'stage' vector
	stage = nc_funs('varget', handles.fname, 'stage', [handles.sliceNumber 0], [1 handles.number_of_points]);
	h = figure;
	trisurf(double(handles.volumes)+1,double(stage(:)),double(handles.x(:)),double(handles.y(:)),'facecolor',[.9 .8 .6],'edgecolor','k');
	axis image, view(90,0), axis vis3d, axis off, zoom_j(1.5)

% -------------------------------------------------------------------------------------
function edit_x_min_CB(hObject, handles)
	dim_funs('xMin', hObject, handles)

% -------------------------------------------------------------------------------------
function edit_x_max_CB(hObject, handles)
	dim_funs('xMax', hObject, handles)

% --------------------------------------------------------------------
function edit_y_min_CB(hObject, handles)
	dim_funs('yMin', hObject, handles)

% --------------------------------------------------------------------
function edit_y_max_CB(hObject, handles)
	dim_funs('yMax', hObject, handles)

% --------------------------------------------------------------------
function edit_x_inc_CB(hObject, handles)
	dim_funs('xInc', hObject, handles)

% --------------------------------------------------------------------
function edit_Ncols_CB(hObject, handles)
	dim_funs('nCols', hObject, handles)

% --------------------------------------------------------------------
function edit_y_inc_CB(hObject, handles)
	dim_funs('yInc', hObject, handles)

% --------------------------------------------------------------------
function edit_Nrows_CB(hObject, handles)
	dim_funs('nRows', hObject, handles)

% --------------------------------------------------------------------
function push_Help_CB(hObject, handles)


% --------------------------------------------------------------------
function push_Cancel_CB(hObject, handles)
	delete(handles.figure1)
	
% --------------------------------------------------------------------
function hh = loc_quiver(hAx,varargin)
%QUIVER Quiver plot.
%   QUIVER(X,Y,U,V) plots velocity vectors as arrows with components (u,v)
%   at the points (x,y).  The matrices X,Y,U,V must all be the same size
%   and contain corresponding position and velocity components (X and Y
%   can also be vectors to specify a uniform grid).  QUIVER automatically
%   scales the arrows to fit within the grid.
%
%   QUIVER(U,V) plots velocity vectors at equally spaced points in
%   the x-y plane.
%
%   QUIVER(U,V,S) or QUIVER(X,Y,U,V,S) automatically scales the 
%   arrows to fit within the grid and then stretches them by S.  Use
%   S=0 to plot the arrows without the automatic scaling.
%
%   H = QUIVER(...) returns a vector of line handles.

	% Arrow head parameters
	alpha = 0.33;		% Size of arrow head relative to the length of the vector
	beta = 0.33;		% Width of the base of the arrow head relative to the length
	autoscale = 1;		% Autoscale if ~= 0 then scale by this.

	nin = nargin - 1;

	% Check numeric input arguments
	if (nin < 4)					% quiver(u,v) or quiver(u,v,s)
		[msg,x,y,u,v] = xyzchk(varargin{1:2});
	else
		[msg,x,y,u,v] = xyzchk(varargin{1:4});
	end
	if ~isempty(msg), error(msg); end

	if (nin == 3 || nin == 5)		% quiver(u,v,s) or quiver(x,y,u,v,s)
		autoscale = varargin{nin};
	end

	% Scalar expand u,v
	if (numel(u) == 1),     u = u(ones(size(x))); end
	if (numel(v) == 1),     v = v(ones(size(u))); end

	if autoscale,
		% Base autoscale value on average spacing in the x and y
		% directions.  Estimate number of points in each direction as
		% either the size of the input arrays or the effective square
		% spacing if x and y are vectors.
		if min(size(x))==1, n=sqrt(numel(x)); m=n; else [m,n]=size(x); end
		delx = diff([min(x(:)) max(x(:))])/n;
		dely = diff([min(y(:)) max(y(:))])/m;
		del = delx.^2 + dely.^2;
		if (del > 0)
			len = sqrt((u.^2 + v.^2)/del);
			maxlen = max(len(:));
		else
			maxlen = 0;
		end
		
		if maxlen > 0
			autoscale = autoscale*0.9 / maxlen;
		else
			autoscale = autoscale*0.9;
		end
		u = u*autoscale; v = v*autoscale;
	end

	% Make velocity vectors
	x = x(:).';		y = y(:).';
	u = u(:).';		v = v(:).';
	uu = [x;x+u;repmat(NaN,size(u))];
	vv = [y;y+v;repmat(NaN,size(u))];

	h1 = line('XData',uu(:), 'YData',vv(:), 'Parent',hAx, 'Color','k');

	% Make arrow heads and plot them
	hu = [x+u-alpha*(u+beta*(v+eps));x+u; ...
		x+u-alpha*(u-beta*(v+eps));repmat(NaN,size(u))];
	hv = [y+v-alpha*(v-beta*(u+eps));y+v; ...
		y+v-alpha*(v+beta*(u+eps));repmat(NaN,size(v))];
	h2 = line('XData',hu(:), 'YData',hv(:), 'Parent',hAx, 'Color','k');

	if (nargout > 0),	hh = [h1;h2];	end


% --- Creates and returns a handle to the GUI figure. 
function view_anuga_LayoutFcn(h1)

set(h1,...
'Position',[520 548 350 330],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','View ANUGA',...
'NumberTitle','off',...
'PaperSize',[20.98404194812 29.67743169791],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 204 331 38],...
'Style','frame',...
'Tag','frame3');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 242 331 35],...
'Style','frame',...
'Tag','frame2');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[299 177 40 15],...
'String','Slice nº',...
'Style','text');

uicontrol('Parent',h1,...
'Position',[10 159 291 17],...
'BackgroundColor',[0.99 0.99 0.993],...
'Call',{@main_uiCB,h1,'slider1_CB'},...
'Enable','inactive',...
'Style','slider',...
'TooltipString','Use this slider to update the displayed data slice',...
'Tag','slider1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_swwName_CB'},...
'HorizontalAlignment','left',...
'Position',[10 288 310 21],...
'Style','edit',...
'TooltipString','Name of .sww file',...
'Tag','edit_swwName');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_swwName_CB'},...
'Position',[320 288 21 21],...
'TooltipString','Browse for a sww file name',...
'Tag','push_swwName');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[10 311 100 17],...
'String','Input sww file',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 179 170 15],...
'String','Scale color to global min/max',...
'Style','checkbox',...
'TooltipString','If checked, image is obtained by scaling between global min/max. Otherwise use Slice''s min/max',...
'Value',0,...
'Tag','checkbox_globalMinMax');

uicontrol('Parent',h1, 'Position',[180 179 100 15],...
'Call',{@main_uiCB,h1,'checkbox_splitDryWet_CB'},...
'String','Split Water/Land',...
'Style','checkbox',...
'TooltipString','If checked, water and land parts of the image are built separately - Nice with shadings.',...
'Value',1,...
'Tag','checkbox_splitDryWet');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_showSlice_CB'},...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','bold',...
'Position',[240 7 100 23],...
'String','Show slice',...
'TooltipString','Extract the slice selected in "Slice nº" and shot it in a Mirone window',...
'Tag','push_showSlice');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 40 331 97],...
'Style','frame',...
'Tag','frame1');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[30 130 121 15],...
'String','Griding Line Geometry',...
'Style','text',...
'Tag','text_GLG');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_x_max_CB'},...
'HorizontalAlignment','left',...
'Position',[126 98 75 21],...
'Style','edit',...
'Tag','edit_x_max');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_y_max_CB'},...
'HorizontalAlignment','left',...
'Position',[126 72 75 21],...
'Style','edit',...
'Tag','edit_y_max');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[18 102 30 15],...
'String','X Dir',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[17 76 30 15],...
'String','Y Dir',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[146 119 41 13],...
'String','Max',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[67 119 41 13],...
'String','Min',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_x_inc_CB'},...
'HorizontalAlignment','left',...
'Position',[206 98 71 21],...
'Style','edit',...
'TooltipString','DX grid spacing',...
'Tag','edit_x_inc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_y_inc_CB'},...
'HorizontalAlignment','left',...
'Position',[206 72 71 21],...
'Style','edit',...
'TooltipString','DY grid spacing',...
'Tag','edit_y_inc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_Ncols_CB'},...
'HorizontalAlignment','left',...
'Position',[282 98 50 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Ncols');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_Nrows_CB'},...
'HorizontalAlignment','left',...
'Position',[282 72 50 21],...
'Style','edit',...
'TooltipString','Number of rows in the grid',...
'Tag','edit_Nrows');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[223 121 41 13],...
'String','Spacing',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[285 121 51 13],...
'String','# of lines',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[0.83137 0.81569 0.78431],...
'Call',{@main_uiCB,h1,'push_Help_CB'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[275 47 61 18],...
'String','?',...
'Tag','push_Help');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_x_min_CB'},...
'HorizontalAlignment','left',...
'Position',[45 98 75 21],...
'Style','edit',...
'Tag','edit_x_min');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_y_min_CB'},...
'HorizontalAlignment','left',...
'Position',[45 73 75 21],...
'Style','edit',...
'Tag','edit_y_min');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_showMesh_CB'},...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','bold',...
'Position',[120 7 100 23],...
'String','Show mesh',...
'TooltipString','Show the mesh triangulation used in simulation',...
'Tag','push_showMesh');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_Cancel_CB'},...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[20 7 80 23],...
'String','Cancel',...
'Tag','push_Cancel');

uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[140 312 200 16],...
'String','Info',...
'Style','text',...
'Tag','text_Info');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_sliceNumber_CB'},...
'Enable','inactive',...
'Position',[300 158 40 20],...
'String','1',...
'Style','edit',...
'Tag','edit_sliceNumber');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'radio_stage_CB'},...
'Position',[17 252 50 15],...
'String','Stage',...
'Style','radiobutton',...
'TooltipString','Plot the "stage" variable',...
'Value',1,...
'Tag','radio_stage');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'radio_xmoment_CB'},...
'Position',[116 252 85 15],...
'String','Xmomentum',...
'Style','radiobutton',...
'TooltipString','Plot the "xmomentum" variable',...
'Tag','radio_xmoment');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'radio_ymoment_CB'},...
'Position',[243 252 80 15],...
'String','Ymomentum',...
'Style','radiobutton',...
'TooltipString','Plot the "ymomentum" variable',...
'Tag','radio_ymoment');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'popup_derivedVar_CB'},...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[111 211 190 22],...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_derivedVar');

uicontrol('Parent',h1, 'Position',[308 211 23 23],...
'Call',{@main_uiCB,h1,'toggle_vel_CB'},...
'Style','togglebutton',...
'TooltipString','Plot arrow field when this button is depressed',...
'Vis', 'off', ...
'Tag','toggle_vel');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'check_derivedVar_CB'},...
'Position',[18 214 80 15],...
'String','Derived var',...
'Style','checkbox',...
'TooltipString','Select a derived quantity from the side popup menu',...
'Tag','check_derivedVar');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[40 269 100 15],...
'String','Primary quantities',...
'Style','text',...
'Tag','text_Pq');

function main_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

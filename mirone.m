function varargout = mirone(varargin)
%	MIRONE, by itself, creates a window bar from which you load a lot of grid/images formats
%	MIRONE(FNAME) opens the file FNAME and displays it
%	H = MIRONE(...) returns the handle to a new mirone window
%
%	mirone('CALLBACK',handles,...) calls the local function named CALLBACK with the given input arguments.

%	Copyright (c) 2004-2016 by J. Luis
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

% $Id: mirone.m 7919 2016-05-31 15:55:41Z j $

	if (nargin > 1 && ischar(varargin{1}))
		if ( ~isempty(strfind(varargin{1},':')) || ~isempty(strfind(varargin{1},filesep)) )
			% Very likely called with a filename with those horrendous stupid blanks
			for (k = 1:nargin-1),	varargin{k}(end+1) = ' ';	end
			varargin = {[varargin{:}]};			h = mirone_OpeningFcn(varargin{:});

			if (nargout),	varargout{1} = h;		end
		else
			gui_Callback = str2func(varargin{1});
			[varargout{1:nargout}] = feval(gui_Callback,varargin{2:end});
		end
	else
		h = mirone_OpeningFcn(varargin{:});
		if (nargout),	varargout{1} = h;		end
	end

% --------------------------------------------------------------------------------------------------
function hObject = mirone_OpeningFcn(varargin)
% PRAGMA SECTION (It's far far from clear when files must be declared here)
%#function uigetfolder_standalone mapproject_m grdproject_m coordinate_system gmtmbgrid_m gmtedit
%#function nearneighbor_m cpt2cmap grdfilter_m grdgradient_m grdsample_m grdtrack_m grdtrend_m 
%#function grdutils scaleto8 bpass3d inv3d rtp3d syn3d igrf_m surface_m
%#function range_change swan tsun2 mansinha_m deform_mansinha deform_okada dim_funs
%----- These are for image
%#function grayto8 grayto16 grayxform imfilter_mex imhistc imlincombc parityscan intlutc ordf
%#function imreconstructmex applylutc bwboundariesmex bwlabel1 bwlabel2 ditherc bwlabelnmex
%----- These are in utils
%#function tabpanelfcn degree2dms dms2degree dec2deg dec_year ivan_the_terrible ddewhite string_token
%#function test_dms text_read double2ascii save_seismicity jd2date guess_file shading_mat 
%#function getline_j new_frame3D histos_seis uirestore_j uisuspend_j
%----- These is for write_gmt_script
%#function draw_scale time_stamp pscoast_options_Mir paint_option w_option
%----- Those are ..., the hell with explanations for what I don't realy understand. They are needed, that's all.
%#function gmtlist_m country_select read_isf choosebox magbarcode listbox_message add_poles animate_seismicity
%#function get_polygon rot_euler datums telha_m find_clusters fft_stuff select_cols uistack_j smoothing_param
%#function patch_meca ui_edit_patch_special bands_list multibandread_j imscroll_j geog2projected_pts
%#function mltable_j iptcheckinput resampsep wgifc telhometro vitrinite edit_line move_obj make_arrow
%#function edit_track_mb save_track_mb houghmex qhullmx writegif mpgwrite cq helpdlg
%#function move2side aguentabar gdal_project gdalwarp_mex poly2mask_fig url2image calc_bonin_euler_pole spline_interp
%#function mat2clip buffer_j PolygonClip trend1d_m akimaspline shake_mex ground_motion wms_tool microlev
%#function write_esri_hdr distmin mag_synthetic image_histo write_gmt_symb mkpj decompress mosaicer
%#function lasreader_mex laszreader_mex escorrega show_manguito travel thresholdit intersections nswing runCB_tintol
%#function usgs_recent_seismicity sat_orbits uisetdate doy
%#function c_cpt2cmap c_grdfilter c_grdinfo c_grdlandmask c_grdproject c_grdread c_grdsample
%#function c_grdtrend c_mapproject c_nearneighbor c_shoredump c_surface popenr

	global gmt_ver;		gmt_ver = 5;	global home_dir;	fsep = filesep;
	toCompile = false;		% To compile set this one to TRUE
	if (toCompile)
		home_dir = cd;
	else
		if (isempty(home_dir))		% First time call. Find out where we are
			home_dir = fileparts(mfilename('fullpath'));			% Get the Mirone home dir and set path
			addpath(home_dir, [home_dir fsep 'src_figs'],[home_dir fsep 'utils']);
			if (exist('OCTAVE_VERSION','builtin') ~= 0)				% This is a repetition of the test later in mirone_uis
				addpath([home_dir fsep 'lib_mex' fsep 'octave' fsep octave_config_info.canonical_host_type]);
			else
				addpath([home_dir fsep 'lib_mex']);
			end
		end
	end

	[hObject,handles,home_dir] = mirone_uis(home_dir);
	pos = get(handles.figure1,'pos');
	if (pos(4) == 32)				% On Win 8 and retina, the bloody OS insists in resetting it to 32
		set(handles.figure1,'Pos',[pos(1), pos(2)+31, pos(3), 1]);
	end

	handles.home_dir = home_dir;
	handles.DefLineThick = 0.5;		% Default line thickness (overwriten by mirone_pref)
	handles.DefLineColor = 'k';		% Default line color (overwriten by mirone_pref)
	handles.DefineMeasureUnit = 'k';% Default measure units to kilometrs (overwriten by mirone_pref)
	handles.grdname = [];			% Contains the name of the current (if it's the case) gmt grid
	handles.nTrack = 0;				% Counter of the number of MB tracks imported
	handles.origFig = [];			% To store the original image copy
	handles.fileName = [];			% To store any input grid/image file name
	handles.image_type = 0;			% 1->grd; 2-> trivial (jpg,png,bmp,etc...); 3->GeoTIFF; 4->DEMs; 20-> white bg
	handles.transparency = [];		% Will hold the word 'alpha' when reading a RGBA image
	handles.computed_grid = 0;		% However, matrices with a gmt header will have this == 1, so that they can be saved
	handles.no_file = 1;			% 0 means a grid is loaded and 1 that it is not (to test when icons can be pushed)
	handles.geog = 1;				% By default grids are assumed to be in geographical coordinates
	handles.swathRatio = 3;			% Default swath width / water depth ratio for multibeam planing
	handles.grdMaxSize = 524288000;	% I use this for limiting the grid size that is stored in RAM (50 Mb)
	handles.firstMBtrack = 1;		% Used for knowing whether to display or not the MB planing info message in "start planing"
	handles.EarthRad = 6371.005076;	% Authalic radius
	handles.maregraphs_count = 0;	% Counter of maregraphs (tsunami modeling)
	handles.Illumin_type = 0;		% 0 -> no illumination; 1 -> grdgradient; 2 -> grdgradient Lambertian; 4 -> Lambertian;
	handles.zoom_state = 0;			% Flag to signal if zoom state is to be re-activated (0 means no re-activation)
	handles.bg_color = [1 1 1];		% (default) Backgoround color used when grid values == NaN
	handles.which_cont = [];		% To store the contour levels (used to not draw repeated contours)
	handles.have_nans = 0;			% Used to know if the grids have NaNs
	handles.is_draped = false;		% Used to know if the image comes from draping
	handles.was_int16 = false;		% Keep track of short int imported grids
	handles.Nodata_int16 = [];		% To store Nodata of short int grids
	handles.ForceInsitu = false;	% Use "insitu" grid importing (transposition)
	handles.DefineEllipsoide = [6378137, 0, 1/298.2572235630];	% Defaults to WGS-84
	handles.path_tmp = [home_dir fsep 'tmp' fsep];
	handles.path_data = [home_dir fsep 'data' fsep];
	handles.last_directories = {handles.path_tmp; home_dir};	% Let it have something existent
	handles.hImg = [];				% To hold the image handle
	handles.firstIllum = 1;			% First illumination will use the displayed image which may have been IP
	handles.flederBurn = 1;			% When build a fleder obj burn eventual coastlines
	handles.flederPlanar = 1;		% Default to planar flder objs (but mirone_pref knows better)
	handles.whichFleder = 1;		% whichFleder = 1 for the free iview4d or 0 for the true thing (fledermaus)
	handles.oldSize = [get(hObject,'Pos'); get(hObject,'Pos')];		% Duplicate so that we store ORIGINAL size
	if (handles.oldSize(1,4) == 0),	handles.oldSize(1,4) = 1;		end
	handles.is_projected = false;	% To keep track if coords are projected or not
	handles.defCoordsIn = 0;		% To use when Load files and have to decide if we need to project
									% 0 -> don't know; -1 -> Coords are already projected; 1 -> geog coords needing project
	try		handles.Projections;			% Use a try/catch since isfield is brain-dead long
	catch,	handles.Projections = nan;		% Set it to something to prevent "unknown field" error
	end
	handles.scale2meanLat = 1;		% Scale geog ref images so that at middle image dx = dy in cartesian units (kind of proj)
	handles.FOpenList = cell(numel(handles.RecentF),1);
	handles.withSliders = true;		% Set Zoom sliders
	handles.validGrid = false;		%
	handles.nLayers = 1;			% If > 1 after reading a netCDF file call aquamoto
	handles.deflation_level = 0;	% If > 0 will create compressed netCDF-4 files
	handles.is_defRegion = false;	% A def region is a particular case to create GMT custom symbols.

	try							% A file named mirone_pref.mat contains the preferences, read them from it
		prf = load([handles.path_data 'mirone_pref.mat']);
		handles.geog = prf.geog;
		handles.grdMaxSize = prf.grdMaxSize;				% 2^20 = 1 Mb
		handles.swathRatio = prf.swathRatio;
		handles.last_directories = prf.directory_list;
		handles.DefLineThick = sscanf(strtok(prf.DefLineThick{1}),'%f');
		% Decode the line color string into the corresponding char (e.g. k,w, etc...)
		if (strcmp(prf.DefLineColor{1},'Black')),	handles.DefLineColor = 'k';
		else										handles.DefLineColor = lower(prf.DefLineColor{1}(1));
		end
		% Decode the Measure unities into a char code (e.g. n, k, m, u from {'nautical miles' 'kilometers' 'meters' 'user'})
		handles.DefineMeasureUnit = prf.DefineMeasureUnit{1}(1);
		handles.DefineEllipsoide = prf.DefineEllipsoide_params;		% Set the default ellipsoide parameters (a,b,f)
		handles.flederBurn = prf.flederBurn;
		handles.flederPlanar = prf.flederPlanar;
		handles.scale2meanLat = prf.scale2meanLat;
		handles.FOpenList = prf.FOpenList;
		handles.whichFleder = prf.whichFleder;
		if (~prf.moveDoubleClick)							% this info is used by UI_EDIT_POLYGON()
			setappdata(handles.axes1,'MovPolyg','extend')	% Move lines with a Shift-click drag-n-drop
		end
		handles.bg_color = prf.nanColor;
		handles.deflation_level = prf.deflation_level;
		try		gmt_ver = prf.gmt_ver;		% Hope that this will be a temporary thing till GMT5 is fully working
		catch,	gmt_ver = 4;
		end
	catch
		% Tell mirone_pref to write up the defaults.
		mirone_pref(handles,'nikles')
		handles = guidata(handles.figure1);					% And need also the updated handles
		gmt_ver = 4;										% Default to GMT4
	end

	j = false(1,numel(handles.last_directories));			% vector for eventual cleaning non-existing dirs
	for (i = 1:numel(handles.last_directories))				% Check that all dirs in last_directories exist
		j(i) = (exist(handles.last_directories{i},'dir') ~= 7);
	end
	handles.last_directories(j) = [];						% clean non-existing directories

	if (isempty(handles.last_directories))					% Don't ever let it be empty
		handles.last_directories = {handles.path_tmp; home_dir};	% Let it have something existent
	end
	if (any(j))					% If any of the old dirs evaporated, update that info in mirone_prefs
		directory_list = handles.last_directories;	isempty(directory_list);	% To shut up MLint
		if (handles.version7),		save([handles.path_data 'mirone_pref.mat'],'directory_list', '-append', '-v6')
		else						save([handles.path_data 'mirone_pref.mat'],'directory_list', '-append')
		end
	end

	handles.work_dir = handles.last_directories{1};
	if (handles.work_dir ~= fsep),		handles.work_dir = [handles.work_dir fsep];		end
	handles.last_dir = handles.last_directories{1};			% Initialize last_dir to work_dir
	setappdata(hObject,'swathRatio',handles.swathRatio);	% I need this for getline_mb

	% Change the MeasureDistance label to the selected (in prefs) units
	set(handles.ToolsMeasureDist,'Label',['Distance in ' handles.DefineMeasureUnit])

	% -------------------------- Find in which mode Mirone was called ----------------------------------
	drv = [];	grd_data_in = 0;	grd_data_interfero = 0;		pal = [];	win_name = 'Mirone';
	if ~isempty(varargin)
		n_argin = nargin;
		if (n_argin == 1 && ischar(varargin{1}))				% Called with a file name as argument
			[pato, fname, EXT] = fileparts(varargin{1});		% Test to check online command input
			if (isempty(pato)),			varargin{1} = [cd fsep fname EXT];	end
			[drv, algures] = aux_funs('findFileType',varargin{1});
			if (ischar(algures)),		varargin{1} = algures;	end 		% File exists but not in Mirone's root dir
			if (ischar(algures) || algures),	handles.fileName = varargin{1};		end		% Can be added by recentFiles
		elseif (isa(varargin{1},'uint8') || isa(varargin{1},'int8') || islogical(varargin{1}))
			% Called with an image as argument and optionaly an struct header (& geog, name, cmap optional fields)
			if (isa(varargin{1},'int8'))		% We cannot represent a int8 image. Do something
				varargin{1} = uint8(cvlib_mex('addS', int16(varargin{1}), 128));	% [-128 127] -> [0 255]
			end
			% Now deal with the case of a eventual multiband ( > than 3 planes) array
			if (size(varargin{1},3) > 3),		aux_funs('toBandsList', handles.figure1, varargin{1}, 'multiband array'),	end

			isReferenced = false;
			if (n_argin == 2 && isa(varargin{2},'struct'))		% An image with coordinates
				tmp = varargin{2};
				handles.head = tmp.head;		X = tmp.X;		Y = tmp.Y;
				handles.image_type = 3;			axis_t = 'xy';
				if (isfield(tmp,'geog')),		handles.geog = tmp.geog;	end % Prevails over the guess in show_image (does it?)
				if (isfield(tmp,'cmap')),		set(handles.figure1,'Colormap',tmp.cmap);	end
				if (isfield(tmp,'name')),		win_name = tmp.name;	end
				if (isfield(tmp,'srsWKT'))
					aux_funs('appP', handles, tmp.srsWKT)			% If we have a WKT proj, store it
					isReferenced = true;
					if (~handles.geog),		handles.is_projected = true;	end		% WEAK LOGIC. SHOULD PARSE WKT TO MAKE SURE
				end
			else
				X = [];			Y = [];			win_name = 'Cropped_image';
				handles.image_type = 2;			handles.geog = 0;		axis_t = 'off';
				handles.head = [1 size(varargin{1}, 2) 1 size(varargin{1}, 1) 0 255 0 1 1];		% Fake a grid reg GMT header
				if (ndims(varargin{1}) == 2),	set(handles.figure1,'Colormap',gray(256));	end
				pal = getappdata(0,'CropedColormap');				% See if we have a colormap to use here
				if (~isempty(pal)),		set(handles.figure1,'Colormap',pal),	rmappdata(0,'CropedColormap'),	end
				setappdata(hObject,'Croped','yes');					% ???
			end

			if (size(varargin{1},3) == 4),	handles.transparency = 'alpha';		end		% Signal show_image that 4th band is transparency

			handles = show_image(handles,win_name,X,Y,varargin{1},0,axis_t,handles.head(7),1);
			if (~isReferenced),		grid_info(handles,[],'iminfo',varargin{1});			% Create a info string
			else					grid_info(handles,tmp.srsWKT,'referenced',varargin{1});
			end
			handles = aux_funs('isProj',handles);				% Check/set about coordinates type
			
		elseif (n_argin == 1 && isa(varargin{1},'struct') && isfield(varargin{1},'ProjectionRefPROJ4'))
			% A GMT5 grid/image structure. (for images we still do not use eventual alpha channel)
			handles.head = [varargin{1}.range varargin{1}.registration varargin{1}.inc];
			if (~isfield(varargin{1}, 'image'))
				Z = varargin{1}.z;			grd_data_in = true;
				if (~isa(Z,'single')),		Z = single(Z);		end
				handles.have_nans = grdutils(Z,'-N');
				X = linspace(varargin{1}.range(1), varargin{1}.range(2), varargin{1}.n_columns);
				Y = linspace(varargin{1}.range(3), varargin{1}.range(4), varargin{1}.n_rows);
			else
				handles.image_type = 3;		isReferenced = false;
				X = [varargin{1}.x(1) varargin{1}.x(end)];		Y = [varargin{1}.y(1) varargin{1}.y(end)];
				handles.geog = aux_funs('guessGeog', [Y Y]);
				ProjectionRefWKT = varargin{1}.ProjectionRefWKT;
				if (isempty(ProjectionRefWKT) && ~isempty(varargin{1}.ProjectionRefPROJ4))
					ProjectionRefWKT = ogrproj(varargin{1}.ProjectionRefPROJ4);
				end
				if (~isempty(ProjectionRefWKT))
					aux_funs('appP', handles, varargin{1}.ProjectionRefWKT)			% If we have a WKT proj, store it
					isReferenced = true;
					if (~handles.geog),		handles.is_projected = true;	end		% WEAK LOGIC. SHOULD PARSE WKT TO MAKE SURE
				end
				if (ndims(varargin{1}.image) == 2),		set(handles.figure1,'Colormap',gray(256));		end
				win_name = 'Nikles';
				if (~isempty(varargin{1}.title)),		win_name = varargin{1}.title;	end
				handles = show_image(handles,win_name,X,Y,varargin{1}.image,0,'xy',varargin{1}.registration,1);
				if (~isReferenced),		grid_info(handles,[],'iminfo',varargin{1}.image);		% Create a info string
				else					grid_info(ProjectionRefWKT, 'referenced', varargin{1}.image);
				end
				handles = aux_funs('isProj',handles);				% Check/set about coordinates type
			end

		elseif (n_argin == 2 && isequal([size(varargin{2},1) size(varargin{2},2)],[1 9]))
			% A matrix with the classic nine elements header vector
			handles.head = varargin{2};
			Z = varargin{1};			grd_data_in = true;
			if (~isa(Z,'single')),		Z = single(Z);		end
			handles.have_nans = grdutils(Z,'-N');
			X = linspace(handles.head(1), handles.head(2), size(Z,2));
			Y = linspace(handles.head(3), handles.head(4), size(Z,1));

		elseif (n_argin < 4 && ~(isa(varargin{1},'uint8') || isa(varargin{1},'int8')))
			% A matrix. Treat it as if it is a gmt grid. No error testing on the grid head descriptor
			Z = varargin{1};			grd_data_in = true;
			if (~isa(Z,'single')),		Z = single(Z);		end
			handles.have_nans = grdutils(Z,'-N');
			if (numel(varargin) == 2 && isa(varargin{2},'struct'))		% An grid with a header
				tmp = varargin{2};
				handles.head = tmp.head;	X = tmp.X;	Y = tmp.Y;
				if (isfield(tmp,'name')),	win_name = tmp.name;	end		% All calls should transmit a name, but ...
				if (isfield(tmp,'cmap')),	pal = tmp.cmap;			end
				if (isfield(tmp,'was_int16'))
					handles.was_int16 = tmp.was_int16;		handles.Nodata_int16 = tmp.Nodata_int16;
				end
				if (isfield(tmp,'srsWKT'))
					grid_info(handles,tmp.srsWKT,'referenced',varargin{1});	% Create a info string
					aux_funs('appP', handles, tmp.srsWKT)					% We have a WKT proj, store it
					handles.is_projected = true;		% WEAK LOGIC. SHOULD PARSE WKT TO MAKE SURE
				elseif (isfield(tmp,'ProjGMT'))			% From geog_calculator. Has opt_J.
					projection_menu(handles, tmp.ProjGMT)
					handles = guidata(hObject);			% Get the updated version changed in the above call
				end
			else
				zz = grdutils(Z,'-L');
				handles.head = [1 size(Z,2) 1 size(Z,1) zz(1) zz(2) 0 1 1];
				X = 1:size(Z,2);			Y = 1:size(Z,1);
			end
		elseif (n_argin == 4 && isnumeric(varargin{1}) && isa(varargin{2},'struct') && ...
				(strcmp(varargin{3},'Deformation') || strcmp(varargin{3},'Interfero')) )
			% A matrix. Treat it as if it'is a gmt grid. No error testing on the grid head descriptor
			% Note: this is a special case of the situation above that will be used to identify this figure
			% as an Okada deformtion data (via its Name). This info is searched by the tsunami modeling option
			if (~isa(varargin{1},'single')),		varargin{1} = single(varargin{1});		end
			handles.head = varargin{2}.head;		X = varargin{2}.X;		Y = varargin{2}.Y;
			Z = varargin{1};
			handles.have_nans = grdutils(Z,'-N');
			if (varargin{3}(1) == 'D')
				grd_data_in = true;					win_name = 'Okada deformation';
				setappdata(hObject,'hFigParent',varargin{4});
			else						% A matrix input containing an interfeogram with cdo == varargin{4}
				grd_data_interfero = true;			win_name = 'Interferogram';
				cdo = varargin{4};
			end
		end
	end

	% The following IF cases deal only with cases where a grid was given in argument
	if (grd_data_in || grd_data_interfero)
		handles.image_type = 1;		handles.computed_grid = 1;	% Signal that this is a computed grid
		if (isempty(pal)),			pal = jet(256);		end
		if (grd_data_interfero)			% Interferogram grid
			pal = load([handles.path_data 'gmt_other_palettes.mat'],'circular');	pal = pal.circular;
			zz = uint8(abs(rem(double(Z),cdo)/cdo)*255);
		else
			zz = scaleto8(Z);
		end
		set(handles.figure1,'Colormap',pal)
		aux_funs('StoreZ',handles,X,Y,Z)			% If grid size is not to big we'll store it
		aux_funs('colormap_bg',handles,Z,pal);
		handles = show_image(handles,win_name,X,Y,zz,1,'xy',handles.head(7));
	end

	% ---------+================+---------- FIGURE VISIBLE HERE ----------+==============+-----------
	set(handles.figure1,'Vis', 'on')
	% ---------+================+---------- FIGURE VISIBLE HERE ----------+==============+-----------

	handles.IAmAMac = strncmp(computer,'MAC',3);
	setappdata(0,'IAmAMac',handles.IAmAMac)
	if (handles.IAmAMac)
		if (isempty(getappdata(0,'have_DYLD')))		% Deal with MacOSX blindness
			%set_gmt(['DYLD_LIBRARY_PATH=' cd ':'],'DYLD_LIBRARY_PATH')		% Prepend current dir to DYLD_LIBRARY_PATH
			DYLD = getenv('DYLD_LIBRARY_PATH');		% Available only on > R14?? versions 
			setenv('DYLD_LIBRARY_PATH',[DYLD ':' cd])
			setappdata(0,'have_DYLD',true)			% Signal to do this only once
		end
		% F. TMW just cannot make things work decently on Macs. Labels are bigger than reserved space
		% And now they screw it even more. The following simply f... on 2012 (and maybe earlier(
		v = version('-release');
		if (str2double(v(1:4)) < 2012)
			set(handles.File,  'Label', ['<HTML><FONT size="3">' get(handles.File, 'Label') '</Font></html>'])
			set(handles.Image, 'Label', ['<HTML><FONT size="3">' get(handles.Image, 'Label') '</Font></html>'])
			set(handles.Tools, 'Label', ['<HTML><FONT size="3">' get(handles.Tools, 'Label') '</Font></html>'])
			set(handles.Draw,  'Label', ['<HTML><FONT size="3">' get(handles.Draw, 'Label') '</Font></html>'])
			set(handles.Geography, 'Label', ['<HTML><FONT size="3">' get(handles.Geography, 'Label') '</Font></html>'])
			set(handles.Plates,'Label', ['<HTML><FONT size="3">' get(handles.Plates, 'Label') '</Font></html>'])
			set(handles.MagGrav,  'Label', ['<HTML><FONT size="3">' get(handles.MagGrav, 'Label') '</Font></html>'])
			set(handles.Seismology, 'Label', ['<HTML><FONT size="3">' get(handles.Seismology, 'Label') '</Font></html>'])
			set(handles.Tsunamis, 'Label', ['<HTML><FONT size="3">' get(handles.Tsunamis, 'Label') '</Font></html>'])
			set(handles.GMT, 'Label', ['<HTML><FONT size="3">' get(handles.GMT, 'Label') '</Font></html>'])
			set(handles.GridTools, 'Label', ['<HTML><FONT size="3">' get(handles.GridTools, 'Label') '</Font></html>'])
			set(handles.Projections, 'Label', ['<HTML><FONT size="3">' get(handles.Projections, 'Label') '</Font></html>'])
			set(handles.Help, 'Label', ['<HTML><FONT size="3">' get(handles.Help, 'Label') '</Font></html>'])
		end
		set(hObject,'DockControls','off')
	end

	%Find out which gmt version is beeing used. 
	info = getappdata(0,'gmt_version');		% See if the info is already there.
	if (isempty(info))
		info = set_gmt(['GMT_USERDIR=' home_dir fsep 'gmt_userdir']);
		setappdata(0,'gmt_version',info);	% Save it so that the next time a new mirone window is opened
	end
	if (~strcmp(info.full, 'y'))
		set([handles.CoastLineFull handles.PBFull handles.RiversFull], 'Enable','off')
	end
	if (~strcmp(info.high, 'y'))
		set([handles.CoastLineHigh handles.PBHigh handles.RiversHigh], 'Enable','off')
	end
	if (~strcmp(info.intermediate, 'y'))
		set([handles.CoastLineInterm handles.PBInterm handles.RiversInterm], 'Enable','off')
	end
	if (~strcmp(info.low, 'y'))
		set([handles.CoastLineLow handles.PBLow handles.RiversLow], 'Enable','off')
	end
	if (~strcmp(info.crude, 'y'))
		set([handles.CoastLineCrude handles.PBCrude handles.RiversCrude], 'Enable','off')
	end

	% Deal with the new (big) troubles introduced by using GMT5.2 that needs to know where to find its own share dir
	if (gmt_ver == 5)		% For GMT5 we have a shity highly police control on sharedir. Use this to cheat it.
		t = set_gmt('GMT5_SHAREDIR', 'whatever');		% Inquire if GMT5_SHAREDIR exists 
		if (isempty(t) || (~(exist(t, 'dir') == 7)))	% Test also that the dir exists
			% If not, set a fake one with minimalist files so that GMT does not complain/errors
			% We have to use GMT5_SHAREDIR and NOT GMT_SHAREDIR because it's the first one checked in gmt_init.c/GMT_set_env()
			set_gmt(['GMT5_SHAREDIR=' home_dir fsep 'gmt_share']);
		end
	end

	% See if we have a TMP dir set by a ENV var that will take precedence over the default one
	t = set_gmt('MIRONE_TMP', 'whatever');			% Inquire if MIRONE_TMP exists
	if (~isempty(t))		% Now check that the dir exists
		if (exist(t, 'dir') == 7)
			handles.path_tmp = [t fsep];
		end
	end

	guidata(hObject, handles);
	tmp.home_dir = home_dir;	tmp.work_dir = handles.work_dir;	tmp.last_dir = handles.last_dir;
	setappdata(0,'MIRONE_DIRS',tmp);		% To access from places where handles.home_dir is unknown (must precede gateLoadFile())
	if (~isempty(drv))
		gateLoadFile(handles,drv,varargin{1});		% load recognized file types
	else
		recentFiles(handles,[]);					% Just make the "Recent files" entry available
	end

	if (~isempty(handles.hImg))						% If we had something in input, check the coordinates type
		handles = aux_funs('isProj', guidata(hObject));
		setAxesDefCoordIn(handles,1);
	end

	set_gmt(['PROJ_LIB=' home_dir fsep 'data' fsep 'proj_lib']);		% For projections with GDAL
	set_gmt(['GDAL_DATA=' home_dir fsep 'data' fsep 'gdal_data']);
	set_gmt(['GEOTIFF_CSV=' home_dir fsep 'data' fsep 'gdal_data']);

 	try		custom_menus(hObject, handles.path_data),	end			% Add any custom menus as commanded by OPTcontrol.txt

% --------------------------------------------------------------------------------------------------
function erro = gateLoadFile(handles,drv,fname)
% Gateway function to load a recognized file type using its name
	erro = 0;
	if (strncmp(drv, 'MB', 2)),	drv_or = drv;	drv = 'MB';		end		% The MB (system) type case are actually several
	switch drv
		case 'gmt',			loadGRID(handles, fname, 'GMT_relatives')
		case 'generic',		FileOpenNewImage_CB(handles, fname);
		case 'geotif',		FileOpenGeoTIFF_CB(handles, 'nikles', fname);
		case 'ecw',			FileOpenGeoTIFF_CB(handles, 'ecw', fname);		% A particular case (includes jp2)
		case 'multiband',	FileOpenGDALmultiBand_CB(handles, 'AVHRR', fname);
		case 'envherd',		FileOpen_ENVI_Erdas_CB(handles, [], fname);
		case 'mat',			FileOpenSession_CB(handles, fname)
		case 'mola',		loadGRID(handles, fname, 'MOLA');
		case 'cpt',			color_palettes(fname);
		case 'dat',			load_xyz(handles, fname);
		case 'ncshape',		load_xyz(handles, fname, drv);
		case 'shp',			DrawImportShape_CB(handles, fname);
		case 'ogr',			DrawImportOGR_CB(handles, fname);
		case 'las',			read_las(handles, fname);
		case 'mgg_gmt',		GeophysicsImportGmtFile_CB(handles, fname);
		case 'ghost',		load_ps(handles, fname);		% Put ghostscript on works
		case 'sww',			aquamoto(fname);
		case 'MB',          load_MB(handles, fname, drv_or);
		case 'dono',		erro = FileOpenGeoTIFF_CB(handles,'dono',fname);	% It means "I don't know"
		otherwise,			erro = 1;
	end
	if (erro),		warndlg(['Sorry but couldn''t figure out what to do with the ' fname ' file'],'Warning'),	end

% --------------------------------------------------------------------------------------------------
function handles = recentFiles(handles, opt)
	if (~isempty(handles.fileName) && ~exist(handles.fileName,'file')),		return,		end		% For example, subdatasets
	jump = false;		N = min(numel(handles.FOpenList), numel(handles.RecentF));
	if (nargin == 1 && ~isempty(handles.fileName))		% Update list
		for (i = 1:N)		% See if new name is already in FOpenList
			if (strcmpi(handles.fileName, handles.FOpenList{i}))
				for (k = i:-1:2)	% Yes, it is. So move it to the top
					handles.FOpenList{k} = handles.FOpenList{k-1};
				end
				handles.FOpenList{1} = handles.fileName;
				jump = true;	break
			end
		end
		if (~jump)				% We got a new name
			for (i = N:-1:2)	% Make room for the new name at top of the list 
				handles.FOpenList{i} = handles.FOpenList{i-1};
			end
			handles.FOpenList{1} = handles.fileName;
		end
	end
	if (isempty(handles.fileName)),	jump = true;	end
	if ( ~jump && (nargin == 1 || (nargin == 2 && ~isempty(opt))) )		% Save only if it worth it
		FOpenList = handles.FOpenList;		fname = [handles.path_data 'mirone_pref.mat'];
		if (~handles.version7),		save(fname,'FOpenList','-append')	% Update the list for "Recent files"
		else						save(fname,'FOpenList','-append', '-v6')
		end
	end
	for (i = 1:N)			% The ishandle test below is crutial when GCPs
		if (isempty(handles.FOpenList{i}) || ~ishandle(handles.RecentF(i))),	continue,	end
		set(handles.RecentF(i), 'Label',handles.FOpenList{i},'Call',{@openRF,i}, 'Vis','on')
	end
	if (~nargout),	guidata(handles.figure1,handles),	end

% --------------------------------------------------------------------------------------------------
function openRF(obj,event,n)
% Open from 'Recent Files'
	handles = guidata(obj);
	[drv, sim] = aux_funs('findFileType',handles.FOpenList{n});
	if (sim)
		gateLoadFile(handles,drv,handles.FOpenList{n});
	else		% File does not exist. Update the FOpenList
		warndlg(['File: ' handles.FOpenList{n} ' no longer exists'],'Warning')
		handles.FOpenList(n) = [];		handles.FOpenList{end+1} = [];	% Delete missing and add one at the end
		delete(handles.RecentF(n)),		guidata(handles.figure1,handles)
	end

% --------------------------------------------------------------------------------------------------
function handles = SetAxesNumericType(handles,event)
% Save original X & Y labels in appdata for easear access when we want to change them
	setappdata(handles.axes1,'XTickOrig',get(handles.axes1,'XTickLabel'))
	setappdata(handles.axes1,'XTickOrigNum',get(handles.axes1,'XTick'))
	setappdata(handles.axes1,'YTickOrig',get(handles.axes1,'YTickLabel'))
	setappdata(handles.axes1,'YTickOrigNum',get(handles.axes1,'YTick'))
	set(handles.axes1, 'FontSize', 9)				% Make this the default
	LFT = 'DegDec';			visibility = 'on';		% For the geog case
	if (~handles.geog),		LFT = 'NotGeog';	visibility = 'off';		end 
	setappdata(handles.axes1,'LabelFormatType',LFT)
	set(handles.LabFormat, 'Vis', visibility)
	if (handles.validGrid), set(handles.PixMode, 'Call', {@PixMode_CB,handles.figure1, true}, 'Vis', 'on')
	else					set(handles.PixMode, 'Vis', 'off')
	end
	set(handles.RCMode, 'Call', {@PixMode_CB,handles.figure1, false})
	set(handles.FancyMode, 'Call', {@FancyMode_CB,handles.figure1})

% --------------------------------------------------------------------------------------------------
function PixMode_CB(hObject, event, hFig, opt)
% Inside each grid cell, which is a pixel on the screen, display only the grid node value
	handles = guidata(hFig);
	if (opt)		% Pixel mode on/off
		if (strcmp(get(hObject,'Checked'),'off'))
			set(hObject,'Checked','on'),	setappdata(hFig,'PixelMode',true)
		else
			set(hObject,'Checked','off'),	setappdata(hFig,'PixelMode',false)
		end
		set(handles.RCMode, 'Checked','off'),		setappdata(hFig,'RCMode',false)		% Put the RowColMode to off
	else			% Row Col mode on/off
		if (strcmp(get(hObject,'Checked'),'off'))
			set(hObject,'Checked','on'),	setappdata(hFig,'RCMode',true)
		else
			set(hObject,'Checked','off'),	setappdata(hFig,'RCMode',false)
		end
		set(handles.PixMode, 'Checked','off'),		setappdata(hFig,'PixelMode',false)	% Put the PixMode to off
	end

% --------------------------------------------------------------------------------------------------
function FancyMode_CB(hObject, evt, hFig)
% Call the function that sets or unsets the FANCY frame depending on the previous state
	handles = guidata(hFig);
	if (strcmp(get(hObject,'Checked'),'off'))
		set(hObject,'Checked','on'),	fancyFrame(handles,'set')
	else
		set(hObject,'Checked','off'),	fancyFrame(handles,'unset')
	end

% --------------------------------------------------------------------------------------------------
function  PlatesAgeLift_CB(handles)
% Apply Parsons & Sclatter relation to compensate sea-bottom age sinking.
% It assumes that bathymetry is the loaded grid (in meters Z up) and age in Ma
	if (handles.no_file),		return,		end
	resp = inputdlg({'Full name (with path) of Age grid:'},'Where is the grid?',[1 80]);
	if (isempty(resp)),		return,		end
	fname = resp{1};
	if (~exist(fname,'file'))
		warndlg('The file name provided does not exist. Bye, Bye','Warning'),	return
	end
	att = gdalread(fname,'-M','-C');
	if (att.GMT_hdr(1) > handles.head(1) || att.GMT_hdr(2) < handles.head(2) || ...
			att.GMT_hdr(3) > handles.head(3) || att.GMT_hdr(4) < handles.head(4))
		errordlg('No way. The Age grid does not cover the limits of the bathymetry grid.','Error'),		return
	end

	[Age, X, Y, srsWKT, miniHandles] = read_grid([], fname, 'GDAL', sprintf('-R%.12f/%.12f/%.12f/%.12f', handles.head(1:4)));
	if (isempty(Age)),	return,		end			% Something bad happened
	lift = single(350 * sqrt(double(Age)));
	[X,Y,Z] = load_grd(handles);
	lift = cvlib_mex('resize', lift, [size(Z,1) size(Z,2)]);
	cvlib_mex('add', lift, Z);
	miniHandles.head(7:9) = handles.head(7:9);		% min/max (5:6) will be updated in GRDdisplay
	GRDdisplay(handles,X,Y,lift,miniHandles.head,[],'AgeLiftedBathymetry',srsWKT)

% --------------------------------------------------------------------------------------------------
function load_ps(handles, fname)
% Load a PostScript file ... via ghostscript and a MEX to access ghost stdout (SHITY QUALITY)
% P6
% # Image generated by Ghostscript (device=ppmraw)
% 612 792
% 255
	P = popenr(['gswin64c -q -r300x300 -sDEVICE=ppmraw -sOutputFile=- ' fname]);
	try
		% The following while's are because we can't fseek a pipe and header size is variable.
		while (popenr(P,1,'char') ~= 10),	end		% Consume First  header line
		while (popenr(P,1,'char') ~= 10),	end		% Consume Second header line
		t = '                    ';
		for (k = 1:20)
			t(k) = popenr(P,1,'char');
			if (t(k) == 10),	break,		end		% EOF reached
		end
		t(k:end) = [];								% Remove excess characters
		ii = find(t == 32);							% Find the dimensions separator (a blank)
		width  = str2double(char( t(1 : ii(1)-1) )');
		height = str2double(char( t(ii(1)+1 : end) )');

		while (popenr(P,1,'char') ~= 10),	end		% Consume Fourth header line
		img = popenr(P, [height width 3], 'char');
	catch
		popenr(P, -1)					% Close the stdin stream
		errordlg(lasterr, 'Error'),		return
	end
	popenr(P, -1)						% Close the stdin stream
	mirone(img);

% --------------------------------------------------------------------------------------------------
function load_MB(handles, fname, frmt)
% Send MB datalist or single MB file to mbimport and get back an Image like that produced by mbswath
	cmd = '';
	fid = fopen(fname,'rt');
	todos = fread(fid,'*char');		fclose(fid);
	if (todos(1) == '#' || todos(1) == '>')
		ind = find(todos == 10);		% Find new lines
		cmd = todos(2:ind(1)-1)';
		%todos = todos(ind(1)+1:end);
	end
	%[nomes, frmt] = strread(todos,'%s %s');

	I = gmtmex(['mbimport -I' fname ' ' cmd]);
	mirone(I)
	if (handles.no_file),	delete(handles.figure1),	end		% Delete the old and empty fig.

% --------------------------------------------------------------------------------------------------
function varargout = ImageCrop_CB(handles, opt, opt2, opt3)
% OPT  is either a handle to a line that may be a rectangle/polygon,
%      a Mx2 or 2xM matrix with the (x y) vector coordinates, OR, if empty
%      calls rubberbandbox to do a interactive croping (called by "Crop Grid")
% OPT2 is a string to direct this function to different operations that
%      apply to the grid and update the image.
% OPT3 contains the interpolation method when OPT2 == 'FillGaps' ('cubic', 'linear' or 'sea')
%      Or the fill value when OPT2 == 'SetConst'
% Note: I won't make the "Drape" option active in the cropped window
%
% VARARGOUT -> If used will hold the result of this function instead of creating a new Fig
%              Currently implemented in cases:
%                  Crop image (opt == hLine), 'CropaWithCoords', 'CropaGrid_pure', 'ROI_Mean', 'ROI_Median', 'ROI_STD'

if (handles.no_file),		return,		end
first_nans = 0;		pal = [];		mask = [];	crop_pol = false;	% Defaults to croping from a rectangle
wasROI = false;		done = false;	invert = false;
if (nargin < 3),	opt2 = [];		end
if (nargin < 4),	opt3 = [];		end

if ~isempty(opt)				% OPT must be a rectangle/polygon handle (the rect may serve many purposes)
	if ((numel(opt) == 1) && ishandle(opt))
		x = get(opt,'XData');	y = get(opt,'YData');
	elseif (numel(opt) == 4)	% Assume that we have a [xmin xmax ymin ymax] BB
		x = [opt(1) opt(1) opt(2) opt(2) opt(1)];
		y = [opt(3) opt(4) opt(4) opt(3) opt(3)];
	else
		if (size(opt,2) > 2),	x = opt(1,1:end);	y = opt(2,1:end);	% Row vectors
		else					x = opt(:,1)';		y = opt(:,2)';		% Were col vectors, make them row for consistency
		end
	end

	xp(1) = min(x);		xp(2) = max(x);
	yp(1) = min(y);		yp(2) = max(y);
	if (numel(x) == 5 && (x(1) == x(end)) && (y(1) == y(end)) && ...		% Test if we have a rectangle
			(x(1) == x(2)) && (x(3) == x(4)) && (y(1) == y(4)) && (y(2) == y(3)) )
		if (strcmp(opt2,'SetConst'))
			[xp, yp] = aux_funs('adjust_rect', handles, xp, yp);	% Adjust such that only inside nodes are selected
		end
	else
		if (xp(1) < handles.head(1) || xp(2) > handles.head(2) || yp(1) < handles.head(3) || yp(2) > handles.head(4))
			% Somewhat rare case where the polygon extends to outside grid/img limits. Must crop it to them.
			P1.x = x(:);	P1.y = y(:);	P1.hole = 0;	P2.hole = 0;
			P2.x = [handles.head(1); handles.head(1); handles.head(2); handles.head(2); handles.head(1)];
			P2.y = [handles.head(3); handles.head(4); handles.head(4); handles.head(3); handles.head(3)];
			outPolyg = PolygonClip(P1, P2);				% Intersection of polygon and map limits
			x = outPolyg.x;		y = outPolyg.y;
			xp(1) = min(x);		xp(2) = max(x);
			yp(1) = min(y);		yp(2) = max(y);
		end
		x_lim = [min(x) max(x)];		y_lim = [min(y) max(y)];
		crop_pol = true;			% Flag that we are croping from a polygon
	end
	rect_crop = [xp(1) yp(1) (xp(2) - xp(1)) (yp(2) - yp(1))];	% Even when from a rectangle (to insure a known order)

	if isempty(opt2)				% Just pure Image croping
		Z_rect = get(handles.hImg,'CData');
		if (isa(Z_rect, 'uint16')),		Z_rect = scaleto8(Z_rect);	end		% In the rare event of "Noe Diluge" images
		limits = getappdata(handles.axes1,'ThisImageLims');
		I = cropimg(limits(1:2),limits(3:4),Z_rect,rect_crop);
		if (crop_pol)				% Shape cropping
			mask = ~(img_fun('roipoly_j',x_lim,y_lim,I,x,y));
			if (ndims(I) == 2),		I(mask) = 0;
			else
				for (k = 1:3)		% Sorry for the clutereness
					tmp = I(:,:,k);		tmp(mask) = 255;	I(:,:,k) = tmp;
				end
				mask = uint8(~mask);	cvlib_mex('CvtScale',mask,255,0);	% NEW. 1st transparency attempt
				I(:,:,4) = mask;
				clear tmp mask
			end
		end
		[m,n] = size(I);

	elseif (strcmp(opt2,'CropaWithCoords'))		% Crop Image with coordinates
		Z_rect = get(handles.hImg,'CData');
		if (isa(Z_rect, 'uint16')),		Z_rect = scaleto8(Z_rect);	end		% In the rare event of "Noe Diluge" images
		[I,r_c] = cropimg(handles.head(1:2),handles.head(3:4),Z_rect,rect_crop,'out_grid');
		if (crop_pol)				% Shape cropping
			mask = ~(img_fun('roipoly_j',x_lim,y_lim,I,x,y));
			if (ndims(I) == 2),		I(mask) = 0;
			else
				for (k = 1:3)
					tmp = I(:,:,k);		tmp(mask) = 255;	I(:,:,k) = tmp;
				end
				mask = uint8(~mask);	cvlib_mex('CvtScale',mask,255,0);	% NEW. 1st transparency attempt
				I(:,:,4) = mask;
				clear tmp mask
			end
		end
		[m,n] = size(I);

	elseif (strcmp(opt2,'ImplantGrid'))				% Read and external grid and implant it in host grid
		[Z_rect, r_c] = transplants(opt, 'grid', true, handles);
		if (isempty(Z_rect)),		return,		end		% User gave up

	else					% Extract the sub-grid inside the rectangle/polygon
		[X,Y,Z,head] = load_grd(handles);
		if isempty(Z),		return,		end		% An error message was already issued
		[Z_rect,r_c] = cropimg(head(1:2),head(3:4),Z,rect_crop,'out_grid');
		if (crop_pol)
			zzz = grdutils(Z_rect,'-L');	z_min = zzz(1);		clear zzz;
			resp = [];
			if (strcmp(opt2,'CropaGrid_pure'))
				resp = inputdlg({'Enter outside polygon value'},'Choose out value',[1 30],{sprintf('%.4f',z_min)});	pause(0.01)
				if isempty(resp),		return,		end
				resp = str2double(resp{1});
			elseif (strcmp(opt2,'ROI_SetConst'))	% Set the polygon in-or-out to cte
				resp = question({'Enter new grid value'},'Replace with cte value',[1 30],'NaN','whatever');
				if isempty(resp),		return,		end
				invert = resp{2};
				resp = str2double(resp{1});
			end
			if (isnan(resp)),	handles.have_nans = 1;	end

			mask = img_fun('roipoly_j', x_lim, y_lim, Z_rect, x, y);
			if (strcmp(opt2,'CropaGrid_pure'))
				Z_rect(~mask) = single(resp);
			elseif (strcmp(opt2,'ROI_SetConst'))
				if (invert),	mask = ~mask;	end			% Mask in the outside instead
				Z_rect(mask) = single(resp);				% Set the mask values to const
				if (invert)
					Z(:) = single(resp);					% Actually, we need to mask the whole outside
				else
					handles.Z_back = Z(r_c(1):r_c(2),r_c(3):r_c(4));	handles.r_c = r_c;		% For the Undo op
				end
				Z(r_c(1):r_c(2),r_c(3):r_c(4)) = Z_rect;
				if (isnan(resp)),		handles.have_nans = 1;	first_nans = 1;		end
			elseif (strcmp(opt2,'ROI_Mean'))
				res = Z_rect(mask);		res(isnan(res)) = [];	res = mean(double(res));
				if (nargout),	varargout{1} = res;		return,		end
				msgbox(sprintf('Mean over ROI = %15g',res),'ROI-Mean'),	return	
			elseif (strcmp(opt2,'ROI_Median'))
				res = Z_rect(mask);		res(isnan(res)) = [];	res = median(res);
				if (nargout),	varargout{1} = res;		return,		end
				msgbox(sprintf('Median over ROI = %15g',res),'ROI-Median'),	return	
			elseif (strcmp(opt2,'ROI_STD'))
				res = Z_rect(mask);		res(isnan(res)) = [];	res = std(double(res));
				if (nargout),	varargout{1} = res;		return,		end
				msgbox(sprintf('STD over ROI = %15g',res),'ROI-STD'),	return	
			elseif (strcmp(opt2,'ROI_MedianFilter'))
				[Z,Z_rect,handles] = roi_filtering(handles, Z, head, Z_rect, r_c, mask);
			elseif (strcmp(opt2,'ROI_SplineSmooth'))
				opt2 = 'SplineSmooth';	% Strip the 'ROI_' part so that we can use the same code as for rectangles
				wasROI = true;			% Signal the SplineSmooth code below that we need to mask result
			elseif (strcmp(opt2,'CropaGrid_histo'))
				Z_rect(~mask) = single(NaN);
			elseif (strcmp(opt2,'ROI_Clip'))
				opt2 = 'Clip';			% Now that we have the mask, make this case == to rectangle clip
			else
				warndlg('Unknown case in ImageCrop','Warning'),		return
			end
		end
		[m,n] = size(Z_rect);
	end
else					% Interactive croping (either Grid or Image)
	if (strcmp(opt2,'CropaGrid'))	% Arrive here when called by "Grid Tools -> Crop Grid"
		[X,Y,Z,head] = load_grd(handles);
		if isempty(Z),		return,		end
		[p1,p2] = rubberbandbox;
		x0 = min(p1(1),p2(1));		y0 = min(p1(2),p2(2));
		dx = abs(p2(1)-p1(1));		dy = abs(p2(2)-p1(2));
		[Z_rect,r_c] = cropimg([head(1) head(2)],[head(3) head(4)],Z,[x0 y0 dx dy],'out_grid');
		X = linspace( head(1) + (r_c(3)-1)*head(8), head(1) + (r_c(4)-1)*head(8), r_c(4) - r_c(3) + 1 );
		Y = linspace( head(3) + (r_c(1)-1)*head(9), head(3) + (r_c(2)-1)*head(9), r_c(2) - r_c(1) + 1 );
		head(1) = X(1);		head(2) = X(end);		head(3) = Y(1);		head(4) = Y(end);
		tit = 'Grid cuted by Mirone';		% Have to change this to reflect the old title
		GRDdisplay(handles,X,Y,Z_rect,head,tit,'Cropped_grid')
		return
	else			% Just a image crop op
		I = cropimg;	[m,n] = size(I);
	end
end

if (isempty(opt2) || strcmp(opt2,'CropaWithCoords'))	% Just pure Image croping
	if (m < 2 || n < 2),	return,	end		% Image too small.
	if (strcmp(get(handles.axes1,'Ydir'),'normal')),	I = flipdim(I,1);	end
	if (ndims(I) == 2)
		pal = get(handles.figure1, 'Colormap');
		if (length(pal) == 64), pal = jet(256);		end		% Risky - This is a patch for "Find Clusters"
		setappdata(0,'CropedColormap',pal);			% indexed image, so I need to save it's colormap
	end
	set(handles.figure1,'pointer','arrow');
	if (~isempty(opt2))
		head(2) = handles.head(1) + (r_c(4)-1)*handles.head(8);
		head(1) = handles.head(1) + (r_c(3)-1)*handles.head(8);
		head(4) = handles.head(3) + (r_c(2)-1)*handles.head(9);
		head(3) = handles.head(3) + (r_c(1)-1)*handles.head(9);
		head(5:9) = [0 255 0 handles.head(8:9)];	tmp.name = 'Cropped_image';
		tmp.head = head;		tmp.geog = handles.geog;			tmp.X = head(1:2);		tmp.Y = head(3:4);
		if (~isempty(pal)),		tmp.cmap = pal;		end
		proj4 = getappdata(handles.figure1,'Proj4');
		if (~isempty(proj4))
			tmp.srsWKT = ogrproj(proj4);
		else
			projWKT = getappdata(handles.figure1,'ProjWKT');
			if (~isempty(projWKT)),		tmp.srsWKT = projWKT;	end
		end
	end
	if (nargout)
		varargout{1} = I;
		if (nargout == 2),		varargout{2} = tmp;	end
	elseif (isempty(opt2))				% Crop without coords
		mirone(I);
	else
		mirone(flipdim(I,1),tmp);		% Crop with coords
	end
	done = true;				% We are done. BYE BYE.

elseif (strncmp(opt2(1:min(length(opt2),9)),'CropaGrid',9))		% Do the operation indicated in opt2(11:end) & return
	curr_opt = opt2(11:end);
	if (~strcmp(curr_opt,'pure'))			% We will need those for all other options
		head(2) = head(1) + (r_c(4)-1)*head(8);			head(1) = head(1) + (r_c(3)-1)*head(8);
		head(4) = head(3) + (r_c(2)-1)*head(9);			head(3) = head(3) + (r_c(1)-1)*head(9);
		if (isa(Z,'single')),	zz = grdutils(Z_rect,'-L');			head(5:6) = [zz(1) zz(2)];
		else					head(5) = double(min(Z_rect(:)));	head(6) = double(max(Z_rect(:)));
		end
		to_func.Z = Z_rect;		to_func.head = head;
	end
	if (strcmp(curr_opt,'pure'))			% PURE means pure CropaGrid
		X = linspace( head(1) + (r_c(3)-1)*head(8), head(1) + (r_c(4)-1)*head(8), r_c(4) - r_c(3) + 1 );
		Y = linspace( head(3) + (r_c(1)-1)*head(9), head(3) + (r_c(2)-1)*head(9), r_c(2) - r_c(1) + 1 );
		head(1) = X(1);		head(2) = X(end);		head(3) = Y(1);		head(4) = Y(end);
		srsWKT = [];
		if (~handles.geog)
			prjInfoStruc = aux_funs('getFigProjInfo',handles);
			if (~isempty(prjInfoStruc.projWKT))	% TODO. Otherwise check if prjInfoStruc.proj4 and convert it to WKT
				srsWKT = prjInfoStruc.projWKT;
			end
		end
		if (~nargout)						% Create a new Fig
			tit = 'Grid cuted by Mirone';	% Have to change this to reflect the old title
			GRDdisplay(handles,X,Y,Z_rect,head,tit,'Cropped_grid',srsWKT)
		else								% Send back the cropped grid to whom asked for it. 
			varargout = {X,Y,Z_rect,head};	% Is not going to be easy to document this
		end
	elseif (strcmp(curr_opt,'histo'))		% HISTO means compute histogram inside the selected rect area
		GridToolsHistogram_CB(guidata(handles.figure1), to_func);
	elseif (strcmp(curr_opt,'power'))		% POWER means compute log10 power spectrum
		GridToolsSectrum_CB(guidata(handles.figure1), 'Power', to_func)
	elseif (strcmp(curr_opt,'autocorr'))	% AUTOCORR means compute the autocorrelation
		GridToolsSectrum_CB(guidata(handles.figure1), 'Autocorr', to_func)
	elseif (strcmp(curr_opt,'fftTools'))	% FFTTOOLS means call the fft_stuff
		GridToolsSectrum_CB(guidata(handles.figure1), 'Allopts', to_func)
	end
	done = true;				% We are done. BYE BYE.

elseif (strcmp(opt2,'FillGaps'))
	if ~any(isnan(Z_rect(:)))	% No gaps
		warndlg('Selected area does not have any voids (NaNs)','Warning'),	 return
	end

	X = linspace( head(1) + (r_c(3)-1)*head(8), head(1) + (r_c(4)-1)*head(8), r_c(4) - r_c(3) + 1 );
	Y = linspace( head(3) + (r_c(1)-1)*head(9), head(3) + (r_c(2)-1)*head(9), r_c(2) - r_c(1) + 1 );
	if (~isempty(opt3) && strcmp(opt3,'surface'))
		opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', X(1), X(end), Y(1), Y(end));
		opt_I = sprintf('-I%.10f/%.10f',head(8),head(9));
	end
	Z_rect = double(Z_rect);	% It has to be
	aa = isnan(Z_rect(:));
	[X,Y] = meshgrid(X,Y);
	ZZ = Z_rect(:);		ZZ(aa) = [];
	XX = X(:);			XX(aa) = [];
	YY = Y(:);			YY(aa) = [];
	if (~isempty(opt3))
		switch opt3
			case 'surface', Z_rect = gmtmbgrid_m(XX,YY,ZZ,opt_R,opt_I,'-T.25', '-Mz');
			case 'cubic',	Z_rect = griddata_j(XX,YY,ZZ,X,Y,'cubic');
			case 'linear',	Z_rect = griddata_j(XX,YY,ZZ,X,Y,'linear');
			case 'sea',		Z_rect(aa) = 0;
		end
	else
		Z_rect = gmtmbgrid_m(XX,YY,ZZ,opt_R,opt_I,'-T.25','-v', '-Mz');
	end
	clear X XX Y YY ZZ;

elseif (strcmp(opt2,'SplineSmooth'))
	X = linspace( head(1) + (r_c(3)-1)*head(8), head(1) + (r_c(4)-1)*head(8), r_c(4) - r_c(3) + 1 );
	Y = linspace( head(3) + (r_c(1)-1)*head(9), head(3) + (r_c(2)-1)*head(9), r_c(2) - r_c(1) + 1 );
	Z_rect = double(Z_rect);	% It has to be
	[pp, p_guess] = spl_fun('csaps',{Y(1:min(m,10)),X(1:min(n,10))},Z_rect(1:min(m,10),1:min(n,10)));% Get a good estimate of p
	prompt = {'Enter smoothing p paramer'};		dlg_title = 'Smoothing parameter input';
	defAns = {sprintf('%.12f',p_guess{1})};		resp = inputdlg(prompt,dlg_title,[1 38],defAns);
	if (isempty(resp)),		return,		end
	resp = str2double(resp{1});
	if (isnan(resp)),		set(handles.figure1,'pointer','arrow'),		return,		end
	pp = spl_fun('csaps',{Y,X},Z_rect,resp);
	Z_rect = spl_fun('fnval',pp,{Y,X});		clear pp;
	handles.Z_back = Z(r_c(1):r_c(2),r_c(3):r_c(4));	handles.r_c = r_c;			% For the Undo op
	if (wasROI)				% Apply the mask and smooth over the mask edges
		Z_rect = smooth_roipoly_edge(head, handles.have_nans, Z, handles.Z_back, Z_rect, r_c, mask, 3);
	end

elseif (strcmp(opt2,'MedianFilter'))
	[Z,Z_rect,handles] = roi_filtering(handles, Z, head, Z_rect, r_c, 'rect', 'no');

elseif (strcmp(opt2,'Clip'))
	handles.Z_back = Z(r_c(1):r_c(2),r_c(3):r_c(4));	handles.r_c = r_c;			% For the Undo op
	[Z_rect,head] = ml_clip(handles, handles.Z_back);
	if (isempty(Z_rect)),		return,		end
	handles.head(5:6) = head(5:6);
	if (crop_pol)			% Means, if ROI_Clip
		Z_rect(~mask) = handles.Z_back(~mask);
		zz = grdutils(Z_rect,'-L');		handles.head(5:6) = [min(zz(1),handles.head(5)) max(zz(2),handles.head(6))];
	end

elseif (strcmp(opt2,'SetConst'))		% Replace grid values inside rect by a cte value
	if (isempty(opt3))
		resp = inputdlg({'Enter new grid value'},'Replace with cte value',[1 30]);	pause(0.01)
		if (isempty(resp)),		return,		end
		resp = str2double(resp);
		if (~isreal(resp)),		resp = NaN;		end		% A 'i' or a 'j' in resp would have caused this
	else
		resp = opt3;
	end
	Z_rect = single(resp);
	handles.Z_back = Z(r_c(1):r_c(2),r_c(3):r_c(4));	handles.r_c = r_c;			% For the Undo op
	if (~handles.have_nans && isnan(resp))				% See if we have new NaNs
		handles.have_nans = 1;		first_nans = 1;
	end
	
elseif (strcmp(opt2,'GetMean'))			% ...
	Z_rect(isnan(Z_rect)) = [];		res = mean(double(Z_rect(:)));
	if (nargout),	varargout{1} = res;		return,		end
	msgbox(sprintf('Mean over ROI = %15g',res),'ROI-Mean')
	done = true;

elseif (strcmp(opt2,'GetMedian'))		% ...
	Z_rect(isnan(Z_rect)) = [];		res = median(Z_rect(:));
	if (nargout),	varargout{1} = res;		return,		end
	msgbox(sprintf('Median over ROI = %15g',res),'ROI-Median')
	done = true;
	
elseif (strcmp(opt2,'GetSTD'))			% ...
	Z_rect(isnan(Z_rect)) = [];		res = std(double(Z_rect(:)));
	if (nargout),	varargout{1} = res;		return,		end
	msgbox(sprintf('STD over ROI = %15g',res),'ROI-STD')
	done = true;

elseif (strcmp(opt2,'ImplantGrid'))		% The first part of this job was done above, where we got Z_rect and r_c
	[X,Y,Z,head] = load_grd(handles);

end

if (done),		return,		end

if (~strcmp(opt2,'MedianFilter') && ~strcmp(opt2,'ROI_SetConst'))		% Otherwise, this was already done in roi_filtering
	if (isa(Z,'single')),		Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
	elseif (isa(Z,'int16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = int16(Z_rect);
	elseif (isa(Z,'uint16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = uint16(Z_rect);
	else						Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
	end
	if (handles.have_nans && ~any(isnan(Z_rect(:))))	% Check that old NaNs had not been erased
		handles.have_nans = grdutils(Z,'-N');
	end
end

if ~isempty(opt2)		% Here we have to update the image in the processed region
	if (isa(Z,'single')),	zz = grdutils(Z,'-L');		z_min = zz(1);		z_max = zz(2);
	else					z_min = double(min(Z(:)));	z_max = double(max(Z(:)));
	end
	img = [];
	if ((abs(z_min - head(5)) > 1e-5 || abs(z_max - head(6)) > 1e-5) && handles.Illumin_type == 0)
		img = scaleto8(Z);		% Z_MIN or Z_MAX have changed. Need to recompute image (but only if no illumin)
	end
	if (first_nans)		% We have NaNs for the first time. Adjust the colormap
		aux_funs('colormap_bg',handles,Z,get(handles.figure1,'Colormap'));
	end
	z_int = uint8(round( ((double(Z_rect) - z_min) / (z_max - z_min))*255 ));
	if (handles.Illumin_type >= 1 && handles.Illumin_type <= 4)
		illumComm = getappdata(handles.figure1,'illumComm');
		z_int = ind2rgb8(z_int,get(handles.figure1,'Colormap'));	% z_int is now RGB
		if (handles.Illumin_type == 1)
			opt_N = sprintf('-Nt1/%.6f/%.6f',handles.grad_sigma, handles.grad_offset);
			if (handles.geog),	R = grdgradient_m(Z_rect,head,'-M',illumComm,opt_N);
			else				R = grdgradient_m(Z_rect,head,illumComm,opt_N);
			end
		else
			R = grdgradient_m(Z_rect,head,illumComm, '-a1');
		end
		z_int = shading_mat(z_int,R,'no_scale');	% and now it is illuminated
	elseif (handles.Illumin_type ~= 0)
		warndlg('Sorry, this operation is not allowed with this shading illumination type','Warning')
		return
	end
	if (isempty(img)),		img = get(handles.hImg,'CData');	end		% If img was not recomputed, get from screen 
	handles.img_back = img(r_c(1):r_c(2),r_c(3):r_c(4),1:end);			% For the undo op
	if (~isempty(mask) && handles.Illumin_type ~= 0)
		mask = repmat(~mask,[1 1 3]);
		z_int(mask) = handles.img_back(mask);
	end
	img(r_c(1):r_c(2),r_c(3):r_c(4),1:end) = z_int;			clear z_int Z_rect R;
	set(handles.hImg,'CData',img)

	head(5) = z_min;	head(6) = z_max;
	handles.computed_grid = 1;		handles.head = head;	%handles.origFig = img;
	setappdata(handles.figure1,'dem_z',Z);
end

% UNDO that works only with these cases
if (~invert && (numel(opt) == 1) && any(strcmp(opt2,{'MedianFilter' 'ROI_MedianFilter' 'SetConst' 'ROI_SetConst' 'SplineSmooth'})))
	cmenuHand = get(opt,'UIContextMenu');
	uimenu(cmenuHand, 'Label', 'Undo', 'Separator','on', 'Callback', {@do_undo,handles.figure1,opt,cmenuHand});
end
guidata(handles.figure1, handles);		set(handles.figure1,'pointer','arrow')

% -----------------------------------------------------------------------------------------
function do_undo(obj,event,hFig,h,img)
	handles = guidata(hFig);
	[X,Y,Z] = load_grd(handles);	% Experimental. No testing for error in loading
	Z(handles.r_c(1):handles.r_c(2),handles.r_c(3):handles.r_c(4)) = handles.Z_back;
	setappdata(handles.figure1,'dem_z',Z);
	handles.origFig = get(handles.hImg,'CData');
	handles.origFig(handles.r_c(1):handles.r_c(2),handles.r_c(3):handles.r_c(4),1:end) = handles.img_back;
	set(handles.hImg,'CData',handles.origFig)
	guidata(handles.figure1,handles)
	delete(findobj(get(h,'UIContextMenu'),'Label', 'Undo'))

% --------------------------------------------------------------------
function ImageResetOrigImg_CB(handles)
	if (handles.no_file || handles.image_type == 20),		return,		end
	set(handles.hImg,'CData', handles.origFig);
	set(handles.figure1,'ColorMap',handles.origCmap)
	handles.Illumin_type = 0;		handles.firstIllum = 1;
	handles.validGrid = handles.validGrid_orig;		handles.was_int16 = handles.was_int16_orig;
	handles.computed_grid = handles.computed_grid_orig;
	set(handles.ImgHist,'checked','off');
	aux_funs('togCheck',get(handles.ImRestore,'UserData'), [handles.ImMod8cor handles.ImMod8gray handles.ImModBW handles.ImModRGB])
	guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function ImageHistEqualize_CB(handles, hObject)
	if (handles.no_file),		return,		end
	zz = get(handles.hImg,'CData');
	if strcmp(get(hObject,'checked'),'off')		% Then equalize
		if (ndims(zz) == 3);
			zz = cvlib_mex('color',zz,'rgb2hls');
			zz(:,:,2) = img_fun('histeq_j',zz(:,:,2));
			%zz(:,:,2) = fcnBPDFHE(zz(:,:,2));
			J = cvlib_mex('color',zz,'hls2rgb');
		else
			J = img_fun('histeq_j',zz);
		end
		set(handles.hImg,'CData', J),	set(hObject,'checked','on')
	else								% Then de-equalize
		set(handles.hImg,'CData', handles.origFig);
		handles.Illumin_type = 0;		set(hObject,'checked','off');
	end
	guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function ImageHistEqualizeGrid_CB(handles, hObject)
	if (~handles.validGrid),	return,		end
	if strcmp(get(hObject,'checked'),'off')		% Then equalize
		[X,Y,Z,head] = load_grd(handles);
		if isempty(Z),		return,		end
		out = grdutils(Z,'-S');			% Get mean and std
		handles.cur_pal = get(handles.figure1, 'Colormap');
		new_pal = cdf2pal(head(5),head(6),out(1),out(2),handles.cur_pal);
		aux_funs('colormap_bg',handles,Z,new_pal);
		handles.hist_grid = 1;		set(hObject,'checked','on')
	else			% Then de-equalize
		set(handles.figure1,'Colormap',handles.cur_pal);	handles.hist_grid = 0;
		handles.Illumin_type = 0;		set(hObject,'checked','off');
	end
	guidata(hObject, handles);

% --------------------------------------------------------------------
function ImageSegment_CB(handles, hObject)
	if (handles.no_file),		return,		end
	rgbIm = get(handles.hImg,'CData');
	if (ndims(rgbIm) == 2),		rgbIm = ind2rgb8(rgbIm, get(handles.figure1, 'Colormap'));	end
	set(handles.figure1,'pointer','watch')

	luvIm = cvlib_mex('CvtScale',single(rgbIm), 1/255); 
	luvIm = cvlib_mex('color',luvIm,'rgb2luv');
	rgbIm = permute(rgbIm,[3 2 1]);
	luvIm = permute(luvIm,[3 2 1]);

	p.steps = 2;		p.SpeedUp = 2;		p.synergistic = true;
	p.SpatialBandWidth = 7;			p.RangeBandWidth = 6.5;
	p.MinimumRegionArea = 30;		p.GradientWindowRadius = 2;
	p.MixtureParameter = .3;		p.EdgeStrengthThreshold = .3;
	[fimage] = edison_wrapper_mex(luvIm, rgbIm, p);

	fimage = permute(fimage, [3 2 1]);
	luvIm = cvlib_mex('color',fimage,'luv2rgb');	clear fimage
	rgbIm = scaleto8(luvIm, -8);
	set(handles.figure1,'pointer','arrow')
	
	if (handles.image_type == 2)
		h = mirone(rgbIm);		set(h,'Name','Color segmented')
	else
		tmp = struct('X',handles.head(1:2), 'Y',handles.head(3:4), 'name','Color segmented', 'geog',handles.geog, 'head',handles.head);
		strWKT = getappdata(handles.figure1,'ProjWKT');
		if (~isempty(strWKT)),	tmp.srsWKT = strWKT;	end
		mirone(rgbIm, tmp)
	end

% --------------------------------------------------------------------
function PanZoom_CB(handles, hObject, opt)
	if (handles.no_file),	set(hObject,'State','off'),		return,		end	
	if (strcmp(get(handles.Tesoura,'State'),'on')),		set(handles.Tesoura,'State','off'),		end		% If Scisors were on

	if (strcmp(opt,'zoom'))
		if ( strcmp(get(hObject,'State'),'on') )
			if (strcmp(get(handles.Mao,'State'),'on'))
				set(handles.Mao,'State','off'),		pan(handles.figure1, 'off');
			end
			zoom_j(handles.figure1, 'on');
		else
			zoom_j(handles.figure1, 'off');
		end
	else		% Pan case
		if ( strcmp(get(hObject,'State'),'on') )
			lims1 = getappdata(handles.axes1,'ThisImageLims');
			lims2 = [get(handles.axes1,'xlim') get(handles.axes1,'ylim')];
			if (isequal(lims1, lims2))
				set(hObject,'State','off')		% If we are in no zoom state there is nothing to pan
				return
			end
			if (strcmp(get(handles.Zoom,'State'),'on'))
				set(handles.Zoom,'State','off'),	zoom_j(handles.figure1, 'off');
			end
			pan(handles.figure1, 'on');
		else
			pan(handles.figure1, 'off');
		end
	end

% --------------------------------------------------------------------
function zoom_state(handles, state)
% Sets the zoom sate to off, or reset it to on if ...
	%if (handles.IAmOctave),		guidata(handles.figure1,handles),	return,		end
	switch state
		case 'off_yes'			% Set zoom permanently off
			if (~handles.no_file),	zoom_j(handles.figure1,'off');		end		% No need to check when first file in
			set(handles.Zoom,'State','off');
			handles.zoom_state = 0;		guidata(handles.figure1,handles)
		case 'maybe_off'		% If zoom was active, keep trace of it
			zoom_j(handles.figure1,'off');		h = handles.Zoom;
			if (strcmp(get(h,'State'),'on')),	handles.zoom_state = 1;
			else								handles.zoom_state = 0;
			end
			set(h,'State','off'),	guidata(handles.figure1,handles)
		case 'maybe_on'			% Check if zoom has to be re-activated
			handles = guidata(handles.figure1);		% Need to get the updated handles
			if (handles.zoom_state)
				zoom_j(handles.figure1,'on');	set(handles.Zoom,'State','on');
			end
	end

% --------------------------------------------------------------------
function hand = FileNewBgFrame_CB(handles, region, imSize, figTitle)
% Create a empty window with a frame selected in bg_region
% However, if REGION was transmited, it is assumed to have [x_min x_max y_min y_max is_geog toDef]
% TODEF	Logical that if true is used to add a new uimenu to call the write def symbol code
% IMSIZE may either be the image size or the Figure title
% HAND	Optional, updated version of HANDLES
	if (nargout),	hand = handles;		end		% In case of a forced 'return' bellow do not error

	if (nargin == 1)
		region = bg_region;		% region contains [x_min x_max y_min y_max is_geog toDef]
		if isempty(region),		return,		end		% User gave up
		if (region(5)),			region(5) = aux_funs('guessGeog',region(1:4));		end		% Refine (can be 2)
		handles.geog = region(5) + 10;			% The +10 instructs show_image to accept this val(and subtracts 10)
	end
	if (nargin <= 2),		imSize = [];		figTitle = 'Mirone Base Map';		end
	if (nargin == 3 && isa(imSize,'char')),		figTitle = imSize;	imSize = [];	end
	if (numel(region) == 6 && region(6))		% Add a new uimenu to call the write def symbol code
		aux_funs('addUI', handles)
		handles.is_defRegion = true;
	end

	if ( any(isnan(region(1:4))) )
		errordlg('The requested region limts is undeterminated (it has NaNs)','Error'),		return
	end
	X = region(1:2);	Y = region(3:4);		handles.head = [X Y 0 255 0];
	if (isempty(imSize) || numel(imSize) ~= 2)
		scrsz = get(0,'ScreenSize');		% Get screen size
		aspect = diff(Y) / diff(X);
		nx = round(scrsz(3)*.75);	ny = round(nx * aspect);
		if (ny > scrsz(4) - 30)
			ny = scrsz(4) - 30;		nx = round(ny / aspect);
		end
	else
		nx = imSize(1);			ny = imSize(2);
	end
	handles.head(8) = diff(X) / (nx - 1);	handles.head(9) = diff(Y) / (ny - 1);
	Z = repmat(uint8(255),ny,nx);
	pal = repmat(handles.bg_color,256,1);	set(handles.figure1,'Colormap',pal);
	handles.image_type = 20;
	handles = show_image(handles,figTitle,X,Y,Z,0,'xy',0,imSize);
	if (handles.is_defRegion)			% def region doesn't use these
		set([handles.Geography handles.Plates handles.MagGrav ...
			handles.Seismology handles.Projections], 'Vis', 'off')
	end
	drawnow			% Otherwise, the damn Java makes a black window until all posterior elements are plotted
	aux_funs('isProj',handles);				% Check about coordinates type
	if (nargout),	hand = handles;		end

% --------------------------------------------------------------------
function FileSaveGMTgrid_CB(handles, opt)
% Save internaly computed grids and GDAL recognized DEM grids into GMT grd grids
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	if (nargin == 1),	opt = [];	end

	[X,Y,Z,head] = load_grd(handles);
	if isempty(Z),	return,		end		% An error message was already issued
	if (~isempty(opt) && strcmp(opt,'Surfer'))
		txt1 = 'Surfer 6 binary grid (*.grd,*.GRD)';	txt2 = 'Select output Surfer 6 grid';
	else			% Internaly computed grid
		tit = 'Grid computed with Mirone';
		txt1 = 'netCDF grid format (*.grd,*.GRD)';		txt2 = 'Select output GMT grid';
	end

	[FileName,PathName] = put_or_get_file(handles,{'*.grd;*.GRD',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.grd');
	if isequal(FileName,0),		return,		end
	f_name = [PathName FileName];

	if (~isempty(opt) && strcmp(opt,'Surfer'))
		fid = fopen(f_name, 'wb');		fwrite(fid,'DSBB','char');
		fwrite(fid,[size(Z,2) size(Z,1)],'int16');		fwrite(fid,head(1:6),'double');
		if (~isa(Z,'single')),			Z = single(Z);	end
		if (handles.have_nans),			Z(isnan(Z)) = 1.701410e38;		end
		fwrite(fid,Z','float32');		fclose(fid);
		return
	end

	if (isappdata(handles.axes1,'DatumProjInfo'))		% If it was a grid imported by gdal, uppdate the title
		DPI = getappdata(handles.axes1,'DatumProjInfo');
		tit = ['Projection: ' DPI.projection ' Datum: ' DPI.datum];
		if (length(tit) > 80),		tit = tit(1:80);	end		% (1:80) otherwise it BOOMs
	end
	% Defaults and srsWKT fishing are set in nc_io
	misc = struct('x_units',[],'y_units',[],'z_units',[],'z_name',[],'desc',[], ...
		'title',tit,'history',[],'srsWKT',[], 'strPROJ4',[]);
	nc_io(f_name, 'w', handles, Z, misc)

% --------------------------------------------------------------------
function File_img2GMT_RGBgrids_CB(handles, opt1, opt2)
% Save image as a triplet of gmt grids - R,G,B.
% OPT1 == 'image' || == [] Capture only the image and not graphical elements
%	(lines, symbols, etc...). The grids nrow & ncol is the same as the image
%	number of lines and pixels.
% OPT1 == 'screen' Does a screen capture that includes all the graphical elements
%	that may have been drawn (lines, symbols, etc...).
% OPT2 == fname. It is used by the write gmt script routine to capture the image and
%	write the image as a triplet of gmt grids with name stem = OPT2.
	if (handles.no_file),	return,		end

	if    (nargin == 1),	opt1 = 'image';		opt2 = [];
	elseif(nargin == 2),	opt2 = [];			end

	str1 = {'*.grd;*.GRD','netCDF int2 grid format (*.grd,*.GRD)'; '*.*', 'All Files (*.*)'};
	if (isempty(opt2))
		[FileName,PathName] = put_or_get_file(handles,str1,'Select output GMT grid','put');
		if isequal(FileName,0),		return,		end
	else		% It means the output file name was transmited in input
		[PathName,FileName] = fileparts(opt2);		PathName = [PathName filesep];
	end

	[PATH,FNAME,EXT] = fileparts([PathName FileName]);
	if isempty(EXT),	EXT = '.grd';		end	
	f_name_r = [PathName FNAME '_r' EXT];	f_name_g = [PathName FNAME '_g' EXT];
	f_name_b = [PathName FNAME '_b' EXT];

	if (strcmp(opt1,'image')),			img = get(handles.hImg,'CData');			% Get image
	elseif (strcmp(opt1,'screen')),		img = snapshot(handles.figure1,'noname');	% Screen capture with resizing option
	else								img = flipdim(imcapture(handles.axes1,'img',0),1);		% Call from write_script
	end

	if (ndims(img) == 2)
		img = ind2rgb8(img,get(handles.figure1,'Colormap'));	% Need this because image is indexed
	end
	flip = false;
	if (~strcmp(get(handles.axes1,'Ydir'),'normal')),	flip = true;	end

	% Defaults and srsWKT fishing are set in nc_io
	if (~flip)
		nc_io(f_name_r, 'w', handles, img(:,:,1));		nc_io(f_name_g, 'w', handles, img(:,:,2))
		nc_io(f_name_b, 'w', handles, img(:,:,3))
	else
		nc_io(f_name_r, 'w', handles, flipud(img(:,:,1)));		nc_io(f_name_g, 'w', handles, flipud(img(:,:,2)))
		nc_io(f_name_b, 'w', handles, flipud(img(:,:,3)))
	end
	set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function FileSaveENCOMgrid_CB(handles)
% Save memory resident grids into the Encom grid format
	if (aux_funs('msg_dlg',14,handles)),	return,		end

	txt1 = 'Encom grid format (*.grd,*.GRD)';	txt2 = 'Select output Encom grid';
	[FileName,PathName] = put_or_get_file(handles,{'*.grd;*.GRD',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.grd');
	if isequal(FileName,0),		return,		end

	[X,Y,Z,head,m,n] = load_grd(handles);
	if isempty(Z),		return,		end		% An error message was already issued
	f_name = [PathName FileName];

	set(handles.figure1,'pointer','watch')
	fid = fopen(f_name,'wb');
	ID = ['Model Vision Grid' repmat(' ',1,80-17)];		fwrite(fid,ID,'80*char');
	ID = ['Mirone_Grid' repmat(' ',1,80-11)];			fwrite(fid,ID,'80*char');
	fwrite(fid,'  NO_REF','char');
	fwrite(fid,'GRIDFPT ZNIL','char');
	fwrite(fid,single(-2e16),'float32');
	fwrite(fid,'ROWS','char');		fwrite(fid,m,'float32');
	fwrite(fid,'COLS','char');		fwrite(fid,n,'float32');
	fwrite(fid,'XORG','char');		fwrite(fid,single(head(1)),'float32');
	fwrite(fid,'YORG','char');		fwrite(fid,single(head(3)),'float32');
	fwrite(fid,'DX  ','char');		fwrite(fid,single(head(8)),'float32');
	fwrite(fid,'DY  ','char');		fwrite(fid,single(head(9)),'float32');
	fwrite(fid,'DEGR','char');		fwrite(fid,single(0),'float32');
	if (handles.have_nans),			Z(isnan(Z)) = -2e16;	end
	Z = (rot90(Z,1));				Z = flipud(Z);		% I cannot do better than this manip
	fwrite(fid,Z,'float32');		fclose(fid);
	set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function out = ExtractProfile_CB(handles, opt, opt2)
% OPT == 'point', 'dynamic', '3D' or 'input'. POINT and 3D, interpolate at the line vertex only
%	3D interpolates along the third dimension (but is not fully 3D interp)
%	'input' means that OPT2 MUST exist and holds a Mx2 matrix with the x,y coordinates.
%	        This option is to be used in a non-interactive mode and normally goes along with 'OUT'
% OUT -> When used it returns a structure with 'x', y', and 'z' members (with obvious meaning)
%        and does not display the result in the 2D profile tool.

	if (nargout),	out = [];	end
	if (handles.no_file),	return,		end
	if (nargin == 1),	opt = '';		end
	point_int = strcmp(opt,'point');
	do_3D = strcmp(opt,'3D');
	track_is_done = false;		do_stack = false;

	[X,Y,Z] = load_grd(handles,'silent');
	if (isempty(Z) && handles.validGrid)					% Grid not in memory error
		warndlg('Grid is not in memory (too big). Bye','Warning'),	return
	elseif (~handles.validGrid && strcmp(opt,'dynamic') && (ndims(get(handles.hImg,'CData')) == 3) )
		errordlg('Extracting a dynamic profile of a RGB image is not suported.','ERROR')
		return
	end

	zoom_state(handles,'maybe_off')
	if ~isempty(getappdata(handles.figure1,'TrackThisLine'))		% Most common use
		hand = getappdata(handles.figure1,'TrackThisLine');
		xp = get(hand,'Xdata');			yp = get(hand,'Ydata');
		rmappdata(handles.figure1,'TrackThisLine')		% Clear it so that the next time it may work when called interactivelly
	elseif ~isempty(getappdata(handles.figure1,'StackTrack'))		% Stacking interpolation
		do_stack = getappdata(handles.figure1,'StackTrack');		% do_stack is in fact the line handle (named for convenience)
		xp = get(do_stack,'Xdata');		yp = get(do_stack,'Ydata');
		rmappdata(handles.figure1,'StackTrack')
	else
		if (strcmp(opt,'dynamic'))
			getline_j(handles.figure1,'dynamic');					% Calls grid_profiler at every ButtonMotion
			hDynProfAx = getappdata(handles.axes1,'dynProfile');	% Get the handle of the dynamic profile Axes
			hLine = getappdata(hDynProfAx,'theLine');
			ud = get(hLine, 'UserData');
			xx = ud(:,1);			yy = ud(:,2);		zz = get(hLine,'YData');
			delete(hDynProfAx);		rmappdata(handles.axes1,'dynProfile');	% We are over with the dynamic tracking. Remove that axes
		else
			if (strcmp(opt,'input'))
				xp = opt2(:,1);			yp = opt2(:,2);				% We are not testing here that nargin == 3 and opt2 is correct
				point_int = true;
			else
				[xp,yp] = getline_j(handles.figure1);
			end
			if (numel(xp) < 2),		zoom_state(handles,'maybe_on'),		return,		end
			[xx, yy, zz] = grid_profiler(handles.figure1, xp, yp, point_int, false, false, do_3D);	% The 'false' is for not 'do_dynamic'
		end
		track_is_done = true;
	end

	if (~track_is_done)			% Otherwise we already know them
		[xx, yy, zz] = grid_profiler(handles.figure1, xp, yp, point_int, track_is_done, do_stack, do_3D);
	end

	if (nargout)				% Than we are over here.
		out = struct('x',xx, 'y',yy, 'z',zz);
		return
	end

	zoom_state(handles,'maybe_on')
	if (do_3D && numel(xx) > 1)				% Save result on file (because it has more than one profile)
		embrulho.x = xx;		embrulho.y = yy;
		embrulho.z = zz;		embrulho.time_z = handles.time_z;
		draw_funs([],'save_xyz',embrulho)
	elseif (point_int && ~do_3D)			% Save result on file
		draw_funs([],'save_xyz',[xx(:) yy(:) zz(:)])
	else									% Disply profile in ecran
		[pato,name,ext] = fileparts(get(handles.figure1,'Name'));
		if (~isempty(ext))
			ext = strtok(ext);				% Remove the "@ ??%" part
		else
			ind = strfind(name, ' @');
			if (~isempty(ind))
				name = name(1:ind(1)-1);	%		"
			end
		end
		if (do_3D)
			ecran(handles,handles.time_z,zz,['3D Profile from ' name ext])
		elseif (~isa(zz,'cell'))
			ecran(handles,xx,yy,zz,['Track from ' name ext])
		else
			% Make 3 call, one for each channel. This is a bit ugly but plotting the 3 curves 
			% on same 'ecran' window would probably open the door to several problems. We'll see.
			cor = {'r' 'g' 'b'};
			for (k = 1:3)
				hFig = ecran(handles,xx,yy,zz{k},['Track from ' name ext ' (' upper(cor{k}) ')']);
				h = findobj(hFig,'type','line');		set(h,'color',cor{k});		% Set color by band
				if (k > 1)				% Cascade the 3 figures so that the user get immediately aware
					pos = get(hFig, 'Pos');		pos(2) = pos(2) - (k-1)*30;		set(hFig, 'Pos', pos);
				end
			end
		end
	end

% --------------------------------------------------------------------
function FileOpen_ENVI_Erdas_CB(handles, opt, opt2)
% This function reads both ENVI or Erdas files. Furthermore, based on the file byte
% type it guesses if we are dealing with a typical grid file (in which case it is
% treated like a native gmt grid) or a raster image file.

	if (nargin == 2)		% Otherwise, OPT2 already contains the File name
		str1 = {'*.img;*.IMG', [opt ' (*.img,*.IMG)']; '*.*', 'All Files (*.*)'};
		[FileName,PathName] = put_or_get_file(handles,str1,['Select ' opt ' File'],'get');
		if isequal(FileName,0),		return,		end
	else
		PathName = [];		FileName = opt2;
	end

	handles.fileName = [PathName FileName];

	att = gdalread(handles.fileName,'-M','-C');
	if (~strcmp(att.Band(1).DataType,'Byte'))			% We have a "grid"
		loadGRID(handles, handles.fileName, 'guess', att);		return
	end
	
	if (~isempty(att.GeoTransform)),	handles.image_type = 3;		opt_U = '-U';
	else								handles.image_type = 2;		opt_U = ' ';
	end

	set(handles.figure1,'pointer','watch')
	zz = gdalread(handles.fileName, opt_U);
	handles.head = att.GMT_hdr;
	if (size(zz,3) > 3)
		if (max(zz(:)) == 1),	zz = logical(zz);	handles.head(6) = 1;	end	% att.Metadata fails to detect NBITS
		FileOpenGDALmultiBand_CB(handles, 'multiband', zz, att)
		recentFiles(guidata(handles.figure1));	% Insert fileName into "Recent Files" & save handles
		return
	end
	X = handles.head(1:2);			Y = handles.head(3:4);
	pal = gray(256);
	if ( strcmp(att.ColorInterp,'Palette') && ~isempty(att.Band(1).ColorMap) )
		pal = att.Band(1).ColorMap.CMap(:,1:3);
	end
	set(handles.figure1,'Colormap',pal)

	handles = show_image(handles,handles.fileName,X,Y,zz,0,'xy',1);
	grid_info(handles,att,'gdal')			% Construct a info message
	handles = aux_funs('isProj',handles);	% Check about coordinates type
	handles = setAxesDefCoordIn(handles,1);	% Sets the value of the axes uicontextmenu that selects whether project or not
	recentFiles(handles);					% Insert fileName into "Recent Files" & save handles

% --------------------------------------------------------------------
function FileOpenNewImage_CB(handles, opt)
	str1 = {'*.jpg', 'JPEG image (*.jpg)'; ...
		'*.png', 'Portable Network Graphics(*.png)';		'*.bmp', 'Windows Bitmap (*.bmp)'; ...
		'*.gif', 'GIF image (*.gif)';						'*.tif', 'Tagged Image File (*.tif)'; ...
		'*.pcx', 'Windows Paintbrush (*.pcx)';				'*.ras', 'SUN rasterfile (*.ras)'; ...
		'*.ppm', 'Portable Pixmap (*.ppm)';					'*.pgm', 'Portable Graymap (*.pgm)'; ...
		'*.raw;*.bin', 'RAW file (*.raw,*.bin)'; 			'*.shade', 'IVS shade File (*.shade)'; ...
		'*.xwd', 'X Windows Dump (*.xwd)';					'*.*', 'All Files (*.*)'};
	if (nargin == 1)
		[FileName,PathName,handles] = put_or_get_file(handles,str1,'Select image format','get');
		if isequal(FileName,0),		return,		end
	else				% Filename was transmited in input
		PathName = [];		FileName = opt;
	end
	handles.fileName = [PathName FileName];

	head_fw = [];			% Used when check for a .*fw registering world file
	[PATH,FNAME,EXT] = fileparts(handles.fileName);
	if (strcmpi(EXT,'.shade'))
		[fid, msg] = fopen(handles.fileName, 'r');
		if (fid < 0)
			set(handles.figure1,'pointer','arrow');		errordlg([handles.fileName ': ' msg],'ERROR'),	return
		end
		fseek(fid, 56, 'bof');					% Seek forward to the image data.
		nm = fread(fid,2,'uint16');				% read n_row & n_col
		n_row = nm(1);		n_col = nm(2);
		fseek(fid, 87, 'bof');					% position the pointer at the end of the header
		nbytes = n_row*n_col*4;					% image is of RGBA type, so it has 4 channels
		I = fread(fid,nbytes,'*uint8');		fclose(fid);
		I = reshape(I, [4 n_row n_col]);	I = permute(I, [3 2 1]);
		I(:,:,1) = flipud(I(:,:,1));		I(:,:,2) = flipud(I(:,:,2));
		I(:,:,3) = flipud(I(:,:,3));		I(:,:,4) = flipud(I(:,:,4));
		I = I(:,:,2:4);							% strip alpha off of I
	elseif (strcmpi(EXT,'.raw') || strcmpi(EXT,'.bin'))
		FileOpenGDALmultiBand_CB(handles, 'RAW', handles.fileName)
		return		% We are done here. Bye Bye.
	else
		info_img = imfinfo(handles.fileName);		% This and att are repeated but not 100%
		if (strcmpi(EXT,'.ppm') || strcmpi(EXT,'.pnm') || strcmpi(EXT,'.pgm')  || strcmpi(EXT,'.pbm') || ...
				strcmpi(EXT,'.pcx') || strcmpi(EXT,'.ras'))		% GDAL screws on these
			try			I = imread(handles.fileName);
			catch,		errordlg(lasterr,'Error'),		return	% It realy may happen
			end
			att.RasterCount = size(I, 3);
		else
			% Try with GDAL and if it fails try again with imread
			try
				[I, att] = gdalread(handles.fileName);
			catch
				try
					I = imread(handles.fileName);
					att.RasterCount = size(I, 3);
				catch
					errordlg(lasterr, 'Error')
				end
			end
		end
		try		handles.transparency = info_img.Transparency;	end
		if (att.RasterCount > 4)
			% Animatted images.	BUT WORK ONLY WITH INDEXED IMAGES, OTHERWISE ... DON'T KNOW WHAT ERROR
			handles.cinemaImgs = I;
			I = I(:,:,1);
			handles.hUIcinanim = cine_animations(handles);
			guidata(handles.figure1,handles)
			ud(3) = att.RasterCount;		ud(2) = 1;		ud(1) = 1;
			set(handles.hUIcinanim(1), 'UserData', ud)
			if (isfield(info_img(1),'ColorTable'))			% Gif images call it 'ColorTable'
				set(handles.hUIcinanim(3), 'UserData', {info_img.ColorTable})
			elseif (isfield(info_img(1),'Colormap') && ~isempty(info_img(1).Colormap))
				set(handles.hUIcinanim(3), 'UserData', {info_img.Colormap})
			end
		end
		[head_fw,err_msg] = tfw_funs('inquire',[size(I,1) size(I,2)],PATH,FNAME,EXT);	% See if we have .*fw file
		if (~isempty(err_msg))
			warndlg(['A registering world file was found but the following error occured: ' err_msg],'Warning')
		end

		if (strcmp(info_img(1).ColorType,'grayscale') || (strcmp(info_img(1).ColorType,'truecolor') && (ndims(I) ~= 3)) )
			set(handles.figure1,'Colormap',gray(256))
		elseif (isfield(info_img(1),'ColorTable'))			% Gif images call it 'ColorTable'
			set(handles.figure1,'Colormap',info_img(1).ColorTable)
		elseif (isfield(info_img(1),'Colormap') && ~isempty(info_img(1).Colormap))
			set(handles.figure1,'Colormap',info_img(1).Colormap)
		end
	end

	if (isempty(head_fw))			% Image is not georeferenced
		handles.head = [1 size(I,2) 1 size(I,1) 0 255 0 1 1];	% Fake a grid reg GMT header
		handles.image_type = 2;		X = [];		Y = [];			ax_dir = 'off';
	else							% Got and decoded a .*fw file
		handles.image_type = 3;		handles.head = head_fw;
		X = handles.head(1:2);		Y = handles.head(3:4);		ax_dir = 'xy';
		I = flipdim(I,1);
	end
	handles = show_image(handles,handles.fileName,X,Y,I,0,ax_dir,0);
	grid_info(handles,handles.fileName,'iminfo');	% Construct a info string
	handles = aux_funs('isProj',handles);			% Check/set about coordinates type
	guidata(handles.figure1,handles)
	if (~isempty(head_fw) && ishandle(handles.Projections))		% Case still unknown to aux_funs('isProj',... 
		set(handles.Projections,'Enable','on')		% We don't know which but at least it is georeferenced
		setappdata(handles.figure1,'ProjGMT','')
	end
	recentFiles(handles);		% Insert fileName into "Recent Files" & save handles

% --------------------------------------------------------------------
function FileOpenWebImage_CB(handles, fname, opt)
% Load a file directly from the web.
% OPT, [optional] informs that we are to get the ClusterMaps Mirone visitors Map
	if (nargin == 1)
		resp = inputdlg({'Enter Image URL'},'Get Image from Web',[1 120]);
		if (isempty(resp)),		return,		end
		fname = resp{1};
	end
	I = gdalread(fname);
	if (isempty(I))
		errordlg('Sorry, but failed to fetch image from the ether. A Proxy problem?', 'Error'),		return
	end
	if (nargin == 3)
		handles.head = [-180 180 -57 78 0 255 0 360/(size(I,2)-1) 135/(size(I,1)-1)];
		X = [-180 180];		Y = [-57 78];		ax_dir = 'xy';		I = flipdim(I,1);
		handles.geog = 1;	handles.image_type = 3;
		if (~handles.no_file)		% When current Mirone fig is not empty, display the visitors map on a new one
			tmp.head = handles.head;	tmp.X = X;		tmp.Y = Y;		tmp.geog = 1;	tmp.name = 'Visitors';
			mirone(I, tmp),		return
		end
	else
		handles.head = [1 size(I,2) 1 size(I,1) 0 255 0 1 1];	% Fake a grid reg GMT header
		handles.image_type = 2;		X = [];		Y = [];			ax_dir = 'off';
	end
	handles = show_image(handles,fname,X,Y,I,0,ax_dir,0);
	handles = aux_funs('isProj',handles);			% Check/set about coordinates type
	guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function FileOpenGDALmultiBand_CB(handles, opt, opt2, opt3)
% Read GDAL files that may be multiband
% OPT2, if present, MUST contain either the full file name OR a multiband array. Currently used to load RAW images
% OPT3, if present, contains the gdalread metadata structure.

	if strcmp(opt,'ENVISAT')
		str1 = {'*.n1;*.N1', 'Envisat (*.n1,*.N1)'; '*.*', 'All Files (*.*)'};
	elseif strcmp(opt,'AVHRR')
		str1 = {'*.n14;*.N14;*n17;*N17', 'NOAA (*.n14,*.N14,*.n17,*.N17)'; '*.*', 'All Files (*.*)'};
	end

	if (nargin == 2)		% Otherwise, OPT2 already contains the File name
		[FileName,PathName] = put_or_get_file(handles,str1,['Select ' opt ' File'],'get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	else
		fname = opt2;
	end

	att.ProjectionRef = [];			X = [];			Y = [];		reseta = false;		ax_dir = 'off';
	reader = 'GDAL';				% this will be used by bands_list to decide which reader to use

	set(handles.figure1,'pointer','watch')
	if (strcmp(opt,'RAW'))
		[I,cmd1,cmd2] = read_FlatFile({fname});
		if (isempty(I)),	set(handles.figure1,'pointer','arrow'),		return,		end
		[att.RasterYSize, att.RasterXSize, n_bands] = size(I);
		bands_inMemory = 1:n_bands;			% Make it a vector
		reader = {cmd1; cmd2};
		handles.head = [1 size(I,2) 1 size(I,1) 0 255 0 1 1];	% Fake GMT header
		handles.image_type = 2;		handles.geog = 0;
	elseif (strcmp(opt,'ENVISAT') || strcmp(opt,'AVHRR'))
		att = gdalread(fname,'-M', '-C');
		bands_inMemory = 10;				% AD-HOC
		opt_U = ' ';
		if (~isempty(att.GeoTransform)),	opt_U = '-U';	end
		opt_B = sprintf('-B1-%d', bands_inMemory);
		I = gdalread(fname,'-S', opt_B, opt_U, '-C');
		n_bands = att.RasterCount;
		if (n_bands == 4)					% We might have a situation here. A Geotiff with transparency
			[P,N,EXT] = fileparts(fname);
			if (strncmpi(EXT, '.tif', 4))	% For tiffs/geotiffs check with matlab
				info_img = imfinfo(fname);	% GDAL doesn't tell us about transparency
				try	handles.transparency = info_img.Transparency;	end
			end
		end
		bands_inMemory = 1:min(n_bands,bands_inMemory);			% Make it a vector
		handles.head = att.GMT_hdr;
		handles.image_type = 2;		handles.geog = 0;
		if (~isempty(att.ProjectionRef))							% Georeferenced image
			X = handles.head(1:2);		Y = handles.head(3:4);
			handles.head(8) = diff(X) / (size(I,2) - 1);		handles.head(9) = diff(Y) / (size(I,1) - 1);
			ax_dir = 'xy';				handles.image_type = 3;
			att.GMT_hdr(8:9) = handles.head(8:9);				% Update the attrib struct
		end
		if (~isempty(att.GCPvalues))	% Save GCPs so that we can plot them and warp the image
			setappdata(handles.figure1,'GCPregImage',att.GCPvalues)
			setappdata(handles.figure1,'fnameGCP',fname)	% Save this to know when GCPs are to be removed
		end													% from appdata. That is donne in show_image()
	elseif (strncmp(opt,'PCA',3) || strncmp(opt,'mul',3))	% Generic multiband file transmited in input
		I = fname;		fname = opt;	nome_tmp = handles.fileName;	opt_bak = opt;
		if (~isempty(nome_tmp)),	fname = nome_tmp;	[pato, opt] = fileparts(fname);		reseta = true;	end
		if (nargin == 4),	att = opt3;		end
		[att.RasterYSize, att.RasterXSize, n_bands] = size(I);	% Use 'att' to be consistent with the other cases
		bands_inMemory = 1:n_bands;
		if (handles.image_type == 3)
			X = handles.head(1:2);		Y = handles.head(3:4);		ax_dir = 'xy';
		end
		reader = [];
	end

	aux_funs('toBandsList', handles.figure1, I, opt, fname, n_bands, bands_inMemory, reader);
	if (reseta),	opt = opt_bak;		end		% Wee need it with a known value for the remaining tests 

	handles.fileName = [];		% Not eligible for automatic re-loading
	handles.was_int16 = false;	handles.computed_grid = 0;
 	if ( any(strcmp(opt, {'ENVISAT' 'AVHRR' 'multiband'})) )
		cmap = att.Band(1).ColorMap;
		if (isempty(cmap)),		set(handles.figure1,'Colormap',jet(256))
		else					set(handles.figure1,'Colormap',cmap(1).CMap(:,1:3))
		end
		handles = recentFiles(handles);			% Insert fileName into "Recent Files"
	end
	if (n_bands == 1),		set(handles.figure1,'Colormap',gray(256)),		end		% Takes precedence over the above
	if (islogical(I)),		I = I(:,:,1);		end		% We don't want to try to compose an RGB out of logicals
	handles = show_image(handles,fname,X,Y,I,0,ax_dir,0);	% It also guidata(...) & reset pointer
	if (isappdata(handles.axes1,'InfoMsg')),	rmappdata(handles.axes1,'InfoMsg'),		end
	if ( any(strcmp(opt, {'ENVISAT' 'AVHRR' 'multiband'})) ),		grid_info(handles,att,'gdal');		end		% Construct a info message
	aux_funs('isProj',handles,1);				% Check/set about coordinates type
	setAxesDefCoordIn(handles);					% Sets the value of the axes uicontextmenu that selects whether project or not

% --------------------------------------------------------------------
function erro = FileOpenGeoTIFF_CB(handles, tipo, opt)
	if (nargin == 2)
		switch lower(tipo)
			case 'geotiff',		str1 = {'*.tif;*.TIF;*.tiff;*.TIFF', 'GeoTiff (*.tif,*.tiff,*.TIF,*.TIFF)'};
			case 'sid',			str1 = {'*.sid;*.SID', 'Mr Sid (*.sid,*.SID)'};
			case 'ecw',			str1 = {'*.ecw;*.ECW', 'ECW (*.ecw,*.ECW)'};
			case 'jp2',			str1 = {'*.jp2;*.JP2', 'Jpeg 2000 (*.jp2,*.JP2)'};
			case 'bil',			str1 = {'*.bil;*.BIL;', 'ESRI BIL (*.bil,*.BIL)'; '*.*', 'All Files (*.*)'};
			case 'bsb',			str1 = {'*.kap;*.KAP;*.nos;*NOS', 'BSB Nautical Chart (*.kap,*.KAP,*.nos,*.NOS)'; '*.*', 'All Files (*.*)'};
			otherwise,			str1 = {'', 'Don''t know (*.*)'};		% Used with the "Try Luck" option
		end
		str1(2,1:2) = {'*.*', 'All Files (*.*)'};
		[FileName,PathName,handles] = put_or_get_file(handles,str1,['Select ' tipo ' file'],'get');
		if isequal(FileName,0),		return,		end
	else				% Filename was transmited in input
		PathName = [];	FileName = opt;
	end
	handles.fileName = [PathName FileName];

	ECWpatch(handles, tipo)					% Check if we need to patch the ECW library memory fragmentation
	erro = 0;			gotHDRcoords = false;		opt_U = '-U';		opt_L = ' ';		fnameBak = handles.fileName;
	try		att = gdalread(handles.fileName,'-M','-C');	% Safety against blind calls from 'try luck'
	catch,	erro = 1;	errordlg(lasterr,'Error'),		return
	end
	att.fname = handles.fileName;						% If hdfread is used, it will need the file name (not eventual dataset name)
	att.subDsName = '';

	if (att.RasterCount == 0 && ~isempty(att.Subdatasets))
		str = strrep(att.Subdatasets, '=', ' ');
		c = false(1, numel(str));
		c3D = false(1, numel(str));						% To keep trace of eventual 3D datasets
		for (k = 2:2:numel(str))						% Seek for non-interesting (params) arrays 
			indF = strfind(str{k}, ']');
			ind = strfind(str{k}(1:indF), 'x');
			if (isempty(ind) || numel(ind) > 2)			% Don't want 1D or > 3D arrays
				c(k) = true;	c(k-1) = true;		continue
			elseif (numel(ind) == 1 && strcmp(str{k}(ind(1)-2:ind(1)-1), '[1'))		% Nor 1D singletons, e.g. [1x1557]
				c(k) = true;	c(k-1) = true;		continue
			elseif (numel(ind) == 2 && ~strcmp(str{k}(ind(1)-2:ind(1)-1), '[1'))	% But exclude 3D singletons, e.g. [1x1557x1557]
				c3D(k) = true;
			end
			if ((indF - ind) == 2),	c(k) = true;	c(k-1) = true;	end		% Don't want arrays with less than 10 (2 char) columns
		end
		if (~all(c) && any(c))		% Remove non-interesting arrays from sight
			str(c) = [];	att.Subdatasets(c) = [];	% Remove them also from the att.Subdatasets
			c3D(c) = [];
		elseif (all(c) && ~any(c3D))
			warndlg('This file has no 2D or 3D subdataset, so nothing I can do, Bye. Perhaps a OC SMI file?','Warning')
			return
		end

		SS = get(0,'ScreenSize');	nChars = 0;
		for (k = 1:numel(str)),		nChars = max(nChars, numel(str{k}));	end
		fac = 1;
		if (handles.IamCompiled),	fac = 1.05;		end
		[s,ok] = listdlg('PromptString',{'This file has subdatasets' 'you have to select one:'}, 'ListSize', ...
				[min(nChars*6*fac, SS(3)-100) min((size(str,1)*20 + 50), 200)], ...
				'Name','DATASET Selection', 'SelectionMode','single', 'ListString',str);	pause(0.01)
		if (~ok),	return,		end						% Uset hit "Cancel"
		if (rem(s,2) == 0),		s = s - 1;		end		% Selection was done over "description" and not the "name"

		ind = strfind(str{s}, ' ');
		FileName = str{s}(ind+1:end);					% First "ind" chars are of the form SUBDATASET_1_NAME=

		% Deal with case of smultaneous existence of one 3D and one or more 2D subdatasets
		if (any(c3D) && c3D(s+1))						% Got a selection of one 3D-SUBDATASET CASE
			handles.hMirFig = [];						% Informs aquamoto/slices that it has to create a new Fig to hold slice
			ind = strfind(FileName, ':');				% We will need this in slices/push_inputName_CB() to ...
			handles.ncVarName = FileName(ind(end)+1:end);	% New member (ncVarName)
			if (any(strncmp(att.Metadata,'NC_GLOBAL#TSU=',13)))
				aquamoto(handles, handles.fileName)
			else
				slices(handles, handles.fileName)
			end
			return
		end

		handles.fileName = FileName;
		ind = strfind(str{s+1}, ' ');
		subDsName = str{s+1}(ind(2)+1:ind(3)-1);		% Get dataset name
		AllSubdatasets = att.Subdatasets;				% Copy this for eventual use in "empilhador"
		att = gdalread(FileName,'-M','-C');				% Try again
		att.subDsName = subDsName;						% Another non-standard
		att.AllSubdatasets = AllSubdatasets;			% and another
	end
	att.fname = fnameBak;								% If hdfread is used, it will need the file name (not eventual dataset name)

	if (att.RasterCount == 0)			% Should never happen given the piece of code above, but ...
		errordlg('Probably a multi-container file. Could not read it since its says that it has no raster bands.','ERROR')
		return
	elseif (att.RasterCount > 3)		% Since it is a multiband file, try luck there
		FileOpenGDALmultiBand_CB(handles,'AVHRR',handles.fileName);		return
	end

	if (att.RasterCount == 2 && strcmp(att.Band(1).DataType,'Int16'))	% Good chances that it's a Radarsat-2 file
		loadRADARSAT(handles, att),		return
	end

	if (~strcmp(att.Band(1).DataType,'Byte'))			% JPK2, for example, may contain DTMs
		tipo = 'guess';			% Probably a temp step during the adaptation period after the OceanColor had f. changed format 
		if (strcmp(att.DriverShortName, 'netCDF') || strcmp(att.DriverShortName, 'HDF5Image')),		tipo = 'ncHDF';		end
		loadGRID(handles,handles.fileName, tipo, att);		return
	end

	if (strncmp(att.DriverShortName, 'HDF4', 4))		% Hmm, and what if it's a new nc/HDF file?
		tmp = att.GMT_hdr(1:4);			% It will change if we find coords in HDF variables not scanned by GDAL
		[head, slope, intercept, base, is_modis, is_linear, is_log, att] = empilhador('getFromMETA', att);
		if (~isequal(att.GMT_hdr(1:4), tmp)),		gotHDRcoords = true;		end
		if (isfield(att, 'hdrModisL2') && ~isempty(att.hdrModisL2)),	opt_L = '-L';	end
	end

	if (isempty(att.GeoTransform) && ~gotHDRcoords),		opt_U = ' ';	end
	Z = gdalread(handles.fileName, opt_U, opt_L);
	handles.head = att.GMT_hdr;

	if (~isempty(att.GeoTransform) || gotHDRcoords)		% Georeferenced image
		X = handles.head(1:2);		Y = handles.head(3:4);
		handles.head(8) = diff(X) / (size(Z,2) - 1);		handles.head(9) = diff(Y) / (size(Z,1) - 1);
		ax_dir = 'xy';				handles.image_type = 3;
		att.GMT_hdr(8:9) = handles.head(8:9);			% Update the attrib struct
	else							% Case of "raw" imagerie
		X = [1 size(Z,2)];		Y = [1 size(Z,1)];	ax_dir = 'off';		handles.image_type = 2;
		handles.head = [X Y 0 255 0 1 1];
		if (opt_L(1) == '-'),	ax_dir = 'xy';		handles.head = head;	end	% A dumb test for MODIS quality flags
	end

	if (strcmpi(att.ColorInterp,'palette'))
		if isempty(att.Band(1).ColorMap),		pal = jet(256);
		else
			try
				pal = att.Band(1).ColorMap.CMap;	pal = pal(:,1:3);		% GDAL creates a Mx4 colormap
				if ((handles.head(6)+1) < size(pal,1) && isequal(pal(handles.head(6)+2,1:end), [0 0 0]) )
					pal(handles.head(6)+2:end,:) = [];		% TEMP!? GDAL forces pals of 256 and that screws some cases
				end
			catch
				warndlg('Figure ColorMap had troubles. Replacing by a default one.','Warning'); pal = jet(256);
			end
		end
	elseif (strcmpi(att.ColorInterp,'gray'))
		pal = repmat( (att.GMT_hdr(5):att.GMT_hdr(6))' / att.GMT_hdr(6), 1, 3);
	elseif (strcmpi(att.ColorInterp,'Undefined') && strcmp(att.Band.DataType, 'Byte'))
		if (att.GMT_hdr(5) ~= att.GMT_hdr(6))		% Strong possibility that this is shit from GDAL. Check it
			att.GMT_hdr(5) = double(min(Z(:)));		att.GMT_hdr(6) = double(max(Z(:)));
		end
		if (att.GMT_hdr(5) ~= att.GMT_hdr(6))		% Right, gdalread/GDAL screwed here
			pal = gray(diff(att.GMT_hdr(5:6))+1);	% This is probably not good enough (think if values are [0 1 3 6 255])
		else
			pal = gray(2);
		end
	else
		pal = gray(256);		% Fall back case
	end

	if (~isempty(att.GCPvalues) && isempty(att.GeoTransform))		% Save GCPs so that we can plot them and warp the image
		setappdata(handles.figure1,'GCPregImage',att.GCPvalues)
		setappdata(handles.figure1,'fnameGCP',handles.fileName)		% Save this to know when GCPs are to be removed
	end																% from appdata. That is donne in show_image()

	set(handles.figure1,'Colormap',pal);
	handles = show_image(handles,att.fname,X,Y,Z,0,ax_dir,0);
	grid_info(handles,att,'gdal')				% Construct a info message
	handles = aux_funs('isProj',handles);		% Check/set about coordinates type
	handles = setAxesDefCoordIn(handles,1);		% Sets the value of the axes uicontextmenu that selects whether project or not
	recentFiles(handles);						% Insert fileName into "Recent Files" & save handles

% --------------------------------------------------------------------
function ECWpatch(handles, tipo)
% This is a tortuous solution to overcome the shity job that the ECW library performs on the memory continuity
% The fact is that the ECW dll chooses to sit rigth in the midle of the largest memory block, thus highly reducing its usability
	if (strcmp(tipo,'jp2') || strcmp(tipo,'ecw'))
		clear gdalread gdalwrite
		set_gmt(['GDAL_DRIVER_PATH=' handles.home_dir]);
	end

% --------------------------------------------------------------------
function FileOpenDEM_CB(handles, opt)
% Files of the following formats are read (well re-directed) here
	tipo = opt;
	switch opt
		case {'GMT' 'Surfer', 'ENCOM', 'GSOFT'}
			str1 = {'*.grd;*.GRD;*.nc;*.NC', 'Grid files (*.grd,*.GRD,*.nc,*.NC)';'*.*', 'All Files (*.*)'};	tipo = 'GMT_relatives';
		case 'SSimg'
			str1 = {'*.img;*.IMG', 'Grid files (*.img,*.IMG)';'*.*', 'All Files (*.*)'};
		case 'MANI'
			str1 = {'*.man;*.MAN', 'Grid files (*.man,*.MAN)';'*.*', 'All Files (*.*)'};	tipo = 'GMT_relatives';
		case 'ArcAscii',		str1 = {'*.asc;*.ASC', 'Arc/Info grid (*.asc,*.ASC)'; '*.*', 'All Files (*.*)'};
		case 'ArcBinary',		str1 = {'*.adf;*.ADF', 'Arc/Info grid (*.adf,*.ADF)'; '*.*', 'All Files (*.*)'};
		case 'DTED',			str1 = {'*.dt0;*.DT0;*.dt1;*.DT1', 'DTED (*.dt0,*.DT0,*.dt1,*.DT1)'; '*.*', 'All Files (*.*)'};
		case 'GTOPO30',			str1 = {'*.dem;*.DEM', 'GTOPO30 DEM (*.dem,*.DEM)'; '*.*', 'All Files (*.*)'};
		case 'GXF',				str1 = {'*.gxf;*.GXF', 'Geosoft GXF (*.gxf,*.GXF)'; '*.*', 'All Files (*.*)'};
		case 'MOLA',			str1 = {'*.img;*.IMG', 'MOLA DEM (*.img,*.IMG)'; '*.*', 'All Files (*.*)'};
		case 'SDTS',			str1 = {'*catd.ddf;*CATD.DDF', 'USGS SDTS DEM (*catd.ddf,*CATD.DDF)'; '*.*', 'All Files (*.*)'};
		case 'SRTM30',			str1 = {'*.srtm;*.SRTM;*.srtm.gz', 'SRTM30 DEM (*.srtm,*.SRTM,*.srtm.gz)'; '*.*', 'All Files (*.*)'};
		case {'SRTM1' 'SRTM3'},	str1 = {'*.hgt;*.HGT;*.hgt.zip', [opt ' DEM (*.hgt,*.HGT,*.hgt.zip)']; '*.*', 'All Files (*.*)'};
		case 'USGS_DEM',		str1 = {'*.dem;*.DEM', 'USGS DEM (*.dem,*.DEM)'; '*.*', 'All Files (*.*)'};
		otherwise
			errordlg(['OOPs, where did this ' opt ' code came in?'],'Error'),	return
	end
	[FileName,PathName,handles] = put_or_get_file(handles,str1,['Select ' opt ' File'],'get');
	if isequal(FileName,0),		return,		end
	fullname = [PathName FileName];
	loadGRID(handles,fullname,tipo)

% --------------------------------------------------------------------
function loadGRID(handles, fullname, tipo, opt)
% This function loads grid files that may contain DEMs or other grid (not images) types
% TIPO	-> 'GMT'		read grids using the read_gmt_type_grids() function
%		   'MOLA_lbl'	To read Mars MOLA .img with a .lbl header file (cannot be compressed)
%		   'MOLA'		To read Mars MOLA .img with a .lbl header file *OR* a V3 PDS .img
%		   'whatever'	Let GDAL guess what to do (it means, any string)
% OPT	-> the "att" attributes structure got from att = gdalread(fname,'-M',...)
 
	if (nargin == 3),	opt = ' ';	end
	[Z, X, Y, srsWKT, handles, att] = read_grid(handles, fullname, tipo, opt);
	if (isempty(Z)),	return,		end
	if (~isempty(srsWKT)),	handles.is_projected = true;	end		% WEAK TEST
	if (~isempty(fullname))
		pato = fileparts(fullname);
		ind1 = strfind(fullname, ':"');		% Check if fulname is a subdataset name as 'NETCDF:"v:\tsu\lagos.nc":bathymetry'
		if (~isempty(ind1))
			ind2 = strfind(fullname, '":');
			if (~isempty(ind2))				% It seams it was a subdatset name
				pato = fileparts(fullname(ind1+2:ind2-1));
			end
		end
		handles.last_dir = pato;
	end

	if (~isempty(att) && ~isempty(att.GCPvalues))					% Save GCPs so that we can plot them and warp the image
		setappdata(handles.figure1,'GCPregImage',att.GCPvalues)
		setappdata(handles.figure1,'fnameGCP',handles.fileName)		% Save this to know when GCPs are to be removed
	end																% from appdata. That is donne in show_image()

	aux_funs('StoreZ',handles,X,Y,Z)				% If grid size is not to big we'll store it
	aux_funs('colormap_bg',handles,Z,jet(256));		% Insert the background color in the palette for arrays that have NaNs
	if (~isa(Z, 'uint8'))
		zz = scaleto8(Z);
		if (isa(Z,'int8') && handles.head(6) < 255)	% Put in a colormap that uses all data range 
			set(handles.figure1,'Colormap',jet(handles.head(6)+1))
		end
	else
		zz = Z;
	end
	handles = show_image(handles,handles.fileName,X,Y,zz,1,'xy',handles.head(7));
	if (isappdata(handles.axes1,'InfoMsg')),	rmappdata(handles.axes1,'InfoMsg'),		end
	if (~isempty(att))
		grid_info(handles,att,'gdal');				% Construct a info message and save proj (if ...)
		handles = setAxesDefCoordIn(handles,1);		% Sets a value on the axes uicontextmenu
	elseif (~isempty(srsWKT))
		aux_funs('appP', handles, srsWKT)			% If we have a WKT proj store it
		handles = setAxesDefCoordIn(handles,1);
	else
		if (isappdata(handles.figure1,'ProjWKT')),	rmappdata(handles.figure1,'ProjWKT'),	end
		if (isappdata(handles.figure1,'ProjGMT')),	rmappdata(handles.figure1,'ProjGMT'),	end
		if (isappdata(handles.figure1,'Proj4')),	rmappdata(handles.figure1,'Proj4'),		end
	end
	handles = aux_funs('isProj',handles);			% Check about coordinates type
	recentFiles(handles);							% Insert fileName into "Recent Files" & save handles
	if (handles.nLayers > 1)						% nLayers was found inside the read_grid() call
		handles.hMirFig = handles.figure1;			% Informs slices of current fig handle but do not save it
		slices(handles, fullname)
	end

% -------------------------------------------------------------------------------------------------
function loadRADARSAT(handles, att)
% Take a geotiff with 2 bands of 2 bytes each and read the real & imaginary components of a RADARSAT-2

	Z = single(gdalread(handles.fileName, '-B1'));
	showRADARSAT(handles, att, 'RADARSAT-Real-component', Z)
	Z = single(gdalread(handles.fileName, '-B2'));
	fname = handles.fileName;		% Make a copy
	h = mirone;		handles = guidata(h);		handles.fileName = fname;
	showRADARSAT(handles, att, 'RADARSAT-Imag-component', Z)

function showRADARSAT(handles, att, fname, Z)
% Plot a RADARSAT-2 band. Unfortunately we still do not allow storing Z as a 16 bits array.

	if (isempty(att.ProjectionRef))								% No georeferenced image
		X = 1:att.RasterXSize;		Y = 1:att.RasterYSize;
		handles.image_type = 1;
		att.GMT_hdr(1:4) = [X(1) X(end) Y(1) Y(end)];			% Make them grid reg
	else
		error('FDeu-se')
	end
	handles.head = att.GMT_hdr;
	handles.geog = true;
	if (~isempty(att.GCPvalues))	% Save GCPs so that we can plot them and warp the image
		setappdata(handles.figure1,'GCPregImage',att.GCPvalues)
		setappdata(handles.figure1,'fnameGCP', fname)			% Save this to know when GCPs are to be removed
	end															% from appdata. That is donne in show_image()
	if (~isempty(att.GCPprojection))
		aux_funs('appP', handles, att.GCPprojection)			% If we have a WKT proj, store it
	end

	zz = grdutils(Z,'-L');
	handles.head(5:6) = [zz(1) zz(2)];
	aux_funs('StoreZ',handles,X,Y,Z)				% If grid size is not to big we'll store it
	aux_funs('colormap_bg',handles,Z,gray(256));
	img = scaleto8(Z, 8, handles.head(5:6));
	handles = show_image(handles, fname, X, Y, img, 1, 'xy', false);

	grid_info(handles,att,'gdal');			% Construct a info message and save proj (if ...)
	setAxesDefCoordIn(handles,1);			% Sets a value on the axes uicontextmenu

% _________________________________________________________________________________________________
% -*-*-*-*-*-*-$-$-$-$-$-$-#-#-#-#-#-#-%-%-%-%-%-%-@-@-@-@-@-@-(-)-(-)-(-)-&-&-&-&-&-&-{-}-{-}-{-}-
function handles = show_image(handles, fname, X, Y, I, validGrid, axis_t, adjust, imSize)
% Show image and set other parameters
% AXIS_T	'xy' or 'off'
% IMSIZE	Optional arg. Use when dx ~= dy

	if (adjust)				% Convert the image limits from pixel reg to grid reg
		[X,Y] = aux_funs('adjust_lims', X, Y, size(I,1), size(I,2));
	end
	if (strcmp(get(handles.GCPtool,'Checked'),'on'))		% Call gcpTool() and return
		if (handles.no_file),	aux_funs('togCheck', handles.GCPtool),	return,		end	% Security check
		handles = gcpTool(handles,axis_t,X,Y,I);	return
	end
	if (nargin < 9),	imSize = [];	end
	if (~isa(handles.head, 'double')),		handles.head = double(handles.head);	end
	if (~isa(X, 'double')),		X = double(X);	Y = double(Y);	end		% Security measure (+ lines above & below). It happened before
	if (handles.head(8) == 0 || handles.head(9) == 0)	% Bad usage of the (Z,struct) mirone input mechanism
		handles.head(8) = diff(handles.head(1:2)) / (size(I,2) - ~handles.head(7));		% Get right values for them
		handles.head(9) = diff(handles.head(3:4)) / (size(I,1) - ~handles.head(7));
	end

	if ((handles.image_type ~= 2 && handles.image_type ~= 20) && (abs(diff(handles.head(8:9))) > 1e-5))	% Check anisotropy
		imSize = handles.head(8) / handles.head(9);		% resizetrue will know what to do with this
	end

	if (~validGrid && handles.validGrid),		aux_funs('cleanGRDappdata',handles);	end
	alpha = 1;
	if (size(I,3) > 3)
		if (strcmp(handles.transparency, 'alpha'))
			alpha = I(:,:,4);
		end
		I(:,:,4:end) = [];		% Make sure I is only MxN or MxNx3
	elseif (size(I,3) == 2),	I(:,:,2) = [];			% (could be otherwise when input from multiband)
	end

	% The following applies only to the case where we have a .def region to create GMT custom symbols
	% The first time it passes here we have a square image and the next code does actually nothing.
	% In the eventuality of loading an image we will force it to have limits [-0.5 0.5] and no deformation.
	% In this case (load an image) we don't call resizetrue() because it changes sizes in a way that we can't allow.
	do_resize = true;
	if (handles.is_defRegion)
		if (size(I,1) > size(I,2))		% Image is taller
			X = [-((size(I,2) / size(I,1)) / 2) ((size(I,2) / size(I,1)) / 2)];		Y = [-0.5 0.5];
		else
			Y = [-((size(I,1) / size(I,2)) / 2) ((size(I,1) / size(I,2)) / 2)];		X = [-0.5 0.5];
		end
		if (handles.image_type ~= 20)
			I = flipdim(I,1);		% Need to because we have YDir = Normal (probably need to test for referenced imgs)
			alpha = flipdim(alpha,1);
			do_resize = false;	% When other than a back-ground image we don't want to call resizetrue
		end
		handles.head(1:4) = [X Y];
		handles.head(8) = diff(handles.head(1:2)) / (size(I,2) - ~handles.head(7));
		handles.head(9) = diff(handles.head(3:4)) / (size(I,1) - ~handles.head(7));
	end

	handles.hImg = image(X,Y,I,'Parent',handles.axes1,'AlphaData',alpha);
	if (handles.is_defRegion && handles.image_type ~= 20)	% Need to reset Lims to [-0.5 0.5] in order to center the image
		set(handles.axes1, 'XLim', [-0.5 0.5], 'YLim', [-0.5 0.5])
	end
	if (handles.is_defRegion)
		set(handles.axes1, 'YDir', 'Normal')	% It was bloody set to reverse by the previous image() call.
	end

	zoom_state(handles,'off_yes')
	if (islogical(I))
		set(handles.hImg,'CDataMapping','scaled');		set(handles.figure1,'ColorMap',gray(16));
	else
		set(handles.hImg,'CDataMapping','direct')
	end
	if (handles.geog < 10),		handles.geog = aux_funs('guessGeog',handles.head(1:4));
	else						handles.geog = handles.geog - 10;		% In this case we believe in the pre-set value
	end
	if (handles.image_type == 2),	handles.geog = 0;	end

	if (do_resize)		% If FALSE image is visible already and we don't want to use resizetrue
		magRatio = resizetrue(handles,imSize,axis_t);			% -------> IMAGE IS VISIBLE HERE. <-------
	end

	handles.origFig = I;			handles.no_file = 0;
	handles.Illumin_type = 0;		handles.validGrid = validGrid;	% Signal that gmt grid opps are allowed
	if (~isempty(fname) && ~handles.is_defRegion)
		set(handles.figure1,'Name',[fname sprintf('  @  %d%%',magRatio)])
	end
	setappdata(handles.axes1,'ThisImageLims',[get(handles.axes1,'XLim') get(handles.axes1,'YLim')])
	handles.oldSize(1,:) = get(handles.figure1,'Pos');		% Save fig size to prevent maximizing
	handles.origCmap = get(handles.figure1,'ColorMap');		% Save original colormap 
	set(handles.ImgHist,'checked','off')
	% Make an extra copy of those to use in "restore" because they may be changed by 'bands_list()'
	handles.validGrid_orig = validGrid;			handles.was_int16_orig = handles.was_int16;
	handles.computed_grid_orig = handles.computed_grid;
	fix_axes_labels(handles)		% Make ALL labels have the same number of decimals.
	handles = SetAxesNumericType(handles);					% Set axes uicontextmenus
	if (handles.image_type ~= 1),	handles.grdname = [];	end
	if (~handles.have_nans),		set(handles.haveNaNs,'Vis','off')	% If no NaNs no need of these
	else							set(handles.haveNaNs,'Vis','on')
	end
	guidata(handles.figure1, handles);

	% Hide uicontrols that are useless to images only. First 4 are button are handles which don't stand hiding(???) 
	st = {'off' 'on'};
	set(handles.Projections,'Vis', st{(handles.image_type ~= 2) + 1})
	set(handles.noVGlist(6:end),'Vis', st{validGrid + 1}),		set(handles.noVGlist(1:5),'Ena', st{validGrid + 1})
	set([handles.Geography handles.MagGrav handles.Seismology handles.Plates],'Vis', st{(handles.image_type ~= 2) + 1})
	set(handles.noAxes,'Vis', st{~strcmp(axis_t,'off') + 1})
	set(handles.toGE,'Enable', st{min(handles.geog,1) + 1})
 	if (handles.is_projected),	set(handles.toGE,'Enable', 'on'),	end		% Regardless of the above line
	set(findobj(handles.Projections,'-depth',1,'Label','GMT project'), 'Vis', st{validGrid + 1})
	if (handles.geog),		set(handles.DrawGeogCirc,'Tooltip','Draw geographical circle')
	else					set(handles.DrawGeogCirc,'Tooltip','Draw circle')
	end

	GCPmemoryVis = 'off';
	if (isappdata(handles.figure1,'GCPregImage'))
		fnameGCP = getappdata(handles.figure1,'fnameGCP');
		if (~strcmp(fnameGCP,fname))		% It means the 'GCPregImage' app must be from a previous file
			rmappdata(handles.figure1,'GCPregImage')
		else
			GCPmemoryVis = 'on';
		end
	end
	set(handles.GCPmemory,'Visible',GCPmemoryVis)

	BL = getappdata(handles.figure1,'BandList');		% We must tell between fakes and true 'BandList'
	if (isempty(BL) || ((numel(BL{end}) == 1) && strcmp(BL{end},'Mirone')))
		if (ndims(I) == 3)		% Some cheating to allow selecting individual bands of a RGB image
			tmp1 = cell(4,2);	tmp2 = cell(4,2);		tmp1{1,1} = 'RGB';		tmp1{1,2} = 'RGB';
			for (i = 1:3)
				tmp1{i+1,1} = ['band' sprintf('%d',i)];		tmp1{i+1,2} = ['banda' sprintf('%d',i)];
				tmp2{i+1,1} = [sprintf('%d',i) sprintf('%d',size(I,1)) 'x' sprintf('%d',size(I,2)) ' BSQ']; tmp2{i+1,2} = i;
			end
			tmp = {['+ ' 'RGB']; I; tmp1; tmp2; ''; 1:3; [size(I,1) size(I,2) 3]; 'Mirone'};
			setappdata(handles.figure1,'BandList',tmp)
			set(findobj(handles.Image,'-depth',1,'Label','Load Bands'), 'Vis','on')
			set(findobj(handles.Image,'-depth',1,'Label','Explore RGB'),'Vis','on')
		elseif (ndims(I) == 2)				% Remove it so it won't try to operate on indexed images
			if (isappdata(handles.figure1,'BandList')),		rmappdata(handles.figure1,'BandList'),	end
			set(findobj(handles.figure1,'Label','Load Bands'), 'Vis','off')
			set(findobj(handles.figure1,'Label','Explore RGB'),'Vis','off')
		end
	end
	if (isappdata(handles.axes1,'DatumProjInfo')),		rmappdata(handles.axes1,'DatumProjInfo'),	end
	% Note that, when it applyies the above are rebuilt with a latter call to grid_info(handles,att,'gdal')

	% Sets the right "Image mode" checkbox on (well, we still miss the BW case)
	if (ndims(I) == 3),		hThis = handles.ImModRGB;
	else					hThis = handles.ImMod8cor;
	end
	aux_funs('togCheck',hThis, [handles.ImMod8cor handles.ImMod8gray handles.ImModBW handles.ImModRGB])
	set(handles.ImRestore,'UserData',hThis)			% We need it when restoring original image

	% See if we have a callback function to be executed here (set for example by TINTOL)
	hfun = getappdata(handles.figure1, 'show_image_RunCB');
	if (~isempty(hfun))
		feval(hfun{1}, hfun{2:end})
	end

% --------------------------------------------------------------------
function fix_axes_labels(handles)
% Fix the annoyingly ugly effect of the varying number of decimals in labels
% Trust in the estimate made by num2str.
	xTick  = get(handles.axes1, 'XTick');
	xLabel = num2str(xTick(:));
	ind    = strfind(xLabel(1,:), '.');
	n_dec  = size(xLabel, 2) - ind;
	xLabel = num2str(xTick(:), sprintf('%%.%df', n_dec));

	yTick  = get(handles.axes1, 'YTick');
	yLabel = num2str(yTick(:));
	ind    = strfind(yLabel(1,:), '.');
	n_dec  = size(yLabel, 2) - ind;
	yLabel = num2str(yTick(:), sprintf('%%.%df', n_dec));
	
	set(handles.axes1, 'XTickLabel', xLabel, 'YTickLabel', yLabel)

% --------------------------------------------------------------------
function ToolsMBplaningStart_CB(handles)
	if (aux_funs('msg_dlg',1,handles)),		return,		end		% Test geog & no_file
	if isempty(getappdata(handles.figure1,'dem_z'))		% Test if the grid is loaded in memory
		warndlg('Grid file is bigger than the declared "Grid Max Size". See "File -> Preferences"','Warning');
		return
	end

	prompt = {['The current value for the swath-width / water depth ratio is:   --> ' sprintf('%g',handles.swathRatio) '  <--']
			'If you want to change it, hit "Cancel", and do it in "File -> Preferences"'; ' '
			'NOTE: this message, once accepted, is shown only once. So, if in the midle'
			'of a session, you decide to change the Swath Ratio don''t forget to do it'
			'there. The current value will be applyied to all new tracks without further questions.'};
	if (handles.firstMBtrack == 1)
		button = questdlg(prompt,'Multibeam planing info','Accept','Cancel','Accept');
		if strcmp(button,'Accept'),		handles.firstMBtrack = 0;
		else							return
		end
	end

	zoom_state(handles,'maybe_off');
	[xp,yp,trackHand,barHand] = getline_mb;
	if (isempty(xp)),	zoom_state(handles,'maybe_on'),		return,		end
	
	% Now we have to join all trackHand's in a single handle
	handles.nTrack = handles.nTrack + 1;	% count to the number of tracks
	delete(trackHand(2:end));		trackHand(2:end) = [];
	% make and set tags strings to tracks and track's Bars
	tagL = sprintf('MBtrack%d',handles.nTrack);				tagB = sprintf('swath_w%d',handles.nTrack);
	set(trackHand,'XData',xp, 'YData',yp, 'Tag',tagL);		set(barHand, 'Tag',tagB)
	setappdata(trackHand,'swathRatio',handles.swathRatio)	% save the swathRatio in line's appdata
	% now set the barHand userdata that contains the vertex order
	for (i = 1:numel(xp)),		set(barHand(i),'Userdata',i),	end
	
	draw_funs(trackHand,'MBtrackUictx')			% Set track's uicontextmenu
	draw_funs(barHand,'MBbarUictx')				% Set track bar's uicontextmenu
	zoom_state(handles,'maybe_on');
	handles.hMBplot = trackHand;		guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function ToolsMBplaningImport_CB(handles)
	if (aux_funs('msg_dlg',1,handles)),		return,		end		% Test geog & no_file
	if (~handles.validGrid)
		errordlg('This operation is deffined only for images derived from DEM grids','Error'),	return
	end
	if isempty(getappdata(handles.figure1,'dem_x'))		% Test if the grid is loaded in memory
		warndlg('Grid file is bigger than the declared "Grid Max Size". See "File -> Preferences"','Warning'),	return
	end
	[FileName,PathName] = put_or_get_file(handles,{'*.dat;*.DAT', 'Data files (*.dat,*.DAT)'},'Select input xy file name','get');
	if isequal(FileName,0),		return,		end
	
	out  = load_xyz(handles,[PathName FileName]);
	resp = inputdlg({'Enter swath-width / water depth ratio'},'Multi-beam planing input',[1 38],{'3'});	pause(0.01);
	resp = str2double(resp{1});
	if (isnan(resp)),	return,		end
	if (iscell(out)),	n_segments = length(out);
	else				n_segments = 1;
	end
	h_line = zeros(1, n_segments);
	for (i = 1:n_segments)
		if (iscell(out)),	xy = out{i}(:,1:2);
		else				xy = out(:,1:2);
		end
		z = abs(bi_linear(getappdata(handles.figure1,'dem_x'),getappdata(handles.figure1,'dem_y'),...
			getappdata(handles.figure1,'dem_z'),xy(:,1),xy(:,2)));
		rad = abs(z) * (resp/2) / 111194.9;		% meters -> degrees
		handles.nTrack = handles.nTrack + 1;	% count the number of imported tracks
		% make tags strings to tracks and track's Bars
		tagL = ['MBtrack' sprintf('%d',handles.nTrack)];	tagB = ['swath_w' sprintf('%d',handles.nTrack)];
		nr = size(xy,1);
		az = zeros(nr,1);		h_circ = zeros(1,nr);
		az(2:nr) = azimuth_geo(xy(1:(nr-1),2), xy(1:(nr-1),1), xy(2:nr,2), xy(2:nr,1));
		az(1) = az(2);				% First and second points have the same azim
		for (k=1:nr)				% NOTE: the Tag is very important for line edition
			[lat1,lon1] = circ_geo(xy(k,2),xy(k,1),rad(k),[az(k)-90-1 az(k)-90+1],3);
			[lat2,lon2] = circ_geo(xy(k,2),xy(k,1),rad(k),[az(k)+90-1 az(k)+90+1],3);
			h_circ(k) = line('XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)], ...
				'Color', [.8 .8 .8],'LineWidth',4,'Tag',tagB,'Userdata',k,'Parent',handles.axes1);
		end
		h_line(i) = line('XData', xy(:,1), 'YData', xy(:,2),'Color','k','LineWidth',1,'Tag',tagL,'Parent',handles.axes1);
		setappdata(h_line(i),'swathRatio',str2double(resp{1}))	% save the swathRatio in line's appdata
		draw_funs(h_line(i),'MBtrackUictx')		% Set track's uicontextmenu
		draw_funs(h_circ,'MBbarUictx')			% Set track bar's uicontextmenu
	end
	handles.hMBplot = h_line(1);	guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function varargout = ImageIllumModel_CB(handles, opt)
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	if (nargin == 1),	opt = 'grdgradient_A';	end

	luz = shading_params(opt);
	if (isempty(luz))
		if (nargout),	varargout{1} = [];		end
		return
	end

	if (luz.illum_model == 1),		[varargout{1:nargout}] = ImageIllum(luz, handles, 'grdgrad_class');	% GMT grdgradient classic
	elseif (luz.illum_model == 2),	[varargout{1:nargout}] = ImageIllum(luz, handles, 'grdgrad_lamb');	% GMT grdgradient Lambertian
	elseif (luz.illum_model == 3),	[varargout{1:nargout}] = ImageIllum(luz, handles, 'grdgrad_peuck');	% GMT grdgradient Peucker
	elseif (luz.illum_model == 4),	[varargout{1:nargout}] = ImageIllum(luz, handles, 'lambertian');	% GMT Lambertian
	elseif (luz.illum_model == 5),	[varargout{1:nargout}] = ImageIllum(luz, handles, 'manip');		% ManipRaster
	elseif (luz.illum_model == 6),	[varargout{1:nargout}] = ImageIllum(luz, handles, 'hill');		% ESRI's hillshade
	else							ImageIllumFalseColor(luz, handles)				% False color
	end
	if (luz.illum_model > 6 && nargout),	varargout{1} = [];		end				% No reflectances here

% --------------------------------------------------------------------
function Reft = ImageIllum(luz, handles, opt)
% OPT ->  Select which of the GMT grdgradient illumination algorithms to use
% Illuminate a DEM file and turn it into a RGB image
% For multiple tries I need to use the original image. Otherwise each attempt would illuminate
% the previously illuminated image. An exception occurs when the image was IP but only for the
% first time, repeated illums will use the original img. Otherwise we would need to make another img copy

	[X,Y,Z,head] = load_grd(handles);	% If needed, load gmt grid again
	if isempty(Z),	return,		end		% An error message was already issued

	OPT_a = '-a1';
	if (sum(handles.bg_color) < 0.01),	OPT_a = ' ';	end		% Near black bg color has a different treatment

	if (strcmp(opt,'grdgrad_class'))		% GMT grdgradient classic illumination
		illumComm = sprintf('-A%.2f',luz.azim);
		if (handles.geog),	[R,offset,sigma] = grdgradient_m(Z,head,'-M',illumComm,'-Nt',OPT_a);
		else				[R,offset,sigma] = grdgradient_m(Z,head,illumComm,'-Nt',OPT_a);
		end
		handles.Illumin_type = 1;
		if (sigma < 1e-6),		sigma = 1e-6;	end		% We cannot let them be zero on sprintf('%.6f',..) somewhere else
		if (offset < 1e-6),		offset = 1e-6;	end
		handles.grad_offset = offset;	handles.grad_sigma = sigma;

	elseif (strcmp(opt,'grdgrad_lamb'))		% GMT grdgradient lambertian illumination
		illumComm = sprintf('-Es%.2f/%.2f',luz.azim, luz.elev);
		R = grdgradient_m(Z,head,illumComm ,OPT_a);
		handles.Illumin_type = 2;
% 		nhs = 256;	ncolors = 256;
% 		R = scaleto8(R);		A = handles.origFig;
% 		IND  = uint16(R+1) + 256*uint16(A) + 1;
% 		cmap = get(handles.figure1,'Colormap');
% 		cmap = cmap(:);
% 		cmap = bsxfun(@times,cmap,linspace(0,1,nhs));
% 		cmap = reshape(cmap,[ncolors 3 nhs]);
% 		cmap = permute(cmap,[3 1 2]);
% 		cmap = reshape(cmap,[ncolors*nhs 3]);
% 		cmapUINT8 = uint8(round(cmap*256));
% 		img = reshape(cmapUINT8(IND(:),:),[size(IND),3]);
% 		%img = ind2rgb8(IND,cmap);		% WTF this one screws (planty of blacks)? It should be the same but more efficient.

	elseif (strcmp(opt,'grdgrad_peuck'))	% GMT grdgradient Peucker illumination
		illumComm = '-Ep';
		R = grdgradient_m(Z,head,illumComm, OPT_a);
		handles.Illumin_type = 3;

	elseif (strcmp(opt,'lambertian'))		% GMT Lambertian lighting illumination
		illumComm = sprintf('-E%g/%g/%g/%g/%g/%g',luz.azim,luz.elev,luz.ambient,luz.diffuse,luz.specular,luz.shine);
		R = grdgradient_m(Z,head,illumComm, OPT_a);
		handles.Illumin_type = 4;

	elseif (strcmp(opt,'manip'))			% GMT5 has this one
		illumComm = sprintf('-Em%.2f/%.2f', luz.azim, luz.elev);
		R = grdgradient_m(Z,head,illumComm, OPT_a);
		handles.Illumin_type = 5;

	else			% ESRI's hillshade
		illumComm = sprintf('-Eh%.2f/%.2f', luz.azim, luz.elev);
		R = grdgradient_m(Z,head,illumComm, OPT_a);
		handles.Illumin_type = 6;
	end
	setappdata(handles.figure1,'illumComm',illumComm);		% Save these for ROI op & write_gmt_script
	setappdata(handles.figure1,'Luz',luz);					% Save these for ROI operations

	if (nargout)	% Send the reflectance back to caller and stop here
		Reft = R;	guidata(handles.figure1, handles);
		return
	end

	if (~handles.firstIllum && handles.is_draped)
		drapInfo = getappdata(handles.hImg, 'drapInfo');
		if (ishandle(drapInfo.hFig))
			handlesSon = guidata(drapInfo.hFig);
			ImageDrape_CB(handlesSon, drapInfo.alpha)
		else
			warndlg('You killed the drapping window (happy?), therefore I cannot continue this game.','Warning')
			set(handles.hImg,'CData',handles.origFig);		% Restore original image
			handles.is_draped = false;
		end
		img = get(handles.hImg,'CData');
	elseif (handles.firstIllum),img = get(handles.hImg,'CData');	handles.firstIllum = 0;
	else						img = handles.origFig;
	end
	if (ndims(img) == 2),		img = ind2rgb8(img,get(handles.figure1,'Colormap'));	end
	mex_illuminate(img,R)		% New. It can now operate insitu too
	clear R

	if (handles.have_nans && ~isequal(handles.bg_color, [1 1 1]) && ~strcmp(OPT_a,' '))	% Non-white or black NaN color requested
		ind = isnan(Z);			bg_color = uint8(handles.bg_color * 255);
		tmp = img(:,:,1);		tmp(ind) = bg_color(1);		img(:,:,1) = tmp;
		tmp = img(:,:,2);		tmp(ind) = bg_color(2);		img(:,:,2) = tmp;
		tmp = img(:,:,3);		tmp(ind) = bg_color(3);		img(:,:,3) = tmp;
	end
	
	set(handles.hImg,'CData',img),		refresh(handles.figure1)		% Crazzy beast does not always update the image !!!!!!!!!
	aux_funs('togCheck',handles.ImModRGB, [handles.ImMod8cor handles.ImMod8gray handles.ImModBW])
	guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function ImageIllumFalseColor(luz, handles)
% Illuminate a grid from 3 different directions and turn it into a RGB image

	[X,Y,Z,head] = load_grd(handles);
	if isempty(Z),		return,		end		% An error message was already issued
	D2R = pi/180;

	if (luz.mercedes_type == 1)
		opt_M = ' ';	opt_N = '-Nt';
		if (handles.geog),		opt_M = '-M';		end
		tmp = grdgradient_m(Z,head, sprintf('-A%f',luz.azim(1)), opt_M, opt_N);		%tmp(tmp < 0) = 0;
		zz(:,:,1) = uint8(cvlib_mex('CvtScale', tmp, 254, 1));
		tmp = grdgradient_m(Z,head, sprintf('-A%f',luz.azim(2)), opt_M, opt_N);		%tmp(tmp < 0) = 0;
		zz(:,:,2) = uint8(cvlib_mex('CvtScale', tmp, 254, 1));
		tmp = grdgradient_m(Z,head, sprintf('-A%f',luz.azim(3)), opt_M, opt_N);		%tmp(tmp < 0) = 0;
		zz(:,:,3) = uint8(cvlib_mex('CvtScale', tmp, 254, 1));
	else
		% For compatibilty reasons. But should be possible to do it with a call to grdgradient_m
		size_amp = luz.ampFactor;
		zz(:,:,1) = shade_manip_raster((luz.azim(1)-90)*D2R, luz.elev*D2R, Z, size_amp);
		zz(:,:,2) = shade_manip_raster((luz.azim(2)-90)*D2R, luz.elev*D2R, Z, size_amp);
		zz(:,:,3) = shade_manip_raster((luz.azim(3)-90)*D2R, luz.elev*D2R, Z, size_amp);
		zz = uint8((254 * zz) + 1);
	end

	if (isappdata(handles.figure1,'illumComm')),	rmappdata(handles.figure1,'illumComm'),		end
	set(handles.hImg,'CData',zz);					handles.Illumin_type = 7;
	aux_funs('togCheck',handles.ImModRGB, [handles.ImMod8cor handles.ImMod8gray handles.ImModBW])
	guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function ImageAnaglyph_CB(handles)
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	[X,Y,Z,head] = load_grd(handles);
	if isempty(Z),	return,		end		% An error message was already issued

	set(handles.figure1,'pointer','watch'),			pause(0.05)
	[ny,nx] = size(Z);		D2R = pi/180;			deg2m = 111194.9;
	z_min = head(5);		z_max = head(6) + 1;

	if (handles.geog)
		p_size = sqrt((head(8) * deg2m) * (head(9) * deg2m * cos ((head(4) + head(3))*0.5 * D2R))); 
	else
		p_size = head(8) * head(9);
	end

	sh = grdgradient_m(Z, head, sprintf('-Em%.2f/%.2f', 0, 20));
	sh = uint8(cvlib_mex('CvtScale', sh, 254, 1));

	str_amp = 2;
	alpha = tan(25 * D2R) * str_amp / p_size;
	decal = fix(2 * alpha * (z_max - z_min));
	decal = max(8, decal + rem(decal,2));			% Make it always even and with a minimum of 8

	left = repmat(uint8(255), 1, nx + decal);	right = left;
	ar = repmat(uint8(0), ny, nx);				ag = ar;	l = 0;	r = l;
	inner_ind = (decal / 2) + 1 : nx + (decal / 2);	% Inner indices that trow away the left-right expansion
	for (i = 1:ny)
		for (j = 1:nx)
			iz = fix(alpha * (double(Z(i,j)) - z_min));
			if (j == 1)
				left(j+iz) = sh(i,j);
				right(decal+j-iz) = sh(i,j);
			else
				for (k = r:decal + j - iz),	right(k) = sh(i,j);		end		% Much faster than vectorized code
				for (k = l:j + iz),			left(k)  = sh(i,j);		end
			end
			l = j + iz;			r = decal + j-iz;
		end
		ar(i,1:end) = left(inner_ind);		ag(i,1:end) = right(inner_ind);
		left(:) = uint8(0);					right(:) = uint8(0);
	end
	zz = cat(3,ar,ag,ag);
	set(handles.figure1,'pointer','arrow')
	tmp.head = head;		tmp.geog = handles.geog;	tmp.name = 'Anaglyph';
	tmp.X = head(1:2);		tmp.Y = head(3:4);
	mirone(zz,tmp);

% --------------------------------------------------------------------
function img = shade_manip_raster(azimuth, elevation, Z, size_amp)
% Almost deprecated function.
	if (nargin == 3),	size_amp = 125;		end

	u1 = sin(azimuth) * cos(elevation);		u2 = -cos(azimuth) * cos(elevation);
	u3 = sin(elevation);
	if (~isa(Z,'double')),	Z = double(Z);	end		% Make sure Z is of double type

	% Derivatives of function with respect to rows and columns
	dZdc = zeros(size(Z));		dZdr = zeros(size(Z));

	% Take forward differences on left and right edges
	dZdr(1,1:end)	= (Z(2,1:end) - Z(1,1:end)) / size_amp;
	dZdr(end,1:end)	= (Z(end,1:end) - Z(end-1,1:end)) / size_amp;
	dZdc(:,1)		= (Z(:,2) - Z(:,1)) / size_amp;
	dZdc(:,end)		= (Z(:,end) - Z(:,end-1)) / size_amp;

	% Take centered differences on interior points
	dZdr(2:end-1,1:end) = (Z(3:end,1:end)-Z(1:end-2,1:end)) / size_amp;
	dZdc(:,2:end-1)     = (Z(:,3:end)-Z(:,1:end-2)) / size_amp;
	img = (dZdr*u1 + dZdc*u2 + 2*u3) ./ (sqrt(dZdr .^ 2 + dZdc .^ 2 + 4));
	img(img < 0) = 0;	img(img > 1) = 1;

% --------------------------------------------------------------------
function ImageLink_CB(handles, opt)
	if (handles.no_file),	return,		end
	hFigs = findobj(0,'type','figure');						% Fish all figures
	if (numel(hFigs) == 1),	return,		end					% No one else arround
	ind = (hFigs - handles.figure1) == 0;
	hFigs(ind) = [];										% Remove current figure from the fished list
	IAmAMir = true(1, numel(hFigs));
	for (k = 1:numel(hFigs))
		if (isempty(getappdata(hFigs(k), 'IAmAMirone'))),	IAmAMir(k) = false;		end
	end
	hFigs = hFigs(IAmAMir);									% Retain only the Mirone figures
	if (isempty(hFigs)),		return,		end
	nomes = get(hFigs,'name');
	if (~isa(nomes,'cell')),	nomes = {nomes};	end
	for (k = 1:numel(hFigs))
		[pato, nomes{k}] = fileparts(nomes{k});
	end
	s = listdlg('PromptString','Link with', 'SelectionMode','single', 'ListString',nomes);
	if (isempty(s)),	return,		end
	linkedFig = hFigs(s);
	setappdata(handles.figure1,'LinkedTo',linkedFig)		% pixval_stsbar will take over from here
	setappdata(linkedFig,'LinkedTo',handles.figure1)

% --------------------------------------------------------------------
function ImageRetroShade_CB(handles)
% Illuminate current figure with the reflectance calculated from parent figure
	if (handles.no_file),	return,		end
	h_f = getappdata(handles.figure1,'hFigParent');		% Get the parent figure handle
	if (~ishandle(h_f))
		msgbox('Parent window no longer exists (you kiled it). Exiting.','Warning');		return
	end
	handParent = guidata(h_f);				% Get the parent handles
	[rect, rect_crop] = aux_funs('rectangle_and', handles.head, handParent.head);	% Get the intersection zone
	if (isempty(rect)),	warndlg('The two images do not overlap.','Warning'),	return,		end
	R = ImageIllumModel_CB(handParent);	% Get parent Reflectance array
	if (isempty(R)),	warndlg('Illumination models >= 5 don''t work here.','Warning'),	return,		end
	[R, r_c.r_c] = cropimg(handParent.head(1:2),handParent.head(3:4),R,rect_crop,'out_grid');

	% If parent image is of different resolution, resize it to fit the son_image resolution
	if ( abs(handParent.head(8) - handles.head(8)) > 1e-8 || abs(handParent.head(9) - handles.head(9)) > 1e-8 )
		nx_parent = round(rect_crop(3) / handles.head(8)) + 1;
		ny_parent = round(rect_crop(4) / handles.head(9)) + 1;
		R = cvlib_mex('resize',R,[ny_parent nx_parent],'bicubic');
	end
	
	img = get(handles.hImg,'CData');
	if (ndims(img) == 2),		img = ind2rgb8(img,get(handles.figure1,'Colormap'));	end
	m = size(img,1);			n = size(img,2);
	if (~isequal(size(R), [m n]))			% Happens when parent doesn't completely overlaps son image
		[img2, r_c] = cropimg(handles.head(1:2),handles.head(3:4),img,rect_crop,'out_grid');
		img2 = shading_mat(img2,R,'no_scale');
		img(r_c(1):r_c(2), r_c(3):r_c(4), 1:end) = img2;
	else
		img = shading_mat(img,R,'no_scale');
	end
	
	set(handles.hImg,'CData',img)
	aux_funs('togCheck',handles.ImModRGB, [handles.ImMod8cor handles.ImMod8gray handles.ImModBW])
	set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function ImageDrape_CB(handles, alfa)
	if (handles.no_file),	return,		end
	if (nargin == 1),	alfa = [];		end				% Transparency will be asked later
	son_img = get(handles.hImg,'CData');				% Get "son" image
	h_f = getappdata(handles.figure1,'hFigParent');		% Get the parent figure handle
	if (~ishandle(h_f))
		msgbox('Parent window no longer exists (you kiled it). Exiting.','Warning');		return
	end
	handParent = guidata(h_f);			% We need the parent handles
	parent_img = get(handParent.hImg,'CData');			parent_was_resized = false;
	y_son = size(son_img,1);			x_son = size(son_img,2);			% Get "son" image dimensions 
	y_parent = size(parent_img,1);		x_parent = size(parent_img,2);		% Get "parent" image dimensions 
	transp = [];

	% Find if image needs to be ud fliped
	if(strcmp(get(handles.axes1,'YDir'),'reverse')),	son_img = flipdim(son_img,1);	end

	% See if son_img comes from a grid with NaNs and if yes make the NaNs transparent
	if (handles.validGrid && handles.have_nans)
		[X,Y,Z] = load_grd(handles);		transp = isnan(Z);
	end

	% See about transparency
	if (isempty(alfa))			% transparency was NOT transmitted as argument to this function
		dlg_title = 'Draping Transparency';		num_lines= [1 38];	defAns = {'0'};
		resp = inputdlg('Use Transparency (0-1)?',dlg_title,num_lines,defAns);		pause(0.01);
		alfa = sscanf(resp{1},'%f');
	end
	if (isempty(alfa) || alfa > 1),		alfa = 1;	end
	if (alfa < 0.01 || isnan(alfa)),	alfa = 0;	end

	son_cm = [];
	if (ndims(son_img) == 2),		son_cm = get(handles.figure1,'Colormap');	end		% Get "son" color map

	blind_drape = false;
	if (handles.image_type == 2 || handles.image_type == 20 || handParent.image_type == 2 || handParent.image_type == 20)
		blind_drape = true;
	end
	if (blind_drape)		% Drape based solely on images sizes
		if (y_son ~= y_parent || x_son ~= x_parent)				% Check if "son" and "parent" images have the same size
			son_img = cvlib_mex('resize',son_img,[y_parent x_parent],'bicubic');
		end
	else					% Drape based on images coords - First find the intersection of the 2 regions 
		P1.x = [handles.head(1) handles.head(1) handles.head(2) handles.head(2) handles.head(1)];	P1.hole = 0;
		P1.y = [handles.head(3) handles.head(4) handles.head(4) handles.head(3) handles.head(3)];
		P2.x = [handParent.head(1) handParent.head(1) handParent.head(2) handParent.head(2) handParent.head(1)];	P2.hole = 0;
		P2.y = [handParent.head(3) handParent.head(4) handParent.head(4) handParent.head(3) handParent.head(3)];
		P3 = PolygonClip(P1, P2, 1);				% Intersection of the two rectangles
		if (isempty(P3)),	warndlg('The two images do not overlap.','Warning'),	return,		end
		rx_min = min(P3.x);			rx_max = max(P3.x);		ry_min = min(P3.y);			ry_max = max(P3.y);
		rect_crop = [rx_min ry_min rx_max-rx_min ry_max-ry_min];

		% If parent image is of lesser resolution (90% lower) resize it to fit the son_image resolution
		if ( (handParent.head(8) > handles.head(8)*1.1) || (handParent.head(9) > handles.head(9)*1.1) )
			x_parent = round(diff(handParent.head(1:2)) / handles.head(8)) + 1;
			y_parent = round(diff(handParent.head(3:4)) / handles.head(9)) + 1;
			parent_img = cvlib_mex('resize',parent_img,[y_parent x_parent],'bicubic');
			parent_was_resized = true;
		end
		[r_c] = cropimg(handParent.head(1:2),handParent.head(3:4),parent_img,rect_crop,'out_ind');
		if (diff(r_c(1:2)) <= 0 || diff(r_c(3:4)) <= 0),	return,		end
		[I,P1.hole] = cropimg(handles.head(1:2),handles.head(3:4),son_img,rect_crop,'out_grid');	% P1.hole to shut up MLint
		son_img = cvlib_mex('resize',I,[diff(r_c(1:2)) diff(r_c(3:4))]+1,'bicubic');

		if (~isempty(transp))		% Have NaNs, make them transparent but some complications arise from interpolation
			[I,P1.hole] = cropimg(handles.head(1:2), handles.head(3:4), transp, rect_crop, 'out_grid');
			transp = cvlib_mex('resize',I,[diff(r_c(1:2)) diff(r_c(3:4))]+1,'bilinear');
			if (ndims(son_img) == 2)
				transp = transp | (son_img == 0);		% Combine color & isnan() interpolations, which are not equal (???)
			else
				transp = transp | ((son_img(:,:,1) == 255) & (son_img(:,:,2) == 255) & (son_img(:,:,3) == 255)); % Not sure why I do this
			end
		end

		% Make sure parent & son are both indexed or true color
		if (ndims(son_img) == 2 && ndims(parent_img) == 3)
			son_img = ind2rgb8(son_img,get(handles.figure1,'Colormap'));
		elseif (ndims(son_img) == 3 && ndims(parent_img) == 2)
			parent_img = ind2rgb8(parent_img,get(h_f,'Colormap'));
		end
		if (alfa)		% We must process this case here
			tmp = parent_img(r_c(1):r_c(2),r_c(3):r_c(4), 1:end);
			cvlib_mex('addweighted',son_img,(1 - alfa),tmp,alfa)	% In-place
			alfa = 0;	% We don't want to repeat the alpha blending below
		end
		if (~isempty(transp))
			for (k = 1:size(son_img,3))
				tmp = son_img(:,:,k);			tmp2 = parent_img(r_c(1):r_c(2),r_c(3):r_c(4),k);
				tmp(transp) = tmp2(transp);		son_img(:,:,k) = tmp;
			end
		end
		parent_img(r_c(1):r_c(2),r_c(3):r_c(4), 1:end) = son_img;
		son_img = parent_img;		% Do this to be compatible with the blind_drape mode
	end

	if (alfa)
		if (ndims(son_img) == 2),		son_img = ind2rgb8(son_img,get(handles.figure1,'Colormap'));	end
		if (ndims(parent_img) == 2),	parent_img = ind2rgb8(parent_img,get(h_f,'Colormap'));			end
		cvlib_mex('addweighted',son_img,(1 - alfa),parent_img,alfa)		% In-place
	end
	if (~parent_was_resized)
		set(handParent.hImg,'CData',son_img);
	else
		imSize = [];
		if ( (handles.image_type ~= 2 && handles.image_type ~= 20) && (abs(diff(handles.head(8:9))) > 1e-4) )	% Check aniso
			imSize = handles.head(8) / handles.head(9);		% resizetrue will know what to do with this
		end
		handParent.hImg = image(handParent.head(1:2),handParent.head(3:4),son_img,'Parent',handParent.axes1);
		resizetrue(handParent,imSize,'xy');
	end

	if (~isempty(son_cm) && (alfa > 0)),	set(h_f,'Colormap',son_cm);		end		% Set "son" colormap to "parent" figure
	% Signal in the parent image handles that it has a draped image
	handParent.is_draped = true;		handParent.Illumin_type = 0;	guidata(h_f,handParent)

	% Store in Parent image the info about this (son) fig handle and transparency. For use in illuminations.
	setappdata(handParent.hImg,'drapInfo',struct('hFig',handles.figure1,'alpha',alfa))

% --------------------------------------------------------------------
function ToolsMeasure_CB(handles, opt)
% OPT = 'LLength' | 'Azim' | 'Area'
	if (handles.no_file),	return,		end
	zoom_state(handles,'maybe_off');
	opt2 = opt;
	if (strcmp(opt,'Area')),	opt2 = 'closed';	end
	[xp,yp] = getline_j(handles.figure1, opt2);		% When opt ~= closed, 2th arg is ignored
	if (numel(xp) < 2),		zoom_state(handles,'maybe_on'),		return,		end
	draw_funs([xp(:) yp(:)],['tell' opt])
	zoom_state(handles,'maybe_on');

% --------------------------------------------------------------------
function ToolsMeasureAreaPerCor_CB(handles)
	if (handles.no_file),	return,		end
	img = get(handles.hImg,'CData');		pal = get(handles.figure1,'ColorMap');
	n = size(pal,1);	ny = size(img,1);	nx = size(img,2);		D2R = pi / 180;
	if (ndims(img) == 3)	% RGB, reject unless it as a "gray RGB"
		ind_m = round(rand(1,5)*(ny-1)+1);
		ind_n = round(rand(1,5)*(nx-1)+1);
		df = diff(double(img(ind_m,ind_n, 1:end)),1,3);
		if (~any(df(:) ~= 0))			% Yeap, a RGB gray
			img(:,:,2:3) = [];			% Remove the non-needed pages
		else
			warndlg('This does not work with true color images (they are just too many)','Warning'),	return
		end
		pal = unique(img(:));
		pal = (0:double(pal(end)))' / 255;		% Pretend this is a colormap
		n = numel(pal);
	end

	if (handles.geog)
		y3L = 1:3:ny;	y3L(end+1) = ny + 1;	% Vector with one every other third line
		counts = zeros(n,1);
		X = linspace(handles.head(1),handles.head(2),nx);
		Y = linspace(handles.head(3),handles.head(4),ny);
		Y(end+1) = Y(end);		% Some cases need one extra value in Y(y3L(i)+1) below
		DX = (X(2) - X(1)) * 111.1949;		DY = (Y(2) - Y(1)) * 111.1949;		% Convert DX|Y to km (at the equator)
		for (i=1:numel(y3L)-1)
			threeL = img(y3L(i):y3L(i+1)-1,1:end);	% Chop 3 lines from image
			counts = counts + imhistc(threeL, n, 0, n) * DY * DX * cos(Y(y3L(i)+1)*D2R);
		end
		ColNames = {'Intensity' 'Area km^2'};
	else
		counts = imhistc(img, n, 0, 255);			% Call MEX file to compute the histogram
		ColNames = {'Intensity' 'N pixels'};
	end
	ind = (counts == 0);
	counts(ind) = [];		pal(ind,1:end) = [];		% Remove zero counts
	tableGUI('array',[round(pal(:,1)*255) counts],'ColWidth',[60 90],'ColNames',ColNames,...
		'FigName','Area measures','RowNumbers','y','MAX_ROWS',20,'modal','');
	double2ascii([handles.path_tmp 'area_per_color.dat'],[round(pal(:,1)*255) counts],'%d\t%g')		% Save the result

% --------------------------------------------------------------------
function DrawLine_CB(handles, opt)
	if (handles.no_file),	return,		end
	if (nargin == 1),		opt = [];	end
	% The following is a necessary patch against two big stupidities.
	% First one is from the user that double-clicked on the line icon
	% Second is from Matlab that doesn't let us test a double-click on a uipushtool
	if ( ~isempty(getappdata(handles.figure1, 'fromGL')) ),		return,		end

	if (strncmp(opt,'free',4)),				[xp,yp] = getline_j(handles.figure1,'freehand');
	elseif (strcmp(opt,'GCPmemory')),		xp = [0 0];			% Jump the manual drawing
	elseif (strcmp(opt,'GCPimport')),		xp = [0 0];			% Jump the manual drawing
	elseif (strcmp(opt,'spline')),			[xp,yp] = getline_j(handles.figure1,'spline');
	else									[xp,yp] = getline_j(handles.figure1);
	end
	n_nodes = numel(xp);					LS = 'none';
	if (n_nodes < 2),	zoom_state(handles,'maybe_on'),		return,		end
	% The polyline Tag is very important to destinguish from MB tracks, which have Tags = MBtrack#
	if (strcmp(opt,'FaultTrace'))		% When this function is used for Okada modeling
		h = line('XData', xp, 'YData', yp,'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,'Tag','FaultTrace');
		% Create empty patches that will contain the surface projection of the fault plane
		hp = zeros(1,numel(xp)-1);
		for (k=1:length(xp)-1),		hp(k) = patch('XData',[], 'YData',[]);		end
		setappdata(h,'PatchHand',hp);
	elseif (strncmp(opt,'GCP',3))
		if (strcmp(opt,'GCPmemory'))
			GCPinMemory = getappdata(handles.figure1,'GCPregImage');
			xp = GCPinMemory(:,1);		yp = GCPinMemory(:,2);	LS = ':';
		elseif (strcmp(opt,'GCPimport'))		% Load GCPs file
			str1 = {'*.dat;*.DAT', 'Data file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
			[FileName,PathName] = put_or_get_file(handles,str1,'Select input GCP file name','get');
			if (isequal(FileName,0)),		return,		end
			xy = load_xyz(handles,[PathName FileName]);
			xp = xy(:,1);				yp = xy(:,2);
			setappdata(handles.figure1,'GCPregImage',xy)
			setappdata(handles.figure1,'fnameGCP',FileName)	% To know when GCPs must be removed from appdata (in show_image)
		end			% else -> opt = 'GCPpline'
		h = line('XData', xp, 'YData', yp,'Color','k','LineWidth',0.5,'LineStyle',LS,'Marker','o',...
			'MarkerFaceColor','y','MarkerSize',4,'Tag','GCPpolyline');
		if (opt(4) == 'm'),		register_img(handles,h,GCPinMemory(:,1:4))		% Set uicontext for img registration
		else					register_img(handles,h)
		end
		return
	elseif (strcmp(opt,'SeismicLine'))
		h = line('XData', xp, 'YData', yp,'Color','y','LineWidth',5,'Tag','SeismicLine');
	else
		if (xp(1) == xp(end) && yp(1) == yp(end))	% If line was close by hiting 'c'
			h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',handles.DefLineColor,...
				'LineWidth',handles.DefLineThick,'Tag','Closedpolygon');
		else
			h = line('XData', xp, 'YData', yp,'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,'Tag','polyline');
		end
	end
	draw_funs(h,'line_uicontext')		% Set lines's uicontextmenu

% --------------------------------------------------------------------
function hand = Draw_CB(handles, tipo, smb)
% Draw text, arrows & symbols
	if (handles.no_file),	return,		end
	if (strcmp(tipo, 'Vector'))
		draw_funs([],'DrawVector')		% Vectors are ploted there
	else
		pt = click_e_point(1,'crosshair');
		if (isempty(pt))
			if (nargout),	hand = [];	end
			zoom_state(handles,'maybe_on');		return
		end
		if (strcmp(tipo, 'Symbol'))
			h = line(pt(1,1),pt(1,2),'Marker',smb,'MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10,'Tag','Symbol');
		elseif (strcmp(tipo, 'Text'))
			h = text(pt(1),pt(2),0,'','Editing','on','VerticalAlignment','baseline','Margin',1);
		end
		draw_funs(h,['Draw' tipo])		% Set uicontextmenu
		if (nargout),	hand = h;	end
	end

% --------------------------------------------------------------------
function DrawClosedPolygon_CB(handles, opt)
% OPT = [] Draw a closed polygon and exit
% OPT = 'EulerTrapezium' Draw a trapezium (only 4 sides => 5 pts). Used to fast compute an Euler pole
% OPT = 'SeismicPolyg' Draw a closed polygon to which special seismicity opts will be added.
% OPT = 'from_ROI' Draw a closed polygon and call the ROI operations window
% OPT = 'rectangle' Draw a rectangle
% OPT = h (where h is a line handle) Calls the ROI operations window to operate on the polygon
%		whose handle is h. This mode is activated from a uicontextmenu that closed polylines share
	if (handles.no_file),	return,		end
	if (nargin == 1),		opt = [];	end
	if (ishandle(opt)),		h_line = opt;	opt = 'from_uicontext';  end
	% The following is a necessary pach agains two big stupidities.
	% First is from the user that double-clicked on the polyline icon
	% Second is from Matlab that doesn't lets us test a double-click on a uipushtool
	if (strcmp(get(handles.figure1,'Pointer'), 'crosshair')),		return,		end

	if ( isempty(opt) || any(strcmp(opt,{'from_ROI' 'EulerTrapezium' 'SeismicPolyg'})) )
		zoom_state(handles,'maybe_off');
		[xp,yp] = getline_j(handles.figure1,'closed');		n_nodes = length(xp);
		if (n_nodes < 4),	return,		end			% 4 because a straight line has 3 vertex (last one repeats)
		if (strcmp(opt,'EulerTrapezium') && n_nodes ~= 5)
			errordlg('OK, I won''t insult you this time. Just RTFM or don''t use this option','Error'),		return
		end
		if (strcmp(opt,'EulerTrapezium')),		tag = opt;
		elseif (strcmp(opt,'SeismicPolyg')),	tag = opt;
		else									tag = 'Closedpolygon';
		end
		zoom_state(handles,'maybe_on');
		xp = xp(:)';	yp = yp(:)';
		% The following Tag is very important to distinguish from MB tracks, which have Tags = MBtrack#
		h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',handles.DefLineColor,...
			'LineWidth',handles.DefLineThick,'Tag',tag);
		if (handles.version7 >= 8.4),	set(h, 'FaceColor', 'w', 'FaceAlpha', 0.005),	end		% F TMW, always breaking things
		draw_funs(h,'line_uicontext')		% Set lines's uicontextmenu
		% We are done in this mode (just draw a closed polygon)
		if (isempty(opt) || any(strcmp(tag,{'EulerTrapezium' 'SeismicPolyg'}))),	return,		end
		% Come here when ROI operations were chosen
		roi_image_operations(handles.axes1,[xp(:),yp(:)])
	elseif (strcmp(opt,'from_uicontext'))
		xp = get(h_line,'XData');	yp = get(h_line,'YData');
		roi_image_operations(handles.axes1,[xp(:),yp(:)])
	elseif (strcmp(opt,'rectangle'))
		zoom_state(handles,'maybe_off');
		[p1,p2,hl] = rubberbandbox(handles.axes1);
		zoom_state(handles,'maybe_on');
		difa = abs(p2 - p1);
		if ((difa(1) < handles.head(7)/4) || (difa(2) < handles.head(8)/4))
			delete(hl),		return			% Don't draw ultra small rectangles
		end
		set(hl,'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick)	% Use defaults LineThick and DefLineColor
		am_I_padded = getappdata(handles.figure1,'SidePadded');
		if (~isempty(am_I_padded))			% TINTOL mode
			draw_funs([],'set_recTsu_uicontext', hl)	% Set uicontextmenus appropriate for grid nesting
		else
			draw_funs(hl,'line_uicontext')	% Set lines's uicontextmenu
		end
	end

% --------------------------------------------------------------------
function DrawEulerPoleCircle_CB(handles)
	if (aux_funs('msg_dlg',1,handles)),		return,		end		% Test geog & no_file
	if (strcmp(get(handles.figure1,'Pointer'), 'crosshair')),	return,		end		% Already drawing something else
	zoom_state(handles,'maybe_off');
	
	out = euler_poles_selector(handles.home_dir);			% The output is a struct with fields: lon lat omega plates model
	if isempty(out),	return,		end % User gave up
	plon = out.lon;		plat = out.lat;
	if (~out.absolute)
		h_circ = uicirclegeo(plon, plat, handles.axes1);
		set(h_circ,'Tag','CircleEuler')		% This is used by draw_funs to allow velocity computations
		s = get(h_circ,'Userdata');
		s.omega = out.omega;
		if ~isempty(out.plates)			% Just in case
			if (~strncmp('absolute',out.plates,8))			% A relative plate model
				s.plates = [out.plates '  -- Model = ' out.model];
			else											% An absolute plate model
				s.plates = [out.plates(end-1:end) ' -- Model = ' out.model ' (Absolute)'];
			end
		else
			s.plates = 'I''m lost';
		end
		set(h_circ,'Userdata',s)
		draw_funs(h_circ,'SessionRestoreCircle')
	else
		pt = click_e_point(1,'crosshair');
		if (isempty(pt)),	return,		end
		lon = pt(1);		lat = pt(2);
		[Vx, Vy] = draw_funs([], 'compute_EulerVel', lat, lon, plat, plon, out.omega, 'Nikles');
		struc_arrow = struct('spacingChanged',[], 'hQuiver', [], 'hAx', handles.axes1);
		hQuiver = draw_funs([], 'loc_quiver', struc_arrow, lon, lat, Vx, Vy);
		x = get(hQuiver,'XData');		y = get(hQuiver,'YData');
		set(hQuiver(1),'XData',[x{1} x{2}], 'YData',[y{1} y{2}])	% Merge the header with the "trunk"
		delete(hQuiver(2));		hQuiver(2) = [];
		set(hQuiver,'Tag','Seta','Userdata',1)
		draw_funs(hQuiver,'line_uicontext')
	end
	zoom_state(handles,'maybe_on');

% --------------------------------------------------------------------
function DrawGeogCircle_CB(handles, opt)
	if (handles.no_file),	return,		end
	if (strcmp(get(handles.figure1,'Pointer'), 'crosshair')),	return,		end		% Already drawing something else
	if (nargin == 1),	opt = [];	end
	zoom_state(handles,'maybe_off');
	if (handles.geog && strcmp(opt,'gcirc'))
		draw_funs([],'DrawGreatCircle')				% All work is done there
	elseif (~handles.geog && strcmp(opt,'gcirc'))
		warndlg('Great Circles are only programed to work with geog coordinates.','Warning')
	elseif (~handles.geog && isempty(opt) || strcmp(opt,'cartCirc'))
		draw_funs([],'DrawCartesianCircle')			% All work is done there
	else
		h_circ = uicirclegeo;						% Draw the circle and associate Button... functions
		draw_funs(h_circ,'SessionRestoreCircle')	% Give uicontext
	end
	zoom_state(handles,'maybe_on');

% --------------------------------------------------------------------
function DrawImportText_CB(handles)
% Read a file with text to plot.
% Minimum content per line is "X Y TEXT", but can also have an option -F like in GMT5 before the TEXT
% Accepted forms of -F are
%	-F[+a<angle>][+j<T|M|B|L|C|R>][+f<size>,<font>,<color>], where +f can be still be incomplete as in:
%			+f<size> || +f<size>,<font> || +f<size>,=,<color>
%			<size> is the font size in points. It can have a 'p' appended as in 12p but that's ignored
%			<font> is a valid Matlab fontName
%			<color> is a 'r/g/b' [0 255] color string

	if (handles.no_file),	return,		end
	str1 = {'*.txt;*.TXT;*.dat;*.DAT', 'Text file (*.txt,*.TXT,*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select input text file name','get');
	if isequal(FileName,0),		return,		end

	fid = fopen([PathName,FileName],'r');
	todos = fread(fid,'*char');		fclose(fid);
	todos = strread(todos','%s','delimiter','\n');
	n_lines = numel(todos);
	str.x(n_lines) = 0;		str.y(n_lines) = 0;		str.name = cell(n_lines,1);		foul = false(n_lines,1);
	for (k = 1:n_lines)
		try
			txt = strrep(todos{k},char(9),' ');		% Get rid of eventual tab characters which print as squares
			[t, r] = strtok(txt);			str.x(k) = str2double(t);
			[t, r] = strtok(r);				str.y(k) = str2double(t);
			str.name{k} = r(2:end);
		catch
			foul(k) = true;
		end
	end
	str.x(foul) = [];			str.y(foul) = [];		str.name(foul) = [];
	if (numel(str.x) == 0),		return,		end			% Only errors. Nothing left

	limits = getappdata(handles.axes1,'ThisImageLims');
	if (handles.is_projected && handles.defCoordsIn > 0)
		tmp = proj2proj_pts(handles,[str.x str.y],'lim',limits,'srcProj4','+proj=longlat');
		str.x = tmp(:,1);		str.y = tmp(:,2);
	end
	% Get rid of Texts that are outside the map limits
	indx = (str.x < limits(1) | str.x > limits(2));		str.x(indx) = [];		str.y(indx) = [];	str.name(indx) = [];
	indx = (str.y < limits(3) | str.y > limits(4));		str.x(indx) = [];		str.y(indx) = [];	str.name(indx) = [];
	if (numel(str.x) == 0),		return,		end			% No texts inside area. Return.

	for (k = 1:numel(str.x))		% Loop over number of text strings
		ang = 0;	fz = 10;		% Default font size and rotation
		fontName = 'Book Antiqua';
		fontColor = 'k';
		HorzAlng = 'center';	VertAlng = 'baseline';
		nc = numel(str.name{k});

		if (nc > 2 && strcmp(str.name{k}(1:2), '-F'))	% A GMT5 type -F option like in pstext
			[t, r] = strtok(str.name{k});
			if (isempty(r)),	r = 'Empty text is nice, isn''t it?';	end		% There was nothing to plot
			str.name{k} = ddewhite(r);
			t = t(3:end);						% Strip the '-F'
			ind = strfind(t, '+');
			ind(end+1) = numel(t) + 1;			% To easy up the algo below
			for (n = 1:numel(ind)-1)
				tok = t(ind(n)+2:ind(n+1)-1);
				switch t(ind(n)+1)
					case 'a'
						ang = str2double(tok);
						if (isnan(ang)),	ang = 0;	end
					case 'j'
						switch upper(tok)
							case {'TL' 'LT'},	HorzAlng = 'left';		VertAlng = 'cap';
							case {'ML' 'LM'},	HorzAlng = 'left';		VertAlng = 'middle';
							case {'BL' 'LB'},	HorzAlng = 'left';		VertAlng = 'baseline';
							case {'TC' 'CT'},	HorzAlng = 'center';	VertAlng = 'cap';
							case {'MC' 'CM'},	HorzAlng = 'center';	VertAlng = 'middle';
							case {'BC' 'CB'},	HorzAlng = 'center';	VertAlng = 'baseline';
							case {'TR' 'RT'},	HorzAlng = 'right';		VertAlng = 'cap';
							case {'MR' 'RM'},	HorzAlng = 'right';		VertAlng = 'middle';
							case {'BR' 'RB'},	HorzAlng = 'right';		VertAlng = 'baseline';
						end
					case 'f'		% Parse for a complete +f<size>,<font>,<color> OR partial +f<size> || +f<size>,<font>
						i = strfind(tok, ',');
						if (isempty(i)),	i = numel(tok) + 1;		end		% Pretend we had one comma
						fz = str2double(tok(1:i(1)-1));				% try a unitless size
						if (isnan(fz)),		fz = str2double(tok(1:i(1)-2));		end		% Try removing last char as in '12p'
						if (isnan(fz)),		fz = 10;	end			% If still fails, fall back to default
						if (numel(i) == 1 && numel(tok) > i)		% -F+f<size>,fontName
							fontName = tok(i(1)+1:end);
						elseif (numel(i) == 2 && numel(tok) > i(2))	% -F+f<size>,fontName,color
							fontName = tok(i(1)+1:i(2)-1);
						else
							continue
						end
						if (fontName == '='),	fontName = 'Book Antiqua';	end
						if (numel(i) == 2)							% -F+fsize,fontName,fontColor
							tok = tok(i(2)+1:end);					% Reuse this var. It's the last token
							is = strfind(tok,'/');
							if (numel(is) ~= 2),	continue,	end		% Not a r/g/b color
							fontColor = reshape(sscanf(tok,'%d/%d/%d') / 255, 1, 3);
						end
				end
			end
		end
		if (str.name{k}(1) == '_')			% If leading field is of the form _numeric, take it as font size
			[t, r] = strtok(str.name{k});		fz = str2double(t(2:end));		str.name{k} = r;
			if (isnan(fz)),		fz = 10;		end		% Bad luck or stupid user error
		end
		h = text(str.x(k),str.y(k),str.name{k},'FontSize', fz, 'FontName',fontName,'FontWeight','bold', ...
			'Rotation',ang, 'Color',fontColor, 'HorizontalAlignment',HorzAlng, 'VerticalAlignment',VertAlng, 'Margin',1);
		draw_funs(h,'DrawText')
	end

% --------------------------------------------------------------------
function DrawImportOGR_CB(handles, fname)
	if (nargin == 1)
		str1 = {'*.shp;*.gpx;*.kml;*.gml;*.dxf', 'Data files (*.shp,*.gpx;*.kml;*.gml;*.dxf)';'*.*', 'All Files (*.*)'};
		[FileName,PathName] = put_or_get_file(handles,str1,'Select file name','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	end
	try			s = ogrread(fname);
	catch,		errordlg([lasterr ' Bad luck. NOT an OGR readable file'],'Error'),		return
	end

	do_project = false;		no_file = handles.no_file;
	theProj = s(1).SRSProj4;

	region = s(1).BoundingBox(:)';
	is_geog = aux_funs('guessGeog',region(1:4));

	% Test if file was correctly read and in case not, try with check_cmop() (profile data)
	if (any(isinf(region)))
		err = check_cmop(fname);
		if (err)
			warndlg(['Failed reading file ' fname ' with OGR'],'Warning')
		elseif (handles.no_file)
			delete(handles.figure1)		% Data is displayed in Ecran and the Mirone fig is of no use
		end
		return
	end

	% If we have a file already need to know about coords compat
	if (~handles.no_file)
		if (isempty(theProj) && handles.geog && ~is_geog)
			errordlg('Your background image is in geographics but the shape file has unknown coords.','ERROR'),	return	
		end
		if (~isempty(theProj))
			if (handles.geog && ~is_geog)
				projStruc.SrcProjWKT = theProj;		do_project = true;		% Inverse projection
			else						% Basicaly a projection from whatever to whatever
				projStruc.SrcProjWKT = theProj;
				prjInfoStruc = aux_funs('getFigProjInfo',handles);
				if (isempty(prjInfoStruc.projWKT) && ~isempty(prjInfoStruc.proj4))	% Give it one more chance
					prjInfoStruc.projWKT = ogrproj(prjInfoStruc.proj4);
				elseif (isempty(prjInfoStruc.projWKT) && ~isempty(prjInfoStruc.projGMT))
					errordlg('Sorry but it''s still not possible to convert between GMT and proj4 strings','Error'), return
				end
				projStruc.DstProjWKT = prjInfoStruc.projWKT;				% If it doesn't exist, BUM
				do_project = true;
			end
		end
		if (do_project)
			tmp = s(1).BoundingBox;
			tmp = ogrproj(tmp, projStruc);		% Reproject all BB right away
			s(1).BoundingBox = tmp;
		end
	else
		handles.geog = is_geog;		guidata(handles.figure1, handles);
	end

	% If we have nothing opened create a background region
	lt = handles.DefLineThick;		lc = handles.DefLineColor;
	if (handles.no_file),	handles = FileNewBgFrame_CB(handles, [region handles.geog]);	lc = 'k';	end

	imgLims = getappdata(handles.axes1,'ThisImageLims');
	reco = aux_funs('rectangle_and', imgLims, s(1).BoundingBox(:)');
	if (isempty(reco))
		warndlg('No data inside display region','Warning'),		return
	end

	nGeoms = numel(s);		h1 = zeros(nGeoms,1);		h2 = zeros(nGeoms,1);
	nParanoia = 1000;		% The name talks. COMPLETELY MATLAB CONDITIONED, I WAS NOT LIKE THAT BEFORE
	for (k = 1:nGeoms)
		is3D = ~isempty(s(k).Z);
		if (isempty(s(k).X)),	continue,		end		% The out struct can have many empty elements
		if (do_project),	ogrproj(s(k).X, s(k).Y, projStruc);		end		% Project into basemap coords
		if ( strncmp(s(k).Type,'Point', 5) || strncmp(s(k).Type,'Line', 4) )
			lsty = {'LineStyle', '-'};
			if (s(k).Type(1) == 'P'),	lsty = {'LineStyle', 'none', 'Marker','o', 'MarkerSize',2, 'MarkerEdgeColor','k'};	end

			h1(k) = line('Xdata',single(s(k).X),'Ydata',single(s(k).Y),'Parent',handles.axes1, ...
				'Color',lc,'LineWidth',lt,'Tag','SHPpolyline',lsty{1:end});
			if (is3D),		set(h1(k),'UserData', single(s(k).Z)),	end

		else				% Polygons
			% Make sure it knows that the Earth is round
			if (handles.geog == 1 && (s(1).BoundingBox(1) < -179.5 || s(1).BoundingBox(2) > 179.5) )
				[s(k).Y, s(k).X] = map_funs('trimwrap', s(k).Y, s(k).X, [-90 90], s(1).BoundingBox(1:2), 'wrap');
			elseif (handles.geog == 2 && (s(1).BoundingBox(1) < 0.5 || s(1).BoundingBox(2) > 359.5) )
				[s(k).Y, s(k).X] = map_funs('trimwrap', s(k).Y, s(k).X, [-90 90], s(1).BoundingBox(1:2), 'wrap');
			end
			h2(k) = patch('XData',s(k).X,'YData', s(k).Y,'FaceColor','none','EdgeColor',lc,'Parent',handles.axes1,'Tag','SHPpolygon');
			if (is3D)								% IT'S IGNORING THE EARTH-IS-ROUND? TEST
				set(h2(k), 'UserData', s(k).Z)		% Fleder can drape it (+ other eventual usages)
			end
			if (nGeoms <= nParanoia),	draw_funs(h2(k),'line_uicontext'),	end
		end
	end

	h1((h1 == 0)) = [];		h2((h2 == 0)) = [];
	if ( isempty(h1) && isempty(h2) ),	warndlg('No data inside display region','Warning'),		return,		end
	if ( ~isempty(h1) )
		draw_funs(h1,'setSHPuictx')				% Set lines's uicontextmenu
	end
	if ( ~isempty(h2) && (nGeoms > nParanoia) )	% With luck, your hardware won't choke to dead with this
		draw_funs(h2,'country_patch')			% nParanoia is an arbitrary number that practice will show dependency
	end											% mostly on hardware, for I don't beleave ML will ever behave decently.

	if (no_file && ~isempty(theProj))			% We need to finish this matter
		aux_funs('appP', handles, theProj)		% If we have a WKT proj store it
		handles = aux_funs('isProj',handles);	% Check about coordinates type
		handles = setAxesDefCoordIn(handles,1);
	end
	recentFiles(handles);						% Insert fileName into "Recent Files" & save handles

% --------------------------------------------------------------------
function DrawImportShape_CB(handles, fname)
	if (nargin == 1)
		str1 = {'*.shp;*.SHP', 'Data files (*.shp,*.SHP)';'*.*', 'All Files (*.*)'};
		[FileName,PathName] = put_or_get_file(handles,str1,'Select shape file name','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	end
	try			[s,t] = mex_shape(fname);
	catch,		errordlg([lasterr ' Quite likely, NOT a shapefile'],'Error'),		return
	end

	if (~(strncmp(t,'Arc',3) || strncmp(t,'Point',5) || strncmp(t,'Polygon',7)))
		errordlg(['Sorry. Dealing with this type of data: ' t ' is not (yet?) supported'],'WarnError'),	return
	end

	theProj = [];		do_project = false;		no_file = handles.no_file;
	[PathName, fnamePRJ] = fileparts(fname);
	if (exist([PathName filesep fnamePRJ '.prj'], 'file'))
		fid = fopen([PathName filesep fnamePRJ '.prj']);
		theProj = fread(fid,inf,'*char')';		fclose(fid);
		if (strncmp(theProj,'+proj',5)),	theProj = ogrproj(theProj);		end		% proj info was a proj4 string
	end

	region = [s(1).BoundingBox(1,1:2) s(1).BoundingBox(2,1:2)];
	is_geog = aux_funs('guessGeog',region(1:4));

	% If we have a file already need to know about coords compat
	if (~handles.no_file)
		if (isempty(theProj) && handles.geog && ~is_geog)
			errordlg('Your background image is in geographics but the shape file has unknown coords.','ERROR'),	return	
		end
		if (~isempty(theProj))
			if (handles.geog && ~is_geog)
				projStruc.SrcProjWKT = theProj;		do_project = true;		% Inverse projection
			else						% Basicaly a projection from whatever to whatever
				projStruc.SrcProjWKT = theProj;
				prjInfoStruc = aux_funs('getFigProjInfo',handles);
				if (isempty(prjInfoStruc.projWKT) && ~isempty(prjInfoStruc.proj4))	% Give it one more chance
					prjInfoStruc.projWKT = ogrproj(prjInfoStruc.proj4);
				elseif (isempty(prjInfoStruc.projWKT) && ~isempty(prjInfoStruc.projGMT))
					errordlg('Sorry but it''s still not possible to convert between GMT and proj4 strings','Error'), return
				end
				projStruc.DstProjWKT = prjInfoStruc.projWKT;				% If it doesn't exist, BUM
				do_project = true;
			end
		end
		if (do_project)
			tmp = [s.BoundingBox];		tmp = tmp(1:2,1:end)';
			tmp = ogrproj(tmp, projStruc);		% Reproject all BB right away
			for (k = 1:numel(s))
				k2 = 2*k;		k1 = k2 - 1;
				s(k).BoundingBox(1,1:2) = tmp(k1:k2,1);
				s(k).BoundingBox(2,1:2) = tmp(k1:k2,2);
			end
		end
	else
		handles.geog = is_geog;		guidata(handles.figure1, handles);
	end

	% If we have nothing opened create a background region
	lt = handles.DefLineThick;		lc = handles.DefLineColor;
	if (handles.no_file),	handles = FileNewBgFrame_CB(handles, [region handles.geog]);	lc = 'k';	end

	nPolygs = length(s);	h = zeros(nPolygs,1);
	imgLims = getappdata(handles.axes1,'ThisImageLims');
	if (strncmp(t,'Arc',3) || strncmp(t,'Point',5))
		is3D = false;		lsty = {'LineStyle', '-'};
		if (t(end) == 'Z'),	is3D = true;	end
		if (t(1) == 'P'),	lsty = {'LineStyle', 'none', 'Marker','o', 'MarkerSize',2, 'MarkerEdgeColor','k'};	end
		for (i = 1:nPolygs)
			reco = aux_funs('rectangle_and', imgLims, [s(i).BoundingBox(1,1:2) s(i).BoundingBox(2,1:2)]);
			if (~isempty(reco))			% It means the polyg BB is at least partially inside
				if (do_project),	ogrproj(s(i).X, s(i).Y, projStruc);		end		% Project into basemap coords
				s(i).X = single(s(i).X);		s(i).Y = single(s(i).Y);
				h(i) = line('Xdata',s(i).X,'Ydata', s(i).Y, 'Parent',handles.axes1,'Color',lc,'LineWidth',lt,'Tag','SHPpolyline',lsty{1:end});
			end
			if (is3D),		set(h(i),'UserData', single(s(i).Z(:)')),	end
		end
		h((h == 0)) = [];				% Those were jumped because they were completely outside map limits
		if (isempty(h)),	warndlg('No data inside display region','Warning'),		return,		end
		if (strncmp(t,'Point',5))
			setappdata(h(1), 'isPoint', true)	% Used in draw_funs/setSHPuictx to NOT set the "Join lines" uictx
		end
		draw_funs(h,'setSHPuictx')		% Set lines's uicontextmenu
	elseif (strncmp(t,'Polygon',7))
		nParanoia = 1000;				% The name talks. COMPLETELY MATLAB CONDITIONED, I WAS NOT LIKE THAT BEFORE
		for (i = 1:nPolygs)
			XMin = s(i).BoundingBox(1,1);		XMax = s(i).BoundingBox(1,2);
			reco = aux_funs('rectangle_and', imgLims, [XMin XMax s(i).BoundingBox(2,1) s(i).BoundingBox(2,2)]);
			if (~isempty(reco))				% It means the polyg BB is at least partially inside
				% Make sure it knows that the Earth is round
				if (do_project),	ogrproj(s(i).X, s(i).Y, projStruc);		end		% Project into basemap coords
				if (handles.geog == 1 && (XMin < -179.5 || XMax > 179.5) )
					[s(i).Y, s(i).X] = map_funs('trimwrap', s(i).Y, s(i).X, [-90 90], [XMin XMax], 'wrap');
				elseif (handles.geog == 2 && (XMin < 0.5 || XMax > 359.5) )
					[s(i).Y, s(i).X] = map_funs('trimwrap', s(i).Y, s(i).X, [-90 90], [XMin XMax], 'wrap');
				end
				h(i) = patch('XData',s(i).X,'YData', s(i).Y,'FaceColor','none','EdgeColor',lc,'Parent',handles.axes1,'Tag','SHPpolygon');
				if ( (numel(t) >= 8) && (t(8) == 'Z') )		% IT'S IGNORING THE EARTH-IS-ROUND? TEST
					set(h(i), 'UserData', s(i).Z(:)')		% Fleder can drape it (+ other eventual usages)
				end
				% With luck, your hardware won't choke to dead with this
				if (nPolygs <= nParanoia),	draw_funs(h(i),'line_uicontext'),	end
			end
		end
		h((h == 0)) = [];					% Those were jumped because thay were completely outside map limits
		if (isempty(h)),	warndlg('No data inside display region','Warning'),		return,		end
		if (nPolygs > nParanoia)			% nParanoia is an arbitrary number that practice will show dependency
			draw_funs(h,'country_patch')	% mostly on hardware, for I don't beleave ML will ever behave decently.
		end
	end

	if (no_file && ~isempty(theProj))			% We need to finish this matter
		aux_funs('appP', handles, theProj)		% If we have a WKT proj store it
		handles = aux_funs('isProj',handles);	% Check about coordinates type
		handles = setAxesDefCoordIn(handles,1);
	end
	recentFiles(handles);						% Insert fileName into "Recent Files" & save handles

% --------------------------------------------------------------------
function GeophysicsImportGmtFile_CB(handles, opt)
% Open a .gmt/.nc(MGD77+) file OR a list of .gmt/.nc files
	if (~handles.no_file && ~handles.geog),		aux_funs('msg_dlg',1,handles);		return,		end
	if (strncmp(opt, 'list', 4))
		if (numel(opt) == 4)			% Else OPT contains already the name of the file with list
			str1 = {'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.DAT,*.txt,*.TXT)';'*.*', 'All Files (*.*)'};
			[FileName,PathName] = put_or_get_file(handles,str1,'Select list file','get');
			if isequal(FileName,0),		return,		end
		else
			FileName = opt(6:end);		PathName = '';
		end
		fid = fopen([PathName FileName]);
		c = fread(fid,inf,'*char');		fclose(fid);
		names = strread(c,'%s','delimiter','\n');	clear c fid;
	elseif (strcmp(opt, 'single'))
		[FileName,PathName] = put_or_get_file(handles,{'*.gmt;*.GMT;*.nc;*.NC', 'gmt files (*.gmt,*.GMT,*.nc,*.NC)'},'Select gmt File','get');
		if isequal(FileName,0),		return,		end
		names = {[PathName FileName]};			% So that below the code is the same as for the list case
	else
		PathName = '';
		if (~isa(opt,'cell')),		names = {opt};		% Filename sent in input
		else						names = opt;		% A list of fnames was sent in
		end
	end

	set(handles.figure1,'Pointer','watch');
	if (handles.no_file)			% We don't have a BG map, so we have to create one
		[track, names, names_ui, vars, x_min, x_max, y_min, y_max] = aux_funs('get_mgg', names, PathName, '-Fxym','-G');
		if (isempty(names)),	return,		end			% Non existing files probably
		handles = FileNewBgFrame_CB(handles, [x_min x_max y_min y_max 1]);		pause(0.05)
		handles.no_file = 0;	handles.geog = 1;	guidata(handles.figure1, handles)
	else
		x_lim = get(handles.axes1,'XLim');		y_lim = get(handles.axes1,'YLim');
		opt_R = sprintf('-R%.6f/%.6f/%.6f/%.6f', x_lim(1), x_lim(2), y_lim(1), y_lim(2));
		[track, names, names_ui, vars] = aux_funs('get_mgg', names, PathName, '-Fxym', '-G', opt_R);
		if (isempty(names)),	return,		end			% Non existing files probably
	end

	% And finaly do the ploting
	colors = rand(numel(track),3);			% Use a random color schema
	%use_aguenta = false;					% Less than 20 tracks don't use aguentabar
% 	if (numel(track) > 20),		aguentabar(0,'title','Plotting tracks'),	use_aguenta = true;		end
	for (k = 1:numel(track))
		if (isempty(track(k).longitude)),	continue,	end		% This track is completely outside the map
		h = line(track(k).longitude,track(k).latitude, 'Parent',handles.axes1,'Linewidth',handles.DefLineThick,'Color',...
			colors(k,:),'Tag',names_ui{k},'Userdata',1);
		setappdata(h,'FullName',names{k})	% Store file name in case the uicontext wants to open it with gmtedit
		setappdata(h,'VarsName',vars(k,:))	% Store field name in case the uicontext wants to open it with gmtedit
		draw_funs(h,'gmtfile',track(k).info)
		ui_edit_polygon(h)
		%if (use_aguenta),	h = aguentabar(k/numel(track));		end
	end
	set(handles.figure1,'Pointer','arrow');
	if (~strcmp(opt, 'list'))				% Insert fileName into "Recent Files" & save handles
		[PATH, FNAME, EXT] = fileparts(names{1});
		handles.fileName = [names{1} EXT];		recentFiles(handles);		
	end

% --------------------------------------------------------------------
function DrawContours_CB(handles, opt)
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	if (nargin == 1),	opt = [];	end
	[X,Y,Z,head] = load_grd(handles);
	if isempty(Z),		return,		end		% An error message was already issued
	set(handles.figure1,'pointer','watch')
	nPtsContour = 0;
	if (isempty(opt))				% Do the automatic contouring
		c = contourc(X,Y,double(Z));
		handles.plotContourLabels = 1;
	elseif (isa(opt,'char'))		% Call the interface contouring GUI
		h_which_cont = findobj(handles.figure1,'Type','line','Tag','contour');		% See if contours were mouse deleted
		if (~isempty(h_which_cont))
			if (numel(h_which_cont) == 1)
				handles.which_cont = get(h_which_cont,'Userdata');
			else
				handles.which_cont = unique(cell2mat(get(h_which_cont,'Userdata')));
			end
		else
			handles.which_cont = [];
		end
		set(handles.figure1,'pointer','arrow'),		contouring(handles.figure1,head,handles.which_cont);
		guidata(handles.figure1, handles)
		return
	elseif (isa(opt,'double'))				% Arrive here from the interface contouring GUI
		if (~isempty(handles.which_cont))	% Do not compute/draw repeated contours
			[c,ib] = setdiff(handles.which_cont,opt);	% Find eventual countours to remove (by GUI deselection)
			for (i = 1:length(c))			% Loop over removing contours (if none, this loop has no effect)
				h = findobj(handles.axes1,'type','line','userdata',handles.which_cont(ib(i)));
				for (j = 1:length(h))		% We can easily have more than one
					labHand = getappdata(h(j),'LabelHands');
					try		delete(labHand),	end	% Delete contour's labels
				end
				delete(h)							% And delete the selected contours
			end
			handles.which_cont(ib) = [];
			guidata(handles.figure1, handles);
			[c,ia,ib] = intersect(handles.which_cont,opt(:));
			opt(ib) = [];							% Remove repeated contours
		end
		if (isempty(opt)),		set(handles.figure1,'pointer','arrow'),		return,		end  % Nothing else to do
		opt(~opt) = eps * 1e-2;				% This can of BUGs sometimes ignores the zero contour
		if (numel(opt) == 1),		opt = [opt opt];	end
		c = contourc(X,Y,double(Z),opt);
		if (isempty(c)),		set(handles.figure1,'pointer','arrow'),		return,		end
		nPtsContour = handles.nPtsContour;
	end

	limit = size(c,2);
	i = 1;		h_cont = [];	cont = [];
	while(i < limit)
		z_level = c(1,i);		npoints = c(2,i);
		if (abs(z_level) < eps),		z_level = 0;	c(1,i) = 0;		end		% Account for another ML BUG
		nexti = i+npoints+1;
		xdata = c(1,i+1:i+npoints);		ydata = c(2,i+1:i+npoints);
		if (nPtsContour && numel(xdata) < nPtsContour),		i = nexti;	continue,	end
		% Create the lines
		cu = line('XData',xdata, 'YData',ydata, 'Parent',handles.axes1, 'LineWidth',handles.DefLineThick, ...
			'Color',handles.DefLineColor, 'userdata',z_level,'Tag','contour');
		h_cont = [h_cont; cu(:)];
		cont = [cont; z_level];
		i = nexti;
	end
	handles.which_cont = unique([handles.which_cont; cont]);
	[zlev, ind] = sort(cont);		clear zlev
	h_cont = h_cont(ind);					% handles are now sorted by level
	h_label = '';							% To not hang below if not clabel (anyway, I think it's not used anymore)
	if (handles.plotContourLabels)
		h_label = clabel_j(c,h_cont);		% Label countours automatically
		set(h_label,'Tag','contour');		% The tag is used in "Delete all contours" to delete also the labels
	end
	drawnow									% Don't wait for the next (SLOW) operations to display the contours
	for (i = 1:numel(h_cont))				% Set convenient uicontexts. One for each contour
		setappdata(h_cont(i),'cont_label',get(h_cont(i),'UserData'))	% Used in write_script
		draw_funs(h_cont(i),'ContourLines',h_label)
	end
	set(handles.figure1,'pointer','arrow')
	guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function FileOpenSession_CB(handles, fname)
	if (nargin == 1)
		str1 = {'*.mat;*.MAT', 'Data files (*.mat,*.MAT)'};
		[FileName,PathName,handles] = put_or_get_file(handles,str1,'Select session file name','get');
		if isequal(FileName,0),		return,		end
	else
		FileName = fname;	PathName = [];
	end
	
	figName = [PathName FileName];		% Used at the end of this function to set up Fig's name (but can change in the middle)

	s = load([PathName FileName]);
	if (isfield(s,'markers') && isfield(s, 'FitLine'))	% Ah, this is ecran session. Send it there and stop here
		h = ecran;
		hands = guidata(h);		set(hands.figure1, 'vis', 'off')
		hands.session_name = [PathName FileName];		% Pass the mat file name in this member
		guidata(hands.figure1, hands);					% We have to save it on source because it's where is going to be read
		cb = get(hands.FileOpenSession, 'Call');
		feval(cb, hands.FileOpenSession, hands)			% Call the ecran's FileOpenSession_CB function
		delete(h)
		return
	end

	if (~isfield(s,'grd_name') || strcmpi(s.grd_name(max(numel(s.grd_name)-3,1):end),'.mat'))	% Otherwise infinite loop below
		s.grd_name = [];
	end

	tala = (~isempty(s.grd_name) && exist(s.grd_name,'file') == 2);		flagIllum = true;	% Illuminate (if it is the case)
	if (~tala && ~isempty(s.grd_name))						% Give user a 2nd chance to tell where the grid is
		[PathName, FileName, EXT] = fileparts(s.grd_name);
		resp = inputdlg({'Full name (with path) of missing grid:'},'Where is the grid?',[1 60],{['.....' filesep FileName EXT]});
		if (~isempty(resp))
			s.grd_name = resp{1};
			tala = exist(s.grd_name,'file');
			if (~tala)
				warndlg('The name provided doesn''t exist either. Give up trying to help you.','Warning')
				s.grd_name = [];			% In this case we need this as empty
			end
		else
			s.grd_name = [];				% In this case we need this as empty
		end
	end

	is_defRegion = false;
	try		is_defRegion = s.is_defRegion;		end
	if (isempty(s.grd_name) || tala == 0)
		if (is_defRegion)		% In this mode we activate the GMT custom symbols drawing tool
			handles = FileNewBgFrame_CB(handles, [-0.5 0.5 -0.5 0.5 0 1], [], 'Whatever');
			handles = guidata(handles.figure1);		% Get the updated version
			handles.geog = 0;		% To play safe
			figName = 'Draw GMT Custom symbol';		% Title bar name for this case
		else
			scrsz = get(0,'ScreenSize');		% Get screen size
			dx = s.map_limits(2) - s.map_limits(1);
			dy = s.map_limits(4) - s.map_limits(3);
			aspect = dy / dx;
			nx = round(scrsz(3)*.75);		ny = round(nx * aspect);
			if (ny > scrsz(4) - 30)
				ny = scrsz(4) - 30;			nx = round(ny / aspect);
			end
			Z = repmat(uint8(255),ny,nx);			% Create a white image
			X = [s.map_limits(1) s.map_limits(2)];	Y = [s.map_limits(3) s.map_limits(4)];
			x_inc = diff(X) / nx;					y_inc = diff(Y) / ny;
			dx2 = x_inc / 2;						dy2 = y_inc / 2;
			X = X + [dx2 -dx2];						Y = Y + [dy2 -dy2];		% Make it such that the pix-reg info = region
			handles.head = [X Y 0 255 0 x_inc y_inc];
			handles.image_type = 20;
			set(handles.figure1,'Colormap', ones( size(get(handles.figure1,'Colormap'),1), 3))
			handles = show_image(handles,'Mirone Base Map',X,Y,Z,0,'xy',1);
		end
	else
		if (is_defRegion)		% In this mode we activate the GMT custom symbols drawing tool
			handles = FileNewBgFrame_CB(handles, [-0.5 0.5 -0.5 0.5 0 1], [], 'Whatever');
			handles.geog = 0;		% To play safe
			figName = 'Draw GMT Custom symbol';		% Title bar name for this case
		end

		drv = aux_funs('findFileType', s.grd_name);
		erro = gateLoadFile(handles, drv, s.grd_name);	% Load and displays the file (or dies)
		if (erro),		return,		end					% Error message already issued
		set(handles.figure1,'Colormap',s.img_pal);
		handles = guidata(handles.figure1);				% Get the updated version
		handles.origCmap = s.img_pal;
	end

	% Have to use a try-catch because of the f. compiler bugs ("exist" won't work)
	try			s.illumComm;
	catch,		s.illumComm = '';		flagIllum = false;
	end
	try			s.illumType;
	catch
		if (numel(strfind(s.illumComm,'/')) == 5),		s.illumType = 4;	end		% Lambertian
		s.illumType = 1;		% Test only one case where this might be otherwise
	end

	if (~isempty(s.grd_name) && ~isempty(s.illumComm) && flagIllum)
		[X,Y,Z,head] = load_grd(handles,'silent');
		handles.Illumin_type = s.illumType;
		if (handles.Illumin_type == 1)
			if (handles.geog),	R = grdgradient_m(Z,head,'-M', s.illumComm,'-Nt');
			else				R = grdgradient_m(Z,head, s.illumComm,'-Nt');
			end
		else
			R = grdgradient_m(Z,head, s.illumComm,'-a1');
		end
		zz = ind2rgb8(get(handles.hImg,'CData'),get(handles.figure1,'ColorMap'));
		zz = shading_mat(zz,R,'no_scale');		set(handles.hImg,'CData',zz)		% and now it is illuminated
		setappdata(handles.figure1,'illumComm', s.illumComm)		% Save the 'illumComm' in this new fig
	end

	if (s.haveMBtrack)				% case of MB tracks
		for (i = 1:length(s.MBtrack))
			h_line = line('Xdata',s.MBtrack(i).x,'Ydata',s.MBtrack(i).y,'Parent',handles.axes1, ...
				'LineWidth',s.MBtrack(i).LineWidth, 'color',s.MBtrack(i).color,'Tag', s.MBtrack(i).tag, ...
				'LineStyle', s.MBtrack(i).LineStyle);
			setappdata(h_line,'swathRatio',s.MBtrack(i).swathRatio)
			draw_funs(h_line,'MBtrackUictx')		% Set track's uicontextmenu
		end
		for (i = 1:length(s.MBbar))	% now their's bars
			h_bar = line('Xdata',s.MBbar(i).x,'Ydata',s.MBbar(i).y,'Parent',handles.axes1,'LineWidth',s.MBbar(i).LineWidth,...
				'color',s.MBbar(i).color,'Tag',s.MBbar(i).tag,'UserData',s.MBbar(i).n_vert, 'LineStyle',s.MBbar(i).LineStyle);
			draw_funs(h_bar,'MBbarUictx')			% Set track bar's uicontextmenu
		end
		handles.hMBplot = h_line;
	end

	if (s.haveCircleGeo)				% case of Geographic circles
		for (i = 1:length(s.CircleGeo))
			h_circ = line('Xdata',s.CircleGeo(i).x,'Ydata',s.CircleGeo(i).y,'Parent',handles.axes1, ...
				'LineWidth',s.CircleGeo(i).LineWidth, 'color',s.CircleGeo(i).color,'Tag',s.CircleGeo(i).tag, ...
				'LineStyle',s.CircleGeo(i).LineStyle);
			setappdata(h_circ,'LonLatRad',s.CircleGeo(i).lon_lat_rad);
			s.CircleGeo(i).ud.hcirc = h_circ;				s.CircleGeo(i).ud.parent = handles.axes1;
			s.CircleGeo(i).ud.h_fig = handles.figure1;		s.CircleGeo(i).ud.h_axes = handles.axes1;
			set(h_circ,'UserData',s.CircleGeo(i).ud,'buttondownfcn','uicirclegeo(''circlemousedown'')')
			draw_funs(h_circ,'SessionRestoreCircle')		% Set circle's uicontextmenu
		end
	end

	if (s.haveCircleCart)				% case of Cartesian circles
		for (i = 1:length(s.CircleCart))
			h_circ = line('Xdata',s.CircleCart(i).x,'Ydata',s.CircleCart(i).y,'Parent',handles.axes1, ...
				'LineWidth',s.CircleCart(i).LineWidth, 'color', s.CircleCart(i).color, 'Tag',s.CircleCart(i).tag, ...
				'LineStyle',s.CircleCart(i).LineStyle);
			setappdata(h_circ,'LonLatRad',s.CircleCart(i).lon_lat_rad);
			x = linspace(-pi,pi,360);
			setappdata(h_circ,'X',cos(x));		setappdata(h_circ,'Y',sin(x))	% Save unit circle coords
			s.CircleCart(i).ud.hcirc = h_circ;
			s.CircleCart(i).ud.parent = handles.axes1;
			draw_funs(h_circ,'SessionRestoreCircleCart')		% Set circle's uicontextmenu
		end
	end

	if (s.havePline)					% case of polylines
		for (i = 1:length(s.Pline))
			h_line = line('Xdata',s.Pline(i).x,'Ydata',s.Pline(i).y,'Parent',handles.axes1,'LineWidth',s.Pline(i).LineWidth,...
				'color',s.Pline(i).color,'Tag',s.Pline(i).tag, 'LineStyle',s.Pline(i).LineStyle);
			if (isfield(s.Pline(i),'Marker') && ~isempty(s.Pline(i).Marker))% New in 21-9-2011
				set(h_line, 'Marker', s.Pline(i).Marker, 'MarkerSize',s.Pline(i).Size, ...
					'MarkerFaceColor',s.Pline(i).FillColor, 'MarkerEdgeColor',s.Pline(i).EdgeColor)
			end		
			if (isfield(s.Pline(i),'appd') && ~isempty(s.Pline(i).appd))	% Need the isfield test for backward compat
				fdnames = fieldnames(s.Pline(i).appd);
				for (fd = 1:numel(fdnames))
					if (strcmp(fdnames{fd}, 'UserData'))		% This is a UserData stored in appdata
						set(h_line, 'UserData', s.Pline(i).appd.(fdnames{fd}))
					else
						setappdata(h_line, fdnames{fd}, s.Pline(i).appd.(fdnames{fd}))
					end
				end
			end
			if (isfield(s.Pline(i),'tag') && strcmp(s.Pline(i).tag, 'NEST'))
				draw_funs([], 'set_recTsu_uicontext', h_line)
			elseif (isfield(s.Pline(i),'LineInfo') && ~isempty(s.Pline(i).LineInfo))
				setappdata(h_line,'LineInfo',s.Pline(i).LineInfo)
				set(h_line,'UserData',1)
				draw_funs(h_line,'isochron',{s.Pline(i).LineInfo})
			elseif (isfield(s.Pline(i),'cont_label') && ~isempty(s.Pline(i).cont_label))
				setappdata(h_line,'cont_label',s.Pline(i).cont_label)		% Used in write_script
				set(h_line,'UserData',s.Pline(i).cont_label)				% Not sure if realy worth save this in UD
				draw_funs(h_line,'ContourLines','')		% The empty ('') is h_lable in DrawContours.
			else
				draw_funs(h_line,'line_uicontext')		% Set lines's uicontextmenu
			end
		end
	end

	if (s.havePlineAsPoints)			% case of polylines as points (markers) only
		for (i = 1:length(s.PlineAsPoints))
			h_line_pt = line('Xdata',s.PlineAsPoints(i).x, 'Ydata',s.PlineAsPoints(i).y,'Parent',handles.axes1, ...
				'LineStyle','none', 'Marker',s.PlineAsPoints(i).Marker, 'MarkerSize',s.PlineAsPoints(i).Size, ...
				'MarkerFaceColor',s.PlineAsPoints(i).FillColor, 'MarkerEdgeColor',s.PlineAsPoints(i).EdgeColor, ...
				'Tag',s.PlineAsPoints(i).tag);
			draw_funs(h_line_pt,'DrawSymbol')		% Set marker's uicontextmenu (tag is very important)
		end
	end

	if (s.haveSymbol)					% case of Symbols (line Markers)
		for (i = 1:length(s.Symbol))
			h_symb = line('Xdata',s.Symbol(i).x, 'Ydata',s.Symbol(i).y, 'Parent',handles.axes1,'Marker', ...
				s.Symbol(i).Marker,'MarkerSize', s.Symbol(i).Size,'MarkerFaceColor', s.Symbol(i).FillColor, ...
				'MarkerEdgeColor', s.Symbol(i).EdgeColor, 'Tag',s.Symbol(i).tag);
			draw_funs(h_symb,'DrawSymbol')			% Set symbol's uicontextmenu
		end
	end

	if (s.haveText)					% case of text strings
		try		s.Texto;			% Compatibility issue (Use a try because of compiler bugs)
		catch
			% Do it this way because compiled version canot tel 'Text' from 'text'
			t = load([PathName FileName],'Text');	s.Texto = t.Text;
		end
		for (i = 1:length(s.Texto))
			if (isempty(s.Texto(i).str)),		continue,	end
			h_text = text(s.Texto(i).pos(1), s.Texto(i).pos(2), s.Texto(i).pos(3), s.Texto(i).str,...
				'Parent',handles.axes1, 'Rotation',s.Texto(i).angle,...
				'FontAngle',s.Texto(i).FontAngle, 'Tag',s.Texto(i).Tag, 'FontWeight',s.Texto(i).FontWeight,...
				'color',s.Texto(i).color, 'FontName',s.Texto(i).FontName, 'FontSize',s.Texto(i).FontSize);
			if (isfield(s.Texto(i),'VerticalAlignment')),	set(h_text,'VerticalAlignment',s.Texto(i).VerticalAlignment),		end
			if (isfield(s.Texto(i),'HorizontalAlignment')),	set(h_text,'HorizontalAlignment',s.Texto(i).HorizontalAlignment),	end
			draw_funs(h_text,'DrawText')		% Set texts's uicontextmenu
		end
	end

	try		s.havePatches;		% 'Try' because compiler BUGs
	catch,	s.havePatches = false;
	end
	if (s.havePatches)			% case of patchs - NOTE, the Tags are currently lost 
		for (i = 1:numel(s.Patches))
			is_telha = false;
			is_arrow = strcmp(s.Patches(i).tag,'Arrow');
			if (strcmp(s.Patches(i).tag,'tapete_R') || strcmp(s.Patches(i).tag,'tapete'))
				s.Patches(i).x = reshape(s.Patches(i).x,4,numel(s.Patches(i).x)/4);
				s.Patches(i).y = reshape(s.Patches(i).y,4,numel(s.Patches(i).y)/4);
				is_telha = true;
			end
			h_patch = patch('XData',s.Patches(i).x, 'YData',s.Patches(i).y, 'Parent',handles.axes1, ...
				'LineWidth',s.Patches(i).LineWidth, 'EdgeColor',s.Patches(i).EdgeColor, 'FaceColor',s.Patches(i).FaceColor,...
				'LineStyle',s.Patches(i).LineStyle, 'Tag', s.Patches(i).tag);
			set(h_patch,'UserData',s.Patches(i).ud)
			if (isfield(s.Patches(i),'appd') && ~isempty(s.Patches(i).appd))	% Need the isfield test for backward compat
				fdnames = fieldnames(s.Patches(i).appd);
				for (fd = 1:numel(fdnames))
					if (~strcmp(fdnames{fd}, 'polygon_data'))
						setappdata(h_patch, fdnames{fd}, s.Patches(i).appd.(fdnames{fd}))
					end
				end
			end
			if (is_telha)
				draw_funs(h_patch,'telhas_patch')		% Set telhas's uicontextmenu
			elseif (is_arrow)
				draw_funs([],'set_vector_uicontext',h_patch)	% Set Arrow's uicontextmenu
			else
				draw_funs(h_patch,'line_uicontext')		% Set patch's uicontextmenu
			end
		end
		try				% 'MecaMag5' is new (20-9-2011), so we must use a try
			if (~isempty(s.MecaMag5)),	setappdata(handles.figure1, 'MecaMag5', s.MecaMag5),	end
		end
	end

	try
		if (s.haveCoasts),	datasets_funs('CoastLines', handles, s.coastUD);  end
		if (s.havePolitic)
			if (iscell(s.politicUD)),		s.politicUD = s.politicUD{1};	end
			datasets_funs('Political',handles, s.politicUD(2), s.politicUD(1));
		end
		if (s.haveRivers)
			if (iscell(s.riversUD)),		s.riversUD = s.riversUD{1};		end
			datasets_funs('Rivers', handles, s.riversUD(2), s.riversUD(1));
		end
	end
	handles.last_dir = fileparts(figName);		% We want it to point to where te .mat file came from.
	guidata(handles.figure1, handles);
	handles.fileName = figName;					% TRICK. To be used only in the next call to recentFiles()
	handles = recentFiles(handles);				% Insert session into "Recent Files" & NOT NOT NOT save handles
	set(handles.figure1, 'Name',figName)
	if (tala == 0 && ~isempty(s.grd_name))		% Only now to not mess with the "current figure"
		warndlg(['The file ' s.grd_name ' doesn''t exists on the directory it was when the session was saved. Put it back there.'],'Warning')
	end

% ------------------------------------------------------------------------------------------------
function FileSaveSession_CB(handles)
	if (handles.image_type == 0),		return,		end
	str1 = {'*.mat;*.MAT', 'Data files (*.mat,*.MAT)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select session file name','put','.mat');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	set(handles.figure1,'pointer','watch')
	grd_name = handles.fileName;	% Use this variable name for compatibility reason
	if (handles.image_type == 20),	grd_name = [];	end		% The idea is if we have only vector data it may have been edited
	img_pal = get(handles.figure1,'Colormap');		illumComm = [];		illumType = handles.Illumin_type;
	map_limits = getappdata(handles.axes1,'ThisImageLims');
	if (handles.validGrid && handles.Illumin_type >= 1 && handles.Illumin_type <= 4)
		illumComm = getappdata(handles.figure1,'illumComm');
	end
	ALLlineHand = findobj(get(handles.axes1,'Child'),'Type','line');
	m = 1;
	haveMBtrack = 0;	havePline = 0;		haveText = 0;	haveSymbol = 0;		haveCircleGeo = 0;
	haveCircleCart = 0; havePlineAsPoints = 0;  havePatches = 0;haveCoasts = 0; havePolitic = 0;	haveRivers = 0;
	MBtrack = [];	MBbar = [];		Pline = [];		Symbol = [];	CircleGeo = [];	MecaMag5 = [];
	CircleCart = [];	PlineAsPoints = [];			coastUD = [];	politicUD = [];	riversUD = [];

	h = findobj(ALLlineHand,'Tag','Symbol');		% case of a Symbol (in fact a line Marker)
	if (~isempty(h))
		nO = numel(h);
		Symbol = struct('x',cell(1,nO), 'y',cell(1,nO), 'Marker',cell(1,nO), 'MarkerSize',cell(1,nO), ...
			'MarkerFaceColor',cell(1,nO), 'MarkerEdgeColor',cell(1,nO), 'tag',cell(1,nO));
		for (i = 1:nO)
			xx = get(h(i),'XData');		yy = get(h(i),'YData');
			Symbol(i).x = xx(:);		Symbol(i).y = yy(:);
			Symbol(i).Marker = get(h(i),'Marker');
			Symbol(i).Size = get(h(i),'MarkerSize');
			Symbol(i).FillColor = get(h(i),'MarkerFaceColor');
			Symbol(i).EdgeColor = get(h(i),'MarkerEdgeColor');
			Symbol(i).tag = get(h(i),'Tag');
		end
		haveSymbol = 1;
		ALLlineHand = setxor(ALLlineHand, h);       % Those are processed, so remove them from handles list
	end
	
	h = findobj(ALLlineHand,'Tag','Pointpolyline');		% Polyline with only markers are particular line cases
	if (~isempty(h))
		nO = numel(h);
		PlineAsPoints = struct('x',cell(1,nO), 'y',cell(1,nO), 'Marker',cell(1,nO), 'MarkerSize',cell(1,nO), ...
			'MarkerFaceColor',cell(1,nO), 'MarkerEdgeColor',cell(1,nO), 'tag',cell(1,nO));
		for (i = 1:nO)
			xx = get(h(i),'XData');			yy = get(h(i),'YData');
			PlineAsPoints(i).x = xx(:);		PlineAsPoints(i).y = yy(:);
			PlineAsPoints(i).Marker = get(h(i),'Marker');
			PlineAsPoints(i).Size = get(h(i),'MarkerSize');
			PlineAsPoints(i).FillColor = get(h(i),'MarkerFaceColor');
			PlineAsPoints(i).EdgeColor = get(h(i),'MarkerEdgeColor');
			PlineAsPoints(i).tag = 'Pointpolyline';
		end
		ALLlineHand = setxor(ALLlineHand, h);		havePlineAsPoints = 1;
	end

	h = findobj(ALLlineHand,'Tag','circleGeo');			% circles are particular line cases
	h = [h; findobj(ALLlineHand,'Tag','CircleEuler')];
	if (~isempty(h))
		nO = numel(h);
		CircleGeo = struct('x',cell(1,nO), 'y',cell(1,nO), 'LineWidth',cell(1,nO), 'LineStyle',cell(1,nO), ...
			'color',cell(1,nO), 'lon_lat_rad',cell(1,nO), 'ud',cell(1,nO), 'tag',cell(1,nO));
		for (i = 1:nO)
			xx = get(h(i),'XData');		yy = get(h(i),'YData');
			CircleGeo(i).x = xx(:);		CircleGeo(i).y = yy(:);
			CircleGeo(i).LineWidth = get(h(i),'LineWidth');
			CircleGeo(i).LineStyle = get(h(i),'LineStyle');
			CircleGeo(i).color = get(h(i),'color');
			CircleGeo(i).tag = get(h(i),'Tag');
			CircleGeo(i).lon_lat_rad = getappdata(h(i),'LonLatRad');
			CircleGeo(i).ud = get(h(i),'UserData');	% UserData contains alot of need info
		end
		ALLlineHand = setxor(ALLlineHand, h);		haveCircleGeo = 1;
	end
	
	h = findobj(ALLlineHand,'Tag','circleCart');
	if (~isempty(h))
		nO = numel(h);
		CircleCart = struct('x',cell(1,nO), 'y',cell(1,nO), 'LineWidth',cell(1,nO), 'LineStyle',cell(1,nO), ...
			'color',cell(1,nO), 'lon_lat_rad',cell(1,nO), 'ud',cell(1,nO), 'tag',cell(1,nO));
		for (i = 1:nO)
			xx = get(h(i),'XData');			yy = get(h(i),'YData');
			CircleCart(i).x = xx(:);		CircleCart(i).y = yy(:);
			CircleCart(i).LineWidth = get(h(i),'LineWidth');
			CircleCart(i).LineStyle = get(h(i),'LineStyle');
			CircleCart(i).color = get(h(i),'color');
			CircleCart(i).tag = 'circleCart';
			CircleCart(i).lon_lat_rad = getappdata(h(i),'LonLatRad');
			CircleCart(i).ud = get(h(i),'UserData');
		end
		ALLlineHand = setxor(ALLlineHand, h);		haveCircleCart = 1;
	end

	h = findobj(ALLlineHand,'Tag','MBtrack');		% case of a MBtrack line (DOESN'T WORK BECAUSE 'MBtrack%d')
	if (~isempty(h))
		nO = numel(h);
		MBtrack = struct('x',cell(1,nO), 'y',cell(1,nO), 'LineWidth',cell(1,nO), 'LineStyle',cell(1,nO), ...
			'color',cell(1,nO), 'swathRatio',cell(1,nO), 'tag',cell(1,nO));
		for (i = 1:nO)
			xx = get(h(i),'XData');		yy = get(h(i),'YData');
			MBtrack(i).x = xx(:);		MBtrack(i).y = yy(:);
			MBtrack(i).LineWidth = get(h(i),'LineWidth');
			MBtrack(i).LineStyle = get(h(i),'LineStyle');
			MBtrack(i).color = get(h(i),'color');
			MBtrack(i).swathRatio = getappdata(h(i),'swathRatio');
			MBtrack(i).tag = 'MBtrack';
		end
		ALLlineHand = setxor(ALLlineHand, h);		haveMBtrack = 1;
		draw_funs(h,'MBtrackUictx')					% Set track's uicontextmenu

		h = findobj(ALLlineHand,'Tag','swath_w');	% case of a MBtrack's bar line
		nO = numel(h);
		MBbar = struct('x',cell(1,nO), 'y',cell(1,nO), 'LineWidth',cell(1,nO), 'LineStyle',cell(1,nO), ...
			'color',cell(1,nO), 'n_vert',cell(1,nO), 'tag',cell(1,nO));
		for (i = 1:nO)
			xx = get(h(i),'XData');		yy = get(h(i),'YData');
			MBbar(i).x = xx(:);			MBbar(i).y = yy(:);
			MBbar(i).LineWidth = get(h(i),'LineWidth');
			MBbar(i).LineStyle = get(h(i),'LineStyle');
			MBbar(i).color = get(h(i),'color');
			MBbar(i).n_vert = get(h(i),'UserData');
			MBbar(i).tag = 'swath_w';
			set(h(i),'Userdata',i)
		end
		ALLlineHand = setxor(ALLlineHand, h);
		draw_funs(h,'MBbarUictx')				% Set track bar's uicontextmenu
	end

	for (i = 1:numel(ALLlineHand))
		tag = get(ALLlineHand(i),'Tag');
		if (strcmp(tag,'CoastLineNetCDF'))
			haveCoasts = 1;		coastUD = get(ALLlineHand(i),'UserData');	
		elseif (strcmp(tag,'PoliticalBoundaries'))
			havePolitic = 1;	politicUD = get(ALLlineHand(i),'UserData');	
		elseif (strcmp(tag,'Rivers'))
			haveRivers = 1;		riversUD = get(ALLlineHand(i),'UserData');	
		else		% for the time beeing, it applyies to simple polylines
			xx = get(ALLlineHand(i),'XData');		yy = get(ALLlineHand(i),'YData');
			Pline(m).x = xx(:);			Pline(m).y = yy(:);
			Pline(m).LineWidth = get(ALLlineHand(i),'LineWidth');
			Pline(m).LineStyle = get(ALLlineHand(i),'LineStyle');
			Pline(m).color = get(ALLlineHand(i),'color');
			Pline(m).tag = tag;
			Marker = get(ALLlineHand(i),'Marker');
			if (Marker(1) ~= 'n')		% Not 'none' so something.
				Pline(m).Marker = Marker;
				Pline(m).Size = get(ALLlineHand(i),'MarkerSize');
				Pline(m).FillColor = get(ALLlineHand(i),'MarkerFaceColor');
				Pline(m).EdgeColor = get(ALLlineHand(i),'MarkerEdgeColor');
			end
			if (strcmp(tag, 'contour'))	% For contours we must pass along the label info
				Pline(m).cont_label = getappdata(ALLlineHand(i),'cont_label');
			elseif (isappdata(ALLlineHand(i),'LineInfo'))
				Pline(m).LineInfo = getappdata(ALLlineHand(i),'LineInfo');
			end
			% NEW at 19-Oct-2014, save all appdatas which will cause that the above two will potentially be repeated
			app = getappdata(ALLlineHand(i));
			ud = get(ALLlineHand(i), 'UserData');
			if (~isempty(ud)),	app.UserData = ud;		end
			if (~isempty(fieldnames(app)))
				if (isfield(app, 'polygon_data'))	% Remove this because it will be restored by ui_edit_polygon()
					app = rmfield(app, 'polygon_data');
				end
				if (~isempty(fieldnames(app)))
					Pline(m).appd = app;
				end
			end
			
			m = m + 1;		havePline = 1;
		end
	end

	% Patches may have associated particular meanings (eg Focal Mecas), but
	% they will loose them here. Maybe in the future I'll do something better.
	ALLpatchHand = findobj(get(handles.axes1,'Child'),'Type','patch');
	ALLpatchHand = ALLpatchHand(end:-1:1);		% Don't know if this always good but respects stack order
	nO = numel(ALLpatchHand);
	Patches = struct('x',cell(1,nO), 'y',cell(1,nO), 'LineWidth',cell(1,nO), 'LineStyle',cell(1,nO), ...
		'EdgeColor',cell(1,nO), 'FaceColor',cell(1,nO), 'ud',cell(1,nO), 'tag',cell(1,nO));
	for (i = 1:nO)
		xx = get(ALLpatchHand(i),'XData');		yy = get(ALLpatchHand(i),'YData');
		Patches(i).x = xx(:);					Patches(i).y = yy(:);
		Patches(i).LineWidth = get(ALLpatchHand(i),'LineWidth');
		Patches(i).LineStyle = get(ALLpatchHand(i),'LineStyle');
		Patches(i).EdgeColor = get(ALLpatchHand(i),'EdgeColor');
		Patches(i).FaceColor = get(ALLpatchHand(i),'FaceColor');
		Patches(i).ud   = get(ALLpatchHand(i),'UserData');
		Patches(i).tag  = get(ALLpatchHand(i),'Tag');
		app = getappdata(ALLpatchHand(i));
		if (~isempty(fieldnames(app)))
			if (isfield(app, 'polygon_data'))	% Remove this because it will be restored by ui_edit_polygon()
				app = rmfield(app, 'polygon_data');
			end
			if (~isempty(fieldnames(app)))
				Patches(i).appd = app;
			end
		end
		havePatches = 1;
	end

	ALLtextHand = findobj(get(handles.axes1,'Child'),'Type','text');
	nO = numel(ALLtextHand);
	Texto = struct('str',cell(1,nO), 'pos',cell(1,nO), 'angle',cell(1,nO), 'color',cell(1,nO), ...
		'FontSize',cell(1,nO), 'HorizontalAlignment',cell(1,nO), 'VerticalAlignment',cell(1,nO));
	for (i = 1:nO)
		Texto(i).str = get(ALLtextHand(i),'String');
		if (isempty(Texto(i).str)),  continue,	end
		Texto(i).pos = get(ALLtextHand(i),'Position');		Texto(i).FontAngle = get(ALLtextHand(i),'FontAngle');
		Texto(i).angle = get(ALLtextHand(i),'Rotation');	Texto(i).Tag = get(ALLtextHand(i),'Tag');
		Texto(i).color = get(ALLtextHand(i),'color');		Texto(i).FontName = get(ALLtextHand(i),'FontName');
		Texto(i).FontSize = get(ALLtextHand(i),'FontSize');	Texto(i).FontWeight = get(ALLtextHand(i),'FontWeight');
		Texto(i).HorizontalAlignment = get(ALLtextHand(i),'HorizontalAlignment');
		Texto(i).VerticalAlignment = get(ALLtextHand(i),'VerticalAlignment');
		haveText = 1;
	end

	if (havePatches)
		MecaMag5 = getappdata(handles.figure1,'MecaMag5');
	end

	IamTINTOL = ~isempty(getappdata(handles.figure1, 'L_UsedByGUIData'));
	is_defRegion = handles.is_defRegion;

	save(fname,'grd_name','img_pal', 'havePline','Pline', 'haveMBtrack', 'MBtrack','MBbar', ...
		'haveText','Texto', 'haveSymbol','Symbol', 'haveCircleGeo','CircleGeo', 'haveCircleCart', ...
		'havePlineAsPoints','PlineAsPoints','CircleCart', 'map_limits', 'havePatches', 'Patches', ...
		'haveCoasts', 'coastUD','havePolitic', 'politicUD','haveRivers', 'riversUD', 'illumComm', ...
		'illumType', 'MecaMag5', 'IamTINTOL', 'is_defRegion', '-v6')
	set(handles.figure1,'pointer','arrow')

	% Trick to shut up stupid MLint warnings
	if (0 && grd_name && img_pal && map_limits && illumComm && haveMBtrack && havePline && haveText && haveSymbol), end
	if (0 && haveCircleCart && havePlineAsPoints && haveCoasts && illumComm && haveMBtrack && illumType), end
	if (0 && havePolitic && haveRivers && MBtrack && MBbar && Symbol && CircleGeo && MecaMag5 && CircleCart), end
	if (0 && PlineAsPoints && coastUD && politicUD && riversUD && haveCircleGeo && IamTINTOL && is_defRegion), end

% --------------------------------------------------------------------
function ImageMapLimits_CB(handles, opt)
% Change the Image OR the Map extents limits by asking it's corner coordinates
% region contains [x_min x_max y_min y_max is_geog] in PIXEL REG MODE
	if (strcmpi(opt,'map'))				% Change the display extentents and go away
		region = bg_region('with_limits', [get(handles.axes1, 'XLim') get(handles.axes1, 'YLim')] );
		if isempty(region),		return,		end
		set(handles.axes1, 'XLim',region(1:2), 'YLim',region(3:4))
		setappdata(handles.axes1,'ThisImageLims',region)
		return
	end

	img = get(handles.hImg,'CData');
	if (strcmp(opt, 'img'))
		region = bg_region('empty');
		if isempty(region),		return,		end
		X = region(1:2);		Y = region(3:4);
		x_inc = diff(X) / size(img,2);			y_inc = diff(Y) / size(img,1);
	else									% Fit to -R-0.5/0.5/-0.5/0.5
		X = [-0.5 0.5];			Y = [-0.5 0.5];
		x_inc = diff(X) / size(img,2);			y_inc = diff(Y) / size(img,1);
		aspect = size(img,1) / size(img,2);
		if (aspect > 1),		X = X / aspect;		% Taller image, contract X
		elseif (aspect < 1),	Y = Y * aspect;		% Wider image, contract Y
		end
		aux_funs('addUI', handles)
	end
	%x_inc = diff(X) / size(img,2);			y_inc = diff(Y) / size(img,1);
	dx2 = x_inc / 2;						dy2 = y_inc / 2;
	X = X + [dx2 -dx2];						Y = Y + [dy2 -dy2];		% X,Y in grid-reg so that the pix-reg info = region
	handles.head(1:4) = [X Y];				handles.head(8:9) = [x_inc y_inc];  
	handles.geog = aux_funs('guessGeog',handles.head(1:4));			% Trust more in this guess than on bg_region()
	handles.fileName = [];					% Not loadable in session
	if (~handles.validGrid),				handles.image_type = 3;
	else
		X = linspace(X(1),X(2),size(img,2));		Y = linspace(Y(1),Y(2),size(img,1));
		setappdata(handles.figure1,'dem_x',X);  	setappdata(handles.figure1,'dem_y',Y);
	end
	
	% Flipud the image if necessary
	if (strcmp(get(handles.axes1,'YDir'),'reverse')),	img = flipdim(img,1);	end
	show_image(handles,'New Limits',X,Y,img,handles.validGrid,'xy',0);
	if (strcmp(opt, 'fit') && aspect ~= 1)	% Since the aspect ratio ~= 1 we must change the display limits
		set(handles.axes1, 'XLim',[-0.5 0.5], 'YLim',[-0.5 0.5])
		setappdata(handles.axes1,'ThisImageLims',[-0.5 0.5 -0.5 0.5])		
	end

% --------------------------------------------------------------------
function GeophysicsSwanPlotStations_CB(handles)
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	zoom_state(handles,'maybe_off');
	pt = click_e_point(1,'crosshair');
	if (isempty(pt)),	return,		end
	handles.maregraphs_count = handles.maregraphs_count + 1;	% Count number of maregraphs
	h = line('XData',pt(1,1),'YData',pt(1,2),'Parent',handles.axes1,'Marker','o','MarkerFaceColor','y',...
			'MarkerEdgeColor','k','MarkerSize',10,'Tag','Maregraph');
	draw_funs(h,'DrawSymbol')			% Set symbol's uicontextmenu
	zoom_state(handles,'maybe_on');		guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function GRDdisplay(handles, X, Y, Z, head, tit, name, srsWKT)
% Show matrix Z in a new window.
	if (isa(Z,'double')),		Z = single(Z);	end
	zz = grdutils(Z,'-L');		% Do not trust that handles.have_nans is always correct
	head(5:6) = double(zz(1:2));
	tmp.head = head;			tmp.X = X;		tmp.Y = Y;		tmp.geog = handles.geog;	tmp.name = name;
	if (nargin == 8 && ~isempty(srsWKT)),		tmp.srsWKT = srsWKT;	end
	set(handles.figure1,'pointer','arrow')
	h = mirone(Z,tmp);			newHand = guidata(h);
	thematic_pal = getappdata(handles.figure1, 'thematic_pal');		% Because of the bloody cpt2cmap crash
	if (~isempty(thematic_pal)),	setappdata(newHand.figure1,'thematic_pal',thematic_pal),	end

% --------------------------------------------------------------------
function FileSaveImgGrdGdal_CB(handles, opt1, opt2)
% OPT1 = DRIVER == GTiff, HFA (erdas), ENVI, ECW, JP2ECW
% OPT2 == grid -> saves the underlaying grid
% ELSE -> do a screen capture
	if (handles.no_file),		return,		end
	if (strcmp(opt2,'grid') && ~handles.validGrid)
		errordlg('You don''t have a Grid loaded, so OBVIOUSLY you cannot save it.','Error');  return
	end
	flip = 0;			% FLIP = 1, when images need to be UD fliped before saving
	switch opt1
		case 'GeoTiff',		str1 = {'*.tif;*.TIF;*.tiff;*.TIFF', 'GeoTiff (*.tif;*.TIF;*.tiff;*.TIFF)'};	driver = 'GTiff';
		case 'Erdas',		str1 = {'*.img;*.IMG', 'Erdas (*.img;*.IMG)'};			driver = 'HFA';
		case 'Envi',		str1 = {'*.img;*.IMG', 'Envi (*.img;*.IMG)'};			driver = 'ENVI';
		case 'ESRI',		str1 = {'*.bil;*.BIL', 'Esri (*.bil;*.BIL)'};			driver = 'EHdr';
		case 'JP2K',		str1 = {'*.jp2;*.JP2', 'Jpeg2000 (*.jp2;*.JP2)'};		driver = 'JP2ECW';
		case 'PDF',			str1 = {'*.pdf;*.PDF', 'GeoPDF (*.pdf;*.PDF)'};			driver = 'GeoPDF';
	end
	[FileName,PathName] = put_or_get_file(handles,str1,['Select ' opt1 ' file name'],'put', str1{1}(2:5));
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	head = handles.head;
	% 'ThisImageLims' contains the limits as seen from the pixel-registration stand-point,
	%  and that is the convention that we will use to save rasters using GDAL.
	imgLims = getappdata(handles.axes1,'ThisImageLims');
	if (strcmp(opt2,'grid'))
		[X,Y,Z,head] = load_grd(handles);
		if isempty(Z),		return,		end		% An error message was already issued
		if (handles.was_int16)
			if (handles.have_nans)		% Restore the original Nodata value, or use -32768 if we don't know it
				Z(isnan(Z(:))) = handles.Nodata_int16;
			end
			Z = int16(Z);
		end
	else										% 'image'
		Z = snapshot(handles.figure1,'noname');		pause(0.01)
		if (isempty(Z)),		return,		end
		head(8) = (imgLims(2)-imgLims(1)) / size(Z,2);
		head(9) = (imgLims(4)-imgLims(3)) / size(Z,1);
	end

	hdr.name = fname;		hdr.driver = driver;	hdr.Geog = handles.geog;
	projWKT = getappdata(handles.figure1,'ProjWKT');
	Proj4 = getappdata(handles.figure1,'Proj4');
	if (~isempty(projWKT)),		hdr.projWKT = projWKT;
	elseif (~isempty(Proj4)),	hdr.projWKT = ogrproj(Proj4);
	end

	try
		hdr.Xinc = head(8);		hdr.Yinc = head(9);
		hdr.ULx = imgLims(1);	hdr.ULy = imgLims(4);
		hdr.Reg = 1;			hdr.Flip = flip;
	catch
		errordlg('Shit, image header was not saved as it should.','Error'),		return
	end
	if ( (strcmp(opt2,'img') || strcmp(opt2,'screen')) && ndims(Z) == 2 )
		hdr.Cmap = get(handles.figure1,'ColorMap');
	end
	if ( strcmp(driver,'GTiff') )
		gcps = getappdata(handles.figure1,'GCPregImage');
		if (~isempty(gcps))
			h = findobj(handles.axes1,'Tag','GCPpolyline');		% If we have them visible, they might have been edited
			if (~isempty(h))
				x = get(h,'XData');		y = get(h,'YData');		gcps(:,1:2) = [x(:) y(:)];
			end
			hdr.gcp = gcps;
		end
	end

	if (ndims(Z) == 2 && strcmp(driver,'HFA') && strcmp(get(handles.hImg, 'CDataMapping'), 'scaled') )
		clim = get(handles.axes1, 'CLim');
		tmp = alloc_mex([size(Z) diff(clim)+1], 'uint8');	n = 1;
		for (k = clim(1):clim(end))
			tmp0 = alloc_mex(size(Z), 'uint8');
			tmp0(Z == k) = 1;
			tmp(:,:, n) = tmp0;				n = n + 1;
		end
		Z = tmp;		hdr.meta = 'NBITS=1';
		if (strcmp(get(handles.axes1,'YDir'),'reverse')),	hdr.Flip = 1;	end
	end
	if (strcmp(driver(1:3), 'JP2')),		ECWpatch(handles,'jp2');	end		% Make sure that the ECW driver is ready
	gdalwrite(Z,hdr)

% --------------------------------------------------------------------
function GridToolsHistogram_CB(handles, opt)
% OPT2 if present is a structure with two fields: opt2.Z (the matrix); opt2.head (its 9 elements header)
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	if (nargin == 1)						% Use entire grid
		[X,Y,Z,head] = load_grd(handles);
		if isempty(Z),  	return,		end
		z_min = head(5);	z_max = head(6);
	else									% Use a subset grid extracted from a rectangular area
		Z = opt.Z;			z_min = opt.head(5);	z_max = opt.head(6);
	end
	binwidth = (z_max - z_min) / 20;	% Default to 20 bins
	resp = inputdlg({'Enter Bin Width (default is 20 bins)'},'Histogram',[1 38],{sprintf('%g',binwidth)});	pause(0.01);
	resp = abs( str2double(resp{1}) );
	if (isnan(resp)),	set(handles.figure1,'pointer','arrow'),		return,		end
	n = round( (z_max - z_min) / resp );
	[n,xout] = histo_m('hist',Z(:),n,[z_min z_max]);
	h = mirone;							% Create a new Mirone figure
	mirone('FileNewBgFrame_CB', guidata(h), [xout(1) xout(end) 0 max(n) 0], [600 600],'Grid Histogram');
	histo_m('bar',xout,n,'hist');
	hand2 = guidata(h);
	set(hand2.axes1, 'XLim',[xout(1) xout(end)], 'YLim', [0 max(n)])	% Have to because histo_m had screwed them
	set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function GridToolsGridMask_CB(handles)
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	if (~handles.have_nans)
		msgbox('This option only works on grids that have NaNs.','Warning'),	return
	end
	[X,Y,Z,head] = load_grd(handles);
	if isempty(Z),		return,		end		% An error message was already issued
	Z(~isnan(Z)) = 1;
	GRDdisplay(handles,X,Y,Z,head,'Mask grid','Mask grid')

% --------------------------------------------------------------------
function GridToolsSectrum_CB(handles, opt1, opt2)
% OPT1 == 'Amplitude'	-> compute amplitude spectrum
% OPT1 == 'Power'		-> compute power spectrum
% OPT1 == 'Autocorr'	-> compute autocorrelation
% OPT1 == 'Allopts'		-> call the fft_stuff window
% OPT2	if present is a structure with two fields: opt2.Z (the matrix); opt2.head (its 9 elements header)
%  -"-	OR a line/patch handle used to do band filtering
	quick = 0;		head = [];
	if (handles.have_nans)
		warndlg('This grid has NaNs. That is not allowed in FFTs','Warning'),	return
	end
	if (nargin == 2)						% Use entire grid
		if (handles.validGrid)
			[X,Y,Z,head] = load_grd(handles);
			if isempty(Z),		return,		end
		else								% IMAGE FFT
			Z = get(handles.hImg, 'CData');
			if (ndims(Z) == 3)
				Z = cvlib_mex('color', Z, 'rgb2gray');
				warndlg('Currently FFT of RGB images is done by converting the RGB to gray.','Warning');
			end
			if (~isempty(Z)),		head = handles.head;	Z = double(Z);		end			% Ghrrrrrr
		end
	else									% Use a (maybe) subset grid extracted from a rectangular area
		if (strcmp(opt1(2:end),'pass'))		% Very special case (no subset)
			h = getappdata(handles.figure1, 'ParentFig');
			if (~ishandle(h))
				errordlg('Too late. You killed the figure with the original data.'),	return
			end
			handles = guidata(h);			% Is the parent figure handles that FFT_STUFF() will need
			if (handles.validGrid)
				[X,Y,Z,head] = load_grd(handles);
				if isempty(Z),	return,		end
			else							% IMAGE FFT
				Z = get(handles.hImg, 'CData');
				if (ndims(Z) == 3),		Z = cvlib_mex('color', Z, 'rgb2gray');		end
				head = handles.head;	Z = double(Z);			% Ghrrrrrr
			end
			handles.XXXhLine = opt2;		% Dumb field to be used only by FFT_STUFF()
			guidata(handles.figure1, handles)
		else
			Z = opt2.Z;		head = opt2.head;	quick = 1;
		end
	end

	if (quick),				set(handles.figure1,'pointer','watch');		end
	if (~isempty(Z)),		fft_stuff(handles.figure1, Z, head, handles.geog, opt1)
	else					fft_stuff
	end
	if (quick),		set(handles.figure1,'pointer','arrow');		end

% --------------------------------------------------------------------
function Zsmoothed = GridToolsSmooth_CB(handles, opt)
% Smooth out a grid using B-splines
%
% Because the smoothing factor is choose interactively, OPT is to contain a multiplicative factor to
% be applied to the automatically estimated 'p' factor. This gives a mechanism to call this
% function from another one and send back the final result without stopping due to any questions.

	if (aux_funs('msg_dlg',14,handles)),		return,		end
	[X,Y,Z,head,m,n] = load_grd(handles,'double');
	if isempty(Z),		return,		end				% An error message was already issued

	[pp, p_guess] = spl_fun('csaps',{Y(1:5),X(1:5)},Z(1:5,1:5));		% Get a good estimate of p
	if (nargin == 1)
		prompt = {'Enter smoothing p paramer'};		dlg_title = 'Smoothing parameter input';
		defAns = {sprintf('%.12f',p_guess{1})};		resp  = inputdlg(prompt,dlg_title,[1 38],defAns);	pause(0.01)
		resp = abs( str2double(resp{1}) );
		if (isnan(resp)),	return,		end
		p_to_use = resp;
	else
		p_to_use = p_guess * opt;
	end

	set(handles.figure1,'pointer','watch')
	Lim = handles.grdMaxSize * 0.06;
	nl0 = round(Lim/(n*16*8) / 2);
	if (rem(nl0,2) ~= 0),	nl0 = nl0 + 1;  end % Don't know why, but I rather have it as a even number
	nl = round(nl0 * 0.8);
	if (rem(nl,2) ~= 0),	nl = nl + 1;	end
	skirt = ceil(nl0 - nl) - 2;

	if (handles.have_nans)
		[Z, indNaNs] = fillGridGaps(handles, Z);	% Fill the gaps with the minimum curvature interpolator
	end

	[ind_s,ind] = tile(m,nl,skirt);				% Get indexes for tiling.
	if (size(ind_s,1) > 1)
		last_row = 1;
		for k = 1:size(ind_s,1)
			tmp1 = (ind_s(k,1):ind_s(k,2));		% Indexes with overlapping zone
			tmp2 = ind(k,1):ind(k,2);			% Indexes of chunks without the overlaping zone
			pp = spl_fun('csaps',{Y(tmp1),X},Z(tmp1, 1:end), p_to_use);
			tmp = spl_fun('fnval',pp,{Y(tmp1),X});
			tmp = tmp(tmp2, 1:end);
			Z(last_row : (last_row + size(tmp,1) - 1), :) = tmp;	% Z is already a copy (is f double)
			last_row = last_row + size(tmp,1);
		end
		clear pp tmp;
	else
		pp = spl_fun('csaps',{Y,X},Z, p_to_use);
		Z = spl_fun('fnval',pp,{Y,X});		clear pp;
	end
	Z = single(Z);
	if (handles.have_nans),			Z(indNaNs) = NaN;		end		% Obey the law of NaN conservation

	if (nargout)
		Zsmoothed = Z;		% No worries, it's not a real copy
	else
		tit = ['Spline smoothed grid. p parameter used = ' sprintf('%.12f', p_to_use)];
		GRDdisplay(handles,X,Y,Z,head,tit,'Spline smoothed grid');
	end

% --------------------------------------------------------------------
function sdgGrid = GridToolsSDG_CB(handles, opt, opt2)
% Compute the Second Derivative in direction of Gradient, which was initially the
% recommended method to determine the FOS (foot of continental slope) for purposes
% of continental shelf extension of maritime states from article 76 of LOS (law of the sea)
%
% Because one normally need to smooth up the grid before computing the derivatives and the
% smoothing factor is choose interactively, OPT2 is to contain a multiplicative factor to
% be applied to the automatically estimated 'p' factor. This gives a mechanism to call this
% function from another one and send back the final result without stopping due to any questions.

	if (aux_funs('msg_dlg',14,handles)),		return,		end
	[X,Y,Z,head,m,n] = load_grd(handles,'double');
	if isempty(Z),		return,		end
	[pp, p_guess] = spl_fun('csaps',{Y(1:5),X(1:5)},Z(1:5,1:5));	% Get a good estimate of p
	if (nargin == 2)
		prompt = {'Enter smoothing p paramer'};		dlg_title = 'Smoothing parameter input';
		defAns = {sprintf('%.12f',p_guess{1})};		resp  = inputdlg(prompt,dlg_title,[1 38],defAns);	pause(0.05)
		if isempty(resp),	return,		end
		p_to_use = str2double(resp{1});
	else
		p_to_use = p_guess * opt2;
	end

	% Apparently the biggest memory monster (pp) obeys roughly to the following relation:
	% n_row x 4 x n_column x 4 * 8. So if I impose a limit "Lim" to this monster, I should
	% be able to find out a chunk height that fits to the above relation. However, 'Lim'
	% should be of the order of handles.grdMaxSize because the tiling is meant to be used only
	% for large grids (relative to the available ram). Also, I don't understand why contiguous
	% chunks don't patch perfectly . From what I could find, it depends heavily on the p parameter.
	% A skirt 20% of chunk's height seams to produce a reasonable result for high p's, but not always perfect.

	set(handles.figure1,'pointer','watch')
	Lim = handles.grdMaxSize * 0.6;
	nl0 = round(Lim / (n*16*8) / 3);			% Divide by 3 because of the auxiliary variables
	if (rem(nl0,2) ~= 0),	nl0 = nl0 + 1;	end	% Don't know why, but I rather have it as a even number
	nl = round(nl0 * 0.8);
	if (rem(nl,2) ~= 0),	nl = nl + 1;	end
	skirt = ceil(nl0 - nl) - 2;

	if (handles.have_nans)
		% Although the SPL toolbox says NaNs and Infs are ignored, the fact is that their presence
		% (well, their removal) makes a hell out of the result with lots of propagation effects. So,
		% we fill the gaps with the minimum curvature interpolator, compute the SDG and reset the
		% NaNs at the end. No tilling needs to be done here since Z is already a copy (double) and
		% and the filling gaps are actually done in-place.
		[Z, indNaNs] = fillGridGaps(handles, Z);
	end

	[ind_s,ind] = tile(m,nl,skirt);				% Get indexes for tiling.
	R = zeros(m,n);
	if (size(ind_s,1) > 1),						% There is still a very strange thing that I don't understand.
		last_row = 1;
		for k = 1:size(ind_s,1)
			tmp1 = (ind_s(k,1):ind_s(k,2));		% Indexes with overlapping zone
			tmp2 = ind(k,1):ind(k,2);			% Indexes of chunks without the overlaping zone
			pp = spl_fun('csaps',{Y(tmp1),X},Z(tmp1, 1:end), p_to_use);
			DfDX = spl_fun('fnder',pp, [1 0]);			% df / dx
			vx = spl_fun('fnval',DfDX,{Y(tmp1),X});		clear DfDX;
			DfDY = spl_fun('fnder',pp, [0 1]);			% df / dy
			vy = spl_fun('fnval',DfDY,{Y(tmp1),X});		clear DfDY;
			D2fDX2 = spl_fun('fnder',pp, [2 0]);		% d^2f / dx^2
			Hxx = spl_fun('fnval',D2fDX2,{Y(tmp1),X});	clear D2fDX2;
			D2fDY2 = spl_fun('fnder',pp, [0 2]);		% d^2f / dy^2
			Hyy = spl_fun('fnval',D2fDY2,{Y(tmp1),X});	clear D2fDY2;
			D2fDXDY = spl_fun('fnder',pp, [1 1]);		clear pp;	% d^2f / (dx dy)
			Hxy = spl_fun('fnval',D2fDXDY,{Y(tmp1),X});	clear D2fDXDY;
			tmp = zeros(ind(k,2),n);
			for (j = 1:n)
				for (i = tmp2)
					v = [vx(i,j) vy(i,j)] / norm([vx(i,j) vy(i,j)]);		% eq(2)
					tmp(i,j) = (v * [Hxx(i,j) Hxy(i,j); Hxy(i,j) Hyy(i,j)] * v') / (v*v');
				end
			end
			tmp = tmp(tmp2, 1:end);
			R(last_row : (last_row + size(tmp,1) - 1), :) = tmp;
			last_row = last_row + size(tmp,1);
		end
	else
		pp = spl_fun('csaps',{Y,X},Z, p_to_use);
		DfDX = spl_fun('fnder',pp, [1 0]);			% df / dx
		vx = spl_fun('fnval',DfDX,{Y,X});			clear DfDX;
		DfDY = spl_fun('fnder',pp, [0 1]);			% df / dy
		vy = spl_fun('fnval',DfDY,{Y,X});			clear DfDY;
		D2fDX2 = spl_fun('fnder',pp, [2 0]);		% d^2f / dx^2
		Hxx = spl_fun('fnval',D2fDX2,{Y,X});		clear D2fDX2;
		D2fDY2 = spl_fun('fnder',pp, [0 2]);		% d^2f / dy^2
		Hyy = spl_fun('fnval',D2fDY2,{Y,X});		clear D2fDY2;
		D2fDXDY = spl_fun('fnder',pp, [1 1]);		clear pp;	% d^2f / (dx dy)
		Hxy = spl_fun('fnval',D2fDXDY,{Y,X});		clear D2fDXDY;
	% 	Z = spl_fun('fnval',pp,{Y,X});
	% 	[vy vx] = gradient_geo(Y,X,Z,'grad');
	% 	Hyy = gradient_geo(Y,X,vy,'gradN');			% d2f/dy2
	% 	Hxx = gradient_geo(Y,X,vx,'gradE');			% d2f/dx2
	% 	Hxy = 2*gradient_geo(Y,X,vx,'gradN');		% d2f/dxdy
		for (j = 1:n)
			for (i = 1:m)
				v = [vx(i,j) vy(i,j)] / norm([vx(i,j) vy(i,j)]);		% eq(2)
				R(i,j) = (v * [Hxx(i,j) Hxy(i,j); Hxy(i,j) Hyy(i,j)] * v') / (v*v');
			end
		end
		clear Hxx Hxy Hyy vx vy;
	end		% end of Tiling

	R = single(R);
	if strcmp(opt,'negative'),		R(R > 0) = 0;
	elseif strcmp(opt,'positive'),	R(R < 0) = 0;
	end
	if (handles.have_nans),			R(indNaNs) = NaN;		end		% Obey the law of NaN conservation
	if (nargout)
		sdgGrid = R;			% No worries, it's not a real copy
	else
		GRDdisplay(handles,X,Y,R,head,'SDG field','SDG field');
	end

% --------------------------------------------------------------------
function [X,Y,slope,head] = GridToolsSlope_CB(handles, opt)
% OPT == 'degrees'	Compute slope in degrees
% OPT == 'percent'	Compute slope in percentage
% OPT == 'aspect'	Compute aspect in degrees
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	[X,Y,Z,head] = load_grd(handles,'double');
	if isempty(Z),		return,		end
	set(handles.figure1,'pointer','watch')
	
	D2R = pi/180;	R2D = 180/pi;	tit = ['Slope in ' opt];
	if (handles.geog)
		if (strcmp(opt,'aspect'))
			slope = gradient_geo(Y,X,Z,'aspect');		tit = 'Terrain aspect in degrees clockwise from North';
		else
			slope = gradient_geo(Y,X,Z,'slope');
		end
	else
		if (strcmp(opt,'aspect'))
			slope = gradient_geo(Y,X,Z,'aspect','cart');	tit = 'Terrain aspect in degrees clockwise from North';
		else
			nz = getnormals(X,Y,Z);			slope = acos(nz);
		end
	end

	if ((strcmp(opt,'percent') && handles.geog)),			slope = 100 * tan(slope*D2R);
	elseif ((strcmp(opt,'percent') && ~handles.geog)),		slope = 100 * tan(slope);
	elseif ((strcmp(opt,'degrees') && ~handles.geog)),		slope = slope * R2D;
	end
	if (nargout == 0),	GRDdisplay(handles,X,Y,slope,head,tit,tit);		end

% --------------------------------------------------------------------
function GridToolsDirDerive_CB(handles, opt)
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	if (nargin == 1),	opt = 'first';		end
	luz = shading_params('dirDerivative');
	if isempty(luz),	return,		end
	azim = (90 - luz.azim) * pi/180;

	[X,Y,Z,head] = load_grd(handles,'double');
	if isempty(Z),		return,		end		% An error message was already issued

	set(handles.figure1,'pointer','watch')
	doGeog = 'geog';
	if (~handles.geog),		doGeog = 'cart';	end
	[gradN, gradE] = gradient_geo(Y,X,Z,'grad',doGeog);				% df/dy & df/dx
	if strcmp(opt,'first')
		Z = gradE * cos(azim) + gradN * sin(azim);
		str = 'First derivative';
	else		% second derivative
		Z = gradient_geo(Y,X,gradN,'gradN',doGeog)*(sin(azim)^2);	clear gradN;		% d2f/dy2
		Z = Z + gradient_geo(Y,X,gradE,'gradE',doGeog)*(cos(azim)^2);					% + d2f/dx2
		Z = Z + 2*gradient_geo(Y,X,gradE,'gradN',doGeog)*(cos(azim) * sin(azim));		% + d2f/dxdy
		str = 'Second derivative';
	end
	GRDdisplay(handles,X,Y,Z,head,str,str);

% --------------------------------------------------------------------
function out = GridToolsFindHoles_CB(handles)
% Find holes in double arrays and draw rectangles arround them
% If an output is required, don't plot and send out the result as a cell array
% where each element is a 5x2 matrix with x y coords starting at LL corner and runing CW

	if (aux_funs('msg_dlg',14,handles)),	return,		end  
	if (~handles.have_nans),	warndlg('This grid has no holes','Warning'),	return,	end
	[X,Y,Z,head,m,n] = load_grd(handles);
	if isempty(Z),		return,		end		% An error message was already issued

	set(handles.figure1,'pointer','watch'),		go_out = false;
	B = img_fun('find_holes',isnan(Z));
	if (nargout),	out = cell(numel(B),1);		go_out = true;		end
	for (i = 1:length(B))
		x_min = min(B{i}(:,2));		x_max = max(B{i}(:,2));
		y_min = min(B{i}(:,1));		y_max = max(B{i}(:,1));
		x_min = max(1,x_min-5);		x_max = min(x_max+5,n);
		y_min = max(1,y_min-5);		y_max = min(y_max+5,m);
		x_min = head(1) + (x_min-1)*head(8);	x_max = head(1) + (x_max-1)*head(8);
		y_min = head(3) + (y_min-1)*head(9);	y_max = head(3) + (y_max-1)*head(9);
		if (go_out)
			out{i} = [x_min y_min; x_min y_max; x_max y_max; x_max y_min; x_min y_min];
		else
			% Draw rectangles arround each hole
			h = line('XData',[x_min x_min x_max x_max x_min], 'YData',[y_min y_max y_max y_min y_min], ...
				'Parent',handles.axes1,'Color',handles.DefLineColor,'LineWidth',1);
			draw_funs(h,'SRTMrect')		% Set uicontexts
		end
	end
	set(handles.figure1,'pointer','arrow');

% --------------------------------------------------------------------
function GridToolsSaveAsSRTM_CB(handles)
% Only grids with the same characteristics as SRTM 3c files are allowed to be saved
	if (handles.no_file),	return,		end
	[X,Y,Z,head] = load_grd(handles);		% No need to test for in-memory Z
	if ( ((head(2)-head(1)) - 1) > 1e-6 || ((head(4)-head(3)) - 1) > 1e-6 )
		errordlg('Grid does not cover a 1 degree square','Error'),	return
	end
	if ( ((abs(head(8)-3/3600) > 1e-6) || (abs(head(9)-3/3600) > 1e-6)) && ...
			((abs(head(8)-1/3600) > 1e-6) || (abs(head(9)-1/3600) > 1e-6)) )
		errordlg('Grid spacing differs from SRTM 1 or 3 arcsec files','Error'),		return
	end
	
	% Build the file name (appended with the '_p' suffix)
	if (sign(head(1)) > 0 ),	w = 'E';
	else						w = 'W';	end
	if (sign(head(3)) > 0 ),	n = 'N';
	else						n = 'S';	end
	name = [n sprintf('%.2d',abs(round(head(3)))) w sprintf('%.3d',abs(round(head(1)))) '_p.hgt'];
	[FileName,PathName] = put_or_get_file(handles, name,'Select SRTM File name', 'put','.hgt');
	if (isequal(FileName,0)),	return,		end
	
	Z(isnan(Z)) = -32768;		Z = int16(rot90(Z,-1));		% Reset eventual NaNs to the SRTM nodata value
	fid = fopen([PathName FileName],'w','b');
	fwrite(fid,Z,'int16');		fclose(fid);
	guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function GridToolsPadd2Const_CB(handles)
% Pad the array to a const value (currently ct = zero) using a Hanning window
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	[X,Y,Z,head,m,n] = load_grd(handles);
	if isempty(Z),		return,		end		% An error message was already issued
	resp  = inputdlg({'Enter number of border lines'},'Skirt width',[1 38],{'10'});		pause(0.01)
	resp = abs( round(str2double(resp{1})) );
	if (isnan(resp)),	return,		end
	n_pad = resp * 2;
	Z = mboard(Z,n,m,n+n_pad,m+n_pad);
	zzz = grdutils(Z,'-L');  head(5) = zzz(1);  head(6) = zzz(2);	clear zzz;
	head(1) = head(1) - n_pad/2 * head(8);		head(2) = head(2) + n_pad/2 * head(8);
	head(3) = head(3) - n_pad/2 * head(9);		head(4) = head(4) + n_pad/2 * head(9);
	X = linspace(head(1),head(2),n+n_pad);		Y = linspace(head(3),head(4),m+n_pad);
	GRDdisplay(handles,X,Y,Z,head,'Padded Grid','Padded Grid')

% --------------------------------------------------------------------
% function GridToolsMesher_CB(handles)
% % 
%  	[v, f] = reduce_qslim(handles, 0.2, 'name','C:\SVN\mironeWC\triLixo.stl','format','STL','type','binary','close',-10000,'verbose','y');
% 
% % 	[v, f] = reduce_qslim(handles, 1, 'close',-5700);
% 	s.type = 'FV';	s.faces = f;	s.vertices = v;
% 	%s.track = [256510.83 4547390.15; 267558.43 4547250.09; 282801.33 4547390.15; 292170.81 4547390.15; 312587.90 4547390.15];
% 	s.hdr = handles.head;
% 	p = xyzokb_m(s,'-C2700');
% % 	p = xyzokb_m(s,'-C2700','-V','-R246819.716/318083.733/4513687.54/4575210.41', '-I1454.36769/1922.58966');
% 	mirone(p)
% 	return
% 	fout.faces = f;		fout.vertices = v;
% 	h=figure; p = patch(fout);
% 	set(p, 'facecolor', [.5 .5 .5])

% --------------------------------------------------------------------
function FileSaveFleder_CB(handles, opt)
% Depending on the OPT value, this function builds either:
% OPT = 'writeAll3' a set of three files: .geo, .dtm, .shade as DMagic would do
% OPT = 'writePlanar' directly build a planar Sonar SD file to be used by Fledermaus
% OPT = 'writeSpherical' directly build a spherical Sonar SD file to be used by Fledermaus
% OPT = 'runPlanar' build a planar .sd file (but don't keep it) and run the viewer
% OPT = 'runSpherical' build a spherical .sd file (but don't keep it) and run the viewer
	if (handles.no_file),	return,		end
	if (nargin == 1),		opt = 'runPlanar';	end
	if ( (strcmp(opt,'writeSpherical') || ~handles.flederPlanar) && ~handles.geog)
		errordlg('Spherical objects are allowed only for geographical grids','Error'),	return
	end

	fname = write_flederFiles(opt, handles);		pause(0.01);
	if (isempty(fname)),	return,		end
	if (fname(end) == 'e'),		comm = [' -scene ' fname ' &'];		% A SCENE file
	else						comm = [' -data ' fname ' &'];		% A SD file
	end
	if (strncmp(opt,'run',3))		% Run the viewer
		if (handles.whichFleder),	fcomm = ['iview4d' comm];		% Free viewer
		else						fcomm = ['fledermaus' comm];	% The real thing
		end
		try
			if (isunix)				% Stupid linux doesn't react to a non-existant iview4d
				resp = unix(fcomm);
				if (resp == 0)
					errordlg('I could not find Fledermaus. Hmmm, do you have it?','Error')
				end
			elseif (ispc)
				s = dos(fcomm);
				% Try again with the 'iview3d'
				if (~s && handles.whichFleder),		fcomm(6) = '3';		dos(fcomm);		end
			else
				errordlg('Unknown platform.','Error'),	return
			end
		catch
			errordlg('I could not find Fledermaus. Hmmm, do you have it?','Error')
		end
		pause(1)
		%builtin('delete',fname);
	end

% --------------------------------------------------------------------
function ImageEdgeDetect_CB(handles, opt)
% OPT	-> 'ppa'     Find ridges
%		-> 'apalpa'  Get the polygons that sorround good data (limit NaN-noNaN)
%		-> 'Vec'
%		-> 'Lines'
%		-> 'Circles' Find circles in image (uses OpenCV and fails a lot)
%		-> 'Rect'    Find rectangles in image (uses OpenCV and fails a lot)
%		-> 'Ras'
%		-> 'SUS'
if (handles.no_file),		return,		end

img = get(handles.hImg,'CData');		% Even when not used, this op cost nothing
dims = [size(img,1) size(img,2)];
set(handles.figure1,'pointer','watch')
if (strcmp(opt,'ppa'))
	[X,Y,Z,handles.head] = load_grd(handles);
	if isempty(Z),		return,		end
	out = grdppa_m(Z,handles.head);
	h_ridge = line(out(1,1:end),out(2,1:end), 'Linewidth',handles.DefLineThick, ...
		'Color',handles.DefLineColor,'Tag','creast_line','Userdata',1);
	multi_segs_str = cell(length(h_ridge),1);	% Just create a set of empty info strings
	draw_funs(h_ridge,'isochron',multi_segs_str)
	set(handles.figure1,'pointer','arrow')
	return
% 	x_lim = get(handles.axes1,'XLim');	y_lim = get(handles.axes1,'YLim');
% 	h_lixo = figure('MenuBar','none');
% 	h_tmp = line('XData',out(1,:),'YData',out(2,:),'Linewidth',0.1);
% 	set(handles.axes1,'XLim',x_lim,'YLim',y_lim)
% 	F = getframe(h_lixo);
% 	img = F.cdata;
% 	[m,n,k] = size(img);
% 	x_inc = (handles.head(2) - handles.head(1)) / n;
% 	y_inc = (handles.head(4) - handles.head(3)) / m;
% 	I = flipud(img(:,:,1));
% 	delete(h_lixo);
elseif (strcmp(opt,'apalpa'))		% Get the polygons that sorround good data
	if (~handles.have_nans)
		warndlg('There are no wholes in this grid.','Chico Clever'),	return
	end
	[X,Y,Z] = load_grd(handles);
	if isempty(Z),		return,		end
	mask = isnan(Z);
	B = img_fun('bwboundaries',mask, 8, 'holes');
	%B = bwbound_unique(B);
	dims = [size(Z,1) size(Z,2)];
	opt = 'Vec';					% Trick to not need to add an extra case in the IF test below

elseif (strcmp(opt,'Vec') || strcmp(opt,'Ras') || strcmp(opt(1:3),'SUS'))
	if (ndims(img) == 3),		img = cvlib_mex('color',img,'rgb2gray');	end
	if (~strcmp(opt(1:3),'SUS'))
		if (~handles.IamCompiled),		img = canny(img);		% Economic version (uses singles & cvlib_mex but crashs in compiled)
		else							img = img_fun('edge',img,'canny');
		end
	else
		img = susan(img,'-e');				% Do SUSAN edge detect
		if (strcmp(opt(4:end),'vec')),	opt = 'Vec';		% This avoids some extra tests later
		else							opt = 'Ras';
		end
	end
	if (strcmp(opt,'Vec'))
		B = img_fun('bwboundaries',img,'noholes');
		B = bwbound_unique(B);
	end
elseif (strcmp(opt,'Lines'))
	%B = cvlib_mex('houghlines2',img,'standard',1,pi/180,100,0,0);
	%BW = cvlib_mex('canny',img,40,200,3);
	if (ndims(img) == 3),			img = cvlib_mex('color',img,'rgb2gray');	end
	if (~handles.IamCompiled),		BW = canny(img);
	else							BW = img_fun('edge',img,'canny');
	end

	[H,T,R] = img_fun('hough',BW);
	P = img_fun('houghpeaks',H,50,'threshold',ceil(0.3*double(max(H(:)))));
	lines = img_fun('houghlines',BW,T,R,P,'FillGap',10,'MinLength',50);
	if (~isfield(lines,'point1')),	set(handles.figure1,'pointer','arrow'),		return,		end
	B = cell(1,length(lines));
	for (k = 1:length(lines))
		B{k} = [lines(k).point1(2:-1:1); lines(k).point2(2:-1:1)];
	end
elseif (strcmp(opt,'Circles'))
	B = cvlib_mex('houghcircles',img);
	if (isempty(B)),	set(handles.figure1,'pointer','arrow'),		return,		end
elseif (strcmp(opt,'Rect'))
	B = cvlib_mex('findrect',img);
	if (isempty(B)),	set(handles.figure1,'pointer','arrow'),		return,		end
end

if (strcmp(opt,'Vec') || strcmp(opt,'Lines') || strcmp(opt,'Rect'))		% Convert the edges found into vector
	x_inc = handles.head(8);		y_inc = handles.head(9);
	x_min = handles.head(1);		y_min = handles.head(3);
	if (handles.head(7))			% Work in pixel registration
		x_min = x_min + x_inc/2;	y_min = y_min + y_inc/2;
	end
	h_edge = zeros(length(B),1);	i = 1;
	for k = 1:length(B)
		bnd = B{k};
		if (strcmp(opt,'Vec') && mask(bnd(1,1),bnd(1,2)) && ...	% Drastic but no better solution to avoid NaNs-in-corners cases
				(min(bnd(:,1)) == 1 || max(bnd(:,1)) == dims(2) || min(bnd(:,1)) == 1 || max(bnd(:,1)) == dims(1)) )
			continue
		end
		if (size(bnd,1) > 4)
			bnd = cvlib_mex('dp', bnd, 0.1);		% Simplify line but only on straight lines
		end
		% Some times we get a BB rectangle, test and ignore it if it's the case
		if (size(bnd,1) < 10)
			ma = max(bnd);
			if (all(abs(min(bnd) - 1)) < 1e-3)
				if ((abs(ma(1) - dims(1)) < 1e-3) && (abs(ma(2) - dims(2)) < 1e-3))
					continue
				end
			end
		end

		y = (bnd(:,1)-1)*y_inc + y_min;
		x = (bnd(:,2)-1)*x_inc + x_min;
		h_edge(i) = line('XData',x, 'YData',y, 'Parent',handles.axes1,'Linewidth',handles.DefLineThick, ...
			'Color',handles.DefLineColor,'Tag','edge_detected','Userdata',i);
		i = i + 1;
		%ellipse_t = fit_ellipse( x,y,handles.axes1 );
	end
	h_edge(h_edge == 0) = [];					% Remove empty handles remaining from pre-declaration
	multi_segs_str = cell(length(h_edge),1);	% Just create a set of empty info strings
	draw_funs(h_edge,'isochron',multi_segs_str);
elseif (strcmp(opt,'Circles'))
	x = linspace(-pi,pi,360);		y = x;
	xx = cos(x);					yy = sin(y);
	%h_circ = line('XData', [], 'YData', []);
	for k = 1:size(B,1)
		x = B(k,1) + B(k,3) * xx;			y = B(k,2) + B(k,3) * yy;
		h_circ = line('XData',x, 'YData',y, 'Parent',handles.axes1, 'Linewidth',handles.DefLineThick, ...
				'Color',handles.DefLineColor,'Userdata',B(k,1:end));
		draw_funs(h_circ,'SessionRestoreCircleCart')	% Give uicontext
		setappdata(h_circ,'LonLatRad',B(k,1:end))
		%setappdata(h_circ,'X',xx);		setappdata(h_circ,'Y',yy);
	end
else							% Display the bw image where the edges are the whites
	setappdata(0,'CropedColormap',gray);
	if (handles.image_type == 2)
		mirone(img);
	else
		tmp.X = handles.head(1:2);	tmp.Y = handles.head(3:4);		tmp.name = 'Edge detected';
		tmp.geog = handles.geog;	tmp.head = [handles.head(1:4) 0 1 handles.head(7:9)];
		mirone(img, tmp);
	end
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DigitalFilt_CB(handles, opt)
	if (handles.no_file),	return,		end
	if (strcmp(opt,'image'))
		digitalFiltering(handles.hImg);
	else		% grid
		[X,Y,Z,handles.head] = load_grd(handles);	% load the grid array here
		if isempty(Z),		return,		end			% An error message was already issued
		[Z, img] = digitalFiltering(handles.hImg,Z,get(handles.figure1,'ColorMap'));
		if (isempty(Z)),	return,		end

		[zzz] = grdutils(Z,'-L');	handles.head(5) = zzz(1);	handles.head(6) = zzz(2);
		setappdata(handles.figure1,'dem_z',Z);
		handles.origFig = img;		guidata(handles.figure1,handles)
	end

% --------------------------------------------------------------------
function RotateTool_CB(handles, opt)
	if (handles.no_file),	return,		end
	if (strcmp(opt,'image'))
		img = rotatetool(get(handles.hImg,'CData'),handles);
		if (ndims(img) == 2),	setappdata(0,'CropedColormap',get(handles.figure1,'ColorMap')),		end
		if(strcmp(get(handles.axes1,'YDir'),'normal')),		img = flipdim(img,1);	end
		if (~isempty(img)),		mirone(img);	end
	else		% grid
		[X,Y,Z] = load_grd(handles);			% load the grid array here
		if isempty(Z),		return,		end		% An error message was already issued
		[newZ, hdr] = rotatetool(Z,handles);
		if (isempty(newZ)),	return,		end
		[ny,nx] = size(newZ);
		X = linspace(hdr(1),hdr(2),nx);		Y = linspace(hdr(3),hdr(4),ny);
		GRDdisplay(handles,X,Y,newZ,hdr,'Rotated grid','Rotated grid')
	end

% --------------------------------------------------------------------
function [Z, indNaNs] = fillGridGaps(handles, Z)
% Inspects the array Z for blobs of NaNs and fill them with minimum curvature interpolation.

	holes = GridToolsFindHoles_CB(handles);
	head = handles.head;
	if (nargout == 2)
		indNaNs = isnan(Z);			% Usefull for restore the NaNs at a later stage
	end
	for (k = 1:numel(holes))
		rect = [holes{k}(1,1) holes{k}(1,2) (holes{k}(3,1) - holes{k}(2,1)) (holes{k}(2,2) - holes{k}(1,2))]; 
		[Z_rect,r_c] = cropimg(head(1:2), head(3:4), Z, rect, 'out_grid');
		XX = linspace( head(1) + (r_c(3)-1)*head(8), head(1) + (r_c(4)-1)*head(8), r_c(4) - r_c(3) + 1 );
		YY = linspace( head(3) + (r_c(1)-1)*head(9), head(3) + (r_c(2)-1)*head(9), r_c(2) - r_c(1) + 1 );
		opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', XX(1), XX(end), YY(1), YY(end));
		opt_I = sprintf('-I%.10f/%.10f',head(8),head(9));
		Z_rect = Z_rect(:);			% 
		ind = ~isnan(Z_rect);		% Finding NaNs in singles should be a little cheaper
		Z_rect = double(Z_rect);	% It has to be
		Z_rect = Z_rect(ind);
		[XX,YY] = meshgrid(XX,YY);
		XX = XX(ind);		YY = YY(ind);
		Z_rect = gmtmbgrid_m(XX, YY, Z_rect, opt_R, opt_I, '-Mz');
		Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
	end

% --------------------------------------------------------------------
function line_levelling(handles)
% A grid based method to do tie-line leveling of magnetic and/or gravity data.
% Use the loaded grid, load a x,y,z data file that can be multi-seg (normally the same used to compute the grid)
% and compute and apply the vertical offset necessary to bring the profile's best fit model to the
% the best fit surface to the loaded grid. The result is stored in a new file with '_leveled' appended to data's name

	[FileName, pato, handles] = put_or_get_file(handles, ...
			{'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select File','get');
	if (isequal(FileName,0)),		return,		end
	fname = [pato FileName];
	[numeric_data, ms_str] = load_xyz([], fname);
	if (isempty(numeric_data)),		return,	end
	if (~isa(numeric_data, 'cell'))
		numeric_data = {numeric_data};		% To go along with the case of a multi-segment file
	end
	[X,Y,Z] = load_grd(handles);			% load the grid array here
	if isempty(Z),		return,		end		% An error message was already issued
	opt_N = '-Nr6';							% Fit a quadratic surface
	newZ = c_grdtrend(Z,handles.head,'-T',opt_N);
	newZ = c_grdtrend(newZ,handles.head,'-D','-N3');
	aux_funs('StoreZ',handles,X,Y,newZ)		% Need to do this because ExtractProfile_CB() will get the grid from there
	[lix,FileName,EXT] = fileparts(FileName);
	fid = fopen([pato FileName '_leveled' EXT],'wt');
	for (k = 1:numel(numeric_data))
		data = numeric_data{k}(:,1:3);
		ll = draw_funs([], 'show_LineLength', [], [], data(:,1:2), 'k');
		rP = [0; cumsum(ll.seg_len)];		% The profile accumulated distance
		zP_fit = trend1d_m([rP data(:,3)],'-Fm','-N4r');	% z of the profile fit
		track  = ExtractProfile_CB(handles, 'input', data(:,1:2));
		diffa  = zP_fit - track.z(:);		% Difference between profile fit and grdtrack on surface fit
		data(:,3) = data(:,3) - diffa;		% The z-leveled profile
		fprintf(fid, '%s\n', ms_str{k});
		fprintf(fid, '%.5f\t%.5f\t%.2f\n', data');
	end
	fclose(fid);
	aux_funs('StoreZ',handles,X,Y,Z)		% Restore original array

% --------------------------------------------------------------------
function TransferB_CB(handles, opt, opt2)
% OPT2 is only used by custom_menus() and holds the directory name to start in
	if (strcmp(opt,'guessType'))
		str = {'*.grd;*.nc;*.tif;*.tiff;*.jpg;*.jp2;*.png;*.gif;*.mat;*.cpt;*.hdf;*.img', ...
				'Files (*.grd,*.nc,*.tif,*.tiff,*.jpg,*.jp2,*.png,*.gif,*.mat,*.cpt,*.hdf,*.img)'; '*.*', 'All Files (*.*)'};
		if (nargin == 3 && ~isempty(opt2)),		handles.last_dir = opt2;		end
		[FileName,PathName,handles] = put_or_get_file(handles,str,'Select file','get');
		if (isequal(FileName,0)),	return,		end				% User gave up
		drv = aux_funs('findFileType',[PathName FileName]);
		if (~isempty(drv)),		gateLoadFile(handles,drv,[PathName FileName]);  end

	elseif (strcmp(opt,'BgMap'))
		out = bg_map(handles.path_data);
		if isempty(out),	return,		end		% User gave up loading the fig tille
		handles.geog = 1;				handles.image_type = 3;		handles.head = out.head;
		show_image(handles,out.imgName,out.X,out.Y,out.img,0,'xy',1,1);

	elseif (strcmp(opt,'NewEmpty'))
		h = mirone;
		newHand = guidata(h);		newHand.last_dir = handles.last_dir;	guidata(h, newHand)
		setappdata(h,'hFigParent',handles.figure1);		% Save the Parent fig handles in this new figure
		if (~handles.no_file)					% Set the Drape option to 'on' in the New window
			set(newHand.ImageDrape,'Vis','on')
		end
		if (handles.validGrid)
			set(newHand.RetroShade,'Vis','on')	% Set the Retro-Illum option to 'on' in the New window
		end
		thematic_pal = getappdata(handles.figure1, 'thematic_pal');		% Because of the bloody cpt2cmap crash
		if (~isempty(thematic_pal)),	setappdata(newHand.figure1,'thematic_pal',thematic_pal),	end

	elseif (strcmp(opt,'scale'))				% Apply a scale factor
		resp = inputdlg({'Enter scale factor'},'Rescale grid',[1 30],{'-1'});	pause(0.01)
		if (isempty(resp)),		return,	end
		scal = str2double(resp{1});
		if (isnan(scal)),	return,		end
		[X,Y,Z] = load_grd(handles);			% load the grid array here
		if isempty(Z),		return,		end		% An error message was already issued
		Z = cvlib_mex('CvtScale',Z, scal);
		head = [handles.head(1:4) handles.head(5:6)*scal handles.head(7:end)];
		if (handles.head(5) > handles.head(6))	% Happens when we change sign
			tmp = handles.head(5);	handles.head(5) = handles.head(6);	handles.head(6) = tmp;
		end
		tmp = struct('X',X, 'Y',Y, 'head',head, 'geog',handles.geog, 'name','Scaled grid');
		projWKT = getappdata(handles.figure1,'ProjWKT');
		if (~isempty(projWKT)),		tmp.srsWKT = projWKT;	end
		mirone(Z, tmp)

	elseif (strcmp(opt,'sinks'))				% Fill sinks
		[X,Y,Z] = load_grd(handles);			% load the grid array here
		if isempty(Z),		return,		end
		Z = img_fun('imfill', Z, 'holes');
		head = [handles.head(1:4) min(Z(:)) handles.head(6) handles.head(7:end)];
		tmp = struct('X',X, 'Y',Y, 'head',head, 'geog',handles.geog, 'name','Sinks filled grid');
		projWKT = getappdata(handles.figure1,'ProjWKT');
		if (~isempty(projWKT)),		tmp.srsWKT = projWKT;	end
		mirone(Z, tmp)

	elseif (strncmp(opt,'Multiscale',5))		% Apply a block-wise operation
		[X,Y,Z] = load_grd(handles);			% load the grid array here
		if isempty(Z),		return,		end
		resp = multiscale;
		if (isempty(resp)),		return,		end
		opt_A = sprintf('-A%d', resp.method);
		opt_W = sprintf('-W%d', resp.size);
		opt_N = sprintf('-N%d', handles.have_nans);
		opt_G = ' ';
		if (handles.geog),	opt_G = '-G';	end

		Z = mirblock(Z, handles.head, opt_A, opt_N, opt_W, opt_G);	% MUST flexiblize mirblock to acept 10 hdr elements
		ind = find(isinf(Z(:)));
		if (~isempty(ind)),	Z(ind) = NaN;	end
		projWKT = getappdata(handles.figure1,'ProjWKT');
		if (strcmp(resp.name, 'Slope')),	resp.name = 'Slope in degrees';		end
		GRDdisplay(handles,X,Y,Z,handles.head,[],resp.name, projWKT);

 	elseif (strcmp(opt,'dump'))					% Show the RAM fragmentation (Windows only)
		dumpmemmex

 	elseif (strcmp(opt,'lasterr'))				% Show last error (standalone only and for debug)
		lstErr = lasterror;
		msgbox(sprintf('Last error message in stack is:\n\n%s\n%s\n',lstErr.message, lstErr.identifier),'Debug message')

 	elseif (strcmp(opt,'sharedir'))				% Show GMT_SHAREDIR env var (standalone only and for debug)
		env4 = set_gmt('GMT_SHAREDIR','lixo');
		env5 = set_gmt('GMT5_SHAREDIR','lixo');
		if (~isempty(env4) && ~isempty(env5))
			msgbox(sprintf('GMT_SHAREDIR = %s\nGMT5_SHAREDIR = %s', env4, env5))
		elseif (~isempty(env4))
			msgbox(sprintf('GMT_SHAREDIR = %s', env4))
		elseif (~isempty(env5))
			msgbox(sprintf('GMT5_SHAREDIR = %s', env5))
		else
			warndlg('No GMT_SHAREDIR variables found. Most likely you do not have GMT installed.', 'Warning')
		end

 	elseif (strncmp(opt,'TransplantGrid',6))	%
		ImageCrop_CB(handles,zeros(5,2),'ImplantGrid')	% The '0' is only to avoid the call to ruberbandbox in ImageCrop_CB

 	elseif (strcmp(opt,'fract'))				% Fractal surf. Have to do it here due to dumb compiler limitations
		if (handles.no_file)					% When called from a virgin figure, use it
			[Z, hdrStruct] = gen_UMF2d;
			handles.computed_grid = 1;			handles.image_type = 1;
			hdrStruct.Z = Z;
			loadGRID(handles,'','IN',hdrStruct)
		else									% For non-virgin figs, create a new one with the result.
			gen_UMF2d;
		end

 	elseif (strcmp(opt,'update'))				% Update via Web the stand-alone version
		dest_fiche = [handles.path_tmp 'apudeita.txt'];		url = 'w3.ualg.pt/~jluis/mirone/updates/v27/';
		dos(['wget "' url 'apudeita.txt' '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
		finfo = dir(dest_fiche);
		if (finfo.bytes == 0)
			builtin('delete',dest_fiche);
			msgbox('This Mirone version is updated to latest.','Nothing New1'),	return
		end
		fid = fopen(dest_fiche,'rt');
		todos = fread(fid,'*char');		fclose(fid);
		[nomes, MD5, V.Vstr] = strread(todos,'%s %s %s');	% In future we will have a use for the version string
		builtin('delete',dest_fiche);	n = 1;		% Remove this one right away
		namedl = cell(1);	MD5_check = cell(1);	% Mostly to shutup MLint
		pato_file = cell(numel(nomes),1);
		ind_all = false(numel(nomes),1);			% To flag the ones truely to be updated later
		for (k = 1:numel(nomes))
			[pato, nome, ext] = fileparts(nomes{k});
			pato_file{k} = pato;					% Store path files to see if that dir must be created by the bat file
			if (exist(nomes{k}, 'file'))			% File exists localy and it's a potential target for update
				localMD5 = CalcMD5(nomes{k},'file');
				if (~strcmp(MD5{k}, localMD5))
					namedl{n} = [url nome ext];		% File name to update with path realtive to Mir root
					MD5_check{n} = MD5{k};
					n = n + 1;		ind_all(k) = true;
				end
			else									% New file. Download for sure.
				namedl{n} = [url nome ext];
				MD5_check{n} = MD5{k};
				n = n + 1;			ind_all(k) = true;
			end
		end
		if (n == 1)									% Nothing new to update
			builtin('delete',dest_fiche);
			msgbox('This Mirone version is updated to latest.','Nothing New2'),	return
		end
		nomes = nomes(ind_all);						% Retain only the original names (that may include a path) to dl
		pato_file = pato_file(ind_all);

		ind = false(1,n-1);			msg = [];
		for (k = 1:n-1)
			[pato, nome, ext] = fileparts(namedl{k});		dest_fiche = [handles.path_tmp nome ext];
			dos(['wget "' url nome ext '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
			localMD5 = CalcMD5(dest_fiche, 'file');
			if (~strcmp(MD5_check{k}, localMD5))	% Troubles in transmission
				ind(k) = true;
			end
		end
		if (any(ind))								% If some file was lost report it and remove them from list
			namedl(~ind) = [];		nomes(~ind) = [];
			msg = cell(numel(find(ind))+1,1);
			msg{1} = 'Failed to download these files:';		msg(2:end) = namedl(ind);
		elseif (all(ind))
			warndlg(sprintf('Failed to download all %d files. Try again.', n-1),'Warning'),		return
		end
		if (~isempty(msg))
			msg{end+1} = '';	msg{end+1} = 'Only the downloaded files will be updated on NEXT Mirone start';
		else
			msg = 'Successfully downloaded all necessary files. They will be updated on NEXT Mirone start';
		end
		msgbox(msg, 'Finished download')
		fid = fopen([handles.path_tmp 'apudeita.bat'],'wt');	% Create the updating batch that will be run by callMir
		fprintf(fid, '@echo off\nREM copy updated files from tmp and place into their destination\n');
		for (k = 1:numel(namedl))
			[pato, nome, ext] = fileparts(namedl{k});
			if (~isempty(pato_file{k}))
				fprintf(fid, 'IF NOT EXIST %s md %s\n', pato_file{k}, pato_file{k});
			end
			fprintf(fid, 'move /Y tmp\\%s\t.\\%s\n', [nome ext], nomes{k});
		end
		fprintf(fid, 'echo Ja ta. Finished update\n');
		fclose(fid);

	elseif (strcmp(opt,'DayNight'))
% 		[sun_params, lon, lat] = solar_params(-7.92, 37.073);
% 		h = patch('XData',lon, 'YData',lat,'FaceColor','none','Parent',handles.axes1,'Tag','DayNight');
% 		draw_funs(h,'line_uicontext')
% 		img = gdalread([handles.path_data 'etopo4.jpg'], '-U');
% 		tmp.X = [-180 180];		tmp.Y = [-90 90];
% 		x_inc = diff(tmp.X) / (size(img,2)-1);
% 		y_inc = diff(tmp.Y) / (size(img,1)-1);
% 		tmp.head = [tmp.X tmp.Y 0 255 0 x_inc y_inc];
% 		h = mirone(img, tmp);
% 		handles.geog = 1;		handles.image_type = 3;		handles.head = out.head;
% 		show_image(handles,out.imgName,out.X,out.Y,out.img,0,'xy',1,1);

%  	elseif (strncmp(opt,'hydro',5))					% ...
% 		if (strcmp(opt(7:end), 'flow'))
% 			hydrology('entry', handles)
% 		elseif (strcmp(opt(7:end), 'basin'))
% 			hydrology('picdrain', handles)
% 		end
	end

% --------------------------------------------------------------------
function Transfer_CB(handles, opt)
	if (handles.no_file),		return,		end

	set(handles.figure1,'pointer','watch')
	img = get(handles.hImg,'CData');
	if (strcmp(opt,'Corners'))
		if (islogical(img)),	img = uint8(img);	end		% cvlib_mex should take care of this case
		corn = cvlib_mex('goodfeatures',img,100,0.05);
		if (handles.image_type ~= 2)
			y = (corn(:,1)-1)*handles.head(9) + handles.head(3);
			x = (corn(:,2)-1)*handles.head(8) + handles.head(1);
		else
			x = corn(:,1);		y = corn(:,2); 
		end
		hLine = line('XData',x,'YData',y,'Parent',handles.axes1,'LineStyle','none','Marker','o', ...
					'MarkerEdgeColor','w','MarkerFaceColor','k','MarkerSize',6,'Tag','corner_detected','Userdata',1);
		multi_segs_str = cell(length(hLine),1);		% Just create a set of empty info strings
		draw_funs(hLine,'isochron',multi_segs_str);

	elseif (strcmp(opt,'toRGB'))
		if (ndims(img) == 3 || isa(img,'logical')),	set(handles.figure1,'pointer','arrow'),		return,	end		% Nothing to do
		img = ind2rgb8(img,get(handles.figure1,'Colormap'));
		set(handles.hImg,'CData', img)
		aux_funs('togCheck', handles.ImModRGB, [handles.ImMod8cor handles.ImMod8gray handles.ImModBW handles.ImModNeg handles.ImModNegBW])

	elseif (strcmp(opt,'8-bit'))
		if (ndims(img) ~= 3),	set(handles.figure1,'pointer','arrow'),		return,		end		% Nothing to do
		resp  = inputdlg({'Number of colors (2-256)'},'Color quantization',[1 30],{'256'});	pause(0.01)
		nColors = round (abs(str2double(resp{1})) );
		if (isnan(nColors)),	set(handles.figure1,'pointer','arrow'),		return,		end
		[img, map] = img_fun( 'rgb2ind', img, max(2, min(nColors, 256)) );	% Ensure we are in the [2-256] int
		set(handles.hImg,'CData', img),			set(handles.figure1,'ColorMap',map)
		aux_funs('togCheck', handles.ImMod8cor, [handles.ImMod8gray handles.ImModBW handles.ImModRGB handles.ImModNeg handles.ImModNegBW])

	elseif (strcmp(opt,'gray'))
		if (ndims(img) == 3)
			img = cvlib_mex('color',img,'rgb2gray');	set(handles.hImg,'CData', img);
		end
		set(handles.figure1,'ColorMap',gray(256))
		aux_funs('togCheck',handles.ImMod8gray, [handles.ImMod8cor handles.ImModBW handles.ImModRGB handles.ImModNeg handles.ImModNegBW])

	elseif (strcmp(opt,'bw'))
		img = img_fun('im2bw',img);
		set(handles.hImg,'CData', img, 'CDataMapping','scaled');
		set(handles.figure1,'ColorMap',gray(256))
		aux_funs('togCheck',handles.ImModBW, [handles.ImMod8cor handles.ImMod8gray handles.ImModRGB handles.ImModNeg handles.ImModNegBW])

	elseif (strcmp(opt,'neg'))
		img = img_fun('imcomplement',img);
		set(handles.hImg,'CData', img)
		aux_funs('togCheck',handles.ImModNeg, [handles.ImMod8cor handles.ImMod8gray handles.ImModRGB handles.ImModBW handles.ImModNegBW])

	elseif (strcmp(opt,'negBW'))
		img = img_fun('im2bw',img_fun('imcomplement',img));
		set(handles.hImg,'CData', img, 'CDataMapping','scaled');
		set(handles.figure1,'ColorMap',gray(256))
		aux_funs('togCheck',handles.ImModNegBW, [handles.ImMod8cor handles.ImMod8gray handles.ImModRGB handles.ImModBW handles.ImModNeg])

	elseif (strcmp(opt,'RGBexp'))
		if (ndims(img) < 3),	return,		end		% Should never happen because non-RGB must have this opt hiden
		P{1} = cvlib_mex('color',img,'rgb2gray');
		P{2} = img(:,:,1);		P{3} = img(:,:,2);		P{4} = img(:,:,3);
		B = cvlib_mex('color',img,'rgb2hsv');
		P{5} = B(:,:,1);		P{6} = B(:,:,2);		P{7} = B(:,:,3);
		B = cvlib_mex('color',img,'rgb2YCrCb');
		P{8} = B(:,:,1);		P{9} = B(:,:,2);		P{10} = B(:,:,3);
		B = cvlib_mex('color',img,'rgb2lab');
		P{11} = B(:,:,1);		P{12} = B(:,:,2);		P{13} = B(:,:,3);
		if (handles.image_type == 2)
			hAx = montage(P);
		else
			hAx = montage(P,'handParent',handles.figure1, 'flipud');
		end
		set(get(hAx(1),'Parent'), 'Name', 'Explore RGB')
		labels = {'Gray' 'Red' 'Green' 'Blue' 'Hue' 'Saturation' 'Value' 'Luminance' 'Red Chrominance' ...
			'Blue Chrominance' 'L*a*b* => L*' 'L*a*b* => a*' 'L*a*b* => b*'};
		cores = {'k' [0.7 0 0] [0 0.7 0] 'b' [1 0 1]*0.6 [1 0 1]*0.6 [1 0 1]*0.6 ...
			[0 1 1]*0.6 [0 1 1]*0.6 [0 1 1]*0.6 [1 1 0]*0.6 [1 1 0]*0.6 [1 1 0]*0.6};
		for (k = 1:13)		% Add labels to individual component images
			text(10, 50, labels{k}, 'Parent', hAx(k), 'FontSize',16, 'FontWeight','bold', 'Color',cores{k}, 'Tag','cor');
		end

	elseif (strcmp(opt,'copyclip'))		% Img and frame capture to ClipBoard
		h = getappdata(handles.figure1,'CoordsStBar');		set(h,'Visible','off');
		imcapture(handles.axes1,'imgAx');					set(h(2:end),'Visible','on')

	elseif (strcmp(opt,'Ctrl-c'))
		h_active = getappdata(handles.figure1,'epActivHand');
		if (h_active)					% We have a line or patch in edit mode. Copy it
			x = get(h_active,'xdata');		y = get(h_active,'ydata');		z = getappdata(h_active,'ZData');
			if (isempty(z)),		mat2clip([x(:) y(:)],8)
			else					mat2clip([x(:) y(:) z(:)],8)
			end
			setappdata(0, 'CtrlCHandle', [h_active handles.axes1])	% Put a handle copy on root's appdata
		end

	elseif (strcmp(opt,'Ctrl-v'))
		h = getappdata(0, 'CtrlCHandle');	% Get what's in this root's appdata
		if (isempty(h) || ~ishandle(h(1)))
			set(handles.figure1,'pointer','arrow'),		return
		end
		draw_funs(h(1), 'Ctrl_v', [h(1) handles.axes1])

% 	elseif (strcmp(opt,'Ctrl-p'))
% 		uiinspect(gco)

	elseif (strncmp(opt,'flip',4))		% LR or UP image flipage. OPT = flipLR or flipUD
		% OPT == 'LR' -> Flips the image left-right. OPT == 'UD' -> Flips the image up-down
		direction = 1;
		if strcmp(opt(5:6),'LR'),	direction = 2;		end		% Flip left-right
		set(handles.hImg,'CData', flipdim(img,direction))

	elseif (strcmp(opt,'KML'))
		[FileName,PathName] = put_or_get_file(handles,{'*.kml', 'KML files (*.kml)'},'Select file','put','.kml');
		if isequal(FileName,0),		set(handles.figure1,'pointer','arrow'),		return,		end			% User gave up
		Z = [];
		if (handles.have_nans)
			[X,Y,Z] = load_grd(handles);
			if isempty(Z),		set(handles.figure1,'pointer','arrow'),		return,		end
		end
		writekml(handles,Z,[PathName FileName])		% Z will be used to setup a alpha channel

	elseif (strncmp(opt,'morph',5))				% Works for either image or grids
		set(handles.figure1,'pointer','arrow')
		structuring_elem(handles.figure1, strcmp(opt(7:end), 'grd'))

	elseif (strcmp(opt,'scatter'))
		str = {'*.dat;*.DAT;*.txt;*.TXT', 'Data file (*.dat,*.DAT,*.txt,*.TXT)'; '*.*', 'All Files (*.*)'};
		[FileName,PathName] = put_or_get_file(handles,str,'Select file','get');
		if isequal(FileName,0),		set(handles.figure1,'pointer','arrow'),		return,		end		% User gave up
		scatter_plot(handles,[PathName,FileName]);

	elseif (strcmp(opt,'print'))
		h = findobj('Type','uicontrol');		set(h,'Visible','off')	% We don't want to print the buttons
		handsStBar = getappdata(handles.figure1,'CoordsStBar');
		set(handsStBar,'Visible','off');		set(handles.figure1,'pointer','arrow')
		ff = findobj(handles.figure1, '-depth',1, 'Style','frame', 'Tag', 'FancyFrame');
		if (~isempty(ff)),		fancyFrame(guidata(handles.figure1),'pset');	end
		if (ispc),		print -v,		else	print;  end
		if (~isempty(ff)),		fancyFrame(handles,'punset');	end
		set(h,'Visible','on');  set(handsStBar(2:end),'Visible','on');

	end
	set(handles.figure1,'pointer','arrow')

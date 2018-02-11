function varargout = plot_composer(varargin)
% Helper window to generate a GMT script that reproduces the Mirone's figure contents

%	Copyright (c) 2004-2018 by J. Luis
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

% $Id: plot_composer.m 10267 2018-02-11 00:43:03Z j $

	handMir = varargin{1};
	if (handMir.no_file)     % Stupid call with nothing loaded on the Mirone window
		return
	end

	hObject = figure('Vis','off');
	plot_composer_LayoutFcn(hObject);
	move2side(hObject,'right');

	load([handMir.path_data 'mirone_icons.mat']);
	hTB = uitoolbar('parent',hObject, 'BusyAction','queue','HandleVisibility','on','Interruptible','on',...
		'Tag','FigureToolBar','Vis','on');
	uipushtool('parent',hTB,'Click','mirone(''Draw_CB'',guidata(gcbo),''Text'')', ...
		'Tag','DrawText','cdata',text_ico,'Tooltip','Insert Text','Sep','on');
	uipushtool('parent',hTB,'Click','mirone(''DrawLine_CB'',guidata(gcbo))', ...
		'Tag','DrawLine','cdata',Mline_ico,'Tooltip','Draw Line');
	uipushtool('parent',hTB,'Click',@draw_rectangle, ...
		'Tag','DrawRect','cdata',rectang_ico,'Tooltip','Draw Rectangle');
	uipushtool('parent',hTB,'Click','mirone(''DrawClosedPolygon_CB'',guidata(gcbo),[])', ...
		'Tag','DrawPolyg','cdata',polygon_ico,'Tooltip','Draw Closed Polygon');
	uipushtool('parent',hTB,'Click','mirone(''Draw_CB'',guidata(gcbo),''Vector'')', ...
		'Tag','DrawArrow','cdata',Marrow_ico,'Tooltip','Draw Arrow');
	uipushtool('parent',hTB,'Click','mirone(''DrawGeogCircle_CB'',guidata(gcbo),''cartCirc'')', ...
		'Tag','DrawGeogCirc','cdata',circ_ico,'Tooltip','Draw circle');
	uitoggletool('parent',hTB,'Click','mirone(''PanZoom_CB'',guidata(gcbo),gcbo,''zoom'')', ...
		'Tag','Zoom','cdata',zoom_ico,'Tooltip','Zooming on/off','Sep','on');

	handles = guihandles(hObject);

	hFigs = findobj(0,'type','figure');		% Fish all figures, but not the one that we are about to create
	hAllFigs = handMir.figure1;
	hFigs = aux_funs('figs_XOR', handles.figure1, hFigs);	% Get all unique Mirone Figs
	if (~isempty(hFigs))
		if (numel(hFigs) > 1)	% The problem here is that hFigs is sorted and I need the handMir.figure1 at the top
			hFigs(hFigs == handMir.figure1) = [];
			hFigs = [handMir.figure1; hFigs(:)];
		end
		hAllFigs = hFigs;
	end

	% Some of these Mirone Figs may be empty. We don't care about them, so wipe them out from list.
	c = false(1,numel(hAllFigs));
	for (k = 1:numel(hAllFigs))
		handTmp = guidata(hAllFigs(k));
		if (handTmp.no_file),	c(k) = true;	end
	end
	hAllFigs(c) = [];

	N_figs = numel(hAllFigs);
	nomes = get(hAllFigs,'name');
	if (~isa(nomes,'cell')),	nomes = {nomes};	end
	for (k = 1:N_figs)
		[pato, nomes{k}] = fileparts(nomes{k});
		ind = strfind(nomes{k}, '@');
		if (~isempty(ind)),		nomes{k}(ind(end)-1:end) = [];	end		% Remove the zooming info
	end
	set(handles.popup_familyPlots, 'Str', nomes)
	grid_figs(handles, N_figs)

	handles.hAllFigs = hAllFigs;	% Those are Fig handles

	sizes_cm = {'A0 (84.1 118.9 cm)'; 'A1 (59.4 84.1 cm)'; 'A2 (42.0 59.4 cm)'; 'A3 (29.70 42.0 cm)'
		'A4 (21.0 29.70 cm)'; 'A5 (14.8 21.0 cm)'; 'B0 (100.05 141.39 cm)'
		'B1 (70.70 100.05 cm)'; 'B2 (50.02 70.70 cm)'; 'B3 (35.35 50.02 cm)'; 'B4 (25.01 35.35 cm)'
		'B5 (17.67 25.01 cm)'; 'archA (22.86 30.48 cm)'; 'archB (30.48 45.72 cm)'; 'archC (45.72 60.96 cm)'
		'archD (60.96 91.44 cm)'; 'archE (91.44 121.92 cm)'; 'flsa (21.59 33.02 cm)'; 'halfletter (13.97 21.59 cm)'
		'note (19.05 25.40 cm)'; 'letter (21.59 27.94 cm)'; 'legal (21.59 35.56 cm)'; '11x17 (27.94 43.18 cm)'
		'ledger (43.18 27.94 cm)'};

	sizes_pt = {'A0  (2380 3368 pt)'; 'A1  (1684 2380 pt)'; 'A2  (1190 1684 pt)'; 'A3  (842 1190 pt)'
		'A4  (595 842 pt)'; 'A5  (421 595 pt)'; 'B0  (2836 4008 pt)'; 'B1  (2004 2836 pt)'; 'B2  (1418 2004 pt)'
		'B3  (1002 1418 pt)'; 'B4  (709 1002 pt)'; 'B5  (501 709 pt)'; 'archA (648 864 pt)'; 'archB (864 1296 pt)'
		'archC (1296 1728 pt)'; 'archD (1728 2592 pt)'; 'archE (2592 3456 pt)'; 'flsa  (612 936 pt)'
		'halfletter (396 612 pt)'; 'note   (540 720 pt)'; 'letter (612 792 pt)'; 'legal  (612 1008 pt)'
		'11x17  (792 1224 pt)'; 'ledger (1224 792 pt)'};

	sizes_in = {'A0 (33.06 46.78 cm)'; 'A1 (23.39 33.06 in)'; 'A2 (16.53 23.39 in)'; 'A3 (11.69 16.53 in)'
		'A4 (8.26 11.69 in)'; 'A5 (5.85 8.26 in)'; 'B0 (39.39 55.67 in)'
		'B1 (27.83 39.39 in)'; 'B2 (19.69 27.83 in)'; 'B3 (13.92 19.69 in)'; 'B4 (9.85 13.92 in)'
		'B5 (6.96 9.85 in)'; 'archA (9.0 12.0 in)'; 'archB (12.0 18.0 in)'; 'archC (18.0 24.0 in)'
		'archD (24.0 36.0 in)'; 'archE (36.0 48.0 in)'; 'flsa (8.5 13.0 in)'; 'halfletter (5.5 8.5 in)'
		'note (7.5 10.0 in)'; 'letter (8.5 11.0 in)'; 'legal (8.5 14.0 in)'; '11x17 (11.0 17.0 in)'
		'ledger (17.0 11.0 in)'};

	paper_cm = [84.1 118.9; 59.40 84.1; 42.0 59.4; 29.7 42.0; 21.0 29.7; 14.8 21.0; ...
		100.05 141.40; 70.7 100.05; ...
		50.02 70.70; 35.35 50.02; 25.01 35.35; 17.67 25.01; 22.86 30.48; 30.48 45.72; ...
		45.72 60.96; 60.96 91.44; 91.44 121.92; 21.59 33.02; 13.97 21.59; 19.05 25.40; ...
		21.59 27.94; 21.59 35.56; 27.94 43.18; 43.18 27.94];

	paper_pt = [2380 3368; 1684 2380; 1190 1684; 842 1190; 595 842; 421 595; 
		2836 4008; 2004 2836; 1418 2004; 1002 1418; 709 1002; 501 709; 648 864; 864 1296; ...
		1296 1728; 1728 2592; 2592 3456; 612 936; 396 612; 540 720; 612 792; 612 1008; 792 1224; 1224 792];

	paper_in = [33.06 46.78; 23.39 33.06; 16.53 23.39; 11.69 16.53; 8.26 11.69; 5.85 8.26;
		39.39 55.67; 27.83 39.39; 19.69 27.83; 13.92 19.69; ...
		9.85 13.92; 6.96 9.85; 9.0 12.0; 12.0 18.0; 18.0 24.0; 24.0 36.0; 36.0 48.0; 8.5 13.0; 5.5 8.5; ...
		7.5 10.0; 8.5 11.0; 8.5 14.0; 11.0 17.0; 17.0 11.0];

	set(handles.popup_paperSize,'String',sizes_cm,'Value',5)

	paper = [paper_cm(5,1) paper_cm(5,2)];			% Set to A4 (x,y)
	set(handles.axes1,'XLim',[0 paper(1)],'YLim',[0 paper(2)], 'DataAspectRatio',[1 1 1]);

	% Set the default projections ofered here
	handles.projGDAL_name = {
			'Linear (meaning, no projection)'; ...
			'Geographical                                   (EPSG:4326)'; ...
			'Equidistant Cilindrical (Plate Caree Sphere)   (EPSG:3786)'; ...
			'World Mercator                                 (EPSG:54004)'; ...
			'Pseudo Mercator                                (EPSG:3857)'; ...
			'ETRS89 / Portugal TM06                         (EPSG:3763)'; ...
			'Datum 73 (Portugal)                            (EPSG:27493)'; ...
			'WGS84 / UTM zone 29N                           (EPSG:32629)'; ...
			'Miller Cylindrical (Sphere)                    (EPSG:53003)'; ...
			'Gall Stereographic (Sphere)                    (EPSG:53016)'; ...
			'Cassini (Sphere)                               (EPSG:53028)'; ...
			'Sinusoidal (Sphere)                            (EPSG:53008)'; ...
			'Mollweide (Sphere)                             (EPSG:53009)'; ...
			'Robinson (Sphere)                              (EPSG:53030)'; ...
			'Eckert IV'; ...
			'Eckert VI'; ...
			'Lambert Conformal Conic'; ...
			'Equidistant Conic'; ...
			'Albers Equal Area'; ...
			'Lambert Equal Area';  ...
			'Bonne Sphere                                   (EPSG:53024)'; ...
			'North Pole Stereographic                       (EPSG:102018)'; ...
			'Van der Grinten I (Sphere)                     (EPSG:53029)'};
	
	handles.projGDAL_pars = {
			''; ...
			'+proj=longlat'; ...
			'+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs'; ...
			'+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'; ...
			'+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs'; ...
			'+proj=tmerc +lat_0=39.66825833333333 +lon_0=-8.133108333333334 +k=1 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'; ...
			'+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111112 +k=1 +x_0=180.598 +y_0=-86.99 +ellps=intl +towgs84=-223.237,110.193,36.649,0,0,0,0 +units=m +no_defs'; ...
			'+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs'; ...
			'+proj=mill +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +R_A +a=6371000 +b=6371000 +units=m +no_defs'; ...
			'+proj=gall +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs'; ...
			'+proj=cass +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs'; ...
			'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs'; ...
			'+proj=moll +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs'; ...
			'+proj=robin +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs'; ...
			'+proj=eck4'; ...
			'+proj=eck6'; ...
			'+proj=lcc +lat_1=20n +lat_2=60n'; ...
			'+proj=eqdc +lat_1=15n +lat_2=75n'; ...
			'+proj=aea +lat_1=20n +lat_2=60n'; ...
			'+proj=laea +lat_1=20n +lat_2=60n'; ...
			'+proj=bonne +lon_0=0 +lat_1=60 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs'; ...
			'+proj=stere +lat_0=90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'; ...
			'+proj=vandg +lon_0=0 +x_0=0 +y_0=0 +R_A +a=6371000 +b=6371000 +units=m +no_defs';};

	set(handles.popup_projections,'String',handles.projGDAL_name)

	handles.sizes_cm = sizes_cm;    handles.paper_cm = paper_cm;
	handles.sizes_pt = sizes_pt;    handles.paper_pt = paper_pt;
	handles.sizes_in = sizes_in;    handles.paper_in = paper_in;
	handles.opt_L = [];
	handles.psc_res = [];			handles.psc_opt_W = [];		% pscoast stuff
	handles.psc_type_p = [];		handles.psc_type_r  = [];	%		"
	handles.opt_psc = [];			% To eventualy hold several of the pscoast options
	handles.scale_set = false;		% To signal that user changed scale
	handles.supported_proj = false;	% Will turn to true when confirmed that we can make a map with GMT
	handles.d_path = handMir.path_data;
	handles.IamCompiled = handMir.IamCompiled;
	handles.path_tmp = handMir.path_tmp;
	handles.which_unit = 'cm';
	handles.last_dir = handMir.last_dir;
	handles.version7 = handMir.version7;
	handles.hRect = zeros(N_figs, 1) * NaN;
	if (strncmp(computer, 'PC', 2)),	handles.script_type = 'bat';
	else,								handles.script_type = 'bash';
	end
	
	% These are for Mirone calls not barf
	handles.no_file = false;
	handles.DefLineColor = 'k';
	handles.DefLineThick = 0.5;
	handles.head = [0 21 0 29.7 0 1 0 0.001 0.001];		% FAKE
	handles.geog = 0;
	handles.image_type = 20;
	handles.validGrid = false;
	handles.Tesoura = [];
	handles.Mao = [];
	setappdata(handles.axes1,'ThisImageLims', [0 21 0 29.7])

	% Loop over the number of figures and do preparatory work
	for (k = 1:N_figs)
		handles = get_img_dims(handles, k);		% In a function so it can later be called on different images
		handles = draw_img_rectangle(handles, k);
	
		% If the caller is 'ecran' we must set several fake handles struct members that exist in Mirone
		this_hMir = guidata(hAllFigs(k));
		if (strcmp(get(this_hMir.axes1, 'UserData'), 'XY'))
			this_hMir.image_type = 20;
			this_hMir.geog = 0;			this_hMir.head = [];
			this_hMir.PalAt = [];		this_hMir.PalIn = [];
			this_hMir.grdname = [];
			this_hMir.path_data = handles.d_path;		% Idiot name confusion
			this_hMir.IamXY = true;
		else
			this_hMir.IamXY = false;
		end
		guidata(this_hMir.figure1, this_hMir)
	end

	% ------------------ Set prefix name based on grid/image name --------------------------------------
	[lixo,figname] = fileparts(get(handMir.figure1, 'Name'));
	set(handles.edit_prefix,'String',strtok(figname))

	% ---------------------------------------------------------------------------------------------------

	% ----------- Pick up the projection initial (and sometimes final) guess ---------------------------
	if (handMir.is_projected)
		prj4 = aux_funs('get_proj_string', handMir.figure1);
		if (~isempty(prj4))
			set(handles.edit_projection,'String', prj4, 'Enable', 'off')
		else
			set(handles.edit_projection,'String', '', 'Enable', 'off')
			warndlg('This file claims to be projected but I can''t find its projection info.', 'Warning')
		end
		set(handles.popup_projections,'String', 'Already projected (can''t change)', 'Enable', 'off')
	elseif (handMir.geog)
		%  Here we should remember the last choice of how to deal with geogs
		set(handles.popup_projections,'Val', 2)
		set(handles.edit_projection,'String', handles.projGDAL_pars{2})
	else
		set(handles.popup_projections,'String', 'Linear (can''t change)', 'Enable', 'off')
		set(handles.edit_projection,'String', '', 'Enable', 'off')
	end
	handles.projection_str = cell(N_figs,1);		% SO FAR MAKE THEM INITIALLY ALL EQUAL
	for (k = 1:N_figs)
		handles.projection_str{k} = get(handles.edit_projection,'String');		% To restore in case of error
	end
	% --------------------------------------------------------------------------------------------------

	% ----------------- Use the directory list from mirone_pref ----------------------------------------
	directory_list = [];
	load([handles.d_path 'mirone_pref.mat']);

	j = false(1,numel(directory_list));					% vector for eventual cleaning non-existing dirs
	if iscell(directory_list)							% When exists a dir list in mirone_pref
		for (i = 1:length(directory_list))
			if (~exist(directory_list{i},'dir')),   j(i) = true;   end
		end
		directory_list(j) = [];							% clean eventual non-existing directories
		%directory_list    = [{handMir.last_dir}; directory_list];
		set(handles.popup_directory_list, 'String', directory_list)
		handles.last_directories = directory_list;
	else												% mirone_pref had no dir list
		handles.last_directories = {handles.last_dir};
		set(handles.popup_directory_list,'String',handles.last_directories)
	end
	% --------------------------------------------------------------------------------------------------

	%------------ Give a Pro look (3D) to the frame boxes  ---------------------------------------------
	new_frame3D(hObject, [handles.text_ps handles.text_out handles.text_dim])
	hFrame = findobj(handles.figure1,'Style','Frame', 'Tag','L');
	% But now we also need to update the UserData of the Right side frames so the resizing works correctly
	for (k = 1:numel(hFrame))
		ud = get(hFrame(k), 'UserData');	set(hFrame(k), 'UserData', ud+1)
	end
	posFig = get(handles.figure1, 'Pos');
	hFrame = findobj(handles.figure1,'Style','Frame', 'Tag','R');
	for (k = 1:numel(hFrame))
		pos = get(hFrame(k), 'Pos');		set(hFrame(k), 'UserData', posFig(3)-pos(1))
	end
	hFrame = findobj(handles.figure1,'Style','Frame', 'Tag','RR');
	for (k = 1:numel(hFrame))
		pos = get(hFrame(k), 'Pos');		set(hFrame(k), 'UserData', posFig(3)-pos(1))
	end
	% --------------------------------------------------------------------------------------------------

	% ------------ Apply inherited projection ----------------------------------------------------------
	guidata(hObject, handles);
	for (k = 1:N_figs)
		update_scales(handles, k)
		handles = guidata(hObject);		% Recover in "this handles" the changes donne in update_scales
	end

	% ------------- Add this figure handle to the carraças list ----------------------------------------
	plugedWin = getappdata(handMir.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handMir.figure1,'dependentFigs',plugedWin);
	% --------------------------------------------------------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% -----------------------------------------------------------------------------------------
function handles = get_img_dims(handles, N, width)
% ...
	if (nargin == 2),	width  = 15;	end			% Default width in cm
	handMir = guidata(handles.hAllFigs(N));
	imgXlim = get(handMir.axes1,'XLim');    imgYlim = get(handMir.axes1,'YLim');
	if (handMir.image_type == 2 || handMir.image_type == 20)		% "trivial" images
		if (~handMir.IamXY)
			[ny,nx,nz] = size(get(handMir.hImg,'CData'));
		else
			nx = 300;	ny = 200;		% Pure INVENTION
		end
		handles.x_min(N) = imgXlim(1);     handles.x_max(N) = imgXlim(2);
		handles.y_min(N) = imgYlim(1);     handles.y_max(N) = imgYlim(2);
	else
		head = handMir.head;
		if (handMir.was_pixreg)			% If we under the hood converted a pix reg to grid reg, undo it now
			dx2 = handMir.head(8) / 2;		dy2 = handMir.head(9) / 2;
			head(1) = head(1)-dx2;			head(2) = head(2)+dx2;
			head(3) = head(3)-dy2;			head(4) = head(4)+dy2;
		end
		if (diff(imgXlim) < diff(head(1:2))*0.9)		% Figure is zoomed
			handles.x_min(N) = imgXlim(1) + head(8)/2;		handles.x_max(N) = imgXlim(2) - head(8)/2;
			handles.y_min(N) = imgYlim(1) + head(9)/2;		handles.y_max(N) = imgYlim(2) - head(9)/2;
		else
			handles.x_min(N) = head(1);    handles.x_max(N) = head(2);
			handles.y_min(N) = head(3);    handles.y_max(N) = head(4);
		end
		nx = round((head(2) - head(1)) / head(8));		% May not be exactly correct but is good enough here
		ny = round((head(4) - head(3)) / head(9));
	end

	ratio  = nx / width;
	height = round(ny) / ratio;
	if (height > 29)				% That is, if height + Y0 nearly outside the A4 page
		while (height > 29)			% Make it approximately fit the page height
			height = round(height * 0.1);
			width  = round(width  * 0.1);
		end
	end

	handles.opt_R{N} = sprintf('-R%.12g/%.12g/%.12g/%.12g', ...
	                   handles.x_min(N), handles.x_max(N), handles.y_min(N), handles.y_max(N));
	handles.width_orig(N) = width;		handles.height_orig(N) = height;	handles.map_width{N} = width;

% -----------------------------------------------------------------------------------------
function handles = draw_img_rectangle(handles, N, X0, Y0)
% ...
	if (nargin == 2)
		X0 = 2.5;			Y0 = 2.0;
		if (N > 1)			% Cascade the figs from Bottom up
			Y0 = Y0 + (N - 1) * 4;
			X0 = X0 + (N - 1) * 0.5;
		end
	end
	rect_x = [X0 X0 X0+handles.width_orig(N) X0+handles.width_orig(N) X0];
	rect_y = [Y0 Y0+handles.height_orig(N) Y0+handles.height_orig(N) Y0 Y0];
	handles.rect_x{N} = rect_x;   handles.rect_y{N} = rect_y;
	handles.map_width{N} = sprintf('%.2g', handles.map_width{N});	% Now it's a string
	set(handles.edit_mapWidth, 'String', handles.map_width{N})		% Fill the width editbox
	set(handles.edit_mapHeight,'String', sprintf('%.2f',handles.height_orig(N)))	% Fill the height editbox

	% ---------- Draw the image rectangle
	if (~ishandle(handles.hRect(N)))			% First time, create it
		handles = create_images(handles, rect_x, rect_y, N);
		h = patch('Parent',handles.axes1,'XData',rect_x,'YData',rect_y,'FaceColor','none', ...
				  'EdgeColor','k', 'LineWidth',0.1, 'Tag','PlotRect');
		if (handles.version7 >= 8.4),	set(h, 'FaceColor', 'w', 'FaceAlpha', 0.005),	end		% F TMW, always breaking things
		cmenuHand = uicontextmenu('Parent',handles.figure1);
		set(h, 'UIContextMenu', cmenuHand);
		uimenu(cmenuHand, 'Label', 'Delete', 'Callback', {@delete_fig, h});
		uimenu(cmenuHand, 'Label', 'Frame settings', 'Callback', {@set_opt_B, h}, 'Sep', 'on');
		set(h, 'ButtonDownFcn', {@popup_familyPlots_CB, h, [], N})		% Must come before the ui_edit_polygon() call
		ui_edit_polygon(h, '')					% Set edition functions
		setappdata(h, 'RunCB', {@update_scales, h})
		setappdata(h, 'opt_B', '-Ba -BWSen')
		setappdata(h, 'My_N', N)				% 
		handles.hRect(N) = h;					% Save the rectangle hand
	else
		set(handles.hRect(N), 'XData',rect_x,'YData',rect_y)
	end

	% -----------------
	function handles = create_images(handles, rect_x, rect_y, N, opt)
		% This is NOT an nested function.
		% Whith 4 argins creates a new image. With 5, updates the size of handles.hAllImgs(N)
		% One of RECT_X or RECT_Y may be empty, case in which only one the axes is updated.
		if (nargin == 4)
			handMir = guidata(handles.hAllFigs(N));
			img = get(handMir.hImg, 'CData');
			n_rows = size(img,1);	n_cols = size(img,2);
			rows_tille = 128;		% Height of the thumb nail 
			fac = rows_tille / n_rows;
			cols_tille = round(n_cols * fac);
			img = cvlib_mex('resize', img, [rows_tille cols_tille], 'nearest');
			if (size(img,3) == 1)		% A indexed image, convert it to RGB to forget about the multi cmaps shit
				img = ind2rgb8(img, get(handMir.figure1,'Colormap'));
			end
			handles.orig_img{N} = img;	% Save a copy for (re)projections
			handles.hAllImgs(N) = image('XData',[rect_x(1) rect_x(3)],'YData',rect_y(1:2),'CData', img);
		else
			if (~isempty(rect_x) && ~isempty(rect_y))
				set(handles.hAllImgs(N), 'XData',[rect_x(1) rect_x(3)],'YData',rect_y(1:2))
			elseif (~isempty(rect_x))
				set(handles.hAllImgs(N), 'XData',[rect_x(1) rect_x(3)])
			elseif (~isempty(rect_y))
				set(handles.hAllImgs(N), 'YData',rect_y(1:2))
			end
		end

% -----------------------------------------------------------------------------------------
function popup_familyPlots_CB(hObject, handles, N)
% ...
	if (nargin == 2)		% Is == 3 when this fun is called by the image's rectangle buttondownfcn
		N = get(hObject, 'Val');
	else
		handles = guidata(hObject);			% In this case hObject is the image's patch rectangle handle
		N = getappdata(hObject, 'My_N');	% This is the true N, that is updated in case some fig is deleted
		set(handles.popup_familyPlots, 'Val', N)
	end
	set(handles.edit_projection,  'Str', handles.projection_str{N})
	set(handles.edit_mapWidth,    'Str', handles.width_orig(N))
	set(handles.edit_mapHeight,   'Str', handles.height_orig(N))
	set(handles.edit_X0,          'Str', handles.X0(N))
	set(handles.edit_Y0,          'Str', handles.Y0(N))
	set(handles.edit_scale,       'Str', handles.scale{N})

% -----------------------------------------------------------------------------------------
function popup_gridFigs_CB(hObject, handles)
% Arrange figs in a grid
	val = get(hObject, 'Value');	str = get(hObject, 'String');	str = str{val};
	N = numel(get(handles.popup_familyPlots, 'Str'));
	X0 = 2.5;	Y0 = 2.0;

	if (strcmp(str, 'cascade'))
		for (n = 1:N)
			y0 = Y0 + (n - 1) * 4;
			x0 = X0 + (n - 1) * 0.5;
			xx = get(handles.hRect(n), 'XData');	yy = get(handles.hRect(n), 'YData');
			xx = x0 + xx - xx(1);		yy = y0 + yy - yy(1);
			set(handles.hRect(n), 'XData',xx, 'YData',yy)
		end
		return
	end

	n_rows = str2double(str(1));	n_cols = str2double(str(3));
	pad = 1.5;	% 1 cm
	iy = 1;
	if (get(handles.radio_Landscape, 'Val')),	iy = 2;		end
	ind = get(handles.popup_paperSize, 'Val');
	f_width = round((handles.paper_cm(ind, iy) - (n_cols-1)*pad - X0) / n_cols);
	row_height = zeros(1, n_rows);
	k = 1;
	for (n = 1:n_rows)
		this_row_height = 0;
		for (m = 1:n_cols)
			if (k <= N)			% Because we may have more cells than figures
				handles = get_img_dims(handles, k, f_width);
				this_row_height = max(this_row_height, handles.height_orig(k));
			end
			k = k + 1;
		end
		row_height(n) = this_row_height;
	end

	% Now that we know the height of each row we can draw the images bounding boxes
	k = 1;
	for (n = 1:n_rows)
		if (n > 1),		Y0 = Y0 + sum(row_height(1:n-1)) + (n-1) * pad;	end
		X0 = 2.5;
		for (m = 1:n_cols)
			if (k <= N)
				X0 = X0 + (m-1) * (f_width + pad);
				handles = draw_img_rectangle(handles, k, X0, Y0);
				update_scales(handles, k)
				handles = guidata(hObject);		% Recover in "this handles" the changes donne in update_scales
			end
			k = k + 1;
		end
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function grid_figs(handles, N_figs)
	if (N_figs > 1)			% Set up the figures grid options
		if     (N_figs == 2),	str = {'2x1'; '1x2'};
		elseif (N_figs == 3),	str = {'3x1'; '1x3'};
		elseif (N_figs == 4),	str = {'2x2'; '4x1'; '1x4'};
		elseif (N_figs == 5),	str = {'2x2'; '3x2'; '2x3'};
		elseif (N_figs == 6),	str = {'3x2'; '2x3'};
		elseif (N_figs == 7),	str = {'4x2'; '2x4'};
		elseif (N_figs == 8),	str = {'4x2'; '2x4'};
		elseif (N_figs == 9),	str = {'3x3'; '5x2'; '2x5'};
		else,					str = {'5x2'; '2x5'};
		end
		str{end+1} = 'cascade';
		set(handles.popup_gridFigs, 'Str', str, 'Vis', 'on')
	else
		set(handles.popup_gridFigs, 'Vis', 'off')
	end

% -----------------------------------------------------------------------------------------
function popup_paperSize_CB(hObject, handles)
	val = get(hObject,'Value');		str = get(hObject,'String');	str = str{val};
	switch str
		case 'cm',			lims = handles.paper_cm(val,1:2);
		case 'in',			lims = handles.paper_in(val,1:2);
		case 'pt',			lims = handles.paper_pt(val,1:2);
	end
	if (get(handles.radio_Portrait,'Value'))
		set(handles.axes1,'XLim',[0 lims(1)],'YLim',[0 lims(2)]);
	else
		set(handles.axes1,'XLim',[0 lims(2)],'YLim',[0 lims(1)]);
	end

% -----------------------------------------------------------------------------------------
function popup_paperUnits_CB(hObject, handles)
	val = get(hObject,'Value');		str = get(hObject, 'String');
	lav = get(handles.popup_paperSize,'Value');
	switch str{val}
		case 'cm'
			set(handles.popup_paperSize,'String',handles.sizes_cm,'Value',lav)
			if (get(handles.radio_Portrait,'Value'))
				set(handles.axes1,'XLim',[0 handles.paper_cm(lav,1)],'YLim',[0 handles.paper_cm(lav,2)])
			else
				set(handles.axes1,'XLim',[0 handles.paper_cm(lav,2)],'YLim',[0 handles.paper_cm(lav,1)])
			end
		case 'in'
			set(handles.popup_paperSize,'String',handles.sizes_in,'Value',lav)
			if (get(handles.radio_Portrait,'Value'))
				set(handles.axes1,'XLim',[0 handles.paper_in(lav,1)],'YLim',[0 handles.paper_in(lav,2)])
			else
				set(handles.axes1,'XLim',[0 handles.paper_in(lav,2)],'YLim',[0 handles.paper_in(lav,1)])
			end
		case 'pt'
			set(handles.popup_paperSize,'String',handles.sizes_pt,'Value',lav)
			if (get(handles.radio_Portrait,'Value'))
				set(handles.axes1,'XLim',[0 handles.paper_pt(lav,1)],'YLim',[0 handles.paper_pt(lav,2)])
			else
				set(handles.axes1,'XLim',[0 handles.paper_pt(lav,2)],'YLim',[0 handles.paper_pt(lav,1)])
			end
	end
	conv_units(handles, str{val})
	handles.which_unit = str{val};
	guidata(hObject,handles);

% -----------------------------------------------------------------------------------------
function conv_units(handles, dest, N)
	if (nargin == 2),	N = 1;	end
	if (strcmp(handles.which_unit,'cm') && strcmp(dest,'in')),		fact = 1  / 2.54;
	elseif (strcmp(handles.which_unit,'cm') && strcmp(dest,'pt')),	fact = 72 / 2.54;
	elseif (strcmp(handles.which_unit,'in') && strcmp(dest,'cm')),	fact = 2.54;
	elseif (strcmp(handles.which_unit,'in') && strcmp(dest,'pt')),	fact = 72;
	elseif (strcmp(handles.which_unit,'pt') && strcmp(dest,'cm')),	fact = 2.54/72;
	elseif (strcmp(handles.which_unit,'pt') && strcmp(dest,'in')),	fact = 1/72;
	end

	set(handles.edit_mapWidth,'String',sprintf('%.2f',str2double(get(handles.edit_mapWidth,'String')) * fact))
	set(handles.edit_mapHeight,'String',sprintf('%.2f',str2double(get(handles.edit_mapHeight,'String')) * fact))
	set(handles.edit_X0,'String',sprintf('%.2f',str2double(get(handles.edit_X0,'String')) * fact))
	set(handles.edit_Y0,'String',sprintf('%.2f',str2double(get(handles.edit_Y0,'String')) * fact))

	xx = get(handles.hRect(N),'XData') * fact;			yy = get(handles.hRect(N),'YData') * fact;
	set(handles.hRect(N),'XData',xx, 'YData',yy);

% -----------------------------------------------------------------------------------------
function radio_Portrait_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_Landscape,'Value',0)

	posR   = get(handles.popup_familyPlots, 'Pos');
	posFig = get(handles.figure1, 'Pos');
	posAx  = get(handles.axes1, 'Pos');
	H = posFig(4) -100 -4;
	posAx(3) = H / sqrt(2);		posAx(4) = H;
	x1 = 25 + (posR(1)-30)/2 - posAx(3)/2;		% 25 is the left margin for the labels and 5 is for right border
	new_pos = [x1 100 posAx(3:4)];
	set(handles.axes1, 'XLim',get(handles.axes1,'YLim'), 'YLim',get(handles.axes1,'XLim'),'Pos', new_pos);

% -----------------------------------------------------------------------------------------
function radio_Landscape_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_Portrait,'Value',0)

	posR   = get(handles.popup_familyPlots, 'Pos');
	posFig = get(handles.figure1, 'Pos');
	posAx  = get(handles.axes1, 'Pos');
	maxW = (posFig(4) -100 -4) * sqrt(2);
	W = min(posR(1) - 30, maxW);	% 25 is the left margin for the labels and 5 is for right border
	%W = posR(1) - 30;			% 25 is the left margin for the labels and 5 is for right border
	posAx(3) = W;		posAx(4) = W / sqrt(2);
	x1 = 25;
	y1 = 100 + (posFig(4) - 100 - 4)/2 - posAx(4)/2;
	new_pos = [x1 y1 posAx(3:4)];
	set(handles.axes1, 'XLim',get(handles.axes1,'YLim'), 'YLim',get(handles.axes1,'XLim'),'Pos', new_pos);

% -----------------------------------------------------------------------------------------
function radio_autoRun_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_writeScript, 'Val', 0)

% -----------------------------------------------------------------------------------------
function radio_writeScript_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_autoRun, 'Val', 0)

% -----------------------------------------------------------------------------------------
function popup_directory_list_CB(hObject, handles)
	val = get(hObject,'Value');		str = get(hObject, 'String');
	% Put the selected field on top of the String list.
	tmp = str(val);			str(val) = [];
	new_str = [tmp; str];	set(hObject,'String',new_str,'Value',1); 

% -----------------------------------------------------------------------------------------
function push_change_dir_CB(hObject, handles)
	pato = handles.last_dir;
	if (strncmp(computer, 'PC', 2))
		work_dir = uigetfolder_win32('Select scripts folder',pato);
	else			% This guy doesn't let to be compiled
		work_dir = uigetdir(pato, 'Select scripts folder');
	end
	if (~isempty(work_dir))
		% Because unique returns a sorted array and I want the just selected dir on top, do this trick
		directory_list = unique([cellstr(work_dir); handles.last_directories]);
		directory_list = [cellstr(work_dir); directory_list];	% Put it on top
		for (k = 2:numel(directory_list))
			if (strcmp(directory_list{k}, work_dir))			% Find the repetition
				break
			end
		end
		directory_list(k) = [];		% This is the repeated entry
		handles.last_directories = directory_list;
		set(handles.popup_directory_list,'String', directory_list)
		guidata(hObject, handles);

		version7 = version;
		V7 = (sscanf(version7(1),'%f') > 6);
		if (~V7)                  % R <= 13
			save([handles.d_path 'mirone_pref.mat'],'directory_list', '-append')
		else
			save([handles.d_path 'mirone_pref.mat'],'directory_list', '-append', '-v6')
		end
	end

% -----------------------------------------------------------------------------------------
function edit_scale_CB(hObject, handles, N)
% ...
	if (nargin == 2),	N = get(handles.popup_familyPlots, 'Value');	end
	handMir = guidata(handles.hAllFigs(N));
	str = get(hObject,'String');
	xx  = get(handles.hRect(N),'XData');	yy = get(handles.hRect(N),'YData');
	if (~handMir.is_projected)
		opt_J = create_opt_J(handles, 1, str);
		in = [handles.x_min(N) handles.y_min(N); handles.x_min(N) handles.y_max(N); ...
		      handles.x_max(N) handles.y_max(N); handles.x_max(N) handles.y_min(N)];
		try
			opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g', handles.x_min(N), handles.x_max(N), handles.y_min(N), handles.y_max(N));
			out = c_mapproject(in,opt_R,opt_J,['-D' handles.which_unit(1)]);
			if (isa(out, 'struct')),	out = out.data;		end		% When GMT5, out is a struct
		catch
			return
		end
		xmax = max(out(:,1));		ymax = max(out(:,2));
		xx(3) = xmax+xx(1);			xx(4) = xmax+xx(1);
		yy(2) = ymax+yy(1);			yy(3) = ymax+yy(1);
	else
		scale = 1 / sscanf(str(3:end),'%f',1);
		dx_prj = handles.x_max(N) - handles.x_min(N);		% It's in projected meters
		dy_prj = handles.y_max(N) - handles.y_min(N);		% 		""
		width  = dx_prj * scale * 100;
		height = dy_prj * scale * 100;
		if (strcmp(handles.which_unit,'in')),     width = width * 2.54;     height = height * 2.54;     end
		if (strcmp(handles.which_unit,'pt')),     width = width * 2.54/72;  height = height * 2.54/72;  end
		xx(3:4) = xx(1) + width;
		yy(2:3) = yy(1) + height;
	end
	set(handles.hRect(N),'XData',xx,'YData',yy)
	handles.scale_set = true;
	guidata(hObject, handles);
	update_scales(handles)

% -----------------------------------------------------------------------------------------
function edit_mapWidth_CB(hObject, handles)
% Set new map width
	N = get(handles.popup_familyPlots, 'Value');
	str = get(hObject,'String');        w = str2double(str);
	if (isnan(w)),     set(hObject,'String',str),		return,		end
	xx = get(handles.hRect(N),'XData');
	xx(3) = xx(2) + w;      xx(4) = xx(1) + w;
	set(handles.hRect(N),'XData',xx)
	update_scales(handles, N)

% -----------------------------------------------------------------------------------------
function edit_mapHeight_CB(hObject, handles)
% Set new map height
	N = get(handles.popup_familyPlots, 'Value');
	str = get(hObject,'String');		h = str2double(str);
	if (isnan(h)),			set(hObject,'String',str),		return,		end
	yy = get(handles.hRect,'YData');
	yy(2) = yy(1) + h;		yy(3) = yy(4) + h;
	set(handles.hRect(N),'YData',yy)
	update_scales(handles, N)

% -----------------------------------------------------------------------------------------
function opt_J = create_opt_J(handles, N, scale)
% Create an option -J strin given the projection selected.
% This can be a classic -J or a new -J<proj4>

	if (isempty(N)),	N = 1;	end
	handMir = guidata(handles.hAllFigs(N));

	prj = handles.projection_str{N};
	if (~isempty(strfind(prj, '+proj=longlat')) || ~isempty(strfind(prj, '+proj=latlong')))
		prj = '-JX';		% +proj=latlong is still giving troubles in GMT
	end

	if (isempty(prj) || strncmp(prj, '-JX', 3))		% A linear image
		if (handMir.IamXY)					% An Ecran figure
			opt_J = sprintf('-JX%s%s/%s%s', get(handles.edit_mapWidth,'Str'), handles.which_unit(1), ...
			                 get(handles.edit_mapHeight,'Str'), handles.which_unit(1));
		else
			if (nargin == 3)		% This is the case used in edit_scale_CB() that needs a true 1:xxxx scale
				opt_J = ['-Jx' scale];
			else
				opt_J = ['-JX' handles.map_width{N} handles.which_unit(1) '/0'];
			end
		end
	else
		if (prj(1) == '+' || strncmpi(prj, 'epsg', 4))
			opt_J = sprintf('-J"%s +width=%s%c"', prj, handles.map_width{N}, handles.which_unit(1));
		else			% Presumably a GMT projection, otherwise error.
			opt_J = [prj '/' handles.map_width{N} handles.which_unit(1)];
			opt_J(3) = lower(opt_J(3));
		end
	end

% ----------------------------------------------------------------------------------------- ***
function update_scales(handles, N)
% Update the scale, sizes or the scale when those were changed.
% Plot a projected mini frame when it's not rectangular.
% This fun is also called as a callback registered in ui_edit_polygon()

	isfrom_RunCB = false;
	if (~isa(handles, 'struct'))		% Than this is a call from the RunCB registered by ui_edit_polygon
		N = getappdata(handles, 'My_N');%
		handles = guidata(handles);		% and handles is actually the rectangle handle.
		isfrom_RunCB = true;
	end
	handMir = guidata(handles.hAllFigs(N));

	xx = get(handles.hRect(N),'XData');		yy = get(handles.hRect(N),'YData');
	handles.X0(N) = xx(1);		handles.Y0(N) = yy(1);		% Save this to use in the push_OK_CB() loop
	set(handles.edit_X0,'String', sprintf('%.2f', handles.X0(N)));
	set(handles.edit_Y0,'String', sprintf('%.2f', handles.Y0(N)));

	new_w = sprintf('%.2f', (xx(3) - xx(2)));
	set(handles.edit_mapWidth,'String',new_w);
	handles.map_width{N} = new_w;

	opt_J = create_opt_J(handles, N);

	if (strncmp(opt_J, '-JX', 3))					% Linear proj has a different treatment
		handles.scale{N} = '1:1';					% Default value when no projected grids
		if (~handMir.IamXY)							% Only rescale if image, not XY plot
			new_x = xx(3) - xx(2);
			new_y = yy(2) - yy(1);
			if (isfrom_RunCB && ~isempty(strfind(handles.projection_str{N}, '+proj')))	% If Geog, maintain the 1:1 aspect ratio
				fac = new_x / (handMir.head(2) - handMir.head(1));		% cm/deg
				new_y = (handMir.head(4) - handMir.head(3)) * fac;
			end
			yy(2) = new_y + yy(1);      yy(3) = yy(2);			% It will become "True" scale
			set(handles.hRect(N),'YData', yy)
			create_images(handles, xx, yy, N, []);				% Update this axes size
			if (isfrom_RunCB)
				ui_edit_polygon(handles.hRect(N), '')			% Get it out of edit mode
			end
			handles.width_orig(N) = new_x;		handles.height_orig(N) = new_y;
			handles.scale{N} = get_scale(handles, new_x, new_y, N);
		end
		set(handles.edit_mapWidth,'String', handles.map_width{N});				% Uppdate map width
		set(handles.edit_mapHeight,'String',sprintf('%.2f', (yy(2) - yy(1))));	% Uppdate map height
		set(handles.edit_scale,'String', handles.scale{N})
		handles.supported_proj = true;

		guidata(handles.figure1, handles)
		return		% We are done
	end

	in = [handles.x_min(N) handles.y_min(N); handles.x_min(N) handles.y_max(N); ...
	      handles.x_max(N) handles.y_max(N); handles.x_max(N) handles.y_min(N)];
	try
		opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g', handles.x_min(N), handles.x_max(N), handles.y_min(N), handles.y_max(N));
		out = c_mapproject(in,opt_R,opt_J,['-D' handles.which_unit(1)]);
		if (isa(out, 'struct')),	out = out.data;		end		% When GMT5, out is a struct
		if (all(out(:) == 0))
			warndlg('Sorry, it is not yet possible to create maps in GMT with this PROJ.4 projection.','Warning')
			handles.supported_proj = false;
		else
			handles.supported_proj = true;
		end
	catch
		handles.supported_proj = false;
	end

	guidata(handles.figure1, handles)
	if (~handles.supported_proj),		return,		end

	new_y = max(out(:,2)) - min(out(:,2));

	yy(2) = new_y + yy(1);      yy(3) = yy(2);				% It will become "True" scale
	set(handles.hRect(N), 'XData', xx, 'YData', yy);
	create_images(handles, xx, yy, N, []);					% Update this image size
	if (isfrom_RunCB)
		ui_edit_polygon(handles.hRect(N), '')				% Get it out of edit mode
	end	
	set(handles.edit_mapWidth, 'String', handles.map_width{N});				% Uppdate map width
	set(handles.edit_mapHeight,'String', sprintf('%.2f', yy(2) - yy(1)));	% Uppdate map height
	handles.width_orig(N) = xx(3) - xx(2);		handles.height_orig(N) = yy(2) - yy(1);

	% ---------- Project the stamp image ---------------------------
	proj_img(handles, handMir, opt_J, N)

	% ----------- Compute scale 1:xxxx -----------------------------
	if (~handles.scale_set)						% If user changed scale, don't compute it here
		xm = (handles.x_min(N) + handles.x_max(N)) / 2;   ym = (handles.y_min(N) + handles.y_max(N)) / 2;
		opt_R = sprintf('-R%f/%f/%f/%f', xm-2, xm+2, ym-2, ym+2);
		in = [handles.x_min(N) handles.y_min(N); handles.x_min(N) handles.y_max(N); ...
		      handles.x_max(N) handles.y_max(N); handles.x_max(N) handles.y_min(N)];

		prj = get(handles.edit_projection,'String');
		if (isempty(prj))
			opt_J = '-Jx1:1';
		else
			if (prj(1) == '+' || strncmpi(prj, 'epsg', 4))
				opt_J = sprintf('-J"%s"/1:1', prj);
			else		% Presumably a GMT projection, otherwise error.
				opt_J = [prj '/1:1'];
				opt_J(3) = lower(opt_J(3));
			end
		end
		out = c_mapproject(in,opt_R,'-C','-F', opt_J);
		if (isa(out, 'struct')),	out = out.data;		end		% When GMT5, out is a struct
		dx_prj = out(4,1) - out(1,1);			% It's in projected meters
		dy_prj = out(2,2) - out(1,2);			% It's in projected meters
		dx_rect = xx(4) - xx(1);				% Is in "cm", "in" or "pt". So convert to "cm"
		dy_rect = yy(2) - yy(1);				% Is in "cm", "in" or "pt". So convert to "cm"
		if (strcmp(handles.which_unit,'in')),     dx_rect = dx_rect * 2.54;     dy_rect = dy_rect * 2.54;     end
		if (strcmp(handles.which_unit,'pt')),     dx_rect = dx_rect * 2.54/72;  dy_rect = dy_rect * 2.54/72;  end
		scale = max(dx_rect/dx_prj/100, dy_rect/dy_prj/100);
		[n,d] = rat(scale,1e-9);
		if (n > 1),		d = d / n;		end
		handles.scale{N} = sprintf('1:%d', d);
		set(handles.edit_scale, 'String', handles.scale{N})
		handles.scale_set = false;
	else
		handles.scale{N} = get(handles.edit_scale, 'String');	% If custom scale, store it too
	end

	guidata(handles.figure1, handles)

	% -------------------------------------------------
	function proj_img(handles, handMir, opt_J, N)
	% ...
		hdrStruct.ULx = handMir.head(1);
		hdrStruct.ULy = handMir.head(4);
		hdrStruct.Xinc = (handMir.head(2) - handMir.head(1)) / (size(handles.orig_img{N},2) - 1);
		hdrStruct.Yinc = (handMir.head(4) - handMir.head(3)) / (size(handles.orig_img{N},1) - 1);
		hdrStruct.DstProjSRS = opt_J(4:end-1);
		hdrStruct.nodata = 255;
		[img, att] = gdalwarp_mex(handles.orig_img{N}, hdrStruct);
		set(handles.hAllImgs(N), 'CData', img)

% -----------------------------------------------------------------------------------------
function scale_str = get_scale(handles, width, height, N)
% For now this only computes the scale when we have a projected grid
	scale_str = '1:1';
	handMir = guidata(handles.hAllFigs(N));
	if (handMir.is_projected)
		dx_prj = handles.x_max(N) - handles.x_min(N);		% It's in projected meters
		dy_prj = handles.y_max(N) - handles.y_min(N);		% It's in projected meters
		scale = max(width/dx_prj/100, height/dy_prj/100);
		[n,d] = rat(scale,1e-9);
		if (n > 1),    d = d / n;      end
		scale_str = sprintf('1:%d', d);
	end

% -----------------------------------------------------------------------------------------
function edit_X0_CB(hObject, handles)
% Set new x origin
	N = get(handles.popup_familyPlots, 'Value');
	str = get(hObject,'String');			x0 = str2double(str);
	if (isnan(x0)),		set(hObject,'String',str),		return,		end
	xx = get(handles.hRect(N),'XData');		xx = xx - xx(1) + x0;
	set(handles.hRect(N),'XData', xx)
	create_images(handles, xx, [], N, []);
	handles.X0(N) = x0;
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function edit_Y0_CB(hObject, handles)
% Set new y origin
	N = get(handles.popup_familyPlots, 'Value');
	str = get(hObject,'String');			y0 = str2double(str);
	if (isnan(y0)),     set(hObject,'String',str),		return,		end
	yy = get(handles.hRect(N),'YData');		yy = yy - yy(1) + y0;
	set(handles.hRect(N),'YData', yy)
	create_images(handles, [], yy, N, []);
	handles.Y0(N) = y0;
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function popup_projections_CB(hObject, handles)
	handMir = guidata(handles.hAllFigs(1));
	if (~handMir.geog)
		warndlg('Only GEOGRAPHIC coordinates can be projected. I don''t do reprojections yet.','Warning')
		return
	end
	prj = handles.projGDAL_pars{get(hObject,'Value')};
	set(handles.edit_projection,'String',prj)
	edit_projection_CB(handles.edit_projection, handles)

% -----------------------------------------------------------------------------------------
function edit_projection_CB(hObject, handles)
	try
		N = get(handles.popup_familyPlots, 'Value');
		prj = get(hObject, 'String');
		handles.projection_str{N} = prj;	% And save a copy
		update_scales(handles, N)			% It saves handles as well	
	catch
		str = get(hObject, 'String');
		errordlg(sprintf('This transformation\n%s\nis not valid. Error:\n%s', str, lasterr), 'ERROR')
		set(hObject, 'String', handles.projection_str{N})
	end

% -----------------------------------------------------------------------------------------
function push_worldCRS_CB(hObject, handles)
	warndlg('Sorry, not yet.','Warning')

% -----------------------------------------------------------------------------------------
function [hLine, res, opt_W, type_p, type_r] = find_psc_stuff(hLine)
% See if we have any pscoast stuff

	haveCoasts = false;     havePolitical = false;  haveRivers = false;
	res = [];           opt_W = [];         type_p = [];        type_r = [];
	h_c = findobj(hLine,'Tag','CoastLineNetCDF');
	if (~isempty(h_c))
		if (length(h_c) > 1),   h_c = h_c(1);     end
		CoastRes    = get(h_c,'UserData');
		LineWidth_c = get(h_c,'LineWidth');
		LineColor_c = get(h_c,'Color');
		LineStyle_c = get(h_c,'LineStyle');
		haveCoasts  = true;
	end
	h_p = findobj(hLine,'Tag','PoliticalBoundaries');
	if (~isempty(h_p))
		if (length(h_p) > 1),   h_p = h_p(1);     end
		zz = get(h_p,'UserData');
		if (iscell(zz)),    zz = zz{1};     end
		PoliticalRes = zz(1);        PoliticalType = zz(2);
		LineWidth_p = get(h_p,'LineWidth');
		LineColor_p = get(h_p,'Color');
		LineStyle_p = get(h_p,'LineStyle');
		havePolitical = true;
	end
	h_r = findobj(hLine,'Tag','Rivers');
	if (~isempty(h_r))
		if (length(h_r) > 1),   h_r = h_r(1);     end
		zz = get(h_r,'UserData');
		if (iscell(zz)),        zz = zz{1};     end
		RiversRes = zz(1);          RiversType = zz(2);
		LineWidth_r = get(h_r,'LineWidth');
		LineColor_r = get(h_r,'Color');
		LineStyle_r = get(h_r,'LineStyle');
		haveRivers = true;
	end
	hLine = setxor(hLine, [h_c; h_p; h_r]);

	if (haveCoasts || havePolitical || haveRivers)
		res_c = '';     res_p = '';     res_r = '';
		if (haveCoasts)
			cor = round(LineColor_c * 255);
			cor = [num2str(cor(1)) '/' num2str(cor(2)) '/' num2str(cor(3))];
			switch CoastRes
				case 'f',   res_c = ['-Df -W' num2str(LineWidth_c) 'p,' cor];
				case 'h',   res_c = ['-Dh -W' num2str(LineWidth_c) 'p,' cor];
				case 'i',   res_c = ['-Di -W' num2str(LineWidth_c) 'p,' cor];
				case 'l',   res_c = ['-Dl -W' num2str(LineWidth_c) 'p,' cor];
				case 'c',   res_c = ['-Dc -W' num2str(LineWidth_c) 'p,' cor];
			end
			if (~strcmp(LineStyle_c,'-'))   % If we have a line style other than solid
				switch LineStyle_c
					case '--',  res_c = [res_c 'ta'];
					case ':',   res_c = [res_c 'to'];
					case '-.',  res_c = [res_c 't10_2_2_5:5'];
				end
			end
		end
		if (havePolitical)
			switch PoliticalRes
				case 'f',   res_p = '-Df';
				case 'h',   res_p = '-Dh';
				case 'i',   res_p = '-Di';
				case 'l',   res_p = '-Dl';
				case 'c',   res_p = '-Dc';
			end
			cor = round(LineColor_p * 255);
			cor = [num2str(cor(1)) '/' num2str(cor(2)) '/' num2str(cor(3))];
			switch PoliticalType
				case '1',   type_p = ['-N1/'  num2str(LineWidth_p) 'p,' cor];
				case '2',   type_p = ['-N2/'  num2str(LineWidth_p) 'p,' cor];
				case '3',   type_p = ['-N3/'  num2str(LineWidth_p) 'p,' cor];
				case 'a',   type_p = ['-Na/'  num2str(LineWidth_p) 'p,' cor];
			end
			if (~strcmp(LineStyle_p,'-'))   % If we have a line style other than solid
				switch LineStyle_p
					case '--',  type_p = [type_p 'ta'];
					case ':',   type_p = [type_p 'to'];
					case '-.',  type_p = [type_p 't10_2_2_5:5'];
				end
			end
		end
		if (haveRivers)
			switch RiversRes
				case 'f',   res_r = '-Df';
				case 'h',   res_r = '-Dh';
				case 'i',   res_r = '-Di';
				case 'l',   res_r = '-Dl';
				case 'c',   res_r = '-Dc';
			end
			cor = round(LineColor_r * 255);
			cor = [num2str(cor(1)) '/' num2str(cor(2)) '/' num2str(cor(3))];
			switch RiversType
				case '1',   type_r = ['-I1/' num2str(LineWidth_r) 'p,' cor];
				case '2',   type_r = ['-I2/' num2str(LineWidth_r) 'p,' cor];
				case '3',   type_r = ['-I3/' num2str(LineWidth_r) 'p,' cor];
				case '4',   type_r = ['-I4/' num2str(LineWidth_r) 'p,' cor];
				case '5',   type_r = ['-I5/' num2str(LineWidth_r) 'p,' cor];
				case '6',   type_r = ['-I6/' num2str(LineWidth_r) 'p,' cor];
				case '7',   type_r = ['-I7/' num2str(LineWidth_r) 'p,' cor];
				case '8',   type_r = ['-I8/' num2str(LineWidth_r) 'p,' cor];
				case '9',   type_r = ['-I9/' num2str(LineWidth_r) 'p,' cor];
				case '10',  type_r = ['-I10/' num2str(LineWidth_r) 'p,' cor];
				case 'a',   type_r = ['-Ia/' num2str(LineWidth_r) 'p,' cor];
				case 'r',   type_r = ['-Ir/' num2str(LineWidth_r) 'p,' cor];
				case 'i',   type_r = ['-Ii/' num2str(LineWidth_r) 'p,' cor];
				case 'c',   type_r = ['-Ic/' num2str(LineWidth_r) 'p,' cor];
			end
			if (~strcmp(LineStyle_r,'-'))   % If we have a line style other than solid
				switch LineStyle_r
					case '--',  type_r = [type_r 'ta'];
					case ':',   type_r = [type_r 'to'];
					case '-.',  type_r = [type_r 't10_2_2_5:5'];
				end
			end
		end
		res = unique([res_c(1:3); res_p; res_r],'rows');  % We don't want repeated resolution strings
		if (size(res,1) > 1)        % Shit, we have mixed resolutions
			res = '-Di';            % TEMPORARY SOLUTION UNTIL I FIND HOW TO FIND THE HIGHEST COMMON RES
		end
		if (~isempty(res_c)),   opt_W = res_c(5:end);
		else,                   opt_W = [];
		end
	end

% ----------------------------------------------------------------------------------------- /////////
function push_OK_CB(hObject, handles)
% ...
	if (~handles.supported_proj)
		errordlg('I told you before that this is not possible (unsupported projection)', 'Error'),	return
	end

	if (~isempty(handles.hAllFigs))					% Otherwise we still can go. Maybe art drawings
		handMir = guidata(handles.hAllFigs(1));
	else
		handMir.version7 = handles.version7;		% We have to ahce at least these to not error bellow
		handMir.IamXY    = false;
		handMir.grdname  = '';
	end

	dest_dir = get(handles.popup_directory_list,'String');
	if (iscell(dest_dir)),		dest_dir = dest_dir{1};		end
	if (dest_dir(end) ~= '\' && dest_dir(end) ~= '/')
		dest_dir = [dest_dir filesep];
	end
	prefix = get(handles.edit_prefix,'String');

	writeScript = (get(handles.radio_writeScript, 'Val') ~= 0);
	do_MEX_fig  = (get(handles.radio_autoRun, 'Val') ~= 0);

	opt_P = '';
	if (get(handles.radio_Portrait,'Value')),	opt_P = ' -P';	end

	sc = handles.script_type;	ellips = 'WGS-84';
	hAlfaPatch = [];

	if (strcmp(sc,'bat')),		comm = 'REM ';		pb = '%';	pf = '%';
	else,						comm = '# ';		pb = '$';	pf = '';
	end
	if (strcmp(ellips,'WGS-84'))     % It is the default, so don't use any
		ellips = '';
	else
		ellips = [' --ELLIPSOID=' ellips];
	end

	opt_len_unit   = '';%' --PROJ_LENGTH_UNIT=point';
	opt_frameWidth = ' --MAP_FRAME_WIDTH=0.15c';
	frmPen = '';
	if (handMir.IamXY),		frmPen = '--MAP_FRAME_PEN=1.25p';	end

	o = 1;					% The MEX script counter. DO NOT USE IT FOR ANYTHING ELSE
	mex_sc = cell(15,2);

	prefix_ddir = [dest_dir prefix];
	pack = struct('comm',comm, 'pb',pb, 'pf',pf, 'do_MEX',do_MEX_fig, 'dest_dir',dest_dir, 'prefix',prefix, ...
		'prefix_ddir',prefix_ddir, 'ellips',ellips, 'opt_len_unit', opt_len_unit, 'RJOK',' -R -J -O -K', 'KORJ','');

	for (N = 1:numel(handles.hAllFigs))

		handMir = guidata(handles.hAllFigs(N));

		opt_annotsize = sprintf(' --FONT_ANNOT_PRIMARY=%dp', get(handMir.axes1, 'FontSize'));

		hLine = findobj(get(handMir.axes1,'Child'),'Type','line','Visible','on');
		% If we have costlines, need to use their (Mirone) settings 
		if (~isempty(findobj(hLine,'Tag','CoastLineNetCDF')) || ...
			~isempty(findobj(hLine,'Tag','Rivers')) || ...
			~isempty(findobj(hLine,'Tag','PoliticalBoundaries')))
			[hLine, handles.psc_res, handles.psc_opt_W, handles.psc_type_p, handles.psc_type_r] = ...
				find_psc_stuff(hLine);
		else
			handles.psc_res = '';
			handles.opt_psc = '';		% Need to reset for every time we don't find coasts in figure
		end
		if (~isempty(handles.psc_res))	% Means that we have coastlines and will use the Mirone settings
			handles.opt_psc = [handles.psc_res ' ' handles.psc_opt_W ' ' handles.psc_type_p ' ' handles.psc_type_r];
		end

		opt_J = create_opt_J(handles, N);
		if (handMir.geog && strncmp(opt_J, '-JX', 3) && opt_J(end) ~= 'd' && ~isempty(handles.psc_res))	% Append the 'd' otherwise pscoast barfs
			opt_J = [opt_J 'd']; %#ok<AGROW>
			ind = strfind(opt_J, '/');
			if (~isempty(ind))				% To turn (e.g) -JX15c/0d into -JX15cd/0d
				opt_J = [opt_J(1:ind(1)-1) 'd' opt_J(ind(1):end)];
			end
		elseif (strncmp(opt_J, '-JX', 3) && isempty(handles.projection_str{N}))		% Linear, has the freedom to deform
			opt_J = [opt_J(1:end-1) get(handles.edit_mapHeight,'String') handles.which_unit(1)];
		end
		pack.KORJ = [' -K -O ' handles.opt_R{N} ' ' opt_J];
		opt_B = getappdata(handles.hRect(N), 'opt_B');
		if (isempty(opt_B))
			opt_B = ' -Blbrt --MAP_FRAME_TYPE=plain --MAP_FRAME_PEN=faint';		% Trick to not plot any frame
		else	% If these were sneaked in via -B deactivate the default ones
			if (strfind(opt_B, '--FONT_ANNOT_PRIMARY')),	opt_annotsize = '';		end
			if (strfind(opt_B, '--MAP_FRAME_WIDTH')),		opt_frameWidth = '';	end
		end

		if (handMir.geog == 1),		opt_deg = '--FORMAT_GEO_MAP=ddd:mm:ss';
		elseif (handMir.geog == 2),	opt_deg = '--FORMAT_GEO_MAP=+ddd:mm:ss';
		else,						opt_deg = '';		% Non geog
		end

		opt_U = '';
		if (N == 1)
			[script, l, saveBind, id_grd, id_cpt] = do_init_script(handles, handMir, opt_B, handles.opt_R{1}, ...
				opt_J, opt_deg, opt_frameWidth, opt_annotsize, frmPen, pack);
			X0 = sprintf('-X%.3g%s', handles.X0(N), handles.which_unit(1));
			Y0 = sprintf('-Y%.3g%s', handles.Y0(N), handles.which_unit(1));
			X0_cum = handles.X0(1);		Y0_cum = handles.Y0(1);
			if (get(handles.check_timeStamp, 'Val')),	opt_U = ' -U';	end
		else			% The N > 1 cases are somewhat different			
			X0 = sprintf('-X%.3g%s', handles.X0(N)-handles.X0(1), handles.which_unit(1));
			Y0 = sprintf('-Y%.3g%s', handles.Y0(N)-handles.Y0(1), handles.which_unit(1));
			X0_cum = X0_cum + handles.X0(N)-handles.X0(1);
			Y0_cum = Y0_cum + handles.Y0(N)-handles.Y0(1);		% The accumulated shifts

			% UPDATE -R, -J and others for this Fig
			script{l} = sprintf('\nset lim=%s', handles.opt_R{N});		l = l + 1; %#ok<AGROW>
			script{l} = sprintf('\nset proj=%s', opt_J);				l = l + 1; %#ok<AGROW>
			script{l} = sprintf('set deg_form=%s', opt_deg);			l = l + 1; %#ok<AGROW>
		end

		% Start writing the GMT commands
		[script, mex_sc, l, o] = do_psbasemap(script, mex_sc, l, o, pack, handles.opt_R{N}, opt_J, opt_B, ...
			X0, Y0, opt_U, opt_P, opt_deg, opt_annotsize, frmPen, opt_frameWidth, N);

		[script, l, o, haveAlfa, used_grd, nameRGB] = do_grdimg(handMir, script, pack, l, o);
		hWait = [];
		if (pack.do_MEX && handMir.image_type ~= 20)
			[mex_sc, o, hWait] = do_grdimg_MEX(handMir, mex_sc, o, pack.KORJ);
		end

		% If we have used a grid, build the GMT palette
		[script, l, o, sc_cpt] = do_palette(handMir, script, l, o, pack, used_grd, id_cpt);

		% Coastlines section
		[script, mex_sc, l, o] = do_pscoast(handles, handMir, script, mex_sc, l, o, pack, handles.opt_R{1}, opt_J);

		% Search for contour lines
		hText = findobj(get(handMir.axes1,'Child'),'Type','text');
		[script, mex_sc, l, o, used_grd, hLine, hText] = do_contour(handMir, script, mex_sc, l, o, ...
																	pack, hLine, hText, used_grd);

		% Search for symbols
		[script, mex_sc, l, o, hLine] = do_symbols(handMir, script, mex_sc, l, o, pack, hLine, opt_J);

		hPatch = findobj(get(handMir.axes1,'Child'),'Type','patch');

		% Search for focal mecanisms
		[script, mex_sc, l, o, hLine, hPatch] = do_meca(handMir, script, mex_sc, l, o, pack, hLine, hPatch, opt_J);

		% Search for "telhas"
		[script, l, hPatch] = do_telhas(script, l, hPatch, pack);

		% Search for countries
		[script, l, haveAlfa, hPatch] = do_countries(handMir, script, l, pack, haveAlfa, hPatch);

		% Search for "Bar graphs"
		[script, mex_sc, l, o, hPatch, xx, yy] = do_bargraphs(handMir, script, mex_sc, l, o, pack, hPatch);

		% Search for "Histograms"
		[script, mex_sc, l, o, hPatch] = do_histograms(handMir, script, mex_sc, l, o, pack, hPatch, xx, yy);

		% Search for closed polygons
		[script, mex_sc, l, o] = do_patches(handMir, script, mex_sc, l, o, pack, hPatch, writeScript);

		% Search for lines associated with GMT custom symbols
		[script, mex_sc, l, o, hLine] = do_custom_symbols(handles, script, mex_sc, l, o, pack, hLine);

		% Search for lines or polylines
		[script, mex_sc, l, o] = do_lines(script, mex_sc, l, o, pack, hLine);

		% Search for text strings
		[script, mex_sc, l, o] = do_text(script, mex_sc, l, o, pack, hText);

		% Search for colorbar
		[script, mex_sc, l, o] = do_colorbar(handles, handMir, script, mex_sc, l, o, pack, sc_cpt);

		% See if we have to do a screen capture to 3 RGB grids
		do_screncapture(handMir, hAlfaPatch, haveAlfa, writeScript, nameRGB)

		% =============================== Search for "MagBars" (XY only) ===============================
		[script, mex_sc, l, o] = do_magbars(handMir, script, mex_sc, l, o, pack, saveBind, opt_J);
		% ==============================================================================================
	end
	
	hText  = findobj(get(handles.axes1,'Child'),'Type','text');
 	hPatch = findobj(get(handles.axes1,'Child'),'Type','patch');
	hLine  = findobj(get(handles.axes1,'Child'),'Type','line','Visible','on');
	hPatch = setxor(hPatch, handles.hRect);		% But setxor screws order by sorting the result

	want_ruler = false;
	%if (isempty(handles.hAllFigs)),	l = 1;	X0_cum = 0;		Y0_cum = 0;	opt_deg = '';	hWait = [];	end
	if (~isempty(hText) || ~isempty(hPatch) || ~isempty(hLine) || want_ruler)
		paper_size = handles.paper_cm(get(handles.popup_paperSize,'Val'), 1:2);
		ix = 1;		iy = 2;
		if (get(handles.radio_Landscape, 'Val')),		ix = 2;		iy = 1;		end
		opt_R = sprintf('-R0/%.3g/0/%.3g', paper_size(ix), paper_size(iy));
		script{l} = sprintf('\nset lim=%s', opt_R);			l = l + 1;
		script{l} = sprintf('set proj=-Jx1c');				l = l + 1;
		if (want_ruler),	ruler = '-Ba5f1WS --MAP_FRAME_TYPE=inside';
		else,				ruler = '-Blbrt';
		end
		if (~pack.do_MEX)
			if (want_ruler)
				script{l} = 'set frame=-Ba5f1WS --MAP_FRAME_TYPE=inside';	l = l + 1;
			else
				script{l} = 'set frame=';	l = l + 1;
			end
		end
		pack.KORJ = ' -K -O -R -Jx1c';
		X0 = sprintf('-X%.2g', -X0_cum);		Y0 = sprintf('-Y%.2g', -Y0_cum);
		[script, mex_sc, l, o] = do_psbasemap(script, mex_sc, l, o, pack, opt_R, '-Jx1c', ruler, ...
			X0, Y0, '', [opt_P ' -Vq'], opt_deg, opt_annotsize, frmPen, opt_frameWidth, 2);

		[script, mex_sc, l, o] = do_text(script, mex_sc, l, o, pack, hText);
		[script, mex_sc, l, o] = do_patches(handMir, script, mex_sc, l, o, pack, hPatch, writeScript);
		[script, mex_sc, l, o, hLine] = do_psimage(handles, script, mex_sc, l, o, pack, hLine);
		[script, mex_sc] = do_lines(script, mex_sc, l, o, pack, hLine);
	end

	if (pack.do_MEX)
		try		gsimage(handles, mex_sc, hWait)
		catch,	errordlg(lasterr, 'Error')
		end
		return
	end

	% ------------------------------------- WRITE THE SCRIPT ---------------------------------------
	% First do some eventual cleaning
	if (~isempty(handMir.grdname) && ~used_grd),	script(id_grd) = [];	end

	% Remove empties at the end of 'script' to not screw the last command patching below
	k = numel(script);
	while (isempty(script{k})),		k = k - 1;		end
	if (k < numel(script)),	script(k+1:end) = [];	end

	if (strcmp(sc,'bat'))
		fid = fopen([prefix_ddir '_mir.' sc],'wt');
	else
		fid = fopen([prefix_ddir '_mir.' sc],'wb');		% This way scripts are directly executable
	end
	for (i = 1:numel(script) - 1)
		fprintf(fid,'%s\n',script{i});
	end

	ind = strfind(script{i+1}, ' -K');
	if (~isempty(ind))
		script{i+1}(ind:ind+2) = [];		% Remove the last '-K'
	end
	fprintf(fid,'%s\n',script{i+1});
	fclose(fid);

	if (get(handles.radio_writeScript, 'Val'))
		msg{1} = ['File ' prefix '_mir.' handles.script_type ' successfuly created in:  ' dest_dir];
		h = msgbox(msg);		% This message self distructs in 4 sec
		pause(4)
		if (ishandle(h)),	delete(h),	end
	end

% ------------------------------------------------------------------------------------------------------------
function [script, l, saveBind, id_grd, id_cpt] = do_init_script(handles, handMir, opt_B, opt_R, opt_J, ...
	opt_deg, opt_frameWidth, opt_annotsize, frmPen, pack)
% Write the script initialization part (variables declaration)

	[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix] = unpack(pack);

	script = cell(40,1);
	id_grd = 0;			% Will force an error if used later and not changed below

	val   = get(handles.popup_paperSize,'Value');
	list  = get(handles.popup_paperSize,'String');
	str   = list{val};			k = strfind(str,' ');
	paper_media = str(1:k(1)-1);

	need_path = false;
	if (~isempty(handMir.grdname))
		[PATH,FNAME,EXT] = fileparts(handMir.grdname);
		just_grd_name = [FNAME EXT];
		ind = 0;
		if (dest_dir(end) == '\' || dest_dir(end) == '/'),	ind = 1;	end
		need_path = true;
		if (strcmp(PATH,dest_dir(1:end-ind))),	need_path = false;		end		% Because PATH has no trailing /
	end

	% ------------ Get size of Rectangle to use as info on the image size ------------------
	units = 'cm';
	val = get(handles.popup_paperUnits, 'Val');
	if (val == 2),		units = 'inch';
	elseif (val == 3),	units = 'points';
	end
	imgDimsInfo = sprintf(' ---- The image area has exactly %s x %s %s (unless you change -R or -J)', ...
		get(handles.edit_mapWidth,'Str'), get(handles.edit_mapHeight,'Str'), units);
	% --------------------------------------------------------------------------------------

	l = 1;
	if (~strcmp(handles.script_type, 'bat'))		% Write a bash script
		script{l} = '#!/bin/bash -f';				l=l+1;
		script{l} = [comm 'Coffeeright Mirone Tec'];l=l+1;
		script{l} = comm;							l=l+1;
		script{l} = [comm ' ---- Projection. You may change it if you know how to'];    l=l+1;
		script{l} = ['proj=' opt_J];		l=l+1;		% Map scale
		script{l} = [comm ' ---- Frame annotations. You may change it if you know how to'];    l=l+1;
		script{l} = ['frame=' opt_B];		l=l+1;      saveBind = l-1;
		script{l} = [comm imgDimsInfo];		l=l+1;
		script{l} = [comm ' ---- Map limits. You may change it if you know how to'];    l=l+1;
		script{l} = ['lim=' opt_R];			l=l+1;
		script{l} = comm;					l=l+1;
		script{l} = [comm ' ---- Longitude annotation style. Use the +ddd:mm:ss form => [0;360] range '];	l=l+1;
		script{l} = ['deg_form=' opt_deg];	l=l+1;
		script{l} = '';						l=l+1;
		script{l} = [comm ' ---- Width (> 0) of map borders for fancy map frame'];    l=l+1;
		script{l} = ['frame_width=' opt_frameWidth];	l=l+1;
		script{l} = [comm ' ---- Annotation font size in points'];		l=l+1;
		script{l} = ['annot_size=' opt_annotsize];      l=l+1;
 		if (~isempty(handMir.grdname))
			if (~need_path)
				script{l} = ['grd=' just_grd_name];		id_grd = l;		l=l+1;
			else
				script{l} = ['grd=' handMir.grdname];	id_grd = l;		l=l+1;
			end
		end
		script{l} = ['cpt=' prefix '.cpt'];		id_cpt = l;   l=l+1;
		script{l} = ['ps=' prefix '.ps'];			l=l+1;
		script{l} = ['gmtset PS_MEDIA=' paper_media]; l=l+1;
	else											% Write a dos batch    
		script{l} = '@echo OFF';					l=l+1;
		script{l} = [comm 'Coffeewrite Mirone Tec'];l=l+1;
		script{l} = comm;							l=l+1;
		script{l} = [comm ' ---- Projection. You may change it if you know how to'];		l=l+1;
		script{l} = ['set proj=' opt_J];			l=l+1;		% Map scale
		script{l} = [comm ' ---- Frame annotations. You may change it if you know how to'];	l=l+1;
		script{l} = ['set frame=' opt_B];			l=l+1;      saveBind = l-1;
		script{l} = [comm imgDimsInfo];				l=l+1;
		script{l} = [comm ' ---- Map limits. You may change it if you know how to'];		l=l+1;
		script{l} = ['set lim=' opt_R];				l=l+1;
		script{l} = comm;							l=l+1;
		script{l} = [comm ' ---- Longitude annotation style. Use the +ddd:mm:ss form => [0;360] range '];	l=l+1;
		script{l} = ['set deg_form=' opt_deg];		l=l+1;
		script{l} = '';								l=l+1;
		script{l} = [comm ' ---- Width (> 0) of map borders for fancy map frame'];    l=l+1;
		script{l} = ['set frame_width=' opt_frameWidth];      l=l+1;
		script{l} = [comm ' ---- Annotation font size in points'];    l=l+1;
		script{l} = ['set annot_size=' opt_annotsize];      l=l+1;
		if (handMir.IamXY)
			script{l} = sprintf('\n%s --- Map frame thickness in points.', comm);		l=l+1;
			script{l} = ['set framePen=' frmPen];	l=l+1;
		end
		script{l} = '';								l=l+1;
		if (~isempty(handMir.grdname))
			if (~need_path)
				script{l} = ['set grd=' just_grd_name];		id_grd = l; l=l+1;
			else
				script{l} = ['set grd=' handMir.grdname];	id_grd = l; l=l+1;
			end
		end
		script{l} = ['set cpt=' prefix '.cpt'];			id_cpt = l; l=l+1;
		script{l} = ['set ps=' prefix '.ps'];			l=l+1;
		script{l} = ['gmtset PS_MEDIA=' paper_media];	l=l+1;
	end

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o] = do_psbasemap(script, mex_sc, l, o, pack, opt_R, opt_J, opt_B, ...
	X0, Y0, opt_U, opt_P, opt_deg, opt_annotsize, frmPen, opt_frameWidth, N)
% ...
	[comm, pb, pf, do_MEX] = unpack(pack);
	if (N == 1)
		script{l} = sprintf('\n%s --- Start by creating the basemap frame.', comm);		l=l+1;
		OK = ' -K  > ';
	else
		script{l} = sprintf('\n%s --- Creating basemap of the %dth image.', comm, N);	l=l+1;
		OK = ' -K -O >> ';
	end

	if (~isempty(frmPen))
		frmPen = [pb 'framePen' pf];		% Only write if it exists
	end
	script{l} = ['gmt psbasemap ' pb 'lim' pf ' ' pb 'proj' pf ' ' pb 'frame' pf ' ' X0 ' ' Y0 opt_U opt_P ...
	             ' ' pb 'deg_form' pf ' ' pb 'annot_size' pf ' ' frmPen ' ' pb 'frame_width' pf OK pb 'ps' pf];
	l = l + 1;
	if (do_MEX)
		mex_sc{o,1} = sprintf('psbasemap %s %s %s %s %s %s %s %s %s %s %s %s', ...
					  opt_R, opt_J, opt_B, X0, Y0, opt_U, opt_P, opt_deg, opt_annotsize, frmPen, ...
					  opt_frameWidth, OK(1:end-3));		% We don't want the '>' or '>>' in 'OK'
		o = o + 1;
	end

% ------------------------------------------------------------------------------------------------------------
function [script, l, o, haveAlfa, used_grd, nameRGB] = do_grdimg(handMir, script, pack, l, o)
% Deal with the grdimage & grdgradient part

	[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir, opt_len_unit] = unpack(pack);
	haveAlfa = false;	used_grd = false;
	nameRGB = [];				% When not empty it means we'll do a screen capture ('image' or to capture transp)
	if (~isempty(handMir.grdname))
		% If renderer == OpenGL, that is interpreted as a transparency request. In that case we need a screen capture
		if (strcmpi(get(handMir.figure1,'Renderer'), 'OpenGL'))
			if (~isempty(findobj(get(handMir.axes1,'Child'),'Type','patch')))		% Extra test
				hP = findobj(get(handMir.axes1,'Child'),'Type','patch');
				for (k = 1:numel(hP))
					% 0.005 is a fake number set in DrawClosedPolygon_CB to trick the R2015 breakage in hitting patches
					if (get(hP(k), 'FaceAlpha') ~= 0.005)	% So if they are all == 0.005 ignore them as transparent
						handMir.Illumin_type = 10;			% Dumb fake value just to force screen capture
						haveAlfa = true;
					end
				end
			end
		end
		if (handMir.nLayers > 1)	% While we could try to read the layer, that's complicated. So, force screen capture
			handMir.Illumin_type = 10;
		end
		if (handMir.Illumin_type > 0 && handMir.Illumin_type <= 4)
			% We have a image illuminated with grdgradient. Rebuild de illumination
			illumComm = getappdata(handMir.figure1,'illumComm');
			opt_M = '';
			if (handMir.Illumin_type == 1 && handMir.geog),		opt_M = ' -M';	end
			if (~isempty(opt_M)),								opt_M = ' -fg';	end
			opt_N = '';
			if (handMir.Illumin_type == 1),		opt_N = ' -Nt';     end
			ind_A = strfind(illumComm, '-A');	ind_N = strfind(opt_N, ' -N');
			if (do_MEX && ~isempty(ind_A) && ~isempty(ind_N))		% Use the new grdimage -I option ... when know GMT version
				illum = sprintf(' -I+a%s+n%s', illumComm(ind_A+2:end), opt_N(4:end));
			else
				name_illum = [prefix '_intens.grd'];
				script{l} = sprintf('\n%s -------- Compute the illumination grid', comm);    l=l+1;
				script{l} = ['gmt grdgradient ' pb 'grd' pf opt_M ' ' illumComm opt_N ' -G' name_illum ellips];    l=l+1;
				illum = [' -I' name_illum];
			end
			have_gmt_illum = true;     used_grd = true;
		elseif (handMir.Illumin_type > 4 || handMir.is_draped)
			% We have a Manip or draping illumination. Here we have to use the R,G,B trick
			nameRGB = [prefix_ddir '_channel'];    name_sc = [prefix '_channel'];
			illum = [nameRGB '_r.grd ' nameRGB '_g.grd ' nameRGB '_b.grd']; % ????
			have_gmt_illum = false;
		else        % We don't have any illumination
			have_gmt_illum = false;
			used_grd = true;
		end

		if (have_gmt_illum)                     % grdimage with illumination
			script{l} = sprintf('\n%s -------- Plot the the base image using grdimage & illumination', comm);    l=l+1;
			script{l} = ['gmt grdimage ' pb 'grd' pf ' -C' pb 'cpt' pf illum ellips RJOK ' >> ' pb 'ps' pf];
			l=l+1;
			used_grd = true;
		elseif (used_grd && ~have_gmt_illum)     % Simple grdimage call
			script{l} = sprintf('\n%s -------- Plot the the base image using grdimage', comm);    l=l+1;
			script{l} = ['gmt grdimage ' pb 'grd' pf ' -C' pb 'cpt' pf ellips RJOK ' >> ' pb 'ps' pf];   l=l+1;
			used_grd = true;
		else                                    % No grd used, use the R,G,B channels
			script{l} = sprintf('\n%s -------- Plot the 3 RGB base images using grdimage', comm);    l=l+1;
			script{l} = ['gmt grdimage ' name_sc '_r.grd ' name_sc '_g.grd ' name_sc '_b.grd' ellips RJOK ' >> ' pb 'ps' pf];
			l=l+1;    
		end
	elseif (handMir.image_type == 20)
		% Do nothing regarding the basemap image (in fact we don't have any image)
	else    % We don't have a grid, so we need to fish the image and save it as R,G,B triplet
		nameRGB = [prefix_ddir '_channel'];    name_sc = [prefix '_channel'];
		script{l} = sprintf('\n%s -------- Plot the 3 RGB base images using grdimage', comm);    l=l+1;
		script{l} = ['gmt grdimage ' name_sc '_r.grd ' name_sc '_g.grd ' name_sc '_b.grd' ellips RJOK ' >> ' pb 'ps' pf];
		l = l + 1;
	end

% ------------------------------------------------------------------------------------------------------------
function [script, l, o, sc_cpt] = do_palette(handMir, script, l, o, pack, used_grd, id_cpt)
% ...
	if (used_grd || strcmp(get(handMir.PalAt,'Check'),'on') || strcmp(get(handMir.PalIn,'Check'),'on') )
		[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir, opt_len_unit] = unpack(pack);
		tmp = cell(261,1);
		pal = get(handMir.figure1,'colormap');
		%Z = getappdata(handMir.figure1,'dem_z');
		% SE Z == [] FAZER QUALQUER COISA
		if (handMir.have_nans),     cor_nan = pal(1,:);     pal = pal(2:end,:);   end     % Remove the bg color

		pal_len = size(pal,1);
		z_min = handMir.head(5);    z_max = handMir.head(6);

		dz = (z_max - z_min) / pal_len;
		tmp{1} = '# Color palette exported by Mirone';
		tmp{2} = '# COLOR_MODEL = RGB';
		cor = round(pal*255);
		for (i = 1:pal_len)
			cor_str = sprintf([num2str(cor(i,1),'%.12g') '\t' num2str(cor(i,2),'%.12g') '\t' num2str(cor(i,3),'%.12g')]);
			z1 = num2str(z_min+dz*(i-1),'%.3f');
			z2 = num2str(z_min+dz*i,'%.3f');
			tmp{i+2} = sprintf([z1 '\t' cor_str '\t' z2 '\t' cor_str]);
		end
		tmp{pal_len+3} = sprintf('F\t255\t255\t255');
		tmp{pal_len+4} = sprintf('B\t0\t0\t0');
		if (handMir.have_nans)
			cor = round(cor_nan*255);
			cor_str = sprintf(['N\t' num2str(cor(1),'%.12g') '\t' num2str(cor(2),'%.12g') '\t' num2str(cor(3),'%.12g')]);
			tmp{pal_len+5} = sprintf(cor_str);
		else
			tmp{pal_len+5} = sprintf('N\t255\t255\t255');
		end
		sc_cpt = [dest_dir prefix '.cpt'];
		fid = fopen(sc_cpt,'wt');
		for (i = 1:pal_len+5),	fprintf(fid,'%s\n',tmp{i});     end
		fclose(fid);
	else        % Remove the cpt declaration. After all we won't go to use it
		sc_cpt = '';
		script(id_cpt) = [];    l = l - 1;
	end

% ------------------------------------------------------------------------------------------------------------
function [mex_sc, o, hWait] = do_grdimg_MEX(handMir, mex_sc, o, KORJ)
% Do the grdimage stuff for MEX
	hWait = aguentabar(0,'title','Computing PDF fig');
	pad = 0;
	if (handMir.geog && isempty(strfind(KORJ, '-JX'))),		pad = 2;	end
	img = get(handMir.hImg, 'CData');
	n_band = size(img,3);
	if (n_band == 1)		% DONT KNOW IF THERE IS A BETTER SOLUTION, BUT img MUST BE RGB
		img = ind2rgb8(img, get(handMir.figure1,'Colormap'));	% img is now RGB
		n_band = 3;
	end
	n_col = size(img,2);
	new_size = [(size(img,1) + 2*pad) (n_col + 2*pad) size(img,3)];
	img2(prod(new_size)) = uint8(0);
	r = 1;		off = n_band - 1;
	n_col_pad = n_col + 2*pad;
	for (row = size(img,1):-1:1)
		k = (n_col_pad * (pad + r - 1) + pad) * n_band + 1;
		for (col = 1:n_col)
			img2(k:k+off) = img(row,col,:);
			k = k + n_band;
		end
		r = r + 1;
	end
	img2 = reshape(img2, new_size);

	I.proj4 = '';	I.wkt = '';		I.range = handMir.head(1:6);	I.inc = handMir.head(8:9);	I.nodata = NaN;
	I.registration = handMir.head(7);	I.title = '';	I.comment = '';	I.command = '';	I.datatype = 'uint8';
	I.x = linspace(handMir.head(1), handMir.head(2), size(img,2));
	I.y = linspace(handMir.head(3), handMir.head(4), size(img,1));
	I.image = img2;	I.x_unit = '';	I.y_unit = '';	I.z_unit = '';	I.colormap = [];	I.alpha = [];
	I.layout = 'BRPa';
	I.pad = pad;
	mex_sc{o,1} = ['grdimage -Dr' KORJ];
	mex_sc{o,2} = I;		o = o + 1;
	%for (k = 1:size(img,3))
	%	img(:,:,k) = reshape(img(:,:,k)', [size(img,1) size(img,2)]);
	%end
	%mex_sc{o,2} = gmt('wrapimage', img, handMir.head);		l = l + 1;
	%mex_sc{o,2}.mem_layout = 'TRBa';
	aguentabar(1/3);

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o] = do_pscoast(handles, handMir, script, mex_sc, l, o, pack, opt_R, opt_J)
% Do wrapping work around the pscoast call

	[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir] = unpack(pack);

	if (~isempty(handles.opt_psc))	% We have pscoast commands
		if (~do_MEX)
			[script, l] = do_pscoast_job(handles, handMir, script, l, o, comm, pb, pf, ellips);
		else
			[script, l, o, mex_sc] = do_pscoast_job(handles, handMir, script, l, o, comm, pb, pf, ellips, mex_sc);
			if (strfind(mex_sc{o-1, 1}, '-J '))		% Because if opt_J = -JX...d we can't just say '-J' to bring it from history
				mex_sc{o-1, 1} = strrep(mex_sc{o-1, 1}, '-J ', [opt_J ' ']);
			end
			if (~strfind(mex_sc{o-1, 1}, '-R '))	% Than it means a new -R was set in do_pscoast_job()
				mex_sc{o, 1} = ['psxy -T' opt_R ' ' opt_J ' -O -K'];	% Bogus command used only to reset the -R & -J
				o = o + 1;
			end
		end
	elseif (~isempty(handles.opt_L))
		script{l} = ['gmt psbasemap ' handles.opt_L RJOK ' >> ' pb 'ps' pf];	l = l + 1;
	end

% ------------------------------------------------------------------------------------------------------------
function [script, l, o, mex_sc] = do_pscoast_job(handles, handMir, script, l, o, comm, pb, pf, ellips, mex_sc)
% Do the actual work of writing a pscoast command

	if (nargin < 10),	mex_sc = '';	end
	script{l} = sprintf('\n%s Plot coastlines', comm);	l = l + 1;
	opt_R = ' -R';		opt_J = ' -J';
	if (~handMir.geog && handMir.is_projected)
		[xy_prj, msg] = geog2projected_pts(handMir, [handles.x_min handles.y_min; handles.x_min handles.y_max; ...
										   handles.x_max handles.y_max; handles.x_max handles.y_min;], ...
										   [get(handMir.axes1,'Xlim') get(handMir.axes1,'Ylim') 1]);
		if (~isempty(msg))
			errordlg(msg, 'Error')	% But we don't stop because of this error
		else
			opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12gr', xy_prj(1,:),xy_prj(3,:));		% Note the -R./././.r construct
		end
		projGMT = getappdata(handMir.figure1,'ProjGMT');
		if (~isempty(projGMT))		% Only simple case
			opt_J = projGMT;
		end

		% Here we need also to use the map scale in the form 1:xxxxx
		escala = get(handles.edit_scale , 'String');
		ind = strfind(script{5}, '-J');
		script{5} = [script{5}(1:ind+1) 'x' escala];	% DANGEROUS. IT RELIES ON THE INDEX 5
	end

	script{l} = ['gmt pscoast ' handles.opt_psc ellips handles.opt_L opt_R opt_J ' -O -K >> ' pb 'ps' pf];	l = l + 1;
	if (~isempty(mex_sc))
		mex_sc{o, 1} = ['pscoast ' handles.opt_psc ellips handles.opt_L opt_R opt_J ' -O -K'];			o = o + 1;
	end

	if (numel(opt_R) > 3)		% We need a trick to reset -R & -J so that the remaining commands can rely on gmt.history
		script{l} = sprintf('\n%s -------- Fake command used only to reset the -R & -J to their script defaults.', comm);    l=l+1;
		script{l} = ['gmt psxy ' pb 'lim' pf ' ' pb 'proj' pf ' -T -O -K >> ' pb 'ps' pf];
		l = l + 1;
	end

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o, used_grd, hLine, hText] = do_contour(handMir, script, mex_sc, l, o, pack, ...
	hLine, hText, used_grd)
% Handle the contour plots

	[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir] = unpack(pack);
	tag = get(hLine,'Tag');
	if (~isempty(tag) && (~isempty(handMir.grdname) || do_MEX))
		h = findobj(hLine,'Tag','contour');
		if (~isempty(h))
			h_label = findobj(hText,'Tag','contour');		% Search for contour labels
			if (~isempty(h_label))
				lab = get(h_label,'UserData');
				if (iscell(lab)),   lab = unique(cat(1,lab{:}));    end
			else
				lab = [];
			end
			conts = zeros(numel(h),1);
			for (i = 1:numel(h))
				conts(i) = getappdata(h(i),'cont_label');
			end
			conts = unique(conts);
			no_anot = setxor(conts,lab);    % Contour levels that are not annotated
			name = [prefix_ddir '_cont.dat'];
			fid = fopen(name,'wt');
			if (isempty(no_anot))           % Annotate all contours
				fprintf(fid,'%.5f\tA\n',conts);
			else                            % Annotate only some contours
				conts = [[lab; no_anot] [ones(length(lab),1)*double('A'); ones(length(no_anot),1)*double('C')]];
				conts = sortrows(conts);
				fprintf(fid,'%.5f\t%c\n',conts');
			end
			fclose(fid);
			if (~isempty(handMir.grdname))
				script{l} = sprintf('\n%s ---- Plot contours', comm);	l=l+1;
				script{l} = ['gmt grdcontour ' pb 'grd' pf ' -C' [prefix '_cont.dat'] ellips RJOK ' >> ' pb 'ps' pf];
				l = l + 1;
			end
			if (do_MEX)
				mex_sc{o,1} = ['grdcontour -C"' dest_dir prefix '_cont.dat' '"' ellips KORJ];
				[X,Y,Z] = load_grd(handMir);
				mex_sc{o,2} = gmt('wrapgrid', Z, handMir.head);
				o = o + 1;
			end
			used_grd = true;
			hLine = setxor(hLine, h);       % h is processed, so remove it from handles list
			hText = setxor(hText, h_label); % same for contour label strings
		end
	end

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o, hLine] = do_symbols(handMir, script, mex_sc, l, o, pack, hLine, opt_J)
% Deal with the symbols plotting

	haveSymbol = false;
	tag = get(hLine,'Tag');
	if (isempty(tag)),	return,		end

	h = findobj(hLine,'Tag','Symbol');
	h = [h; findobj(hLine,'Tag','City_major')];
	h = [h; findobj(hLine,'Tag','City_other')];
	h = [h; findobj(hLine,'Tag','volcano')];
	h = [h; findobj(hLine,'Tag','hotspot')];
	h = [h; findobj(hLine,'Tag','Earthquakes')];
	h = [h; findobj(hLine,'Tag','DSDP')];
	h = [h; findobj(hLine,'Tag','ODP')];

	% Search for points as Markers (that is, line with no line - just symbols on vertices)
	h_shit = get(hLine,'LineStyle');
	h_num_shit = find(strcmp('none', h_shit));
	if (h_num_shit)
		id = ismember(h, hLine(h_num_shit));		% Many, if not all, can be repeated
		h(id) = [];									% This will remove repeted elements
		h = [h; hLine(h_num_shit)];
	end

	if (~isempty(h))
		symbols = get_symbols(h);
		haveSymbol = true;
		hLine = setxor(hLine, h);			% h is processed, so remove it from handles list
	end

	if (~haveSymbol),	return,		end

	[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir, opt_len_unit] = unpack(pack);
	ns = numel(symbols.x);
	name = [prefix_ddir '_symb.dat'];	name_sc = [prefix '_symb.dat'];
	fc = symbols.FillColor{1};			ec = symbols.EdgeColor{1};
	opt_G = '';			opt_W = '';
	if (~ischar(fc)),	opt_G = sprintf(' -G%d/%d/%d', round(fc * 255));	end
	if (ischar(ec) && strcmp(ec, 'auto'))			% WRONG. Should be line's 'Color' property
		opt_W = ' -W1p';
	elseif (~ischar(ec))
		opt_W = sprintf(' -W1p,%d/%d/%d', round(ec * 255));
	end

	if (ns > 1 && numel(symbols.Size) == 1)			% We have the same symbol repeated ns times
		fid = fopen(name,'wt');
		fprintf(fid,'%.5f\t%.5f\n',[symbols.x{:}; symbols.y{:}]);
		script{l} = sprintf('\n%s ---- Plot symbols', comm);    l = l + 1;
		script{l} = ['gmt psxy ' name_sc ' -S' symbols.Marker num2str(symbols.Size{1}) 'p' opt_G ...
					 opt_W ellips RJOK ' >> ' pb 'ps' pf];		l = l + 1;
		fclose(fid);
		if (do_MEX)
			mex_sc{o,1} = ['psxy -S' symbols.Marker num2str(symbols.Size{1}) 'p' opt_G opt_W ellips KORJ];
			mex_sc{o,2} = [symbols.x{:}; symbols.y{:}];			o = o + 1;
		end
	elseif (ns == 1 && numel(symbols.Size) == 1)	% We have only one symbol
		script{l} = sprintf('\n%s  ---- Plot symbol', comm);		l=l+1;
		script{l} = [sprintf('echo %.6f\t%.6f',symbols.x{1},symbols.y{1}) ' | ' ...
					'gmt psxy -S' symbols.Marker num2str(symbols.Size{1}) 'p' opt_G ...
					opt_W ellips RJOK ' >> ' pb 'ps' pf];		l = l + 1;
		if (do_MEX)
			mex_sc{o,1} = ['psxy -S' symbols.Marker num2str(symbols.Size{1}) 'p' opt_G opt_W ellips KORJ];
			mex_sc{o,2} = [symbols.x{1}; symbols.y{1}];			o = o + 1;
		end
	else								% We have ns different symbols
		m = zeros(ns,1);
		for (i = 1:ns)
			m(i) = size(symbols.x{i},2);
		end
		n = find(m ~= 1);
		if (~isempty(n))				% We have a mixed scenario. Individual as well as group symbols
			script = write_group_symb(prefix,prefix_ddir,comm,pb,pf,ellips,symbols,n,script, opt_J);
			symbols.x(n) = [];			symbols.y(n) = [];  % Clear processed symbols
			symbols.FillColor(n) = [];	symbols.EdgeColor(n) = [];
			symbols.Size(n) = [];		symbols.Marker(n,:) = [];
			l = numel(script) + 1;
			ns = ns - numel(n);
		end
		fid = fopen(name,'wt');
		for (i = 1:ns)
			fc = symbols.FillColor{i};			ec = symbols.EdgeColor{i};
			opt_G = '';			opt_W = '';
			if (~ischar(fc)),	opt_G = sprintf(' -G%d/%d/%d', round(fc * 255));	end
			if (ischar(ec) && strcmp(ec, 'auto'))			% WRONG. Should be line's 'Color' property
				opt_W = ' -W1p';
			elseif (~ischar(ec))
				opt_W = sprintf(' -W1p,%d/%d/%d', round(ec * 255));
			end
			fprintf(fid,'>%s\n',[opt_G opt_W]);
			fprintf(fid,'%.5f\t%.5f\t%.0f\t%s\n',symbols.x{i},symbols.y{i},symbols.Size{i},symbols.Marker(i,:));
		end
		fclose(fid);
		script{l} = ' ';                        	l = l + 1;
		script{l} = [comm ' ---- Plot symbols'];    l = l + 1;
		script{l} = ['gmt psxy ' name_sc ellips opt_len_unit RJOK ' -S >> ' pb 'ps' pf];		l = l + 1;
		if (do_MEX)
			mex_sc{o,1} = ['psxy -S "' dest_dir name_sc '"' ellips opt_len_unit KORJ];		o = o + 1;
		end
	end

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o, hLine, hPatch] = do_meca(handMir, script, mex_sc, l, o, pack, hLine, ...
	hPatch, opt_J)
% Handle focal mechanisms plotting

	if (~isempty(hPatch))
		focHand = findobj(hPatch,'Tag','FocalMeca');
		if (~isempty(focHand))
			[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir, opt_len_unit] = unpack(pack);
			% First deal with the 'line anchors'
			focHandAnchor = findobj(hLine,'Tag','FocalMecaAnchor');   % Handles of the line anchors
			x = get(focHandAnchor,'XData');			y = get(focHandAnchor,'YData');
			if (iscell(x)),		x = cell2mat(x);	y = cell2mat(y);    end
			id_anch = find(diff(x,1,2));
			if (isempty(id_anch)),		id_anch = find(diff(y,1,2));	end		% Rare cases were movement was vertical

			psmeca_line = cell(numel(focHand),1);
			for (k = 1:numel(focHand)),		psmeca_line{k} = getappdata(focHand(k),'psmeca_com');	end
			psmeca_line = cat(1,psmeca_line{:});    % This also get us rid of empty cell fields.

			if (isequal(psmeca_line(1,1:2), [x(end,1) y(end,1)]))	% Don't know why but both vars have reverse order
				x = x(end:-1:1, :);		y = y(end:-1:1, :);
			end
			n_cols = size(psmeca_line,2);
			if (n_cols == 10 || n_cols == 14),		with_label = 1;
			else,									with_label = 0;
			end
			name = [prefix_ddir '_meca.dat'];   name_sc = [prefix '_meca.dat'];     opt_C = '';
			fid = fopen(name,'wt');
			if (n_cols == 9 || n_cols == 10)		% Aki & Richard convention
				% If beach-bals are not ploted at their origin update the ploting coords columns
				if (~isempty(id_anch))
					psmeca_line(:,8) = x(:,2);		psmeca_line(:,9) = y(:,2);     opt_C = ' -C';
				end
				opt_S = ['-Sa' getappdata(handMir.figure1,'MecaMag5') 'c'];
				format = '%.4f\t%.4f\t%.1f\t%.0f\t%.0f\t%.0f\t%.1f\t%.4f\t%.4f';
				for (k=1:size(psmeca_line,1))
					fprintf(fid,format,psmeca_line(k,1:9));
					if (with_label),	fprintf(fid,'\t%s\n',num2str(psmeca_line(k,10)));
					else,				fprintf(fid,'\n');
					end
				end
			elseif (n_cols == 13 || n_cols == 14)	% CMT convention
				% If beach-bals are not ploted at their origin update the ploting coords columns
				if (~isempty(id_anch))
					psmeca_line(:,12) = x(:,2);		psmeca_line(:,13) = y(:,2);     opt_C = ' -C';
				end
				psmeca_line(:,11) = psmeca_line(:,11) + 7;		% psmeca uses Moment in Dyn-cm
				opt_S = ['-Sc' getappdata(handMir.figure1,'MecaMag5') 'c'];
				format = '%.4f\t%.4f\t%.1f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.2f\t%d\t%.4f\t%.4f';           
				for (k=1:size(psmeca_line,1))
					fprintf(fid,format,psmeca_line(k,1:13));
					if (with_label),	fprintf(fid,'\t%s\n',num2str(psmeca_line(k,14)));
					else,				fprintf(fid,'\n');
					end
				end
			end
			fclose(fid);
			script{l} = sprintf('\n%s ---- Plot Focal Mechanisms', comm);   l=l+1;
			script{l} = ['gmt psmeca ' opt_S opt_C ' ' name_sc ellips RJOK ' >> ' pb 'ps' pf];		l = l + 1;
			if (do_MEX)
				mex_sc{o,1} = ['psmeca "' dest_dir name_sc '"' opt_S opt_C ' ' ellips KORJ];	o = o + 1;
			end
			hPatch = setxor(hPatch, focHand);		% focHand is processed, so remove it from handles list
			hLine  = setxor(hLine, focHandAnchor);	%       iden
		end
	end

% ------------------------------------------------------------------------------------------------------------
function [script, l, hPatch] = do_telhas(script, l, hPatch, pack)
% Handles the "Telhas", but telhas is no mex

	if (~isempty(hPatch))
		hTelhas = findobj(hPatch,'Tag','tapete');
		if (~isempty(hTelhas))
			[comm, pb, pf, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir] = unpack(pack);
			tmp = findobj(hPatch,'Tag','tapete_R');
			hPatch = setxor(hPatch, tmp);       % Remove the reverse "telhas" 
			n_tapetes = length(hTelhas);
			for (i=1:n_tapetes)
				saved = get(hTelhas(i),'UserData');
				name = [prefix_ddir sprintf('_telha_%d.dat',i)];
				name_sc = [prefix sprintf('_telha_%d.dat',i)];
				if (~isempty(saved))
					fid = fopen(name,'wt');
					fprintf(fid,'%.5f\t%.5f\n',saved.line');
					fclose(fid);
					script{l} = sprintf('\n%s ---- Plot telhas. NOTE: THIS IS NOT A GMT PROGRAM', comm);   l=l+1;
					script{l} = ['telha ' name_sc ' ' saved.opt_E ' ' saved.opt_I ' ',...
					             saved.opt_N ' ' saved.opt_T ' -Blixo.dat'];
%					mex_sc{o,1} = sprintf('telha		% BUT TELHA IS NOT A MEX
					l = l + 1;
					script{l} = ['psxy lixo.dat ' ellips RJOK ' -L >> ' pb 'ps' pf];
%					mex_sc{o,1} = ['psxy lixo.dat -L ' ellips KORJ];
					l = l + 1;
				end
			end
			hPatch = setxor(hPatch, hTelhas);       % hTelhas is processed, so remove it from handles list
		end
	end

% ------------------------------------------------------------------------------------------------------------
function [script, l, haveAlfa, hPatch] = do_countries(handMir, script, l, pack, haveAlfa, hPatch)
% Deal with country plots but this have to be replaced by a pscoast call

	if (~isempty(hPatch))
		% First see about these still remaining patch transparency
		% But if both write script and auto PDF and transparencies, the script will probably be wrong
		[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir] = unpack(pack);
		if (~do_MEX)
			[hPatch, hAlfaPatch] = findTransparents(hPatch);
			if (isempty(hAlfaPatch)),		haveAlfa = false;		end			% An extra test
		end

		hAtlas = findobj(hPatch,'Tag','Atlas');
		if (~isempty(hAtlas))
			n_cts = length(hAtlas);
			if (n_cts > 1)                  % We have multiple countries
				ct_names = cell(n_cts,1);   % To hold the country names
				for (k=1:n_cts)             % Loop over all countries found and store theyr names
					ct_names{k} = get(hAtlas(k),'UserData');
				end
			else                            % We have only one country
				ct_names = {get(hAtlas,'UserData')};
			end
			ct_names = unique(ct_names);    % Many countries have several polygons (e.g. islands).
			name = [prefix_ddir '_country_names.txt'];
			fid = fopen(name,'wt');
			fprintf(fid,'%s\n',ct_names{:});        fclose(fid);
			script{l} = sprintf('\n%s ---- Plot countries. NOTE: THIS IS NOT A GMT PROGRAM', comm);   l=l+1;
			ct_with_pato = getappdata(handMir.figure1,'AtlasResolution');
			script{l} = [cd filesep 'country_extract -P' name ' ' ct_with_pato ' -C | ',...
					'gmt psxy -W0.5p ' ellips RJOK ' >> ' pb 'ps' pf];
			l = l + 1;
			hPatch = setxor(hPatch, hAtlas);       % AtlasHand is processed, so remove it from handles list
		end
	end

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o, hPatch, xx, yy] = do_bargraphs(handMir, script, mex_sc, l, o, pack, hPatch)
% ...
	xx = [];	yy = [];
	if (~isempty(hPatch))
		thisHand = findobj(hPatch,'Tag','BarGraph');
		if (~isempty(thisHand))
			[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir, opt_len_unit] = unpack(pack);
			xx = get(thisHand,'XData');     yy = get(thisHand,'YData');
			bar_W = xx(4,1) - xx(1,1);
			xx = (xx(1,:) + xx(4,:)) / 2;		% Coordinates at the bars half width
			yy = yy(2, :);

			name = [prefix_ddir '_bars.dat'];	name_sc = [prefix '_bars.dat'];
			FillColor = get(thisHand,'FaceColor');		LineWidth = get(thisHand,'LineWidth');
			EdgeColor = get(thisHand,'EdgeColor');
			cor_fill = sprintf('%d/%d/%d', round(FillColor * 255));
			cor_edge = sprintf('%d/%d/%d', round(EdgeColor * 255));
			cor = [' -G' cor_fill ' -W' num2str(LineWidth) 'p,' cor_edge];
			fid = fopen(name,'wt');
			fprintf(fid,'%f\t%f\n',[xx; yy]);
			fclose(fid);
			script{l} = sprintf('\n%s ---- Plot Bar graph', comm);		l = l + 1;
			opt_S = sprintf(' -Sb%fu', bar_W);
			script{l} = ['gmt psxy ' name_sc opt_len_unit opt_S RJOK cor ' >> ' pb 'ps' pf];	l = l + 1;
			if (do_MEX)
				mex_sc{o,1} = ['psxy ' opt_len_unit opt_S KORJ cor];
				mex_sc{o,2} = [xx(:) yy(:)];		o = o + 1;
			end
			hPatch = setxor(hPatch, thisHand);       % thisHand is processed, so remove it from handles list
		end
	end
	
% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o, hPatch] = do_histograms(handMir, script, mex_sc, l, o, pack, hPatch, xx, yy)
% ...
	if (~isempty(hPatch))
		thisHand = findobj(hPatch,'Tag','Histogram');
		if (~isempty(thisHand))
			[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir, opt_len_unit] = unpack(pack);
			xh = get(thisHand,'XData');			% 2xN matrix
			bar_W = xh(4,1) - xh(1,1);

			name = [prefix_ddir '_histo.dat'];	name_sc = [prefix '_histo.dat'];
			FillColor = get(thisHand,'FaceColor');		LineWidth = get(thisHand,'LineWidth');
			EdgeColor = get(thisHand,'EdgeColor');
			y = getappdata(thisHand, 'xy');	% Have to store it here in order that GMT can reconstruct the command
			fid = fopen(name,'wt');
			fprintf(fid,'%f\n', y);
			fclose(fid);
			script{l} = sprintf('\n%s ---- Plot Histogram', comm);   l=l+1;
			opt_G = sprintf(' -G%d/%d/%d', round(FillColor * 255));
			opt_W = sprintf(' -W%g', bar_W);
			opt_L = sprintf(' -L%gp,%d/%d/%d', LineWidth, round(EdgeColor * 255));
			script{l} = ['gmt pshistogram ' name_sc opt_len_unit opt_G opt_W opt_L RJOK ' -F >> ' pb 'ps' pf];
			l = l + 1;
			if (do_MEX)
				mex_sc{o,1} = ['pshistogram -F ' opt_len_unit opt_G opt_W opt_L KORJ];
				mex_sc{o,2} = [xx(:) yy(:)];		o = o + 1;
			end
			hPatch = setxor(hPatch, thisHand);       % thisHand is processed, so remove it from handles list
		end
	end

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o] = do_patches(handMir, script, mex_sc, l, o, pack, hPatch, writeScript)
% ...
	if (~isempty(hPatch))
		[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir, opt_len_unit] = unpack(pack);
		xx = get(hPatch,'XData');     yy = get(hPatch,'YData');
		n_patch = numel(hPatch);
		%LineStyle = get(ALLpatchHand,'LineStyle');
		LineWidth = get(hPatch,'LineWidth');
		if (iscell(LineWidth)),     LineWidth = cat(1,LineWidth{:});     end
		EdgeColor = get(hPatch,'EdgeColor');
		if (iscell(EdgeColor)),     EdgeColor = cat(1,EdgeColor{:});     end
		FillColor = get(hPatch,'FaceColor');
		if (iscell(FillColor))
			if (handMir.version7 >= 8.4)	% Must check if have to undo a trick to workarround a R2015 bug
				for (k = 1:n_patch)
					if (get(hPatch(k), 'FaceAlpha') == 0.005)
						FillColor{k} = 'none';
					end
				end
			end
			resp = strcmp('none',FillColor);
			if (~any(resp))
				FillColor = cat(1,FillColor{:});
			else
				for (i = 1:numel(resp))					% Signal down if this is a non colored polygon
					if (resp(i))
						FillColor{i} = [-1 -1 -1];		% This is a non-colored one
					end
				end
				FillColor = cat(1,FillColor{:});
			end
		else				% We have only one patch
			xx = num2cell(xx,1);   yy = num2cell(yy,1);	% Make it a cell for reducing the head-hakes
			resp = strcmp('none',FillColor);
			if (resp || get(hPatch, 'FaceAlpha') == 0.005)	% The 0.005 is the flag to workaround R2015 bug
				FillColor = [-1 -1 -1];					% Signal down that this is a non colored polygon
			end
		end
		name = [prefix_ddir '_patch.dat'];		name_sc = [prefix '_patch.dat'];
		if (writeScript),	fid = fopen(name,'wt');		end
		for (i = 1:n_patch)
			cor_edge = sprintf('%d/%d/%d', round(EdgeColor(i,1:3) * 255));
			if (FillColor(i,1) >= 0)		% Color filled polygon
				cor_fill = sprintf('%d/%d/%d', round(FillColor(i,1:3) * 255));
				mlt_comm = ['> -G' cor_fill ' -W' num2str(LineWidth(i)) 'p,' cor_edge];
			else							% No filling color
				mlt_comm = ['> -G- -W' num2str(LineWidth(i)) 'p,' cor_edge];
			end

			if (writeScript)
				if (any(isnan(xx{i})))      % If we have NaNs we need to split into segments
					[latcells,loncells] = aux_funs('polysplit', yy{i}(:),xx{i}(:));
					for (j = 1:numel(loncells))
						fprintf(fid,'%s\n',mlt_comm);
						fprintf(fid,'%.5f\t%.5f\n',[loncells{j}(:)'; latcells{j}(:)']);
					end
				else
					fprintf(fid,'%s\n',mlt_comm);
					fprintf(fid,'%.5f\t%.5f\n',[xx{i}(:)'; yy{i}(:)']);
				end
			end
			if (do_MEX)
				transp = get(hPatch(i), 'FaceAlpha');		opt_t = '';
				if (transp > 0.005),	opt_t = sprintf(' -t%d', round((1-transp) * 100));	end		% Have transparency
		 		mex_sc{o,1} = ['psxy ' mlt_comm(3:end) ellips opt_len_unit opt_t KORJ];
				mex_sc{o,2} = [xx{i}(:) yy{i}(:)];		o = o + 1;
			end
		end
		if (writeScript),	fclose(fid);	end
		script{l} = sprintf('\n%s ---- Plot closed AND colored polygons', comm);		l = l + 1;
		script{l} = ['gmt psxy ' name_sc ellips opt_len_unit RJOK ' >> ' pb 'ps' pf];		l = l + 1;
	end

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o, hLine] = do_psimage(handles, script, mex_sc, l, o, pack, hLine)
% Deal with rectangles that might have an assiciated PS (like the custom PS symbols)
	if (isempty(hLine)),	return,		end			% Some stray shit
	[comm, pb, pf, do_MEX, ellips, RJOK, KORJ] = unpack(pack);
	c = false(1, numel(hLine));
	for (i = 1:numel(hLine))
		cs_fname = getappdata(hLine(i), 'cust_symb');
		if (isempty(cs_fname)),		continue,	end			% Just a generic element.
		if (i == 1)
			script{l} = sprintf('\n%s ---- Plot GMT a PS with psimage', comm);  l = l + 1;
		end
		xx = get(hLine(i),'XData');		yy = get(hLine(i),'YData');
		unit = handles.which_unit(1);
		script{l} = sprintf('gmt psimage %s %s -Dx%.4g%c/%.4g%c+w%.4g%c >> %sps%s', cs_fname, RJOK, xx(1), ...
		                    unit, yy(1), unit, xx(3)-xx(2), unit, pb, pf);		l = l + 1;
		if (do_MEX)
			mex_sc{o,1} = sprintf('psimage %s %s -Dx%.4g%c/%.4g%c+w%.4g%c', cs_fname, KORJ, xx(1), ...
			                      unit, yy(1), unit, xx(3)-xx(2), unit);		o = o + 1;
		end
		c(i) = true;
	end
	hLine(c) = [];			% Delete these since we don't want to plot them

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o, hLine] = do_custom_symbols(handles, script, mex_sc, l, o, pack, hLine)
% ...
	if (isempty(hLine)),	return,		end			% Some stray shit
	[comm, pb, pf, do_MEX, ellips, RJOK, KORJ] = unpack(pack);
	c = false(1, numel(hLine));
	for (i = 1:numel(hLine))
		cs_fname = getappdata(hLine(i), 'cust_symb');
		if (isempty(cs_fname)),		continue,	end			% Just a regular element.

		if (i == 1)
			script{l} = sprintf('\n%s ---- Plot GMT custom symbols', comm);  l = l + 1;
		end
		[PATH,FNAME] = fileparts(cs_fname);
		cs_fname = [PATH filesep FNAME];		% Must remove the extension
		xx = get(hLine(i),'XData');		yy = get(hLine(i),'YData');
		x_max = max(xx);		x_min = min(xx);
		y_max = max(yy);		y_min = min(yy);
		sym_width = x_max - x_min;
		x0 = x_min + sym_width/2;	y0 = y_min + (y_max - y_min)/2;
		m_width = handles.x_max - handles.x_min;
		w = sym_width / m_width * str2double(handles.map_width);	% OK, and if handles.which_unit(1) is not cm?
		script{l} = sprintf('echo %0.10g %0.10g | gmt psxy %s -Sk%s/%f%c >> %sps%s', ...
							x0,y0, RJOK, cs_fname, w, handles.which_unit(1), pb, pf);		l = l + 1;
		if (do_MEX)
			mex_sc{o,1} = sprintf('psxy %s -Sk"%s"/%f%c', KORJ, cs_fname, w, handles.which_unit(1));
			mex_sc{o,2} = [x0 y0];		o = o + 1;
		end
		c(i) = true;
	end
	hLine(c) = [];		% Delete these since we don't want to plot them

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o, hLine] = do_lines(script, mex_sc, l, o, pack, hLine)
% ...
	if (isempty(hLine)),	return,		end			% Some stray shit
	[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir, opt_len_unit] = unpack(pack);
	xx = get(hLine,'XData');		yy = get(hLine,'YData');
	if (~iscell(xx))				% We have only one line
		xx = num2cell(xx(:),1);		yy = num2cell(yy(:),1);
	end

	script{l} = sprintf('\n%s ---- Plot lines', comm);  l=l+1;
	LineStyle = get(hLine,'LineStyle');
	[LineStyle,LineStyle_gmt] = lineStyle2num(LineStyle);
	LineWidth = get(hLine,'LineWidth');
	if (iscell(LineWidth)),     LineWidth = cat(1,LineWidth{:});    end
	LineColor = get(hLine,'Color');
	if (iscell(LineColor)),     LineColor = cat(1,LineColor{:});    end
	[b,m] = sortrows([LineWidth LineColor LineStyle]);
	m = m(end:-1:1);			% Revert order because I want thicker lines ploted first --- WHY???
	xx = xx(m);     yy = yy(m);
	LineWidth = LineWidth(m,:);     LineColor = LineColor(m,:);
	LineStyle = LineStyle(m);       LineStyle_gmt = LineStyle_gmt(m,:);
	[b,m] = unique([LineWidth LineColor LineStyle],'rows');   % reuse b,m
	m = m(end:-1:1);			% OK, now we have to put it back in ascending order        
	m = [0; m];					% use this first index to help file creation algo
	for (i = 1:length(m)-1)
		name = sprintf('%s_line_%d.dat', prefix_ddir, i);
		name_sc = sprintf('%s_line_%d.dat', prefix, i);
		fid = fopen(name,'wt');
		for (j = m(i)+1:m(i+1))
			if (any(isnan(xx{j})))          % If we have NaNs we need to split into segments
				[latcells,loncells] = aux_funs('polysplit', yy{j}(:),xx{j}(:));
				for (k = 1:numel(loncells))
					fprintf(fid,'>\n');
					fprintf(fid,'%.5f\t%.5f\n',[loncells{k}(:)'; latcells{k}(:)']);
				end
			else
				fprintf(fid,'>\n');
				fprintf(fid,'%.5f\t%.5f\n',[xx{j}(:) yy{j}(:)]');
			end
		end
		fclose(fid);
		cor = round(LineColor(j,:) * 255);
		cor = [num2str(cor(1)) '/' num2str(cor(2)) '/' num2str(cor(3))];
		script{l} = ['gmt psxy ' name_sc ellips ' -W' num2str(LineWidth(j)) 'p,' ...
					 cor LineStyle_gmt{j} opt_len_unit RJOK ' >> ' pb 'ps' pf];	l = l + 1;
		if (do_MEX)
			mex_sc{o,1} = ['psxy "' dest_dir name_sc '"' ellips ' -W' num2str(LineWidth(j)) 'p,' cor LineStyle_gmt{j} opt_len_unit KORJ];
			o = o + 1;
		end
	end

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o, hText] = do_text(script, mex_sc, l, o, pack, hText)
% ...
	if (isempty(hText)),	return,		end

	[comm, pb, pf, do_MEX, ellips, RJOK, KORJ] = unpack(pack);
	pos = get(hText,'Position');      %font = get(ALLtextHand,'FontName');
	fsize = get(hText,'FontSize');    fcolor = get(hText,'Color');
	% Convert to string right away
	if (isa(fsize, 'cell'))
		for (k = 1:numel(fsize))
			fsize{k} = sprintf('%d', fsize{k});
		end
	else
		fsize = sprintf('%d', fsize);
	end

	% Find the Hor/Vert alignment
	HA = get(hText, 'HorizontalAlignment');
	VA = get(hText, 'VerticalAlignment');
	if (isa(HA,'cell'))		% Get only the first char of each row
		HA_ = char(zeros(numel(HA),1));
		VA_ = char(zeros(numel(HA),1));
		for (k = 1:numel(HA))
			HA_(k) = upper(HA{k}(1));
			VA_(k) = upper(VA{k}(1));
		end
		HA = HA_;	VA = VA_;
	else
		HA = upper(HA(1));		VA = upper(VA(1));
	end
	ind = (HA ~= 'L') & (HA ~= 'C') & (HA ~= 'R');		% Equivalent to GMT's 'LCR'
	if (any(ind))
		HA(ind) = 'L';		% 'Others' get 'L(eft)'
	end
	ind = (VA ~= 'T') & (VA ~= 'M') & (VA ~= 'M');		% Equivalent to GMT's 'TMB'
	if (any(ind))
		VA(ind) = 'B';		% 'Others' get 'B(ottom)'
	end
	HV = [HA VA];			% Hor/Ver justification code

	fcolor_s = {''};
	if (isnumeric(fcolor))
		fcolor = round(fcolor * 255);
		if (numel(fcolor) == 1 && fcolor ~= 0)
			fcolor_s{1} = sprintf('%d/%d/%d', fcolor, fcolor, fcolor);
		elseif (~isequal(fcolor, [0 0 0]))
			fcolor_s{1} = sprintf('%d/%d/%d', fcolor(1:3));
		end
	elseif (ischar(fcolor))		% Shit, we have to decode the color letter
		switch fcolor
			case 'w',		fcolor_s{1} = '255/255/255';
			case 'y',		fcolor_s{1} = '255/255/0';
			case 'c',		fcolor_s{1} = '0/255/255';
			case 'r',		fcolor_s{1} = '255/0/0';
			case 'g',		fcolor_s{1} = '0/255/0';
			case 'b',		fcolor_s{1} = '0/0/255';
		end
	elseif (iscell(fcolor))			% Double shit, we have to convert a Mx3 cell matrix into texts
		tmp = cell2mat(fcolor) * 255;
		fcolor_s = cell(size(tmp,1),1);
		for (m = 1:size(tmp,1))
			if (isequal(tmp(m,:), [0 0 0])),	fcolor_s{m} = '';	continue,	end
			fcolor_s{m} = sprintf('%d/%d/%d', tmp(m,1:3));
		end
	end
	str = get(hText,'String');		angle = get(hText,'Rotation');
	if (~iscell(pos))				% Make them cells for author's mental sanity
		pos = {pos};				fsize = {fsize};		angle = {angle};
		str = {str};				%font = {font};
	end

	n_text = numel(str);
	script{l} = sprintf('\n%s ---- Plot text strings', comm);   l=l+1;    
	for (i = 1:n_text)
		% Quick and dirty patch for when opt_G is a cell of cells and it than crash below on sprintf
		if (~isempty(fcolor_s{i})),	fsize{i} = [fsize{i} ',' fcolor_s{i}];	end
		opt_F = sprintf('-F+f%s+a%g+j%s', fsize{i}, angle{i}, HV);
		script{l} = sprintf('echo %.5f %.5f %s | gmt pstext %s %s %s >> %sps%s', pos{i}(1), pos{i}(2), str{i}(1,:), opt_F, ellips, RJOK, pb, pf);
		l = l + 1;
		if (do_MEX)
			mex_sc{o,1} = sprintf('pstext %s %s %s', ellips, opt_F, KORJ);
			mex_sc{o,2} = struct('data',pos{i}(1:2), 'text',str{i}(1,:));
			o = o + 1;
		end
		this_nLines = size(str{i},1);		% How many lines has this text element?
		if (this_nLines > 1)				% More than one. So try to estimate each line Pos from a simple calculus
			ext = get(hText(i), 'Extent');
			for (k = 2:this_nLines)
				yPos = pos{i}(2) - (k - 1) * (ext(4) / this_nLines);
				script{l} = sprintf('echo %.5f %.5f %s | gmt pstext %s %s %s >> %sps%s', pos{i}(1), yPos, str{i}(k,:), opt_F, ellips, RJOK, pb, pf);
				l = l + 1;
				if (do_MEX)
					mex_sc{o,1} = sprintf('pstext %s %s %s', ellips, opt_F, KORJ);
					mex_sc{o,2} = struct('data', [pos{i}(1) yPos], 'text',str{i}(k,:));
					o = o + 1;
				end
			end
		end
	end

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o] = do_colorbar(handles, handMir, script, mex_sc, l, o, pack, sc_cpt)
% Handle psscale plottings

	if (strcmp(get(handMir.PalAt,'Check'),'on') || strcmp(get(handMir.PalIn,'Check'),'on'))
		[comm, pb, pf, do_MEX] = unpack(pack);
		if (strcmp(get(handMir.PalAt,'Check'),'on')),	axHandle = get(handMir.PalAt,'UserData');
		else,											axHandle = get(handMir.PalIn,'UserData');
		end

		if (isa(axHandle, 'cell')),		axHandle = axHandle{1};		end		% After the ML great breakage it may be a cell
		axUnits = get(axHandle(1), 'Units');		set(axHandle(1), 'Units', 'pixels');
		posCB = get(axHandle(1),'pos');				set(axHandle(1), 'Units', axUnits);
		axUnits = get(handMir.axes1, 'Units');		set(handMir.axes1, 'Units', 'pixels');
		posAx = get(handMir.axes1,'pos');			set(handMir.axes1, 'Units', axUnits);

		mapW = str2double(get(handles.edit_mapWidth,'String'));
		mapH = str2double(get(handles.edit_mapHeight,'String'));
		cbH = posCB(4) / posAx(4) * mapH;       % Estimate the colorbar height like this
		marg = 0.3;     cbW = 0.5;    % Margin between Image and colorbar (in cm) and colorbar width
		unitC = handles.which_unit(1);

		if (handles.which_unit(1) == 'i')
			marg = marg / 2.54;			cbW = cbW / 2.54;
		end
		if (handles.which_unit(1) == 'p')
			marg = marg / 2.54 * 72;	cbW = cbW / 2.54 * 72;
		end
		YTick = get(axHandle(1),'YTick');		bInt = YTick(2) - YTick(1);		% To use in -B option
		opt_D = sprintf(' -D%.2f%c/%.2f%c/%.2f%c/%.2f%c',mapW+marg,unitC, cbH/2,unitC, cbH,unitC, cbW,unitC);
		script{l} = sprintf('\n%s ---- Plot colorbar ---', comm);   l=l+1;
		script{l} = ['gmt psscale' opt_D ' -S -C' pb 'cpt' pf ' -B' num2str(bInt) ' -O -K >> ' pb 'ps' pf];
		l = l + 1;
		if (do_MEX)
			mex_sc{o,1} = sprintf('psscale %s -S -C"%s" -B%g -O -K', opt_D, sc_cpt, bInt);	o = o + 1;
		end
	end

% ------------------------------------------------------------------------------------------------------------
function do_screncapture(handMir, hAlfaPatch, haveAlfa, do_writeScript, nameRGB)
% ...
	if (~isempty(nameRGB) && ~haveAlfa && do_writeScript)
        mirone('File_img2GMT_RGBgrids_CB', handMir, 'image', nameRGB)
	elseif (~isempty(nameRGB) && haveAlfa && do_writeScript)
		% Here we'll hide everything except the patches with transparency
		ALLlineHand = findobj(get(handMir.axes1,'Child'),'Type','line');
		ALLpatchHand = findobj(get(handMir.axes1,'Child'),'Type','patch');
		ALLtextHand = findobj(get(handMir.axes1,'Child'),'Type','text');
		set(ALLlineHand, 'Vis', 'off');		set(ALLpatchHand, 'Vis', 'off');	set(ALLtextHand, 'Vis', 'off')
		set(hAlfaPatch, 'Vis', 'on')					% Only semi-transparent ones are visible now
		try
			refresh(handMir.figure1)		% F... Matlab OpenGL driver has more bugs than a dead rat
			if (isempty(ALLlineHand) && isempty(ALLpatchHand) && isempty(ALLtextHand))	% No need to SC because image is clean
	        	mirone('File_img2GMT_RGBgrids_CB', handMir, 'image', nameRGB)
			else
				mirone('File_img2GMT_RGBgrids_CB', handMir, 'fromWS', nameRGB)
			end
		end
		% Make everybody visible again
		set(ALLlineHand, 'Vis', 'on');		set(ALLpatchHand, 'Vis', 'on');	set(ALLtextHand, 'Vis', 'on')
	end

% ------------------------------------------------------------------------------------------------------------
function [script, mex_sc, l, o] = do_magbars(handMir, script, mex_sc, l, o, pack, saveBind, opt_J)
% ...
	if (handMir.IamXY && strcmp(get(handMir.axes2, 'Vis'), 'on'))
		[comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir] = unpack(pack);
 		hMagBar = findobj(handMir.axes2, 'type', 'patch');
		if (~isempty(hMagBar))
			name = sprintf('%s_magbar.dat', prefix_ddir);
			name_sc = sprintf('%s_magbar.dat', prefix);
			xx = get(hMagBar,'XData');     yy = get(hMagBar,'YData');
			cor = get(hMagBar, 'FaceVertexCData');
			fid = fopen(name,'wt');
			for (i = 1:size(xx,2))
				if (cor(i,1) == 0)
					fprintf(fid, '> -G0\n');
					fprintf(fid, '%.5f\t%.1f\n',[xx(:,i) yy(:,i)]');
				end
			end
			fclose(fid);
		end
		opt_R = sprintf(' -R%.12g/%.12g/0/1', get(handMir.axes2,'xlim'));	% We need new limits here (different axes)
		i = strfind(opt_J, '/');
		Y0 = sprintf(' -Y%.12gc', str2double(opt_J(i+1:end-1))+0.1);	% Vertical offset equal to frame height + 0.1 cm
		opt_J = sprintf(' %s0.6c', opt_J(1:i));
		if (strcmpi(script{saveBind}(end), 'n')),	script{saveBind}(end) = [];		end		% Don't want top frame line
		script{l} = sprintf('\n%s ---- Plot the magnetic reversals bars (positives only)', comm);   l=l+1;
		script{l} = ['gmt psxy ' name_sc opt_R opt_J Y0 ' -O -K >> ' pb 'ps' pf];
		if (do_MEX)
			mex_sc{o,1} = ['psxy ' dest_dir name_sc opt_R opt_J Y0 ' -O -K'];
		end
	end

% ------------------------------------------------------------------------------------------------------------
function [comm, pb, pf, do_MEX, ellips, RJOK, KORJ, dest_dir, prefix, prefix_ddir, opt_len_unit] = unpack(pack)
%
	comm = pack.comm;
	pb = pack.pb;
	pf = pack.pf;
	do_MEX = pack.do_MEX;
	ellips = pack.ellips;
	RJOK = pack.RJOK;
	KORJ  = pack.KORJ;
	dest_dir = pack.dest_dir;
	prefix = pack.prefix;
	prefix_ddir = pack.prefix_ddir;
	opt_len_unit = pack.opt_len_unit;

% ------------------------------------------------------------------------------------------------------------
function gsimage(handles, script, hWait)
% ...
	% Start by deleting the empties
	c = false(size(script,1), 1);
	for (k = 1:size(script,1))
		if (isempty(script{k,1})),	c(k) = true;	end
	end
	script(c,:) = [];
	script{end,1} = strrep(script{end,1}, ' -K', '');		% And remove last -K

	for (k = 1:size(script,1))
		if (isempty(script{k,2}))
			gmtmex(script{k,1});
		else
			gmtmex(script{k,1}, script{k,2});
			if (strncmp(script{k,1}, 'grdim', 5))
				aguentabar(2/3);
			end
		end
	end
	if (~isempty(hWait)),	aguentabar(1);		end

	val = get(handles.popup_figFormat, 'Val');
	str = get(handles.popup_figFormat, 'Str');	EXT = str{val};
	if (strcmp(EXT, 'pdf')),		opt_T = '-Tf';
	elseif (strcmp(EXT, 'png')),	opt_T = '-Tg';
	elseif (strcmp(EXT, 'jpg')),	opt_T = '-Tj';
	elseif (strcmp(EXT, 'ps')),		opt_T = '-Te';
	end

	fname = [handles.path_tmp 'auto.' EXT];
	if (get(handles.check_trimWhite, 'Val'))
		gmtmex(['psconvert = ' opt_T ' -A0.5p -F' fname]);
	else
		gmtmex(['psconvert = ' opt_T ' -F' fname]);
	end
	if (handles.IamCompiled)
		win_open_mex(handles.path_tmp, ['auto.' EXT]);
	else
		if ispc
			winopen(fname);
		elseif strncmp(computer,'MAC',3) 
			unix(['open "' fname '" &']);
		else
			warndlg(sprintf('Automatically opening of PDF files not implemented in Linux\nOpen it youself in\n%s', fname), 'Warning')
		end
	end

% ------------------------------------------------------------------------------------------------------------
function symbol = get_symbols(hand)
	xx = get(hand,'XData');     yy = get(hand,'YData');
	if (~iscell(xx))
		xx = num2cell(xx,1);   yy = num2cell(yy,1);   % Make it a cell for reducing the head-hakes
	end
	symbol.x = xx(:);       symbol.y = yy(:);
	symbol.Marker = get(hand,'Marker');
	zz = get(hand,'MarkerSize');
	if (~iscell(zz)),   symbol.Size = num2cell(zz,1);
	else,				symbol.Size = zz;
	end
	zz = get(hand,'MarkerFaceColor');
	if (~iscell(zz)),	symbol.FillColor = num2cell(zz(:),1);
	else,				symbol.FillColor = zz;
	end
	zz = get(hand,'MarkerEdgeColor');
	if (~iscell(zz)),	symbol.EdgeColor = num2cell(zz(:),1);
	else,				symbol.EdgeColor = zz;
	end

	symbol.Marker = char(symbol.Marker);
	symbol.Marker = symbol.Marker(:,1);

	symbol.Marker(symbol.Marker == '^') = 't';
	symbol.Marker(symbol.Marker == '>') = 't';      % not in GMT
	symbol.Marker(symbol.Marker == '<') = 't';      % not in GMT
	symbol.Marker(symbol.Marker == 'v') = 'i';
	symbol.Marker(symbol.Marker == '.') = 'p';
	symbol.Marker(symbol.Marker == 'd') = 'd';
	symbol.Marker(symbol.Marker == 'o') = 'c';
	symbol.Marker(symbol.Marker == '+') = 'x';      % not in GMT
	symbol.Marker(symbol.Marker == 'x') = 'x';
	symbol.Marker(symbol.Marker == 's') = 's';
	symbol.Marker(symbol.Marker == '*') = 'a';
	symbol.Marker(symbol.Marker == 'p') = 'a';      % not in GMT
	symbol.Marker(symbol.Marker == 'h') = 'a';      % not in GMT

% ------------------------------------------------------------------------------------------------------------
function [LineStyle_num,LineStyle_gmt] = lineStyle2num(LineStyle)
	if (~iscell(LineStyle)),    LineStyle = {LineStyle};    end
	lt = {'-'; '--'; ':'; '-.'};
	LineStyle_num = strrep(LineStyle,lt{4},'4');
	LineStyle_num = strrep(LineStyle_num,lt{3},'3');
	LineStyle_num = strrep(LineStyle_num,lt{2},'2');
	LineStyle_num = strrep(LineStyle_num,lt{1},'1');
	tmp = LineStyle_num;
	LineStyle_num = str2num(cat(1,LineStyle_num{:}));
	% Convert to GMT linestyles
	tmp = strrep(tmp,'4',',.-');
	tmp = strrep(tmp,'3',',.');
	tmp = strrep(tmp,'2',',-');
	LineStyle_gmt = strrep(tmp,'1','');

% ------------------------------------------------------------------------------------------------------------
function script = write_group_symb(prefix,prefix_ddir,comm,pb,pf,ellips,symbols,n,script, opt_J)
% Write a group symbol to file, and uppdate the "script"
	l = numel(script) + 1;
	for (i = 1:numel(n))
		name = sprintf('%s_symb_%d.dat', prefix_ddir, i);
		name_sc = sprintf('%s_symb_%d.dat', prefix, i);
		fid = fopen(name,'wt');
		fc = symbols.FillColor{n(i)};		ec = symbols.EdgeColor{n(i)};
		if (ischar(fc)),	opt_G = '';
		else,				opt_G = sprintf(' -G%d/%d/%d', round(fc * 255));
		end
		if (ischar(ec))
			if (strcmp(ec, 'none')),	opt_W = '';
			else,						opt_W = ' -W1p';		% 'auto'. WRONG. Should be line's 'Color' property
			end
		else
			opt_W = sprintf(' -W1p,%d/%d/%d', round(ec * 255));
		end
		fprintf(fid,'%.5f\t%.5f\n',[symbols.x{n(i)}; symbols.y{n(i)}]);
		script{l} = ' ';					l=l+1;
		script{l} = [comm 'Plot symbols'];	l=l+1;
		script{l} = ['psxy ' name_sc ' -S' symbols.Marker(n(i)) num2str(symbols.Size{n(i)}) 'p' opt_G ...
                opt_W ellips opt_J ' -R -O -K >> ' pb 'ps' pf];    l=l+1;
		fclose(fid);
	end

% ------------------------------------------------------------------------------------------------------------
function [ALLpatchHand, hAlfaPatch] = findTransparents(ALLpatchHand)
% Find patches which have a level of transparency > 0.05
	ind = false(1,numel(ALLpatchHand));
	for (k = 1:numel(ALLpatchHand))
		fa = get(ALLpatchHand(k),'FaceAlpha');
		if (fa < 0.95 && fa > 0.05)			% Patch has transparency
			ind(k) = true;
		end
	end
	hAlfaPatch = ALLpatchHand(ind);			% Split the transparent and non-transparent
	ALLpatchHand(ind) = [];

% ------------------------------------------------------------------------------------------------------------
function set_opt_B(obj, evt, hPatch)
% Ask for a Frame Settings, an -B option for now, and store as an appdata in hPatch
	old_B = getappdata(hPatch, 'opt_B');
	opt_B = inputdlg({'Enter frame settings in the form of a GMT -B option'},'Frame settings',[1 60],{old_B});
	if (~isempty(opt_B)),	setappdata(hPatch, 'opt_B', opt_B{1}),	end

% ------------------------------------------------------------------------------------------------------------
function delete_fig(obj, evt, hPatch)
% Delete a patch fig and do the corresponding house cleaning
	handles = guidata(hPatch);
	N = get(handles.popup_familyPlots, 'Val');
	str = get(handles.popup_familyPlots, 'Str');
	if (numel(str) == 1)
		delete(handles.figure1)		% BYE BYE
		return
		%set(handles.popup_familyPlots, 'Vis', 'off');	% Since we have no more --- NOT YET READY
	end
	str(N) = [];
	set(handles.popup_familyPlots, 'Val', 1, 'Str', str)
	grid_figs(handles, numel(str))

	% Update the N stored in appdatas of all hRects that have a number higher than the one being killed
	for (k = N+1:numel(str)+1)
		setappdata(handles.hRect(k), 'My_N', k-1)
	end

	handles.hAllFigs(N) = [];
	delete(handles.hRect(N));	handles.hRect(N) = [];

	handles.x_min(N) = [];		handles.x_max(N) = [];	handles.y_min(N) = [];	handles.y_max(N) = [];
	handles.X0(N) = [];			handles.Y0(N) = [];		handles.scale{N} = [];	handles.opt_R(N) = [];
	handles.width_orig(N)  = [];		handles.height_orig(N) = [];	handles.map_width(N) = [];
	handles.projection_str(N) = [];
	delete(handles.hAllImgs(N)),		handles.hAllImgs(N) = [];	handles.orig_img{N} = [];
	guidata(handles.figure1, handles)

% ------------------------------------------------------------------------------------------------------------
function draw_rectangle(hObj, evt)
	handles = guidata(hObj);
	[p1,p2,hLine] = rubberbandbox(handles.axes1);
	difa = abs(p2 - p1);
	if ((difa(1) < handles.head(7)/4) || (difa(2) < handles.head(8)/4))
		delete(hLine),		return			% Don't draw ultra small rectangles
	end
	set(hLine,'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick)
	draw_funs([],'set_rect_uictx_PC', hLine)	% Set lines's uicontextmenu

% ------------------------------------------------------------------------------------------------------------
function figure1_ResizeFcn(hObj, evt)
% Move the right column of uicontrols such that they stay at the same absolute distance from Fig right side
	handles = guidata(hObj);
	posFig = get(handles.figure1, 'Pos');
	if (isequal(posFig(3:4), [800 653]))		% Used to avoid this function to run when on fig creation
		return
	end
	allHands = [getappdata(handles.figure1, 'allHands'); findobj(handles.figure1,'Style','Frame')];
	for (k = 1:numel(allHands))			% Movall all uicontrols on right column
		ud = get(allHands(k), 'UserData');		x1 = posFig(3) - ud;
		pos = get(allHands(k), 'Pos');			pos(1) = x1;
		set(allHands(k), 'Pos', pos)
	end
	sqrt2 = 1.414213562373095;		% sqrt(2)
	posAx = get(handles.axes1, 'Pos');
	posR  = get(handles.popup_familyPlots, 'Pos');
	if (get(handles.radio_Portrait, 'Val'))
		H = min(posFig(4) - posAx(2) - 4, (posR(1) - 30) * sqrt2);
		W = H / sqrt2;
		posAx(1) = 25 + (posR(1)-30)/2 - posAx(3)/2;		posAx(2) = 100;		
	else
		maxW = (posFig(4) -100 -4) * sqrt2;
		W = min(posR(1) - 30, maxW);
		H = W / sqrt2;
		posAx(1) = 25;
		posAx(2) = (posFig(4) - posAx(4))/2 + 48;	% = 100 + (posFig(4) - 100 - 4)/2 - posAx(4)/2;
	end
	posAx(3) = W;		posAx(4) = H;
	set(handles.axes1, 'Pos', posAx);
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------------
function plot_composer_LayoutFcn(hFig)
figW = 800;
h = zeros(31,1);
set(hFig, 'Position',[520 131 800 653],...
	'PaperUnits','centimeters',...
	'ResizeFcn',@figure1_ResizeFcn,...
	'Color',get(0,'factoryUicontrolBackgroundColor'),...
	'MenuBar','none',...
	'Name','Plot Composer',...
	'NumberTitle','off',...
	'DoubleBuffer','on',...
	'PaperPosition',[0.6345175 6.345175 20.30456 15.22842],...
	'PaperSize',[21.573595 27.91877],...
	'RendererMode','manual',...
	'Resize','on',...
	'HandleVisibility','callback',...
	'Tag','figure1');

axes('Parent',hFig, 'Units','pixels', 'Position',[99 100 382 543],...
	'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
	'Tag','axes1');

uicontrol('Parent',hFig, 'Position',[570 499 222 91],  'Style','frame', 'UserData', figW-570);
uicontrol('Parent',hFig, 'Position',[570 349 222 131], 'Style','frame', 'UserData', figW-570);
uicontrol('Parent',hFig, 'Position',[570 199 222 131], 'Style','frame', 'UserData', figW-570);

h(1) = uicontrol('Parent',hFig, 'Position',[576 529 71 15],...
	'Callback',@plot_composer_uiCB,...
	'String','Portrait',...
	'Style','radiobutton',...
	'Value',1,...
	'UserData', figW-576, ...
	'Tag','radio_Portrait');

h(2) = uicontrol('Parent',hFig, 'Position',[570 608 170 22],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'String','',...
	'Style','popupmenu',...
	'TooltipString','Select the Mirone window to plot',...
	'Value',1,...
	'UserData', figW-570, ...
	'Tag','popup_familyPlots');

h(3) = uicontrol('Parent',hFig, 'Position',[738 608 54 22],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'String','',...
	'Style','popupmenu',...
	'TooltipString','Choose how to arrange the figures in a grid',...
	'Visible','off',...
	'UserData', figW-738, ...
	'Tag','popup_gridFigs');

h(4) = uicontrol('Parent',hFig, 'Position',[576 553 161 22],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'String','A4 595 842',...
	'Style','popupmenu',...
	'Value',1,...
	'UserData', figW-576, ...
	'Tag','popup_paperSize');

h(5) = uicontrol('Parent',hFig, 'Position',[735 553 50 22],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'String',{'cm'; 'in'; 'pt'},...
	'Style','popupmenu',...
	'Value',1,...
	'UserData', figW-735, ...
	'Tag','popup_paperUnits');

h(6) = uicontrol('Parent',hFig, 'Position',[658 529 75 15],...
	'Callback',@plot_composer_uiCB,...
	'String','Landscape',...
	'Style','radiobutton',...
	'UserData', figW-658, ...
	'Tag','radio_Landscape');

h(7) = uicontrol('Parent',hFig, 'Position',[655 207 61 21],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'String','2.5',...
	'Style','edit',...
	'TooltipString','Plot X origin',...
	'UserData', figW-655, ...
	'Tag','edit_X0');

h(8) = uicontrol('Parent',hFig, 'Position',[623 507 120 15],...
	'String','Trim white borders',...
	'Style','checkbox',...
	'TooltipString','Remove all border spaces arround the figure.',...
	'Value',1,...
	'UserData', figW-623, ...
	'Tag','check_trimWhite');

h(9) = uicontrol('Parent',hFig, 'Position',[576 448 65 15],...
	'Callback',@plot_composer_uiCB,...
	'String','Auto run',...
	'Style','radiobutton',...
	'TooltipString','Run GMT commands in the backgroud to recreate this figure.',...
	'Value',1,...
	'UserData', figW-576, ...
	'Tag','radio_autoRun');

h(10) = uicontrol('Parent',hFig, 'Position',[644 445 50 21],...
	'BackgroundColor',[1 1 1],...
	'String',{'pdf'; 'jpg'; 'png'; 'ps'},...
	'Style','popupmenu',...
	'Tooltip','Select format for the auto run figure',...
	'Value',1,...
	'UserData', figW-644, ...
	'Tag','popup_figFormat');

h(11) = uicontrol('Parent',hFig, 'Position',[576 417 105 15],...
	'Callback',@plot_composer_uiCB,...
	'String','Write GMT script',...
	'Style','radiobutton',...
	'Tooltip','Save a GMT scrit that will reproduce this plot',...
	'UserData', figW-576, ...
	'Tag','radio_writeScript');

h(12) = uicontrol('Parent',hFig,'Position',[576 388 191 22],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'String','Nikles',...
	'Style','popupmenu',...
	'TooltipString','Save script and files in this directory',...
	'Value',1,...
	'UserData', figW-576, ...
	'Tag','popup_directory_list');

h(13) = uicontrol('Parent',hFig, 'Position',[765 390 21 21],...
	'Callback',@plot_composer_uiCB,...
	'FontWeight','bold',...
	'String','...',...
	'UserData', figW-765, ...
	'Tag','push_change_dir');

h(14) = uicontrol('Parent',hFig, 'Position',[686 360 100 21],...
	'BackgroundColor',[1 1 1],...
	'String','',...
	'Style','edit',...
	'TooltipString','Script and files name prefix',...
	'UserData', figW-686, ...
	'Tag','edit_prefix');

h(15) = uicontrol('Parent',hFig, 'Position',[653 292 131 21],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'String','',...
	'Style','edit',...
	'TooltipString','Aproximate map scale',...
	'UserData', figW-653, ...
	'Tag','edit_scale');

h(16) = uicontrol('Parent',hFig, 'Position',[655 248 60 21],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'String','',...
	'Style','edit',...
	'TooltipString','Map width',...
	'UserData', figW-655, ...
	'Tag','edit_mapWidth');

h(17) = uicontrol('Parent',hFig, 'Position',[725 248 60 21],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'String','',...
	'Style','edit',...
	'TooltipString','Map height',...
	'UserData', figW-725, ...
	'Tag','edit_mapHeight');

h(18) = uicontrol('Parent',hFig, 'Position',[725 207 61 21],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'String','2.5',...
	'Style','edit',...
	'Tooltip','Plot Y origin',...
	'UserData', figW-725, ...
	'Tag','edit_Y0');

h(19) = uicontrol('Parent',hFig, 'Position',[600 252 51 15],...
	'HorizontalAlignment','left',...
	'String','Map Dims',...
	'Style','text',...
	'UserData', figW-600, ...
	'Tag','text1');

h(20) = uicontrol('Parent',hFig, 'Position',[671 271 35 16], 'String','Width', 'Style','text',...
	'UserData', figW-671, ...
	'Tag','text2');

h(21) = uicontrol('Parent',hFig, 'Position',[665 228 41 15],...
	'HorizontalAlignment','left',...
	'String','X origin',...
	'Style','text',...
	'UserData', figW-665, ...
	'Tag','text3');

h(22) = uicontrol('Parent',hFig, 'Position',[735 228 41 15],...
	'HorizontalAlignment','left',...
	'String','Y origin',...
	'Style','text',...
	'UserData', figW-735, ...
	'Tag','text4');

h(23) = uicontrol('Parent',hFig, 'Position',[625 363 58 15],...
	'HorizontalAlignment','left',...
	'String','Name prefix',...
	'Style','text',...
	'UserData', figW-625, ...
	'Tag','text5');

h(24) = uicontrol('Parent',hFig, 'Position',[575 295 76 15],...
	'HorizontalAlignment','left',...
	'String','Map scale (apr)',...
	'Style','text',...
	'UserData', figW-575, ...
	'Tag','text6');

uicontrol('Parent',hFig, 'Position',[151 45 331 22],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'String',' ',...
	'Style','popupmenu',...
	'Value',1,...
	'Tag','popup_projections');

uicontrol('Parent',hFig, 'Position',[491 46 120 23],...
	'Callback',@plot_composer_uiCB,...
	'FontSize',8,...
	'FontWeight','normal',...
	'String','CRS of the world',...
	'Tooltip','Select a Coordinate Reference System from a World list',...
	'Tag','push_worldCRS');

uicontrol('Parent',hFig, 'Position',[10 47 140 16],...
	'FontName','Helvetica',...
	'HorizontalAlignment','center',...
	'String','Recent Coord Ref Systems',...
	'Style','text',...
	'Tag','text8');

uicontrol('Parent',hFig, 'Position',[5 5 606 40],...
	'BackgroundColor',[1 1 1],...
	'Callback',@plot_composer_uiCB,...
	'HorizontalAlignment','left',...
	'Max',3,...
	'String','',...
	'Style','edit',...
	'Tooltip','This is the Proj4 definition string that describes the plot coordinate system. Blank, defaults to geogs',...
	'Tag','edit_projection');

h(25) = uicontrol('Parent',hFig, 'Position',[734 271 40 16], 'String','Height', 'Style','text',...
	'UserData', figW-734, ...
	'Tag','text9');

h(26) = uicontrol('Parent',hFig, 'Position',[600 210 51 15],...
	'HorizontalAlignment','left',...
	'String','Map origin',...
	'Style','text',...
	'UserData', figW-600, ...
	'Tag','text10');

h(27) = uicontrol('Parent',hFig, 'Position',[636 582 70 15],...
	'FontSize',9,...
	'String','Paper size',...
	'Style','text',...
	'UserData', figW-636, ...
	'Tag','text_ps');

h(28) = uicontrol('Parent',hFig, 'Position',[642 472 50 15],...
	'FontSize',9,...
	'String','Output',...
	'Style','text',...
	'UserData', figW-642, ...
	'Tag','text_out');

uicontrol('Parent',hFig, 'Position',[690 127 70 15],...
	'String','Scale Bar',...
	'Style','checkbox',...
	'Visible', 'off', ...
	'Tag','check_scaleBar');

h(29) = uicontrol('Parent',hFig,'Position',[691 105 80 15],...
	'String','Time stamp',...
	'Style','checkbox',...
	'Tooltip', 'Plot a GMT time stamp in lower left corner', ...
	'UserData', figW-691, ...
	'Tag','check_timeStamp');

h(30) = uicontrol('Parent',hFig, 'Position',[639 322 80 15],...
	'FontSize',9,...
	'String','Dimensions',...
	'Style','text',...
	'UserData', figW-639, ...
	'Tag','text_dim');

h(31) = uicontrol('Parent',hFig, 'Position',[680 20 100 24],...
	'Callback',@plot_composer_uiCB,...
	'FontSize',9,...
	'FontWeight','bold',...
	'String','OK',...
	'UserData', figW-680, ...
	'Tag','push_OK');

setappdata(hFig, 'allHands', h)

function plot_composer_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

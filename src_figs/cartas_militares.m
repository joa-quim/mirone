function varargout = cartas_militares(varargin)
% Load GIF files with the Portuguese "Cartas Militares" and display them georeferenced in Mirone
%
% Georeferencing is acomplished in one of 3 ways (by that order):
% 1. A corresponding Ozi Explorer .map file is found on the same directory as the GIF file
% 2. A .gfw or .wrl (with contents like in a .tfw file) file is found on the same directory as the GIF file
% 3. Internal georeferencing based on the known coordinates (in "Coordenadas Militares") of each image
% 
% Thirth method should be prefered (more accurate, I believe)

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

	hObject = figure('Visible','off');
	cartas_militares_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	if (numel(varargin) > 0)
		handMir = varargin{1};
        path_data = handMir.path_data;
	else
	    path_data = [cd filesep 'data' filesep];
	end

	% Load the directory list stored in mirone_pref
	load([path_data 'mirone_pref.mat'], 'directory_list');
	set(handles.popup_directory_list,'String', directory_list)
	handles.last_directories = directory_list;
	handles.files_dir = handles.last_directories{1};

	% ===== UGLY PATCH TO IMPLEMENT THE LIDAR PT MOOSAIC TOOL ====
	if (numel(varargin) == 2)
		handles.path_data = path_data;
		lidarPT(handles);
		set(hObject,'Vis','on')
		guidata(hObject, handles);
		return
	end

    naos = cell(19,1);				% By Columns
    naos{1} = [1:19 22:58];
    naos{2} = [1:19 28:58];
    naos{3} = [1:16 30:58];
    naos{4} = [4:16 35:58];
    naos{5} = [1 43:47 57:58];
    naos{6} = [1 58];
    naos{7} = [1 58];
    naos{8} = 1;
    naos{9} = [56 58];
    naos{10} = 56:58;
    naos{11} = 56:58;
    naos{12} = [1 2 56:58];
    naos{13} = [1:9 29:30 56:58];
    naos{14} = [1:11 16:20 26:30 57:58];
    naos{15} = [1:12 15:22 25:31 36 57:58];
    naos{16} = [1:46 57:58];
    naos{17} = [1:47 57:58];
    naos{18} = [1:49 54:58];
    naos{19} = 1:58;

    atlas_file = [path_data 'countries_dp5.bin'];
    paises.ct = country_select(atlas_file,'-Pportugal');
	% Clean up the empty fields in the ct struct (given I could not do it at mex level)
	id = false(numel(paises.ct),1);
	for (k = 1:numel(paises.ct))
        if (isempty(paises.ct(k).Country)),     id(k) = true;    end
	end
	paises.ct(id) = [];

	xMap_min = 72000;		yMap_min = 0;

	set(handles.axes1,'XLim',[0 (19*16000)] + xMap_min + [-1000 1000])
	set(handles.axes1,'YLim',[0 (58*10000)] + yMap_min + [-1000 1000])
	set(handles.axes1,'XTick',[],'YTick',[], 'DataAspectRatio', [1 1 1])

	projStruc.DstProjWKT = ogrproj('+proj=tmerc +lat_0=39.668258333333333 +lon_0=-8.133108611111111 +k=1.0 +x_0=200000 +y_0=300000 +ellps=intl');
	xy_prj = ogrproj([paises.ct.Country(1,:)' paises.ct.Country(2,:)']+0, projStruc);
	line('XData',xy_prj(:,1),'YData',xy_prj(:,2),'Parent',handles.axes1);

	set(hObject,'Vis','on')
	handles.mapa = cell(58,18);
	nr = 1;                     % counter on the number of rectangles
	for (m=58:-1:1)
		yp = [(m-1) (m-1) m m] * 10000 + yMap_min;
		for (n=1:19)            % Loop over columns
 			if (find(naos{n} == m)),    continue,	end
			xp = [(n-1) n n (n-1)] * 16000 + xMap_min; 
			h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',[.7 .7 .7],'LineWidth',1,'Parent',handles.axes1);
			set(h,'ButtonDownFcn',{@bdnTile,hObject},'UserData',[m n])
			handles.mapa{m,n} = sprintf('%.3d',nr);
			text(xp(1)+8000,yp(4)-5000,sprintf('%d',nr),'HorizontalAlignment','center', 'FontSize',5,'HitTest','off','Parent',handles.axes1)
			nr = nr + 1;
		end
		if (~rem(m,5)),			pause(0.01),	end
	end

	% Exceptions to the pure numeric name tile. For example '401A' ---------------------------------------------
    ABs = cell(19,2);
	ABs{1,1} = 22;			ABs{1,2} = {'401A'};
	ABs{2,1} = 19;			ABs{2,2} = {'441B'};
	ABs{3,1} = 30;			ABs{3,2} = {'306B'};
	ABs{4,1} = [4 12 35:37];		ABs{4,2} = {'583A' '515A' '248B' '238A' '227B'};
	ABs{5,1} = 43;			ABs{5,2} = {'162A'};
	ABs{13,1}= [29 30 56];	ABs{13,2}= {'325A' '315A' '9A'};
	ABs{14,1}= [16 20 30];	ABs{14,2}= {'483A' '441A' '315B'};
	ABs{15,1}= [12 31 36];	ABs{15,2}= {'525A' '306A' '248A'};
	ABs{16,1}= [38 40 46];	ABs{16,2}= {'227A' '205A' '142A'};
	ABs{18,1}= 49;			ABs{18,2}= {'108A'};
	ABs{19,1}= 52;			ABs{19,2}= {'67A'};

	for (n=1:19)
		xp = [(n-1) n n (n-1)] * 16000 + xMap_min; 
		for (k = 1:numel(ABs{n,1}))
			m = ABs{n,1}(k);
			yp = [(m-1) (m-1) m m] * 10000 + yMap_min;
			h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',[.7 .7 .7],'LineWidth',1,'Parent',handles.axes1);
			set(h,'ButtonDownFcn',{@bdnTile,hObject},'UserData',{ABs{n,2}{k},[m n]})	% Uf, complicated this UD
			text(xp(1)+8000,yp(4)-5000, ABs{n,2}{k}, 'HorizontalAlignment','center', 'FontSize',5,'HitTest','off','Parent',handles.axes1)			
		end
	end
	% -------------------------------------------------------------------------------------------------------------

	% Exceptions TYPE II ------------------------------------------------------------------------------------------
	xp = [76000 83659 83659 76000];		yp = [270000 270000 280125 280125];
	h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',[.7 .7 .7],'LineWidth',1,'Parent',handles.axes1);
	tile.name = '325B';			tile.coords_x = xp;		tile.coords_y = yp;
	set(h,'ButtonDownFcn',{@bdnTile,hObject},'UserData',tile)		% Use a structure in UD to destinguish from the other 2 types
	text(xp(1)+4000,yp(4)-5000, '325B', 'HorizontalAlignment','center', 'FontSize',5,'HitTest','off','Parent',handles.axes1)			

	xp = [96344 104000 104000 96344];		yp = [270000 270000 280000 280000];
	h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',[.7 .7 .7],'LineWidth',1,'Parent',handles.axes1);
	tile.name = '325C';			tile.coords_x = xp;		tile.coords_y = yp;
	set(h,'ButtonDownFcn',{@bdnTile,hObject},'UserData',tile)		% Use a structure in UD to destinguish from the other 2 types

	xp = [9 10 10 9] * 16000 + xMap_min;		yp = [-617 -617 10000 10000];
	h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',[.7 .7 .7],'LineWidth',1,'Parent',handles.axes1);
	tile.name = '611';			tile.coords_x = xp;		tile.coords_y = yp;
	set(h,'ButtonDownFcn',{@bdnTile,hObject},'UserData',tile)		% Use a structure in UD to destinguish from the other 2 types
	% -------------------------------------------------------------------------------------------------------------
	
	% Tooltip
	msg = sprintf(['When you want to get the files directly from Web\n\n' ...
		'As for example at:\nhttp://www.civil.ist.utl.pt/~ruif/HRF/cartas/Cartas%%20Militares%%20GIF_25000/']);
	set(handles.radio_inWeb,'Tooltip',msg)

	% Initiliaze with this site. This also means that it doesn't recall changes between usage sessions
	set(handles.edit_forWeb,'String', 'http://www.civil.ist.utl.pt/~ruif/HRF/cartas/Cartas%20Militares%20GIF_25000/')

	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end

% -----------------------------------------------------------------------------------------
function lidarPT(handles)
% This function ends the configuration work nomaly done in main and load/display
% a background image plus the lidar tiles matrix, created as a single patch object.

	set(handles.radio_inWeb,'Enable','off')
	set(handles.push_lidarMosaico,'Vis','on')

	load([handles.path_data 'lidarPT.mat'])
	x_min = double(mosaico(1,1));	x_max = double(mosaico(1,2));
	y_min = double(mosaico(1,3));	y_max = double(mosaico(1,4));
	x_inc = 1600;		y_inc = 1000;

	[I, att] = gdalread([handles.path_data 'PTimg_lidar.jp2']);
	X = att.GMT_hdr(1:2);		Y = att.GMT_hdr(3:4);
	image(X,Y,flipdim(I,1), 'Parent',handles.axes1)
	x_lim = [x_min x_max] + [-35 30]*x_inc;
	y_lim = [y_min y_max] + [-2 2]*y_inc;
	set(handles.axes1,'XLim', x_lim, 'YLim', y_lim, 'YDir','normal')
	set(handles.axes1,'XTick',[],'YTick',[], 'DataAspectRatio', [1 1 1])
	addHelpLegend(handles)

	nColsQ = (x_max - x_min) / x_inc + 1;		% Number of cells along X in our PT
	X = x_min:x_inc:x_max;
	Y = y_min:y_inc:y_max;
	nCols = numel(X);
	nRows = numel(Y);
	sz =  nCols*nRows;

	n_quadric = size(mosaico,1)-1;				% Number of actual sub-grids in the LIDAR2011 PT
	ind(n_quadric) = int32(0);
	cell_yes = false(nCols * (nRows-1),1);	% ...
	for (k = 2:n_quadric+1)
		col = round((double(mosaico(k,1)) - x_min) / x_inc) + 1;
		row = round((double(mosaico(k,3)) - y_min) / y_inc ) + 1;
		c = int32((row - 1) * nColsQ + col);
		cell_yes(c) = true;
		ind(k-1) = c;
	end

	[Y, X] = meshgrid (Y, X);
	X = reshape(X, sz, 1); 
	Y = reshape(Y, sz, 1);
	vert = [X Y];
	clear X Y
	faces = zeros(nCols*(nRows-1),4);
	j = 1;
	for (i = 1:nCols * (nRows-1))
		if (rem(i,nCols) ~= 0 && cell_yes(i))
			faces(j,:) = [i, i+1, i+1+nCols, i+nCols];             
			j = j + 1;
		end
	end
	faces(j:end,:) = [];

	FV.vertices = vert;		FV.faces = faces;
	h = patch(FV,'FaceColor','none','Parent',handles.axes1, 'EdgeColor','y', 'Tag','LidarMatrix');
	set(h,'ButtonDownFcn',{@bdnLidar,handles.figure1},'UserData',[x_lim y_lim x_inc y_inc])
	setappdata(h, 'infoBB', [x_min x_max y_min y_max nColsQ])
	setappdata(h, 'infoIND', ind)
	setappdata(h, 'infoNOME', nomes)
	sldT = 9;		% Slider thickness
	H = 50;			% Distance from bottom of Fig
	setSliders(handles, sldT, H);

% -----------------------------------------------------------------------------------------
function addHelpLegend(handles)
% Create another axes with the help legend

	figPos = get(handles.figure1, 'pos');
	pos = [figPos(3)/2-10 figPos(4)-198 figPos(3)/2-8 190];
	hAx = axes('Parent',handles.figure1,...
			'Units','pixels', 'Pos',pos, 'Tag','axes2','XTick',[],'YTick',[],'HandleVisibility','off');
	txt = sprintf(['Seleccionar o dir dos dados\n\n',...
		'Para zoom, use as teclas +/-\nUse as setas para deslocar\nhorizontal e verticalmente\n\n' ...
		'Clique num quadradinho para\nseleccionar. Mais que um\nselecciona um rectângulo\n\n' ...
		'No fim clica "Faz mosaico"']);
	text(0.02,0.98,txt,'Parent', hAx, 'VerticalAlignment','Top')

% -----------------------------------------------------------------------------------------
function bdnLidar(obj,evt,hFig)
% Create a red patch in the cell of tha main data matrix clicked by user

	handles = guidata(hFig);
	hQ = gcbo;
	pt = get(handles.axes1, 'CurrentPoint');		% Current click in map coordinates
	ud = get(hQ,'UserData');						% = [x_lim y_lim x_inc y_inc]
	col = fix((pt(1,1) - ud(1)) / ud(5));			% col number starting at 0 from axes (image) origin
	row = fix((pt(1,2) - ud(3)) / ud(6));			% rown ----"----
	x1 = ud(1) + col * ud(5);	x2 = x1 + ud(5);	% coordinates of the clicked rectangle
	y1 = ud(3) + row * ud(6);	y2 = y1 + ud(6);
	xp = [x1 x1 x2 x2];			yp = [y1 y2 y2 y1];
	h = patch('XData',xp,'YData',yp,'FaceColor','r','Parent',handles.axes1,'Tag','oneTile');

	[info.nome, info.row_col] = findInnerLidarTiles(hQ, x1, y1);
	set(h, 'UserData', info)
	set(h,'ButtonDownFcn','delete(gco)')			% Kill me when I'm clicked

% -----------------------------------------------------------------------------------------
function [nome, row_col] = findInnerLidarTiles(hQ, x, y)
% Get the grid name associated to the cell where X,Y falls and its ROW,COL address in the data matrix

	nome = [];	row_col = [];
	ud = get(hQ,'UserData');						% = [x_lim y_lim x_inc y_inc]
	infoBB = getappdata(hQ, 'infoBB');				% infoBB holds the true data limits, [x_min x_max y_min y_max nColsQ]
	col = round((x - infoBB(1)) / ud(5)) + 1;		% col number starting at 1 from data origin
	row = round((y - infoBB(3)) / ud(6)) + 1;		% rown ----"----
	ind = int32((row - 1) * infoBB(5) + col);		% row major linear index of clicked cell (one based)
	infoIND = getappdata(hQ, 'infoIND');
	infoNOME = getappdata(hQ, 'infoNOME');
	k = find(infoIND == ind);						% Find the index of the 'ind' linear index
	if (isempty(k)),	return,		end
	nome = infoNOME{k+1};							% + 1 because first name is empty (the global BB)
	row_col = [row col]; 

% -----------------------------------------------------------------------------------------
function push_lidarMosaico_CB(hObject, handles)
% Check which tiles are selected and build a mosaic that embraces them all (even the inner non-selected ones)

	hQ = findobj(handles.axes1,'Type','patch','Tag','LidarMatrix');
	t = findobj(handles.axes1,'Type','patch','Tag','oneTile');
	if (isempty(t)),	return,		end
	row_min = 10000;	row_max = 0;
	col_min = 10000;	col_max = 0;
	z_min   = 10000;	z_max   = -1000;
	for (k = 1:numel(t))						% Find our corners
		ud = get(t(k), 'UserData');
		row_min = min(row_min, ud.row_col(1));
		row_max = max(row_max, ud.row_col(1));
		col_min = min(col_min, ud.row_col(2));
		col_max = max(col_max, ud.row_col(2));
	end

	infoBB = getappdata(hQ, 'infoBB');			% infoBB holds the true data limits, [x_min x_max y_min y_max nColsQ]

	pato = get(handles.popup_directory_list,'String');
	pato = pato{1};
	if (pato(end) ~= filesep),		pato = [pato filesep];		end
	%pato = 'C:\dat\LiDAR2011\laz\';

	% Start filling our mosaic
	nx = (col_max - col_min + 1) * 1600 / 2 + 1;% 2 is the grid step and + 1 because of grid registration
	ny = (row_max - row_min + 1) * 1000 / 2 + 1;
	Z = alloc_mex(ny,nx,'single',nan);
	X0 = infoBB(1) + (col_min - 1) * 1600;
	Y0 = infoBB(3) + (row_min - 1) * 1000;
	got_one = false;
	for (m = 1:(row_max - row_min + 1))
		y = Y0 + (m - 1) * 1000 + 100;			% The 100 is just put the point well inside this cell
		for (n = 1:(col_max - col_min + 1))
			x = X0 + (n - 1) * 1600 + 100;
			nome = findInnerLidarTiles(hQ, x, y);
			if (isempty(nome)),		continue,	end
			fname = [pato nome '-mis_orto.laz'];
			if (exist(fname,'file') ~= 2)
				fname = [pato nome '-mis_orto.xyz'];		% Try the XYZ version
				if (exist(fname,'file') ~= 2)
					fname = [pato nome '-mis_orto.asc'];	% Try the ASC version
				end
			end
			if (exist(fname,'file') ~= 2),	continue,	end	% Should give a warning here
			z = readLidarTile(fname);
			got_one = true;
			zzz = grdutils(z,'-L');		z_min = min(zzz(1), z_min);		z_max = max(zzz(2), z_max);
			Z((m-1)*500+1:m*500+1, (n-1)*800+1:n*800+1) = z;
		end
	end

	if (~got_one)
		warndlg('Nickles de LIDAR grelhas no directório indicado. Bye Bye','Warning'),	return
	end
	X1 = infoBB(1) + col_max * 1600;
	Y1 = infoBB(3) + row_max * 1000;
	meta.X = linspace(X0, X1, nx);		meta.Y = linspace(Y0, Y1, ny);
	meta.name = 'Mosaico LIDAR2011';	meta.geog = 0;
	meta.head = [X0 X1 Y0 Y1 z_min z_max 0 2 2];
	meta.srsWKT = ['PROJCS["ETRS89_Portugal_TM06",GEOGCS["GCS_ETRS89",DATUM["D_ETRS_1989",' ...
			'SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],' ...
			'PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",39.66825833333333],' ...
			'PARAMETER["central_meridian",-8.133108333333334],PARAMETER["scale_factor",1],'...
			'PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'];
	mirone(Z, meta)

% -----------------------------------------------------------------------------------------
function z = readLidarTile(fname)
% Read a LIDAR_PT tile. Check if FNAME is a LASZ file

	[PATO, NAME, EXT] = fileparts(fname);
	if (strcmpi(EXT,'.laz'))
		xyz = laszreader_mex(fname);
		z = single(xyz(3,:));
		clear xyz
		z = reshape(z,801,501)';
	else
		z = gdalread(fname);
	end
	z(z == -999) = NaN;
	z = flipud(z);

% -----------------------------------------------------------------------------------------
function bdnTile(obj,event,hFig)
% Do what ever it must to get a referenced 1:25000 image
	handles = guidata(hFig);
	ud = get(gcbo,'UserData');
	proj = [];		knowLimits = false;

	if (isa(ud,'cell'))				% Exception TYPE I (tile names with inbeded letters)
		tile_name = ud{1};
		m = ud{2}(1);		n = ud{2}(2);		% We still need these for computing the corner coordinates
	elseif (isa(ud,'struct'))		% Exception TYPE II (size of tiles not 10x16 km)
		tile_name = ud.name;
		tmp.X = ud.coords_x(1:2);	tmp.Y = ud.coords_y(2:3);
		knowLimits = true;
	else							% Most of the cases
		m = ud(1);			n = ud(2);
		tile_name = handles.mapa{m,n};
	end

	if (get(handles.radio_inloco, 'Value'))
		pato = get(handles.popup_directory_list,'String');
		pato = pato{1};
		if (pato(end) ~= filesep),		pato = [pato filesep];		end
		fname = [pato tile_name '.gif'];
		if (exist(fname,'file') ~= 2)
			fname = [pato tile_name '.sid'];		% Try in Mr. Sid format (MUST ADD A TEST IF MRSID DRIVER EXISTS)
		end
		if (exist(fname,'file') ~= 2)
			h = text(100000,50000, ['FILE  ' tile_name '.GIF  NOT FOUND IN THERE'], 'HorizontalAlignment','left', 'FontSize',25, 'Rotation', 65, 'Color','r');
			pause(2),	delete(h)
			return
		end
		fnameRef = [pato tile_name '.map'];
	else
		url = get(handles.edit_forWeb,'String');
		if (url(end) ~= '/'),			url = [url '/'];		end
		if (~strncmp(url,'http://',7)),	url = ['http://' url];	end
		fname = [url tile_name '.GIF'];
		fnameRef = [];			% We need to find a way to fetch the corresponding .map file
	end
		
	set(handles.figure1,'pointer','watch'),		pause(0.1)
	[img,att] = gdalread(fname, '-U');
	if (isempty(img))
		errordlg('There was an error reading the image. We got nothing.','Error')
		return
	end

	if (~isempty(fnameRef) && exist(fnameRef,'file') == 2)			% Found a Ozi .map reference file
		[gcp, proj] = getProjFromOzi(fnameRef);
		% Compute and apply the affine transformation
		trans = AffineTransform(gcp(:,1:2),gcp(:,3:4));
		x_pt = [1; size(img,2)];    y_pt = [1; size(img,1)];
		X1 = [x_pt y_pt ones(size(x_pt,1),1)];
		U1 = X1 * trans;
		tmp.X = U1(:,1)';		tmp.Y = U1(:,2)';
		if (tmp.Y(1) > tmp.Y(2))			% They just love that stupid origin in the UL corner
			yy = tmp.Y(1);		tmp.Y(1) = tmp.Y(2);	tmp.Y(2) = yy;
		end
	elseif (~isempty(fnameRef) && ~isempty(att.GeoTransform))		% We had a world file
		tmp.X = [att.Corners.LL(1) att.Corners.LR(1)];
		tmp.Y = [att.Corners.LL(2) att.Corners.UL(2)];
	else
		if (~knowLimits)		% Otherwise we already know them
			tmp.X = [(n-1) n] * 16000 + 72000;
			tmp.Y = [(m-1) m] * 10000;
		end
		% Set the proj string for the "Coordenadas militares, datum Lisboa"
		proj = '+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111111 +k=1.0 +x_0=200000 +y_0=300000 +ellps=intl +towgs84=-304.046,-60.576,103.640,0,0,0,0';
	end

	tmp.geog = 0;       tmp.name = fname;
	x_inc = diff(tmp.X) / (size(img,2) - 1);
	y_inc = diff(tmp.Y) / (size(img,1) - 1);
	tmp.head = [tmp.X tmp.Y att.GMT_hdr(5:6) 0 x_inc y_inc];
	
	try
		if (ndims(img) ~= 3)	% The Sided version are rgb
			tmp.cmap = att.Band.ColorMap.CMap(:,1:3);
		end
	catch
		errordlg(['Bad colormap in file ' fname],'ERROR')
		set(handles.figure1,'pointer','arrow')
		return
	end
	if (~isempty(proj))			% If we know the projection, convert it to the WKT form
		tmp.srsWKT = ogrproj(proj);
	end
	mirone(img, tmp);
	set(handles.figure1,'pointer','arrow')

% ----------------------------------------------------------------------------
function popup_directory_list_CB(hObject, handles, opt)
% OPT is used by push_change_dir (just to save code)

	if (nargin == 2),	opt = [];   end
	if isempty(opt)
		val = get(hObject,'Value');     str = get(hObject, 'String');
		% Put the selected field on top of the String list. This is necessary because the "OK" button will
		tmp = str(val);         str(val) = [];
		new_str = [tmp; str];   set(hObject,'String',new_str); 
		set(hObject,'Value',1)
		if iscell(tmp),			new_dir = tmp{1};
		elseif ischar(tmp),		new_dir = tmp;
		else                    return        % ???
		end
	else
		new_dir = opt;
	end
	handles.files_dir = new_dir;
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function push_change_dir_CB(hObject, handles)
	if (strcmp(computer, 'PCWIN'))
		work_dir = uigetfolder_win32('Select a directory', cd);
	else            % This guy doesn't let to be compiled
		work_dir = uigetdir(cd, 'Select a directory');
	end
	if (isempty(work_dir)),		return,		end
	handles.last_directories = [cellstr(work_dir); handles.last_directories];
	set(handles.popup_directory_list,'String',handles.last_directories)
	guidata(hObject, handles);
	popup_directory_list_CB(handles.popup_directory_list, handles, work_dir)

% ----------------------------------------------------------------------------
function radio_inloco_CB(hObject, handles)
	if (~get(hObject,'Val')),		set(hObject,'Val',1),	return,		end
	set(handles.radio_inWeb,'Val',0)
	set(handles.edit_forWeb,'Vis','off')
	set([handles.popup_directory_list handles.push_change_dir], 'Vis','on')

% ----------------------------------------------------------------------------
function radio_inWeb_CB(hObject, handles)
	if (~get(hObject,'Val')),		set(hObject,'Val',1),	return,		end
	set(handles.radio_inloco,'Val',0)
	set(handles.edit_forWeb,'Vis','on')
	set([handles.popup_directory_list handles.push_change_dir], 'Vis','off')

% ----------------------------------------------------------------------------
function [gcp, proj] = getProjFromOzi(fname)
	% Decode relevant info in a .map Ozi reference file 
	fid = fopen(fname);
	c = fread(fid,inf,'*char');
	fclose(fid);
	map = strread(c,'%s','delimiter','\n');   clear c fid;
	
	% Now get the GCPs
	ind = strfind(map{10},',');
	% We are decoding strings like this one
	%Point01,xy,  337,  336,in, deg,    ,        ,N,    ,        ,W, grid,   ,     185000,     579000,N
	gcp = zeros(30,4);
	k = 1;
	while ( map{k+9}(ind(3)-1) ~= ' ')		% While we have a pixel<->coord point pair
		gcp(k,1) = str2double( map{k+9}((ind(2)+1) : (ind(3)-1)) ); 
		gcp(k,2) = str2double( map{k+9}((ind(3)+1) : (ind(4)-1)) ); 
		gcp(k,3) = str2double( map{k+9}((ind(14)+1) : (ind(15)-1)) ); 
		gcp(k,4) = str2double( map{k+9}((ind(15)+1) : (ind(16)-1)) ); 
		k = k + 1;
	end

	% Remove non-used pre-allocated gcp
	gcp = gcp(1:(k-1),:);

	% Get mapp projection - THE ONLY TWO CASES THAT I SAW. OBVIOUSLY THIS SHOULD NOT BE APPLYED TO OTHER OZI FILES
	ind = strfind(map{9},',');
	proj = '';
	if ( strncmp(map{9}( ind(1)+1:ind(2)-1 ),'Transverse',10) )		% Transverse Mercator
		proj = '+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111111 +k=1.0 +x_0=200000 +y_0=300000 +ellps=intl +towgs84=-302.581,-61.360,103.047,0,0,0,0';
	elseif ( strncmp(map{9}( ind(1)+1:ind(2)-1 ),'(UTM)',5) )		% (UTM) Universal Transverse Mercator
		proj = '+proj=utm +zone=29 +k=0.9996 +ellps=intl +towgs84=-85.858,-108.681,-120.361,0,0,0,0';
	end

% -----------------------------------------------------------------------------------------
function trans = AffineTransform(uv,xy)
% For an affine transformation:
%                     [ A D 0 ]
% [u v 1] = [x y 1] * [ B E 0 ]
%                     [ C F 1 ]
% There are 6 unknowns: A,B,C,D,E,F
% Another way to write this is:
%                   [ A D ]
% [u v] = [x y 1] * [ B E ]
%                   [ C F ]
% Rewriting the above matrix equation:
% U = X * T, where T = reshape([A B C D E F],3,2)
%
% With 3 or more correspondence points we can solve for T,
% T = X\U which gives us the first 2 columns of T, and
% we know the third column must be [0 0 1]'.

	K = 3;      M = size(xy,1);     X = [xy ones(M,1)];
	U = uv;         % just solve for the first two columns of T

	% We know that X * T = U
	if rank(X) >= K
		Tinv = X \ U;
	else
		errordlg(sprintf('At least %d non-collinear points needed to infer %s transform.',K,'affine'),'Error');
	end

	Tinv(:,3) = [0 0 1]';       % add third column
	trans = inv(Tinv);
	trans(:,3) = [0 0 1]';

%--------------------------------------------------------------------------
function setSliders(handles, sldT, H)
% Create a pair of sliders and register them to use in 'imscroll_j'

	uBak = get(handles.axes1, 'Units');		set(handles.axes1, 'Units', 'pixels')
	axPos = get(handles.axes1,'Pos');		set(handles.axes1, 'Units', uBak)

    hSliders = getappdata(handles.axes1,'SliderAxes');
	if (isempty(hSliders))
		sliderVer = uicontrol('Units','pixels','Style','slider','Parent',handles.figure1, 'Visible','on',...
            'Pos',[axPos(1)+axPos(3)-sldT axPos(2) sldT axPos(4)+1],'Background',[.9 .9 .9]);
		sliderHor = uicontrol('Units','pixels','Style','slider','Parent',handles.figure1, 'Visible','on',...
            'Pos',[axPos(1) H axPos(3)+1 sldT],'Background',[.95 .95 .95]);
		set(sliderHor,'Min',0,'Max',1,'Value',0,'Tag','HOR','Callback',{@slider_Cb,handles.axes1,'SetSliderHor'})
		set(sliderVer,'Min',0,'Max',1,'Value',0,'Tag','VER','Callback',{@slider_Cb,handles.axes1,'SetSliderVer'})
		% Register the sliders in the axe's appdata
		drawnow			% Needed because otherwise - BUG BUG - next line would change the VER slider thickness !!!!!
		setappdata(handles.axes1,'SliderAxes',[sliderHor sliderVer])
		imscroll_j(handles.axes1,'ZoomSetSliders')              % ...
	else			% We have them already. They just need to be updated
        set(hSliders(1), 'Pos',[axPos(1)+axPos(3)+1 axPos(2) sldT axPos(4)+1],'Vis','off')
        set(hSliders(2), 'Pos',[axPos(1) H-1 axPos(3)+1 sldT],'Vis','off')
	end

% -----------------------------------------------------------------------------------------
function slider_Cb(obj,evt,ax,opt)
	imscroll_j(ax,opt)

% -----------------------------------------------------------------------------
function figure1_KeyPressFcn(hObj, event)
	handles = guidata(hObj);
	if isequal(get(hObj,'CurrentCharacter'),'+')
		zoom_j(hObj,2,[]);
	elseif isequal(get(hObj,'CurrentCharacter'),'-')
		zoom_j(hObj,0.5,[]);
	end
	hSliders = getappdata(handles.axes1,'SliderAxes');
	if (~isempty(hSliders) && strcmp( get(hSliders(1),'Vis'),'on' ) )	% If (1) is visible so is the other
		CK = get(hObj,'CurrentKey');
		if (strcmp(CK,'rightarrow') || strcmp(CK,'leftarrow'))
			SS = get(hSliders(1),'SliderStep');			val = get(hSliders(1),'Value');
			if (CK(1) == 'r'),		newVal = min(val + SS(1), 1);	% I know that imscroll_j sliders are [0 1]
			else					newVal = max(0, val - SS(1));
			end
			set(hSliders(1),'Value', newVal)
			imscroll_j(handles.axes1,'SetSliderHor')
		elseif (strcmp(CK,'uparrow') || strcmp(CK,'downarrow'))
			SS = get(hSliders(2),'SliderStep');			val = get(hSliders(2),'Value');
			if (CK(1) == 'u'),		newVal = min(val + SS(1), 1);	% I know that imscroll_j sliders are [0 1]
			else					newVal = max(0, val - SS(1));
			end
			set(hSliders(2),'Value', newVal)
			imscroll_j(handles.axes1,'SetSliderVer')
		end
	end

% --- Creates and returns a handle to the GUI figure. 
function cartas_militares_LayoutFcn(h1)

	rootUnits = get(0, 'Units');			set(0, 'Units', 'pixels');
	screenSize = get(0, 'ScreenSize');      set(0, 'Units', rootUnits);
	FigHeight = screenSize(4) - 65;			FigWidth = round(FigHeight * 0.5);
	AxHeight  = FigHeight - 60;
	pos = [520 113 FigWidth FigHeight];

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'DoubleBuffer','on',...
'MenuBar','none',...
'Name','Cartas Militares',...
'NumberTitle','off',...
'Pos',pos,...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','Call',...
'Tag','figure1');

	pos = [0 60 FigWidth AxHeight];

axes('Parent',h1,...
'Units','pixels',...
'Box','on',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'Pos',pos,...
'Tag','axes1');

uicontrol('Parent',h1, 'Pos',[3 3 FigWidth-40 22],...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Enter web site address where we can find the image files',...
'Tag','edit_forWeb',...
'Visible','off');

uicontrol('Parent',h1, 'Pos',[10 3 FigWidth-40 22],...
'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'Style','popupmenu',...
'String', {' '},...
'Tooltip','Select the directory where the image files reside',...
'Value',1,...
'Tag','popup_directory_list');

uicontrol('Parent',h1, 'Pos',[FigWidth+10-41 3 21 23],...
'Call',@main_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'Tooltip','Select a different directory',...
'Tag','push_change_dir');

uicontrol('Parent',h1, 'Pos',[10 30 61 15],...
'Call',@main_uiCB,...
'FontName','Helvetica',...
'String','In loco',...
'Style','radiobutton',...
'Tooltip','When you have the files on disk ',...
'Value',1,...
'Tag','radio_inloco');

uicontrol('Parent',h1, 'Pos',[130 30 71 15],...
'Call',@main_uiCB,...
'FontName','Helvetica',...
'String','In Web',...
'Style','radiobutton',...
'Tooltip','When you want to get the files from Web',...
'Tag','radio_inWeb');

uicontrol('Parent',h1, 'Position',[FigWidth-110 26 100 21],...
'Call',@main_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','Faz Mosaico',...
'Vis','off',...
'Tag','push_lidarMosaico');

function main_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

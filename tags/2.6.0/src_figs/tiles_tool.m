function varargout = tiles_tool(varargin)
% Get tiles and mosaic them from Virtual Earth (and others) web tile servers 

%	Copyright (c) 2004-2014 by J. Luis
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

	hObject = figure('Vis','off');
	tiles_tool_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	if (numel(varargin) > 0)
		handMir = varargin{1};
        handles.path_data = handMir.path_data;
	else
        handles.path_data = [pwd filesep 'data' filesep];
	end

	% -------------- Import icons -----------------------------------------------
	load([handles.path_data 'mirone_icons.mat'],'zoom_ico','tools_ico','help_ico','um_ico','dois_ico','tres_ico','Mplay_ico');

	handles.uistbar = zeros(1,8);
	h_tb = uitoolbar('parent',hObject,'Clipping', 'on', 'BusyAction','queue','HandleVisibility','on',...
       'Interruptible','on','Tag','FigureToolBar','Visible','on');
	handles.uistbar(1) = uitoggletool('parent',h_tb,'Click',{@click_whatkind_CB,1},'cdata',um_ico, 'Tooltip','Only Satellite Images', 'State', 'on');
	handles.uistbar(2) = uitoggletool('parent',h_tb,'Click',{@click_whatkind_CB,2},'cdata',dois_ico, 'Tooltip','Only Map');
	handles.uistbar(3) = uitoggletool('parent',h_tb,'Click',{@click_whatkind_CB,3},'cdata',tres_ico, 'Tooltip','Map over Satellite');
	handles.uistbar(4) = uipushtool('parent',h_tb,'Click',@click_source_CB,'cdata',tools_ico, 'Tooltip','Select Source servers','Sep','on');
	handles.uistbar(5) = uitoggletool('parent',h_tb,'Click',@click_zoom_CB,'cdata',zoom_ico, 'Tooltip','Zoom','Sep','on');
	ico = icon_ancora;		ico(~ico) = NaN;		ico = single(ico / 255);
	handles.uistbar(6) = uipushtool('parent',h_tb,'Click',@click_anchor_CB,'cdata',ico, 'Tooltip','Set zoom anchor point','Sep','on');
	handles.uistbar(7) = uipushtool('parent',h_tb,'Click',@click_MOSAIC_e_GO_CB,'cdata',Mplay_ico, 'Tooltip','Get necessary images and build Mosaic','Sep','on');
	handles.uistbar(8) = uipushtool('parent',h_tb,'Click',@click_help_CB,'cdata',help_ico, 'Tooltip','Help','Sep','on');

	handles.proxy = [];			handles.port = [];
	handles.proxyPort = [];		handles.patchHandles = [];		handles.hImgZoomed = [];
	handles.whatkind = {'aerial' 'road', 'hybryd'};		% The 3 possible imge types
	handles.slected_whatkind = 'aerial';

	w_map = flipdim(imread([handles.path_data 'etopo2.jpg']),1);
	handles.hImg = image([-180 180],[-90 90],w_map,'Parent',handles.axes1);
	set(handles.axes1, 'XLim',[-180 180], 'YLim',[-90 90],'YDir','normal')	% Oblige limits to be what we want

	% --------------------- Read the cache directory list from mirone_pref ----------------------
	load([handles.path_data 'mirone_pref.mat']);
	try			cacheDirs = cacheTilesDir;	% Try if we already have a cache directory store in prefs
	catch,		cacheDirs = [];
	end

	if ( ~isempty(cacheDirs) )			% We have cache dir(s), but make sure it exists
		ind = true(1,numel(cacheDirs));
		for (k = 1:numel(ind))					% dir(s) are stored in a cell array
			if ( exist(cacheDirs{k}, 'dir') == 7)
				ind(k) = false;
			end
		end
		cacheDirs(ind) = [];			% Remove eventual non existent dirs
		if ( isempty(cacheDirs) ),		cacheDirs = {''};		val = 1;		% It can't be an empty var
		else							cacheDirs = [{''}; cacheDirs];		val = 2;
		end
		set(handles.popup_directory_list,'String',cacheDirs, 'Val', val)
	end
	% ------------------------------------------------------------------------------------------

	% ---------- Read the tilesServers.txt file and swalow its contents ------------------------
	% NOTE: an important part of this section is futuristic code
	handles.serversImageOut = 'http://a0.ortho.tiles.virtualearth.net/tiles/a';		% DEFAULT TO VE
	handles.serversRoadOut = 'http://r0.ortho.tiles.virtualearth.net/tiles/r';
	handles.serversHybridOut = 'http://h0.ortho.tiles.virtualearth.net/tiles/h';
	handles.serversQuadkeyOut = {'0' '1'; '2' '3'};
	handles.serversOrder = ones(3,1);
	
	fid = fopen([handles.path_data 'tilesServers.txt']);
	if (fid > 0)				% Otherwise, non-existent file, revert to defaults
		c = (fread(fid,'*char'))';      fclose(fid);
		servers = strread(c,'%s','delimiter','\n');   clear c fid;
		% Remove comment and also eventual empty lines
		m = numel(servers);		c = false(m,1);
		for (k = 1:m)
			if ( isempty(servers{k}) || servers{k}(1) == '#' ),		c(k) = true;	end
		end
		servers(c) = [];			n_servers = numel(servers);

		handles.servers_image  = cell(1,n_servers);		handles.servers_road = cell(1,n_servers);
		handles.servers_hybrid = cell(1,n_servers);		handles.servers_quadkey = cell(1,n_servers);
		%handles.serversOrder = ones(n_servers, 1);
		for (k = 1:n_servers)					% Loop over number of servers (not that many, we hope)
			[tok,rem] = strtok(servers{k});
			handles.servers_image{k} = tok;		% Images server
			if ( ~isempty(rem) )				% Maps server
				[tok,rem] = strtok(rem);		handles.servers_road{k} = tok;
			end
			if ( ~isempty(rem) )				% Hybrids server
				[tok,rem] = strtok(rem);		handles.servers_hybrid{k} = tok;
			end
			if ( ~isempty(rem) )				% Quadtree keeword
				tok = strtok(rem);				handles.servers_quadkey{k} = tok;
			end
		end
	end
	% -------------------------------------------------------------------------------------------

	% ------------- Some tooltips strings -------------------------
	str = sprintf('Return image in its original Mercator coordinates\n NOTE: THIS IS THE CORRECT THING TO DO');
	set(handles.radio_mercator, 'Tooltip', str)
	str = sprintf('Return image in approximate geographical coordinates\n WARNING: THIS IS IS ONLY A CRUDE APPROXIMATION');
	set(handles.radio_geogs, 'Tooltip', str)

	% ----------------- Finish sliders configurations ------------------------------------------
	st = [1 1] / 23;		% 23 is the currently maximum as far as i know
	set(handles.slider_zoomFactor,'Min',0,'Max',23,'Val',0,'SliderStep',st) 
	
	set(handles.HOR,'Call',{@slider_Cb,handles.axes1,'SetSliderHor'})
	set(handles.VER,'Call',{@slider_Cb,handles.axes1,'SetSliderVer'})
	% Register the sliders in the axe's appdata
	setappdata(handles.axes1,'SliderAxes',[handles.HOR handles.VER])

	% Choose default command line output for tiles_tool
	handles.output = hObject;
	guidata(hObject, handles);

	set(hObject,'Vis','on');
	if (nargout), 	varargout{1} = hObject;		end

% -----------------------------------------------------------------------------------------
function click_whatkind_CB(obj, evt, opt)
	if ( strcmp(get(obj, 'State'), 'off') )			% Don't let deselect a toggled button
		set(obj,'State', 'on'),		return
	end
	handles = guidata(obj);
	all3 = [1 2 3];
	ind = ([opt opt opt] ~= all3);
	other_two = all3(ind);
	set(handles.uistbar(other_two), 'State', 'off')		% Set the other two buttons to deselected
	handles.slected_whatkind = handles.whatkind{opt};	% Store currently selection
	guidata(handles.figure1, handles);

% -----------------------------------------------------------------------------------------
function click_zoom_CB(hObject, evt)
	handles = guidata(hObject);
	if ( strcmp(get(hObject,'State'), 'on') )
		zoom_j(handles.figure1, 'on')
	else
		zoom_j(handles.figure1, 'off')
	end

% -----------------------------------------------------------------------------------------
function click_source_CB(hObject, evt)
% Get the response of the "Tiles servers" window
	handles = guidata(hObject);
	order = handles.serversOrder;		% Order of last call
	out = tiles_servers(handles.servers_image, handles.servers_road, handles.servers_hybrid, order);
	% out.order contains the order that the out.whatkind{i} has in the tilesServers file
	% For example [1 2 1] means that Image and Hybrid are to be fetched from the servers specified
	% on the first line of tilesServers.txt whilst Maps server is the second entry of that file
	if ( ~isempty(out) )
		handles.serversImageOut = out.whatkind{1};
		handles.serversRoadOut = out.whatkind{2};
		handles.serversHybridOut = out.whatkind{3};
		handles.serversOrder = out.order;
		guidata(handles.figure1, handles)
	end

% -----------------------------------------------------------------------------------------
function popup_directory_list_CB(hObject, handles, opt)
% OPT, used by push_change_dir, is char array

	if (nargin == 2)    opt = [];   end
	if ( ~isempty(opt) )				% Add a new entry to the cache dir list. Otherwise, just normal popup functioning.
		contents = get(hObject, 'String');
		if ( numel(contents) == 1 ),	rest = [];
		else							rest = contents(2:end);
		end

		cacheTilesDir = [{opt}; rest];			% Also the var that will be saved in 'mirone_pref'
		set(hObject, 'String', [{''}; cacheTilesDir], 'Val', 2)			% Empty, <=> no cache, always on top
		save([handles.path_data 'mirone_pref.mat'],'cacheTilesDir', '-append', '-v6')		% Update the prefs file
	end

% -----------------------------------------------------------------------------------------
function push_change_dir_CB(hObject, handles)
	if (strcmp(computer, 'PCWIN'))
		cache_dir = uigetfolder_win32('Select a directory', cd);
	else            % This guy doesn't let to be compiled
		cache_dir = uigetdir(cd, 'Select a directory');
	end
	if ~isempty(cache_dir)
		popup_directory_list_CB(handles.popup_directory_list, handles, cache_dir)
	end

% -----------------------------------------------------------------------------------------
function push_clearCacheInfo_CB(hObject, handles)
	cacheTilesDir = [];
	save([handles.path_data 'mirone_pref.mat'],'cacheTilesDir', '-append', '-v6')		% Update the prefs file

% -----------------------------------------------------------------------------------------
function click_anchor_CB(hObject, event)
	handles = guidata(hObject);
	pt = click_e_point(1,'crosshair');
	hAnchor = findobj(handles.axes1,'Tag','Anchor');
	if (isempty(hAnchor))
		hAnchor = line(pt(1,1),pt(1,2),'Marker','p','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',12,'Tag','Anchor');
	else
		set(hAnchor,'XData',pt(1,1) ,'YData',pt(1,2))
	end
	ui_edit_polygon(hAnchor)				% Set edition functions

% -----------------------------------------------------------------------------------------
function click_MOSAIC_e_GO_CB(hObject, event)
% Test if everything is ok and call mosaicing function with (which outputs to Mirone)

	handles = guidata(hObject);
	ind = get(handles.patchHandles, 'UserData');
	if (isempty(ind)),		return,		end			% No patches ploted
	ind = logical(cat(1,ind{:}));
	if ( ~any(ind) ),		return,		end			% No patch selected
	
	% ---------------- Have cache info? -----------------
	val = get(handles.popup_directory_list,'Value');
	contents = get(handles.popup_directory_list, 'String');
	str = contents{val};		cacheDir = [];
	if ( ~isempty(str) ),		cacheDir = str;		end

	tiles_bb = handles.tiles_bb(ind,:);
	lon = [min(tiles_bb(:,1)) max(tiles_bb(:,2))] + [1 -1] * 1e-6;		% eps is for not getting neighboring tiles
	lat = [min(tiles_bb(:,3)) max(tiles_bb(:,4))] + [1 -1] * 1e-6;
	zoomLevel = get(handles.slider_zoomFactor, 'Val') + 1;
	addLevel = get(handles.popup_addLevel, 'Val') - 1;					% Case user wants to download higher than Slider's level
	zoomLevel = zoomLevel + addLevel;
	if (zoomLevel > 18 && handles.slected_whatkind(1) == 'r')			% Roads have 18 as zoom level limit (hybrid???)
		zoomLevel = 18;
	end

	[whatkind, source_PN, source_PV] = get_kind(handles);

	if ( get(handles.radio_geogs, 'Val') )		% Somewhat wrong but practical
		url2image('callmir',lon,lat, zoomLevel, 'lonlat', 'yes', 'cache', cacheDir, ...
			'what',whatkind, source_PN, source_PV, 'verbose','yes');
	else
		url2image('callmir',lon,lat, zoomLevel, 'cache', cacheDir, 'what',whatkind, source_PN, source_PV, 'verbose','yes');
	end

% -----------------------------------------------------------------------------------------
function [whatkind, source_PN, source_PV] = get_kind(handles)
% Get the type of imagery selected by the (1) (2) (3) buttons

	whatkind = handles.slected_whatkind;
	source_PN = 'treta';		source_PV = [];		% Dumb value used to default to VE 
	if ( whatkind(1) == 'a' && isempty(strfind(handles.serversImageOut, 'virtualearth')) )
		server = handles.serversImageOut;		quadkey = handles.servers_quadkey(handles.serversOrder(1));
		if ( ~strncmp(server, 'http', 4) ),		server = ['http://' server];	end
		source_PN = 'source';					source_PV = {server; quadkey{1}};
	elseif ( whatkind(1) == 'r' && isempty(strfind(handles.serversRoadOut, 'virtualearth')) )
		server = handles.serversRoadOut;		quadkey = handles.servers_quadkey(handles.serversOrder(2));
		if ( ~strncmp(server, 'http', 4) ),		server = ['http://' server];	end
		source_PN = 'source';					source_PV = {server; quadkey{1}};
	elseif ( whatkind(1) == 'h' && isempty(strfind(handles.serversHybridOut, 'virtualearth')) )
		server = handles.serversHybridOut;		quadkey = handles.servers_quadkey(handles.serversOrder(3));
		if ( ~strncmp(server, 'http', 4) ),		server = ['http://' server];	end
		source_PN = 'source';					source_PV = {server; quadkey{1}};
	end

% -----------------------------------------------------------------------------------------
function region2tiles(handles,lon,lat,zoomFactor)
% Get tiles BoundingBoxs and plot patches for each tile of the region delimited by LON, LAT
	% Maybe one day we can devise a strategy to reuse the ones that already exist, but meanwhile ... 
	if (~isempty(handles.patchHandles))
		try		delete(handles.patchHandles),	handles.patchHandles = [];	end
	end
	if ( zoomFactor == 1 ),		return,		end

	url = url2image('tile2url', lon, lat, zoomFactor,'quadonly',1);
	[lims, tiles_bb]  = url2image('quadcoord', url);
	[m,n] = size(url);		hp = zeros(m, n);

	k = 1;
	for (j = 1:n)			% col
		xp = [tiles_bb(k,1) tiles_bb(k,1) tiles_bb(k,2) tiles_bb(k,2)];
		for (i = 1:m)
			yp = [tiles_bb(k,3) tiles_bb(k,4) tiles_bb(k,4) tiles_bb(k,3)];
			hp(i,j) = patch('XData',xp,'YData',yp,'Parent',handles.axes1,'FaceColor','none','EdgeColor','k', ...
				'UserData',0, 'ButtonDownFcn',@bdn_tile, 'HitTest', 'on');
			k = k + 1;
		end
	end

	handles.tiles_bb = tiles_bb;
	handles.patchHandles = hp;
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function bdn_tile(obj,eventdata)
	stat = get(gcbo,'UserData');
	if (~stat),		set(gcbo,'FaceColor','y','UserData',1)        % If not selected
	else			set(gcbo,'FaceColor','none','UserData',0)
	end
	refresh

% -------------------------------------------------------------------------------------
function slider_zoomFactor_CB(hObject, handles)
%
	zoomLevel = round(get(hObject,'Value')) + 1;
	lon = get(handles.axes1,'XLim')+[1e-6 -1e-6];	lat = get(handles.axes1,'YLim');
	lat(1) = max(lat(1), -85);			lat(2) = min(lat(2), 85);

	hAnchor = findobj(handles.axes1,'Tag','Anchor');	% See if we have zooming anchor point
	anchor_pt = [];
	if (~isempty(hAnchor))
		anchor_pt = [get(hAnchor,'XData') get(hAnchor,'YData')];	% Yes, so send this to the zoom fun
	end

	nXpatch = getPixel(lon, zoomLevel);
	if (nXpatch > 16)
		zoom_j(handles.figure1, 2, anchor_pt)
		lon = get(handles.axes1,'XLim');	lat = get(handles.axes1,'YLim');
		lon = lon + [-.2 +.2]*diff(lon);	lat = lat + [-.2 +.2]*diff(lat);	% Add 20% on each side to create squares beyound visible
		lon(1) = max(lon(1), -180);			lon(2) = min(lon(2), 180);
		lat(1) = max(lat(1), -85);			lat(2) = min(lat(2), 85);
	elseif (nXpatch < 15)
		zoom_j(handles.figure1, 0.5, anchor_pt)
		lon = get(handles.axes1,'XLim');	lat = get(handles.axes1,'YLim');
		lon = lon + [-.2 +.2]*diff(lon);	lat = lat + [-.2 +.2]*diff(lat);
		lon(1) = max(lon(1), -180);			lon(2) = min(lon(2), 180);
		lat(1) = max(lat(1), -85);			lat(2) = min(lat(2), 85);
	end

	if (zoomLevel > 8)				% At higher zoom levels the bg image is too poor, so get a better one
		bgZoomLevel = zoomLevel - 3;
		val = get(handles.popup_directory_list,'Value');			% Have cache info?
		contents = get(handles.popup_directory_list, 'String');
		str = contents{val};		cacheDir = [];
		if ( ~isempty(str) ),		cacheDir = str;		end
		[whatkind, src_PN, src_PV] = get_kind(handles);			% Decide based on the (1)(2)(3) buttons
		[img, hdr] = ...
			url2image('tile2img',lon,lat, bgZoomLevel, 'cache', cacheDir, 'what',whatkind, src_PN, src_PV,'lonlat', 'yes', 'verbose','y');
		h = image('XData',hdr.X, 'YData',hdr.Y, 'CData',img, 'Parent', handles.axes1);
		if (ishandle(handles.hImgZoomed)),		delete(handles.hImgZoomed),		end		% Delete the previous zoomed image
		handles.hImgZoomed = h;
		guidata(handles.figure1, handles)
	elseif (ishandle(handles.hImgZoomed))
		delete(handles.hImgZoomed)
	end

	region2tiles(handles,lon,lat,zoomLevel)
	set(handles.text_zoomFactor, 'String', zoomLevel)
	if (~isempty(hAnchor))					% If we have one, make sure it's always on top
		Parent = get(hAnchor,{'Parent'});
		Children = allchild(Parent{1});
		Children = setxor(Children, hAnchor);
		if (ishandle(handles.hImgZoomed))
			Children = setxor(Children, handles.hImgZoomed);
			Children = [Children; handles.hImgZoomed];
		end
		Children = setxor(Children, handles.hImg);
		Children = [hAnchor; Children; handles.hImg];
		set(handles.axes1,'Children',Children);
	end
	
% ----------------------------------------------------------------------------	
function nXpatch = getPixel(lon, zoomL)
% Compute the number of tiles from LON(1) to LON(2) at zomm level ZOOML.
% We don't do it for LAT because it would have to be done in isometric latitudes
% x is the fractional number of the 256 bins counting from origin
	pixPerDeg = 2^(zoomL - 1) / 360;
	x = (lon + 180) * pixPerDeg;
	nXpatch = fix(x(2)) - fix(x(1));

% -------------------------------------------------------------------------------------
function slider_Cb(obj,evt,hAxes,opt)
% Control side image sliders
	imscroll_j(hAxes,opt)

% -------------------------------------------------------------------------------------
function radio_mercator_CB(hObject, handles)
	if (get(hObject, 'Val')),		set(handles.radio_geogs, 'Val', 0)
	else							set(hObject, 'Val', 1)
	end

% -------------------------------------------------------------------------------------
function radio_geogs_CB(hObject, handles)
	if (get(hObject, 'Val')),		set(handles.radio_mercator, 'Val', 0)
	else							set(hObject, 'Val', 1)
	end

% -------------------------------------------------------------------------------------
function check_proxy_CB(hObject, handles)
	if ( get(hObject, 'Val') )
		set([handles.edit_proxy handles.edit_port handles.text_port], 'Enable', 'on')
		ind = strfind(handles.proxyPort, ':');
		if ( ~isempty(ind) )		% We already have a full proxy adress
			set_gmt(['http_proxy=' handles.proxyPort]);
		end
	else
		set([handles.edit_proxy handles.edit_port handles.text_port], 'Enable', 'off')
		set_gmt('http_proxy=');
	end
	guidata(handles.figure1, handles)		% as a request to change the proxy settings

% -------------------------------------------------------------------------------------
function edit_proxy_CB(hObject, handles)
	handles.proxy = get(hObject, 'String');
	handles.proxyPort = [handles.proxy ':' handles.port];
	if (~isempty(handles.port))
		set_gmt(['http_proxy=' handles.proxyPort]);
	end
	guidata(handles.figure1, handles)
	
% -------------------------------------------------------------------------------------
function edit_port_CB(hObject, handles)
	handles.port = get(hObject, 'String');
	handles.proxyPort = [handles.proxy ':' handles.port];
	if (~isempty(handles.port))
		set_gmt(['http_proxy=' handles.proxyPort]);
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function click_help_CB(hObject, eventdata)
msg{1} = ['Tool to select several tile images from internet (or local cache) and mosaic them in a single image. ' ...
		'Image resolution is selected with the "Zoom Level" slider, which will create a mesh overlain over the background ' ...
		'image map. Clicking on the mesh thus created turns individual rectangles to yelow, indicating that they are selected. ' ...
		'There are 3 types of images: Satelitte only (button with icon "1"); Maps (button with icon "2"); or Hybrid ' ...
		'(the other button). Hit the the button with the green arrow to fetch data and build the mosaic.'];
msg{2} = ' ';
msg{3} = ['Tiles images are obtained, by default, from MS servers - the same that serve the live.com maps (Virtual Earth). ' ...
		'However it is possible to declare other servers. To do it, follow the instruction on top of the "tilesServers.txt" ' ...
		'that resides on the "data" subdirectory of Mirone''s instalation. If you have more than one set of servers, you can ' ...
		'select among them by clicking on the button with hammers icon. '];
msg{4} = ' ';
msg{5} = ['CACHE: the cache directory is very important since imgas will be searched there before trying to download them ' ...
		'from internet. The cache directory structure can have one of the two following forms.'];
%msg{5} = ' ';
msg{6} = sprintf(['1) NASA World Wind type cache. There is a plugin for WW that allows downloading of Virtual Earth tiles. ' ...
	'In such cases, tiles are stored in a directory tree typicaly like:\n' ...
	'  C:\\Program Files\\NASA\\World Wind 1.4\\Cache\\Virtual Earth\n' ...
	'So, if you have WW installed with the VE plugin, give the above adress to be used as a search ' ...
	'path to VE files before trying to dowload them.\n' ...
'2) C:\\whatever\\you\\want\\to\\call\\this\\path\\cache\n' ...
	'In this case 3 further subdirectories are appended to the CACHE dir (see also below) to acomodate the ''what'' propertie.\n' ...
		'IF ''what'' == ''aerial'' cache = [cache ''kh'']\n' ...
		'IF ''what'' == ''road''   cache = [cache ''mt'']\n' ...
		'IF ''what'' == ''hybrid'' cache = [cache ''tt'']\n' ...
	'Furthermore, a subsequent directory is still appended based on the ZOOM level required.\n' ...
	'Example of an absolute dir of a 12 zoom level aerial request:\n' ...
		' C:/lixo/cache/kh/12\n\n' ...
'Summary. CACHE is a base name directory of which subdirectories are assumed to exist (or ' ...
'created if they do not) to hold progressive refinement zoom level images.\n' ...
'Please follow EXACTLY one of the two possible forms as explained above. ... Otherwise cache is ignored.']); 

helpdlg(msg,'Help on Tiles tool');

%-------------------------------------------------------------------------------
function img = icon_ancora
%
img(:,:,1) = [...
 0   0   0   0   0   0   3   7   7   3   0   0   0   0   0 0
 0   0   0   0   0   0   7   0   0   7   0   0   0   0   0 0
 0   0   0   0   0   0   7   0   0   7   0   0   0   0   0 0
 0   0   0   0   0   0   3   7   7   3   0   0   0   0   0 0
 0   0   0   0  35  23  52 103 103  52  23  35   0   0   0 0
 0   0   0   0  17  35  70  87  87  70  35  17   0   0   0 0
 0   0   0   0   0   0   0  28  28   0   0   0   0   0   0 0
 0   0   0   0   0   0   0  78  72   0   0   0   0   0   0 0
 0   0   0   0   0   0   0  64  58   0   0   0   0   0   0 0
 1  39  31  26   0   0   0  50  45   0   0   0  25  16  39 0
 1  30  88  77  27   0   0  39  33   0   0  27  36  17  16 1
 0  22  68  53  26   4   0  28  24   0   3  13  32  17  11 0
 0   4  28  72  53  35  30  30  30  17  17  17  17  15   2 0
 0   0   5  33  46  69  83  42  41  34  32  17  12   3   0 0
 0   0   0   2  11  23  15  17  17  15  11   4   1   0   0 0
 0   0   0   0   5   7   8   9   9   8   6   4   0   0   0 0 ];

img(:,:,2) = [...
 0   0   0   0   0   0  37  84  84  37   0   0   0   0   0 0
 0   0   0   0   0   0  84   0   0  84   0   0   0   0   0 0
 0   0   0   0   0   0  84   0   0  84   0   0   0   0   0 0
 0   0   0   0   0   0  37  84  84  37   0   0   0   0   0 0
 0   0   0   0 128 116 141 202 202 141 116 128   0   0   0 0
 0   0   0   0  61 128 165 183 183 165 128  61   0   0   0 0
 0   0   0   0   0   0   0 124 124   0   0   0   0   0   0 0
 0   0   0   0   0   0   0 166 161   0   0   0   0   0   0 0
 0   0   0   0   0   0   0 155 149   0   0   0   0   0   0 0
 2 136 129 108   0   0   0 143 139   0   0   0 106 107 136 2
 4 128 196 184 116   0   0 133 129   0   0 114 144 112 108 3
 0  89 174 157  66  13   0 124 121   0  12  52 135 112  70 0
 0  13 118 178 157 137 132 132 132 112 112 112 112  98  10 0
 0   0  19 104 149 175 188 155 152 138 134 112  76  15   0 0
 0   0   0   2  40  92  99 112 112  99  72  29   2   0   0 0
 0   0   0   0   7   9  11  12  12  11   8   6   0   0   0 0 ];

img(:,:,3) = [...
 0   0   0   0   0   0  66 150 150  66   0   0   0   0   0 0
 0   0   0   0   0   0 150   0   0 150   0   0   0   0   0 0
 0   0   0   0   0   0 150   0   0 150   0   0   0   0   0 0
 0   0   0   0   0   0  66 150 150  66   0   0   0   0   0 0
 0   0   0   0 206 200 216 240 240 216 200 206   0   0   0 0
 0   0   0   0  99 206 221 228 228 221 206  99   0   0   0 0
 0   0   0   0   0   0   0 204 204   0   0   0   0   0   0 0
 0   0   0   0   0   0   0 223 220   0   0   0   0   0   0 0
 0   0   0   0   0   0   0 218 215   0   0   0   0   0   0 0
 3 209 200 165   0   0   0 213 211   0   0   0 165 187 209 3
 5 202 239 234 179   0   0 208 206   0   0 178 217 196 189 5
 0 135 229 220  86  19   0 204 203   0  18  81 211 196 122 0
 0  19 184 231 220 211 208 208 208 196 196 196 196 171  17 0
 0   0  27 149 217 229 235 224 223 213 210 196 134  24   0 0
 0   0   0   3  58 139 173 196 196 173 127  51   2   0   0 0
 0   0   0   0   8  11  14  15  15  13  10   7   0   0   1 0 ];

% -------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObj, eventdata)
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
function tiles_tool_LayoutFcn(h1)

set(h1, 'Position',[320 30 940 520],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Tiles Tool',...
'NumberTitle','off',...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

axes('Parent',h1, 'Units','pixels', 'FontSize',8, 'Position',[30 75 900 450], 'Tag','axes1');

uicontrol('Parent',h1, 'Position',[90 4 281 22],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_tool_uiCB,...
'Style','popupmenu',...
'String',{''},...
'Tooltip','Select a cache directory where to search/save tiles files',...
'Value',1,...
'Tag','popup_directory_list');

uicontrol('Parent',h1, 'Position',[370 5 20 21],...
'Call',@tiles_tool_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'Tooltip','Select a different directory',...
'Tag','push_change_dir');

uicontrol('Parent',h1, 'Position',[398 5 23 21],...
'Call',@tiles_tool_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','C',...
'Tooltip','Clear cache info (NOT the cache itself)',...
'Tag','push_clearCacheInfo');

uicontrol('Parent',h1, 'Position',[450 30 41 19],...
'Tooltip','When dowloading, add this level to the slider level value',...
'String',{'0'; '+1'; '+2'; '+3'; '+4'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_addLevel');

uicontrol('Parent',h1, 'Position',[90 30 300 14],...
'BackgroundColor',[0.96 0.96 0.96],...
'Call',@tiles_tool_uiCB,...
'Max',23,...
'Style','slider',...
'Tag','slider_zoomFactor');

uicontrol('Parent',h1, 'Position',[398 29 20 15],...
'FontSize',10,...
'FontWeight','demi',...
'HorizontalAlignment','left',...
'String','1',...
'Style','text',...
'Tag','text_zoomFactor');

uicontrol('Parent',h1, 'Position',[30 48 900 10],...
'Style','slider',...
'Tag','HOR',...
'Visible','off');

uicontrol('Parent',h1, 'Position',[931 74 10 450],...
'Style','slider',...
'Tag','VER',...
'Visible','off');

uicontrol('Parent',h1, 'Position',[10 5 80 17],...
'String','Cache directory',...
'Style','text');

uicontrol('Parent',h1, 'Position',[550 30 110 15],...
'Call',@tiles_tool_uiCB,...
'String','Mercator (correct)',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_mercator');

uicontrol('Parent',h1, 'Position',[550 8 100 15],...
'Call',@tiles_tool_uiCB,...
'String','Geogs (approx)',...
'Style','radiobutton',...
'Tag','radio_geogs');

uicontrol('Parent',h1, 'Position',[18 28 70 17],...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Zoom Level',...
'Style','text');

uicontrol('Parent',h1, 'Position',[782 27 55 15],...
'Call',@tiles_tool_uiCB,...
'String','proxy?',...
'Style','checkbox',...
'Tag','check_proxy');

uicontrol('Parent',h1, 'Position',[838 25 91 20],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_tool_uiCB,...
'Enable','off',...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','proxy adress here',...
'Tag','edit_proxy');

uicontrol('Parent',h1, 'Position',[839 3 51 20],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_tool_uiCB,...
'Enable','off',...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','port here',...
'Tag','edit_port');

uicontrol('Parent',h1,'Position',[808 7 28 15],...
'Enable','off',...
'String','Port',...
'Style','text',...
'Tag','text_port');

function tiles_tool_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

% ----------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------
function varargout = tiles_servers(varargin)
% ... 

	hObject = figure('Vis','off');
	tiles_servers_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject, 'center')

	if (numel(varargin) < 3)
		delete(hObject)
		error('tiles_servers: input must have 3 arguments')
	end

	handles.servers_image  = varargin{1};
	handles.servers_road   = varargin{2};
	handles.servers_hybrid = varargin{3};
	order = varargin{4};

	% ------- Fill edit and popup menus with just recieved info -----------
	set(handles.edit_aerial,  'String', handles.servers_image{order(1)} )
	set(handles.edit_road,    'String', handles.servers_road{order(2)} )
	set(handles.edit_hybrid,  'String', handles.servers_hybrid{order(3)} )
	set(handles.popup_aerial, 'String', handles.servers_image, 'Val', order(1))
	set(handles.popup_road,   'String', handles.servers_road,  'Val', order(2))
	set(handles.popup_hybrid, 'String', handles.servers_hybrid,'Val', order(3))
	handles.aerial = handles.servers_image{order(1)};					% Initializations
	handles.road   = handles.servers_road{order(2)};
	handles.hybrid = handles.servers_hybrid{order(3)};
	handles.aerial_ind = order(1);		handles.road_ind = order(2);	handles.hybrid_ind = order(3);
	% ---------------------------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Vis','on');

	% UIWAIT makes tiles_servers wait for user response
	uiwait(handles.figure1);
	handles = guidata(hObject);
	varargout{1} = handles.output;
	if ( ishandle(handles.figure1) ),	delete(handles.figure1),	end

% -----------------------------------------------------------------------------------------
function edit_aerial_CB(hObject, handles, opt)
	if (nargin == 3),		handles.aerial = opt;		% Called by popupmenu
	else					handles.aerial = get(hObject, 'String');
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function edit_road_CB(hObject, handles, opt)
	if (nargin == 3),		handles.road = opt;			% Called by popupmenu
	else					handles.road = get(hObject, 'String');
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function edit_hybrid_CB(hObject, handles, opt)
	if (nargin == 3),		handles.hybrid = opt;		% Called by popupmenu
	else					handles.hybrid = get(hObject, 'String');
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function popup_aerial_CB(hObject, handles)
% Get selected server and transmit info to corresponding edit box (which will do the rest) 
	val = get(hObject,'Value');
	contents = get(hObject,'String');		server = contents{val};
	set(handles.edit_aerial, 'String', server)
	handles.aerial_ind = val;
	edit_aerial_CB(hObject, handles, server)

% -----------------------------------------------------------------------------------------
function popup_road_CB(hObject, handles)
% Get selected server and transmit info to corresponding edit box (which will do the rest) 
	val = get(hObject,'Value');
	contents = get(hObject,'String');		server = contents{val};
	set(handles.edit_road, 'String', server)
	handles.road_ind = val;
	edit_road_CB(hObject, handles, server)

% -----------------------------------------------------------------------------------------
function popup_hybrid_CB(hObject, handles)
% Get selected server and transmit info to corresponding edit box (which will do the rest) 
	val = get(hObject,'Value');
	contents = get(hObject,'String');		server = contents{val};
	set(handles.edit_hybrid, 'String', server)
	handles.hybrid_ind = val;
	edit_hybrid_CB(hObject, handles, server)

% -----------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
	out.whatkind = {handles.aerial; handles.road; handles.hybrid};		% Selected server for "whatkind"
	out.order = [handles.aerial_ind; handles.road_ind; handles.hybrid_ind];		% Indices order to recover correspondent quadkeey
	handles.output = out;
	guidata(handles.figure1, handles)
	uiresume(handles.figure1);

% -----------------------------------------------------------------------------------------
function push_cancel_CB(hObject, handles)
	handles.output = [];
	guidata(hObject,handles)
	uiresume(handles.figure1);

% -----------------------------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure_servers_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	if (exist('OCTAVE_VERSION','builtin'))		% To know if we are running under Octave
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	else
		do_uiresume = strcmp(get(handles.figure1, 'waitstatus'), 'waiting');
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		handles.output = [];		% User gave up, return nothing
		guidata(handles.figure1, handles);	uiresume(handles.figure1);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

% -----------------------------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure_servers_KeyPressFcn(hObject, eventdata)
	handles = guidata(hObject);
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles.output = [];		% User said no by hitting escape
		guidata(hObject, handles);	uiresume(handles.figure1);
	end

% --- Creates and returns a handle to the GUI figure. 
function tiles_servers_LayoutFcn(h1)

set(h1, 'Position',[520 541 550 175],...
'CloseRequestFcn',@figure_servers_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure_servers_KeyPressFcn,...
'MenuBar','none',...
'Name','Tiles servers',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 133 251 21],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'HorizontalAlignment','left',...
'Style','edit', 'Tag','edit_aerial');

uicontrol('Parent',h1, 'Position',[10 83 251 21],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'HorizontalAlignment','left',...
'Style','edit', 'Tag','edit_road');

uicontrol('Parent',h1, 'Position',[10 33 251 21],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_hybrid');

uicontrol('Parent',h1, 'Position',[10 55 110 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Hybrid request from:',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 106 110 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Map request from:',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 157 110 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Image request from:',...
'Style','text');

uicontrol('Parent',h1, 'Position',[269 133 271 22],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'Style','popupmenu',...
'Tooltip','Read from ''tileServers.txt'' on ''data'' directory',...
'Value',1,...
'Tag','popup_aerial');

uicontrol('Parent',h1, 'Position',[269 83 271 22],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'Style','popupmenu',...
'Tooltip','Read from ''tileServers.txt'' on ''data'' directory',...
'Value',1,...
'Tag','popup_road');

uicontrol('Parent',h1, 'Position',[269 33 271 22],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'Style','popupmenu',...
'Tooltip','Read from ''tileServers.txt'' on ''data'' directory',...
'Value',1,...
'Tag','popup_hybrid');

uicontrol('Parent',h1, 'Position',[320 157 150 16],...
'FontName','Helvetica',...
'String','Alernatives imported from file',...
'Tooltip','Read from ''tileServers.txt'' on ''data'' directory',...
'Style','text');

uicontrol('Parent',h1, 'Position',[381 3 70 23],...
'Call',@tiles_servers_uiCB,...
'FontName','Helvetica',...
'FontWeight','bold',...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1, 'Position',[471 3 70 23],...
'Call',@tiles_servers_uiCB,...
'FontName','Helvetica',...
'FontWeight','bold',...
'String','Cancel',...
'Tag','push_cancel');

function tiles_servers_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

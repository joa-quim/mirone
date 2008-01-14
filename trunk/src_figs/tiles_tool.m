function varargout = tiles_tool(varargin)
% M-File changed by desGUIDE 
 
	hObject = figure('Tag','figure1','Visible','off');
	tiles_tool_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'east')

	if (numel(varargin) > 0)
		handMir = varargin{1};
        handles.path_data = handMir.path_data;
	else
        handles.path_data = [pwd filesep 'data' filesep];
	end

	% -------------- Import icons -----------------------------------------------
	load([handles.path_data 'mirone_icons.mat'],'zoom_ico','tools_ico','help_ico','um_ico','dois_ico','tres_ico','Mplay_ico');

	handles.uistbar = zeros(1,7);
	h_tb = uitoolbar('parent',hObject,'Clipping', 'on', 'BusyAction','queue','HandleVisibility','on',...
       'Interruptible','on','Tag','FigureToolBar','Visible','on');
	handles.uistbar(1) = uitoggletool('parent',h_tb,'Click',{@click_whatkind_CB,1},'cdata',um_ico, 'Tooltip','Only Satellite Images', 'State', 'on');
	handles.uistbar(2) = uitoggletool('parent',h_tb,'Click',{@click_whatkind_CB,2},'cdata',dois_ico, 'Tooltip','Only Map');
	handles.uistbar(3) = uitoggletool('parent',h_tb,'Click',{@click_whatkind_CB,3},'cdata',tres_ico, 'Tooltip','Map over Satellite');
	handles.uistbar(4) = uipushtool('parent',h_tb,'Click',@click_source_CB,'cdata',tools_ico, 'Tooltip','Select Source servers','Sep','on');
	handles.uistbar(5) = uitoggletool('parent',h_tb,'Click',@click_zoom_CB,'cdata',zoom_ico, 'Tooltip','Zoom','Sep','on');
	handles.uistbar(6) = uipushtool('parent',h_tb,'Click',@click_MOSAIC_e_GO_CB,'cdata',Mplay_ico, 'Tooltip','Get necessary images and build Mosaic','Sep','on');
	handles.uistbar(7) = uipushtool('parent',h_tb,'Click',@click_help_CB,'cdata',help_ico, 'Tooltip','Help','Sep','on');

	handles.proxy = [];			handles.port = [];
	handles.proxyPort = [];		handles.patchHandles = [];
	handles.doProxySetting = false;
	handles.whatkind = {'aerial' 'road', 'hybryd'};		% The 3 possible imge types
	handles.slected_whatkind = 'aerial';

	w_map = flipdim(imread([handles.path_data 'etopo2.jpg']),1);
	image([-180 180],[-90 90],w_map,'Parent',handles.axes1);
	set(handles.axes1, 'XLim',[-180 180], 'YLim',[-90 90],'YDir','normal')	% Oblige limits to be what we want

	load([handles.path_data 'm_coasts.mat']);
	h_boundaries = line('XData',ncst(:,1),'YData',ncst(:,2),'Parent',handles.axes1,'Linewidth',1,'Color','w');

	% --------------------- Read the cache directory list from mirone_pref ----------------------
	load([handles.path_data 'mirone_pref.mat']);
	try			cacheDirs = cacheTilesDir;	% Try if we already have a cache directory store in prefs
	catch		cacheDirs = [];
	end

	if ( ~isempty(cacheDirs) )			% We have cache dir(s), but make sure it exists
		ind = true(1,numel(cacheDirs));
		for (k = 1:numel(ind))					% dir(s) are stored in a cell array
			if ( exist(cacheDirs{k}, 'dir') == 7)
				ind(k) = false;
			end
		end
		cacheDirs(ind) = [];			% Remove eventual non existent dirs
		set(handles.popup_directory_list,'String',cacheDirs)
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
			if ( servers{k}(1) == '#' || isempty(servers{k}) ),		c(k) = true;	end
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
				[tok,rem] = strtok(rem);		handles.servers_quadkey{k} = tok;
			end
		end
	end
	% -------------------------------------------------------------------------------------------

	% ------------- Some tooltips strings -------------------------
	str = sprintf('Return image in its original Mercator coordinates\n NOTE: THIS IS THE CORRECT THING TO DO');
	set(handles.radio_mercator, 'TooltipString', str)
	str = sprintf('Return image in approximate geographical coordinates\n WARNING: THIS IS IS ONLY A CRUDE APPROXIMATION');
	set(handles.radio_geogs, 'TooltipString', str)

	% ----------------- Finish sliders configurations ------------------------------------------
	st = [1 1] / 23;		% 23 is the currently maximum as far as i know
	set(handles.slider_zoomFactor,'Min',0,'Max',23,'Val',0,'SliderStep',st) 
	
	set(handles.HOR,'Callback',{@slider_Cb,handles.axes1,'SetSliderHor'})
	set(handles.VER,'Callback',{@slider_Cb,handles.axes1,'SetSliderVer'})
	% Register the sliders in the axe's appdata
	setappdata(handles.axes1,'SliderAxes',[handles.HOR handles.VER])

	% Choose default command line output for tiles_tool
	handles.output = hObject;
	guidata(hObject, handles);

	set(hObject,'Visible','on');
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
% 		zoom_j(handles.figure1, 'on', {@meshtiles,handles.axes1})
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
function popup_directory_list_Callback(hObject, eventdata, handles, opt)
% OPT is used by pushbutton_change_dir (just to save code)

	if (nargin == 3)    opt = [];   end
	val = get(hObject,'Value');     str = get(hObject, 'String');
	if isempty(opt)
		% Put the selected field on top of the String list.
		tmp = str(val);         str(val) = [];
		cacheTilesDir = [tmp; str];
	else
		cacheTilesDir = [opt; str];
		if (~iscell(cacheTilesDir)),	cacheTilesDir = {cacheTilesDir};	end
		save([handles.path_data 'mirone_pref.mat'],'cacheTilesDir', '-append')		% Update the prefs file
	end
	set(handles.popup_directory_list, 'String',cacheTilesDir ,'Value',1); 

% -----------------------------------------------------------------------------------------
function pushbutton_change_dir_Callback(hObject, eventdata, handles)
	%cache_dir = uigetdir;        % This guy doesn't let to be compiled
	cache_dir = uigetfolder_standalone;
	if ~isempty(cache_dir)
		popup_directory_list_Callback(handles.popup_directory_list, eventdata, handles, cache_dir)
	end

% -----------------------------------------------------------------------------------------
function click_MOSAIC_e_GO_CB(hObject, eventdata)
% Test if everything is ok and call mosaicing function with (which outputs to Mirone)

	handles = guidata(hObject);
	ind = get(handles.patchHandles, 'UserData');
	if (isempty(ind)),		return,		end			% No patches ploted
	ind = logical(cat(1,ind{:}));
	if ( ~any(ind) ),		return,		end			% No patch selected
	
	% ---------------- Have cache info? -----------------
	val = get(handles.popup_directory_list,'Value');
	str = get(handles.popup_directory_list, 'String');
	if ( ~isempty(str) ),		cacheDir = str{val};
	else						cacheDir = [];
	end

	% --------------- Proxy settings --------------------
	if (~isempty(handles.proxyPort))
		ind = strfind(handles.proxyPort, ':');
		if ( isempty(ind) && get(handles.check_proxy,'Val') )		% either adress or port are missing
			warndlg('Your proxy settings are incomplete. Ignoring','Warning')
		elseif (handles.doProxySetting)
			set_gmt(['http_proxy=' handles.proxyPort]);
			handles.doProxySetting = false;
			guidata(handles.figure1, handles)
		end
	end
	% ----------------------------------------------------

	tiles_bb = handles.tiles_bb(ind,:);
	lon = [min(tiles_bb(:,1)) max(tiles_bb(:,2))] + [1 -1] * 1e-6;		% eps is for not getting neighboring tiles
	lat = [min(tiles_bb(:,3)) max(tiles_bb(:,4))] + [1 -1] * 1e-6;
	zoomLevel = get(handles.slider_zoomFactor, 'Val') + 1;
	if (zoomLevel > 18 && handles.slected_whatkind(1) == 'r')			% Roads have 18 as zoom level limit (hybrid???)
		zoomLevel = 18;
	end

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

	if ( get(handles.radio_geogs, 'Val') )		% Somewhat wrong but practical
		url2image('callmir',lon,lat, zoomLevel, 'lonlat', 'yes', 'cache', cacheDir, ...
			'what',whatkind, source_PN, source_PV, 'verbose','yes');
	else
		url2image('callmir',lon,lat, zoomLevel, 'cache', cacheDir, 'what',whatkind, source_PN, source_PV, 'verbose','yes');
	end

% -----------------------------------------------------------------------------------------
function region2tiles(handles,lon,lat,zoomFactor)
% Get tiles BoundingBoxs and plot patches for each tile of the region delimited by LON, LAT
	% Maybe one day we can devise a strategy to reuse the ones that already exist, but meanwhile ... 
	if (~isempty(handles.patchHandles))
		try		delete(handles.patchHandles),	handles.patchHandles = [];	end
	end
	if ( zoomFactor == 1 ),		return,		end

	%profile on
	[url, lonT, latT] = url2image('tile2url', lon, lat, zoomFactor,'quadonly',1);
	%profile viewer
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
function bdn_tile(obj,eventdata,handles)
	stat = get(gcbo,'UserData');
	if ~stat        % If not selected
		set(gcbo,'FaceColor','y','UserData',1),			refresh
	else
		set(gcbo,'FaceColor','none','UserData',0),		refresh
	end

% -------------------------------------------------------------------------------------
function slider_zoomFactor_Callback(hObject, eventdata, handles)
%
	zoomLevel = round(get(hObject,'Value')) + 1;

	lon = get(handles.axes1,'XLim');	lat = get(handles.axes1,'YLim');
	lat(1) = max(lat(1), -85);			lat(2) = min(lat(2), 85);

	nXpatch = getPixel(lon, zoomLevel);
	while (nXpatch > 20 || nXpatch > 20)
		zoom_j(handles.figure1, 2)
		lon = get(handles.axes1,'XLim');	lat = get(handles.axes1,'YLim');
		lat(1) = max(lat(1), -85);			lat(2) = min(lat(2), 85);
		nXpatch = getPixel(lon, zoomLevel);
	end

	region2tiles(handles,lon,lat,zoomLevel)
	set(handles.text_zoomFactor, 'String', zoomLevel)
	
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

% -----------------------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
	delete(handles.figure1);

% -------------------------------------------------------------------------------------
function radio_mercator_Callback(hObject, eventdata, handles)
	if (get(hObject, 'Val'))
		set(handles.radio_geogs, 'Val', 0)
	else
		set(hObject, 'Val', 1)
	end

% -------------------------------------------------------------------------------------
function radio_geogs_Callback(hObject, eventdata, handles)
	if (get(hObject, 'Val'))
		set(handles.radio_mercator, 'Val', 0)
	else
		set(hObject, 'Val', 1)
	end

% -------------------------------------------------------------------------------------
function check_proxy_Callback(hObject, eventdata, handles)
	if ( get(hObject, 'Val') )
		set([handles.edit_proxy handles.edit_port handles.text_port], 'Enable', 'on')
	else
		set([handles.edit_proxy handles.edit_port handles.text_port], 'Enable', 'off')
	end
	handles.doProxySetting = true;			% If the checkbox state was changed, interpret it
	guidata(handles.figure1, handles)		% as a request to change the proxy settings

% -------------------------------------------------------------------------------------
function edit_proxy_Callback(hObject, eventdata, handles)
	handles.proxy = get(hObject, 'String');
	handles.proxyPort = [handles.proxy ':' handles.port];
	if (~isempty(handles.port)),		handles.doProxySetting = true;	end
	guidata(handles.figure1, handles)

% -------------------------------------------------------------------------------------
function edit_port_Callback(hObject, eventdata, handles)
	handles.port = get(hObject, 'String');
	handles.proxyPort = [handles.proxy ':' handles.port];
	if (~isempty(handles.proxy)),		handles.doProxySetting = true;	end
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
		'However it is ppossible to declare other servers. To do it, follow the instruction on top of the "tilesServers.txt" ' ...
		'that resides on the "data" subdirectory of Mirone''s instalation. If you have more than one set of servers, you can ' ...
		'selection among them by clicking on the button with hammers icon. '];
msg{4} = ' ';
msg{3} = ['CACHE: the cache directory is very important since imgas will be searched there before trying to download them ' ...
		'from internet. The cache directory structure can have one of the two following forms.'];
%msg{5} = ' ';
msg{5} = sprintf(['1) NASA World Wind type cache. There is a plugin for WW that allows downloading of Virtual Earth tiles. ' ...
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
		' C:/whatever/you/want/to/call/this/path/cache/kh/12\n\n' ...
'Summary. CACHE is a base name directory of which subdirectories are assumed to exist (or ' ...
'created if they do not) to hold proguessive refinement zoom level images.\n' ...
'Please follow EXACTLY one of the two possible forms as explained above. ... Otherwise cache is ignored.']); 

helpdlg(msg,'Help on Tiles tool');

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
function tiles_tool_LayoutFcn(h1);

set(h1, 'Position',[520 365 765 435],...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
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

axes('Parent',h1,...
'Units','pixels',...
'FontSize',8,...
'Position',[30 75 721 361],...
'Tag','axes1');

uicontrol('Parent',h1,...
'Callback',{@tiles_tool_uicallback,h1,'check_proxy_Callback'},...
'Position',[570 27 72 15],...
'String','Use proxy',...
'Style','checkbox',...
'Tag','check_proxy');

uicontrol('Parent',h1, 'Position',[90 4 281 22],...
'BackgroundColor',[1 1 1],...
'Callback',{@tiles_tool_uicallback,h1,'popup_directory_list_Callback'},...
'Style','popupmenu',...
'TooltipString','Select a cache directory where to search/save tiles files',...
'Value',1,...
'Tag','popup_directory_list');

uicontrol('Parent',h1, 'Position',[371 5 18 21],...
'Callback',{@tiles_tool_uicallback,h1,'pushbutton_change_dir_Callback'},...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'TooltipString','Select a different directory',...
'Tag','pushbutton_change_dir');

uicontrol('Parent',h1, 'Position',[90 30 281 14],...
'BackgroundColor',[0.96 0.96 0.96],...
'Callback',{@tiles_tool_uicallback,h1,'slider_zoomFactor_Callback'},...
'Max',23,...
'Style','slider',...
'Tag','slider_zoomFactor');

uicontrol('Parent',h1, 'Position',[377 29 20 15],...
'FontSize',10,...
'FontWeight','demi',...
'HorizontalAlignment','left',...
'String','1',...
'Style','text',...
'Tag','text_zoomFactor');

uicontrol('Parent',h1, 'Position',[30 50 721 8],...
'Style','slider',...
'Tag','HOR',...
'Visible','off');

uicontrol('Parent',h1, 'Position',[754 74 9 361],...
'Style','slider',...
'Tag','VER',...
'Visible','off');

uicontrol('Parent',h1, 'Position',[10 5 80 17],...
'String','Cache directory',...
'Style','text');

uicontrol('Parent',h1, 'Position',[426 30 110 15],...
'Callback',{@tiles_tool_uicallback,h1,'radio_mercator_Callback'},...
'String','Mercator (correct)',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_mercator');

uicontrol('Parent',h1, 'Position',[426 10 130 15],...
'Callback',{@tiles_tool_uicallback,h1,'radio_geogs_Callback'},...
'String','Geographics (approx)',...
'Style','radiobutton',...
'Tag','radio_geogs');

uicontrol('Parent',h1, 'Position',[18 28 70 17],...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Zoom Level',...
'Style','text');

uicontrol('Parent',h1, 'Position',[640 25 111 20],...
'BackgroundColor',[1 1 1],...
'Callback',{@tiles_tool_uicallback,h1,'edit_proxy_Callback'},...
'Enable','off',...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','proxy adress here',...
'Tag','edit_proxy');

uicontrol('Parent',h1, 'Position',[641 3 51 20],...
'BackgroundColor',[1 1 1],...
'Callback',{@tiles_tool_uicallback,h1,'edit_port_Callback'},...
'Enable','off',...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','port here',...
'Tag','edit_port');

uicontrol('Parent',h1,'Position',[610 7 30 15],...
'Enable','off',...
'String','Port',...
'Style','text',...
'Tag','text_port');

function tiles_tool_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

% ----------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------
function varargout = tiles_servers(varargin)
% M-File changed by desGUIDE 

	hObject = figure('Tag','figure1','Visible','off');
	tiles_servers_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject, 'center')

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
	handles.aerial = handles.servers_image{1}(order(1));					% Initializations
	handles.road   = handles.servers_road{1}(order(2));
	handles.hybrid = handles.servers_hybrid{1}(order(3));
	handles.aerial_ind = 1;		handles.road_ind = 1;	handles.hybrid_ind = 1;
	% ---------------------------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');

	% UIWAIT makes tiles_servers wait for user response
	uiwait(handles.figure1);
	handles = guidata(hObject);
	varargout{1} = handles.output;
	if ( ishandle(handles.figure1) ),	delete(handles.figure1),	end

% -----------------------------------------------------------------------------------------
function edit_aerial_Callback(hObject, eventdata, handles, opt)
	if (nargin == 4),		handles.aerial = opt;		% Called by popupmenu
	else					handles.aerial = get(hObject, 'String');
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function edit_road_Callback(hObject, eventdata, handles, opt)
	if (nargin == 4),		handles.road = opt;			% Called by popupmenu
	else					handles.road = get(hObject, 'String');
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function edit_hybrid_Callback(hObject, eventdata, handles, opt)
	if (nargin == 4),		handles.hybrid = opt;		% Called by popupmenu
	else					handles.hybrid = get(hObject, 'String');
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function popup_aerial_Callback(hObject, eventdata, handles)
% Get selected server and transmit info to corresponding edit box (which will do the rest) 
	val = get(hObject,'Value');
	contents = get(hObject,'String');		server = contents{val};
	set(handles.edit_aerial, 'String', server)
	handles.aerial_ind = val;
	edit_aerial_Callback(hObject, [], handles, server)

% -----------------------------------------------------------------------------------------
function popup_road_Callback(hObject, eventdata, handles)
% Get selected server and transmit info to corresponding edit box (which will do the rest) 
	val = get(hObject,'Value');
	contents = get(hObject,'String');		server = contents{val};
	set(handles.edit_road, 'String', server)
	handles.road_ind = val;
	edit_road_Callback(hObject, [], handles, server)

% -----------------------------------------------------------------------------------------
function popup_hybrid_Callback(hObject, eventdata, handles)
% Get selected server and transmit info to corresponding edit box (which will do the rest) 
	val = get(hObject,'Value');
	contents = get(hObject,'String');		server = contents{val};
	set(handles.edit_hybrid, 'String', server)
	handles.hybrid_ind = val;
	edit_hybrid_Callback(hObject, [], handles, server)

% -----------------------------------------------------------------------------------------
function push_OK_Callback(hObject, eventdata, handles)
	out.whatkind = {handles.aerial; handles.road; handles.hybrid};		% Selected server for "whatkind"
	out.order = [handles.aerial_ind; handles.road_ind; handles.hybrid_ind];		% Indices order to recover correspondent quadkeey
	handles.output = out;           guidata(hObject,handles)
	guidata(handles.figure1, handles)
	uiresume(handles.figure1);

% -----------------------------------------------------------------------------------------
function push_cancel_Callback(hObject, eventdata, handles)
	handles.output = [];
	guidata(hObject,handles)
	uiresume(handles.figure1);

% -----------------------------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure_servers_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
		% The GUI is still in UIWAIT, UIRESUME
		handles.output = [];		% User gave up, return nothing
		guidata(hObject, handles);	uiresume(handles.figure1);
	else    % The GUI is no longer waiting, just close it
		handles.output = [];		% User gave up, return nothing
		guidata(hObject, handles);	delete(handles.figure1);
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
function tiles_servers_LayoutFcn(h1);

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
'Callback',{@tiles_servers_uicallback,h1,'edit_aerial_Callback'},...
'HorizontalAlignment','left',...
'Style','edit', 'Tag','edit_aerial');

uicontrol('Parent',h1, 'Position',[10 83 251 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@tiles_servers_uicallback,h1,'edit_road_Callback'},...
'HorizontalAlignment','left',...
'Style','edit', 'Tag','edit_road');

uicontrol('Parent',h1, 'Position',[10 33 251 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@tiles_servers_uicallback,h1,'edit_hybrid_Callback'},...
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
'Callback',{@tiles_servers_uicallback,h1,'popup_aerial_Callback'},...
'Style','popupmenu',...
'TooltipString','Read from ''tileServers.txt'' on ''data'' directory',...
'Value',1,...
'Tag','popup_aerial');

uicontrol('Parent',h1, 'Position',[269 83 271 22],...
'BackgroundColor',[1 1 1],...
'Callback',{@tiles_servers_uicallback,h1,'popup_road_Callback'},...
'Style','popupmenu',...
'TooltipString','Read from ''tileServers.txt'' on ''data'' directory',...
'Value',1,...
'Tag','popup_road');

uicontrol('Parent',h1, 'Position',[269 33 271 22],...
'BackgroundColor',[1 1 1],...
'Callback',{@tiles_servers_uicallback,h1,'popup_hybrid_Callback'},...
'Style','popupmenu',...
'TooltipString','Read from ''tileServers.txt'' on ''data'' directory',...
'Value',1,...
'Tag','popup_hybrid');

uicontrol('Parent',h1, 'Position',[320 157 150 16],...
'FontName','Helvetica',...
'String','Alernatives imported from file',...
'TooltipString','Read from ''tileServers.txt'' on ''data'' directory',...
'Style','text');

uicontrol('Parent',h1, 'Position',[381 3 70 23],...
'Callback',{@tiles_servers_uicallback,h1,'push_OK_Callback'},...
'FontName','Helvetica',...
'FontWeight','bold',...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1, 'Position',[471 3 70 23],...
'Callback',{@tiles_servers_uicallback,h1,'push_cancel_Callback'},...
'FontName','Helvetica',...
'FontWeight','bold',...
'String','Cancel',...
'Tag','push_cancel');

function tiles_servers_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

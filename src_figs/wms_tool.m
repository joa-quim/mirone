function varargout = wms_tool(varargin)
% Helper tool to construct an URL to download an image from a WMS
%
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

	if (ishandle(varargin{1}))				% An handle to a rectangle
		x = get(varargin{1}, 'XData');
		y = get(varargin{1}, 'YData');
		region = double([min(x) max(x) min(y) max(y)]);
	elseif (numel(varargin{1}) == 4)		% A [xmin xmax ymin ymax] vector
		region = varargin{1};
	else
		errordlg('WMS_TOOL: bad input argument','Error')
		return
	end

	hObject = figure('Vis','off');
	wms_tool_LayoutFcn(hObject);
	handles = guihandles(hObject);

	handles.proxy = [];			handles.port = [];		handles.proxyPort = [];
	handles.region = region;

	data = clock;
	handles.date = sprintf('%d-%.2d-%.2d', data(1:3));
	set(handles.edit_date, 'String', sprintf('%.2d-%.2d-%d',data(3:-1:1)))

	% Compute image size assuming the default resolution of 8 arc-sec/pix
	W = round(diff(handles.region(1:2)) / 8 * 3600);
	H = round(diff(handles.region(3:4)) / 8 * 3600);
	str = sprintf('Try to download a [%d x %d] image', H,W);
	set(handles.push_OK, 'ToolTip', str)
	set(handles.edit_res, 'ToolTip', sprintf('Resolution in arc-sec/pixel. Not lower than 1.\nAt this res\n%s',str))
	
	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% ------------------------------------------------------------------------------------
function edit_res_CB(hObject, handles)
	res = str2double(get(hObject, 'String'));
	if (isnan(res) || res < 1)
		set(hObject, 'String', '8')
		res = 8;
	end
	% Compute image size for this resolution and update tooltips
	W = round(diff(handles.region(1:2)) / res * 3600);
	H = round(diff(handles.region(3:4)) / res * 3600);
	str = sprintf('Try to download a [%d x %d] image', H,W);
	set(handles.push_OK, 'ToolTip', str)
	set(handles.edit_res, 'ToolTip', sprintf('Resolution in arc-sec/pixel. Not lower than 1.\nAt this res we''ll %s',str))

% ------------------------------------------------------------------------------------
function edit_date_CB(hObject, handles)
	data = get(hObject, 'String');
	ind = strfind(data,'-');
	erro = false;
	if (numel(ind) ~= 2)							erro = true;
	elseif (numel(data(1:ind(1)-1)) ~= 2)			erro = true;
	elseif (numel(data(ind(1)+1:ind(2)-1)) ~= 2)	erro = true;
	elseif (numel(data(ind(2)+1:end)) ~= 4)			erro = true;
	end
	if (erro)
		errordlg('Wrongly formated date. Please enter date as dd-mm-yyyy','Error')
		set(hObject, 'String', handles.date)
		return
	end

	%Ok, now that date is in the way I want it. Convert it into the way WMS uses (blheak!!)
	handles.date = [data(7:10) data(3:6) data(1:2)];
	guidata(handles.figure1, handles)

% ------------------------------------------------------------------------------------
function popup_servers_CB(hObject, handles)
% Do nothing for now. Maybe in a near future

% ------------------------------------------------------------------------------------
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

% ------------------------------------------------------------------------------------
function edit_proxy_CB(hObject, handles)
	handles.proxy = get(hObject, 'String');
	handles.proxyPort = [handles.proxy ':' handles.port];
	if (~isempty(handles.port))
		set_gmt(['http_proxy=' handles.proxyPort]);
	end
	guidata(handles.figure1, handles)

% ------------------------------------------------------------------------------------
function edit_port_CB(hObject, handles)
	handles.port = get(hObject, 'String');
	handles.proxyPort = [handles.proxy ':' handles.port];
	if (~isempty(handles.port))
		set_gmt(['http_proxy=' handles.proxyPort]);
	end
	guidata(handles.figure1, handles)

% ------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...
	res = str2double(get(handles.edit_res, 'String'));
	W = round(diff(handles.region(1:2)) / res * 3600);
	H = round(diff(handles.region(3:4)) / res * 3600);
	if (min(W,H) < 512)		% Don't know the true min size
		errordlg(sprintf('Requested image size is too small [%d x %d].', H,W),'Error')
		return
	end
	
	URL = sprintf(['http://onearth.jpl.nasa.gov/wms.cgi?TIME=%s&SERVICE=WMS&LAYERS=daily_planet' ...
				'&EXCEPTIONS=application/vnd.ogc.se_xml&FORMAT=image/jpeg&TRANSPARENT=FALSE' ...
				'&HEIGHT=%d&BGCOLOR=0xFFFFFF&REQUEST=GetMap&WIDTH=%d&BBOX=%.4f,%.4f,%.4f,%.4f' ...
				'&STYLES=&SRS=EPSG:4326&VERSION=1.1.1'], handles.date, H, W, handles.region([1 3 2 4]));
	
	set(handles.figure1, 'pointer', 'watch')
	if (false)		% --> FORCE/NOT USE OF WGET 
		dest_fiche = 'lixogrr';
		if (ispc),		dos(['wget "' URL '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
		else			unix(['wget ''' URL ''' -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
		end
		set(handles.figure1, 'pointer', 'arrow')
		finfo = dir(dest_fiche);
		if (finfo.bytes < 1024)					% Delete the file anyway because it exists but is empty
			str = 'Failed to download image: ';
			if (min(W,H) > 5000)	str = [str '(possibly because it was too big) '];		end
			errordlg([str URL],'Error')
			builtin('delete',dest_fiche);
			return
		end
		img = gdalread(dest_fiche,'-U');
		if (size(img,3) == 4),		img(:,:,4) = [];	end		% We don't deal yet with alpha-maps in images
		builtin('delete',dest_fiche);
	else
		img = gdalread(URL,'-U');
		set(handles.figure1, 'pointer', 'arrow')
		if (isempty(img))
			str = 'Failed to download image: ';
			if (min(W,H) > 5000)	str = [str '(possibly because it was too big) '];		end
			errordlg([str URL],'Error')
			return
		end
	end

	struc.X = handles.region(1:2);		struc.Y = handles.region(3:4);
	struc.geog = 1;						struc.name = 'WMS image';
	struc.head = [handles.region 0 125 0 res/3600 res/3600];
	mirone(img, struc)

% ------------------------------------------------------------------------------------------
function wms_tool_LayoutFcn(h1)

set(h1, 'Position',[500 450 361 171],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','WMS give-me',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[91 100 51 22],...
'BackgroundColor',[1 1 1],...
'Call',@wms_tool_uiCB,...
'String','8',...
'Style','edit',...
'Tag','edit_res');

uicontrol('Parent',h1, 'Position',[0 103 90 15],...
'HorizontalAlignment','right', 'String','Image resolution',...
'Style','text');

uicontrol('Parent',h1, 'Position',[183 104 65 15],...
'HorizontalAlignment','right', 'String','Image date',...
'Style','text');

uicontrol('Parent',h1, 'Position',[250 100 101 22],...
'BackgroundColor',[1 1 1],...
'Call',@wms_tool_uiCB,...
'String','',...
'Style','edit',...
'Tooltip','Date in dd-mm-yyyy',...
'Tag','edit_date');

uicontrol('Parent',h1, 'Position',[11 50 201 21],...
'BackgroundColor',[1 1 1],...
'Call',@wms_tool_uiCB,...
'String','Daily planet (MODIS from JPL)',...
'Style','popupmenu',...
'Tooltip','Currently no other choice',...
'Value',1,...
'Tag','popup_servers');

uicontrol('Parent',h1, 'Position',[11 71 90 15],...
'HorizontalAlignment','left', 'String','WMS server',...
'Style','text');

uicontrol('Parent',h1, 'Position',[30 143 310 18],...
'FontAngle','italic',...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','bold',...
'String','Get referenced images from Web Map Server(s)',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 13 55 15],...
'Call',@wms_tool_uiCB,...
'String','proxy?',...
'Style','checkbox',...
'Tooltip','Check if you are Proxy-dependent',...
'Tag','check_proxy');

uicontrol('Parent',h1, 'Position',[66 11 91 20],...
'BackgroundColor',[1 1 1],...
'Call',@wms_tool_uiCB,...
'Enable','off',...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tooltip','proxy adress here',...
'Tag','edit_proxy');

uicontrol('Parent',h1, 'Position',[190 11 51 20],...
'BackgroundColor',[1 1 1],...
'Call',@wms_tool_uiCB,...
'Enable','off',...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tooltip','port here',...
'Tag','edit_port');

uicontrol('Parent',h1, 'Position',[161 14 28 15],...
'Enable','off',...
'HorizontalAlignment','right',...
'String','Port',...
'Style','text',...
'Tag','text_port');

uicontrol('Parent',h1, 'Position',[281 10 71 21],...
'Call',@wms_tool_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','OK',...
'Tag','push_OK');


function wms_tool_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

function varargout = travel(varargin)
% Tool to load a figure, reference it in Geogs and send it GE

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

	hObject = figure('Vis','off');
	travel_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center');

	handles.IPcmap = [];        % Initialize implanting image colormap to empty
	handles.ptLon = [];
	handles.ptLat = [];
	handles.no_file = false;	% writekml checks for these field
	handles.have_nans = false;
	handles.image_type = 2;

	mir_dirs = getappdata(0,'MIRONE_DIRS');
	if (~isempty(mir_dirs))
		handles.home_dir = mir_dirs.home_dir;		% Start in values
		handles.work_dir = mir_dirs.work_dir;
		handles.last_dir = mir_dirs.last_dir;
	else
		handles.home_dir = cd;		handles.work_dir = cd;		handles.last_dir = cd;
	end
	handles.path_tmp = [handles.home_dir filesep 'tmp' filesep];

	set(hObject,'Vis','on');
	if (nargout),       varargout{1} = hObject;    end
	guidata(handles.figure1, handles)

% ---------------------------------------------------------------------------------------
function push_loadImg_CB(hObject, handles)
	[FileName,PathName] = uigetfile({ ...
		'*.jpg', 'JPEG image (*.jpg)'; ...
		'*.tif', 'Tagged Image File (*.tif)'; ...
		'*.png', 'Portable Network Graphics(*.png)'; ...
		'*.bmp', 'Windows Bitmap (*.bmp)'; ...
		'*.pcx', 'Windows Paintbrush (*.pcx)'; ...
		'*.*', 'All Files (*.*)'}, ...
		'Select Image to Transplant');
	if isequal(FileName,0);     return;     end
	[PATH,FNAME,EXT] = fileparts([PathName FileName]);
	if (strcmpi(EXT,'.jpg') || strcmpi(EXT,'.jpeg') || strcmpi(EXT,'.bmp'))
		% I think those formats don't have colormaps and gdal is MUCH faster than imread
		I = gdalread([PathName FileName]);      handles.IPcmap = [];
	else
		[I,handles.IPcmap] = imread([PathName FileName]);
	end
	[m,n,k] = size(I);		handles.IPsize = [m n k];
	axes(handles.axes1);	handles.hImg = image(I);           axis off;

	% Compute implanting image aspect ratio and set axes 'PlotBoxAspectRatio' to it
	ip_aspect = m / n;
	set(handles.axes1,'PlotBoxAspectRatio',[1 ip_aspect 1])

	guidata(handles.figure1, handles)

% ---------------------------------------------------------------------------------------
function edit_ptLon_CB(hObject, handles)
	lon = get(hObject,'String');
	if (isnan(str2double(lon)))
		set(hObject,'String','')
		handles.ptLon = [];
		guidata(handles.figure1, handles),		return
	end
	val = test_dms(lon);
	lon = 0;
	if str2double(val{1}) > 0
		for i = 1:length(val),  lon = lon + str2double(val{i}) / (60^(i-1));    end
	else
		for i = 1:length(val),  lon = lon - abs(str2double(val{i})) / (60^(i-1));   end
	end
	handles.ptLon = lon;
	guidata(handles.figure1, handles)

% ---------------------------------------------------------------------------------------
function edit_ptLat_CB(hObject, handles)
	lat = get(hObject,'String');
	if (isnan(str2double(lat)))
		set(hObject,'String','')
		handles.ptLat = [];
		guidata(handles.figure1, handles),		return
	end
	val = test_dms(lat);
	lat = 0;
	if str2double(val{1}) > 0
		for i = 1:length(val),  lat = lat + str2double(val{i}) / (60^(i-1));    end
	else
		for i = 1:length(val),  lat = lat - abs(str2double(val{i})) / (60^(i-1));   end
	end
	handles.ptLat = lat;
	guidata(handles.figure1, handles)

% ---------------------------------------------------------------------------------------
function edit_height_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),		set(hObject,'String',''),	return,		end

% ---------------------------------------------------------------------------------------
function edit_width_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),		set(hObject,'String',''),	return,		end

% ---------------------------------------------------------------------------------------
function edit_picWidth_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),		set(hObject,'String',''),	return,		end

% ---------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...

	if (isempty(handles.ptLon))
		errordlg('I don''t know the Longitude of visiting point','Error'),	return
	end
	if (isempty(handles.ptLat))
		errordlg('I don''t know the Latitude of visiting point','Error'),	return
	end

	img_pbar = handles.IPsize(1) / handles.IPsize(2);	% PlotBoxAspectRatio

	width = str2double(get(handles.edit_picWidth,'Str'));
	if (isnan(width))		% Right, try size in loco
		width  = str2double(get(handles.edit_width,'Str'));
		height = str2double(get(handles.edit_height,'Str'));
		if (isnan(width) || isnan(height))
			errordlg('Singularities are not allowed. Neither picture width or size on the ground provided.','Error'), return
		end
	else
		height = round(width * img_pbar);
	end
	M_PER_DEG = 111195;		% Good enough approximation here of meter per degree in a spherical Earth
	latHeight = height / M_PER_DEG;
	lonWidth  = width / M_PER_DEG / cos(handles.ptLat * pi / 180);
	west  = handles.ptLon - lonWidth/2;			east  = handles.ptLon + lonWidth/2;
	south = handles.ptLat - latHeight/2;		north = handles.ptLat + latHeight/2;

	img = flipdim(get(handles.hImg,'CData'),1);		% Flipud because of the direct/reverse YDir shit
	handles.hImg = image([west east],[south north],img,'Parent',handles.axes1);
	set(handles.axes1,'PlotBoxAspectRatio',[1 img_pbar 1])
	handles.head = [west east south north];
	guidata(handles.figure1, handles)
	writekml(handles);


% ---------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
        delete(hObject);
	end

% --- Creates and returns a handle to the GUI figure. 
function travel_LayoutFcn(h1)

set(h1, 'Position',[520 430 316 370],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Mirone Travel',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'Visible','on');

axes('Parent',h1, ...
'Units','pixels',...
'Position',[4 139 307 228],...
'Box','on',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'Tag','axes1');

uicontrol('Parent',h1, 'Position',[100 108 51 22],...
'Call',@travel_uiCB,...
'Style','edit',...
'BackgroundColor',[1 1 1],...
'TooltipString','Width in meters that the picture have at its destination',...
'Tag','edit_picWidth');

uicontrol('Parent',h1, 'Position',[200 107 111 23],...
'Call',@travel_uiCB,...
'String','Load Image',...
'Tag','push_loadImg');

uicontrol('Parent',h1, 'Position',[45 81 161 14],...
'String','Coordinates at destination',...
'Style','text');

uicontrol('Parent',h1, 'Position',[3 63 25 14],...
'HorizontalAlignment','right',...
'String','Lon',...
'Style','text');

uicontrol('Parent',h1, 'Position',[30 59 81 22],...
'Call',@travel_uiCB,...
'Style','edit',...
'BackgroundColor',[1 1 1],...
'TooltipString','Longitude of landing point',...
'Tag','edit_ptLon');

uicontrol('Parent',h1, 'Position',[123 63 25 14],...
'HorizontalAlignment','right',...
'String','Lat',...
'Style','text');

uicontrol('Parent',h1, 'Position',[150 59 81 22],...
'Call',@travel_uiCB,...
'Style','edit',...
'BackgroundColor',[1 1 1],...
'TooltipString','Latitude  of landing point',...
'Tag','edit_ptLat');

uicontrol('Parent',h1, 'Position',[4 33 190 14],...
'String','OR - Size of picture in loco (meters)',...
'Style','text');

uicontrol('Parent',h1, 'Position',[3 13 35 14],...
'HorizontalAlignment','right',...
'String','Width',...
'Style','text');

uicontrol('Parent',h1, 'Position',[93 12 35 14],...
'HorizontalAlignment','right',...
'String','Height',...
'Style','text');

uicontrol('Parent',h1, 'Position',[40 9 51 22],...
'Call',@travel_uiCB,...
'Style','edit',...
'BackgroundColor',[1 1 1],...
'TooltipString','Width of the picture on ground in meters',...
'Tag','edit_width');

uicontrol('Parent',h1, 'Position',[130 9 51 22],...
'BackgroundColor',[1 1 1],...
'Call',@travel_uiCB,...
'Style','edit',...
'BackgroundColor',[1 1 1],...
'TooltipString','Height of the picture on ground in meters',...
'Tag','edit_height');

uicontrol('Parent',h1, 'Position',[3 113 95 14],...
'HorizontalAlignment','right',...
'String','Picture width (m)',...
'Style','text');

uicontrol('Parent',h1, 'Position',[200 8 111 23],...
'Call',@travel_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'String','Take me There',...
'Tag','push_OK');

function travel_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

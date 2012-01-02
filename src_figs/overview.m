function varargout = overview(varargin)
% This GUI builds an overview of grids readed by GDAL
% The selected grid is sub-sampled, during the
% reading stage, in order to obtain an overview of
% about 200 rows/columns. This way, arbitrarely
% large grids may be previewed, but the quality of
% the preview will may degrade with grid size.

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

if (nargin >= 4 && ischar(varargin{1}))
    gui_CB = str2func(varargin{1});
    feval(gui_CB,varargin{2:end})
else
    h = overview_OpeningFcn(varargin{:});
    if (nargout)    varargout{1} = h;   end
end

% --- Executes just before overview is made visible.
function hObject = overview_OpeningFcn(varargin)
	hObject = overview_LayoutFcn;
	handles = guihandles(hObject);
	move2side(hObject,'center');

	if (~isempty(varargin))
		handles.ForceInsitu = varargin{1}.ForceInsitu;
		handles.home_dir = varargin{1}.home_dir;
		handles.last_dir = varargin{1}.last_dir;
		handles.work_dir = varargin{1}.work_dir;
		handles.grdMaxSize = varargin{1}.grdMaxSize;
	else
		% Default values for those
		handles.home_dir = cd;
		handles.last_dir = cd;
		handles.work_dir = cd;
		handles.ForceInsitu = 0;
		handles.grdMaxSize = 1e15;
	end

	handles.path_tmp = [handles.home_dir filesep 'tmp' filesep];
	handles.hImg = [];
	handles.name_uncomp = [];
	handles.was_int16 = 0;
	handles.Nodata_int16 = [];

	load([handles.home_dir filesep 'data' filesep 'mirone_icons.mat'],'rectang_ico','info_ico','boreas_ico','help_ico','zoom_ico');

	h_toolbar = uitoolbar('parent',hObject,'Clipping', 'on', 'BusyAction','queue','HandleVisibility','on',...
		'Interruptible','on','Tag','FigureToolBar','Visible','on');
	uipushtool('parent',h_toolbar,'Click',@rectang_clickedCB,'Tag','rectang','cdata',rectang_ico,...
		'Tooltip','Draw Rectangle','Sep','on');
	uipushtool('parent',h_toolbar,'Click',@info_clickedCB,'Tag','info','Tooltip','Grid Info','cdata',info_ico);
	uipushtool('parent',h_toolbar,'Click',@all_clickedCB,'Tooltip','Swallow entire grid','cdata',boreas_ico);
	uitoggletool('parent',h_toolbar,'Click',@zoom_clickedCB,'Tag','zoom','Tooltip','Zoom','cdata',zoom_ico);
	uipushtool('parent',h_toolbar,'Click',@help_clickedCB,'Tag','help','Tooltip','Help','cdata',help_ico);

	logo = gdalread(['data' filesep 'logo.png']);
	image('parent',handles.axes1,'CData',flipdim(logo,1));

	set(hObject,'Visible','on','HandleVisibility','callback');
	guidata(hObject, handles);

% --------------------------------------------------------------------------------------------------
function rectang_clickedCB(obj,eventdata)
	handles = guidata(obj);     % get handles
	if (isempty(handles.hImg))		return,		end
try
	hl = findobj(handles.figure1,'Type','patch');
	if (~isempty(hl))		return,		end     % We already have a rectangle so don't need another.
	[p1,p2,hl] = rubberbandbox;
	% Round rectangle coords to the nearest node
	p1(1) = round((p1(1) - handles.head_orig(1)) / handles.head_orig(8)) * handles.head_orig(8) + handles.head_orig(1);
	p2(1) = round((p2(1) - handles.head_orig(1)) / handles.head_orig(8)) * handles.head_orig(8) + handles.head_orig(1);
	p1(2) = round((p1(2) - handles.head_orig(3)) / handles.head_orig(9)) * handles.head_orig(9) + handles.head_orig(3);
	p2(2) = round((p2(2) - handles.head_orig(3)) / handles.head_orig(9)) * handles.head_orig(9) + handles.head_orig(3);
	delete(hl);     % Given that we cannot change a line into a patch, the easy way is to delete and reborn it
	hl = patch('XData',[p1(1) p1(1) p2(1) p2(1) p1(1)],'YData',[p1(2) p2(2) p2(2) p1(2) p1(2)],...
		'FaceColor','none','EdgeColor','k','LineWidth',1);

	cmenuHand = uicontextmenu;
	set(hl, 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
	ui_edit_polygon(hl)    % Set edition functions
	uimenu(cmenuHand, 'Label', 'Rectangle limits', 'Separator','on', 'Callback', {@rectangle_limits,hl});
	uimenu(cmenuHand, 'Label', 'Crop Grid', 'Callback', {@CropaGrid,hl});
catch       % Don't know why but uisuspend sometimes breaks
	set(handles.figure1,'Pointer','arrow');
end

% --------------------------------------------------------------------------------------------------
function info_clickedCB(obj,eventdata)
	handles = guidata(obj);     % get handles
	if (isempty(handles.hImg))		return,		end
	w = getappdata(handles.axes1,'InfoMsg');
	if (~isempty(w))
		message_win('create',w, 'figname','Image Info');
	else
		message_win('create','Hoops, I lost the info', 'figname','Error');
	end

% --------------------------------------------------------------------------------------------------
function help_clickedCB(obj,eventdata)
str = sprintf(['This tool is meant to be used on very large grids.\n'...
	'The selected grid is sub-sampled, during the\n'...
	'reading stage, in order to obtain an overview of\n'...
	'about 200 rows/columns. This way, arbitrarily\n'...
	'large grids may be previewed, but the quality of\n'...
	'the preview will degrade with grid size.\n\n'...
	'Click on the rectangle icon to draw a rectangle\n'...
	'on the image. Double clicking the rectangle lets\n'...
	'you edit it. A precise control of the rect size is\n'...
	'available by a right-click and selecting "rectangle\n'...
	'-limits". Selecting "Crop Grid" extracts the region\n'...
	'inside the rectangle at the grid''s full resolution.']);
message_win('create',str, 'figname','Help')

% --------------------------------------------------------------------------------------------------
function all_clickedCB(obj,eventdata)
	handles = guidata(obj);			% get handles
	if (isempty(handles.hImg))		return,		end
	CropaGrid(obj,eventdata,[],'all')

% --------------------------------------------------------------------------------------------------
function zoom_clickedCB(obj,eventdata)
	handles = guidata(obj);     % get handles
	if (isempty(handles.hImg))
		set(obj,'State','off');    return
	end
	if strcmp(get(obj,'State'),'on')
		zoom_j('on');
	else
		zoom_j('off');
	end

% -----------------------------------------------------------------------------------------
function rectangle_limits(obj,eventdata,h)
	% Change the Rectangle's limits by asking it's corner coordinates
	x = get(h,'XData');   y = get(h,'YData');
	handles = guidata(obj);

	region = bg_region('with_limits',[x(1) x(3) y(1) y(3)]);
	if isempty(region),    return;  end     % User gave up
	x_min = region(1);      x_max = region(2);
	y_min = region(3);      y_max = region(4);

	% Make sure -R & -I are exactly compatible
	x_min = handles.head(1) + round( (x_min - handles.head(1)) / handles.head(8) ) * handles.head(8);
	x_max = handles.head(1) + round( (x_max - handles.head(1)) / handles.head(8) ) * handles.head(8);
	y_min = handles.head(3) + round( (y_min - handles.head(3)) / handles.head(9) ) * handles.head(9);
	y_max = handles.head(3) + round( (y_max - handles.head(3)) / handles.head(9) ) * handles.head(9);

	set(h, 'XData', [x_min,x_min,x_max,x_max,x_min], 'YData', [y_min,y_max,y_max,y_min,y_min]);

% --------------------------------------------------------------------
function FileOpen_DEM_CB(hObject, handles, opt)
% Files of the following formats are read (well directed) here
	switch opt
		case 'netCDF'
			str1 = {'*.grd;*.GRD;*.nc;*.NC', 'Grid files (*.grd,*.GRD,*.nc,*.NC)';'*.*', 'All Files (*.*)'};
		case 'ArcBinary'
			str1 = {'*.adf;*.ADF', 'Arc/Info grid (*.adf,*.ADF)'; '*.*', 'All Files (*.*)'};
		case 'DTED'
			str1 = {'*.dt0;*.DT0;*.dt1;*.DT1', 'DTED (*.dt0,*.DT0,*.dt1,*.DT1)'; '*.*', 'All Files (*.*)'};
		case 'ESRI_hdr'
			str1 = {'*.bil;*.BIL;', 'ESRI BIL (*.bil,*.BIL)'; '*.*', 'All Files (*.*)'};
		case (strcmp(opt,'ENVI'))
			str1 = {'*.img;*.IMG', 'ENVI (*.img,*.IMG)'; '*.*', 'All Files (*.*)'};
		case (strcmp(opt,'ERDAS'))
			str1 = {'*.img;*.IMG', 'Erdas (*.img,*.IMG)'; '*.*', 'All Files (*.*)'};
		case 'GTOPO30'
			str1 = {'*.dem;*.DEM', 'GTOPO30 DEM (*.dem,*.DEM)'; '*.*', 'All Files (*.*)'};
		case 'SRTM30'
			str1 = {'*.srtm;*.SRTM;*.srtm.gz', 'SRTM30 DEM (*.srtm,*.SRTM,*.srtm.gz)'; '*.*', 'All Files (*.*)'};
		case 'GeoTiff_DEM'
			str1 = {'*.tif;*.TIF;*.tiff;*.TIFF', 'GeoTiff DEM(*.tif,*.tiff,*.TIF,*.TIFF)'; '*.*', 'All Files (*.*)'};
		case {'MOLA_lbl' 'MOLA'}
			str1 = {'*.img;*.IMG', 'MOLA DEM (*.img,*.IMG)'; '*.*', 'All Files (*.*)'};
		case 'SDTS'
			str1 = {'*catd.ddf;*CATD.DDF', 'USGS SDTS DEM (*catd.ddf,*CATD.DDF)'; '*.*', 'All Files (*.*)'};
		otherwise
			str1 = {'*.grd;*.nc;*.tif;*.tiff;*.jpg;*.jp2;*.png;*.hdf', ...
					'Files (*.grd,*.nc,*.tif,*.tiff,*.jpg,*.jp2,*.png,*.hdf)'; '*.*', 'All Files (*.*)'};
			opt = '';
	end

	handles.name_uncomp = [];
	str2 = ['Select ' opt ' File'];
	[FileName,PathName,handles] = put_or_get_file(handles,str1,str2,'get');
	if isequal(FileName,0)		return,		end	
	fname = [PathName FileName];
	fname_t = fname;
	[PATH,FNAME,EXT] = fileparts(fname);

	if (strcmp(opt, 'SRTM30'))
		[name_hdr,comp_type] = write_esri_hdr(fname,'SRTM30');
		if (~isempty(comp_type))
			fname_t = decompress(fname,'warn');			% Name with compression ext removed
			if (isempty(fname_t))		return,		end	% Error message already issued.
			fname = fname_t;
			handles.name_uncomp = fname;				% We need a copy of the decompressed file name for removing
		end
	elseif (strcmpi(EXT,'.zip'))
		fname_t = ['/vsizip/' PathName FileName filesep FileName(1:end-4)];
	elseif (strcmpi(EXT,'.gz'))
		fname_t = ['/vsigzip/' PathName FileName];
	end
	
	if (strcmp(opt, 'MOLA_lbl'))		% This means that a .img MOLA file with a .lbl cannot be cmpressed
		fname = [PATH filesep FNAME '.lbl'];
		if (~exist(fname, 'file'))
			errordlg(['ERROR: Could not find format descriptive file: ' fname],'Error'),	return
		end
		fname_t = fname;
	end

	att =  gdalread(fname_t,'-M','-C');		% Use fname_t because fle may be /vsizip | /vsigzip
	opt_P = sprintf('-P%d',min(round(att.RasterXSize / 200), round(att.RasterYSize / 200)));

	[Z, X, Y, srsWKT, handles, att] = read_grid(handles, fname, 'OVR', opt_P);
	if (~strcmp(att.Band(1).DataType,'Byte'))
		if (~isa(Z,'single')),		Z = single(Z);		end
		if ( ~isempty(att.Band(1).NoDataValue) && ~isnan(att.Band(1).NoDataValue) && att.Band(1).NoDataValue ~= 0 )
			% Do this because some formats (e.g MOLA PDS v3) are so dumb that they declare a NoDataValue
			% and than don't use it !!!!!!
			if (att.Band(1).NoDataValue < 0)
				ind = (Z <= single(att.Band(1).NoDataValue));
			else
				ind = (Z >= single(att.Band(1).NoDataValue));
			end
			if (any(ind(:)))
				Z(ind) = NaN;		handles.have_nans = 1;
			end
		end
		Z = scaleto8(Z);
	end
	handles.head_orig = handles.head;

	handles.hImg = image(X,Y,Z);
	set(handles.axes1,'YDir','normal')
	Resize1(handles.axes1, handles.hImg)
	handles.fname = fname;
	grid_info(handles,att,'gdal');				% Construct a info message and save proj (if ...)
	set(handles.figure1,'pointer','arrow')
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function CropaGrid(obj,eventdata,h,opt)
% Select the region to be extracted from grid.
% If a fourth argument (OPT) is transmited, then the entire grid is read in
	handles = guidata(obj);     % get handles

	if (nargin == 3)
		x = get(h,'XData');		y = get(h,'YData');
		if (min(x) < handles.head_orig(1) || max(x) > handles.head_orig(2) || ...
				min(y) < handles.head_orig(3) || max(y) > handles.head_orig(4))
			errordlg('Selected region is partially outside the grid domain','Error')
			return
		end
		opt_R = sprintf('-R%.15f/%.15f/%.15f/%.15f', x(1),x(4),y(1),y(2));
	else
		opt_R = ' ';
	end

	[Z, X, Y, srsWKT, handles] = read_grid(handles, handles.fname, 'OVR', opt_R);

	tmp.head = handles.head;	tmp.X = X;		tmp.Y = Y;
	if (~isempty(srsWKT))		tmp.srsWKT = srsWKT;	end
	tmp.name = 'Subregion from Overview';
	tmp.was_int16 = handles.was_int16;
	tmp.Nodata_int16 = handles.Nodata_int16;
	mirone(Z,tmp);

%--------------------------------------------
function Resize1(axHandle, imHandle)
% Resize figure containing a single axes object with a single image.

set(axHandle,'Units','normalized','Position',[0 0 1 1]) % Don't realy understand why, but I need this

imageWidth  = size(get(imHandle, 'CData'), 2);
imageHeight = size(get(imHandle, 'CData'), 1);

% Don't try to handle the degenerate case.
if (imageWidth * imageHeight == 0),    return;  end

% What are the screen dimensions
screenSize = get(0, 'ScreenSize');      screenWidth = screenSize(3);    screenHeight = screenSize(4);
if ((screenWidth <= 1) || (screenHeight <= 1))
    screenWidth = Inf;    screenHeight = Inf;
end

% For small images, compute the minimum side as 60% of largest of the screen dimensions
% Except in the case of croped images, where 512 is enough for the pushbuttons
% LeastImageSide = fix(max([screenWidth screenHeight] * 0.6));
LeastImageSide = 400;

if (imageWidth < LeastImageSide && imageHeight < LeastImageSide)   % Augment very small images
    while (imageWidth < LeastImageSide && imageHeight < LeastImageSide)
        imageWidth = imageWidth*1.05;  imageHeight = imageHeight*1.05;
    end
    imageWidth = fix(imageWidth);   imageHeight = fix(imageHeight);
    % Large aspect ratio figures may still need to have their size adjusted
    if (imageWidth < 512)
        while (imageWidth < LeastImageSide && imageHeight < screenHeight-50)
            imageWidth = imageWidth*1.05;  imageHeight = imageHeight*1.05;
        end
    end
end

axUnits = get(axHandle, 'Units');
set(axHandle, 'Units', 'pixels');
axPos = get(axHandle, 'Position');

figHandle = get(axHandle, 'Parent');
figUnits = get(figHandle, 'Units');
rootUnits = get(0, 'Units');
set(figHandle, 'Units', 'pixels');
set(0, 'Units', 'pixels');

figLeftBorder = 10;  % assume left figure decorations are 10 pixels
figRightBorder = 10;
figBottomBorder = 30;
figTopBorder = 50;
figTopBorder = figTopBorder + 30;  % scribe hack

minFigWidth = 128; % don't try to display a figure smaller than this.
minFigHeight = 128;

% What are the gutter sizes?
figPos = get(figHandle, 'Position');
gutterLeft = max(axPos(1) - 1, 0);

scale = 100;    done = 0;
defAxesPos = get(0,'DefaultAxesPosition');
nonzeroGutters = (gutterLeft > 0);
while (~done)
    if (nonzeroGutters)
        gutterWidth = round((1 - defAxesPos(3)) * imageWidth / defAxesPos(3));
        gutterHeight = round((1 - defAxesPos(4)) * imageHeight / defAxesPos(4));
        newFigWidth = imageWidth + gutterWidth;
        newFigHeight = imageHeight + gutterHeight;
    else
        newFigWidth = imageWidth;        newFigHeight = imageHeight;
    end
    if (((newFigWidth + figLeftBorder + figRightBorder) > screenWidth) || ...
                ((newFigHeight + figBottomBorder + figTopBorder) > screenHeight))
        scale = 9 * scale / 10;
        imageWidth  = round(imageWidth * scale / 100);
        imageHeight = round(imageHeight * scale / 100);
    else
        done = 1;
    end
end

% ---------------------------------------------
h_Xlabel = get(axHandle,'Xlabel');  h_Ylabel = get(axHandle,'Ylabel');
units_save = get(h_Xlabel,'units');
set(h_Xlabel,'units','pixels');     set(h_Ylabel,'units','pixels');
Xlabel_pos = get(h_Xlabel,'pos');   Ylabel_pos = get(h_Ylabel,'pos');

if (abs(Ylabel_pos(1)) < 20)    % Stupid hack, but there is a bug somewhere
    Ylabel_pos(1) = 25;
end

setappdata(axHandle,'Backup_LabelPos',[Xlabel_pos Ylabel_pos])

y_margin = abs(Xlabel_pos(2))+get(h_Xlabel,'Margin') + 10;  % Devera conter a altura em pixeis do Xlabel
x_margin = abs(Ylabel_pos(1))+get(h_Ylabel,'Margin') + 10;  % Devera conter a largura em pixeis do Ylabel
if strcmp(get(axHandle,'Visible'),'off')               % No Labels, give only a 20 pixels margin to account for Status bar
    x_margin = 0;   y_margin = 0;
end
set(h_Xlabel,'units',units_save);     set(h_Ylabel,'units',units_save);

newFigWidth = max(newFigWidth, minFigWidth) + x_margin;
newFigHeight = max(newFigHeight, minFigHeight) + y_margin;

figPos(1) = max(1, figPos(1) - floor((newFigWidth - figPos(3))/2));
figPos(2) = max(1, figPos(2) - floor((newFigHeight - figPos(4))/2));
figPos(3) = newFigWidth;
figPos(4) = newFigHeight;

% Figure out where to place the axes object in the resized figure
gutterWidth = figPos(3) - imageWidth;
gutterHeight = figPos(4) - imageHeight;
gutterLeft = floor(gutterWidth/2) + x_margin/2 - 1;
gutterBottom = floor(gutterHeight/2) + y_margin/2 - 1;

axPos(1) = gutterLeft*0.8 + 1;  axPos(2) = gutterBottom - 5;
axPos(3) = imageWidth;      axPos(4) = imageHeight;

% Force the window to be in the "north" position. 75 is the height of the blue Bar + ...
figPos(2) = screenHeight - figPos(4) - 75;
set(figHandle, 'Position', figPos);     set(axHandle, 'Position', axPos);

% Restore the units
drawnow;		% necessary to work around HG bug   -SLE
set(figHandle, 'Units', figUnits);
set(axHandle, 'Units', axUnits);
set(0, 'Units', rootUnits);

% --------------------------------------------------------------------
function figure1_DeleteFcn(hObject, eventdata)
% Before deleting the fig see if there was files left to remove
% (files that resulted from the uncompression)
	handles = guidata(hObject);
	if (~isempty(handles.name_uncomp))
		try
			delete(handles.name_uncomp);
		end
	end
	delete(hObject)

% --------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
if isequal(get(hObject,'CurrentKey'),'escape')
    % Before deleting the fig see if there was files left to remove
    % (files that resulted from the uncompression)
	handles = guidata(hObject);
	if (~isempty(handles.name_uncomp))
		try
			delete(handles.name_uncomp);
		end
	end
	delete(hObject)
end

% --------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function h1 = overview_LayoutFcn()

h1 = figure('Position',[520 560 244 240],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'IntegerHandle','off',...
'KeyPressFcn',@figure1_KeyPressFcn,...
'DeleteFcn',@figure1_DeleteFcn,...
'MenuBar','none',...
'Name','Overview',...
'NumberTitle','off',...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Colormap',jet(256),...
'Visible','off',...
'Tag','figure1');

h2 = axes('Parent',h1,...
'Units','pixels',...
'Position',[1 1 244 240],...
'Tag','axes1',...
'Visible','off');

h4 = get(h2,'xlabel');
set(h4,'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[0.497916666666667 -0.0979166666666667 1.00005459937205],...
'VerticalAlignment','cap',...
'HandleVisibility','off',...
'Visible','off');

h5 = get(h2,'ylabel');
set(h5,'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[-0.11875 0.49375 1.00005459937205],...
'Rotation',90,...
'VerticalAlignment','bottom',...
'HandleVisibility','off',...
'Visible','off');

h7 = uimenu('Parent',h1,'Label','File Open');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''netCDF'')',...
'Label','netCDF (GMT) grid');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''ArcBinary'')',...
'Label','Arc/Info Binary grid');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''DTED'')', 'Label','DTED');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''ENVI'')', 'Label','ENVI');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''ERDAS'')', 'Label','ERDAS');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''ESRI_hdr'')', 'Label','ESRI BIL');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''GTOPO30'')', 'Label','GTOPO30');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''GeoTiff_DEM'')', 'Label','GeoTiff DEM');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''MOLA_lbl'')', 'Label','MOLA (+ .lbl)');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''MOLA'')', 'Label','MOLA (v3 PDS)');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''SDTS'')', 'Label','SDTS');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''SRTM30'')', 'Label','SRTM30');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_CB'',gcbo,guidata(gcbo),''Other'')',...
'Sep','on',...
'Label','Others');

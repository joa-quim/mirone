function varargout = vitrinite(varargin)
% Helper window to assist carbons Age guessing black magic

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

% $Id: vitrinite.m 3608 2012-07-12 16:37:13Z j $

	hObject = figure('Tag','figure1','Visible','off');
	vitrinite_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(handles.figure1,'center')

	mir_dirs = getappdata(0,'MIRONE_DIRS');
	if (~isempty(mir_dirs))
		handles.home_dir = mir_dirs.home_dir;		% Start in values
		handles.work_dir = mir_dirs.work_dir;
		handles.last_dir = mir_dirs.last_dir;
	else
		handles.home_dir = cd;		handles.work_dir = cd;		handles.last_dir = cd;
	end

	% Import icons
	load([handles.home_dir filesep 'data' filesep 'mirone_icons.mat'],'Mfopen_ico','zoom_ico');

	hTB = uitoolbar('parent',handles.figure1, 'BusyAction','queue','HandleVisibility','on',...
	   'Interruptible','on','Tag','FigureToolBar','Visible','on');
	uipushtool('parent',hTB,'Click',@clicked_loadMarkersImg_CB, ...
	   'cdata',Mfopen_ico,'Tooltip','Open Marker images');
	uipushtool('parent',hTB,'Click',@clicked_loadSampleImg_CB, ...
	   'cdata',Mfopen_ico,'Tooltip','Open Sample image','Sep','on');
	uitoggletool('parent',hTB,'Click','zoom', ...
	   'cdata',zoom_ico,'Tooltip','Zooming on/off','Sep','on');

	handles.hImgMarker = [];
	handles.hImgSample = [];
	handles.countMark = 0;         % Counter of number of loaded image markers
	handles.tiePoints = [];
	handles.countTiePoints = 0;
	handles.haveStatusBar = 0;
	handles.IamCallibrated = 0;
	handles.geog = 0;

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% --------------------------------------------------------------------
function clicked_loadMarkersImg_CB(hObject, eventdata)
	handles = guidata(hObject);
	str = { '*.jpg', 'JPEG image (*.jpg)'; ...
			'*.png', 'Portable Network Graphics(*.png)'; ...
			'*.bmp', 'Windows Bitmap (*.bmp)'; ...
			'*.gif', 'GIF image (*.gif)'; ...
			'*.pcx', 'Windows Paintbrush (*.pcx)'; ...
			'*.ras', 'SUN rasterfile (*.ras)'; ...
			'*.tif', 'Tagged Image File (*.tif)'; ...
			'*.*', 'All Files (*.*)'};
	[FileName,PathName, handles] = put_or_get_file(handles,str,'Select image','get');
	if isequal(FileName,0),		return,		end

	fname = [PathName FileName];
	try             % Use a try because if the name was given via edit box it may be wrong
		[img,cmap] = imread(fname);
		handles.ImgMarker{handles.countMark+1} = img;    % Make a copy
		handles.countMark = handles.countMark + 1;
	catch
		errordlg(['Error: -> ' fname ' does not exist or is not a valid image file'],'Error')
		return
	end
	handles.head = [1 size(img,2) 1 size(img,1) 0 255 0 1 1];	% Fake a grid reg GMT header

	if (ndims(img) == 3)
		errordlg('Error: True color images are not allowed.','Error');        return
	end
	[m,n] = size(img);
	if (isempty(cmap))
		set(handles.figure1,'Colormap', gray(256));
	elseif (~isempty(cmap))
		set(handles.figure1,'Colormap', cmap);
	end

	% Compute image aspect ratio and set axes 'PlotBoxAspectRatio' to it
	aspect = m / n;
	handles.hImgMarker = image(img,'Parent',handles.axes2);
	handles.imgAspect{handles.countMark} = [1 aspect 1];    % Store image aspect ratio
	set(handles.axes2,'PlotBoxAspectRatio',[1 aspect 1],'Visible','off')

	% Update the Image Markers popup
	[PATH,FNAME,EXT] = fileparts(fname);
	if (handles.countMark == 1)
		set(handles.popup_markersImgs,'String',[FNAME EXT])
	else
		str = get(handles.popup_markersImgs,'String');
		if (iscell(str))
			str{end+1} = [FNAME EXT];
		else
			str = {str; [FNAME EXT]};
		end
		set(handles.popup_markersImgs,'String',str,'Value',handles.countMark)	% Set it the visible one
	end

	% Try to fish the reflectance factor from image name
	sep = strfind(FNAME,'_');
	if (~isempty(sep))
		reflectance = str2double(FNAME(sep(end)+1:end));
		if (~isnan(reflectance) && reflectance > 0.1 && reflectance < 10)
			handles.reflectance{handles.countMark} = reflectance;				% Yes, we got it
		else
			handles.reflectance{handles.countMark} = [];
		end
	else
		handles.reflectance{handles.countMark} = [];
	end

	% If we got a reflectance put it right away in the the editbox
	if (~isempty(handles.reflectance{handles.countMark}))
		set(handles.edit_reflectance,'String',reflectance)
	end

	handles.is_projected = 0;
	if (handles.countMark >= 1 && ~isempty(handles.hImgSample))
		if (~handles.haveStatusBar)
			createStatusBar(handles);            handles.haveStatusBar = 1;
		end
	end
	guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function clicked_loadSampleImg_CB(hObject, eventdata)
	handles = guidata(hObject);
	str = { '*.jpg', 'JPEG image (*.jpg)'; ...
			'*.png', 'Portable Network Graphics(*.png)'; ...
			'*.bmp', 'Windows Bitmap (*.bmp)'; ...
			'*.gif', 'GIF image (*.gif)'; ...
			'*.pcx', 'Windows Paintbrush (*.pcx)'; ...
			'*.ras', 'SUN rasterfile (*.ras)'; ...
			'*.tif', 'Tagged Image File (*.tif)'; ...
			'*.*', 'All Files (*.*)'};
	[FileName,PathName, handles] = put_or_get_file(handles,str,'Select image','get');
	if isequal(FileName,0),		return,		end

    fname = [PathName FileName];
	try             % Use a try because if the name was given via edit box it may be wrong
        [handles.ImgSample,cmap] = imread(fname);
	catch
        errordlg(['Error: -> ' fname ' does not exist or is not a valid image file'],'Error')
        return
	end
	handles.head = [1 size(handles.ImgSample,2) 1 size(handles.ImgSample,1) 0 255 0 1 1];	% Fake a grid reg GMT header

    if (ndims(handles.ImgSample) == 3)
        errordlg('Error: True color images are not allowed.','Error')
        return
    end
	[m,n] = size(handles.ImgSample);
	if (isempty(cmap))
        set(handles.figure1,'Colormap', gray(256));
	elseif (~isempty(cmap))
        set(handles.figure1,'Colormap', cmap);
	end
	
	% Compute image aspect ratio and set axes 'PlotBoxAspectRatio' to it
	aspect = m / n;
    if (isempty(handles.hImgSample))    % First time a sample image is loaded
	    handles.hImgSample = image(handles.ImgSample,'Parent',handles.axes1);
    else
        set(handles.hImgSample,'CData',handles.ImgSample)
    end
	set(handles.axes1,'PlotBoxAspectRatio',[1 aspect 1],'Visible','off')
	
	handles.is_projected = 0;
    if (~handles.IamCallibrated)
		if (handles.countMark >= 1 && ~handles.haveStatusBar)
			createStatusBar(handles);            handles.haveStatusBar = 1;
		end
    else        % A previous callibrated image already exists. Need to callibrate this one
       push_callibrate_CB(handles.push_callibrate, [], handles) 
    end
	guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function popup_markersImgs_CB(hObject, handles)
	if (handles.countMark < 2)		return,		end      % Too soon
	set(handles.edit_reflectance,'String',handles.reflectance{get(hObject,'Value')})
	resetImg(handles,get(hObject,'Value'))

% --------------------------------------------------------------------
function resetImg(handles, number)
% Update image marker accordingly to what was selected in the popup
	set(handles.hImgMarker,'CData',handles.ImgMarker{number})
	set(handles.axes2,'PlotBoxAspectRatio',handles.imgAspect{number})    

% --------------------------------------------------------------------
function push_callibrate_CB(hObject, handles)
% ... 
	if (isempty(handles.hImgSample))
		errordlg('Calibrate what? Your eyes? Load the sample image.','ERROR');return
	end
	if (handles.countTiePoints < 2)					% Not yet possible to calibrate
		errordlg('No I won''t. I need at least two tie points to do that.','ERROR');return
	end

	x = sort(cell2mat(handles.tiePoints));			% Marker pixel values
	Y = sort(cell2mat(handles.reflectance));		% Marker reflectance values
	%     x = [0.0 x];								% Don't let extrapolation pass to negative values
	%     Y = [0.15 Y];
	Y_bak = Y;

	yi = single(interp1(x,Y,0:255,'linear','extrap'));
	Z = yi(handles.ImgSample);
	X = 1:size(handles.ImgSample,2);
	Y = 1:size(handles.ImgSample,1);
	head = [1 X(end) 1 Y(end) 0 255 0 1 1];
	setappdata(handles.figure1,'dem_z',Z);  setappdata(handles.figure1,'dem_x',X);
	setappdata(handles.figure1,'dem_y',Y);
	set(handles.push_getAvgReflec,'Enable','on')
	handles.IamCallibrated = 1;
	handles.callibCurv = yi;                % Save the callibration curve

	% Display the callibration curve
	try             % Use a try because on the first time they don't yet exist
		delete(handles.callLine)
		delete(handles.callPts)
		%delete(handles.errorEnvelop)
	end
	%     x = double(x);      yi = double(yi);
	%     y1 = [Y_bak+[0, yi(round( cat(2,handles.sigma{:}) ))]];
	%     y2 = [Y_bak-[0, yi(round( cat(2,handles.sigma{:}) ))]];
	%     y = [y1 y2(end:-1:1)];
	%     handles.errorEnvelop = patch('XData',[x x(end:-1:1)],'YData', y,'FaceColor',[.9 .9 .9],'Parent',handles.axes3);
	handles.callLine = line('XData',0:255,'YData',yi,'Parent',handles.axes3);
	handles.callPts = line('XData',x,'YData',Y_bak,'Parent',handles.axes3,'LineStyle','none',...
		'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',7);

	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function push_getTiePoint_CB(hObject, handles)
% Get a tie point from the current Marker image
	[x,y,but]  = click_e_point(1,'crosshair');
	if (but ~= 1),   return;     end
	params.Point = [x y];	params.Tolerance = 10;		params.Connect = 4;

	img = get(handles.hImgMarker,'CData');					% Get the image
	dumb = repmat(img,[1 1 3]);				% Because something strange in cvlib_mex for gray case.
	[dumb,mask] = cvlib_mex('floodfill',dumb,params);
	clear dumb;

	% TESTAR SE TENHO REFLECTANCE PARA ESTE PONTO
	tiePoint = round(mean2(img(mask)));
	whichMarker = get(handles.popup_markersImgs,'Value');	% We need to know which tie point is this
	if (numel(handles.tiePoints) < whichMarker)
		handles.tiePoints{whichMarker} = tiePoint;
		handles.countTiePoints = handles.countTiePoints + 1;
	else
		if (isempty(handles.tiePoints{whichMarker}))
			handles.tiePoints{whichMarker} = tiePoint;
			handles.countTiePoints = handles.countTiePoints + 1;
		else
			% New value for an already existing tiePoint. Do not increment countTiePoints
			handles.tiePoints{whichMarker} = tiePoint;
		end
	end
	%	handles.sigma{whichMarker} = std(double(img(mask)));
	guidata(handles.figure1,handles)

%--------------------------------------------------------------------------
function push_getAvgReflec_CB(hObject, handles)
% Get the average reflectance of the clicked shape
	if (~handles.IamCallibrated)
		errordlg('You need to callibrate the sample image first.','ERROR'); return
	end
	[x,y,but]  = click_e_point(1,'crosshair');
	if (but ~= 1),   return;     end
	params.Point = [x y];    params.Tolerance = 5;    params.Connect = 8;
	hAx = get(handles.figure1,'CurrentAxes');
	img = get(findobj(hAx,'Type','Image'),'CData');                   % Get the image
	if (ndims(img) == 2)	dumb = repmat(img,[1 1 3]);
	else					dumb = img;
	end

	[dumb,mask] = cvlib_mex('floodfill',dumb,params);
	clear dumb;

	avgPix = round(mean2(img(mask)));
	avgReflect = double(handles.callibCurv(avgPix));
	hText = text(x,y,sprintf('%.2f',avgReflect),'Fontsize',8,'Parent',hAx,'HorizontalAlignment','center');

	% Set a uicontext with the "Deleting" option
	cmenuHand = uicontextmenu;		set(hText, 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Delete', 'Call', 'delete(gco)');

%--------------------------------------------------------------------------
function edit_reflectance_CB(hObject, handles)
% Will have code when we let enter the reflectance value here
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0.1 || xx > 10)
		set(hObject, ''),		return
	end
	whichMarker = get(handles.popup_markersImgs,'Value');	% Need to know which tie point this refers to
	handles.reflectance{whichMarker} = xx;
	guidata(handles.figure1, handles)

%--------------------------------------------------------------------------
function createStatusBar(handles)
% simulates a box at the bottom of the figure
	figPos = get(handles.figure1,'Pos');
	H = 22;
	sbPos(1) = 1;               sbPos(2) = 2;
	sbPos(3) = figPos(3)-2;     sbPos(4) = H-1;
	h = axes('Parent',handles.figure1,'Box','off','Visible','off','Tag','sbAxes','Units','Pixels',...
        'Position',sbPos,'XLim',[0 sbPos(3)],'YLim',[0 H-1]);
	hFieldFrame = createframe(h,[1 (figPos(3) - 1)],H);
	setappdata(handles.figure1,'CoordsStBar',[h hFieldFrame]);  % Save it for use in ...
	set(hFieldFrame,'Visible','on')
	set(h,'HandleVisibility','off')
	pixval_stsbar(handles.figure1);

%--------------------------------------------------------------------------
function hFrame = createframe(ah,fieldPos,H)
% Creates a virtual panel surrounding the field starting at fieldPos(1) and
% ending end fieldPos(2) pixels. ah is the sb's handle (axes).
% It returns a handle array designating the frame.
	
	from = fieldPos(1);     to = fieldPos(2);
	% col = rgb2hsv(get(fh,'Color'));       % fh was the figure's handle
	% lightColor = col;   lightColor(2) = 0.5*lightColor(2);  lightColor(3) = 0.9; lightColor = hsv2rgb(lightColor);
	% darkColor = col;    darkColor(3) = 0.4;  darkColor = hsv2rgb(darkColor);
	% This is the result of the above. I just don't want the extra burden of compiling those routines
	lightColor = [0.9 0.89150943396226 0.87452830188679];
	darkColor  = [0.4 0.39245283018868 0.37735849056604];
	
	hFrame(1) = line([from to],[H-2 H-2],'Color',darkColor,'Visible','off','Tag','Sts_T','parent',ah);    % Top line
	hFrame(2) = line([from from],[1 H-2],'Color',darkColor,'Visible','off','Tag','Sts_L','parent',ah);    % Left line
	hFrame(3) = line([from+1 to-1],[1 1],'Color',lightColor,'Visible','off','Tag','Sts_B','parent',ah);   % Bottom line
	hFrame(4) = line([to-1 to-1],[1 H-2],'Color',lightColor,'Visible','off','Tag','Sts_R','parent',ah);   % Right line

% -------------------------------------------------------------------------------------
function y = mean2(x)
%MEAN2 Compute mean of matrix elements.
	y = sum(x(:)) / numel(x);

% -------------------------------------------------------------------------------------
function y = median2(x)
%MEDIAN2 Compute median of matrix elements.
	y = median(x(:));


% --- Creates and returns a handle to the GUI figure. 
function vitrinite_LayoutFcn(h1)

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Vitrinite',...
'NumberTitle','off',...
'Position',[520 268 780 532],...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','Callback',...
'Tag','figure1');

h2 = axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Position',[400 183 371 331],...
'XTick', [], 'YTick', [],...
'Tag','axes1');

h3 = get(h2,'title');

set(h3,'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[0.5 1.01963746223565 1.00005459937205],...
'VerticalAlignment','bottom',...
'HandleVisibility','off');

h7 = axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Position',[10 183 371 331],...
'XTick', [], 'YTick', [],...
'Tag','axes2');

h8 = get(h7,'title');

set(h8,'Parent',h7,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[0.5 1.01963746223565 1.00005459937205],...
'VerticalAlignment','bottom',...
'HandleVisibility','off');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'popup_markersImgs_CB'},...
'Position',[10 134 241 22],...
'Style','popupmenu',...
'String',{''},...
'TooltipString','Loaded Marker Images',...
'Value',1,...
'Tag','popup_markersImgs');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_callibrate_CB'},...
'FontSize',9,...
'Position',[520 132 91 23],...
'String','Calibrate',...
'Tag','push_callibrate');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_getTiePoint_CB'},...
'FontSize',9,...
'Position',[410 132 91 23],...
'String','Get Tie Point',...
'Tag','push_getTiePoint');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_getAvgReflec_CB'},...
'Enable','inactive',...
'FontSize',9,...
'Position',[634 132 137 24],...
'String','Get average reflectance',...
'Tag','push_getAvgReflec');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_reflectance_CB'},...
'Position',[308 134 71 21],...
'Style','edit',...
'TooltipString','Reflectance of current marker.',...
'Tag','edit_reflectance');

uicontrol('Parent',h1,...
'FontSize',9,...
'Position',[305 157 78 16],...
'String','Reflectance',...
'Style','text',...
'Tag','text3');

h18 = axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[127.5 4 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Position',[20 41 741 81],...
'XColor',get(0,'defaultaxesXColor'),...
'XLim',[0 255],...
'XLimMode','manual',...
'YColor',get(0,'defaultaxesYColor'),...
'YLim',[0 8],...
'YLimMode','manual',...
'HandleVisibility','off',...
'Tag','axes3');

h19 = get(h18,'title');

set(h19,'Parent',h18,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[127.5 8.64197530864197 1.00005459937205],...
'VerticalAlignment','bottom',...
'HandleVisibility','off');

uicontrol('Parent',h1,...
'FontSize',9,...
'Position',[10 158 131 16],...
'String','Current Marker Image',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[23 101 108 16],...
'String','Callibration Curve',...
'Style','text');

uicontrol('Parent',h1,...
'FontSize',9,...
'Position',[138 515 131 16],...
'String','Marker Image',...
'Style','text');

uicontrol('Parent',h1,...
'FontSize',9,...
'Position',[510 515 131 16],...
'String','Sample Image',...
'Style','text');

function main_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

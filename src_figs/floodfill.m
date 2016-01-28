function varargout = floodfill(varargin)
% Helper Fig to do color segmentation or painting like the magick wand

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

% $Id: floodfill.m 4457 2014-05-06 01:40:33Z j $

	if (isempty(varargin))
		errordlg('FLOODFILL: wrong number of arguments.','Error'),	return
	end

	hObject = figure('Tag','figure1','Visible','off');
	floodfill_LayoutFcn(hObject);
	handles = guihandles(hObject);

    handles.hCallingFig = varargin{1};
    handMir = guidata(handles.hCallingFig);
	if (handMir.no_file)
		errordlg('You didn''t even load a file. What are you expecting then?','ERROR')
		delete(hObject);    return
	end

	% Import icons
	load ([handMir.path_data 'mirone_icons.mat'],'lapis_ico','trincha_ico','pipeta_ico','balde_ico',...
		'circ_ico','rectang_ico','ellipse_ico');

	h_toolbar = uitoolbar('parent',hObject, 'BusyAction','queue','HandleVisibility','on',...
		'Interruptible','on','Tag','FigureToolBar','Visible','on');
	uipushtool('parent',h_toolbar,'Click',{@line_clickedcallback,'pencil'}, ...
		'cdata',lapis_ico,'TooltipString','Pencil');
	uipushtool('parent',h_toolbar,'Click',{@line_clickedcallback,'paintbrush'}, ...
		'cdata',trincha_ico,'TooltipString','Paintbrush');
	uipushtool('parent',h_toolbar,'Click',@pipeta_clickedcallback,'cdata',pipeta_ico,'TooltipString','Color Picker');
	uipushtool('parent',h_toolbar,'Click',@flood_clickedcallback,'cdata',balde_ico,'TooltipString','Floodfill');
	uipushtool('parent',h_toolbar,'Click',{@shape_clickedcallback,'circ'},'cdata',circ_ico,'TooltipString','Circle');
	uipushtool('parent',h_toolbar,'Click',{@shape_clickedcallback,'rect'},'cdata',rectang_ico,'TooltipString','Rectangle');
	uipushtool('parent',h_toolbar,'Click',{@shape_clickedcallback,'ellipse'},'cdata',ellipse_ico,'TooltipString','Ellipse');

	% Try to position this figure glued to the right of calling figure
	posThis = get(hObject,'Pos');
	posParent = get(handles.hCallingFig,'Pos');
	ecran = get(0,'ScreenSize');
	xLL = posParent(1) + posParent(3) + 6;
	xLR = xLL + posThis(3);
	if (xLR > ecran(3))         % If figure is partially out, bring totally into screen
		xLL = ecran(3) - posThis(3);
	end
	yLL = (posParent(2) + posParent(4)/2) - posThis(4) / 2;
	set(hObject,'Pos',[xLL yLL posThis(3:4)])

	handles.hCallingAxes = get(handles.hCallingFig,'CurrentAxes');
	handles.hImage = handMir.hImg;
	img = get(handles.hImage,'CData');
	handles.origFig = img;			% Make a copy of the original image
	handles.must_restore_img = false;
	handles.img_was_changed  = false;
	if (ndims(img) == 2)			% Ghrr, we need to convert to RGB because the floodfill here works bad with indexed
		pal = get(handles.hCallingFig,'Colormap');
		img = ind2rgb8(img, pal);
		set(handles.hImage,'CData', img);
		handles.must_restore_img = true;	% At the end (in closerequest), restore original image
	end
	handles.imgSize = size(img);

	handles.IAmAMirone = getappdata(handles.hCallingFig,'IAmAMirone');
	if (~isempty(handles.IAmAMirone))
		handlesMir = guidata(handles.hCallingFig);
		handles.head = handlesMir.head;
		handles.DefLineThick = handlesMir.DefLineThick;
		handles.DefLineColor = handlesMir.DefLineColor;
		handles.image_type   = handlesMir.image_type;
	else
		handles.IAmAMirone = 0;
		% Build a GMT type header info
		handles.head    = [get(handles.hCallingAxes,'xlim') get(handles.hCallingAxes,'ylim') 0 0 1];
		handles.head(8) = (handles.head(2)-handles.head(1)) / size(img,2);
		handles.head(9) = (handles.head(4)-handles.head(3)) / size(img,1);

		handles.DefLineThick = 1;
		handles.DefLineColor = [0 0 0];
	end

	% Initialize some vars
	handles.connect = 4;
	handles.tol = 20;
	handles.randColor = 1;
	handles.fillColor = [255 255 255];
	handles.useDilation = 1;            % Dilate mask before using it in picking shapes
	handles.colorSegment = 1;           % Default to do color segmentation
	%handles.colorModel = 'YCrCb';       % When color segmentation, default to this color space
	handles.colorModel = [];            % Do color segmentation in RGB
	handles.minPts = 50;                % When digitizing, dont create polygs with less than this pts
	handles.udCount = 1;                % Unique incremental identifier to set in UserData of polygons
	handles.single_poly = 0;            % IF == 1 -> all detected polygons are drawn in a single multi-polygon
	handles.bg_color = 0;               % Background color (or color number in cmap)
	handles.lineWidth = 3;
	handles.elemSquare = 1;

	set(handles.listbox_lineWidth,'String',1:99,'Val',2)
	set(handles.slider_tolerance,'Value',handles.tol)

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, [handles.text_Paint handles.text_DS])
	%------------- END Pro look (3D) -----------------------------------------------------

	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end
	set(hObject,'Visible','on');

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handMir.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handMir.figure1,'dependentFigs',plugedWin);

% --------------------------------------------------------------------
function line_clickedcallback(hObject, eventdata, opt)
    handles = guidata(hObject);     % get handles
    figure(handles.hCallingFig)
    
    state = uisuspend_j(handles.hCallingFig);        % Remember initial figure state
    if (strcmp(opt,'pencil'))
        set(handles.hCallingFig,'Pointer', 'custom','PointerShapeCData',getPointer('pencil'),'PointerShapeHotSpot',[14 2])
    else
        set(handles.hCallingFig,'Pointer', 'custom','PointerShapeCData',getPointer('brush'),'PointerShapeHotSpot',[14 2])
    end
    w = waitforbuttonpress;
    if (w == 0),    paintFirstButtonDown(handles,state,opt)       % A mouse click
    else            set(handles.hCallingFig,'Pointer', 'arrow');
    end

% -------------------
function paintFirstButtonDown(handles,state,opt)
	if (strcmp(opt,'pencil')),      lineThick = 1;
	else                            lineThick = handles.lineWidth;
	end
	lineType = 8;       % Default to 8 connectivity
	if (get(handles.checkbox_AA,'Value')),      lineType = 16;      end
	pt = get(handles.hCallingAxes, 'CurrentPoint');
	setappdata(handles.figure1,'prev_pt',pt(1,1:2))
	set(handles.hCallingFig,'WindowButtonMotionFcn',{@wbm_line,handles,lineThick,lineType},...
		'WindowButtonDownFcn',{@wbd_paint,handles.hCallingFig,state});
	if (~handles.elemSquare && lineThick > 1)       % Use a round element
		set(handles.hCallingFig,'WindowButtonMotionFcn',{@wbm_circ,handles,lineThick,lineType})
	end

	handles.img_was_changed = true;
	guidata(handles.figure1, handles)

% -------------------
function wbm_line(obj,eventdata,handles,lineThick,lineType)
    % Draw the line using a square element
    pt = get(handles.hCallingAxes, 'CurrentPoint');
    prev_pt = getappdata(handles.figure1,'prev_pt');
    setappdata(handles.figure1,'prev_pt',pt(1,1:2))
    [x,y] = getpixcoords(handles,[prev_pt(1) pt(1,1)],[prev_pt(2) pt(1,2)]);
    x = round(x);       y = round(y);
    if (~insideRect(handles,x(2),y(2))),      return;     end
    win_dx = abs(diff(x)) + 4;         win_dy = abs(diff(y)) + 4;
    
    % Notice that the (2) index denotes the current point
    win_left = max(x(2)-win_dx,1);      win_right = min(x(2)+win_dx,handles.imgSize(2));
    win_top  = max(y(2)-win_dy,1);      win_bot   = min(y(2)+win_dy,handles.imgSize(1));
    if (win_top > win_bot)  % Buble sort
        fds = win_top;        win_top = win_bot;        win_bot = fds;
    end
       
    r = win_top:win_bot;    c = win_left:win_right;
    img = get(handles.hImage,'CData');
    img_s = img(r,c,:);     % Extract a sub-image
    x = x - c(1);           % We need the PTs coords relative to the subimage
    y = y - r(1);
    cvlib_mex('line',img_s,[x(1) y(1)],[x(2) y(2)],handles.fillColor,lineThick,lineType)
    img(r,c,:) = img_s;     % Reimplant the painted sub-image
    set(handles.hImage,'CData',img);

% -------------------
function wbm_circ(obj,eventdata,handles,lineThick,lineType)
    % Draw the line using a circular element
    pt = get(handles.hCallingAxes, 'CurrentPoint');
    [x,y] = getpixcoords(handles,pt(1,1),pt(1,2));
    if (~insideRect(handles,x,y)),      return;     end
    x = round(x);       y = round(y);
    win_left = max(x-lineThick-1,1);      win_right = min(x+lineThick+1,handles.imgSize(2));
    win_top  = max(y-lineThick-1,1);      win_bot   = min(y+lineThick+1,handles.imgSize(1));
    if (win_top > win_bot)
        fds = win_top;        win_top = win_bot;        win_bot = fds;
    end

    r = win_top:win_bot;    c = win_left:win_right;
    img = get(handles.hImage,'CData');
    img_s = img(r,c,:);     % Extract a sub-image
    x = x - c(1);           y = y - r(1);      % We need the PTs coords relative to the subimage
    cvlib_mex('circle',img_s,[x y],lineThick,handles.fillColor,-1,lineType)
    img(r,c,:) = img_s;     % Reimplant the painted sub-image
    set(handles.hImage,'CData',img);

% -------------------
function wbd_paint(obj,eventdata,hCallFig,state)
	uirestore_j(state, 'nochildren');		% Restore the figure's initial state

% --------------------------------------------------------------------
function pipeta_clickedcallback(hObject, eventdata)
% Pick one color from image and make it the default painting one
	handles = guidata(hObject);     % get handles
	figure(handles.hCallingFig)
	%pal = get(handles.hCallingFig,'Colormap');
	set(handles.hCallingFig,'Pointer', 'custom','PointerShapeCData',getPointer('pipeta'),'PointerShapeHotSpot',[15 1])
	w = waitforbuttonpress;
	if (w == 0)       % A mouse click
		pt = get(handles.hCallingAxes, 'CurrentPoint');
		[c,r] = getpixcoords(handles,pt(1,1),pt(1,2));
		c = round(c);       r = round(r);
		if (~insideRect(handles,c,r))
			set(handles.hCallingFig,'Pointer', 'arrow');            return;
		end
		img = get(handles.hImage,'CData');
		fillColor = double(img(r,c,:));
		handles.fillColor = reshape(fillColor,1,numel(fillColor));
		set(handles.toggle_currColor,'BackgroundColor',handles.fillColor / 255)
		set(handles.hCallingFig,'Pointer', 'arrow');
		guidata(handles.figure1,handles)
	end
    
% -------------------------------------------------------------------------------------
function flood_clickedcallback(hObject, eventdata)
% ...
	handles = guidata(hObject);     % get handles
	figure(handles.hCallingFig)
	set(handles.hCallingFig,'Pointer', 'custom','PointerShapeCData',getPointer('bucket'),'PointerShapeHotSpot',[16 15])
	[params,but] = prepareParams(handles,getPointer('bucket'),[16 15]);
	if (isempty(params) || but ~= 1 || ~insideRect(handles,params.Point(1),params.Point(2)))
		set(handles.hCallingFig,'Pointer', 'arrow');
		return
	end
	while (but == 1)
		img = get(handles.hImage,'CData');
		if (insideRect(handles,params.Point(1),params.Point(2)))
			img = cvlib_mex('floodfill',img,params);
			set(handles.hImage,'CData', img);
			handles.img_was_changed = true;
			guidata(handles.figure1, handles)
		end
		[x,y,but] = click_e_point(1,getPointer('bucket'),[16 15]);  % Get next point
		[x,y] = getpixcoords(handles,x,y);
		params.Point = [x y];
	end
	set(handles.hCallingFig,'Pointer', 'arrow')

% --------------------------------------------------------------------
function shape_clickedcallback(hObject, eventdata, opt)
	handles = guidata(hObject);     % get handles
	state = uisuspend_j(handles.hCallingFig);        % Remember initial figure state
	set(handles.hCallingFig,'Pointer', 'crosshair')
	figure(handles.hCallingFig)

	w = waitforbuttonpress;
	if (w == 0),	ShapeFirstButtonDown(handles,state,opt)       % A mouse click
	else			set(handles.hCallingFig,'Pointer', 'arrow');
	end

% -------------------
function ShapeFirstButtonDown(handles,state,opt)
	handles = guidata(handles.figure1);     % Update handles
	lineType = 8;       % Default to 8 connectivity
	if (get(handles.checkbox_AA,'Value')),      lineType = 16;      end
	pt = get(handles.hCallingAxes, 'CurrentPoint');
	[x,y] = getpixcoords(handles,pt(1,1),pt(1,2));
	x = round(x);    y = round(y);
	if (~insideRect(handles,x,y)),  set(handles.hCallingFig,'Pointer', 'arrow');    return;     end
	img = get(handles.hImage,'CData');
	lineThick = handles.lineWidth;
	if (get(handles.checkbox_filled,'Value'))
		lineThick = -lineThick;
	end
	if (strcmp(opt,'circ'))
		set(handles.hCallingFig,'WindowButtonMotionFcn',{@wbm_circle,handles,[x y],img,lineThick,lineType})
	elseif (strcmp(opt,'rect'))
		set(handles.hCallingFig,'WindowButtonMotionFcn',{@wbm_rectangle,handles,[x y],img,lineThick,lineType})
	elseif (strcmp(opt,'ellipse'))
		set(handles.hCallingFig,'WindowButtonMotionFcn',{@wbm_ellipse,handles,[x y],img,lineThick,lineType})
	end
	set(handles.hCallingFig,'WindowButtonDownFcn',{@wbd_paint,handles.hCallingFig,state})

	handles.img_was_changed = true;
	guidata(handles.figure1, handles)

% -------------------
function wbm_circle(obj,eventdata,handles,first_pt,IMG,lineThick,lineType)
    % Draw a circle
    pt = get(handles.hCallingAxes, 'CurrentPoint');
    x = round( getPixel_coords(handles.imgSize(2),handles.head(1:2),pt(1,1)) );
    y = round( getPixel_coords(handles.imgSize(1),handles.head(3:4),pt(1,2)) );
    if (~insideRect(handles,x,y)),      return;     end
    dx = abs(x - first_pt(1));          dy = abs(y - first_pt(2));
    rad = round(sqrt(dx*dx + dy*dy));   dt = round(abs(lineThick)/2) + 2;

    win_left = max(first_pt(1) - rad - dt,1);     win_right = min(first_pt(1) + rad + dt,handles.imgSize(2));
    win_top = max(first_pt(2) - rad - dt,1);      win_bot = min(first_pt(2) + rad + dt,handles.imgSize(1));
    
    r = win_top:win_bot;    c = win_left:win_right;
    img_s = IMG(r,c,:);     % Extract a sub-image
    x = first_pt(1) - c(1);           y = first_pt(2) - r(1);      % We need the PTs coords relative to the subimage
    cvlib_mex('circle',img_s,[x y],rad,handles.fillColor,lineThick,lineType)
    IMG(r,c,:) = img_s;     % Reimplant the painted sub-image
    set(handles.hImage,'CData',IMG);
    
% -------------------
function wbm_ellipse(obj,eventdata,handles,first_pt,IMG,lineThick,lineType)
    % Draw an ellipse
    pt = get(handles.hCallingAxes, 'CurrentPoint');
    x = round( getPixel_coords(handles.imgSize(2),handles.head(1:2),pt(1,1)) );
    y = round( getPixel_coords(handles.imgSize(1),handles.head(3:4),pt(1,2)) );
    if (~insideRect(handles,x,y)),      return;     end
    dx = abs(x - first_pt(1));          dy = abs(y - first_pt(2));
    dt = round(abs(lineThick)/2) + 2;

    win_left = max(first_pt(1) - dx - dt,1);     win_right = min(first_pt(1) + dx + dt,handles.imgSize(2));
    win_top = max(first_pt(2) - dy - dt,1);      win_bot = min(first_pt(2) + dy + dt,handles.imgSize(1));
    
    r = win_top:win_bot;    c = win_left:win_right;
    img_s = IMG(r,c,:);     % Extract a sub-image
    x = first_pt(1) - c(1);           y = first_pt(2) - r(1);      % We need the PTs coords relative to the subimage
    box.center = [x y];
    box.size = [2*dx 2*dy];
    cvlib_mex('eBox',img_s,box,handles.fillColor,lineThick,lineType)
    IMG(r,c,:) = img_s;     % Reimplant the painted sub-image
    set(handles.hImage,'CData',IMG);
    
% -------------------
function wbm_rectangle(obj,eventdata,handles,first_pt,IMG,lineThick,lineType)
    % Draw a circle
    pt = get(handles.hCallingAxes, 'CurrentPoint');
    x = round( getPixel_coords(handles.imgSize(2),handles.head(1:2),pt(1,1)) );
    y = round( getPixel_coords(handles.imgSize(1),handles.head(3:4),pt(1,2)) );
    if (~insideRect(handles,x,y)),      return;     end
    dt = round(abs(lineThick)/2) + 2;
    win_dx = abs(x - first_pt(1)) + dt + 2;         win_dy = abs(y - first_pt(2)) + dt + 2;
    
    win_left = max(x-win_dx,1);      win_right = min(x+win_dx,handles.imgSize(2));
    win_top  = max(y-win_dy,1);      win_bot   = min(y+win_dy,handles.imgSize(1));
    if (win_top > win_bot)
        fds = win_top;        win_top = win_bot;        win_bot = fds;
    end
    
    r = win_top:win_bot;    c = win_left:win_right;
    img_s = IMG(r,c,:);     % Extract a sub-image
    x = x - c(1);                       y = y - r(1);      % We need the PTs coords relative to the subimage
    first_pt(1) = first_pt(1) - c(1);   first_pt(2) = first_pt(2) - r(1);
    cvlib_mex('rectangle',img_s,[first_pt(1) first_pt(2)],[x y],handles.fillColor,lineThick,lineType)
    IMG(r,c,:) = img_s;     % Reimplant the painted sub-image
    set(handles.hImage,'CData',IMG);

% -------------------------------------------------------------------------------------
function slider_tolerance_CB(hObject, handles)
    handles.tol = round(get(hObject,'Value'));
    set(handles.text_tol,'String',['Tolerance = ' num2str(handles.tol)])
    guidata(handles.figure1,handles)
	handles.img_was_changed
        
% -------------------------------------------------------------------------------------
function [params,but] = prepareParams(handles, opt, opt2)
% Prepare the params structure that is to be transmited to cvlib_mex (gets also the point)
% OPT, if provided and contains a 16x16 array, is taken as a pointer
% OPT2, if provided, must be a 2 element vector with 'PointerShapeHotSpot'
% BUT returns which button has been pressed.
	if (nargin == 1),   opt = [];   opt2 = [];   end
	if (nargin == 2),   opt2 = [];  end
	if (~handles.randColor)        % That is, Cte color
		if (~isempty(handles.fillColor))
			params.FillColor = handles.fillColor;
		else
			errordlg('I don''t have yet a filling color. You must select one first.','Error')
			params = [];        return
		end
	end
	figure(handles.hCallingFig)         % Bring the figure containing image forward
	poin = 'crosshair';
	if (numel(opt) == 256),        poin = opt;    end
	[x,y,but]  = click_e_point(1,poin,opt2);
	[x,y] = getpixcoords(handles,x,y);
	params.Point = [x y];
	params.Tolerance = handles.tol;
	params.Connect = handles.connect;

% -------------------------------------------------------------------------------------
function [x,y] = getpixcoords(handles,x,y)
% Convert x,y to pixel coordinates (they are not when the image has other coordinates)
	if (handles.head(7))                % Image is pixel registered
		X = [handles.head(1) handles.head(2)] + [handles.head(8) -handles.head(8)]/2;
		Y = [handles.head(3) handles.head(4)] + [handles.head(9) -handles.head(9)]/2;
	else                                % Image is grid registered
		X = [handles.head(1) handles.head(2)];
		Y = [handles.head(3) handles.head(4)];
	end
	x = getPixel_coords(handles.imgSize(2),X,x);
	y = getPixel_coords(handles.imgSize(1),Y,y);

% --------------------------------------------------------------------
function res = insideRect(handles,x,y)
% Check if the point x,y in PIXELS is inside the rectangle RECT
% RECT = [1 handles.imgSize(2) 1 handles.imgSize(1)]
	res = ( x >= 1 && x <= handles.imgSize(2) && y >= 1 && y <= handles.imgSize(1) );

% -------------------------------------------------------------------------------------
function radio_randColor_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),	return,	end
	set(handles.radio_cteColor,'Value',0)
	handles.randColor = 1;
	set(findobj(handles.figure1,'Style','toggle'),'Enable','Inactive')
	set(handles.push_moreColors,'Enable','Inactive')
	set(handles.toggle_currColor,'BackgroundColor','w')
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function radio_cteColor_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),	return,	end
	set(handles.radio_randColor,'Value',0)
	handles.randColor = 0;
	set(findobj(handles.figure1,'Style','toggle'),'Enable','on');
	set(handles.push_moreColors,'Enable','on')
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function radio_fourConn_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),	return,	end
	set(handles.radio_eightConn,'Value',0)
	handles.connect = 4;
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function radio_eightConn_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),	return,	end
	set(handles.radio_fourConn,'Value',0)
	handles.connect = 8;
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function toggle00_CB(hObject, eventdata, handles)
% All color toggle end up here
	toggleColors(hObject,handles)

% -------------------------------------------------------------------------------------
function toggleColors(hCurr,handles)
% hCurr is the handle of the current slected toggle color. Get its color
% and assign it to toggle_currColor.
	set(hCurr,'Value',1);               % Reset it to pressed state
	set(handles.toggle_currColor,'BackgroundColor',get(hCurr,'BackgroundColor'))
	handles.fillColor = round(get(hCurr,'BackgroundColor')*255);
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function push_moreColors_CB(hObject, handles)
	c = uisetcolor;
	if (length(c) > 1)          % That is, if a color was selected
		handles.fillColor = round(c*255);
		set(handles.toggle_currColor,'BackgroundColor',c)
		guidata(handles.figure1,handles)
	end

% -------------------------------------------------------------------------------------
function push_pickSingle_CB(hObject, handles)
	[params,but] = prepareParams(handles);
	if (isempty(params) || but ~= 1),   return;     end
	img = get(handles.hImage,'CData');              % Get the image
	if (~get(handles.check_mahal,'Val'))
		[dumb,mask] = cvlib_mex('floodfill',img,params);
		clear dumb		% Humm, shouldn't we take care of this inside cvlib_mex?
	else
		r = round(params.Point(2));	c = round(params.Point(1));
		r_min = max(r-3,0);			r_max = min(r+3, size(img,1));
		c_min = max(c-3,0);			c_max = min(c+3, size(img,2));
		I = img(r_min:r_max, c_min:c_max, :);
		I = reshape(I, size(I,1)*size(I,2), 3);
		[C, m] = covmatrix(I);
		T = handles.tol;
		mask = colorseg(img, T, m, C);
	end
	if (get(handles.checkbox_useDilation,'Value'))
		mask  = img_fun('bwmorph',mask,'dilate');
	end
	if (handles.colorSegment)
		if (ndims(img) == 3)
			mask = repmat(mask,[1 1 3]);
		end
		img(~mask) = handles.bg_color;
		if (handles.image_type == 2 || handles.image_type == 20)
			h = mirone(img);
			set(h,'ColorMap',get(handles.hCallingFig,'ColorMap'),'Name','Color segmentation')
		else
			tmp.X = handles.head(1:2);  tmp.Y = handles.head(3:4);  tmp.head = handles.head;
			tmp.name = 'Color segmentation';
			mirone(img,tmp);
		end
	else
		digitize(handles, mask)
	end

% -------------------------------------------------------------------------------------
function digitize(handles,img)
% IMG is binary mask. Digitize its outer contours

	B = img_fun('bwboundaries',img,'noholes');

	for (k = 1:numel(B))
		boundary = B{k};
		if (numel(boundary) > 4)
			boundary = cvlib_mex('dp', boundary, 0.5);		% Simplify line
			B{k} = boundary;
		end
	end

	x_inc = handles.head(8);    y_inc = handles.head(9);
	x_min = handles.head(1);    y_min = handles.head(3);
	if (handles.head(7))							% Work in grid registration
		x_min = x_min + x_inc/2;
		y_min = y_min + y_inc/2;
	end

	% Renew the catch of line thickness/color because they may have been reset meanwhile.
	if (~isempty(handles.IAmAMirone))
		handMir = guidata(handles.hCallingFig);
		lineThickness = handMir.DefLineThick;
		lineColor = handMir.DefLineColor;
	else
		lineThickness = handles.DefLineThick;
		lineColor = handles.DefLineColor;
	end

	if (handles.single_poly)						% Draw a single polygon
		% Add NaNs to the end of each polygon
		nElem = zeros(length(B)+1,1);
		for k = 1:length(B)
			B{k}(end+1,1:2) = [NaN NaN];
			nElem(k+1) = size(B{k},1);
		end
		soma = cumsum(nElem);

		x = zeros(soma(end),1);     y = x;
		for k = 1:length(B)
			y(soma(k)+1:soma(k+1)) = (B{k}(:,1)-1) * y_inc + y_min;
			x(soma(k)+1:soma(k+1)) = (B{k}(:,2)-1) * x_inc + x_min;
		end

		h_edge = line(x, y, 'Linewidth',lineThickness, 'Color',lineColor, ...
				'Tag','shape_detected','Userdata',handles.udCount);

		multi_segs_str = cell(length(h_edge),1);	% Just create a set of empty info strings
		draw_funs(h_edge,'isochron',multi_segs_str);
		handles.udCount = handles.udCount + 1;
	else											% Draw separate polygons
		for k = 1:numel(B)
			x = (B{k}(:,2)-1) * x_inc + x_min;
			y = (B{k}(:,1)-1) * y_inc + y_min;
			h_edge = line(x, y, 'Linewidth',lineThickness, 'Color',lineColor, 'Tag','shape_detected');
			draw_funs(h_edge,'line_uicontext')		% Set lines's uicontextmenu
		end
	end
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function push_pickMultiple_CB(hObject, handles)
% I left the code that deals with processing in Lab & HSV color models but removed
% those options from the GUI. This is for the case that I change my mind and decide
% to reintroduce it. For the time beeing, RGB seams to work better.
	[params, but] = prepareParams(handles);
	if (isempty(params) || but ~= 1),   return;     end
	img = get(handles.hImage,'CData');              % Get the image
	nColors = 0;
	mask = false([size(img,1) size(img,2)]);        % Initialize the mask
	while (but == 1)
		nColors = nColors + 1;
		[dumb,mask(:,:,nColors)] = cvlib_mex('floodfill',img,params);
		[x,y,but] = click_e_point(1,'crosshair');  % Get next point
		[x,y] = getpixcoords(handles,x,y);
		params.Point = [x y];
	end

	if (~isempty(handles.colorModel))
		if (ndims(img) == 2)
			img = ind2rgb8(img,get(handles.hCallingFig,'ColorMap'));
		end
		img = cvlib_mex('color',img,['rgb2' handles.colorModel]);
	end

	if (isempty(handles.colorModel) && ndims(img) == 3)        % That is, RGB model
		a = img(:,:,1);		b = img(:,:,2);		c = img(:,:,3);
		cm = zeros(nColors, 6);
		im_out = repmat(uint8(handles.bg_color), size(img));
		for count = 1:nColors
			cm(count,1) = mean2(a(mask(:,:,count)));
			cm(count,1:2) = [max(cm(count,1) - handles.tol, 0) min(cm(count,1) + handles.tol, 255)];
			cm(count,3)   = mean2(b(mask(:,:,count)));
			cm(count,3:4) = [max(cm(count,3) - handles.tol, 0) min(cm(count,3) + handles.tol, 255)];
			cm(count,5)   = mean2(c(mask(:,:,count)));
			cm(count,5:6) = [max(cm(count,5) - handles.tol, 0) min(cm(count,5) + handles.tol, 255)];
			tmp = (a >= cm(count,1) & a <= cm(count,2) & b >= cm(count,3) & b <= cm(count,4) & ...
				c >= cm(count,5) & c <= cm(count,6));
			tmp  = img_fun('bwmorph',tmp,'clean');      % Get rid of isolated pixels
			if (get(handles.checkbox_useDilation,'Value'))
				tmp  = img_fun('bwmorph',tmp,'dilate');
				%tmp  = cvlib_mex('dilate',tmp);
				%tmp  = cvlib_mex('morphologyex',tmp,'close');
			end
			if (handles.colorSegment)			% Do color segmentation. One figure for each color
				t = im_out(:,:,1);		t(tmp) = a(tmp);	im_out(:,:,1) = t;
				t = im_out(:,:,2);		t(tmp) = b(tmp);	im_out(:,:,2) = t;
				t = im_out(:,:,3);		t(tmp) = c(tmp);	im_out(:,:,3) = t;
			else								% Create contours from mask
				digitize(handles, tmp)
			end
		end
		if (handles.colorSegment && handles.IAmAMirone)
			if (handles.image_type == 2 || handles.image_type == 20)
				h = mirone(im_out);
				set(h,'Name','Color segmentation')
			else
				hdr.X = handles.head(1:2);  hdr.Y = handles.head(3:4);  hdr.head = handles.head;
				hdr.name = sprintf('Color segmentation n? %d', count);
				mirone(im_out,hdr);
			end
		elseif (handles.colorSegment)
			figure;		image(im_out);
		end

	elseif (isempty(handles.colorModel) && ndims(img) == 2)        % That is, indexed image & RGB model
		cm = zeros(nColors, 2);
		im_out = repmat(uint8(handles.bg_color), size(img));
		for count = 1:nColors
			cm(count,1) = mean2(img(mask(:,:,count)));
			cm(count,1:2) = [max(cm(count,1) - handles.tol, 0) min(cm(count,1) + handles.tol, 255)];
			tmp = (img >= cm(count,1) & img <= cm(count,2));
			tmp  = img_fun('bwmorph',tmp,'clean');      % Get rid of isolated pixels
			if (get(handles.checkbox_useDilation,'Value'))
				tmp  = img_fun('bwmorph',tmp,'dilate');
			end
			if (handles.colorSegment)			% Do color segmentation. One figure for each color
				im_out(tmp) = img(tmp);
			else								% Create contours from mask
				digitize(handles, tmp)
			end
		end
		if (handles.colorSegment && handles.IAmAMirone)
			if (handles.image_type == 2 || handles.image_type == 20)
				h = mirone(im_out);
				set(h,'ColorMap',get(handles.hCallingFig,'ColorMap'),'Name','Color segmentation')
			else
				hdr.X = handles.head(1:2);  hdr.Y = handles.head(3:4);  hdr.head = handles.head;
				hdr.name = sprintf('Color segmentation n? %d', count);
				hdr.cmap = get(handles.hCallingFig,'ColorMap');
				mirone(im_out, hdr);
			end
		elseif (handles.colorSegment)
			h = figure;     image(im_out);
			set(h,'ColorMap',get(handles.hCallingFig,'ColorMap'),'Name','Color segmentation')
		end

	else										% Either YCrCb or Lab model were used
		a = img(:,:,2);
		b = img(:,:,3);
		cm = zeros(nColors, 4);
		for count = 1:nColors
			cm(count,1) = min(min(a(mask(:,:,count))));
			cm(count,2) = max(max(a(mask(:,:,count))));
			cm(count,3) = min(min(b(mask(:,:,count))));
			cm(count,4) = max(max(b(mask(:,:,count))));
			tmp = (a >= cm(count,1) & a <= cm(count,2) & b >= cm(count,3) & b <= cm(count,4));
			%tmp = bwareaopen(tmp,20,4);
			if (handles.colorSegment)           % Do color segmentation. One figure for each color
				tmp = repmat(tmp,[1 1 3]);
				img(~tmp) = handles.bg_color;
				if (handles.IAmAMirone)
					if (handles.image_type == 2 || handles.image_type == 20)
						mirone(img);
					else
						tmp.X = handles.head(1:2);  tmp.Y = handles.head(3:4);  tmp.head = handles.head;
						tmp.name = sprintf('Color segmentation n? %d', count);
						mirone(img,tmp);
					end
				else
					figure;		image(img);
				end
			else                                % Create contours from mask
				digitize(handles,tmp)
			end
		end    
	end

% -------------------------------------------------------------------------------------
% function radio_YCrCb_CB(hObject, handles)
%     if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
%     set(handles.radio_Lab,'Value',0)
%     set(handles.radio_RGB,'Value',0)
%     handles.colorModel = 'YCrCb';
%     guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
% function radio_Lab_CB(hObject, handles)
%     if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
%     set(handles.radio_YCrCb,'Value',0)
%     set(handles.radio_RGB,'Value',0)
%     handles.colorModel = 'lab';
%     guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
% function radio_RGB_CB(hObject, handles)
%     % This the only one know (it seams to do better than the others)
%     if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
%     set(handles.radio_YCrCb,'Value',0)
%     set(handles.radio_Lab,'Value',0)
%     handles.colorModel = '';
%     guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function radio_colorSegment_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_digitize,'Value',0)
	handles.colorSegment = 1;
	set(handles.edit_minPts,'Visible','off')
	set(handles.text_minPts,'Visible','off')
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function radio_digitize_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_colorSegment,'Value',0)
	handles.colorSegment = 0;
	set(handles.edit_minPts,'Visible','on')
	set(handles.text_minPts,'Visible','on')
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function edit_minPts_CB(hObject, handles)
	xx = round( str2double(get(hObject,'String')) );
	if (isnan(xx)),     set(hObject,'String','50');     return;     end
	handles.minPts = xx;
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function push_restoreImg_CB(hObject, handles)
    set(handles.hImage,'CData',handles.origFig)

%-----------------------------------------------------------------------------------------
function pix_coords = getPixel_coords(img_length, XData, axes_coord)
% Convert coordinates from axes (real coords) to image (pixel) coordinates.
% IMG_LENGTH is the image width (n_columns)
% XDATA is the image's [x_min x_max] in axes coordinates
% AXES_COORD is the (x,y) coordinate of the point(s) to be converted

	slope = (img_length - 1) / (XData(end) - XData(1));
	if ((XData(1) == 1) && (slope == 1))
		pix_coords = axes_coord;
	else
		pix_coords = slope * (axes_coord - XData(1)) + 1;
	end

% -------------------------------------------------------------------------------------
function y = mean2(x)
% Compute mean of matrix elements.
	y = sum(x(:)) / numel(x);

% -------------------------------------------------------------------------------------
function listbox_lineWidth_CB(hObject, handles)
    handles.lineWidth = get(hObject,'Value');
    guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------
function radio_square_CB(hObject, handles)
    if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
    handles.elemSquare = 1;
    set(handles.radio_round,'Value',0)
    guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------
function radio_round_CB(hObject, handles)
    if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
    handles.elemSquare = 0;
    set(handles.radio_square,'Value',0)
    guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------
function check_mahal_CB(hObject, handles)
	if (get(hObject,'Value')),		set(handles.push_pickMultiple,'Enable', 'off')
	else							set(handles.push_pickMultiple,'Enable', 'on')
	end

% ---------------------------------------------------------------------
function I = colorseg(f, T, m, C)
%COLORSEG Performs segmentation of a color image.
%	Minimalist code from the DIPUM version
%
%   S = COLORSEG(F, T, M, C) performs segmentation of
%   color image F using the Mahalanobis distance as a measure of
%   similarity. C is the 3-by-3 covariance matrix of the sample color
%   vectors of the class of interest. See function covmatrix for the
%   computation of C and M. 
%
%   I is the segmented image (a binary matrix) in which 0s denote the
%   background. 

	% Preliminaries. Recall that varargin is a cell array.
	if ( (ndims(f) ~= 3) || (size(f, 3) ~= 3) ),	error('Input image must be RGB.');	end
	M = size(f, 1);		N = size(f, 2);
	% Convert f to vector format using function imstack2vectors.
	f = reshape(f, M*N, size(f, 3));
	
	f = double(f);
	% Initialize I as a column vector.  It will be reshaped later into an image.
	I = false(M*N, 1); 
	m = m(:)'; % Make sure that m is a row vector.

	D = mahalanobis(f, C, m);

	% D is a vector of size MN-by-1 containing the distance computations
	% from all the color pixels to vector m. Find the distances <= T.
	I(D <= T) = true;
	% Reshape I into an M-by-N image.
	I = reshape(I, M, N);

% -------------------------------------------------------------------------------------
function [C, m] = covmatrix(X)
%COVMATRIX Computes the covariance matrix of a vector population.
%   [C, M] = COVMATRIX(X) computes the covariance matrix C and the
%   mean vector M of a vector population organized as the rows of
%   matrix X. C is of size N-by-N and M is of size N-by-1, where N is
%   the dimension of the vectors (the number of columns of X).

	K = size(X,1);
	X = double(X);
	% Compute an unbiased estimate of m.
	m = sum(X, 1)/K;
	% Subtract the mean from each row of X.
	X = X - m(ones(K, 1), :);
	% Compute an unbiased estimate of C. Note that the product is
	% X'*X because the vectors are rows of X.	
	C = (X'*X)/(K - 1);
	m = m';		% Convert to a column vector.	 

% ---------------------------------------------------------------------
function d = mahalanobis(varargin)
%MAHALANOBIS Computes the Mahalanobis distance.
%   D = MAHALANOBIS(Y, X) computes the Mahalanobis distance between
%   each vector in Y to the mean (centroid) of the vectors in X, and
%   outputs the result in vector D, whose length is size(Y, 1).  The
%   vectors in X and Y are assumed to be organized as rows.  The
%   input data can be real of complex. The outputs are real quantities.

	param = varargin; % Keep in mind that param is a cell array.
	Y = param{1};
	ny = size(Y, 1); % Number of vectors in Y.
	
	if numel(param) == 2
		X = param{2};
		% Compute the mean vector and covariance matrix of the vectors in X.
		[Cx, mx] = covmatrix(X);
	elseif (numel(param) == 3)						% Cov. matrix and mean vector provided.
		Cx = param{2};
		mx = param{3};
	else 
		error('Wrong number of inputs.')
	end
	if (size(mx,1) > 1),	mx = mx(:)';	end		% Make sure that mx is a row vector.
	
	% Subtract the mean vector from each vector in Y.
	Yc = Y - mx(ones(ny, 1), :);	
	
	% Compute the Mahalanobis distances.
	d = real(sum(Yc / Cx .* conj(Yc), 2));

% -------------------------------------------------------------------------------------------------------
function p = getPointer(opt)

if (strcmp(opt,'pipeta'))
	p = [...
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	1	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	1
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	1	1	1
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	NaN	1	1	1	1	1
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	1	1	1	1	1
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	1	1	1	1	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	1	1	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	1	1	1	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	1	1	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	1	NaN	NaN	1	1	1	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	1	NaN	1	1	1	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	1	NaN	1	1	1	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	1	NaN	1	1	1	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	1	1	1	1	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		1	1	1	1	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN];
        
elseif (strcmp(opt,'pencil'))
	p = [...
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	1	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	1	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	1	NaN	NaN	NaN	1	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	1	NaN	1	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	1	1	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	1	1	1	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	1	1	1	1	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	1	1	1	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN];
        
elseif (strcmp(opt,'brush'))
	p = [...
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	1	1
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	1	NaN	1	1	1	1	1
		NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	1	1	1	1	1	NaN
		NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	1	1	NaN	NaN
		NaN	NaN	NaN	1	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN
		1	1	1	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	1	NaN	NaN
		NaN	1	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN
		NaN	NaN	1	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN
		NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	1	1	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN];
elseif (strcmp(opt,'bucket'))
	p = [...
		NaN	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	1	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	1	1	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	NaN	1	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	NaN	1	NaN	NaN	1	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	NaN	1	NaN	NaN	NaN	1	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN
		NaN	NaN	1	NaN	NaN	NaN	1	NaN	1	NaN	1	NaN	NaN	NaN	NaN	NaN
		NaN	1	NaN	NaN	NaN	NaN	NaN	1	1	NaN	NaN	1	NaN	NaN	NaN	NaN
		NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	NaN	NaN
		NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	1	1	1	NaN
		NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	1	1	NaN
		NaN	NaN	NaN	NaN	1	NaN	NaN	NaN	1	NaN	NaN	NaN	NaN	1	1	NaN
		NaN	NaN	NaN	NaN	NaN	1	1	1	NaN	NaN	NaN	NaN	NaN	1	1	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	NaN
		NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1	1	NaN];    
end

% -----------------------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
% Executes when user closes figure
	handles = guidata(hObject);
	if (handles.must_restore_img && ~handles.img_was_changed)	% Restore the original indexed image
		set(handles.hImage, 'CData', handles.origFig)			% but only when image was not changed
	end
	delete(handles.figure1)

% --- Creates and returns a handle to the GUI figure. 
function floodfill_LayoutFcn(h1)

set(h1,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Flood Fill',...
'NumberTitle','off',...
'Position',[520 390 185 366],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[4 30 177 110],'Style','frame');
uicontrol('Parent',h1,'Position',[4 209 177 153],'Style','frame');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@floodfill_uiCB,...
'Max',255,...
'Position',[4 171 176 16],...
'Style','slider',...
'SliderStep',[0.00390625 0.1],...
'TooltipString','Color detection equal pixel color +/- tolerance',...
'Tag','slider_tolerance');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[12 105 160 21],...
'String','Pick single shape',...
'TooltipString','Pick up the body''s shape with the selected color',...
'Tag','push_pickSingle');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[12 277 70 15],...
'String','4 conn',...
'Style','radiobutton',...
'TooltipString','Floodfill with a connectivity of 4',...
'Value',1,...
'Tag','radio_fourConn');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[80 277 90 15],...
'String','8 connectivity',...
'Style','radiobutton',...
'TooltipString','Floodfill with a connectivity of 8',...
'Tag','radio_eightConn');

uicontrol('Parent',h1,'BackgroundColor',[1 0.501960784313725 0.501960784313725],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[13 315 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,'BackgroundColor',[1 0 0],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[29 315 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,'BackgroundColor',[0.501960784313725 0.250980392156863 0.250980392156863],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[45 315 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,'BackgroundColor',[0.250980392156863 0 0],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[61 315 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,'BackgroundColor',[1 1 0.501960784313725],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[77 315 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'BackgroundColor',[1 0.501960784313725 0],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[93 315 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'BackgroundColor',[0.501960784313725 0.501960784313725 0],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[109 315 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'BackgroundColor',[0 1 0],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[125 315 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'BackgroundColor',[0 0.501960784313725 0],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[13 299 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'BackgroundColor',[0 0.501960784313725 0.501960784313725],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[29 299 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[45 299 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'BackgroundColor',[0 0 1],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[61 299 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'BackgroundColor',[0 1 1],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[77 299 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'BackgroundColor',[1 0.501960784313725 0.752941176470588],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[93 299 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'BackgroundColor',[0.250980392156863 0 0.501960784313725],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[109 299 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@floodfill_uicallback,h1,'toggle00_CB'},...
'Enable','inactive',...
'Position',[125 299 15 15],...
'Style','togglebutton',...
'Value',1);

uicontrol('Parent',h1,...
'FontSize',9,...
'Position',[69 353 40 15],...
'String','Paint',...
'Style','text','Tag','text_Paint');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[12 339 110 15],...
'String','Random colors',...
'Style','radiobutton',...
'TooltipString','Repainting will use a randomly selected color',...
'Value',1,...
'Tag','radio_randColor');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[115 339 65 15],...
'String','Cte color',...
'Style','radiobutton',...
'TooltipString','Repainting will use a selected color',...
'Tag','radio_cteColor');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Enable','inactive',...
'Position',[143 314 30 19],...
'String','More',...
'TooltipString','Chose other color',...
'Tag','push_moreColors');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'Position',[143 299 30 16],...
'Style','togglebutton',...
'TooltipString','Current selected color',...
'Value',1,...
'HandleVisibility','off',...
'Tag','toggle_currColor');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[12 77 160 21],...
'String','Pick multiple shapes',...
'TooltipString','Find out all bounding polygons that share the selected color',...
'Tag','push_pickMultiple');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[12 55 140 15],...
'String','Color segmentation',...
'Style','radiobutton',...
'TooltipString','Create a separate figure with shapes colored as you selected',...
'Value',1,...
'Tag','radio_colorSegment');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[12 35 65 15],...
'String','Digitize',...
'Style','radiobutton',...
'TooltipString','Detect bounding polygon to the colored selected body',...
'Tag','radio_digitize');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@floodfill_uiCB,...
'Position',[140 35 33 20],...
'String','50',...
'Style','edit',...
'TooltipString','Shapes with less than this number of vertex won''t be ploted',...
'Tag','edit_minPts',...
'Visible','off');

uicontrol('Parent',h1,...
'Position',[118 38 20 15],...
'String','Min',...
'Style','text',...
'Tag','text_minPts',...
'Visible','off');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[48 187 90 14],'String','Tolerance = 20',...
'Style','text','Tag','text_tol');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[37 4 111 21],'String','Restore image',...
'Tag','push_restoreImg');

uicontrol('Parent',h1,...
'Position',[4 152 95 15],'String','Use Dilation',...
'Style','checkbox',...
'TooltipString','Use dilation operation to find better limits between neighboring shapes',...
'Value',1,...
'Tag','checkbox_useDilation');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[95 152 92 15],'String','Mahalanobis',...
'Style','checkbox',...
'TooltipString','Do color segmentation using Mahalanobis distance',...
'Value',0,...
'Tag','check_mahal');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@floodfill_uiCB,...
'Position',[10 238 40 32],'String',{'LineThickness'},...
'Style','listbox',...
'TooltipString','Line thickness in pixels',...
'Value',1,...
'Tag','listbox_lineWidth');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[80 255 95 15],'String','Square element',...
'Style','radiobutton',...
'TooltipString','Use square elements in line drawings',...
'Value',1,...
'Tag','radio_square');

uicontrol('Parent',h1,...
'Call',@floodfill_uiCB,...
'Position',[80 235 91 15],...
'String','Round element',...
'Style','radiobutton',...
'TooltipString','Use round elements in line drawings',...
'Tag','radio_round');

uicontrol('Parent',h1,'Position',[10 215 80 15],...
'String','Antialiasing','Style','checkbox',...
'TooltipString','Use anti aliasing in line drwaings',...
'Tag','checkbox_AA');

uicontrol('Parent',h1,'Position',[90 215 100 15],...
'String','Filled forms','Style','checkbox',...
'TooltipString','When drawing circles, ellipses or rectangles draw them filled',...
'Tag','checkbox_filled');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[49 131 90 16],...
'String','Digit / Segment',...
'Style','text','Tag','text_DS');

function floodfill_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,[],guidata(h1));

function floodfill_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));


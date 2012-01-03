function varargout = mpaint(varargin)
% Works i similar way to the Windows Paint program (except the ROI image menipulations)
%
% USAGE:
%   MPAINT(HANDLE) operates in the single image contained in HANDLE.
%   HANDLE can be a handle to a figure, axes, or image object.
%   For the floodfill painting, every left-click paint a group of connected pixels.
%   (see below the 'connectivity' and 'tolerance'). A right click stops the painting.
%   For all other operations a left-click starts the painting, which continues for
%   every mouth movement and stops whith any mouse click.
%
%   H = MPAINT(HANDLE); returns the handles of the MPAINT figure
%   
%   Specifying connectivity:
%   By default, MPAINT uses 4-connected neighbors but you can select
%   8-connected in the apropriate radio buttons
%
%   Specifying tolerance:
%   Tolerance is used only when the doing painting operations (when using the bucket icon)
%   The "Tolerance" slider selects the maximal lower/higher brightness/color difference
%   between the currently hited pixel and its neighbors. Pixels within the -> hit pixel
%   +- |tolerance| are considered to belong to the repainted domain.
%
%   NOTE: MPAINT is an GUI interface to some functions in the CVLIB_MEX
%   library, so you must have it installed to run this function
%
% EXAMPLES:
%   img = imread('j.jpg')
%   h = image(img);                 % display the image
%   hFig = mpaint(h);               % Click at will
%   B = getimage(hFig);             % Get your work of art

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
 
hObject = figure('Tag','figure1','Visible','off');
mpaint_LayoutFcn(hObject);
handles = guihandles(hObject);

if (~isempty(varargin))
    handles.hCallingFig = varargin{1};
else
    handles.hCallingFig = gcf;        % Useless with Mirone figures because they are hiden to gcf
end

if (~ishandle(handles.hCallingFig))
    errordlg('Error, first argument is not a valid handle','ERROR');
    delete(hObject);    return
end

if ( strcmp(get(handles.hCallingFig,'type'),'axes') || strcmp(get(handles.hCallingFig,'type'),'figure'))
    handles.hImage = findobj(handles.hCallingFig,'type','image');
elseif ( strcmp(get(handles.hCallingFig,'type'),'image') )
    handles.hImage = handles.hCallingFig;
end

if (isempty(handles.hImage))
    errordlg('Error, no image found','ERROR');
    delete(hObject);    return
end

% Ok now we will fish the tue fig handle (before it could have been an img handle for example)
handles.hCallingFig = get(get(handles.hImage,'Parent'),'Parent');

handMir = guidata(handles.hCallingFig);
 
% Import icons
load ([handMir.path_data 'mirone_icons.mat'],'lapis_ico','trincha_ico','pipeta_ico','balde_ico',...
    'circ_ico','rectang_ico','ellipse_ico','R_ico');

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
uipushtool('parent',h_toolbar,'Click',{@restore_clickedcallback},'cdata',R_ico,'Tooltip','Restore Image','Sep','on');

% Add this figure handle to the carra?as list
plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
plugedWin = [plugedWin hObject];
setappdata(handles.hCallingFig,'dependentFigs',plugedWin);

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
handles.XData = get(handles.hImage,'XData');
handles.YData = get(handles.hImage,'YData');
handles.origImg = get(handles.hImage,'CData');        % Make a copy of the original image
handles.imgSize = size(handles.origImg);
handles.head = handMir.head;

handles.connect = 4;
handles.tol = 20;
handles.fillColor = [255 255 255];
handles.lineWidth = 2;
handles.elemSquare = 1;

set(handles.listbox_lineWidth,'String',1:99,'Val',2)
set(handles.slider_tolerance,'Value',handles.tol)

% Choose default command line output for mpaint_export
if (nargout),   varargout{1} = hObject;     end
guidata(hObject, handles);
set(hObject,'Visible','on');

% --------------------------------------------------------------------
function restore_clickedcallback(hObject, eventdata)
    % Restore original (well its RGB version if it was indexed) image
    handles = guidata(hObject);     % get handles
    set(handles.hImage,'CData',handles.origImg)

% --------------------------------------------------------------------
function shape_clickedcallback(hObject, eventdata, opt)
    handles = guidata(hObject);     % get handles
    state = uisuspend_j(handles.hCallingFig);        % Remember initial figure state
    set(handles.hCallingFig,'Pointer', 'crosshair')
    figure(handles.hCallingFig)
    
    w = waitforbuttonpress;
    if (w == 0),    ShapeFirstButtonDown(handles,state,opt)       % A mouse click
    else            set(handles.hCallingFig,'Pointer', 'arrow');
    end

% -------------------
function ShapeFirstButtonDown(handles,state,opt)
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

% --------------------------------------------------------------------
function line_clickedcallback(hObject, eventdata, opt)
    handles = guidata(hObject);     % get handles
    figure(handles.hCallingFig)
    
    state = uisuspend_j(handles.hCallingFig);        % Remember initial figure state
    if (strcmp(opt,'pencil'))
        set(handles.hCallingFig,'Pointer', 'custom','PointerShapeCData',getPointer('pencil'),'PointerShapeHotSpot',[14 2])
    elseif (strcmp(opt,'paintbrush'))
        set(handles.hCallingFig,'Pointer', 'custom','PointerShapeCData',getPointer('brush'),'PointerShapeHotSpot',[14 2])
    end
    w = waitforbuttonpress;
    if (w == 0)       % A mouse click
        paintFirstButtonDown(handles,state,opt)
    else
        set(handles.hCallingFig,'Pointer', 'arrow');
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

% -------------------
function wbm_line(obj,eventdata,handles,lineThick,lineType)
    % Draw the line using a square element
    pt = get(handles.hCallingAxes, 'CurrentPoint');
    prev_pt = getappdata(handles.figure1,'prev_pt');
    setappdata(handles.figure1,'prev_pt',pt(1,1:2))
    x = round( getPixel_coords(handles.imgSize(2),handles.XData,[prev_pt(1) pt(1,1)]) );
    y = round( getPixel_coords(handles.imgSize(1),handles.YData,[prev_pt(2) pt(1,2)]) );
    if (~insideRect(handles,x(2),y(2))),      return;     end
    win_dx = abs(diff(x)) + 4;         win_dy = abs(diff(y)) + 4;
    
    % Notice that the (2) index denotes the current point
    win_left = max(x(2)-win_dx,1);      win_right = min(x(2)+win_dx,handles.imgSize(2));
    win_top  = max(y(2)-win_dy,1);      win_bot   = min(y(2)+win_dy,handles.imgSize(1));
    if (win_top > win_bot)
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
    x = round( getPixel_coords(handles.imgSize(2),handles.XData,pt(1,1)) );
    y = round( getPixel_coords(handles.imgSize(1),handles.YData,pt(1,2)) );
    if (~insideRect(handles,x,y)),      return;     end
    
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
function wbm_circle(obj,eventdata,handles,first_pt,IMG,lineThick,lineType)
    % Draw a circle
    pt = get(handles.hCallingAxes, 'CurrentPoint');
    x = round( getPixel_coords(handles.imgSize(2),handles.XData,pt(1,1)) );
    y = round( getPixel_coords(handles.imgSize(1),handles.YData,pt(1,2)) );
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
    x = round( getPixel_coords(handles.imgSize(2),handles.XData,pt(1,1)) );
    y = round( getPixel_coords(handles.imgSize(1),handles.YData,pt(1,2)) );
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
    x = round( getPixel_coords(handles.imgSize(2),handles.XData,pt(1,1)) );
    y = round( getPixel_coords(handles.imgSize(1),handles.YData,pt(1,2)) );
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
    
% -------------------
function wbd_paint(obj,eventdata,hCallFig,state)
    uirestore_j(state, 'nochildren');		% Restore the figure's initial state

% --------------------------------------------------------------------
function pipeta_clickedcallback(hObject, eventdata)
    % Pick one color from image and make it the default painting one
    handles = guidata(hObject);     % get handles
    figure(handles.hCallingFig)
    pal = get(handles.hCallingFig,'Colormap');
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
        if (ndims(img) == 2)            % Here we have to permanently change the image type to RGB
            img = ind2rgb8(img,get(handles.hCallingFig,'ColorMap'));
            handles.origImg = img;      % Update the copy of the original image
            set(handles.hImage,'CData', img); 
        end
        fillColor = double(img(r,c,:));
        handles.fillColor = reshape(fillColor,1,numel(fillColor));
        set(handles.toggle_currColor,'BackgroundColor',handles.fillColor / 255)
        set(handles.hCallingFig,'Pointer', 'arrow');
        guidata(handles.figure1,handles)
    end

% -------------------------------------------------------------------------------------
function flood_clickedcallback(hObject, eventdata)
    handles = guidata(hObject);     % get handles
    figure(handles.hCallingFig)
    set(handles.hCallingFig,'Pointer', 'custom','PointerShapeCData',getPointer('bucket'),'PointerShapeHotSpot',[16 15])
    [params,but] = prepareParams(handles,getPointer('bucket'),[16 15]);
    if (isempty(params) || but ~= 1 || ~insideRect(handles,params.Point(1),params.Point(2)))
        set(handles.hCallingFig,'Pointer', 'arrow');        return;
    end
    while (but == 1)
        img = get(handles.hImage,'CData');
        if (ndims(img) == 2)            % Here we have to permanently change the image type to RGB
            img = ind2rgb8(img,get(handles.hCallingFig,'ColorMap'));
            handles.origImg = img;      % Update the copy of the original image
            guidata(handles.figure1,handles)
        end
        img = cvlib_mex('floodfill',img,params);
        set(handles.hImage,'CData', img); 
        [x,y,but] = click_e_point(1,getPointer('bucket'),[16 15]);  % Get next point
        [x,y] = getpixcoords(handles,x,y);
        params.Point = [x y];
    end
    set(handles.hCallingFig,'Pointer', 'arrow')

% -------------------------------------------------------------------------------------
function [params,but] = prepareParams(handles, opt, opt2)
% Prepare the params structure that is to be transmited to cvlib_mex (gets also the point)
% OPT, if provided and contains a 16x16 array, is taken as a pointer
% OPT2, if provided, must be a 2 element vector with 'PointerShapeHotSpot'
% BUT returns which button has been pressed.
    if (nargin == 1),   opt = [];   opt2 = [];   end
    if (nargin == 2),   opt2 = [];  end
    params.FillColor = handles.fillColor;
    figure(handles.hCallingFig)         % Bring the figure containing image forward
    but = 1;                            % Initialize to a left click
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
    [m n k] = size(get(handles.hImage,'CData'));
    x = getPixel_coords(n,X,x);
    y = getPixel_coords(m,Y,y);

% --------------------------------------------------------------------
function res = insideRect(handles,x,y)
    % Check if the point x,y in PIXELS is inside the rectangle RECT
    % RECT = [1 handles.imgSize(2) 1 handles.imgSize(1)]
    res = ( x >= 1 && x <= handles.imgSize(2) && y >= 1 && y <= handles.imgSize(1) );

% -------------------------------------------------------------------------------------
function toggle00_Callback(hObject, eventdata, handles)
    % All color toggle end up here
    toggleColors(hObject,handles)

% -------------------------------------------------------------------------------------
function toggleColors(hCurr,handles)
% hCurr is the handle of the current slected toggle color. Get its color
% and assign it to toggle_currColor.
    set(hCurr,'Value',1);               % Reset it to pressed state
    set(handles.toggle_currColor,'BackgroundColor',get(hCurr,'BackgroundColor'))
    handles.fillColor = round(get(hCurr,'BackgroundColor')*255);
    img = get(handles.hImage,'CData');
    if (ndims(img) == 2)            % Here we have to permanently change the image type to RGB
        img = ind2rgb8(img,get(handles.hCallingFig,'ColorMap'));
        handles.origImg = img;      % Update the copy of the original image
        set(handles.hImage,'CData', img); 
    end
    guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function pushbutton_moreColors_Callback(hObject, eventdata, handles)
    c = uisetcolor;
    if (length(c) > 1)          % That is, if a color was selected
        handles.fillColor = round(c*255);
        set(handles.toggle_currColor,'BackgroundColor',c)
        img = get(handles.hImage,'CData');
        if (ndims(img) == 2)            % Here we have to permanently change the image type to RGB
            img = ind2rgb8(img,get(handles.hCallingFig,'ColorMap'));
            handles.origImg = img;      % Update the copy of the original image
            set(handles.hImage,'CData', img); 
        end        
        guidata(handles.figure1,handles)
    end

% ---------------------------------------------------------------------------
function slider_tolerance_Callback(hObject, eventdata, handles)
    handles.tol = round(get(hObject,'Value'));
    set(handles.text_tol,'String',['Tolerance = ' num2str(handles.tol)])
    guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
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
function listbox_lineWidth_Callback(hObject, eventdata, handles)
    handles.lineWidth = get(hObject,'Value');
    guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------
function radio_square_Callback(hObject, eventdata, handles)
    if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
    handles.elemSquare = 1;
    set(handles.radio_round,'Value',0)
    guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------
function radio_round_Callback(hObject, eventdata, handles)
    if (get(hObject,'Value'))
        handles.elemSquare = 0;
        set(handles.radio_square,'Value',0)
        guidata(handles.figure1, handles);
    else
        set(hObject,'Value',1)
    end
  
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

% --- Creates and returns a handle to the GUI figure. 
function mpaint_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','mpaint',...
'NumberTitle','off',...
'Position',[520 682 229 118],...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0.501960784313725 0.501960784313725],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[60 98 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0 0],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[76 98 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton2');

uicontrol('Parent',h1,...
'BackgroundColor',[0.501960784313725 0.250980392156863 0.250980392156863],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[92 98 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton3');

uicontrol('Parent',h1,...
'BackgroundColor',[0.250980392156863 0 0],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[108 98 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton4');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 0.501960784313725],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[124 98 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton5');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0.501960784313725 0],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[140 98 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton6');

uicontrol('Parent',h1,...
'BackgroundColor',[0.501960784313725 0.501960784313725 0],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[156 98 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton7');

uicontrol('Parent',h1,...
'BackgroundColor',[0 1 0],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[172 98 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton8');

uicontrol('Parent',h1,...
'BackgroundColor',[0 0.501960784313725 0],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[60 82 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton9');

uicontrol('Parent',h1,...
'BackgroundColor',[0 0.501960784313725 0.501960784313725],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[76 82 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton10');

uicontrol('Parent',h1,...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[92 82 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton11');

uicontrol('Parent',h1,...
'BackgroundColor',[0 0 1],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[108 82 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton12');

uicontrol('Parent',h1,...
'BackgroundColor',[0 1 1],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[124 82 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton13');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0.501960784313725 0.752941176470588],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[140 82 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton14');

uicontrol('Parent',h1,...
'BackgroundColor',[0.250980392156863 0 0.501960784313725],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[156 82 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton15');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@mpaint_uicallback,h1,'toggle00_Callback'},...
'Position',[172 82 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','togglebutton16');

uicontrol('Parent',h1,...
'Callback',{@mpaint_uicallback,h1,'pushbutton_moreColors_Callback'},...
'Position',[190 96 30 19],...
'String','More',...
'TooltipString','Chose other color',...
'Tag','pushbutton_moreColors');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'Position',[190 82 30 15],...
'Style','togglebutton',...
'TooltipString','Current selected color',...
'Value',1,...
'HandleVisibility','off',...
'Tag','toggle_currColor');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@mpaint_uicallback,h1,'listbox_lineWidth_Callback'},...
'Position',[10 82 40 32],...
'String',{'Listbox'},...
'Style','listbox',...
'TooltipString','Line thickness in pixels',...
'Value',1,...
'Tag','listbox_lineWidth');

uicontrol('Parent',h1,...
'Position',[10 58 120 15],...
'String','Enable antialiasing',...
'Style','checkbox',...
'TooltipString','Use anti aliasing in line drwaings',...
'Tag','checkbox_AA');

uicontrol('Parent',h1,...
'Callback',{@mpaint_uicallback,h1,'radio_square_Callback'},...
'Position',[10 38 120 15],...
'String','Square element',...
'Style','radiobutton',...
'TooltipString','Use square elements in line drawings',...
'Value',1,...
'Tag','radio_square');

uicontrol('Parent',h1,...
'Callback',{@mpaint_uicallback,h1,'radio_round_Callback'},...
'Position',[125 38 100 15],...
'String','Round element',...
'Style','radiobutton',...
'TooltipString','Use round elements in line drawings',...
'Tag','radio_round');

uicontrol('Parent',h1,...
'Position',[125 58 95 15],...
'String','Filled forms',...
'Style','checkbox',...
'TooltipString','Draw filled circles, ellipses or rectangles',...
'Tag','checkbox_filled');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@mpaint_uicallback,h1,'slider_tolerance_Callback'},...
'Max',255,...
'Position',[10 4 210 15],...
'Style','slider',...
'SliderStep',[0.00390625 0.1],...
'TooltipString','Color detection equal pixel color +/- tolerance',...
'Tag','slider_tolerance');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[82 19 90 14],...
'String','Tolerance = 20',...
'Style','text',...
'Tag','text_tol');

function mpaint_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

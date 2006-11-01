function varargout = floodFill(varargin)
% M-File changed by desGUIDE 
 
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to floodFill_export (see VARARGIN) 
 
hObject = figure('Tag','figure1','Visible','off');
handles = guihandles(hObject);
guidata(hObject, handles);
floodFill_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

if (~isempty(varargin))
    handles.hCallingFig = varargin{1};
else
    handles.hCallingFig = gcf;        % Useless with Mirone figures because they are hiden to gcf
end

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
handles.hImage = findobj(handles.hCallingFig,'Type','image');
img = get(handles.hImage,'CData');
handles.origFig = img;        % Make a copy of the original image

handles.IAmAMirone = getappdata(handles.hCallingFig,'IAmAMirone');
if (~isempty(handles.IAmAMirone))
    handlesMir = guidata(handles.hCallingFig);
    handles.head = handlesMir.head;
    %handles.origFig = handlesMir.origFig;
    handles.DefLineThick = handlesMir.DefLineThick;
    handles.DefLineColor = handlesMir.DefLineColor;
    handles.image_type = handlesMir.image_type;
else
    handles.IAmAMirone = 0;
    % Build a GMT type header info
    handles.head = [get(handles.hCallingAxes,'xlim') get(handles.hCallingAxes,'ylim') 0 0 1];
    handles.head(8) = (handles.head(2)-handles.head(1)) / size(img,2);
    handles.head(9) = (handles.head(4)-handles.head(3)) / size(img,1);
    
    handles.DefLineThick = 1;
    handles.DefLineColor = [0 0 0];
end

% Initialize some vars
handles.connect = 4;
handles.tol = 20;
handles.randColor = 1;
handles.fillColor = [];
handles.useDilation = 1;            % Dilate mask before using it in picking shapes
handles.colorSegment = 1;           % Default to do color segmentation
%handles.colorModel = 'YCrCb';       % When color segmentation, default to this color space
handles.colorModel = [];            % Do color segmentation in RGB
handles.minPts = 50;                % When digitizing, dont create polygs with less than this pts
handles.udCount = 1;                % Unique incremental identifier to set in UserData of polygons
handles.single_poly = 0;            % IF == 1 -> all detected polygons are drawn in a single multi-polygon
handles.bg_color = 0;               % Background color (or color number in cmap) 
set(handles.slider_tolerance,'Value',handles.tol)

handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes floodFill_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = floodFill_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = floodFill_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------------------
function slider_tolerance_Callback(hObject, eventdata, handles)
    handles.tol = round(get(hObject,'Value'));
    set(handles.text_tol,'String',['Tolerance = ' num2str(handles.tol)])
    guidata(handles.figure1,handles)
    
% -------------------------------------------------------------------------------------
function pushbutton_repaintSingle_Callback(hObject, eventdata, handles)
    [params,but] = prepareParams(handles, 1);
    if (isempty(params) || but ~= 1),   return;     end
    img = get(handles.hImage,'CData');
    if (ndims(img) == 2)            % Here we have to permanently change the image type to RGB
        img = ind2rgb8(img,get(handles.hCallingFig,'ColorMap'));
        handles.origFig = img;      % Update the copy of the original image
        guidata(handles.figure1,handles)
    end
    img = cvlib_mex('floodfill',img,params);
    set(handles.hImage,'CData', img);
    
% -------------------------------------------------------------------------------------
function pushbutton_repaintMultiple_Callback(hObject, eventdata, handles)
    [params,but] = prepareParams(handles, 1);
    if (isempty(params) || but ~= 1),   return;     end
    while (but == 1)
        img = get(handles.hImage,'CData');
        if (ndims(img) == 2)            % Here we have to permanently change the image type to RGB
            img = ind2rgb8(img,get(handles.hCallingFig,'ColorMap'));
            handles.origFig = img;      % Update the copy of the original image
            guidata(handles.figure1,handles)
        end
        img = cvlib_mex('floodfill',img,params);
        set(handles.hImage,'CData', img); 
        [x,y,but] = ginput_pointer(1,'crosshair');  % Get next point
        params.Point = [x y];
    end

% -------------------------------------------------------------------------------------
function [params,but] = prepareParams(handles, opt)
% Prepare the params structure that is to be transmited to cvlib_mex (gets also the point)
% OPT is optional and ment to be used (its enough if ~= []) when painting with Cte color
% BUT returns which button has been pressed.
    if (nargin == 1),   opt = [];   end
    if (~handles.randColor && ~isempty(opt))        % That is, Cte color
        if (~isempty(handles.fillColor))
            params.FillColor = handles.fillColor;
        else
            errordlg('I don''t have yet a filling color. You must select one first.','Error')
            params = [];
            return
        end
    end
    figure(handles.hCallingFig)         % Bring the figure containing image forward
    but = 1;                            % Initialize to a left click
    [x,y,but]  = ginput_pointer(1,'crosshair');
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
    x = localAxes2pix(n,X,x);
    y = localAxes2pix(m,Y,y);

% -------------------------------------------------------------------------------------
function radio_randColor_Callback(hObject, eventdata, handles)
    if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
    set(handles.radio_cteColor,'Value',0)
    handles.randColor = 1;
    set(findobj(handles.figure1,'Style','toggle'),'Enable','Inactive')
    set(handles.pushbutton_moreColors,'Enable','Inactive')
    set(handles.toggle_currColor,'BackgroundColor','w')
    guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function radio_cteColor_Callback(hObject, eventdata, handles)
    if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
    set(handles.radio_randColor,'Value',0)
    handles.randColor = 0;
    set(findobj(handles.figure1,'Style','toggle'),'Enable','on');
    set(handles.pushbutton_moreColors,'Enable','on')
    guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function radio_fourConn_Callback(hObject, eventdata, handles)
    if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
    set(handles.radio_eightConn,'Value',0)
    handles.connect = 4;
    guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function radio_eightConn_Callback(hObject, eventdata, handles)
    if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
    set(handles.radio_fourConn,'Value',0)
    handles.connect = 8;
    guidata(handles.figure1,handles)

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
    guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function pushbutton_moreColors_Callback(hObject, eventdata, handles)
    c = uisetcolor;
    if (length(c) > 1)          % That is, if a color was selected
        handles.fillColor = round(c*255);
        set(handles.toggle_currColor,'BackgroundColor',c)
        guidata(handles.figure1,handles)
    end

% -------------------------------------------------------------------------------------
function pushbutton_pickSingle_Callback(hObject, eventdata, handles)
    [params,but] = prepareParams(handles, []);
    if (isempty(params) || but ~= 1),   return;     end
    img = get(handles.hImage,'CData');              % Get the image
    [dumb,mask] = cvlib_mex('floodfill',img,params);
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
        digitize(handles,mask)
    end

% -------------------------------------------------------------------------------------
function digitize(handles,img)
    % IMG is binary mask. Digitize its outer contours

    B = img_fun('bwboundaries',img,'noholes');
    
    % Remove short polygons, and reduce pts along rows & cols
    j = false(1,length(B));
	for k = 1:length(B)
		if (length(B{k}) < handles.minPts)
            j(k) = 1;
        else
            df = diff(B{k}(:,1));
            id = df == 0;
            id = id & [false; id(1:end-1)];
            if (any(id)),       B{k}(id,:) = [];   end
            df = diff(B{k}(:,2));
            id = df == 0;
            id = id & [false; id(1:end-1)];
            if (any(id)),       B{k}(id,:) = [];   end
        end
    end
    B(j) = [];
    
    x_inc = handles.head(8);    y_inc = handles.head(9);
    x_min = handles.head(1);    y_min = handles.head(3);
    if (handles.head(7))            % Work in grid registration
        x_min = x_min + x_inc/2;
        y_min = y_min + y_inc/2;
    end

    if (handles.single_poly)                    % Draw a single polygon
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
			
        h_edge = line(x, y,'Linewidth',handles.DefLineThick,'Color',handles.DefLineColor, ...
                'Tag','shape_detected','Userdata',handles.udCount);
            
		multi_segs_str = cell(length(h_edge),1);    % Just create a set of empty info strings
		draw_funs(h_edge,'isochron',multi_segs_str);
        handles.udCount = handles.udCount + 1;
    else                                        % Draw separate polygons
		for k = 1:length(B)
            x = (B{k}(:,2)-1) * x_inc + x_min;
            y = (B{k}(:,1)-1) * y_inc + y_min;
            
            h_edge = line(x, y,'Linewidth',handles.DefLineThick,'Color',handles.DefLineColor, ...
                    'Tag','shape_detected');
            draw_funs(h_edge,'line_uicontext')      % Set lines's uicontextmenu
        end    
    end
    guidata(handles.figure1,handles)
    
% -------------------------------------------------------------------------------------
function pushbutton_pickMultiple_Callback(hObject, eventdata, handles)
    % I left the code that deals with processing in Lab & HSV color models but removed
    % those options from the GUI. This is for the case that I change my mind and decide
    % to reintroduce it. For the time beeing, RGB seams to work better.
    [params,but] = prepareParams(handles, 1);
    if (isempty(params) || but ~= 1),   return;     end
    img = get(handles.hImage,'CData');              % Get the image
    nColors = 0;
    mask = false([size(img,1) size(img,2)]);        % Initialize the mask
    while (but == 1)
        nColors = nColors + 1;
        [dumb,mask(:,:,nColors)] = cvlib_mex('floodfill',img,params);
        [x,y,but] = ginput_pointer(1,'crosshair');  % Get next point
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
        a = img(:,:,1);
        b = img(:,:,2);
        c = img(:,:,3);
        cm = zeros(nColors, 6);
        for count = 1:nColors
            cm(count,1) = mean2(a(mask(:,:,count)));
            %cm(count,1) = a(round(params.Point(2)),round(params.Point(1)));
            cm(count,1:2) = [max(cm(count,1) - handles.tol, 0) min(cm(count,1) + handles.tol, 255)];
            cm(count,3) = mean2(b(mask(:,:,count)));
            %cm(count,3) = b(round(params.Point(2)),round(params.Point(1)));
            cm(count,3:4) = [max(cm(count,3) - handles.tol, 0) min(cm(count,3) + handles.tol, 255)];
            cm(count,5) = mean2(c(mask(:,:,count)));
            %cm(count,5) = c(round(params.Point(2)),round(params.Point(1)));
            cm(count,5:6) = [max(cm(count,5) - handles.tol, 0) min(cm(count,5) + handles.tol, 255)];
            tmp = (a >= cm(count,1) & a <= cm(count,2) & b >= cm(count,3) & b <= cm(count,4) & ...
                c >= cm(count,5) & c <= cm(count,6));
            tmp  = img_fun('bwmorph',tmp,'clean');      % Get rid of isolated pixels
            if (get(handles.checkbox_useDilation,'Value'))
                tmp  = img_fun('bwmorph',tmp,'dilate');
                %tmp  = cvlib_mex('dilate',tmp);
                %tmp  = cvlib_mex('morphologyex',tmp,'close');
            end
            if (handles.colorSegment)           % Do color segmentation. One figure for each color
                tmp = repmat(tmp,[1 1 3]);
                img(~tmp) = handles.bg_color;
                if (handles.IAmAMirone)
                    if (handles.image_type == 2 || handles.image_type == 20)
                        h = mirone(img);
                        set(h,'Name','Color segmentation')
                    else
                        tmp.X = handles.head(1:2);  tmp.Y = handles.head(3:4);  tmp.head = handles.head;
                        tmp.name = ['Color segmentation nº ' num2str(1)];
                        mirone(img,tmp);
                    end
                else
        	        figure; image(img);
                end
            else                                % Create contours from mask
                digitize(handles,tmp)
            end
        end
    elseif (isempty(handles.colorModel) && ndims(img) == 2)        % That is, indexed image & RGB model
        cm = zeros(nColors, 2);
        for count = 1:nColors
            cm(count,1) = mean2(img(mask(:,:,count)));
            %cm(count,1) = img(round(params.Point(2)),round(params.Point(1)));
            cm(count,1:2) = [max(cm(count,1) - handles.tol, 0) min(cm(count,1) + handles.tol, 255)];
            tmp = (img >= cm(count,1) & img <= cm(count,2));
            tmp  = img_fun('bwmorph',tmp,'clean');      % Get rid of isolated pixels
            if (get(handles.checkbox_useDilation,'Value'))
                tmp  = img_fun('bwmorph',tmp,'dilate');
            end
            if (handles.colorSegment)           % Do color segmentation. One figure for each color
                img(~tmp) = handles.bg_color;
                if (handles.IAmAMirone)
                    if (handles.image_type == 2 || handles.image_type == 20)
                        h = mirone(img);
                        set(h,'ColorMap',get(handles.hCallingFig,'ColorMap'),'Name','Color segmentation')
                    else
                        tmp.X = handles.head(1:2);  tmp.Y = handles.head(3:4);  tmp.head = handles.head;
                        tmp.name = ['Color segmentation nº ' num2str(1)];
                        tmp.cmap = get(handles.hCallingFig,'ColorMap');
                        mirone(img,tmp);
                    end
                else
        	        h = figure;     image(img);
                    set(h,'ColorMap',get(handles.hCallingFig,'ColorMap'),'Name','Color segmentation')
                end
            else                                % Create contours from mask
                digitize(handles,tmp)
            end
        end
    else                                        % Either YCrCb or Lab model were used
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
                        tmp.name = ['Color segmentation nº ' num2str(1)];
                        mirone(img,tmp);
                    end
                else
        	        figure; image(img);
                end
            else                                % Create contours from mask
                digitize(handles,tmp)
            end
        end    
    end
    
% -------------------------------------------------------------------------------------
% function radio_YCrCb_Callback(hObject, eventdata, handles)
%     if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
%     set(handles.radio_Lab,'Value',0)
%     set(handles.radio_RGB,'Value',0)
%     handles.colorModel = 'YCrCb';
%     guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
% function radio_Lab_Callback(hObject, eventdata, handles)
%     if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
%     set(handles.radio_YCrCb,'Value',0)
%     set(handles.radio_RGB,'Value',0)
%     handles.colorModel = 'lab';
%     guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
% function radio_RGB_Callback(hObject, eventdata, handles)
%     % This the only one know (it seams to do better than the others)
%     if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
%     set(handles.radio_YCrCb,'Value',0)
%     set(handles.radio_Lab,'Value',0)
%     handles.colorModel = '';
%     guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function radio_colorSegment_Callback(hObject, eventdata, handles)
    if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
    set(handles.radio_digitize,'Value',0)
    handles.colorSegment = 1;
    set(handles.edit_minPts,'Visible','off')
    set(handles.text_minPts,'Visible','off')
    guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function radio_digitize_Callback(hObject, eventdata, handles)
    if (~get(hObject,'Value')),      set(hObject,'Value',1);   return;     end
    set(handles.radio_colorSegment,'Value',0)
    handles.colorSegment = 0;
    set(handles.edit_minPts,'Visible','on')
    set(handles.text_minPts,'Visible','on')
    guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function edit_minPts_Callback(hObject, eventdata, handles)
    xx = round( str2double(get(hObject,'String')) );
    if (isnan(xx)),     set(hObject,'String','50');     return;     end
    handles.minPts = xx;
    guidata(handles.figure1,handles)


% -------------------------------------------------------------------------------------
function push_restoreImg_Callback(hObject, eventdata, handles)
    set(handles.hImage,'CData',handles.origFig)

% -------------------------------------------------------------------------------------
function pixelx = localAxes2pix(dim, x, axesx)
%   Convert axes coordinates to pixel coordinates.
%   PIXELX = AXES2PIX(DIM, X, AXESX) converts axes coordinates
%   (as returned by get(gca, 'CurrentPoint'), for example) into
%   pixel coordinates.  X should be the vector returned by
%   X = get(image_handle, 'XData') (or 'YData').  DIM is the
%   number of image columns for the x coordinate, or the number
%   of image rows for the y coordinate.
	xfirst = x(1);
	xlast = x(max(size(x)));
	
	if (dim == 1)
        pixelx = axesx - xfirst + 1;
        return;
	end
	xslope = (dim - 1) / (xlast - xfirst);
	if ((xslope == 1) & (xfirst == 1))
        pixelx = axesx;
	else
        pixelx = xslope * (axesx - xfirst) + 1;
	end

% -------------------------------------------------------------------------------------
function y = mean2(x)
%MEAN2 Compute mean of matrix elements.
y = sum(x(:)) / numel(x);

% -------------------------------------------------------------------------------------
function checkbox_useDilation_Callback(hObject, eventdata, handles)
% It does nothing


% --- Creates and returns a handle to the GUI figure. 
function floodFill_LayoutFcn(h1,handles);

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Flood Fill',...
'NumberTitle','off',...
'Position',[520 390 185 395],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

h2 = uicontrol('Parent',h1,...
'Position',[4 35 177 110],...
'Style','frame','Tag','frame2');

uicontrol('Parent',h1,...
'Position',[4 222 176 165],...
'Style','frame','Tag','frame1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@floodFill_uicallback,h1,'slider_tolerance_Callback'},...
'Max',255,...
'Position',[4 181 176 16],...
'Style','slider',...
'SliderStep',[0.00390625 0.1],...
'TooltipString','Color detection equal pixel color +/- tolerance',...
'Tag','slider_tolerance');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'pushbutton_pickSingle_Callback'},...
'Position',[12 113 160 23],...
'String','Pick single shape',...
'TooltipString','Pick up the body''s shape with the selected color',...
'Tag','pushbutton_pickSingle');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'pushbutton_repaintSingle_Callback'},...
'Position',[12 258 160 23],...
'String','Paint once',...
'TooltipString','Repaint single body',...
'Tag','pushbutton_repaintSingle');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'radio_fourConn_Callback'},...
'Position',[12 302 90 15],...
'String','4 connectivity',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_fourConn');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'radio_eightConn_Callback'},...
'Position',[12 286 90 15],...
'String','8 connectivity',...
'Style','radiobutton',...
'Tag','radio_eightConn');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0.501960784313725 0.501960784313725],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[13 340 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle11');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0 0],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[29 340 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle12');

uicontrol('Parent',h1,...
'BackgroundColor',[0.501960784313725 0.250980392156863 0.250980392156863],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[45 340 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle13');

uicontrol('Parent',h1,...
'BackgroundColor',[0.250980392156863 0 0],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[61 340 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle14');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 0.501960784313725],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[77 340 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle15');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0.501960784313725 0],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[93 340 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle16');

uicontrol('Parent',h1,...
'BackgroundColor',[0.501960784313725 0.501960784313725 0],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[109 340 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle17');

uicontrol('Parent',h1,...
'BackgroundColor',[0 1 0],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[125 340 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle18');

uicontrol('Parent',h1,...
'BackgroundColor',[0 0.501960784313725 0],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[13 324 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle21');

uicontrol('Parent',h1,...
'BackgroundColor',[0 0.501960784313725 0.501960784313725],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[29 324 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle22');

uicontrol('Parent',h1,...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[45 324 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle23');

uicontrol('Parent',h1,...
'BackgroundColor',[0 0 1],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[61 324 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle24');

uicontrol('Parent',h1,...
'BackgroundColor',[0 1 1],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[77 324 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle25');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0.501960784313725 0.752941176470588],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[93 324 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle26');

uicontrol('Parent',h1,...
'BackgroundColor',[0.250980392156863 0 0.501960784313725],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[109 324 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle27');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@floodFill_uicallback,h1,'toggle00_Callback'},...
'Enable','inactive',...
'Position',[125 324 15 15],...
'Style','togglebutton',...
'Value',1,...
'Tag','toggle28');

uicontrol('Parent',h1,...
'FontSize',9,...
'Position',[54 379 60 15],...
'String','Flood fill',...
'Style','text','Tag','text1');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'radio_randColor_Callback'},...
'Position',[12 363 90 15],...
'String','Random colors',...
'Style','radiobutton',...
'TooltipString','Repainting will use a randomly selected color',...
'Value',1,...
'Tag','radio_randColor');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'radio_cteColor_Callback'},...
'Position',[115 363 61 15],...
'String','Cte color',...
'Style','radiobutton',...
'TooltipString','Repainting will use a selected color',...
'Tag','radio_cteColor');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'pushbutton_moreColors_Callback'},...
'Enable','inactive',...
'Position',[143 338 30 19],...
'String','More',...
'TooltipString','Chose other color',...
'Tag','pushbutton_moreColors');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'Position',[145 324 25 15],...
'Style','togglebutton',...
'TooltipString','Current selected color',...
'Value',1,...
'HandleVisibility','off',...
'Tag','toggle_currColor');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'pushbutton_repaintMultiple_Callback'},...
'Position',[12 230 160 23],...
'String','Paint multiple',...
'TooltipString','Each mouse click will repaint the hit shape',...
'Tag','pushbutton_repaintMultiple');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'pushbutton_pickMultiple_Callback'},...
'Position',[12 81 160 23],...
'String','Pick multiple shapes',...
'TooltipString','Find out all bounding polygons that share the selected color',...
'Tag','pushbutton_pickMultiple');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'radio_colorSegment_Callback'},...
'Position',[12 59 115 15],...
'String','Color segmentation',...
'Style','radiobutton',...
'TooltipString','Create a separate figure with shapes colored as you selected',...
'Value',1,...
'Tag','radio_colorSegment');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'radio_digitize_Callback'},...
'Position',[12 42 55 15],...
'String','Digitize',...
'Style','radiobutton',...
'TooltipString','Detect bounding polygon to the colored selected body',...
'Tag','radio_digitize');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@floodFill_uicallback,h1,'edit_minPts_Callback'},...
'Position',[140 40 33 20],...
'String','50',...
'Style','edit',...
'TooltipString','Shapes with less than this number of vertex won''t be ploted',...
'Tag','edit_minPts',...
'Visible','off');

uicontrol('Parent',h1,...
'Position',[118 43 20 15],...
'String','Min',...
'Style','text',...
'Tag','text_minPts',...
'Visible','off');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[57 197 90 15],...
'String','Tolerance = 20',...
'Style','text','Tag','text_tol');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'push_restoreImg_Callback'},...
'Position',[37 5 111 23],...
'String','Restore image',...
'Tag','push_restoreImg');

uicontrol('Parent',h1,...
'Callback',{@floodFill_uicallback,h1,'checkbox_useDilation_Callback'},...
'Position',[4 158 80 15],...
'String','Use Dilation',...
'Style','checkbox',...
'TooltipString','Use dilation operation to find better limits between neighboring shapes',...
'Value',1,...
'Tag','checkbox_useDilation');


function floodFill_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

function varargout = vitrinite(varargin)
% M-File changed by desGUIDE 
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vitrinite_export (see VARARGIN)

hObject = figure('Tag','figure1','Visible','off');
handles = guihandles(hObject);
guidata(hObject, handles);
vitrinite_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

movegui(handles.figure1,'north')

% Import icons
load (['data' filesep 'mirone_icons.mat'],'Mfopen_ico','zoom_ico');

h_toolbar = uitoolbar('parent',handles.figure1, 'BusyAction','queue','HandleVisibility','on',...
   'Interruptible','on','Tag','FigureToolBar','Visible','on');
uipushtool('parent',h_toolbar,'Click',@clicked_loadMarkersImg_Callback, ...
   'cdata',Mfopen_ico,'TooltipString','Open Marker images');
uipushtool('parent',h_toolbar,'Click',@clicked_loadSampleImg_Callback, ...
   'cdata',Mfopen_ico,'TooltipString','Open Sample image','Sep','on');
uitoggletool('parent',h_toolbar,'Click','zoom', ...
   'cdata',zoom_ico,'TooltipString','Zooming on/off','Sep','on');


handles.hImgMarker = [];
handles.hImgSample = [];
handles.countMark = 0;         % Counter of number of loaded image markers
handles.tiePoints = [];
handles.countTiePoints = 0;
handles.haveStatusBar = 0;
handles.IamCallibrated = 0;
handles.geog = 0;

% Choose default command line output for vitrinite_export
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes vitrinite_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = vitrinite_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = vitrinite_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function clicked_loadMarkersImg_Callback(hObject, eventdata)
    handles = guidata(hObject);
	[FileName,PathName] = uigetfile({ ...
        '*.jpg', 'JPEG image (*.jpg)'; ...
        '*.png', 'Portable Network Graphics(*.png)'; ...
        '*.bmp', 'Windows Bitmap (*.bmp)'; ...
        '*.gif', 'GIF image (*.gif)'; ...
        '*.pcx', 'Windows Paintbrush (*.pcx)'; ...
        '*.ras', 'SUN rasterfile (*.ras)'; ...
        '*.tif', 'Tagged Image File (*.tif)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Select image format');
	if isequal(FileName,0);     return;     end
    fname = [PathName FileName];
	try             % Use a try because if the name was given via edit box it may be wrong
        [img,cmap] = imread(fname);
        handles.ImgMarker{handles.countMark+1} = img;    % Make a copy
        handles.countMark = handles.countMark + 1;
	catch
        errordlg(['Error: -> ' fname ' does not exist or is not a valid image file'],'Error')
        return
	end

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
        set(handles.popup_markersImgs,'String',str,'Value',handles.countMark)    % Set it the visible one
    end
    
    % Try to fish the reflectance factor from image name
    sep = findstr(FNAME,'_');
    if (~isempty(sep))
        reflectance = str2double(FNAME(sep(end)+1:end));
        if (~isnan(reflectance) && reflectance > 0.1 && reflectance < 10)
            handles.reflectance{handles.countMark} = reflectance;               % Yes, we got it
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
	
	if (handles.countMark >= 1 && ~isempty(handles.hImgSample))
        if (~handles.haveStatusBar)
            createStatusBar(handles);            handles.haveStatusBar = 1;
        end
	end
	guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function clicked_loadSampleImg_Callback(hObject, eventdata)
    handles = guidata(hObject);
	[FileName,PathName] = uigetfile({ ...
        '*.jpg', 'JPEG image (*.jpg)'; ...
        '*.png', 'Portable Network Graphics(*.png)'; ...
        '*.bmp', 'Windows Bitmap (*.bmp)'; ...
        '*.gif', 'GIF image (*.gif)'; ...
        '*.pcx', 'Windows Paintbrush (*.pcx)'; ...
        '*.ras', 'SUN rasterfile (*.ras)'; ...
        '*.tif', 'Tagged Image File (*.tif)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Select image format');
	if isequal(FileName,0);     return;     end
    fname = [PathName FileName];
	try             % Use a try because if the name was given via edit box it may be wrong
        [handles.ImgSample,cmap] = imread(fname);
	catch
        errordlg(['Error: -> ' fname ' does not exist or is not a valid image file'],'Error')
        return
	end

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
	
    if (~handles.IamCallibrated)
	    if (handles.countMark >= 1 && ~handles.haveStatusBar)
            createStatusBar(handles);            handles.haveStatusBar = 1;
        end
    else        % A previous callibrated image already exists. Need to callibrate this one
       pushbutton_callibrate_Callback(handles.pushbutton_callibrate, [], handles) 
    end
	guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function popup_markersImgs_Callback(hObject, eventdata, handles)
    if (handles.countMark < 2),    return;     end      % Too soon
    set(handles.edit_reflectance,'String',handles.reflectance{get(hObject,'Value')})
    resetImg(handles,get(hObject,'Value'))

% --------------------------------------------------------------------
function resetImg(handles, number)
    % Update image merker accordingly to what was selected in the popup
	set(handles.hImgMarker,'CData',handles.ImgMarker{number})
	set(handles.axes2,'PlotBoxAspectRatio',handles.imgAspect{number})    

% --------------------------------------------------------------------
function pushbutton_callibrate_Callback(hObject, eventdata, handles)
    % 
    if (isempty(handles.hImgSample))
        errordlg('Calibrate what? Your eyes? Load the sample image.','ERROR');return
    end
    if (handles.countTiePoints < 2)                   % Not yet possible to calibrate
        errordlg('No I won''t. I need at least two tie points to do that.','ERROR');return
    end
    
    x = sort(cell2mat(handles.tiePoints));            % Marker pixel values
    Y = sort(cell2mat(handles.reflectance));          % Marker reflectance values
%     x = [0.0 x];                                  % Don't let extrapolation pass to negative values
%     Y = [0.15 Y];
    x = x;
    Y = Y;
    Y_bak = Y;
    
    yi = single(interp1(x,Y,0:255,'linear','extrap'));
    Z = yi(handles.ImgSample);
    X = 1:size(handles.ImgSample,2);
    Y = 1:size(handles.ImgSample,1);
    head = [1 X(end) 1 Y(end) 0 255 0 1 1];
    setappdata(handles.figure1,'dem_z',Z);  setappdata(handles.figure1,'dem_x',X);
    setappdata(handles.figure1,'dem_y',Y);  setappdata(handles.figure1,'GMThead',head);
    set(handles.pushbutton_getAvgReflec,'Enable','on')
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
function pushbutton_getTiePoint_Callback(hObject, eventdata, handles)
    % Get a tie point from the current Marker image
    [x,y,but]  = ginput_pointer(1,'crosshair');
    if (but ~= 1),   return;     end
    params.Point = [x y];    params.Tolerance = 10;    params.Connect = 4;
    
    img = get(handles.hImgMarker,'CData');                   % Get the image
    [dumb,mask] = cvlib_mex('floodfill',img,params);
    
    % TESTAR SE TENHO REFLECTANCE PARA ESTE PONTO
    tiePoint = round(mean2(img(mask)));
    whichMarker = get(handles.popup_markersImgs,'Value');   % We need to know which tie point is this
    handles.tiePoints{whichMarker} = tiePoint;
    handles.sigma{whichMarker} = std2(img(mask));
    handles.countTiePoints = handles.countTiePoints + 1;
	guidata(handles.figure1,handles)

%--------------------------------------------------------------------------
function pushbutton_getAvgReflec_Callback(hObject, eventdata, handles)
    % Get the average reflectance of the clicked shape
    if (~handles.IamCallibrated)
        errordlg('You need to callibrate the sample image first.','ERROR'); return
    end
    [x,y,but]  = ginput_pointer(1,'crosshair');
    if (but ~= 1),   return;     end
    params.Point = [x y];    params.Tolerance = 5;    params.Connect = 8;
    hAx = get(handles.figure1,'CurrentAxes');
    img = get(findobj(hAx,'Type','Image'),'CData');                   % Get the image
    
    [dumb,mask] = cvlib_mex('floodfill',img,params);
    avgPix = round(mean2(img(mask)));
    avgReflect = double(handles.callibCurv(avgPix));
    hText = text(x,y,sprintf('%.2f',avgReflect),'Fontsize',8,'Parent',hAx,'HorizontalAlignment','center');

    % Set a uicontext with the "Deleting" option
    cmenuHand = uicontextmenu;      set(hText, 'UIContextMenu', cmenuHand);
    uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');

%--------------------------------------------------------------------------
function edit_reflectance_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of edit_reflectance as text
%        str2double(get(hObject,'String')) returns contents of edit_reflectance as a double

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
function vitrinite_LayoutFcn(h1,handles);

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Vitrinite',...
'NumberTitle','off',...
'Position',[520 268 780 532],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
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

h12 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@vitrinite_uicallback,h1,'popup_markersImgs_Callback'},...
'Position',[10 134 241 22],...
'Style','popupmenu',...
'String',{''},...
'TooltipString','Loaded Marker Images',...
'Value',1,...
'Tag','popup_markersImgs');

h13 = uicontrol('Parent',h1,...
'Callback',{@vitrinite_uicallback,h1,'pushbutton_callibrate_Callback'},...
'FontSize',9,...
'Position',[520 132 91 23],...
'String','Calibrate',...
'Tag','pushbutton_callibrate');

h14 = uicontrol('Parent',h1,...
'Callback',{@vitrinite_uicallback,h1,'pushbutton_getTiePoint_Callback'},...
'FontSize',9,...
'Position',[410 132 91 23],...
'String','Get Tie Point',...
'Tag','pushbutton_getTiePoint');

h15 = uicontrol('Parent',h1,...
'Callback',{@vitrinite_uicallback,h1,'pushbutton_getAvgReflec_Callback'},...
'Enable','inactive',...
'FontSize',9,...
'Position',[634 132 137 24],...
'String','Get average reflectance',...
'Tag','pushbutton_getAvgReflec');

h16 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@vitrinite_uicallback,h1,'edit_reflectance_Callback'},...
'Position',[308 134 71 21],...
'Style','edit',...
'TooltipString','Reflectance of current marker.',...
'Tag','edit_reflectance');

h17 = uicontrol('Parent',h1,...
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

h23 = uicontrol('Parent',h1,...
'FontSize',9,...
'Position',[10 158 131 16],...
'String','Current Marker Image',...
'Style','text',...
'Tag','text4');

h24 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[23 101 108 16],...
'String','Callibration Curve',...
'Style','text',...
'Tag','text5');

h25 = uicontrol('Parent',h1,...
'FontSize',9,...
'Position',[138 515 131 16],...
'String','Marker Image',...
'Style','text',...
'Tag','text6');

h26 = uicontrol('Parent',h1,...
'FontSize',9,...
'Position',[510 515 131 16],...
'String','Sample Image',...
'Style','text',...
'Tag','text7');

function vitrinite_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

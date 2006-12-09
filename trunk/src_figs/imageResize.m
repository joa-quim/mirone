function varargout = imageResize(varargin)
% M-File changed by desGUIDE 

hObject = figure('Tag','figure1','Visible','off');
handles = guihandles(hObject);
guidata(hObject, handles);
imageResize_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

if (~isempty(varargin))
    handles.hCallingFig = varargin{1};
else
    delete(hObject)
    return
end

% Get the Mirone handles. We need it here
handlesMir = guidata(handles.hCallingFig);
if (handlesMir.no_file)
    errordlg('You didn''t even load a file. What are you expecting then?','Error')
    delete(hObject);    return
end

plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
plugedWin = [plugedWin hObject];            % Add this figure handle to the carraças list
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
handles.hImage = findobj(handles.hCallingFig,'Type','image');
handles.imgSize   = size(get(handles.hImage,'CData'));

handles.head = handlesMir.head;
handles.DefLineThick = handlesMir.DefLineThick;
handles.DefLineColor = handlesMir.DefLineColor;
handles.geog = handlesMir.geog;
handles.image_type = handlesMir.image_type;
if (numel(handles.imgSize) == 2)
    handles.parentCmap = get(handlesMir.figure1,'Colormap');
end

handles.pixWidth  = handles.imgSize(2);
handles.pixHeight = handles.imgSize(1);
handles.resolution = 72;
handles.unitFact  = 2.54;
handles.resolutionFact = 1 / handles.resolution * handles.unitFact;
handles.constrainProp = 1;
handles.intepMethod = 'nearest';
handles.isPercent = 0;

% Fill in the necessary uicontrols
set(handles.edit_pixHeight,'String',handles.imgSize(1))
set(handles.edit_pixWidth,'String',handles.imgSize(2))
set(handles.edit_docResolution,'String',handles.resolution)
handles = pix2size(handles,'w');           % Fill the "Document" edits
handles = pix2size(handles,'h');

%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
bgcolor = get(0,'DefaultUicontrolBackgroundColor');
framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
%set(0,'Units','pixels');    set(hObject,'Units','pixels')    % Pixels are easier to reason with
h_f = findobj(hObject,'Style','Frame');
for i=1:length(h_f)
    frame_size = get(h_f(i),'Position');
    f_bgc = get(h_f(i),'BackgroundColor');
    usr_d = get(h_f(i),'UserData');
    if abs(f_bgc(1)-bgcolor(1)) > 0.01           % When the frame's background color is not the default's
        frame3D(hObject,frame_size,framecolor,f_bgc,usr_d)
    else
        frame3D(hObject,frame_size,framecolor,'',usr_d)
        delete(h_f(i))
    end
end
% Recopy the text fields on top of previously created frames (uistack is to slow)
h_t = [handles.text_docSize handles.text_pixDim];
for i=1:length(h_t)
    usr_d = get(h_t(i),'UserData');
    t_size = get(h_t(i),'Position');   t_str = get(h_t(i),'String');    fw = get(h_t(i),'FontWeight');
    bgc = get (h_t(i),'BackgroundColor');   fgc = get (h_t(i),'ForegroundColor');
    t_just = get(h_t(i),'HorizontalAlignment');     t_tag = get (h_t(i),'Tag');
    h = uicontrol('Parent',hObject, 'Style','text', 'Position',t_size,'String',t_str,'Tag',t_tag,...
        'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw,...
        'UserData',usr_d,'HorizontalAlignment',t_just);
end
handles.text_pixDim = h;        % We only need this one
%------------- END Pro look (3D) -------------------------------------------------------

handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes imageResize_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% ----------------------------------------------------------------------------------
set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = imageResize_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = imageResize_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ----------------------------------------------------------------------------------
function hand = pix2size(handles,opt)
    % Convert from pixels to "Document Size" unities
    if ( opt == 'w' )       % Width
        handles.docWidth = handles.pixWidth * handles.resolutionFact;
        set(handles.edit_docWidth,'String',handles.docWidth)
    else                    % Height
        handles.docHeight = handles.pixHeight * handles.resolutionFact;
        set(handles.edit_docHeight,'String',handles.docHeight)
    end
    guidata(handles.figure1, handles);
    if (nargout),   hand = handles;     end

% ----------------------------------------------------------------------------------
function size2pix(handles,opt)
    % Convert from "Document Size" to pixels unities
    if ( opt == 'w' )       % Width
        handles.pixWidth = handles.docWidth / handles.resolutionFact;
        set(handles.edit_pixWidth,'String',round(handles.pixWidth))
    else                    % Height
        handles.pixHeight = handles.docHeight / handles.resolutionFact;
        set(handles.edit_pixHeight,'String',round(handles.pixHeight))
    end
    guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------------
function edit_pixWidth_Callback(hObject, eventdata, handles)
    xx = round( str2double(get(hObject,'String')) );        p = xx / 100;
    if (isnan(xx)),     set(hObject,'String',handles.pixWidth);    return;      end
    if (handles.isPercent),     xx = handles.imgSize(2) * p;       end
    pixWidthOld = handles.pixWidth;
    handles.pixWidth  = xx;
    pix2size(handles,'w')           % Update the "Document" edit
    if (handles.constrainProp)      % Recompute height
        if (~handles.isPercent)
            handles.pixHeight = handles.pixHeight * handles.pixWidth / pixWidthOld;
            set(handles.edit_pixHeight,'String',round(handles.pixHeight))
        else
            handles.pixHeight = handles.imgSize(1) * p;
            set(handles.edit_pixHeight,'String',round(p*100))
        end
        pix2size(handles,'h')       % Update the "Document" edit
    end
    guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------------
function edit_pixHeight_Callback(hObject, eventdata, handles)
    xx = round( str2double(get(hObject,'String')) );        p = xx / 100;
    if (handles.isPercent),     xx = round(handles.imgSize(1) * xx / 100);      end
    if (isnan(xx)),     set(hObject,'String',handles.pixHeight);    return;     end
    pixHeightOld = handles.pixHeight;
    handles.pixHeight  = xx;
    pix2size(handles,'h')           % Update the "Document" edit
    if (handles.constrainProp)      % Recompute height
        if (~handles.isPercent)
            handles.pixWidth = handles.pixWidth * handles.pixHeight / pixHeightOld;
            set(handles.edit_pixWidth,'String',round(handles.pixWidth))
        else
            handles.pixWidth = handles.imgSize(2) * p;
            set(handles.edit_pixWidth,'String',round(p*100))
        end
        pix2size(handles,'w')       % Update the "Document" edit
    end
    guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------------
function edit_docWidth_Callback(hObject, eventdata, handles)
    xx = round( str2double(get(hObject,'String')) );
    if (isnan(xx)),     set(hObject,'String',handles.docWidth);    return;     end
    docWidthOld = handles.docWidth;
    handles.docWidth  = xx;
    size2pix(handles,'w')           % Update the "Document" edit
    if (handles.constrainProp)      % Recompute height
        handles.docHeight = handles.docHeight * handles.docWidth / docWidthOld;
        set(handles.edit_docHeight,'String',handles.docHeight)
        size2pix(handles,'h')       % Update the "Document" edit
    end

% ----------------------------------------------------------------------------------
function edit_docHeight_Callback(hObject, eventdata, handles)
    xx = round( str2double(get(hObject,'String')) );
    if (isnan(xx)),     set(hObject,'String',handles.docHeight);    return;     end
    docHeightOld = handles.docHeight;
    handles.docHeight  = xx;
    size2pix(handles,'h')           % Update the "Pixel" edit
    if (handles.constrainProp)      % Recompute height
        handles.docWidth = handles.docWidth * handles.docHeight / docHeightOld;
        set(handles.edit_docWidth,'String',handles.docWidth)
        size2pix(handles,'w')       % Update the "Pixel" edit
    end

% ----------------------------------------------------------------------------------
function edit_docResolution_Callback(hObject, eventdata, handles)
    % Note that even if the resolution is set to pixels/cm we store in pixels/inch
    xx = str2double(get(hObject,'String'));
    if (isnan(xx)),     set(hObject,'String',handles.resolution);    return;     end
    if (get(handles.popup_docResolution,'Value') == 1)     % pixel/inch
        handles.resolution  = round(xx);
    else
        handles.resolution  = round(xx*2.54);
    end
    handles.resolutionFact = 1 / handles.resolution * handles.unitFact;
    guidata(handles.figure1, handles);
    size2pix(handles,'w')           % Update the "Pixel" edits
    size2pix(handles,'h')

% ----------------------------------------------------------------------------------
function popup_pixWidth_Callback(hObject, eventdata, handles)
    % Change between pixels and percent, but only on the display
    val = get(hObject,'Value');
    switch val
        case 1
            w = handles.pixWidth;            h = handles.pixHeight;
            handles.isPercent = 0;
        case 2
            w = handles.pixWidth / handles.imgSize(2) * 100;
            h = handles.pixHeight / handles.imgSize(1) * 100;
            handles.isPercent = 1;
    end
    set(handles.edit_pixWidth,'String',round(w))
    set(handles.edit_pixHeight,'String',round(h))
    set(handles.popup_pixHeight,'Value',val)
    guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------------
function popup_pixHeight_Callback(hObject, eventdata, handles)
    % Let the popup_pixWidth do the work
    set(handles.popup_pixWidth,'Value',get(hObject,'Value'))
    popup_pixWidth_Callback(handles.popup_pixWidth, [], handles)

% ----------------------------------------------------------------------------------
function popup_docWidth_Callback(hObject, eventdata, handles)
    switch get(hObject,'Value')
        case 1,     handles.unitFact  = 2.54;
        case 2,     handles.unitFact  = 254;
        case 3,     handles.unitFact  = 1;
        case 4,     handles.unitFact  = 1 / 72;
    end
    set(handles.popup_docHeight,'Value',get(hObject,'Value'))
    handles.resolutionFact = 1 / handles.resolution * handles.unitFact;
    guidata(handles.figure1,handles)
    pix2size(handles,'w')
    pix2size(handles,'h')
    
% ----------------------------------------------------------------------------------
function popup_docHeight_Callback(hObject, eventdata, handles)
    % Let the popup_docWidth do the work
    set(handles.popup_docWidth,'Value',get(hObject,'Value'))
    popup_pixWidth_Callback(handles.popup_docWidth, [], handles)

% ----------------------------------------------------------------------------------
function popup_docResolution_Callback(hObject, eventdata, handles)
    % In fact we never change the resolution unites out of DPI
    switch get(hObject,'Value')
        case 1,     resolution = handles.resolution;
        case 2,     resolution = handles.resolution / 2.54;
    end
    set(handles.edit_docResolution,'String',resolution)

% ----------------------------------------------------------------------------------
function checkbox_constProportions_Callback(hObject, eventdata, handles)
    if (get(hObject,'Value'))
        handles.constrainProp = 1;
    else
        handles.constrainProp = 0;
    end
    guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------------
function popup_resampMethod_Callback(hObject, eventdata, handles)
    % Change the interpolation method
    switch get(hObject,'Value')
        case 1,    handles.intepMethod = 'nearest';
        case 2,    handles.intepMethod = 'bilinear';
        case 3,    handles.intepMethod = 'bicubic';
        case 4,    handles.intepMethod = 'area';
    end
    guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
    % Do the job
    try
        img = get(handles.hImage,'CData');
    catch
        errordlg('Figure no longer exists. Why did you kill it?','Error');  return
    end
    set(handles.figure1,'pointer','watch');    set(handles.hCallingFig,'pointer','watch')
    %img = img_fun('imresize',img,round([handles.pixHeight handles.pixWidth]),handles.intepMethod);
    img = cvlib_mex('resize',img,round([handles.pixHeight handles.pixWidth]),handles.intepMethod);
    set(handles.figure1,'pointer','arrow');    set(handles.hCallingFig,'pointer','arrow')
    
    if (handles.image_type == 2)
        if (ndims(img) == 2)
            setappdata(0,'CropedColormap',handles.parentCmap);
        end
        h = mirone(img);
        set(h, 'Name', 'Resized Image')
    else
        tmp.head = [handles.head(1:6) 0];
        tmp.head(8) = diff(handles.head(1:2)) / (round(handles.pixWidth) - ~handles.head(7));
        tmp.head(9) = diff(handles.head(3:4)) / (round(handles.pixHeight) - ~handles.head(7));
        tmp.X = handles.head(1:2);      tmp.Y = handles.head(3:4);
        tmp.name = 'Resized Image';
        if (ndims(img) == 2)
            tmp.cmap = handles.parentCmap;
        end
        mirone(img,tmp);
    end
    
% ----------------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
    delete(handles.figure1)

% --- Creates and returns a handle to the GUI figure. 
function imageResize_LayoutFcn(h1,handles);

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',[0.831372549019608 0.815686274509804 0.784313725490196],...
'MenuBar','none',...
'Name','Resize Image',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',[20.98404194812 29.67743169791],...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',[520 539 328 261],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');


uicontrol('Parent',h1,'Position',[7 177 241 75],'Style','frame','Tag','frame1');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@imageResize_uicallback,h1,'edit_pixWidth_Callback'},...
'Position',[79 214 71 21],...
'Style','edit',...
'TooltipString','Describes the width in pixels','Tag','edit_pixWidth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@imageResize_uicallback,h1,'edit_pixHeight_Callback'},...
'Position',[79 187 71 21],...
'Style','edit',...
'TooltipString','Describes the height in pixels','Tag','edit_pixHeight');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[29 219 41 15],...
'String','Width',...
'Style','text',...
'TooltipString','Describes the width in pixels','Tag','text1');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[29 190 41 18],...
'String','Height',...
'Style','text',...
'TooltipString','Describes the height in pixels','Tag','text2');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@imageResize_uicallback,h1,'popup_pixWidth_Callback'},...
'Position',[158 214 70 22],...
'String',{  'pixels'; 'percent' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_pixWidth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@imageResize_uicallback,h1,'popup_pixHeight_Callback'},...
'Position',[158 187 70 22],...
'String',{  'pixels'; 'percent' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_pixHeight');

uicontrol('Parent',h1,'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[28 242 181 17],...
'String',' Pixel Dimensions:',...
'Style','text','Tag','text_pixDim');

uicontrol('Parent',h1,'Position',[7 60 241 100],'Style','frame','Tag','frame2');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@imageResize_uicallback,h1,'edit_docWidth_Callback'},...
'Position',[79 122 71 21],...
'Style','edit',...
'TooltipString','Set the document width','Tag','edit_docWidth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@imageResize_uicallback,h1,'edit_docHeight_Callback'},...
'Position',[79 95 71 21],...
'Style','edit',...
'TooltipString','Set the document height','Tag','edit_docHeight');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[30 127 41 15],...
'String','Width',...
'Style','text',...
'TooltipString','Set the document width','Tag','text4');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[30 98 41 18],...
'String','Height',...
'Style','text',...
'TooltipString','Set the document height','Tag','text5');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@imageResize_uicallback,h1,'popup_docWidth_Callback'},...
'Position',[159 122 81 22],...
'String',{  'cm'; 'mm'; 'inch'; 'points' },...
'Style','popupmenu',...
'Value',1,'Tag','popup_docWidth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@imageResize_uicallback,h1,'popup_docHeight_Callback'},...
'Position',[159 95 81 22],...
'String',{  'cm'; 'mm'; 'inch'; 'points' },...
'Style','popupmenu',...
'Value',1,'Tag','popup_docHeight');

uicontrol('Parent',h1,'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[29 150 101 17],...
'String',' Document Size:',...
'Style','text','Tag','text_docSize');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@imageResize_uicallback,h1,'edit_docResolution_Callback'},...
'Position',[79 69 71 21],...
'Style','edit',...
'TooltipString','Set the document resolution','Tag','edit_docResolution');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[13 71 65 18],...
'String','Resolution',...
'Style','text',...
'TooltipString','Set the document resolution','Tag','text7');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@imageResize_uicallback,h1,'popup_docResolution_Callback'},...
'Position',[158 68 83 22],...
'String',{  'pixels/inch'; 'pixels/cm' },...
'Style','popupmenu',...
'TooltipString','Set the document resolution',...
'Value',1,'Tag','popup_docResolution');

uicontrol('Parent',h1,...
'Callback',{@imageResize_uicallback,h1,'checkbox_constProportions_Callback'},...
'FontSize',9,'Position',[7 36 151 16],...
'String','Constrain Proportions',...
'Style','checkbox',...
'TooltipString','Constrain aspect ratio',...
'Value',1,'Tag','checkbox_constProportions');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@imageResize_uicallback,h1,'popup_resampMethod_Callback'},...
'Position',[127 10 121 22],...
'String',{  'Nearest Neighbor'; 'Bilinear'; 'Bicubic'; 'Area Relation' },...
'Style','popupmenu',...
'Value',1,'Tag','popup_resampMethod');

uicontrol('Parent',h1,...
'Callback',{@imageResize_uicallback,h1,'pushbutton_OK_Callback'},...
'FontSize',10,...
'Position',[257 229 66 23],...
'String','OK','Tag','pushbutton_OK');

uicontrol('Parent',h1,...
'Callback',{@imageResize_uicallback,h1,'pushbutton_cancel_Callback'},...
'FontSize',10,...
'Position',[257 189 66 23],...
'String','Cancel','Tag','pushbutton_cancel');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[7 11 116 18],...
'String','Resample Method:',...
'Style','text','Tag','text8');

function imageResize_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

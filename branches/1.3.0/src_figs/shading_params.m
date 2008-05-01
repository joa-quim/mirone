function varargout = shading_params(varargin)
% M-File changed by desGUIDE 

%	Copyright (c) 2004-2006 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

hObject = figure('Tag','figure1','Visible','off');
handles = guihandles(hObject);
guidata(hObject, handles);
shading_params_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
 
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% varargin   command line arguments to shading_params_export (see VARARGIN)

global home_dir
movegui(hObject,'center');               % Reposition the window on screen

% Case when this function was called directly
if isempty(home_dir),   f_path = [pwd filesep 'data' filesep];
else                    f_path = [home_dir filesep 'data' filesep];   end

handles.mercedes = 0;       % Flag for False Color option  (when set to 1)
handles.dirDerivative = 0;  % Flag for Directional Derivative option

if (length(varargin) >= 1)
    if (~strcmp(varargin{1},'dirDerivative'))
        load([f_path 'mirone_icons.mat'],'um_ico','dois_ico','tres_ico','quatro_ico','cinco_ico','seis_ico','sete_ico');
		h_toolbar = uitoolbar('parent',hObject,'Clipping', 'on', 'BusyAction','queue','HandleVisibility','on',...
           'Interruptible','on','Tag','FigureToolBar','Visible','on');
		handles.ui_grdgrad_A = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'grdgradient_A'},...
            'TooltipString','GMT grdgradient classic', 'CData',um_ico);
		handles.ui_grdgrad_E1 = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'grdgradient_E1'},...
            'TooltipString','GMT grdgradient Lambertian', 'CData',dois_ico);
		handles.ui_grdgrad_E2 = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'grdgradient_E2'},...
            'TooltipString','GMT grdgradient Peucker', 'CData',tres_ico);
		handles.ui_lambert = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'lambertian'},...
            'TooltipString','Lambertian with lighting', 'CData',quatro_ico);
		handles.ui_color = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'color'},...
            'TooltipString','Color (Manip-Raster)', 'CData',cinco_ico);
		handles.ui_gray = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'gray'},...
            'TooltipString','Gray (Manip-Raster)', 'CData',seis_ico);
		handles.ui_falseColor = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'mercedes'},...
            'TooltipString','False color', 'CData',sete_ico);
        
        % The following is for use in toggle_uis(...). It's easier there to deal with
        % numbers, but in other places its preferable to have names. So we have a
        % duplicate. Attention, the order in .ui_tools must reproduce its declarations
        handles.ui_tools = [handles.ui_grdgrad_A handles.ui_grdgrad_E1 handles.ui_grdgrad_E2 ...
            handles.ui_lambert handles.ui_color handles.ui_gray handles.ui_falseColor];
    else
        % With this option the better is realy not to show the rest of the window elements
        pos_f = get(hObject,'Position');    % Original Fig size
        pos_a = get(handles.axes1,'Position');    pos_a(1) = 0;
        set(handles.axes2,'Visible','off')
        pos_textAzim = get(handles.text_azim,'Position');       % Get and change text_azim position
        set(handles.text_azim,'Position',[2 pos_textAzim(2) pos_textAzim(3)-2 pos_textAzim(4)])
        pos_azim = get(handles.edit_azim,'Position');       % Get and change edit_azim position
        set(handles.edit_azim,'Position',pos_azim+[-14 0 0 0])
        pos_ok = get(handles.pushbutton_OK,'Position');     pos_ok(1) = 90; pos_ok(3) = 40;    pos_ok(4) = 20;
        set(handles.pushbutton_OK,'Position',pos_ok)
        set(hObject,'Position',[pos_f(1) pos_f(2) pos_a(1)+pos_a(3)+14 pos_f(4)],'Name','')   % New figure's size
        handles.dirDerivative = 1;
    end
else
    errordlg('Unknown Illumination option','Error')
end

% Import background image
astrolabio = imread([f_path 'astrolabio.jpg']);
image(astrolabio,'parent',handles.axes1);

pos = get(handles.axes1,'Position');
set(handles.axes1,'Visible','off')

% Draw everything that may be needed for all options. Later, depending on the
% option selected, only the allowed features will be made visible
x0 = pos(3)/2;      y0 = pos(4)/2;      radius = pos(3)/2;
h_line(1) = line('parent',handles.axes1,'XData',[x0 x0],'YData',[y0 0],'Color','r','Tag','red','LineWidth',3,'Userdata',radius);
if (handles.dirDerivative == 0)         % Otherwise there is no point in creating those
	x1 = x0 + radius * cos(30*pi/180);      y1 = y0 + radius * sin(30*pi/180);
	h_line(2) = line('parent',handles.axes1,'XData',[x0 x1],'YData',[y0 y1],'Color','g','Tag','green','LineWidth',3,'Visible','off');
	x1 = x0 + radius * cos(150*pi/180);     y1 = y0 + radius * sin(150*pi/180);
	h_line(3) = line('parent',handles.axes1,'XData',[x0 x1],'YData',[y0 y1],'Color','b','Tag','blue','LineWidth',3,'Visible','off');
	set(h_line,'Userdata',radius)        % save radius of circumscribed circle (image is square)
	% Now draw, on axes2, a quarter of circle and a line
	t = 0:0.02:pi/2;    x = [0 cos(t) 0];     y = [0 sin(t) 0];
	line('parent',handles.axes2,'XData',x,'YData',y,'HitTest','off','Color','k','LineWidth',1);
	h_line(4) = line('parent',handles.axes2,'XData',[0 cos(30*pi/180)],'YData',[0 sin(30*pi/180)],'Color','k','LineWidth',3,'Visible','off');
	set(h_line(4),'Tag','Elev','Userdata',1)        % save radius of circumscribed circle
end

switch varargin{1}
    case 'grdgradient_A'
        set(handles.ui_grdgrad_A,'State','on');
    case 'dirDerivative'
        % Do nothing here. Just to account for this option            
    otherwise
        errordlg('Unknown Illumination option','Error')
end

handles.h_line = h_line;
guidata(hObject, handles);
show_needed(hObject,[],varargin{1})
set(hObject,'WindowButtonDownFcn',{@ButtonDown,h_line,handles});

% Choose default command line output for shading_params_export
handles.output = hObject;
guidata(hObject, handles);

set(hObject,'Visible','on');
% UIWAIT makes shading_params_export wait for user response (see UIRESUME)
uiwait(handles.figure1);

handles = guidata(hObject);
out = shading_params_OutputFcn(hObject, [], handles);
varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = shading_params_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;
% The figure can be deleted now
delete(handles.figure1);

% -----------------------------------------------------------------------------------------
function show_needed(obj,eventdata,opt)
handles = guidata(obj);         % Get handles
h_all = handles.h_line;
handles.mercedes = 0;
if (strncmp(opt,'grdgradient',11))
    set(handles.edit_elev,'Enable','off');          set(handles.edit_azim,'Visible','on')
    set(handles.edit_ambient,'Enable','off');       set(handles.edit_diffuse,'Enable','off')
    set(handles.edit_specular,'Enable','off');      set(handles.edit_shine,'Enable','off')
    set(handles.edit_azimR,'Visible','off');        set(handles.edit_azimG,'Visible','off')
    set(handles.edit_azimB,'Visible','off');        set(handles.text_elev,'Enable','on');
    set(handles.text_ambient,'Enable','off');       set(handles.text_diffuse,'Enable','off');
    set(handles.text_reflection,'Enable','off');    set(handles.text_shine,'Enable','off');
    if (strcmp(opt(12:end),'_A'))
        set(handles.edit_azim,'Enable','on');       set(handles.text_azim,'Enable','on');
        set(h_all(1),'Visible','on');               set(h_all(2:4),'Visible','off')
        toggle_uis(handles,1);                      set(handles.figure1,'Name','GMT grdgradient')
    elseif (strcmp(opt(12:end),'_E1'))
        set(h_all([1 4]),'Visible','on');           set(h_all(2:3),'Visible','off')
        set(handles.edit_azim,'Enable','on');       set(handles.text_azim,'Enable','on');
        set(handles.edit_elev,'Enable','on');       toggle_uis(handles,2);
        set(handles.figure1,'Name','GMT grdgradient - Lambertian')
    else        % _E2
        set(handles.edit_azim,'Enable','off');      set(handles.text_azim,'Enable','off');
        set(handles.text_elev,'Enable','off');
        set(h_all(1:4),'Visible','off');            toggle_uis(handles,3);
        set(handles.figure1,'Name','GMT grdgradient - Peucker')
    end
elseif (strcmp(opt,'color') | strcmp(opt,'gray'))
    set(handles.edit_elev,'Enable','on');           set(handles.edit_azim,'Visible','on')
    set(handles.edit_ambient,'Enable','off');       set(handles.edit_diffuse,'Enable','off')
    set(handles.edit_specular,'Enable','off');      set(handles.edit_shine,'Enable','off')
    set(handles.edit_azimR,'Visible','off');        set(handles.edit_azimG,'Visible','off')
    set(handles.edit_azimB,'Visible','off');
    set(handles.text_elev,'Enable','on');
    set(handles.text_ambient,'Enable','off');       set(handles.text_diffuse,'Enable','off');
    set(handles.text_reflection,'Enable','off');    set(handles.text_shine,'Enable','off');
    set(handles.edit_azim,'Enable','on');           set(handles.text_azim,'Enable','on');
    set(h_all(1),'Visible','on');                   set(h_all(4),'Visible','on')
    set(h_all(2:3),'Visible','off')
    if (strcmp(opt,'color'))
        toggle_uis(handles,5);                      set(handles.figure1,'Name','Color')
    else
        toggle_uis(handles,6);                      set(handles.figure1,'Name','Gray')
    end
elseif (strcmp(opt,'lambertian'))
    set(handles.edit_elev,'Enable','on');           set(handles.edit_azim,'Visible','on')
    set(handles.edit_ambient,'Enable','on');        set(handles.edit_diffuse,'Enable','on')
    set(handles.edit_specular,'Enable','on');       set(handles.edit_shine,'Enable','on')
    set(handles.edit_azimR,'Visible','off');        set(handles.edit_azimG,'Visible','off')
    set(handles.edit_azimB,'Visible','off');
    set(handles.text_elev,'Enable','on');
    set(handles.text_ambient,'Enable','on');        set(handles.text_diffuse,'Enable','on');
    set(handles.text_reflection,'Enable','on');     set(handles.text_shine,'Enable','on');
    set(handles.edit_azim,'Enable','on');           set(handles.text_azim,'Enable','on');
    set(h_all(1),'Visible','on');                   set(h_all(4),'Visible','on')
    set(h_all(2:3),'Visible','off')
    toggle_uis(handles,4)
    set(handles.figure1,'Name','Lambertian lighting')
elseif (strcmp(opt,'mercedes'))
    set(handles.edit_elev,'Enable','on');           set(handles.edit_azim,'Visible','off')
    set(handles.edit_ambient,'Enable','off');       set(handles.edit_diffuse,'Enable','off')
    set(handles.edit_specular,'Enable','off');      set(handles.edit_shine,'Enable','off')
    set(handles.edit_azimR,'Visible','on');         set(handles.edit_azimG,'Visible','on')
    set(handles.edit_azimB,'Visible','on');
    set(handles.text_elev,'Enable','on');
    set(handles.text_ambient,'Enable','off');       set(handles.text_diffuse,'Enable','off');
    set(handles.text_reflection,'Enable','off');    set(handles.text_shine,'Enable','off');
    set(h_all(1:4),'Visible','on')
    handles.mercedes = 1;
    toggle_uis(handles,7)
    set(handles.figure1,'Name','False color')
elseif (strcmp(opt,'dirDerivative'))            % This for good because this function won't be called again
    %set(h_all(1),'Visible','on');                   %set(h_all(2:4),'Visible','off')
    set(handles.edit_azimR,'Visible','off');        set(handles.edit_azimG,'Visible','off')
    set(handles.edit_azimB,'Visible','off');        set(handles.text_elev,'Enable','off');
    set(handles.figure1,'Name','Azim')
end
guidata(obj,handles)

% -----------------------------------------------------------------------------------------
function toggle_uis(handles,ui)
% Do not let more the one uitoggletool be on the state of pushed
n = 1:length(handles.ui_tools);
n(n == ui) = [];        % Remove current ui index
set(handles.ui_tools(n),'State','off');

% -----------------------------------------------------------------------------------------
function ButtonDown(obj,eventdata,h_all,handles)
% It could be cleverer.
pt = get(gca, 'CurrentPoint');
x_lim = get(gca,'xlim');      y_lim = get(gca,'ylim');
% check if x,y is inside of axis
if ~((pt(1,1)>=x_lim(1)) & (pt(1,1)<=x_lim(2)) & (pt(1,2)>=y_lim(1)) & (pt(1,2)<=y_lim(2)))    % outside axis limits
    return
end
if any(h_all == gco)
    h = h_all(h_all == gco);    % When more than one line handle exists, find only the selected one
    set(gcf,'WindowButtonMotionFcn',{@ButtonMotion,h,handles},'WindowButtonUpFcn',{@ButtonUp,h_all,handles},...
        'Pointer', 'crosshair');
else
    return;
end

% -----------------------------------------------------------------------------------------
function ButtonMotion(obj,eventdata,h,handles)
selectionType = get(gcf, 'SelectionType');
pt = get(gca, 'CurrentPoint');
if strcmp(selectionType, 'normal')      % right-cick
    xx = get(h,'XData');    yy = get(h,'YData');
    theta = cart2pol(pt(1,1)-xx(1),pt(1,2)-yy(1));
    radius = get(h,'Userdata');
    x2 = xx(1) + radius * cos(theta);      y2 = yy(1) + radius * sin(theta);
    if strcmp(get(h,'Tag'),'Elev') & (theta >= 0 & theta <= pi/2)   % Elevation line
        set(h,'XData',[xx(1) x2],'YData',[yy(1) y2]);
        set(handles.edit_elev,'String',num2str(fix(theta *180/pi)) )
    elseif ~strcmp(get(h,'Tag'),'Elev')     % Azimuth line(s)
        set(h,'XData',[xx(1) x2],'YData',[yy(1) y2]);
    
        % NOTE to if I ever want to reuse this code. Normally ang_2pi should be = pi/2 - (pi*.....)
        % for the normal y origin at bottm left corner. However, due to the stupid habit of using y=0
        % at top left corner when dealing with images, to get an azimuth angle we have to do like following. 

        % truncate angles into [-pi pi] range
        ang_2pi = pi/2 + ( pi*((abs(theta)/pi) - 2*ceil(((abs(theta)/pi)-1)/2)) * sign(theta) );
        epsilon = -1e-7;        %  Allow points near zero to remain there
        indx = find(ang_2pi < epsilon);
        %  Shift the points in the [-pi 0] range to [pi 2pi] range
        if ~isempty(indx);  ang_2pi(indx) = ang_2pi(indx) + 2*pi;  end;
        if strcmp(get(h,'Tag'),'red')
            if (~handles.mercedes)
                set(handles.edit_azim,'String',num2str(fix(ang_2pi *180/pi)) )
            else
                set(handles.edit_azimR,'String',num2str(fix(ang_2pi *180/pi)) )
            end
        elseif strcmp(get(h,'Tag'),'green')
            set(handles.edit_azimG,'String',num2str(fix(ang_2pi *180/pi)) )
        elseif strcmp(get(h,'Tag'),'blue')
            set(handles.edit_azimB,'String',num2str(fix(ang_2pi *180/pi)) )
        end
    end
end

% -----------------------------------------------------------------------------------------
function ButtonUp(obj,eventdata,h,handles)
set(gcf,'WindowButtonMotionFcn','','WindowButtonDownFcn',{@ButtonDown,h,handles},'WindowButtonUpFcn','');
set(gcf,'Pointer', 'arrow')

% -----------------------------------------------------------------------------------------
function edit_azim_Callback(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------------
function edit_elev_Callback(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------------
function edit_ambient_Callback(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------------
function edit_diffuse_Callback(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------------
function edit_specular_Callback(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------------
function edit_shine_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function edit_azimR_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function edit_azimG_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function edit_azimB_Callback(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
if (handles.mercedes == 0)
    out.azim = str2num(get(handles.edit_azim,'String'));
else
    out.azim(1) = str2num(get(handles.edit_azimR,'String'));
    out.azim(2) = str2num(get(handles.edit_azimG,'String'));
    out.azim(3) = str2num(get(handles.edit_azimB,'String'));
end
out.elev = str2num(get(handles.edit_elev,'String'));
out.ambient = str2num(get(handles.edit_ambient,'String'));
out.diffuse = str2num(get(handles.edit_diffuse,'String'));
out.specular = str2num(get(handles.edit_specular,'String'));
out.shine = str2num(get(handles.edit_shine,'String'));

if (handles.dirDerivative == 0)
	% Find out which illumination model has been choosen. 
	% This is needed in Mirone to know what to do with the out vars
	if (strcmp(get(handles.ui_grdgrad_A,'State'),'on'))
        out.illum_model = 1;
	elseif (strcmp(get(handles.ui_grdgrad_E1,'State'),'on'))
        out.illum_model = 2;
	elseif (strcmp(get(handles.ui_grdgrad_E2,'State'),'on'))
        out.illum_model = 3;
	elseif (strcmp(get(handles.ui_lambert,'State'),'on'))
        out.illum_model = 4;
	elseif (strcmp(get(handles.ui_color,'State'),'on'))
        out.illum_model = 5;
	elseif (strcmp(get(handles.ui_gray,'State'),'on'))
        out.illum_model = 6;
	elseif (strcmp(get(handles.ui_falseColor,'State'),'on'))
        out.illum_model = 7;
	else
        errordlg('Uknown illumination model.','Error')
        out = [];
	end
end

handles.output = out;
guidata(hObject,handles);
uiresume(handles.figure1);

% --------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% Hint: delete(hObject) closes the figure
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    handles.output = [];        % User gave up, return nothing
    guidata(hObject, handles);
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    handles.output = [];        % User gave up, return nothing
    guidata(hObject, handles);
    delete(handles.figure1);
end

% --------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% Check for "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = [];    % User said no by hitting escape
    guidata(hObject, handles);
    uiresume(handles.figure1);
end


% --- Creates and returns a handle to the GUI figure. 
function shading_params_LayoutFcn(h1,handles);

set(h1,...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','shading_params',...
'NumberTitle','off',...
'Position',[520 400 320 150],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

h2 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@shading_params_uicallback,h1,'edit_azim_Callback'},...
'Position',[56 6 34 18],...
'String','0',...
'Style','edit',...
'TooltipString','Azimuth direction',...
'Tag','edit_azim');

h3 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@shading_params_uicallback,h1,'edit_elev_Callback'},...
'Position',[131 6 30 18],...
'String','30',...
'Style','edit',...
'TooltipString','Elevation light direction',...
'Tag','edit_elev');

h4 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@shading_params_uicallback,h1,'edit_ambient_Callback'},...
'Position',[186 102 30 18],...
'String','.55',...
'Style','edit',...
'Tag','edit_ambient');

h5 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@shading_params_uicallback,h1,'edit_diffuse_Callback'},...
'Position',[186 79 30 18],...
'String','.6',...
'Style','edit',...
'Tag','edit_diffuse');

h6 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@shading_params_uicallback,h1,'edit_specular_Callback'},...
'Position',[186 56 30 18],...
'String','.4',...
'Style','edit',...
'Tag','edit_specular');

h7 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@shading_params_uicallback,h1,'edit_shine_Callback'},...
'Position',[186 33 30 18],...
'String','10',...
'Style','edit',...
'Tag','edit_shine');

h8 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[219 102 61 16],...
'String','Ambient light',...
'Style','text',...
'Tag','text_ambient');

h9 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[219 78 83 17],...
'String','Diffuse reflection',...
'Style','text',...
'Tag','text_diffuse');

h10 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[219 57 91 15],...
'String','Specular reflection',...
'Style','text',...
'Tag','text_reflection');

h11 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[219 35 77 15],...
'String','Specular shine',...
'Style','text',...
'Tag','text_shine');

h12 = uicontrol('Parent',h1,...
'Callback',{@shading_params_uicallback,h1,'pushbutton_OK_Callback'},...
'Position',[229 6 66 23],...
'String','OK',...
'Tag','pushbutton_OK');

h13 = uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[13 7 42 16],...
'String','Azimuth',...
'Style','text',...
'Tag','text_azim');

h14 = uicontrol('Parent',h1,...
'FontSize',9,...
'Position',[123 104 50 16],...
'String','Elevation',...
'Style','text',...
'Tag','text_elev');

h15 = axes('Parent',h1,'Units','pixels','Position',[16 29 91 91],...
'Tag','axes1','Visible','off');

h20 = axes('Parent',h1,...
'Units','pixels',...
'Position',[126 49 51 51],...
'Tag','axes2',...
'Visible','off');

h25 = uicontrol('Parent',h1,...
'BackgroundColor',[1 0 0],...
'Callback',{@shading_params_uicallback,h1,'edit_azimR_Callback'},...
'Position',[10 6 34 18],...
'String','0',...
'Style','edit',...
'TooltipString','Red component azimuth',...
'Tag','edit_azimR');

h26 = uicontrol('Parent',h1,...
'BackgroundColor',[0 1 0],...
'Callback',{@shading_params_uicallback,h1,'edit_azimG_Callback'},...
'Position',[44 6 34 18],...
'String','120',...
'Style','edit',...
'TooltipString','Green component azimuth',...
'Tag','edit_azimG');

h27 = uicontrol('Parent',h1,...
'BackgroundColor',[0 0 1],...
'Callback',{@shading_params_uicallback,h1,'edit_azimB_Callback'},...
'Position',[79 6 34 18],...
'String','240',...
'Style','edit',...
'TooltipString','Blue component azimuth',...
'Tag','edit_azimB');

function shading_params_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

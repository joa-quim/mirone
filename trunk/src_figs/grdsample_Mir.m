function varargout = grdsample_Mir(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to grdsample_Mir (see VARARGIN)
%
%   The output is a structure with one or more of the following fields (or empty):
%   out.opt_R   -> only if grid limits are changed
%   out.opt_N   -> allways (even if user only changes x_inc/y_inc)
%   out.opt_Q   -> if bilinear was choosen
%   out.opt_L   -> if boundary conditions were selected

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
grdsample_Mir_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

movegui(hObject,'center')

if ~isempty(varargin)
    handles.h_calling_fig = varargin{1};
else
    % Tenho de prever este caso (mas como?)
end

handles.x_min = [];             handles.x_max = [];
handles.y_min = [];             handles.y_max = [];
handles.x_inc = [];             handles.y_inc = [];
handles.dms_xinc = 0;           handles.dms_yinc = 0;
handles.bd_cond = [];

%-----------
% Fill in the grid limits boxes with calling fig values and save some limiting value
head = getappdata(handles.h_calling_fig,'GMThead');
if isempty(head)
    % Fazer qualquer coisa (mas o que? aqui e dificil)
end
str = ddewhite(num2str(head(1),'%.10f'),'0');    set(handles.edit_x_min,'String',str)
str = ddewhite(num2str(head(2),'%.10f'),'0');    set(handles.edit_x_max,'String',str)
str = ddewhite(num2str(head(3),'%.10f'),'0');    set(handles.edit_y_min,'String',str)
str = ddewhite(num2str(head(4),'%.10f'),'0');    set(handles.edit_y_max,'String',str)
handles.x_min = head(1);            handles.x_max = head(2);
handles.y_min = head(3);            handles.y_max = head(4);
handles.x_min_or = head(1);         handles.x_max_or = head(2);
handles.y_min_or = head(3);         handles.y_max_or = head(4);
handles.one_or_zero = head(7);
[m,n] = size(getappdata(handles.h_calling_fig,'dem_z'));
handles.nr_or = m;                  handles.nc_or = n;;

% Fill in the x,y_inc and nrow,ncol boxes
set(handles.edit_Nrows,'String',num2str(m))
set(handles.edit_Ncols,'String',num2str(n))
% Compute default xinc, yinc based on map limits
yinc = (head(4) - head(3)) / (m-1);   xinc = (head(2) - head(1)) / (n-1);
set(handles.edit_y_inc,'String',num2str(yinc,10))
set(handles.edit_x_inc,'String',num2str(xinc,10))
handles.x_inc = xinc;    handles.y_inc = yinc;
%----------------

% Give a Pro look (3D) to the frame boxes 
bgcolor = get(0,'DefaultUicontrolBackgroundColor');
framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
set(hObject,'Units','pixels')    % Pixels are easier to reason with
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
h_t = findobj(hObject,'Style','Text');
for i=1:length(h_t)
    usr_d = get(h_t(i),'UserData');
    t_size = get(h_t(i),'Position');   t_str = get(h_t(i),'String');    fw = get(h_t(i),'FontWeight');
    bgc = get (h_t(i),'BackgroundColor');   fgc = get (h_t(i),'ForegroundColor');
    uicontrol('Parent',hObject, 'Style','text', 'Position',t_size,'String',t_str, ...
        'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw,'UserData',usr_d);
end
delete(h_t)

% Choose default command line output for grdsample_Mir_export
handles.output = hObject;
guidata(hObject, handles);

set(hObject,'Visible','on');
% UIWAIT makes grdsample_Mir_export wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --------------------------------------------------------------------
handles = guidata(hObject);
out = grdsample_Mir_OutputFcn(hObject, [], handles);
varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = grdsample_Mir_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;
% The figure can be deleted now
delete(handles.figure1);

% --------------------------------------------------------------------
function edit_x_min_Callback(hObject, eventdata, handles)
x_min_or = handles.x_min_or;
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)            % when dd:mm or dd:mm:ss was given
    x_min = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   x_min = x_min + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   x_min = x_min - abs(str2double(val{i})) / (60^(i-1));   end
    end
    if (x_min < x_min_or);  set(hObject,'String',num2str(x_min_or,6));    return;     end; 
    handles.x_min = x_min;
    if ~isempty(handles.x_max) & x_min >= handles.x_max
        errordlg('West Longitude >= East Longitude ','Error in Longitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    end
    nc = get(handles.edit_Ncols,'String');
    if ~isempty(handles.x_max) & ~isempty(nc)       % x_max and ncols boxes are filled
        % Compute Ncols, but first must recompute x_inc
        x_inc = ivan_the_terrible((handles.x_max - x_min),round(abs(str2double(nc))),1);
        xx = floor((handles.x_max - str2double(xx)) / (str2double(get(handles.edit_x_inc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_x_inc,'String',num2str(x_inc,8))
        guidata(hObject, handles);
    elseif ~isempty(handles.x_max)      % x_max box is filled but ncol is not, so put to the default (100)
        x_inc = ivan_the_terrible((handles.x_max - x_min),100,1);
        set(handles.edit_x_inc,'String',num2str(x_inc,8))
        set(handles.edit_Ncols,'String','100')
        guidata(hObject, handles);
    end
else                % box is empty, so clear also x_inc and ncols
    set(handles.edit_x_inc,'String','');     set(handles.edit_Ncols,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% --------------------------------------------------------------------
function edit_x_max_Callback(hObject, eventdata, handles)
x_max_or = handles.x_max_or;
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    x_max = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   x_max = x_max + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   x_max = x_max - abs(str2double(val{i})) / (60^(i-1));   end
    end
    if (x_max > x_max_or);  set(hObject,'String',num2str(x_max_or,6));    return;     end; 
    handles.x_max = x_max;
    if ~isempty(handles.x_min) & x_max <= handles.x_min 
        errordlg('East Longitude <= West Longitude','Error in Longitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    end
    nc = get(handles.edit_Ncols,'String');
    if ~isempty(handles.x_min) & ~isempty(nc)       % x_max and ncols boxes are filled
        % Compute Ncols, but first must recompute x_inc
        x_inc = ivan_the_terrible((x_max - handles.x_min),round(abs(str2double(nc))),1);
        xx = floor((handles.x_min - str2double(xx)) / (str2double(get(handles.edit_x_inc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_x_inc,'String',num2str(x_inc,8))
        guidata(hObject, handles);    
    elseif ~isempty(handles.x_min)      % x_min box is filled but ncol is not, so put to the default (100)
        x_inc = ivan_the_terrible((x_max - handles.x_min),100,1);
        set(handles.edit_x_inc,'String',num2str(x_inc,8))
        set(handles.edit_Ncols,'String','100')
        guidata(hObject, handles);
    end
else                % box is empty, so clear also x_inc and ncols
    set(handles.edit_x_inc,'String','');     set(handles.edit_Ncols,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% --------------------------------------------------------------------
function edit_y_min_Callback(hObject, eventdata, handles)
% Read value either in decimal or in the dd:mm or dd_mm:ss formats and do some tests
y_min_or = handles.y_min_or;
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    y_min = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   y_min = y_min + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   y_min = y_min - abs(str2double(val{i})) / (60^(i-1));   end
    end
    if (y_min < y_min_or);  set(hObject,'String',num2str(y_min_or,6));    return;     end; 
    handles.y_min = y_min;
    if ~isempty(handles.y_max) & y_min >= handles.y_max
        errordlg('South Latitude >= North Latitude','Error in Latitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    end
    nr = get(handles.edit_Nrows,'String');
    if ~isempty(handles.y_max) & ~isempty(nr)       % y_max and nrows boxes are filled
        % Compute Nrowss, but first must recompute y_inc
        y_inc = ivan_the_terrible((handles.y_max - y_min),round(abs(str2double(nr))),1);
        xx = floor((handles.y_max - str2double(xx)) / (str2double(get(handles.edit_y_inc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_y_inc,'String',num2str(y_inc,8))
        guidata(hObject, handles);
    elseif ~isempty(handles.y_max)      % y_max box is filled but nrows is not, so put to the default (100)
        y_inc = ivan_the_terrible((handles.y_max - y_min),100,1);
        set(handles.edit_y_inc,'String',num2str(y_inc,8))
        set(handles.edit_Nrows,'String','100')
        guidata(hObject, handles);
    end
else                % box is empty, so clear also y_inc and nrows
    set(handles.edit_y_inc,'String','');     set(handles.edit_Nrows,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% --------------------------------------------------------------------
function edit_y_max_Callback(hObject, eventdata, handles)
y_max_or = handles.y_max_or;
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    y_max = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   y_max = y_max + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   y_max = y_max - abs(str2double(val{i})) / (60^(i-1));   end
    end
    if (y_max > y_max_or);  set(hObject,'String',num2str(y_max_or,6));    return;     end; 
    handles.y_max = y_max;
    if ~isempty(handles.y_min) & y_max <= handles.y_min 
        errordlg('North Latitude <= South Latitude','Error in Latitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    end
    nr = get(handles.edit_Nrows,'String');
    if ~isempty(handles.y_min) & ~isempty(nr)       % y_min and nrows boxes are filled
        % Compute Nrows, but first must recompute y_inc
        y_inc = ivan_the_terrible((y_max - handles.y_min),round(abs(str2double(nr))),1);
        xx = floor((handles.y_min - str2double(xx)) / (str2double(get(handles.edit_y_inc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_y_inc,'String',num2str(y_inc,8))
        guidata(hObject, handles);
    elseif ~isempty(handles.y_min)      % y_min box is filled but nrows is not, so put to the default (100)
        y_inc = ivan_the_terrible((y_max - handles.y_min),100,1);
        set(handles.edit_y_inc,'String',num2str(y_inc,8))
        set(handles.edit_Nrows,'String','100')
        guidata(hObject, handles);
    end
else                % This box is empty, so clear also y_inc and nrows
    set(handles.edit_y_inc,'String','');     set(handles.edit_Nrows,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% --------------------------------------------------------------------
function edit_x_inc_Callback(hObject, eventdata, handles)
dms = 0;
xx = get(hObject,'String');     val = test_dms(xx);
if isempty(val)
    set(hObject, 'String', '');    return
end
% If it survived then ...
if length(val) > 1    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
x_inc = 0;
for i = 1:length(val)   x_inc = x_inc + str2double(val{i}) / (60^(i-1));    end
if ~isempty(handles.x_min) & ~isempty(handles.x_max)
    % Make whatever x_inc given compatible with GMT_grd_RI_verify
    x_inc = ivan_the_terrible((handles.x_max - handles.x_min), x_inc,2);
    if ~dms         % case of decimal unities
        set(hObject,'String',num2str(x_inc,8))
        ncol = floor((handles.x_max - handles.x_min) / x_inc + 0.5) + handles.one_or_zero;
    else            % inc was in dd:mm or dd:mm:ss format
        ncol = floor((handles.x_max - handles.x_min) / x_inc + 0.5) + handles.one_or_zero;
        ddmm = dec2deg(x_inc);
        set(hObject,'String',ddmm)
    end
    set(handles.edit_Ncols,'String',num2str(ncol))
end
handles.dms_xinc = dms;
guidata(hObject, handles);
if isempty(get(handles.edit_y_inc,'String'))     set(handles.edit_y_inc,'String',xx);    end

% --------------------------------------------------------------------
function edit_Ncols_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if ~isempty(get(handles.edit_x_min,'String')) & ~isempty(get(handles.edit_x_max,'String')) & ...
        ~isempty(get(handles.edit_x_inc,'String')) & ~isempty(xx)
    x_inc = ivan_the_terrible((handles.x_max - handles.x_min),round(abs(str2double(xx))),1);
    if handles.dms_xinc        % x_inc was given in dd:mm:ss format
        ddmm = dec2deg(x_inc);
        set(handles.edit_x_inc,'String',ddmm)
    else                    % x_inc was given in decimal format
        set(handles.edit_x_inc,'String',num2str(x_inc,8));
    end
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function edit_y_inc_Callback(hObject, eventdata, handles)
dms = 0;
xx = get(hObject,'String');     val = test_dms(xx);
if isempty(val)
    set(hObject, 'String', '');    return
end
% If it survived then ...
if length(val) > 1    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
y_inc = 0;
for i = 1:length(val)   y_inc = y_inc + str2double(val{i}) / (60^(i-1));    end
if ~isempty(handles.y_min) & ~isempty(handles.y_max)
    % Make whatever y_inc given compatible with GMT_grd_RI_verify
    y_inc = ivan_the_terrible((handles.y_max - handles.y_min), y_inc,2);
    if ~dms         % case of decimal unities
        set(hObject,'String',num2str(y_inc,8))
        nrow = floor((handles.y_max - handles.y_min) / y_inc + 0.5) + handles.one_or_zero;
    else            % inc was in dd:mm or dd:mm:ss format
        nrow = floor((handles.y_max - handles.y_min) / y_inc + 0.5) + handles.one_or_zero;
        ddmm = dec2deg(y_inc);
        set(hObject,'String',ddmm)
    end
    set(handles.edit_Nrows,'String',num2str(nrow))
end
handles.dms_yinc = dms;
guidata(hObject, handles);

% --------------------------------------------------------------------
function edit_Nrows_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if ~isempty(get(handles.edit_y_min,'String')) & ~isempty(get(handles.edit_y_max,'String')) & ...
        ~isempty(get(handles.edit_y_inc,'String')) & ~isempty(xx)
    y_inc = ivan_the_terrible((handles.y_max - handles.y_min),round(abs(str2double(xx))),1);
    if handles.dms_yinc        % y_inc was given in dd:mm:ss format
        ddmm = dec2deg(y_inc);
        set(handles.edit_y_inc,'String',ddmm)
    else                    % y_inc was given in decimal format
        set(handles.edit_y_inc,'String',num2str(y_inc,8));
    end
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function pushbutton_Help_R_F_T_Callback(hObject, eventdata, handles)
message = {'Min and Max, of "X Direction" and "Y Direction" specify the Region of'
    'interest. To specify boundaries in degrees and minutes [and seconds],'
    'use the dd:mm[:ss.xx] format.'
    '"Spacing" sets the grid size for grid output. You may choose different'
    'spacings for X and Y. Also here you can use the dd:mm[:ss.xx] format.'
    'In "#of lines" it is offered the easyeast way of controling the grid'
    'dimensions (lines & columns).'};
helpdlg(message,'Help on Grid Line Geometry');

% --------------------------------------------------------------------
function checkbox_Option_Q_Callback(hObject, eventdata, handles)
% No need for coding

% --------------------------------------------------------------------
function popup_BoundaryCondition_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case ' ',        handles.bd_cond = [];
    case '',         handles.bd_cond = [];
    case 'x',        handles.bd_cond = ['-Lx'];
    case 'y',        handles.bd_cond = ['-Ly'];
    case 'xy',       handles.bd_cond = ['-Lxy']; 
    case 'g',        handles.bd_cond = ['-Lg'];
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function pushbutton_Help_L_Callback(hObject, eventdata, handles)
message = {'Boundary condition flag may be "x" or "y" or "xy" indicating data is periodic'
           'in range of x or y or both set by the grids limits in the above boxes,'
           'or flag may be "g" indicating geographical conditions (x and y may be'
           'lon and lat). [Default is no boundary conditions].'};
helpdlg(message,'Help -L option');

% --------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
handles.output = '';        % User gave up, return nothing
guidata(hObject, handles);  uiresume(handles.figure1);

% --------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
out = [];   n_set = 0;
x_min = get(handles.edit_x_min,'String');   x_max = get(handles.edit_x_max,'String');
y_min = get(handles.edit_y_min,'String');   y_max = get(handles.edit_y_max,'String');
if isempty(x_min) | isempty(x_max) | isempty(y_min) | isempty(y_max)
    errordlg('One or more grid limits are empty. Open your yes.','Error');    return
end

nx = str2double(get(handles.edit_Ncols,'String'));
ny = str2double(get(handles.edit_Nrows,'String'));
if (isnan(nx) | isnan(ny))      % I think this was already tested, but ...
    errordlg('One (or two) of the grid dimensions are not valid. Do your best.','Error');   return
end

if (nx ~= handles.nc_or | ny ~= handles.nr_or)
    out.opt_N = ['-N' get(handles.edit_Ncols,'String') '/' get(handles.edit_Nrows,'String')];
    n_set = 1;
end

if (get(handles.checkbox_Option_Q,'Value'))
    out.opt_Q = '-Q';
end

if (~isempty(handles.bd_cond))
    out.opt_L = handles.bd_cond;
end

% See if grid limits were changed
if ( (abs(handles.x_min-handles.x_min_or) > 1e-4) | (abs(handles.x_max-handles.x_max_or) > 1e-4) | ...
        (abs(handles.y_min-handles.y_min_or) > 1e-4) | (abs(handles.y_max-handles.y_max_or) > 1e-4))
    out.opt_R = ['-R' num2str(handles.x_min,6) '/' num2str(handles.x_max,6) '/' ...
        num2str(handles.y_min,6) '/' num2str(handles.y_max,6)];
    out.x_min = handles.x_min;    out.x_max = handles.x_max;
    out.y_min = handles.y_min;    out.y_max = handles.y_max;
end

if (~n_set)
    errordlg('You haven''t select anything usefull to do.','Chico Clever');   return
end

handles.output = out;
guidata(hObject, handles);    uiresume(handles.figure1);

% --------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    handles.output = [];        % User gave up, return nothing
    guidata(hObject, handles);    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    handles.output = [];        % User gave up, return nothing
    guidata(hObject, handles);    delete(handles.figure1);
end

% --------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% Check for "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = '';    % User said no by hitting escape
    guidata(hObject, handles);    uiresume(handles.figure1);
end


% --- Creates and returns a handle to the GUI figure. 
function grdsample_Mir_LayoutFcn(h1,handles);

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','grdsample',...
'NumberTitle','off',...
'Position',[265.768111202607 398.94409531652 441 158],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

h2 = uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 75 421 75],...
'String',{  '' },...
'Style','frame',...
'Tag','frame2');

h3 = uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[30 141 121 15],...
'String','Griding Line Geometry',...
'Style','text',...
'Tag','text8');

h4 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdsample_Mir_uicallback,h1,'edit_x_min_Callback'},...
'HorizontalAlignment','left',...
'Position',[77 111 80 21],...
'Style','edit',...
'Tag','edit_x_min');

h5 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdsample_Mir_uicallback,h1,'edit_x_max_Callback'},...
'HorizontalAlignment','left',...
'Position',[163 111 80 21],...
'Style','edit',...
'Tag','edit_x_max');

h6 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdsample_Mir_uicallback,h1,'edit_y_min_Callback'},...
'HorizontalAlignment','left',...
'Position',[77 85 80 21],...
'Style','edit',...
'Tag','edit_y_min');

h7 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdsample_Mir_uicallback,h1,'edit_y_max_Callback'},...
'HorizontalAlignment','left',...
'Position',[163 85 80 21],...
'Style','edit',...
'Tag','edit_y_max');

h8 = uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[22 115 55 15],...
'String','X Direction',...
'Style','text',...
'Tag','text2');

h9 = uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[21 89 55 15],...
'String','Y Direction',...
'Style','text',...
'Tag','text3');

h10 = uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[183 132 41 13],...
'String','Max',...
'Style','text',...
'Tag','text4');

h11 = uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[99 132 41 13],...
'String','Min',...
'Style','text',...
'Tag','text5');

h12 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdsample_Mir_uicallback,h1,'edit_x_inc_Callback'},...
'HorizontalAlignment','left',...
'Position',[248 111 71 21],...
'Style','edit',...
'TooltipString','DX grid spacing',...
'Tag','edit_x_inc');

h13 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdsample_Mir_uicallback,h1,'edit_y_inc_Callback'},...
'HorizontalAlignment','left',...
'Position',[248 85 71 21],...
'Style','edit',...
'TooltipString','DY grid spacing',...
'Tag','edit_y_inc');

h14 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdsample_Mir_uicallback,h1,'edit_Ncols_Callback'},...
'HorizontalAlignment','left',...
'Position',[324 111 65 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Ncols');

h15 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdsample_Mir_uicallback,h1,'edit_Nrows_Callback'},...
'HorizontalAlignment','left',...
'Position',[324 85 65 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Nrows');

h16 = uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[265 133 41 13],...
'String','Spacing',...
'Style','text',...
'Tag','text6');

h17 = uicontrol('Parent',h1,...
'CData',[],...
'Enable','inactive',...
'Position',[332 133 51 13],...
'String','# of lines',...
'Style','text',...
'Tag','text7',...
'UserData',[]);

h18 = uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Callback',{@grdsample_Mir_uicallback,h1,'pushbutton_Help_R_F_T_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[400 84 21 48],...
'String','?',...
'Tag','pushbutton_Help_R_F_T');

h19 = uicontrol('Parent',h1,...
'Callback',{@grdsample_Mir_uicallback,h1,'checkbox_Option_Q_Callback'},...
'Position',[221 42 114 15],...
'String','Bilinear interpolation',...
'Style','checkbox',...
'TooltipString','Use bilinear rather than bicubic interpolation',...
'Tag','checkbox_Option_Q');

h20 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdsample_Mir_uicallback,h1,'popup_BoundaryCondition_Callback'},...
'HorizontalAlignment','right',...
'Position',[114 39 47 21],...
'String',{  ' '; 'x'; 'y'; 'xy'; 'g' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_BoundaryCondition');

h21 = uicontrol('Parent',h1,...
'Callback',{@grdsample_Mir_uicallback,h1,'pushbutton_Help_L_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[172 38 21 23],...
'String','?',...
'Tag','pushbutton_Help_L');

h22 = uicontrol('Parent',h1,...
'Callback',{@grdsample_Mir_uicallback,h1,'pushbutton_cancel_Callback'},...
'Position',[295 8 66 23],...
'String','Cancel',...
'Tag','pushbutton_cancel');

h23 = uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[17 41 95 16],...
'String','Boundary condition',...
'Style','text',...
'Tag','text10');

h24 = uicontrol('Parent',h1,...
'Callback',{@grdsample_Mir_uicallback,h1,'pushbutton_OK_Callback'},...
'Position',[365 8 66 23],...
'String','OK',...
'Tag','pushbutton_OK');

function grdsample_Mir_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

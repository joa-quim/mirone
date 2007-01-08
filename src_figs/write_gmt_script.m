function varargout = write_gmt_script(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to write_gmt_script (see VARARGIN)

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

%# mex

hObject = figure('Tag','figure1','Visible','off');
handles = guihandles(hObject);
guidata(hObject, handles);
write_gmt_script_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
movegui(hObject,'center');

sizes_cm = {'A0 (83.96 118.82 cm)'; 'A1 (59.41 83.96 cm)'; 'A2 (41.98 59.41 cm)'; 'A3 (29.70 41.98 cm)'
    'A4 (20.99 29.70 cm)'; 'A5 (14.85 20.99 cm)'; 'A6 (10.48 14.85 cm)'; 'A7 (7.41 10.48 cm)'
    'A8 (5.22 7.41 cm)'; 'A9 (3.70 5.22 cm)'; 'A10 (2.61 3.70 cm)'; 'B0 (100.05 141.39 cm)'
    'B1 (70.70 100.05 cm)'; 'B2 (50.02 70.70 cm)'; 'B3 (35.35 50.02 cm)'; 'B4 (25.01 35.35 cm)'
    'B5 (17.67 25.01 cm)'; 'archA (22.86 30.48 cm)'; 'archB (30.48 45.72 cm)'; 'archC (45.72 60.96 cm)'
    'archD (60.96 91.44 cm)'; 'archE (91.44 121.92 cm)'; 'flsa (21.59 33.02 cm)'; 'halfletter (13.97 21.59 cm)'
    'note (19.05 25.40 cm)'; 'letter (21.59 27.94 cm)'; 'legal (21.59 35.56 cm)'; '11x17 (27.94 43.18 cm)'
    'ledger (43.18 27.94 cm)'};

sizes_pt = {'A0  (2380 3368 pt)'; 'A1  (1684 2380 pt)'; 'A2  (1190 1684 pt)'; 'A3  (842 1190 pt)'
    'A4  (595 842 pt)'; 'A5  (421 595 pt)'; 'A6  (297 421 pt)'; 'A7  (210 297 pt)'; 'A8  (148 210 pt)'
    'A9  (105 148 pt)'; 'A10 (74 105 pt)'; 'B0  (2836 4008 pt)'; 'B1  (2004 2836 pt)'; 'B2  (1418 2004 pt)'
    'B3  (1002 1418 pt)'; 'B4  (709 1002 pt)'; 'B5  (501 709 pt)'; 'archA (648 864 pt)'; 'archB (864 1296 pt)'
    'archC (1296 1728 pt)'; 'archD (1728 2592 pt)'; 'archE (2592 3456 pt)'; 'flsa  (612 936 pt)'
    'halfletter (396 612 pt)'; 'note   (540 720 pt)'; 'letter (612 792 pt)'; 'legal  (612 1008 pt)'
    '11x17  (792 1224 pt)'; 'ledger (1224 792 pt)'};

sizes_in = {'A0 (33.06 46.78 cm)'; 'A1 (23.39 33.06 in)'; 'A2 (16.53 23.39 in)'; 'A3 (11.69 16.53 in)'
    'A4 (8.26 11.69 in)'; 'A5 (5.85 8.26 in)'; 'A6 (4.13 5.85 in)'; 'A7 (2.92 4.13 in)'
    'A8 (2.06 2.92 in)'; 'A9 (1.46 2.06 in)'; 'A10 (1.03 1.46 in)'; 'B0 (39.39 55.67 in)'
    'B1 (27.83 39.39 in)'; 'B2 (19.69 27.83 in)'; 'B3 (13.92 19.69 in)'; 'B4 (9.85 13.92 in)'
    'B5 (6.96 9.85 in)'; 'archA (9.0 12.0 in)'; 'archB (12.0 18.0 in)'; 'archC (18.0 24.0 in)'
    'archD (24.0 36.0 in)'; 'archE (36.0 48.0 in)'; 'flsa (8.5 13.0 in)'; 'halfletter (5.5 8.5 in)'
    'note (7.5 10.0 in)'; 'letter (8.5 11.0 in)'; 'legal (8.5 14.0 in)'; '11x17 (11.0 17.0 in)'
    'ledger (17.0 11.0 in)'};

paper_cm = [83.96 118.82; 59.41 83.96; 41.98 59.41; 29.70 41.98; 20.99 29.70; 14.85 20.99; ...
        10.48 14.85; 7.41 10.48; 5.22 7.41; 3.70 5.22; 2.61 3.70; 100.05 141.40; 70.70 100.05; ...
        50.02 70.70; 35.35 50.02; 25.01 35.35; 17.67 25.01; 22.86 30.48; 30.48 45.72; ...
        45.72 60.96; 60.96 91.44; 91.44 121.92; 21.59 33.02; 13.97 21.59; 19.05 25.40; ...
        21.59 27.94; 21.59 35.56; 27.94 43.18; 43.18 27.94];

paper_pt = [2380 3368; 1684 2380; 1190 1684; 842 1190; 595 842; 421 595; 297 421; 210 297; 148 210; ...
    105 148; 74 105; 2836 4008; 2004 2836; 1418 2004; 1002 1418; 709 1002; 501 709; 648 864; 864 1296; ...
    1296 1728; 1728 2592; 2592 3456; 612 936; 396 612; 540 720; 612 792; 612 1008; 792 1224; 1224 792];

paper_in = [33.06 46.78; 23.39 33.06; 16.53 23.39; 11.69 16.53; 8.26 11.69; 5.85 8.26; 4.13 5.85; ...
    2.92 4.13; 2.06 2.92; 1.46 2.06; 1.03 1.46; 39.39 55.67; 27.83 39.39; 19.69 27.83; 13.92 19.69; ...
    9.85 13.92; 6.96 9.85; 9.0 12.0; 12.0 18.0; 18.0 24.0; 24.0 36.0; 36.0 48.0; 8.5 13.0; 5.5 8.5; ...
    7.5 10.0; 8.5 11.0; 8.5 14.0; 11.0 17.0; 17.0 11.0];

set(handles.popup_PaperSize,'String',sizes_cm,'Value',5)

handles.sizes_cm = sizes_cm;    handles.paper_cm = paper_cm;
handles.sizes_pt = sizes_pt;    handles.paper_pt = paper_pt;
handles.sizes_in = sizes_in;    handles.paper_in = paper_in;
handles.opt_L = [];             handles.opt_U = [];
handles.opt_psc = [];           % To eventualy hold several of the pscoast options
handles.scale_set = 0;          % To signal that user changed scale

mirone_handles = varargin{1};
if (mirone_handles.no_file)     % Stupid call with nothing loaded on the Mirone window
    delete(hObject);    return
end
handles.script_type = varargin{2};
if (strcmp(handles.script_type,'bat'))
    set(hObject,'Name','Write GMT batch')
    set(handles.pushbutton_OK,'String','Write batch')
else
    set(hObject,'Name','Write GMT script')
end

% Add this figure handle to the carraças list
plugedWin = getappdata(mirone_handles.figure1,'dependentFigs');
plugedWin = [plugedWin hObject];
setappdata(mirone_handles.figure1,'dependentFigs',plugedWin);

if (mirone_handles.image_type == 1 || mirone_handles.image_type == 3 || mirone_handles.image_type == 4)
    head = mirone_handles.head;
    handles.x_min = head(1);    handles.x_max = head(2);
    handles.y_min = head(3);    handles.y_max = head(4);
    nx = round((head(2) - head(1)) / head(8));          % May not be exactly correct but is good enough here
    ny = round((head(4) - head(3)) / head(9));
elseif (mirone_handles.image_type == 2 | mirone_handles.image_type == 20)     % "trivial" images
    [ny,nx,nz] = size(get(mirone_handles.hImg,'CData'));
    zz1 = get(mirone_handles.axes1,'XLim');    zz2 = get(mirone_handles.axes1,'YLim');
    handles.x_min = zz1(1);     handles.x_max = zz1(2);
    handles.y_min = zz2(1);     handles.y_max = zz2(2);
    clear nz zz1 zz2;
end    
width  = 15;                    % Default starting width in cm
fac_15 = nx / width;
height = round(ny) / fac_15;
if (height > 27)                % That is, if height + Y0 nearly outside the A4 page
    while (height > 27)         % Make it approximately fit the page height
        height = round(height * 0.1);
        width  = round(width * 0.1);
    end
end
handles.opt_R = ['-R' num2str(handles.x_min,'%.9g') '/' num2str(handles.x_max,'%.9g') ...
        '/' num2str(handles.y_min,'%.9g') '/' num2str(handles.y_max,'%.9g')];

handles.mirone_handles = mirone_handles;
handles.width_or = width;  handles.height_or = height;
handles.scale = width;
handles.which_unit = 'cm';
handles.d_path = mirone_handles.path_data;

% Compute image aspect ratio and set axes 'PlotBoxAspectRatio' to it
handles.paper = [paper_cm(5,1) paper_cm(5,2)];         % Set to A4 (x,y)
handles.paper_aspect = handles.paper(2)/handles.paper(1);
% set(handles.axes1,'XLim',[0 handles.paper(1)],'YLim',[0 handles.paper(2)], ...
%         'PlotBoxAspectRatio',[1 handles.paper_aspect 1]);
set(handles.axes1,'XLim',[0 handles.paper(1)],'YLim',[0 handles.paper(2)], 'DataAspectRatio',[1 1 1]);
X0 = 2.5;               Y0 = 2.5;       % This is the GMT default's plot origin in cm
rect_x = [X0 X0 X0+width X0+width X0];
rect_y = [Y0 Y0+height Y0+height Y0 Y0];
handles.rect_x = rect_x;   handles.rect_y = rect_y;
handles.scale = num2str(width,'%.2g');
set(handles.edit_mapWidth,'String',num2str(width,'%.2f'))       % Fill the width editbox
set(handles.edit_mapHeight,'String',num2str(height,'%.2f'))     % Fill the height editbox

% ---------- Draw the grid->image rectangle
h = line('XData',rect_x,'YData',rect_y, 'Color','k','LineWidth',.5,'Tag','PlotRect');
% h = patch('XData',rect_x,'YData',rect_y,'FaceColor','w','EdgeColor','k','LineWidth',.5,'Tag','PlotRect');
% set(h,'ButtonDownFcn',{@move_rectangle,handles,h})

% ---------- Set line uicontexts
cmenuHand = uicontextmenu;
set(h, 'UIContextMenu', cmenuHand);
cb_moveRECT = {@move_rectangle,handles,h};
uimenu(cmenuHand, 'Label', 'Move rectangle', 'Callback', cb_moveRECT);
ui_edit_polygon(h)    % Set edition functions

handles.hand_rect = h;      % Save the rectangle hand
handles.hand_frame_proj = [];

coord_system_script = [];
directory_list = [];
load([handles.d_path 'mirone_pref.mat']);
% Check that the coord_system structure has no errors. If it has, load the default value.
% The result is used to update the handles structure.
[handles,handles.coord_system_script] = check_coord_system(handles,coord_system_script,'_script');

handles.all_ellipsoides = DefineEllipsoide;     % This is already in mirone_prefs

%----------- Recall previous settings stored in mirone_pref -------------------
handles.h_txt_info = findobj(hObject,'Tag','text_ProjDescription');
handles.txt_info_pos = get(handles.h_txt_info,'Position');
if (~mirone_handles.geog & iscell(handles.proj_info_txt_script))     % Need to do this because grid is not geog
    if (length(handles.proj_info_txt_script) == 5)      handles.proj_info_txt_script(3) = [];   end
    k = strfind(handles.proj_info_txt_script{1},'->');
    handles.proj_info_txt_script{1} = [handles.proj_info_txt_script{1}(1:k+1)];
    k = strfind(handles.proj_info_txt_script{2},'->');
    handles.proj_info_txt_script{2} = [handles.proj_info_txt_script{2}(1:k+2) '   Linear'];
    % Also remove the Ellipsoid info
    k = strfind(handles.proj_info_txt_script{end-1},'->');
    handles.proj_info_txt_script{end-1} = [handles.proj_info_txt_script{end-1}(1:k+2) '   NA'];
    % Change the -J... to -JX...
    k = strfind(handles.proj_info_txt_script{end},'-J');
    handles.proj_info_txt_script{end} = [handles.proj_info_txt_script{end}(1:k+2) 'X' width handles.which_unit(1)];
    handles.coord_system_script.projection = ['-JX' width handles.which_unit(1)];
end
set(handles.h_txt_info,'String',handles.proj_info_txt_script,'Position',handles.txt_info_pos)
handles.all_datums = datums;    % datums is a function in utils

% ---------- Split the scale from the projection string
tmp = handles.coord_system_script.projection;
if (~isempty(tmp))
	if (length(tmp) == 4 & strcmp(tmp(3),'m'))      % Simple Mercator has the form "-Jm1"
        tmp = tmp(1:end-1);
	elseif (length(tmp) == 3 & strcmp(upper(tmp(3)),'X'))  % Linear proj has the form "-JX"
        handles.opt_J_no_scale = [tmp(1:2) upper(tmp(3))]; % Save this
	else                                            % All other should terminate as "-J.../1"
        tmp = tmp(1:end-2);
        handles.opt_J_no_scale = [tmp(1:2) upper(tmp(3)) tmp(4:end)];           % Save this
	end
else
    handles.opt_J_no_scale = '-JX10';               % Use this default
end
%opt_J = [tmp(1:2) upper(tmp(3)) tmp(4:end) '/' handles.scale handles.which_unit(1)];
handles.curr_datum = handles.all_datums{handles.coord_system_script.datum_val,2};   % Save this
%handles.opt_J_no_scale = [tmp(1:2) upper(tmp(3)) tmp(4:end)];               % And this


% ----------- Use the directory list from mirone_pref
j = logical(zeros(1,length(directory_list)));           % vector for eventual cleaning non-existing dirs

if iscell(directory_list)                               % When exists a dir list in mirone_pref
    for i = 1:length(directory_list)
        if ~exist(directory_list{i},'dir'),   j(i) = 1;   end
    end
    directory_list(j) = [];                             % clean eventual non-existing directories
    if ~isempty(directory_list)                         % If there is one left
        set(handles.popup_directory_list,'String',directory_list)
        handles.last_directories = directory_list;
    else
        set(handles.popup_directory_list,'String',pwd)
        handles.last_directories = cellstr(pwd);
    end
else                                                    % mirone_pref had no dir list
    handles.last_directories = cellstr(pwd);
    set(handles.popup_directory_list,'String',handles.last_directories)
end

% --------- Set prefix name based on month and day numbers
prefix = clock;
prefix = ['mir' num2str(prefix(3)) '-' num2str(prefix(2))];
set(handles.edit_prefix,'String',prefix)

if (~mirone_handles.geog)   % Non geogs don't use scale bars
    set(handles.togglebutton_Option_L,'Visible','off')
    set(findobj(hObject,'Style','text','Tag','text_MapScale'), 'Visible','off');
end

% ---------- See if we have pscoast stuff
ALLlineHand = findobj(get(mirone_handles.axes1,'Child'),'Type','line');
handles.psc_res = [];   handles.psc_opt_W = [];     handles.psc_type_p = [];    handles.psc_type_r  = [];
if (~isempty(findobj(ALLlineHand,'Tag','CoastLineNetCDF')) | ~isempty(findobj(ALLlineHand,'Tag','Rivers')) ...
        | ~isempty(findobj(ALLlineHand,'Tag','PoliticalBoundaries')) )
	[handles.ALLlineHand, handles.psc_res, handles.psc_opt_W, handles.psc_type_p, handles.psc_type_r] = ...
        find_psc_stuff(ALLlineHand);
else
    set(handles.pushbutton_coastLines,'Visible', 'off')
    handles.ALLlineHand = ALLlineHand;
end;    clear ALLlineHand;

%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
bgcolor = get(0,'DefaultUicontrolBackgroundColor');
framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
set(0,'Units','pixels');    set(hObject,'Units','pixels')    % Pixels are easier to reason with
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
% ------------ END Pro look (3D) -------------------------------------------------------

% ------------ Apply inherited projection
guidata(hObject, handles);
pushbutton_uppdate_Callback(handles.pushbutton_uppdate, [], handles)
handles = guidata(hObject);     % Recover in "this handles" the changes donne in pushbutton_uppdate

% Choose default command line output for write_gmt_script_export
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes write_gmt_script_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = write_gmt_script_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = write_gmt_script_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;

% -----------------------------------------------------------------------------------
function popup_PaperSize_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
switch handles.which_unit
    case 'cm'
        lims = handles.paper_cm(val,1:2);
    case 'in'
        lims = handles.paper_in(val,1:2);
    case 'pt'
        lims = handles.paper_pt(val,1:2);
end
if (get(handles.radiobutton_P,'Value'))
    set(handles.axes1,'XLim',[0 lims(1)],'YLim',[0 lims(2)]);
else
    set(handles.axes1,'XLim',[0 lims(2)],'YLim',[0 lims(1)]);
end

% -----------------------------------------------------------------------------------
function radiobutton_P_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton_P
if (get(hObject,'Value'))
    img_size_x = get(handles.axes1,'YLim');     % Just swap x & y
    img_size_y = get(handles.axes1,'XLim');
    set(handles.axes1,'XLim',img_size_x,'YLim',img_size_y);
    %set(handles.axes1,'XLim',img_size_x,'YLim',img_size_y,'PlotBoxAspectRatio',[1/handles.paper_aspect 1 1]);
    set(handles.radiobutton_L,'Value',0)
else
    set(hObject,'Value',1)
end

% -----------------------------------------------------------------------------------
function radiobutton_L_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    img_size_x = get(handles.axes1,'YLim');     % Just swap x & y
    img_size_y = get(handles.axes1,'XLim');
    pbar = [1 1/handles.paper_aspect 1];
    set(handles.axes1,'XLim',img_size_x,'YLim',img_size_y);
    %set(handles.axes1,'XLim',img_size_x,'YLim',img_size_y,'PlotBoxAspectRatio',[1 1/handles.paper_aspect 1]);
    set(handles.radiobutton_P,'Value',0)
else
    set(hObject,'Value',1)
end

% -----------------------------------------------------------------------------------
function radiobutton_cm_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    val = get(handles.popup_PaperSize,'Value');
    set(handles.popup_PaperSize,'String',handles.sizes_cm,'Value',val)
    if (get(handles.radiobutton_P,'Value'))
        set(handles.axes1,'XLim',[0 handles.paper_cm(val,1)],'YLim',[0 handles.paper_cm(val,2)])
    else
        set(handles.axes1,'XLim',[0 handles.paper_cm(val,2)],'YLim',[0 handles.paper_cm(val,1)])
    end    
    conv_units(handles,'cm')
    set(handles.radiobutton_in,'Value',0)
    set(handles.radiobutton_pt,'Value',0)
    handles.which_unit = 'cm';      guidata(hObject,handles);
else
    set(hObject,'Value',1)
end

% -----------------------------------------------------------------------------------
function radiobutton_in_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    val = get(handles.popup_PaperSize,'Value');
    set(handles.popup_PaperSize,'String',handles.sizes_in,'Value',val)
    if (get(handles.radiobutton_P,'Value'))
        set(handles.axes1,'XLim',[0 handles.paper_in(val,1)],'YLim',[0 handles.paper_in(val,2)])
    else
        set(handles.axes1,'XLim',[0 handles.paper_in(val,2)],'YLim',[0 handles.paper_in(val,1)])
    end    
    conv_units(handles,'in')
    set(handles.radiobutton_cm,'Value',0)
    set(handles.radiobutton_pt,'Value',0)
    handles.which_unit = 'in';      guidata(hObject,handles);
else
    set(hObject,'Value',1)
end

% -----------------------------------------------------------------------------------
function radiobutton_pt_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    val = get(handles.popup_PaperSize,'Value');
    set(handles.popup_PaperSize,'String',handles.sizes_pt,'Value',val)
    if (get(handles.radiobutton_P,'Value'))
        set(handles.axes1,'XLim',[0 handles.paper_pt(val,1)],'YLim',[0 handles.paper_pt(val,2)])
    else
        set(handles.axes1,'XLim',[0 handles.paper_pt(val,2)],'YLim',[0 handles.paper_pt(val,1)])
    end    
    conv_units(handles,'pt')
    set(handles.radiobutton_cm,'Value',0)
    set(handles.radiobutton_in,'Value',0)
    handles.which_unit = 'pt';      guidata(hObject,handles);
else
    set(hObject,'Value',1)
end

% -----------------------------------------------------------------------------------
function conv_units(handles,dest)
xx = get(handles.hand_rect,'XData');        yy = get(handles.hand_rect,'YData');
xf = get(handles.hand_frame_proj,'XData');  yf = get(handles.hand_frame_proj,'YData');
if (strcmp(handles.which_unit,'cm') & strcmp(dest,'in'))
    xx = xx / 2.54;     yy = yy / 2.54;
    xf = xf / 2.54;     yf = yf / 2.54;
    set(handles.edit_mapWidth,'String',num2str( str2num(get(handles.edit_mapWidth,'String'))/2.54,'%.2f' ))
    set(handles.edit_mapHeight,'String',num2str( str2num(get(handles.edit_mapHeight,'String'))/2.54,'%.2f' ))
    set(handles.edit_X0,'String',num2str( str2num(get(handles.edit_X0,'String'))/2.54,'%.2f' ))
    set(handles.edit_Y0,'String',num2str( str2num(get(handles.edit_Y0,'String'))/2.54,'%.2f' ))
elseif (strcmp(handles.which_unit,'cm') & strcmp(dest,'pt'))
    xx = xx * 72/2.54;  yy = yy * 72/2.54;
    xf = xf * 72/2.54;  yf = yf * 72/2.54;
    set(handles.edit_mapWidth,'String',num2str( str2num(get(handles.edit_mapWidth,'String'))* 72/2.54,'%.2f' ))
    set(handles.edit_mapHeight,'String',num2str( str2num(get(handles.edit_mapHeight,'String'))* 72/2.54,'%.2f' ))
    set(handles.edit_X0,'String',num2str( str2num(get(handles.edit_X0,'String'))* 72/2.54,'%.2f' ))
    set(handles.edit_Y0,'String',num2str( str2num(get(handles.edit_Y0,'String'))* 72/2.54,'%.2f' ))
elseif (strcmp(handles.which_unit,'in') & strcmp(dest,'cm'))
    xx = xx * 2.54;     yy = yy * 2.54;
    xf = xf * 2.54;     yf = yf * 2.54;
    set(handles.edit_mapWidth,'String',num2str( str2num(get(handles.edit_mapWidth,'String'))*2.54,'%.2f' ))
    set(handles.edit_mapHeight,'String',num2str( str2num(get(handles.edit_mapHeight,'String'))*2.54,'%.2f' ))
    set(handles.edit_X0,'String',num2str( str2num(get(handles.edit_X0,'String'))*2.54,'%.2f' ))
    set(handles.edit_Y0,'String',num2str( str2num(get(handles.edit_Y0,'String'))*2.54,'%.2f' ))
elseif (strcmp(handles.which_unit,'in') & strcmp(dest,'pt'))
    xx = xx * 72;       yy = yy * 72;
    xf = xf * 72;       yf = yf * 72;
    set(handles.edit_mapWidth,'String',num2str( str2num(get(handles.edit_mapWidth,'String'))*72,'%.2f' ))
    set(handles.edit_mapHeight,'String',num2str( str2num(get(handles.edit_mapHeight,'String'))*72,'%.2f' ))
    set(handles.edit_X0,'String',num2str( str2num(get(handles.edit_X0,'String'))*72,'%.2f' ))
    set(handles.edit_Y0,'String',num2str( str2num(get(handles.edit_Y0,'String'))*72,'%.2f' ))
elseif (strcmp(handles.which_unit,'pt') & strcmp(dest,'cm'))
    xx = xx * 2.54/72;  yy = yy * 2.54/72;
    xf = xf * 2.54/72;  yf = yf * 2.54/72;
    set(handles.edit_mapWidth,'String',num2str( str2num(get(handles.edit_mapWidth,'String'))*2.54/72,'%.2f' ))
    set(handles.edit_mapHeight,'String',num2str( str2num(get(handles.edit_mapHeight,'String'))*2.54/72,'%.2f' ))
    set(handles.edit_X0,'String',num2str( str2num(get(handles.edit_X0,'String'))*2.54/72,'%.2f' ))
    set(handles.edit_Y0,'String',num2str( str2num(get(handles.edit_Y0,'String'))*2.54/72,'%.2f' ))
elseif (strcmp(handles.which_unit,'pt') & strcmp(dest,'in'))
    xx = xx / 72;       yy = yy / 72;
    xf = xf / 72;       yf = yf / 72;
    set(handles.edit_mapWidth,'String',num2str( str2num(get(handles.edit_mapWidth,'String'))/72,'%.2f' ))
    set(handles.edit_mapHeight,'String',num2str( str2num(get(handles.edit_mapHeight,'String'))/72,'%.2f' ))
    set(handles.edit_X0,'String',num2str( str2num(get(handles.edit_X0,'String'))/72,'%.2f' ))
    set(handles.edit_Y0,'String',num2str( str2num(get(handles.edit_Y0,'String'))/72,'%.2f' ))
end
set(handles.hand_rect,'XData',xx);      set(handles.hand_rect,'YData',yy);
if (~isempty(handles.hand_frame_proj))
    set(handles.hand_frame_proj,'XData',xf);      set(handles.hand_frame_proj,'YData',yf);
end

% Also uppdate the projection info text
str = get(handles.h_txt_info,'String');
try
	k = strfind(str{end},'/');
	if (~isempty(k))
        new_w = get(handles.edit_mapWidth,'String');
        str{end} = [str{end}(1:k(end)) new_w handles.which_unit(1)];
        set(handles.h_txt_info,'String',str)
	end
end

% -----------------------------------------------------------------------------------
function radiobutton_setWidth_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_setHeight,'Value',0)
else
    set(hObject,'Value',1)
end

% -----------------------------------------------------------------------------------
function radiobutton_setHeight_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_setWidth,'Value',0)
else
    set(hObject,'Value',1)
end

% -----------------------------------------------------------------------------------
function radiobutton_180_180_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_0_360,'Value',0)
else
    set(hObject,'Value',1)
end

% -----------------------------------------------------------------------------------
function radiobutton_0_360_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_180_180,'Value',0)
else
    set(hObject,'Value',1)
end

% -----------------------------------------------------------------------------------------
function move_rectangle(obj,eventdata,handles,h)
handles = guidata(handles.figure1);     % Uppdate handles
hf = handles.hand_frame_proj;
state = uisuspend(gcf);     % Remember initial figure state
x = get(h,'XData');     y = get(h,'YData');
if (~isempty(hf))
    xf = get(hf,'XData');   yf = get(hf,'YData');
    frame = [xf(:) yf(:)];  old_orig = [x(1) y(1)];
else
    frame = [];         old_orig = [];
end
WidthHeight = [x(3)-x(1) y(3)-y(1)];
set(gcf,'WindowButtonMotionFcn',{@wbm_MoveRectangle,handles,h,WidthHeight,hf,frame,old_orig},'WindowButtonDownFcn', ...
    {@wbd_MoveRectangle,handles,h,state}, 'Pointer', 'crosshair');

function wbm_MoveRectangle(obj,eventdata,handles,h,WidthHeight,hf,frame,old_orig)
pt = get(gca, 'CurrentPoint');
x = [pt(1,1) pt(1,1) pt(1,1)+WidthHeight(1) pt(1,1)+WidthHeight(1) pt(1,1)];
y = [pt(1,2) pt(1,2)+WidthHeight(2) pt(1,2)+WidthHeight(2) pt(1,2) pt(1,2)];
set(h, 'XData', x, 'YData', y);
set(handles.edit_X0,'String',num2str(pt(1,1),'%.2f'))
set(handles.edit_Y0,'String',num2str(pt(1,2),'%.2f'))
if (~isempty(hf))
    set(hf,'XData', frame(:,1)+pt(1,1)-old_orig(1), 'YData', frame(:,2)+pt(1,2)-old_orig(2));
end

function wbd_MoveRectangle(obj,eventdata,handles,h,state)
% check if x,y is inside of axis
pt = get(gca, 'CurrentPoint');  x = pt(1,1);    y = pt(1,2);
x_lim = get(gca,'xlim');      y_lim = get(gca,'ylim');
if (x<x_lim(1)) | (x>x_lim(2)) | (y<y_lim(1)) | (y>y_lim(2));   return; end
set(gcf,'WindowButtonMotionFcn','','WindowButtonDownFcn','', 'Pointer', 'arrow');
uirestore_j(state);           % Restore the figure's initial state
%set(handles.pushbutton_uppdate,'Visible','on')

% -----------------------------------------------------------------------------------
function pushbutton_uppdate_Callback(hObject, eventdata, handles)
xx = get(handles.hand_rect,'XData');
yy = get(handles.hand_rect,'YData');
set(handles.edit_X0,'String',num2str(xx(1),'%.1f'));
set(handles.edit_Y0,'String',num2str(yy(1),'%.1f'));

if (strcmp(handles.opt_J_no_scale(1:3),'-JX'))      % Linear proj has a different treatment
	scale_x = (xx(3) - xx(2)) / handles.width_or;
	scale_y = (yy(2) - yy(1)) / handles.height_or;
	new_y = handles.height_or * scale_x;
	new_x = handles.width_or * scale_y;
	if (get(handles.radiobutton_setWidth,'Value'))
        yy(2) = new_y + yy(1);      yy(3) = new_y + yy(1);
	elseif (get(handles.radiobutton_setHeight,'Value'))
        xx(3) = new_x + xx(2);      xx(4) = new_x + xx(1);
	else
        yy(2) = new_y + yy(1);      yy(3) = new_y + yy(1);  % It will become "True" scale
	end
	set(handles.hand_rect, 'XData', xx, 'YData', yy);
	set(handles.edit_mapWidth,'String',num2str((xx(3) - xx(2)),'%.2f'));    % Uppdate map width
	set(handles.edit_mapHeight,'String',num2str((yy(2) - yy(1)),'%.2f'));   % Uppdate map height
    set(handles.edit_scale,'String','1:1')
    % Also update the projection info text
    str = get(handles.h_txt_info,'String');
    try
        new_w = get(handles.edit_mapWidth,'String');
        str{end} = ['J<options> ->  -JX' new_w handles.which_unit(1)];
        set(handles.h_txt_info,'String',str)
    end
    return      % We are donne
end

new_w = num2str((xx(3) - xx(2)),'%.2f');
set(handles.edit_mapWidth,'String',new_w);
opt_J = [handles.opt_J_no_scale '/' new_w handles.which_unit(1)];
in = [handles.x_min handles.y_min; handles.x_min handles.y_max; handles.x_max handles.y_max; handles.x_max handles.y_min];
try
    opt_R = ['-R' num2str(handles.x_min,'%.6f') '/' num2str(handles.x_max,'%.6f') '/' ...
            num2str(handles.y_min,'%.6f') '/' num2str(handles.y_max,'%.6f')];
    out = mapproject_m(in,opt_R,opt_J,['-D' handles.which_unit(1)]);
catch
    return
end

% scale_prj = abs(out(2,2) - out(1,2)) / abs(out(4,1) - out(1,1));    % = dy / dx     NA TA SEMPRE CERTO
% scale_x = (xx(3) - xx(2)) / handles.width_or;
% scale_y = (yy(2) - yy(1)) / handles.height_or;
% new_y = handles.height_or * scale_x * scale_prj;
% new_x = handles.width_or * scale_y / scale_prj;

new_x = max(out(:,1)) - min(out(:,1));
new_y = max(out(:,2)) - min(out(:,2));

if (get(handles.radiobutton_setWidth,'Value'))
    yy(2) = new_y + yy(1);      yy(3) = new_y + yy(1);
elseif (get(handles.radiobutton_setHeight,'Value'))
    xx(3) = new_x + xx(2);      xx(4) = new_x + xx(1);
else
    yy(2) = new_y + yy(1);      yy(3) = new_y + yy(1);  % It will become "True" scale
end
set(handles.hand_rect, 'XData', xx, 'YData', yy);
set(handles.edit_mapWidth,'String',num2str((xx(3) - xx(2)),'%.2f'));    % Uppdate map width
set(handles.edit_mapHeight,'String',num2str((yy(2) - yy(1)),'%.2f'));   % Uppdate map height
handles.scale = num2str((xx(3) - xx(2)),'%.2f');

% --- Compute a projected mini frame
n = 21;
xf = linspace(handles.x_min,handles.x_max,n)';
yf = linspace(handles.y_min,handles.y_max,n)';
in = [repmat(xf(1),n,1) yf; xf(2:end) repmat(yf(end),n-1,1); repmat(xf(end),n-1,1) yf(end-1:-1:1); ...
        xf(end:-1:1) repmat(yf(1),n,1)];
%out_f = mapproject_m(in,opt_R,'-C','-F',[handles.opt_J_no_scale '/' handles.scale handles.which_unit(1)]);
out_f = mapproject_m(in,opt_R,opt_J,['-D' handles.which_unit(1)]);
new_x = xx(1) + out_f(:,1);
new_y = yy(1) + out_f(:,2);

% Draw it if it's not a rectangle
if ~(out_f(1,1) == out_f(n,1) & out_f(n,2) == out_f(2*n-1,2))
    if (isempty(handles.hand_frame_proj))      % First time. Creat it.
        handles.hand_frame_proj = line('XData',new_x,'YData',new_y, 'Color','r','LineWidth',.5,'Tag','PlotFrameProj');
        uistack(handles.hand_frame_proj, 'down')
    else
        set(handles.hand_frame_proj, 'XData', new_x, 'YData', new_y);
    end
else    % It is a rectangle
    if (~isempty(handles.hand_frame_proj))  % If we have a previous red frame, delete it
        delete(handles.hand_frame_proj);  handles.hand_frame_proj = [];
    end
end

% ----------- Compute scale 1:xxxx
if (~handles.scale_set)     % If user changed scale, don't compute it here
	xm = (handles.x_min + handles.x_max) / 2;   ym = (handles.y_min + handles.y_max) / 2;
	opt_R = ['-R' num2str(xm-2) '/' num2str(xm+2) '/' num2str(ym-2) '/' num2str(ym+2)];
	in = [handles.x_min handles.y_min; handles.x_min handles.y_max; handles.x_max handles.y_max; handles.x_max handles.y_min];
	opt_J = [handles.opt_J_no_scale '/1'];
	out = mapproject_m(in,opt_R,'-C','-F',opt_J);
	dx_prj = out(4,1) - out(1,1);   % It's in projected meters
	dy_prj = out(2,2) - out(1,2);   % It's in projected meters
	dx_rect = xx(4) - xx(1);        % Is in "cm", "in" or "pt". So convert to "cm"
	dy_rect = yy(2) - yy(1);        % Is in "cm", "in" or "pt". So convert to "cm"
	if (strcmp(handles.which_unit,'in'))      dx_rect = dx_rect * 2.54;     dy_rect = dy_rect * 2.54;     end
	if (strcmp(handles.which_unit,'pt'))      dx_rect = dx_rect * 2.54/72;  dy_rect = dy_rect * 2.54/72;  end
	scale = max(dx_rect/dx_prj/100, dy_rect/dy_prj/100);
	[n,d] = rat(scale,1e-9);
	if (n > 1)     d = d / n;      end
	set(handles.edit_scale,'String',['1:' num2str(d)])
	handles.scale_set = 0;
end

% ------------ Also uppdate the projection info text
str = get(handles.h_txt_info,'String');
try
    new_w = get(handles.edit_mapWidth,'String');
	k = strfind(str{end},'/');
	if (~isempty(k))
        str{end} = [str{end}(1:k(end)) new_w handles.which_unit(1)];
    else        % We have a -J without any '/'. Linear projection
        k = strfind(str{end},'-J');
        str{end} = [str{end}(1:k+2) new_w handles.which_unit(1)];
	end
    set(handles.h_txt_info,'String',str)
end
%zoom on
guidata(hObject,handles)
%set(hObject,'Visible','off')

% -----------------------------------------------------------------------------------
function edit_X0_Callback(hObject, eventdata, handles)
% Set new x origin
str = get(hObject,'String');        x0 = str2double(str);
xx = get(handles.hand_rect,'XData');
if (isnan(x0))      set(hObject,'String',str);      return;     end
set(handles.hand_rect,'XData',xx - xx(1) + x0)
if (~isempty(handles.hand_frame_proj))
    set(handles.hand_frame_proj,'XData',get(handles.hand_frame_proj,'XData') - xx(1) + x0)
end

% -----------------------------------------------------------------------------------
function edit_Y0_Callback(hObject, eventdata, handles)
% Set new y origin
str = get(hObject,'String');        y0 = str2double(str);
yy = get(handles.hand_rect,'YData');
if (isnan(y0))      set(hObject,'String',str);      return;     end
set(handles.hand_rect,'YData',yy - yy(1) + y0)
if (~isempty(handles.hand_frame_proj))
    set(handles.hand_frame_proj,'YData',get(handles.hand_frame_proj,'YData') - yy(1) + y0)
end

% -----------------------------------------------------------------------------------
function edit_mapWidth_Callback(hObject, eventdata, handles)
% Set new map width
str = get(hObject,'String');        w = str2double(str);
xx = get(handles.hand_rect,'XData');
if (isnan(w))      set(hObject,'String',str);      return;     end
xx(3) = xx(2) + w;      xx(4) = xx(1) + w;
set(handles.hand_rect,'XData',xx)
pushbutton_uppdate_Callback(handles.pushbutton_uppdate, eventdata, handles)

% -----------------------------------------------------------------------------------
function edit_mapHeight_Callback(hObject, eventdata, handles)
% Set new map height
str = get(hObject,'String');        h = str2double(str);
yy = get(handles.hand_rect,'YData');
if (isnan(h))      set(hObject,'String',str);      return;     end
yy(2) = yy(1) + h;      yy(3) = yy(4) + h;
set(handles.hand_rect,'YData',yy)
pushbutton_uppdate_Callback(handles.pushbutton_uppdate, eventdata, handles)

% -----------------------------------------------------------------------------------------
function popup_directory_list_Callback(hObject, eventdata, handles)
if (nargin == 3)    opt = [];   end
val = get(hObject,'Value');     str = get(hObject, 'String');
% Put the selected field on top of the String list.
tmp = str(val);         str(val) = [];
new_str = [tmp; str];   set(hObject,'String',new_str,'Value',1); 

% -----------------------------------------------------------------------------------------
function pushbutton_change_dir_Callback(hObject, eventdata, handles)
if (ispc)
    contents = get(handles.popup_directory_list,'String');
    if (iscell(contents))
        pato = contents{1};         % Start at default's directory
    else
        pato = contents;            % Start at default's directory
    end
    work_dir = uigetfolder_standalone('Select scripts folder',pato);
else            % This guy doesn't let to be compiled
    work_dir = uigetdir;
end
if ~isempty(work_dir)
    handles.last_directories = [cellstr(work_dir); handles.last_directories];
    set(handles.popup_directory_list,'String',handles.last_directories)
    guidata(hObject, handles);
end

% -----------------------------------------------------------------------------------------
function edit_prefix_Callback(hObject, eventdata, handles)
% Nothing to do

% -----------------------------------------------------------------------------------
function edit_scale_Callback(hObject, eventdata, handles)
str = get(hObject,'String');
k = strfind(str,':');
xx = get(handles.hand_rect,'XData');    yy = get(handles.hand_rect,'YData');
new_w = num2str((xx(3) - xx(2)),'%.2f');
opt_J = [handles.opt_J_no_scale '/' str];   opt_J(3) = lower(opt_J(3));
in = [handles.x_min handles.y_min; handles.x_min handles.y_max; handles.x_max handles.y_max; handles.x_max handles.y_min];
try
    opt_R = ['-R' num2str(handles.x_min,'%.6f') '/' num2str(handles.x_max,'%.6f') '/' ...
            num2str(handles.y_min,'%.6f') '/' num2str(handles.y_max,'%.6f')];
    out = mapproject_m(in,opt_R,opt_J,['-D' handles.which_unit(1)]);
catch    return;    end
xmax = max(out(:,1));   ymax = max(out(:,2));
xx(3) = xmax+xx(1);     xx(4) = xmax+xx(1);
yy(2) = ymax+yy(1);     yy(3) = ymax+yy(1);
set(handles.hand_rect,'XData',xx,'YData',yy)
handles.scale_set = 1;
guidata(hObject, handles);
pushbutton_uppdate_Callback(handles.pushbutton_uppdate, eventdata, handles)

% -----------------------------------------------------------------------------------
function pushbutton_mapProjections_Callback(hObject, eventdata, handles)
if (~handles.mirone_handles.geog)
    msg = ['I will tell you a secret. A map projection is an operation where GEOGRAPHIC ' ...
            'coordinates (representing a nearly spherical surface) are transformed into ' ...
            'planar (flat) coordinates. You got the message?'];
    warndlg(msg,'Warning');     return
end
fname = [handles.d_path 'mirone_pref.mat'];
coord_system_script = coordinate_system(handles.coord_system_script,handles.all_datums,[]);
if (isempty(coord_system_script)),   return;     end
handles.coord_system_script = coord_system_script;

% Split the scale from the projection string
tmp = coord_system_script.projection;
if (length(tmp) == 4 & strcmp(tmp(3),'m'))     % Simple Mercator cames in the form "-Jm1"
    tmp = tmp(1:end-1);
elseif (length(tmp) == 3 & strcmp(tmp(3),'x')) % Linear proj cames in the form "-Jx"
    tmp = tmp;
else                                            % All other should terminate as "-J.../1"
    tmp = tmp(1:end-2);
end
xx = get(handles.hand_rect,'XData');        yy = get(handles.hand_rect,'YData');
handles.scale = num2str( (xx(3) - xx(2)),'%.2g');
if (length(tmp) > 3)
    opt_J = [tmp(1:2) upper(tmp(3)) tmp(4:end) '/' handles.scale handles.which_unit(1)];
    handles.opt_J_no_scale = [tmp(1:2) upper(tmp(3)) tmp(4:end)];           % Save this
else        % Linear projections
    opt_J = [tmp(1:2) upper(tmp(3)) handles.scale handles.which_unit(1)];
    handles.opt_J_no_scale = [tmp(1:2) upper(tmp(3))];                      % Save this
end
set(handles.edit_mapWidth,'String',handles.scale)
handles.curr_datum = handles.all_datums{coord_system_script.datum_val,2};   % Save this

string = {['System   -> ' coord_system_script.SysName];...
        ['Projection -> ' coord_system_script.ProjName];...
        ['Ellipsoid  ->  ' handles.all_datums{coord_system_script.datum_val,2}];...
        ['J<options> ->  ' opt_J]};
[outstring,newpos] = textwrap(handles.h_txt_info,string);
pos = handles.txt_info_pos;
pos(4) = newpos(4);
set(handles.h_txt_info,'String',outstring,'Position',[pos(1),pos(2),pos(3),pos(4)])
coord_system_script.proj_info_txt = outstring;
coord_system_script.proj_info_pos = pos;

save(fname,'coord_system_script','-append');      % Update mirone_pref
guidata(hObject,handles)
pushbutton_uppdate_Callback(handles.pushbutton_uppdate, eventdata, handles)

% ----------------------------------------------------------------------------------------
function togglebutton_Option_L_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    xx = draw_scale;
    if (~isempty(xx))
        set(handles.checkbox_removeOptionL,'Visible','on','Value',0)
        handles.opt_L = [' ' xx];
    else
        set(hObject,'Value',0)
    end
else
    handles.opt_L = [];
    set(handles.checkbox_removeOptionL,'Visible','off')
end
guidata(hObject, handles);

% ----------------------------------------------------------------------------------------
function togglebutton_Option_U_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    xx = time_stamp;
    if (~isempty(xx))
        set(handles.checkbox_removeOptionU,'Visible','on','Value',0)
        handles.opt_U = [' ' xx];
    else
        set(hObject,'Value',0)
    end
else
    handles.opt_U = [];
    set(handles.checkbox_removeOptionU,'Visible','off')
end
guidata(hObject, handles);

% ----------------------------------------------------------------------------------------
function checkbox_removeOptionL_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.togglebutton_Option_L,'Value',0)    % Reset the associated togglebutton to off
    set(hObject,'Visible','off')                    % Hide this button
    handles.opt_L = [];
    guidata(hObject, handles);
end

% -----------------------------------------------------------------------------------------
function checkbox_removeOptionU_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.togglebutton_Option_U,'Value',0)    % Reset the associated togglebutton to off
    set(hObject,'Visible','off')                    % Hide this button
    handles.opt_U = [];
    guidata(hObject, handles);
end

% -----------------------------------------------------------------------------------
function pushbutton_coastLines_Callback(hObject, eventdata, handles)
handles.opt_psc = pscoast_options_Mir(handles.mirone_handles, handles.psc_res, handles.psc_opt_W, ...
    handles.psc_type_p, handles.psc_type_r);
    % Testar se ha -N e -I repetidas
guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
% Here we transmit the: -J<...>, paper name, files prefix, X0, Y0 and destination directory
if (~strcmp(handles.opt_J_no_scale(1:3),'-JX'))
    opt_J = [handles.opt_J_no_scale '/' handles.scale handles.which_unit(1)];
else        % Linear projection
    opt_J = [handles.opt_J_no_scale handles.scale handles.which_unit(1)];
end
val = get(handles.popup_PaperSize,'Value');
list = get(handles.popup_PaperSize,'String');
str = list{val};        k = strfind(str,' ');
paper = str(1:k(1)-1);
d_dir = get(handles.popup_directory_list,'String');
if (iscell(d_dir))      d_dir = d_dir{1};    end
prefix = get(handles.edit_prefix,'String');

X0 = get(handles.edit_X0,'String');     Y0 = get(handles.edit_Y0,'String');
X0 = ['-X' X0 handles.which_unit(1)];   Y0 = ['-Y' Y0 handles.which_unit(1)];

if (get(handles.radiobutton_180_180,'Value'))   % [-180;180] range
    opt_deg = '--PLOT_DEGREE_FORMAT=ddd:mm:ss';
else                                            % [0;360] range
    opt_deg = '--PLOT_DEGREE_FORMAT=+ddd:mm:ss';
end

%try
    % Before calling the write script routine we have to find if we have any pscoast stuff
    if (isempty(handles.opt_psc))       % Means that the pscoast_options was not used
        if (~isempty(handles.psc_res))  % Means that we have coastlines and will use the Mirone settings
            handles.opt_psc = [handles.psc_res ' ' handles.psc_opt_W ' ' handles.psc_type_p ' ' handles.psc_type_r];
        end
    end
    if (get(handles.radiobutton_P,'Value'))     opt_P = ' -P';
    else                                        opt_P = '';
    end
    out_msg = build_write_script(handles.mirone_handles,handles.ALLlineHand,handles.opt_R, opt_J, d_dir, prefix, ...
        paper, X0, Y0, handles.script_type, handles.curr_datum, handles.opt_L, handles.opt_U, opt_P, ...
        handles.opt_psc, opt_deg);
    msg{1} = ['File ' prefix '_mir.' handles.script_type ' successufuly created in:  ' d_dir];
    if (out_msg)
        msg{2} = [];
        msg{3} = 'WARNING: Read the important message on the header of the script';
    end
    msgbox(msg);
%catch
    %errordlg(['An unknown error occured while writing the ' handles.script_type ' file'],'Error');
%end

% --------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1)

%-------------------------------------------------------------------------------------
function [handles,out] = check_coord_system(handles,coord_system,side)
% Currently the coord_system structure must contain the following fields:
%              group_val
%             system_val
%              datum_val
%        cilindrical_val
%          azimuthal_val
%              conic_val
%       miscelaneous_val
%               ProjName    % Projection Name (used in the text info)
%       map_scale_factor
%        MeasureUnit_val    % Currently selected measure unitie
%      DegreeFormat1_val    % Degree format value for point conversions
%      DegreeFormat2_val    % Degree format value for file conversions
%           system_FE_FN    % False eastings/northings
%                is_geog    % Signals if we have geogs coordinates
%             projection    % -J string
%     ProjParameterValue    % projection parameters
%          proj_info_txt
%          proj_info_pos    % This one has to be dealt inside the opening function (don't know here a def value)

if (isempty(coord_system))   % If it doesn't exist, create an empty one
    coord_system = struct([]);
end

% If any of those is missing, assign it a default value
if (~isfield(coord_system,'group_val'))     out.group_val = 1;
else        out.group_val = coord_system.group_val;    end
if (~isfield(coord_system,'system_val'))    out.system_val = 1;
else        out.system_val = coord_system.system_val;    end
if (~isfield(coord_system,'datum_val'))     out.datum_val = 221;   % Default to wgs84
else        out.datum_val = coord_system.datum_val;    end
if (~isfield(coord_system,'cilindrical_val'))   out.cilindrical_val = 1;
else        out.cilindrical_val = coord_system.cilindrical_val;    end
if (~isfield(coord_system,'azimuthal_val'))     out.azimuthal_val = 1;
else        out.azimuthal_val = coord_system.azimuthal_val;    end
if (~isfield(coord_system,'conic_val'))         out.conic_val = 1;
else        out.conic_val = coord_system.conic_val;    end
if (~isfield(coord_system,'miscelaneous_val'))  out.miscelaneous_val = 1;
else        out.miscelaneous_val = coord_system.miscelaneous_val;    end
if (~isfield(coord_system,'ProjName'))          out.ProjName = 'Unknown';
else        out.ProjName = coord_system.ProjName;    end
if (~isfield(coord_system,'map_scale_factor'))  out.map_scale_factor = [];
else        out.map_scale_factor = coord_system.map_scale_factor;    end
if (~isfield(coord_system,'system_FE_FN'))      out.system_FE_FN = [];
else        out.system_FE_FN = coord_system.system_FE_FN;    end
if (~isfield(coord_system,'projection'))        out.projection = [];
else        out.projection = coord_system.projection;    end
if (~isfield(coord_system,'ProjParameterValue')) out.ProjParameterValue = [];
else        out.ProjParameterValue = coord_system.ProjParameterValue;    end
if (~isfield(coord_system,'proj_info_txt'))     out.proj_info_txt = 'Nikles';
else        out.proj_info_txt = coord_system.proj_info_txt;    end
if (~isfield(coord_system,'MeasureUnit_val'))   out.MeasureUnit_val = 1;
else        out.MeasureUnit_val = coord_system.MeasureUnit_val;    end
if (~isfield(coord_system,'DegreeFormat1_val'))   out.DegreeFormat1_val = 1;
else        out.DegreeFormat1_val = coord_system.DegreeFormat1_val;    end
if (~isfield(coord_system,'DegreeFormat2_val'))   out.DegreeFormat2_val = 1;
else        out.DegreeFormat2_val = coord_system.DegreeFormat2_val;    end
if (~isfield(coord_system,'is_geog'))   out.is_geog = 1;
else        out.is_geog = coord_system.is_geog;    end

% This is my solution to cat 2 structures. There must be a clever way.
hand_cell = struct2cell(handles);       % Convert handles struct to cell
out_cell = struct2cell(out);            % Convert coord_system struct to cell
both_cell = [hand_cell; out_cell];      % Cat them
names_hand = fieldnames(handles);       % Get handles field names
names_out = fieldnames(out);            % Get coord_system field names
for (i=1:length(names_out))             % Append the "side" to the coord_system field names
    names_out{i} = [names_out{i} side];
end
both_names = [names_hand; names_out];   % Cat the handles and the coord_system field names
handles = cell2struct(both_cell,both_names,1);  % Finaly, rebuild the handles structure.

% ----------------------------------------------------------------------------------
function [ALLlineHand, res, opt_W, type_p, type_r] = find_psc_stuff(ALLlineHand)
% See if we have any pscoast stuff
haveSymbol = 0;     haveCoasts = 0;     havePolitical = 0;  haveRivers = 0;
res = [];           opt_W = [];         type_p = [];        type_r = [];
h_c = findobj(ALLlineHand,'Tag','CoastLineNetCDF');
if (~isempty(h_c))
    if (length(h_c) > 1)    h_c = h_c(1);     end
    CoastRes = get(h_c,'UserData');
    LineWidth_c = get(h_c,'LineWidth');
    LineColor_c = get(h_c,'Color');
    LineStyle_c = get(h_c,'LineStyle');
    haveCoasts = 1;
end
h_p = findobj(ALLlineHand,'Tag','PoliticalBoundaries');
if (~isempty(h_p))
    if (length(h_p) > 1)    h_p = h_p(1);     end
    zz = get(h_p,'UserData');
    if (iscell(zz))     zz = zz{1};     end
    PoliticalRes = zz(1);        PoliticalType = zz(2);
    LineWidth_p = get(h_p,'LineWidth');
    LineColor_p = get(h_p,'Color');
    LineStyle_p = get(h_p,'LineStyle');
    havePolitical = 1;
end
h_r = findobj(ALLlineHand,'Tag','Rivers');
if (~isempty(h_r))
    if (length(h_r) > 1)    h_r = h_r(1);     end
    zz = get(h_r,'UserData');
    if (iscell(zz))         zz = zz{1};     end
    RiversRes = zz(1);          RiversType = zz(2);
    LineWidth_r = get(h_r,'LineWidth');
    LineColor_r = get(h_r,'Color');
    LineStyle_r = get(h_r,'LineStyle');
    haveRivers = 1;
end
ALLlineHand = setxor(ALLlineHand, [h_c; h_p; h_r]);

if (haveCoasts | havePolitical | haveRivers)
    res_c = '';     res_p = '';     res_r = '';
    if (haveCoasts)
        cor = round(LineColor_c * 255);
        cor = [num2str(cor(1)) '/' num2str(cor(2)) '/' num2str(cor(3))];
        switch CoastRes
            case 'f',   res_c = ['-Df -W' num2str(LineWidth_c) 'p/' cor];
            case 'h',   res_c = ['-Dh -W' num2str(LineWidth_c) 'p/' cor];
            case 'i',   res_c = ['-Di -W' num2str(LineWidth_c) 'p/' cor];
            case 'l',   res_c = ['-Dl -W' num2str(LineWidth_c) 'p/' cor];
            case 'c',   res_c = ['-Dc -W' num2str(LineWidth_c) 'p/' cor];
        end
        if (~strcmp(LineStyle_c,'-'))   % If we have a line style other than solid
            switch LineStyle_c
                case '--',  res_c = [res_c 'ta'];
                case ':',   res_c = [res_c 'to'];
                case '-.',  res_c = [res_c 't10_2_2_5:5'];
            end
        end
    end
    if (havePolitical)
        switch PoliticalRes
            case 'f',   res_p = ['-Df'];
            case 'h',   res_p = ['-Dh'];
            case 'i',   res_p = ['-Di'];
            case 'l',   res_p = ['-Dl'];
            case 'c',   res_p = ['-Dc'];
        end
        cor = round(LineColor_p * 255);
        cor = [num2str(cor(1)) '/' num2str(cor(2)) '/' num2str(cor(3))];
        switch PoliticalType
            case '1',   type_p = ['-N1/'  num2str(LineWidth_p) 'p/' cor];
            case '2',   type_p = ['-N2/'  num2str(LineWidth_p) 'p/' cor];
            case '3',   type_p = ['-N3/'  num2str(LineWidth_p) 'p/' cor];
            case 'a',   type_p = ['-Na/'  num2str(LineWidth_p) 'p/' cor];
        end
        if (~strcmp(LineStyle_p,'-'))   % If we have a line style other than solid
            switch LineStyle_p
                case '--',  type_p = [type_p 'ta'];
                case ':',   type_p = [type_p 'to'];
                case '-.',  type_p = [type_p 't10_2_2_5:5'];
            end
        end
    end
    if (haveRivers)
        switch RiversRes
            case 'f',   res_r = ['-Df'];
            case 'h',   res_r = ['-Dh'];
            case 'i',   res_r = ['-Di'];
            case 'l',   res_r = ['-Dl'];
            case 'c',   res_r = ['-Dc'];
        end
        cor = round(LineColor_r * 255);
        cor = [num2str(cor(1)) '/' num2str(cor(2)) '/' num2str(cor(3))];
        switch RiversType
            case '1',   type_r = ['-I1/' num2str(LineWidth_r) 'p/' cor];
            case '2',   type_r = ['-I2/' num2str(LineWidth_r) 'p/' cor];
            case '3',   type_r = ['-I3/' num2str(LineWidth_r) 'p/' cor];
            case '4',   type_r = ['-I4/' num2str(LineWidth_r) 'p/' cor];
            case '5',   type_r = ['-I5/' num2str(LineWidth_r) 'p/' cor];
            case '6',   type_r = ['-I6/' num2str(LineWidth_r) 'p/' cor];
            case '7',   type_r = ['-I7/' num2str(LineWidth_r) 'p/' cor];
            case '8',   type_r = ['-I8/' num2str(LineWidth_r) 'p/' cor];
            case '9',   type_r = ['-I9/' num2str(LineWidth_r) 'p/' cor];
            case '10',  type_r = ['-I10/' num2str(LineWidth_r) 'p/' cor];
            case 'a',   type_r = ['-Ia/' num2str(LineWidth_r) 'p/' cor];
            case 'r',   type_r = ['-Ir/' num2str(LineWidth_r) 'p/' cor];
            case 'i',   type_r = ['-Ii/' num2str(LineWidth_r) 'p/' cor];
            case 'c',   type_r = ['-Ic/' num2str(LineWidth_r) 'p/' cor];
        end
            if (~strcmp(LineStyle_r,'-'))   % If we have a line style other than solid
            switch LineStyle_r
                case '--',  type_r = [type_r 'ta'];
                case ':',   type_r = [type_r 'to'];
                case '-.',  type_r = [type_r 't10_2_2_5:5'];
            end
        end
end
    res = unique([res_c(1:3); res_p; res_r],'rows');  % We don't want repeated resolution strings
    if (size(res,1) > 1)        % Shit, we have mixed resolutions
        res = '-Di';            % TEMPORARY SOLUTION UNTIL I FIND HOW TO FIND THE HIGHEST COMMON RES
    end
    if (~isempty(res_c))    opt_W = res_c(5:end);
    else                    opt_W = [];     end
end

% --------------------------------------------------------------------------------------------------------
function out_msg = build_write_script(handles_mirone, ALLlineHand, opt_R, opt_J, dest_dir, ...
    prefix, paper, X0, Y0, sc, ellips, opt_L, opt_U, opt_P, opt_psc, opt_deg)
% This function do most of the hard work in finding the script components. The pscoast stuff is
% worked out by the "find_psc_stuff" function.

if (isempty(opt_psc))   have_psc = 0;       % We do not have any pscoast commands
else                    have_psc = 1;       end

if (~strcmp(paper,'A4'))    paper_media = paper;
else                        paper_media = [];   end
if (strcmp(sc,'bat'))
    comm = 'REM ';      pb = '%';   pf = '%';
else
    comm = '# ';        pb = '$';   pf = '';
end
if (strcmp(ellips,'WGS-84'))     % It is the default, so don't use any
    ellips = [];
else
    ellips = [' --ELLIPSOID=' ellips];
end

% ------------ Some (maybe) needed vars ------------------------------------------------------------------
haveSymbol = 0;     i_c = [];   used_grd = 0;   sc_cpt = [dest_dir filesep prefix '.cpt'];  out_msg = 0;
need_path = 0;      used_countries = 0;
script = cell(16,1);
if (~isempty(handles_mirone.grdname))
    [PATH,FNAME,EXT] = fileparts(handles_mirone.grdname);
    just_grd_name = [FNAME EXT];
    if (strcmp(PATH,dest_dir))      need_path = 0;
    else                            need_path = 1;  end
    clear PATH FNAME EXT;
else
    need_path = 0;
end
grd_name = handles_mirone.grdname;

% -------------------- Build -B string -------------------------------------------------
try
	h_axes = findobj(handles_mirone.figure1,'Type','Axes');
	Bx = get(h_axes,'XTick');      d_Bx = diff(Bx);
	By = get(h_axes,'YTick');      d_By = diff(By);
	opt_B = ['-B' num2str(d_Bx(1)) '/' num2str(d_By(1))];
	clear h_axes Bx By d_Bx d_By;
catch
    opt_B = '-B1000000';    % invented value
end
% --------------------------------------------------------------------------------------

l = 1;
if (~strcmp(sc,'bat'))                          % Write a csh script
    script{l} = ['#!/bin/csh -f'];              l=l+1;
	script{l} = comm;                           l=l+1;
	script{l} = [comm 'Coffeeright Mirone Tec'];l=l+1;
	script{l} = comm;                           l=l+1;
	script{l} = [comm ' ---- Projection. You may change it if you know how to'];    l=l+1;
	script{l} = ['set proj = ' opt_J];          l=l+1;      % Map scale
	script{l} = [comm ' ---- Frame annotations. You may change it if you know how to'];    l=l+1;
	script{l} = ['set frm = ' opt_B];           l=l+1;
	script{l} = [comm ' ---- Map limits. You may change it if you know how to'];    l=l+1;
	script{l} = ['set lim = ' opt_R];           l=l+1;
	script{l} = comm;                           l=l+1;
	script{l} = [comm ' ---- Longitude annotation style. The +ddd:mm:ss form => [0;360] range '];    l=l+1;
	script{l} = ['set deg_form=' opt_deg];      l=l+1;
	script{l} = '';                             l=l+1;
    prefix_ddir = [dest_dir filesep prefix];    % Add destination dir to the name prefix
    if (~isempty(grd_name))
        if (~need_path)
    	    script{l} = ['set grd = ' just_grd_name];   id_grd = l; l=l+1;
        else
    	    script{l} = ['set grd = ' grd_name];        id_grd = l; l=l+1;
        end
    end
    script{l} = ['set cpt = ' prefix '.cpt']; id_cpt = l;   l=l+1;
	script{l} = ['set ps = ' prefix '.ps'];     l=l+1;
    if (~isempty(paper_media))
    	script{l} = [comm ' We are not using A4'];  l=l+1;
        script{l} = ['gmtset PAPER_MEDIA=' paper_media]; l=l+1;
    	script{l} = comm;                       l=l+1;        
    end
else                                            % Write a dos batch    
	script{l} = ['@echo OFF'];                  l=l+1;
	script{l} = [comm 'Coffeewrite Mirone Tec'];l=l+1;
	script{l} = comm;                           l=l+1;
	script{l} = [comm ' ---- Projection. You may change it if you know how to'];    l=l+1;
	script{l} = ['set proj=' opt_J];            l=l+1;      % Map scale
	script{l} = [comm ' ---- Frame annotations. You may change it if you know how to'];    l=l+1;
	script{l} = ['set frm=' opt_B];             l=l+1;
	script{l} = [comm ' ---- Map limits. You may change it if you know how to'];    l=l+1;
	script{l} = ['set lim=' opt_R];             l=l+1;
	script{l} = comm;                           l=l+1;
	script{l} = [comm ' ---- Longitude annotation style. The +ddd:mm:ss form => [0;360] range '];    l=l+1;
	script{l} = ['set deg_form=' opt_deg];      l=l+1;
	script{l} = '';                             l=l+1;
    prefix_ddir = [dest_dir filesep prefix];    % Add destination dir to the name prefix
    if (~isempty(grd_name))
        if (~need_path)
	        script{l} = ['set grd=' just_grd_name];     id_grd = l; l=l+1;
        else
	        script{l} = ['set grd=' grd_name];          id_grd = l; l=l+1;
        end
    end
    script{l} = ['set cpt=' prefix '.cpt'];     id_cpt = l; l=l+1;
	script{l} = ['set ps=' prefix '.ps'];       l=l+1;
    if (~isempty(paper_media))
    	script{l} = [comm ' ---- We are not using A4'];  l=l+1;
        script{l} = ['gmtset PAPER_MEDIA=' paper_media]; l=l+1;
    	script{l} = comm;                       l=l+1;        
    end
end

% ------------ Start writing GMT commands --------------------------------
script{l} = [' '];                              l=l+1;
script{l} = [comm '-------- Start by creating the basemap frame'];  l=l+1;
script{l} = ['psbasemap ' pb 'lim' pf ' ' pb 'proj' pf ' ' pb 'frm' pf ' ' X0 ' ' Y0 opt_U opt_P ' ' pb 'deg_form' pf ' -K > ' pb 'ps' pf];
l=l+1;
if (~isempty(grd_name))
    if ( handles_mirone.Illumin_type > 0 && handles_mirone.Illumin_type <= 4 )
        % We have a image illuminated with grdgradient. Rebuild de illumination
        illumComm = getappdata(handles_mirone.figure1,'illumComm');
        opt_M = '';
        if (handles_mirone.Illumin_type == 1 && handles_mirone.geog),   opt_M = ' -M';     end
        opt_N = '';
        if (handles_mirone.Illumin_type == 1),      opt_N = ' -Nt';     end
        name_illum = [prefix '_intens.grd=cf'];
        script{l} = [''];                          l=l+1;
        script{l} = [comm '-------- Compute the illumination grid'];    l=l+1;
        script{l} = ['grdgradient ' pb 'grd' pf opt_M ' ' illumComm opt_N ' -G' name_illum ellips];    l=l+1;
        have_gmt_illum = 1;     used_grd = 1;
        illum = [' -I' name_illum];
    elseif ( handles_mirone.Illumin_type > 4 || handles_mirone.is_draped )
        % We have a Manip or draping illumination. Here we have to use the R,G,B trick
        name = [prefix_ddir '_channel'];    name_sc = [prefix '_channel'];
        mirone('File_img2GMT_RGBgrids_CB',gcbo,[],handles_mirone,'image',name)
        illum = [name '_r.grd ' name '_g.grd ' name '_b.grd']; % ????
        have_gmt_illum = 0;
    else        % We don't have any illumination
        have_gmt_illum = 0;
        used_grd = 1;
    end
    script{l} = [' '];                          l=l+1;
    if (have_gmt_illum)                     % grdimage with illumination
        script{l} = [comm '-------- Plot the the base image using grdimage & illumination'];    l=l+1;
        script{l} = ['grdimage ' pb 'grd' pf ' -R -J -C' pb 'cpt' pf illum ellips opt_L ' -O -K >> ' pb 'ps' pf];
        l=l+1;
        used_grd = 1;
    elseif (used_grd & ~have_gmt_illum)     % Simple grdimage call
        script{l} = [comm '-------- Plot the the base image using grdimage'];    l=l+1;
        script{l} = ['grdimage ' pb 'grd' pf ' -R -J -C' pb 'cpt' pf ellips opt_L ' -O -K >> ' pb 'ps' pf];   l=l+1;
        used_grd = 1;
    else                                    % No grd used, use the R,G,B channels
        script{l} = [comm '-------- Plot the 3 RGB base images using grdimage'];    l=l+1;
        script{l} = ['grdimage ' name_sc '_r.grd ' name_sc '_g.grd ' name_sc '_b.grd' ellips ' -R -J -O -K >> ' pb 'ps' pf];
        l=l+1;    
    end
    clear grd_name have_gmt_illum illum;
elseif (handles_mirone.image_type == 20)
    % Do nothing regarding the basemap image (in fact we don't have any image)
else    % We don't have a grid, so we need to fish the image and save it as R,G,B triplet
    name = [prefix_ddir '_channel'];    name_sc = [prefix '_channel'];
    mirone('File_img2GMT_RGBgrids_CB',gcbo,[],handles_mirone,'image',name) 
    script{l} = [' '];
    script{l} = [comm '-------- Plot the 3 RGB base images using grdimage'];    l=l+1;
    script{l} = ['grdimage ' name_sc '_r.grd ' name_sc '_g.grd ' name_sc '_b.grd' ellips ' -R -J -O -K >> ' pb 'ps' pf];
    l=l+1;    
end

% ------------ If we have used a GMT grid file build the GMT palette -----------------------
if (used_grd)
    tmp = cell(261,1);
    pal = colormap(handles_mirone.axes1);
    Z = getappdata(handles_mirone.figure1,'dem_z');
    % SE Z == [] FAZER QUALQUER COISA
    if (handles_mirone.have_nans)      cor_nan = pal(1,:);     pal = pal(2:end,:);   end     % Remove the bg color
    
    pal_len = size(pal,1);
    z_min = handles_mirone.head(5);    z_max = handles_mirone.head(6);
    
    dz = (z_max - z_min) / pal_len;
    tmp{1} = ['# Color palette exported by Mirone'];
    tmp{2} = ['# COLOR_MODEL = RGB'];
    cor = round(pal*255);
	for i=1:pal_len
        cor_str = sprintf([num2str(cor(i,1),'%.12g') '\t' num2str(cor(i,2),'%.12g') '\t' num2str(cor(i,3),'%.12g')]);
        z1 = num2str(z_min+dz*(i-1),'%.3f');
        z2 = num2str(z_min+dz*i,'%.3f');
        tmp{i+2} = sprintf([z1 '\t' cor_str '\t' z2 '\t' cor_str]);
	end
    tmp{pal_len+3} = sprintf('F\t255\t255\t255');
    tmp{pal_len+4} = sprintf('B\t0\t0\t0');
	if (handles_mirone.have_nans)
        cor = round(cor_nan*255);
        cor_str = sprintf(['N\t' num2str(cor(1),'%.12g') '\t' num2str(cor(2),'%.12g') '\t' num2str(cor(3),'%.12g')]);
        tmp{pal_len+5} = sprintf(cor_str);
	else
        tmp{pal_len+5} = sprintf('N\t255\t255\t255');
	end
    sc_cpt = [dest_dir filesep prefix '.cpt'];
	fid = fopen(sc_cpt,'wt');
	for (i=1:pal_len+5)    fprintf(fid,'%s\n',tmp{i});     end
	fclose(fid);
    clear tmp z_min z_max pal_len pal cor cor_str fid Z dz z1 z2
else        % Remove the cpt declaration. After all we won't go to use it
    script(id_cpt) = [];    l=l-1;
end

% ------------- Coastlines section -----------------------------------
if (have_psc)       % We have pscoast commands
    script{l} = [' '];                      l=l+1;
    script{l} = [comm 'Plot coastlines'];   l=l+1;
    script{l} = ['pscoast ' opt_psc ellips ' -R -J -O -K >> ' pb 'ps' pf];    l=l+1;
end

% ------------- Search for contour lines ----------------------------------------------------
ALLtextHand = findobj(get(handles_mirone.axes1,'Child'),'Type','text');
% % % If we have focal mecanisms with labels, remove their handles right away
% % h = findobj(ALLtextHand,'Tag','TextMeca');                                  % I'M NOT SURE ON THIS ONE
% % if (~isempty(h))    ALLtextHand = setxor(ALLtextHand, h);   end

tag = get(ALLlineHand,'Tag');
if (~isempty(tag) & ~isempty(handles_mirone.grdname))
    h = findobj(ALLlineHand,'Tag','contour');
    if (~isempty(h))
        h_label = findobj(ALLtextHand,'Tag','contour');   % Search for contour labels
        if (~isempty(h_label))
            lab = get(h_label,'UserData');
            if (iscell(lab))    lab = unique(cat(1,lab{:}));    end
        else
            lab = [];
        end
        conts = zeros(length(h),1);
        for i=1:length(h)
            conts(i) = getappdata(h(i),'cont_label');
        end
        conts = unique(conts);
        no_anot = setxor(conts,lab);    % Contour levels that are not annotated
        name = [prefix_ddir '_cont.dat'];
        fid = fopen(name,'wt');
        if (isempty(no_anot))           % Annotate all contours
            fprintf(fid,'%.5f\tA\n',conts);
        else                            % Annotate only some contours
            conts = [[lab; no_anot] [ones(length(lab),1)*double('A'); ones(length(no_anot),1)*double('C')]];
            conts = sortrows(conts);
            fprintf(fid,'%.5f\t%c\n',conts');
        end
        fclose(fid);
        script{l} = [' '];                                  l=l+1;
        script{l} = [comm ' ---- Plot contours'];           l=l+1;
        script{l} = ['grdcontour ' pb 'grd' pf ' -R -J -C' [prefix '_cont.dat'] ellips ' -O -K >> ' pb 'ps' pf];
        l=l+1;
        used_grd = 1;
        ALLlineHand = setxor(ALLlineHand, h);       % h is processed, so remove it from handles list
        ALLtextHand = setxor(ALLtextHand, h_label); % same for contour label strings
        clear h conts fid name no_anot lab;
    end
end

% ------------- Search for symbols -----------------------------------------------------------    
tag = get(ALLlineHand,'Tag');
if (~isempty(tag))
    h = findobj(ALLlineHand,'Tag','Symbol');
    h = [h; findobj(ALLlineHand,'Tag','City_major')];
    h = [h; findobj(ALLlineHand,'Tag','City_other')];
    h = [h; findobj(ALLlineHand,'Tag','volcano')];
    h = [h; findobj(ALLlineHand,'Tag','hotspot')];
    h = [h; findobj(ALLlineHand,'Tag','Earthquakes')];
    h = [h; findobj(ALLlineHand,'Tag','DSDP')];
    h = [h; findobj(ALLlineHand,'Tag','ODP')];
    
    % Search for points as Markers (that is, line with no line - just symbols on vertices)
	h_shit = get(ALLlineHand,'LineStyle');
	h_num_shit = strmatch('none',h_shit);
	if (h_num_shit)
        h = [h; ALLlineHand(h_num_shit)];
	end
    clear h_shit h_num_shit;
    
    if (~isempty(h))
        symbols = get_symbols(h);
        haveSymbol = 1;
        ALLlineHand = setxor(ALLlineHand, h);    % h is processed, so remove it from handles list
    end
    clear h;
end

if (haveSymbol)
    ns = length(symbols.x);
    name = [prefix_ddir '_symb.dat'];    name_sc = [prefix '_symb.dat'];
    if (ns > 1 & length(symbols.Size) == 1)      % We have the same symbol repeated ns times
    	fid = fopen(name,'wt');
        cor_fill = round(symbols.FillColor{1} * 255);
        cor_fill = [num2str(cor_fill(1)) '/' num2str(cor_fill(2)) '/' num2str(cor_fill(3))];
        cor_edge = round(symbols.EdgeColor{1} * 255);
        cor_edge = [num2str(cor_edge(1)) '/' num2str(cor_edge(2)) '/' num2str(cor_edge(3))];
        fprintf(fid,'%.5f\t%.5f\n',[symbols.x{:}; symbols.y{:}]);
		script{l} = [' '];                                  l=l+1;
		script{l} = [comm ' ---- Plot symbols'];            l=l+1;
		script{l} = ['psxy ' name_sc ' -S' symbols.Marker num2str(symbols.Size{1}) 'p' ' -G' cor_fill ...
                ' -W1/' cor_edge ellips ' -R -J -O -K >> ' pb 'ps' pf];    l=l+1;
    	fclose(fid);
    elseif (ns == 1 & length(symbols.Size) == 1)      % We have only one symbol
        cor_fill = round(symbols.FillColor{1} * 255);
        cor_fill = [num2str(cor_fill(1)) '/' num2str(cor_fill(2)) '/' num2str(cor_fill(3))];
        cor_edge = round(symbols.EdgeColor{1} * 255);
        cor_edge = [num2str(cor_edge(1)) '/' num2str(cor_edge(2)) '/' num2str(cor_edge(3))];
		script{l} = [' '];                                  l=l+1;
		script{l} = [comm ' ---- Plot symbol'];             l=l+1;
		script{l} = ['echo ' num2str(symbols.x{1},'%.5f') ' ' num2str(symbols.y{1},'%.5f') ' | ' ...
                    'psxy -S' symbols.Marker num2str(symbols.Size{1}) 'p' ' -G' cor_fill ...
                    ' -W1/' cor_edge ellips ' -R -J -O -K >> ' pb 'ps' pf];    l=l+1;        
    else                                % We have ns different symbols
        m = zeros(ns,1);
        for i=1:ns
            m(i) = size(symbols.x{i},2);
        end
        n = find(m ~= 1);
        if (~isempty(n))                % We have a mixed scenario. Individual as well as group symbols
            script = write_group_symb(prefix,prefix_ddir,comm,pb,pf,ellips,symbols,n,script);
            symbols.x(n) = [];          symbols.y(n) = [];  % Clear processed symbols
            symbols.FillColor(n) = [];  symbols.EdgeColor(n) = [];
            symbols.Size(n) = [];       symbols.Marker(n,:) = [];
            l = length(script) + 1;
            ns = ns - length(n);
        end
        clear m n;
    	fid = fopen(name,'wt');
        for i=1:ns
            cor_fill = round(symbols.FillColor{i} * 255);
            cor_fill = [num2str(cor_fill(1)) '/' num2str(cor_fill(2)) '/' num2str(cor_fill(3))];
            cor_edge = round(symbols.EdgeColor{i} * 255);
            cor_edge = [num2str(cor_edge(1)) '/' num2str(cor_edge(2)) '/' num2str(cor_edge(3))];
            fprintf(fid,'%s\n',['>' ' -G' cor_fill ' -W1/' cor_edge]);
            fprintf(fid,'%.5f\t%.5f\t%.0f\t%s\n',symbols.x{i},symbols.y{i},symbols.Size{i},symbols.Marker(i,:));
        end
		script{l} = [' '];                                  l=l+1;
		script{l} = [comm ' ---- Plot symbols'];            l=l+1;
		script{l} = ['psxy ' name_sc ellips ' -S -R -J --MEASURE_UNIT=point -M -O -K >> ' pb 'ps' pf];    l=l+1;
    	fclose(fid);
    end
    clear cor_fill cor_edge ns symbols haveSymbol name name_sc;
end
% ------------------------------------------------------------------------------------------------

% ------------- Search for focal mecanisms ----------------------------
ALLpatchHand = findobj(get(handles_mirone.axes1,'Child'),'Type','patch');
if (~isempty(ALLpatchHand))
    focHand = findobj(ALLpatchHand,'Tag','FocalMeca');
    if (~isempty(focHand))
        % First deal with the 'line anchors'
        focHandAnchor = findobj(ALLlineHand,'Tag','FocalMecaAnchor');   % Handles of the line anchors
        x = get(focHandAnchor,'XData');         y = get(focHandAnchor,'YData');
        if (iscell(x))      x = cell2mat(x);    y = cell2mat(y);    end
        id_anch = find(diff(x,1,2));
        
        psmeca_line = cell(length(focHand),1);
        for(k=1:length(focHand))
            psmeca_line{k} = getappdata(focHand(k),'psmeca_com');
        end
        psmeca_line = cat(1,psmeca_line{:});    % This also get us rid of empty cell fields.
                
        n_cols = size(psmeca_line,2);
        if (n_cols == 10 | n_cols == 14)
            with_label = 1;
        else
            with_label = 0;
        end
        name = [prefix_ddir '_meca.dat'];   name_sc = [prefix '_meca.dat'];     opt_C = '';
        fid = fopen(name,'wt');
        if (n_cols == 9 | n_cols == 10)         % Aki & Richard convention
            % If beach-bals are not ploted at their origin update the ploting coords columns
            if (~isempty(id_anch))
                psmeca_line(:,8) = x(:,2);  psmeca_line(:,9) = y(:,2);     opt_C = '-C';
            end
            opt_S = ['-Sa' getappdata(handles_mirone.figure1,'MecaMag5') 'c'];
            format = ['%.4f\t%.4f\t%.1f\t%.0f\t%.0f\t%.0f\t%.1f\t%.4f\t%.4f'];
            for (k=1:size(psmeca_line,1))
                fprintf(fid,format,[psmeca_line(k,1:9)]);
                if (with_label)
                    fprintf(fid,'\t%s\n',num2str(psmeca_line(k,10)));
                else
                    fprintf(fid,'\n');
                end
            end
        elseif (n_cols == 13 | n_cols == 14)    % CMT convention
            % If beach-bals are not ploted at their origin update the ploting coords columns
            if (~isempty(id_anch))
                psmeca_line(:,12) = x(:,2); psmeca_line(:,13) = y(:,2);     opt_C = '-C';
            end
            psmeca_line(:,11) = psmeca_line(:,11) + 7;      % psmeca uses Moment in Dyn-cm
            opt_S = ['-Sc' getappdata(handles_mirone.figure1,'MecaMag5') 'c'];
            format = ['%.4f\t%.4f\t%.1f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.2f\t%d\t%.4f\t%.4f'];
            for (k=1:size(psmeca_line,1))
                fprintf(fid,format,[psmeca_line(k,1:13)]);
                if (with_label)
                    fprintf(fid,'\t%s\n',num2str(psmeca_line(k,14)));
                else
                    fprintf(fid,'\n');
                end
            end
        end
        fclose(fid);
    	script{l} = [' '];              l=l+1;
        script{l} = [comm ' ---- Plot Focal Mechanisms'];   l=l+1;
        script{l} = ['psmeca ' opt_S ' ' opt_C ' ' name_sc ellips ' -R -J -O -K >> ' pb 'ps' pf];    l=l+1;
        ALLpatchHand = setxor(ALLpatchHand, focHand);       % focHand is processed, so remove it from handles list
        ALLlineHand  = setxor(ALLlineHand, focHandAnchor);  %       iden
        clear focHand name name_sc psmeca_line with_label n_cols id_anch opt_S opt_C
    end
end
% -------------------------------------------------------------------------------------------------------

% ------------- Search for countries ----------------------------
if (~isempty(ALLpatchHand))
    AtlasHand = findobj(ALLpatchHand,'Tag','Atlas');
    if (~isempty(AtlasHand))
        used_countries = 1;     need_path = 1;
        n_cts = length(AtlasHand);
        if (n_cts > 1)                  % We have multiple countries
            ct_names = cell(n_cts,1);   % To hold the country names
            for (k=1:n_cts)             % Loop over all countries found and store theyr names
                ct_names{k} = get(AtlasHand(k),'UserData');
            end
        else                            % We have only one country
            ct_names = {get(AtlasHand(k),'UserData')};
        end
        ct_names = unique(ct_names);    % Many countries have several polygons (e.g. islands).
        name = [prefix_ddir '_country_names.txt'];   name_sc = [prefix '_country_names.txt'];
        fid = fopen(name,'wt');
        fprintf(fid,'%s\n',ct_names{:});        fclose(fid);
		script{l} = [' '];              l=l+1;
        script{l} = [comm ' ---- Plot countries. NOTE: THIS IS NOT A GMT PROGRAM'];   l=l+1;
        ct_with_pato = getappdata(handles_mirone.figure1,'AtlasResolution');
		script{l} = [pwd filesep 'country_extract -P' name ' ' ct_with_pato ' -C | ',...
                'psxy ' ellips ' -R -J -W0.5p -M -O -K >> ' pb 'ps' pf];    l=l+1;
        ALLpatchHand = setxor(ALLpatchHand, AtlasHand);       % AtlasHand is processed, so remove it from handles list
        clear AtlasHand name name_sc n_cts ct_with_pato
    end
end
% -------------------------------------------------------------------------------------------------------

% ------------- Search for "telhas" ---------------------------------
if (~isempty(ALLpatchHand))
    TelhasHand = findobj(ALLpatchHand,'Tag','tapete');
    if (~isempty(TelhasHand))
        tmp = findobj(ALLpatchHand,'Tag','tapete_R');
        ALLpatchHand = setxor(ALLpatchHand, tmp);       % Remove the reverse "telhas" 
        n_tapetes = length(TelhasHand);
        for (i=1:n_tapetes)
            saved = get(TelhasHand(i),'UserData');
            name = [prefix_ddir '_telha_' num2str(i) '.dat'];
            name_sc = [prefix '_telha_' num2str(i) '.dat'];
            if (~isempty(saved))
                fid = fopen(name,'wt');
                fprintf(fid,'%.5f\t%.5f\n',saved.line');
                fclose(fid);
			    script{l} = [' '];              l=l+1;
                script{l} = [comm ' ---- Plot telhas. NOTE: THIS IS NOT A GMT PROGRAM'];   l=l+1;
			    script{l} = ['telha ' name_sc ' ' saved.opt_E ' ' saved.opt_I ' ',...
                    saved.opt_N ' ' saved.opt_T ' -Blixo.dat'];     l=l+1;
                script{l} = ['psxy lixo.dat ' ellips ' -R -J -M -L -O -K >> ' pb 'ps' pf];    l=l+1;
                fclose(fid);
                ALLpatchHand = setxor(ALLpatchHand, TelhasHand);       % TelhasHand is processed, so remove it from handles list
            end
        end
        %ALLpatchHand = setxor(ALLpatchHand, TelhasHand);       % TelhasHand is processed, so remove it from handles list
        clear TelhasHand name name_sc n_tapetes tmp 
    end
end

% -------------------------------------------------------------------------------------------------------

% ------------- Search for closed polygons ----------------------------
if (~isempty(ALLpatchHand))
    xx = get(ALLpatchHand,'XData');     yy = get(ALLpatchHand,'YData');
    n_patch = length(ALLpatchHand);
    LineStyle = get(ALLpatchHand,'LineStyle');
    LineWidth = get(ALLpatchHand,'LineWidth');
    if (iscell(LineWidth))      LineWidth = cat(1,LineWidth{:});     end
    EdgeColor = get(ALLpatchHand,'EdgeColor');
    if (iscell(EdgeColor))      EdgeColor = cat(1,EdgeColor{:});     end
    FillColor = get(ALLpatchHand,'FaceColor');
    if (iscell(FillColor))
        resp = strmatch('none',char(FillColor{:}));
        if (isempty(resp))
            FillColor = cat(1,FillColor{:});
        else
            for (i=1:length(resp))                  % Signal down that this is a non colored polygon
                FillColor{resp(i)} = [-1 -1 -1];    % FDS it worked. I wonder why FillColor{resp} = repmat([-1 -1 -1],length(resp),1); DOESN'T
            end
            FillColor = cat(1,FillColor{:});
        end
    else                % We have only one patch
        xx = num2cell(xx,1);   yy = num2cell(yy,1);   % Make it a cell for reducing the head-hakes
        resp = strmatch('none',FillColor);
        if (~isempty(resp))
            FillColor = [-1 -1 -1];                 % Signal down that this is a non colored polygon
        end
    end
    name = [prefix_ddir '_patch.dat'];    name_sc = [prefix '_patch.dat'];
	fid = fopen(name,'wt');
    for (i=1:n_patch)
        cor_edge = round(EdgeColor(i,1:3) * 255);
        cor_edge = [num2str(cor_edge(1)) '/' num2str(cor_edge(2)) '/' num2str(cor_edge(3))];
        cor_fill = round(FillColor(i,1:3) * 255);
        if (cor_fill(1) >= 0)       % Color filled polygon
            cor_fill = [num2str(cor_fill(1)) '/' num2str(cor_fill(2)) '/' num2str(cor_fill(3))];
            mlt_comm = ['>' ' -G' cor_fill ' -W' num2str(LineWidth(i)) 'p/' cor_edge];
        else                        % No filling color
            mlt_comm = ['> -W' num2str(LineWidth(i)) 'p/' cor_edge];
        end
        
        if (any(isnan(xx{i})))      % If we have NaNs we need to split into segments
            [latcells,loncells] = polysplit(yy{i}(:),xx{i}(:));
            for (j=1:numel(loncells))
                fprintf(fid,'%s\n',mlt_comm);
                fprintf(fid,'%.5f\t%.5f\n',[loncells{j}(:)'; latcells{j}(:)']);
            end
        else
            fprintf(fid,'%s\n',mlt_comm);
            fprintf(fid,'%.5f\t%.5f\n',[xx{i}(:)'; yy{i}(:)']);
        end
    end
	fclose(fid);
	script{l} = [' '];              l=l+1;
    script{l} = [comm ' ---- Plot closed AND colored polygons'];   l=l+1;
	script{l} = ['psxy ' name_sc ellips ' -R -J --MEASURE_UNIT=point -M -O -K >> ' pb 'ps' pf];    l=l+1;
    clear ALLpatchHand name name_sc n_patch xx yy LineStyle LineWidth EdgeColor FillColor cor_edge resp
end

% ------------- Search for lines or polylines ----------------------------
if (~isempty(ALLlineHand))      % OK, now the only left line handles must be, plines, mb-tracks, etc
    xx = get(ALLlineHand,'XData');     yy = get(ALLlineHand,'YData');
    if (~iscell(xx))            % We have only one line
        xx = num2cell(xx(:),1);   yy = num2cell(yy(:),1);
    end
    n_lin = length(xx);
	script{l} = [' '];                      l=l+1;
    script{l} = [comm ' ---- Plot lines'];  l=l+1;
    if (n_lin > 0)     % We have more than one line.         E SENAO?
        LineStyle = get(ALLlineHand,'LineStyle');
        [LineStyle,LineStyle_gmt] = lineStyle2num(LineStyle);
        LineWidth = get(ALLlineHand,'LineWidth');
        if (iscell(LineWidth)),     LineWidth = cat(1,LineWidth{:});    end
        LineColor = get(ALLlineHand,'Color');
        if (iscell(LineColor)),     LineColor = cat(1,LineColor{:});    end
        [b,m] = sortrows([LineWidth LineColor LineStyle]);
        m = m(end:-1:1);            % Revert order because I want thicker lines ploted first
        xx = xx(m);     yy = yy(m);
        LineWidth = LineWidth(m,:);     LineColor = LineColor(m,:);
        LineStyle = LineStyle(m);       LineStyle_gmt = LineStyle_gmt(m,:);
        [b,m] = unique([LineWidth LineColor LineStyle],'rows');   % reuse b,m
        m = m(end:-1:1);            % OK, now we have to put it back in ascending order        
        reps = setxor(m,1:n_lin);   % Find repeated (e.g. same color & line thickness)
        m = [0; m];                 % use this first index to help file creation algo
        for i=1:length(m)-1
            name = [prefix_ddir '_line_' num2str(i) '.dat'];
            name_sc = [prefix '_line_' num2str(i) '.dat'];
            fid = fopen(name,'wt');
            for j=m(i)+1:m(i+1)
                if (any(isnan(xx{j})))          % If we have NaNs we need to split into segments
                    [latcells,loncells] = polysplit(yy{j}(:),xx{j}(:));
                    for (k=1:numel(loncells))
                        fprintf(fid,'%s\n','>');
                        fprintf(fid,'%.5f\t%.5f\n',[loncells{k}(:)'; latcells{k}(:)']);
                    end
                else
                    fprintf(fid,'%s\n','>');
                    fprintf(fid,'%.5f\t%.5f\n',[xx{j}(:) yy{j}(:)]');
                end
            end
            fclose(fid);
            cor = round(LineColor(j,:) * 255);
            cor = [num2str(cor(1)) '/' num2str(cor(2)) '/' num2str(cor(3))];
            script{l} = ['psxy ' name_sc ellips ' -R -J -W' num2str(LineWidth(j)) 'p,' ...
                    cor LineStyle_gmt{j} ' --MEASURE_UNIT=point -M -O -K >> ' pb 'ps' pf];
            l=l+1;
        end
    end
    clear xx yy cor fid m name name_sc reps LineStyle LineWidth LineColor
end

% ------------- Search for text strings ---------------------------------------------
if (~isempty(ALLtextHand))          % ALLtextHand was found above in the search for contours
	if (~isempty(ALLtextHand))      % We (still) have text fields
        pos = get(ALLtextHand,'Position');      font = get(ALLtextHand,'FontName');
        fsize = get(ALLtextHand,'FontSize');    fcolor = get(ALLtextHand,'Color');
        if (isnumeric(fcolor))
            fcolor = round(fcolor * 255);
            if (numel(fcolor) == 1)
                opt_G = {[' -G' num2str(fcolor)]};
            else
                opt_G = {[' -G' num2str(fcolor(1)) '/' num2str(fcolor(2)) '/' num2str(fcolor(3))]};
            end
        elseif (ischar(fcolor))     % Shit, we have to decode the color letter
            switch fcolor
                case 'w',       opt_G = {' -G255'};
                case 'k',       opt_G = {' -G0'};
                case 'y',       opt_G = {' -G255/255/0'};
                case 'c',       opt_G = {' -G0/255/255'};
                case 'r',       opt_G = {' -G255/0/0'};
                case 'g',       opt_G = {' -G0/255/0'};
                case 'b',       opt_G = {' -G0/0/255'};
                otherwise,      opt_G = {''};
            end
        elseif (iscell(fcolor))     % Double shit, we have to convert a Mx3 cell matrix into texts
            tmp = cell2mat(fcolor) * 255;
            opt_G = cell(size(tmp,1),1);
            for (m = 1:size(tmp,1))
                opt_G{m} = [' -G' num2str(tmp(1)) '/' num2str(tmp(2)) '/' num2str(tmp(3))];
            end
        else
            opt_G = {''};
        end
        str = get(ALLtextHand,'String');      angle = get(ALLtextHand,'Rotation');
        if (~iscell(pos))           % Make them cells for author's mental sanity
            pos = num2cell(pos(:),1);   fsize = num2cell(fsize,1);      angle = num2cell(angle,1);
            str = {str};                font = {font};
        end
        n_text = length(str);
		script{l} = [' '];                     l=l+1;
        script{l} = [comm ' ---- Plot text strings'];   l=l+1;    
        for (i=1:n_text)
            script{l} = ['echo ' num2str(pos{i}(1),'%.5f') ' ' num2str(pos{i}(2),'%.5f') ' ' num2str(fsize{i}) ' ' ...
                num2str(angle{i}) ' 4 LB ' str{i} ' | pstext' ellips opt_G{i} ' -R -J -O -K >> ' pb 'ps' pf];
            l=l+1;
        end
	end
end

%if (~isempty(ALLlineHand))     % Perdemos linhas

% -------------- Write the script ---------------------------------------------
% First do some eventual cleaning
if (~isempty(handles_mirone.grdname) & ~used_grd)          script(id_grd) = [];        end;
if (ispc & (used_grd | used_countries) & need_path & ~strcmp(sc,'bat'))
    tmp = cell(7,1);
    tmp{1} = [comm 'If you see this message is because you choosed to generate a c-shell'];
    tmp{2} = [comm 'script in a windows running machine. Notice that I have no means to'];
    tmp{3} = [comm 'guess on what imullation schema (e.g. SFU, cygwin, etc..) you intend'];
    tmp{4} = [comm 'to run this script. The point is that they use different file names'];
    tmp{5} = [comm 'mapping. While cygwin accepts the c:\somewhere\somefile, SFU wants'];
    tmp{6} = [comm '/dev/fs/C/somewhere/somefile. So it''s your responsability to set'];
    if (used_grd & ~used_countries)
        tmp{7} = [comm 'the $grd variable with the correct path.'];
    elseif (~used_grd & used_countries)
        tmp{7} = [comm 'the paths correctly in the "country_select" command line.'];
    else            % Both cases
        tmp{7} = [comm 'the $grd variable with the correct path. And the same'];
        tmp{9} = [comm 'for the paths in the "country_select" command line.'];
    end
    tmp{end+1} = '';
    script = [script(1:4); tmp; script(5:end)];
    out_msg = 1;
end
fid = fopen([prefix_ddir '_mir.' sc],'wt');
for i = 1:length(script)-1
    fprintf(fid,'%s\n',script{i});
end
if (strcmp(sc,'bat'))   cut = 11;
else                    cut = 10;    end
last = [script{i+1}(1:end-cut) ' >> ' pb 'ps' pf];    % Remove the last '-K'
fprintf(fid,'%s\n',last);
fclose(fid);

% ----------------------------------------------------------------------------------
function symbol = get_symbols(hand)
xx = get(hand,'XData');     yy = get(hand,'YData');
if (~iscell(xx))
    xx = num2cell(xx,1);   yy = num2cell(yy,1);   % Make it a cell for reducing the head-hakes
end
symbol.x = xx(:);       symbol.y = yy(:);
symbol.Marker = get(hand,'Marker');
zz = get(hand,'MarkerSize');
if (~iscell(zz))    symbol.Size = num2cell(zz,1);
else                symbol.Size = zz;       end
zz = get(hand,'MarkerFaceColor');
if (~iscell(zz))    symbol.FillColor = num2cell(zz(:),1);
else                symbol.FillColor = zz;  end
zz = get(hand,'MarkerEdgeColor');
if (~iscell(zz))    symbol.EdgeColor = num2cell(zz(:),1);
else                symbol.EdgeColor = zz;  end

symbol.Marker = char(symbol.Marker);
symbol.Marker = symbol.Marker(:,1);

symbol.Marker(symbol.Marker == '^') = 't';
symbol.Marker(symbol.Marker == '>') = 't';      % not in GMT
symbol.Marker(symbol.Marker == '<') = 't';      % not in GMT
symbol.Marker(symbol.Marker == 'v') = 'i';
symbol.Marker(symbol.Marker == '.') = 'p';
symbol.Marker(symbol.Marker == 'd') = 'd';
symbol.Marker(symbol.Marker == 'o') = 'c';
symbol.Marker(symbol.Marker == '+') = 'x';      % not in GMT
symbol.Marker(symbol.Marker == 'x') = 'x';
symbol.Marker(symbol.Marker == 's') = 's';
symbol.Marker(symbol.Marker == '*') = 'a';
symbol.Marker(symbol.Marker == 'p') = 'a';      % not in GMT
symbol.Marker(symbol.Marker == 'h') = 'a';      % not in GMT

% ----------------------------------------------------------------------------------
function [LineStyle_num,LineStyle_gmt] = lineStyle2num(LineStyle)
if (~iscell(LineStyle))     LineStyle = {LineStyle};    end
lt = {'-'; '--'; ':'; '-.'};
LineStyle_num = strrep(LineStyle,lt{4},'4');
LineStyle_num = strrep(LineStyle_num,lt{3},'3');
LineStyle_num = strrep(LineStyle_num,lt{2},'2');
LineStyle_num = strrep(LineStyle_num,lt{1},'1');
tmp = LineStyle_num;
LineStyle_num = str2num(cat(1,LineStyle_num{:}));
% Convert to GMT linestyles
tmp = strrep(tmp,'4',',.-');
tmp = strrep(tmp,'3',',.');
tmp = strrep(tmp,'2',',-');
LineStyle_gmt = strrep(tmp,'1','');

% --------------------------------------------------------------------
function script = write_group_symb(prefix,prefix_ddir,comm,pb,pf,ellips,symbols,n,script)
% Write a group symbol to file, and uppdate the "script"
l = length(script) + 1;
for i=1:length(n)
    name = [prefix_ddir '_symb_' num2str(i) '.dat'];
    name_sc = [prefix '_symb_' num2str(i) '.dat'];
	fid = fopen(name,'wt');
	cor_fill = round(symbols.FillColor{n(i)} * 255);
	cor_fill = [num2str(cor_fill(1)) '/' num2str(cor_fill(2)) '/' num2str(cor_fill(3))];
	cor_edge = round(symbols.EdgeColor{n(i)} * 255);
	cor_edge = [num2str(cor_edge(1)) '/' num2str(cor_edge(2)) '/' num2str(cor_edge(3))];
	fprintf(fid,'%.5f\t%.5f\n',[symbols.x{n(i)}; symbols.y{n(i)}]);
	script{l} = [' '];                                  l=l+1;
	script{l} = [comm 'Plot symbols'];                  l=l+1;
	script{l} = ['psxy ' name_sc ' -S' symbols.Marker(n(i)) num2str(symbols.Size{n(i)}) 'p' ' -G' cor_fill ...
            ' -W1/' cor_edge ellips ' -R -J -O -K >> ' pb 'ps' pf];    l=l+1;
	fclose(fid);
end

% --------------------------------------------------------------------------------
function [latcells,loncells] = polysplit(lat,lon)
%POLYSPLIT Extract segments of NaN-delimited polygon vectors to cell arrays
%
%   [LATCELLS,LONCELLS] = POLYSPLIT(LAT,LON) returns the NaN-delimited
%   segments of the vectors LAT and LON as N-by-1 cell arrays with one
%   polygon segment per cell.  LAT and LON must be the same size and have
%   identically-placed NaNs.  The polygon segments are column vectors if
%   LAT and LON are column vectors, and row vectors otherwise.

% Copyright 1996-2006 The MathWorks, Inc.
% $Revision: 1.4.4.5 $    $Date: 2006/05/24 03:35:26 $

% checkinput(lat,{'numeric'},{'real','vector'},mfilename,'LAT',1);
% checkinput(lon,{'numeric'},{'real','vector'},mfilename,'LON',2);
[lat, lon] = removeExtraNanSeparators(lat, lon);

% Find NaN locations.
indx = find(isnan(lat(:)));

% Simulate the trailing NaN if it's missing.
if ~isempty(lat) && ~isnan(lat(end))
    indx(end+1,1) = numel(lat) + 1;
end

%  Extract each segment into pre-allocated N-by-1 cell arrays, where N is
%  the number of polygon segments.  (Add a leading zero to the indx array
%  to make indexing work for the first segment.)
N = numel(indx);
latcells = cell(N,1);
loncells = cell(N,1);
indx = [0; indx];
for k = 1:N
    iStart = indx(k)   + 1;
    iEnd   = indx(k+1) - 1;
    latcells{k} = lat(iStart:iEnd);
    loncells{k} = lon(iStart:iEnd);
end

% --------------------------------------------------------------------------------
function [xdata, ydata, zdata] = removeExtraNanSeparators(xdata, ydata, zdata)
%removeExtraNanSeparators  Clean up NaN separators in polygons and lines
%
%   [XDATA, YDATA] = removeExtraNanSeparators(XDATA, YDATA) removes NaNs
%   from the vectors XDATA and YDATA, leaving only isolated NaN separators.
%   If present, one or more leading NaNs are removed entirely.  If present,
%   a single trailing NaN is preserved.  NaNs are removed, but never added,
%   so if the input lacks a trailing NaN, so will the output.  XDATA and
%   YDATA must match in size and have identical NaN locations.
%
%   [XDATA, YDATA, ZDATA] = removeExtraNanSeparators(XDATA, YDATA, ZDATA)
%   removes NaNs from the vectors XDATA, YDATA, and ZDATA, leaving only
%   isolated NaN separators and optionally, consistent with the input, a
%   single trailing NaN.
%
%   Examples
%   --------
%   xin = [NaN NaN 1:3 NaN 4:5 NaN NaN NaN 6:9 NaN NaN]
%   yin = xin;
%   [xout, yout] = removeExtraNanSeparators(xin, yin);
%   xout
%
%   xin = [NaN 1:3 NaN NaN 4:5 NaN NaN NaN 6:9]'
%   yin = xin;
%   zin = xin;
%   [xout, yout, zout] = removeExtraNanSeparators(xin, yin, zin);
%   xout

% Copyright 2005-2006 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2006/06/15 20:11:13 $

% if nargin < 3
%     if ~isequal(isnan(xdata), isnan(ydata))
%         eid = sprintf('%s:%s:inconsistentXY', getcomp, mfilename);
%         error(eid,'XDATA and YDATA mismatch in size or NaN locations.')
%     end
% else
%     if ~isequal(isnan(xdata), isnan(ydata), isnan(zdata))
%         eid = sprintf('%s:%s:inconsistentXYZ', getcomp, mfilename);
%         error(eid,'XDATA, YDATA (or ZDATA) mismatch in size or NaN locations.')
%     end
% end

p = find(isnan(xdata(:)'));     % Determing the positions of each NaN.

% Determine the position of each NaN that is not the final element in a sequence of contiguous NaNs.
q = p(diff(p) == 1);

% If there's a leading sequence of NaNs (a sequence starting with a NaN in
% position 1), determine the position of each NaN in this sequence.
if isempty(p),      r = [];
else                r = find((p - (1:numel(p))) == 0);
end

% Determine the position of each excess NaN.
if isempty(r),      s = q;
else                s = [r q(q > r(end))];
end

% Remove the excess NaNs.
xdata(s) = [];      ydata(s) = [];
if (nargin >= 3),   zdata(s) = [];  end

% ---------------------- Creates and returns a handle to the GUI figure. 
function write_gmt_script_LayoutFcn(h1,handles);

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','write_gmt_script',...
'NumberTitle','off',...
'Position',[520 340 561 460],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

h2 = uicontrol('Parent',h1,'Position',[30 9 205 121],'Style','frame','Tag','frame1');

h3 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@write_gmt_script_uicallback,h1,'popup_PaperSize_Callback'},...
'Position',[40 98 181 22],...
'String','A4 595 842',...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_PaperSize');

h4 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'radiobutton_P_Callback'},...
'Position',[40 42 71 15],...
'String','Portrait',...
'Style','radiobutton',...
'Value',1,...
'Tag','radiobutton_P');

h5 = axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Position',[30 169 281 281],...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','axes1');

h10 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'radiobutton_L_Callback'},...
'Position',[40 18 72 15],...
'String','Landscape',...
'Style','radiobutton',...
'Tag','radiobutton_L');

h11 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'radiobutton_setWidth_Callback'},...
'Position',[135 42 95 15],...
'String','Set map width',...
'Style','radiobutton',...
'TooltipString','Check this to force Map width = rectangle width',...
'Value',1,...
'Tag','radiobutton_setWidth');

h12 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'radiobutton_setHeight_Callback'},...
'Position',[135 18 95 15],...
'String','Set map height',...
'Style','radiobutton',...
'TooltipString','Check this to force Map height = rectangle height',...
'Tag','radiobutton_setHeight');

uicontrol('Parent',h1, ...
'Callback',{@write_gmt_script_uicallback,h1,'radiobutton_180_180_Callback'},...
'Position',[340 150 95 15],...
'String','[-180 180]', 'Style','radiobutton',...
'TooltipString','Plot longitudes in the [-180;180] range',...
'Value',1,...
'Tag','radiobutton_180_180');

uicontrol('Parent',h1, ...
'Callback',{@write_gmt_script_uicallback,h1,'radiobutton_0_360_Callback'},...
'Position',[340 130 95 15],...
'String','[0 360]', 'Style','radiobutton',...
'TooltipString','Plot longitudes in the [0;360] range',...
'Tag','radiobutton_0_360');

h13 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'pushbutton_uppdate_Callback'},...
'Position',[340 60 68 33],...
'String','Update', 'FontWeight','bold',...
'TooltipString','Update the rectangle size to allowed dimensions according to projection and side constraints',...
'Tag','pushbutton_uppdate');

h14 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@write_gmt_script_uicallback,h1,'edit_X0_Callback'},...
'Position',[480 89 61 21],...
'String','2.5',...
'Style','edit',...
'TooltipString','Plot X origin',...
'Tag','edit_X0');

h15 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@write_gmt_script_uicallback,h1,'edit_Y0_Callback'},...
'Position',[480 59 61 21],...
'String','2.5',...
'Style','edit',...
'TooltipString','Plot Y origin',...
'Tag','edit_Y0');

h16 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@write_gmt_script_uicallback,h1,'edit_mapWidth_Callback'},...
'Position',[480 149 60 21],...
'Style','edit',...
'TooltipString','Map width',...
'Tag','edit_mapWidth');

h17 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'pushbutton_mapProjections_Callback'},...
'Position',[400 367 95 23],...
'String','Map projection',...
'TooltipString','Select/Change the map projection',...
'Tag','pushbutton_mapProjections');

h18 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'radiobutton_cm_Callback'},...
'Position',[40 75 40 15],...
'String','cm',...
'Style','radiobutton',...
'TooltipString','Show paper size in centimeters',...
'Value',1,...
'Tag','radiobutton_cm');

h19 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'radiobutton_in_Callback'},...
'Position',[116 75 40 15],...
'String','in',...
'Style','radiobutton',...
'TooltipString','Show paper size in inches',...
'Tag','radiobutton_in');

h20 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'radiobutton_pt_Callback'},...
'Position',[179 75 40 15],...
'String','pt',...
'Style','radiobutton',...
'TooltipString','Show paper size in points',...
'Tag','radiobutton_pt');

h21 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'pushbutton_OK_Callback'},...
'FontSize',9,...
'FontWeight','bold',...
'Position',[460 8 82 24],...
'String','Write script',...
'Tag','pushbutton_OK');

h22 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@write_gmt_script_uicallback,h1,'edit_mapHeight_Callback'},...
'Position',[480 119 60 21],...
'Style','edit',...
'TooltipString','Map height',...
'Tag','edit_mapHeight');

h23 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@write_gmt_script_uicallback,h1,'edit_scale_Callback'},...
'Position',[410 179 131 21],...
'Style','edit',...
'TooltipString','Aproximate map scale',...
'Tag','edit_scale');

h24 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@write_gmt_script_uicallback,h1,'popup_directory_list_Callback'},...
'Position',[330 240 191 22],...
'Style','popupmenu',...
'TooltipString','Save script and files in this directory',...
'Value',1,...
'Tag','popup_directory_list');

h25 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@write_gmt_script_uicallback,h1,'edit_prefix_Callback'},...
'Position',[430 209 111 21],...
'Style','edit',...
'TooltipString','Script and files name prefix',...
'Tag','edit_prefix');

h26 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[426 152 51 15],...
'String','Map width','Style','text','Tag','text2');

h27 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[426 122 53 15],...
'String','Map height','Style','text','Tag','text3');

h28 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[426 92 41 15],...
'String','X origin','Style','text','Tag','text4');

h29 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[426 62 41 15],...
'String','Y origin','Style','text','Tag','text5');

h30 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'pushbutton_change_dir_Callback'},...
'FontWeight','bold','Position',[520 241 21 21],...
'String','...',...
'Tag','pushbutton_change_dir');

h31 = uicontrol('Parent',h1,'Position',[120 9 3 51],'Style','frame','Tag','frame2');

h32 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[369 212 58 15],...
'String','Name prefix','Style','text','Tag','text6');

h33 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[331 182 76 15],...
'String','Map scale (apr)','Style','text','Tag','text8');

h34 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'togglebutton_Option_L_Callback'},...
'Position',[260 95 66 34],...
'Style','togglebutton',...
'TooltipString','Call a window to draw a map scale bar',...
'SelectionHighlight','off',...
'Tag','togglebutton_Option_L');

h35 = uicontrol('Parent',h1,'Enable','inactive','HorizontalAlignment','left',...
'Position',[269 100 48 27],...
'String',{'    Map'; 'scale bar'},...
'Style','text',...
'Tag','text_MapScale');

h36 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'togglebutton_Option_U_Callback'},...
'Position',[260 26 66 34],...
'Style','togglebutton',...
'TooltipString','Call a window to draw a Time Stamp',...
'SelectionHighlight','off',...
'Tag','togglebutton_Option_U');

h37 = uicontrol('Parent',h1,'Enable','inactive','HorizontalAlignment','left',...
'Position',[265 30 57 28],'String',{'Time Stamp'; '& signature' },...
'Style','text','Tag','text11');

h38 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'checkbox_removeOptionL_Callback'},...
'Position',[265 30 66 15],...
'String','Remove',...
'Style','checkbox',...
'TooltipString','Remove the scale bar',...
'Tag','checkbox_removeOptionL',...
'Visible','off');

h39 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'checkbox_removeOptionU_Callback'},...
'Position',[260 10 66 15],...
'String','Remove',...
'Style','checkbox',...
'TooltipString','Remove the Time Stamp',...
'Tag','checkbox_removeOptionU',...
'Visible','off');

h40 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'pushbutton_coastLines_Callback'},...
'Position',[360 408 181 23],...
'String','Apply finer control to coast lines',...
'Tag','pushbutton_coastLines');

h41 = uicontrol('Parent',h1,...
'Callback',{@write_gmt_script_uicallback,h1,'pushbutton_cancel_Callback'},...
'FontSize',9,...
'Position',[370 8 82 24],...
'String','Cancel',...
'Tag','pushbutton_cancel');

h42 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[340 281 201 71],...
'String','Nikles',...
'Style','text',...
'Tag','text_ProjDescription');

function write_gmt_script_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

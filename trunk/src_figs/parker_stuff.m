function varargout = parker_stuff(varargin)
% M-File changed by desGUIDE 
% Last Modified 31-May-2005

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
parker_stuff_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

global home_dir
movegui(hObject,'center')
% Case when this function was called directly
if isempty(home_dir),   home_dir = pwd;     end

if isempty(home_dir)        % Case when this function was called directly
    handles.path_data = ['data' filesep];
else
    handles.path_data = [home_dir filesep 'data' filesep];
end

% Import icons
load([handles.path_data 'mirone_icons.mat'],'Mfopen_ico');
set(handles.pushbutton_BatGrid,'CData',Mfopen_ico)
set(handles.pushbutton_SourceGrid,'CData',Mfopen_ico)
clear Mfopen_ico;

handles.h_calling_fig = [];     % Handles to the calling figure
handles.geog = 0;               % Set this as default
handles.is_redPole = 0;         % Flag to signal that options are for use in reduction to the pole
handles.Z_bat = [];
handles.Z_src = [];
handles.zobs = 0;               % Default observation level
handles.data = [];              % Date has no default value
last_dir = [];
no_igrf  = 0;                   % In reduction to the pole we don't need to compute the IGRF

if (~isempty(varargin))
    if (strcmp(varargin{1},'parker_direct'))
        h_txt = findobj(hObject,'Tag','text_FieldMag');
        set(h_txt,'String','Mag')
        set(handles.edit_SourceGrid,'TooltipString','Enter Magnetization grid name')
        set(handles.edit_wlong,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
        set(handles.edit_wshort,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
        set(handles.edit_sdec,'String','0')     % Default is geocentric dipole
        set(handles.edit_sdip,'String','0')
        set(handles.checkbox_CenterDipole,'Value',1)
        set(handles.checkbox_fieldIsRTP,'Visible','off')
        set(hObject,'Name','Parker Direct')
        handles.what_parker = 'direct';
        handles.h_calling_fig = varargin{2};
    elseif (strcmp(varargin{1},'parker_inverse'))
        set(handles.edit_sdec,'String','0')     % Default is geocentric dipole
        set(handles.edit_sdip,'String','0')
        set(handles.checkbox_CenterDipole,'Value',1)
        set(handles.checkbox_fieldIsRTP,'Value',0)
        set(hObject,'Name','Parker Inverse')
        handles.what_parker = 'inverse';        
        handles.h_calling_fig = varargin{2};
    elseif (strcmp(varargin{1},'redPole'))
        no_igrf = 1;
        set(handles.edit_BatGrid,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
        set(handles.pushbutton_BatGrid,'Enable','inactive')
        set(handles.edit_wlong,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
        set(handles.edit_wshort,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
        set(handles.edit_date,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
        set(handles.checkbox_CenterDipole,'Visible','off')
        set(handles.checkbox_fieldIsRTP,'Visible','off')
        h_txt = findobj(hObject,'Tag','text_wshort');
        set(h_txt,'Enable','off')
        h_txt = findobj(hObject,'Tag','text_wlong');
        set(h_txt,'Enable','off')
        % Reuse some edit boxes, but we have to change their text
        h_txt = findobj(hObject,'Tag','text_Level');
        set(h_txt,'String','Field dip')
        set(handles.edit_zobs,'String','','TooltipString','Inclination of the magnetic field.')
        h_txt = findobj(hObject,'Tag','text_Thickness');
        set(h_txt,'String','Field dec')
        set(handles.edit_thickness,'String','','TooltipString','Declination of the magnetic field.')
        set(hObject,'Name','Reduction to the Pole')
        handles.what_parker = 'redPole';
        handles.h_calling_fig = varargin{2};
    else        % defaults to "direct"
        handles.what_parker = 'direct';
    end
else
    % Else defaults to "direct" but without knowing h_calling_fig
    handles.what_parker = 'direct';
end

if (~no_igrf)
    handles.start_stop_epoch = [1900 2010];     % ISTO TEM DE SER AUTOMATIZADO
end

% first 2 cols from table III of Singleton's paper on fft
nlist = {64,72,75,80,81,90,96,100,108,120,125,128,135,144,150,160,162,180,192,200,...
    216,225,240,243,250,256,270,288,300,320,324,360,375,384,400,405,432,450,480,...
    486,500,512,540,576,600,625,640,648,675,720,729,750,768,800,810,864,900,960,...
    972,1000,1024,1080,1125,1152,1200,1215,1250,1280,1296,1350,1440,1458,1500,...
    1536,1600,1620,1728,1800,1875}';
set(handles.listbox_nnx,'String',nlist)
set(handles.listbox_nny,'String',nlist)

% Set upt some useful tooltips
str = sprintf(['The default value is the number of rows in the grid\n',...
    'However, for reducing border effects you may want to apply\n',...
    'a skirt to the grid. For that, select a value from the side\n',...
    'listbox. Extra points will be padded by mirroiong the west side.']);
set(handles.edit_Nrows,'TooltipString',str)
str = sprintf(['The default value is the number of cols in the grid\n',...
    'However, for reducing border effects you may want to apply\n',...
    'a skirt to the grid. For that, select a value from the side\n',...
    'listbox. Extra points will be padded by mirroiong the south side.']);
set(handles.edit_Ncols,'TooltipString',str)

str = sprintf('Good FFT numbers for padding the grid');
set(handles.listbox_nnx,'TooltipString',str)
set(handles.listbox_nny,'TooltipString',str)

% Give a Pro look (3D) to the frame boxes 
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

% Choose default command line output for parker_stuff_export
handles.output = hObject;
guidata(hObject, handles);
set(hObject,'Visible','on');
% UIWAIT makes parker_stuff_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = parker_stuff_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = parker_stuff_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure

% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------------------------------
function edit_BatGrid_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if isempty(fname)
    handles.Z_bat = [];    return;
end
% Let the pushbutton_BatGrid_Callback do all the work
parker_stuff('pushbutton_BatGrid_Callback',gcbo,[],guidata(gcbo),fname)

% -------------------------------------------------------------------------------------------------
function pushbutton_BatGrid_Callback(hObject, eventdata, handles, opt)
if (nargin == 4)    fname = opt;    end

if (isempty(opt))       % Otherwise 'opt' already transmited the file name.
    if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
        cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
        last_dir = cfig_handles.last_dir;
        home = cfig_handles.home_dir;
    else
        last_dir = [];
    end

    if (~isempty(last_dir)),    cd(last_dir);   end
    [FileName,PathName] = uigetfile({'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'},'Select GMT grid');
    pause(0.01);
    if (~isempty(last_dir)),    cd(home);   end
    if isequal(FileName,0);     return;     end
    fname = [PathName FileName];
end

% Because GMT and Surfer share the .grd extension, find out which kind grid we are dealing with
[fid, msg] = fopen(fname, 'r');
if (fid < 0)
    errordlg([PathName FileName ': ' msg],'ERROR'); return
end
ID = fread(fid,4,'*char');
ID = strread(ID,'%s');      fclose(fid);
if strcmp(ID,'DSBB') | strcmp(ID,'DSRB') 
    fname = [fname '=6'];
elseif strcmp(ID,'DSAA')
    warndlg('I don''t know and do not intend to learn how to read ASCII Surfer grids.','Warning')
    return
else        % It must (we hope) be a gmt grid
end
[X,Y,handles.Z_bat,handles.head_bat] = grdread_m(fname);
% See if Source/Mag grid is already loaded and, if yes, if they are compatible
if (~isempty(get(handles.edit_SourceGrid,'String')))
    if ( abs(handles.head_bat(1) - handles.head_src(1)) > 1e-4 | abs(handles.head_bat(2) - handles.head_src(2)) > 1e-4 |...
            abs(handles.head_bat(3) - handles.head_src(3)) > 1e-4 | abs(handles.head_bat(4) - handles.head_src(4)) > 1e-4)
        errordlg('Error: Bathymetry & Source grids do not cover the same region','Error');  return
    elseif(abs(handles.head_bat(8) - handles.head_src(8)) > 1e-6 | abs(handles.head_bat(9) - handles.head_src(9)) > 1e-6)
        errordlg('Error: Bathymetry & Source grids do not have the same size.','Error');     return
    end
end
% Try to guess if bat is in meters. If yes convert to km
if (abs(handles.head_bat(6) - handles.head_bat(5)) > 15)
    grdutils(handles.Z_bat,['-M' num2str(0.001,'%.3f')]);
    handles.head_bat(5) = handles.head_bat(5) / 1000;
    handles.head_bat(6) = handles.head_bat(6) / 1000;
end
set(handles.edit_BatGrid,'String',fname)
guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_SourceGrid_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if isempty(fname)
    handles.Z_src = [];    return;
end
% Let the pushbutton_SourceGrid_Callback do all the work
parker_stuff('pushbutton_InputFile_Callback',gcbo,[],guidata(gcbo),fname)

% -------------------------------------------------------------------------------------------------
function pushbutton_SourceGrid_Callback(hObject, eventdata, handles,opt)
if (nargin == 4)    fname = opt;    end

if (isempty(opt))       % Otherwise 'opt' already transmited the file name.
    if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
        cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
        last_dir = cfig_handles.last_dir;
        home = cfig_handles.home_dir;
    else
        last_dir = [];
    end

    if (~isempty(last_dir)),    cd(last_dir);   end
    [FileName,PathName] = uigetfile({'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'},'Select GMT grid');
    pause(0.01);
    if (~isempty(last_dir)),    cd(home);   end
    if isequal(FileName,0);     return;     end
    fname = [PathName FileName];
end

% Because GMT and Surfer share the .grd extension, find out which kind grid we are dealing with
[fid, msg] = fopen(fname, 'r');
if fid < 0
    errordlg([PathName FileName ': ' msg],'ERROR'); return
end
ID = fread(fid,4,'*char');
ID = strread(ID,'%s');
if strcmp(ID,'DSBB') | strcmp(ID,'DSRB') 
    fname = [fname '=6'];
elseif strcmp(ID,'DSAA')
    warndlg('I don''t know and do not intend to learn how to read ASCII Surfer grids.','Warning')
    return
end
[handles.X,handles.Y,handles.Z_src,handles.head_src] = grdread_m(fname);
% See if Bat grid is already loaded and, if yes, if both grids are compatible
if (~isempty(get(handles.edit_BatGrid,'String')))
    if ( abs(handles.head_bat(1) - handles.head_src(1)) > 1e-4 | abs(handles.head_bat(2) - handles.head_src(2)) > 1e-4 |...
            abs(handles.head_bat(3) - handles.head_src(3)) > 1e-4 | abs(handles.head_bat(4) - handles.head_src(4)) > 1e-4)
        errordlg('Error: Bathymetry & Source grids do not cover the same region','Error');  return
    elseif(abs(handles.head_bat(8) - handles.head_src(8)) > 1e-6 | abs(handles.head_bat(9) - handles.head_src(9)) > 1e-6)
        errordlg('Error: Bathymetry & Source grids do not have the same size.','Error');     return
    end
end

[handles.orig_nrows,handles.orig_ncols] = size(handles.Z_src);
handles.nrows = handles.orig_nrows;             % Make them equal by default
handles.ncols = handles.orig_ncols;
set(handles.edit_SourceGrid,'String',fname)
set(handles.edit_Nrows,'string',num2str(handles.nrows))
set(handles.edit_Ncols,'string',num2str(handles.ncols))

% Try to guess if grid is in geogs (restricting to 90 degrees spaning is more than enough as a test)
if (abs(handles.head_src(2)-handles.head_src(1)) < 90 | abs(handles.head_src(4)-handles.head_src(3)) < 90)
    handles.geog = 1;   % We probably have a geog grid
    set(handles.checkbox_geog,'Value',1)
end

% The easeast way of not leting the user screw things by selecting a nnx and/or nny lower
% than nx or ny is to delete the forbiden numbers from the listboxes 
contents = get(handles.listbox_nny,'String');
%xx=str2num(cat(1,contents{:}));
[m,n] = size(contents);
for (i=1:m)    xx(i) = str2num(contents{i});    end
ind = find(xx > handles.nrows);
nlist = num2cell(xx(ind),1);
set(handles.listbox_nny,'String',nlist)
ind = find(xx > handles.ncols);
nlist = num2cell(xx(ind),1);
set(handles.listbox_nnx,'String',nlist)
handles.rlon = (handles.head_src(2) + handles.head_src(1)) / 2;
handles.rlat = (handles.head_src(4) + handles.head_src(3)) / 2;

if (handles.geog)
    [sclat,sclon] = scltln(handles.rlat);
    dx = handles.head_src(8) * sclon;
    dy = handles.head_src(9) * sclat;
    handles.scaled_dx = dx;     handles.scaled_dy = dy;
else
    handles.scaled_dx = handles.head_src(8);     handles.scaled_dy = handles.head_src(9);
end

% Compute the wshort & wlong default values (only for the 'inverse' case)
if (strcmp(handles.what_parker,'inverse'))
    if (handles.geog)
        wshort = max(dx*2, dy*2);
        wlong = max(dx*handles.edit_Ncols, dy*handles.edit_Nrows);
    else
        wshort = max(2*handles.head_src(8),2*handles.head_src(9));
        wlong = max(handles.edit_Ncols*handles.head_src(8),handles.edit_Nrows*handles.head_src(9));
    end
    wlong = max(wlong,150);     % Beter use this as the wlong default 
    set(handles.edit_wshort,'string',num2str(wshort,'%.1f'))
    set(handles.edit_wlong,'string',num2str(wlong,'%.0f'))
end
guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function checkbox_mirror_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.edit_Ncols,'Enble','off')
    set(handles.listbox_nny,'Enble','off')
    set(handles.edit_Nrows,'Enble','off')
    set(handles.listbox_nnx,'Enble','off')
else
    set(handles.edit_Ncols,'Enble','on')
    set(handles.listbox_nny,'Enble','on')
    set(handles.edit_Nrows,'Enble','on')
    set(handles.listbox_nnx,'Enble','on')
end

% -------------------------------------------------------------------------------------------------
function edit_date_Callback(hObject, eventdata, handles)
if (strcmp(handles.what_parker,'redPole'))
    return;     % In this mode the box serves only to store a value
end
if (get(handles.checkbox_fieldIsRTP,'Value'))
    set(hObject,'String','')
    return
end
xx = str2double(get(hObject,'String'));
if (xx < handles.start_stop_epoch(1) | xx > handles.start_stop_epoch(2))
    errordlg('Date outside the current IGRF model limits','Error');
    handles.date = [];
    set(hObject,'String','');       return
else
    handles.date = xx;
    elev = str2double(get(handles.edit_zobs,'String'));
    try     % I use a try here because some vars may not have been yet defined
        %out = igrf(handles.rlat, handles.rlon, elev, xx, handles.igrf_coefs, handles.start_stop_epoch);
        out = igrf_m(handles.rlon, handles.rlat, elev, xx);
        if (~get(handles.checkbox_CenterDipole,'Value'))
            set(handles.edit_sdec,'String',num2str(out(6)))
            set(handles.edit_sdip,'String',num2str(out(7)))
        end
    end
end
guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_zobs_Callback(hObject, eventdata, handles)
zobs = str2double(get(hObject,'String'));
if (isnan(zobs))
    set(hObject,'String','0');      return;
else
    handles.zobs = zobs;    
end
guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_thickness_Callback(hObject, eventdata, handles)
if (strcmp(handles.what_parker,'redPole'))
    return;     % In this mode the box serves only to store a value
end
xx = str2double(get(hObject,'String'));
if (xx <= 0)
    errordlg('You must be dreaming. What is a layer with zero or negative thickness?','Chico Clever');
    set(hObject,'String','');       return
else
    handles.thick = xx;
end
guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_sdip_Callback(hObject, eventdata, handles)
if (get(handles.checkbox_fieldIsRTP,'Value'))
    set(hObject,'String','90');    return        % Don't let it be changed if Field is RTP
end
if (~strcmp(handles.what_parker,'redPole'))
	if (get(handles.checkbox_CenterDipole,'Value'))     % sdip = 0 for centered dipole
        set(hObject,'String','0');    return
	end
end
xx = str2double(get(hObject,'String'));
if (xx < -90 | xx > 90)
    errordlg('Inlinations are restricted to the [-90;+90] interval.','Error');
    set(hObject,'String','');       return
else
    handles.sdip = xx;
end
guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_sdec_Callback(hObject, eventdata, handles)
if (get(handles.checkbox_fieldIsRTP,'Value'))
    set(hObject,'String','0');    return        % Don't let it be changed if Field is RTP
end
if (~strcmp(handles.what_parker,'redPole'))
	if (get(handles.checkbox_CenterDipole,'Value'))     % sdec = 0 for centered dipole
        set(hObject,'String','0');    return
	end
end
xx = str2double(get(hObject,'String'));
if (xx < -90 | xx > 90)
    errordlg('Declinations are restricted to the [-90;+90] interval.','Error');
    set(hObject,'String','');       return
else
    handles.sdec = xx;
end
guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_wshort_Callback(hObject, eventdata, handles)
% Nothing to do. The "Compute" button will read this field

% -------------------------------------------------------------------------------------------------
function edit_wlong_Callback(hObject, eventdata, handles)
% Nothing to do. The "Compute" button will read this field

% -------------------------------------------------------------------------------------------------
function listbox_nnx_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns listbox_nnx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_nnx
contents = get(hObject,'String');
nnx = str2double(contents{get(hObject,'Value')});
set(handles.edit_Ncols,'String',num2str(nnx))
handles.ncols = nnx;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function listbox_nny_Callback(hObject, eventdata, handles)
contents = get(hObject,'String');
nny = str2double(contents{get(hObject,'Value')});
set(handles.edit_Nrows,'String',num2str(nny))
handles.nrows = nny;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_Ncols_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isempty(get(hObject,'String')))
    try,    set(hObject,'String',num2str(handles.ncols));   return;     end
end
if (xx < handles.cols)
    set(hObject,'String',num2str(handles.ncols));    return;
end
handles.ncols = xx;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_Nrows_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx))
    try,    set(hObject,'String',num2str(handles.nrows));   return;     end
end
if (xx < handles.nrows)
    set(hObject,'String',num2str(handles.nrows));   return;
end
handles.nrows = xx;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function checkbox_geog_Callback(hObject, eventdata, handles)
% Nothing to do

% -------------------------------------------------------------------------------------------------
function checkbox_fieldIsRTP_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.edit_sdip,'String','90');   set(handles.edit_sdec,'String','0')
    set(handles.edit_date,'String','');     set(handles.checkbox_CenterDipole,'Value',0)
else
    set(handles.edit_sdip,'String','0');   set(handles.edit_sdec,'String','0')
    set(handles.checkbox_CenterDipole,'Value',1)
end

% -------------------------------------------------------------------------------------------------
function checkbox_CenterDipole_Callback(hObject, eventdata, handles)
if (get(handles.checkbox_fieldIsRTP,'Value'))
    set(hObject,'Value',0)
end

% -------------------------------------------------------------------------------------------------
function pushbutton_compute_Callback(hObject, eventdata, handles)
% Before asking the apropriate function to do the work we have to ... TEST
% Source grid, Nrows & Ncols are common to all options. So test them first
if (isempty(handles.Z_src))
    errordlg('You didn''t give me a Source grid (Field or Magnetization). What do you want me to do?','Chico Clever')
    return
end
nrows = str2double(get(handles.edit_Nrows,'String'));
ncols = str2double(get(handles.edit_Ncols,'String'));
if (isempty(get(handles.edit_Nrows,'String')) | isempty(get(handles.edit_Ncols,'String')))
    errordlg('One or both of grid size dimensions are empty. What have you done?','Error')
    return
end

% Now those that are common to direct/inverse cases
if (strcmp(handles.what_parker,'direct') | strcmp(handles.what_parker,'inverse'))
    if (isempty(handles.Z_bat))
        errordlg('Must give me a grid with the bathymetry','Error');    return
    end
    date = str2double(get(handles.edit_date,'String'));
    if (isnan(date) & ~get(handles.checkbox_fieldIsRTP,'Value'))
        errordlg('I need to know the year of the survey (see Date box)','Error');   return
    end
    thick = str2double(get(handles.edit_thickness,'String'));
    if (isnan(thick))
        errordlg('I need to know the thickness of the magnetic layer (see Thickness box)','Error');   return
    end
    zobs = str2double(get(handles.edit_zobs,'String'));
    if (isnan(zobs))
        errordlg('I need to know the level of the observation of the survey (see Level box)','Error');   return
    end
    dx = handles.scaled_dx;     dy = handles.scaled_dy;
end
% Get Mag/Field Dec & Dip
sdec = str2double(get(handles.edit_sdec,'String'));
sdip = str2double(get(handles.edit_sdip,'String'));
new_nx = str2double(get(handles.edit_Ncols,'String'));
new_ny = str2double(get(handles.edit_Nrows,'String'));

switch handles.what_parker
    case 'direct'
        if (handles.orig_ncols < handles.ncols | handles.orig_nrows < handles.nrows)      % Padding was asked for
            if (get(handles.checkbox_mirror,'Value'))   % Do mirror
                h = mboard(handles.Z_bat,handles.orig_ncols,handles.orig_nrows);
                f = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows);
                f3d = syn3d(double(f),double(h),handles.rlat,handles.rlon,date,zobs,thick,0,dx,dy,sdip,sdec);
                f3d = f3d(1:handles.orig_nrows,1:handles.orig_ncols);   % Remove the mirror
            else
                [h,band] = mboard(handles.Z_bat,handles.orig_ncols,handles.orig_nrows,new_nx,new_ny);
                f = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows,new_nx,new_ny);
                f3d = syn3d(double(f),double(h),handles.rlat,handles.rlon,date,zobs,thick,0,dx,dy,sdip,sdec);
                m1 = band(1)+1;     m2 = m1 + handles.orig_nrows - 1;
                n1 = band(3)+1;     n2 = n1 + handles.orig_ncols - 1;
                f3d = f3d(m1:m2,n1:n2);         % Remove the padding skirt
            end
        else
            f3d = syn3d(double(handles.Z_src),double(handles.Z_bat),handles.rlat,handles.rlon, ...
                date,zobs,thick,0,dx,dy,sdip,sdec);
        end
        z_min = min(f3d(:));    z_max = max(f3d(:));
        tmp.head = [handles.head_src(1) handles.head_src(2) handles.head_src(3) handles.head_src(4) ...
                z_min z_max 0 handles.head_src(8) handles.head_src(9)];
        tmp.X = handles.X;      tmp.Y = handles.Y;      tmp.name = 'Magnetic Anomaly (nT)';
        new_window = mirone(single(f3d),tmp);
    case 'inverse'
        if (get(handles.checkbox_fieldIsRTP,'Value'))   % Case of a already RTP Field
            handles.rlat = 90;      handles.rlon = 0;
            sdip = 90;              sdec = 0;       date = 2000;
        end
        ws = str2double(get(handles.edit_wshort,'String'));
        wl = str2double(get(handles.edit_wlong,'String'));
        if (handles.orig_ncols < handles.ncols | handles.orig_nrows < handles.nrows)      % Padding was asked for
            if (get(handles.checkbox_mirror,'Value'))   % Do mirror
                h = mboard(handles.Z_bat,handles.orig_ncols,handles.orig_nrows);
                f = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows);
                m3d = inv3d(double(f),double(h),wl,ws,handles.rlat,handles.rlon,date,zobs,thick,0,dx,dy,sdip,sdec);
                m3d = m3d(1:handles.orig_nrows,1:handles.orig_ncols);   % Remove the mirror
            else
                [h,band] = mboard(handles.Z_bat,handles.orig_ncols,handles.orig_nrows,new_nx,new_ny);
                f = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows,new_nx,new_ny);
                m3d = inv3d(double(f),double(h),wl,ws,handles.rlat,handles.rlon,date,zobs,thick,0,dx,dy,sdip,sdec);
                m1 = band(1)+1;     m2 = m1 + handles.orig_nrows - 1;
                n1 = band(3)+1;     n2 = n1 + handles.orig_ncols - 1;
                m3d = m3d(m1:m2,n1:n2);         % Remove the padding skirt
            end
        else
            m3d = inv3d(double(handles.Z_src),double(handles.Z_bat),wl,ws,handles.rlat,handles.rlon, ...
                date,zobs,thick,0,dx,dy,sdip,sdec);
        end
        z_min = min(m3d(:));        z_max = max(m3d(:));
        tmp.head = [handles.head_src(1) handles.head_src(2) handles.head_src(3) handles.head_src(4) ...
                z_min z_max 0 handles.head_src(8) handles.head_src(9)];
        tmp.X = handles.X;      tmp.Y = handles.Y;      tmp.name = 'Magnetization (A/m^2)';
        new_window = mirone(single(m3d),tmp);
    case 'redPole'
        incl_fld = str2double(get(handles.edit_zobs,'String'));
        decl_fld = str2double(get(handles.edit_thickness,'String'));
        incl_mag = str2double(get(handles.edit_sdip,'String'));
        decl_mag = str2double(get(handles.edit_sdec,'String'));
        if (isnan(incl_fld) | isnan(decl_fld))
            errordlg('You need to give me valid magnetic field Inclination and Declination.','Error');  return;
        end
        if (isnan(incl_mag) | isnan(decl_mag))
            errordlg('You need to give me valid magnetization Inclination and Declination.','Error');  return;
        end
        if (handles.orig_ncols < handles.ncols | handles.orig_nrows < handles.nrows)      % Padding was asked for
            if (get(handles.checkbox_mirror,'Value'))   % Do mirror
                f = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows);
                f = rtp3d(double(f),incl_fld,decl_fld,incl_mag,decl_mag);
                f = f(1:handles.orig_nrows,1:handles.orig_ncols);   % Remove the mirror
            else
                [f,band] = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows,new_nx,new_ny);
                f = rtp3d(double(f),incl_fld,decl_fld,incl_mag,decl_mag);
                m1 = band(1)+1;     m2 = m1 + handles.orig_nrows - 1;
                n1 = band(3)+1;     n2 = n1 + handles.orig_ncols - 1;
                f = f(m1:m2,n1:n2);         % Remove the padding skirt
            end
        else
            f = rtp3d(double(handles.Z_src),incl_fld,decl_fld,incl_mag,decl_mag);
        end
        z_min = min(f(:));    z_max = max(f(:));
        tmp.head = [handles.head_src(1) handles.head_src(2) handles.head_src(3) handles.head_src(4) ...
                z_min z_max 0 handles.head_src(8) handles.head_src(9)];
        tmp.X = handles.X;      tmp.Y = handles.Y;      tmp.name = 'Reduction to the Pole anomaly (nT)';
        new_window = mirone(single(f),tmp);
end

% -------------------------------------------------------------------------------------------------
function [sclat,sclon] = scltln(orlat)
% Routine to determine lat-lon scales, km/deg, for ellipsoids
% of revolution,  using equations of:
%       Snyder, J.P., 1987, Map Projections -- A Working Manual,
%       USGS Professional Paper 1395, Washington DC, 383p. cf. pp 24-25.
%
% Currently, this is hard-wired for the WGS-84 ellipsoid.
%
% The basic equations are:
% 	sclat = a * (1-e*e)    /  (1 - e*e * sin(orlat)*sin(orlat))**(1.5)
%	sclon = a * cos(orlat) /  (1 - e*e * sin(orlat)*sin(orlat))**(0.5)
%
% where:    a  is the equatorial radius
%           b  is the polar radius
%           e  is the eccentricity
%           f  is the flattening
% Also:
%	e*e = 1. - b*b/a*a
%	f   = 1. - b/a
%
% Dan Scheirer, 21 May 1991

% These constants belong to the: WGS, 1984 ellipsoid (gmt_defaults.h)
a = 6378.137;   b = 6356.7521;

% Now, do the calculations...
e2 = 1 - (b*b)/(a*a);
sinlat = sin(orlat*pi/180);
denom  = sqrt(1 - e2 * sinlat * sinlat);
sclat = (pi/180) * a * (1 - e2)  / denom / denom / denom;
sclon = (pi/180) * a * cos(orlat*pi/180) / denom;


% --- Creates and returns a handle to the GUI figure. 
function parker_stuff_LayoutFcn(h1,handles);

set(h1, 'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','parker_stuff',...
'NumberTitle','off',...
'Position',[520 498 400 302],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

h2 = uicontrol('Parent',h1,'Position',[10 11 301 71],'String',{''},'Style','frame','Tag','frame3');
h3 = uicontrol('Parent',h1,'Position',[10 101 381 81],'String',{''},'Style','frame','Tag','frame2');
h4 = uicontrol('Parent',h1,'Position',[10 201 381 91],'String',{''},'Style','frame','Tag','frame1');

h5 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback4,h1,[],'edit_BatGrid_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 232 311 21],...
'Style','edit',...
'TooltipString','Enter bathymetry grid name (km +ve up)',...
'Tag','edit_BatGrid');

h6 = uicontrol('Parent',h1,...
'Callback',{@parker_stuff_uicallback4,h1,[],'pushbutton_BatGrid_Callback'},...
'CData',[],...
'Position',[361 231 23 23],...
'Tag','pushbutton_BatGrid',...
'UserData',[]);

h7 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback4,h1,[],'edit_SourceGrid_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 262 311 21],...
'Style','edit',...
'TooltipString','Enter magnetic field grid name',...
'Tag','edit_SourceGrid');

h8 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[17 236 31 15],...
'String','Bat','Style','text','Tag','text1');

h9 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[18 265 31 15],...
'String','Field','Style','text','Tag','text_FieldMag');

h10 = uicontrol('Parent',h1,...
'Callback',{@parker_stuff_uicallback4,h1,[],'pushbutton_SourceGrid_Callback'},...
'Position',[361 261 23 23],...
'Tag','pushbutton_SourceGrid');

h11 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback,h1,'edit_date_Callback'},...
'Position',[263 53 41 21],...
'Style','edit',...
'TooltipString','Decimal year for IGRF field calculation',...
'Tag','edit_date');

h12 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback,h1,'edit_thickness_Callback'},...
'Position',[162 53 41 21],...
'String','0.5',...
'Style','edit',...
'TooltipString','Thickness of magnetic source layer (km)',...
'Tag','edit_thickness');

h13 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback,h1,'edit_zobs_Callback'},...
'Position',[55 53 41 21],...
'String','0',...
'Style','edit',...
'TooltipString','Observation level above sealevel (km +ve up)',...
'Tag','edit_zobs');

h14 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback,h1,'edit_sdip_Callback'},...
'Position',[55 18 41 21],...
'Style','edit',...
'TooltipString','inclination of magnetization',...
'Tag','edit_sdip');

h15 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback,h1,'edit_sdec_Callback'},...
'Position',[161 18 41 21],...
'Style','edit',...
'TooltipString','declination of magnetization',...
'Tag','edit_sdec');

h16 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback,h1,'edit_wshort_Callback'},...
'Position',[335 150 41 21],...
'Style','edit',...
'TooltipString','filter short wavelength cutoff (km)',...
'Tag','edit_wshort');

h17 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback,h1,'edit_wlong_Callback'},...
'Position',[335 110 41 21],...
'Style','edit',...
'TooltipString','filter long wavelength cutoff (km)',...
'Tag','edit_wlong');

h18 = uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[229 56 31 15],...
'String','Date','Style','text','Tag','text_Date');

h19 = uicontrol('Parent',h1,'Position',[110 56 49 15],'String','Thickness','Style','text','Tag','text_Thickness');

h20 = uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[12 56 40 15],...
'String','Level','Style','text','Tag','text_Level');

h21 = uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[107 21 51 15],...
'String','Mag dec','Style','text','Tag','text_MagDec');

h22 = uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[12 22 40 15],...
'String','Mag dip','Style','text','Tag','text_MagDip');

h23 = uicontrol('Parent',h1,'Position',[291 153 41 15],'String','Wshort','Style','text','Tag','text_wshort');

h24 = uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[299 114 31 15],...
'String','Wlong','Style','text','Tag','text_wlong');

h25 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback,h1,'listbox_nny_Callback'},...
'Position',[126 102 51 78],'Style','listbox','Value',1,'Tag','listbox_nny');

h26 = uicontrol('Parent',h1,...
'Callback',{@parker_stuff_uicallback,h1,'pushbutton_compute_Callback'},...
'FontWeight','bold',...
'Position',[326 11 65 23],...
'String','Compute',...
'Tag','pushbutton_compute');

h27 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback,h1,'listbox_nnx_Callback'},...
'Position',[232 102 51 78],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_nnx');

h28 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback,h1,'edit_Nrows_Callback'},...
'Position',[86 135 41 21],...
'Style','edit',...
'TooltipString','Number of grid rows',...
'Tag','edit_Nrows');

h29 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@parker_stuff_uicallback,h1,'edit_Ncols_Callback'},...
'Position',[192 135 41 21],...
'Style','edit',...
'TooltipString','Number of grid columns',...
'Tag','edit_Ncols');

h30 = uicontrol('Parent',h1,'Position',[87 159 39 15],'String','# Rows','Style','text','Tag','text11');
h31 = uicontrol('Parent',h1,'Position',[192 159 39 15],'String','# Cols','Style','text','Tag','text12');

h32 = uicontrol('Parent',h1,...
'Callback',{@parker_stuff_uicallback,h1,'checkbox_geog_Callback'},...
'Position',[50 210 117 15],...
'String','Geographic coords?',...
'Style','checkbox',...
'TooltipString','Are the grids in geographical coordinates?',...
'Tag','checkbox_geog');

h33 = uicontrol('Parent',h1,...
'Callback',{@parker_stuff_uicallback,h1,'checkbox_CenterDipole_Callback'},...
'Position',[227 19 77 15],...
'String','Geocentric?',...
'Style','checkbox',...
'TooltipString','Check this to assume geocentric dipole hypothesis',...
'Tag','checkbox_CenterDipole');

h34 = uicontrol('Parent',h1,...
'Callback',{@parker_stuff_uicallback,h1,'checkbox_mirror_Callback'},...
'Position',[19 139 50 15],...
'String','Mirror',...
'Style','checkbox',...
'TooltipString','Check this to Mirror the grid before FFT',...
'Tag','checkbox_mirror');

h35 = uicontrol('Parent',h1,...
'Callback',{@parker_stuff_uicallback,h1,'checkbox_fieldIsRTP_Callback'},...
'Position',[225 210 117 15],...
'String','Field is already RTP',...
'Style','checkbox',...
'TooltipString','Check this box if the anomalous field is already Reduced To the Pole',...
'Tag','checkbox_fieldIsRTP');

function parker_stuff_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

function parker_stuff_uicallback4(hObject, eventdata, h1, opt, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1),opt);

function varargout = overview(varargin)
% This GUI builds an overview of grids readed by GDAL
% The selected grid is sub-sampled, during the
% reading stage, in order to obtain an overview of
% about 200 rows/columns. This way, arbitrarely
% large grids may be previewed, but the quality of
% the preview will may degrade with grid size.

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

if (nargin >= 4 & isstr(varargin{1}))
    gui_Callback = str2func(varargin{1});
    feval(gui_Callback,varargin{2:end})
else
    h = overview_OpeningFcn(varargin{:});
    if (nargout)    varargout{1} = h;   end
end

% --- Executes just before overview is made visible.
function hObject = overview_OpeningFcn(varargin)
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to overview (see VARARGIN)
hObject = overview_LayoutFcn;
handles = guihandles(hObject);
movegui(hObject,'north');

% Default values for those
handles.home_dir = pwd;
handles.last_dir = handles.home_dir;
handles.work_dir = handles.home_dir;
handles.ForceInsitu = 0;

if (length(varargin) > 0)
    handles.ForceInsitu = varargin{1}.ForceInsitu;
    handles.home_dir = varargin{1}.home_dir;
    handles.last_dir = varargin{1}.last_dir;
    handles.work_dir = varargin{1}.work_dir;
end

handles.path_tmp = [handles.home_dir filesep 'tmp' filesep];
handles.icon_img = [];
handles.name_uncomp = [];
handles.was_int16 = 0;
handles.Nodata_int16 = [];

load([handles.home_dir filesep 'data' filesep 'mirone_icons.mat'],'rectang_ico','info_ico','boreas_ico','help_ico','zoom_ico');

h_toolbar = uitoolbar('parent',hObject,'Clipping', 'on', 'BusyAction','queue','HandleVisibility','on',...
   'Interruptible','on','Tag','FigureToolBar','Visible','on');
uipushtool('parent',h_toolbar,'Click',@rectang_clickedcallback,'Tag','rectang','cdata',rectang_ico,...
   'TooltipString','Draw Rectangle','Separator','on');
uipushtool('parent',h_toolbar,'Click',@info_clickedcallback,'Tag','info','TooltipString','Grid Info','cdata',info_ico);
uipushtool('parent',h_toolbar,'Click',@all_clickedcallback,'TooltipString','Swallow entire grid','cdata',boreas_ico);
uitoggletool('parent',h_toolbar,'Click',@zoom_clickedcallback,'Tag','zoom','TooltipString','Zoom','cdata',zoom_ico);
uipushtool('parent',h_toolbar,'Click',@help_clickedcallback,'Tag','help','TooltipString','Help','cdata',help_ico);

logo = gdalread(['data' filesep 'logo.png']);
image('parent',handles.axes1,'CData',flipdim(logo,1));
%set(handles.axes1,'YDir','reverse');

set(hObject,'Visible','on','HandleVisibility','callback');

% Choose default command line output for overview
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------------------------------------
function rectang_clickedcallback(obj,eventdata)
handles = guidata(obj);     % get handles
if (isempty(handles.icon_img))      return;     end
try
    hl = findobj(handles.figure1,'Type','patch');
    if (~isempty(hl))       return;     end     % We already have a rectangle so don't need another.
    [p1,p2,hl] = rubberbandbox;
    % Round rectangle coords to the nearest node
    p1(1) = round((p1(1) - handles.head_orig(1)) / handles.head_orig(8)) * handles.head_orig(8) + handles.head_orig(1);
    p2(1) = round((p2(1) - handles.head_orig(1)) / handles.head_orig(8)) * handles.head_orig(8) + handles.head_orig(1);
    p1(2) = round((p1(2) - handles.head_orig(3)) / handles.head_orig(9)) * handles.head_orig(9) + handles.head_orig(3);
    p2(2) = round((p2(2) - handles.head_orig(3)) / handles.head_orig(9)) * handles.head_orig(9) + handles.head_orig(3);
    delete(hl);     % Given that we cannot change a line into a patch, the easy way is to delete and reborn it
    hl = patch('XData',[p1(1) p1(1) p2(1) p2(1) p1(1)],'YData',[p1(2) p2(2) p2(2) p1(2) p1(2)],...
        'FaceColor','none','EdgeColor','k','LineWidth',1);

    cmenuHand = uicontextmenu;
    set(hl, 'UIContextMenu', cmenuHand);
    uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
    ui_edit_polygon(hl)    % Set edition functions
    uimenu(cmenuHand, 'Label', 'Rectangle limits', 'Separator','on', 'Callback', {@rectangle_limits,hl});
    uimenu(cmenuHand, 'Label', 'Crop Grid', 'Callback', {@CropaGrid,hl});
catch       % Don't know why but uisuspend sometimes breaks
    set(handles.figure1,'Pointer','arrow');
end

% --------------------------------------------------------------------------------------------------
function info_clickedcallback(obj,eventdata)
handles = guidata(obj);     % get handles
if (isempty(handles.icon_img))      return;     end
if (handles.image_type == 1)
    info_m(['grdinfo ' handles.fname]);
elseif (handles.image_type == 2)
    w{1} = ['   Xmin:  ' num2str(handles.head(1)) '    Xmax: ' num2str(handles.head(2))];
    w{2} = ['   Ymin:  ' num2str(handles.head(3)) '    Ymax: ' num2str(handles.head(4))];
    w{3} = ['   Zmin:  ' num2str(handles.head(5)) '    Zmax: ' num2str(handles.head(6))];
    w{4} = ['   Xinc:  ' num2str(handles.head(8)) '    Yinc: ' num2str(handles.head(9))];
    one_or_zero = ~(handles.head(7) == 1);      % To give correct nx,ny with either grid or pixel registration
    nx = round((handles.head(2) - handles.head(1))/handles.head(8) + one_or_zero);
    ny = round((handles.head(4) - handles.head(3))/abs(handles.head(9)) + one_or_zero);
    w{5} = ['   nx:  ' num2str(nx) '    ny: ' num2str(ny)];
    msgbox(w,'Image Info');
else
    w{1} = handles.Hdr.Driver;
    w{2} = ['Width:  ' num2str(handles.Hdr.Dim_nx) '    Height:  ' num2str(handles.Hdr.Dim_ny)];
    w{3} = ['Projection:  ' handles.Hdr.projection];
    w{4} = ['Datum:  ' handles.Hdr.datum];
    w{5} = ['Ellipsoide:  ' handles.Hdr.ellipsoid];
    w{6} = ['Pizel Size:  (' num2str(handles.Hdr.pixel_x) ',' num2str(handles.Hdr.pixel_y) ')'];
    w{7} = 'Projected corner coordinates';
    w{8} = ['   Xmin:  ' num2str(handles.Hdr.LL_prj_xmin) '    Xmax: ' num2str(handles.Hdr.UR_prj_xmax)];
    w{9} = ['   Ymin:  ' num2str(handles.Hdr.LL_prj_ymin) '    Ymax: ' num2str(handles.Hdr.UR_prj_ymax)];
    if ~isempty(handles.Hdr.UL_geo_xmin)
        w{10} = 'Geographical corner coordinates';
        w{11} = ['   Xmin:  ' num2str(handles.Hdr.LL_geo_xmin) '    Xmax: ' num2str(handles.Hdr.LR_geo_xmax)];
        w{12} = ['   Ymin:  ' num2str(handles.Hdr.LL_geo_ymin) '    Ymax: ' num2str(handles.Hdr.UR_geo_ymax)];
    end
    try,        w{13} = ['Zmin:  ' handles.Hdr.Zmin '   Zmax: ' handles.Hdr.Zmax];    end
    w{14} = ['Color Type:  ' handles.Hdr.Cmap];
    w{15} = [];
    [PATH,FNAME,EXT] = fileparts(handles.fname);
    
    w{16} = ['NOTE: The file "' handles.path_tmp FNAME '.info" contains more detailed information about this Image'];
    msgbox(w,'Image Info');
end

% --------------------------------------------------------------------------------------------------
function help_clickedcallback(obj,eventdata)
str = sprintf(['This tool is meant to be used on very large grids.\n'...
    'The selected grid is sub-sampled, during the\n'...
    'reading stage, in order to obtain an overview of\n'...
    'about 200 rows/columns. This way, arbitrarily\n'...
    'large grids may be previewed, but the quality of\n'...
    'the preview will degrade with grid size.\n\n'...
    'Click on the rectangle icon to draw a rectangle\n'...
    'on the image. Double clicking the rectangle lets\n'...
    'you edit it. A precise control of the rect size is\n'...
    'available by a right-click and selecting "rectangle\n'...
    '-limits". Selecting "Crop Grid" extracts the region\n'...
    'inside the rectangle at the grid''s full resolution.']);
helpdlg(str,'Help')

% --------------------------------------------------------------------------------------------------
function all_clickedcallback(obj,eventdata)
handles = guidata(obj);     % get handles
if (isempty(handles.icon_img))      return;     end
CropaGrid(obj,eventdata,[],'all')

% --------------------------------------------------------------------------------------------------
function zoom_clickedcallback(obj,eventdata)
handles = guidata(obj);     % get handles
if (isempty(handles.icon_img))
    set(obj,'State','off');    return;
end
if strcmp(get(obj,'State'),'on')
    zoom_j('on');
else
    zoom_j('off');
end

% -----------------------------------------------------------------------------------------
function rectangle_limits(obj,eventdata,h)
% Change the Rectangle's limits by asking it's corner coordinates
x = get(h,'XData');   y = get(h,'YData');

region = bg_region('with_limits',[x(1) x(3) y(1) y(3)]);
if isempty(region),    return;  end     % User gave up
x_min = region(1);      x_max = region(2);
y_min = region(3);      y_max = region(4);
set(h, 'XData', [x_min,x_min,x_max,x_max,x_min], 'YData', [y_min,y_max,y_max,y_min,y_min]);

% --------------------------------------------------------------------
function FileOpenArcGrid_Callback(hObject, eventdata, handles, opt)
if (strcmp(opt,'ascii'))
    str1 = {'*.grd;*.GRD', 'Arc/Info grid (*.grd,*.GRD)'; '*.*', 'All Files (*.*)'};
    str2 = 'Select Arc/Info ASCII grid File';
    handles.type = 'ArcAscii';
else        % Binary
    str1 = {'*.adf;*.ADF', 'Arc/Info grid (*.adf,*.ADF)'; '*.*', 'All Files (*.*)'};
    str2 = 'Select Arc/Info Binary grid File';
    handles.type = 'ArcBinary';
end

[FileName,PathName] = put_or_get_file(handles,str1,str2,'get');
if isequal(FileName,0);     return;     end
fname = [PathName FileName];
att =  gdalread(fname,'-M','-C');
m = att.RasterXSize;            n = att.RasterYSize;
handles.head = att.GMT_hdr;     handles.head_orig = handles.head;

% The following serves only for storing file infos
str = ['gdalinfo ' fname];
[s,w] = mat_lyies(str,[handles.path_tmp FileName '.info']);     % FileName may retain the comp ext
if ~(isequal(s,0))                  % An error has occured
    msg = ['Error getting file info (gdalinfo).', sprintf('\n\n'), 'Error message was:', sprintf('\n'), w];
    h = errordlg(msg,'Error');    movegui(h,'north');     return
end
set(handles.figure1,'pointer','watch')
file_info = dataread('file',[handles.path_tmp FileName '.info'],'%s','delimiter','\n','whitespace','');
Hdr = read_gdal_info(file_info);
Hdr.LL_prj_xmin = num2str(att.Corners.LL(1));   Hdr.LL_prj_ymin = num2str(att.Corners.LL(2));
Hdr.UR_prj_xmax = num2str(att.Corners.UR(1));   Hdr.UR_prj_ymax = num2str(att.Corners.UR(2));
Hdr.Zmin = num2str(handles.head(5));            Hdr.Zmax = num2str(handles.head(6));

jump = min(round(m / 200), round(n / 200));
opt_P = ['-P' num2str(jump)];
Z =  gdalread(fname,'-U','-S',opt_P);      [m,n] = size(Z);
X = linspace(handles.head(1),handles.head(2),n);
Y = linspace(handles.head(3),handles.head(4),m);
colormap(handles.axes1,jet(256));       handles.icon_img = image(X,Y,Z);
set(handles.axes1,'YDir','normal')
Resize1(handles.axes1, handles.icon_img)
handles.fname = fname;
handles.image_type = 4;
handles.Hdr = Hdr;
set(handles.figure1,'pointer','arrow')
guidata(hObject,handles)

% --------------------------------------------------------------------
function FileOpenGMTgrid_Callback(hObject, eventdata, handles)
% Read a GMT grid and convert it into an indexed image
if (isunix)
    warndlg('This option works only under windows. The culprit is netCDF that doesn''t build a shared library in other OS than Win (long life to *nix wonders))','Warning')
    return
end
str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select GMT grid','get');
if isequal(FileName,0);     return;     end
fname = [PathName FileName];

% Because GMT and Surfer share the .grd extension, find out which kind grid we are dealing with
[fid, msg] = fopen(fname, 'r');
if (fid < 0)    errordlg([PathName FileName ': ' msg],'ERROR'); return
end
ID = fread(fid,4,'*char');      ID = strread(ID,'%s');      fclose(fid);
if (strcmp(ID,'DSBB') | strcmp(ID,'DSAA') | strcmp(ID,'DSRB'))
    warndlg('Reading Surfer grids is not supported here.','Warning')
    return
elseif strcmp(ID,'Mode')        % Model Vision Grid
    warndlg('Reading ENCOM grids is not supported here.','Warning')
    return
else        % It must (we hope) be a gmt grid
end

set(handles.figure1,'pointer','watch')
head = grdinfo_m(fname,'silent');       % Use grdinfo because gdal gives wrong info on GMT grids
if (head(7))                            % Convert to grid registration
        head(1) = head(1) + head(8) / 2;        head(2) = head(2) - head(8) / 2;
        head(3) = head(1) + head(9) / 2;        head(4) = head(4) - head(9) / 2;
        head(7) = 0;
end
att =  gdalread(fname,'-M','-C');
is_coards = false;
if (att.RasterCount == 0)       % A GMT new format or a generic netcdf grid. Try luck.
    try
        ds = att.Subdatasets{1};
        ind = strfind(ds, '=');
        if (isempty(ind))
            errordlg('Whoops. Could not find a subdataset description. Quiting.','ERROR');  return
        end
        fname = ds(ind+1:end);
        att =  gdalread(fname,'-M','-C');
        is_coards = true;
    catch
        set(handles.figure1,'pointer','arrow')
        w{1} = 'Sorry. Something screw up with this NETCDF (new GMT version grid?) file';
        w{2} = lasterr;
        errordlg(w,'ERROR');        return
    end
end
m = att.RasterXSize;            n = att.RasterYSize;
% handles.head = att.GMT_hdr;
% handles.head(8) = head(8);              % Correct the wrong x_inc/y_inc info given by gdal
% handles.head(9) = head(9);

handles.head = head;
handles.head_orig = handles.head;
jump = min(round(m / 200), round(n / 200));
opt_P = ['-P' num2str(jump)];
if (~is_coards)
    Z =  gdalread(fname,'-U','-S',opt_P);
else    % In this case GDAL seams incapable of determine band's min/max and therefore no scaling
    Z =  gdalread(fname,'-U',opt_P);
    Z = scaleto8(Z);
end
[m,n] = size(Z);
X = linspace(handles.head(1),handles.head(2),n);
Y = linspace(handles.head(3),handles.head(4),m);
colormap(handles.axes1,jet(256));       handles.icon_img = image(X,Y,Z);
set(handles.axes1,'YDir','normal')
Resize1(handles.axes1, handles.icon_img)
handles.fname = fname;                  handles.type = 'GMT';
handles.image_type = 1;
set(handles.figure1,'pointer','arrow')
guidata(hObject,handles)

% --------------------------------------------------------------------
function FileOpen_DEM_Callback(hObject, eventdata, handles, opt)
% Files of the following formats are read (well directed) here
switch opt
    case 'DTED'
        str1 = {'*.dt0;*.DT0;*.dt1;*.DT1', 'DTED (*.dt0,*.DT0,*.dt1,*.DT1)'; '*.*', 'All Files (*.*)'};
        handles.type = 'DTED';          opt_C = ' ';
    case 'ESRI_hdr'
	    str1 = {'*.bil;*.BIL;', 'ESRI BIL (*.bil,*.BIL)'; '*.*', 'All Files (*.*)'};
        handles.type = 'ESRI_hdr';      opt_C = ' ';
    case 'GTOPO30'
        str1 = {'*.dem;*.DEM', 'GTOPO30 DEM (*.dem,*.DEM)'; '*.*', 'All Files (*.*)'};
        handles.type = 'GTOPO30';       opt_C = ' ';
    case 'GeoTiff_DEM'
        str1 = {'*.tif;*.TIF;*.tiff;*.TIFF', 'GeoTiff DEM(*.tif,*.tiff,*.TIF,*.TIFF)'; '*.*', 'All Files (*.*)'};
        handles.type = 'GeoTiff_DEM';   opt_C = ' ';
    case 'USGS_DEM'
        str1 = {'*.dem;*.DEM', 'USGS DEM (*.dem,*.DEM)'; '*.*', 'All Files (*.*)'};
        handles.type = 'USGS_DEM';      opt_C = '-C';
    case 'SDTS'
        str1 = {'*catd.ddf;*CATD.DDF', 'USGS SDTS DEM (*catd.ddf,*CATD.DDF)'; '*.*', 'All Files (*.*)'};
        handles.type = 'SDTS';          opt_C = '-C';
    case 'GXF'
        str1 = {'*.gxf;*.GXF', 'Geosoft GXF (*.gxf,*.GXF)'; '*.*', 'All Files (*.*)'};
        handles.type = 'GXF';           opt_C = ' ';
    otherwise
        return
end
str2 = ['Select ' opt ' File'];
[FileName,PathName] = put_or_get_file(handles,str1,str2,'get');
if isequal(FileName,0);     return;     end
fname = [PathName FileName];

att =  gdalread(fname,'-M','-C');
m = att.RasterXSize;            n = att.RasterYSize;
handles.head = att.GMT_hdr;     handles.head_orig = handles.head;

% The following serves only for storing file infos
str = ['gdalinfo ' fname];
[s,w] = mat_lyies(str,[handles.path_tmp FileName '.info']);     % FileName may retain the comp ext
if ~(isequal(s,0))                  % An error has occured
    msg = ['Error getting file info (gdalinfo).', sprintf('\n\n'), 'Error message was:', sprintf('\n'), w];
    h = errordlg(msg,'Error');    movegui(h,'north');     return
end
set(handles.figure1,'pointer','watch')
file_info = dataread('file',[handles.path_tmp FileName '.info'],'%s','delimiter','\n','whitespace','');
Hdr = read_gdal_info(file_info);
Hdr.LL_prj_xmin = num2str(att.Corners.LL(1));   Hdr.LL_prj_ymin = num2str(att.Corners.LL(2));
Hdr.UR_prj_xmax = num2str(att.Corners.UR(1));   Hdr.UR_prj_ymax = num2str(att.Corners.UR(2));
Hdr.Zmin = num2str(handles.head(5));            Hdr.Zmax = num2str(handles.head(6));

jump = min(round(m / 200), round(n / 200));
opt_P = ['-P' num2str(jump)];
Z =  gdalread(fname,'-U','-S',opt_P);      [m,n] = size(Z);
X = linspace(handles.head(1),handles.head(2),n);
Y = linspace(handles.head(3),handles.head(4),m);
colormap(handles.axes1,jet(256));           handles.icon_img = image(X,Y,Z);
set(handles.axes1,'YDir','normal')
Resize1(handles.axes1, handles.icon_img)
handles.fname = fname;
handles.image_type = 4;
handles.Hdr = Hdr;
set(handles.figure1,'pointer','arrow')
guidata(hObject,handles)

% --------------------------------------------------------------------
function FileOpenSRTM_Callback(hObject, eventdata, handles, opt)
if (nargin == 3)    opt = 'SRTM3';   end
if (strcmp(opt,'SRTM30'))
	str1 = {'*.srtm;*.SRTM;*.srtm.gz', 'SRTM30 DEM (*.srtm,*.SRTM,*.srtm.gz)'; '*.*', 'All Files (*.*)'};
else
	str1 = {'*.hgt;*.HGT;*.hgt.zip', [opt ' DEM (*.hgt,*.HGT,*.hgt.zip)']; '*.*', 'All Files (*.*)'};
end
[FileName,PathName] = put_or_get_file(handles,str1,'Select SRTM DEM File','get');
if isequal(FileName,0);     return;     end

fname = [PathName FileName];
tipo = opt;     handles.type = opt;
name_uncomp = [];

[name_hdr,comp_type] = write_ESRI_hdr(fname,tipo);
if (~isempty(comp_type))
    fname = decompress(fname,'warn');           % Name with compression ext removed
    name_uncomp = fname;        % Here we need a copy of the decompressed file name for removing
    if (isempty(fname))   return;     end;      % Error message already issued.
end        

set(handles.figure1,'pointer','watch')
att =  gdalread(fname,'-M','-C');
m = att.RasterXSize;            n = att.RasterYSize;
handles.head = att.GMT_hdr;     handles.head_orig = handles.head;

jump = min(round(m / 200), round(n / 200));
opt_P = ['-P' num2str(jump)];
Z =  gdalread(fname,'-U','-S',opt_P,'-C');      [m,n] = size(Z);
X = linspace(handles.head(1),handles.head(2),n);
Y = linspace(handles.head(3),handles.head(4),m);
colormap(handles.axes1,jet(256));           handles.icon_img = image(X,Y,Z);
set(handles.axes1,'YDir','normal')
Resize1(handles.axes1, handles.icon_img)
handles.fname = fname;
handles.image_type = 2;             % Flags that info is minimal
handles.name_uncomp = name_uncomp;  % If not empty this will signal read_DEMs to use the uncompressed file
handles.name_hdr = name_hdr;
set(handles.figure1,'pointer','arrow')

guidata(hObject,handles)

% --------------------------------------------------------------------
function FileOpenMOLA_Callback(hObject, eventdata, handles)
str1 = {'*.img;*.IMG', 'MOLA DEM (*.img,*.IMG)'; '*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select MOLA DEM File','get');
if isequal(FileName,0);     return;     end

handles.type = 'MOLA';
name_uncomp = [];

[PATH,FNAME,EXT] = fileparts([PathName FileName]);
fname = [PATH filesep FNAME '.lbl'];
fp = fopen(fname,'rt');
if (fp < 0)
    errordlg(['ERROR: Could not find format descriptive file: ' fname],'Error');  return
end
s = strread(fread(fp,'*char').','%s','delimiter','\n');

LINES = findcell('LINES', s);
if (isempty(LINES))   error = 1;  end
[t,r] = strtok(s{LINES.cn},'=');                n_lines = str2num(r(3:end));

LINE_SAMPLES = findcell('LINE_SAMPLES', s);
if (isempty(LINE_SAMPLES))   error = 1;  end
[t,r] = strtok(s{LINE_SAMPLES.cn},'=');         n_samples = str2num(r(3:end));

SAMPLE_BITS = findcell('SAMPLE_BITS', s);
if (isempty(SAMPLE_BITS))   error = 1;  end
[t,r] = strtok(s{SAMPLE_BITS.cn},'=');          n_bits = str2num(r(3:end));

CENTER_LATITUDE = findcell('CENTER_LATITUDE', s);
if (isempty(CENTER_LATITUDE))   error = 1;  end
[t,r] = strtok(s{CENTER_LATITUDE.cn},'=');      lat0 = str2num(strtok(strtok(r,'=')));

CENTER_LONGITUDE = findcell('CENTER_LONGITUDE', s);
if (isempty(CENTER_LONGITUDE))   error = 1;  end
[t,r] = strtok(s{CENTER_LONGITUDE.cn},'=');     lon0 = str2num(strtok(strtok(r,'=')));

LINE_FIRST_PIXEL = findcell('LINE_FIRST_PIXEL', s);
if (isempty(LINE_FIRST_PIXEL))   error = 1;  end
[t,r] = strtok(s{LINE_FIRST_PIXEL.cn},'=');     line_first_pix = str2num(r(3:end));

LINE_LAST_PIXEL = findcell('LINE_LAST_PIXEL', s);
if (isempty(LINE_LAST_PIXEL))   error = 1;  end
[t,r] = strtok(s{LINE_LAST_PIXEL.cn},'=');      line_last_pix = str2num(r(3:end));

SAMPLE_FIRST_PIXEL = findcell('SAMPLE_FIRST_PIXEL', s);
if (isempty(SAMPLE_FIRST_PIXEL))   error = 1;  end
[t,r] = strtok(s{SAMPLE_FIRST_PIXEL.cn},'=');   sample_first_pix = str2num(r(3:end));

SAMPLE_LAST_PIXEL = findcell('SAMPLE_LAST_PIXEL', s);
if (isempty(SAMPLE_LAST_PIXEL))   error = 1;  end
[t,r] = strtok(s{SAMPLE_LAST_PIXEL.cn},'=');    sample_last_pix = str2num(r(3:end));

MAP_RESOLUTION = findcell('MAP_RESOLUTION', s);
if (isempty(MAP_RESOLUTION))   error = 1;  end
[t,r] = strtok(s{MAP_RESOLUTION.cn},'=');       res = str2num(strtok(strtok(r,'=')));

LINE_PROJECTION_OFFSET = findcell('LINE_PROJECTION_OFFSET', s);
if (isempty(LINE_PROJECTION_OFFSET))   error = 1;  end
[t,r] = strtok(s{LINE_PROJECTION_OFFSET.cn},'=');line_off = str2num(r(3:end));

SAMPLE_PROJECTION_OFFSET = findcell('SAMPLE_PROJECTION_OFFSET', s);
if (isempty(SAMPLE_PROJECTION_OFFSET))   error = 1;  end
[t,r] = strtok(s{SAMPLE_PROJECTION_OFFSET.cn},'=');sample_off = str2num(r(3:end));

limits = [(sample_first_pix - sample_off)/res+lon0 (sample_last_pix - sample_off)/res+lon0 ...
    (line_first_pix - line_off)/res+lat0 (line_last_pix - line_off)/res+lat0];

fname = [PathName FileName];

opt = [limits n_samples n_lines 1/res];    % Must inform write_ESRI_hdr what to write
tipo = [opt(1) opt(4) opt(5) opt(6) opt(7) opt(7) -99999];

[name_hdr,comp_type] = write_ESRI_hdr(fname,tipo);
if (~isempty(comp_type))
    fname = decompress(fname,'warn');           % Name with compression ext removed
    name_uncomp = fname;        % Here we need a copy of the decompressed file name for removing
    if (isempty(fname))   return;     end;      % Error message already issued.
end

set(handles.figure1,'pointer','watch')
att =  gdalread(fname,'-M','-C');
m = att.RasterXSize;            n = att.RasterYSize;
handles.head = att.GMT_hdr;     handles.head_orig = handles.head;

jump = min(round(m / 200), round(n / 200));
opt_P = ['-P' num2str(jump)];
Z =  gdalread(fname,'-U','-S',opt_P);      [m,n] = size(Z);
X = linspace(handles.head(1),handles.head(2),n);
Y = linspace(handles.head(3),handles.head(4),m);
colormap(handles.axes1,jet(256));           handles.icon_img = image(X,Y,Z);
set(handles.axes1,'YDir','normal')
Resize1(handles.axes1, handles.icon_img)
handles.fname = fname;
handles.image_type = 2;             % Flags that info is minimal
handles.name_uncomp = name_uncomp;  % If not empty this will signal read_DEMs to use the uncompressed file
handles.name_hdr = name_hdr;
set(handles.figure1,'pointer','arrow')
guidata(hObject,handles)

% --------------------------------------------------------------------
function FileOpen_ENVI_Erdas_Callback(hObject, eventdata, handles, opt)
% This function reads both ENVI or Erdas files. Furthermore, based on the file byte
% type it guesses if we are dealing with a typical grid file (in which case it is
% treated like a native gmt grid) or a raster image file.

if (strcmp(opt,'ENVI'))
    str1 = {'*.img;*.IMG', 'ENVI (*.img,*.IMG)'; '*.*', 'All Files (*.*)'};
    str2 = 'Select ENVI file';
    handles.type = 'ENVI';
elseif (strcmp(opt,'Erdas'))
    str1 = {'*.img;*.IMG', 'Erdas (*.img,*.IMG)'; '*.*', 'All Files (*.*)'};
    str2 = 'Select Erdas file';
    handles.type = 'Erdas';
else    % Stupid error
    return
end
if (handles.ForceInsitu)    opt_I = '-I';   % Use only in desperate cases.
else                        opt_I = ' ';    end
handles.was_int16 = 0;  % To make sure that it wasnt left = 1 from a previous use.

[FileName,PathName] = put_or_get_file(handles,str1,str2,'get');
if isequal(FileName,0);     return;     end
fname = [PathName FileName];

str = ['gdalinfo ' fname];
[s,w] = mat_lyies(str,[handles.path_tmp FileName '.info']);

if ~(isequal(s,0))                  % An error as occured
    msg = ['Error getting file info (gdalinfo).', sprintf('\n\n'), 'Error message was:', sprintf('\n'), w];
    h = errordlg(msg,'Error');    movegui(h,'north');     return
end

file_info = dataread('file',[handles.path_tmp FileName '.info'],'%s','delimiter','\n','whitespace','');
Hdr = read_gdal_info(file_info);

if (strcmp(Hdr.Type,'Byte') | ~isempty(Hdr.CTable))     % We have a raster image
    warndlg('Reading ENVI or Erdas images are not yet supported (only DEMs).','Warning')
    return
end

set(handles.figure1,'pointer','watch')
att =  gdalread(fname,'-M','-C');
m = att.RasterXSize;            n = att.RasterYSize;
handles.head = att.GMT_hdr;     handles.head_orig = handles.head;

Hdr.LL_prj_xmin = num2str(att.Corners.LL(1));   Hdr.LL_prj_ymin = num2str(att.Corners.LL(2));
Hdr.UR_prj_xmax = num2str(att.Corners.UR(1));   Hdr.UR_prj_ymax = num2str(att.Corners.UR(2));
Hdr.Zmin = num2str(handles.head(5));            Hdr.Zmax = num2str(handles.head(6));

jump = min(round(m / 200), round(n / 200));
opt_P = ['-P' num2str(jump)];
Z =  gdalread(fname,'-U','-S',opt_P);      [m,n] = size(Z);
X = linspace(handles.head(1),handles.head(2),n);
Y = linspace(handles.head(3),handles.head(4),m);
colormap(handles.axes1,jet(256));           handles.icon_img = image(X,Y,Z);
set(handles.axes1,'YDir','normal')
Resize1(handles.axes1, handles.icon_img)
handles.fname = fname;
handles.image_type = 4;
set(handles.figure1,'pointer','watch')
guidata(hObject,handles)

% -----------------------------------------------------------------------------------------
function CropaGrid(obj,eventdata,h,opt)
% Select the region to be extracted from grid.
% If a fourth argument (OPT) is transmited, then the entire grid is read in
handles = guidata(obj);     % get handles

if (nargin == 3)
	x = get(h,'XData');     y = get(h,'YData');
	if (min(x) < handles.head_orig(1) | max(x) > handles.head_orig(2) | ...
            min(y) < handles.head_orig(3) | max(y) > handles.head_orig(4))
        errordlg('Selected region is partially outside the grid domain','Error')
        return
	end
	opt_R = ['-R' num2str(x(1),'%.15f') '/' num2str(x(4),'%.15f') '/' num2str(y(1),'%.15f') '/' num2str(y(2),'%.15f')];
else
    opt_R = ' ';
end

[Z,head,Hdr] = read_DEMs(handles,opt_R);

tmp.head = head;                        [m,n] = size(Z);
tmp.X = linspace(head(1),head(2),n);    tmp.Y = linspace(head(3),head(4),m);
tmp.name = 'Subsample from Overview';
tmp.was_int16 = handles.was_int16;
tmp.Nodata_int16 = handles.Nodata_int16;
new_window = mirone(Z,tmp);

% --------------------------------------------------------------------
function [Z,head,Hdr] = read_DEMs(handles,opt_R)
% This function loads grid files that are may contain DEMs or other grid (not images) types
% OPT is used when reading MOLA files, OPT = [x_min x_max y_min y_max n_cols n_rows grid_inc]

Hdr = [];
fname = handles.fname;
tipo = handles.type;

set(handles.figure1,'pointer','watch')
if (handles.ForceInsitu)    opt_I = '-I';   % Use only in desperate cases.
else                        opt_I = ' ';    end
handles.was_int16 = 0;      % To make sure that it wasnt left = 1 from a previous use.
if (strcmp(tipo,'GMT'))
    [Z,att] =  gdalread(fname,'-U',opt_I,opt_R,'-C');   Z = single(Z);
    head = att.GMT_hdr;
    % Need to recompute x_inc/y_inc due to a Gdal bug (Not in new versions, but it doesn't harm)
    head(8) = (head(2) - head(1)) / (size(Z,2) - ~head(7));
    head(9) = (head(4) - head(3)) / (size(Z,1) - ~head(7));
    
    if (strcmp(att.Band.DataType,'Int16'))
        handles.was_int16 = 1;
        handles.Nodata_int16 = att.Band.NoDataValue;
    end
elseif ( strcmp(tipo,'SRTM30') | strcmp(tipo,'SRTM3') | strcmp(tipo,'SRTM1') | strcmp(tipo,'MOLA') )    
    name_uncomp = [];
    if (~isempty(handles.name_uncomp))      % That is, if file was compressed
        fname = handles.name_uncomp;
        name_uncomp = handles.name_uncomp;
    end

    [Z,att] =  gdalread(fname,'-U','-C',opt_I,opt_R);   Z = single(Z);
    Z(Z <= single(att.Band.NoDataValue)) = NaN;
    head = att.GMT_hdr;
    handles.was_int16 = 1;      handles.Nodata_int16 = att.Band.NoDataValue;
elseif ( strcmp(tipo,'USGS_DEM') | strcmp(tipo,'GTOPO30') | strcmp(tipo,'DTED') | strcmp(tipo,'SDTS') | ...
        strcmp(tipo,'GeoTiff_DEM') | strcmp(tipo,'ArcAscii')  | strcmp(tipo,'ArcBinary') | ...
        strcmp(tipo,'GXF') | strcmp(tipo,'ENVI') | strcmp(tipo,'Erdas') | strcmp(tipo,'ESRI_hdr'))
    
    [Z,att] =  gdalread(fname,'-U','-C',opt_I,opt_R);   Z = single(Z);
    head = att.GMT_hdr;
    if (~isempty(att.Band.NoDataValue))   Z(Z <= single(att.Band.NoDataValue)) = NaN;    end
    
    if (strcmp(att.Band.DataType,'Int16'))
        handles.was_int16 = 1;
        handles.Nodata_int16 = att.Band.NoDataValue;
    end
end

handles.head = head;
set(handles.figure1,'pointer','arrow')
guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function [FileName,PathName] = put_or_get_file(handles,str1,str2,tipo)
% Use this function to select input or output filename
if (strcmp(tipo,'get'))
    cd(handles.last_dir)
    [FileName,PathName] = uigetfile(str1,str2);
    if (PathName ~= 0),         handles.last_dir = PathName;    end
elseif (strcmp(tipo,'put'))
    cd(handles.work_dir)
    [FileName,PathName] = uiputfile(str1,str2);
    if (PathName ~= 0),         handles.last_dir = PathName;    end
end
pause(0.01)
cd(handles.home_dir);       % allways go home to avoid troubles
guidata(handles.figure1,handles)

%--------------------------------------------
function Resize1(axHandle, imHandle)
% Resize figure containing a single axes object with a single image.

set(axHandle,'Units','normalized','Position',[0 0 1 1]) % Don't realy understand why, but I need this

imageWidth  = size(get(imHandle, 'CData'), 2);
imageHeight = size(get(imHandle, 'CData'), 1);

% Don't try to handle the degenerate case.
if (imageWidth * imageHeight == 0),    return;  end

% What are the screen dimensions
screenSize = get(0, 'ScreenSize');      screenWidth = screenSize(3);    screenHeight = screenSize(4);
if ((screenWidth <= 1) | (screenHeight <= 1))
    screenWidth = Inf;    screenHeight = Inf;
end

% For small images, compute the minimum side as 60% of largest of the screen dimensions
% Except in the case of croped images, where 512 is enough for the pushbuttons
% LeastImageSide = fix(max([screenWidth screenHeight] * 0.6));
LeastImageSide = 400;

if (imageWidth < LeastImageSide & imageHeight < LeastImageSide)   % Augment very small images
    while (imageWidth < LeastImageSide & imageHeight < LeastImageSide)
        imageWidth = imageWidth*1.05;  imageHeight = imageHeight*1.05;
    end
    imageWidth = fix(imageWidth);   imageHeight = fix(imageHeight);
    % Large aspect ratio figures may still need to have their size adjusted
    if (imageWidth < 512)
        while (imageWidth < LeastImageSide & imageHeight < screenHeight-50)
            imageWidth = imageWidth*1.05;  imageHeight = imageHeight*1.05;
        end
    end
end

axUnits = get(axHandle, 'Units');
set(axHandle, 'Units', 'pixels');
axPos = get(axHandle, 'Position');

figHandle = get(axHandle, 'Parent');
figUnits = get(figHandle, 'Units');
rootUnits = get(0, 'Units');
set(figHandle, 'Units', 'pixels');
set(0, 'Units', 'pixels');

figLeftBorder = 10;  % assume left figure decorations are 10 pixels
figRightBorder = 10;
figBottomBorder = 30;
figTopBorder = 50;
figTopBorder = figTopBorder + 30;  % scribe hack

minFigWidth = 128; % don't try to display a figure smaller than this.
minFigHeight = 128;

% What are the gutter sizes?
figPos = get(figHandle, 'Position');
gutterLeft = max(axPos(1) - 1, 0);
gutterRight = max(figPos(3) - (axPos(1) + axPos(3)) + 1, 0);
gutterBottom = max(axPos(2) - 1, 0);
gutterTop = max(figPos(4) - (axPos(2) + axPos(4)) + 1, 0);

scale = 100;    done = 0;
defAxesPos = get(0,'DefaultAxesPosition');
nonzeroGutters = (gutterLeft > 0);
while (~done)
    if (nonzeroGutters)
        gutterWidth = round((1 - defAxesPos(3)) * imageWidth / defAxesPos(3));
        gutterHeight = round((1 - defAxesPos(4)) * imageHeight / defAxesPos(4));
        newFigWidth = imageWidth + gutterWidth;
        newFigHeight = imageHeight + gutterHeight;
    else
        newFigWidth = imageWidth;        newFigHeight = imageHeight;
    end
    if (((newFigWidth + figLeftBorder + figRightBorder) > screenWidth) | ...
                ((newFigHeight + figBottomBorder + figTopBorder) > screenHeight))
        scale = 9 * scale / 10;
        imageWidth  = round(imageWidth * scale / 100);
        imageHeight = round(imageHeight * scale / 100);
    else
        done = 1;
    end
end

% ---------------------------------------------
h_Xlabel = get(axHandle,'Xlabel');  h_Ylabel = get(axHandle,'Ylabel');
units_save = get(h_Xlabel,'units');
set(h_Xlabel,'units','pixels');     set(h_Ylabel,'units','pixels');
Xlabel_pos = get(h_Xlabel,'pos');   Ylabel_pos = get(h_Ylabel,'pos');

if (abs(Ylabel_pos(1)) < 20)    % Stupid hack, but there is a bug somewhere
    Ylabel_pos(1) = 25;
end

setappdata(axHandle,'Backup_LabelPos',[Xlabel_pos Ylabel_pos])

y_margin = abs(Xlabel_pos(2))+get(h_Xlabel,'Margin') + 10;  % Devera conter a altura em pixeis do Xlabel
x_margin = abs(Ylabel_pos(1))+get(h_Ylabel,'Margin') + 10;  % Devera conter a largura em pixeis do Ylabel
if strcmp(get(axHandle,'Visible'),'off')               % No Labels, give only a 20 pixels margin to account for Status bar
    x_margin = 0;   y_margin = 0;
end
set(h_Xlabel,'units',units_save);     set(h_Ylabel,'units',units_save);

newFigWidth = max(newFigWidth, minFigWidth) + x_margin;
newFigHeight = max(newFigHeight, minFigHeight) + y_margin;

figPos(1) = max(1, figPos(1) - floor((newFigWidth - figPos(3))/2));
figPos(2) = max(1, figPos(2) - floor((newFigHeight - figPos(4))/2));
figPos(3) = newFigWidth;
figPos(4) = newFigHeight;

% Figure out where to place the axes object in the resized figure
gutterWidth = figPos(3) - imageWidth;
gutterHeight = figPos(4) - imageHeight;
gutterLeft = floor(gutterWidth/2) + x_margin/2 - 1;
gutterBottom = floor(gutterHeight/2) + y_margin/2 - 1;

axPos(1) = gutterLeft*0.8 + 1;  axPos(2) = gutterBottom - 5;
axPos(3) = imageWidth;      axPos(4) = imageHeight;

% Force the window to be in the "north" position. 75 is the height of the blue Bar + ...
figPos(2) = screenHeight - figPos(4) - 75;
set(figHandle, 'Position', figPos);     set(axHandle, 'Position', axPos);

% Restore the units
drawnow;  % necessary to work around HG bug   -SLE
set(figHandle, 'Units', figUnits);
set(axHandle, 'Units', axUnits);
set(0, 'Units', rootUnits);

% --------------------------------------------------------------------
function figure1_DeleteFcn(hObject, eventdata, handles)
% Before deleting the fig see if there was files left to remove
% (files that resulted from the uncompression)
try         % This bludy stupid sometimes (and it's realy when it pleases) doesn't know handles
	if (~isempty(handles.name_uncomp))
        try
            delete(handles.name_uncomp);    delete(handles.name_hdr);
        end
	end
end
delete(hObject)

% --------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')
    % Before deleting the fig see if there was files left to remove
    % (files that resulted from the uncompression)
    try         % This bludy stupid sometimes (and it's realy when it pleases) doesn't know handles
		if (~isempty(handles.name_uncomp))
            try
                delete(handles.name_uncomp);    delete(handles.name_hdr);
            end
		end
    end
	delete(hObject)
end

% --------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function h1 = overview_LayoutFcn()

h1 = figure(...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'IntegerHandle','off',...
'KeyPressFcn','overview(''figure1_KeyPressFcn'',gcbo,[],guidata(gcbo))',...
'MenuBar','none',...
'Name','overview',...
'NumberTitle','off',...
'Position',[520 560 244 240],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'DeleteFcn','overview(''figure1_DeleteFcn'',gcbo,[],guidata(gcbo))',...
'Visible','off',...
'Tag','figure1',...
'UserData',[]);

h2 = axes('Parent',h1,...
'Units','pixels',...
'Position',[1 1 240 240],...
'Tag','axes1',...
'Visible','off');

h4 = get(h2,'xlabel');
set(h4,'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[0.497916666666667 -0.0979166666666667 1.00005459937205],...
'VerticalAlignment','cap',...
'HandleVisibility','off',...
'Visible','off');

h5 = get(h2,'ylabel');
set(h5,'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[-0.11875 0.49375 1.00005459937205],...
'Rotation',90,...
'VerticalAlignment','bottom',...
'HandleVisibility','off',...
'Visible','off');

h7 = uimenu('Parent',h1,'Label','File Open','Tag','FileOpen');

h8 = uimenu('Parent',h7,...
'Callback','overview(''FileOpenGMTgrid_Callback'',gcbo,[],guidata(gcbo))',...
'Label','GMT grid',...
'Tag','FileOpenGMTgrid');

h9 = uimenu('Parent',h7,...
'Callback','overview(''FileOpenArcGrid_Callback'',gcbo,[],guidata(gcbo),''binary'')',...
'Label','Arc/Info Binary Grid',...
'Tag','FileOpenArcGrid');

h10 = uimenu('Parent',h7,...
'Callback','overview(''FileOpenArcGrid_Callback'',gcbo,[],guidata(gcbo),''ascii'')',...
'Label','Arc/Info Ascii Grid',...
'Tag','FileOpenArcGrid');

h11 = uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_Callback'',gcbo,[],guidata(gcbo),''DTED'')',...
'Label','DTED',...
'Tag','FileOpen_DEM');

uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_Callback'',gcbo,[],guidata(gcbo),''ESRI_hdr'')',...
'Label','ESRI BIL',...
'Tag','FileOpen_DEM');

h12 = uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_Callback'',gcbo,[],guidata(gcbo),''GTOPO30'')',...
'Label','GTOPO30',...
'Tag','FileOpen_DEM');

h13 = uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_Callback'',gcbo,[],guidata(gcbo),''GeoTiff_DEM'')',...
'Label','GeoTiff DEM',...
'Tag','FileOpen_DEM');

h14 = uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_Callback'',gcbo,[],guidata(gcbo),''GXF'')',...
'Label','GXF',...
'Tag','FileOpen_DEM');

h15 = uimenu('Parent',h7,...
'Callback','overview(''FileOpenMOLA_Callback'',gcbo,[],guidata(gcbo))',...
'Label','MOLA DEM',...
'Tag','FileOpenMOLA');

h16 = uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_Callback'',gcbo,[],guidata(gcbo),''SDTS'')',...
'Label','SDTS',...
'Tag','FileOpen_DEM');

h17 = uimenu('Parent',h7,...
'Callback','overview(''FileOpenSRTM_Callback'',gcbo,[],guidata(gcbo),''SRTM1'')',...
'Label','SRTM 1 arc sec',...
'Tag','FileOpenSRTM');

h18 = uimenu('Parent',h7,...
'Callback','overview(''FileOpenSRTM_Callback'',gcbo,[],guidata(gcbo),''SRTM3'')',...
'Label','SRTM 3 arc sec',...
'Tag','FileOpenSRTM');

h19 = uimenu('Parent',h7,...
'Callback','overview(''FileOpenSRTM_Callback'',gcbo,[],guidata(gcbo),''SRTM30'')',...
'Label','SRTM30',...
'Tag','FileOpenSRTM');

h20 = uimenu('Parent',h7,...
'Callback','overview(''FileOpen_DEM_Callback'',gcbo,[],guidata(gcbo),''USGS_DEM'')',...
'Label','USGS DEM',...
'Tag','FileOpen_DEM');

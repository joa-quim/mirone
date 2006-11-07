function varargout = mirone(varargin)
%   MIRONE, by itself, creates a window bar from which you load a lot of grid/images formats
%   MIRONE(GMT_GRID) opens the grid (must be a GMT grid) and displays it
%   H = MIRONE returns the handle to a new mirone window
%
%   mirone('CALLBACK',hObject,eventData,handles,...) calls the local
%   function named CALLBACK in mirone.m with the given input arguments.

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

if (nargin >= 4 && ischar(varargin{1}))
    gui_Callback = str2func(varargin{1});
    feval(gui_Callback,varargin{2:end});
    %if (nargout),   varargout{1} = get(0,'CurrentFigure');     end
    if (nargout),   varargout{1} = varargin{4}.figure1;    end
else
    h = mirone_OpeningFcn(varargin{:});
    if (nargout),   varargout{1} = h;       end
end

function hObject = mirone_OpeningFcn(varargin)
% varargin   command line arguments to mirone (see VARARGIN)
[hObject,version7,IamCompiled] = mirone_uis;
handles = guihandles(hObject);

% PRAGMA SECTION (It's far far from clear when files must be declared here)
%#function uigetfolder_standalone mapproject_m grdproject_m coordinate_system surface_m
%#function nearneighbor_m cpt2cmap grdfilter_m grdgradient_m grdsample_m grdtrack_m grdtrend_m 
%#function grdutils scaleto8
%----- These are in mag
%#function bpass3d inv3d nskew rtp3d syn3d igrf_m
%----- These are for okada tsunamis
%#function range_change swan tsun2 mansinha_m
%----- These are for image
%#function grayto8 grayto16 grayxform imfilter_mex imhistc imlincombc parityscan uintlutc ordf
%#function imreconstructmex applylutc bwboundariesmex bwlabel1 bwlabel2
%----- These are for griddata_j
%#function qhullmx
%----- These are in utils
%#function tabpanelfcn degree2dms dms2degree dec2deg dec_year ivan_the_terrible ddewhite string_token
%#function test_dms text_read double2ascii save_seismicity jd2date trimpatch
%#function run_and_report_error guess_file shading_mat getline_j frame3D histos_seis
%----- These is for ecran
%#function smoothing_param
%----- These is for write_gmt_script
%#function draw_scale time_stamp pscoast_options_Mir paint_option w_option
%----- Those are ..., the hell with explanations for what I don't realy understand. They are needed, that's all.
%#function gmtlist_m country_select read_isf choosebox MagBarCode listbox_message add_poles animate_seismicity
%#function get_polygon rot_euler datums telha_m find_clusters fft_stuff select_cols uistack_j 
%#function patch_meca ui_edit_patch_special bands_list multibandread_j imscroll_j transform_fun
%#function load_defFilters mltable_j iptcheckinput resampsep intmax wgifc

global home_dir;    home_dir = pwd;

% Reposition the window on screen
if (version7),  pos = get(hObject,'Pos');   set(hObject,'Pos',[pos(1:3) 1]);    end
movegui(hObject,'north');
set(0,'CurrentFigure',hObject)      % Due to a R2006a incredible BUG

% The addpath command cannot be compiled, so the paths bellow have be added manually
% to the matlab path. However, in the distribution version the next line(s) must be
% uncommented in order that the user doesn't have to know that those need to be added.
% addpath([home_dir filesep 'src_figs']);
% addpath([home_dir filesep 'lib_mex']);
% addpath([home_dir filesep 'utils']);
% addpath([home_dir filesep 'mag']);

handles.home_dir = home_dir;
handles.version7 = version7;% % If == 1 => R14 or latter
handles.IamCompiled = IamCompiled; % If == 1 than we know that we are dealing with a compiled (V3) version
handles.DefLineThick = 1;   % Default line thickness (overwriten by mirone_pref)
handles.DefLineColor = 'k'; % Default line color (overwriten by mirone_pref)
handles.DefineMeasureUnit = 'k'; % Default measure units to kilometrs (overwriten by mirone_pref)
handles.grdname = [];       % Contains the name of the current (if it's the case) gmt grid
handles.h_MBplot  = [];     % Handles to multi-beam tracks
handles.nTrack = 0;         % Counter of the number of MB tracks imported
handles.hist = 0;           % 0 -> linear colormap; 1 -> histogram equalized
handles.origFig = [];       % To store the original image copy
handles.DEM_name = [];      % DEM file name for eventual GMT grid conversion
handles.grdformat = 0;      % Flag to signal the format of the imported gmt grid (defaults to dafault)
handles.image_type = 0;     % Image type. 1->grd; 2-> trivial (jpg,png,bmp,etc...); 3->GeoTIFF; 4->DEMs; 5->
handles.computed_grid = 0;  % However, matrices with a gmt header will have this == 1, so that they can be saved
handles.no_file = 1;        % 0 means a grid is loaded and 1 that it is not (to test when icons can be pushed)
handles.BgMap = 0;          % 1 means a Background image map was loaded
handles.geog = 1;           % By default grids are assumed to be in geographical coordinates
handles.swathRatio = 3;     % Default swath width / water depth ratio for multibeam planing
handles.grdMaxSize = 20971520;   % I use this for limiting the grid size that is stored in RAM (20 Mb)
handles.firstMBtrack = 1;   % Used for knowing whether to display or not the MB planing info message in "start planing"
handles.EarthRad = 6371;    % Authalic radius
handles.is_grid_reg = 0;    % Due to the unknown type of DEMs loaded with GDAL, I set this as default.
handles.calling_figure = [];% When Mirone is called with an array in argument, this contains the calling figure handle
handles.maregraphs_count = 0; % Counter of maregraphs (tsunami modeling)
handles.Illumin_type = 0;   % 0 -> no illumination; 1 -> grdgradient; 2 -> grdgradient Lambertian; 4 -> Lambertian;
handles.out_in_NewWindow = 1;   % 1 -> every new computed grid is outputed in a new window. 0 -> grids are saved on disk
handles.zoom_state = 0;     % Flag to signal if zoom state is to be re-activated (0 means no re-activation)
handles.bg_color = [1 1 1]; % Backgoround color used when grid values == NaN
handles.which_cont = [];    % To store the contour levels (used to not draw repeated contours)
handles.have_nans = [];     % Used to know if the grids have NaNs
handles.is_draped = 0;      % Used to know if the image comes from draping
handles.was_int16 = 0;      % Keep track of short int imported grids
handles.saveAsInt16 = 0;    % Save gmt grids in the format #8
handles.Nodata_int16 = [];  % To store Nodata of short int grids
handles.ForceInsitu = 0;    % Use "insitu" grid importing (transposition)
handles.DefineEllipsoide = [6378137, 0, 1/298.2572235630];    % Defaults to WGS-84
handles.path_tmp = [home_dir filesep 'tmp' filesep];
handles.path_data = [home_dir filesep 'data' filesep];
handles.last_directories = {handles.path_tmp; home_dir};    % Let it have something existent
handles.grd_img = [];       % Stupid name which will hold the image handle
handles.imgName = [];       % To hold the full path of loaded images
handles.firstIllum = 1;     % First illumination will use the displayed image which may have been IP
handles.oldSize = get(hObject,'Pos');  
if (handles.oldSize(4) == 0),      handles.oldSize(4) = 1;    end
setappdata(hObject,'VirginSize',handles.oldSize)   % To use inside resizetrue

try     % A file named mirone_pref.mat (under data) contains the preferences, read them from it
    load([handles.path_data 'mirone_pref.mat']);
    handles.geog = geog;
    handles.grdMaxSize = grdMaxSize * 2^20;                     % 2^20 = 1 Mb
    handles.swathRatio = swathRatio;
    handles.last_directories = directory_list;
    handles.DefLineThick = str2double(DefLineThick{1}(1));
    handles.out_in_NewWindow = out_in_NewWindow;
    handles.saveAsInt16 = saveAsInt16;
    % Decode the line color string into the corresponding char (e.g. k,w, etc...)
    if (strcmp(DefLineColor{1},'Black')),   handles.DefLineColor = 'k';
    else                                    handles.DefLineColor = lower(DefLineColor{1}(1));   end
    % Decode the Measure unites into a char code (e.g. n, k, m, u from {'nautical miles' 'kilometers' 'meters' 'user'})
    handles.DefineMeasureUnit = DefineMeasureUnit{1}(1);
    handles.DefineEllipsoide = DefineEllipsoide_params;     % Set the default ellipsoide parameters (a,b,f)
end

j = false(1,length(handles.last_directories));          % vector for eventual cleaning non-existing dirs
for (i = 1:length(handles.last_directories))            % Check that all dirs in last_directories exist
    try         cd(handles.last_directories{i});        % NOTE. I don't use 'exist' anymore because
    catch       j(i) = 1;       cd(home_dir);           % the stupid compiler allways return something > 0
    end
end
handles.last_directories(j) = [];                       % clean non-existing directories
cd(home_dir);               % Need to come back home because it was probably somewere out there

if (isempty(handles.last_directories))              % Don't never let it be empty
    handles.last_directories = {handles.path_tmp; home_dir};    % Let it have something existent
end
handles.work_dir = handles.last_directories{1};     %
handles.last_dir = handles.last_directories{1};     % Initialize last_dir to work_dir

setappdata(hObject,'DefLineThick',handles.DefLineThick)     % Save this for accessability in draw_funs
setappdata(hObject,'DefLineColor',handles.DefLineColor)     % Save this for accessability in draw_funs
setappdata(hObject,'swathRatio',handles.swathRatio);        % I need this for getline_mb
setappdata(hObject,'ValidGrid',0)        % Flag to signal draw_funs that "Crop grid" isn't possible

% Change the MeasureDistance label to the selected (in prefs) unites
set(handles.ToolsMeasureDist,'Label',['Distance in ' handles.DefineMeasureUnit])
set(handles.ImageDrape,'Enable','off');             % Set the Drape option to it's default value (off)

% Detect in which mode Mirone was called
grd_fname_in = 0;   grd_data_in = 0;    grd_data_interfero_in = 0;  grd_data_deformation_in = 0;
n_argin = length(varargin);
if ~isempty(varargin)
    if (n_argin == 1 && ischar(varargin{1}))            % Called with a gmt grid name as argument
        grd_fname = varargin{1};        grd_fname_in = 1;
    elseif ( isa(varargin{1},'uint8') || isa(varargin{1},'logical') )
        % Called with an image as argument and optionaly an struct header (& geog, name, cmap optional fields)
        dims = size(varargin{1});
        if ( n_argin == 2 && isstruct(varargin{2}) )            % An image with coordinates
            tmp = varargin{2};
            handles.head = tmp.head;        X = tmp.X;      Y = tmp.Y;
            handles.image_type = 3;         axis_t = 'xy';
            if (isfield(tmp,'geog')),       handles.geog = tmp.geog;
            else                            handles.geog = guessGeog(handles.head(1:4));
            end
            if (isfield(tmp,'cmap')),       set(handles.figure1,'Colormap',tmp.cmap);   end
            if (isfield(tmp,'name')),       win_name = tmp.name;
            else                            win_name = 'Mirone';
            end
        else
            X = [];         Y = [];         win_name = 'Cropped Image';
            handles.image_type = 2;         handles.geog = 0;       axis_t = 'off';
            handles.head = [1 dims(2) 1 dims(1) 0 255 0 1 1];   % Fake a grid reg GMT header
            pal = getappdata(0,'CropedColormap');               % See if we have a colormap to use here
            if (~isempty(pal)),     set(handles.figure1,'Colormap',pal);    rmappdata(0,'CropedColormap');  end
            setappdata(hObject,'Croped','yes');             % ???
        end
        handles = show_image(handles,win_name,X,Y,varargin{1},0,axis_t,handles.head(7),1);
        grid_info(handles,[],'iminfo',varargin{1});         % Contruct a info string
        if (isa(varargin{1},'logical'))
            set(handles.grd_img,'CDataMapping','scaled');   set(handles.figure1,'ColorMap',gray(255));
        end        
    elseif (length(varargin) == 2 && isnumeric(varargin{1}) && isstruct(varargin{2}))
        % A matrix. Treat it as if it is a gmt grid. No error testing on the grid head descriptor
        grd_data_in = 1;
        Z = varargin{1};    tmp = varargin{2};
        handles.have_nans = grdutils(Z,'-N');
        head = tmp.head;    X = tmp.X;  Y = tmp.Y;
        try     name = tmp.name;    % All calls should transmit a name, but if they do not ...
        catch   name = [];  end
        if (isfield(tmp,'was_int16'))
            handles.was_int16 = tmp.was_int16;      handles.Nodata_int16 = tmp.Nodata_int16;
        end
        clear tmp;
    elseif ( length(varargin) == 4 && isnumeric(varargin{1}) && isstruct(varargin{2}) && ...
            strcmp(varargin{3},'Deformation') && ishandle(varargin{4}) )
        % A matrix. Treat it as if it'is a gmt grid. No error testing on the grid head descriptor
        % Note: this is a special case of the situation above that will be used to identify this figure
        % as an Okada deformtion data (via its Name). This info is searched by the tsunami modeling option
        grd_data_deformation_in = 1;
        Z = varargin{1};    tmp = varargin{2};
        head = tmp.head;    X = tmp.X;  Y = tmp.Y;  clear tmp;
        handles.calling_figure = varargin{4};
    elseif ( length(varargin) == 4 && isnumeric(varargin{1}) && isstruct(varargin{2}) && ...
            strcmp(varargin{3},'Interfero') && isnumeric(varargin{4}) )
        % A matrix input containing an interfeogram with cdo == varargin{4}
        grd_data_interfero_in = 1;
        Z = varargin{1};    tmp = varargin{2};      cdo = varargin{4};
        head = tmp.head;    X = tmp.X;  Y = tmp.Y;  clear tmp;
    end
else        % Called with no arguments
    set(hObject,'Name','Mirone')
end

% The following IF cases deals only with cases where a grid or a gmt grid name was given in argument
if (grd_fname_in || grd_data_in || grd_data_interfero_in || grd_data_deformation_in)
    if (grd_fname_in)
        win_name = grd_fname;
        handles.grdname = grd_fname;
        [X,Y,Z,head] = grdread_m(grd_fname,'single');
        handles.have_nans = grdutils(Z,'-N');
    elseif (grd_data_in)    % That is, grid as a matrix (the header was already parsed above)
        handles.grdname = '';
        handles.computed_grid = 1;          % Signal that this is a computed gmt grid
        if (~isempty(name)),    win_name = name;
        else                    win_name = 'Mirone';
        end
    elseif (grd_data_deformation_in)        % A deformation (Okada) array
        win_name = 'Okada deformation';
        handles.grdname = '';
        handles.computed_grid = 1;          % Signal that this is a computed gmt grid
    elseif (grd_data_interfero_in)
        win_name = 'Interferogram';
        handles.grdname = '';
        handles.computed_grid = 1;          % Signal that this is a computed gmt grid
    end
    handles.image_type = 1;
    handles.geog = guessGeog(head(1:4));
    handles.head = head;            % Save header info
    if (grd_data_interfero_in)      % Interferogram grid
        load([handles.path_data 'gmt_other_palettes.mat'],'circular');
        pal = circular;
        zz = uint8(abs(rem(double(Z),cdo)/cdo)*255);
    else                            % Not an interferogram grid
        zz = scaleto8(Z);       pal = jet(256);
    end
    % If grid size is not to big I'll store it
    if (numel(Z)*4 <= handles.grdMaxSize)
        if (~isa(Z,'single')),  setappdata(hObject,'dem_z',single(Z));
        else                    setappdata(hObject,'dem_z',Z);  end
        setappdata(hObject,'dem_x',X);  setappdata(hObject,'dem_y',Y);
        setappdata(hObject,'GMThead',head);
    end
    set(hObject, 'Units', 'pixels');    setappdata(hObject,'Zmin_max',[head(5) head(6)])
    aux_funs('colormap_bg',handles,Z,pal);
    handles = show_image(handles,win_name,X,Y,zz,1,'xy',head(7));
end

% Set some accelerators
set(findobj(hObject,'Tag','ToolsMBplaningEdit'), 'Accelerator','e');
set(findobj(hObject,'Tag','DrawEditLine'), 'Accelerator','l');

%Find out which gmt version is beeing used. 
info = getappdata(0,'gmt_version');     % See if the info is already there.
if (isempty(info))
    info = test_gmt;
    setappdata(0,'gmt_version',info);   % Save it so that the next time a new mirone window is opened
end                                     % (but within the same session) we don't have to run test_gmt.

if (info.version(1) ~= '0')             % we have a full GMT installation, find out which coastlines are installed.
    if (info.full ~= 'y')
        set([handles.DatasetsCoastLineFull handles.DatasetsPBFull handles.DatasetsRiversFull], 'Enable','off')
    end
    if (info.high ~= 'y')
        set([handles.DatasetsCoastLineHigh handles.DatasetsPBHigh handles.DatasetsRiversHigh], 'Enable','off')
    end    
else                                    % we don't have a full GMT, so use our own minimalist init files
    info = test_gmt(['GMT_USERDIR=' home_dir filesep 'gmt_userdir']);  % Set the GMT_USERDIR env and recheck for coastline files
    if (info.full ~= 'y')
        set([handles.DatasetsCoastLineFull handles.DatasetsPBFull handles.DatasetsRiversFull], 'Enable','off')
    end
    if (info.high ~= 'y')
        set([handles.DatasetsCoastLineHigh handles.DatasetsPBHigh handles.DatasetsRiversHigh], 'Enable','off')
    end
    if (info.intermediate ~= 'y')
        set([handles.DatasetsCoastLineInterm handles.DatasetsPBInterm handles.DatasetsRiversInterm], 'Enable','off')
    end
    if (info.low ~= 'y')
        set([handles.DatasetsCoastLineLow handles.DatasetsPBLow handles.DatasetsRiversLow], 'Enable','off')
    end
    if (info.crude ~= 'y')
        set([handles.DatasetsCoastLineCrude handles.DatasetsPBCrude handles.DatasetsRiversCrude], 'Enable','off')
    end
end

guidata(hObject, handles);
set(hObject,'Visible','on','HandleVisibility','callback');  limpa(handles);
setappdata(hObject,'IAmAMirone',1);         % Use this appdata to identify Mirone figures
%setappdata(handles.axes1,'ProjWKT',geogWKT) % The Geog WGS84 string in WKT format

% --------------------------------------------------------------------------------------------------
function SetAxesNumericType(handles,eventdata,hand1)
if (nargin == 3)        % This happens when this function was called trough mirone by an
    handles = hand1;    % an external program. Then hand1 must contain the gcf handles.
end
% Save original X & Y labels in appdata for easear access when we want to change them
h_axes = findobj(handles.figure1,'Type','Axes');
setappdata(handles.figure1,'XTickOrig',get(h_axes,'XTickLabel'))
setappdata(handles.figure1,'YTickOrig',get(h_axes,'YTickLabel'))
set(h_axes, 'FontSize', 9)     % Make this the default
if (handles.geog),      setappdata(handles.figure1,'LabelFormatType','DegDec')
else                    setappdata(handles.figure1,'LabelFormatType','NotGeog')
end
Parent = get(h_axes,'Parent');

figure(handles.figure1)                 % Allways bring main figure to stack's top because of uicontextmenu
cmenu_axes = uicontextmenu('Parent',Parent);    % Ned to set Parent for when Mirone is called from M_GMT
set(h_axes, 'UIContextMenu', cmenu_axes);
if (handles.geog)
    uimenu(cmenu_axes, 'Label', 'Label Format -> DD.xx', 'Callback', 'draw_funs([],''ChangeAxesLabels'',''ToDegDec'')');
    uimenu(cmenu_axes, 'Label', 'Label Format -> DD MM', 'Callback', 'draw_funs([],''ChangeAxesLabels'',''ToDegMin'')');
    uimenu(cmenu_axes, 'Label', 'Label Format -> DD MM.xx', 'Callback', 'draw_funs([],''ChangeAxesLabels'',''ToDegMinDec'')');
    uimenu(cmenu_axes, 'Label', 'Label Format -> DD MM SS', 'Callback', 'draw_funs([],''ChangeAxesLabels'',''ToDegMinSec'')');
    uimenu(cmenu_axes, 'Label', 'Label Format -> DD MM SS.x', 'Callback', 'draw_funs([],''ChangeAxesLabels'',''ToDegMinSecDec'')');
end
itemFS = uimenu(cmenu_axes, 'Label', 'Label Font Size', 'Separator','on');
uimenu(itemFS, 'Label', '7   pt', 'Callback', 'set(gca, ''FontSize'', 7)');
uimenu(itemFS, 'Label', '8   pt', 'Callback', 'set(gca, ''FontSize'', 8)');
uimenu(itemFS, 'Label', '9   pt', 'Callback', 'set(gca, ''FontSize'', 9)');
uimenu(itemFS, 'Label', '10 pt', 'Callback', 'set(gca, ''FontSize'', 10)');
uimenu(cmenu_axes, 'Label', 'Grid on/off', 'Callback', 'grid', 'Separator','on');
if (handles.ValidGrid)
    uimenu(cmenu_axes, 'Label', 'Pixel mode on/off', 'Callback', {@PixMode_callback,handles.figure1},'Separator','on');
end

h_warning = findobj('Type','figure','Name','Warning');
if ~isempty(h_warning),     figure(h_warning);   end     % If a warning message exists, bring it forward

% --------------------------------------------------------------------
function PixMode_callback(hObject, eventdata, hFig)
% Inside each grid cell, which is a pixel in the screen, display only the grid node value
if (strcmp(get(hObject,'Checked'),'off'))
    set(hObject,'Checked','on');    setappdata(hFig,'PixelMode',1)
else
    set(hObject,'Checked','off');   setappdata(hFig,'PixelMode',0)
end

% --------------------------------------------------------------------
function ImageCrop_Callback(hObject, eventdata, handles, opt, opt2, opt3)
% OPT is either a handle to a line that may be a rectangle/polygon, OR, if empty
%   calls rubberbandbox to do a interactive croping (called by "Crop Grid")
% OPT2 is a string to direct this function to different operations that
%   apply to the grid and update the image.
% OPT3 contains the interpolation method when OPT2 == 'FillGaps' ('cubic', 'linear' or 'sea')
% Note: I won't make the "Drape" option active in the cropped window
% NAO TA ACABADA. FALTA PROGRAMAR OS CASE ILUMINACAO == 2,3,4   (PT)

set(handles.figure1,'pointer','watch')
crop_pol = 0;   % Defaults to croping from a rectangle
if (nargin < 5),    opt2 = [];      end
if (nargin < 6),    opt3 = [];      end
if ~isempty(opt)        % OPT must be a rectangle/polygon handle (the rect may serve many purposes)
    x = get(opt,'XData');   y = get(opt,'YData');
    if ~( (x(1) == x(end)) && (y(1) == y(end)) && length(x) == 5 && ...
            (x(1) == x(2)) && (x(3) == x(4)) && (y(1) == y(4)) && (y(2) == y(3)) )
        xp(1) = min(x);     xp(2) = max(x);
        yp(1) = min(y);     yp(2) = max(y);
        rect_crop = [xp(1) yp(1) (xp(2) - xp(1)) (yp(2) - yp(1))];
        crop_pol = 1;       % Flag that we are croping from a polygon
    else
        rect_crop = [x(1) y(1) (x(3)-x(2)) (y(2)-y(1))];
    end
    if isempty(opt2)                            % Just pure Image croping
        Z_rect = get(findobj(handles.figure1,'Type','image'),'CData');
        limits = getappdata(handles.figure1,'ThisImageLims');
        I = cropimg(limits(1:2),limits(3:4),Z_rect,rect_crop);
        [m,n] = size(I);
    elseif (strcmp(opt2,'CropaWithCoords'))     % Crop Image with coordinates
        Z_rect = get(findobj(handles.figure1,'Type','image'),'CData');
        [I,r_c] = cropimg([handles.head(1:2)],[handles.head(3:4)],Z_rect,rect_crop,'out_grid');
        [m,n] = size(I);
    else                    % Extract the sub-grid inside the rectangle/polygon
        [X,Y,Z,head] = load_grd(handles);
        if isempty(Z),  set(handles.figure1,'pointer','arrow');    return;     end;    % An error message was already issued
        [Z_rect,r_c] = cropimg([head(1) head(2)],[head(3) head(4)],Z,rect_crop,'out_grid');
        if (crop_pol)
            [zzz] = grdutils(Z_rect,'-L');  z_min = zzz(1);     clear zzz;
            if (strcmp(opt2,'CropaGrid_pure'))
                defAns = {num2str(z_min)};
                resp  = inputdlg({'Enter outside polygon value'},'Choose out value',[1 30],defAns);    pause(0.01)
                if isempty(resp);    set(handles.figure1,'pointer','arrow');    return;     end
            elseif (strcmp(opt2,'ROI_SetConst'))    % Set the polygon interiour to cte
                resp  = inputdlg({'Enter new grid value'},'Replace with cte value',[1 30]);    pause(0.01)
                if isempty(resp);    set(handles.figure1,'pointer','arrow');    return;     end
            end
            x_lim = [min(x) max(x)];        y_lim = [min(y) max(y)];
            mask = img_fun('roipoly_j',x_lim,y_lim,double(Z_rect),x,y);
            if (strcmp(opt2,'CropaGrid_pure'))
                Z_rect(~mask) = single(str2double(resp));
            elseif (strcmp(opt2,'ROI_SetConst'))
                Z_rect(mask) = single(str2double(resp));    % Set the mask values to const
                Z(r_c(1):r_c(2),r_c(3):r_c(4)) = Z_rect;
                if (isnan(str2double(resp))),  handles.have_nans = 1;  end
            elseif (strcmp(opt2,'ROI_MedianFilter'))
                [Z,Z_rect,handles] = roi_filtering(handles, Z, head, Z_rect, r_c, mask);
            elseif (strcmp(opt2,'CropaGrid_histo'))
                Z_rect(~mask) = single(NaN);
            else
                warndlg('Unknown case in ImageCrop','Warning');     return
            end
        end
        [m,n] = size(Z_rect);
    end
else                    % Interactive croping (either Grid or Image)
    if (strcmp(opt2,'CropaGrid'))   % Arrive here when called by "Extract Region"
        [X,Y,Z,head] = load_grd(handles);
        if isempty(Z),  set(handles.figure1,'pointer','arrow');    return;     end;
        [p1,p2] = rubberbandbox;
        x0 = min(p1(1),p2(1));      y0 = min(p1(2),p2(2));
        dx = abs(p2(1)-p1(1));      dy = abs(p2(2)-p1(2));
        [Z_rect,r_c] = cropimg([head(1) head(2)],[head(3) head(4)],Z,[x0 y0 dx dy],'out_grid');
        X = (head(1) + (r_c(3)-1)*head(8)):head(8):(head(1) + (r_c(4)-1)*head(8));
        Y = (head(3) + (r_c(1)-1)*head(9)):head(9):(head(3) + (r_c(2)-1)*head(9));
        head(1) = X(1);     head(2) = X(end);       head(3) = Y(1);     head(4) = Y(end);
        tit = 'Grid cuted by Mirone';      % Have to change this to reflect the old title
        GRD_save_or_display(handles,X,Y,Z_rect,head,tit,'Croped grid')
        return
    else            % Just a image crop op
        I = cropimg;     [m,n] = size(I);
    end
end

first_nans = 0;     pal = [];
if (isempty(opt2) || strcmp(opt2,'CropaWithCoords'))   % Just pure Image croping
    if (m < 2 || n < 2),  set(handles.figure1,'pointer','arrow');    return;     end;    % Image too small. Probably a user bad mouse control
    if (handles.BgMap == 1)     % Image needs to be re-fliped ud (it was previously flipud)
        I = flipdim(I,1);
    end
    if (strcmp(get(handles.axes1,'Ydir'),'normal')),    I = flipdim(I,1);    end
    if (ndims(I) == 2)
        pal = get(handles.figure1, 'Colormap');
        if (length(pal) == 64), pal = jet(256);     end     % Risky - This is a patch for "Find Clusters"
        setappdata(0,'CropedColormap',pal);         % indexed image, so I need to save it's colormap
    end
    set(handles.figure1,'pointer','arrow');
    if (isempty(opt2))
        mirone(I);
    else
        head(2) = handles.head(1) + (r_c(4)-1)*handles.head(8);     head(1) = handles.head(1) + (r_c(3)-1)*handles.head(8);
        head(4) = handles.head(3) + (r_c(2)-1)*handles.head(9);     head(3) = handles.head(3) + (r_c(1)-1)*handles.head(9);
        head(5) = 0;            head(6) = 255;     head(7) = 0;     head(8:9) = handles.head(8:9);  tmp.name = 'Croped Image';
        tmp.head = head;        tmp.geog = handles.geog;            tmp.X = head(1:2);    tmp.Y = head(3:4);
        if (~isempty(pal)),     tmp.cmap = pal;     end
        mirone(flipdim(I,1),tmp);
    end
    return;                 % We are done. BYE BYE.
elseif ( strncmp(opt2(1:min(length(opt2),9)),'CropaGrid',9) )       % Do the operatio indicated in opt2(11:end) & return
    curr_opt = opt2(11:end);
    if (~strcmp(curr_opt,'pure'))           % We will need those for all other options
        head(2) = head(1) + (r_c(4)-1)*head(8);         head(1) = head(1) + (r_c(3)-1)*head(8);
        head(4) = head(3) + (r_c(2)-1)*head(9);         head(3) = head(3) + (r_c(1)-1)*head(9);
        head(5) = double(min(Z_rect(:)));               head(6) = double(max(Z_rect(:)));
        to_func.Z = Z_rect;                             to_func.head = head;
    end
    if (strcmp(curr_opt,'pure'))            % PURE means pure CropaGrid
        X = (head(1) + (r_c(3)-1)*head(8)):head(8):(head(1) + (r_c(4)-1)*head(8));
        Y = (head(3) + (r_c(1)-1)*head(9)):head(9):(head(3) + (r_c(2)-1)*head(9));
        head(1) = X(1);     head(2) = X(end);       head(3) = Y(1);     head(4) = Y(end);
        tit = 'Grid cuted by Mirone';       % Have to change this to reflect the old title
        GRD_save_or_display(handles,X,Y,Z_rect,head,tit,'Croped grid')
    elseif (strcmp(curr_opt,'histo'))       % HISTO means compute histogram inside the selected rect area
        GridToolsHistogram_Callback([], [], guidata(handles.figure1), to_func);
    elseif (strcmp(curr_opt,'power'))       % POWER means compute log10 power spectrum
        GridToolsSectrum_Callback([], [], guidata(handles.figure1), 'Power', to_func)
    elseif (strcmp(curr_opt,'autocorr'))    % AUTOCORR means compute the autocorrelation
        GridToolsSectrum_Callback([], [], guidata(handles.figure1), 'Autocorr', to_func)
    elseif (strcmp(curr_opt,'fftTools'))    % FFTTOOLS means call the fft_stuff
        GridToolsSectrum_Callback([], [], guidata(handles.figure1), 'Allopts', to_func)
    end
    return
elseif (strcmp(opt2,'FillGaps'))
    if ~any(isnan(Z_rect(:)))    % No gaps
        set(handles.figure1,'pointer','arrow');    warndlg('Selected area does not have any voids (NaNs)','Warning');   return;
    else
        X = (head(1) + (r_c(3)-1)*head(8)):head(8):(head(1) + (r_c(4)-1)*head(8));
        Y = (head(3) + (r_c(1)-1)*head(9)):head(9):(head(3) + (r_c(2)-1)*head(9));
        if (~isempty(opt3) && strcmp(opt3,'surface'))
            opt_R = ['-R' num2str(X(1),10) '/' num2str(X(end),10) '/' num2str(Y(1),10) '/' num2str(Y(end),10)];
            opt_I = ['-I' num2str(head(8),10) '/' num2str(head(9),10)];
        end
        Z_rect = double(Z_rect);      % It has to be
        aa = isnan(Z_rect(:));
        [X,Y] = meshgrid(X,Y);
        ZZ = Z_rect(:);     ZZ(aa) = [];
        XX = X(:);          XX(aa) = [];
        YY = Y(:);          YY(aa) = [];
        if (~isempty(opt3))
            switch opt3
                case 'surface', Z_rect = surface_m(XX,YY,ZZ,opt_R,opt_I,'-T.25');
                case 'cubic',   Z_rect = griddata_j(XX,YY,ZZ,X,Y,'cubic');
                case 'linear',  Z_rect = griddata_j(XX,YY,ZZ,X,Y,'linear');
                case 'sea',     Z_rect(aa) = 0;  clear X XX Y YY ZZ;
            end
        else
            Z_rect = surface_m(XX,YY,ZZ,opt_R,opt_I,'-T.25','-v');
        end
        clear X XX Y YY ZZ;
        Z_rect = single(Z_rect);      % Revert it
        Z(r_c(1):r_c(2),r_c(3):r_c(4)) = Z_rect;
    end
elseif (strcmp(opt2,'SplineSmooth'))
    X = (head(1) + (r_c(3)-1)*head(8)):head(8):(head(1) + (r_c(4)-1)*head(8));
    Y = (head(3) + (r_c(1)-1)*head(9)):head(9):(head(3) + (r_c(2)-1)*head(9));
    Z_rect = double(Z_rect);      % It has to be
    [pp p_guess] = spl_fun('csaps',{Y(1:min(m,10)),X(1:min(n,10))},Z_rect(1:min(m,10),1:min(n,10)));% Get a good estimate of p
    prompt = {'Enter smoothing p paramer'};     dlg_title = 'Smoothing parameter input';
    defAns = {num2str(p_guess{1},12)};     resp  = inputdlg(prompt,dlg_title,[1 38],defAns);
    pause(0.01)
    if isempty(resp);    set(handles.figure1,'pointer','arrow');    return;     end
    pp = spl_fun('csaps',{Y,X},Z_rect,str2double(resp{1}));
    Z_rect = single(spl_fun('fnval',pp,{Y,X}));    clear pp;
    Z(r_c(1):r_c(2),r_c(3):r_c(4)) = Z_rect;
elseif (strcmp(opt2,'MedianFilter'))
    [Z,Z_rect,handles] = roi_filtering(handles, Z, head, Z_rect, r_c, 'rect', 'no');
elseif (strcmp(opt2,'SetConst'))        % Replace grid values inside rect by a cte value
    prompt = {'Enter new grid value'};     dlg_title = 'Replace with cte value';
    resp  = inputdlg(prompt,dlg_title,[1 30]);    pause(0.01)
    if isempty(resp);    set(handles.figure1,'pointer','arrow');    return;     end
    Z_rect = repmat(single(str2double(resp)),m,n);
    Z(r_c(1):r_c(2),r_c(3):r_c(4)) = Z_rect;
    if (~handles.have_nans && isnan(str2double(resp)))       % See if we have new NaNs
        handles.have_nans = 1;      first_nans = 1;
    elseif (handles.have_nans && ~isnan(str2double(resp)))   % Check that old NaNs had not been erased
        handles.have_nans = grdutils(Z_rect,'-N');
    end
end

if ~isempty(opt2)       % Here we have to update the image in the processed region
    X = (head(1) + (r_c(3)-1)*head(8)):head(8):(head(1) + (r_c(4)-1)*head(8));
    Y = (head(3) + (r_c(1)-1)*head(9)):head(9):(head(3) + (r_c(2)-1)*head(9));
    h_img = findobj(handles.axes1,'Type','image');
    [zzz] = grdutils(Z,'-L');       z_min = zzz(1);     z_max = zzz(2);     clear zzz;      img = [];
    if ( (abs(z_min - head(5)) > 1e-5 || abs(z_max - head(6)) > 1e-5) && handles.Illumin_type == 0 )
        img = scaleto8(Z);              % Z_MIN or Z_MAX have changed. Need to recompute image (but only if no illumin)
    end
    if (handles.Illumin_type == 0)
        z_int = uint8(round( ((double(Z_rect) - z_min) / (z_max - z_min))*255 ));
    elseif ( handles.Illumin_type == 1 || handles.Illumin_type == 2 || handles.Illumin_type == 4 )
        luz = getappdata(handles.figure1,'Luz');
        z_int = uint8(round( ((double(Z_rect) - z_min) / (z_max - z_min))*255 ));
        z_int = ind2rgb8(z_int,get(handles.figure1,'Colormap'));    % z_int is now 3D
        if (handles.Illumin_type == 1)
            head_tmp = [X(1) X(end) Y(1) Y(end) head(5:9)];
            opt_N = ['-Nt1/' num2str(handles.grad_sigma) '/' num2str(handles.grad_offset)];
            if (handles.geog),  R = grdgradient_m(Z_rect,head_tmp,'-M',['-A' num2str(luz.azim)],opt_N);
            else                R = grdgradient_m(Z_rect,head_tmp,['-A' num2str(luz.azim)],opt_N); end
            z_int = shading_mat(z_int,R,'no_scale');% and now it is illuminated
        elseif (handles.Illumin_type == 2)          % grdgradient lambertian
            head_tmp = [X(1) X(end) Y(1) Y(end) head(5:9)];
            R = grdgradient_m(Z_rect,head_tmp,['-Es' num2str(luz.azim) '/' num2str(luz.elev)]);
            z_int = shading_mat(z_int,R,'no_scale');
        else                                        % Lambertian lighting illum
           R = grdgradient_m(Z,head,['-E' num2str(luz.azim) '/' num2str(luz.elev) '/' num2str(luz.ambient),...
                   '/' num2str(luz.diffuse) '/' num2str(luz.specular) '/' num2str(luz.shine)]);
            z_int = shading_mat(z_int,R,'no_scale');% and now it is illuminated
        end
    end
    if (isempty(img))  img = get(h_img,'CData');    end     % That is, if img was not recomputed, get from screen 
    handles.img_back = img(r_c(1):r_c(2),r_c(3):r_c(4),:);  % For the undo op
    img(r_c(1):r_c(2),r_c(3):r_c(4),:) = z_int;
    clear z_int Z_rect R;
    if (first_nans)       % We have NaNs for the first time. Adjust the colormap
        aux_funs('colormap_bg',handles,Z,get(handles.figure1,'Colormap'));
    end
    set(h_img,'CData',img)

    head(5) = z_min;     head(6) = z_max;
    handles.computed_grid = 1;              handles.head = head;    handles.origFig = img;
    setappdata(handles.figure1,'dem_z',Z);  setappdata(handles.figure1,'GMThead',head);
    setappdata(handles.figure1,'Zmin_max',[head(5) head(6)])
    guidata(handles.figure1, handles);      set(handles.figure1,'pointer','arrow')
end

% Experimental UNDO
if (~isempty(opt2) && strcmp(opt2,'MedianFilter')),       set_undo(handles,opt);  end
if (~isempty(opt2) && strcmp(opt2,'ROI_MedianFilter')),   set_undo(handles,opt);  end

% --------------------------------------------------------------------
function set_undo(handles,h)
% Experimental UNDO  that works only with the "Median Filter" option
cmenuHand = get(h,'UIContextMenu');
cb_undo = {@do_undo,handles,h,cmenuHand};
uimenu(cmenuHand, 'Label', 'Undo', 'Separator','on', 'Callback', cb_undo);

% -----------------------------------------------------------------------------------------
function do_undo(obj,eventdata,handles,h,img)
[X,Y,Z,head] = load_grd(handles);   % Experimental. No testing for error in loading
Z(handles.r_c(1):handles.r_c(2),handles.r_c(3):handles.r_c(4)) = handles.Z_back;
[zzz] = grdutils(Z,'-L');       z_min = zzz(1);     z_max = zzz(2);     clear zzz;
head(5) = z_min;                head(6) = z_max;
setappdata(handles.figure1,'dem_z',Z);    setappdata(handles.figure1,'GMThead',head);
h_img = findobj(handles.axes1,'Type','image');
handles.origFig = get(h_img,'CData');
handles.origFig(handles.r_c(1):handles.r_c(2),handles.r_c(3):handles.r_c(4),:) = handles.img_back;
set(h_img,'CData',handles.origFig)
guidata(handles.figure1,handles)
ui_u = findobj(get(h,'UIContextMenu'),'Label', 'Undo');
set(ui_u,'Visible','off')       % If I delete it we get an error

% --------------------------------------------------------------------
function ImageFlip_Callback(hObject, eventdata, handles, opt)
% OPT == 'LR'   Flips the image left-right
% OPT == 'UD'   Flips the image up-down
if (handles.no_file == 1),    return;      end
h_img = findobj(handles.figure1,'Type','image');    img = get(h_img,'CData');
if length(size(img)) == 3               % RGB image
    if strcmp(opt,'LR'),    img = flipdim(img,2);
    else                    img = flipdim(img,1);   end   % Flip UD
else
    if strcmp(opt,'LR'),    img = fliplr(img);
    else                    img = flipud(img);      end
end
set(h_img,'CData', img);
if ~isempty(handles.origFig),   handles.origFig = img;  end
guidata(hObject, handles);

% --------------------------------------------------------------------
function ImageResetOrigImg_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),    return;      end
h_img = findobj(handles.axes1,'Type','image');
try         % In some cases (e.g. histograms, countries) img may not exist
	set(h_img,'CData', handles.origFig);
    set(handles.figure1,'ColorMap',handles.origCmap)
	handles.hist = 0;   handles.Illumin_type = 0;   handles.firstIllum = 1;
    handles.ValidGrid = handles.ValidGrid_orig;     handles.was_int16 = handles.was_int16_orig;
    handles.computed_grid = handles.computed_grid_orig;
    set(handles.ImgHist,'checked','off');
	guidata(hObject, handles);
end

% --------------------------------------------------------------------
function ImageHistEqualize_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),     return;      end

zz = get(handles.grd_img,'CData');
if strcmp(get(hObject,'checked'),'off')     % Then equalize
    if (ndims(zz) == 3);
        zz = cvlib_mex('color',zz,'rgb2hsv');
        zz(:,:,3) = img_fun('histeq_j',zz(:,:,3));
        J = cvlib_mex('color',zz,'hsv2rgb');
    else
        J = img_fun('histeq_j',zz);
    end
    set(handles.grd_img,'CData', J);    handles.hist = 1;       set(hObject,'checked','on')
elseif ~isempty(handles.origFig)            % Then de-equalize
    set(handles.grd_img,'CData', handles.origFig);    handles.hist = 0;
    handles.Illumin_type = 0;       set(hObject,'checked','off');
else
    msgbox('Sorry. To save memory I didn''t make an image copy. You''ll have to reload the file again.','Warning')
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function ImageHistEqualizeGrid_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),      return;      end
if (~handles.ValidGrid),        return,     end
set(handles.figure1,'pointer','watch')

if strcmp(get(hObject,'checked'),'off')     % Then equalize
    [X,Y,Z,head] = load_grd(handles);
    if isempty(Z),   set(handles.figure1,'pointer','watch');   return;     end;
    out = grdutils(Z,'-S');     % Get mean and std
    handles.cur_pal = get(handles.figure1, 'Colormap');
    new_pal = cdf2pal(head(5),head(6),out(1),out(2),handles.cur_pal);
    aux_funs('colormap_bg',handles,Z,new_pal);
    handles.hist_grid = 1;          set(hObject,'checked','on')
else            % Then de-equalize
    set(handles.figure1,'Colormap',handles.cur_pal);      handles.hist_grid = 0;
    handles.Illumin_type = 0;       set(hObject,'checked','off');
end
guidata(hObject, handles);  set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function ImageShowPalette_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),     return;      end
h = findobj(handles.figure1,'Type','image');
if (ndims(get(h,'CData')) == 3)
    msgbox('True color images do not use color palettes.','Warning');    return
end
cmap = get(handles.figure1,'Colormap');
if (handles.image_type == 1)
    dz = (handles.head(6) - handles.head(5)) / 10;
    show_palette([1 length(cmap)], [handles.head(5) handles.head(6)], dz, cmap)
elseif (handles.image_type > 2 || handles.image_type <= 5)
    mm =  getappdata(handles.figure1,'Zmin_max');
    dz = diff(mm) / 10;
    show_palette([1 length(cmap)], [mm(1) mm(2)], dz, cmap)
else
    msgbox('Sorry, Palette not available for this Image type.')
end

% --------------------------------------------------------------------
function PanZoom_Callback(hObject,event,handles,opt)
if (handles.no_file == 1),    return;      end
if (strcmp(opt,'zoom'))
	if strcmp(get(hObject,'State'),'on')
        zoom_j('on');    h = findobj(handles.figure1,'Tag','Mao');
        if (strcmp(get(h,'State'),'on'))
            set(h,'State','off');   pan('off');
        end
	else
        zoom_j('off');
	end
else        % Pan case
	if strcmp(get(hObject,'State'),'on')
        pan('on');    h = findobj(handles.figure1,'Tag','Zoom');
        if (strcmp(get(h,'State'),'on'))
            set(h,'State','off');   zoom_j('off');
        end
	else
        pan('off');
	end
end

% --------------------------------------------------------------------
function zoom_state(handles, state)
% Sets the zoom sate to off, or reset it to on if ...
switch state
    case 'off_yes'          % Set zoom permanently off
        zoom_j('off');
        set(findobj(handles.figure1,'Tag','Zoom'),'State','off');
        handles.zoom_state = 0;     guidata(handles.figure1,handles)
    case 'maybe_off'        % If zoom was active, keep trace of it
        zoom_j('off');        h = findobj(handles.figure1,'Tag','Zoom');
        if (strcmp(get(h,'State'),'on')),   handles.zoom_state = 1;
        else                                handles.zoom_state = 0;    end
        set(h,'State','off');       guidata(handles.figure1,handles)
    case 'maybe_on'         % Check if zoom has to be re-activated
        handles = guidata(handles.figure1);  % Need to get the updated handles
        if (handles.zoom_state)
            zoom_j('on');     set(findobj(handles.figure1,'Tag','Zoom'),'State','on');
        end
end

% --------------------------------------------------------------------
function FilePreferences_Callback(hObject, eventdata, handles)
% mirone_pref deals with sizes in MB, but here we need it in bytes
max = fix(handles.grdMaxSize / (2^20));
out = mirone_pref(handles.geog, max, handles.swathRatio, handles.ForceInsitu);
if isempty(out),    return;     end
handles.geog = out.geog;        handles.swathRatio = out.swathRatio;
handles.grdMaxSize = out.grdMaxSize * 2^20;
handles.work_dir = out.last_dir;        handles.last_dir = out.last_dir;
handles.DefLineThick = out.DefLineThick;            % Default line thickness
handles.DefLineColor = out.DefLineColor;            % Default line color
handles.DefineMeasureUnit = out.DefineMeasureUnit;  % Default measure unit (.'n', .'k', .'m', .'u')
handles.DefineEllipsoide = out.DefineEllipsoide;    % Default ellipsoide
handles.out_in_NewWindow = out.out_in_NewWindow;    % Choose what to do with new computed grids
handles.saveAsInt16 = out.saveAsInt16;
handles.ForceInsitu = out.ForceInsitu;
setappdata(handles.figure1,'swathRatio',handles.swathRatio);    % I need this in getline_mb
setappdata(handles.figure1,'DefLineThick',handles.DefLineThick)     % Save this for accessability in draw_funs
setappdata(handles.figure1,'DefLineColor',handles.DefLineColor)     % Save this for accessability in draw_funs
guidata(hObject, handles);

% --------------------------------------------------------------------
function FileNewBgFrame_Callback(hObject, eventdata, handles, region, imSize)
% Create a empty window with a frame selected in bg_region
% However, if REGION was transmited, it is assumed to have  [x_min x_max y_min y_max is_geog]
if (nargin == 3)
    region = bg_region;     % region contains [x_min x_max y_min y_max is_geog]
    if isempty(region),     return;     end     % User gave up
end

X = region(1:2);      Y = region(3:4);      handles.geog = guessGeog(region(1:4));
handles.head = [X Y 0 255 0];

scrsz = get(0,'ScreenSize');         % Get screen size
dx = X(2) - X(1);   dy = Y(2) - Y(1);
aspect = dy / dx;
nx = round(scrsz(3)*.75);   ny = round(nx * aspect);
if (ny > scrsz(4) - 30)
    ny = scrsz(4) - 30;
    nx = round(ny / aspect);
end
handles.head(8) = diff(X) / (nx - 1);    handles.head(9) = diff(Y) / (ny - 1);
Z = repmat(uint8(255),ny,nx);
pal = repmat(handles.bg_color,256,1);    set(handles.figure1,'Colormap',pal);
handles.image_type = 20;     handles.grdname = [];
if (nargin <= 4),   imSize = [];    end
show_image(handles,'Mirone Base Map',X,Y,Z,0,'xy',0,imSize);

% --------------------------------------------------------------------
function FileNewBgMap_Callback(hObject, eventdata, handles)
% The output of bg_regionMapTilled is a structure with the following fields:
% out.img   -> the image                % out.X         -> image's X limits
% out.Y     -> image's Y limits         % out.imgName   -> image's full path and name
out = bg_region_map_tilled;
if isempty(out),    return;     end     % User gave up loading the fig tille
handles.imgName = out.imgName;          handles.grdname = [];
handles.geog = 1;                       handles.image_type = 2;
handles.head(1:2) = out.X;              handles.head(3:4) = out.Y;      handles.head(5:7) = [0 255 0];
handles.head(8) = diff(out.X) / (size(out.img,2)-1);    handles.head(9) = diff(out.Y) / (size(out.img,1)-1);
handles = show_image(handles,out.imgName,out.X,out.Y,out.img,0,'xy',0,1);
handles.BgMap = 1;      guidata(hObject,handles);

% --------------------------------------------------------------------
function FileNewEmpty_Callback(hObject, eventdata, handles)
setappdata(0,'parent_gcf',handles.figure1);
h = mirone;
set(findobj(h,'Tag','ImageDrape'),'Enable','on')        % Set the Drape option to 'on' in the New window 

% --------------------------------------------------------------------
function FileExport_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),    return;      end
set(handles.figure1,'pointer','watch');
h = findobj('Type','uicontrol');    set(h,'Visible','off')
handsStBar = getappdata(handles.figure1,'CoordsStBar');
set(handsStBar,'Visible','off');
cd(handles.work_dir)
filemenufcn(handles.figure1,'FileExport')
cd(handles.home_dir);       % allways go home to avoid troubles
set(h,'Visible','on');  set(handsStBar,'Visible','on');
set(handles.figure1,'pointer','arrow');

% --------------------------------------------------------------------
function SaveGenericImage(handles)
if (handles.no_file == 1),    return;      end
str1 = {'*.jpg', 'JPEG image (*.jpg)'; ...
    '*.bmp', 'Windows Bitmap (*.bmp)'; ...
    '*.hdf', 'Hieralchical Data Format (*.hdf)'; ...
    '*.gif', 'GIF image (*.gif)'; ...
    '*.pcx', 'Windows Paintbrush (*.pcx)'; ...
    '*.png', 'Portable Network Graphics(*.png)'; ...
    '*.ras', 'SUN rasterfile (*.ras)'; ...
    '*.raw', 'Raw RGB format (*.raw)'; ...
    '*.tif', 'Tagged Image File (*.tif)'; ...
    '*.xwd', 'X Windows Dump (*.xwd)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select image format','put');
if isequal(FileName,0);     return;     end

[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if isempty(EXT) && ~isempty(FileName)
    msgbox('Sorry, but you have to give the filename extention (consequence of Matlab bugs!)','Error'); return
end
set(handles.figure1,'pointer','watch')
img = get(handles.grd_img,'CData');
if (isa(img,'logical')),      img = uint8(bitshift(uint16(img),8));   end
if (strcmp(get(handles.axes1,'Ydir'),'normal')),   img = flipdim(img,1);     end
if (strcmpi(EXT,'.jpg') || strcmpi(EXT,'.jpeg'))
    if isa(img, 'uint8' && ndims(img) == 2)
        img = ind2rgb8(img,get(handles.figure1,'Colormap'));
    end
    imwrite(img,[PathName FileName],'Quality',100);
elseif (strcmpi(EXT,'.raw'))
    img = get(handles.grd_img,'CData');
    if isa(img, 'uint8')        % Not yet programed for other image types
        if (ndims(img) == 2),   pix = ind2rgb8(img,get(handles.figure1,'Colormap'));
        else                    pix = img;      end
        fid = fopen([PathName FileName],'wb');
        [nl,nc,np] = size(pix);                 l = 1;      pix1 = repmat(uint8(0),nl*nc*3,1);
        if (strcmp(get(handles.axes1,'Ydir'),'normal'));    m = nl:-1:1;
        else                                                m = 1:nl;        end
        for i=m
            for j=1:nc
                for k=1:3;  pix1(l) = pix(i,j,k);    l = l + 1;     end
            end
        end
        fwrite(fid,pix1,'uint8');
        %fwrite(fid,permute(pix, [3 2 1]),'uint8');
        fclose(fid);
    else
        msgbox('Sorry: RAW format can only be used with byte type images.','Warning'); return
    end
elseif strcmpi(EXT,'.gif')              % Non existent in < R14
    if (ndims(img) == 3)
        errordlg('GIF format does not support saving RGB images','Error');
        set(handles.figure1,'pointer','arrow');     return
    end
    writegif(img,get(handles.figure1,'Colormap'),[PathName FileName]);    
else        % All other image formats
    if isa(img, 'uint8')
        try
            imwrite(img,get(handles.figure1,'Colormap'),[PathName FileName]);
        catch       % For example RGB images canot be saved as pcx
            errordlg('Error writing image file','Error');   lasterror
        end
    else
        try        imwrite(img,[PathName FileName]);
        catch      errordlg('Error writing image file','Error');    lasterror
        end
    end
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function FileSaveGMTgrid_Callback(hObject, eventdata, handles,opt)
% Save internaly computed grids and GDAL recognized DEM grids into GMT grd grids
if (aux_funs('msg_dlg',14,handles));     return;      end
if (nargin == 3),   opt = [];   end

[X,Y,Z,head] = load_grd(handles);
if isempty(Z),   return;     end;    % An error message was already issued
if (~isempty(opt) && strcmp(opt,'Surfer'))
    cdf_format = '=6';
    tit = ' ';      % Have to change this to reflect the old title
    txt1 = 'Surfer 6 binary grid (*.grd,*.GRD)';    txt2 = 'Select output Surfer 6 grid';
else           % Internaly computed grid    
    if (handles.was_int16 && handles.saveAsInt16)
        cdf_format = '=ns';
    else
        cdf_format = '=nf';
    end
    tit = 'Grid computed inside Mirone';
    txt1 = 'netCDF grid format (*.grd,*.GRD)';      txt2 = 'Select output GMT grid';
end

[FileName,PathName] = put_or_get_file(handles,{'*.grd;*.GRD',txt1; '*.*', 'All Files (*.*)'},txt2,'put');
if isequal(FileName,0);     return;     end

set(handles.figure1,'pointer','watch')
[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if isempty(EXT),    f_name = [PathName FNAME '.grd' cdf_format];
else                f_name = [PathName FNAME EXT cdf_format];       end

% If it was a grid imported by gdal, uppdate the title
if (isappdata(handles.axes1,'DatumProjInfo'))
    DPI = getappdata(handles.axes1,'DatumProjInfo');
    tit = ['Projection: ' DPI.projection ' Datum: ' DPI.datum];
    if (length(tit) > 80),      tit = tit(1:80);    end     % (1:80) otherwise it BOOMs
end
grdwrite_m(Z,head,f_name,tit);      set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function File_img2GMT_RGBgrids_Callback(hObject, eventdata, handles, opt1, opt2, opt3)
% Save image as a triplet of gmt grids - R,G,B.
% OPT1 == 'image' || == [] Capture only the image and not graphical elements
%   (lines, symbols, etc...). The grids nrow & ncol is the same as the image
%   number of lines and pixels.
% OPT1 == 'screen' Does a screen capture that includes all the graphical elements
%   that may have been drawn (lines, symbols, etc...). On the other hand I still
%   don't know how to control the number of lines and pixels.
% OPT2 == 'geotiff' Uses the choice made by OPT1 and writes a GeoTIFF file
% OPT3 == fname. It is used by the write gmt script routine to capture the image and
%   write the image as a triplet of gmt grids with name stem = OPT3. NOTE: must have OPT1 = 'image'
if (handles.no_file == 1),      aux_funs('msg_dlg',1,handles);     return;      end

if    (nargin == 3),   opt1 = 'image';     opt2 = [];  opt3 = [];
elseif(nargin == 4),   opt2 = [];          opt3 = [];
elseif(nargin == 5),   opt3 = [];          end

if (isempty(opt2))
    str1 = {'*.grd;*.GRD','netCDF int2 grid format (*.grd,*.GRD)'; '*.*', 'All Files (*.*)'};
    str2 = 'Select output GMT grid';
else
    str1 = {'*.tiff', 'GeoTiff (*.tiff)'};    str2 = 'Select Geotiff file name';
end
if (isempty(opt3))
    [FileName,PathName] = put_or_get_file(handles,str1,str2,'put');
    if isequal(FileName,0);     return;     end
else        % It means the output file name was transmited in input
    [PathName,FileName] = fileparts(opt3);
    PathName = [PathName filesep];
end

set(handles.figure1,'pointer','watch')
fmt = '=ns';

[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if isempty(EXT)
    f_name_r = [PathName FNAME ['_r.grd' fmt]];     f_name_g = [PathName FNAME ['_g.grd' fmt]];
    f_name_b = [PathName FNAME ['_b.grd' fmt]];     f_name_gray = [PathName FNAME ['_gray.grd' fmt]];
else
    f_name_r = [PathName FNAME '_r' EXT fmt];       f_name_g = [PathName FNAME '_g' EXT fmt];
    f_name_b = [PathName FNAME '_b' EXT fmt];       f_name_gray = [PathName FNAME '_gray' EXT fmt];
end

tit = ' ';
if (isappdata(handles.axes1,'DatumProjInfo'))
    DPI = getappdata(handles.axes1,'DatumProjInfo');
    tit = ['Projection: ' DPI.projection ' Datum: ' DPI.datum];
end
D = handles.head(1:7);

if (handles.image_type == 1)                % GMT grid
    tit = 'GMT grid converted to 3 RGB grids';
elseif (handles.image_type == 2)            % Generic formats
    tit = handles.imgName;
end
if (length(tit) > 80),  tit = tit(1:80);    end     % (1:80) otherwise it BOOMs

if (strcmp(opt1,'image'))
    img = get(findobj(handles.axes1,'Type','image'),'CData');   % Get image
    if (ndims(img) == 2 && isempty(opt2))
        img = ind2rgb8(img,get(handles.figure1,'Colormap'));    % Need this because image is indexed
    end
    %if (isempty(opt2)),    img = flipdim(img,1);      end      % That, is no GeoTIFF, so image needs to be ud fliped
    if (~strcmp(get(handles.axes1,'Ydir'),'normal') || ~isempty(opt2)),    img = flipdim(img,1);      end
elseif (strcmp(opt1,'screen'))                      % Do a screen capture
    resizetrue(handles.figure1,'screen_capture');
    img = screen_capture(handles.figure1,'lixo_screen.jpg');    % Only the .jpg really matters
    resizetrue(handles.figure1,'after_screen_capture');
    img = aux_funs('strip_bg_color',handles,img);
    if (isempty(opt2)),     img = flipdim(img,1);       end      % That, is no GeoTIFF. Image needs to be ud fliped
end

if (isempty(opt2))              % That is, no GeoTIFF
    grdwrite_m(single(img(:,:,1)),D,f_name_r,tit)
    grdwrite_m(single(img(:,:,2)),D,f_name_g,tit)
    grdwrite_m(single(img(:,:,3)),D,f_name_b,tit)
else                            % GeoTIFF
%     if isempty(EXT),    FileName = [FileName '.tiff'];  end
%     write_geotiff(handles,img,D,[PathName FileName])
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function FileSaveENCOMgrid_Callback(hObject, eventdata, handles)
% Save memory resident grids into the Encom grid format
if (aux_funs('msg_dlg',14,handles));     return;      end

txt1 = 'Encom grid format (*.grd,*.GRD)';   txt2 = 'Select output Encom grid';
[FileName,PathName] = put_or_get_file(handles,{'*.grd;*.GRD',txt1; '*.*', 'All Files (*.*)'},txt2,'put');
if isequal(FileName,0);     return;     end

[X,Y,Z,head,m,n] = load_grd(handles);
if isempty(Z),   return;     end;    % An error message was already issued

set(handles.figure1,'pointer','watch')
[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if isempty(EXT),    f_name = [PathName FNAME '.grd'];
else                f_name = [PathName FNAME EXT];       end
fid = fopen(f_name,'wb');
ID = ['Model Vision Grid' repmat(' ',1,80-17)];     fwrite(fid,ID,'80*char');   % Fk stupid ML oblyged me to this
ID = ['Mirone Grid (Title unknown here)' repmat(' ',1,80-32)];     fwrite(fid,ID,'80*char');
fwrite(fid,'  NO_REF','char');
fwrite(fid,'GRIDFPT ZNIL','char');
fwrite(fid,single(-2e16),'float32');
fwrite(fid,'ROWS','char');      fwrite(fid,m,'float32');
fwrite(fid,'COLS','char');      fwrite(fid,n,'float32');
fwrite(fid,'XORG','char');      fwrite(fid,single(head(1)),'float32');
fwrite(fid,'YORG','char');      fwrite(fid,single(head(3)),'float32');
fwrite(fid,'DX  ','char');      fwrite(fid,single(head(8)),'float32');
fwrite(fid,'DY  ','char');      fwrite(fid,single(head(9)),'float32');
fwrite(fid,'DEGR','char');      fwrite(fid,single(0),'float32');
if (handles.have_nans),      Z(isnan(Z)) = -2e16;    end
Z = (rot90(Z,1));               Z = flipud(Z);      % I cannot do better than this manip
fwrite(fid,Z,'float32');        fclose(fid);
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function FilePrintSetup_Callback(hObject, eventdata, handles)
print -dsetup

% --------------------------------------------------------------------
function FilePrint_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),    return;      end
h = findobj('Type','uicontrol');
set(h,'Visible','off')              % We don't want to print the buttons
handsStBar = getappdata(handles.figure1,'CoordsStBar');
set(handsStBar,'Visible','off');
if (ispc),      print -v
else            print;  end
set(h,'Visible','on');  set(handsStBar,'Visible','on');

% --------------------------------------------------------------------
function ImageInfo_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),     return;      end
if (handles.image_type == 1 || handles.image_type == 4)
    [X,Y,Z,head] = load_grd(handles);
    grid_info(handles,X,Y,head,Z)
else
    grid_info(handles,[],[],[],[])
end

% --------------------------------------------------------------------
function ExtractProfile_Callback(hObject, eventdata, handles, opt)
if (aux_funs('msg_dlg',100,handles));     return;      end  % Test only for no_file.
point_int = 0;                                  % Default to "profile" interpolation
if (nargin == 4),   point_int = 1;     end      % Interpolate at the line vertex only
[X,Y,Z,head] = load_grd(handles,'silent');
if (isempty(Z) && (handles.image_type == 1 || handles.image_type == 4))     % Grid not in memory error
    errordlg('Grid was not on memory. Increase "Grid max size" and start over again.','ERROR'); return
elseif (isempty(Z) && ndims(get(handles.grd_img,'CData')) == 2)
    Z = get(handles.grd_img,'CData');
    img_lims = getappdata(handles.figure1,'ThisImageLims'); % Get limits and correct them for the pix reg problem
    x_inc = (img_lims(2)-img_lims(1)) / size(Z,2);      y_inc = (img_lims(4)-img_lims(3)) / size(Z,1);
    img_lims = img_lims + [x_inc -x_inc y_inc -y_inc]/2;    % Remember that the Image is ALWAYS pix reg
    X = linspace(img_lims(1),img_lims(2),size(Z,2));    Y = linspace(img_lims(3),img_lims(4),size(Z,1));
    head = [X(1) X(end) Y(1) Y(end) double(min(Z(:))) double(max(Z(:))) 1 X(2)-X(1) Y(2)-Y(1)];
elseif (isempty(Z))
    errordlg('Extracting profile of a RGB image is not suported.','ERROR');     return
end

zoom_state(handles,'maybe_off')
if ~isempty(getappdata(handles.figure1,'TrackThisLine'))
    hand = getappdata(handles.figure1,'TrackThisLine');
    xp = get(hand,'Xdata');     yp = get(hand,'Ydata');
    rmappdata(handles.figure1,'TrackThisLine')      % Clear it so that the next time it may work when called interactivelly
else
    [xp,yp] = getline_j(handles.figure1);
end
n_nodes = length(xp);
if (n_nodes < 2),       zoom_state(handles,'maybe_on');      return;     end
if (~point_int)         % Profile interp
	dx = X(2) - X(1);   dy = Y(2) - Y(1);
	% Here I don't realy know what is the good increment for interpolation, so I'll
	% interpolate at the half grid spacing for each dimension.
	% Construct the vectors with the points where to interpolate the profile
	xx = [];     yy = [];
	for i=1:n_nodes-1
        n_int = round( max( abs(xp(i+1)-xp(i))/(dx/2), abs(yp(i+1)-yp(i))/(dy/2) ) );         % find ...
        xx = [xx linspace(xp(i),xp(i+1),n_int)];     yy = [yy linspace(yp(i),yp(i+1),n_int)];  % at nodes, values are repeated
	end
else                    % Interpolation at line vetex
    xx = xp;    yy = yp;
end

% Interpolate
if ~isempty(getappdata(handles.figure1,'dem_x'))        % Grid is in memory
    if (~getappdata(handles.figure1,'PixelMode'))       % Interpolation mode
        zz = grdtrack_m(Z,head,[xx' yy'],'-Z')';        % It uses less memory (and non doubles)
        %zz = interp2(X,Y,double(Z),xx,yy,'*cubic');
    else                                                % NEARNEIGBOR mode
        [rows,cols] = size(Z);
        rp = aux_funs('axes2pix',rows, get(handles.grd_img,'YData'),yy);
        cp = aux_funs('axes2pix',cols, get(handles.grd_img,'XData'),xx);
        r = min(rows, max(1, round(rp)));   c = min(cols, max(1, round(cp)));
        rc = (c - 1) * rows + r;
        zz = double(Z(rc));
    end
else                                % grid was loaded here (big according to preferences), so interp linearly
    zz = bi_linear(X,Y,Z,xx,yy);
end

zoom_state(handles,'maybe_on')
if (~point_int)         % Disply profile in ecran
    [pato,name,ext] = fileparts(get(handles.figure1,'Name'));
    ecran('Image',xx,yy,zz,['Track from ' name ext])
else                    % Save result on file
    draw_funs([],'save_xyz',[xx(:) yy(:) zz(:)])
end

% --------------------------------------------------------------------
function FileOpenOverview_Callback(hObject, eventdata, handles)
tmp.home_dir = handles.home_dir;        tmp.last_dir = handles.last_dir;
tmp.work_dir = handles.work_dir;        tmp.ForceInsitu = handles.ForceInsitu;
overview(tmp)

% --------------------------------------------------------------------
function FileOpenBSB_Callback(hObject, eventdata, handles)
str1 = {'*.kap;*.KAP;*.nos;*NOS', 'BSB Nautical Chart (*.kap,*.KAP,*.nos,*.NOS)'; '*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select BSB Nautical Chart File','get');
if isequal(FileName,0);     return;     end

set(handles.figure1,'pointer','watch')
att = gdalread([PathName FileName],'-M','-C');
if (~isempty(att.ProjectionRef)),   Z = gdalread([PathName FileName],'-U','-C');
else                                Z = gdalread([PathName FileName],'-C');
end
handles.head = att.GMT_hdr;
X = handles.head(1:2);      Y = handles.head(3:4);
if (~isempty(att.Band(1).ColorMap)),    pal = att.Band(1).ColorMap.CMap(:,1:3);
else                                    pal = jet(256);     end

handles.geog = guessGeog([X Y]);
handles.image_type = 3;
handles.grdname = [];   set(handles.figure1,'Colormap',pal);
aux_funs('cleanGRDappdata',handles);        % Remove eventual grid stuff variables from appdata
show_image(handles,[PathName FileName],X,Y,Z,0,'xy',1);
grid_info(handles,att,'gdal')           % Construct a info message

% --------------------------------------------------------------------
function FileOpen_ENVI_Erdas_Callback(hObject, eventdata, handles, opt)
% This function reads both ENVI or Erdas files. Furthermore, based on the file byte
% type it guesses if we are dealing with a typical grid file (in which case it is
% treated like a native gmt grid) or a raster image file.

str1 = {'*.img;*.IMG', [opt ' (*.img,*.IMG)']; '*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,['Select ' opt ' File'],'get');
if isequal(FileName,0);     return;     end

if (handles.ForceInsitu),   opt_I = '-I';   % Use only in desperate cases.
else                        opt_I = ' ';    end
handles.was_int16 = 0;      % To make sure that it wasnt left = 1 from a previous use.

att = gdalread([PathName FileName],'-M','-C');
if ((att.RasterXSize * att.RasterYSize * 4) > handles.grdMaxSize)
    if ( strcmp(yes_or_no('title','Warning'),'Yes')),  return;     end      % Advise accepted
end

set(handles.figure1,'pointer','watch')
if (strcmp(att.Band(1).DataType,'Byte') || ~isempty(att.Band(1).ColorMap))     % We have a raster image
    if (strcmp(att.Band(1).DataType,'Byte'))
        zz = gdalread([PathName FileName],'-U');
    else
        zz = gdalread([PathName FileName],'-S','-U');
    end
    head = att.GMT_hdr;             [m,n,k] = size(zz);
    X = [head(1) head(2)];          Y = [head(3) head(4)];
    handles.geog = guessGeog([X Y]);
    if strcmp(att.ColorInterp,'gray'),          pal = gray(256);
    elseif strcmp(att.ColorInterp,'Palette')
        if (~isempty(att.Band(1).ColorMap)),    pal = att.Band(1).ColorMap.CMap(:,1:3);
        else                                    pal = jet(256);     end
    else                                        pal = jet(256);
    end
    aux_funs('cleanGRDappdata',handles);        % Remove eventual grid stuff variables from appdata
    handles.image_type = 3;
    ValidGrid = 0;                              % Signal that grid opps are not allowed
    set(handles.figure1,'Colormap',pal)
else
    Z = gdalread([PathName FileName],'-U','-C',opt_I);       Z = single(Z);
    head = att.GMT_hdr;
    if (~isempty(att.Band(1).NoDataValue)),  Z(Z <= single(att.Band(1).NoDataValue)) = NaN;    end
    handles.have_nans = grdutils(Z,'-N');    
    zz = scaleto8(Z);       [m,n] = size(Z);
    X = linspace(head(1),head(2),n);  Y = linspace(head(3),head(4),m);  % Need this for image
    handles.geog = guessGeog(head(1:4));

    grd_size = m * n * 4;
    if (grd_size <= handles.grdMaxSize)         % Save grid if it is not to big
        setappdata(handles.figure1,'dem_z',Z);  setappdata(handles.figure1,'dem_x',X);
        setappdata(handles.figure1,'dem_y',Y);  setappdata(handles.figure1,'GMThead',head);
    end

    handles.DEM_name = [PathName FileName];     % Save DEM file name for eventual writing of a GMT grid
    handles.image_type = 4;
    ValidGrid = 1;                              % Signal that grid opps are allowed
    setappdata(handles.figure1,'Zmin_max',[head(5) head(6)])
    aux_funs('colormap_bg',handles,Z,jet(256));
end

handles.head = head;    handles.grdname = [];
show_image(handles,[PathName FileName],X,Y,zz,ValidGrid,'xy',1);
grid_info(handles,att,'gdal')           % Construct a info message

% --------------------------------------------------------------------
function FileOpenNewImage_Callback(hObject, eventdata, handles)
str1 = {'*.jpg', 'JPEG image (*.jpg)'; ...
    '*.png', 'Portable Network Graphics(*.png)'; ...
    '*.bmp', 'Windows Bitmap (*.bmp)'; ...
    '*.gif', 'GIF image (*.gif)'; ...
    '*.tif', 'Tagged Image File (*.tif)'; ...
    '*.pcx', 'Windows Paintbrush (*.pcx)'; ...
    '*.ras', 'SUN rasterfile (*.ras)'; ...
    '*.hdf', 'Hieralchical Data Format (*.hdf)'; ...
    '*.ppm', 'Portable Pixmap (*.ppm)'; ...
    '*.pgm', 'Portable Graymap (*.pgm)'; ...
    '*.raw;*.bin', 'RAW file (*.raw,*.bin)'; ...
    '*.shade', 'IVS shade File (*.shade)'; ...
    '*.xwd', 'X Windows Dump (*.xwd)'; ...
    '*.*', 'All Files (*.*)'};
[FileName,PathName,handles] = put_or_get_file(handles,str1,'Select image format','get');
if isequal(FileName,0);     return;     end

set(handles.figure1,'pointer','watch')
[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if (strcmpi(EXT,'.shade'))
    [fid, msg] = fopen([PathName FileName], 'r');
    if (fid < 0)
        set(handles.figure1,'pointer','arrow');    errordlg([PathName FileName ': ' msg],'ERROR');   return
    end
    fseek(fid, 56, 'bof');                  % Seek forward to the image data.
    nm = fread(fid,2,'uint16');             % read n_row & n_col
    n_row = nm(1);      n_col = nm(2);
    fseek(fid, 87, 'bof');                  % position the pointer at the end of the header
    nbytes = n_row*n_col*4;                 % image is of RGBA type, so it has 4 channels
    I = fread(fid,nbytes,'*uint8');     fclose(fid);
    I = reshape(I, [4 n_row n_col]);    I = permute(I, [3 2 1]);
    I(:,:,1) = flipud(I(:,:,1));        I(:,:,2) = flipud(I(:,:,2));
    I(:,:,3) = flipud(I(:,:,3));        I(:,:,4) = flipud(I(:,:,4));
    I = I(:,:,2:4);                         % strip alpha off of I
elseif (strcmpi(EXT,'.raw') || strcmpi(EXT,'.bin'))
    FileOpenGDALmultiBand_Callback(hObject, [], handles, 'RAW', [PathName FileName])
    return      % We are done here. Bye Bye.
else
    info_img = imfinfo([PathName FileName]);
    if ( any(strcmpi(EXT,{'.tif' '.tiff'})) && strcmp(info_img.Compression,'LZW'))
        % If Tiffs are LZW, imread R13 is not able to read them 
        I = gdalread([PathName FileName]);
    else
        [I,map] = imread([PathName FileName]);
    end
    if (strcmp(info_img.ColorType,'grayscale'))
        set(handles.figure1,'Colormap',gray(256))
    elseif (isfield(info_img,'ColorTable'))         % Gif images call it 'ColorTable'
        set(handles.figure1,'Colormap',info_img.ColorTable)
    elseif (isfield(info_img,'Colormap') && ~isempty(info_img.Colormap))
        set(handles.figure1,'Colormap',info_img.Colormap)
    end
end

handles.image_type = 2;     handles.grdname = [];       handles.geog = 0;    % None of this image types is coordinated
handles.head = [1 size(I,2) 1 size(I,1) 0 255 0 1 1];   % Fake a grid reg GMT header
aux_funs('cleanGRDappdata',handles);            % Remove eventual grid stuff variables from appdata
handles = show_image(handles,[PathName FileName],[],[],I,0,'off',0);
grid_info(handles,[PathName FileName],'iminfo');% Contruct a info string
guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function FileOpenGDALmultiBand_Callback(hObject, eventdata, handles, opt, opt2)
% Read GDAL files that may be multiband
% OPT2, if present, MUST contain the full file name. Currently used to load RAW images

if strcmp(opt,'ENVISAT')
    str1 = {'*.n1;*.N1', 'Envisat (*.n1,*.N1)'; '*.*', 'All Files (*.*)'};
elseif strcmp(opt,'AVHRR')
    str1 = {'*.n14;*.N14;*n17;*N17', 'NOAA (*.n14,*.N14,*.n17,*.N17)'; '*.*', 'All Files (*.*)'};
end

if (nargin == 4)    % Otherwise, OPT2 already contains the File name
    [FileName,PathName] = put_or_get_file(handles,str1,['Select ' opt ' File'],'get');
    if isequal(FileName,0),     return;     end
	fname = [PathName FileName];
else
    fname = opt2;
end

att.ProjectionRef = [];
reader = 'GDAL';                % this will be used by bands_list to decide which reader to use

set(handles.figure1,'pointer','watch')
if (strcmp(opt,'RAW'))          % 
    [I,cmd1,cmd2] = read_FlatFile({fname});
    if (isempty(I)),    set(handles.figure1,'pointer','arrow');     return;     end
    [att.RasterYSize,att.RasterXSize,n_bands] = size(I);
    bands_inMemory = 1:n_bands;         % Make it a vector
    reader = {cmd1; cmd2};
    handles.head = [1 size(I,2) 1 size(I,1) 0 255 0 1 1];   % Fake GMT header
elseif (strcmp(opt,'ENVISAT') || strcmp(opt,'AVHRR'))    % 
    bands_inMemory = 3;
    opt_B = ['-B1-' num2str(bands_inMemory)];
    [I,att] = gdalread(fname,'-S', opt_B,'-C');
    n_bands = att.RasterCount;
    bands_inMemory = 1:min(n_bands,bands_inMemory);      % Make it a vector
    handles.head = att.GMT_hdr;
end

tmp1 = cell(n_bands+1,2);    tmp2 = cell(n_bands+1,2);
tmp1{1,1} = opt;    tmp1{1,2} = opt;
for (i = 1:n_bands)
    tmp1{i+1,1} = ['band' sprintf('%d',i)];
    tmp1{i+1,2} = ['banda' sprintf('%d',i)];    % TEMP
    tmp2{i+1,1} = [num2str(i) 'x1 bip'];        % TEMP
    tmp2{i+1,2} = i;
end
tmp = {['+ ' opt]; I; tmp1; tmp2; fname; bands_inMemory; [att.RasterYSize att.RasterXSize n_bands]; reader};

setappdata(handles.figure1,'BandList',tmp)
handles.image_type = 2;
handles.was_int16 = 0;  handles.computed_grid = 0;
handles.grdname = [];   handles.geog = 0;       % None of this image types is coordinated (nor geog nor anything else)
aux_funs('cleanGRDappdata',handles);            % Remove eventual grid stuff variables from appdata
if (isempty(att.Band(1).ColorMap))
    set(handles.figure1,'Colormap',jet(256))
else
    set(handles.figure1,'Colormap',att.Band(1).ColorMap.CMap)
end
show_image(handles,fname,[],[],I,0,'off',0);    % It also guidata(...) & reset pointer
if (isappdata(handles.axes1,'InfoMsg')),    rmappdata(handles.axes1,'InfoMsg');     end
if (~strcmp(opt,'RAW')),    grid_info(handles,att,'gdal')           % Construct a info message
end

% --------------------------------------------------------------------
function FileOpenGeoTIFF_Callback(hObject, eventdata, handles, tipo)
switch lower(tipo)
    case 'geotiff',     str1 = {'*.tif;*.TIF;*.tiff;*.TIFF', 'GeoTiff (*.tif,*.tiff,*.TIF,*.TIFF)'};
    case 'sid',         str1 = {'*.sid;*.SID', 'GeoTiff (*.sid,*.SID)'};
    case 'ecw',         str1 = {'*.ecw;*.ECW', 'GeoTiff (*.ecw,*.ECW)'};
    case 'jp2',         str1 = {'*.jp2;*.JP2', 'GeoTiff (*.jp2,*.JP2)'};
    otherwise,          str1 = {'', 'Don''t know (*.*)'};       % Used with the "Try Luck" option
end
str1(2,1:2) = {'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,['Select ' tipo ' file'],'get');
if isequal(FileName,0);     return;     end

att = gdalread([PathName FileName],'-M','-C');
set(handles.figure1,'pointer','watch')

if (att.RasterCount > 3)            % Since it is a multiband file, try luck there
    FileOpenGDALmultiBand_Callback([],[],handles,'AVHRR',[PathName FileName]);
    return
end

if ~(strcmp(att.Band(1).DataType,'Byte'))       % JPK2, for example, may contain DTMs
    fullname{1} = PathName;     fullname{2} = FileName;
    read_DEMs(handles,fullname,'JP2_DEM');    return
end

Z = gdalread([PathName FileName],'-U');
handles.head = att.GMT_hdr;

if (~isempty(att.GeoTransform))     % Georeferenced image
    X = handles.head(1:2);      Y = handles.head(3:4);
    handles.head(8) = diff(X) / (size(Z,2) - 1);      handles.head(9) = diff(Y) / (size(Z,1) - 1);
    ax_dir = 'xy';              handles.image_type = 3;
    att.GMT_hdr(8:9) = handles.head(8:9);   % Update the attrib struct
else                                % Case of "raw" imagerie
    X = [1 size(Z,2)];      Y = [1 size(Z,1)];   ax_dir = 'off';     handles.image_type = 2;
    Z = flipdim(Z,1);       handles.head = [X Y 0 255 0 1 1];
end

if (strcmpi(att.ColorInterp,'palette'))
    if isempty(att.Band(1).ColorMap),       pal = jet(256);
    else
        try     pal = att.Band(1).ColorMap.CMap;     pal = pal(:,1:3);       % GDAL creates a Mx4 colormap
        catch   warndlg('Figure ColorMap had troubles. Replacing by a default one.','Warning'); pal = jet(256);
        end
    end
elseif (strcmpi(att.ColorInterp,'gray'))
    pal = repmat( (att.GMT_hdr(5):att.GMT_hdr(6))' / 255, 1, 3);
else                            pal = jet(256);
end

handles.geog = guessGeog([X Y]);
handles.grdname = [];   set(handles.figure1,'Colormap',pal);
aux_funs('cleanGRDappdata',handles);        % Remove eventual grid stuff variables from appdata
show_image(handles,[PathName FileName],X,Y,Z,0,ax_dir,0);
grid_info(handles,att,'gdal')               % Construct a info message

% --------------------------------------------------------------------
function FileOpenDEM_Callback(hObject, eventdata, handles, opt)
% Files of the following formats are read (well re-directed) here
tipo = opt;
switch opt
    case {'GMT' 'Surfer', 'ENCOM'}
        str1 = {'*.grd;*.GRD;*.nc;*.NC', 'Grid files (*.grd,*.GRD,*.nc,*.NC)';'*.*', 'All Files (*.*)'};    tipo = 'GMT_relatives';
    case 'MANI'
        str1 = {'*.man;*.MAN', 'Grid files (*.man,*.MAN)';'*.*', 'All Files (*.*)'};    tipo = 'GMT_relatives';
    case 'ArcAscii'
        str1 = {'*.grd;*.GRD', 'Arc/Info grid (*.grd,*.GRD)'; '*.*', 'All Files (*.*)'};
    case 'ArcBinary'
        str1 = {'*.adf;*.ADF', 'Arc/Info grid (*.adf,*.ADF)'; '*.*', 'All Files (*.*)'};
    case 'DTED'
        str1 = {'*.dt0;*.DT0;*.dt1;*.DT1', 'DTED (*.dt0,*.DT0,*.dt1,*.DT1)'; '*.*', 'All Files (*.*)'};
    case 'ESRI_hdr'
	    str1 = {'*.bil;*.BIL;', 'ESRI BIL (*.bil,*.BIL)'; '*.*', 'All Files (*.*)'};
    case 'GTOPO30'
        str1 = {'*.dem;*.DEM', 'GTOPO30 DEM (*.dem,*.DEM)'; '*.*', 'All Files (*.*)'};
    case 'GeoTiff_DEM'
        str1 = {'*.tif;*.TIF;*.tiff;*.TIFF', 'GeoTiff DEM(*.tif,*.tiff,*.TIF,*.TIFF)'; '*.*', 'All Files (*.*)'};
    case 'GXF'
        str1 = {'*.gxf;*.GXF', 'Geosoft GXF (*.gxf,*.GXF)'; '*.*', 'All Files (*.*)'};
    case 'SDTS'
        str1 = {'*catd.ddf;*CATD.DDF', 'USGS SDTS DEM (*catd.ddf,*CATD.DDF)'; '*.*', 'All Files (*.*)'};
    case 'SRTM30'
    	str1 = {'*.srtm;*.SRTM;*.srtm.gz', 'SRTM30 DEM (*.srtm,*.SRTM,*.srtm.gz)'; '*.*', 'All Files (*.*)'};
    case {'SRTM1' 'SRTM3'}
    	str1 = {'*.hgt;*.HGT;*.hgt.zip', [opt ' DEM (*.hgt,*.HGT,*.hgt.zip)']; '*.*', 'All Files (*.*)'};
    case 'USGS_DEM'
        str1 = {'*.dem;*.DEM', 'USGS DEM (*.dem,*.DEM)'; '*.*', 'All Files (*.*)'};
    otherwise
        errordlg(['OOPs, where did this ' opt ' code came in?'],'Error');   return
end
[FileName,PathName] = put_or_get_file(handles,str1,['Select ' opt ' File'],'get');
if isequal(FileName,0);     return;     end
fullname{1} = PathName;     fullname{2} = FileName;
read_DEMs(handles,fullname,tipo)

% --------------------------------------------------------------------
function FileOpenMOLA_Callback(hObject, eventdata, handles)
str1 = {'*.img;*.IMG', 'MOLA DEM (*.img,*.IMG)'; '*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select MOLA DEM File','get');
if isequal(FileName,0),     return;     end

type = 'MOLA';      error = 0;
[PATH,FNAME] = fileparts([PathName FileName]);
fname = [PATH filesep FNAME '.lbl'];
fp = fopen(fname,'rt');
if (fp < 0)
    errordlg(['ERROR: Could not find format descriptive file: ' fname],'Error');  return
end
s = strread(fread(fp,'*char').','%s','delimiter','\n');

LINES = findcell('LINES', s);                   if (isempty(LINES)),  error = 1;  end
[t,r] = strtok(s{LINES.cn},'=');                n_lines = str2double(r(3:end));

LINE_SAMPLES = findcell('LINE_SAMPLES', s);     if (isempty(LINE_SAMPLES)),  error = 1;  end
[t,r] = strtok(s{LINE_SAMPLES.cn},'=');         n_samples = str2double(r(3:end));

SAMPLE_BITS = findcell('SAMPLE_BITS', s);       if (isempty(SAMPLE_BITS)),  error = 1;  end
[t,r] = strtok(s{SAMPLE_BITS.cn},'=');          n_bits = str2double(r(3:end));

CENTER_LATITUDE = findcell('CENTER_LATITUDE', s);       if (isempty(CENTER_LATITUDE)),  error = 1;  end
[t,r] = strtok(s{CENTER_LATITUDE.cn},'=');      lat0 = str2double(strtok(strtok(r,'=')));

CENTER_LONGITUDE = findcell('CENTER_LONGITUDE', s);     if (isempty(CENTER_LONGITUDE)),  error = 1;  end
[t,r] = strtok(s{CENTER_LONGITUDE.cn},'=');     lon0 = str2double(strtok(strtok(r,'=')));

LINE_FIRST_PIXEL = findcell('LINE_FIRST_PIXEL', s);     if (isempty(LINE_FIRST_PIXEL)),  error = 1;  end
[t,r] = strtok(s{LINE_FIRST_PIXEL.cn},'=');     line_first_pix = str2double(r(3:end));

LINE_LAST_PIXEL = findcell('LINE_LAST_PIXEL', s);       if (isempty(LINE_LAST_PIXEL)),  error = 1;  end
[t,r] = strtok(s{LINE_LAST_PIXEL.cn},'=');      line_last_pix = str2double(r(3:end));

SAMPLE_FIRST_PIXEL = findcell('SAMPLE_FIRST_PIXEL', s); if (isempty(SAMPLE_FIRST_PIXEL)),  error = 1;  end
[t,r] = strtok(s{SAMPLE_FIRST_PIXEL.cn},'=');   sample_first_pix = str2double(r(3:end));

SAMPLE_LAST_PIXEL = findcell('SAMPLE_LAST_PIXEL', s);   if (isempty(SAMPLE_LAST_PIXEL)),  error = 1;  end
[t,r] = strtok(s{SAMPLE_LAST_PIXEL.cn},'=');    sample_last_pix = str2double(r(3:end));

MAP_RESOLUTION = findcell('MAP_RESOLUTION', s);         if (isempty(MAP_RESOLUTION)),  error = 1;  end
[t,r] = strtok(s{MAP_RESOLUTION.cn},'=');       res = str2double(strtok(strtok(r,'=')));

LINE_PROJECTION_OFFSET = findcell('LINE_PROJECTION_OFFSET', s);
if (isempty(LINE_PROJECTION_OFFSET)),  error = 1;  end
[t,r] = strtok(s{LINE_PROJECTION_OFFSET.cn},'=');line_off = str2double(r(3:end));

SAMPLE_PROJECTION_OFFSET = findcell('SAMPLE_PROJECTION_OFFSET', s);
if (isempty(SAMPLE_PROJECTION_OFFSET)),  error = 1;  end
[t,r] = strtok(s{SAMPLE_PROJECTION_OFFSET.cn},'=');sample_off = str2double(r(3:end));

if (error),     errordlg('An error has occured in header file parsing','Error');    return;     end
limits = [(sample_first_pix - sample_off)/res+lon0 (sample_last_pix - sample_off)/res+lon0 ...
    (line_first_pix - line_off)/res+lat0 (line_last_pix - line_off)/res+lat0];

fullname{1} = PathName;     fullname{2} = FileName;
read_DEMs(handles,fullname,type,[limits n_samples n_lines 1/res])

% --------------------------------------------------------------------
function read_DEMs(handles,fullname,tipo,opt)
% This function loads grid files that may contain DEMs or other grid (not images) types
% OPT is used when reading MOLA files, OPT = [x_min x_max y_min y_max n_cols n_rows grid_inc]

att = [];     % If a file is read by gdal, this won't be empty at the end of this function
set(handles.figure1,'pointer','watch')
PathName = fullname{1}; FileName = fullname{2};
if (handles.ForceInsitu),   opt_I = '-I';   % Use only in desperate cases.
else                        opt_I = ' ';    end
handles.was_int16 = 0;      % To make sure that it wasnt left = 1 from a previous use.
if (~isempty(getappdata(handles.figure1,'dem_z')))   % Free memory that may be necessary bellow
    rmappdata(handles.figure1,'dem_z')
end

if (strcmp(tipo,'GMT_relatives'))   % The reading is done by the read_gmt_type_grids function
    [handles,X,Y,Z,head] = read_gmt_type_grids(handles,[PathName FileName]);
    if (isempty(X)),    set(handles.figure1,'pointer','arrow');     return;     end
    zz = scaleto8(Z);
elseif ( strcmp(tipo,'SRTM30') || strcmp(tipo,'SRTM3') || strcmp(tipo,'SRTM1') || strcmp(tipo,'MOLA') )
    fname = [PathName FileName];        name_uncomp = [];
    tipo1 = tipo;
    if (strcmp(tipo,'MOLA'))    % Must inform write_ESRI_hdr what to write
        tipo1 = [opt(1) opt(4) opt(5) opt(6) opt(7) opt(7) -99999];
    end
    [name_hdr,comp_type] = write_ESRI_hdr(fname,tipo1);
    if (~isempty(comp_type))
        fname = decompress(fname,'warn');           % Named with compression ext removed
        name_uncomp = fname;        % Here we need a copy of the decompressed file name for removing
        if (isempty(fname)),  return;     end;      % Error message already issued.
    end

    [Z,att] =  gdalread(fname,'-U','-C',opt_I);   Z = single(Z);
    head = att.GMT_hdr;
    Z(Z <= single(att.Band(1).NoDataValue)) = NaN;     handles.have_nans = grdutils(Z,'-N');
    zz = scaleto8(Z);           [m,n] = size(Z);
    X = linspace(head(1),head(2),n);  Y = linspace(head(3),head(4),m);  % Need this for image
    handles.image_type = 1;     handles.computed_grid = 1;  handles.geog = 1;
    handles.was_int16 = 1;      handles.grdname = [];       handles.Nodata_int16 = att.Band(1).NoDataValue;
    delete(name_hdr);
    if (~isempty(name_uncomp)),     delete(name_uncomp);    end
elseif ( strcmp(tipo,'USGS_DEM') || strcmp(tipo,'GTOPO30') || strcmp(tipo,'DTED') || strcmp(tipo,'SDTS') || ...
        strcmp(tipo,'GeoTiff_DEM') || strcmp(tipo,'ArcAscii')  || strcmp(tipo,'ArcBinary') || ...
        strcmp(tipo,'GXF') || strcmp(tipo,'ESRI_hdr') || strcmp(tipo,'JP2_DEM'))
    att = gdalread([PathName FileName],'-M','-C');
    if ((att.RasterXSize * att.RasterYSize * 4) > handles.grdMaxSize)
        if ( strcmp(yes_or_no('title','Warning'),'Yes')),  return;     end      % Advise accepted
    end
    
    Z =  gdalread([PathName FileName],'-U',opt_I,'-C');     Z = single(Z);
    head = att.GMT_hdr;
    if (~isempty(att.Band(1).NoDataValue)),  Z(Z <= single(att.Band(1).NoDataValue)) = NaN;    end
    handles.have_nans = grdutils(Z,'-N');
    if (isequal(head(5:6),[0 0]))   % It happens with GeoTiff_DEM
        head(5) = double(min(Z(:)));        head(6) = double(max(Z(:)));
    end
    zz = scaleto8(Z);    [m,n] = size(Z);
    X = linspace(head(1),head(2),n);  Y = linspace(head(3),head(4),m);  % Need this for image
    
    handles.DEM_name = [PathName FileName]; % Save DEM file name for eventual writing of a GMT grid
    if (strcmp(att.Band(1).DataType,'Int16')), handles.was_int16 = 1;  end
    handles.image_type = 4;                 handles.grdname = [];   handles.Nodata_int16 = att.Band(1).NoDataValue;
end

% Guess if file is in geogs
switch tipo
    case {'GTOPO30' 'SRTM30' 'SRTM1' 'SRTM3'},    handles.geog = 1;
    otherwise
        handles.geog = guessGeog(head(1:4));    % A bit risky this guess
end

% If grid size is not to big I'll store it
if (numel(zz)*4 <= handles.grdMaxSize)
    setappdata(handles.figure1,'dem_z',Z);  setappdata(handles.figure1,'dem_x',X);
    setappdata(handles.figure1,'dem_y',Y);  setappdata(handles.figure1,'GMThead',head);
end

handles.head = head;
setappdata(handles.figure1,'Zmin_max',[head(5) head(6)])
aux_funs('colormap_bg',handles,Z,jet(256));
show_image(handles,[PathName FileName],X,Y,zz,1,'xy',head(7));
if (isappdata(handles.axes1,'InfoMsg')),    rmappdata(handles.axes1,'InfoMsg');     end
if (~isempty(att)),    grid_info(handles,att,'gdal');   end            % Construct a info message

% --------------------------------------------------------------------
function handles = show_image(handles,fname,X,Y,I,ValidGrid,axis_t,adjust,imSize)
% Show image and set other parameters
if (adjust)         % Convert the image limits from pixel reg to grid reg
    [m,n,k] = size(I);  [X,Y] = aux_funs('adjust_lims',X,Y,m,n);
end
if (nargin < 9),        imSize = [];    end
if (ValidGrid && ~isempty(X)),  dx = X(2) - X(1);       dy = Y(2) - Y(1);
else                            dx = 0;                 dy = 0;
end
if (isempty(imSize) && (abs(dx - dy) > 1e-4))       % Check for grid node spacing anisotropy
    imSize = dx / dy;                               % resizetrue will know what to do with this
end
if (strcmp(get(findobj(handles.figure1,'Tag','GCPtool'),'Checked'),'on'))
%     handles = gcpTool_old(handles,axis_t,X,Y,I);
    handles = gcpTool(handles,axis_t,X,Y,I);
    return
end

handles.grd_img = image(X,Y,I,'Parent',handles.axes1);
zoom_state(handles,'off_yes');      set(handles.grd_img,'CDataMapping','direct')
if (strcmp(axis_t,'xy')),           set(handles.axes1,'XDir','normal','YDir','normal')
elseif (strcmp(axis_t,'off')),      set(handles.axes1,'Visible','off')
else    warndlg('Warning: Unknown axes setting in show_image','Warning')
end
resizetrue(handles.figure1,imSize);
handles.origFig = I;            handles.no_file = 0;            handles.BgMap = 0;
handles.Illumin_type = 0;       handles.ValidGrid = ValidGrid;  % Signal that gmt grid opps are allowed
setappdata(handles.figure1,'ValidGrid',ValidGrid)   % Flag to signal draw_funs if "Crop grid" is possible
set(handles.figure1,'Name',fname)
SetAxesNumericType(handles);                        % Set axes uicontextmenus
setappdata(handles.figure1,'ThisImageLims',[get(handles.axes1,'XLim') get(handles.axes1,'YLim')])
setappdata(handles.axes1,'ThisImageLims',[get(handles.axes1,'XLim') get(handles.axes1,'YLim')])     % Used in pan and somewhere else
setappdata(handles.figure1,'GMThead',handles.head);
handles.oldSize = get(handles.figure1,'Pos');       % Save fig size to prevent maximizing
handles.origCmap = get(handles.figure1,'ColorMap'); % Save original colormap 
set(findobj('Tag','ImgHist'),'checked','off');  set(findobj('Tag','ImgHistGrd'),'checked','off');
% Make an extra copy of those to use in "restore" because they may be changed by 'bands_list()'
handles.ValidGrid_orig = ValidGrid;         handles.was_int16_orig = handles.was_int16;
handles.computed_grid_orig = handles.computed_grid;
guidata(handles.figure1, handles);              set(handles.figure1,'pointer','arrow')

if(~ValidGrid)      % Delete uicontrols that are useless to images only
    delete(findobj(handles.figure1,'Tag','ImgHistGrd'));    delete(findobj(handles.figure1,'Tag','Illuminate'))
    delete(findobj(handles.figure1,'Tag','GridTools'));     delete(findobj(handles.figure1,'Tag','Contours_a'))
    delete(findobj(handles.figure1,'Tag','Contours_i'));
    set(findobj(handles.figure1,'Tag','SaveGMTgrid'),'Enable','off');
end

if (ndims(I) == 3)          % Some cheating to allow selecting individual bands of a RGB image
    tmp1 = cell(4,2);   tmp2 = cell(4,2);    tmp1{1,1} = 'RGB';     tmp1{1,2} = 'RGB';
    for (i = 1:3)
        tmp1{i+1,1} = ['band' sprintf('%d',i)];     tmp1{i+1,2} = ['banda' sprintf('%d',i)];
        tmp2{i+1,1} = [num2str(i) num2str(size(I,1)) 'x' num2str(size(I,2)) ' BSQ'];        tmp2{i+1,2} = i;
    end
    tmp = {['+ ' 'RGB']; I; tmp1; tmp2; ''; 1:3; [size(I,1) size(I,2) 3]; 'Mirone'};
    setappdata(handles.figure1,'BandList',tmp)
    set(findobj(handles.figure1,'Label','Load Bands'),'Enable','on')
elseif (ndims(I) == 2)      % Remove it so it won't try to operate on indexed images
    %if (~isempty(getappdata(handles.figure1,'BandList'))),      rmappdata(handles.figure1,'BandList');  end
    if (isappdata(handles.figure1,'BandList')),     rmappdata(handles.figure1,'BandList');  end
    set(findobj(handles.figure1,'Label','Load Bands'),'Enable','off')
end
if (isappdata(handles.axes1,'ProjWKT')),            rmappdata(handles.axes1,'ProjWKT');         end
if (isappdata(handles.axes1,'DatumProjInfo')),      rmappdata(handles.axes1,'DatumProjInfo');   end

% --------------------------------------------------------------------
function ToolsMBplaningStart_Callback(hObject, eventdata, handles)
if (handles.image_type ~= 1),   aux_funs('msg_dlg',2,handles);     return;      end
if (handles.geog ~= 1),         aux_funs('msg_dlg',3,handles);     return;      end
if isempty(getappdata(handles.figure1,'dem_x'))     % Test if the grid is loaded in memory
    warndlg('Grid file is bigger than the declared "Grid Max Size". See "File -> Preferences"','Warning');
    return
end

prompt = {['The current value for the swath-width / water depth ratio is:   --> ' num2str(handles.swathRatio) '  <--']
        'If you want to change it, hit "Cancel", and do it in "File -> Preferences"'
        ' '
        'NOTE: this message, once accepted, is shown only once. So, if in the midle'
        'of a session, you decide to change the Swath Ratio don''t forget to do it'
        'there. The current value will be applyied to all new tracks without further questions.'};
if (handles.firstMBtrack == 1)
    button = questdlg(prompt,'Multibeam planing info','Accept','Cancel','Accept');
    if strcmp(button,'Accept'),     handles.firstMBtrack = 0;
    else                            return;    end
end

zoom_state(handles,'maybe_off');
[xp,yp,trackHand,barHand] = getline_mb;
if (isempty(xp)),   zoom_state(handles,'maybe_on');    return;  end

% Now we have to join all trackHand's in a single handle
handles.nTrack = handles.nTrack + 1;    % count to the number of tracks
delete(trackHand(2:end));       trackHand(2:end) = [];
% make and set tags strings to tracks and track's Bars
tagL = ['MBtrack' num2str(handles.nTrack)];         tagB = ['swath_w' num2str(handles.nTrack)];
set(trackHand,'XData',xp, 'YData',yp, 'Tag',tagL);  set(barHand, 'Tag',tagB)
setappdata(trackHand,'swathRatio',handles.swathRatio)  % save the swathRatio in line's appdata
% now set the barHand userdata that contains the vertex order
for i = 1:length(xp),     set(barHand(i),'Userdata',i);     end

draw_funs(trackHand,'MBtrack_uicontext')        % Set track's uicontextmenu
draw_funs(barHand,'MBbar_uicontext')            % Set track bar's uicontextmenu

zoom_state(handles,'maybe_on');
handles.h_MBplot = trackHand;       guidata(hObject, handles);
str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
str2 = 'Select output xy file name';
[FileName,PathName] = put_or_get_file(handles,str1,str2,'put');
if isequal(FileName,0);     return;     end
double2ascii([PathName FileName],[xp yp],'%f\t%f');

% --------------------------------------------------------------------
function ToolsMBplaningImport_Callback(hObject, eventdata, handles)
if (handles.image_type ~= 1),   aux_funs('msg_dlg',2,handles);     return;      end
if (handles.geog ~= 1),         aux_funs('msg_dlg',3,handles);     return;      end
if isempty(getappdata(handles.figure1,'dem_x'))     % Test if the grid is loaded in memory
    warndlg('Grid file is bigger than the declared "Grid Max Size". See "File -> Preferences"','Warning');  return
end
str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select input xy file name','get');
if isequal(FileName,0);     return;     end

out = draw_funs([PathName FileName],'ImportLine');
prompt = {'Enter swath-width / water depth ratio'};
resp  = inputdlg(prompt,'Multi-beam planing input',[1 38],{'3'});     pause(0.01);
if isempty(resp);   return;     end
if (iscell(out)),   n_segments = length(out);
else                n_segments = 1;     end
hold on
for (i = 1:n_segments)
    if (iscell(out)),   xy = out{i}(:,1:2);
    else                xy = out(:,1:2);    end
    z = abs(bi_linear(getappdata(handles.figure1,'dem_x'),getappdata(handles.figure1,'dem_y'),...
        getappdata(handles.figure1,'dem_z'),xy(:,1),xy(:,2)));
    rad = abs(z) * (str2double(resp{1})/2) / 111194.9; % meters -> degrees
    handles.nTrack = handles.nTrack + 1;    % count the number of imported tracks
    % make tags strings to tracks and track's Bars
    tagL = ['MBtrack' num2str(handles.nTrack)];     tagB = ['swath_w' num2str(handles.nTrack)];
    nr = size(xy,1);
    az = zeros(nr,1);    h_circ = zeros(1,nr);
    az(2:nr) = azimuth_geo(xy(1:(nr-1),2),xy(1:(nr-1),1),xy(2:nr,2),xy(2:nr,1));
    az(1) = az(2);              % First and second points have the same azim
    for (k=1:nr)                % NOTE: the Tag is very important for line edition
        [lat1,lon1] = circ_geo(xy(k,2),xy(k,1),rad(k),[az(k)-90-1 az(k)-90+1],3);
        [lat2,lon2] = circ_geo(xy(k,2),xy(k,1),rad(k),[az(k)+90-1 az(k)+90+1],3);
        h_circ(k) = line('XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)], ...
            'Color', [.8 .8 .8],'LineWidth',4,'Tag',tagB,'Userdata',k);
    end        
    h_line(i) = line('XData', xy(:,1), 'YData', xy(:,2),'Color','k','LineWidth',1,'Tag',tagL);
    setappdata(h_line(i),'swathRatio',str2double(resp{1}))  % save the swathRatio in line's appdata
    draw_funs(h_line(i),'MBtrack_uicontext')    % Set track's uicontextmenu
    draw_funs(h_circ,'MBbar_uicontext')         % Set track bar's uicontextmenu
end
hold off    
handles.h_MBplot = h_line(1);   guidata(hObject, handles);

% --------------------------------------------------------------------
function ToolsMBplaningEdit_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),      return;         end
if ~isempty(handles.h_MBplot),  edit_track_mb;  end

% --------------------------------------------------------------------
function ToolsMBplaningSave_Callback(hObject, eventdata, handles, opt)
% Note. This function is also used to save polylines as asked from draw_funs
if (handles.no_file == 1),      return;     end
if (nargin == 3),               opt = [];   end

if (ishandle(opt)),     x = get(opt,'XData');   y = get(opt,'YData');
else                    [x,y] = save_track_mb;  end
if isempty(x);          warndlg('Empty line; exiting','Warning');    return;  end

str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select output xy file name','put');
if isequal(FileName,0);     return;     end

[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if isempty(EXT),    f_name = [PathName FNAME '.dat'];
else                f_name = [PathName FNAME EXT];       end

% Save data with a format determined by axes format
labelType = getappdata(handles.figure1,'LabelFormatType');  % find the exes label format
if isempty(labelType),      labelType = ' ';        end     % untempered matlab axes labels
switch labelType
    case {' ','DegDec','NotGeog'}
        double2ascii(f_name,[x(:) y(:)],'%f\t%f');
    case 'DegMin'
        out_x = degree2dms(x,'DDMM',0,'numeric');        out_y = degree2dms(y,'DDMM',0,'numeric');
        double2ascii(f_name,[out_x.dd(:) out_x.mm(:) out_y.dd(:) out_y.mm(:)],'%4d %02d\t%4d %02d');
    case 'DegMinDec'        % I'm writing the minutes with a precision of 2 decimals
        out_x = degree2dms(x,'DDMM.x',2,'numeric');      out_y = degree2dms(y,'DDMM.x',2,'numeric');
        double2ascii(f_name,[out_x.dd(:) out_x.mm(:) out_y.dd(:) out_y.mm(:)],'%4d %02.2f\t%4d %02.2f');
    case 'DegMinSec'
        out_x = degree2dms(x,'DDMMSS',0,'numeric');      out_y = degree2dms(y,'DDMMSS',0,'numeric');
        double2ascii(f_name,[out_x.dd(:) out_x.mm(:) out_x.ss(:) ...
                out_y.dd(:) out_y.mm(:) out_y.ss(:)],'%4d %02d %02d\t%4d %02d %02d');
    case 'DegMinSecDec'     % I'm writing the seconds with a precision of 2 decimals
        out_x = degree2dms(x,'DDMMSS',2,'numeric');      out_y = degree2dms(y,'DDMMSS',2,'numeric');
        double2ascii(f_name,[out_x.dd(:) out_x.mm(:) out_x.ss(:) ...
                out_y.dd(:) out_y.mm(:) out_y.ss(:)],'%4d %02d %02.2f\t%4d %02d %02.2f');
end

% --------------------------------------------------------------------
function ToolsMBplaningDelete_Callback(hObject, eventdata, handles)
% The input argument must be "1" for deleting or "2" for computing the track length
if (handles.no_file == 1),    return;      end
save_track_mb(1);

% --------------------------------------------------------------------
function ToolsMBplaningTrackLength_Callback(hObject, eventdata, handles)
% The input argument must be "1" for deleting or "2" for computing the track length
if (handles.no_file == 1),    return;      end
dist = save_track_mb(2);
if (isempty(dist)),     warndlg('Empty line; exiting','Warning');    return;  end
msgbox(['Track length = ' num2str(dist,4) ' NM'])

% --------------------------------------------------------------------
function ImageIlluminationModel_Callback(hObject, eventdata, handles, opt)
if (aux_funs('msg_dlg',14,handles));     return;      end
if (nargin == 3),   opt = 'grdgradient_A';   end

luz = shading_params(opt);
pause(0.01)         % Give time to the shading_params window be deleted
if (isempty(luz)),  return;     end

if (luz.illum_model == 1)       % GMT grdgradient classic
    ImageIlluminateLambertian(hObject, luz, handles, 'grdgradient_class')
elseif (luz.illum_model == 2)   % GMT grdgradient Lambertian
    ImageIlluminateLambertian(hObject, luz, handles, 'grdgradient_lamb')
elseif (luz.illum_model == 3)   % GMT grdgradient Peucker
    ImageIlluminateLambertian(hObject, luz, handles, 'grdgradient_peuck')
elseif (luz.illum_model == 4)   % Lambertian illumination
    ImageIlluminateLambertian(hObject, luz, handles, 'lambertian')
elseif (luz.illum_model == 5)   % ManipRaster color illumination algo
    ImageIlluminateGray(hObject, luz, handles, 'color')
elseif (luz.illum_model == 6)   % ManipRaster gray illumination algo
    ImageIlluminateGray(hObject, luz, handles, 'gray')
else                            % False color illumination
    ImageIlluminateFalseColor(hObject, luz, handles)
end

% --------------------------------------------------------------------
function ImageIlluminateLambertian(hObject, luz, handles, opt)
% OPT == 'lambertian'   Use the Lambertian illumination algorithm
% OPT == whatever       Use the GMT grdgradient illumination algorithm
% Illuminate a DEM file and turn it into a RGB image
% For multiple tryies I need to use the original image. Otherwise each attempt would illuminate
% the previously illuminated image. An exception occurs when the image was IP but only for the
% first time, repeated illums will use the original img. Otherwise we would need to make another img copy

[X,Y,Z,head] = load_grd(handles);   % If needed, load gmt grid again
if isempty(Z),   return;     end;   % An error message was already issued

set(handles.figure1,'pointer','watch')

if (handles.firstIllum),    img = get(handles.grd_img,'CData');     handles.firstIllum = 0;
else                        img = handles.origFig;      end
if (isempty(img)),          img = get(handles.grd_img,'CData');     end             % No copy in memory
if (ndims(img) == 2),       img = ind2rgb8(img,get(handles.figure1,'Colormap'));    end    % Image is 2D   

setappdata(handles.figure1,'Xmin_max',[X(1) X(end)]);   setappdata(handles.figure1,'Ymin_max',[Y(1) Y(end)]);
setappdata(handles.figure1,'Zmin_max',[head(5) head(6)])

if (strcmp(opt,'grdgradient_class'))    % GMT grdgradient classic illumination
    if (handles.geog),  [R,offset,sigma] = grdgradient_m(Z,head,'-M',['-A' num2str(luz.azim)],'-Nt');
    else                [R,offset,sigma] = grdgradient_m(Z,head,['-A' num2str(luz.azim)],'-Nt'); end
    img = shading_mat(img,R,'no_scale');
    handles.Illumin_type = 1;
    handles.grad_offset = offset;   handles.grad_sigma = sigma;
elseif (strcmp(opt,'grdgradient_lamb')) % GMT grdgradient lambertian illumination
    R = grdgradient_m(single(Z),head,['-Es' num2str(luz.azim) '/' num2str(luz.elev)]);
    img = shading_mat(img,R,'no_scale');
    handles.Illumin_type = 2;
elseif (strcmp(opt,'grdgradient_peuck'))% GMT grdgradient Peucker illumination
    R = grdgradient_m(single(Z),head,'-Ep');
    img = shading_mat(img,double(R),'no_scale');
    handles.Illumin_type = 0;           % Should be 3 but its not programmed
elseif (strcmp(opt,'lambertian'))       % Lambertian lighting illumination
    R = grdgradient_m(Z,head,['-E' num2str(luz.azim) '/' num2str(luz.elev) '/' num2str(luz.ambient),...
           '/' num2str(luz.diffuse) '/' num2str(luz.specular) '/' num2str(luz.shine)]);
    img = shading_mat(img,R,'no_scale');           handles.Illumin_type = 4;
end

setappdata(handles.figure1,'Luz',luz);   % Save these for region operations
set(handles.grd_img,'CData',img)
guidata(hObject, handles);              set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function ImageIlluminateGray(hObject, luz, handles, color)
% Illuminate a previously open GMT or DEM file and turn it into a RGB image
% For multiple tryies I need to use the original image. Otherwise each attempt would illuminate
% the previously illuminated image. An exception occurs when the image was IP but only for the
% first time, repeated illums will use the original img. Otherwise we would need to make another img copy

[X,Y,Z,head,m] = load_grd(handles);
if isempty(Z),        return;     end;    % An error message was already issued
D2R = pi/180;

set(handles.figure1,'pointer','watch')

% Tiling
[ind_s,ind] = tile(m,400,4);      % shade_manip_raster "only consumes" 3 times Z grid size
if size(ind_s,1) > 1
    img = [];
    for i = 1:size(ind_s,1)
        tmp1 = (ind_s(i,1):ind_s(i,2));     % Indexes with overlapping zone
        tmp2 = ind(i,1):ind(i,2);           % Indexes of chunks without the overlaping zone
        tmp = shade_manip_raster((luz.azim-90)*D2R,luz.elev*D2R,Z(tmp1,:));
        img = [img; tmp(tmp2,:)];
    end
else
    img = shade_manip_raster((luz.azim-90)*D2R,luz.elev*D2R,Z);     % [0-1] matrix (reflectance)
end

if (strcmp(color,'color'))
	if (handles.firstIllum),    img1 = get(handles.grd_img,'CData');     handles.firstIllum = 0;
	else                        img1 = handles.origFig;      end
	if (isempty(img)),          img1 = get(handles.grd_img,'CData');     end             % No copy in memory
	if (ndims(img) == 2),       img1 = ind2rgb8(img1,get(handles.figure1,'Colormap'));    end    % Image is 2D   
    img = shading_mat(img1,img);
else
    img = uint8((254 * img) + 1);   % Need to convert the reflectance matrix into a gray indexed image
end

handles.Illumin_type = 0;       % This is temporary. It should be either 5 (color) or 6 (gray)
handles.no_file = 0;
set(handles.grd_img,'CData',img);
if (strcmp(color,'gray')),  set(handles.figure1,'Colormap',gray(256));    end     % Only gray images need colormap seting
set(handles.figure1,'pointer','arrow');     guidata(hObject, handles);

% --------------------------------------------------------------------
function ImageIlluminateFalseColor(hObject, luz, handles)
% Illuminate a GMT grid from 3 different directions and turn it into a RGB image

[X,Y,Z,head,m] = load_grd(handles);
if isempty(Z),        return;     end;    % An error message was already issued
set(handles.figure1,'pointer','watch')
D2R = pi/180;

% Tiling
[ind_s,ind] = tile(m,400,4);      % shade_manip_raster "only consumes" 3 times Z grid size
if size(ind_s,1) > 1
    zz1 = uint8([]);    zz2 = uint8([]);    zz3 = uint8([]);
    for i = 1:size(ind_s,1)
        tmp1 = (ind_s(i,1):ind_s(i,2));     % Indexes with overlapping zone
        tmp2 = ind(i,1):ind(i,2);           % Indexes of chunks without the overlaping zone
        tmp_1 = shade_manip_raster((luz.azim(1)-90)*D2R,luz.elev*D2R,Z(tmp1,:));
        tmp_2 = shade_manip_raster((luz.azim(2)-90)*D2R,luz.elev*D2R,Z(tmp1,:));
        tmp_3 = shade_manip_raster((luz.azim(3)-90)*D2R,luz.elev*D2R,Z(tmp1,:));
        zz1 = [zz1; tmp_1(tmp2,:)];
        zz2 = [zz2; tmp_2(tmp2,:)];
        zz3 = [zz3; tmp_3(tmp2,:)];
    end
    zz(:,:,1) = zz1;    zz(:,:,2) = zz2;    zz(:,:,3) = zz3;
else
    zz(:,:,1) = shade_manip_raster((luz.azim(1)-90)*D2R,luz.elev*D2R,Z);
    zz(:,:,2) = shade_manip_raster((luz.azim(2)-90)*D2R,luz.elev*D2R,Z);
    zz(:,:,3) = shade_manip_raster((luz.azim(3)-90)*D2R,luz.elev*D2R,Z);
end

handles.Illumin_type = 0;       % This is temporary. It should be 5
set(handles.grd_img,'CData',zz);
set(handles.figure1,'pointer','arrow');     guidata(hObject, handles);

% --------------------------------------------------------------------
function ImageAnaglyph_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',14,handles));     return;      end
[X,Y,Z,head] = load_grd(handles);
if isempty(Z),   return;     end;    % An error message was already issued

set(handles.figure1,'pointer','watch')
[ny,nx] = size(Z);    D2R = pi/180;         deg2m = 111194.9;
x_min = head(1);      y_min = head(3);      x_max = head(2);     y_max = head(4);
z_min = head(5);      z_max = head(6) + 1;  x_inc = head(8);      y_inc = head(9);

if (handles.geog == 1)
    p_size = sqrt((x_inc * deg2m) * (y_inc * deg2m * cos ((y_max + y_min)*0.5 * D2R))); 
else     p_size = x_inc * y_inc;    end

azimuth	= -90 * D2R;     elevation = 20 * D2R;

% Tiling
[m,n] = size(Z);
[ind_s,ind] = tile(m,400,4);      % shade_manip_raster "only consumes" 3 times Z grid size
if (size(ind_s,1) > 1)
    sh = [];
    for i = 1:size(ind_s,1)
        tmp1 = (ind_s(i,1):ind_s(i,2));     % Indexes with overlapping zone
        tmp2 = ind(i,1):ind(i,2);           % Indexes of chunks without the overlaping zone
        tmp = shade_manip_raster(azimuth,elevation,Z(tmp1,:));
        sh = [sh; tmp(tmp2,:)];
    end
else
    sh = shade_manip_raster(azimuth,elevation,Z);
end
sh = uint8((254 * sh) + 1);

str_amp = 2 ;
alpha = tan(25 * D2R) * str_amp / p_size;
decal = 1 + fix(2 * alpha * (z_max - z_min));
ana_header.nx = nx + decal;
ana_header.x_min = x_min - x_inc * (decal / 2);
ana_header.x_max = x_max + x_inc * (decal / 2);

left = repmat(uint8(255), 1, ana_header.nx);    right = left;
ar = repmat(uint8(0), ny, ana_header.nx);       ag = ar;    l = 0;  r = l;
for i=1:ny
    for j=1:nx
        iz=fix(alpha * (double(Z(i,j)) - z_min));
        if (j == 1)
            left(j+iz) = sh(i,j);
            right(decal+j-iz) = sh(i,j);
        else
            for (k=r:decal + j - iz),   right(k) = sh(i,j);  end
            for (k=l:j + iz),           left(k)  = sh(i,j);  end
        end
        l = j + iz;     r = decal + j-iz;
    end
    ar(i,:) = left;      ag(i,:) = right;
    left  = repmat(uint8(0), 1, ana_header.nx);     right = repmat(uint8(0), 1, ana_header.nx);
end
zz(:,:,1) = ar;    zz(:,:,2) = ag;  zz(:,:,3) = ag; 
show_image(handles,'Anaglyph',X,Y,zz,1,'xy',0);

% --------------------------------------------------------------------
function img = shade_manip_raster(azimuth,elevation,Z)
size_amp = 125;
u1 = sin(azimuth) * cos(elevation);     u2 = -cos(azimuth) * cos(elevation);
u3 = -sin(elevation);
if (~isa(Z,'double')),  Z = double(Z);  end     % Make sure Z is of double type

% Derivatives of function with respect to rows and columns
dZdc = zeros(size(Z));      dZdr = zeros(size(Z));

% Take forward differences on left and right edges
dZdr(1,:) = Z(2,:) - Z(1,:);      dZdr(end,:) = Z(end,:) - Z(end-1,:);
dZdc(:,1) = Z(:,2) - Z(:,1);      dZdc(:,end) = Z(:,end) - Z(:,end-1);

% Take centered differences on interior points
dZdr(2:end-1,:) = (Z(3:end,:)-Z(1:end-2,:)) / size_amp;
dZdc(:,2:end-1) = (Z(:,3:end)-Z(:,1:end-2)) / size_amp;
img  = (dZdr*u1 + dZdc*u2 - 2*u3) ./ (sqrt(dZdr .^ 2 + dZdc .^ 2 + 4));
img(img < 0) = 0;   img(img > 1) = 1;

% --------------------------------------------------------------------
function ImageDrape_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),    return;      end
son_img = get(findobj(handles.axes1,'Type','image'),'CData');     % Get "son" image
h_f = getappdata(0,'parent_gcf');                       % Get the parent figure handle
try
    h_a = findobj(h_f,'Type','axes');                   % Get the parent figure axes handle
catch
    msgbox('Parent window no longer exists (you kiled it). Exiting.','Warning');
    set(handles.ImageDrape,'Enable','off');             % Set the Drape option to it's default value (off)
    return
end
h_parent_img = findobj(get(h_a,'Children'),'Type','image');     % Get parent image handle
parent_img = get(findobj(h_f,'Type','image'),'CData');
[y_son x_son z] = size(son_img);                        % Get "son" image dimensions 
[y_parent x_parent z] = size(parent_img);               % Get "parent" image dimensions

% Find if image needs to be ud fliped
flip = 0;
if(strcmp(get(handles.axes1,'YDir'),'reverse')),    flip = 1;   end
% See about transparency
dlg_title = 'Draping Transparency';        num_lines= [1 38];  defAns = {'0'};
resp  = inputdlg('Use Transparency (0-1)?',dlg_title,num_lines,defAns);     pause(0.01);
if (isempty(resp) || str2double(resp{1}) < 0.01),       no_alfa = 1;
elseif (str2double(resp{1}) > 1),       no_alfa = 0;    alfa = 1;
else                                    no_alfa = 0;    alfa = str2double(resp{1}); end

son_cm = [];
if (ndims(son_img) == 2),       son_cm = get(handles.figure1,'Colormap');   end       % Get "son" color map
if(strcmp(get(handles.axes1,'YDir'),'reverse')),    son_img = flipdim(son_img,1);   end
if (y_son ~= y_parent || x_son ~= x_parent)              % Check if "son" and "parent" images have the same size
    son_img = cvlib_mex('resize',son_img,[y_parent x_parent],'bicubic');
end
if (~no_alfa)
    if (ndims(son_img) == 2),       son_img = ind2rgb8(son_img,get(handles.figure1,'Colormap'));    end
    if (ndims(parent_img) == 2),    parent_img = ind2rgb8(parent_img,get(h_f,'Colormap'));      end
    %son_img = uint8(double(parent_img) * alfa + double(son_img) * (1 - alfa));
    cvlib_mex('addweighted',son_img,(1 - alfa),parent_img,alfa)     % In-place
end
set(h_parent_img,'CData',son_img);
if (~isempty(son_cm) && no_alfa),    set(h_f,'Colormap',son_cm);     end  % Set "son" colormap to "parent" figure


% if (y_son ~= y_parent || x_son ~= x_parent)              % Check if "son" and "parent" images have the same size
%     head = [1 x_son 1 y_son 0 255 0 1 1];
%     opt_N = ['-N' num2str(x_parent) '/' num2str(y_parent)]; % option for grdsample
%     if (ndims(son_img) == 3)                            % RGB image
%         %for (i=1:3),    ZI(:,:,i) = grdsample_m(single(son_img(:,:,i)),head,opt_N);    end
%         ZI = cvlib_mex('resize',son_img,[y_parent x_parent],'bicubic');
%         if (flip),      ZI = flipdim(ZI,1);   end
%         set(h_parent_img,'CData',uint8(ZI))
%     else                                                % Indexed image
%         son_cm = get(handles.figure1,'Colormap');       % Get "son" color map
%         %ZI = grdsample_m(single(son_img),head,opt_N);
%         ZI = cvlib_mex('resize',son_img,[y_parent x_parent],'bicubic');
%         if (flip),  set(h_parent_img,'CData',flipud(uint8(ZI)))
%         else        set(h_parent_img,'CData',uint8(ZI));    end
%         set(h_f,'Colormap',son_cm)                      % Set "son" colormap to "parent" figure
%     end
% else                                                    % "son" and "parent" images have the same size
%     if (ndims(son_img) == 3)                            % RGB image
%         % Fliping is not allways true (e.g. if image comes from a
%         % grd illumination). In those cases, this is wrong.
%         if (flip),  son_img = flipdim(son_img,1);   end
%         nans_son = isnan(getappdata(handles.figure1,'dem_z'));
% 		dlg_title = 'Draping Transparency';        num_lines= [1 38];  defAns = {'0'};
% 		resp  = inputdlg('Use Transparency (0-1)?',dlg_title,num_lines,defAns);     pause(0.01);
% 		if (isempty(resp) || str2double(resp{1}) < 0.01),       no_alfa = 1;
% 		elseif (str2double(resp{1}) > 1),       no_alfa = 0;    alfa = 1;
%         else                                    no_alfa = 0;    alfa = str2double(resp{1}); end
%         if (ndims(parent_img) == 1),    parent_img = ind2rgb8(parent_img,get(h_f,'Colormap'));   end
%         if (any(nans_son(:)))
%             [m,n,k] = size(son_img);    id = repmat(reshape(nans_son,m,n),[1 1 3]);
%             son_img(id) = parent_img(id);
%         end
%         if (~no_alfa)       % Apply transparency
%             son_img = uint8(double(parent_img) * alfa + double(son_img) * (1 - alfa));
%         end
%         set(h_parent_img,'CData',son_img)
%     else                                                % Indexed image
%         if (flip),  set(h_parent_img,'CData',flipud(son_img));
%         else        set(h_parent_img,'CData',son_img);      end
%         set(h_f,'Colormap',get(handles.figure1,'Colormap'))     % Set "son" colormap to "parent" figure
%     end
% end
% Signal in the parent image handles that it has a draped image
handles = guidata(h_f);     handles.is_draped = 1;      handles.Illumin_type = 0;
guidata(h_f,handles)

% --------------------------------------------------------------------
function ImageColorPalettes_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),   color_palettes;     return;    end
if (ndims(get(handles.grd_img,'CData')) == 3)
    warndlg('True color images do not use color palettes. So you cannot change it.','Warning'); return
end
color_palettes(handles.figure1,handles.home_dir,handles.work_dir)

% --------------------------------------------------------------------
function ToolsMeasureDist_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),    return;      end
zoom_state(handles,'maybe_off');
[xp,yp] = getline_j(handles.figure1);  n_nodes = length(xp);
if (n_nodes < 2);     zoom_state(handles,'maybe_on');   return;     end
if (handles.geog)
    lat_i = yp(1:end-1);   lat_f = yp(2:end);
    lon_i = xp(1:end-1);   lon_f = xp(2:end);
    ellips = handles.DefineEllipsoide;      tmp = [];
    for (i=1:length(lat_i))
        s = vdist(lat_i(i),lon_i(i),lat_f(i),lon_f(i),ellips);
        tmp = [tmp s];
    end
    switch handles.DefineMeasureUnit
        case 'n',       scale = 1852;   str_unit = ' NM';           % Nautical miles
        case 'k',       scale = 1000;   str_unit = ' kilometers';   % Kilometers
        case {'m','u'}, scale = 1;      str_unit = ' meters';       % Meters or user unites
    end
    dist = sum(tmp) / scale;
else
    dx = diff(xp);      dy = diff(yp);
    dist = sum(sqrt(dx.*dx + dy.*dy));
    str_unit = ' map units';
end
zoom_state(handles,'maybe_on');
msgbox(['Distance = ' num2str(dist,'%.5f') str_unit])

% --------------------------------------------------------------------
function ToolsMeasureArea_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),    return;      end
zoom_state(handles,'maybe_off');
[xp,yp] = getline_j(handles.figure1,'closed');    n_nodes = length(xp);
if (n_nodes < 2),   zoom_state(handles,'maybe_on');     return;     end
if (handles.geog)
    area = area_geo(yp, xp);        % Area is reported on the unit sphere
    area = area * 4 * pi * handles.EarthRad^2;
    unit = ' km^2'; 
else
    area = polyarea(xp,yp);   % Area is reported in map user unites
    unit = ' map units ^2'; 
end
zoom_state(handles,'maybe_on');
msgbox(['Area = ' num2str(area) unit])

% --------------------------------------------------------------------
function ToolsMeasureAzimuth_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),    return;      end
zoom_state(handles,'maybe_off');
msg = [];
[xp,yp] = getline_j(handles.figure1);      n_nodes = length(xp);
if (n_nodes < 2),   zoom_state(handles,'maybe_on');     return;     end
if (handles.geog)
    az = azimuth_geo(yp(1:end-1), xp(1:end-1), yp(2:end), xp(2:end));
else
    dx = diff(xp);   dy = diff(yp);
    angs = atan2(dy,dx) * 180/pi;       % and convert to degrees
    az = (90 - angs);                   % convert to azim (cw from north)
    ind = find(az < 0);
    az(ind) = 360 + az(ind);
end

for (i = 1:length(az))
    msg = [msg; {['Azimuth' num2str(i) '  =  ' num2str(az(i),'%3.1f') '  degrees']}];
end
zoom_state(handles,'maybe_on');
msgbox(msg,'Azimuth(s)')

% --------------------------------------------------------------------
function [X,Y,Z,head,m,n] = load_grd(handles, opt)
% Load a GMT grid either from memory, or re-read it again if it is to
% big to fit in it (biger than handles.grdMaxSize)
if (handles.image_type == 20),  return;     end     % Fake image with bg_color. Nothing loadable
if (nargin == 1),   opt = ' ';  end
X = getappdata(handles.figure1,'dem_x');    Y = getappdata(handles.figure1,'dem_y');
Z = getappdata(handles.figure1,'dem_z');    head = getappdata(handles.figure1,'GMThead');
if (handles.ForceInsitu),   opt_I = '-I';   % Use only in desperate cases.
else                        opt_I = ' ';    end
if (nargout > 4),           [m,n] = size(Z);end
err_msg = 0;

if (isempty(X) && handles.image_type == 1 && ~handles.computed_grid)
    [X,Y,Z,head] = grdread_m(handles.grdname,'single',opt_I);
    if (nargout == 6),          [m,n] = size(Z);end
elseif (isempty(X) && handles.image_type == 4 && ~handles.computed_grid)
    Z = [];     % Check for this to detect a (re)-loading error
    err_msg = 'Grid was not on memory. Increase "Grid max size" and start over again.';
elseif (isempty(X) && isempty(handles.grdname))
    Z = [];     % Check for this to detect a (re)-loading error
    err_msg = 'Grid could not be reloaded. You probably need to increase "Grid max size"';
end
if (err_msg & ~strcmp(opt,'silent')),    errordlg(err_msg,'ERROR');    end

% --------------------------------------------------------------------
function DrawLine_Callback(hObject, eventdata, handles, opt)
if (handles.no_file == 1),  return;     end
if (nargin == 3),           opt = [];   end
% The following is a necessary pach agains two big stupidities.
% First is from the user that double-clicked on the line icon
% Second is from Matlab that doesn't lets us test a double-click on a uipushtool
if ( ~isempty(getappdata(handles.figure1, 'FromGetLine_j')) )
    return
end
zoom_state(handles,'maybe_off');


if (~isempty(opt) && strcmp(opt,'freehand')),   [xp,yp] = getline_j(handles.figure1,'freehand');
else                                            [xp,yp] = getline_j(handles.figure1);    end
n_nodes = length(xp);
if (n_nodes < 2),   zoom_state(handles,'maybe_on');     return;     end
% The polyline Tag is very important to destinguish from MB tracks, which have Tags = MBtrack#
if (~isempty(opt) && strcmp(opt,'FaultTrace'))    % When this function is used for Okada modeling
    if (handles.image_type ~= 1),   aux_funs('msg_dlg',2,handles);     zoom_state(handles,'maybe_on'); return;      end
    lineHand = line('XData', xp, 'YData', yp,'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,'Tag','FaultTrace');
    % Create empty patches that will contain the surface projection of the fault plane
    hp = zeros(length(xp)-1);
    for (k=1:length(xp)-1), hp(k) = patch('XData', [], 'YData',[]);    end
    setappdata(lineHand,'PatchHand',hp);
else
    if (xp(1) == xp(end) && yp(1) == yp(end))    % If line was close by hiting 'c'
        lineHand = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',handles.DefLineColor,...
        'LineWidth',handles.DefLineThick,'Tag','Closedpolygon');
    else
        lineHand = line('XData', xp, 'YData', yp,'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,'Tag','polyline');
    end
end
draw_funs(lineHand,'line_uicontext')        % Set lines's uicontextmenu
zoom_state(handles,'maybe_on');

% --------------------------------------------------------------------
function DrawVector_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
zoom_state(handles,'maybe_off');
draw_funs([],'DrawVector')      % Vectors are drawn there
zoom_state(handles,'maybe_on');

% --------------------------------------------------------------------
function DrawClosedPolygon_Callback(hObject, eventdata, handles, opt)
% OPT = [] Draw a closed polygon and exit
% OPT = 'EulerTrapezium' Draw a trapezium (only 4 sides => 5 pts). Used to fast compute an Euler pole
% OPT = 'SeismicityPolygon' Draw a closed polygon to which special seismicity opts will be added.
% OPT = 'from_ROI' Draw a closed polygon and call the ROI operations window
% OPT = 'rectangle' Draw a rectangle
% OPT = h (where h is a line handle) Calls the ROI operations window to operate on the polygon
%       whose handle is h. This mode is activated from a uicontextmenu that closed polylines share
if (handles.no_file == 1),      return;      end
if (nargin == 3),   opt = [];   end
if (ishandle(opt)), h_line = opt;   opt = 'from_uicontext';  end
% The following is a necessary pach agains two big stupidities.
% First is from the user that double-clicked on the polyline icon
% Second is from Matlab that doesn't lets us test a double-click on a uipushtool
if ( ~isempty(getappdata(handles.figure1, 'FromGetLine_j')) )
    return
end

if ( isempty(opt) || any(strcmp(opt,{'from_ROI' 'EulerTrapezium' 'SeismicityPolygon'})) )
    zoom_state(handles,'maybe_off');
    [xp,yp] = getline_j(handles.figure1,'closed');      n_nodes = length(xp);
    if (n_nodes < 4),    return;     end        % 4 because a straight line has 3 vertex (last one repeats)
    if (strcmp(opt,'EulerTrapezium') && n_nodes ~= 5)
        errordlg('OK, I won''t insult you this time. Just RTFM or don''t use this option','Error');     return;
    end
    if (strcmp(opt,'EulerTrapezium')),  tag = 'EulerTrapezium';
    elseif (strcmp(opt,'SeismicityPolygon')),  tag = 'SeismicityPolygon';
    else                                tag = 'Closedpolygon';      end
    zoom_state(handles,'maybe_on');
    xp = xp(:)';    yp = yp(:)';
    % The following Tag is very important to distinguish from MB tracks, which have Tags = MBtrack#
    lineHand = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',handles.DefLineColor,...
        'LineWidth',handles.DefLineThick,'Tag',tag);
    draw_funs(lineHand,'line_uicontext')      % Set lines's uicontextmenu
    if (isempty(opt) || any(strcmp(tag,{'EulerTrapezium' 'SeismicityPolygon'}))),  return;     end      % We are done in this mode (just draw a closed polygon)
    % If we come here it's because ROI operations were chosen
    roi_image_operations(handles.axes1,[xp(:),yp(:)])
elseif (strcmp(opt,'from_uicontext'))
    xp = get(h_line,'XData');   yp = get(h_line,'YData');
    roi_image_operations(handles.axes1,[xp(:),yp(:)])
elseif (strcmp(opt,'rectangle'))
    zoom_state(handles,'maybe_off');
    [p1,p2,hl] = rubberbandbox;
    set(hl,'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick)    % Use defaults LineThick and DefLineColor
    draw_funs(hl,'line_uicontext')       % Set lines's uicontextmenu
    zoom_state(handles,'maybe_on');
end

% --------------------------------------------------------------------
function DrawEulerPoleCircle_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
zoom_state(handles,'maybe_off');

out = euler_poles_selector;         % The output is a struct with fields: lon lat omega plates model
if isempty(out),    return;     end % User gave up
lon = out.lon;      lat = out.lat;
h_circ = uicirclegeo(lon,lat);
set(h_circ,'Tag','CircleEuler')     % This is used by draw_funs to allow velocity computations
if ~isempty(out)
    s = get(h_circ,'Userdata');
    s.omega = out.omega;
    if ~isempty(out.plates)     % Just in case
        if isempty(strmatch('absolute',out.plates))     % A relative plate model
            s.plates = [out.plates '  -- Model = ' out.model];
        else                                            % An absolute plate model
            s.plates = [out.plates(end-1:end) ' -- Model = ' out.model ' (Absolute)'];
        end
    else
        s.plates = 'I''m lost';
    end
    set(h_circ,'Userdata',s)
end
draw_funs(h_circ,'SessionRestoreCircle')
zoom_state(handles,'maybe_on');

% --------------------------------------------------------------------
function DrawGeographicalCircle_Callback(hObject, eventdata, handles, opt)
% If OPT exist (no matter what it contains) draw a great circle arc instead of a circle
if (handles.no_file == 1),      return;      end
if (nargin == 3),   opt = [];   end
zoom_state(handles,'maybe_off');
if (handles.geog && strcmp(opt,'gcirc'))
    draw_funs([],'DrawGreatCircle')             % All work is done there
elseif (handles.geog ~= 1 && strcmp(opt,'gcirc'))
    warndlg('Great Circles are only programed to work with geog coordinates.','Warning')
elseif (handles.geog ~= 1 && isempty(opt) || strcmp(opt,'cartCirc'))
    draw_funs([],'DrawCartesianCircle')         % All work is done there
else
    h_circ = uicirclegeo;                       % Draw the circle and associate Button... functions
    draw_funs(h_circ,'SessionRestoreCircle')    % Give uicontext
end
zoom_state(handles,'maybe_on');

% --------------------------------------------------------------------
function DrawImportLine_Callback(hObject, eventdata, handles, opt)
% OPT is a string with either "AsLine", or "AsPoint", or "AsMaregraph", or "FaultTrace"
if (handles.no_file == 1),      aux_funs('msg_dlg',1,handles);     return;      end
str1 = {'*.dat;*.DAT', 'Data file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select input xy file name','get');
if isequal(FileName,0);     return;     end
draw_funs([PathName FileName],'ImportLine',opt)

% --------------------------------------------------------------------
function DrawImportText_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),  return;      end
str1 = {'*.txt;*.TXT;*.dat;*.DAT', 'Text file (*.txt,*.TXT,*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select input text file name','get');
if isequal(FileName,0);     return;     end

fid = fopen([PathName,FileName],'r');
todos = fread(fid,'*char');     fclose(fid);
[str.x str.y str.name] = strread(todos,'%f %f %s');     % Note: text.name is a cell array of chars
limits = getappdata(handles.figure1,'ThisImageLims');   xx = limits(1:2);   yy = limits(3:4);
% Get rid of Texts that are outside the map limits
indx = find(str.x < xx(1) | str.x > xx(2));
str.x(indx) = [];      str.y(indx) = [];      str.name(indx) = [];
indx = find(str.y < yy(1) | str.y > yy(2));
str.x(indx) = [];      str.y(indx) = [];      str.name(indx) = [];

n_str = length(str.x);
if (n_str == 0),   return;     end             % No texts inside area. Return.

str.name = strrep(str.name,'_',' ');          % Replace '_' by ' '
h = text(str.x(1),str.y(1),str.name{1},'FontName','Book Antiqua','FontWeight','bold');
draw_funs(h,'DrawText')
h = text(str.x(2),str.y(2),str.name{2},'FontSize',18,'FontWeight','bold');
draw_funs(h,'DrawText')

% --------------------------------------------------------------------
function DrawImportShape_Callback(hObject, eventdata, handles)
str1 = {'*.shp;*.SHP', 'Data files (*.shp,*.SHP)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select shape file name','get');
if isequal(FileName,0);     return;     end
[s,t] = mex_shape([PathName FileName]);

lt = getappdata(handles.figure1,'DefLineThick');    lc = getappdata(handles.figure1,'DefLineColor');
region = [s(1).BoundingBox(1,1:2) s(1).BoundingBox(2,1:2)];
is_geog = guessGeog(region(1:4));

if (handles.geog && ~is_geog)
    errordlg('Error. Your background image is in geographics but the shape file is not','ERROR');   return
elseif (handles.geog == 0 && is_geog == 1 && ~isempty(handles.grd_img))
    warndlg('WARNING: Your background image is not in geographics but the shape file seams to be. Probable mistake','Warning')
else
    handles.geog = is_geog;     guidata(handles.figure1, handles);
end

if (handles.no_file)
    FileNewBgFrame_Callback(handles.figure1,[],handles, [region handles.geog]);
end

h = zeros(length(s),1);
imgLims = getappdata(handles.axes1,'ThisImageLims');
if (strcmp(t,'Arc'))
    for i = 1:length(s)
        out = aux_funs('insideRect',imgLims,[s(i).BoundingBox(1,1) s(i).BoundingBox(2,1); s(i).BoundingBox(1,1) s(i).BoundingBox(2,2); ...
            s(i).BoundingBox(1,2) s(i).BoundingBox(2,2); s(i).BoundingBox(1,2) s(i).BoundingBox(2,1)]);
        if (any(out))       % It means the polyg BB is at least partially inside
            h(i) = line('Xdata',s(i).X,'Ydata',s(i).Y,'Parent',handles.axes1,'Color',lc,'LineWidth',lt,'Tag','SHPpolyline');
        end
        h((h == 0)) = [];   % Those were jumped because thay were completely outside map limits
    end
    draw_funs(h,'SHP_uicontext')        % Set lines's uicontextmenu
elseif (strcmp(t,'Polygon'))
    for i = 1:length(s)
        out = aux_funs('insideRect',imgLims,[s(i).BoundingBox(1,1) s(i).BoundingBox(2,1); s(i).BoundingBox(1,1) s(i).BoundingBox(2,2); ...
            s(i).BoundingBox(1,2) s(i).BoundingBox(2,2); s(i).BoundingBox(1,2) s(i).BoundingBox(2,1)]);
        if (any(out))       % It means the polyg BB is at least partially inside
            h(i) = patch('XData',s(i).X,'YData', s(i).Y,'FaceColor','none','EdgeColor',handles.DefLineColor,'Tag','SHPpolygon');
        end
        h((h == 0)) = [];   % Those were jumped because thay were completely outside map limits
        %draw_funs(h(i),'line_uicontext')           % In future do this when they are not to many many
    end
    draw_funs(h,'country_patch')
end

% --------------------------------------------------------------------
function GeophysicsImportFaultFile_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),      aux_funs('msg_dlg',1,handles);     return;      end
str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select input xy file name','get');
if isequal(FileName,0);     return;     end
[numeric_data,multi_segs_str] = text_read([PathName FileName],NaN,NaN);

% fac = 1 / 6371 * 180 / pi;      % Desenrasque pro Miguel
% for (k=1:size(numeric_data,1))
%     lon1 = numeric_data(k,1);    lat1 = numeric_data(k,2);
%     [lat2(k),lon2(k)] = circ_geo(lat1,lon1,numeric_data(k,3)*fac,numeric_data(k,6)+0,1);
% end
% numeric_data(:,1) = lon2';
% numeric_data(:,2) = lat2';
% for (k=1:size(numeric_data,1))
%     if (k == 1),    lon1 = numeric_data(k,1);    lat1 = numeric_data(k,2);
%     else            lon1 = x(4);    lat1 = y(4);
%     end
%     [lat2,lon2] = circ_geo(lat1,lon1,numeric_data(k,3)*fac,numeric_data(k,6)+180,1);
%     lat = [lat1;lat2];      lon = [lon1;lon2];
%     %lineHand1(k) = plot(lon,lat,'Color','r','Tag','FaultTrace');
%     lineHand1(k) = line('XData',lon,'YData',lat,'Color','r','Tag','FaultTrace');
%     draw_funs(lineHand1(k),'line_uicontext')     % Set lines's uicontextmenu
%     depth{k} = numeric_data(k,5);      slip{k} = numeric_data(k,9);
%     rake{k} = numeric_data(k,8);       strike{k} = numeric_data(k,6)+180;
%     dip{k} = numeric_data(k,7);
%     % Draw patches
%     rng_p = numeric_data(k,4) * cos(dip{k}*pi/180) * fac;
%     [lat1,lon1] = circ_geo(lat(1),lon(1),rng_p,strike{k}+90,1);
%     [lat2,lon2] = circ_geo(lat(2),lon(2),rng_p,strike{k}+90,1);
%     x = [lon(1) lon(2) lon2 lon1];    y = [lat(1) lat(2) lat2 lat1];
%     hp = patch('XData',x,'YData',y,'FaceColor',[.8 .8 .8]);
%     setappdata(lineHand1(k),'PatchHand',hp);
% end
% deform_okada(handles.figure1,handles.geog,lineHand1,strike,depth,slip,dip,rake);

nf = 15;     nseg = 30;
rng = repmat(7.5 / 6371 * 180 / pi, nseg, 1);
hold on
for (k=1:nf)
    tmpx = numeric_data((k-1)*nseg+1:k*nseg,2);       tmpy = numeric_data((k-1)*nseg+1:k*nseg,1);
    [lat1,lon1] = circ_geo(tmpy,tmpx,rng,numeric_data((k-1)*nseg+1:k*nseg,6)+180,1);
    [lat2,lon2] = circ_geo(tmpy(end),tmpx(end),rng(end),numeric_data((k-1)*nseg+1,6),1);
    lat = [lat1;lat2];      lon = [lon1;lon2];
    lineHand1(k) = plot(lon,lat,'Color','r','Tag','FaultTrace');
    draw_funs(lineHand1(k),'line_uicontext')     % Set lines's uicontextmenu
    %Lon. Lat. depth slip rake strike dip
    depth{k} = numeric_data((k-1)*nseg+1:k*nseg,3);      slip{k} = numeric_data((k-1)*nseg+1:k*nseg,4);
    rake{k} = numeric_data((k-1)*nseg+1:k*nseg,5);       strike{k} = numeric_data((k-1)*nseg+1:k*nseg,6);
    dip{k} = numeric_data((k-1)*nseg+1:k*nseg,7);
    for (i=1:nseg)      % Draw patches
        rng_p = 12 * cos(dip{k}(i)*pi/180) / 6371 * 180 / pi;
        [lat1,lon1] = circ_geo(lat(i),lon(i),rng_p,strike{k}(i)+90,1);
        [lat2,lon2] = circ_geo(lat(i+1),lon(i+1),rng_p,strike{k}(i)+90,1);
        x = [lon(i) lon(i+1) lon2 lon1];    y = [lat(i) lat(i+1) lat2 lat1];
        hp(i) = patch('XData',x,'YData',y,'FaceColor',[.8 .8 .8]);
    end
    setappdata(lineHand1(k),'PatchHand',hp);
end
hold off
deform_okada(handles.figure1,handles.geog,lineHand1,strike,depth,slip,dip,rake);

% --------------------------------------------------------------------
function GeophysicsImportGmtFile_Callback(hObject, eventdata, handles)
% Open a .gmt file
if (handles.no_file == 0 && handles.geog ~= 1),  aux_funs('msg_dlg',3,handles);     return;      end

[FileName,PathName] = put_or_get_file(handles,{'*.gmt;*.GMT', 'gmt files (*.gmt,*.GMT)'},'Select gmt File','get');
if isequal(FileName,0);     return;     end
[PATH,FNAME] = fileparts([PathName FileName]);
f_name = [PATH filesep FNAME];          % Rip the .gmt extension
set(handles.figure1,'Pointer','watch')
track = gmtlist_m(f_name,'-Fxy','-G');
if (handles.no_file)        % We don't have a BG map, so we have to create one
    y_min = min(track.latitude);    y_max = max(track.latitude);
    x_min = min(track.longitude);   x_max = max(track.longitude);
    FileNewBgFrame_Callback(handles.figure1,[],handles, [x_min x_max y_min y_max 1])
else                        % Get rid of points that are outside the map limits
    [track.longitude,track.latitude] = aux_funs('in_map_region',handles,track.longitude,track.latitude,0.1,[]);
    if (isempty(track.longitude))       % There was nothing left
        set(handles.figure1,'Pointer','arrow');     return
    end
end
h = line(track.longitude,track.latitude,'Linewidth',handles.DefLineThick,'Color',handles.DefLineColor,'Tag',FileName,'Userdata',1);
setappdata(h,'FullName',[PathName FileName])    % Store file name in case the uicontext wants to open it with gmtedit
draw_funs(h,'gmtfile',track.info)
set(handles.figure1,'Pointer','arrow')

% --------------------------------------------------------------------
function GeophysicsImportGmtFileList_Callback(hObject, eventdata, handles)
% Open a list of .gmt files
if (handles.no_file == 0 && handles.geog ~= 1),   aux_funs('msg_dlg',3,handles);     return;      end
str1 = {'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.DAT,*.txt,*.TXT)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select list file','get');
if isequal(FileName,0);     return;     end
fid = fopen([PathName FileName]);
c = char(fread(fid))';      fclose(fid);
names = strread(c,'%s','delimiter','\n');   clear c fid;
% Test if any file in the list does not exist
m = length(names);      nao = zeros(m,1);
for (k=1:m),    nao(k) = exist(names{k});   end
id = find(nao ~= 2);
if (~isempty(id))
    msgbox(names{id},'FILES NOT FOUND');
    names{id} = [];             % Remove them from the list
    if (isempty(names)),    return;     end;        % empty list
end
names_ui = names;               % For 'Tagging' lines purpose
for (k=1:length(names))         % Rip the .gmt extension
    [PATH,FNAME] = fileparts(names{k});
    names{k} = FNAME;       names_ui{k} = FNAME;
    if (~isempty(PATH))     % File names in the list have a path
        names{k} = [PATH filesep names{k}];
    else                    % They do not have a path, but we need it. So prepend the list-file path
        names{k} = [PathName names{k}];
    end
end

set(handles.figure1,'Pointer','watch');
if (handles.no_file)        % We don't have a BG map, so we have to create one
    track = gmtlist_m(names{1:end},'-Fxy','-G');
    len_t = length(track);  x_min = zeros(1,len_t); x_max = x_min;  y_min = x_min;  y_max = x_min;
    for (k=1:length(track))
        x_min(k) = min(track(k).longitude);     x_max(k) = max(track(k).longitude);
        y_min(k) = min(track(k).latitude);      y_max(k) = max(track(k).latitude);
    end
    x_min = min(x_min);    x_max = max(x_max);
    y_min = min(y_min);    y_max = max(y_max);
    FileNewBgFrame_Callback(handles.figure1,[],handles, [x_min x_max y_min y_max 1])
else
    x_lim = get(handles.axes1,'XLim');    y_lim = get(handles.axes1,'YLim');
    opt_R = ['-R' num2str(x_lim(1)) '/' num2str(x_lim(2)) '/' num2str(y_lim(1)) '/' num2str(y_lim(2))];
    track = gmtlist_m(names{1:end},'-Fxy','-G',opt_R);
end

% And finaly do the ploting
colors = rand(length(track),3);         % Use a random color schema
for (k=1:length(track))
    id0 = (track(k).longitude == 0);    % I must change gmtlist to do this
    track(k).longitude(id0) = [];       track(k).latitude(id0) = [];
    if (isempty(track(k).longitude)),   continue;   end     % This track is completely outside the map
	h = line(track(k).longitude,track(k).latitude,'Linewidth',handles.DefLineThick,'Color',...
        colors(k,:),'Tag',names_ui{k},'Userdata',1);
	setappdata(h,'FullName',names{k})    % Store file name in case the uicontext wants to open it with gmtedit
	draw_funs(h,'gmtfile',track(k).info)
end
set(handles.figure1,'Pointer','arrow');

% --------------------------------------------------------------------
function DrawText_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),    return;      end
zoom_state(handles,'maybe_off');
pt = ginput_pointer(1,'crosshair');
textHand = text(pt(1),pt(2),0,'','Editing','on');
draw_funs(textHand,'DrawText')          % Set text's uicontextmenu
zoom_state(handles,'maybe_on');

% --------------------------------------------------------------------
function DrawSymbol_Callback(hObject, eventdata, handles, opt)
if (handles.no_file == 1),    return;      end
switch opt
    case 'circle',   smb = 'o';
    case 'square',   smb = 's';
    case 'triangle', smb = '^';
    case 'star',     smb = 'p';
    case 'cross',    smb = 'x';
end
zoom_state(handles,'maybe_off');
pt = ginput_pointer(1,'crosshair');
symbHand = line(pt(1,1),pt(1,2),'Marker',smb,'MarkerFaceColor','y',...
        'MarkerEdgeColor','k','MarkerSize',10,'Tag','Symbol');
draw_funs(symbHand,'DrawSymbol')          % Set symbol's uicontextmenu
zoom_state(handles,'maybe_on');

% --------------------------------------------------------------------
function DrawContours_Callback(hObject, eventdata, handles, opt)
if (aux_funs('msg_dlg',14,handles));     return;      end
if (nargin == 3),   opt = [];   end
set(handles.figure1,'pointer','watch')
[X,Y,Z,head] = load_grd(handles);
if isempty(Z),  return;     end;    % An error message was already issued
if (isempty(opt))               % Do the automatic contouring
    c = contourc(X,Y,double(Z));
elseif (isa(opt,'char'))        % Call the interface contouring GUI
    h_which_cont = findobj(handles.figure1,'Type','line','Tag','contour');     % See if contours were mouse deleted
    if (~isempty(h_which_cont))
        handles.which_cont = unique(cell2mat(get(h_which_cont,'Userdata')));
    else
        handles.which_cont = [];
    end
    set(handles.figure1,'pointer','arrow');     contouring(handles.figure1,head,handles.which_cont);
    guidata(handles.figure1, handles);
    return
elseif (isa(opt,'double'))      % Arrive here from the interface contouring GUI
    if (~isempty(handles.which_cont))   % Do not compute/draw repeated contours
        [c,ib] = setdiff(handles.which_cont,opt);   % Find eventual countours to remove (by GUI deselection)
        for (i = 1:length(c))   % Loop over removing contours (if none, this loop has no effect)
            h = findobj(handles.axes1,'type','line','userdata',handles.which_cont(ib(i)));
            for (j = 1:length(h))   % We can easily have more than one
                labHand = getappdata(h(j),'LabelHands');
                try    delete(labHand);   end  % Delete contour's labels
            end
            delete(h)                           % And delete the selected contours
        end
        handles.which_cont(ib) = [];
        guidata(handles.figure1, handles);
        [c,ia,ib] = intersect(handles.which_cont,opt(:));
        opt(ib) = [];                   % Remove repeated contours
    end    
    if (isempty(opt)),      set(handles.figure1,'pointer','arrow');    return;     end  % Nothing else to do
    if (length(opt) == 1),  opt = [opt opt];    end
    c = contourc(X,Y,double(Z),opt);
    if (isempty(c)),        set(handles.figure1,'pointer','arrow');    return;     end
end
limit = size(c,2);
i = 1;      h_cont = [];      cont = [];
while(i < limit)
    z_level = c(1,i);    npoints = c(2,i);
    nexti = i+npoints+1;
    xdata = c(1,i+1:i+npoints);     ydata = c(2,i+1:i+npoints);
    % Create the lines
    cu = line('XData',xdata,'YData',ydata,'LineWidth',1,'Color','k','userdata',z_level,'Tag','contour');
    h_cont = [h_cont; cu(:)];
    cont = [cont; z_level];
    i = nexti;
end
handles.which_cont = unique([handles.which_cont; cont]);
[zlev, ind] = sort(cont);       clear zlev
h_cont = h_cont(ind);           % handles are now sorted by level
h_label = clabel_j(c,h_cont);   % Label countours automatically
set(h_label,'Tag','contour');   % The tag is used in "Delete all contours" to delete also the labels
for i = 1:length(h_cont)        % Set convenient uicontexts. One for each contour
    setappdata(h_cont(i),'cont_label',get(h_cont(i),'UserData'))
    draw_funs(h_cont(i),'ContourLines',h_label)
end
set(handles.figure1,'pointer','arrow')
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function FileOpenSession_Callback(hObject, eventdata, handles)
str1 = {'*.mat;*.MAT', 'Data files (*.mat,*.MAT)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select session file name','get');
if isequal(FileName,0);     return;     end

set(handles.figure1,'pointer','watch')
load([PathName FileName])

if isempty(grd_name)
    if isempty(map_limits)      % Neither grid or map_limits (???)
        set(handles.figure1,'pointer','arrow');    errordlg('Unknown error in "Open Session": Quiting','Error');
        return
    end
    scrsz = get(0,'ScreenSize');         % Get screen size
    dx = map_limits(2) - map_limits(1);   dy = map_limits(4) - map_limits(3);
    aspect = dy / dx;
    nx = round(scrsz(3)*.75);       ny = round(nx * aspect);
    if (ny > scrsz(4) - 30)
        ny = scrsz(4) - 30;         nx = round(ny / aspect);
    end
    Z = repmat(uint8(255),ny,nx);       % Create a white image
    X = [map_limits(1) map_limits(2)];      Y = [map_limits(3) map_limits(4)];
    x_inc = diff(X) / nx;                   y_inc = diff(Y) / ny;
    dx2 = x_inc / 2;                        dy2 = y_inc / 2;
    X = X + [dx2 -dx2];                     Y = Y + [dy2 -dy2];     % Make it such that the pix-reg info = region
    handles.head = [X Y 0 255 0 x_inc y_inc];    
    handles.image_type = 20;
    if (~isempty(img_pal)),     set(handles.figure1,'Colormap',img_pal)   % Patch for a stupid users use
    else                        set(handles.figure1,'Colormap', ones( size(get(handles.figure1,'Colormap'),1), 3));
    end
    handles = show_image(handles,'Mirone Base Map',X,Y,Z,0,'xy',1);
elseif ~strcmp(handles.grdname,grd_name)    % Session's grid is different than current grid, so we have to open "grd_name"
    handles.image_type = 1;
    [X,Y,Z,head] = grdread_m(grd_name,'single');
    handles.geog = guessGeog(head(1:4));
    handles.grdname = grd_name;     handles.head = head;     % Save header info (needed ?)
    zz = scaleto8(Z);
    handles.have_nans = grdutils(Z,'-N');
    % If grid size is not to big I'll store it in memory
    if (numel(zz)*4 <= handles.grdMaxSize)
        setappdata(handles.figure1,'dem_z',Z);  setappdata(handles.figure1,'dem_x',X);
        setappdata(handles.figure1,'dem_y',Y);  setappdata(handles.figure1,'GMThead',head);
    end
    
    try         aux_funs('colormap_bg',handles,Z,img_pal);
    catch       aux_funs('colormap_bg',handles,Z,jet(256));
    end
    handles = show_image(handles,grd_name,X,Y,zz,1,'xy',head(7));
else
    err = 0;
    try         soma = sum(img_pal - get(handles.figure1,'Colormap'));
    catch       err = 1;    soma = 0;        end
    if (soma > 0.1 | err == 1)      % grid is the same but colormap is different
        set(handles.figure1,'Colormap',img_pal)
    end
end 
h_figs = findobj('Type','figure');
if (length(h_figs) > 1)                                 % Large images generate a warning message figure
    h_kill = findobj(h_figs,'Tag','Msgbox_Warning');    % that has to be killed, otherwise the drawing
    delete(h_kill)                                      % elements are put on it and not on the image.
end

hold on
if (haveMBtrack)                % case of MB tracks
    for i=1:length(MBtrack)
        h_line = line(MBtrack(i).x,MBtrack(i).y,'LineWidth',MBtrack(i).LineWidth,'color',...
            MBtrack(i).color,'Tag',MBtrack(i).tag, 'LineStyle',MBtrack(i).LineStyle);
        setappdata(h_line,'swathRatio',MBtrack(i).swathRatio)
        draw_funs(h_line,'MBtrack_uicontext')       % Set track's uicontextmenu
    end
    for i=1:length(MBbar)       % now their's bars
        h_bar = line(MBbar(i).x,MBbar(i).y,'LineWidth',MBbar(i).LineWidth,'color',MBbar(i).color,...
            'Tag',MBbar(i).tag,'UserData',MBbar(i).n_vert, 'LineStyle',MBbar(i).LineStyle);
        draw_funs(h_bar,'MBbar_uicontext')         % Set track bar's uicontextmenu
    end
    handles.h_MBplot = h_line;
end
if (haveCircleGeo)              % case of Geographic circles
    for i=1:length(CircleGeo)
        h_circ = line(CircleGeo(i).x,CircleGeo(i).y,'LineWidth',CircleGeo(i).LineWidth,'color',...
            CircleGeo(i).color,'Tag',CircleGeo(i).tag, 'LineStyle',CircleGeo(i).LineStyle);
        setappdata(h_circ,'LonLatRad',CircleGeo(i).lon_lat_rad);
        CircleGeo(i).ud.hcirc = h_circ;
        CircleGeo(i).ud.parent = get(h_circ,'parent');
        set(h_circ,'UserData',CircleGeo(i).ud,'buttondownfcn','uicirclegeo(''circlemousedown'')')
        draw_funs(h_circ,'SessionRestoreCircle')       % Set circle's uicontextmenu
    end
end
if (haveCircleCart)             % case of Cartesian circles
    for i=1:length(CircleCart)
        h_circ = line(CircleCart(i).x,CircleCart(i).y,'LineWidth',CircleCart(i).LineWidth,'color',...
            CircleCart(i).color,'Tag',CircleCart(i).tag, 'LineStyle',CircleCart(i).LineStyle);
        setappdata(h_circ,'LonLatRad',CircleCart(i).lon_lat_rad);
        CircleCart(i).ud.hcirc = h_circ;
        CircleCart(i).ud.parent = get(h_circ,'parent');
        draw_funs(h_circ,'SessionRestoreCircleCart')       % Set circle's uicontextmenu
    end
end
if (havePline)                  % case of polylines
    for i=1:length(Pline)
        h_line = line(Pline(i).x,Pline(i).y,'LineWidth',Pline(i).LineWidth,'color',...
            Pline(i).color,'Tag',Pline(i).tag, 'LineStyle',Pline(i).LineStyle);
        draw_funs(h_line,'line_uicontext')       % Set lines's uicontextmenu
    end
end
if (havePlineAsPoints)          % case of polylines as points (markers) only
    for i=1:length(PlineAsPoints)
        h_line_pt = plot(PlineAsPoints(i).x, PlineAsPoints(i).y, 'LineStyle','none', ...
            'Marker',PlineAsPoints(i).Marker, 'MarkerSize',PlineAsPoints(i).Size, ...
            'MarkerFaceColor',PlineAsPoints(i).FillColor, ...
            'MarkerEdgeColor',PlineAsPoints(i).EdgeColor, 'Tag',PlineAsPoints(i).tag);
        draw_funs(h_line_pt,'DrawSymbol')        % Set marker's uicontextmenu (tag is very important)        
    end
end
if (haveSymbol)                 % case of Symbols (line Markers)
    for i=1:length(Symbol)
        h_symb = line(Symbol(i).x,Symbol(i).y,'Marker',Symbol(i).Marker,'MarkerSize',Symbol(i).Size,...
            'MarkerFaceColor',Symbol(i).FillColor, 'MarkerEdgeColor',Symbol(i).EdgeColor, 'Tag',Symbol(i).tag);
        draw_funs(h_symb,'DrawSymbol')          % Set symbol's uicontextmenu
    end
end
if (haveText)                   % case of text strings
    for i=1:length(Texto)
        if (isempty(Texto(i).str)),  continue;   end
        h_text = text(Texto(i).pos(1),Texto(i).pos(2),Texto(i).pos(3), Texto(i).str, 'Rotation',Texto(i).angle,...
            'FontAngle',Texto(i).FontAngle, 'Tag',Texto(i).Tag, 'FontWeight',Texto(i).FontWeight,...
            'color',Texto(i).color, 'FontName',Texto(i).FontName, 'FontSize',Texto(i).FontSize);
        if (isfield(Texto(i),'VerticalAlignment')),      set(h_text,'VerticalAlignment',Texto(i).VerticalAlignment);        end
        if (isfield(Texto(i),'HorizontalAlignment')),    set(h_text,'HorizontalAlignment',Texto(i).HorizontalAlignment);    end
        draw_funs(h_text,'DrawText')       % Set texts's uicontextmenu
    end
end
if (havePatches)                % case of patchs - NOTE, the Tags are currently lost
    for (i=1:length(Patches))
        try                     % We never know with those guys, so it's better to play safe
            is_telha = 0;
            if (strcmp(Patches(i).tag,'tapete_R') || strcmp(Patches(i).tag,'tapete'))
                Patches(i).x = reshape(Patches(i).x,4,length(Patches(i).x)/4);
                Patches(i).y = reshape(Patches(i).y,4,length(Patches(i).y)/4);
                is_telha = 1;
            end
%             if (Patches(i).FaceColor(1) == 1)
%                 Patches(i).tag = 'tapete';
%             else
%                 Patches(i).tag = 'tapete_R';
%             end
            h_patch = patch('XData',Patches(i).x, 'YData',Patches(i).y, 'LineWidth',Patches(i).LineWidth,...
                'EdgeColor',Patches(i).EdgeColor, 'FaceColor',Patches(i).FaceColor,...
                'LineStyle',Patches(i).LineStyle, 'Tag', Patches(i).tag);
            if (is_telha)
                set(h_patch,'UserData',Patches(i).ud)
                draw_funs(h_patch,'telhas_patch')       % Set telhas's uicontextmenu
            else
                draw_funs(h_patch,'line_uicontext')     % Set patch's uicontextmenu
            end
        end
    end
end
if (haveCoasts),    DatasetsCoastLineNetCDF_Callback([], [], handles,coastUD);  end
if (havePolitic)
    if (iscell(politicUD)),     politicUD = politicUD{1};     end
    DatasetsPoliticalBoundaries_Callback([], [], handles,politicUD(2),politicUD(1));
end
if (haveRivers)
    if (iscell(riversUD)),      riversUD = riversUD{1};     end
    DatasetsRivers_Callback([], [], handles,riversUD(2),riversUD(1));
end
guidata(hObject, handles);      hold off
set(handles.figure1,'pointer','arrow')

% ------------------------------------------------------------------------------------------------
function FileSaveSession_Callback(hObject, eventdata, handles)
if (handles.image_type == 0),   return;      end
str1 = {'*.mat;*.MAT', 'Data files (*.mat,*.MAT)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select session file name','put');
if isequal(FileName,0);     return;     end

set(handles.figure1,'pointer','watch')
[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if isempty(EXT),    fname = [PathName FNAME '.mat'];
else                fname = [PathName FNAME EXT];       end

grd_name = handles.grdname;     img_pal = get(handles.figure1,'Colormap');     map_limits = [];
if isempty(grd_name)        % Even if map does not come from a gmt grid allow at least the elements recovery
    map_limits = getappdata(handles.figure1,'ThisImageLims');
    img_pal = [];           % We will recreate a white bg image
end
ALLlineHand = findobj(get(handles.axes1,'Child'),'Type','line');
j = 1;  k = 1;  m = 1;  n = 1;  cg = 1; cc = 1; pp = 1;
haveMBtrack = 0;    havePline = 0;          haveText = 0;   haveSymbol = 0; haveCircleGeo = 0;
haveCircleCart = 0; havePlineAsPoints = 0;  havePatches = 0;haveCoasts = 0; havePolitic = 0;    haveRivers = 0;
MBtrack = [];   MBbar = [];     Pline = [];     Symbol = [];    Texto = [];     CircleGeo = [];
CircleCart = [];    PlineAsPoints = [];         Patches = [];   coastUD = [];   politicUD = []; riversUD = [];
for i = 1:length(ALLlineHand)
    tag = get(ALLlineHand(i),'Tag');
    if (strcmp(tag,'MBtrack'))       % case of a MBtrack line
        xx = get(ALLlineHand(i),'XData');     yy = get(ALLlineHand(i),'YData');
        MBtrack(j).x = xx(:);       MBtrack(j).y = yy(:);
        MBtrack(j).LineWidth = get(ALLlineHand(i),'LineWidth');
        MBtrack(j).LineStyle = get(ALLlineHand(i),'LineStyle');
        MBtrack(j).color = get(ALLlineHand(i),'color');
        MBtrack(j).tag = tag;       MBtrack(j).swathRatio = getappdata(ALLlineHand(i),'swathRatio');
        j = j + 1;      haveMBtrack = 1;
    elseif (strcmp(tag,'swath_w'))   % case of a MBtrack's bar line
        xx = get(ALLlineHand(i),'XData');     yy = get(ALLlineHand(i),'YData');
        MBbar(k).x = xx(:);         MBbar(k).y = yy(:);
        MBbar(k).LineWidth = get(ALLlineHand(i),'LineWidth');
        MBbar(k).LineStyle = get(ALLlineHand(i),'LineStyle');
        MBbar(k).color = get(ALLlineHand(i),'color');
        MBbar(k).tag = tag;
        MBbar(k).n_vert = get(ALLlineHand(i),'UserData');
        k = k + 1;
    elseif (strcmp(tag,'Symbol'))    % case of a Symbol (in fact a line Marker)
        xx = get(ALLlineHand(i),'XData');     yy = get(ALLlineHand(i),'YData');
        Symbol(n).x = xx(:);          Symbol(n).y = yy(:);
        Symbol(n).Marker = get(ALLlineHand(i),'Marker');
        Symbol(n).Size = get(ALLlineHand(i),'MarkerSize');
        Symbol(n).FillColor = get(ALLlineHand(i),'MarkerFaceColor');
        Symbol(n).EdgeColor = get(ALLlineHand(i),'MarkerEdgeColor');
        Symbol(n).tag = get(ALLlineHand(i),'Tag');
        n = n + 1;      haveSymbol = 1;
    elseif (strcmp(tag,'circleGeo') || strcmp(tag,'CircleEuler'))   % circles are particular line cases
        xx = get(ALLlineHand(i),'XData');     yy = get(ALLlineHand(i),'YData');
        CircleGeo(cg).x = xx(:);      CircleGeo(cg).y = yy(:);
        CircleGeo(cg).LineWidth = get(ALLlineHand(i),'LineWidth');
        CircleGeo(cg).LineStyle = get(ALLlineHand(i),'LineStyle');
        CircleGeo(cg).color = get(ALLlineHand(i),'color');
        CircleGeo(cg).tag = tag;
        CircleGeo(cg).lon_lat_rad = getappdata(ALLlineHand(i),'LonLatRad');
        CircleGeo(cg).ud = get(ALLlineHand(i),'UserData');   % UserData contains alot of need info
        cg = cg + 1;      haveCircleGeo = 1;
    elseif (strcmp(tag,'circleCart'))   % circles are particular line cases
        xx = get(ALLlineHand(i),'XData');     yy = get(ALLlineHand(i),'YData');
        CircleCart(cc).x = xx(:);      CircleCart(cc).y = yy(:);
        CircleCart(cc).LineWidth = get(ALLlineHand(i),'LineWidth');
        CircleCart(cc).LineStyle = get(ALLlineHand(i),'LineStyle');
        CircleCart(cc).color = get(ALLlineHand(i),'color');
        CircleCart(cc).tag = tag;
        CircleCart(cc).lon_lat_rad = getappdata(ALLlineHand(i),'LonLatRad');
        CircleCart(cc).ud = get(ALLlineHand(i),'UserData');
        cc = cc + 1;      haveCircleCart = 1;
%     elseif ~isempty(strfind(tag,'Arrow'))   % Arrows are tricky - they have two handles (line and head)
%         xx = get(ALLlineHand(i),'XData');     yy = get(ALLlineHand(i),'YData');
    elseif (strcmp(tag,'Pointpolyline'))   % Polyline with only markers are particular line cases
        xx = get(ALLlineHand(i),'XData');     yy = get(ALLlineHand(i),'YData');
        PlineAsPoints(pp).x = xx(:);      PlineAsPoints(pp).y = yy(:);
        PlineAsPoints(pp).Marker = get(ALLlineHand(i),'Marker');
        PlineAsPoints(pp).Size = get(ALLlineHand(i),'MarkerSize');
        PlineAsPoints(pp).FillColor = get(ALLlineHand(i),'MarkerFaceColor');
        PlineAsPoints(pp).EdgeColor = get(ALLlineHand(i),'MarkerEdgeColor');
        PlineAsPoints(pp).tag = tag;
        pp = pp + 1;      havePlineAsPoints = 1;
    elseif (strcmp(tag,'CoastLineNetCDF'))
        haveCoasts = 1;     coastUD = get(ALLlineHand(i),'UserData');
    elseif (strcmp(tag,'PoliticalBoundaries'))
        havePolitic = 1;    politicUD = get(ALLlineHand(i),'UserData');
    elseif (strcmp(tag,'Rivers'))
        haveRivers = 1;     riversUD = get(ALLlineHand(i),'UserData');
    else        % for the time beeing, it applyies to simple polylines
        xx = get(ALLlineHand(i),'XData');     yy = get(ALLlineHand(i),'YData');
        Pline(m).x = xx(:);         Pline(m).y = yy(:);
        Pline(m).LineWidth = get(ALLlineHand(i),'LineWidth');
        Pline(m).LineStyle = get(ALLlineHand(i),'LineStyle');
        Pline(m).color = get(ALLlineHand(i),'color');
        Pline(m).tag = tag;
        m = m + 1;      havePline = 1;
    end
end

% Patches may have associated particular meanings (eg Focal Mecas), but
% they will loose them here. Maybe in the future I'll do something better.
j = 1;
ALLpatchHand = findobj(get(handles.axes1,'Child'),'Type','patch');
for (i = 1:length(ALLpatchHand))
    tag = get(ALLpatchHand(i),'Tag');
    xx = get(ALLpatchHand(i),'XData');     yy = get(ALLpatchHand(i),'YData');
    Patches(j).x = xx(:);           Patches(j).y = yy(:);
    Patches(j).LineWidth = get(ALLpatchHand(i),'LineWidth');
    Patches(j).LineStyle = get(ALLpatchHand(i),'LineStyle');
    Patches(j).EdgeColor = get(ALLpatchHand(i),'EdgeColor');
    Patches(j).FaceColor = get(ALLpatchHand(i),'FaceColor');
    Patches(j).ud = get(ALLpatchHand(i),'UserData');
    Patches(j).tag = tag;
    j = j + 1;      havePatches = 1;
end

ALLtextHand = findobj(get(handles.axes1,'Child'),'Type','text');
for (i = 1:length(ALLtextHand))
    Texto(n).str = get(ALLtextHand(i),'String');
    if (isempty(Texto(n).str)),  continue;   end
    Texto(n).pos = get(ALLtextHand(i),'Position');       Texto(n).FontAngle = get(ALLtextHand(i),'FontAngle');
    Texto(n).angle = get(ALLtextHand(i),'Rotation');     Texto(n).Tag = get(ALLtextHand(i),'Tag');
    Texto(n).color = get(ALLtextHand(i),'color');        Texto(n).FontName = get(ALLtextHand(i),'FontName');
    Texto(n).FontSize = get(ALLtextHand(i),'FontSize');  Texto(n).FontWeight = get(ALLtextHand(i),'FontWeight');
    Texto(n).HorizontalAlignment = get(ALLtextHand(i),'HorizontalAlignment');
    Texto(n).VerticalAlignment = get(ALLtextHand(i),'VerticalAlignment');
    haveText = 1;   n = n + 1;
end
save(fname,'grd_name','img_pal', 'havePline','Pline', 'haveMBtrack', 'MBtrack','MBbar', ...
    'haveText','Texto', 'haveSymbol','Symbol', 'haveCircleGeo','CircleGeo', 'haveCircleCart', ...
    'havePlineAsPoints','PlineAsPoints','CircleCart', 'map_limits', 'havePatches', 'Patches', ...
    'haveCoasts', 'coastUD','havePolitic', 'politicUD','haveRivers', 'riversUD')
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function FileSaveGMT_script_Callback(hObject, eventdata, handles, opt)
if (handles.no_file == 1),      return;      end
write_gmt_script(handles, opt)

% --------------------------------------------------------------------
function DatasetsEarthquakes_Callback(hObject, eventdata, handles, opt)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
if (nargin == 3),   opt = [];   end         % When using the internal seismicity file
earthquakes(handles.figure1,handles.axes1,opt);

% --------------------------------------------------------------------
function GeophysicsFocMec_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
focal_meca(handles.figure1,handles.axes1);

% --------------------------------------------------------------------
function DatasetsHotspots_Callback(hObject, eventdata, handles)
% Read hotspot.dat which has 4 columns (lon lat name age)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
fid = fopen([handles.path_data 'hotspots.dat'],'r');
tline = fgetl(fid);             % Jump the header line
todos = fread(fid,'*char');     fclose(fid);
[hot.x hot.y hot.name hot.age] = strread(todos,'%f %f %s %f');     % Note: hot.name is a cell array of chars
clear todos;

% Get rid of Fogspots that are outside the map limits
[x,y,indx,indy] = aux_funs('in_map_region',handles,hot.x,hot.y,0,[]);
hot.name(indx) = [];   hot.age(indx) = [];
hot.name(indy) = [];   hot.age(indy) = [];
n_hot = length(x);    h_hotspot = zeros(1,n_hot);
for (i = 1:n_hot)
    h_hotspot(i) = line(x(i),y(i),'Marker','p','MarkerFaceColor','r',...
        'MarkerEdgeColor','k','MarkerSize',10,'Tag','hotspot','Userdata',i);
end
draw_funs(h_hotspot,'hotspot',hot)

% --------------------------------------------------------------------
function DatasetsVolcanoes_Callback(hObject, eventdata, handles)
% Read volcanoes.dat which has 6 columns (lat lon name ...)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
fid = fopen([handles.path_data 'volcanoes.dat'],'r');
todos = fread(fid,'*char');
[volc.y volc.x volc.name region volc.desc volc.dating] = strread(todos,'%f %f %s %s %s %s');
fclose(fid);    clear region todos

% Get rid of Volcanoes that are outside the map limits
[x,y,indx,indy] = aux_funs('in_map_region',handles,volc.x,volc.y,0,[]);
volc.name(indx) = [];       volc.desc(indx) = [];       volc.dating(indx) = [];
volc.name(indy) = [];       volc.desc(indy) = [];       volc.dating(indy) = [];
n_volc = length(x);    h_volc = zeros(1,n_volc);
for (i = 1:n_volc)
    h_volc(i) = line(x(i),y(i),'Marker','^','MarkerFaceColor','y',...
        'MarkerEdgeColor','k','MarkerSize',8,'Tag','volcano','Userdata',i);
end
draw_funs(h_volc,'volcano',volc)

% --------------------------------------------------------------------
function DatasetsTides_Callback(hObject, eventdata, handles)
load([handles.path_data 't_xtide.mat']);
% Get rid of Tide stations that are outside the map limits
[x,y] = aux_funs('in_map_region',handles,xharm.longitude,xharm.latitude,0,[]);
h_tides = line(x,y,'Marker','^','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',6,...
    'LineStyle','none','Tag','TideStation');
draw_funs(h_tides,'TideStation',[])

% --------------------------------------------------------------------
function DatasetsIsochrons_Callback(hObject, eventdata, handles, opt)
% Read multisegment isochrons.dat which has 3 columns (lat lon id)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
if (nargin == 4)        % Read a ascii multi-segment with info file
    str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
    [FileName,PathName] = put_or_get_file(handles,str1,'Select File','get');
    if (FileName == 0),     return;     end
    tag = 'Unnamed';        fname = [PathName FileName];
else
    tag = 'isochron';       fname = [handles.path_data 'isochrons.dat'];
end
xx = get(handles.axes1,'Xlim');           yy = get(handles.axes1,'Ylim');
set(handles.figure1,'pointer','watch')
[bin,n_column,multi_seg,n_headers] = guess_file(fname);
if (n_column == 1 && multi_seg == 0)        % Take it as a file names list
    fid = fopen(fname);
    c = char(fread(fid))';      fclose(fid);
    names = strread(c,'%s','delimiter','\n');   clear c fid;
else
    names = {fname};
end

if (nargin == 3),   ix = 2;     iy = 1;
else                ix = 1;     iy = 2;     end
tol = 0.5;

for (k=1:length(names))
    fname = names{k};
    j = strfind(fname,filesep);
    if (isempty(j)),    fname = [PathName fname];   end
    [numeric_data,multi_segs_str] = text_read(fname,NaN,NaN,'>');
	n_isoc = 0;     n_segments = length(numeric_data);
	h_isoc = ones(n_segments,1)*NaN;   % This is the maximum we can have
	n_clear = false(n_segments,1);
	for i=1:n_segments
        % Get rid of points that are outside the map limits
        [tmpx,tmpy] = aux_funs('in_map_region',handles,numeric_data{i}(:,ix),numeric_data{i}(:,iy),tol,[xx yy]);
        if (~isempty(tmpx))
            n_isoc = n_isoc + 1;
            h_isoc(i) = line(tmpx,tmpy,'Linewidth',handles.DefLineThick,'Color',handles.DefLineColor,'Tag',tag,'Userdata',n_isoc);
        else
            n_clear(i) = 1;             % Store indexes for clearing vanished segments info
        end
	end
	multi_segs_str(n_clear) = [];       % Clear the unused info
	
	ind = isnan(h_isoc);    h_isoc(ind) = [];      % Clear unused rows in h_isoc (due to over-dimensioning)
	draw_funs(h_isoc,'isochron',multi_segs_str)
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DatasetsPlateBound_PB_All_Callback(hObject, eventdata, handles)
% Read and plot the of the modified (by me) Peter Bird's Plate Boundaries
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
set(handles.figure1,'pointer','watch')
load([handles.path_data 'PB_boundaries.mat'])

% ------------------
% Get rid of boundary segments that are outside the map limits
xx = get(handles.axes1,'Xlim');      yy = get(handles.axes1,'Ylim');
tol = 0.5;
% ------------------ OTF class
n = length(OTF);    k = [];
for i = 1:n
    ind = find(OTF(i).x_otf < xx(1)-tol | OTF(i).x_otf > xx(2)+tol);
    OTF(i).x_otf(ind) = [];     OTF(i).y_otf(ind) = [];
    if isempty(OTF(i).x_otf),   k = [k i];  end         % k is a counter to erase out-of-map segments
end;    OTF(k) = [];
n = length(OTF);    k = [];
for i = 1:n
    ind = find(OTF(i).y_otf < yy(1)-tol | OTF(i).y_otf > yy(2)+tol);
    OTF(i).x_otf(ind) = [];     OTF(i).y_otf(ind) = [];
    if isempty(OTF(i).x_otf),   k = [k i];  end
end;    OTF(k) = [];
% ------------------ OSR class
n = length(OSR);    k = [];
for i = 1:n
    ind = find(OSR(i).x_osr < xx(1)-tol | OSR(i).x_osr > xx(2)+tol);
    OSR(i).x_osr(ind) = [];     OSR(i).y_osr(ind) = [];
    if isempty(OSR(i).x_osr),   k = [k i];  end
end;    OSR(k) = [];
n = length(OSR);    k = [];
for i = 1:n
    ind = find(OSR(i).y_osr < yy(1)-tol | OSR(i).y_osr > yy(2)+tol);
    OSR(i).x_osr(ind) = [];     OSR(i).y_osr(ind) = [];
    if isempty(OSR(i).x_osr),   k = [k i];  end
end;    OSR(k) = [];
% ------------------ CRB class
n = length(CRB);    k = [];
for i = 1:n
    ind = find(CRB(i).x_crb < xx(1)-tol | CRB(i).x_crb > xx(2)+tol);
    CRB(i).x_crb(ind) = [];     CRB(i).y_crb(ind) = [];
    if isempty(CRB(i).x_crb),   k = [k i];  end
end;    CRB(k) = [];
n = length(CRB);    k = [];
for i = 1:n
    ind = find(CRB(i).y_crb < yy(1)-tol | CRB(i).y_crb > yy(2)+tol);
    CRB(i).x_crb(ind) = [];     CRB(i).y_crb(ind) = [];
    if isempty(CRB(i).x_crb),   k = [k i];  end
end;    CRB(k) = [];
% ------------------ CTF class
n = length(CTF);    k = [];
for i = 1:n
    ind = find(CTF(i).x_ctf < xx(1)-tol | CTF(i).x_ctf > xx(2)+tol);
    CTF(i).x_ctf(ind) = [];     CTF(i).y_ctf(ind) = [];
    if isempty(CTF(i).x_ctf),   k = [k i];  end
end;    CTF(k) = [];
n = length(CTF);    k = [];
for i = 1:n
    ind = find(CTF(i).y_ctf < yy(1)-tol | CTF(i).y_ctf > yy(2)+tol);
    CTF(i).x_ctf(ind) = [];     CTF(i).y_ctf(ind) = [];
    if isempty(CTF(i).x_ctf),   k = [k i];  end
end;    CTF(k) = [];
% ------------------ CCB class
n = length(CCB);    k = [];
for i = 1:n
    ind = find(CCB(i).x_ccb < xx(1)-tol | CCB(i).x_ccb > xx(2)+tol);
    CCB(i).x_ccb(ind) = [];     CCB(i).y_ccb(ind) = [];
    if isempty(CCB(i).x_ccb),   k = [k i];  end
end;    CCB(k) = [];
n = length(CCB);    k = [];
for i = 1:n
    ind = find(CCB(i).y_ccb < yy(1)-tol | CCB(i).y_ccb > yy(2)+tol);
    CCB(i).x_ccb(ind) = [];     CCB(i).y_ccb(ind) = [];
    if isempty(CCB(i).x_ccb),   k = [k i];  end
end;    CCB(k) = [];
% ------------------ OCB class
n = length(OCB);    k = [];
for i = 1:n
    ind = find(OCB(i).x_ocb < xx(1)-tol | OCB(i).x_ocb > xx(2)+tol);
    OCB(i).x_ocb(ind) = [];     OCB(i).y_ocb(ind) = [];
    if isempty(OCB(i).x_ocb),   k = [k i];  end
end;    OCB(k) = [];
n = length(OCB);    k = [];
for i = 1:n
    ind = find(OCB(i).y_ocb < yy(1)-tol | OCB(i).y_ocb > yy(2)+tol);
    OCB(i).x_ocb(ind) = [];     OCB(i).y_ocb(ind) = [];
    if isempty(OCB(i).x_ocb),   k = [k i];  end
end;    OCB(k) = [];
% ------------------ SUB class
n = length(SUB);    k = [];
for i = 1:n
    ind = find(SUB(i).x_sub < xx(1)-tol | SUB(i).x_sub > xx(2)+tol);
    SUB(i).x_sub(ind) = [];     SUB(i).y_sub(ind) = [];
    if isempty(SUB(i).x_sub),   k = [k i];  end
end;    SUB(k) = [];
n = length(SUB);    k = [];
for i = 1:n
    ind = find(SUB(i).y_sub < yy(1)-tol | SUB(i).y_sub > yy(2)+tol);
    SUB(i).x_sub(ind) = [];     SUB(i).y_sub(ind) = [];
    if isempty(SUB(i).x_sub),   k = [k i];  end
end;    SUB(k) = [];

% ------------------ Finally do the ploting ------------------------------------
% Plot the OSR class
n = length(OSR);    h_PB_All_OSR = zeros(n,1);
for i = 1:n
    line(OSR(i).x_osr,OSR(i).y_osr,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_OSR(i) = line(OSR(i).x_osr,OSR(i).y_osr,'Linewidth',2,'Color','r','Tag','PB_All','Userdata',i);
end
% Plot the OTF class
n = length(OTF);    h_PB_All_OTF = zeros(n,1);
for i = 1:n
    line(OTF(i).x_otf,OTF(i).y_otf,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_OTF(i) = line(OTF(i).x_otf,OTF(i).y_otf,'Linewidth',2,'Color','g','Tag','PB_All','Userdata',i);
end
% Plot the CRB class
n = length(CRB);    h_PB_All_CRB = zeros(n,1);
for i = 1:n
    line(CRB(i).x_crb,CRB(i).y_crb,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_CRB(i) = line(CRB(i).x_crb,CRB(i).y_crb,'Linewidth',2,'Color','b','Tag','PB_All','Userdata',i);
end
% Plot the CTF class
n = length(CTF);    h_PB_All_CTF = zeros(n,1);
for i = 1:n
    line(CTF(i).x_ctf,CTF(i).y_ctf,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_CTF(i) = line(CTF(i).x_ctf,CTF(i).y_ctf,'Linewidth',2,'Color','y','Tag','PB_All','Userdata',i);
end
% Plot the CCB class
n = length(CCB);    h_PB_All_CCB = zeros(n,1);
for i = 1:n
    line(CCB(i).x_ccb,CCB(i).y_ccb,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_CCB(i) = line(CCB(i).x_ccb,CCB(i).y_ccb,'Linewidth',2,'Color','m','Tag','PB_All','Userdata',i);
end
% Plot the OCB class
n = length(OCB);    h_PB_All_OCB = zeros(n,1);
for i = 1:n
    line(OCB(i).x_ocb,OCB(i).y_ocb,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_OCB(i) = line(OCB(i).x_ocb,OCB(i).y_ocb,'Linewidth',2,'Color','c','Tag','PB_All','Userdata',i);
end
% Plot the SUB class
n = length(SUB);    h_PB_All_SUB = zeros(n,1);
for i = 1:n
    line(SUB(i).x_sub,SUB(i).y_sub,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_SUB(i) = line(SUB(i).x_sub,SUB(i).y_sub,'Linewidth',2,'Color','c','Tag','PB_All','Userdata',i);
end

% Join all line handles into a single variable
h.OSR = h_PB_All_OSR;    h.OTF = h_PB_All_OTF;    h.CRB = h_PB_All_CRB;    h.CTF = h_PB_All_CTF;
h.CCB = h_PB_All_CCB;    h.OCB = h_PB_All_OCB;    h.SUB = h_PB_All_SUB;
% Join all data into a single variable
data.OSR = OSR;    data.OTF = OTF;    data.CRB = CRB;    data.CTF = CTF;
data.CCB = CCB;    data.OCB = OCB;    data.SUB = SUB;
draw_funs(h,'PlateBound_All_PB',data);      set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DatasetsPlateBoundAfrica_Callback(hObject, eventdata, handles, opt)
% This function is currently not used
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
switch opt
    case 'Africa',  file = 'PB2002_AF.dat';
    case 'Eurasia', file = 'PB2002_EU.dat';
    case 'NorthAmerica',  file = 'PB2002_NA.dat';
end
fid = fopen([handles.path_data file],'r');
todos = fread(fid,'*char');     [x y] = strread(todos,'%f %f');
fclose(fid);    clear todos
line(x,y,'Linewidth',3,'Tag','Africa','Color','k');

% --------------------------------------------------------------------
function DatasetsODP_DSDP_Callback(hObject, eventdata, handles,opt)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
set(handles.figure1,'pointer','watch')
fid = fopen([handles.path_data 'DSDP_ODP.dat'],'r');
todos = fread(fid,'*char');
[ODP.x ODP.y zz ODP.leg ODP.site ODP.z ODP.penetration] = strread(todos,'%f %f %s %s %s %s %s');
fclose(fid);    clear todos zz

% Get rid of Sites that are outside the map limits
[ODP.x,ODP.y,indx,indy] = aux_funs('in_map_region',handles,ODP.x,ODP.y,0,[]);

ODP.leg(indx) = [];     ODP.site(indx) = [];    ODP.z(indx) = [];   ODP.penetration(indx) = [];
ODP.leg(indy) = [];     ODP.site(indy) = [];    ODP.z(indy) = [];   ODP.penetration(indy) = [];

% If there no sites left, return
if isempty(ODP.x)
    set(handles.figure1,'pointer','arrow');    msgbox('Warning: There are no sites inside this area.','Warning');    return;
end

% Find where in file is the separation of DSDP from ODP legs
ind = find(str2double(ODP.leg) >= 100);
if ~isempty(ind),   ind = ind(1);   end
if (strcmp(opt,'ODP'))      % If only ODP sites were asked remove DSDP from data structure
    ODP.x(1:ind-1) = [];    ODP.y(1:ind-1) = [];    ODP.z(1:ind-1) = [];
    ODP.leg(1:ind-1) = [];  ODP.site(1:ind-1) = []; ODP.penetration(1:ind-1) = [];
elseif (strcmp(opt,'DSDP'))
    ODP.x(ind:end) = [];    ODP.y(ind:end) = [];    ODP.z(ind:end) = [];
    ODP.leg(ind:end) = [];  ODP.site(ind:end) = []; ODP.penetration(ind:end) = [];
end

n_sites = length(ODP.x);    h_sites = zeros(n_sites,1);
if (strcmp(opt,'DSDP'))
    if (n_sites == 0)           % If there are no sites, give a warning and exit
        set(handles.figure1,'pointer','arrow');        msgbox('Warning: There are no DSDP sites inside this area.','Warning');    return;
    end
    for i = 1:n_sites
        h_sites(i) = line(ODP.x(i),ODP.y(i),'Marker','o','MarkerFaceColor','g',...
            'MarkerEdgeColor','k','MarkerSize',8,'Tag','DSDP','Userdata',i);
    end
    draw_funs(h_sites,'ODP',ODP)
elseif (strcmp(opt,'ODP'))
    if (n_sites == 0)           % If there are no sites, give a warning and exit
        set(handles.figure1,'pointer','arrow');        msgbox('Warning: There are no ODP sites inside this area.','Warning');    return;
    end
    for i = 1:n_sites
        h_sites(i) = line(ODP.x(i),ODP.y(i),'Marker','o','MarkerFaceColor','r',...
            'MarkerEdgeColor','k','MarkerSize',8,'Tag','ODP','Userdata',i);
    end
    draw_funs(h_sites,'ODP',ODP)
else
    h_sites = zeros(length(1:ind-1),1);
    for i = 1:ind-1
        h_sites(i) = line(ODP.x(i),ODP.y(i),'Marker','o','MarkerFaceColor','g',...
            'MarkerEdgeColor','k','MarkerSize',8,'Tag','DSDP','Userdata',i);
    end
    draw_funs(h_sites,'ODP',ODP)
    h_sites = zeros(length(ind:n_sites),1);
    for (i = 1:length(ind:n_sites))
        j = i + ind - 1;
        h_sites(i) = line(ODP.x(j),ODP.y(j),'Marker','o','MarkerFaceColor','r',...
            'MarkerEdgeColor','k','MarkerSize',8,'Tag','ODP','Userdata',j);
    end
    draw_funs(h_sites,'ODP',ODP)
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DatasetsCoastLineNetCDF_Callback(hObject, eventdata, handles, res)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
if (nargin == 3),   res = 'l';  end
set(handles.figure1,'pointer','watch')

lon = get(handles.axes1,'Xlim');      lat = get(handles.axes1,'Ylim');
lt = handles.DefLineThick;  lc = handles.DefLineColor;
opt_R = ['-R' num2str(lon(1)) '/' num2str(lon(2)) '/' num2str(lat(1)) '/' num2str(lat(2))];

switch res
    case 'c',        opt_res = '-Dc';        pad = 1.0;
    case 'l',        opt_res = '-Dl';        pad = 0.5;
    case 'i',        opt_res = '-Di';        pad = 0.1;
    case 'h',        opt_res = '-Dh';        pad = 0.02;
    case 'f',        opt_res = '-Df';        pad = 0.005;
end
coast = shoredump(opt_R,opt_res,'-A1/1/1');

% Get rid of data that are outside the map limits
%lon(1) = max(lon(1)-pad, -179.99);      lon(2) = min(lon(2)+pad, 179.99);
lon = lon - [pad -pad];     lat = lat - [pad -pad];
indx = (coast(1,:) < lon(1) | coast(1,:) > lon(2));
coast(:,indx) = [];
indx = (coast(2,:) < lat(1) | coast(2,:) > lat(2));
coast(:,indx) = [];

if (~all(isnan(coast(:))))
	h_coast = line(coast(1,:),coast(2,:),'Linewidth',lt,'Color',lc,'Tag','CoastLineNetCDF','UserData',opt_res(3));
	draw_funs(h_coast,'Coastline_uicontext')    % Set line's uicontextmenu
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DatasetsPoliticalBoundaries_Callback(hObject, eventdata, handles, type, res)
% TYPE is: '1' -> National Boundaries
%          '2' -> State Boundaries
%          '3' -> Marine Boundaries
%          'a' -> All Boundaries
% RES is:  'c' or 'l' or 'i' or 'h' or 'f' (gmt database resolution)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file

set(handles.figure1,'pointer','watch')
lon = get(handles.axes1,'Xlim');      lat = get(handles.axes1,'Ylim');
opt_R = ['-R' num2str(lon(1)) '/' num2str(lon(2)) '/' num2str(lat(1)) '/' num2str(lat(2))];

switch type
    case '1',        opt_N = '-N1';
    case '2',        opt_N = '-N2';
    case '3',        opt_N = '-N3';
    case 'a',        opt_N = '-Na';
end

switch res
    case 'c',        opt_res = '-Dc';        pad = 1;
    case 'l',        opt_res = '-Dl';        pad = 0.5;
    case 'i',        opt_res = '-Di';        pad = 0.1;
    case 'h',        opt_res = '-Dh';        pad = 0.05;
    case 'f',        opt_res = '-Df';        pad = 0.01;
end
boundaries = shoredump(opt_R,opt_N,opt_res);

% Get rid of data that are outside the map limits
lon = lon - [pad -pad];     lat = lat - [pad -pad];
indx = (boundaries(1,:) < lon(1) | boundaries(1,:) > lon(2));
boundaries(:,indx) = [];
indx = (boundaries(2,:) < lat(1) | boundaries(2,:) > lat(2));
boundaries(:,indx) = [];

if (~all(isnan(boundaries(:))))
	h_boundaries = line(boundaries(1,:),boundaries(2,:),'Linewidth',handles.DefLineThick,...
        'Color',handles.DefLineColor,'Tag','PoliticalBoundaries', 'UserData',[opt_res(3) opt_N(3)]);
	draw_funs(h_boundaries,'Coastline_uicontext')    % Set line's uicontextmenu
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DatasetsRivers_Callback(hObject, eventdata, handles, type, res)
% TYPE is: '1' -> Permanent major rivers;           '2' -> Additional major rivers
%          '3' -> Additional rivers                 '4' -> Minor rivers
%          '5' -> Intermittent rivers - major       '6' -> Intermittent rivers - additional
%          '7' -> Intermittent rivers - minor       '8' -> Major canals
%          '9' -> Minor canals
%          'a' -> All rivers and canals (1-10)      'r' -> All permanent rivers (1-4)
%          'i' -> All intermittent rivers (5-7)
% RES is:  'c' or 'l' or 'i' or 'h' or 'f' (gmt database resolution)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file

set(handles.figure1,'pointer','watch')
lon = get(handles.axes1,'Xlim');      lat = get(handles.axes1,'Ylim');
opt_R = ['-R' num2str(lon(1)) '/' num2str(lon(2)) '/' num2str(lat(1)) '/' num2str(lat(2))];

switch type
    case '1',        opt_I = '-I1';         case '2',        opt_I = '-I2';
    case '3',        opt_I = '-I3';
    case '5',        opt_I = '-I5';         case '6',        opt_I = '-I6';
    case '7',        opt_I = '-I7';
    case 'a',        opt_I = '-Ia';
    case 'r',        opt_I = '-Ir';         case 'i',        opt_I = '-Ii';
end

switch res
    case 'c',        opt_res = '-Dc';        pad = 1;
    case 'l',        opt_res = '-Dl';        pad = 0.5;
    case 'i',        opt_res = '-Di';        pad = 0.1;
    case 'h',        opt_res = '-Dh';        pad = 0.05;
    case 'f',        opt_res = '-Df';        pad = 0.01;
end
rivers = shoredump(opt_R,opt_I,opt_res);

% Get rid of data that are outside the map limits
lon = lon - [pad -pad];     lat = lat - [pad -pad];
indx = (rivers(1,:) < lon(1) | rivers(1,:) > lon(2));
rivers(:,indx) = [];
indx = (rivers(2,:) < lat(1) | rivers(2,:) > lat(2));
rivers(:,indx) = [];

if (~all(isnan(rivers(:))))
	h_rivers = line(rivers(1,:),rivers(2,:),'Linewidth',handles.DefLineThick,'Color',handles.DefLineColor,...
        'Tag','Rivers', 'UserData',[opt_res(3) opt_I(3:end)]);
	draw_funs(h_rivers,'Coastline_uicontext')    % Set line's uicontextmenu
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DatasetsAtlas_Callback(hObject, eventdata, handles)
% Plot countries as patches with pre-set colors. If the axes was not yet
% created, the atlas GUI will take care of that as well.
if (handles.geog && ~handles.no_file)   % If a file is loaded and is not geog
    errordlg('This operation is currently possible only for geographic type data','Error');     return;
end
atlas(handles.figure1,handles.axes1);

% --------------------------------------------------------------------
function DatasetsCities_Callback(hObject, eventdata, handles,opt)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
if strcmp(opt,'major')
    fid = fopen([handles.path_data 'wcity_major.dat'],'r');
    tag = 'City_major';
elseif strcmp(opt,'other')
    fid = fopen([handles.path_data 'wcity.dat'],'r');
    tag = 'City_other';
end
todos = fread(fid,'*char');     fclose(fid);
[city.x city.y city.name] = strread(todos,'%f %f %s');     % Note: city.name is a cell array of chars
% Get rid of Cities that are outside the map limits
[x,y,indx,indy] = aux_funs('in_map_region',handles,city.x,city.y,0,[]);
city.name(indx) = [];       city.name(indy) = [];
n_city = length(x);

if (n_city == 0),   return;     end     % No cities inside area. Return.
h_city = line(x,y,'LineStyle','none','Marker','o','MarkerFaceColor','k',...
    'MarkerEdgeColor','w','MarkerSize',6,'Tag',tag);
draw_funs(h_city,'DrawSymbol')                  % Set symbol's uicontextmenu

% Estimate the text position shift in order that it doesn't fall over the city symbol 
pos = get(handles.figure1,'Position');
x_lim = get(handles.axes1,'xlim');
z1 = 7 / pos(3);
dx = z1 * (x_lim(2) - x_lim(1));

city.name = strrep(city.name,'_',' ');          % Replace '_' by ' '
textHand = zeros(1,n_city);
for i = 1:n_city                                % Plot the City names
    textHand(i) = text(x(i)+dx,y(i),0,city.name{i},'Tag',tag);
    draw_funs(textHand(i),'DrawText')           % Set text's uicontextmenu
end

% --------------------------------------------------------------------
function ImageMapLimits_Callback(hObject, eventdata, handles)
% Change the Image limits by asking it's corner coordinates
region = bg_region('empty');     % region contains [x_min x_max y_min y_max is_geog]
if isempty(region),     return;     end     % User gave up
img = get(findobj(handles.axes1,'Type','image'),'CData');
X = region(1:2);    Y = region(3:4);        handles.geog = guessGeog(region(1:4));
x_inc = diff(X) / size(img,2);              y_inc = diff(Y) / size(img,1);
dx2 = x_inc / 2;                            dy2 = y_inc / 2;
X = X + [dx2 -dx2];                         Y = Y + [dy2 -dy2];     % Make it such that the pix-reg info = region
handles.head = [X Y 0 255 0 x_inc y_inc];   handles.image_type = 3;

% Flipud the image if necessary
if (strcmp(get(handles.axes1,'YDir'),'reverse')),    img = flipdim(img,1);      end
show_image(handles,'New Limits',X,Y,img,0,'xy',0,1);

% --------------------------------------------------------------------
function GeophysicsTTT_Callback(hObject, eventdata, handles,opt)
if (aux_funs('msg_dlg',14,handles));     return;      end
if (nargin == 3),   opt = [];   end
if isempty(opt)     % Plot the source
    pt = ginput_pointer(1,'crosshair');
    symbHand = line(pt(1,1),pt(1,2),'Marker','o','MarkerFaceColor','y',...
            'MarkerEdgeColor','k','MarkerSize',10,'Tag','TTT');
    draw_funs(symbHand,'DrawSymbol')            % Set symbol's uicontextmenu
elseif (strcmp(opt,'load'))
    str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)'};
    [FileName,PathName] = put_or_get_file(handles,str1,'Select input xy_time file name','get');
    if isequal(FileName,0);     return;     end
    out = draw_funs([PathName FileName],'ImportLine');
    if (size(out,2) ~= 3)
        errordlg('Wrong choice. For using this option the file MUST have 3 columns (position and time).','Error'); return
    end
    symbHand = line(out(:,1),out(:,2),'Marker','o','MarkerFaceColor','y','linestyle','none',...
        'MarkerEdgeColor','k','MarkerSize',10,'Tag','TTT','UserData',out(:,3));
    setappdata(symbHand,'TTTimes',out)
    draw_funs(symbHand,'DrawSymbol')            % Set symbol's uicontextmenu    
else                % Compute
    h_src = findobj('Type','line','Tag','TTT');
    if (isempty(h_src)),    errordlg('Yes I compute, but WHAT?','Error');   return;     end
    if (length(h_src) > 1)
        errordlg('More than one source found. This is not allowed neither in single source or Ray tracing modes. Be modest.','Error');     return
    end
    xy_t = getappdata(h_src,'TTTimes');             % See if we have anything in appdata
    if (isempty(xy_t)),     single_src = 1;         % Single source mode
    else                    single_src = 0; end     % Ray tracing mode
    [X,Y,Z,head] = load_grd(handles);
    if isempty(Z),   return;     end;    % An error message was already issued
    if (handles.have_nans)
        errordlg('Bathymetry grid cannot have NaNs','Error');   return;
    end
    h_info = [head(1) head(2) head(3) head(4) head(8) head(9)];
    if (single_src)         % Single source mode
        xx = get(h_src,'XData');    yy = get(h_src,'YData');
        set(handles.figure1,'pointer','watch');
        tt = wave_travel_time(double(Z)+0,h_info,[xx yy],handles.geog);     % Adding 0 is important (pointers, mex, etc)
        tit = 'Tsunami Travel Times';
    else                    % Find ray tracing solution
        h = waitbar(0,'Hold on your camels: computing solution');
        xx = get(h_src,'XData');    yy = get(h_src,'YData');    tempo = get(h_src,'UserData');
        if (length(xx) ~= length(tempo))            % Some(s) station(s) has been killed
            [c,ia,ib] = setxor(xx(:),xy_t(:,1));    % Find which
            tempo(ib) = [];                         % Remove the corresponding times
        end
        tmp = 0;
        for (i = 1:length(xx))
            tt = wave_travel_time(double(Z)+0,h_info,[xx(i) yy(i)],handles.geog);     % 0 because of pointers, mex, etc
            tmp = tmp + (tt - tempo(i)) .^2;
            waitbar(i/length(xx))
        end
        tt = sqrt(tmp/size(xy_t,1));        tit = 'Ray tracing solution';        close(h)
    end
    head(5) = min(min(tt));    head(6) = max(max(tt));
    GRD_save_or_display(handles,X,Y,tt,head,tit,tit)
end

% --------------------------------------------------------------------
function GeophysicsSwanCompute_Callback(hObject, eventdata, handles)
%#function waitbar

small = 1e-6;   % Used in postion comparation. It places the accuracy at sub-meter level
Z_bat = [];     head_bat = [];      Z_src = [];     head_src = [];
I_am_in_bat = 0;        % Used for trying to guess if the call was made from the bthymetry
I_am_in_source = 0;     % figure, or the tsunami source figure.
% The following are flags to signal how much information is already known before
% calling "swan_options". They will be set to 1 by the guessing code. 
% I don't try to guess a "deform_only".
nothing_in_memory = 0;      bat_and_deform_with_maregs = 0;
bat_and_deform = 0;         bat_with_maregs = 0;            bat_only = 0;

% Well, for writing a movie the user would have a chance to select a directory, but thats too complicated to advance
if (strcmp(handles.last_dir,handles.home_dir))
    errordlg('Writting on Mirone''s installation directory is not allowed. Change the working (see File -> preferences)','Error')
    return
end

if (handles.no_file == 1 || handles.image_type ~= 1 || length(findobj('Type','figure')) == 1)
    if (handles.image_type == 1)
        % This is a risk because nothing guaranties that the grid contains bathymetry data
        [X,Y,Z_bat,head_bat] = load_grd(handles);   Z_bat = double(Z_bat);
        if isempty(Z_bat),  return;     end;    % An error message was already issued
        bat_only = 1;
    end
else
    if (isempty(handles.calling_figure))    % I am in the bathymetry figure
        I_am_in_bat = 1;    I_am_in_source = 0;
        h_src = findobj('Type','figure','Name','Okada deformation');
        if ishandle(h_src)
            Z_src = getappdata(h_src,'dem_z');   Z_src = double(Z_src);
            head_src = getappdata(h_src,'GMThead');
        else
            I_am_in_source = 0;
            %errordlg('You need a valid deformation Figure to run this option.','Error');    return;
        end
    else        % I (hope) am in the tsunami source (deformation) figure
        I_am_in_bat = 0;    I_am_in_source = 1;
        Z_bat  = getappdata(handles.calling_figure,'dem_z');   Z_bat = double(Z_bat);
        head_bat  = getappdata(handles.calling_figure,'GMThead');
    end

    if (I_am_in_bat)
        [X,Y,Z_bat,head_bat] = load_grd(handles);       Z_bat = double(Z_bat);
        if isempty(Z_bat),    return;     end;    % An error message was already issued
        % Now, if they where no errors, we can compare Z_bat and Z_src to see if they are compatible
        if (numel(Z_bat) ~= numel(Z_src))   % Shit
            if (isempty(Z_src))     % Why?
                errordlg('There was an error in fishing the deformation data.','Error');    return;
            end
            errordlg('Bathymetry and deformation grids have not the same size.','Error');   return;
        end
    elseif (I_am_in_source)
        [X,Y,Z_src,head_src] = load_grd(handles);       Z_src = double(Z_src);
        if isempty(Z_src),   return;     end;    % An error message was already issued
        if (numel(Z_bat) ~= numel(Z_src))   % Shit
            if (isempty(Z_bat))     % Why?
                errordlg('There was an error in fishing the bathymetry data.','Error');     return;
            end
            errordlg('Bathymetry and deformation grids have not the same size.','Error');   return;
        end
    end

    % If everything went well, we know that grids are of the same size but we still
    % have to test that the region covered is the same.
    if (I_am_in_bat &&  I_am_in_source)       % If bat & source figures test that they are compatible.
        if ( abs(head_src(1) - head_bat(1)) > small || abs(head_src(2) - head_bat(2)) > small || ...
               abs(head_src(3) - head_bat(3)) > small || abs(head_src(4) - head_bat(4)) > small )
            errordlg('Bathymetry & Source grids do not cover the same region','Error');     return
        end
        bat_and_deform = 1;
    end
end

% See if we have maregraphs (if they exist, that is interpreted as a computation request)
h_mareg = findobj('Type','line','Tag','Maregraph');
mareg_pos = [];
if (~isempty(h_mareg))
    side = zeros(1,length(h_mareg));
    for (i=1:length(h_mareg)),      side(i) = isappdata(h_mareg(i),'Side');    end
    if (any(side) & ~all(side))
        errordlg('ERROR: you cannot mix individual stations with stations on grid borders.','ERROR');   return
    end
    if (all(side))      % We are in the stations on grid borders case
        side = cell(1,length(h_mareg));
        for (i=1:length(h_mareg)),      side{i} = getappdata(h_mareg(i),'Side');    end
        side = strrep(side,'W','1');    side = strrep(side,'S','2');
        side = strrep(side,'E','3');    side = strrep(side,'N','4');
        [i,j] = sort(side);             % Now we can sort them to find the correct WSEN order
        h_mareg = h_mareg(j);           % And reorder the handles to follow the expected order
    end
    for i=1:length(h_mareg)
        mareg_pos = [mareg_pos; [get(h_mareg(i),'XData')' get(h_mareg(i),'YData')']];
    end
    if (bat_and_deform)
        bat_and_deform_with_maregs = 1;
    else                % OK, this also implies that we have a possibly valid bat file
        bat_with_maregs = 1;
    end
end

if (bat_and_deform_with_maregs)
    out = swan_options('bat_and_deform_with_maregs',mareg_pos,handles.figure1);
elseif (bat_and_deform && ~bat_and_deform_with_maregs)
    out = swan_options('bat_and_deform',handles.figure1);
elseif (bat_with_maregs)
    out = swan_options('bat_with_maregs',Z_bat,head_bat,mareg_pos,handles.figure1);
elseif (bat_only)
    out = swan_options('bat_only',Z_bat,head_bat,handles.figure1);
else
    out = swan_options;
    nothing_in_memory = 1;
end
pause(0.05);        % Give time to swan_options window to die

if (isempty(out)),  return;     end

if (nothing_in_memory)
    Z_bat = out.grid_Z_bat;     head_bat = out.grid_head_bat;
    Z_src = out.grid_Z_src;     head_src = out.grid_head_src;
end
if (bat_with_maregs || bat_only)
    Z_src = out.grid_Z_src;     head_src = out.grid_head_src;
end

extra_args1 = ' ';   extra_args2 = ' ';   extra_args3 = ' ';
extra_args4 = ' ';   extra_args5 = ' ';
if (isfield(out,'maregraph_xy') && isfield(out,'maregraph_data_name'))
    extra_args1 = ['-O' out.maregraph_data_name];
end
if (isfield(out,'opt_M')),    extra_args2 = out.opt_M;   end
if (isfield(out,'opt_N')),    extra_args3 = out.opt_N;   end
if (isfield(out,'opt_m')),    extra_args4 = out.opt_m;   end
if (isfield(out,'opt_G')),    extra_args5 = out.opt_G;   end

if (isempty(Z_bat) || isempty(head_bat) || isempty(Z_src) || isempty(head_src))
    errordlg('ERROR: one or more of the bat/source variables are empty where they souldn''t be.','Error');
    return;
end

Z_bat = double(Z_bat);      Z_src = double(Z_src);      % make sure they are both doubles

if (isfield(out,'maregraph_xy'))    % Ask for computation of maregraphs
    if (isfield(out,'opt_m'))       % Movie option
        movie = swan(Z_bat, head_bat, Z_src, head_src, out.params, out.maregraph_xy, extra_args1, ...
            extra_args2,extra_args3,extra_args5,'-f');
    else
        swan(Z_bat, head_bat, Z_src, head_src, out.params, out.maregraph_xy, extra_args1, ...
            extra_args2,extra_args3,extra_args5);
    end    
else                                % Compute grids or movie
    if (isfield(out,'opt_m'))       % Movie option
        movie = swan(Z_bat, head_bat, Z_src, head_src, out.params, extra_args2,extra_args3,extra_args5,'-f');
    else
        swan(Z_bat, head_bat, Z_src, head_src, out.params, extra_args2,extra_args3,extra_args5);
    end
end
pause(0.01);        % Give time to waitbar window to die
if (isfield(out,'opt_m')),   do_movie(handles,movie,'swan');   end

% --------------------------------------------------------------------
function GeophysicsSwanPlotStations_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',14,handles));     return;      end
zoom_state(handles,'maybe_off');
pt = ginput_pointer(1,'crosshair');
if (isempty(pt)),   return;  end
handles.maregraphs_count = handles.maregraphs_count + 1;    % Count number of maregraphs
%tag = ['Maregraph' num2str(handles.maregraphs_count)];
symbHand = line(pt(1,1),pt(1,2),'Marker','o','MarkerFaceColor','y',...
        'MarkerEdgeColor','k','MarkerSize',10,'Tag','Maregraph');
draw_funs(symbHand,'DrawSymbol')          % Set symbol's uicontextmenu
zoom_state(handles,'maybe_on');
guidata(hObject, handles)

% --------------------------------------------------------------------
function GeophysicsSwanGridBorderStations_Callback(hObject, eventdata, handles)
% Get the limits of smaller grid and plot them as a rectangle (individual lines) on current fig
if (aux_funs('msg_dlg',14,handles));     return;      end
small = 1e-5;       % Used in relative origin comparation.
str_R = [];         str_I =[];
adjust_w = 0;       adjust_s = 0;       % Flags used when the finer grid has to be adjusted
adjust_x_inc = 0;   adjust_y_inc = 0;   % in order to fit correctly with coarser grid
in_w_to_e = 0;      in_s_to_n = 0;

[X,Y,Z,head] = load_grd(handles);
if isempty(Z),   return;     end;    % An error message was already issued
str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select GMT grid','get');
if isequal(FileName,0);     return;     end

if (handles.grdformat > 1), FileName = [FileName '=' num2str(handles.grdformat)];     end
D = grdinfo_m([PathName FileName],'silent');

% Verify if, on the overlapping zone, the nodes of the larger grid cuincide with nodes of the smaler
xoff_w = abs(head(1) - D(1));   xoff_e = abs(head(2) - D(2));
yoff_s = abs(head(3) - D(3));   yoff_n = abs(head(4) - D(4));
dx_w = xoff_w / D(8);           dx_e = xoff_e / D(8);
dy_s = yoff_s / D(9);           dy_n = yoff_n / D(9);
if ((dx_w - fix(dx_w)) > small)     % Need to adjust at the west border
    adjust_w = 1;
    if ((dx_w - fix(dx_w)) > head(8)/2) % The adjustment is to the east
        in_w_to_e = 1;
    end
end
if ((dy_s - fix(dy_s)) > small)     % Need to adjust at the south border
    adjust_s = 1;
    if ((dy_s - fix(dy_s)) > head(9)/2) % The adjustment is to the north
        in_s_to_n = 1;
    end
end
% We also have to test if the grid increments of both grids are multiples
inc_x_ratio = head(8) / D(8);    inc_y_ratio = head(9) / D(9);
if ((inc_x_ratio - fix(inc_x_ratio)) > 1e-3)
    adjust_x_inc = 1;
    new_x_inc = (D(2) - D(1)) / (fix(inc_x_ratio) - 1);
end
if ((inc_y_ratio - fix(inc_y_ratio)) > 1e-3)
    adjust_y_inc = 1;
    new_y_inc = (D(4) - D(3)) / (fix(inc_y_ratio) - 1);
end

if (adjust_w || adjust_s || adjust_x_inc || adjust_y_inc)
    warndlg('The fine grid doesn''t fit well within the larger grid. Trying to fix it.','Warning')
    if (adjust_w && adjust_s) % We have to adjust all borders
        if (in_w_to_e)      % Move x origin to the east
            new_w = fix(dx_w) * (head(8) + 1);
            new_e = new_w + D(2) - D(1);
        else                % Move x origin to the west
            new_w = fix(dx_w) * head(8);
            new_e = new_w + D(2) - D(1);
        end
        if (in_s_to_n)      % Move y origin to the north
            new_s = fix(dy_s) * (head(9) + 1);
            new_n = new_s + D(4) - D(3);
        else                % Move y origin to the south
            new_s = fix(dy_s) * head(9);
            new_n = new_s + D(4) - D(3);
        end
        str_R = [' -R' num2str(new_w,10) '/' num2str(new_s,10) '/' num2str(new_e,10) '/' num2str(new_n,10)];
    elseif (adjust_w)        % Need to adjust only est and west borders
        if (in_w_to_e)      % Move x origin to the east
            new_w = fix(dx_w) * (head(8) + 1);
            new_e = new_w + D(2) - D(1);
        else                % Move x origin to the west
            new_w = fix(dx_w) * head(8);
            new_e = new_w + D(2) - D(1);
        end
        str_R = [' -R' num2str(new_w,10) '/' num2str(D(3),10) '/' num2str(new_e,10) '/' num2str(D(4),10)];
    elseif (adjust_s)        % Need to adjust only south and north borders
        if (in_s_to_n)      % Move y origin to the north
            new_s = fix(dy_s) * (head(9) + 1);
            new_n = new_s + D(4) - D(3);
        else                % Move y origin to the south
            new_s = fix(dy_s) * head(9);
            new_n = new_s + D(4) - D(3);
        end
        str_R = [' -R' num2str(D(1),10) '/' num2str(new_s,10) '/' num2str(D(2),10) '/' num2str(new_n,10)];
    else
        errordlg('Asneira desconhecida (Unknown error).','Error')
    end
end

if (adjust_x_inc)       % Nao vou testar/usar o y_inc. A malha tem (tera?) de ser quadrada
    str_I = [' -I' num2str(new_x_inc,10)];
end

if (~isempty(str_R))        % Finer grid needs adjustment. Do it and return.
    str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
    [fName,pName] = put_or_get_file(handles,str1,'Select GMT grid','get');
    if (isempty(fName)),     return;     end
    fName = [fName '=6'];       % I want them in Surfer format
    str = ['grdsample ' [PathName FileName] str_R ' -G' [pName fName]];
    if (~isempty(str_I)),   str = [str str_I];    end
    if isunix,      s = unix(str);
    elseif ispc,    s = dos(str);
    else            errordlg('Unknown platform.','Error');    end
    if ~(isequal(s,0))                  % An error as occured
        errordlg('Error running grdsample. Finer grid was not adjusted.','Error')
    end
    return
end

% The following Tag is very important to destinguish from MB tracks, which have Tags = MBtrack#
% West border
ny = fix((D(4) - D(3)) / head(9)) + 1;
x1 = repmat(D(1),1,ny);
y1 = linspace(D(3),D(3)+(ny-1)*head(9),ny);     % y_min -> y_max
hold on;
lineHand_w = plot(x1,y1,'-o','Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,...
    'MarkerEdgeColor','w','MarkerFaceColor','k', 'MarkerSize',4, 'Tag','Maregraph');
draw_funs(lineHand_w,'line_uicontext')        % Set lines's uicontextmenu
% South border
nx = fix((D(2) - D(1)) / head(8)) + 1;
x2 = linspace(D(1),D(1)+(nx-1)*head(8),nx);     % x_min -> x_max
y2 = repmat(D(3),1,nx);
lineHand_s = plot(x2,y2,'-o','Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,...
    'MarkerEdgeColor','w','MarkerFaceColor','k', 'MarkerSize',4, 'Tag','Maregraph');
draw_funs(lineHand_s,'line_uicontext')        % Set lines's uicontextmenu
% East border
x3 = repmat(D(2),1,ny);
lineHand_e = plot(x3,y1,'-o','Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,...
    'MarkerEdgeColor','w','MarkerFaceColor','k', 'MarkerSize',4, 'Tag','Maregraph');
draw_funs(lineHand_e,'line_uicontext')        % Set lines's uicontextmenu
% North border
y3 = repmat(D(4),1,nx);
lineHand_n = plot(x2,y3,'-o','Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,...
    'MarkerEdgeColor','w','MarkerFaceColor','k', 'MarkerSize',4, 'Tag','Maregraph');
draw_funs(lineHand_n,'line_uicontext')        % Set lines's uicontextmenu
hold off;

% Now we have to compute the index of the maregraphs positions on the finer grid borders
ind.x = inc_x_ratio * (0:length(x2)-1);
ind.y = inc_y_ratio * (0:length(y1)-1);
% Save the index in the corresponding line handles. This way, the user may delete the
% line maregraphs that he is not interested in. The remaining maregraphs will be fished
% out in swan_options, where a tsun2.dat file will be created
set(lineHand_w,'UserData',ind);     setappdata(lineHand_w,'Side','W')
set(lineHand_s,'UserData',ind);     setappdata(lineHand_s,'Side','S')
set(lineHand_e,'UserData',ind);     setappdata(lineHand_e,'Side','E')
set(lineHand_n,'UserData',ind);     setappdata(lineHand_n,'Side','N')

% --------------------------------------------------------------------
function GeophysicsTsun2_Callback(hObject, eventdata, handles, opt)
if (aux_funs('msg_dlg',14,handles));     return;      end
if (nargin == 3), opt = [];     end

if (strcmp(opt,'write_params'))     % Just do what it says and return
    % See if we have lines of maregraphs corresponding to the finer grid edges
    h_mareg = findobj('Type','line','Tag','Maregraph');
    side = [];
    if (~isempty(h_mareg))
        for i=1:length(h_mareg)
            side(i) = getappdata(h_mareg(i),'Side');
            ind(i) = get(h_mareg(i),'UserData');
        end
    else
        errordlg('You don''t have any maregraphs. Tsun2 needs to be feed with water height trhough maregraphs.','Error')
        return
    end

    if (length(side) > 2)
        errordlg('Tsunamis cannot arrive at more than two edges of the finer grid. (If you are not convinced, think a bit more)','Error')
        return
    elseif (length(side) == 2 && strcmp(side(1),'W') && strcmp(side(1),'E'))
        errordlg('This is completly idiot. The wave arrives on the West & East borders?','Chico Clever')
        return    
    elseif (length(side) == 2 && strcmp(side(1),'S') && strcmp(side(1),'N'))
        errordlg('This is completly idiot. The wave arrives on the North & South borders?','Chico Clever')
        return    
    elseif (isempty(side))
        errordlg('These maregraphs are not of the correct type to use with the tsun2 code.','Error')
        return
    end

    str1 = {'*.par', 'params file (*.par)';'*.*', 'All Files (*.*)'};
    [FileName,PathName] = put_or_get_file(handles,str1,'Select tsun2 params file','put');
    if isequal(FileName,0);     return;     end
    fid = fopen([PathName FileName],'wt');
    fprintf(fid,'%s\n%s\n','# Mirone generated tsun2 parameter file','#');

    % OK, we have only one or two edges. Lets find out whitch are they
    if (length(side) == 1 && strcmp(opt,'write_params'))     % Single edge
        switch side
            case 'W',       fprintf(fid,'%s\n','W');     fprintf(fid,'%d ',ind.y);
            case 'S',       fprintf(fid,'%s\n','S');     fprintf(fid,'%d ',ind.x);
            case 'E',       fprintf(fid,'%s\n','E');     fprintf(fid,'%d ',ind.y);
            case 'N',       fprintf(fid,'%s\n','N');     fprintf(fid,'%d ',ind.x);
        end
    elseif (length(side) == 2 && strcmp(opt,'write_params'))
        if (strcmp(char(side(1)),'S') && strcmp(char(side(2)),'W'))
            fprintf(fid,'%s %s\n','W','S');
            fprintf(fid,'%d ',ind(2).y);   fprintf(fid,'\n');      fprintf(fid,'%d ',ind(1).x);
        elseif (strcmp(char(side(1)),'E') && strcmp(char(side(2)),'S'))
            fprintf(fid,'%s %s\n','S','E');
            fprintf(fid,'%d ',ind(2).x);   fprintf(fid,'\n');      fprintf(fid,'%d ',ind(1).y);
        elseif (strcmp(char(side(1)),'N') && strcmp(char(side(2)),'E'))
            fprintf(fid,'%s %s\n','E','N');
            fprintf(fid,'%d ',ind(2).y);   fprintf(fid,'\n');      fprintf(fid,'%d ',ind(1).x);
        elseif (strcmp(char(side(1)),'N') && strcmp(char(side(2)),'W'))
            fprintf(fid,'%s %s\n','W','N');
            fprintf(fid,'%d ',ind(2).y);   fprintf(fid,'\n');      fprintf(fid,'%d ',ind(1).x);
        end
    end
    fclose(fid);
    return
end         % END of write_params

out = swan_options('Tsun2',handles.figure1);
if isempty(out);    return;     end
pause(0.05);        % Give time to swan_options window to die

extra_args1 = ' ';   extra_args2 = ' ';   extra_args3 = ' ';
extra_args4 = ' ';   extra_args5 = ' ';   extra_args6 = ' ';   extra_args7 = ' ';
if (isfield(out,'maregraph_xy') && isfield(out,'params_file_name'))
    extra_args1 = ['-P' out.params_file_name];
end
if (isfield(out,'opt_M')),  extra_args2 = out.opt_M;   end
if (isfield(out,'opt_N')),  extra_args3 = out.opt_N;   end
if (isfield(out,'opt_m')),  extra_args4 = out.opt_m;   end
if (isfield(out,'opt_G')),  extra_args5 = out.opt_G;   end
if (isfield(out,'opt_I')),   extra_args6 = out.opt_I;   end

[X,Y,Z_bat,head_bat] = load_grd(handles);   Z_bat = double(Z_bat);
if (isempty(Z_bat) || isempty(head_bat))
    errordlg('ERROR: one or more of the bat variables are empty where they souldn''t be.','Error');
    return;
end
if ( abs(head_bat(8) - head_bat(9)) > 1e-3 )
    warndlg('Grid cells are not square. I don''t know the effect of this.','Warning')
end

if (isfield(out,'maregraph_xy'))    % Ask for computation of maregraphs (WRONG - This is not an option, but will be in future)
    dt1 = diff(out.maregraph_xy(:,1));    dt2 = diff(dt1);      t0 = out.maregraph_xy(1,1);
    if any(dt2 ~= 0)        % First column doesn't have the time
        warndlg('The maregraph file does not have a time increment. I will assume it is 1 second.','SEVERE WARNING')
        extra_args6 = ['-I' num2str(head_bat(8)) '/1'];
    else                    % First column of maregs file has the time. We don't want it
        dt = dt1(1);        % This is the time increment to be used as option to tsun2
        out.maregraph_xy(:,1) = [];
        cfl = head_bat(8) / sqrt(abs(head_bat(5))*9.8);
        if (cfl < dt)
            warndlg('Your wave propagates faster than the maregraph time increment. Expect divergent results.','SEVERE WARNING')
        end
        extra_args6 = ['-I' num2str(head_bat(8)) '/' num2str(dt)];
        if (t0 > 1),    extra_args7 = ['-J' num2str(t0)];   end     % If we start computin at a latter time
    end
    if (isfield(out,'opt_m'))       % Movie option
        movie = tsun2(Z_bat, head_bat, out.maregraph_xy, extra_args1, ...
            extra_args2,extra_args3,extra_args5,extra_args6,extra_args7,'-f');
    else
        tsun2(Z_bat, head_bat, out.maregraph_xy, extra_args1,extra_args2,extra_args3,extra_args5,extra_args6,extra_args7);
    end    
else                                % Compute grids or movie
    if (isfield(out,'opt_m'))       % Movie option
        movie = tsun2(Z_bat, head_bat, extra_args2,extra_args3,extra_args5,extra_args6,'-f');
    else
        tsun2(Z_bat, head_bat, extra_args2,extra_args3,extra_args5,extra_args6);
    end
end
pause(0.01);        % Give time to waitbar window to die

if (isfield(out,'opt_m')),   do_movie(handles,movie,'tsun2');   end

%------------------------------------------------------------------------------------------
function do_movie(handles,movie,opt)
if (strcmp(opt,'swan')),    is_swan = 1;
else                        is_swan = 0;    end
[m,n,k] = size(movie);
opt_E = '-E60/30/0.55/0.6/0.4/10';      % Should be "controlable"
n_crop = 0;    % Should be a variabe. It's used to crop borders and therefore hide grid borders reflections
[X,Y,Z0,head] = load_grd(handles);
if isempty(Z0),   return;     end;    % An error message was already issued
if (n_crop)
    Z0 = double(Z0(n_crop+1:m-n_crop, n_crop+1:n-n_crop));
else
    Z0 = double(Z0);
end
h = waitbar(0,'Wait again (computing movie)');

cmap = flipud(hot(256));
cmap(1,:) = [0 0 0.8];      cmap(2,:) = [0 0 0.8];
idx0 = (Z0 <= 0);           Z0(idx0) = 0;
if (is_swan)
    head(5) = 0;
    R0 = grdgradient_m(Z0,head,opt_E);
    R0 = flipud(R0);    idx0 = flipud(idx0);
else    clear idx0;     end
Z0 = flipud(Z0);

i = 1;
while (i <= k)
    Z = movie(:,:,1);
    movie(:,:,1) = [];      % Free this page (we don't need it anymore)
    if (n_crop),    Z = double(Z(n_crop+1:m-n_crop, n_crop+1:n-n_crop));
    else            Z = double(Z);    end
    z_max = max(max(Z));  z_min = min(min(Z));
    head(5) = z_min;      head(6) = z_max;
    R = grdgradient_m(Z,head,opt_E);
    if (is_swan)
        Z(~idx0) = Z0(~idx0);        R(~idx0) = R0(~idx0);
        z_max = head(6);
    end
    
    dif_z = abs(Z - Z0);    idx = (dif_z > 0.05);
    if (n_crop),    img = repmat(uint8(0),m-2*n_crop,n-2*n_crop);
    else            img = repmat(uint8(0),m,n);     end
    img(~idx) = uint8(round( (Z(~idx) / z_max)*254 + 1 ));
    img(idx) = 1;
    idx = (Z == 0);    img(idx) = 1;
    
    img = ind2rgb8(img,cmap);
    img = shading_mat(img,R);
    M(i) = im2frame(img);
    waitbar(i/k)
    i = i + 1;
end
close(h)

str1 = {'*.avi;*.AVI', 'avi files (*.avi,*.AVI)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select movie name','put');
if isequal(FileName,0);     return;     end
set(handles.figure1,'pointer','watch')
movie2avi_j(M,[PathName FileName],'compression','none','fps',5)
set(handles.figure1,'pointer','arrow');

% --------------------------------------------------------------------
function GeophysicsGeogCalculator_Callback(hObject, eventdata, handles)
% Call the geographic calculator. There is no output when called like this
geog_calculator(handles.figure1);

% --------------------------------------------------------------------
function GeophysicsEulerStuff_Callback(hObject, eventdata, handles, opt)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
if (strcmp(opt,'CompEulerPole'))
    compute_euler(handles.figure1)      % Launch a new GUI to do it
elseif (strcmp(opt,'AdjustEulerPole'))
    manual_pole_adjust(handles.figure1) %       "
end

% --------------------------------------------------------------------
function GRD_save_or_display(handles,X,Y,Z,head,tit,name)
% Choose what to do with Z. Show it in a new window or save it as a gmt grid.
% TIT is used only when saving a grid
if (nargin < 7),    name = [];  end
if (handles.out_in_NewWindow)
    if (isa(Z,'double')),   Z = single(Z);  end
    head(5) = double(min(Z(:)));        head(6) = double(max(Z(:)));
    tmp.head = [head(1) head(2) head(3) head(4) head(5) head(6) head(7) head(8) head(9)];
    tmp.X = X;    tmp.Y = Y;    tmp.name = name;
    mirone(Z,tmp);
else
    str1 = {'*.grd;*.GRD','netCDF grid format (*.grd,*.GRD)'; '*.*', 'All Files (*.*)'};
    [FileName,PathName] = put_or_get_file(handles,str1,'Select output GMT grid','put');
    if isequal(FileName,0);     return;     end

    [PATH,FNAME,EXT] = fileparts([PathName FileName]);
    if (handles.grdformat == 6)       % Preserve the original format
        if isempty(EXT),    f_name = [PathName FNAME '.grd=' num2str(handles.grdformat)];
        else                f_name = [PathName FNAME EXT '=' num2str(handles.grdformat)];       end
    else
        if isempty(EXT),    f_name = [PathName FNAME '.grd=nf'];
        else                f_name = [PathName FNAME EXT '=nf'];       end    
    end
    grdwrite_m(Z,head,f_name,tit);     pause(0.01);
end
set(handles.figure1,'pointer','arrow');

% --------------------------------------------------------------------
function FileSaveImgGrdGdal_Callback(hObject, eventdata, handles, opt1, opt2)
% OPT1 = DRIVER == GTiff, HFA (erdas), ENVI, ECW, JP2ECW
% OPT2 == img  -> saves the displayed image;        OPT2 == screen  -> do a screen capture
% OPT2 == grid -> saves the underlaying grid
if (handles.no_file == 1),    return;      end
if (strcmp(opt2,'grid') && ~handles.ValidGrid)
    errordlg('You don''t have a Grid loaded, so OBVIOUSLY you cannot save it.','Error');  return
end
flip = 0;           % FLIP = 1, when images need to be UD fliped before saving
switch opt1
    case 'generic',     SaveGenericImage(handles);  return
    case 'GeoTiff',     str1 = {'*.tif;*.TIF;*.tiff;*.TIFF', 'GeoTiff (*.tif;*.TIF;*.tiff;*.TIFF)'};      driver = 'GTiff';
    case 'Erdas',       str1 = {'*.img;*.IMG', 'Erdas (*.img;*.IMG)'};          driver = 'HFA';
    case 'Envi',        str1 = {'*.img;*.IMG', 'Envi (*.img;*.IMG)'};           driver = 'ENVI';
    case 'ESRI',        str1 = {'*.bil;*.BIL', 'Envi (*.bil;*.BIL)'};           driver = 'EHdr';
    case 'JP2K',        str1 = {'*.jp2;*.JP2', 'Jpeg2000 (*.jp2;*.JP2)'};       driver = 'JP2ECW';
end
str2 = ['Select ' opt1 ' file name'];
[FileName,PathName] = put_or_get_file(handles,str1,str2,'put');
if isequal(FileName,0);     return;     end
[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if isempty(EXT),    FileName = [FileName str1{1}(2:5)];  end
fname = [PathName FileName];

head = handles.head;
% 'ThisImageLims' contains the limits as seen from the pixel-registration stand-point,
%  and that is the convention that we will use to save rasters using GDAL.
imgLims = getappdata(handles.figure1,'ThisImageLims');
if (strcmp(opt2,'grid'))
    [X,Y,Z,head,m,n] = load_grd(handles);
    if isempty(Z),   return;     end;    % An error message was already issued
    if (handles.was_int16)
        if (handles.have_nans)      % Restore the original Nodata value, or use -32768 if we don't know it
            Z(isnan(Z(:))) = handles.Nodata_int16;
        end
        Z = int16(Z);
    end
elseif (strcmp(opt2,'screen'))          % Do a screen capture
    resizetrue(handles.figure1,'screen_capture');
    img = screen_capture(handles.figure1,'lixo_screen.jpg'); % Only the .jpg really matters
    resizetrue(handles.figure1,'after_screen_capture');
    img = aux_funs('strip_bg_color',handles,img);
    head(8) = (imgLims(2)-imgLims(1)) / size(img,2);
    head(9) = (imgLims(4)-imgLims(3)) / size(img,1);
    flip = 1;
else                                    % 'image'
    img = get(handles.grd_img,'CData');
end

hdr.name = fname;       hdr.driver = driver;    hdr.Geog = handles.geog;
try
	hdr.Xinc = head(8);     hdr.Yinc = head(9);
	hdr.ULx = imgLims(1);   hdr.ULy = imgLims(4);
	hdr.Reg = 1;            hdr.Flip = flip;
catch
    errordlg('Shit, image header was not saved as it should.','Error');     return
end
if ( (strcmp(opt2,'img') || strcmp(opt2,'screen')) && ndims(img) == 2 )
    hdr.Cmap = get(handles.figure1,'ColorMap');
end

% hdr.DstProjSRS = '+proj=utm +zone=29';
% [Z2,hdr2]=gdalvirtual(Z,hdr);
% head = hdr2.GMT_hdr;
% X = linspace(head(1),head(2),size(Z2,2));  Y = linspace(head(3),head(4),size(Z2,1));
% tmp.head = head;    tmp.X = X;      tmp.Y = Y;
% tmp.name = 'projected';
% mirone(Z2,tmp);
% return

if (strcmp(opt2,'grid')),   gdalwrite(Z,hdr)
else                        gdalwrite(img,hdr)
end

% --------------------------------------------------------------------
function [FileName,PathName,handles] = put_or_get_file(handles,str1,str2,type)
% Use this function to select input or output filename
if (strcmp(type,'get'))
    cd(handles.last_dir)
    [FileName,PathName] = uigetfile(str1,str2);
    if (PathName ~= 0),         handles.last_dir = PathName;    end
elseif (strcmp(type,'put'))
    cd(handles.work_dir)
    [FileName,PathName] = uiputfile(str1,str2);
    if (PathName ~= 0),         handles.last_dir = PathName;    end
end
pause(0.01);        cd(handles.home_dir);       % allways go home
if (~isempty(strfind([PathName FileName],' ')))
    errordlg('If you had RTFM you should know that names (path included) with white spaces are totaly FORBIDEN here.','ERROR')
    FileName = 0;   % Make calling routine quit
end
guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function GridToolsGrdfilter_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',14,handles));     return;      end
[X,Y,Z] = load_grd(handles);
if isempty(Z),  return;     end;    % An error message was already issued
grdfilter_Mir(handles.figure1,handles.geog,Z);

% --------------------------------------------------------------------
function GridToolsGrdsample_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',14,handles));     return;      end
out = grdsample_Mir(handles.figure1);
if (isempty(out)),      return;     end
[X,Y,Z,head] = load_grd(handles);
if isempty(Z),  return;     end;    % An error message was already issued

set(handles.figure1,'pointer','watch')
arg1 = ' ';     arg2 = ' ';     arg3 = ' ';     arg4 = ' ';
if (isfield(out,'opt_R')),   arg1 = out.opt_R;   end
if (isfield(out,'opt_N')),   arg2 = out.opt_N;   end
if (isfield(out,'opt_Q')),   arg3 = out.opt_Q;   end
if (isfield(out,'opt_L')),   arg4 = out.opt_L;   end
newZ = grdsample_m(Z,head,arg1,arg2,arg3,arg4);
[zzz] = grdutils(newZ,'-L');    z_min = zzz(1);     z_max = zzz(2);
[ny,nx] = size(newZ);
if (strcmp(arg1,' '))       % Grid limits did not change
    x_inc = (head(2) - head(1)) / (nx - ~head(7));
    y_inc = (head(4) - head(3)) / (ny - ~head(7));
    new_head = [head(1) head(2) head(3) head(4) z_min z_max head(7) x_inc y_inc];
    X = linspace(head(1),head(2),nx);       Y = linspace(head(3),head(4),ny);
else
    x_inc = (out.x_min - out.x_max) / (nx - ~head(7));
    y_inc = (out.y_min - out.y_max) / (ny - ~head(7));
    new_head = [out.x_min out.x_max out.y_min out.y_max z_min z_max head(7) x_inc y_inc];    
    X = linspace(out.x_min,out.x_max,nx);   Y = linspace(out.y_min,out.y_max,ny);
end
GRD_save_or_display(handles,X,Y,newZ,new_head,'Resampled grid','Resampled grid')

% --------------------------------------------------------------------
function GridToolsGrdgradient_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',14,handles));     return;      end
out = grdgradient_Mir(handles.figure1);      pause(0.01);
if (isempty(out)),      return;     end
[X,Y,Z,head] = load_grd(handles);
if isempty(Z),  return;     end;    % An error message was already issued

set(handles.figure1,'pointer','watch')
arg1 = ' ';     arg2 = ' ';     arg3 = ' ';     arg4 = ' ';     arg5 = ' ';
if (isfield(out,'opt_A')),   arg1 = out.opt_A;   end
if (isfield(out,'opt_D')),   arg3 = out.opt_D;   end
if (isfield(out,'opt_N')),   arg2 = out.opt_N;   end
if (isfield(out,'opt_L')),   arg4 = out.opt_L;   end
if (handles.geog),           arg5 = '-M';        end
newZ = grdgradient_m(single(-double(Z)),head,arg1,arg2,arg3,arg4,arg5);     % should be clever
head(5) = double(min(newZ(:)));      head(6) = double(max(newZ(:)));
GRD_save_or_display(handles,X,Y,newZ,head,'Gradient grid','Gradient grid')

% --------------------------------------------------------------------
function GridToolsGrdtrend_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',14,handles));     return;      end
out = grdtrend_Mir(handles.figure1);      pause(0.01);
if (isempty(out)),      return;     end
[X,Y,Z,head] = load_grd(handles);
if isempty(Z),  return;     end;    % An error message was already issued
set(handles.figure1,'pointer','watch')
newZ = grdtrend_m(Z,head,out.opt_what,out.opt_N);
head(5) = double(min(newZ(:)));      head(6) = double(max(newZ(:)));
GRD_save_or_display(handles,X,Y,newZ,head,'Grdtrend grid','Grdtrend grid')

% --------------------------------------------------------------------
function GridToolsGrdproject_Callback(hObject, eventdata, handles)
% Call the geographic calculator in the grid only mode.
if (aux_funs('msg_dlg',14,handles));     return;      end
[X,Y,Z,head] = load_grd(handles);
if isempty(Z),          return;     end;    % An error message was already issued
out = geog_calculator(Z,head,handles.figure1,'onlyGrid');
if (isempty(out)),      return;     end     % User just gave up
Z = out.Z;      head = out.head;
head(5) = double(min(Z(:)));        head(6) = double(max(Z(:)));
GRD_save_or_display(handles,X,Y,Z,head,'Grdprojected grid','Grdprojected grid')

% --------------------------------------------------------------------
function GridToolsGridClip_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',14,handles));     return;      end
[X,Y,Z,head] = load_grd(handles);    %  load gmt grid
if isempty(Z),   return;     end;    % An error message was already issued
out = ml_clip(head(5),head(6));
if isempty(out),    set(handles.figure1,'pointer','arrow');    return;     end;
if ~isempty(out{2}),    Z(Z > out{1}) = out{2};     end;     % Clip above
if ~isempty(out{4}),    Z(Z < out{3}) = out{4};     end;     % Clip below
GRD_save_or_display(handles,X,Y,Z,head,'Cliped grid','Cliped grid')

% --------------------------------------------------------------------
function GridToolsCropGrid_Callback(hObject, eventdata, handles,opt)
if (aux_funs('msg_dlg',14,handles));     return;      end
ImageCrop_Callback([],[],handles,[],'CropaGrid')

% --------------------------------------------------------------------
function GridToolsHistogram_Callback(hObject, eventdata, handles, opt)
% OPT2 if present is a structure with two fields: opt2.Z (the matrix); opt2.head (its 9 elements header)
if (aux_funs('msg_dlg',14,handles));     return;      end
if (nargin == 3)                        % Use entire grid
    [X,Y,Z,head] = load_grd(handles);
    if isempty(Z),  return;     end;    % An error message was already issued
    z_min = head(5);        z_max = head(6);
else                                    % Use a subset grid extracted from a rectangular area
    Z = opt.Z;              z_min = opt.head(5);      z_max = opt.head(6);
end
binwidth = (z_max - z_min) / 20;    % Default to 20 bins
resp = inputdlg({'Enter Bin Width (default is 20 bins)'},'Histogram',[1 38],{num2str(binwidth)});     pause(0.01);
if isempty(resp);    set(handles.figure1,'pointer','arrow');     return;     end
n = round( (z_max - z_min) / str2double(resp{1}) );
[n,xout] = histo_m('hist',Z(:),n,[z_min z_max]);
h = mirone;                         % Create a new Mirone figure
mirone('FileNewBgFrame_Callback',h,[],guidata(h), [xout(1) xout(end) 0 max(n) 0], [600 600]);
set(h,'Name','Grid Histogram');       axes(get(h,'CurrentAxes'));
histo_m('bar',xout,n,'hist');
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function GridToolsGridMask_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',14,handles));     return;      end
if (~handles.have_nans)
    msgbox('This option only works on grids that have NaNs.','Warning');    return;
end
[X,Y,Z,head] = load_grd(handles);    %  load gmt grid
if isempty(Z),   return;     end;    % An error message was already issued
Z(~isnan(Z)) = 1;
GRD_save_or_display(handles,X,Y,Z,head,'Mask grid','Mask grid')

% --------------------------------------------------------------------
function GridToolsSectrum_Callback(hObject, eventdata, handles, opt1, opt2)
% OPT1 == 'Amplitude'   -> compute amplitude spectrum
% OPT1 == 'Power'       -> compute power spectrum
% OPT1 == 'Autocorr'    -> compute autocorrelation
% OPT1 == 'Allopts'     -> call the fft_stuff window
% OPT2 if present is a structure with two fields: opt2.Z (the matrix); opt2.head (its 9 elements header)
if (aux_funs('msg_dlg',14,handles));     return;      end
if (nargin == 4)                        % Use entire grid
    [X,Y,Z,head] = load_grd(handles);   quick = 0;
    if isempty(Z),  return;     end;    % An error message was already issued
else                                    % Use a subset grid extracted from a rectangular area
    Z = opt2.Z;             head = opt2.head;       quick = 1;
end
if (handles.have_nans)
    warndlg('This grid has NaNs. That is not allowed in FFTs','Warning');    return;
end

if (quick),     set(handles.figure1,'pointer','watch'); end
fft_stuff(handles.figure1, Z, head, handles.geog, opt1);
if (quick),     set(handles.figure1,'pointer','arrow'); end

% --------------------------------------------------------------------
function GridToolsSmooth_Callback(hObject, eventdata, handles)
if (aux_funs('msg_dlg',14,handles));     return;      end
[X,Y,Z,head,m,n] = load_grd(handles);       %  If gmt grid is not in memory, so read it again
if isempty(Z),      return;     end;        % An error message was already issued
if (~isa(Z,'double')),  Z = double(Z);  end;            % Make sure Z is of double type

[pp p_guess] = spl_fun('csaps',{Y(1:5),X(1:5)},Z(1:5,1:5));   % Get a good estimate of p
prompt = {'Enter smoothing p paramer'};     dlg_title = 'Smoothing parameter input';
defAns = {num2str(p_guess{1},12)};     resp  = inputdlg(prompt,dlg_title,[1 38],defAns);
pause(0.01)
if isempty(resp);    return;     end

set(handles.figure1,'pointer','watch')
Lim = handles.grdMaxSize*.6;
nl0 = round(Lim/(n*16*8));
if (rem(nl0,2) ~= 0), nl0 = nl0 + 1;  end   % Don't know why, but I rather have it as a even number
nl = round(nl0*.8);
if (rem(nl,2) ~= 0),    nl = nl + 1;  end
skirt = ceil(nl0 - nl) - 2;

[ind_s,ind] = tile(m,nl,skirt);             % Get indexes for tiling.
if size(ind_s,1) > 1                        % There is still a very strange thing that I don't understand.
    Zs = [];
    for k = 1:size(ind_s,1)
        tmp1 = (ind_s(k,1):ind_s(k,2));     % Indexes with overlapping zone
        tmp2 = ind(k,1):ind(k,2);           % Indexes of chunks without the overlaping zone
        pp = spl_fun('csaps',{Y(tmp1),X},Z(tmp1,:),str2double(resp{1}));
        tmp = spl_fun('fnval',pp,{Y(tmp1),X});
        Zs = [Zs; tmp(tmp2,:)];
    end
    clear pp tmp;
    Z = Zs;     clear Zs;
else
    pp = spl_fun('csaps',{Y,X},Z,str2double(resp{1}));
    Z = spl_fun('fnval',pp,{Y,X});    clear pp;
end

tit = ['Spline smoothed grid. p parameter used = ' str2double(resp{1})];
GRD_save_or_display(handles,X,Y,Z,head,tit,'Spline smoothed grid');

% --------------------------------------------------------------------
function GridToolsSDG_Callback(hObject, eventdata, handles, opt)
if (aux_funs('msg_dlg',14,handles));     return;      end

[X,Y,Z,head] = load_grd(handles);           [m,n] = size(Z);
if isempty(Z),   return;     end;    % An error message was already issued
if (~isa(Z,'double')),  Z = double(Z);  end;            % Make sure Z is of double type
[pp p_guess] = spl_fun('csaps',{Y(1:5),X(1:5)},Z(1:5,1:5));   % Get a good estimate of p
prompt = {'Enter smoothing p paramer'};     dlg_title = 'Smoothing parameter input';
defAns = {num2str(p_guess{1},12)};     resp  = inputdlg(prompt,dlg_title,[1 38],defAns);
pause(0.01)
if isempty(resp),     return;     end

% Apparently the biggest memory monster (pp) obeys roughly to the following relation:
% n_row x 4 x n_column x 4 * 8. So if I impose a limit "Lim" to this monster, I should
% be able to find out a chunk height that fits to the above relation. However, 'Lim'
% should be of the order of handles.grdMaxSize because the tiling is meant to be used only
% for large grids (relative to the available ram). Also, I don't understand why contiguous
% chunks don't patch perfectly (very easealy seen with the Lambertian illumination). From
% what I could find, it depends havely on the p parameter. A skirt 20% of the chunk height
% seams to produce a reasonable result for high p's, but not allways perfect.

set(handles.figure1,'pointer','watch')
Lim = handles.grdMaxSize*.5;
nl0 = round(Lim/(n*16*8));
if rem(nl0,2) ~= 0, nl0 = nl0 + 1;  end     % Don't know why, but I rather have it as a even number
nl = round(nl0*.8);
if rem(nl,2) ~= 0, nl = nl + 1;  end
skirt = ceil(nl0 - nl) - 2;

[ind_s,ind] = tile(m,nl,skirt);           % Get indexes for tiling.
if size(ind_s,1) > 1                        % There is still a very strange thing that I don't understand.
    R = [];
    for k = 1:size(ind_s,1)
        tmp1 = (ind_s(k,1):ind_s(k,2));     % Indexes with overlapping zone
        tmp2 = ind(k,1):ind(k,2);           % Indexes of chunks without the overlaping zone
        pp = spl_fun('csaps',{Y(tmp1),X},Z(tmp1,:),str2double(resp{1}));
        DfDX = spl_fun('fnder',pp, [1 0]);            % df / dx
        vx = spl_fun('fnval',DfDX,{Y(tmp1),X});       clear DfDX;
        DfDY = spl_fun('fnder',pp, [0 1]);            % df / dy
        vy = spl_fun('fnval',DfDY,{Y(tmp1),X});       clear DfDY;
        D2fDX2 = spl_fun('fnder',pp, [2 0]);          % d^2f / dx^2
        Hxx = spl_fun('fnval',D2fDX2,{Y(tmp1),X});    clear D2fDX2;
        D2fDY2 = spl_fun('fnder',pp, [0 2]);          % d^2f / dy^2
        Hyy = spl_fun('fnval',D2fDY2,{Y(tmp1),X});    clear D2fDY2;
        D2fDXDY = spl_fun('fnder',pp, [1 1]);         clear pp;    % d^2f / (dx dy)
        Hxy = spl_fun('fnval',D2fDXDY,{Y(tmp1),X});   clear D2fDXDY;
        tmp = zeros(ind(k,2),n);
        for i = tmp2
            for j=1:n
                v = [vx(i,j) vy(i,j)] / norm([vx(i,j) vy(i,j)]);        % eq(2)       
                tmp(i,j) = (v * [Hxx(i,j) Hxy(i,j); Hxy(i,j) Hyy(i,j)] * v') / (v*v');
            end
        end
        R = [R; tmp(tmp2,:)];
    end
else
    pp = spl_fun('csaps',{Y,X},Z,str2double(resp{1}));     clear Z;
    DfDX = spl_fun('fnder',pp, [1 0]);        % df / dx
    vx = spl_fun('fnval',DfDX,{Y,X});         clear DfDX;
    DfDY = spl_fun('fnder',pp, [0 1]);        % df / dy
    vy = spl_fun('fnval',DfDY,{Y,X});         clear DfDY;
    D2fDX2 = spl_fun('fnder',pp, [2 0]);      % d^2f / dx^2
    Hxx = spl_fun('fnval',D2fDX2,{Y,X});      clear D2fDX2;
    D2fDY2 = spl_fun('fnder',pp, [0 2]);      % d^2f / dy^2
    Hyy = spl_fun('fnval',D2fDY2,{Y,X});      clear D2fDY2;
    D2fDXDY = spl_fun('fnder',pp, [1 1]);     clear pp;    % d^2f / (dx dy)
    Hxy = spl_fun('fnval',D2fDXDY,{Y,X});     clear D2fDXDY;
    m=length(Y); n=length(X);    R = zeros(m,n);
    for i = 1:m
        for j=1:n
            v = [vx(i,j) vy(i,j)] / norm([vx(i,j) vy(i,j)]);        % eq(2)       
            R(i,j) = (v * [Hxx(i,j) Hxy(i,j); Hxy(i,j) Hyy(i,j)] * v') / (v*v');
        end
    end
    clear Hxx Hxy Hyy vx vy;
end     % end of Tiling

if strcmp(opt,'negative'),      R(R > 0) = 0;
elseif strcmp(opt,'positive'),  R(R < 0) = 0;
end
GRD_save_or_display(handles,X,Y,R,head,'SDG field','SDG field');

% --------------------------------------------------------------------
function GridToolsSlope_Callback(hObject, eventdata, handles, opt)
% OPT == 'degrees'  Compute a DEM slope in degrees
% OPT == 'percent'  Compute a DEM slope in percentage
% OPT == 'aspect'   Compute a DEM aspect in degrees
if (aux_funs('msg_dlg',14,handles));     return;      end
[X,Y,Z,head] = load_grd(handles);
if isempty(Z),      return;     end;            % An error message was already issued
set(handles.figure1,'pointer','watch')
if (~isa(Z,'double')),  Z = double(Z);  end;    % Make sure Z is of double type

[m,n] = size(Z);       D2R = pi/180;   R2D = 180/pi;
if (handles.geog == 1)
    if (strcmp(opt,'aspect'))
        slope = gradient_geo(Y,X,Z,'aspect');   % As you can see "slope" is in fact "aspect"
        tit = 'Terrain aspect in degrees clockwise from North';
    else
        slope = gradient_geo(Y,X,Z,'slope');
        tit   = ['Slope in ' opt];
    end
else
    if (strcmp(opt,'aspect'))
        set(handles.figure1,'pointer','arrow');    warndlg('Not yet programed for cartesian grids','Warning');  return
    end
    [nx,ny] = meshgrid(X,Y);
    [nx,ny,nz] = getnormals(nx,ny,Z);
    vert = [0 0 1];     slope = zeros(m,n);
    for i=1:m
        for j=1:n
            slope(i,j) = acos(sum(vert.*[nx(i,j) ny(i,j) nz(i,j)]));
        end
    end
    tit   = ['Slope in ' opt];
end

if (strcmp(opt,'percent') && handles.geog == 1)
    slope = 100 * tan(slope*D2R);
elseif (strcmp(opt,'percent') && handles.geog == 0)
    slope = 100 * tan(slope);
elseif (strcmp(opt,'degrees') && handles.geog == 0)
    slope = slope * R2D;
end
GRD_save_or_display(handles,X,Y,slope,head,tit,tit);

% --------------------------------------------------------------------
function GridToolsDirectionalDerivative_Callback(hObject, eventdata, handles, opt)
if (aux_funs('msg_dlg',14,handles));     return;      end
if (nargin == 3),   opt = 'first';     end
luz = shading_params('dirDerivative');
pause(0.01)         % Give time to the azimuth window be deleted
if isempty(luz),   return;     end
azim = (90 - luz.azim) * pi/180;

[X,Y,Z,head] = load_grd(handles);
if isempty(Z),   return;     end;    % An error message was already issued
if (~isa(Z,'double')),  Z = double(Z);  end;            % Make sure Z is of double type

set(handles.figure1,'pointer','watch')
if (handles.geog == 1)
    [gradN gradE] = gradient_geo(Y,X,Z,'grad');             % df/dy & df/dx
    if strcmp(opt,'first')
        Z = gradE * cos(azim) + gradN * sin(azim);
        str = 'First derivative';
    else        % second derivative
        Z = gradient_geo(Y,X,gradN,'gradN')*(sin(azim)^2);   clear gradN;       % d2f/dy2
        Z = Z + gradient_geo(Y,X,gradE,'gradE')*(cos(azim)^2);                  % + d2f/dx2
        Z = Z + 2*gradient_geo(Y,X,gradE,'gradN')*(cos(azim) * sin(azim));      % + d2f/dxdy
        str = 'Second derivative';
    end
else
    set(handles.figure1,'pointer','arrow');     warndlg('Not yet programed for cartesian grids','Warning');   return
end
GRD_save_or_display(handles,X,Y,Z,head,str,str);

% --------------------------------------------------------------------
function GridToolsSRTM_mosaic_Callback(hObject, eventdata, handles, opt)
% OPT is either == [] (SRTM3c) or == 'srtm1', or == 'srtm30'
if (nargin == 3),   opt = [];   end
if (isempty(opt)),  srtm_tool
else                srtm_tool(opt); end

% --------------------------------------------------------------------
function GridToolsFindHoles_Callback(hObject, eventdata, handles)
% Find holes in double arrays and draw rectangles arround them
if (aux_funs('msg_dlg',14,handles));     return;      end
[X,Y,Z,head,m,n] = load_grd(handles);
if isempty(Z),      return;     end;    % An error message was already issued
a = isnan(Z);
if (~any(a(:))),    warndlg('This grid has no holes','Warning');    return;  end

set(handles.figure1,'pointer','watch')
B = img_fun('find_holes',a);
% Draw rectangles arround each hole
hold on
lt = getappdata(handles.figure1,'DefLineThick');    lc = getappdata(handles.figure1,'DefLineColor');
for i=1:length(B)
    x_min = min(B{i}(:,2));    x_max = max(B{i}(:,2));
    y_min = min(B{i}(:,1));    y_max = max(B{i}(:,1));
    if (x_min > 5),    x_min = x_min - 5;
    else               x_min = 1;      end
    if (x_max < n-5),  x_max = x_max + 5;
    else               x_max = n;      end
    if (y_min > 5),    y_min = y_min - 5;
    else               y_min = 1;      end
    if (y_max < m-5),  y_max = y_max + 5;
    else               y_max = m;      end
    x_min = head(1) + (x_min-1)*head(8);    x_max = head(1) + (x_max-1)*head(8);
    y_min = head(3) + (y_min-1)*head(9);    y_max = head(3) + (y_max-1)*head(9);
    h = plot([x_min x_min x_max x_max x_min],[y_min y_max y_max y_min y_min],'Color',lc,'LineWidth',lt);
    draw_funs(h,'SRTM_rectangle')       % Set uicontexts
end
hold off;   set(handles.figure1,'pointer','arrow');

% --------------------------------------------------------------------
function GridToolsSaveAsSRTM_Callback(hObject, eventdata, handles)
% Only grids with the same characteristics as SRTM 3c files are allowed to be saved
if (handles.no_file == 1),     return;      end
[X,Y,Z,head] = load_grd(handles);       % No need to test for in-memory Z
if ( ((head(2)-head(1)) - 1) > 1e-6 || ((head(4)-head(3)) - 1) > 1e-6 )
    errordlg('Grid does not cover a 1 degree square','Error');  return
end
if ( ((abs(head(8)-3/3600) > 1e-6) || (abs(head(9)-3/3600) > 1e-6)) && ...
        ((abs(head(8)-1/3600) > 1e-6) || (abs(head(9)-1/3600) > 1e-6)) )
    errordlg('Grid spacing differs from SRTM 1 or 3 arcsec files','Error');    return
end

% Build the file name (appended with the '_p' suffix)
if (sign(head(1)) > 0 ),    w = 'E';
else                        w = 'W';    end
if (sign(head(3)) > 0 ),    n = 'N';
else                        n = 'S';    end
name = [n num2str(abs(round(head(3))),'%.2d') w num2str(abs(round(head(1))),'%.3d') '_p.hgt'];
cd(handles.work_dir)
[FileName,PathName] = uiputfile(name,'Select SRTM File name');
if (PathName ~= 0),         handles.last_dir = PathName;    end
pause(0.01)
if isequal(FileName,0);     return;     end
cd(handles.home_dir);       % allways go home to avoid troubles
[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if (isempty(EXT) && ~isempty(FileName))    FileName = [FileName '.hgt'];     end

Z(isnan(Z)) = -32768;  Z = int16(rot90(Z,-1));      % Reset eventual NaNs to the SRTM nodata value
fid = fopen([PathName FileName],'w','b');
fwrite(fid,Z,'int16');  fclose(fid);
guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function GridToolsPadd2Const_Callback(hObject, eventdata, handles)
% Pad the array to a const value (currently ct = zero) using a Hanning window
if (aux_funs('msg_dlg',14,handles));     return;      end
[X,Y,Z,head,m,n] = load_grd(handles);
if isempty(Z),   return;     end;    % An error message was already issued
resp  = inputdlg({'Enter number of border lines'},'Skirt width',[1 38],{'10'});
pause(0.01)
if isempty(resp);   return;     end
n_pad = str2double(resp{1}) * 2;
Z = mboard(Z,n,m,n+n_pad,m+n_pad);
zzz = grdutils(Z,'-L');  z_min = zzz(1);    z_max = zzz(2);     clear zzz;
head(1) = head(1) - n_pad/2 * head(8);      head(2) = head(2) + n_pad/2 * head(8);
head(3) = head(3) - n_pad/2 * head(9);      head(4) = head(4) + n_pad/2 * head(9);
head(5) = z_min;                            head(6) = z_max;
X = linspace(head(1),head(2),n+n_pad);      Y = linspace(head(3),head(4),m+n_pad);
GRD_save_or_display(handles,X,Y,Z,head,'Padded Grid','Padded Grid')

% --------------------------------------------------------------------
function FileSaveFlederSD_Callback(hObject, eventdata, handles, opt)
% Depending of the OPT value, this function builds either:
% OPT = 'writePlanarSD' directly build a planar Sonar SD file to be used by Fledermaus
% OPT = 'writeSphericalSD' directly build a spherical Sonar SD file to be used by Fledermaus
% OPT = 'runPlanarSD' build a planar .sd file (but don't keep it) and run the viewer
% OPT = 'runSphericalSD' build a spherical .sd file (but don't keep it) and run the viewer
if (aux_funs('msg_dlg',14,handles));     return;      end
if (nargin == 3),   opt = 'runPlanarSD';   end
if ( (strcmp(opt,'writeSphericalSD') || strcmp(opt,'runSphericalSD')) && ~handles.geog)
    errordlg('Spherical objects are allowed only for geographical grids','Error');    return
end

[X,Y,Z,head] = load_grd(handles);
if isempty(Z),  return;     end;    % An error message was already issued

str1 = {'*.sd;*.SD', 'Fledermaus object (*.sd,*.SD)'};
if (strcmp(opt,'writePlanarSD') || strcmp(opt,'writeSphericalSD'))
    [FileName,PathName] = put_or_get_file(handles,str1,'Name of Fledermaus object','put');
    if isequal(FileName,0);     return;     end
    fname = [PathName FileName];
else
    fname = [handles.path_tmp 'lixoSD.sd'];
end

set(handles.figure1,'pointer','watch')
[PATH,FNAME,EXT] = fileparts(fname);
if isempty(EXT),    fname = [fname '.sd'];  end
fid = fopen(fname,'wb');
if (strcmp(opt,'writePlanarSD') || strcmp(opt,'runPlanarSD'))
    write_flederFiles('main_SD',fid,'writePlanarSD',handles.figure1,handles.axes1,Z,head(1:6))
else
    write_flederFiles('main_SD',fid,'writeSphericalSD',handles.figure1,handles.axes1,Z,head(1:6))
end
write_flederFiles('line_or_points',fid,handles.figure1,handles.axes1,Z,head(1:6))
write_flederFiles('eof',fid)        % Write EOF block and close the file
set(handles.figure1,'pointer','arrow')

if (strcmp(opt,'runPlanarSD') || strcmp(opt,'runSphericalSD'))   % Run the viewer and remove the tmp .sd file
    try
		if (isunix)         % Stupid linux doesn't react to a non-existant iview3d
            resp = unix(['iview3d -data ' fname ' &']);
            if (resp == 0)
                errordlg('I could not find the Free Fledermaus viewer. Have you install it?','Error')
            end
		elseif (ispc),  dos(['iview3d -data ' fname ' &']);
		else            errordlg('Unknown platform.','Error');  return;
		end
    catch
        errordlg('I could not find the Free Fledermaus viewer. Did you install it?','Error')
    end
    delete(fname);
end

% --------------------------------------------------------------------
function figure1_ResizeFcn(hObject, eventdata, handles)
screen = get(0,'ScreenSize');
if (~isempty(handles))
	pos = get(handles.figure1,'Pos');
    if (screen(3) == pos(3))        % Do not allow figure miximizing
        set(handles.figure1,'Pos',handles.oldSize)
    else
        handles.oldSize = pos;
        guidata(handles.figure1,handles)
    end
end

% --------------------------------------------------------------------
function ImageEdgeDetect_Callback(hObject, eventdata, handles, opt)
if (handles.no_file == 1),      return;      end

if (~strcmp(opt,'ppa')),        img = get(handles.grd_img,'CData'); end
if (handles.image_type == 1 || handles.image_type == 4)
    [X,Y,Z,handles.head] = load_grd(handles);
end
set(handles.figure1,'pointer','watch');
if (strcmp(opt,'ppa'))
    if ~(handles.image_type == 1 || handles.image_type == 4),   aux_funs('msg_dlg',2,handles);     return;      end
    out = grdppa_m(Z,handles.head);
    h_ridge = line(out(1,:),out(2,:),'Linewidth',handles.DefLineThick,'Color',handles.DefLineColor,'Tag','creast_line','Userdata',1);
    multi_segs_str = cell(length(h_ridge),1);    % Just create a set of empty info strings
    draw_funs(h_ridge,'isochron',multi_segs_str)
    set(handles.figure1,'pointer','arrow')
    return
%     x_lim = get(handles.axes1,'XLim');    y_lim = get(handles.axes1,'YLim');
%     h_lixo = figure('MenuBar','none');
%     h_tmp = line('XData',out(1,:),'YData',out(2,:),'Linewidth',0.1);
%     set(handles.axes1,'XLim',x_lim,'YLim',y_lim)
%     F = getframe(h_lixo);
%     img = F.cdata;
%     [m,n,k] = size(img);
%     x_inc = (handles.head(2) - handles.head(1)) / n;
%     y_inc = (handles.head(4) - handles.head(3)) / m;
%     I = flipud(img(:,:,1));
%     delete(h_lixo);
elseif (strcmp(opt,'Vec') || strcmp(opt,'Ras') || strcmp(opt(1:3),'SUS'))
    if (ndims(img) == 3),   img = cvlib_mex('color',img,'rgb2gray');      end
    if (~strcmp(opt(1:3),'SUS'))
        %img = img_fun('edge',img,'canny');
        img = cvlib_mex('canny',img,40,200,3);
    else
        img = susan(img,'-e');              % Do SUSAN edge detect
        if (strcmp(opt(4:end),'vec')),      opt = 'Vec';        % This avoids some extra tests later
        else                                opt = 'Ras';
        end
    end
    B = img_fun('bwboundaries',img,'noholes');
elseif (strcmp(opt,'Lines'))
    %B = cvlib_mex('houghlines2',img);       % If img == 3D cvlib_mex will take care
    %B = cvlib_mex('houghlines2',img,'standard',1,pi/180,100,0,0);
    %B = cvlib_mex('contours',img);
    %
    BW = cvlib_mex('canny',img,40,200,3);
        
    [H,T,R] = img_fun('hough',BW);
    P = img_fun('houghpeaks',H,50,'threshold',ceil(0.3*double(max(H(:)))));
    lines = img_fun('houghlines',BW,T,R,P,'FillGap',10,'MinLength',50);
    if (~isfield(lines,'point1')),      set(handles.figure1,'pointer','arrow');     return;   end
    for k = 1:length(lines)
        B{k} = [lines(k).point1(2:-1:1); lines(k).point2(2:-1:1)];
    end
elseif (strcmp(opt,'Circles'))
    B = cvlib_mex('houghcircles',img);
    if (isempty(B))
        set(handles.figure1,'pointer','arrow');     return
    end
end

if (strcmp(opt,'Vec') || strcmp(opt,'Lines'))          % Convert the edges found into vector
    x_inc = handles.head(8);    y_inc = handles.head(9);
    x_min = handles.head(1);    y_min = handles.head(3);
    if (handles.head(7))            % Work in grid registration
        x_min = x_min + x_inc/2;    y_min = y_min + y_inc/2;
    end
	h_edge = zeros(length(B),1);    i = 1;
	for k = 1:length(B)
		boundary = B{k};
		if (length(boundary) < 40 && strcmp(opt,'Vec')),   continue;   end
		y = (boundary(:,1)-1)*y_inc + y_min;
		x = (boundary(:,2)-1)*x_inc + x_min;
		%x = x(1:fix(end/2));  y = y(1:fix(end/2));    % This because the stupid has many duplicated points
		h_edge(i) = line(x, y,'Linewidth',handles.DefLineThick,'Color',handles.DefLineColor,'Tag','edge_detected','Userdata',i);
		i = i + 1;
        %ellipse_t = fit_ellipse( x,y,handles.axes1 );
	end
	
	h_edge(h_edge == 0) = [];                   % Remove empty handles remaining from pre-declaration
	multi_segs_str = cell(length(h_edge),1);    % Just create a set of empty info strings
	draw_funs(h_edge,'isochron',multi_segs_str);
elseif (strcmp(opt,'Circles'))
    x = linspace(-pi,pi,360);       y = x;
    x = cos(x);                     y = sin(y);
    lt = getappdata(get(0,'CurrentFigure'),'DefLineThick');    lc = getappdata(get(0,'CurrentFigure'),'DefLineColor');
    h_circ = line('XData', [], 'YData', []);
    for k = 1:size(B,1)
        x = B(k,1) + B(k,3) * x;           y = B(k,2) + B(k,3) * y;
        set(h_circ, 'XData', x, 'YData', y,'Color',lc,'LineWidth',lt,'Userdata',B(k,:));
        draw_funs(h_circ,'SessionRestoreCircleCart')    % Give uicontext
    end
else                            % Display the bw image where the edges are the whites    
    if (strcmp(get(handles.axes1,'Ydir'),'normal')),    img = flipdim(img,1);   end
    setappdata(0,'CropedColormap',gray);
    mirone(img);
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function ImageMovieFromList_Callback(hObject, eventdata, handles)
str1 = {'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.DAT,*.txt,*.TXT)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select file with image list','get');
if isequal(FileName,0);     return;     end
fid = fopen([PathName FileName]);
c = char(fread(fid))';      fclose(fid);
names = strread(c,'%s','delimiter','\n');   clear c fid;
m = length(names);
for (k=1:m)
    j = strfind(names{k},filesep);
    if (isempty(j)), names{k} = [PathName names{k}]; end
    img = imread(names{k});    M(k) = im2frame(img);
end

str1 = {'*.avi;*.AVI', 'avi files (*.avi,*.AVI)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select movie name','put');
if isequal(FileName,0);     return;     end
movie2avi_j(M,[PathName FileName],'compression','none','fps',5)

% --------------------------------------------------------------------
% function ML_Surface_Callback(hObject, eventdata, handles)
% if (handles.no_file == 1)     return;      end
% if (aux_funs('msg_dlg',14,handles));     return;      end
% 
% [X,Y,Z,head] = load_grd(handles);
% if isempty(Z),      return;     end;
% set(handles.figure1,'pointer','watch')
% surface(X,Y,double(Z))
% set(handles.figure1,'pointer','arrow')
% view(-35,45)

% % --------------------------------------------------------------------
% function ToolsMathSurfs_Callback(hObject, eventdata, handles)
% % http://www.winosi.onlinehome.de/Gallery_t12.htm
% %X = -1:0.01:1;
% X = linspace(-pi,pi);
% Y = X;
% [X,Y] = meshgrid(X,Y);
% %Z = cos(sqrt(X.^2 + Y.^2) + pi/2);
% %Z = sqrt(abs(1-abs(1-abs(1-abs(X.*Y+0.5)))));       % Carves
% Z = cos(atan(X./Y).*sign(Y)*8)/4 .* sin(sqrt(X.^2 + Y.^2)*3);    % Rose (erro)
% %Z = (1-cos(X)).*(1-cos(Y));     % Mamas
% %Z = cos(log(sqrt(X.^2 + Y.^2)/30+0.001));       % Peak
% %Z = abs(sin(4*X)+sin(4*Y))/4 .* cos(sqrt(X.^2 + Y.^2)*16);  % Ripples
% %Z = cos(X.^2 + Y.^2) ./ exp((X.^2 + Y.^2)/4);       % Blupp
% %s = sqrt(X.^2 + Y.^2);    Z = sin(1.25*s)+0.48*sin(3.75*s);     % Hat
% %Z = cos(X .* Y);        % Waves around a Cross
% %Z = cos(X) .* cos(Y);   % A ghost
% %Z = -(1-cos(X)).*(1-cos(Y));     % A Crypt
% h = figure;
% surface(X,Y,Z)
% view(-35,45)

% --------------------------------------------------------------------
function DigitalFilt_Callback(hObject, eventdata, handles, opt)
if (handles.no_file == 1),     return;      end
if (strcmp(opt,'image'))
    digitalFiltering(handles.grd_img);
else        % grid
    if (aux_funs('msg_dlg',14,handles));     return;      end
	[X,Y,Z,head] = load_grd(handles);   % load the grid array here
	if isempty(Z),      return;     end;    % An error message was already issued
	[Z, img] = digitalFiltering(handles.grd_img,Z,get(handles.figure1,'ColorMap'));
	if (isempty(Z)),    return;     end
	
	head(5) = double(min(min(Z)));    head(6) = double(max(max(Z)));
	setappdata(handles.figure1,'dem_z',Z);    setappdata(handles.figure1,'GMThead',head);
	setappdata(handles.figure1,'Zmin_max',[head(5) head(6)])
	handles.origFig = img;
	handles.head = head;
	guidata(handles.figure1,handles)
end

% --------------------------------------------------------------------
function RotateTool_Callback(hObject, eventdata, handles, opt)
if (handles.no_file == 1),     return;      end
if (strcmp(opt,'image'))
    img = rotatetool(get(handles.grd_img,'CData'),get(handles.figure1,'ColorMap'));
    if (ndims(img) == 2),   setappdata(0,'CropedColormap',get(handles.figure1,'ColorMap'));     end
    mirone(img);
else        % grid
    if (aux_funs('msg_dlg',14,handles));     return;      end
	[X,Y,Z,head] = load_grd(handles);   % load the grid array here
	if isempty(Z),      return;     end;    % An error message was already issued
	[newZ, hdr] = rotatetool(Z,head,get(handles.figure1,'ColorMap'));
	if (isempty(newZ)),    return;     end
    [ny,nx] = size(newZ);
    X = linspace(hdr(1),hdr(2),nx);       Y = linspace(hdr(3),hdr(4),ny);
    GRD_save_or_display(handles,X,Y,newZ,hdr,'Rotated grid','Rotated grid')
end

% --------------------------------------------------------------------
function ImageEnhance_Callback(hObject, eventdata, handles, opt)
if (handles.no_file == 1),     return;      end
if (strcmp(opt,'1'))    % The extended Adjust Contrast tool
    image_enhance(handles.figure1);
else                    % The reborn old imadjdemo tool
    image_adjust(handles.figure1);
end

% --------------------------------------------------------------------
function ImageGCPtool_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),     return;      end
if (~strcmp(get(hObject,'Checked'),'on'))
    set(hObject,'Checked','on')
else
    set(hObject,'Checked','off')
end

% ----------------------------------------------------------------------------------
function ImageGCPpline_Callback(hObject, eventdata, handles)
if (handles.no_file == 1),  return;     end

[xp,yp] = getline_j(handles.figure1);
n_nodes = length(xp);
if (n_nodes < 2),     return;     end
lineHand = line('XData', xp, 'YData', yp,'Color','k','LineWidth',0.5,'LineStyle',':','Marker','o',...
    'MarkerFaceColor','y','MarkerSize',7,'Tag','GCPpolyline');
register_img(handles,lineHand)

% --------------------------------------------------------------------
function geog = guessGeog(lims)
    geog = double( ( (lims(1) >= -180 && lims(2) <= 180) || (lims(1) >= 0 && lims(2) <= 360) )...
        && (lims(3) >= -90 || lims(4) <= 90) );

% --------------------------------------------------------------------
function Transfer_Callback(hObject, eventdata, handles, opt)
if (handles.no_file == 1),      return;      end

if (strcmp(opt,'Shape')),       floodFill(handles.figure1);     return;     end
if (strcmp(opt,'ImgResize')),   imageResize(handles.figure1);   return;     end
set(handles.figure1,'pointer','watch')
img = get(handles.grd_img,'CData');
if (strcmp(opt,'Corners'))
    corn = cvlib_mex('goodfeatures',img,100,0.05);
    hold on
    lineHand = plot(corn(:,1),corn(:,2),'ko','MarkerEdgeColor','w','MarkerFaceColor','k', ...
                    'MarkerSize',6,'Tag','corner_detected','Userdata',1);
    hold off
	multi_segs_str = cell(length(lineHand),1);    % Just create a set of empty info strings
	draw_funs(lineHand,'isochron',multi_segs_str);
elseif (strcmp(opt,'gray'))
    if (ndims(img) == 3)
        img = cvlib_mex('color',img,'rgb2gray');
        set(handles.grd_img,'CData', img);
    end
    set(handles.figure1,'ColorMap',gray(256))
elseif (strcmp(opt,'bw'))
    img = img_fun('im2bw',img);                             % UGLY resorces consumption by IPTbx
    set(handles.grd_img,'CData', img, 'CDataMapping','scaled');
    set(handles.figure1,'ColorMap',gray(256))
elseif ( strcmp(opt,'isGDAL') || (strcmp(opt,'isGMT')) )    % Test if GDAL or GMT are able to read this file
    str = {'*.*', 'All Files (*.*)'};                       % Just a strings to the eye
    if (strcmp(opt,'isGMT'))
        str = {'*.grd;*.GRD;*.nc;*.NC', 'Grid files (*.grd,*.GRD,*.nc,*.NC)'; '*.*', 'All Files (*.*)'};
    end
    [FileName,PathName] = put_or_get_file(handles,str,['Select ' opt(3:end) ' file'],'get');
    if isequal(FileName,0);     set(handles.figure1,'pointer','arrow');     return;     end             % User gave up
    if (strcmp(opt,'isGDAL'))                               % Again, ...
        [z,att] = gdalread([PathName FileName]);            % See if GDAL can do it ...
    else
        str = ['grdinfo ' [PathName FileName]];             % or GMT ...
        [s,att] = mat_lyies(str,[handles.path_tmp FileName '.info']);
        if ~(isequal(s,0)),     att = [];       end         % File could not be read
    end
    if (isempty(att)),      msg = ['Sorry, ' opt(3:end) ' Was not able to read the given file.'];
    else                    msg = ['Yeap, ' opt(3:end) ' can read the given file (which is not the same thing as ', ...
                                    'beeing able to correctly decode it).'];
    end
    msgbox(msg,'Test result')
end
set(handles.figure1,'pointer','arrow')

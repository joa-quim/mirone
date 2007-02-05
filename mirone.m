function varargout = mirone(varargin)
%   MIRONE, by itself, creates a window bar from which you load a lot of grid/images formats
%   MIRONE(FNAME) opens the file FNAME and displays it
%   H = MIRONE(...) returns the handle to a new mirone window
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
    if (nargout),   varargout{1} = varargin{4}.figure1;    end
else
    h = mirone_OpeningFcn(varargin{:});
    if (nargout),   varargout{1} = h;       end
end

function hObject = mirone_OpeningFcn(varargin)
[hObject,version7,IamCompiled] = mirone_uis;
handles = guihandles(hObject);

% PRAGMA SECTION (It's far far from clear when files must be declared here)
%#function uigetfolder_standalone mapproject_m grdproject_m coordinate_system surface_m
%#function nearneighbor_m cpt2cmap grdfilter_m grdgradient_m grdsample_m grdtrack_m grdtrend_m 
%#function grdutils scaleto8 waitbar bpass3d inv3d nskew rtp3d syn3d igrf_m
%----- These are for okada tsunamis
%#function range_change swan tsun2 mansinha_m
%----- These are for image
%#function grayto8 grayto16 grayxform imfilter_mex imhistc imlincombc parityscan uintlutc ordf
%#function imreconstructmex applylutc bwboundariesmex bwlabel1 bwlabel2
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
%#function patch_meca ui_edit_patch_special bands_list multibandread_j imscroll_j transform_fun iptchecknargin
%#function load_defFilters mltable_j iptcheckinput resampsep intmax wgifc telhometro vitrinite edit_line
%#function edit_track_mb save_track_mb houghmex qhullmx uisuspend_fig uirestore_fig

global home_dir;    home_dir = pwd;

if (version7),  pos = get(hObject,'Pos');   set(hObject,'Pos',[pos(1:3) 1]);    end
movegui(hObject,'north');           % Reposition the window on screen
set(0,'CurrentFigure',hObject)      % Due to a R2006a incredible BUG
set(hObject,'Visible','on','HandleVisibility','callback');

% The addpath command cannot be compiled, so the paths bellow have be added manually
% to the matlab path. However, in the distribution version the next line(s) must be
% uncommented in order that the user doesn't have to know that those need to be added.
% addpath([home_dir filesep 'src_figs']);
% addpath([home_dir filesep 'lib_mex']);
% addpath([home_dir filesep 'utils']);

handles.home_dir = home_dir;
handles.version7 = version7;% If == 1 => R14 or latter
handles.IamCompiled = IamCompiled; % If == 1 than we know that we are dealing with a compiled (V3) version
handles.DefLineThick = 1;   % Default line thickness (overwriten by mirone_pref)
handles.DefLineColor = 'k'; % Default line color (overwriten by mirone_pref)
handles.DefineMeasureUnit = 'k'; % Default measure units to kilometrs (overwriten by mirone_pref)
handles.grdname = [];       % Contains the name of the current (if it's the case) gmt grid
handles.h_MBplot  = [];     % Handles to multi-beam tracks
handles.nTrack = 0;         % Counter of the number of MB tracks imported
handles.hist = 0;           % 0 -> linear colormap; 1 -> histogram equalized
handles.origFig = [];       % To store the original image copy
handles.fileName = [];      % To store any input grid/image file name
handles.grdformat = 0;      % Flag to signal the format of the imported gmt grid (defaults to dafault)
handles.image_type = 0;     % Image type. 1->grd; 2-> trivial (jpg,png,bmp,etc...); 3->GeoTIFF; 4->DEMs; 5->
handles.computed_grid = 0;  % However, matrices with a gmt header will have this == 1, so that they can be saved
handles.no_file = 1;        % 0 means a grid is loaded and 1 that it is not (to test when icons can be pushed)
handles.geog = 1;           % By default grids are assumed to be in geographical coordinates
handles.swathRatio = 3;     % Default swath width / water depth ratio for multibeam planing
handles.grdMaxSize = 20971520;   % I use this for limiting the grid size that is stored in RAM (20 Mb)
handles.firstMBtrack = 1;   % Used for knowing whether to display or not the MB planing info message in "start planing"
handles.EarthRad = 6371;    % Authalic radius
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
handles.hImg = [];          % New name of the image handle that wiil eventually replace the above one
handles.firstIllum = 1;     % First illumination will use the displayed image which may have been IP
handles.flederBurn = 1;     % When build a fleder obj burn eventual coastlines
handles.flederPlanar = 1;   % Default to planar flder objs (but mirone_pref knows better)
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
    handles.flederBurn = flederBurn;
    handles.flederPlanar = flederPlanar;
end

j = false(1,length(handles.last_directories));          % vector for eventual cleaning non-existing dirs
for (i = 1:length(handles.last_directories))            % Check that all dirs in last_directories exist
    try         cd(handles.last_directories{i});        % NOTE. I don't use 'exist' anymore because
    catch       j(i) = 1;       cd(home_dir);           % the stupid compiler allways return something > 0
    end
end
handles.last_directories(j) = [];                       % clean non-existing directories
cd(home_dir);               % Come back home because it was probably somewhere out there

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

% Detect in which mode Mirone was called
drv = [];   grd_data_in = 0;    grd_data_interfero_in = 0;  grd_data_deform_in = 0;     win_name = 'Mirone';
if ~isempty(varargin)
    if (length(varargin) == 1 && ischar(varargin{1}))               % Called with a file name as argument
        [pato, fname, EXT] = fileparts(varargin{1});                % Test to check online command input
        if (isempty(pato)),     varargin{1} = [handles.home_dir filesep fname EXT];     end
        drv = aux_funs('findFileType',varargin{1});
    elseif ( isa(varargin{1},'uint8') || isa(varargin{1},'logical') )
        % Called with an image as argument and optionaly an struct header (& geog, name, cmap optional fields)
        dims = size(varargin{1});
        if ( length(varargin) == 2 && isstruct(varargin{2}) )       % An image with coordinates
            tmp = varargin{2};
            handles.head = tmp.head;        X = tmp.X;      Y = tmp.Y;
            handles.image_type = 3;         axis_t = 'xy';
            if (isfield(tmp,'geog')),       handles.geog = tmp.geog;    end % Prevails over the guess in show_image
            if (isfield(tmp,'cmap')),       set(handles.figure1,'Colormap',tmp.cmap);   end
            if (isfield(tmp,'name')),       win_name = tmp.name;    end
        else
            X = [];         Y = [];         win_name = 'Cropped Image';
            handles.image_type = 2;         handles.geog = 0;       axis_t = 'off';
            handles.head = [1 dims(2) 1 dims(1) 0 255 0 1 1];       % Fake a grid reg GMT header
            if (ndims(varargin{1}) == 2),   set(handles.figure1,'Colormap',gray(256));  end
            pal = getappdata(0,'CropedColormap');                   % See if we have a colormap to use here
            if (~isempty(pal)),     set(handles.figure1,'Colormap',pal);    rmappdata(0,'CropedColormap');  end
            setappdata(hObject,'Croped','yes');                 % ???
        end
        handles = show_image(handles,win_name,X,Y,varargin{1},0,axis_t,handles.head(7),1);
        grid_info(handles,[],'iminfo',varargin{1});         % Contruct a info string
        if (isa(varargin{1},'logical'))
            set(handles.hImg,'CDataMapping','scaled');   set(handles.figure1,'ColorMap',gray(256));
        end        
    elseif (length(varargin) == 2 && isnumeric(varargin{1}) && isstruct(varargin{2}))
        % A matrix. Treat it as if it is a gmt grid. No error testing on the grid head descriptor
        grd_data_in = 1;
        Z = varargin{1};            tmp = varargin{2};
        handles.have_nans = grdutils(Z,'-N');
        handles.head = tmp.head;    X = tmp.X;  Y = tmp.Y;
        if (isfield(tmp,'name')),   win_name = tmp.name;    end    % All calls should transmit a name, but ...
        if (isfield(tmp,'was_int16'))
            handles.was_int16 = tmp.was_int16;      handles.Nodata_int16 = tmp.Nodata_int16;
        end
        clear tmp;
    elseif ( length(varargin) == 4 && isnumeric(varargin{1}) && isstruct(varargin{2}) && ...
            strcmp(varargin{3},'Deformation') && ishandle(varargin{4}) )
        % A matrix. Treat it as if it'is a gmt grid. No error testing on the grid head descriptor
        % Note: this is a special case of the situation above that will be used to identify this figure
        % as an Okada deformtion data (via its Name). This info is searched by the tsunami modeling option
        grd_data_deform_in = 1;
        Z = varargin{1};            tmp = varargin{2};
        handles.head = tmp.head;    X = tmp.X;  Y = tmp.Y;  clear tmp;
        handles.calling_figure = varargin{4};
        win_name = 'Okada deformation';
    elseif ( length(varargin) == 4 && isnumeric(varargin{1}) && isstruct(varargin{2}) && ...
            strcmp(varargin{3},'Interfero') && isnumeric(varargin{4}) )
        % A matrix input containing an interfeogram with cdo == varargin{4}
        grd_data_interfero_in = 1;
        Z = varargin{1};            tmp = varargin{2};      cdo = varargin{4};
        handles.head = tmp.head;    X = tmp.X;  Y = tmp.Y;  clear tmp;
        win_name = 'Interferogram';
    end
end

% The following IF cases deal only with cases where a grid was given in argument
if (grd_data_in || grd_data_interfero_in || grd_data_deform_in)
    handles.image_type = 1;    handles.computed_grid = 1;      % Signal that this is a computed grid
    if (grd_data_interfero_in)      % Interferogram grid
        load([handles.path_data 'gmt_other_palettes.mat'],'circular');
        pal = circular;
        zz = uint8(abs(rem(double(Z),cdo)/cdo)*255);
    else                            % Not an interferogram grid
        zz = scaleto8(Z);       pal = jet(256);
    end
    set(handles.figure1,'Colormap',pal)
    aux_funs('StoreZ',handles,X,Y,Z)     % If grid size is not to big we'll store it
    set(hObject, 'Units', 'pixels');    setappdata(hObject,'Zmin_max',handles.head(5:6))
    aux_funs('colormap_bg',handles,Z,pal);
    handles = show_image(handles,win_name,X,Y,zz,1,'xy',handles.head(7));
end

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

guidata(hObject, handles);  limpa(handles);
%setappdata(handles.axes1,'ProjWKT',geogWKT) % The Geog WGS84 string in WKT format
if (~isempty(drv)),     gateLoadFile(handles,drv,varargin{1});  end
set(handles.ctrLine, 'Accelerator','l');

% --------------------------------------------------------------------------------------------------
function erro = gateLoadFile(handles,drv,fname)
    % Gateway function to load a recognized file type using its name
    erro = 0;
    switch drv
        case 'gmt',         read_DEMs(handles,{[],fname},'GMT_relatives')
        case 'generic',     FileOpenNewImage_CB([], [], handles, fname);
        case 'geotif',      FileOpenGeoTIFF_CB([], [], handles, 'nikles', fname)
        case 'multiband',   FileOpenGDALmultiBand_CB([], [], handles, 'AVHRR', fname)
        case 'envherd',     FileOpen_ENVI_Erdas_CB([], [], handles, [], fname)
        otherwise
            warndlg(['Sorry but couldn''t figure out what to do with the ' fname ' file'],'Warning')
            erro = 1;
    end

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
cmenu_axes = uicontextmenu('Parent',Parent);    % Need to set Parent for when Mirone is called from M_GMT
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
    uimenu(cmenu_axes, 'Label', 'Pixel mode on/off', 'Callback', {@PixMode_CB,handles.figure1},'Separator','on');
end

h_warning = findobj('Type','figure','Name','Warning');
if ~isempty(h_warning),     figure(h_warning);   end     % If a warning message exists, bring it forward

% --------------------------------------------------------------------
function PixMode_CB(hObject, eventdata, hFig)
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
% NAO TA ACABADA. FALTA PROGRAMAR OS CASE ILUMINACAO == 5,6   (PT)

if (handles.no_file),       return;      end
set(handles.figure1,'pointer','watch')
first_nans = 0;     pal = [];       mask = [];      crop_pol = 0;     % Defaults to croping from a rectangle
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
        [I,r_c] = cropimg(handles.head(1:2),handles.head(3:4),Z_rect,rect_crop,'out_grid');
        [m,n] = size(I);
    else                    % Extract the sub-grid inside the rectangle/polygon
        [X,Y,Z,head] = load_grd(handles);
        if isempty(Z),  set(handles.figure1,'pointer','arrow');    return;     end;    % An error message was already issued
        [Z_rect,r_c] = cropimg(head(1:2),head(3:4),Z,rect_crop,'out_grid');
        if (crop_pol)
            zzz = grdutils(Z_rect,'-L');  z_min = zzz(1);     clear zzz;
            if (strcmp(opt2,'CropaGrid_pure'))
                defAns = {sprintf('%.4f',z_min)};
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
                handles.Z_back = Z(r_c(1):r_c(2),r_c(3):r_c(4));    handles.r_c = r_c;
                Z(r_c(1):r_c(2),r_c(3):r_c(4)) = Z_rect;
                if (isnan(str2double(resp))),  handles.have_nans = 1;  first_nans = 1;  end
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
    if (strcmp(opt2,'CropaGrid'))   % Arrive here when called by "Grid Tools -> Crop Grid"
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

if (isempty(opt2) || strcmp(opt2,'CropaWithCoords'))   % Just pure Image croping
    if (m < 2 || n < 2),  set(handles.figure1,'pointer','arrow');    return;     end;    % Image too small. Probably a user bad mouse control
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
        zzz = grdutils(Z_rect,'-L');                    head(5) = zzz(1);       head(6) = zzz(2);
        to_func.Z = Z_rect;                             to_func.head = head;
    end
    if (strcmp(curr_opt,'pure'))            % PURE means pure CropaGrid
        X = (head(1) + (r_c(3)-1)*head(8)):head(8):(head(1) + (r_c(4)-1)*head(8));
        Y = (head(3) + (r_c(1)-1)*head(9)):head(9):(head(3) + (r_c(2)-1)*head(9));
        head(1) = X(1);     head(2) = X(end);       head(3) = Y(1);     head(4) = Y(end);
        tit = 'Grid cuted by Mirone';       % Have to change this to reflect the old title
        GRD_save_or_display(handles,X,Y,Z_rect,head,tit,'Croped grid')
    elseif (strcmp(curr_opt,'histo'))       % HISTO means compute histogram inside the selected rect area
        GridToolsHistogram_CB([], [], guidata(handles.figure1), to_func);
    elseif (strcmp(curr_opt,'power'))       % POWER means compute log10 power spectrum
        GridToolsSectrum_CB([], [], guidata(handles.figure1), 'Power', to_func)
    elseif (strcmp(curr_opt,'autocorr'))    % AUTOCORR means compute the autocorrelation
        GridToolsSectrum_CB([], [], guidata(handles.figure1), 'Autocorr', to_func)
    elseif (strcmp(curr_opt,'fftTools'))    % FFTTOOLS means call the fft_stuff
        GridToolsSectrum_CB([], [], guidata(handles.figure1), 'Allopts', to_func)
    end
    return
elseif (strcmp(opt2,'FillGaps'))
    if ~any(isnan(Z_rect(:)))    % No gaps
        set(handles.figure1,'pointer','arrow');    warndlg('Selected area does not have any voids (NaNs)','Warning');   return;
    else
        X = (head(1) + (r_c(3)-1)*head(8)):head(8):(head(1) + (r_c(4)-1)*head(8));
        Y = (head(3) + (r_c(1)-1)*head(9)):head(9):(head(3) + (r_c(2)-1)*head(9));
        if (~isempty(opt3) && strcmp(opt3,'surface'))
            opt_R = ['-R' sprintf('%.10f',X(1)) '/' sprintf('%.10f',X(end)) '/' sprintf('%.10f',Y(1)) '/' sprintf('%.10f',Y(end))];
            opt_I = ['-I' sprintf('%.10f',head(8)) '/' sprintf('%.10f',head(9))];
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
    defAns = {sprintf('%.12f',p_guess{1})};     resp  = inputdlg(prompt,dlg_title,[1 38],defAns);
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
    handles.Z_back = Z(r_c(1):r_c(2),r_c(3):r_c(4));        handles.r_c = r_c;
    Z(r_c(1):r_c(2),r_c(3):r_c(4)) = Z_rect;
    if (~handles.have_nans && isnan(str2double(resp)))      % See if we have new NaNs
        handles.have_nans = 1;      first_nans = 1;
    elseif (handles.have_nans && ~isnan(str2double(resp)))  % Check that old NaNs had not been erased
        handles.have_nans = grdutils(Z_rect,'-N');
    end
end

if ~isempty(opt2)       % Here we have to update the image in the processed region
    X = (head(1) + (r_c(3)-1)*head(8)):head(8):(head(1) + (r_c(4)-1)*head(8));
    Y = (head(3) + (r_c(1)-1)*head(9)):head(9):(head(3) + (r_c(2)-1)*head(9));
    [zzz] = grdutils(Z,'-L');       z_min = zzz(1);     z_max = zzz(2);     clear zzz;      img = [];
    if ( (abs(z_min - head(5)) > 1e-5 || abs(z_max - head(6)) > 1e-5) && handles.Illumin_type == 0 )
        img = scaleto8(Z);              % Z_MIN or Z_MAX have changed. Need to recompute image (but only if no illumin)
    end
    if (first_nans)       % We have NaNs for the first time. Adjust the colormap
        aux_funs('colormap_bg',handles,Z,get(handles.figure1,'Colormap'));
    end
    z_int = uint8(round( ((double(Z_rect) - z_min) / (z_max - z_min))*255 ));
    if ( handles.Illumin_type == 0)     % Nothing to do in particular
    elseif ( handles.Illumin_type >= 1 && handles.Illumin_type <= 4 )
        illumComm = getappdata(handles.figure1,'illumComm');
        z_int = ind2rgb8(z_int,get(handles.figure1,'Colormap'));    % z_int is now 3D
        head_tmp = [X(1) X(end) Y(1) Y(end) head(5:9)];
        if (handles.Illumin_type == 1)
            opt_N = ['-Nt1/' sprintf('%.4f',handles.grad_sigma) '/' sprintf('%.4f',handles.grad_offset)];
            if (handles.geog),  R = grdgradient_m(Z_rect,head_tmp,'-M',illumComm,opt_N);
            else                R = grdgradient_m(Z_rect,head_tmp,illumComm,opt_N); end
        else
            R = grdgradient_m(Z_rect,head_tmp,illumComm);
        end
        z_int = shading_mat(z_int,R,'no_scale');    % and now it is illuminated
    else
        warndlg('Sorry, this operation is not allowed with this shading illumination type','Warning')
        set(handles.figure1,'pointer','arrow');     return
    end
    if (isempty(img))  img = get(handles.hImg,'CData');      end     % If img was not recomputed, get from screen 
    handles.img_back = img(r_c(1):r_c(2),r_c(3):r_c(4),:);      % For the undo op
    if (~isempty(mask) && handles.Illumin_type ~= 0)
        mask = repmat(~mask,[1 1 3]);
        z_int(mask) = handles.img_back(mask);
    end
    img(r_c(1):r_c(2),r_c(3):r_c(4),:) = z_int;         clear z_int Z_rect R;
    set(handles.hImg,'CData',img)

    head(5) = z_min;     head(6) = z_max;
    handles.computed_grid = 1;              handles.head = head;    handles.origFig = img;
    setappdata(handles.figure1,'dem_z',Z);  setappdata(handles.figure1,'GMThead',head);
    setappdata(handles.figure1,'Zmin_max',[head(5) head(6)])
end
guidata(handles.figure1, handles);          set(handles.figure1,'pointer','arrow')

% Experimental UNDO
if strmatch(opt2,{'MedianFilter' 'ROI_MedianFilter' 'SetConst' 'ROI_SetConst'}),      set_undo(handles,opt);  end

% --------------------------------------------------------------------
function set_undo(handles,h)
	% Experimental UNDO  that works only with the "Median Filter" option
	cmenuHand = get(h,'UIContextMenu');     cb_undo = {@do_undo,handles,h,cmenuHand};
	uimenu(cmenuHand, 'Label', 'Undo', 'Separator','on', 'Callback', cb_undo);

% -----------------------------------------------------------------------------------------
function do_undo(obj,eventdata,handles,h,img)
	[X,Y,Z,head] = load_grd(handles);   % Experimental. No testing for error in loading
	Z(handles.r_c(1):handles.r_c(2),handles.r_c(3):handles.r_c(4)) = handles.Z_back;
	[zzz] = grdutils(Z,'-L');       z_min = zzz(1);     z_max = zzz(2);     clear zzz;
	head(5) = z_min;                head(6) = z_max;
	setappdata(handles.figure1,'dem_z',Z);    setappdata(handles.figure1,'GMThead',head);
	handles.origFig = get(handles.hImg,'CData');
	handles.origFig(handles.r_c(1):handles.r_c(2),handles.r_c(3):handles.r_c(4),:) = handles.img_back;
	set(handles.hImg,'CData',handles.origFig)
	guidata(handles.figure1,handles)
	ui_u = findobj(get(h,'UIContextMenu'),'Label', 'Undo');
	set(ui_u,'Visible','off')       % If I delete it we get an error

% --------------------------------------------------------------------
function ImageFlip_CB(hObject, eventdata, handles, opt)
	% OPT == 'LR'   Flips the image left-right
	% OPT == 'UD'   Flips the image up-down
	if (handles.no_file == 1),    return;      end
	img = get(handles.hImg,'CData');
	if strcmp(opt,'LR'),    img = flipdim(img,2);
	else                    img = flipdim(img,1);
	end
	set(handles.hImg,'CData', img);

% --------------------------------------------------------------------
function ImageResetOrigImg_CB(hObject, eventdata, handles)
	if (handles.no_file == 1),    return;      end
	try         % In some cases (e.g. histograms, countries) img may not exist
		set(handles.hImg,'CData', handles.origFig);
        set(handles.figure1,'ColorMap',handles.origCmap)
		handles.hist = 0;   handles.Illumin_type = 0;   handles.firstIllum = 1;
        handles.ValidGrid = handles.ValidGrid_orig;     handles.was_int16 = handles.was_int16_orig;
        handles.computed_grid = handles.computed_grid_orig;
        set(handles.ImgHist,'checked','off');
		guidata(hObject, handles);
	end

% --------------------------------------------------------------------
function ImageHistEqualize_CB(hObject, eventdata, handles)
	if (handles.no_file == 1),     return;      end
	
	zz = get(handles.hImg,'CData');
	if strcmp(get(hObject,'checked'),'off')     % Then equalize
        if (ndims(zz) == 3);
            zz = cvlib_mex('color',zz,'rgb2hsv');
            zz(:,:,3) = img_fun('histeq_j',zz(:,:,3));
            J = cvlib_mex('color',zz,'hsv2rgb');
        else
            J = img_fun('histeq_j',zz);
        end
        set(handles.hImg,'CData', J);    handles.hist = 1;       set(hObject,'checked','on')
	elseif ~isempty(handles.origFig)            % Then de-equalize
        set(handles.hImg,'CData', handles.origFig);    handles.hist = 0;
        handles.Illumin_type = 0;       set(hObject,'checked','off');
	else
        msgbox('Sorry. To save memory I didn''t make an image copy. You''ll have to reload the file again.','Warning')
	end
	guidata(hObject, handles);

% --------------------------------------------------------------------
function ImageHistEqualizeGrid_CB(hObject, eventdata, handles)
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
function ImageShowPalette_CB(hObject, eventdata, handles)
	if (handles.no_file == 1),     return;      end
	if (ndims(get(handles.hImg,'CData')) == 3)
        msgbox('True color images do not use color palettes.','Warning');    return
	end
	cmap = get(handles.figure1,'Colormap');     dz = inf;
	if (handles.ValidGrid),     dz = diff(handles.head(6:7)) / 10;      end
    show_palette([1 length(cmap)], [handles.head(5) handles.head(6)], dz, cmap)

% --------------------------------------------------------------------
function PanZoom_CB(hObject,event,handles,opt)
if (handles.no_file == 1),    set(hObject,'State','off');   return;      end
if (strcmp(get(handles.Tesoura,'State'),'on'))  % If Scisors were on
    set(handles.Tesoura,'State','off')
end
if (strcmp(opt,'zoom'))
	if strcmp(get(hObject,'State'),'on')
        zoom_j('on');
        if (strcmp(get(handles.Mao,'State'),'on'))
            set(handles.Mao,'State','off');   pan('off');
        end
	else
        zoom_j('off');
	end
else        % Pan case
	if strcmp(get(hObject,'State'),'on')
        pan('on');
        if (strcmp(get(handles.Zoom,'State'),'on'))
            set(handles.Zoom,'State','off');   zoom_j('off');
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
function FilePreferences_CB(hObject, eventdata, handles)
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
    handles.flederBurn = out.flederBurn;
    handles.flederPlanar = out.flederPlanar;
	setappdata(handles.figure1,'swathRatio',handles.swathRatio);    % I need this in getline_mb
	setappdata(handles.figure1,'DefLineThick',handles.DefLineThick)     % Save this for accessability in draw_funs
	setappdata(handles.figure1,'DefLineColor',handles.DefLineColor)     % Save this for accessability in draw_funs
	guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function FileNewBgFrame_CB(hObject, eventdata, handles, region, imSize)
% Create a empty window with a frame selected in bg_region
% However, if REGION was transmited, it is assumed to have  [x_min x_max y_min y_max is_geog]
if (nargin == 3)
    region = bg_region;     % region contains [x_min x_max y_min y_max is_geog]
    if isempty(region),     return;     end     % User gave up
end

X = region(1:2);      Y = region(3:4);
handles.head = [X Y 0 255 0];

scrsz = get(0,'ScreenSize');         % Get screen size
aspect = diff(Y) / diff(X);
nx = round(scrsz(3)*.75);   ny = round(nx * aspect);
if (ny > scrsz(4) - 30)
    ny = scrsz(4) - 30;     nx = round(ny / aspect);
end
handles.head(8) = diff(X) / (nx - 1);    handles.head(9) = diff(Y) / (ny - 1);
Z = repmat(uint8(255),ny,nx);
pal = repmat(handles.bg_color,256,1);    set(handles.figure1,'Colormap',pal);
handles.image_type = 20;
if (nargin <= 4),   imSize = [];    end
show_image(handles,'Mirone Base Map',X,Y,Z,0,'xy',0,imSize);

% --------------------------------------------------------------------
function FileNewBgMap_CB(hObject, eventdata, handles)
	% The output of bg_regionMapTilled is a structure with the following fields:
	% out.img   -> the image                % out.X         -> image's X limits
	% out.Y     -> image's Y limits         % out.imgName   -> image's full path and name
	out = bg_region_map_tilled;
	if isempty(out),    return;     end     % User gave up loading the fig tille
	handles.geog = 1;                       handles.image_type = 3;
	handles.head(1:2) = out.X;              handles.head(3:4) = out.Y;      handles.head(5:7) = [0 255 0];
	handles.head(8) = diff(out.X) / (size(out.img,2)-1);    handles.head(9) = diff(out.Y) / (size(out.img,1)-1);
	handles = show_image(handles,out.imgName,out.X,out.Y,out.img,0,'xy',0,1);
	guidata(hObject,handles);

% --------------------------------------------------------------------
function FileNewEmpty_CB(hObject, eventdata, handles)
	setappdata(0,'parent_gcf',handles.figure1);
	h = mirone;
	set(findobj(h,'Tag','ImageDrape'),'Enable','on')        % Set the Drape option to 'on' in the New window 
    
% --------------------------------------------------------------------
function FileSaveGMTgrid_CB(hObject, eventdata, handles, opt)
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
    cdf_format = '=nf';
    if (handles.was_int16 && handles.saveAsInt16),  cdf_format = '=ns';    end
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
function File_img2GMT_RGBgrids_CB(hObject, eventdata, handles, opt1, opt2)
% Save image as a triplet of gmt grids - R,G,B.
% OPT1 == 'image' || == [] Capture only the image and not graphical elements
%   (lines, symbols, etc...). The grids nrow & ncol is the same as the image
%   number of lines and pixels.
% OPT1 == 'screen' Does a screen capture that includes all the graphical elements
%   that may have been drawn (lines, symbols, etc...). On the other hand I still
%   don't know how to control the number of lines and pixels.
% OPT2 == fname. It is used by the write gmt script routine to capture the image and
%   write the image as a triplet of gmt grids with name stem = OPT2.
if (handles.no_file == 1),      aux_funs('msg_dlg',1,handles);     return;      end

if    (nargin == 3),    opt1 = 'image';     opt2 = [];
elseif(nargin == 4),    opt2 = [];          end

str1 = {'*.grd;*.GRD','netCDF int2 grid format (*.grd,*.GRD)'; '*.*', 'All Files (*.*)'};
if (isempty(opt2))
    [FileName,PathName] = put_or_get_file(handles,str1,'Select output GMT grid','put');
    if isequal(FileName,0);     return;     end
else        % It means the output file name was transmited in input
    [PathName,FileName] = fileparts(opt2);
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

if (handles.image_type == 1)                        % GMT grid
    tit = 'GMT grid converted to 3 RGB grids';
elseif (handles.image_type == 2)                    % Generic formats
    tit = get(handles.figure1,'Name');
end
if (length(tit) > 80),  tit = tit(1:80);    end     % (1:80) otherwise it BOOMs

if (strcmp(opt1,'image')),  img = get(handles.hImg,'CData');            % Get image
else                        img = snapshot(handles.figure1,'noname');   % Screen capture
end

if (ndims(img) == 2)
    img = ind2rgb8(img,get(handles.figure1,'Colormap'));    % Need this because image is indexed
end
if (~strcmp(get(handles.axes1,'Ydir'),'normal')),    img = flipdim(img,1);      end

grdwrite_m(single(img(:,:,1)),D,f_name_r,tit)
grdwrite_m(single(img(:,:,2)),D,f_name_g,tit)
grdwrite_m(single(img(:,:,3)),D,f_name_b,tit)
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function FileSaveENCOMgrid_CB(hObject, eventdata, handles)
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
function FilePrint_CB(hObject, eventdata, handles)
	if (handles.no_file == 1),    return;      end
	h = findobj('Type','uicontrol');
	set(h,'Visible','off')              % We don't want to print the buttons
	handsStBar = getappdata(handles.figure1,'CoordsStBar');
	set(handsStBar,'Visible','off');
	if (ispc),      print -v
	else            print;  end
	set(h,'Visible','on');  set(handsStBar(2:end),'Visible','on');

% --------------------------------------------------------------------
function ExtractProfile_CB(hObject, eventdata, handles, opt)
if (handles.no_file == 1),     return;      end
point_int = 0;                                  % Default to "profile" interpolation
if (nargin == 4),   point_int = 1;     end      % Interpolate at the line vertex only
[X,Y,Z,head] = load_grd(handles,'silent');
if (isempty(Z) && (handles.image_type == 1 || handles.image_type == 4))     % Grid not in memory error
    errordlg('Grid was not on memory. Increase "Grid max size" and start over again.','ERROR'); return
elseif (isempty(Z) && ndims(get(handles.hImg,'CData')) == 2)
    Z = get(handles.hImg,'CData');
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
        rp = aux_funs('axes2pix',rows, get(handles.hImg,'YData'),yy);
        cp = aux_funs('axes2pix',cols, get(handles.hImg,'XData'),xx);
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
function FileOpenBSB_CB(hObject, eventdata, handles)
str1 = {'*.kap;*.KAP;*.nos;*NOS', 'BSB Nautical Chart (*.kap,*.KAP,*.nos,*.NOS)'; '*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select BSB Nautical Chart File','get');
if isequal(FileName,0);     return;     end

set(handles.figure1,'pointer','watch');     handles.fileName = [PathName FileName];
att = gdalread(handles.fileName,'-M','-C');
if (~isempty(att.ProjectionRef)),   Z = gdalread(handles.fileName,'-U','-C');
else                                Z = gdalread(handles.fileName,'-C');
end
handles.head = att.GMT_hdr;
X = handles.head(1:2);      Y = handles.head(3:4);
if (~isempty(att.Band(1).ColorMap)),    pal = att.Band(1).ColorMap.CMap(:,1:3);
else                                    pal = jet(256);     end

handles.image_type = 3;
set(handles.figure1,'Colormap',pal);
aux_funs('cleanGRDappdata',handles);        % Remove eventual grid stuff variables from appdata
show_image(handles,handles.fileName,X,Y,Z,0,'xy',1);
grid_info(handles,att,'gdal')           % Construct a info message

% --------------------------------------------------------------------
function FileOpen_ENVI_Erdas_CB(hObject, eventdata, handles, opt, opt2)
% This function reads both ENVI or Erdas files. Furthermore, based on the file byte
% type it guesses if we are dealing with a typical grid file (in which case it is
% treated like a native gmt grid) or a raster image file.

if (nargin == 4)    % Otherwise, OPT2 already contains the File name
    str1 = {'*.img;*.IMG', [opt ' (*.img,*.IMG)']; '*.*', 'All Files (*.*)'};
    [FileName,PathName] = put_or_get_file(handles,str1,['Select ' opt ' File'],'get');
    if isequal(FileName,0);     return;     end
else
    PathName = [];      FileName = opt2;
end

if (handles.ForceInsitu),   opt_I = '-I';   % Use only in desperate cases.
else                        opt_I = ' ';    end
handles.was_int16 = 0;      % To make sure that it wasnt left = 1 from a previous use.
handles.fileName = [PathName FileName];

att = gdalread(handles.fileName,'-M','-C');
if ((att.RasterXSize * att.RasterYSize * 4) > handles.grdMaxSize)
    if ( strcmp(yes_or_no('title','Warning'),'Yes')),  return;     end      % Advise accepted
end

set(handles.figure1,'pointer','watch')
if (strcmp(att.Band(1).DataType,'Byte') || ~isempty(att.Band(1).ColorMap))     % We have a raster image
    if (strcmp(att.Band(1).DataType,'Byte'))
        zz = gdalread(handles.fileName,'-U');
    else
        zz = gdalread(handles.fileName,'-S','-U');
    end
    head = att.GMT_hdr;             [m,n,k] = size(zz);
    X = [head(1) head(2)];          Y = [head(3) head(4)];
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
    Z = gdalread(handles.fileName,'-U','-C',opt_I);       Z = single(Z);
    head = att.GMT_hdr;
    if (~isempty(att.Band(1).NoDataValue)),  Z(Z <= single(att.Band(1).NoDataValue)) = NaN;    end
    handles.have_nans = grdutils(Z,'-N');    
    zz = scaleto8(Z);       [m,n] = size(Z);
    X = linspace(head(1),head(2),n);  Y = linspace(head(3),head(4),m);  % Need this for image

    aux_funs('StoreZ',handles,X,Y,Z)            % If grid size is not to big we'll store it
    handles.image_type = 4;
    ValidGrid = 1;                              % Signal that grid opps are allowed
    setappdata(handles.figure1,'Zmin_max',[head(5) head(6)])
    aux_funs('colormap_bg',handles,Z,jet(256));
end

handles.head = head;
show_image(handles,handles.fileName,X,Y,zz,ValidGrid,'xy',1);
grid_info(handles,att,'gdal')           % Construct a info message

% --------------------------------------------------------------------
function FileOpenNewImage_CB(hObject, eventdata, handles, opt)
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
if (nargin == 3)
    [FileName,PathName,handles] = put_or_get_file(handles,str1,'Select image format','get');
    if isequal(FileName,0);     return;     end
else                % Filename was transmited in input
    PathName = [];      FileName = opt;
end
handles.fileName = [PathName FileName];

set(handles.figure1,'pointer','watch')
[PATH,FNAME,EXT] = fileparts(handles.fileName);
if (strcmpi(EXT,'.shade'))
    [fid, msg] = fopen(handles.fileName, 'r');
    if (fid < 0)
        set(handles.figure1,'pointer','arrow');    errordlg([handles.fileName ': ' msg],'ERROR');   return
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
    FileOpenGDALmultiBand_CB(hObject, [], handles, 'RAW', handles.fileName)
    return      % We are done here. Bye Bye.
else
    info_img = imfinfo(handles.fileName);
    if ( any(strcmpi(EXT,{'.tif' '.tiff'})) && strcmp(info_img.Compression,'LZW'))
        % If Tiffs are LZW, imread R13 is not able to read them 
        I = gdalread(handles.fileName);
    elseif (strcmpi(EXT,'.gif') && handles.IamCompiled)   % Try with GDAL because, ML uses f... java
        [I,att] = gdalread(handles.fileName);
        info_img.ColorType = att.Band(1).ColorMap;
    else
        try,        [I,map] = imread(handles.fileName);
        catch       errordlg(lasterr,'Error');      return % It realy may happen
        end
    end
    if (strcmp(info_img.ColorType,'grayscale'))
        set(handles.figure1,'Colormap',gray(256))
    elseif (isfield(info_img,'ColorTable'))         % Gif images call it 'ColorTable'
        set(handles.figure1,'Colormap',info_img.ColorTable)
    elseif (isfield(info_img,'Colormap') && ~isempty(info_img.Colormap))
        set(handles.figure1,'Colormap',info_img.Colormap)
    end
end

handles.image_type = 2;     handles.geog = 0;   % None of this image types is coordinated
handles.head = [1 size(I,2) 1 size(I,1) 0 255 0 1 1];   % Fake a grid reg GMT header
aux_funs('cleanGRDappdata',handles);            % Remove eventual grid stuff variables from appdata
handles = show_image(handles,handles.fileName,[],[],I,0,'off',0);
grid_info(handles,handles.fileName,'iminfo');    % Construct a info string
guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function FileOpenGDALmultiBand_CB(hObject, eventdata, handles, opt, opt2)
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
    opt_B = sprintf('%s%d','-B1-',bands_inMemory);
    [I,att] = gdalread(fname,'-S', opt_B,'-C');
    n_bands = att.RasterCount;
    bands_inMemory = 1:min(n_bands,bands_inMemory);      % Make it a vector
    handles.head = att.GMT_hdr;
    if (~isempty(att.GCPvalues))    % Save GCPs so that we can plot them and warp the image
        setappdata(handles.figure1,'GCPregImage',att.GCPvalues)
        setappdata(handles.figure1,'fnameGCP',fname)    % Save this to know when GCPs are to be removed
    end                                                 % from appdata. That is donne in show_image()
end

tmp1 = cell(n_bands+1,2);    tmp2 = cell(n_bands+1,2);
tmp1{1,1} = opt;    tmp1{1,2} = opt;
for (i = 1:n_bands)
    tmp1{i+1,1} = ['band' sprintf('%d',i)];
    tmp1{i+1,2} = ['banda' sprintf('%d',i)];    % TEMP
    tmp2{i+1,1} = [sprintf('%d',i) 'x1 bip'];        % TEMP
    tmp2{i+1,2} = i;
end
tmp = {['+ ' opt]; I; tmp1; tmp2; fname; bands_inMemory; [att.RasterYSize att.RasterXSize n_bands]; reader};
if (n_bands > 3),       I(:,:,4:end) = [];      % Now I can only be MxN or MxNx3
elseif (n_bands == 2)   I(:,:,2) = [];
end

setappdata(handles.figure1,'BandList',tmp)
handles.image_type = 2;     handles.fileName = [];  % No eligible for automatic re-loading
handles.was_int16 = 0;      handles.computed_grid = 0;
handles.geog = 0;           % None of this image types is coordinated (nor geog nor anything else)
aux_funs('cleanGRDappdata',handles);        % Remove eventual grid stuff variables from appdata
if (~strcmp(opt,'RAW'))
    if (isempty(att.Band(1).ColorMap)),     set(handles.figure1,'Colormap',jet(256))
    else                                    set(handles.figure1,'Colormap',att.Band(1).ColorMap.CMap)
    end
end
if (n_bands == 1),      set(handles.figure1,'Colormap',gray(256));   end    % Takes precedence over the above
show_image(handles,fname,[],[],I,0,'off',0);    % It also guidata(...) & reset pointer
if (isappdata(handles.axes1,'InfoMsg')),    rmappdata(handles.axes1,'InfoMsg');     end
if (~strcmp(opt,'RAW')),    grid_info(handles,att,'gdal');      end     % Construct a info message

% --------------------------------------------------------------------
function FileOpenGeoTIFF_CB(hObject, eventdata, handles, tipo, opt)
switch lower(tipo)
    case 'geotiff',     str1 = {'*.tif;*.TIF;*.tiff;*.TIFF', 'GeoTiff (*.tif,*.tiff,*.TIF,*.TIFF)'};
    case 'sid',         str1 = {'*.sid;*.SID', 'GeoTiff (*.sid,*.SID)'};
    case 'ecw',         str1 = {'*.ecw;*.ECW', 'GeoTiff (*.ecw,*.ECW)'};
    case 'jp2',         str1 = {'*.jp2;*.JP2', 'GeoTiff (*.jp2,*.JP2)'};
    otherwise,          str1 = {'', 'Don''t know (*.*)'};       % Used with the "Try Luck" option
end
if (nargin == 4)
    str1(2,1:2) = {'*.*', 'All Files (*.*)'};
    [FileName,PathName] = put_or_get_file(handles,str1,['Select ' tipo ' file'],'get');
    if isequal(FileName,0);     return;     end
else                % Filename was transmited in input
    PathName = [];      FileName = opt;
end
handles.fileName = [PathName FileName];

att = gdalread(handles.fileName,'-M','-C');
set(handles.figure1,'pointer','watch')

if (att.RasterCount > 3)            % Since it is a multiband file, try luck there
    FileOpenGDALmultiBand_CB([],[],handles,'AVHRR',handles.fileName);
    return
end

if ~(strcmp(att.Band(1).DataType,'Byte'))       % JPK2, for example, may contain DTMs
    fullname{1} = PathName;     fullname{2} = FileName;
    read_DEMs(handles,fullname,'JP2_DEM');    return
end

Z = gdalread(handles.fileName,'-U');
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

set(handles.figure1,'Colormap',pal);
aux_funs('cleanGRDappdata',handles);        % Remove eventual grid stuff variables from appdata
show_image(handles,handles.fileName,X,Y,Z,0,ax_dir,0);
grid_info(handles,att,'gdal')               % Construct a info message

% --------------------------------------------------------------------
function FileOpenDEM_CB(hObject, eventdata, handles, opt)
% Files of the following formats are read (well re-directed) here
tipo = opt;
switch opt
    case {'GMT' 'Surfer', 'ENCOM'}
        str1 = {'*.grd;*.GRD;*.nc;*.NC', 'Grid files (*.grd,*.GRD,*.nc,*.NC)';'*.*', 'All Files (*.*)'};    tipo = 'GMT_relatives';
    case 'MANI'
        str1 = {'*.man;*.MAN', 'Grid files (*.man,*.MAN)';'*.*', 'All Files (*.*)'};    tipo = 'GMT_relatives';
    case 'ArcAscii',        str1 = {'*.grd;*.GRD', 'Arc/Info grid (*.grd,*.GRD)'; '*.*', 'All Files (*.*)'};
    case 'ArcBinary',       str1 = {'*.adf;*.ADF', 'Arc/Info grid (*.adf,*.ADF)'; '*.*', 'All Files (*.*)'};
    case 'DTED',            str1 = {'*.dt0;*.DT0;*.dt1;*.DT1', 'DTED (*.dt0,*.DT0,*.dt1,*.DT1)'; '*.*', 'All Files (*.*)'};
    case 'ESRI_hdr',	    str1 = {'*.bil;*.BIL;', 'ESRI BIL (*.bil,*.BIL)'; '*.*', 'All Files (*.*)'};
    case 'GTOPO30',         str1 = {'*.dem;*.DEM', 'GTOPO30 DEM (*.dem,*.DEM)'; '*.*', 'All Files (*.*)'};
    case 'GeoTiff_DEM',     str1 = {'*.tif;*.TIF;*.tiff;*.TIFF', 'GeoTiff DEM(*.tif,*.tiff,*.TIF,*.TIFF)'; '*.*', 'All Files (*.*)'};
    case 'GXF',             str1 = {'*.gxf;*.GXF', 'Geosoft GXF (*.gxf,*.GXF)'; '*.*', 'All Files (*.*)'};
    case 'SDTS',            str1 = {'*catd.ddf;*CATD.DDF', 'USGS SDTS DEM (*catd.ddf,*CATD.DDF)'; '*.*', 'All Files (*.*)'};
    case 'SRTM30',    	    str1 = {'*.srtm;*.SRTM;*.srtm.gz', 'SRTM30 DEM (*.srtm,*.SRTM,*.srtm.gz)'; '*.*', 'All Files (*.*)'};
    case {'SRTM1' 'SRTM3'}, str1 = {'*.hgt;*.HGT;*.hgt.zip', [opt ' DEM (*.hgt,*.HGT,*.hgt.zip)']; '*.*', 'All Files (*.*)'};
    case 'USGS_DEM',        str1 = {'*.dem;*.DEM', 'USGS DEM (*.dem,*.DEM)'; '*.*', 'All Files (*.*)'};
    otherwise
        errordlg(['OOPs, where did this ' opt ' code came in?'],'Error');   return
end
[FileName,PathName] = put_or_get_file(handles,str1,['Select ' opt ' File'],'get');
if isequal(FileName,0);     return;     end
fullname{1} = PathName;     fullname{2} = FileName;
read_DEMs(handles,fullname,tipo)

% --------------------------------------------------------------------
function FileOpenMOLA_CB(hObject, eventdata, handles)
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
handles.fileName = [fullname{1} fullname{2}];       % To store any input grid/image file name
if (handles.ForceInsitu),   opt_I = '-I';           % Use only in desperate cases.
else                        opt_I = ' ';    end
handles.was_int16 = 0;      % To make sure that it wasnt left = 1 from a previous use.
if (~isempty(getappdata(handles.figure1,'dem_z')))   % Free memory that may be necessary bellow
    rmappdata(handles.figure1,'dem_z')
end

if (strcmp(tipo,'GMT_relatives'))   % The reading is done by the read_gmt_type_grids function
    [handles,X,Y,Z,head] = read_gmt_type_grids(handles,handles.fileName);
    if (isempty(X)),    set(handles.figure1,'pointer','arrow');     return;     end
    zz = scaleto8(Z);
elseif ( strcmp(tipo,'SRTM30') || strcmp(tipo,'SRTM3') || strcmp(tipo,'SRTM1') || strcmp(tipo,'MOLA') )
    name_uncomp = [];    tipo1 = tipo;
    if (strcmp(tipo,'MOLA'))    % Must inform write_ESRI_hdr what to write
        tipo1 = [opt(1) opt(4) opt(5) opt(6) opt(7) opt(7) -99999];
    end
    [name_hdr,comp_type] = write_ESRI_hdr(handles.fileName,tipo1);
    if (~isempty(comp_type))
        fname = decompress(handles.fileName,'warn');    % Named with compression ext removed
        name_uncomp = fname;        % Here we need a copy of the decompressed file name for removing
        if (isempty(fname)),  return;     end;          % Error message already issued.
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
    att = gdalread(handles.fileName,'-M','-C');
    if ((att.RasterXSize * att.RasterYSize * 4) > handles.grdMaxSize)
        if ( strcmp(yes_or_no('title','Warning'),'Yes')),  return;     end      % Advise accepted
    end
    
    Z =  gdalread(handles.fileName,'-U',opt_I,'-C');     Z = single(Z);
    head = att.GMT_hdr;
    if (~isempty(att.Band(1).NoDataValue)),  Z(Z <= single(att.Band(1).NoDataValue)) = NaN;    end
    handles.have_nans = grdutils(Z,'-N');
    if (isequal(head(5:6),[0 0]))   % It happens with GeoTiff_DEM
        zzz = grdutils(Z,'-L');             head(5) = zzz(1);       head(6) = zzz(2);
    end
    zz = scaleto8(Z);    [m,n] = size(Z);
    X = linspace(head(1),head(2),n);  Y = linspace(head(3),head(4),m);  % Need this for image
    
    if (strcmp(att.Band(1).DataType,'Int16')), handles.was_int16 = 1;  end
    handles.image_type = 4;     handles.Nodata_int16 = att.Band(1).NoDataValue;
end

aux_funs('StoreZ',handles,X,Y,Z)    % If grid size is not to big we'll store it
handles.head = head;
setappdata(handles.figure1,'Zmin_max',[head(5) head(6)])
aux_funs('colormap_bg',handles,Z,jet(256));
show_image(handles,handles.fileName,X,Y,zz,1,'xy',head(7));
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
    handles = gcpTool(handles,axis_t,X,Y,I);
    return
end

handles.hImg = image(X,Y,I,'Parent',handles.axes1);      handles.hImg = handles.hImg;
zoom_state(handles,'off_yes');      set(handles.hImg,'CDataMapping','direct')
if (strcmp(axis_t,'xy')),           set(handles.axes1,'XDir','normal','YDir','normal')
elseif (strcmp(axis_t,'off')),      set(handles.axes1,'Visible','off')
else    warndlg('Warning: Unknown axes setting in show_image','Warning')
end
resizetrue(handles.figure1,imSize);
handles.origFig = I;            handles.no_file = 0;
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
handles.geog = aux_funs('guessGeog',handles.head(1:4));
if (handles.image_type ~= 1),   handles.grdname = [];   end
guidata(handles.figure1, handles);              set(handles.figure1,'pointer','arrow')

if (~ValidGrid)      % Hide uicontrols that are useless to images only
    set(handles.ImgHistGrd,'Visible','off');    set(handles.Illuminate,'Visible','off')
    set(handles.Contours_a,'Visible','off');    set(handles.Contours_i,'Visible','off')
    set(handles.MBplan,'Visible','off');        set(handles.GridTools,'Visible','off');
    set(handles.SaveGMTgrid,'Enable','off');
else
    set(handles.ImgHistGrd,'Visible','on');     set(handles.Illuminate,'Visible','on')
    set(handles.Contours_a,'Visible','on');     set(handles.Contours_i,'Visible','on')
    set(handles.MBplan,'Visible','on');         set(handles.GridTools,'Visible','on');
    set(handles.SaveGMTgrid,'Enable','on');
end
if (strcmp(axis_t,'off')),      set(handles.noAxes,'Visible','off');
else                            set(handles.noAxes,'Visible','on');
end

GCPmemoryVis = 'off';
if (isappdata(handles.figure1,'GCPregImage'))
    fnameGCP = getappdata(handles.figure1,'fnameGCP');
    if (~strcmp(fnameGCP,fname))        % It means the 'GCPregImage' app must be from a previous file
        rmappdata(handles.figure1,'GCPregImage')
    else
        GCPmemoryVis = 'on';
    end
end
set(handles.GCPmemory,'Visible',GCPmemoryVis)

BL = getappdata(handles.figure1,'BandList');        % We must tell between fakes and true 'BandList'
if ( isempty(BL) || ((numel(BL{end}) == 1) && strcmp(BL{end},'Mirone')) )
	if (ndims(I) == 3)          % Some cheating to allow selecting individual bands of a RGB image
        tmp1 = cell(4,2);   tmp2 = cell(4,2);    tmp1{1,1} = 'RGB';     tmp1{1,2} = 'RGB';
        for (i = 1:3)
            tmp1{i+1,1} = ['band' sprintf('%d',i)];     tmp1{i+1,2} = ['banda' sprintf('%d',i)];
            tmp2{i+1,1} = [sprintf('%d',i) sprintf('%d',size(I,1)) 'x' sprintf('%d',size(I,2)) ' BSQ']; tmp2{i+1,2} = i;
        end
        tmp = {['+ ' 'RGB']; I; tmp1; tmp2; ''; 1:3; [size(I,1) size(I,2) 3]; 'Mirone'};
        setappdata(handles.figure1,'BandList',tmp)
        set(findobj(handles.figure1,'Label','Load Bands'),'Enable','on')
	elseif (ndims(I) == 2)      % Remove it so it won't try to operate on indexed images
        if (isappdata(handles.figure1,'BandList')),     rmappdata(handles.figure1,'BandList');  end
        set(findobj(handles.figure1,'Label','Load Bands'),'Enable','off')
	end
end
if (isappdata(handles.axes1,'ProjWKT')),            rmappdata(handles.axes1,'ProjWKT');         end
if (isappdata(handles.axes1,'DatumProjInfo')),      rmappdata(handles.axes1,'DatumProjInfo');   end
% Note that, when it applyies the above are rebuilt with a latter call to grid_info(handles,att,'gdal')

% --------------------------------------------------------------------
function ToolsMBplaningStart_CB(hObject, eventdata, handles)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
if isempty(getappdata(handles.figure1,'dem_z'))     % Test if the grid is loaded in memory
    warndlg('Grid file is bigger than the declared "Grid Max Size". See "File -> Preferences"','Warning');
    return
end

prompt = {['The current value for the swath-width / water depth ratio is:   --> ' sprintf('%g',handles.swathRatio) '  <--']
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
tagL = ['MBtrack' sprintf('%d',handles.nTrack)];    tagB = ['swath_w' sprintf('%d',handles.nTrack)];
set(trackHand,'XData',xp, 'YData',yp, 'Tag',tagL);  set(barHand, 'Tag',tagB)
setappdata(trackHand,'swathRatio',handles.swathRatio)  % save the swathRatio in line's appdata
% now set the barHand userdata that contains the vertex order
for i = 1:length(xp),     set(barHand(i),'Userdata',i);     end

draw_funs(trackHand,'MBtrack_uicontext')        % Set track's uicontextmenu
draw_funs(barHand,'MBbar_uicontext')            % Set track bar's uicontextmenu

zoom_state(handles,'maybe_on');
handles.h_MBplot = trackHand;       guidata(handles.figure1, handles);
str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select output xy file name','put');
if isequal(FileName,0);     return;     end
double2ascii([PathName FileName],[xp yp],'%f\t%f');

% --------------------------------------------------------------------
function ToolsMBplaningImport_CB(hObject, eventdata, handles)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
if (~handles.ValidGrid)
    errordlg('This operation is deffined only for images derived from DEM grids','Error');  return
end
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
    tagL = ['MBtrack' sprintf('%d',handles.nTrack)];     tagB = ['swath_w' sprintf('%d',handles.nTrack)];
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
function ImageIlluminationModel_CB(hObject, eventdata, handles, opt)
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
elseif (luz.illum_model == 4)   % GMT Lambertian illumination
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
% OPT ->  Select which of the GMT grdgradient illumination algorithms to use
% Illuminate a DEM file and turn it into a RGB image
% For multiple tryies I need to use the original image. Otherwise each attempt would illuminate
% the previously illuminated image. An exception occurs when the image was IP but only for the
% first time, repeated illums will use the original img. Otherwise we would need to make another img copy

[X,Y,Z,head] = load_grd(handles);   % If needed, load gmt grid again
if isempty(Z),   return;     end;   % An error message was already issued

set(handles.figure1,'pointer','watch')

if (handles.firstIllum),    img = get(handles.hImg,'CData');     handles.firstIllum = 0;
else                        img = handles.origFig;      end
if (isempty(img)),          img = get(handles.hImg,'CData');     end             % No copy in memory
if (ndims(img) == 2),       img = ind2rgb8(img,get(handles.figure1,'Colormap'));    end    % Image is 2D   

setappdata(handles.figure1,'Xmin_max',[X(1) X(end)]);   setappdata(handles.figure1,'Ymin_max',[Y(1) Y(end)]);
setappdata(handles.figure1,'Zmin_max',[head(5) head(6)])

if (strcmp(opt,'grdgradient_class'))    % GMT grdgradient classic illumination
    illumComm = ['-A' sprintf('%.2f',luz.azim)];
    if (handles.geog),  [R,offset,sigma] = grdgradient_m(Z,head,'-M',illumComm,'-Nt');
    else                [R,offset,sigma] = grdgradient_m(Z,head,illumComm,'-Nt'); end
    handles.Illumin_type = 1;
    handles.grad_offset = offset;   handles.grad_sigma = sigma;
elseif (strcmp(opt,'grdgradient_lamb')) % GMT grdgradient lambertian illumination
    illumComm = ['-Es' sprintf('%.2f',luz.azim) '/' sprintf('%.2f',luz.elev)];
    R = grdgradient_m(Z,head,illumComm);
    handles.Illumin_type = 2;
elseif (strcmp(opt,'grdgradient_peuck'))% GMT grdgradient Peucker illumination
    illumComm = '-Ep';
    R = grdgradient_m(Z,head,illumComm);
    handles.Illumin_type = 3;
elseif (strcmp(opt,'lambertian'))       % GMT Lambertian lighting illumination
    illumComm = ['-E' sprintf('%g',luz.azim) '/' sprintf('%g',luz.elev) '/' sprintf('%g',luz.ambient),...
           '/' sprintf('%g',luz.diffuse) '/' sprintf('%g',luz.specular) '/' sprintf('%g',luz.shine)];
    R = grdgradient_m(Z,head,illumComm);
    handles.Illumin_type = 4;
end
img = shading_mat(img,R,'no_scale');

setappdata(handles.figure1,'illumComm',illumComm);   % Save these for region op & write_gmt_script
setappdata(handles.figure1,'Luz',luz);   % Save these for region operations
set(handles.hImg,'CData',img)
guidata(handles.figure1, handles);              set(handles.figure1,'pointer','arrow')

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
	if (handles.firstIllum),    img1 = get(handles.hImg,'CData');     handles.firstIllum = 0;
	else                        img1 = handles.origFig;      end
	if (isempty(img1)),         img1 = get(handles.hImg,'CData');     end             % No copy in memory
	if (ndims(img1) == 2),      img1 = ind2rgb8(img1,get(handles.figure1,'Colormap'));    end    % Image is 2D   
    img = shading_mat(img1,img);
else
    img = uint8((254 * img) + 1);   % Need to convert the reflectance matrix into a gray indexed image
end

if (isappdata(handles.figure1,'illumComm')),        rmappdata(handles.figure1,'illumComm');   end
handles.no_file = 0;        handles.Illumin_type = 5;
set(handles.hImg,'CData',img);
if (strcmp(color,'gray')),  set(handles.figure1,'Colormap',gray(256));    end
set(handles.figure1,'pointer','arrow');     guidata(handles.figure1, handles);

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

if (isappdata(handles.figure1,'illumComm')),        rmappdata(handles.figure1,'illumComm');   end
handles.Illumin_type = 7;       % This is temporary. It should be 7
set(handles.hImg,'CData',zz);
set(handles.figure1,'pointer','arrow');     guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function ImageAnaglyph_CB(hObject, eventdata, handles)
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
function ImageDrape_CB(hObject, eventdata, handles)
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

% Signal in the parent image handles that it has a draped image
handles = guidata(h_f);     handles.is_draped = 1;      handles.Illumin_type = 0;
guidata(h_f,handles)

% --------------------------------------------------------------------
function ImageColorPalettes_CB(hObject, eventdata, handles)
if (handles.no_file == 1),   color_palettes;     return;    end
if (ndims(get(handles.hImg,'CData')) == 3)
    warndlg('True color images do not use color palettes. So you cannot change it.','Warning'); return
end
color_palettes(handles.figure1,handles.home_dir,handles.work_dir)

% --------------------------------------------------------------------
function ToolsMeasureDist_CB(hObject, eventdata, handles)
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
msgbox(['Distance = ' sprintf('%.5f',dist) str_unit])

% --------------------------------------------------------------------
function ToolsMeasureArea_CB(hObject, eventdata, handles)
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
msgbox(['Area = ' sprintf('%g',area) unit])

% --------------------------------------------------------------------
function ToolsMeasureAzimuth_CB(hObject, eventdata, handles)
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
    msg = [msg; {['Azimuth' sprintf('%d',i) '  =  ' sprintf('%3.1f',az(i)) '  degrees']}];
end
zoom_state(handles,'maybe_on');
msgbox(msg,'Azimuth(s)')

% --------------------------------------------------------------------
function [X,Y,Z,head,m,n] = load_grd(handles, opt)
% Load a GMT grid either from memory, or re-read it again if it is to
% big to fit in it (biger than handles.grdMaxSize)
if (handles.image_type == 20),  X=[];Y=[];Z=[];head=[];   return;     end     % Fake image with bg_color. Nothing loadable
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
function DrawLine_CB(hObject, eventdata, handles, opt)
if (handles.no_file == 1),  return;     end
if (nargin == 3),           opt = [];   end
% The following is a necessary patch against two big stupidities.
% First is from the user that double-clicked on the line icon
% Second is from Matlab that doesn't lets us test a double-click on a uipushtool
if ( ~isempty(getappdata(handles.figure1, 'FromGetLine_j')) )
    return
end
zoom_state(handles,'maybe_off');

if (~isempty(opt) && strcmp(opt,'freehand')),   [xp,yp] = getline_j(handles.figure1,'freehand');
elseif (strcmp(opt,'GCPmemory')),   xp = [0 0];    % Jump the manual drawing
else                                            [xp,yp] = getline_j(handles.figure1);    end
n_nodes = length(xp);
if (n_nodes < 2),   zoom_state(handles,'maybe_on');     return;     end
% The polyline Tag is very important to destinguish from MB tracks, which have Tags = MBtrack#
if (~isempty(opt) && strcmp(opt,'FaultTrace'))    % When this function is used for Okada modeling
    if (~handles.ValidGrid),    aux_funs('msg_dlg',2,handles);     zoom_state(handles,'maybe_on'); return;      end
    lineHand = line('XData', xp, 'YData', yp,'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,'Tag','FaultTrace');
    % Create empty patches that will contain the surface projection of the fault plane
    hp = zeros(length(xp)-1);
    for (k=1:length(xp)-1), hp(k) = patch('XData', [], 'YData',[]);    end
    setappdata(lineHand,'PatchHand',hp);
elseif (strncmp(opt,'GCP',3))
    if (strcmp(opt,'GCPmemory'))
        GCPinMemory = getappdata(handles.figure1,'GCPregImage');
        xp = GCPinMemory(:,1);    yp = GCPinMemory(:,2);
    end         % else -> opt = 'GCPpline'
    lineHand = line('XData', xp, 'YData', yp,'Color','k','LineWidth',0.5,'LineStyle',':','Marker','o',...
        'MarkerFaceColor','y','MarkerSize',5,'Tag','GCPpolyline');
    if (opt(4) == 'm'),     register_img(handles,lineHand,GCPinMemory(:,1:4))      % Set uicontext for img registration
    else                    register_img(handles,lineHand)
    end
    return
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
function DrawVector_CB(hObject, eventdata, handles)
	if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
	zoom_state(handles,'maybe_off');
	draw_funs([],'DrawVector')      % Vectors are drawn there
	zoom_state(handles,'maybe_on');

% --------------------------------------------------------------------
function DrawClosedPolygon_CB(hObject, eventdata, handles, opt)
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
    [p1,p2,hl] = rubberbandbox(handles.axes1);
    set(hl,'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick)    % Use defaults LineThick and DefLineColor
    draw_funs(hl,'line_uicontext')       % Set lines's uicontextmenu
    zoom_state(handles,'maybe_on');
end

% --------------------------------------------------------------------
function DrawEulerPoleCircle_CB(hObject, eventdata, handles)
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
function DrawGeographicalCircle_CB(hObject, eventdata, handles, opt)
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
function DrawImportLine_CB(hObject, eventdata, handles, opt)
	% OPT is a string with either "AsLine", or "AsPoint", or "AsMaregraph", or "FaultTrace"
	if (handles.no_file == 1),  return;      end
	str1 = {'*.dat;*.DAT', 'Data file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select input xy file name','get');
	if isequal(FileName,0);     return;     end
	draw_funs([PathName FileName],'ImportLine',opt)

% --------------------------------------------------------------------
function DrawImportText_CB(hObject, eventdata, handles)
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
function DrawImportShape_CB(hObject, eventdata, handles)
str1 = {'*.shp;*.SHP', 'Data files (*.shp,*.SHP)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select shape file name','get');
if isequal(FileName,0);     return;     end
try         [s,t] = mex_shape([PathName FileName]);
catch       errordlg([lasterr ' Most probably, NOT a shapefile'],'Error');      return
end

lt = getappdata(handles.figure1,'DefLineThick');    lc = getappdata(handles.figure1,'DefLineColor');
region = [s(1).BoundingBox(1,1:2) s(1).BoundingBox(2,1:2)];
is_geog = aux_funs('guessGeog',region(1:4));

if (~handles.no_file && handles.geog && ~is_geog)
    errordlg('Error. Your background image is in geographics but the shape file is not','ERROR');   return
elseif (~handles.no_file && handles.geog == 0 && is_geog == 1 && ~isempty(handles.hImg))
    warndlg('WARNING: Your background image is not in geographics but the shape file seams to be. Probable mistake','Warning')
else
    handles.geog = is_geog;     guidata(handles.figure1, handles);
end

% If we have no background region
if (handles.no_file),   FileNewBgFrame_CB(handles.figure1,[],handles, [region handles.geog]);   end

nPolygs = length(s);        h = zeros(nPolygs,1);
imgLims = getappdata(handles.axes1,'ThisImageLims');
if (strcmp(t,'Arc'))
    for i = 1:nPolygs
        out = aux_funs('insideRect',imgLims,[s(i).BoundingBox(1,1) s(i).BoundingBox(2,1); s(i).BoundingBox(1,1) s(i).BoundingBox(2,2); ...
            s(i).BoundingBox(1,2) s(i).BoundingBox(2,2); s(i).BoundingBox(1,2) s(i).BoundingBox(2,1)]);
        if (any(out))       % It means the polyg BB is at least partially inside
            h(i) = line('Xdata',s(i).X,'Ydata',s(i).Y,'Parent',handles.axes1,'Color',lc,'LineWidth',lt,'Tag','SHPpolyline');
        end
        h((h == 0)) = [];   % Those were jumped because thay were completely outside map limits
    end
    draw_funs(h,'SHP_uicontext')            % Set lines's uicontextmenu
elseif (strcmp(t,'Polygon'))
    nParanoia = 1000;                       % The name talks. COMPLETELY MATLAB CONDITIONED, I WAS NOT LIKE THAT BEFORE
    for i = 1:nPolygs
        out = aux_funs('insideRect',imgLims,[s(i).BoundingBox(1,1) s(i).BoundingBox(2,1); s(i).BoundingBox(1,1) s(i).BoundingBox(2,2); ...
            s(i).BoundingBox(1,2) s(i).BoundingBox(2,2); s(i).BoundingBox(1,2) s(i).BoundingBox(2,1)]);
        if (any(out))                       % It means the polyg BB is at least partially inside
            h(i) = patch('XData',s(i).X,'YData', s(i).Y,'FaceColor','none','EdgeColor',handles.DefLineColor,'Tag','SHPpolygon');
        end
        if ((h(i) ~= 0) && nPolygs <= nParanoia)           % With luck, your hardware won't choke to dead with this
            draw_funs(h(i),'line_uicontext')
        end
    end
    if (nPolygs > nParanoia)                % nParanoia is an arbitrary number that practice will show dependency
        h((h == 0)) = [];                   % Those were jumped because thay were completely outside map limits
        draw_funs(h,'country_patch')        % mostly on hardware, for I don't beleave ML will ever behave decently.
    end
end

% --------------------------------------------------------------------
function GeophysicsImportFaultFile_CB(hObject, eventdata, handles)
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
function GeophysicsImportGmtFile_CB(hObject, eventdata, handles)
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
    FileNewBgFrame_CB(handles.figure1,[],handles, [x_min x_max y_min y_max 1])
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
function GeophysicsImportGmtFileList_CB(hObject, eventdata, handles)
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
for (k=1:m),    nao(k) = exist(names{k},'file');   end
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
    FileNewBgFrame_CB(handles.figure1,[],handles, [x_min x_max y_min y_max 1])
else
    x_lim = get(handles.axes1,'XLim');    y_lim = get(handles.axes1,'YLim');
    opt_R = ['-R' sprintf('%.4f',x_lim(1)) '/' sprintf('%.4f',x_lim(2)) '/' sprintf('%.4f',y_lim(1)) '/' sprintf('%.4f',y_lim(2))];
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
function DrawText_CB(hObject, eventdata, handles)
	if (handles.no_file == 1),    return;      end
	zoom_state(handles,'maybe_off');
	pt = ginput_pointer(1,'crosshair');
	textHand = text(pt(1),pt(2),0,'','Editing','on');
	draw_funs(textHand,'DrawText')          % Set text's uicontextmenu
	zoom_state(handles,'maybe_on');

% --------------------------------------------------------------------
function DrawSymbol_CB(hObject, eventdata, handles, opt)
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
function FileOpenSession_CB(hObject, eventdata, handles)
str1 = {'*.mat;*.MAT', 'Data files (*.mat,*.MAT)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select session file name','get');
if isequal(FileName,0);     return;     end

set(handles.figure1,'pointer','watch')
load([PathName FileName])

tala = exist(grd_name,'file');      flagIllum = true;       % Illuminate (if it is the case)
if (isempty(grd_name) || tala == 0)
    scrsz = get(0,'ScreenSize');         % Get screen size
    dx = map_limits(2) - map_limits(1);   dy = map_limits(4) - map_limits(3);
    aspect = dy / dx;
    nx = round(scrsz(3)*.6);        ny = round(nx * aspect);
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
    set(handles.figure1,'Colormap', ones( size(get(handles.figure1,'Colormap'),1), 3))
    handles = show_image(handles,'Mirone Base Map',X,Y,Z,0,'xy',1);
elseif ( ~handles.no_file && strcmp(handles.fileName,grd_name) )    % Currently loaded background is the same as in session
    set(handles.figure1,'Colormap',img_pal)     % Not harmfull anyway
    flagIllum = false;                          % Do not try to illuminate again (if it is the case)
else
    drv = aux_funs('findFileType',grd_name);
    erro = gateLoadFile(handles,drv,grd_name);  % It loads the file (or dies)
    if (erro),      set(handles.figure1,'pointer','arrow');     return;     end     % Error message already issued
    handles = guidata(handles.figure1);     % Get the updated version
end

try         % The "illumComm" variable only exists after 0.97 version and should be ~= [] only with grids
    if (~isempty(illumComm) && flagIllum)
        [X,Y,Z,head] = load_grd(handles,'silent');
        handles.Illumin_type = illumType;
        if (handles.Illumin_type == 1)
            if (handles.geog),  R = grdgradient_m(Z,head,'-M',illumComm,'-Nt');
            else                R = grdgradient_m(Z,head,illumComm,'-Nt');      end
        else
            R = grdgradient_m(Z,head,illumComm);
        end
        zz = ind2rgb8(get(handles.hImg,'CData'),get(handles.figure1,'ColorMap'));
        zz = shading_mat(zz,R,'no_scale');      set(handles.hImg,'CData',zz)         % and now it is illuminated
        setappdata(handles.figure1,'illumComm',illumComm)       % Save the 'illumComm' in this new fig
    end
end

if (haveMBtrack)                % case of MB tracks
    for i=1:length(MBtrack)
        h_line = line('Xdata',MBtrack(i).x,'Ydata',MBtrack(i).y,'Parent',handles.axes1,'LineWidth',MBtrack(i).LineWidth,...
            'color',MBtrack(i).color,'Tag',MBtrack(i).tag, 'LineStyle',MBtrack(i).LineStyle);
        setappdata(h_line,'swathRatio',MBtrack(i).swathRatio)
        draw_funs(h_line,'MBtrack_uicontext')       % Set track's uicontextmenu
    end
    for i=1:length(MBbar)       % now their's bars
        h_bar = line('Xdata',MBbar(i).x,'Ydata',MBbar(i).y,'Parent',handles.axes1,'LineWidth',MBbar(i).LineWidth,...
            'color',MBbar(i).color,'Tag',MBbar(i).tag,'UserData',MBbar(i).n_vert, 'LineStyle',MBbar(i).LineStyle);
        draw_funs(h_bar,'MBbar_uicontext')         % Set track bar's uicontextmenu
    end
    handles.h_MBplot = h_line;
end
if (haveCircleGeo)              % case of Geographic circles
    for i=1:length(CircleGeo)
        h_circ = line('Xdata',CircleGeo(i).x,'Ydata',CircleGeo(i).y,'Parent',handles.axes1,'LineWidth',CircleGeo(i).LineWidth,...
            'color',CircleGeo(i).color,'Tag',CircleGeo(i).tag, 'LineStyle',CircleGeo(i).LineStyle);
        setappdata(h_circ,'LonLatRad',CircleGeo(i).lon_lat_rad);
        CircleGeo(i).ud.hcirc = h_circ;                 CircleGeo(i).ud.parent = handles.axes1;
        CircleGeo(i).ud.h_fig = handles.figure1;        CircleGeo(i).ud.h_axes = handles.axes1;
        set(h_circ,'UserData',CircleGeo(i).ud,'buttondownfcn','uicirclegeo(''circlemousedown'')')
        draw_funs(h_circ,'SessionRestoreCircle')       % Set circle's uicontextmenu
    end
end
if (haveCircleCart)             % case of Cartesian circles
    for i=1:length(CircleCart)
        h_circ = line('Xdata',CircleCart(i).x,'Ydata',CircleCart(i).y,'Parent',handles.axes1,'LineWidth',CircleCart(i).LineWidth,...
            'color',CircleCart(i).color,'Tag',CircleCart(i).tag, 'LineStyle',CircleCart(i).LineStyle);
        setappdata(h_circ,'LonLatRad',CircleCart(i).lon_lat_rad);
        x = linspace(-pi,pi,360);
	    setappdata(h_circ,'X',cos(x));       setappdata(h_circ,'Y',sin(x))    % Save unit circle coords
        CircleCart(i).ud.hcirc = h_circ;
        CircleCart(i).ud.parent = handles.axes1;
        draw_funs(h_circ,'SessionRestoreCircleCart')       % Set circle's uicontextmenu
    end
end
if (havePline)                  % case of polylines
    for i=1:length(Pline)
        h_line = line('Xdata',Pline(i).x,'Ydata',Pline(i).y,'Parent',handles.axes1,'LineWidth',Pline(i).LineWidth,...
            'color',Pline(i).color,'Tag',Pline(i).tag, 'LineStyle',Pline(i).LineStyle);
        if (isfield(Pline(i),'LineInfo') && ~isempty(Pline(i).LineInfo))
            setappdata(h_line,'LineInfo',Pline(i).LineInfo)
            set(h_line,'UserData',1)
            draw_funs(h_line,'isochron',{Pline(i).LineInfo})
        else
            draw_funs(h_line,'line_uicontext')       % Set lines's uicontextmenu
        end
    end
end
if (havePlineAsPoints)          % case of polylines as points (markers) only
    for i=1:length(PlineAsPoints)
        h_line_pt = line('Xdata',PlineAsPoints(i).x, 'Ydata',PlineAsPoints(i).y,'Parent',handles.axes1, 'LineStyle','none', ...
            'Marker',PlineAsPoints(i).Marker, 'MarkerSize',PlineAsPoints(i).Size, ...
            'MarkerFaceColor',PlineAsPoints(i).FillColor, ...
            'MarkerEdgeColor',PlineAsPoints(i).EdgeColor, 'Tag',PlineAsPoints(i).tag);
        draw_funs(h_line_pt,'DrawSymbol')        % Set marker's uicontextmenu (tag is very important)        
    end
end
if (haveSymbol)                 % case of Symbols (line Markers)
    for i=1:length(Symbol)
        h_symb = line('Xdata',Symbol(i).x,'Ydata',Symbol(i).y,'Parent',handles.axes1,'Marker',Symbol(i).Marker,'MarkerSize',...
            Symbol(i).Size,'MarkerFaceColor',Symbol(i).FillColor, 'MarkerEdgeColor',Symbol(i).EdgeColor, 'Tag',Symbol(i).tag);
        draw_funs(h_symb,'DrawSymbol')          % Set symbol's uicontextmenu
    end
end
if (haveText)                   % case of text strings
    for i=1:length(Texto)
        if (isempty(Texto(i).str)),  continue;   end
        h_text = text(Texto(i).pos(1),Texto(i).pos(2),Texto(i).pos(3), Texto(i).str,...
            'Parent',handles.axes1, 'Rotation',Texto(i).angle,...
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
%             if (Patches(i).FaceColor(1) == 1),    Patches(i).tag = 'tapete';
%             else                                  Patches(i).tag = 'tapete_R';
%             end
            h_patch = patch('XData',Patches(i).x, 'YData',Patches(i).y, 'Parent',handles.axes1,'LineWidth',Patches(i).LineWidth,...
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
if (haveCoasts),    DatasetsCoastLineNetCDF_CB([], [], handles,coastUD);  end
if (havePolitic)
    if (iscell(politicUD)),     politicUD = politicUD{1};     end
    DatasetsPoliticalBoundaries_CB([], [], handles,politicUD(2),politicUD(1));
end
if (haveRivers)
    if (iscell(riversUD)),      riversUD = riversUD{1};     end
    DatasetsRivers_CB([], [], handles,riversUD(2),riversUD(1));
end
guidata(hObject, handles);
set(handles.figure1,'pointer','arrow','Name',[PathName FileName])
if (tala == 0 && ~isempty(grd_name))    % Only now to not mess with the "current figure"
    warndlg(['The file ' grd_name ' doesn''t exists on the directory it was when the session was saved. Put it back there.'],'Warning')
end

% ------------------------------------------------------------------------------------------------
function FileSaveSession_CB(hObject, eventdata, handles)
if (handles.image_type == 0),   return;      end
str1 = {'*.mat;*.MAT', 'Data files (*.mat,*.MAT)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select session file name','put');
if isequal(FileName,0);     return;     end

set(handles.figure1,'pointer','watch')
[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if isempty(EXT),    fname = [PathName FNAME '.mat'];
else                fname = [PathName FNAME EXT];       end

grd_name = handles.fileName;    % Use this variable name for compatibility reason
img_pal = get(handles.figure1,'Colormap');     illumComm = [];      illumType = handles.Illumin_type;
map_limits = getappdata(handles.figure1,'ThisImageLims');
if (handles.ValidGrid && handles.Illumin_type >= 1 && handles.Illumin_type <= 4)
    illumComm = getappdata(handles.figure1,'illumComm');
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
        if (isappdata(ALLlineHand(i),'LineInfo'))
            Pline(m).LineInfo = getappdata(ALLlineHand(i),'LineInfo');
        end
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
    'haveCoasts', 'coastUD','havePolitic', 'politicUD','haveRivers', 'riversUD', 'illumComm', 'illumType')
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DatasetsCoastLineNetCDF_CB(hObject, eventdata, handles, res)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
if (nargin == 3),   res = 'l';  end
set(handles.figure1,'pointer','watch')

lon = get(handles.axes1,'Xlim');      lat = get(handles.axes1,'Ylim');
opt_R = ['-R' sprintf('%.4f',lon(1)) '/' sprintf('%.4f',lon(2)) '/' sprintf('%.4f',lat(1)) '/' sprintf('%.4f',lat(2))];

switch res
    case 'c',        opt_res = '-Dc';        pad = 2.0;
    case 'l',        opt_res = '-Dl';        pad = 0.5;
    case 'i',        opt_res = '-Di';        pad = 0.1;
    case 'h',        opt_res = '-Dh';        pad = 0.02;
    case 'f',        opt_res = '-Df';        pad = 0.005;
end
coast = shoredump(opt_R,opt_res,'-A1/1/1');

% Get rid of data that are outside the map limits
lon = lon - [pad -pad];     lat = lat - [pad -pad];
indx = (coast(1,:) < lon(1) | coast(1,:) > lon(2));
coast(:,indx) = [];
indx = (coast(2,:) < lat(1) | coast(2,:) > lat(2));
coast(:,indx) = [];
coast = single(coast);

if (~all(isnan(coast(:))))
	h = line('XData',coast(1,:),'YData',coast(2,:),'Parent',handles.axes1,'Linewidth',handles.DefLineThick,...
        'Color',handles.DefLineColor,'Tag','CoastLineNetCDF','UserData',opt_res(3));
	draw_funs(h,'Coastline_uicontext')    % Set line's uicontextmenu
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DatasetsPoliticalBoundaries_CB(hObject, eventdata, handles, type, res)
% TYPE is: '1' -> National Boundaries
%          '2' -> State Boundaries
%          '3' -> Marine Boundaries
%          'a' -> All Boundaries
% RES is:  'c' or 'l' or 'i' or 'h' or 'f' (gmt database resolution)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file

set(handles.figure1,'pointer','watch')
lon = get(handles.axes1,'Xlim');      lat = get(handles.axes1,'Ylim');
opt_R = ['-R' sprintf('%.4f',lon(1)) '/' sprintf('%.4f',lon(2)) '/' sprintf('%.4f',lat(1)) '/' sprintf('%.4f',lat(2))];

switch type
    case '1',        opt_N = '-N1';
    case '2',        opt_N = '-N2';
    case '3',        opt_N = '-N3';
    case 'a',        opt_N = '-Na';
end

switch res
    case 'c',        opt_res = '-Dc';        pad = 2;
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
	h = line('XData',boundaries(1,:),'YData',boundaries(2,:),'Parent',handles.axes1,'Linewidth',handles.DefLineThick,...
        'Color',handles.DefLineColor,'Tag','PoliticalBoundaries', 'UserData',[opt_res(3) opt_N(3)]);
	draw_funs(h,'Coastline_uicontext')    % Set line's uicontextmenu
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DatasetsRivers_CB(hObject, eventdata, handles, type, res)
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
opt_R = ['-R' sprintf('%.4f',lon(1)) '/' sprintf('%.4f',lon(2)) '/' sprintf('%.4f',lat(1)) '/' sprintf('%.4f',lat(2))];

switch type
    case '1',        opt_I = '-I1';         case '2',        opt_I = '-I2';
    case '3',        opt_I = '-I3';
    case '5',        opt_I = '-I5';         case '6',        opt_I = '-I6';
    case '7',        opt_I = '-I7';
    case 'a',        opt_I = '-Ia';
    case 'r',        opt_I = '-Ir';         case 'i',        opt_I = '-Ii';
end

switch res
    case 'c',        opt_res = '-Dc';        pad = 2;
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
	h = line('XData',rivers(1,:),'YData',rivers(2,:),'Parent',handles.axes1,'Linewidth',handles.DefLineThick,...
        'Color',handles.DefLineColor,'Tag','Rivers', 'UserData',[opt_res(3) opt_I(3:end)]);
	draw_funs(h,'Coastline_uicontext')    % Set line's uicontextmenu
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function ImageMapLimits_CB(hObject, eventdata, handles)
% Change the Image limits by asking it's corner coordinates
region = bg_region('empty');     % region contains [x_min x_max y_min y_max is_geog]
if isempty(region),     return;     end     % User gave up
img = get(findobj(handles.axes1,'Type','image'),'CData');
X = region(1:2);                            Y = region(3:4);
x_inc = diff(X) / size(img,2);              y_inc = diff(Y) / size(img,1);
dx2 = x_inc / 2;                            dy2 = y_inc / 2;
X = X + [dx2 -dx2];                         Y = Y + [dy2 -dy2];     % Make it such that the pix-reg info = region
handles.head = [X Y 0 255 0 x_inc y_inc];   handles.image_type = 3; handles.fileName = [];  % Not loadable in session

% Flipud the image if necessary
if (strcmp(get(handles.axes1,'YDir'),'reverse')),    img = flipdim(img,1);      end
show_image(handles,'New Limits',X,Y,img,0,'xy',0,1);

% --------------------------------------------------------------------
function GeophysicsTTT_CB(hObject, eventdata, handles,opt)
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
    if isempty(Z),   return;    end;                % An error message was already issued
    if (handles.have_nans),     errordlg('Bathymetry grid cannot have NaNs','Error');   return;    end
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
function GeophysicsSwanCompute_CB(hObject, eventdata, handles)

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
    if (I_am_in_bat && I_am_in_source)       % If bat & source figures test that they are compatible.
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
if (~handles.IamCompiled),  swan_hand = @swan;          % To stop once for all with the bloody version mixing
else                        swan_hand = @swan_sem_wbar;
end

if (isfield(out,'maregraph_xy'))    % Ask for computation of maregraphs
    if (isfield(out,'opt_m'))       % Movie option
        tmovie = feval(swan_hand,Z_bat, head_bat, Z_src, head_src, out.params, out.maregraph_xy, extra_args1, ...
                extra_args2,extra_args3,extra_args5,'-f');
    else
        feval(swan_hand,Z_bat, head_bat, Z_src, head_src, out.params, out.maregraph_xy, extra_args1, ...
                extra_args2,extra_args3,extra_args5);
    end    
else                                % Compute grids or movie
    if (isfield(out,'opt_m'))       % Movie option
        tmovie = feval(swan_hand,Z_bat, head_bat, Z_src, head_src, out.params, extra_args2,extra_args3,extra_args5,'-f');
    else
        feval(swan_hand,Z_bat, head_bat, Z_src, head_src, out.params, extra_args2,extra_args3,extra_args5);
    end
end
pause(0.01);        % Give time to waitbar window to die
if (isfield(out,'opt_m')),   do_movie(handles,tmovie,'swan');   end

% --------------------------------------------------------------------
function GeophysicsSwanPlotStations_CB(hObject, eventdata, handles)
	if (aux_funs('msg_dlg',14,handles));     return;      end
	zoom_state(handles,'maybe_off');
	pt = ginput_pointer(1,'crosshair');
	if (isempty(pt)),   return;  end
	handles.maregraphs_count = handles.maregraphs_count + 1;    % Count number of maregraphs
	%tag = ['Maregraph' num2str(handles.maregraphs_count)];
	symbHand = line(pt(1,1),pt(1,2),'Marker','o','MarkerFaceColor','y',...
            'MarkerEdgeColor','k','MarkerSize',10,'Tag','Maregraph');
	draw_funs(symbHand,'DrawSymbol')          % Set symbol's uicontextmenu
	zoom_state(handles,'maybe_on');	    guidata(hObject, handles)

% --------------------------------------------------------------------
function GeophysicsSwanGridBorderStations_CB(hObject, eventdata, handles)
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

if (handles.grdformat > 1), FileName = [FileName '=' sprintf('%d',handles.grdformat)];     end
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
        str_R = [' -R' sprintf('%.10f',new_w) '/' sprintf('%.10f',new_s) '/' sprintf('%.10f',new_e) '/' sprintf('%.10f',new_n)];
    elseif (adjust_w)        % Need to adjust only est and west borders
        if (in_w_to_e)      % Move x origin to the east
            new_w = fix(dx_w) * (head(8) + 1);
            new_e = new_w + D(2) - D(1);
        else                % Move x origin to the west
            new_w = fix(dx_w) * head(8);
            new_e = new_w + D(2) - D(1);
        end
        str_R = [' -R' sprintf('%.10f',new_w) '/' sprintf('%.10f',D(3)) '/' sprintf('%.10f',new_e) '/' sprintf('%.10f',D(4))];
    elseif (adjust_s)        % Need to adjust only south and north borders
        if (in_s_to_n)      % Move y origin to the north
            new_s = fix(dy_s) * (head(9) + 1);
            new_n = new_s + D(4) - D(3);
        else                % Move y origin to the south
            new_s = fix(dy_s) * head(9);
            new_n = new_s + D(4) - D(3);
        end
        str_R = [' -R' sprintf('%.10f',D(1)) '/' sprintf('%.10f',new_s) '/' sprintf('%.10f',D(2)) '/' sprintf('%.10f',new_n)];
    else
        errordlg('Asneira desconhecida (Unknown error).','Error')
    end
end

if (adjust_x_inc)       % Nao vou testar/usar o y_inc. A malha tem (tera?) de ser quadrada
    str_I = [' -I' sprintf('%.10f',new_x_inc)];
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
function GeophysicsTsun2_CB(hObject, eventdata, handles, opt)
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
if (isfield(out,'opt_I')),  extra_args6 = out.opt_I;   end

[X,Y,Z_bat,head_bat] = load_grd(handles);   Z_bat = double(Z_bat);
if (isempty(Z_bat) || isempty(head_bat))
    errordlg('ERROR: one or more of the bat variables are empty where they souldn''t be.','Error');
    return;
end
if ( abs(head_bat(8) - head_bat(9)) > 1e-3 )
    warndlg('Grid cells are not square. I don''t know the effect of this.','Warning')
end
if (~handles.IamCompiled),  tsun2_hand = @tsun2;          % To stop once for all with the bloody version mixing
else                        tsun2_hand = @tsun2_sem_wbar;
end

if (isfield(out,'maregraph_xy'))    % Ask for computation of maregraphs (WRONG - This is not an option, but will be in future)
    dt1 = diff(out.maregraph_xy(:,1));    dt2 = diff(dt1);      t0 = out.maregraph_xy(1,1);
    if any(abs(dt2) > 100*eps)      % First column doesn't have the time (eps due to very small dts rounding errors)
        warndlg('The maregraph file does not have a time increment. I will assume it is 1 second.','SEVERE WARNING')
        extra_args6 = ['-I' num2str(head_bat(8)) '/1'];
    else                    % First column of maregs file has the time. We don't want it
        dt = dt1(1);        % This is the time increment to be used as option to tsun2
        out.maregraph_xy(:,1) = [];
        cfl = head_bat(8) / sqrt(abs(head_bat(5))*9.8);
        if (cfl <= dt)
            warndlg('Your wave propagates faster than the maregraph time increment. Expect divergent results.','SEVERE WARNING')
        end
        extra_args6 = ['-I' num2str(head_bat(8)) '/' num2str(dt)];
        if (t0 > 1),    extra_args7 = ['-J' num2str(t0)];   end     % If we start computing at a latter time
    end
    if (isfield(out,'opt_m'))       % Movie option
        tmovie = feval(tsun2_hand,Z_bat, head_bat, out.maregraph_xy, extra_args1, ...
            extra_args2,extra_args3,extra_args5,extra_args6,extra_args7,'-f');
    else
        feval(tsun2_hand,Z_bat, head_bat, out.maregraph_xy, extra_args1,extra_args2,extra_args3,extra_args5,extra_args6,extra_args7);
    end    
else                                % Compute grids or movie
    if (isfield(out,'opt_m'))       % Movie option
        tmovie = feval(tsun2_hand,Z_bat, head_bat, extra_args2,extra_args3,extra_args5,extra_args6,'-f');
    else
        feval(tsun2_hand,Z_bat, head_bat, extra_args2,extra_args3,extra_args5,extra_args6);
    end
end
pause(0.01);        % Give time to waitbar window to die

if (isfield(out,'opt_m')),   do_movie(handles,tmovie,'tsun2');   end

%------------------------------------------------------------------------------------------
function do_movie(handles,tmovie,opt)
if (strcmp(opt,'swan')),    is_swan = 1;
else                        is_swan = 0;    end
[m,n,k] = size(tmovie);
opt_E = '-E60/30/0.55/0.6/0.4/10';      % Should be "controlable"
n_crop = 0;    % Should be a variabe. It's used to crop borders and therefore hide grid borders reflections
[X,Y,Z0,head] = load_grd(handles);
if isempty(Z0),   return;     end;    % An error message was already issued
if (n_crop),    Z0 = double(Z0(n_crop+1:m-n_crop, n_crop+1:n-n_crop));
else            Z0 = double(Z0);
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
    Z = tmovie(:,:,1);
    tmovie(:,:,1) = [];      % Free this page (we don't need it anymore)
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
    
    img = ind2rgb8(img,cmap);   img = shading_mat(img,R);
    M(i) = im2frame(img);
    waitbar(i/k);               i = i + 1;
end
close(h)

str1 = {'*.avi;*.AVI', 'avi files (*.avi,*.AVI)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select movie name','put');
if isequal(FileName,0);     return;     end
set(handles.figure1,'pointer','watch')
movie2avi_j(M,[PathName FileName],'compression','none','fps',5)
set(handles.figure1,'pointer','arrow');

% --------------------------------------------------------------------
function GRD_save_or_display(handles,X,Y,Z,head,tit,name)
% Choose what to do with Z. Show it in a new window or save it as a gmt grid.
% TIT is used only when saving a grid
if (nargin < 7),    name = [];  end
if (handles.out_in_NewWindow)
    if (isa(Z,'double')),   Z = single(Z);  end
    zzz = grdutils(Z,'-L');        head(5) = zzz(1);       head(6) = zzz(2);
    tmp.head = [head(1) head(2) head(3) head(4) head(5) head(6) head(7) head(8) head(9)];
    tmp.X = X;    tmp.Y = Y;    tmp.name = name;
    mirone(Z,tmp);
else
    str1 = {'*.grd;*.GRD','netCDF grid format (*.grd,*.GRD)'; '*.*', 'All Files (*.*)'};
    [FileName,PathName] = put_or_get_file(handles,str1,'Select output GMT grid','put');
    if isequal(FileName,0);     return;     end

    [PATH,FNAME,EXT] = fileparts([PathName FileName]);
    if (handles.grdformat == 6)       % Preserve the original format
        if isempty(EXT),    f_name = [PathName FNAME '.grd=' sprintf('%d',handles.grdformat)];
        else                f_name = [PathName FNAME EXT '=' sprintf('%d',handles.grdformat)];       end
    else
        if isempty(EXT),    f_name = [PathName FNAME '.grd=nf'];
        else                f_name = [PathName FNAME EXT '=nf'];       end    
    end
    grdwrite_m(Z,head,f_name,tit);     pause(0.01);
end
set(handles.figure1,'pointer','arrow');

% --------------------------------------------------------------------
function FileSaveImgGrdGdal_CB(hObject, eventdata, handles, opt1, opt2)
% OPT1 = DRIVER == GTiff, HFA (erdas), ENVI, ECW, JP2ECW
% OPT2 == grid -> saves the underlaying grid
% ELSE -> do a screen capture
if (handles.no_file == 1),    return;      end
if (strcmp(opt2,'grid') && ~handles.ValidGrid)
    errordlg('You don''t have a Grid loaded, so OBVIOUSLY you cannot save it.','Error');  return
end
flip = 0;           % FLIP = 1, when images need to be UD fliped before saving
switch opt1
    case 'GeoTiff',     str1 = {'*.tif;*.TIF;*.tiff;*.TIFF', 'GeoTiff (*.tif;*.TIF;*.tiff;*.TIFF)'};      driver = 'GTiff';
    case 'Erdas',       str1 = {'*.img;*.IMG', 'Erdas (*.img;*.IMG)'};          driver = 'HFA';
    case 'Envi',        str1 = {'*.img;*.IMG', 'Envi (*.img;*.IMG)'};           driver = 'ENVI';
    case 'ESRI',        str1 = {'*.bil;*.BIL', 'Envi (*.bil;*.BIL)'};           driver = 'EHdr';
    case 'JP2K',        str1 = {'*.jp2;*.JP2', 'Jpeg2000 (*.jp2;*.JP2)'};       driver = 'JP2ECW';
end
[FileName,PathName] = put_or_get_file(handles,str1,['Select ' opt1 ' file name'],'put');
if isequal(FileName,0);     return;     end
[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if isempty(EXT),    FileName = [FileName str1{1}(2:5)];  end
fname = [PathName FileName];

head = handles.head;
% 'ThisImageLims' contains the limits as seen from the pixel-registration stand-point,
%  and that is the convention that we will use to save rasters using GDAL.
imgLims = getappdata(handles.figure1,'ThisImageLims');
if (strcmp(opt2,'grid'))
    [X,Y,Z,head] = load_grd(handles);
    if isempty(Z),   return;     end;    % An error message was already issued
    if (handles.was_int16)
        if (handles.have_nans)      % Restore the original Nodata value, or use -32768 if we don't know it
            Z(isnan(Z(:))) = handles.Nodata_int16;
        end
        Z = int16(Z);
    end
else                                    % 'image'
    img = snapshot(handles.figure1,'noname');       pause(0.01)
    if (isempty(img)),      return;     end
    head(8) = (imgLims(2)-imgLims(1)) / size(img,2);
    head(9) = (imgLims(4)-imgLims(3)) / size(img,1);
    flip = 0;
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
	elseif (strcmp(type,'put'))
        cd(handles.work_dir)
        [FileName,PathName] = uiputfile(str1,str2);
	end
    if (PathName ~= 0),         handles.last_dir = PathName;    end
	pause(0.01);        cd(handles.home_dir);       % allways go home
	if (~isempty(strfind([PathName FileName],' ')))
        errordlg('If you had RTFM you should know that names (path included) with white spaces are totaly FORBIDEN here.','ERROR')
        FileName = 0;   % Make calling routine quit
	end
	guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function GridToolsGrdfilter_CB(hObject, eventdata, handles)
	if (aux_funs('msg_dlg',14,handles));     return;      end
	[X,Y,Z] = load_grd(handles);
	if isempty(Z),  return;     end     % An error message was already issued
	grdfilter_Mir(handles.figure1,handles.geog,Z);

% --------------------------------------------------------------------
function GridToolsGrdsample_CB(hObject, eventdata, handles)
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
function GridToolsGrdGrad_CB(hObject, eventdata, handles)
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
	[zzz] = grdutils(newZ,'-L');    head(5) = zzz(1);     head(6) = zzz(2);
	GRD_save_or_display(handles,X,Y,newZ,head,'Gradient grid','Gradient grid')

% --------------------------------------------------------------------
function GridToolsGrdtrend_CB(hObject, eventdata, handles)
	if (aux_funs('msg_dlg',14,handles));     return;      end
	out = grdtrend_Mir(handles.figure1);      pause(0.01);
	if (isempty(out)),      return;     end
	[X,Y,Z,head] = load_grd(handles);
	if isempty(Z),  return;     end;    % An error message was already issued
	set(handles.figure1,'pointer','watch')
	newZ = grdtrend_m(Z,head,out.opt_what,out.opt_N);
	[zzz] = grdutils(newZ,'-L');    head(5) = zzz(1);     head(6) = zzz(2);
	GRD_save_or_display(handles,X,Y,newZ,head,'Grdtrend grid','Grdtrend grid')

% --------------------------------------------------------------------
function GridToolsGrdproject_CB(hObject, eventdata, handles)
	% Call the geographic calculator in the grid only mode.
	if (aux_funs('msg_dlg',14,handles));     return;      end
	[X,Y,Z,head] = load_grd(handles);
	if isempty(Z),          return;     end     % An error message was already issued
	out = geog_calculator(Z,head,handles.figure1,'onlyGrid');
	if (isempty(out)),      return;     end     % User just gave up
	Z = out.Z;      head = out.head;
	[zzz] = grdutils(Z,'-L');    head(5) = zzz(1);     head(6) = zzz(2);
	GRD_save_or_display(handles,X,Y,Z,head,'Grdprojected grid','Grdprojected grid')

% --------------------------------------------------------------------
function GridToolsGridClip_CB(hObject, eventdata, handles)
	if (aux_funs('msg_dlg',14,handles));     return;      end
	[X,Y,Z,head] = load_grd(handles);    %  load gmt grid
	if isempty(Z),   return;     end     % An error message was already issued
	out = ml_clip(head(5),head(6));
	if isempty(out),        return;     end;
	if ~isempty(out{2}),    Z(Z > out{1}) = out{2};     end;     % Clip above
	if ~isempty(out{4}),    Z(Z < out{3}) = out{4};     end;     % Clip below
	GRD_save_or_display(handles,X,Y,Z,head,'Cliped grid','Cliped grid')

% --------------------------------------------------------------------
function GridToolsHistogram_CB(hObject, eventdata, handles, opt)
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
resp = inputdlg({'Enter Bin Width (default is 20 bins)'},'Histogram',[1 38],{sprintf('%g',binwidth)});     pause(0.01);
if isempty(resp);    set(handles.figure1,'pointer','arrow');     return;     end
n = round( (z_max - z_min) / str2double(resp{1}) );
[n,xout] = histo_m('hist',Z(:),n,[z_min z_max]);
h = mirone;                         % Create a new Mirone figure
mirone('FileNewBgFrame_CB',h,[],guidata(h), [xout(1) xout(end) 0 max(n) 0], [600 600]);
set(h,'Name','Grid Histogram');       axes(get(h,'CurrentAxes'));
histo_m('bar',xout,n,'hist');
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function GridToolsGridMask_CB(hObject, eventdata, handles)
	if (aux_funs('msg_dlg',14,handles));     return;      end
	if (~handles.have_nans)
        msgbox('This option only works on grids that have NaNs.','Warning');    return;
	end
	[X,Y,Z,head] = load_grd(handles);    %  load gmt grid
	if isempty(Z),   return;     end;    % An error message was already issued
	Z(~isnan(Z)) = 1;
	GRD_save_or_display(handles,X,Y,Z,head,'Mask grid','Mask grid')

% --------------------------------------------------------------------
function GridToolsSectrum_CB(hObject, eventdata, handles, opt1, opt2)
% OPT1 == 'Amplitude'   -> compute amplitude spectrum
% OPT1 == 'Power'       -> compute power spectrum
% OPT1 == 'Autocorr'    -> compute autocorrelation
% OPT1 == 'Allopts'     -> call the fft_stuff window
% OPT2 if present is a structure with two fields: opt2.Z (the matrix); opt2.head (its 9 elements header)
if (aux_funs('msg_dlg',14,handles));     return;      end
if (handles.have_nans)
    warndlg('This grid has NaNs. That is not allowed in FFTs','Warning');    return;
end
if (nargin == 4)                        % Use entire grid
    [X,Y,Z,head] = load_grd(handles);   quick = 0;
    if isempty(Z),  return;     end;    % An error message was already issued
else                                    % Use a subset grid extracted from a rectangular area
    Z = opt2.Z;             head = opt2.head;       quick = 1;
end

if (quick),     set(handles.figure1,'pointer','watch'); end
fft_stuff(handles.figure1, Z, head, handles.geog, opt1);
if (quick),     set(handles.figure1,'pointer','arrow'); end

% --------------------------------------------------------------------
function GridToolsSmooth_CB(hObject, eventdata, handles)
if (aux_funs('msg_dlg',14,handles));     return;      end
[X,Y,Z,head,m,n] = load_grd(handles);       %  If gmt grid is not in memory, so read it again
if isempty(Z),      return;     end;        % An error message was already issued
if (~isa(Z,'double')),  Z = double(Z);  end;            % Make sure Z is of double type

[pp p_guess] = spl_fun('csaps',{Y(1:5),X(1:5)},Z(1:5,1:5));   % Get a good estimate of p
prompt = {'Enter smoothing p paramer'};     dlg_title = 'Smoothing parameter input';
defAns = {sprintf('%.12f',p_guess{1})};     resp  = inputdlg(prompt,dlg_title,[1 38],defAns);
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
function GridToolsSDG_CB(hObject, eventdata, handles, opt)
if (aux_funs('msg_dlg',14,handles));     return;      end

[X,Y,Z,head] = load_grd(handles);           [m,n] = size(Z);
if isempty(Z),   return;     end;    % An error message was already issued
if (~isa(Z,'double')),  Z = double(Z);  end;            % Make sure Z is of double type
[pp p_guess] = spl_fun('csaps',{Y(1:5),X(1:5)},Z(1:5,1:5));   % Get a good estimate of p
prompt = {'Enter smoothing p paramer'};     dlg_title = 'Smoothing parameter input';
defAns = {sprintf('%.12f',p_guess{1})};     resp  = inputdlg(prompt,dlg_title,[1 38],defAns);
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
function GridToolsSlope_CB(hObject, eventdata, handles, opt)
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
function GridToolsDirDerive_CB(hObject, eventdata, handles, opt)
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
function GridToolsFindHoles_CB(hObject, eventdata, handles)
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
function GridToolsSaveAsSRTM_CB(hObject, eventdata, handles)
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
name = [n sprintf('%.2d',abs(round(head(3)))) w sprintf('%.3d',abs(round(head(1)))) '_p.hgt'];
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
function GridToolsPadd2Const_CB(hObject, eventdata, handles)
% Pad the array to a const value (currently ct = zero) using a Hanning window
if (aux_funs('msg_dlg',14,handles));     return;      end
[X,Y,Z,head,m,n] = load_grd(handles);
if isempty(Z),   return;     end;    % An error message was already issued
resp  = inputdlg({'Enter number of border lines'},'Skirt width',[1 38],{'10'});
pause(0.01)
if isempty(resp);   return;     end
n_pad = str2double(resp{1}) * 2;
Z = mboard(Z,n,m,n+n_pad,m+n_pad);
zzz = grdutils(Z,'-L');  head(5) = zzz(1);  head(6) = zzz(2);     clear zzz;
head(1) = head(1) - n_pad/2 * head(8);      head(2) = head(2) + n_pad/2 * head(8);
head(3) = head(3) - n_pad/2 * head(9);      head(4) = head(4) + n_pad/2 * head(9);
X = linspace(head(1),head(2),n+n_pad);      Y = linspace(head(3),head(4),m+n_pad);
GRD_save_or_display(handles,X,Y,Z,head,'Padded Grid','Padded Grid')

% --------------------------------------------------------------------
function FileSaveFlederSD_CB(hObject, eventdata, handles, opt)
% Depending of the OPT value, this function builds either:
% OPT = 'writePlanarSD' directly build a planar Sonar SD file to be used by Fledermaus
% OPT = 'writeSphericalSD' directly build a spherical Sonar SD file to be used by Fledermaus
% OPT = 'runPlanarSD' build a planar .sd file (but don't keep it) and run the viewer
% OPT = 'runSphericalSD' build a spherical .sd file (but don't keep it) and run the viewer
% Note: 20-1-07 planar or spherical obj may be set on preferences. 
if (aux_funs('msg_dlg',14,handles));     return;      end
if (nargin == 3),   opt = 'runPlanarSD';   end
if ( (strcmp(opt,'writeSphericalSD') || ~handles.flederPlanar) && ~handles.geog)
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
if (strcmp(opt,'writePlanarSD') || handles.flederPlanar)
    write_flederFiles('main_SD',fid,'writePlanarSD',handles.figure1,handles.axes1,Z,head(1:6),handles.flederBurn)
else
    write_flederFiles('main_SD',fid,'writeSphericalSD',handles.figure1,handles.axes1,Z,head(1:6),handles.flederBurn)
end
if (handles.flederBurn ~= 2)       % If not screen capture, see if there are lines &pts to flederize
    write_flederFiles('line_or_points',fid,handles.figure1,handles.axes1,Z,head(1:6),handles.flederBurn)
end
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

% -----------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata, handles)
    h = getappdata(handles.figure1,'dependentFigs');
    delete(handles.figure1);        delete(h(ishandle(h)))  % Delete also any eventual 'carraças'
    
% --------------------------------------------------------------------
function ImageEdgeDetect_CB(hObject, eventdata, handles, opt)
if (handles.no_file == 1),      return;      end

if (~strcmp(opt,'ppa')),        img = get(handles.hImg,'CData'); end
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
function ImageMovieFromList_CB(hObject, eventdata, handles)
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
function DigitalFilt_CB(hObject, eventdata, handles, opt)
if (handles.no_file == 1),     return;      end
if (strcmp(opt,'image'))
    digitalFiltering(handles.hImg);
else        % grid
    if (aux_funs('msg_dlg',14,handles));     return;      end
	[X,Y,Z,head] = load_grd(handles);   % load the grid array here
	if isempty(Z),      return;     end;    % An error message was already issued
	[Z, img] = digitalFiltering(handles.hImg,Z,get(handles.figure1,'ColorMap'));
	if (isempty(Z)),    return;     end
	
    [zzz] = grdutils(Z,'-L');    head(5) = zzz(1);     head(6) = zzz(2);
	setappdata(handles.figure1,'dem_z',Z);    setappdata(handles.figure1,'GMThead',head);
	setappdata(handles.figure1,'Zmin_max',[head(5) head(6)])
	handles.head = head;        handles.origFig = img;
	guidata(handles.figure1,handles)
end

% --------------------------------------------------------------------
function RotateTool_CB(hObject, eventdata, handles, opt)
if (handles.no_file == 1),     return;      end
if (strcmp(opt,'image'))
    img = rotatetool(get(handles.hImg,'CData'),get(handles.figure1,'ColorMap'));
    if (ndims(img) == 2),   setappdata(0,'CropedColormap',get(handles.figure1,'ColorMap'));     end
    mirone(img);
else        % grid
    if (aux_funs('msg_dlg',14,handles));     return;      end
	[X,Y,Z,head] = load_grd(handles);       % load the grid array here
	if isempty(Z),      return;     end     % An error message was already issued
	[newZ, hdr] = rotatetool(Z,head,get(handles.figure1,'ColorMap'));
	if (isempty(newZ)),    return;     end
    [ny,nx] = size(newZ);
    X = linspace(hdr(1),hdr(2),nx);       Y = linspace(hdr(3),hdr(4),ny);
    GRD_save_or_display(handles,X,Y,newZ,hdr,'Rotated grid','Rotated grid')
end

% --------------------------------------------------------------------
function ImageGCPtool_CB(hObject, eventdata, handles)
	if (handles.no_file == 1),     return;      end
	if (~strcmp(get(hObject,'Checked'),'on'))
        set(hObject,'Checked','on')
	else
        set(hObject,'Checked','off')
	end

% --------------------------------------------------------------------
function Transfer_CB(dumb1, dumb2, handles, opt)
if (any(strcmp({'Corners' 'gray' 'bw'},opt)) && handles.no_file),      return;      end

if (strcmp(opt,'Shape')),       floodFill(handles.figure1);     return;     end
set(handles.figure1,'pointer','watch')
img = get(handles.hImg,'CData');
if (strcmp(opt,'Corners'))
    corn = cvlib_mex('goodfeatures',img,100,0.05);
    lineHand = line('XData',corn(:,1),'YData',corn(:,2),'Parent',handles.axes1,'ko','MarkerEdgeColor','w', ...
                    'MarkerFaceColor','k','MarkerSize',6,'Tag','corner_detected','Userdata',1);
	multi_segs_str = cell(length(lineHand),1);    % Just create a set of empty info strings
	draw_funs(lineHand,'isochron',multi_segs_str);
elseif (strcmp(opt,'gray'))
    if (ndims(img) == 3)
        img = cvlib_mex('color',img,'rgb2gray');
        set(handles.hImg,'CData', img);
    end
    set(handles.figure1,'ColorMap',gray(256))
elseif (strcmp(opt,'bw'))
    img = img_fun('im2bw',img);
    set(handles.hImg,'CData', img, 'CDataMapping','scaled');
    set(handles.figure1,'ColorMap',gray(256))
elseif (strcmp(opt,'guessType'))
    str = {'*.grd;*.nc;*.tif;*.tiff;*.jpg;*.jp2;*.png', 'Files (*.grd,*.nc,*.tif,*.tiff,*.pjg,*.jp2,*.png)'; '*.*', 'All Files (*.*)'};
    [FileName,PathName] = put_or_get_file(handles,str,'Select file','get');
    if isequal(FileName,0);     set(handles.figure1,'pointer','arrow');     return;     end             % User gave up
    drv = aux_funs('findFileType',[PathName FileName]);
    if (~isempty(drv)),         gateLoadFile(handles,drv,[PathName FileName]);  end
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
% elseif (strcmp(opt,'Paint'))
%     mpaint(handles.figure1);
end
set(handles.figure1,'pointer','arrow')

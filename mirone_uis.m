function [H1,handles,home_dir] = mirone_uis(home_dir)
% --- Creates and returns a handle to the GUI MIRONE figure. 
%#function pan resetplotview igrf_options rally_plater plate_calculator gmtedit ecran snapshot
%#function about_box parker_stuff plate_calculator euler_stuff grid_calculator tableGUI
%#function datasets_funs earthquakes manual_pole_adjust compute_euler focal_meca srtm_tool atlas
%#function image_enhance image_adjust datasets_funs write_gmt_script vitrinite telhometro mpaint
%#function imcapture filter_funs overview imageresize classificationfig tfw_funs mirone_pref
%#function griding_mir grdfilter_mir grdsample_mir grdtrend_mir grdgradient_mir ml_clip show_palette 
%#function geog_calculator color_palettes diluvio fault_models tsu_funs mk_movie_from_list
%#function mxgridtrimesh aquamoto tiles_tool empilhador grdlandmask_win grdlandmask_m escadeirar
%#function run_cmd line_operations world_is_not_round_enough cartas_militares ice_m magbarcode

	% The following test will tell us if we are using the compiled or the ML version
	try
        %dumb=evalin('base','who');
        dumb = which('mirone');
        IamCompiled = 0;
	catch
        IamCompiled = 1;
	end

	% Import icons and fetch home_dir if compiled and called by extension association
	% Here is what will happen. When called by windows extension association the 'home_dir'
	% will contain the path of the file and not of the Mirone installation, and it will fail.
	% In the 'catch' branch we check if the MIRONE_HOME environment variable exists. If yes
	% its value, the correct home_dir, is returned back to mirone.m
	try
        load ([home_dir filesep 'data' filesep 'mirone_icons.mat']);
	catch
        if (IamCompiled && ispc)
            home_dir = winqueryreg('HKEY_CURRENT_USER', 'Environment', 'MIRONE_HOME');
            load ([home_dir filesep 'data' filesep 'mirone_icons.mat']);
        end
	end

pos = [520 758 581 21];     % R13 honest figure dimension
H1 = figure('PaperUnits','centimeters',...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'ResizeFcn',@figure1_ResizeFcn,...
'KeyPressFcn',@figure1_KeyPressFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','Mirone 1.4.0b',...
'NumberTitle','off',...
'PaperPositionMode','auto',...
'PaperSize',[20.98404194812 29.67743169791],...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',pos,...
'HandleVisibility','callback',...
'Tag','figure1',...
'Visible','off');

setappdata(H1,'IAmAMirone',1)           % Use this appdata to identify Mirone figures
setappdata(H1,'PixelMode',0)            % Default

% Detect which matlab version is beeing used. For the moment I'm only interested to know if R13 or >= R14
version7 = version;
if (double(version7(1)) > 54),      version7 = 1;
else                                version7 = 0;
end

hVG = zeros(1,17);		kv = 5;		% hVG will contain the handles of "not valid grid" uis to hide when they are not usable
h_toolbar = uitoolbar('parent',H1, 'BusyAction','queue','HandleVisibility','on','Interruptible','on',...
	'Tag','FigureToolBar','Visible','on');
uipushtool('parent',h_toolbar,'Click','mirone(''TransferB_CB'',guidata(gcbo),''NewEmpty'')', ...
	'Tag','NewFigure','cdata',Mfnew_ico,'Tooltip','Open New figure');
uipushtool('parent',h_toolbar,'Click','mirone(''TransferB_CB'',guidata(gcbo),''guessType'')', ...
	'Tag','ImportKnownTypes','cdata',Mfopen_ico,'Tooltip','Load recognized file types');
hVG(1) = uipushtool('parent',h_toolbar,'Click','mirone(''FileSaveGMTgrid_CB'',guidata(gcbo))', ...
	'Tag','SaveGMTgrid','cdata',Mfsave_ico,'Tooltip','Save GMT grid');
uipushtool('parent',h_toolbar,'Click','mirone_pref(guidata(gcbo))', ...
	'Tag','Preferences','cdata',tools_ico,'Tooltip','Preferences');
uipushtool('parent',h_toolbar,'Click','mirone(''Transfer_CB'',guidata(gcbo),''print'')', ...
	'Tag','Print','cdata',Mprint_ico,'Tooltip','Print image');
uipushtool('parent',h_toolbar,'Click','mirone(''Draw_CB'',guidata(gcbo),''Text'')', ...
	'Tag','DrawText','cdata',text_ico,'Tooltip','Insert Text','Sep','on');
uipushtool('parent',h_toolbar,'Click','mirone(''DrawGeographicalCircle_CB'',guidata(gcbo))', ...
	'Tag','DrawGeogCirc','cdata',circ_ico,'Tooltip','Draw geographical circle');
uipushtool('parent',h_toolbar,'Click','mirone(''DrawLine_CB'',guidata(gcbo))', ...
	'Tag','DrawLine','cdata',Mline_ico,'Tooltip','Draw Line');
uipushtool('parent',h_toolbar,'Click','mirone(''DrawClosedPolygon_CB'',guidata(gcbo),''rectangle'')', ...
	'Tag','DrawRect','cdata',rectang_ico,'Tooltip','Draw Rectangle');
uipushtool('parent',h_toolbar,'Click','mirone(''DrawClosedPolygon_CB'',guidata(gcbo),[])', ...
	'Tag','DrawPolyg','cdata',polygon_ico,'Tooltip','Draw Closed Polygon');
uipushtool('parent',h_toolbar,'Click','mirone(''Draw_CB'',guidata(gcbo),''Vector'')', ...
	'Tag','DrawArrow','cdata',Marrow_ico,'Tooltip','Draw Arrow');
uitoggletool('parent',h_toolbar,'Click','mirone(''PanZoom_CB'',guidata(gcbo),gcbo,''zoom'')', ...
	'Tag','Zoom','cdata',zoom_ico,'Tooltip','Zooming on/off','Sep','on');
uitoggletool('parent',h_toolbar,'Click','mirone(''PanZoom_CB'',guidata(gcbo),gcbo,''pan'')', ...
	'Tag','Mao','cdata',mao,'Tooltip','Pan');
uitoggletool('parent',h_toolbar,'Click','draw_funs(gcbo,''deleteObj'')', ...
	'Tag','Tesoura','cdata',cut_ico,'Tooltip','Delete objects');
uipushtool('parent',h_toolbar,'Click','color_palettes(guidata(gcbo))', ...
	'Tag','ColorPal','cdata',color_ico,'Tooltip','Color Palettes');
hVG(2) = uipushtool('parent',h_toolbar,'Click','mirone(''ImageIlluminationModel_CB'',guidata(gcbo),''grdgradient_A'')', ...
	'Tag','Shading','cdata',shade2_ico,'Tooltip','Shaded illumination','Sep','on');
hVG(3) = uipushtool('parent',h_toolbar,'Click','mirone(''ImageAnaglyph_CB'',guidata(gcbo))', ...
	'Tag','Anaglyph','cdata',anaglyph_ico,'Tooltip','Anaglyph');
hVG(4) = uipushtool('parent',h_toolbar,'Click','mirone(''ToolsMBplaningStart_CB'',guidata(gcbo))', ...
	'Tag','MBplaning','cdata',MB_ico,'Tooltip','Multi-beam planing');
uipushtool('parent',h_toolbar,'Click','mirone(''FileSaveFleder_CB'',guidata(gcbo),''runPlanar'')', ...  
	'Tag','FlederPlanar','cdata',olho_ico,'Tooltip','Run Fleder 3D Viewer');
uipushtool('parent',h_toolbar,'Click','writekml(guidata(gcbo))', ...  
	'Tag','toGE','cdata',GE_ico,'Tooltip','See it in Google Earth');
uipushtool('parent',h_toolbar,'Click',@refresca, 'Tag','Refresh','cdata',refresh_ico,'Tooltip','Refresh','Sep','on');
uipushtool('parent',h_toolbar,'Click','grid_info(guidata(gcbo))','Tag','ImageInfo','cdata',info_ico,'Tooltip','Image info');

h_axes = axes('Parent',H1,'Units','pixels','Position',[60 0 50 10],'Tag','axes1','Visible','off');
cmenu_axes = uicontextmenu('Parent',H1);
set(h_axes, 'UIContextMenu', cmenu_axes);
uimenu(cmenu_axes, 'Label', 'Label Format -> DD.xx', 'Call', 'draw_funs([],''ChngAxLabels'',''ToDegDec'')','Tag','LabFormat');
uimenu(cmenu_axes, 'Label', 'Label Format -> DD MM', 'Call', 'draw_funs([],''ChngAxLabels'',''ToDegMin'')','Tag','LabFormat');
uimenu(cmenu_axes, 'Label', 'Label Format -> DD MM.xx', 'Call', 'draw_funs([],''ChngAxLabels'',''ToDegMinDec'')','Tag','LabFormat');
uimenu(cmenu_axes, 'Label', 'Label Format -> DD MM SS', 'Call', 'draw_funs([],''ChngAxLabels'',''ToDegMinSec'')','Tag','LabFormat');
uimenu(cmenu_axes, 'Label', 'Label Format -> DD MM SS.x', 'Call', 'draw_funs([],''ChngAxLabels'',''ToDegMinSecDec'')','Tag','LabFormat');
itemFS = uimenu(cmenu_axes, 'Label', 'Label Font Size', 'Separator','on');
uimenu(itemFS, 'Label', '7   pt', 'Call', 'set(gca, ''FontSize'', 7)');
uimenu(itemFS, 'Label', '8   pt', 'Call', 'set(gca, ''FontSize'', 8)');
uimenu(itemFS, 'Label', '9   pt', 'Call', 'set(gca, ''FontSize'', 9)');
uimenu(itemFS, 'Label', '10 pt', 'Call', 'set(gca, ''FontSize'', 10)');
uimenu(cmenu_axes, 'Label', 'Grid on/off', 'Call', 'grid', 'Separator','on');
uimenu(cmenu_axes, 'Label', 'Row-Col mode on/off', 'Tag','RCMode', 'Separator','on');
uimenu(cmenu_axes, 'Label', 'Pixel mode on/off', 'Tag','PixMode');
% Those ones are manipulated in setAxesDefCoordIn()
uimenu(cmenu_axes, 'Label', 'Load in projected coords', 'Checked','on', 'Vis','off','Tag','hAxMenuLF');
uimenu(cmenu_axes, 'Label', 'Display projected coords', 'Vis','off','Tag','hAxMenuDM');

hFL = uimenu('Parent',H1,'Label','File','Tag','File');
uimenu('Parent',hFL,'Call','mirone_pref(guidata(gcbo))','Label','Preferences');
uimenu('Parent',hFL,'Callback','mirone(''TransferB_CB'',guidata(gcbo),''NewEmpty'')','Label','New empty window','Sep','on');

h = uimenu('Parent',hFL,'Label','Background window');
uimenu('Parent',h,'Callback','mirone(''TransferB_CB'',guidata(gcbo),''BgMap'')','Label','Map');
uimenu('Parent',h,'Callback','mirone(''FileNewBgFrame_CB'', guidata(gcbo))','Label','Frame');

if (~IamCompiled)
    uimenu('Parent',hFL,'Label','Open a 2D array .mat file','Call','InOut2WS(guidata(gcbo),''loadmat'')');
end

h2 = uimenu('Parent',hFL,'Label','Open Grid/Image','Sep','on','Tag','OpenGI');
uimenu('Parent',h2,'Call','mirone(''FileOpenNewImage_CB'',guidata(gcbo));','Label','Images -> Generic Formats');
uimenu('Parent',h2,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''GMT'')','Label','GMT Grid','Sep','on');
uimenu('Parent',h2,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''Surfer'')','Label','Surfer 6/7 grid');
uimenu('Parent',h2,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''ENCOM'')','Label','Encom grid');
uimenu('Parent',h2,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''MANI'')','Label','Mani grid');
uimenu('Parent',h2,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''ArcAscii'')','Label','Arc/Info ASCII Grid','Sep','on');
uimenu('Parent',h2,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''ArcBinary'')','Label','Arc/Info Binary Grid');
uimenu('Parent',h2,'Call','mirone(''FileOpenGeoTIFF_CB'',guidata(gcbo),''bsb'')','Label','BSB Nautical Chart Format');
uimenu('Parent',h2,'Call','mirone(''FileOpen_ENVI_Erdas_CB'',guidata(gcbo),''ENVI'')','Label','ENVI Raster');
uimenu('Parent',h2,'Call','mirone(''FileOpen_ENVI_Erdas_CB'',guidata(gcbo),''ERDAS'')','Label','Erdas (.img)');
uimenu('Parent',h2,'Call','mirone(''FileOpenGeoTIFF_CB'',guidata(gcbo),''bil'')','Label','ESRI BIL');
uimenu('Parent',h2,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''GXF'')','Label','Geosoft GXF');
uimenu('Parent',h2,'Call','mirone(''FileOpenGeoTIFF_CB'',guidata(gcbo),''geotiff'')','Label','GeoTIFF');
uimenu('Parent',h2,'Call','mirone(''FileOpenGeoTIFF_CB'',guidata(gcbo),''ecw'')','Label','ECW');
uimenu('Parent',h2,'Call','mirone(''FileOpenGeoTIFF_CB'',guidata(gcbo),''sid'')','Label','MrSID');
uimenu('Parent',h2,'Call','mirone(''FileOpenGeoTIFF_CB'',guidata(gcbo),''jp2'')','Label','JPEG2000');
uimenu('Parent',h2,'Call','mirone(''FileOpenGDALmultiBand_CB'',guidata(gcbo),''ENVISAT'')','Label','ENVISAT','Sep','on');
uimenu('Parent',h2,'Call','mirone(''FileOpenGDALmultiBand_CB'',guidata(gcbo),''AVHRR'')','Label','AVHRR');
uimenu('Parent',h2,'Call','mirone(''FileOpenGeoTIFF_CB'',guidata(gcbo),''UNKNOWN'')','Label','Try Luck with GDAL');

h = uimenu('Parent',h2,'Label','Digital Elevation','Sep','on');
uimenu('Parent',h,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''DTED'')','Label','DTED');
uimenu('Parent',h,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''GTOPO30'')','Label','GTOPO30');
uimenu('Parent',h,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''GeoTiff_DEM'')','Label','GeoTIFF DEM');
uimenu('Parent',h,'Call','mirone(''FileOpenMOLA_CB'',guidata(gcbo))','Label','MOLA DEM');
uimenu('Parent',h,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''SRTM1'')','Label','SRTM 1 arcsec');
uimenu('Parent',h,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''SRTM3'')','Label','SRTM 3 arcsec');
uimenu('Parent',h,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''SRTM30'')','Label','SRTM30');
uimenu('Parent',h,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''USGS_DEM'')','Label','USGS DEM');
uimenu('Parent',h,'Call','mirone(''FileOpenDEM_CB'',guidata(gcbo),''SDTS'')','Label','USGS SDTS DEM');

uimenu('Parent',hFL,'Callback','overview(guidata(gcbo))','Label','Open Overview Tool');
uimenu('Parent',hFL,'Call','mirone(''FileOpenSession_CB'',guidata(gcbo))','Label','Open Session');

% ----------------------- Save Images section
h = uimenu('Parent',hFL,'Label','Save Image As...','Sep','on');

uimenu('Parent',h,'Call','snapshot(gcf)', 'Label','Generic Formats');
uimenu('Parent',h,'Call','mirone(''FileSaveImgGrdGdal_CB'',guidata(gcbo),''GeoTiff'',''img'')', 'Label','GeoTiff','Sep','on');
uimenu('Parent',h,'Call','mirone(''FileSaveImgGrdGdal_CB'',guidata(gcbo),''Erdas'',''img'')', 'Label','Erdas Imagine');
uimenu('Parent',h,'Call','mirone(''FileSaveImgGrdGdal_CB'',guidata(gcbo),''Envi'',''img'')', 'Label','Envi .hdr Labeled');
uimenu('Parent',h,'Call','mirone(''FileSaveImgGrdGdal_CB'',guidata(gcbo),''ESRI'',''img'')', 'Label','ESRI .hdr Labeled');
uimenu('Parent',h,'Call','mirone(''FileSaveImgGrdGdal_CB'',guidata(gcbo),''JP2K'',''img'')', 'Label','JPEG2000');
uimenu('Parent',h,'Call','tfw_funs(''write'',guidata(gcbo))', 'Label','Save .tfw file','Tag','saveTFW','Sep','on');
uimenu('Parent',hFL,'Call','snapshot(gcf,''frame'')', 'Label','Export ...','Tag','noAxes');

if (strncmp(computer,'PC',2))
    h = uimenu('Parent',hFL,'Label','Copy to Clipboard','Tag','CopyClip');
    uimenu('Parent',h,'Call','imcapture(gca,''img'');','Label','Image only')
    uimenu('Parent',h,'Call','mirone(''Transfer_CB'',guidata(gcbo),''copyclip'')', 'Label','Image and frame','Tag','noAxes')
	uimenu('Parent',h,'Call','mirone(''Transfer_CB'',guidata(gcbo),''Ctrl-c'')','Label','Copy active line', 'Accel','c', 'Sep','on')
	uimenu('Parent',h,'Call','mirone(''Transfer_CB'',guidata(gcbo),''Ctrl-v'')','Label','Paste line', 'Accel','v')
end

uimenu('Parent',hFL,'Call','mirone(''Transfer_CB'',guidata(gcbo),''KML'')', 'Label','Export to GoogleEarth','Sep','on');

% ----------------------- Save Grids section --------------------- IT ALSO STARTS the hVG(kv)
h = uimenu('Parent',hFL,'Label','Save Grid As...','Sep','on');		hVG(kv) = h;	kv = kv + 1;
uimenu('Parent',h,'Call','mirone(''FileSaveGMTgrid_CB'',guidata(gcbo))','Label','GMT grid');
uimenu('Parent',h,'Call','mirone(''FileSaveGMTgrid_CB'',guidata(gcbo),''Surfer'')','Label','Surfer 6 grid');
uimenu('Parent',h,'Call','mirone(''FileSaveENCOMgrid_CB'',guidata(gcbo))','Label','Encom grid');
uimenu('Parent',h,'Call','mirone(''FileSaveImgGrdGdal_CB'',guidata(gcbo),''GeoTiff'',''grid'')', 'Label','GeoTiff','Sep','on');
uimenu('Parent',h,'Call','mirone(''FileSaveImgGrdGdal_CB'',guidata(gcbo),''Erdas'',''grid'')', 'Label','Erdas Imagine');
uimenu('Parent',h,'Call','mirone(''FileSaveImgGrdGdal_CB'',guidata(gcbo),''Envi'',''grid'')', 'Label','Envi .hdr Labelled');
uimenu('Parent',h,'Call','mirone(''FileSaveImgGrdGdal_CB'',guidata(gcbo),''ESRI'',''grid'')', 'Label','ESRI .hdr Labelled');
uimenu('Parent',h,'Call','mirone(''FileSaveImgGrdGdal_CB'',guidata(gcbo),''JP2K'',''grid'')', 'Label','JPEG2000');

h = uimenu('Parent',hFL,'Label','Save As 3 GMT grids (R,G,B)');
uimenu('Parent',h,'Call','mirone(''File_img2GMT_RGBgrids_CB'',guidata(gcbo))','Label','Image only');
uimenu('Parent',h,'Call','mirone(''File_img2GMT_RGBgrids_CB'',guidata(gcbo),''screen'')','Label','Screen capture');

h = uimenu('Parent',hFL,'Label','Save GMT script','Sep','on');
uimenu('Parent',h,'Call','write_gmt_script(guidata(gcbo),''csh'')','Label','csh script');
uimenu('Parent',h,'Call','write_gmt_script(guidata(gcbo),''bat'')','Label','dos batch');

h = uimenu('Parent',hFL,'Label','Save As Fledermaus Objects');
hVG(kv) = uimenu('Parent',h,'Call','mirone(''FileSaveFleder_CB'',guidata(gcbo),''writeSpherical'')',...
'Label','Spherical Fledermaus Obj');		kv = kv + 1;
uimenu('Parent',h,'Call','mirone(''FileSaveFleder_CB'',guidata(gcbo),''writePlanar'')','Label','Planar Fledermaus Obj');
hVG(kv) = uimenu('Parent',h,'Call','mirone(''FileSaveFleder_CB'',guidata(gcbo),''writeAll3'')', ...
'Label','Dmagic .dtm, .geo, .shade files', 'Sep','on');		kv = kv + 1;

uimenu('Parent',hFL,'Call','mirone(''FileSaveSession_CB'',guidata(gcbo))','Label','Save Session');

if (~IamCompiled)
    uimenu('Parent',hFL,'Label','Grid/Image -> Workspace','Sep','on','Call','InOut2WS(guidata(gcbo),''direct'')');
    hVG(kv) = uimenu('Parent',hFL,'Label','Workspace -> Grid','Call','InOut2WS(guidata(gcbo),''GRID_inverse'')');	kv = kv + 1;
    uimenu('Parent',hFL,'Label','Workspace -> Image','Call','InOut2WS(guidata(gcbo),''IMG_inverse'')');
    uimenu('Parent',hFL,'Label','Clear Workspace','Call','InOut2WS(guidata(gcbo),''clear'')');
end

h = uimenu('Parent',hFL,'Label','Recent Files','Tag','RecentFiles','Sep','on');
for (i=1:6),    uimenu('Parent',h,'Vis','off','Tag','RecentF');   end

uimenu('Parent',hFL,'Call','print -dsetup','Label','Print Setup','Sep','on');
uimenu('Parent',hFL,'Call','mirone(''Transfer_CB'',guidata(gcbo),''print'')','Label','Print...');

% --------------------------- IMAGE MENU ------------------------------------
hIM = uimenu('Parent',H1,'Label','Image','Tag','Image');

h = uimenu('Parent',hIM,'Label','Color Palettes');
uimenu('Parent',h,'Call','color_palettes(guidata(gcbo))','Label','Change Palette');
uimenu('Parent',h,'Call','show_palette(guidata(gcbo),''At'')','Label','Show At side','Tag','PalAt','Sep','on');
uimenu('Parent',h,'Call','show_palette(guidata(gcbo),''In'')','Label','Show Inside','Tag','PalIn');
uimenu('Parent',h,'Call','show_palette(guidata(gcbo),''Float'')','Label','Show Floating');
uimenu('Parent',hIM,'Call','mirone(''ImageCrop_CB'',guidata(gcbo),[])','Label','Crop (interactive)','Sep','on');

h = uimenu('Parent',hIM,'Label','Flip');
uimenu('Parent',h,'Call','mirone(''Transfer_CB'',guidata(gcbo),''flipUD'')','Label','Flip Up-Down');
uimenu('Parent',h,'Call','mirone(''Transfer_CB'',guidata(gcbo),''flipLR'')','Label','Flip Left-Right');

uimenu('Parent',hIM,'Call','mirone(''ImageResetOrigImg_CB'',guidata(gcbo))','Label','Restore Original Image','Tag','ImRestore');
uimenu('Parent',hIM,'Call','mirone(''ImageHistEqualize_CB'',guidata(gcbo),gcbo)','Label','Histogram Equalization (image)','Tag','ImgHist');

hVG(kv) = uimenu('Parent',hIM,'Call','mirone(''ImageHistEqualizeGrid_CB'',guidata(gcbo),gcbo)',...
'Label','Histogram Equalization (grid)','Tag','ImgHistGrd');		kv = kv + 1;

h = uimenu('Parent',hIM,'Label','Illuminate','Tag','Illuminate');	hVG(kv) = h;	kv = kv + 1;
uimenu('Parent',h,'Call','mirone(''ImageIlluminationModel_CB'',guidata(gcbo),''grdgradient_A'')','Label','GMT grdgradient');
uimenu('Parent',h,'Call','mirone(''ImageIlluminationModel_CB'',guidata(gcbo),''falseColor'')','Label','False color');

hVG(kv) = uimenu('Parent',hIM,'Call','mirone(''ImageAnaglyph_CB'',guidata(gcbo))','Label','Anaglyph'); kv = kv + 1;
uimenu('Parent',hIM,'Call','mirone(''ImageDrape_CB'',guidata(gcbo))','Label','Drape','Tag','ImageDrape','Enable','off');
uimenu('Parent',hIM,'Call','mirone(''ImageMapLimits_CB'',guidata(gcbo))','Label','Map Limits');

h = uimenu('Parent',hIM,'Label','Edge detect (Canny)','Sep','on');
uimenu('Parent',h,'Call','mirone(''ImageEdgeDetect_CB'',guidata(gcbo),''Vec'')','Label','Vector');
uimenu('Parent',h,'Call','mirone(''ImageEdgeDetect_CB'',guidata(gcbo),''Ras'')','Label','Raster');

h = uimenu('Parent',hIM,'Label','Edge detect (SUSAN)');
uimenu('Parent',h,'Call','mirone(''ImageEdgeDetect_CB'',guidata(gcbo),''SUSvec'')','Label','Vector');
uimenu('Parent',h,'Call','mirone(''ImageEdgeDetect_CB'',guidata(gcbo),''SUSras'')','Label','Raster');

h = uimenu('Parent',hIM,'Label','Features detection');
uimenu('Parent',h,'Call','mirone(''ImageEdgeDetect_CB'',guidata(gcbo),''Lines'')','Label','Lines');
uimenu('Parent',h,'Call','mirone(''ImageEdgeDetect_CB'',guidata(gcbo),''Circles'')','Label','Circles');
uimenu('Parent',h,'Call','mirone(''Transfer_CB'',guidata(gcbo),''Corners'')','Label','Good features to track');

uimenu('Parent',hIM,'Call','mirone(''Transfer_CB'',guidata(gcbo),''Shape'')','Label','Shape detector')
uimenu('Parent',hIM,'Call','mpaint(gcf)','Label','Paint Brush');
uimenu('Parent',hIM,'Call','mirone(''ImageSegment_CB'',guidata(gcbo))','Label','Image segmentation');
uimenu('Parent',hIM,'Call','classificationfig(gcf);','Label','K-means classification');
uimenu('Parent',hIM,'Call','imageresize(gcf)','Label','Image resize','Sep','on');
uimenu('Parent',hIM,'Call','mirone(''RotateTool_CB'',guidata(gcbo),''image'')','Label','Image rotation');

h = uimenu('Parent',hIM,'Label','Image mode');
uimenu('Parent',h,'Call','mirone(''Transfer_CB'',guidata(gcbo),''toRGB'')','Label','RGB truecolor','Tag','ImModRGB');
uimenu('Parent',h,'Call','mirone(''Transfer_CB'',guidata(gcbo),''8-bit'')','Label','8-bit color','Tag','ImMod8cor');
uimenu('Parent',h,'Call','mirone(''Transfer_CB'',guidata(gcbo),''gray'')','Label','Gray scale','Tag','ImMod8gray');
uimenu('Parent',h,'Call','mirone(''Transfer_CB'',guidata(gcbo),''bw'')','Label','Black and White','Tag','ImModBW');
uimenu('Parent',h,'Call','mirone(''ImageResetOrigImg_CB'',guidata(gcbo))','Label','Original Image','Sep','on');

% ------------ Image filters _______ TO BE CONTINUED
h = uimenu('Parent',hIM,'Label','Filters','Sep','on');
uimenu('Parent',h,'Call','filter_funs(guidata(gcbo),''SUSAN'');','Label','Smooth (SUSAN)');
uimenu('Parent',h,'Call','filter_funs(guidata(gcbo),''Median'');','Label','Median (3x3)');
uimenu('Parent',h,'Call','filter_funs(guidata(gcbo),''Adaptive'');','Label','Adaptive Median (7x7)');
uimenu('Parent',h,'Call','filter_funs(guidata(gcbo),''STD'');','Label','STD (3x3)');
uimenu('Parent',h,'Call','filter_funs(guidata(gcbo),''Min'');','Label','Min (3x3)');
uimenu('Parent',h,'Call','filter_funs(guidata(gcbo),''Max'');','Label','Max (3x3)');
uimenu('Parent',h,'Call','filter_funs(guidata(gcbo),''range'');','Label','Range (3x3)');
uimenu('Parent',hIM,'Call','mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''Allopts'')','Label','FFT Spectrum');

uimenu('Parent',hIM,'Call','mirone(''DigitalFilt_CB'',guidata(gcbo),''image'')','Label','Digital Filtering Tool','Sep','on');
uimenu('Parent',hIM,'Call','image_enhance(gcf)','Label','Image Enhance (1 - Indexed and RGB)');
uimenu('Parent',hIM,'Call','image_adjust(gcf)','Label', 'Image Enhance (2 - Indexed only)');
uimenu('Parent',hIM,'Call','ice_m(gcf,''space'',''rgb'')','Label','Image Color Editor (Indexed and RGB)');
uimenu('Parent',hIM,'Call','aux_funs(''togCheck'',gcbo)','Label','Activate Image-to-Image/Map GCP Tool','Tag','GCPtool','Sep','on');
uimenu('Parent',hIM,'Call','mirone(''DrawLine_CB'',guidata(gcbo),''GCPpline'')','Label','Register Image (Draw GCP points)');
uimenu('Parent',hIM,'Call','mirone(''DrawLine_CB'',guidata(gcbo),''GCPmemory'')',...
	'Tag','GCPmemory','Label','Register Image (Plot in memory GCPs)','Visible','off');  % To GDAL imported file with GCPs
uimenu('Parent',hIM,'Callback','bands_list(gcf)','Label','Load Bands','Sep','on');
%uimenu('Parent',hIM,'Callback','grid_calculator(gcf)','Label','Bands Arithmetic');
uimenu('Parent',hIM,'Call','mk_movie_from_list(guidata(gcbo))', 'Label','Make movie from image list','Sep','on');

% --------------------------- TOOLS MENU ------------------------------------
hTL = uimenu('Parent',H1,'Label','Tools','Tag','Tools');
h = uimenu('Parent',hTL,'Label','Extract Profile');
uimenu('Parent',h,'Call','mirone(''ExtractProfile_CB'',guidata(gcbo))','Label','Static');
uimenu('Parent',h,'Call','mirone(''ExtractProfile_CB'',guidata(gcbo),''dynamic'')','Label','Dynamic');

h = uimenu('Parent',hTL,'Label','Measure','Sep','on');
uimenu('Parent',h,'Call','mirone(''ToolsMeasure_CB'',guidata(gcbo),''LLength'')','Label','Distance','Tag','ToolsMeasureDist');
uimenu('Parent',h,'Call','mirone(''ToolsMeasure_CB'',guidata(gcbo),''Azim'')','Label','Azimuth');
uimenu('Parent',h,'Call','mirone(''ToolsMeasure_CB'',guidata(gcbo),''Area'')','Label','Area');
uimenu('Parent',h,'Call','mirone(''ToolsMeasureAreaPerCor_CB'',guidata(gcbo))','Label','Area per color');

h = uimenu('Parent',hTL,'Label','Multi-beam planing','Tag','MBplan','Sep','on');		hVG(kv) = h;	kv = kv + 1;
uimenu('Parent',h,'Call','mirone(''ToolsMBplaningStart_CB'',guidata(gcbo))','Label','Start track');
uimenu('Parent',h,'Call','mirone(''ToolsMBplaningImport_CB'',guidata(gcbo))','Label','Import track');

uimenu('Parent',hTL,'Call','mirone(''DrawClosedPolygon_CB'',guidata(gcbo),''from_ROI'')','Label','Region-Of-Interest','Sep','on');
uimenu('Parent',hTL,'Call','mirone(''TransferB_CB'',guidata(gcbo),''isGMT'')','Label','Does GMT know this file?','Sep','on');
uimenu('Parent',hTL,'Call','mirone(''TransferB_CB'',guidata(gcbo),''isGDAL'')','Label','Does GDAL know this file?');
uimenu('Parent',hTL,'Call','ecran','Label','X,Y grapher','Sep','on');
hVG(kv) = uimenu('Parent',hTL,'Call','gmtedit','Label','gmtedit');	kv = kv + 1;
uimenu('Parent',hTL,'Call','rally_plater','Label','Rally Plater');
uimenu('Parent',hTL,'Label','entry_vtr','Sep','on');
uimenu('Parent',hTL,'Call','aquamoto(guidata(gcbo))','Label','Aquamoto Viewer','Sep','on');
uimenu('Parent',hTL,'Call','empilhador(guidata(gcbo))','Label','Empilhador');
uimenu('Parent',hTL,'Call','tiles_tool(guidata(gcbo))','Label','Tiling Tool','Sep','on');
hVG(kv) = uimenu('Parent',hTL,'Call','diluvio(guidata(gcbo))','Label','Noe Diluge','Sep','on');		kv = kv + 1;
uimenu('Parent',hTL,'Call','cartas_militares(guidata(gcbo))','Label','Cartas Militares','Sep','on')
uimenu('Parent',hTL,'Call','world_is_not_round_enough(guidata(gcbo))','Label','World is not (round) enough','Sep','on');
uimenu('Parent',hTL,'Call','run_cmd(guidata(gcbo))','Label','Run ML Command','Sep','on');
uimenu('Parent',hTL,'Call','line_operations(guidata(gcbo))','Label','Line Operations','Tag','lineOP');
% uimenu('Parent',hTL,'Call','shape_tool(gcf)','Label','Limiares','Sep','on');
% uimenu('Parent',hTL,'Call','autofaults(guidata(gcbo))','Label','Auto falhas','Sep','on');

% --------------------------- DRAW MENU ------------------------------------
hDR = uimenu('Parent',H1,'Label','Draw','Tag','Draw');
uimenu('Parent',hDR,'Call','mirone(''DrawLine_CB'',guidata(gcbo))','Label','Draw line','Tag','ctrLine', 'Accelerator','l');
uimenu('Parent',hDR,'Call','mirone(''DrawLine_CB'',guidata(gcbo),''spline'')','Label','Draw interpolating spline');
uimenu('Parent',hDR,'Call','mirone(''DrawLine_CB'',guidata(gcbo),''freehand'')','Label','Freehand draw');
uimenu('Parent',hDR,'Call','mirone(''DrawClosedPolygon_CB'',guidata(gcbo),[])','Label','Draw closed polygon');
uimenu('Parent',hDR,'Call','mirone(''DrawClosedPolygon_CB'',guidata(gcbo),''EulerTrapezium'')','Label','Draw Euler trapezium');
uimenu('Parent',hDR,'Call','mirone(''DrawImportLine_CB'',guidata(gcbo),''AsLine'')','Label','Import line','Sep','on');
uimenu('Parent',hDR,'Call','mirone(''DrawImportLine_CB'',guidata(gcbo),''AsPoint'')','Label','Import points');
uimenu('Parent',hDR,'Call','mirone(''Transfer_CB'',guidata(gcbo),''scatter'')','Label','Import scaled symbols');
uimenu('Parent',hDR,'Call','mirone(''DrawImportText_CB'',guidata(gcbo))','Label','Import text');
uimenu('Parent',hDR,'Call','mirone(''DrawImportShape_CB'',guidata(gcbo))','Label','Import shape file');

h = uimenu('Parent',hDR,'Label','Draw circle','Sep','on');
uimenu('Parent',h,'Call','mirone(''DrawGeographicalCircle_CB'',guidata(gcbo))','Label','Geographical circle');
uimenu('Parent',h,'Call','mirone(''DrawGeographicalCircle_CB'',guidata(gcbo),''gcirc'')','Label','Great circle arc');
uimenu('Parent',h,'Call','mirone(''DrawEulerPoleCircle_CB'',guidata(gcbo))','Label','Circle about an Euler pole');
uimenu('Parent',h,'Call','mirone(''DrawGeographicalCircle_CB'',guidata(gcbo),''cartCirc'')','Label','Cartesian circle');
uimenu('Parent',hDR,'Call','mirone(''DrawClosedPolygon_CB'',guidata(gcbo),''rectangle'')','Label','Draw rectangle');
uimenu('Parent',hDR,'Call','mirone(''Draw_CB'',guidata(gcbo),''Vector'')','Label','Draw vector');

h = uimenu('Parent',hDR,'Label','Draw symbol');
uimenu('Parent',h,'Call','mirone(''Draw_CB'',guidata(gcbo),''Symbol'',''o'')','Label','Circle');
uimenu('Parent',h,'Call','mirone(''Draw_CB'',guidata(gcbo),''Symbol'',''s'')','Label','Square');
uimenu('Parent',h,'Call','mirone(''Draw_CB'',guidata(gcbo),''Symbol'',''^'')','Label','Triangle');
uimenu('Parent',h,'Call','mirone(''Draw_CB'',guidata(gcbo),''Symbol'',''p'')','Label','Star');
uimenu('Parent',h,'Call','mirone(''Draw_CB'',guidata(gcbo),''Symbol'',''x'')','Label','Cross');

uimenu('Parent',hDR,'Call','mirone(''Draw_CB'',guidata(gcbo),''Text'')','Label','Insert text');
hVG(kv) = uimenu('Parent',hDR,'Callback','mirone(''DrawContours_CB'',guidata(gcbo))',...
'Label','Contours (automatic)','Sep','on','Tag','Contours_a');		kv = kv + 1;

hVG(kv) = uimenu('Parent',hDR,...
'Call','mirone(''DrawContours_CB'',guidata(gcbo),''gui'')','Label','Contours','Tag','Contours_i');	kv = kv + 1;

% --------------------------- DATASETS MENU ------------------------------------
hDS = uimenu('Parent',H1,'Label','Datasets','Tag','Datasets');
h = uimenu('Parent',hDS,'Label','Draw coastline','Tag','VoidDatasetsCoastLine');

uimenu('Parent',h,'Call','datasets_funs(''CoastLines'',guidata(gcbo),''c'')','Label','Crude resolution','Tag','CoastLineCrude');
uimenu('Parent',h,'Call','datasets_funs(''CoastLines'',guidata(gcbo),''l'')','Label','Low resolution','Tag','CoastLineLow');
uimenu('Parent',h,'Call','datasets_funs(''CoastLines'',guidata(gcbo),''i'')','Label','Intermediate resolution','Tag','CoastLineInterm');
uimenu('Parent',h,'Call','datasets_funs(''CoastLines'',guidata(gcbo),''h'')','Label','High resolution','Tag','CoastLineHigh');
uimenu('Parent',h,'Call','datasets_funs(''CoastLines'',guidata(gcbo),''f'')','Label','Full resolution','Tag','CoastLineFull');

h2 = uimenu('Parent',hDS,'Label','Draw political boundaries','Tag','VoidDatasetsPB');
h = uimenu('Parent',h2,'Label','National boundaries');

uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''1'',''c'')','Label','Crude resolution','Tag','PBCrude');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''1'',''l'')','Label','Low resolution','Tag','PBLow');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''1'',''i'')','Label','Intermediate resolution','Tag','PBInterm');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''1'',''h'')','Label','High resolution','Tag','PBHigh');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''1'',''f'')','Label','Full resolution','Tag','PBFull');

h = uimenu('Parent',h2,'Label','State boundaries (US)');

uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''2'',''c'')','Label','Crude resolution','Tag','PBCrude');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''2'',''l'')','Label','Low resolution','Tag','PBLow');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''2'',''i'')','Label','Intermediate resolution','Tag','PBInterm');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''2'',''h'')','Label','High resolution','Tag','PBHigh');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''2'',''f'')','Label','Full resolution','Tag','PBFull');

h = uimenu('Parent',h2,'Label','Marine boundaries');

uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''3'',''c'')','Label','Crude resolution','Tag','PBCrude');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''3'',''l'')','Label','Low resolution','Tag','PBLow');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''3'',''i'')','Label','Intermediate resolution','Tag','PBInterm');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''3'',''h'')','Label','High resolution','Tag','PBHigh');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''3'',''f'')','Label','Full resolution','Tag','PBFull');

h = uimenu('Parent',h2,'Label','All boundaries');

uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''a'',''c'')','Label','Crude resolution','Tag','PBCrude');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''a'',''l'')','Label','Low resolution','Tag','PBLow');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''a'',''i'')','Label','Intermediate resolution','Tag','PBInterm');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''a'',''h'')','Label','High resolution','Tag','PBHigh');
uimenu('Parent',h,'Call','datasets_funs(''Political'',guidata(gcbo),''a'',''f'')','Label','Full resolution','Tag','PBFull');

h2 = uimenu('Parent',hDS,'Label','Draw rivers','Tag','VoidDatasetsRivers');
h = uimenu('Parent',h2,'Label','Permanent major rivers');

uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''1'',''c'')','Label','Crude resolution','Tag','RiversCrude');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''1'',''l'')','Label','Low resolution','Tag','RiversLow');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''1'',''i'')','Label','Intermediate resolution','Tag','RiversInterm');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''1'',''h'')','Label','High resolution','Tag','RiversHigh');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''1'',''f'')','Label','Full resolution','Tag','RiversFull');

h = uimenu('Parent',h2,'Label','Additional major rivers');

uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''2'',''c'')','Label','Crude resolution','Tag','RiversCrude');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''2'',''l'')','Label','Low resolution','Tag','RiversLow');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''2'',''i'')','Label','Intermediate resolution','Tag','RiversInterm');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''2'',''h'')','Label','High resolution','Tag','RiversHigh');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''2'',''f'')','Label','Full resolution','Tag','RiversFull');

h = uimenu('Parent',h2,'Label','Additional rivers');

uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''3'',''c'')','Label','Crude resolution','Tag','RiversCrude');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''3'',''l'')','Label','Low resolution','Tag','RiversLow');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''3'',''i'')','Label','Intermediate resolution','Tag','RiversInterm');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''3'',''h'')','Label','High resolution','Tag','RiversHigh');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''3'',''f'')','Label','Full resolution','Tag','RiversFull');

h = uimenu('Parent',h2,'Label','Intermittent rivers - major');

uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''5'',''c'')','Label','Crude resolution','Tag','RiversCrude');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''5'',''l'')','Label','Low resolution','Tag','RiversLow');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''5'',''i'')','Label','Intermediate resolution','Tag','RiversInterm');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''5'',''h'')','Label','High resolution','Tag','RiversHigh');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''5'',''f'')','Label','Full resolution','Tag','RiversFull');

h = uimenu('Parent',h2,'Label','Intermittent rivers - additional');

uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''6'',''c'')','Label','Crude resolution','Tag','RiversCrude');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''6'',''l'')','Label','Low resolution','Tag','RiversLow');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''6'',''i'')','Label','Intermediate resolution','Tag','RiversInterm');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''6'',''h'')','Label','High resolution','Tag','RiversHigh');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''6'',''f'')','Label','Full resolution','Tag','RiversFull');

h = uimenu('Parent',h2,'Label','Intermittent rivers - minor');

uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''7'',''c'')','Label','Crude resolution','Tag','RiversCrude');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''7'',''l'')','Label','Low resolution','Tag','RiversLow');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''7'',''i'')','Label','Intermediate resolution','Tag','RiversInterm');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''7'',''h'')','Label','High resolution','Tag','RiversHigh');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''7'',''f'')','Label','Full resolution','Tag','RiversFull');

h = uimenu('Parent',h2,'Label','All rivers and canals');

uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''a'',''c'')','Label','Crude resolution','Tag','RiversCrude');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''a'',''l'')','Label','Low resolution','Tag','RiversLow');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''a'',''i'')','Label','Intermediate resolution','Tag','RiversInterm');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''a'',''h'')','Label','High resolution','Tag','RiversHigh');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''a'',''f'')','Label','Full resolution','Tag','RiversFull');

h = uimenu('Parent',h2,'Label','All permanent rivers');

uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''r'',''c'')','Label','Crude resolution','Tag','RiversCrude');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''r'',''l'')','Label','Low resolution','Tag','RiversLow');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''r'',''i'')','Label','Intermediate resolution','Tag','RiversInterm');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''r'',''h'')','Label','High resolution','Tag','RiversHigh');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''r'',''f'')','Label','Full resolution','Tag','RiversFull');

h = uimenu('Parent',h2,'Label','All intermittent rivers');

uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''i'',''c'')','Label','Crude resolution','Tag','RiversCrude');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''i'',''l'')','Label','Low resolution','Tag','RiversLow');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''i'',''i'')','Label','Intermediate resolution','Tag','RiversInterm');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''i'',''h'')','Label','High resolution','Tag','RiversHigh');
uimenu('Parent',h,'Call','datasets_funs(''Rivers'',guidata(gcbo),''i'',''f'')','Label','Full resolution','Tag','RiversFull');

uimenu('Parent',hDS,'Call','earthquakes(gcf);','Label','Global seismicity','Sep','on');
uimenu('Parent',hDS,'Call','datasets_funs(''Hotspots'',guidata(gcbo))','Label','Hotspot locations');
uimenu('Parent',hDS,'Call','datasets_funs(''Isochrons'',guidata(gcbo))','Label','Magnetic isochrons');
uimenu('Parent',hDS,'Call','datasets_funs(''Volcanoes'',guidata(gcbo))','Label','Volcanoes');
uimenu('Parent',hDS,'Call','datasets_funs(''Tides'',guidata(gcbo))','Label','Tide Stations');
uimenu('Parent',hDS,'Call','datasets_funs(''Plate'',guidata(gcbo))','Label','Plate boundaries');

h = uimenu('Parent',hDS,'Label','Cities');
uimenu('Parent',h,'Call','datasets_funs(''Cities'',guidata(gcbo),''major'')','Label','Major cities');
uimenu('Parent',h,'Call','datasets_funs(''Cities'',guidata(gcbo),''other'')','Label','Other cities');

h = uimenu('Parent',hDS,'Label','ODP/DSDP sites');
uimenu('Parent',h,'Call','datasets_funs(''ODP'',guidata(gcbo),''ODP'')','Label','ODP');
uimenu('Parent',h,'Call','datasets_funs(''ODP'',guidata(gcbo),''DSDP'')','Label','DSDP');
uimenu('Parent',h,'Call','datasets_funs(''ODP'',guidata(gcbo),''ODP_DSDP'')','Label','ODP and DSDP');

uimenu('Parent',hDS,'Call','magbarcode','Label','Mgnetic Bar Code');
uimenu('Parent',hDS,'Call','atlas(guidata(gcbo))','Label','Atlas','Tag','Atlas','Sep','on');
uimenu('Parent',hDS,'Call','datasets_funs(''Isochrons'',guidata(gcbo),[])','Label','External db','Sep','on');
% uimenu('Parent',hDS,'Call','datasets_funs(''GTiles'', guidata(gcbo))','Label','GTiles Map','Sep','on');

% --------------------------- GEOPHYSICS MENU ------------------------------------
hGP = uimenu('Parent',H1,'Label','Geophysics','Tag','Geophysics');
h = uimenu('Parent',hGP,'Label','Elastic deformation');
uimenu('Parent',h,'Call','mirone(''DrawLine_CB'',guidata(gcbo),''FaultTrace'')','Label','Draw Fault');
uimenu('Parent',h,'Call','mirone(''DrawImportLine_CB'',guidata(gcbo),''FaultTrace'')','Label','Import Trace Fault');
uimenu('Parent',h,'Call','fault_models(guidata(gcbo))', 'Label','Import Model Slip','Sep', 'on');

h = uimenu('Parent',hGP,'Label','Tsunami Travel Time','Sep','on');
uimenu('Parent',h,'Call','tsu_funs(''TTT'',guidata(gcbo))','Label','Plot source');
uimenu('Parent',h,'Call','tsu_funs(''TTT'',guidata(gcbo),''load'')','Label','Import maregrphs and time');
uimenu('Parent',h,'Call','tsu_funs(''TTT'',guidata(gcbo),''compute'')','Label','Compute','Sep','on');

h = uimenu('Parent',hGP,'Label','Swan');
uimenu('Parent',h,'Call','tsu_funs(''SwanCompute'',guidata(gcbo))','Label','Compute');
uimenu('Parent',h,'Call','mirone(''DrawImportLine_CB'',guidata(gcbo),''AsMaregraph'')','Label','Import Stations','Sep','on');
uimenu('Parent',h,'Call','mirone(''GeophysicsSwanPlotStations_CB'',guidata(gcbo))','Label','Plot Stations');
uimenu('Parent',h,'Call','tsu_funs(''SwanGridBorder'',guidata(gcbo))','Label','Stations on grid borders');

h = uimenu('Parent',hGP,'Label','Tsun2');
uimenu('Parent',h,'Call','tsu_funs(''Tsun2'',guidata(gcbo),''compute'')','Label','Compute');
uimenu('Parent',h,'Call','tsu_funs(''Tsun2'',guidata(gcbo),''write_params'')', 'Label','Write params file');

uimenu('Parent',hGP,'Call','igrf_options(guidata(gcbo))','Label','IGRF calculator','Sep','on');
uimenu('Parent',hGP,'Call','parker_stuff(''parker_direct'',gcf)','Label','Parker Direct');
uimenu('Parent',hGP,'Call','parker_stuff(''parker_inverse'',gcf)','Label','Parker Inversion');
uimenu('Parent',hGP,'Call','parker_stuff(''redPole'',gcf)','Label','Reduction to the Pole');
uimenu('Parent',hGP,'Call','plate_calculator','Label','Plate calculator','Sep','on');
uimenu('Parent',hGP,'Call','geog_calculator(guidata(gcbo))','Label','Geographic calculator');
uimenu('Parent',hGP,'Call','euler_stuff(gcf)','Label','Euler rotations','Sep','on');
uimenu('Parent',hGP,'Call','compute_euler(gcf)','Label','Compute Euler pole');
uimenu('Parent',hGP,'Call','manual_pole_adjust(gcf)','Label','Manual adjust Euler pole');

uimenu('Parent',hGP,'Label','entry_tl');

h = uimenu('Parent',hGP,'Label','Seismicity','Sep','on');
uimenu('Parent',h,'Callback','earthquakes(gcf,''external'');','Label','Epicenters');
uimenu('Parent',h,'Callback','focal_meca(gcf)','Label','Focal mechanisms');

h = uimenu('Parent',hGP,'Label','Import *.gmt files(s)','Sep','on');
uimenu('Parent',h,'Call','mirone(''GeophysicsImportGmtFile_CB'',guidata(gcbo),''single'')','Label','Single *.gmt file');
uimenu('Parent',h,'Call','mirone(''GeophysicsImportGmtFile_CB'',guidata(gcbo),''list'')','Label','List of files');

% --------------------------- GRID TOOLS MENU ------------------------------------
hGT = uimenu('Parent',H1,'Label','Grid Tools','Tag','GridTools');		hVG(kv) = hGT;

sep = 'off';
if (~IamCompiled)
    uimenu('Parent',hGT,'Call','grid_calculator(gcf)','Label','grid calculator');	sep = 'on';
end

uimenu('Parent',hGT,'Call','grdfilter_mir(guidata(gcbo))','Label','grdfilter','Sep',sep);
uimenu('Parent',hGT,'Call','grdgradient_mir(guidata(gcbo))','Label','grdgradient');
uimenu('Parent',hGT,'Call','grdsample_mir(guidata(gcbo))','Label','grdsample');
uimenu('Parent',hGT,'Call','grdtrend_mir(guidata(gcbo))','Label','grdtrend');
uimenu('Parent',hGT,'Call','grdlandmask_win(guidata(gcbo))','Label','grdlandmask');
uimenu('Parent',hGT,'Call','geog_calculator(guidata(gcbo),''onlyGrid'')','Label','grdproject');
uimenu('Parent',hGT,'Call','mirone(''DigitalFilt_CB'',guidata(gcbo),''grid'')','Label','Digital filtering Tool','Sep','on');
uimenu('Parent',hGT,'Call','ml_clip(guidata(gcbo))','Label','Clip Grid');
uimenu('Parent',hGT,'Call','mirone(''ImageCrop_CB'',guidata(gcbo),[],''CropaGrid'')','Label','Crop Grid');
uimenu('Parent',hGT,'Call','mirone(''RotateTool_CB'',guidata(gcbo),''grid'')','Label','Rotate Grid');
uimenu('Parent',hGT,'Call','mirone(''GridToolsHistogram_CB'',guidata(gcbo))','Label','Histogram');
uimenu('Parent',hGT,'Call','mirone(''GridToolsGridMask_CB'',guidata(gcbo))','Label','Write Mask', 'Tag','haveNaNs');
uimenu('Parent',hGT,'Call','inpaint_nans(guidata(gcbo))','Label','Inpaint NaNs', 'Tag','haveNaNs');

h = uimenu('Parent',hGT,'Label','Hammer grid','Sep','on');
uimenu('Parent',h,'Call','escadeirar(guidata(gcbo))','Label','Rice-field Grid');
uimenu('Parent',h,'Call','mirone(''TransferB_CB'',guidata(gcbo),''scale'')','Label','Rescale');
uimenu('Parent',h,'Call','mirone(''GridToolsPadd2Const_CB'',guidata(gcbo))','Label','Padd to zero');

uimenu('Parent',hGT,'Call','mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''Allopts'')','Label','FFT tool','Sep','on');

h = uimenu('Parent',hGT,'Label','Spectrum');
uimenu('Parent',h,'Call','mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''Amplitude'')','Label','Amplitude spectrum');
uimenu('Parent',h,'Call','mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''Power'')','Label','Power spectrum');
uimenu('Parent',h,'Call','mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''Autocorr'')','Label','Autocorrelation');

uimenu('Parent',hGT,'Call','mirone(''GridToolsSmooth_CB'',guidata(gcbo))','Label','Spline Smooth','Sep','on');

h = uimenu('Parent',hGT,'Label','SDG');
uimenu('Parent',h,'Call','mirone(''GridToolsSDG_CB'',guidata(gcbo),''positive'')','Label','Positive');
uimenu('Parent',h,'Call','mirone(''GridToolsSDG_CB'',guidata(gcbo),''negative'')','Label','Negative');
uimenu('Parent',h,'Call','mirone(''GridToolsSDG_CB'',guidata(gcbo),[])','Label','Both');

h2 = uimenu('Parent',hGT,'Label','Terrain Modeling','Sep','on');
h = uimenu('Parent',h2,'Label','Slope');
uimenu('Parent',h,'Call','mirone(''GridToolsSlope_CB'',guidata(gcbo),''degrees'')','Label','In degrees');
uimenu('Parent',h,'Call','mirone(''GridToolsSlope_CB'',guidata(gcbo),''percent'')','Label','In percentage');

uimenu('Parent',h2,'Call','mirone(''GridToolsSlope_CB'',guidata(gcbo),''aspect'')','Label','Aspect');
h = uimenu('Parent',h2,'Label','Directional derivative');
uimenu('Parent',h,'Call','mirone(''GridToolsDirDerive_CB'',guidata(gcbo),''first'')','Label','First derivative');
uimenu('Parent',h,'Call','mirone(''GridToolsDirDerive_CB'',guidata(gcbo),''second'')','Label','Second derivative');

h = uimenu('Parent',hGT,'Label','Interpolate','Sep','on');
uimenu('Parent',h,'Callback','griding_mir(gcf,''surface'');', 'Label','Minimum curvature');
uimenu('Parent',h,'Callback','griding_mir(gcf,''nearneighbor'');', 'Label','Near neighbor');

h2 = uimenu('Parent',hGT,'Label','SRTM tools','Sep','on');

h = uimenu('Parent',h2,'Label','SRTM mosaic');
uimenu('Parent',h,'Callback','srtm_tool','Label','SRTM 3sec');
uimenu('Parent',h,'Callback','srtm_tool(''srtm1'')','Label','SRTM 1sec');
uimenu('Parent',h,'Callback','srtm_tool(''srtm30'')','Label','SRTM30');

uimenu('Parent',h2,'Call','mirone(''GridToolsFindHoles_CB'',guidata(gcbo))','Label','Find holes');
uimenu('Parent',h2,'Call','mirone(''GridToolsSaveAsSRTM_CB'',guidata(gcbo))','Label','Save as SRTM');

uimenu('Parent',hGT,'Call','mirone(''ImageEdgeDetect_CB'',guidata(gcbo),''ppa'')',...
'Label','Extract ridges/valleys','Sep','on');

% --------------------------- PROJECTIONS MENU -----------------------------
h = uimenu('Parent',H1,'Label','Projections','Tag','Projections');
projection_menu(H1, h, home_dir);
uimenu('Parent',h,'Label','-- REPROJECT --','HitTest','off','Sep','on');
uimenu('Parent',h,'Call','gdal_project(guidata(gcbo))','Label','GDAL project');
uimenu('Parent',h,'Call','geog_calculator(guidata(gcbo),''onlyGrid'')','Label','GMT project');

% --------------------------- HELP MENU ------------------------------------
h = uimenu('Parent',H1,'Label','Help','Tag','Help');
uimenu('Parent',h,'Call','aux_funs(''help'',guidata(gcbo))','Label','Mirone Help (v1.3.0)');
uimenu('Parent',h, 'Call', @showGDALdrivers,'Label','List GDAL formats','Sep','on')
uimenu('Parent',h, 'Call','about_box(guidata(gcbo),''Mirone Last modified at 14 Nov 2008'',''1.4.0b'')','Label','About','Sep','on');

% --------------------------- Build HANDLES and finish things here
	handles = guihandles(H1);
	handles.version7 = version7;			% If == 1 => R14 or latter
	handles.IamCompiled = IamCompiled;		% If == 1 than we know that we are dealing with a compiled (V3) version
	if (version7),  set(H1,'Pos',[pos(1:3) 1]);    end     % Adjust for the > R13 bugginess
	handles.RecentF = handles.RecentF(6:-1:1);  % Inverse creation order so that newest files show on top of the list
	handles.noVGlist = hVG;					% List of ui handles that will show when "not valid grid"
	movegui(H1,'north');					% Reposition the window on screen
	set(0,'CurrentFigure',H1)				% Due to a R2006a incredible BUG
	set(H1,'Visible','on');

% --------------------------------------------------------------------------------------------------
% We need this function also when the pixval_stsbar got stucked
function refresca(obj,eventdata)
	hFig = get(0,'CurrentFigure');
    set(hFig,'Pointer','arrow');
	%set(hFig,'Renderer','painters', 'RendererMode','auto')
	refresh(hFig)

% --------------------------------------------------------------------------------------------------
function showGDALdrivers(hObj,event)
	att   = gdalread('','-M');
	long  = {att.Driver.DriverLongName}';
	short = {att.Driver.DriverShortName}';
	list  = cat(2,short,long);
    tableGUI('array',list,'ColWidth',[60 220],'ColNames',{'Short' 'Long Format Name'},...
        'FigName','Potentialy Available GDAL formats','RowNumbers','y','MAX_ROWS',20,'modal','');

% -----------------------------------------------------------------------------
function figure1_KeyPressFcn(hObj, event)
	handles = guidata(hObj);
	if (handles.no_file),	return,		end
	if isequal(get(hObj,'CurrentCharacter'),'+')
		zoom_j(hObj,2,[]);
	elseif isequal(get(hObj,'CurrentCharacter'),'-')
		zoom_j(hObj,0.5,[]);
	end
	hSliders = getappdata(handles.axes1,'SliderAxes');
	if (~isempty(hSliders) && strcmp( get(hSliders(1),'Vis'),'on' ) )	% If (1) is visible so is the other
		CK = get(hObj,'CurrentKey');
		if (strcmp(CK,'rightarrow') || strcmp(CK,'leftarrow'))
			SS = get(hSliders(1),'SliderStep');			val = get(hSliders(1),'Value');
			if (CK(1) == 'r'),		newVal = min(val + SS(1), 1);	% I know that imscroll_j sliders are [0 1]
			else					newVal = max(0, val - SS(1));
			end
			set(hSliders(1),'Value', newVal)
			imscroll_j(handles.axes1,'SetSliderHor')
		elseif (strcmp(CK,'uparrow') || strcmp(CK,'downarrow'))
			SS = get(hSliders(2),'SliderStep');			val = get(hSliders(2),'Value');
			if (CK(1) == 'u'),		newVal = min(val + SS(1), 1);	% I know that imscroll_j sliders are [0 1]
			else					newVal = max(0, val - SS(1));
			end
			set(hSliders(2),'Value', newVal)
			imscroll_j(handles.axes1,'SetSliderVer')
		end
	end

% --------------------------------------------------------------------
function figure1_ResizeFcn(hObj, event)
	handles = guidata(hObj);
	if (isempty(handles)),      return,     end
	screen = get(0,'ScreenSize');	    pos = get(handles.figure1,'Pos');
	if (screen(3) == pos(3))            % Do not allow figure miximizing
		set(handles.figure1,'Pos',handles.oldSize)
	else
		hSliders = getappdata(handles.axes1,'SliderAxes');
		if (~isempty(hSliders))		% Reposition the vertical slider so that its always out of image
			oldUnit = get(handles.axes1, 'Units');	set(handles.axes1, 'Units', 'pixels');
			axPos = get(handles.axes1,'Pos');		sldT = 7;		% from sldT in resizetrue
			set(hSliders(2), 'Pos',[axPos(1)+axPos(3)+1 axPos(2) sldT axPos(4)+1])
			set(handles.axes1, 'Units', oldUnit);
		end
		handles.oldSize = pos;          guidata(handles.figure1,handles)
	end

% -----------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObj, event)
	handles = guidata(hObj);
	try		h = getappdata(handles.figure1,'dependentFigs');
	catch	delete(gcf),	return
	end
	delete(handles.figure1);		delete(h(ishandle(h)))      % Delete also any eventual 'carraças'
	FOpenList = handles.FOpenList;		fname = [handles.path_data 'mirone_pref.mat'];
	if (~handles.version7),    	save(fname,'FOpenList','-append')   % Update the list for "Recent files"
	else    	                save(fname,'FOpenList','-append', '-v6')
	end

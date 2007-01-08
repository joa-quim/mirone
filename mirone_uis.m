function [h1,version7,IamCompiled] = mirone_uis()
% --- Creates and returns a handle to the GUI MIRONE figure. 
%#function pan resetplotview igrf_options rally_plater plate_calculator gmtedit ecran snapshot
%#function about_box parker_stuff plate_calculator euler_stuff grid_calculator tableGUI imageResize
%#function datasets_funs earthquakes manual_pole_adjust compute_euler focal_meca srtm_tool atlas
%#function image_enhance image_adjust datasets_funs InOut2WS write_gmt_script vitrinite telhometro imcapture

h1 = figure('PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn','mirone(''figure1_CloseRequestFcn'',gcbo,[],guidata(gcbo))',...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','Mirone',...
'NumberTitle','off',...
'PaperPositionMode','auto',...
'PaperSize',[20.98404194812 29.67743169791],...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',[520 758 581 21],...
'ResizeFcn','mirone(''figure1_ResizeFcn'',gcbo,[],guidata(gcbo))',...
'Tag','figure1',...
'Visible','off');
setappdata(h1,'IAmAMirone',1);         % Use this appdata to identify Mirone figures

% Detect which matlab version is beeing used. For the moment I'm only interested to know if R13 or >= R14
version7 = version;
if (str2double(version7(1)) > 6),   version7 = 1;
else                                version7 = 0;
end

% The following test will tell us if we are using the compiled or the ML version
try
    %dumb=evalin('base','who');
    dumb=which('mirone');
    IamCompiled = 0;
catch
    IamCompiled = 1;
end

% Import icons
load (['data' filesep 'mirone_icons.mat']);

h_toolbar = uitoolbar('parent',h1, 'BusyAction','queue','HandleVisibility','on',...
   'Interruptible','on','Tag','FigureToolBar','Visible','on');
uipushtool('parent',h_toolbar,'Click','mirone(''FileNewEmpty_CB'',[],[],guidata(gcbo))', ...
   'Tag','NewFigure','cdata',Mfnew_ico,'TooltipString','Open New figure');
uipushtool('parent',h_toolbar,'Click','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''GMT'')', ...
   'Tag','ImportGMTgrid','cdata',Mfopen_ico,'TooltipString','Load GMT grid');
uipushtool('parent',h_toolbar,'Click','mirone(''FileSaveGMTgrid_CB'',[],[],guidata(gcbo))', ...
   'Tag','SaveGMTgrid','cdata',Mfsave_ico,'TooltipString','Save GMT grid');
uipushtool('parent',h_toolbar,'Click','mirone(''FilePreferences_CB'',[],[],guidata(gcbo))', ...
   'Tag','Preferences','cdata',tools_ico,'TooltipString','Preferences');
uipushtool('parent',h_toolbar,'Click','mirone(''FilePrint_CB'',[],[],guidata(gcbo))', ...
   'Tag','Print','cdata',Mprint_ico,'TooltipString','Print image');
uipushtool('parent',h_toolbar,'Click','mirone(''DrawText_CB'',[],[],guidata(gcbo))', ...
   'Tag','DrawText','cdata',text_ico,'TooltipString','Insert Text','Sep','on');
uipushtool('parent',h_toolbar,'Click','mirone(''DrawGeographicalCircle_CB'',[],[],guidata(gcbo))', ...
   'Tag','DrawGeogCirc','cdata',circ_ico,'TooltipString','Draw geographical circle');
uipushtool('parent',h_toolbar,'Click','mirone(''DrawLine_CB'',[],[],guidata(gcbo))', ...
   'Tag','DrawLine','cdata',Mline_ico,'TooltipString','Draw Line');
uipushtool('parent',h_toolbar,'Click','mirone(''DrawClosedPolygon_CB'',[],[],guidata(gcbo),''rectangle'')', ...
   'Tag','DrawRect','cdata',rectang_ico,'TooltipString','Draw Rectangle');
uipushtool('parent',h_toolbar,'Click','mirone(''DrawClosedPolygon_CB'',[],[],guidata(gcbo),[])', ...
   'Tag','DrawPolyg','cdata',polygon_ico,'TooltipString','Draw Closed Polygon');
uipushtool('parent',h_toolbar,'Click','mirone(''DrawVector_CB'',[],[],guidata(gcbo))', ...
   'Tag','DrawArrow','cdata',Marrow_ico,'TooltipString','Draw Arrow');
uitoggletool('parent',h_toolbar,'Click','mirone(''PanZoom_CB'',gcbo,[],guidata(gcbo),''zoom'')', ...
   'Tag','Zoom','cdata',zoom_ico,'TooltipString','Zooming on/off','Sep','on');
uitoggletool('parent',h_toolbar,'Click','mirone(''PanZoom_CB'',gcbo,[],guidata(gcbo),''pan'')', ...
    'cdata',mao,'Tag','Mao','TooltipString','Pan');
uitoggletool('parent',h_toolbar,'Click','draw_funs([],''DeleteObj'')', ...
   'Tag','Tesoura','cdata',cut_ico,'TooltipString','Delete objects');
uipushtool('parent',h_toolbar,'Click','mirone(''ImageColorPalettes_CB'',[],[],guidata(gcbo))', ...
   'Tag','ColorPal','cdata',color_ico,'TooltipString','Color Palettes');
uipushtool('parent',h_toolbar,'Click','mirone(''ImageIlluminationModel_CB'',[],[],guidata(gcbo),''grdgradient_A'')', ...
   'Tag','Shading','cdata',shade2_ico,'TooltipString','Shaded illumination','Sep','on');
uipushtool('parent',h_toolbar,'Click','mirone(''ImageAnaglyph_CB'',[],[],guidata(gcbo))', ...
   'Tag','Anaglyph','cdata',anaglyph_ico,'TooltipString','Anaglyph');
uipushtool('parent',h_toolbar,'Click','mirone(''GridToolsSlope_CB'',[],[],guidata(gcbo),''degrees'')', ...
   'Tag','TerrainMod','cdata',terrain_ico,'TooltipString','Compute Slope');
uipushtool('parent',h_toolbar,'Click','mirone(''ToolsMBplaningStart_CB'',[],[],guidata(gcbo))', ...
   'Tag','MBplaning','cdata',MB_ico,'TooltipString','Multi-beam planing');
uipushtool('parent',h_toolbar,'Click','mirone(''FileSaveFlederSD_CB'',[],[],guidata(gcbo),''runPlanarSD'')', ...
   'Tag','FlederPlanar','cdata',olho_ico,'TooltipString','Run Fleder 3D Viewer');
uipushtool('parent',h_toolbar,'Click',@refresca, 'Tag','Refresh','cdata',refresh_ico,...
   'TooltipString','Refresh','Sep','on');
uipushtool('parent',h_toolbar,'Click','grid_info(guidata(gcbo))','Tag','ImageInfo','cdata',info_ico,'TooltipString','Image info');

axes('Parent',h1,'Units','pixels','Position',[60 0 50 10],'Tag','axes1','Visible','off');

h2 = uimenu('Parent',h1,'Label','File','Tag','File');
uimenu('Parent',h2,'Call','mirone(''FilePreferences_CB'',[],[],guidata(gcbo))','Label','Preferences');
uimenu('Parent',h2,'Callback','mirone(''FileNewEmpty_CB'',[],[],guidata(gcbo))',...
'Label','New empty window','Sep','on');

h5 = uimenu('Parent',h2,'Label','Background window');
uimenu('Parent',h5,'Callback','mirone(''FileNewBgMap_CB'',gcbo,[],guidata(gcbo))','Label','Map');
uimenu('Parent',h5,'Callback','mirone(''FileNewBgFrame_CB'',gcbo,[],guidata(gcbo))','Label','Frame');

if (~IamCompiled)
    uimenu('Parent',h2,'Label','Open a 2D array .mat file','Call','InOut2WS(guidata(gcbo),''loadmat'')','Tag','LoadMAT');
end

h8 = uimenu('Parent',h2,'Label','Open Grid/Image','Sep','on','Tag','OpenGI');
uimenu('Parent',h8,'Call','mirone(''FileOpenNewImage_CB'',gcbo,[],guidata(gcbo));','Label','Images -> Generic Formats');
uimenu('Parent',h8,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''GMT'')','Label','GMT Grid','Sep','on');
uimenu('Parent',h8,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''Surfer'')','Label','Surfer 6/7 grid');
uimenu('Parent',h8,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''ENCOM'')','Label','Encom grid');
uimenu('Parent',h8,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''MANI'')','Label','Mani grid');
uimenu('Parent',h8,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''ArcAscii'')','Label','Arc/Info ASCII Grid','Sep','on');
uimenu('Parent',h8,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''ArcBinary'')','Label','Arc/Info Binary Grid');
uimenu('Parent',h8,'Call','mirone(''FileOpenBSB_CB'',gcbo,[],guidata(gcbo))','Label','BSB Nautical Chart Format');
uimenu('Parent',h8,'Call','mirone(''FileOpen_ENVI_Erdas_CB'',gcbo,[],guidata(gcbo),''ENVI'')','Label','ENVI Raster');
uimenu('Parent',h8,'Call','mirone(''FileOpen_ENVI_Erdas_CB'',gcbo,[],guidata(gcbo),''ERDAS'')','Label','Erdas (.img)');
uimenu('Parent',h8,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''GXF'')','Label','Geosoft GXF');
uimenu('Parent',h8,'Call','mirone(''FileOpenGeoTIFF_CB'',gcbo,[],guidata(gcbo),''geotiff'')','Label','GeoTIFF');
uimenu('Parent',h8,'Call','mirone(''FileOpenGeoTIFF_CB'',gcbo,[],guidata(gcbo),''ecw'')','Label','ECW');
uimenu('Parent',h8,'Call','mirone(''FileOpenGeoTIFF_CB'',gcbo,[],guidata(gcbo),''sid'')','Label','MrSID');
uimenu('Parent',h8,'Call','mirone(''FileOpenGeoTIFF_CB'',gcbo,[],guidata(gcbo),''jp2'')','Label','JPEG2000');
uimenu('Parent',h8,'Call','mirone(''FileOpenGDALmultiBand_CB'',gcbo,[],guidata(gcbo),''ENVISAT'')','Label','ENVISAT','Sep','on');
uimenu('Parent',h8,'Call','mirone(''FileOpenGDALmultiBand_CB'',gcbo,[],guidata(gcbo),''AVHRR'')','Label','AVHRR');
uimenu('Parent',h8,'Call','mirone(''FileOpenGeoTIFF_CB'',gcbo,[],guidata(gcbo),''UNKNOWN'')','Label','Try Luck with GDAL');

h19 = uimenu('Parent',h8,'Label','Digital Elevation','Sep','on');
uimenu('Parent',h19,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''DTED'')','Label','DTED');
uimenu('Parent',h19,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''ESRI_hdr'')','Label','ESRI BIL');
uimenu('Parent',h19,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''GTOPO30'')','Label','GTOPO30');
uimenu('Parent',h19,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''GeoTiff_DEM'')','Label','GeoTIFF DEM');
uimenu('Parent',h19,'Call','mirone(''FileOpenMOLA_CB'',gcbo,[],guidata(gcbo))','Label','MOLA DEM');
uimenu('Parent',h19,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''SRTM1'')','Label','SRTM 1 arcsec');
uimenu('Parent',h19,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''SRTM3'')','Label','SRTM 3 arcsec');
uimenu('Parent',h19,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''SRTM30'')','Label','SRTM30');
uimenu('Parent',h19,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''USGS_DEM'')','Label','USGS DEM');
uimenu('Parent',h19,'Call','mirone(''FileOpenDEM_CB'',[],[],guidata(gcbo),''SDTS'')','Label','USGS SDTS DEM');

uimenu('Parent',h2,'Callback','overview(guidata(gcbo))','Label','Open Overview Tool');
uimenu('Parent',h2,'Call','mirone(''FileOpenSession_CB'',gcbo,[],guidata(gcbo))','Label','Open Session');

%uimenu('Parent',h2,'Call','snapshot(gcf)','Label','Snapchuta');
%uimenu('Parent',h2,'Call','shape_tool(gcf)','Label','Limiares');
% ----------------------- Save Images section
h9 = uimenu('Parent',h2,'Label','Save Image As...','Sep','on');

uimenu('Parent',h9,'Call','snapshot(gcf)', 'Label','Generic Formats');

uimenu('Parent',h9,...
'Call','mirone(''FileSaveImgGrdGdal_CB'',gcbo,[],guidata(gcbo),''GeoTiff'',''img'')', 'Label','GeoTiff','Sep','on');

uimenu('Parent',h9,...
'Call','mirone(''FileSaveImgGrdGdal_CB'',gcbo,[],guidata(gcbo),''Erdas'',''img'')', 'Label','Erdas Imagine');

uimenu('Parent',h9,...
'Call','mirone(''FileSaveImgGrdGdal_CB'',gcbo,[],guidata(gcbo),''Envi'',''img'')', 'Label','Envi .hdr Labelled');

uimenu('Parent',h9,...
'Call','mirone(''FileSaveImgGrdGdal_CB'',gcbo,[],guidata(gcbo),''ESRI'',''img'')', 'Label','ESRI .hdr Labelled');

uimenu('Parent',h9,...
'Call','mirone(''FileSaveImgGrdGdal_CB'',gcbo,[],guidata(gcbo),''JP2K'',''img'')', 'Label','JPEG2000');

uimenu('Parent',h2,'Call','snapshot(gcf,''frame'')', 'Label','Export ...','Tag','noAxes');

if (strncmp(computer,'PC',2))
    h = uimenu('Parent',h2,'Label','Copy to Clipboard','Tag','CopyClip');
    uimenu('Parent',h,'Callback','imcapture(gca,''img'');','Label','Image only')
    uimenu('Parent',h,'Callback','imcapture(gca,''imgAx'');','Label','Image and frame','Tag','noAxes')
end

% ----------------------- Save Grids section
h9 = uimenu('Parent',h2,'Label','Save Grid As...','Sep','on');
uimenu('Parent',h9,'Call','mirone(''FileSaveGMTgrid_CB'',[],[],guidata(gcbo))','Label','GMT grid');
uimenu('Parent',h9,'Call','mirone(''FileSaveGMTgrid_CB'',[],[],guidata(gcbo),''Surfer'')','Label','Surfer 6 grid');
uimenu('Parent',h9,'Call','mirone(''FileSaveENCOMgrid_CB'',gcbo,[],guidata(gcbo))','Label','Encom grid');

uimenu('Parent',h9,...
'Call','mirone(''FileSaveImgGrdGdal_CB'',gcbo,[],guidata(gcbo),''GeoTiff'',''grid'')', 'Label','GeoTiff','Sep','on');

uimenu('Parent',h9,...
'Call','mirone(''FileSaveImgGrdGdal_CB'',gcbo,[],guidata(gcbo),''Erdas'',''grid'')', 'Label','Erdas Imagine');

uimenu('Parent',h9,...
'Call','mirone(''FileSaveImgGrdGdal_CB'',gcbo,[],guidata(gcbo),''Envi'',''grid'')', 'Label','Envi .hdr Labelled');

uimenu('Parent',h9,...
'Call','mirone(''FileSaveImgGrdGdal_CB'',gcbo,[],guidata(gcbo),''ESRI'',''grid'')', 'Label','ESRI .hdr Labelled');

uimenu('Parent',h9,'Call','mirone(''FileSaveImgGrdGdal_CB'',gcbo,[],guidata(gcbo),''JP2K'',''grid'')', 'Label','JPEG2000');

h37 = uimenu('Parent',h2,'Label','Save As 3 GMT grids (R,G,B)');
uimenu('Parent',h37,'Call','mirone(''File_img2GMT_RGBgrids_CB'',gcbo,[],guidata(gcbo))','Label','Image only');
uimenu('Parent',h37,'Call','mirone(''File_img2GMT_RGBgrids_CB'',gcbo,[],guidata(gcbo),''screen'')','Label','Screen capture');

h33 = uimenu('Parent',h2,'Label','Save GMT script','Sep','on');
uimenu('Parent',h33,'Call','write_gmt_script(guidata(gcbo),''csh'')','Label','csh script');
uimenu('Parent',h33,'Call','write_gmt_script(guidata(gcbo),''bat'')','Label','dos batch');

h45 = uimenu('Parent',h2,'Label','Save As Fledermaus Objects');
uimenu('Parent',h45,'Call','mirone(''FileSaveFlederSD_CB'',[],[],guidata(gcbo),''writeSphericalSD'')',...
'Label','Spherical Fledermaus Obj');

uimenu('Parent',h45,'Call','mirone(''FileSaveFlederSD_CB'',[],[],guidata(gcbo),''writePlanarSD'')',...
'Label','Planar Fledermaus Obj');

uimenu('Parent',h2,'Call','mirone(''FileSaveSession_CB'',gcbo,[],guidata(gcbo))','Label','Save Session');

if (~IamCompiled)
    uimenu('Parent',h2,'Label','Grid/Image -> Workspace','Sep','on','Callback','InOut2WS(guidata(gcbo),''direct'')');
    uimenu('Parent',h2,'Label','Workspace -> Grid','Callback','InOut2WS(guidata(gcbo),''GRID_inverse'')');
    uimenu('Parent',h2,'Label','Workspace -> Image','Callback','InOut2WS(guidata(gcbo),''IMG_inverse'')');
    uimenu('Parent',h2,'Label','Clear Workspace','Callback','InOut2WS(guidata(gcbo),''clear'')');
end

uimenu('Parent',h2,'Call','print -dsetup','Label','Print Setup','Sep','on');
uimenu('Parent',h2,'Call','mirone(''FilePrint_CB'',[],[],guidata(gcbo))','Label','Print...');

h54 = uimenu('Parent',h1,'Label','Image','Tag','Image');
uimenu('Parent',h54,'Call','mirone(''ImageColorPalettes_CB'',[],[],guidata(gcbo))','Label','Color Palettes');
uimenu('Parent',h54,'Call','mirone(''ImageShowPalette_CB'',gcbo,[],guidata(gcbo))','Label','Show Palette');
uimenu('Parent',h54,'Call','mirone(''ImageCrop_Callback'',gcbo,[],guidata(gcbo),[])','Label','Crop','Sep','on');

h59 = uimenu('Parent',h54,'Label','Flip');
uimenu('Parent',h59,'Call','mirone(''ImageFlip_CB'',gcbo,[],guidata(gcbo),''UD'')','Label','Flip Up-Down');
uimenu('Parent',h59,'Call','mirone(''ImageFlip_CB'',gcbo,[],guidata(gcbo),''LR'')','Label','Flip Left-Right');

uimenu('Parent',h54,'Call','mirone(''ImageResetOrigImg_CB'',gcbo,[],guidata(gcbo))','Label','Restore Original Image');

uimenu('Parent',h54,'Call','mirone(''ImageHistEqualize_CB'',gcbo,[],guidata(gcbo))',...
'Label','Histogram Equalization (image)','Tag','ImgHist');

uimenu('Parent',h54,'Call','mirone(''ImageHistEqualizeGrid_CB'',gcbo,[],guidata(gcbo))',...
'Label','Histogram Equalization (grid)','Tag','ImgHistGrd');

h64 = uimenu('Parent',h54,'Label','Illuminate','Tag','Illuminate');

uimenu('Parent',h64,'Call','mirone(''ImageIlluminationModel_CB'',[],[],guidata(gcbo),''grdgradient_A'')',...
'Label','GMT grdgradient');

uimenu('Parent',h64,'Call','mirone(''ImageIlluminationModel_CB'',[],[],guidata(gcbo),''falseColor'')',...
'Label','False color');

uimenu('Parent',h54,'Call','mirone(''ImageAnaglyph_CB'',[],[],guidata(gcbo))','Label','Anaglyph');
uimenu('Parent',h54,'Call','mirone(''ImageDrape_CB'',gcbo,[],guidata(gcbo))','Label','Drape','Tag','ImageDrape','Enable','off');
uimenu('Parent',h54,'Call','mirone(''ImageMapLimits_CB'',gcbo,[],guidata(gcbo))','Label','Map Limits');

h9 = uimenu('Parent',h54,'Label','Edge detect (Canny)','Sep','on');
uimenu('Parent',h9,'Call','mirone(''ImageEdgeDetect_CB'',gcbo,[],guidata(gcbo),''Vec'')','Label','Vector');
uimenu('Parent',h9,'Call','mirone(''ImageEdgeDetect_CB'',gcbo,[],guidata(gcbo),''Ras'')','Label','Raster');

h9 = uimenu('Parent',h54,'Label','Edge detect (SUSAN)');
uimenu('Parent',h9,'Call','mirone(''ImageEdgeDetect_CB'',gcbo,[],guidata(gcbo),''SUSvec'')','Label','Vector');
uimenu('Parent',h9,'Call','mirone(''ImageEdgeDetect_CB'',gcbo,[],guidata(gcbo),''SUSras'')','Label','Raster');

h9 = uimenu('Parent',h54,'Label','Features detection');
uimenu('Parent',h9,'Call','mirone(''ImageEdgeDetect_CB'',gcbo,[],guidata(gcbo),''Lines'')','Label','Lines');
uimenu('Parent',h9,'Call','mirone(''ImageEdgeDetect_CB'',gcbo,[],guidata(gcbo),''Circles'')','Label','Circles');
uimenu('Parent',h9,'Call','mirone(''Transfer_CB'',gcbo,[],guidata(gcbo),''Corners'')','Label','Good features to track');

uimenu('Parent',h54,'Label','entry_sh');
%uimenu('Parent',h54,'Call','mirone(''Transfer_CB'',gcbo,[],guidata(gcbo),''Paint'')','Label','PaintBrosh');
uimenu('Parent',h54,'Call','imageResize(gcf)','Label','Image resize','Sep','on');
uimenu('Parent',h54,'Call','mirone(''RotateTool_CB'',gcbo,[],guidata(gcbo),''image'')','Label','Image rotation');

h9 = uimenu('Parent',h54,'Label','Image mode');
uimenu('Parent',h9,'Call','mirone(''Transfer_CB'',gcbo,[],guidata(gcbo),''gray'')','Label','Gray scale');
uimenu('Parent',h9,'Call','mirone(''Transfer_CB'',gcbo,[],guidata(gcbo),''bw'')','Label','Black and White');
uimenu('Parent',h9,'Call','mirone(''ImageResetOrigImg_CB'',gcbo,[],guidata(gcbo))','Label','Original Image');

% % ------------ Image filters _______ TO BE CONTINUED
% h = uimenu('Parent',h54,'Label','Filters','Sep','on');
% uimenu('Parent',h,'Callback','filter_funs(gcf,''smSUSAN'')','Label','Smooth (SUSAN)');

uimenu('Parent',h54,'Call','mirone(''DigitalFilt_CB'',gcbo,[],guidata(gcbo),''image'')',...
'Label','Digital filtering Tool','Sep','on');

uimenu('Parent',h54,'Call','image_enhance(gcf)','Label','Image Enhance (1 - Indexed and RGB)');
uimenu('Parent',h54,'Call','image_adjust(gcf)','Label','Image Enhance (2 - Indexed only)');
uimenu('Parent',h54,'Call','mirone(''ImageGCPtool_CB'',gcbo,[],guidata(gcbo))',...
    'Label','Activate Image-to-Image/Map GCP Tool','Tag','GCPtool','Sep','on');

uimenu('Parent',h54,'Call','mirone(''DrawLine_CB'',[],[],guidata(gcbo),''GCPpline'')','Label','Register Image (Draw GCP points)');
uimenu('Parent',h54,'Callback','bands_list(gcf)','Label','Load Bands','Sep','on');
%uimenu('Parent',h54,'Callback','grid_calculator(gcf)','Label','Bands Arithmetic');

uimenu('Parent',h54,'Call','mirone(''ImageMovieFromList_CB'',gcbo,[],guidata(gcbo))',...
'Label','Make movie from image list','Sep','on');

h76 = uimenu('Parent',h1,'Label','Tools','Tag','Tools');
uimenu('Parent',h76,'Call','mirone(''ExtractProfile_CB'',gcbo,[],guidata(gcbo))','Label','Extract Profile');

h78 = uimenu('Parent',h76,'Label','Measure','Sep','on');
uimenu('Parent',h78,'Call','mirone(''ToolsMeasureDist_CB'',gcbo,[],guidata(gcbo))','Label','Distance','Tag','ToolsMeasureDist');
uimenu('Parent',h78,'Call','mirone(''ToolsMeasureAzimuth_CB'',gcbo,[],guidata(gcbo))','Label','Azimuth');
uimenu('Parent',h78,'Call','mirone(''ToolsMeasureArea_CB'',gcbo,[],guidata(gcbo))','Label','Area');

h82 = uimenu('Parent',h76,'Label','Multi-beam planing','Tag','MBplan','Sep','on');
uimenu('Parent',h82,'Call','mirone(''ToolsMBplaningStart_CB'',[],[],guidata(gcbo))','Label','Start track');
uimenu('Parent',h82,'Call','mirone(''ToolsMBplaningImport_CB'',gcbo,[],guidata(gcbo))','Label','Import track');

uimenu('Parent',h76,'Call','mirone(''DrawClosedPolygon_CB'',[],[],guidata(gcbo),''from_ROI'')',...
'Label','Region-Of-Interest','Sep','on');

uimenu('Parent',h76,'Call','mirone(''Transfer_CB'',gcbo,[],guidata(gcbo),''isGMT'')','Label','Does GMT know this file?','Sep','on');
uimenu('Parent',h76,'Call','mirone(''Transfer_CB'',gcbo,[],guidata(gcbo),''isGDAL'')','Label','Does GDAL know this file?');
uimenu('Parent',h76,'Call','ecran','Label','X,Y grapher','Sep','on');
uimenu('Parent',h76,'Call','gmtedit','Label','gmtedit');
uimenu('Parent',h76,'Call','rally_plater','Label','Rally Plater');
uimenu('Parent',h76,'Label','entry_vtr','Sep','on');

h92 = uimenu('Parent',h1,'Label','Draw','Tag','Draw');
uimenu('Parent',h92,'Call','mirone(''DrawLine_CB'',[],[],guidata(gcbo))','Label','Draw line','Tag','ctrLine');
uimenu('Parent',h92,'Call','mirone(''DrawLine_CB'',[],[],guidata(gcbo),''freehand'')','Label','Freehand draw');
uimenu('Parent',h92,'Call','mirone(''DrawClosedPolygon_CB'',[],[],guidata(gcbo),[])','Label','Draw closed polygon');
uimenu('Parent',h92,'Call','mirone(''DrawClosedPolygon_CB'',[],[],guidata(gcbo),''EulerTrapezium'')','Label','Draw Euler trapezium');
uimenu('Parent',h92,'Call','mirone(''DrawImportLine_CB'',gcbo,[],guidata(gcbo),''AsLine'')','Label','Import line','Sep','on');
uimenu('Parent',h92,'Call','mirone(''DrawImportLine_CB'',gcbo,[],guidata(gcbo),''AsPoint'')','Label','Import points');
uimenu('Parent',h92,'Call','mirone(''DrawImportText_CB'',gcbo,[],guidata(gcbo))','Label','Import text');
uimenu('Parent',h92,'Call','mirone(''DrawImportShape_CB'',gcbo,[],guidata(gcbo))','Label','Import shape file');

h9 = uimenu('Parent',h92,'Label','Draw circle','Sep','on');
uimenu('Parent',h9,'Call','mirone(''DrawGeographicalCircle_CB'',[],[],guidata(gcbo))','Label','Geographical circle');
uimenu('Parent',h9,'Call','mirone(''DrawGeographicalCircle_CB'',[],[],guidata(gcbo),''gcirc'')','Label','Great circle arc');
uimenu('Parent',h9,'Call','mirone(''DrawEulerPoleCircle_CB'',gcbo,[],guidata(gcbo))','Label','Circle about an Euler pole');
uimenu('Parent',h9,'Call','mirone(''DrawGeographicalCircle_CB'',[],[],guidata(gcbo),''cartCirc'')','Label','Cartesian circle');

uimenu('Parent',h92,'Call','mirone(''DrawClosedPolygon_CB'',[],[],guidata(gcbo),''rectangle'')','Label','Draw rectangle');
uimenu('Parent',h92,'Call','mirone(''DrawVector_CB'',[],[],guidata(gcbo))','Label','Draw vector');

h105 = uimenu('Parent',h92,'Label','Draw symbol');
uimenu('Parent',h105,'Call','mirone(''DrawSymbol_CB'',gcbo,[],guidata(gcbo),''circle'')','Label','Circle');
uimenu('Parent',h105,'Call','mirone(''DrawSymbol_CB'',gcbo,[],guidata(gcbo),''square'')','Label','Square');
uimenu('Parent',h105,'Call','mirone(''DrawSymbol_CB'',gcbo,[],guidata(gcbo),''triangle'')','Label','Triangle');
uimenu('Parent',h105,'Call','mirone(''DrawSymbol_CB'',gcbo,[],guidata(gcbo),''star'')','Label','Star');
uimenu('Parent',h105,'Call','mirone(''DrawSymbol_CB'',gcbo,[],guidata(gcbo),''cross'')','Label','Cross');

uimenu('Parent',h92,'Call','mirone(''DrawText_CB'',[],[],guidata(gcbo))','Label','Insert text');

uimenu('Parent',h92,'Callback','mirone(''DrawContours_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Contours (automatic)','Sep','on','Tag','Contours_a');

uimenu('Parent',h92,...
'Callback','mirone(''DrawContours_Callback'',gcbo,[],guidata(gcbo),''gui'')','Label','Contours','Tag','Contours_i');

h113 = uimenu('Parent',h1,'Label','Datasets','Tag','Datasets');
h114 = uimenu('Parent',h113,'Label','Draw coastline','Tag','VoidDatasetsCoastLine');

uimenu('Parent',h114,'Call','mirone(''DatasetsCoastLineNetCDF_CB'',gcbo,[],guidata(gcbo),''c'')',...
'Label','Crude resolution','Tag','DatasetsCoastLineCrude');

uimenu('Parent',h114,'Call','mirone(''DatasetsCoastLineNetCDF_CB'',gcbo,[],guidata(gcbo),''l'')',...
'Label','Low resolution','Tag','DatasetsCoastLineLow');

uimenu('Parent',h114,'Call','mirone(''DatasetsCoastLineNetCDF_CB'',gcbo,[],guidata(gcbo),''i'')',...
'Label','Intermediate resolution','Tag','DatasetsCoastLineInterm');

uimenu('Parent',h114,'Call','mirone(''DatasetsCoastLineNetCDF_CB'',gcbo,[],guidata(gcbo),''h'')',...
'Label','High resolution','Tag','DatasetsCoastLineHigh');

uimenu('Parent',h114,'Call','mirone(''DatasetsCoastLineNetCDF_CB'',gcbo,[],guidata(gcbo),''f'')',...
'Label','Full resolution','Tag','DatasetsCoastLineFull');

h121 = uimenu('Parent',h113,'Label','Draw political boundaries','Tag','VoidDatasetsPB');
h122 = uimenu('Parent',h121,'Label','National boundaries');

uimenu('Parent',h122,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''1'',''c'')',...
'Label','Crude resolution','Tag','DatasetsPBCrude');

uimenu('Parent',h122,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''1'',''l'')',...
'Label','Low resolution','Tag','DatasetsPBLow');

uimenu('Parent',h122,'Callback','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''1'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsPBInterm');

uimenu('Parent',h122,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''1'',''h'')',...
'Label','High resolution','Tag','DatasetsPBHigh');

uimenu('Parent',h122,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''1'',''f'')',...
'Label','Full resolution','Tag','DatasetsPBFull');

h128 = uimenu('Parent',h121,'Label','State boundaries (US)');

uimenu('Parent',h128,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''2'',''c'')',...
'Label','Crude resolution','Tag','DatasetsPBCrude');

uimenu('Parent',h128,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''2'',''l'')',...
'Label','Low resolution','Tag','DatasetsPBLow');

uimenu('Parent',h128,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''2'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsPBInterm');

uimenu('Parent',h128,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''2'',''h'')',...
'Label','High resolution','Tag','DatasetsPBHigh');

uimenu('Parent',h128,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''2'',''f'')',...
'Label','Full resolution','Tag','DatasetsPBFull');

h134 = uimenu('Parent',h121,'Label','Marine boundaries');

uimenu('Parent',h134,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''3'',''c'')',...
'Label','Crude resolution','Tag','DatasetsPBCrude');

uimenu('Parent',h134,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''3'',''l'')',...
'Label','Low resolution','Tag','DatasetsPBLow');

uimenu('Parent',h134,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''3'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsPBInterm');

uimenu('Parent',h134,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''3'',''h'')',...
'Label','High resolution','Tag','DatasetsPBHigh');

uimenu('Parent',h134,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''3'',''f'')',...
'Label','Full resolution','Tag','DatasetsPBFull');

h140 = uimenu('Parent',h121,'Label','All boundaries');

uimenu('Parent',h140,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''a'',''c'')',...
'Label','Crude resolution','Tag','DatasetsPBCrude');

uimenu('Parent',h140,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''a'',''l'')',...
'Label','Low resolution','Tag','DatasetsPBLow');

uimenu('Parent',h140,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''a'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsPBInterm');

uimenu('Parent',h140,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''a'',''h'')',...
'Label','High resolution','Tag','DatasetsPBHigh');

uimenu('Parent',h140,'Call','mirone(''DatasetsPoliticalBoundaries_CB'',gcbo,[],guidata(gcbo),''a'',''f'')',...
'Label','Full resolution','Tag','DatasetsPBFull');

h146 = uimenu('Parent',h113,'Label','Draw rivers','Tag','VoidDatasetsRivers');
h147 = uimenu('Parent',h146,'Label','Permanent major rivers');

uimenu('Parent',h147,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''1'',''c'')',...
'Label','Crude resolution','Tag','DatasetsRiversCrude');

uimenu('Parent',h147,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''1'',''l'')',...
'Label','Low resolution','Tag','DatasetsRiversLow');

uimenu('Parent',h147,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''1'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsRiversInterm');

uimenu('Parent',h147,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''1'',''h'')',...
'Label','High resolution','Tag','DatasetsRiversHigh');

uimenu('Parent',h147,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''1'',''f'')',...
'Label','Full resolution','Tag','DatasetsRiversFull');

h153 = uimenu('Parent',h146,'Label','Additional major rivers');

uimenu('Parent',h153,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''2'',''c'')',...
'Label','Crude resolution','Tag','DatasetsRiversCrude');

uimenu('Parent',h153,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''2'',''l'')',...
'Label','Low resolution','Tag','DatasetsRiversLow');

uimenu('Parent',h153,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''2'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsRiversInterm');

uimenu('Parent',h153,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''2'',''h'')',...
'Label','High resolution','Tag','DatasetsRiversHigh');

uimenu('Parent',h153,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''2'',''f'')',...
'Label','Full resolution','Tag','DatasetsRiversFull');

h159 = uimenu('Parent',h146,'Label','Additional rivers');

uimenu('Parent',h159,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''3'',''c'')',...
'Label','Crude resolution','Tag','DatasetsRiversCrude');

uimenu('Parent',h159,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''3'',''l'')',...
'Label','Low resolution','Tag','DatasetsRiversLow');

uimenu('Parent',h159,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''3'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsRiversInterm');

uimenu('Parent',h159,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''3'',''h'')',...
'Label','High resolution','Tag','DatasetsRiversHigh');

uimenu('Parent',h159,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''3'',''f'')',...
'Label','Full resolution','Tag','DatasetsRiversFull');

h171 = uimenu('Parent',h146,'Label','Intermittent rivers - major');

uimenu('Parent',h171,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''5'',''c'')',...
'Label','Crude resolution','Tag','DatasetsRiversCrude');

uimenu('Parent',h171,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''5'',''l'')',...
'Label','Low resolution','Tag','DatasetsRiversLow');

uimenu('Parent',h171,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''5'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsRiversInterm');

uimenu('Parent',h171,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''5'',''h'')',...
'Label','High resolution','Tag','DatasetsRiversHigh');

uimenu('Parent',h171,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''5'',''f'')',...
'Label','Full resolution','Tag','DatasetsRiversFull');

h177 = uimenu('Parent',h146,'Label','Intermittent rivers - additional');

uimenu('Parent',h177,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''6'',''c'')',...
'Label','Crude resolution','Tag','DatasetsRiversCrude');

uimenu('Parent',h177,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''6'',''l'')',...
'Label','Low resolution','Tag','DatasetsRiversLow');

uimenu('Parent',h177,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''6'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsRiversInterm');

uimenu('Parent',h177,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''6'',''h'')',...
'Label','High resolution','Tag','DatasetsRiversHigh');

uimenu('Parent',h177,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''6'',''f'')',...
'Label','Full resolution','Tag','DatasetsRiversFull');

h183 = uimenu('Parent',h146,'Label','Intermittent rivers - minor');

uimenu('Parent',h183,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''7'',''c'')',...
'Label','Crude resolution','Tag','DatasetsRiversCrude');

uimenu('Parent',h183,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''7'',''l'')',...
'Label','Low resolution','Tag','DatasetsRiversLow');

uimenu('Parent',h183,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''7'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsRiversInterm');

uimenu('Parent',h183,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''7'',''h'')',...
'Label','High resolution','Tag','DatasetsRiversHigh');

uimenu('Parent',h183,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''7'',''f'')',...
'Label','Full resolution','Tag','DatasetsRiversFull');

h207 = uimenu('Parent',h146,'Label','All rivers and canals');

uimenu('Parent',h207,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''a'',''c'')',...
'Label','Crude resolution','Tag','DatasetsRiversCrude');

uimenu('Parent',h207,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''a'',''l'')',...
'Label','Low resolution','Tag','DatasetsRiversLow');

uimenu('Parent',h207,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''a'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsRiversInterm');

uimenu('Parent',h207,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''a'',''h'')',...
'Label','High resolution','Tag','DatasetsRiversHigh');

uimenu('Parent',h207,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''a'',''f'')',...
'Label','Full resolution','Tag','DatasetsRiversFull');

h213 = uimenu('Parent',h146,'Label','All permanent rivers');

uimenu('Parent',h213,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''r'',''c'')',...
'Label','Crude resolution','Tag','DatasetsRiversCrude');

uimenu('Parent',h213,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''r'',''l'')',...
'Label','Low resolution','Tag','DatasetsRiversLow');

uimenu('Parent',h213,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''r'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsRiversInterm');

uimenu('Parent',h213,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''r'',''h'')',...
'Label','High resolution','Tag','DatasetsRiversHigh');

uimenu('Parent',h213,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''r'',''f'')',...
'Label','Full resolution','Tag','DatasetsRiversFull');

h219 = uimenu('Parent',h146,'Label','All intermittent rivers');

uimenu('Parent',h219,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''i'',''c'')',...
'Label','Crude resolution','Tag','DatasetsRiversCrude');

uimenu('Parent',h219,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''i'',''l'')',...
'Label','Low resolution','Tag','DatasetsRiversLow');

uimenu('Parent',h219,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''i'',''i'')',...
'Label','Intermediate resolution','Tag','DatasetsRiversInterm');

uimenu('Parent',h219,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''i'',''h'')',...
'Label','High resolution','Tag','DatasetsRiversHigh');

uimenu('Parent',h219,'Call','mirone(''DatasetsRivers_CB'',gcbo,[],guidata(gcbo),''i'',''f'')',...
'Label','Full resolution','Tag','DatasetsRiversFull');

uimenu('Parent',h113,'Call','earthquakes(gcf);','Label','Global seismicity','Sep','on');
uimenu('Parent',h113,'Call','datasets_funs(''Hotspots'',guidata(gcbo))','Label','Hotspot locations');
uimenu('Parent',h113,'Call','datasets_funs(''Isochrons'',guidata(gcbo))','Label','Magnetic isochrons');
uimenu('Parent',h113,'Call','datasets_funs(''Volcanoes'',guidata(gcbo))','Label','Volcanoes');
uimenu('Parent',h113,'Call','datasets_funs(''Tides'',guidata(gcbo))','Label','Tide Stations');
uimenu('Parent',h113,'Call','datasets_funs(''Plate'',guidata(gcbo))','Label','Plate boundaries');

h237 = uimenu('Parent',h113,'Label','Cities');
uimenu('Parent',h237,'Call','datasets_funs(''Cities'',guidata(gcbo),''major'')','Label','Major cities');
uimenu('Parent',h237,'Call','datasets_funs(''Cities'',guidata(gcbo),''other'')','Label','Other cities');

h240 = uimenu('Parent',h113,'Label','ODP/DSDP sites');
uimenu('Parent',h240,'Call','datasets_funs(''ODP'',guidata(gcbo),''ODP'')','Label','ODP');
uimenu('Parent',h240,'Call','datasets_funs(''ODP'',guidata(gcbo),''DSDP'')','Label','DSDP');
uimenu('Parent',h240,'Call','datasets_funs(''ODP'',guidata(gcbo),''ODP_DSDP'')','Label','ODP and DSDP');

uimenu('Parent',h113,'Call','draw_funs([],''MagneticBarCode'')','Label','Mgnetic Bar Code');
uimenu('Parent',h113,'Call','atlas(guidata(gcbo))','Label','Atlas','Sep','on');
uimenu('Parent',h113,'Call','datasets_funs(''Isochrons'',guidata(gcbo),[])','Label','External db','Sep','on');

h247 = uimenu('Parent',h1,'Label','Geophysics','Tag','Geophysics');
h248 = uimenu('Parent',h247,'Label','Elastic deformation');
uimenu('Parent',h248,'Call','mirone(''DrawLine_CB'',[],[],guidata(gcbo),''FaultTrace'')','Label','Draw Fault');
uimenu('Parent',h248,'Call','mirone(''DrawImportLine_CB'',gcbo,[],guidata(gcbo),''FaultTrace'')',...
'Label','Import Trace Fault');

uimenu('Parent',h248,'Call','mirone(''GeophysicsImportFaultFile_CB'',gcbo,[],guidata(gcbo))',...
'Label','Import Fault File');

h252 = uimenu('Parent',h247,'Label','Tsunami Travel Time','Sep','on');
uimenu('Parent',h252,'Call','mirone(''GeophysicsTTT_CB'',gcbo,[],guidata(gcbo))','Label','Plot source');
uimenu('Parent',h252,'Call','mirone(''GeophysicsTTT_CB'',gcbo,[],guidata(gcbo),''load'')','Label','Import maregrphs and time');
uimenu('Parent',h252,'Call','mirone(''GeophysicsTTT_CB'',gcbo,[],guidata(gcbo),''compute'')','Label','Compute');

h255 = uimenu('Parent',h247,'Label','Swan');
uimenu('Parent',h255,'Call','mirone(''GeophysicsSwanCompute_CB'',gcbo,[],guidata(gcbo))','Label','Compute');
uimenu('Parent',h255,'Call','mirone(''DrawImportLine_CB'',gcbo,[],guidata(gcbo),''AsMaregraph'')',...
'Label','Import Stations','Sep','on');

uimenu('Parent',h255,'Call','mirone(''GeophysicsSwanPlotStations_CB'',gcbo,[],guidata(gcbo))','Label','Plot Stations');

uimenu('Parent',h255,'Call','mirone(''GeophysicsSwanGridBorderStations_CB'',gcbo,[],guidata(gcbo))',...
'Label','Stations on grid borders');

h260 = uimenu('Parent',h247,'Label','Tsun2');
uimenu('Parent',h260,'Call','mirone(''GeophysicsTsun2_CB'',gcbo,[],guidata(gcbo),''compute'')','Label','Compute');
uimenu('Parent',h260,'Call','mirone(''GeophysicsTsun2_CB'',gcbo,[],guidata(gcbo),''write_params'')',...
'Label','Write params file');

uimenu('Parent',h247,'Callback','igrf_options(gcf)','Label','IGRF calculator','Sep','on');
uimenu('Parent',h247,'Callback','parker_stuff(''parker_direct'',gcf)','Label','Parker Direct');
uimenu('Parent',h247,'Callback','parker_stuff(''parker_inverse'',gcf)','Label','Parker Inversion');
uimenu('Parent',h247,'Callback','parker_stuff(''redPole'',gcf)','Label','Reduction to the Pole');
uimenu('Parent',h247,'Callback','plate_calculator','Label','Plate calculator','Sep','on');
uimenu('Parent',h247,'Callback','geog_calculator(gcf)','Label','Geographic calculator');
uimenu('Parent',h247,'Callback','euler_stuff(gcf)','Label','Euler rotations','Sep','on');
uimenu('Parent',h247,'Callback','compute_euler(gcf)','Label','Compute Euler pole');
uimenu('Parent',h247,'Callback','manual_pole_adjust(gcf)','Label','Manual adjust Euler pole');

uimenu('Parent',h247,'Label','entry_tl');

h272 = uimenu('Parent',h247,'Label','Seismicity','Sep','on');
uimenu('Parent',h272,'Callback','earthquakes(gcf,''external'');','Label','Epicenters');
uimenu('Parent',h272,'Callback','focal_meca(gcf)','Label','Focal mechanisms');

h275 = uimenu('Parent',h247,'Label','Import *.gmt files(s)','Sep','on');
uimenu('Parent',h275,'Call','mirone(''GeophysicsImportGmtFile_CB'',gcbo,[],guidata(gcbo))','Label','Single *.gmt file');
uimenu('Parent',h275,'Call','mirone(''GeophysicsImportGmtFileList_CB'',gcbo,[],guidata(gcbo))','Label','List of files');

h278 = uimenu('Parent',h1,'Label','Grid Tools','Tag','GridTools');

if (~IamCompiled)
    uimenu('Parent',h278,'Call','grid_calculator(gcf)','Label','grid calculator');
end

uimenu('Parent',h278,'Call','mirone(''GridToolsGrdfilter_CB'',gcbo,[],guidata(gcbo))','Label','grdfilter','Sep','on');
uimenu('Parent',h278,'Call','mirone(''GridToolsGrdGrad_CB'',gcbo,[],guidata(gcbo))','Label','grdgradient');
uimenu('Parent',h278,'Call','mirone(''GridToolsGrdproject_CB'',gcbo,[],guidata(gcbo))','Label','grdproject');
uimenu('Parent',h278,'Call','mirone(''GridToolsGrdsample_CB'',gcbo,[],guidata(gcbo))','Label','grdsample');
uimenu('Parent',h278,'Call','mirone(''GridToolsGrdtrend_CB'',gcbo,[],guidata(gcbo))','Label','grdtrend');
uimenu('Parent',h278,'Call','mirone(''DigitalFilt_CB'',gcbo,[],guidata(gcbo),''grid'')',...
'Label','Digital filtering Tool','Sep','on');

uimenu('Parent',h278,'Call','mirone(''GridToolsGridClip_CB'',gcbo,[],guidata(gcbo))','Label','Clip Grid');

uimenu('Parent',h278,...
'Callback','mirone(''ImageCrop_Callback'',gcbo,[],guidata(gcbo),[],''CropaGrid'')','Label','Crop Grid');

uimenu('Parent',h278,'Call','mirone(''RotateTool_CB'',gcbo,[],guidata(gcbo),''grid'')','Label','Rotate Grid');
uimenu('Parent',h278,'Call','mirone(''GridToolsHistogram_CB'',gcbo,[],guidata(gcbo))','Label','Histogram');
uimenu('Parent',h278,'Call','mirone(''GridToolsGridMask_CB'',gcbo,[],guidata(gcbo))','Label','Write Mask');
uimenu('Parent',h278,'Call','mirone(''GridToolsSectrum_CB'',gcbo,[],guidata(gcbo), ''Allopts'')','Label','FFT tool','Sep','on');

h9 = uimenu('Parent',h278,'Label','Spectrum');
uimenu('Parent',h9,'Call','mirone(''GridToolsSectrum_CB'',gcbo,[],guidata(gcbo), ''Amplitude'')','Label','Amplitude spectrum');
uimenu('Parent',h9,'Call','mirone(''GridToolsSectrum_CB'',gcbo,[],guidata(gcbo), ''Power'')','Label','Power spectrum');

uimenu('Parent',h278,'Call','mirone(''GridToolsSectrum_CB'',gcbo,[],guidata(gcbo), ''Autocorr'')','Label','Autocorrelation');
uimenu('Parent',h278,'Call','mirone(''GridToolsSmooth_CB'',gcbo,[],guidata(gcbo))','Label','Spline Smooth','Sep','on');

h293 = uimenu('Parent',h278,'Label','SDG');
uimenu('Parent',h293,'Call','mirone(''GridToolsSDG_CB'',gcbo,[],guidata(gcbo),''positive'')','Label','Positive');
uimenu('Parent',h293,'Call','mirone(''GridToolsSDG_CB'',gcbo,[],guidata(gcbo),''negative'')','Label','Negative');
uimenu('Parent',h293,'Call','mirone(''GridToolsSDG_CB'',gcbo,[],guidata(gcbo),[])','Label','Both');

h297 = uimenu('Parent',h278,'Label','Terrain Modeling','Sep','on');
h298 = uimenu('Parent',h297,'Label','Slope');
uimenu('Parent',h298,'Call','mirone(''GridToolsSlope_CB'',[],[],guidata(gcbo),''degrees'')','Label','In degrees');
uimenu('Parent',h298,'Call','mirone(''GridToolsSlope_CB'',[],[],guidata(gcbo),''percent'')','Label','In percentage');
uimenu('Parent',h297,'Call','mirone(''GridToolsSlope_CB'',[],[],guidata(gcbo),''aspect'')','Label','Aspect');

h302 = uimenu('Parent',h297,'Label','Directional derivative');
uimenu('Parent',h302,'Call','mirone(''GridToolsDirDerive_CB'',gcbo,[],guidata(gcbo),''first'')','Label','First derivative');
uimenu('Parent',h302,'Call','mirone(''GridToolsDirDerive_CB'',gcbo,[],guidata(gcbo),''second'')','Label','Second derivative');

h305 = uimenu('Parent',h278,'Label','Interpolate','Sep','on');
uimenu('Parent',h305,'Callback','griding_Mir(gcf,''surface'');', 'Label','Minimum curvature');
uimenu('Parent',h305,'Callback','griding_Mir(gcf,''nearneighbor'');', 'Label','Near neighbor');

h308 = uimenu('Parent',h278,'Label','SRTM tools','Sep','on');

h309 = uimenu('Parent',h308,'Label','SRTM mosaic');
uimenu('Parent',h309,'Callback','srtm_tool','Label','SRTM 3sec');
uimenu('Parent',h309,'Callback','srtm_tool(''srtm1'')','Label','SRTM 1sec');
uimenu('Parent',h309,'Callback','srtm_tool(''srtm30'')','Label','SRTM30');

uimenu('Parent',h308,'Call','mirone(''GridToolsFindHoles_CB'',gcbo,[],guidata(gcbo))','Label','Find holes');
uimenu('Parent',h308,'Call','mirone(''GridToolsSaveAsSRTM_CB'',gcbo,[],guidata(gcbo))','Label','Save as SRTM');

uimenu('Parent',h278,'Call','mirone(''GridToolsPadd2Const_CB'',gcbo,[],guidata(gcbo))',...
'Label','Padd to zero','Sep','on');

uimenu('Parent',h278,'Call','mirone(''ImageEdgeDetect_CB'',gcbo,[],guidata(gcbo),''ppa'')',...
'Label','Extract ridges/valleys','Sep','on');

h9 = uimenu('Parent',h1,'Label','Help','Tag','Help');
uimenu('Parent',h9,'Call','aux_funs(''help'',guidata(gcbo))','Label','Mirone Help (v7)');
uimenu('Parent',h9, 'Call', @showGDALdrivers,'Label','List GDAL formats','Sep','on')
uimenu('Parent',h9,...
'Call','about_box([''Mirone_Last_modified_at_08_Jan_2007''],''Mirone'')','Label','About','Sep','on');

% --------------------------------------------------------------------------------------------------
% We need this function also when the pixval_stsbar got stucked
function refresca(obj,eventdata)
set(get(0,'CurrentFigure'),'Pointer','arrow');     refresh

% --------------------------------------------------------------------------------------------------
function showGDALdrivers(hObj,event)
	att   = gdalread('','-M');
	long  = cellstr(strvcat(att.Driver.DriverLongName));
	short = cellstr(strvcat(att.Driver.DriverShortName));
	list  = cat(2,short,long);
    tableGUI('array',list,'ColWidth',[60 220],'ColNames',{'Short' 'Long Format Name'},...
        'FigName','Potentialy Available GDAL formats','RowNumbers','y','MAX_ROWS',20,'modal','');
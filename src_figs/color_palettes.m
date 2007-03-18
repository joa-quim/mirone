function varargout = color_palettes(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to color_palettes

%	Copyright (c) 2004-2007 by J. Luis
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
color_palettes_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
movegui(hObject,'center')

handles.z_min = [];     handles.z_max = [];     handles.z_min_orig = [];    handles.z_max_orig = [];
if (length(varargin) >= 1 && isstruct(varargin{1}))      % varargin{1} must be the Mirone handles
    handles.hCallingFig = varargin{1}.figure1;
    handles.home_dir = varargin{1}.home_dir;
    handles.work_dir = varargin{1}.work_dir;
    Z = getappdata(handles.hCallingFig,'dem_z');
    if (~isempty(Z))
        handles.have_nans = varargin{1}.have_nans;
        handles.z_min_orig = varargin{1}.head(5);
        handles.z_max_orig = varargin{1}.head(6);
        set(handles.edit_Zmin,'String',handles.z_min)
        set(handles.edit_Zmax,'String',handles.z_max)
    else                        % File was too big to stay on memory
        handles.have_nans = 0;
    end
    % Add this figure handle to the carraças list
    plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
    plugedWin = [plugedWin hObject];
    setappdata(handles.hCallingFig,'dependentFigs',plugedWin);
else
    set(handles.OptionsAutoApply,'checked','off','Enable','off')   % Prevent trying to update an unexisting figure cmap
    set(handles.OptionsApply,'Enable','off')
    set(findobj(hObject,'Tag','FileSavePaletteGrid'),'Visible','off')   % It makes no sense here
    set(handles.edit_Zmax,'Enable','off');    set(handles.edit_Zmin,'Enable','off')
    set(handles.text_MinZ,'Enable','off');    set(handles.text_MaxZ,'Enable','off')
    handles.hCallingFig = [];
    handles.have_nans = 0;
    handles.z_min = [];
    handles.z_max = [];
    handles.home_dir = pwd;     handles.work_dir = handles.home_dir;
end
handles.d_path = [handles.home_dir filesep 'data' filesep];

handles.pal_top = [];       % will contain new colormaps as changed by the the sliders
handles.pal_bot = [];
handles.bg_color = [1 1 1];
handles.z_intervals = [];

%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
% This stupid doesn't allow a frame in the background of an axes. I tried
% everything but the axes is allways behind the frame. So the trick will be to
% change it's size here (it was to small in guide) and let frame3D do the job.
posf_b = get(handles.frame_bot,'Position');
posf_t = get(handles.frame_top,'Position');
set(handles.frame_top,'Position',[posf_t(1) posf_t(2) posf_b(3) posf_t(4)]);
bgcolor = get(0,'DefaultUicontrolBackgroundColor');
framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
h_f = [handles.frame_top handles.frame_bot];
for i=1:numel(h_f)
    frame_size = get(h_f(i),'Position');
    frame3D(hObject,frame_size,framecolor,'',[])
    delete(h_f(i))
end
%------------- END Pro look (3D) -------------------------------------------------------

handles.h_txt_cZ = findobj(hObject,'Style','Text','Tag','text_currentZ');
handles.txt_cZ_pos = get(handles.h_txt_cZ,'Position');

% Show the current colormap in axes
if strcmp(get(handles.OptionsAutoApply,'checked'),'on')
    cmap = get(handles.hCallingFig,'Colormap');
else
    cmap = colormap(jet(256));
end
handles.cmap = cmap;
handles.cmap_original = cmap;
colormap(cmap);      I = 1:length(cmap);
h_img = image(I,'CDataMapping','direct');   set(gca,'YTick',[],'XTick',[]);
set(h_img,'ButtonDownFcn',{@bdn_pal,handles})

set(handles.slider_Top,'max',length(cmap),'min',1,'value',length(cmap))
set(handles.slider_Bottom,'max',length(cmap),'min',1,'value',1)
handles.no_slider = 1;

% Generate lists of available color palettes
palsML = {'ML -- autumn' 'ML -- bone' 'ML -- colorcube' 'ML -- cool' 'ML -- copper' ...
    'ML -- flag' 'ML -- gray' 'ML -- hot' 'ML -- hsv' 'ML -- jet' 'ML -- lines' 'ML -- pink' ...
    'ML -- prism' 'ML -- summer' 'ML -- winter'};
palsGMT = {'GMT -- drywet' 'GMT -- gebco' 'GMT -- globe' 'GMT -- rainbow' ...
    'GMT -- haxby' 'GMT -- no_green' 'GMT -- ocean' 'GMT -- polar' 'GMT -- red2green' ...
    'GMT -- sealand' 'GMT -- seis' 'GMT -- split' 'GMT -- topo' 'GMT -- wysiwyg'};
palsA = {'mag' 'ArcEnCiel' 'circular' 'ChromaDepth' 'Mer' 'MetalChaud' 'Paysage' 'RougeVert' ...
    'Sbm' 'Sismique' 'Terre' 'Terre_Mer' 'Tubulare' 'Tubulare_inv' 'atlas' 'bvr_180' 'bvr_90' ...
    'bvr_clair' 'bvr_sombre' 'pente_90' 'rainbow_hist'};
palsCAR = {'CAR -- Blue' 'CAR -- Carnation' 'CAR -- Cyan' 'CAR -- Desert' 'CAR -- Earth' 'CAR -- Green' ...
    'CAR -- HotMetal' 'CAR -- Jelly' 'CAR -- Magenta' 'CAR -- MorningGlory' 'CAR -- Mustard' ...
    'CAR -- Ocean' 'CAR -- OceanLight' 'CAR -- Olive' 'CAR -- Oysters' 'CAR -- Pumpkin' 'CAR -- Red' ...
    'CAR -- Rose' 'CAR -- Saturn' 'CAR -- Seafloor' 'CAR -- Space' 'CAR -- SuperNova' ...
    'CAR -- Topographic' 'CAR -- TrackLine' 'CAR -- Yellow' 'CAR -- colors10'};
palsGIMP = {'Abstract_1' 'Abstract_2' 'Abstract_3' 'Aneurism' 'Blinds' 'Browns' 'Brushed_Aluminium' ...
    'Burning_Paper' 'Burning_Transparency' 'Caribbean_Blues' 'CD' 'CD_Half' 'Cold_Steel' ...
    'Deep_Sea' 'Flare_Glow_Angular_1' 'Flare_Glow_Radial_2' 'Flare_Radial_102' 'Flare_Radial_103' ...
    'Flare_Rays_Size_1' 'Four_bars' 'Full_saturation_spectrum_CCW' 'Full_saturation_spectrum_CW' ...
    'Golden' 'Greens' 'Horizon_1' 'Incandescent' 'Land_1' 'Land_and_Sea' 'Metallic_Something' ...
    'Nauseating_Headache' 'Pastels' 'Pastel_Rainbow' 'Purples' 'Radial_Eyeball_Blue' ...
    'Radial_Eyeball_Brown' 'Radial_Eyeball_Green' 'Radial_Rainbow_Hoop' 'Rounded_edge' ...
    'Shadows_1' 'Shadows_2' 'Shadows_3' 'Skyline' 'Skyline_polluted' 'Sunrise' ...
    'Tropical_Colors' 'Wood_1' 'Wood_2' 'Yellow_Contrast' 'Yellow_Orange'};

handles.palsML = palsML;        handles.palsGMT = palsGMT;
handles.palsA = palsA;          handles.palsCAR = palsCAR;
handles.palsGIMP = palsGIMP;
%pals = {'Current' palsML{:} palsGMT{:} palsA{:}};
pals = {'Current' palsML{:}};
set(handles.listbox1,'String',pals);

% Choose default command line output for color_palettes_export
handles.output = hObject;
guidata(hObject, handles);
set(hObject,'Visible','on');

% UIWAIT makes color_palettes_export wait for user response (see UIRESUME)
if (nargout)
    set(handles.push_retColorMap,'Visible','on')
    set(handles.edit_Zmax,'Visible','off')      % The other on the row are already hiden by the pushbutton
    uiwait(handles.figure1);
    handles = guidata(handles.figure1);     % Get updated version
    if (handles.killed)          % Don't try to output eventual non-existing variables
        varargout{1} = [];
    else
        varargout{1} = get(handles.figure1,'Colormap');
    end
    delete(handles.figure1);        % The figure can be deleted now
end

% -----------------------------------------------------------------------------------
function listbox1_Callback(hObject, eventdata, handles)
contents = get(hObject,'String');       pal = contents{get(hObject,'Value')};
if (length(pal) > 8 & strcmp(pal(1:8),'Imported'))
    pal = pal(1:8);
end
switch pal
    case 'Current'
        pal = handles.cmap_original;
    case 'Imported'
        pal = handles.imported_cmap;
    case 'ML -- autumn'
        pal = autumn(256);
    case 'ML -- bone'
        pal = bone(256);
    case 'ML -- colorcube'
        pal = colorcube(256);
    case 'ML -- cool'
        pal = cool(256);
    case 'ML -- copper'
        pal = copper(256);
    case 'ML -- flag'
        pal = flag(256);
    case 'ML -- gray'
        pal = gray(256);
    case 'ML -- hot'
        pal = hot(256);
    case 'ML -- hsv'
        pal = hsv(256);
    case 'ML -- jet'
        pal = jet(256);
    case 'ML -- lines'
        pal = lines(256);
    case 'ML -- pink'
        pal = pink(256);
    case 'ML -- prism'
        pal = prism(256);
    case 'ML -- summer'
        pal = summer(256);
    case 'ML -- white'
        pal = white(256);
    case 'ML -- winter'
        pal = winter(256);
    case 'GMT -- drywet'
        load([handles.d_path 'gmt_other_palettes.mat'],'drywet');        pal = drywet;
    case 'GMT -- gebco'
        load([handles.d_path 'gmt_other_palettes.mat'],'gebco');         pal = gebco;
    case 'GMT -- globe'
        load([handles.d_path 'gmt_other_palettes.mat'],'globo');         pal = globo;
    case 'GMT -- gray'
        load([handles.d_path 'gmt_other_palettes.mat'],'gray');          pal = gray;
    case 'GMT -- haxby'
        load([handles.d_path 'gmt_other_palettes.mat'],'haxby');         pal = haxby;
    case 'GMT -- no_green'
        load([handles.d_path 'gmt_other_palettes.mat'],'no_green');      pal = no_green;
    case 'GMT -- ocean'
        load([handles.d_path 'gmt_other_palettes.mat'],'ocean');         pal = ocean;
    case 'GMT -- polar'
        load([handles.d_path 'gmt_other_palettes.mat'],'polar');         pal = polar;
    case 'GMT -- rainbow'
        load([handles.d_path 'gmt_other_palettes.mat'],'rainbow');       pal = rainbow;
    case 'GMT -- red2green'
        load([handles.d_path 'gmt_other_palettes.mat'],'red2green');     pal = red2green;
    case 'GMT -- sealand'
        load([handles.d_path 'gmt_other_palettes.mat'],'sealand');       pal = sealand;
    case 'GMT -- seis'
        load([handles.d_path 'gmt_other_palettes.mat'],'seis');          pal = seis;
    case 'GMT -- split'
        load([handles.d_path 'gmt_other_palettes.mat'],'split');         pal = split;
    case 'GMT -- topo'
        load([handles.d_path 'gmt_other_palettes.mat'],'topo');          pal = topo;
    case 'GMT -- wysiwyg'
        load([handles.d_path 'gmt_other_palettes.mat'],'wysiwyg');       pal = wysiwyg;
    case 'mag'
        load([handles.d_path 'gmt_other_palettes.mat'],'mag');           pal = mag;
    case 'ArcEnCiel'
        load([handles.d_path 'gmt_other_palettes.mat'],'ArcEnCiel');     pal = ArcEnCiel;
    case 'circular'
        load([handles.d_path 'gmt_other_palettes.mat'],'circular');      pal = circular;
    case 'ChromaDepth'
        load([handles.d_path 'gmt_other_palettes.mat'],'ChromaDepth');   pal = ChromaDepth;
    case 'Mer'
        load([handles.d_path 'gmt_other_palettes.mat'],'Mer');           pal = Mer;
    case 'MetalChaud'
        load([handles.d_path 'gmt_other_palettes.mat'],'MetalChaud');    pal = MetalChaud;
    case 'Paysage'
        load([handles.d_path 'gmt_other_palettes.mat'],'Paysage');       pal = Paysage;
    case 'RougeVert'
        load([handles.d_path 'gmt_other_palettes.mat'],'RougeVert');     pal = RougeVert;
    case 'Sbm'
        load([handles.d_path 'gmt_other_palettes.mat'],'Sbm');           pal = Sbm;
    case 'Sismique'
        load([handles.d_path 'gmt_other_palettes.mat'],'Sismique');      pal = Sismique;
    case 'Terre'
        load([handles.d_path 'gmt_other_palettes.mat'],'Terre');         pal = Terre;
    case 'Terre_Mer'
        load([handles.d_path 'gmt_other_palettes.mat'],'Terre_Mer');     pal = Terre_Mer;
    case 'Tubulare'
        load([handles.d_path 'gmt_other_palettes.mat'],'Tubulare');      pal = Tubulare;
    case 'Tubulare_inv'
        load([handles.d_path 'gmt_other_palettes.mat'],'Tubulare_inv');  pal = Tubulare_inv;
    case 'atlas'
        load([handles.d_path 'gmt_other_palettes.mat'],'atlas');         pal = atlas;
    case 'bvr_180'
        load([handles.d_path 'gmt_other_palettes.mat'],'bvr_180');       pal = bvr_180;
    case 'bvr_90'
        load([handles.d_path 'gmt_other_palettes.mat'],'bvr_90');        pal = bvr_90;
    case 'bvr_clair'
        load([handles.d_path 'gmt_other_palettes.mat'],'bvr_clair');     pal = bvr_clair;
    case 'bvr_sombre'
        load([handles.d_path 'gmt_other_palettes.mat'],'bvr_sombre');    pal = bvr_sombre;
    case 'pente_90'
        load([handles.d_path 'gmt_other_palettes.mat'],'pente_90');      pal = pente_90;
    case 'rainbow_hist'
        load([handles.d_path 'gmt_other_palettes.mat'],'rainbow_hist');  pal = rainbow_hist;
    case 'CAR -- Blue'
        load([handles.d_path 'caris256.mat'],'Blue');                    pal = Blue/255;
    case 'CAR -- Carnation'
        load([handles.d_path 'caris256.mat'],'Carnation');               pal = Carnation/255;
    case 'CAR -- Cyan'
        load([handles.d_path 'caris256.mat'],'Cyan');                    pal = Cyan/255;
    case 'CAR -- Desert'
        load([handles.d_path 'caris256.mat'],'Desert');                  pal = Desert/255;
    case 'CAR -- Earth'
        load([handles.d_path 'caris256.mat'],'Earth');                   pal = Earth/255;
    case 'CAR -- Green'
        load([handles.d_path 'caris256.mat'],'Green');                   pal = Green/255;
    case 'CAR -- HotMetal'
        load([handles.d_path 'caris256.mat'],'HotMetal');                pal = HotMetal/255;
    case 'CAR -- Jelly'
        load([handles.d_path 'caris256.mat'],'Jelly');                   pal = Jelly/255;
    case 'CAR -- Magenta'
        load([handles.d_path 'caris256.mat'],'Magenta');                 pal = Magenta/255;
    case 'CAR -- MorningGlory'
        load([handles.d_path 'caris256.mat'],'MorningGlory');            pal = MorningGlory/255;
    case 'CAR -- Mustard'
        load([handles.d_path 'caris256.mat'],'Mustard');                 pal = Mustard/255;
    case 'CAR -- Ocean'
        load([handles.d_path 'caris256.mat'],'Ocean');                   pal = Ocean/255;
    case 'CAR -- OceanLight'
        load([handles.d_path 'caris256.mat'],'OceanLight');              pal = OceanLight/255;
    case 'CAR -- Olive'
        load([handles.d_path 'caris256.mat'],'Olive');                   pal = Olive/255;
    case 'CAR -- Oysters'
        load([handles.d_path 'caris256.mat'],'Oysters');                 pal = Oysters/255;
    case 'CAR -- Pumpkin'
        load([handles.d_path 'caris256.mat'],'Pumpkin');                 pal = Pumpkin/255;
    case 'CAR -- Red'
        load([handles.d_path 'caris256.mat'],'Red');                     pal = Red/255;
    case 'CAR -- Rose'
        load([handles.d_path 'caris256.mat'],'Rose');                    pal = Rose/255;
    case 'CAR -- Saturn'
        load([handles.d_path 'caris256.mat'],'Saturn');                  pal = Saturn/255;
    case 'CAR -- Seafloor'
        load([handles.d_path 'caris256.mat'],'Seafloor');                pal = Seafloor/255;
    case 'CAR -- Space'
        load([handles.d_path 'caris256.mat'],'Space');                   pal = Space/255;
    case 'CAR -- SuperNova'
        load([handles.d_path 'caris256.mat'],'SuperNova');               pal = SuperNova/255;
    case 'CAR -- Topographic'
        load([handles.d_path 'caris256.mat'],'Topographic');             pal = Topographic/255;
    case 'CAR -- TrackLine'
        load([handles.d_path 'caris256.mat'],'TrackLine');               pal = TrackLine/255;
    case 'CAR -- Yellow'
        load([handles.d_path 'caris256.mat'],'Yellow');                  pal = Yellow/255;
    case 'CAR -- colors10'
        load([handles.d_path 'caris256.mat'],'colors10');                pal = colors10/255;
    case 'Abstract_1'
        load([handles.d_path 'gimp256.mat'],'Abstract_1');               pal = Abstract_1/255;
    case 'Abstract_2'
        load([handles.d_path 'gimp256.mat'],'Abstract_2');               pal = Abstract_2/255;
    case 'Abstract_3'
        load([handles.d_path 'gimp256.mat'],'Abstract_3');               pal = Abstract_3/255;
    case 'Aneurism'
        load([handles.d_path 'gimp256.mat'],'Aneurism');                 pal = Aneurism/255;
    case 'Blinds'
        load([handles.d_path 'gimp256.mat'],'Blinds');                   pal = Blinds/255;
    case 'Browns'
        load([handles.d_path 'gimp256.mat'],'Browns');                   pal = Browns/255;
    case 'Brushed_Aluminium'
        load([handles.d_path 'gimp256.mat'],'Brushed_Aluminium');        pal = Brushed_Aluminium/255;
    case 'Burning_Paper'
        load([handles.d_path 'gimp256.mat'],'Burning_Paper');            pal = Burning_Paper/255;
    case 'Burning_Transparency'
        load([handles.d_path 'gimp256.mat'],'Burning_Transparency');     pal = Burning_Transparency/255;
    case 'Caribbean_Blues'
        load([handles.d_path 'gimp256.mat'],'Caribbean_Blues');          pal = Caribbean_Blues/255;
    case 'CD'
        load([handles.d_path 'gimp256.mat'],'CD');                       pal = CD/255;
    case 'CD_Half'
        load([handles.d_path 'gimp256.mat'],'CD_Half');                  pal = CD_Half/255;
    case 'Cold_Steel'
        load([handles.d_path 'gimp256.mat'],'Cold_Steel');               pal = Cold_Steel/255;
    case 'Deep_Sea'
        load([handles.d_path 'gimp256.mat'],'Deep_Sea');                 pal = Deep_Sea/255;
    case 'Flare_Glow_Angular_1'
        load([handles.d_path 'gimp256.mat'],'Flare_Glow_Angular_1');     pal = Flare_Glow_Angular_1/255;
    case 'Flare_Glow_Radial_2'
        load([handles.d_path 'gimp256.mat'],'Flare_Glow_Radial_2');      pal = Flare_Glow_Radial_2/255;
    case 'Flare_Radial_102'
        load([handles.d_path 'gimp256.mat'],'Flare_Radial_102');         pal = Flare_Radial_102/255;
    case 'Flare_Radial_103'
        load([handles.d_path 'gimp256.mat'],'Flare_Radial_103');         pal = Flare_Radial_103/255;
    case 'Flare_Rays_Size_1'
        load([handles.d_path 'gimp256.mat'],'Flare_Rays_Size_1');        pal = Flare_Rays_Size_1/255;
    case 'Four_bars'
        load([handles.d_path 'gimp256.mat'],'Four_bars');                pal = Four_bars/255;
    case 'Full_saturation_spectrum_CCW'
        load([handles.d_path 'gimp256.mat'],'Full_saturation_spectrum_CCW');    pal = Full_saturation_spectrum_CCW/255;
    case 'Full_saturation_spectrum_CW'
        load([handles.d_path 'gimp256.mat'],'Full_saturation_spectrum_CW');     pal = Full_saturation_spectrum_CW/255;
    case 'Golden'
        load([handles.d_path 'gimp256.mat'],'Golden');                   pal = Golden/255;
    case 'Greens'
        load([handles.d_path 'gimp256.mat'],'Greens');                   pal = Greens/255;
    case 'Horizon_1'
        load([handles.d_path 'gimp256.mat'],'Horizon_1');                pal = Horizon_1/255;
    case 'Incandescent'
        load([handles.d_path 'gimp256.mat'],'Incandescent');             pal = Incandescent/255;
    case 'Land_1'
        load([handles.d_path 'gimp256.mat'],'Land_1');                   pal = Land_1/255;
    case 'Land_and_Sea'
        load([handles.d_path 'gimp256.mat'],'Land_and_Sea');             pal = Land_and_Sea/255;
    case 'Metallic_Something'
        load([handles.d_path 'gimp256.mat'],'Metallic_Something');       pal = Metallic_Something/255;
    case 'Nauseating_Headache'
        load([handles.d_path 'gimp256.mat'],'Nauseating_Headache');      pal = Nauseating_Headache/255;
    case 'Pastels'
        load([handles.d_path 'gimp256.mat'],'Pastels');                  pal = Pastels/255;
    case 'Pastel_Rainbow'
        load([handles.d_path 'gimp256.mat'],'Pastel_Rainbow');           pal = Pastel_Rainbow/255;
    case 'Purples'
        load([handles.d_path 'gimp256.mat'],'Purples');                  pal = Purples/255;
    case 'Radial_Eyeball_Blue'
        load([handles.d_path 'gimp256.mat'],'Radial_Eyeball_Blue');      pal = Radial_Eyeball_Blue/255;
    case 'Radial_Eyeball_Brown'
        load([handles.d_path 'gimp256.mat'],'Radial_Eyeball_Brown');     pal = Radial_Eyeball_Brown/255;
    case 'Radial_Eyeball_Green'
        load([handles.d_path 'gimp256.mat'],'Radial_Eyeball_Green');     pal = Radial_Eyeball_Green/255;
    case 'Radial_Rainbow_Hoop'
        load([handles.d_path 'gimp256.mat'],'Radial_Rainbow_Hoop');      pal = Radial_Rainbow_Hoop/255;
    case 'Rounded_edge'
        load([handles.d_path 'gimp256.mat'],'Rounded_edge');             pal = Rounded_edge/255;
    case 'Shadows_1'
        load([handles.d_path 'gimp256.mat'],'Shadows_1');                pal = Shadows_1/255;
    case 'Shadows_2'
        load([handles.d_path 'gimp256.mat'],'Shadows_2');                pal = Shadows_2/255;
    case 'Shadows_3'
        load([handles.d_path 'gimp256.mat'],'Shadows_3');                pal = Shadows_3/255;
    case 'Skyline'
        load([handles.d_path 'gimp256.mat'],'Skyline');                  pal = Skyline/255;
    case 'Skyline_polluted'
        load([handles.d_path 'gimp256.mat'],'Skyline_polluted');         pal = Skyline_polluted/255;
    case 'Sunrise'
        load([handles.d_path 'gimp256.mat'],'Sunrise');                  pal = Sunrise/255;
    case 'Tropical_Colors'
        load([handles.d_path 'gimp256.mat'],'Tropical_Colors');          pal = Tropical_Colors/255;
    case 'Wood_1'
        load([handles.d_path 'gimp256.mat'],'Wood_1');                   pal = Wood_1/255;
    case 'Wood_2'
        load([handles.d_path 'gimp256.mat'],'Wood_2');                   pal = Wood_2/255;
    case 'Yellow_Contrast'
        load([handles.d_path 'gimp256.mat'],'Yellow_Contrast');          pal = Yellow_Contrast/255;
    case 'Yellow_Orange'
        load([handles.d_path 'gimp256.mat'],'Yellow_Orange');            pal = Yellow_Orange/255;
end
% Search eventual color markers and delete them
hp = findobj(handles.axes1,'Tag','Picos');
if (~isempty(hp)),      delete(hp);     end
handles.no_slider = 1;
guidata(hObject, handles);
change_cmap(handles,pal);

% -----------------------------------------------------------------------------------
function slider_Bottom_CreateFcn(hObject, eventdata, handles)
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'BackgroundColor',[.9 .9 .9]);

% -----------------------------------------------------------------------------------
% --- Executes on slider movement.
function slider_Bottom_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = round(get(hObject,'Value'));

if ~isempty(handles.pal_top)        % The other slider has been activated
    val_t = round(get(handles.slider_Top,'Value'));
    cmap = handles.pal_top;         % Make a copy of the full current colormap
else
    cmap = handles.cmap;            % Make a copy of the full current colormap
    val_t = length(handles.cmap);
end

for i = 1:val,       cmap(i,:) = handles.cmap(1,:);     end
%yi = interp1(handles.cmap,linspace(1,length(cmap),length(cmap)-round(val)));
%cmap(round(val)+1:end,:) = yi(:,:);
if (val < val_t)
    yi = interp1(handles.cmap,linspace(1,length(cmap),val_t-val+1));
    cmap(val:val_t,:) = yi(:,:);
elseif (val > val_t)
    yi = interp1(handles.cmap,linspace(1,length(cmap),val-val_t+1));
    yi(1:end,:) = yi(end:-1:1,:);       % Invert the colormap
    cmap(val_t:val,:) = yi(:,:);
    for i = val:length(cmap),   cmap(i,:) = handles.cmap(1,:);   end
    for i = 1:val_t,            cmap(i,:) = handles.cmap(end,:); end
end

handles.pal_bot = cmap;     handles.no_slider = 0;      guidata(hObject, handles);
change_cmap(handles,cmap)

% -----------------------------------------------------------------------------------
function slider_Top_CreateFcn(hObject, eventdata, handles)
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'BackgroundColor',[.9 .9 .9]);

% -----------------------------------------------------------------------------------
% --- Executes on slider movement.
function slider_Top_Callback(hObject, eventdata, handles)
val = round(get(hObject,'Value'));

if ~isempty(handles.pal_bot)        % The other slider has been activated
    val_b = round(get(handles.slider_Bottom,'Value'));
    cmap = handles.pal_bot;         % Make a copy of the full current colormap
else
    cmap = handles.cmap;            % Make a copy of the full current colormap
    val_b = length(handles.cmap);
end

for i = round(val):length(cmap),       cmap(i,:) = handles.cmap(end,:);     end
if (val > val_b)
    yi = interp1(handles.cmap,linspace(1,length(cmap),val-val_b+1));
    cmap(val_b:val,:) = yi(:,:);
elseif (val < val_b)
    yi = interp1(handles.cmap,linspace(1,length(cmap),val_b-val+1));
    yi(1:end,:) = yi(end:-1:1,:);       % Invert the colormap
    cmap(val:val_b,:) = yi(:,:);
    for i = val_b:length(cmap),   cmap(i,:) = handles.cmap(1,:);   end
    for i = 1:val,            cmap(i,:) = handles.cmap(end,:); end
end

handles.pal_top = cmap;     handles.no_slider = 0;      guidata(hObject, handles);
change_cmap(handles,cmap)

% -----------------------------------------------------------------------------------
function change_cmap(handles,pal)
% Change the Image's colormap to 'pal'
if strcmp(get(handles.OptionsAutoApply,'checked'),'on') % otherwise we are just playing with color_palettes alone
    del = 0;
    h = findobj(handles.hCallingFig,'Type','image');
    if (length(size(get(h,'CData'))) > 2)
        h = handles.figure1;
        waitfor(msgbox('True color images do not use color palettes. So you cannot change it.','Warning','warn'))
        del = 1;
    end
    if (del),   delete(h);     return;     end
end

if (~isempty(handles.z_intervals))  % GMT palette with Z levels
    len_Pal = length(pal);    
    z_grd = linspace(handles.z_min_orig,handles.z_max_orig,len_Pal)';
    z_pal = [handles.z_intervals(:,1); handles.z_intervals(end,2)];
    len_zPal = length(z_pal);
    % Interpolate the grid levels into the palette Z levels
    ind_pal = interp1(z_pal,linspace(0,1,len_zPal),z_grd,'linear','extrap');
    ind_pal = round(ind_pal * len_Pal);
    % Map the old pal indices into the new ones
    pal = interp1(linspace(0,len_zPal-1,length(pal)),pal,ind_pal,'linear','extrap');
    handles.z_intervals = [];       % Reset for the case of a new GMT palette import
end

if (handles.z_min ~= handles.z_min_orig | handles.z_max ~= handles.z_max_orig)
% %     %index = fix((C-cmin)/(cmax-cmin)*m)+1
% %     index = fix((double(getappdata(handles.hCallingFig,'dem_z'))- ...
% %         double(handles.z_min))/(double(handles.z_max)-double(handles.z_min))*length(pal))+1;
% %     index(index < 1) = 1;
%     bg_slot = 1;    end_slot = length(pal);
%     if (handles.z_min < handles.z_min_orig)
%         pct = (handles.z_min_orig - handles.z_min) / (handles.z_max_orig - handles.z_min_orig);
%         bg_slot = round(pct * length(pal));
%     elseif (handles.z_min > handles.z_min_orig)
%         pct = (handles.z_min - handles.z_min_orig) / (handles.z_max_orig - handles.z_min_orig);
%         bg_slot = round(pct * length(pal));        
%     end
%     if (handles.z_max > handles.z_max_orig)
%         pct = (handles.z_max - handles.z_max_orig) / (handles.z_max_orig - handles.z_min_orig);
%         end_slot = length(pal) - round(pct * length(pal));
%     elseif (handles.z_max < handles.z_max_orig)
%         pct = (handles.z_max_orig - handles.z_max) / (handles.z_max_orig - handles.z_min_orig);
%         end_slot = length(pal) - round(pct * length(pal));        
%     end
%     pal = pal(bg_slot:end_slot,:);
    len_Pal = length(pal);    
    z_grd = linspace(handles.z_min_orig,handles.z_max_orig,len_Pal)';
    z_pal = linspace(handles.z_min,handles.z_max,len_Pal)';
    len_zPal = length(z_pal);
    % Interpolate the grid levels into the User selected extrema
    ind_pal = interp1(z_pal,linspace(0,1,len_zPal),z_grd,'linear','extrap');
    ind_pal = round(ind_pal * len_Pal);
    % Map the old pal indices into the new ones
    pal = interp1(linspace(0,len_zPal-1,length(pal)),pal,ind_pal,'linear','extrap');
end

pal(pal > 1) = 1; % Sometimes interpolation gives values that are out of [0,1] range...
pal(pal < 0) = 0;
set(handles.figure1,'Colormap',pal);

if strcmp(get(handles.OptionsAutoApply,'checked'),'on')
    if (handles.have_nans)      pal_bg = [handles.bg_color; pal];
    else                        pal_bg = pal;   end
    set(handles.hCallingFig,'Colormap',pal_bg)         % Change the image colormap
end

if (handles.no_slider == 1)         % A new colormap was loaded
    set(handles.slider_Bottom,'max',length(pal),'min',1,'value',1)
    set(handles.slider_Top,'max',length(pal),'min',1,'value',length(pal))
    handles.pal_top = [];    handles.pal_bot = [];
    handles.cmap = pal;
end
guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function FileSavePalette_Callback(hObject, eventdata, handles, opt)
% OPT == [] writes the current cmap as a descrete GMT palette, but with 256 colors (in fact, a continuous cpt)
% OPT == 'master_disc' writes the current cmap as a descrete GMT palette with 16 colors
% OPT == 'master_cont' writes the current cmap as a continuous GMT palette with 16 colors
if (nargin == 3),   opt = [];   end

% Get directory history
if (~isempty(handles.hCallingFig))
    hand_parent = guidata(handles.hCallingFig);
    handles.last_dir = hand_parent.last_dir;
    handles.home_dir = hand_parent.home_dir;
else
    handles.last_dir = pwd;
    handles.home_dir = handles.last_dir;
end

cd(handles.last_dir)
[FileName,PathName] = uiputfile({'*.cpt', 'GMT color palette (*.cpt)'},'Select CPT File name');
if isequal(FileName,0);     return;     end
pause(0.01)

[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if (isempty(EXT) && ~isempty(FileName)),   FileName = [FileName '.cpt'];     end
cd(handles.home_dir);       % allways go home to avoid troubles

pal = get(handles.figure1,'Colormap');
if (handles.have_nans),   pal = pal(2:end,:);   end     % Remove the bg color
if (strcmp(opt,'master_disc') || strcmp(opt,'master_cont'))
    pal_len = 16;
    dz = 1 / pal_len;
    x = linspace(0,1,pal_len)';
    pal = interp1(linspace(0,1,length(pal)),pal,x);
    z_min = 0;
else        % current cmap as a descrete GMT palette, but with 256 colors
    pal_len = size(pal,1);
    dz = (handles.z_max - handles.z_min) / pal_len;
    z_min = handles.z_min;
end

tmp{1} = '# Color palette exported by Mirone';
tmp{2} = '# COLOR_MODEL = RGB';
if (isempty(opt) || strcmp(opt,'master_disc'))   % Current cmap or descrete GMT master palette with 16 colors
	for i=1:pal_len
        cor = round(pal(i,:)*255);
        cor_str = sprintf([num2str(cor(1),'%.1d') '\t' num2str(cor(2),'%.1d') '\t' num2str(cor(3),'%.1d')]);
        z1 = num2str(z_min+dz*(i-1),'%.3f');
        z2 = num2str(z_min+dz*i,'%.3f');
        tmp{i+2} = sprintf([z1 '\t' cor_str '\t' z2 '\t' cor_str]);
	end
else                % Continuous GMT master palette with 16 colors
    for i=1:15
        cor_s = round(pal(i,:)*255);
        cor_e = round(pal(i+1,:)*255);
        cor_s_str = sprintf([num2str(cor_s(1),'%.1d') '\t' num2str(cor_s(2),'%.1d') '\t' num2str(cor_s(3),'%.1d')]);
        cor_e_str = sprintf([num2str(cor_e(1),'%.1d') '\t' num2str(cor_e(2),'%.1d') '\t' num2str(cor_e(3),'%.1d')]);
        z1 = num2str(z_min+dz*(i-1),'%.3f');
        z2 = num2str(z_min+dz*i,'%.3f');
        tmp{i+2} = sprintf([z1 '\t' cor_s_str '\t' z2 '\t' cor_e_str]);
    end
end

tmp{pal_len+3} = sprintf('F\t255\t255\t255');
tmp{pal_len+4} = sprintf('B\t0\t0\t0');

if (handles.have_nans)
    cor = round(pal(1,:)*255);
    cor_str = sprintf(['N\t' num2str(cor(1),'%.1d') '\t' num2str(cor(2),'%.1d') '\t' num2str(cor(3),'%.1d')]);
    tmp{pal_len+5} = sprintf(cor_str);
else
    tmp{pal_len+5} = sprintf('N\t128\t128\t128');
end

fid = fopen([PathName FileName],'wt');
for (i=1:pal_len+5),   fprintf(fid,'%s\n',tmp{i});     end
fclose(fid);

% --------------------------------------------------------------------
function OptionsDiscretizePalette_Callback(hObject, eventdata, handles, opt)
if (nargin == 3),   opt = '16';     end
n_color = str2double(opt);
pal = get(handles.figure1,'Colormap');
if (handles.have_nans),   pal = pal(2:end,:);   end     % Remove the bg color

xi = linspace(0,1,n_color)';
x = linspace(0,1,length(pal));
pal = interp1(x,pal,xi);
% Now we have to replicate the n_colors until we get 256 because uint8 has that many colors
pal256 = [];    n_rep = 256 / n_color;
try         % Wrap it in a try-catch for the case it fails
	for i=1:n_color
        pal256 = [pal256; repmat(pal(i,:),n_rep,1)];
	end
catch
    pal256 = get(handles.figure1,'Colormap');
end
change_cmap(handles,pal256)
set(handles.listbox1,'Value',1)     % Now we have a new "Current cmap"

% --------------------------------------------------------------------
function FileExit_Callback(hObject, eventdata, handles)
    delete(handles.figure1)

% --------------------------------------------------------------------
function OptionsApply_Callback(hObject, eventdata, handles)
if ~isempty(handles.hCallingFig)
    cmap = get(handles.figure1,'Colormap');
    if strcmp(get(handles.OptionsAutoApply,'checked'),'on')
        change_cmap(handles,cmap)
    else
        set(handles.OptionsAutoApply,'checked','on')    % temporarly set it to on
        change_cmap(handles,cmap)                   % in order that change_cmap
        set(handles.OptionsAutoApply,'checked','off')   % may work.
    end
else
    msgbox('Apply what?','Chico Clever')
end

% --------------------------------------------------------------------
function OptionsAutoApply_Callback(hObject, eventdata, handles)
if strcmp(get(handles.OptionsAutoApply,'checked'),'on')
    set(handles.OptionsAutoApply,'checked','off')
else
    if ~isempty(handles.hCallingFig)
        set(handles.OptionsAutoApply,'checked','on')
        change_cmap(handles,get(handles.figure1,'Colormap'))
    else
        set(handles.OptionsAutoApply,'checked','off')
        msgbox('Auto Apply where?','Chico Clever')
    end
end

% --------------------------------------------------------------------
function FileReadPalette_Callback(hObject, eventdata, handles, opt)
if (nargin == 3),   opt = [];   end
str1 = {'*.cpt;*.CPT', 'CPT files (*.cpt,*.CPT)';'*.*', 'All Files (*.*)'};
cd(handles.work_dir);
[FileName,PathName] = uigetfile(str1,'Select CPT file');
cd(handles.home_dir);
if isequal(FileName,0);     return;     end
fname = [PathName FileName];
try
    if (~isempty(opt))  % Use the cpt Z levels as well
        [cmap,handles.z_intervals] = cpt2cmap(['-C' fname]);
    else                % Use only the cpt colors
        cmap = cpt2cmap(['-C' fname]);
    end
catch
    errordlg('There was an error reading the CPT file.','Error')
    return
end
handles.cmap = cmap;        handles.imported_cmap = cmap;
handles.no_slider = 1;
list = get(handles.listbox1,'String');
list{1} = ['Imported (' FileName ')'];
set(handles.listbox1,'String',list,'Value',1);
guidata(gcbo,handles)
change_cmap(handles,handles.cmap);

% --------------------------------------------------------------------
function [z_min,z_max] = min_max_single(Z)
	% Compute the min/max of single precision Z arrays. I need this due to (another) Matlab
	% bug that gives wrong results when the Z (single) array has NaNs. Ouput are doubles.
	z_min = double(min(Z(~isnan(Z(:)))));   z_max = double(max(Z(~isnan(Z(:)))));

% --------------------------------------------------------------------
function push_retColorMap_Callback(hObject, eventdata, handles)
    handles.killed = false;
    guidata(hObject, handles);    uiresume(handles.figure1);

% --------------------------------------------------------------------
function radiobutton_ML_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    des = [handles.radiobutton_GMT handles.radiobutton_MR handles.radiobutton_CAR handles.radiobutton_GIMP];
    set(des,'Value', 0)
    pals = get(handles.listbox1,'String');
    set(handles.listbox1,'String',{pals{1} handles.palsML{:}},'Value',1)
else
    set(hObject,'Value',1)
end

% --------------------------------------------------------------------
function radiobutton_GMT_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    des = [handles.radiobutton_ML handles.radiobutton_MR handles.radiobutton_CAR handles.radiobutton_GIMP];
    set(des,'Value', 0)
    pals = get(handles.listbox1,'String');
    set(handles.listbox1,'String',{pals{1} handles.palsGMT{:}},'Value',1)
else
    set(hObject,'Value',1)
end

% --------------------------------------------------------------------
function radiobutton_MR_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    des = [handles.radiobutton_ML handles.radiobutton_GMT handles.radiobutton_CAR handles.radiobutton_GIMP];
    set(des,'Value', 0)
    pals = get(handles.listbox1,'String');
    set(handles.listbox1,'String',{pals{1} handles.palsA{:}},'Value',1)
else
    set(hObject,'Value',1)
end

% --------------------------------------------------------------------
function radiobutton_CAR_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    des = [handles.radiobutton_ML handles.radiobutton_GMT handles.radiobutton_MR handles.radiobutton_GIMP];
    set(des,'Value', 0)
    pals = get(handles.listbox1,'String');
    set(handles.listbox1,'String',{pals{1} handles.palsCAR{:}},'Value',1)
else
    set(hObject,'Value',1)
end

% --------------------------------------------------------------------
function radiobutton_GIMP_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    des = [ handles.radiobutton_ML handles.radiobutton_GMT handles.radiobutton_MR handles.radiobutton_CAR];
    set(des,'Value', 0)
    pals = get(handles.listbox1,'String');
    set(handles.listbox1,'String',{pals{1} handles.palsGIMP{:}},'Value',1)
else
    set(hObject,'Value',1)
end

% -----------------------------------------------------------------------------------------
function bdn_pal(obj,eventdata,handles)
stype = get(handles.figure1,'selectiontype');
if (~strcmp(stype,'open'))
    return
else
	p = get(handles.axes1,'currentpoint');  px = p(1,1);
    x = get(handles.axes1,'XLim');
	y = get(handles.axes1,'YLim');
	yh = (y(2) - y(1)) / 4;     % Stick height
	yv = yh / 3;                % Stick point height
	xw = 4;                     % Stick width
	pico = [px y(2); px-xw y(2)-yv; px-xw y(2)-yh; px+xw y(2)-yh; px+xw y(2)-yv; px y(2)];
	h = patch(pico(:,1),pico(:,2),'k');
    % Set Z grid value corresponding to the picked palette value
	z_inc = (handles.z_max_orig - handles.z_min_orig) / length(handles.cmap);
    ind_c = round(px - x(1) / (x(2)-x(1)) * length(handles.cmap)) + 1;
	z_cur = handles.z_min_orig + ind_c * z_inc;
	set(handles.h_txt_cZ,'String',['Z = ' num2str(z_cur,'%.3f')],'Visible','on')
	set(h,'Tag','Picos','ButtonDownFcn',{@bdn_pico,handles,h,x,y,yh,yv,xw})
end

% -----------------------------------------------------------------------------------------
function bdn_pico(obj,eventdata,handles,h,x,y,yh,yv,xw)
handles = guidata(handles.figure1);       % We may need to have an updated version of handles
stype = get(handles.figure1,'selectiontype');
if strcmp(stype,'open')
    delete(h)
    set(handles.h_txt_cZ,'String','','Visible','off')   % Hide pico-Z-value
    set(handles.figure1,'WindowButtonMotionFcn','','WindowButtonUpFcn','')
else
	p = get(handles.axes1,'currentpoint');  px = p(1,1);
    ind_c = round(px - x(1) / (x(2)-x(1)) * length(handles.cmap)) + 1;
    set(h,'UserData',ind_c)
    hp = findobj(handles.axes1,'Tag','Picos');
    inds = get(hp,'UserData');
    if (iscell(inds))
        inds = cat(1,inds{:});
    end
    inds = sort([1; inds; length(handles.cmap)]);
    inds = unique(inds);
    ii = find(inds == ind_c);
    ind_b = inds(ii-1);     % Index of previous marker (or 1 if it doesn't exist)
    ind_a = inds(ii+1);     % Index of following marker (or 256 if it doesn't exist)
    set(handles.figure1,'WindowButtonMotionFcn',{@wbm_pico,handles,h,x,y,yh,yv,xw,ind_c,ind_b,ind_a}, ...
        'WindowButtonUpFcn',{@wbu_pico,handles});
end

% -----------------------------------------------------------------------------------------
function wbm_pico(obj,eventdata,handles,h,x,y,yh,yv,xw,ind_old,ind_b,ind_a)
% ind_old is the cmap index of starting marker position
% ind_b is the index of the precedent marker (or 1 if it doesn't exist)
% ind_a is the index of the following marker (or ncolors if it doesn't exist)
p = get(handles.axes1,'currentpoint');  px = p(1,1);
cmap = handles.cmap;
nc = length(cmap);
ind_c = round(px - x(1) / (x(2)-x(1)) * nc) + 1;
if (ind_c < 1 || ind_c > nc),    return;     end     % Don't report outside values
if ((ind_c < ind_b) || (ind_c > ind_a)),   return;     end   % Do not let markers cross each other
%disp([num2str(ind_c) ' ' num2str(ind_b) ' ' num2str(ind_a) ' ' num2str(nc)])

% Uppdate pico-Z-value
z_inc = (handles.z_max_orig - handles.z_min_orig) / (nc - 1);
z_cur = handles.z_min_orig + (ind_c - 1) * z_inc;
set(handles.h_txt_cZ,'String',['Z = ' num2str(z_cur,'%.3f')],'Position',handles.txt_cZ_pos)

nl = ind_old - ind_b + 1;
nu = ind_c - ind_b + 1;

new_cmap_l = interp1(linspace(0,1,nl),cmap(ind_b:ind_b+nl-1,:),linspace(0,1,nu));
new_cmap_u = interp1(linspace(0,1,ind_a-ind_old),cmap(ind_old+1:ind_a,:),linspace(0,1,ind_a-ind_c));

pico = [px y(2); px-xw y(2)-yv; px-xw y(2)-yh; px+xw y(2)-yh; px+xw y(2)-yv; px y(2)];
set(h,'XData',pico(:,1),'YData',pico(:,2),'UserData',ind_c)

cmap = [cmap(1:ind_b-1,:); new_cmap_l; new_cmap_u; cmap(ind_a+1:nc,:)];
change_cmap(handles,cmap)

% -----------------------------------------------------------------------------------------
function wbu_pico(obj,eventdata,handles)
handles = guidata(handles.figure1);       % We may need to have an updated version of handles
handles.cmap = get(handles.figure1,'Colormap');
guidata(handles.figure1,handles)
set(handles.figure1,'WindowButtonMotionFcn','','WindowButtonUpFcn','')

% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
msg = sprintf(['All visible features are so obvious that they don''t need any help.\n' ...
        'There are, however, hiden and very powerful ones:\n\n' ...
        'Double clicking over the colorbar to insert color markers. Drag these\n' ...
        'markers to modify the color map. A second double click deletes the marker.\n\n' ...
        'The Min Z and Max Z fields may be used to impose a fixed colormap between\n' ...
        'those values. Alernatively, import a GMT palette with the "Use Z values."\n' ...
        'In that case, the grid values will be maped into the palette defined Z levels.\n\n' ...
        'The "Finish and return the Color Map" button returns the selected color map\n' ...
        'in output. This button is visible when the form cmap = color_palettes(); is used,\n\n' ...
        'NOTE: If you have previously drawn contours into your image, you can\n' ...
        'easily fit the colors with the contours.']);
msgbox(msg,'Help')

% --------------------------------------------------------------------
function edit_Zmin_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),     set(hObject,'String',num2str(handles.z_min));   return; end
	handles.z_min = xx;
	guidata(hObject,handles)
	change_cmap(handles,get(handles.figure1,'Colormap'))

% --------------------------------------------------------------------
function edit_Zmax_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),     set(hObject,'String',num2str(handles.z_max));   return; end
	handles.z_max = xx;
	guidata(hObject,handles)
	change_cmap(handles,get(handles.figure1,'Colormap'))

% The following is because the stupid dumb compiler refuses to find these functions
% --------------------------------------------------------------------
function c = autumn(m)
    r = (0:m-1)'/max(m-1,1);    c = [ones(m,1) r zeros(m,1)];

% --------------------------------------------------------------------
function b = bone(m)
    b = (7*gray(m) + fliplr(hot(m)))/8;

% --------------------------------------------------------------------
function c = cool(m)
    r = (0:m-1)'/max(m-1,1);    c = [r 1-r ones(m,1)];

% --------------------------------------------------------------------
function c = copper(m)
    c = min(1,gray(m)*diag([1.2500 0.7812 0.4975]));

% --------------------------------------------------------------------
function map = flag(m)
	f = [1 0 0; 1 1 1; 0 0 1; 0 0 0];
	% Generate m/4 vertically stacked copies of f with Kronecker product.
	e = ones(ceil(m/4),1);
	map = kron(e,f);	map = map(1:m,:);

% --------------------------------------------------------------------
function map = lines(n)
    c = get(0,'defaultaxescolororder');
    map = c(rem(0:n-1,size(c,1))+1,:);

% --------------------------------------------------------------------
function p = pink(m)
    p = sqrt((2*gray(m) + hot(m))/3);

% --------------------------------------------------------------------
function map = prism(m)
	% R = [red; orange; yellow; green; blue; violet]
	R = [1 0 0; 1 1/2 0; 1 1 0; 0 1 0; 0 0 1; 2/3 0 1];
	% Generate m/6 vertically stacked copies of r with Kronecker product.
	e = ones(ceil(m/6),1);
	R = kron(e,R);
	R = R(1:m,:);
	map = R;

% --------------------------------------------------------------------
function c = summer(m)
    r = (0:m-1)'/max(m-1,1);    c = [r .5+r/2 .4*ones(m,1)];

% --------------------------------------------------------------------
function c = winter(m)
    r = (0:m-1)'/max(m-1,1);    c = [zeros(m,1) r .5+(1-r)/2];

% --------------------------------------------------------------------
function g = gray(m)
    g = (0:m-1)'/max(m-1,1);    g = [g g g];

% --------------------------------------------------------------------
function map = colorcube(m)
% 1: find biggest cube that can fit in less than the amount of space
%    available (if the perfect cube exactly fits, then drop down the
%    blue resolution by one, because we need extra room for the higher
%    resolution pure color ramps and gray ramp, and the eye is least
%    sensitive to blue).  But don't drop the blue resolution if it
%    is currently two - because we can't go lower than that...

nrgsteps = fix(m^(1/3)+eps);
extra = m-nrgsteps^3;
if (extra == 0) && (nrgsteps > 2)
  nbsteps = nrgsteps - 1;
  extra = m-(nrgsteps^2 * nbsteps);
else
  nbsteps = nrgsteps;
end

% 2: create the colormap consisting of this cube:

rgstep = 1/(nrgsteps-1);
bstep  = 1/(nbsteps-1);
[r,g,b]=meshgrid(0:rgstep:1,0:rgstep:1,0:bstep:1);
map = [r(:) g(:) b(:)];

% 3: remove gray points from white to black (ones where 3 elements are equal values):

diffmap = diff(map')';
summap = sum(abs(diffmap),2);
notgrays = find(summap ~= 0);
map = map(notgrays,:);

% 4: remove pure colors (ones with two elements zero):

summap = [sum(map(:,[1 2]),2) sum(map(:,[2 3]),2) sum(map(:,[1 3]),2)];
map = map(find(min(summap,[],2) ~= 0),:);

% 5: find out how many slots are left (saving one for black)

remlen = m - size(map,1) - 1;

% 6: divide by four, and put in the biggest r, g, b and gray
%    lines that can fit in the remaining slots.  If evenly divisible,
%    each line will have same length.  If not, red/green/blue will
%    have floor(length), and gray will have the extra.

rgbnsteps = floor(remlen / 4);
knsteps   = remlen - 3*rgbnsteps;

rgbstep = 1/(rgbnsteps);
kstep   = 1/(knsteps  );

rgbramp = (rgbstep:rgbstep:1)';
rgbzero = zeros(length(rgbramp), 1);
kramp   = (kstep:kstep:1)';

map = [map                  % cube minus r, g, b and gray ramps
    rgbramp rgbzero rgbzero % red ramp
    rgbzero rgbramp rgbzero % green ramp
    rgbzero rgbzero rgbramp % blue ramp
    0       0       0       % black
    kramp   kramp   kramp]; % gray ramp

% ------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata, handles)
    handles = guidata(handles.figure1);
	if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
        % The GUI is still in UIWAIT, us UIRESUME
        handles.killed = true;      % User gave up, return nothing
        guidata(hObject, handles);    uiresume(handles.figure1);
	else
        % The GUI is no longer waiting, just close it
        delete(handles.figure1)
	end
    
% --------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function color_palettes_LayoutFcn(h1,handles)

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'MenuBar','none',...
'Name','Color Palettes',...
'NumberTitle','off',...
'Position',[520 500 300 388],...
'Resize','off',...
'DoubleBuffer','on',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[5 201 21 128], 'Style','frame', 'Tag','frame_top');
uicontrol('Parent',h1, 'Position',[5 8 285 187], 'Style','frame', 'Tag','frame_bot');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@color_palettes_uicallback,h1,'listbox1_Callback'},...
'HorizontalAlignment','left',...
'Position',[12 14 271 161],...
'String','Color Tables',...
'Style','listbox',...
'Value',1,...
'Tag','listbox1');

uicontrol('Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',{@color_palettes_uicallback,h1,'slider_Bottom_Callback'},...
'CData',[],...
'Position',[12 254 271 16],...
'String',{  '' },...
'Style','slider',...
'Tag','slider_Bottom',...
'UserData',[]);

axes('Parent',h1, 'Units','pixels',...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Position',[12 274 271 50],...
'XColor',get(0,'defaultaxesXColor'),...
'XLim',get(0,'defaultaxesXLim'),...
'XLimMode','manual',...
'XTick',[],...
'XTickMode','manual',...
'YColor',get(0,'defaultaxesYColor'),...
'YTick',[],...
'YTickMode','manual',...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','axes1');

uicontrol('Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',{@color_palettes_uicallback,h1,'slider_Top_Callback'},...
'Position',[12 221 271 16],...
'String',{  '' },...
'Style','slider',...
'Tag','slider_Top');

uicontrol('Parent',h1,...
'Position',[16 238 71 15],...
'String','Stretch Bottom',...
'Style','text',...
'Tag','text2');

uicontrol('Parent',h1,...
'Position',[14 206 60 14],...
'String','Stretch Top',...
'Style','text',...
'Tag','text3');

h14 = uimenu('Parent',h1,'Label','File','Tag','VoidFile');
h15 = uimenu('Parent',h14,'Label','Read GMT palette','Tag','VoidFileReadPalette');

uimenu('Parent',h15,...
'Callback',{@color_palettes_uicallback4,h1,[],'FileReadPalette_Callback'},...
'Label','As master',...
'Tag','FileReadPalette_Master');

uimenu('Parent',h15,...
'Callback',{@color_palettes_uicallback4,h1,'z_levels','FileReadPalette_Callback'},...
'Label','Use Z levels',...
'Tag','FileReadPalette_Zlevels');

h18 = uimenu('Parent',h14,'Label','Save as GMT palette','Tag','VoidFileSavePalette');

uimenu('Parent',h18,...
'Callback',{@color_palettes_uicallback4,h1,[],'FileSavePalette_Callback'},...
'Label','Grid limits with 256 colors',...
'Tag','FileSavePaletteGrid');

uimenu('Parent',h18,...
'Callback',{@color_palettes_uicallback4,h1,'master_disc','FileSavePalette_Callback'},...
'Label','descrete master cpt (16 colors) ',...
'Tag','FileSavePaletteMasterDisc');

uimenu('Parent',h18,...
'Callback',{@color_palettes_uicallback4,h1,'master_cont','FileSavePalette_Callback'},...
'Label','continuous master cpt (16 colors) ',...
'Tag','FileSavePaletteMasterCont');

uimenu('Parent',h14,...
'Callback',{@color_palettes_uicallback,h1,'FileExit_Callback'},...
'Label','Exit','Separator','on','Tag','FileExit');

h23 = uimenu('Parent',h1,'Label','Options','Tag','VoidOptions');

uimenu('Parent',h23,...
'Callback',{@color_palettes_uicallback,h1,'OptionsApply_Callback'},...
'Label','Apply','Tag','OptionsApply');

uimenu('Parent',h23,...
'Callback',{@color_palettes_uicallback,h1,'OptionsAutoApply_Callback'},...
'Checked','on','Label','Auto Apply','Tag','OptionsAutoApply');

h26 = uimenu('Parent',h23,'Label','Discretize Palette','Separator','on','Tag','VoidDiscretizePalette');

uimenu('Parent',h26,...
'Callback',{@color_palettes_uicallback4,h1,'8','OptionsDiscretizePalette_Callback'},...
'Label','8 colors',...
'Tag','OptionsDiscretizePalette8');

uimenu('Parent',h26,...
'Callback',{@color_palettes_uicallback4,h1,'16','OptionsDiscretizePalette_Callback'},...
'Label','16 colors',...
'Tag','OptionsDiscretizePalette16');

uimenu('Parent',h26,...
'Callback',{@color_palettes_uicallback4,h1,'32','OptionsDiscretizePalette_Callback'},...
'Label','32 colors',...
'Tag','OptionsDiscretizePalette32');

uimenu('Parent',h26,...
'Callback',{@color_palettes_uicallback4,h1,'64','OptionsDiscretizePalette_Callback'},...
'Label','64 colors',...
'Tag','OptionsDiscretizePalette64');

uicontrol('Parent',h1,...
'Callback',{@color_palettes_uicallback,h1,'radiobutton_ML_Callback'},...
'Position',[12 176 35 15], 'String','ML', 'Style','radiobutton', 'Tag','radiobutton_ML');

uicontrol('Parent',h1,...
'Callback',{@color_palettes_uicallback,h1,'radiobutton_GMT_Callback'},...
'Position',[60 176 45 15], 'String','GMT', 'Style','radiobutton', 'Tag','radiobutton_GMT');

uicontrol('Parent',h1,...
'Callback',{@color_palettes_uicallback,h1,'radiobutton_MR_Callback'},...
'Position',[113 176 40 15], 'String','MR', 'Style','radiobutton', 'Tag','radiobutton_MR');

uicontrol('Parent',h1,...
'Callback',{@color_palettes_uicallback,h1,'radiobutton_CAR_Callback'},...
'Position',[163 176 40 15],...
'String','CAR',...
'Style','radiobutton',...
'Tag','radiobutton_CAR');

uicontrol('Parent',h1,...
'Callback',{@color_palettes_uicallback,h1,'radiobutton_GIMP_Callback'},...
'Position',[219 176 45 15],...
'String','GIMP',...
'Style','radiobutton',...
'Tag','radiobutton_GIMP');

uimenu('Parent',h1,...
'Callback',{@color_palettes_uicallback,h1,'Help_Callback'},...
'Label','Help', 'Tag','Help');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@color_palettes_uicallback,h1,'edit_Zmin_Callback'},...
'Position',[50 357 80 21],...
'Style','edit',...
'TooltipString','Use a different value to set a fixed color minimum',...
'Tag','edit_Zmin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@color_palettes_uicallback,h1,'edit_Zmax_Callback'},...
'Position',[200 357 80 21],...
'Style','edit',...
'TooltipString','Use a different value to set a fixed color miximum',...
'Tag','edit_Zmax');

uicontrol('Parent',h1,'Position',[160 359 40 15],...
'String','Max Z','Style','text','Tag','text_MaxZ');

uicontrol('Parent',h1,'Position',[10 359 40 15],...
'String','Min Z','Style','text','Tag','text_MinZ');

uicontrol('Parent',h1,...
'FontSize',10, 'HorizontalAlignment','left',...
'Position',[100 333 85 16],...
'Style','text', 'Tag','text_currentZ');

uicontrol('Parent',h1,...
'Callback',{@color_palettes_uicallback,h1,'push_retColorMap_Callback'},...
'FontName','Helvetica','FontSize',9,...
'Position',[10 356 201 23],...
'String','Finish and return the Color Map',...
'Visible','off',...
'Tag','push_retColorMap');

function color_palettes_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

function color_palettes_uicallback4(hObject, eventdata, h1, opt, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1),opt);

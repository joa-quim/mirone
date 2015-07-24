function varargout = color_palettes(varargin)
% Helper window to select/modify a color palette

%	Copyright (c) 2004-2014 by J. Luis
%
% 	This program is part of Mirone and is free software; you can redistribute
% 	it and/or modify it under the terms of the GNU Lesser General Public
% 	License as published by the Free Software Foundation; either
% 	version 2.1 of the License, or any later version.
% 
% 	This program is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% 	Lesser General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id$

	hObject = figure('Vis','off');
	color_palettes_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	handles.z_min = [];     handles.z_max = [];     handles.z_min_orig = [];    handles.z_max_orig = [];	cmap = [];
	handles.have_nans = 0;	handles.hCallingFig = [];	later_ReadPalette = false;
	handles.home_dir = [];
	handles.IAmOctave = (exist('OCTAVE_VERSION','builtin') ~= 0);	% To know if we are running under Octave
	handles.pal_top = [];			% will contain new colormaps as changed by the the sliders
	handles.pal_bot = [];
	handles.bg_color = [1 1 1];
	handles.z_intervals = [];
	handles.thematic = false;		% Some thematic pals will use a pre-set handles.z_intervals
	handles.hinge = false;			% Thematic pals may have a hinge point
	handles.check_custom_pal = true;% Check once the OPTcontrol.txt file for custom CPTs to apear in 'Thematic'
	handles.custom_thematic_name = [];% Will contain eventual CPT names
	handles.txt_cZ_pos = get(handles.h_txt_cZ,'Pos');

	if (nargin == 1 && isstruct(varargin{1}))
		handMir = varargin{1};
		handles.hCallingFig = handMir.figure1;
		if (~handMir.no_file)			% We have something on the Mirone window
			if (ndims(get(handMir.hImg,'CData')) == 3)
				warndlg('True color images do not use color palettes. So you cannot change it.','Warning');
				delete(hObject);		return
			end
			Z = getappdata(handMir.figure1,'dem_z');
			if (~isempty(Z))
				handles.have_nans = handMir.have_nans;
				if (handMir.head(5) == 0 && handMir.head(6) == 0)	% Happens for example with the stacks
					zzz = grdutils(Z,'-L');
					handles.z_min = zzz(1);			handles.z_max = zzz(2);
					handles.z_min_orig = zzz(1);	handles.z_max_orig = zzz(2);
				else
					handles.z_min = handMir.head(5);		handles.z_max = handMir.head(6);
					handles.z_min_orig = handMir.head(5);	handles.z_max_orig = handMir.head(6);
				end
				set(handles.edit_Zmin,'String',handles.z_min_orig)
				set(handles.edit_Zmax,'String',handles.z_max_orig)
			end
		else
			set(handles.OptionsAutoApply,'checked','off','Enable','off')	% Prevent trying to update an unexisting fig cmap
			set(handles.OptionsApply,'Enable','off')
			set([handles.edit_Zmax handles.edit_Zmin handles.text_MinZ handles.text_MaxZ],'Enable','off');
		end
		handles.home_dir = handMir.home_dir;
		handles.work_dir = handMir.work_dir;
		handles.last_dir = handMir.last_dir;
		% Add this figure handle to the carraças list
		plugedWin = getappdata(handMir.figure1,'dependentFigs');
		plugedWin = [plugedWin hObject];
		setappdata(handMir.figure1,'dependentFigs',plugedWin);

		% See if we have a 'thematic_pal' copy from a previous incarnation
		thematic_pal = getappdata(handles.hCallingFig, 'thematic_pal');
		if (~isempty(thematic_pal))
			handles.custom_thematic_name = thematic_pal{1};
			handles.custom_thematic_pal  = thematic_pal{2};
			handles.check_custom_pal = false;
		end
	elseif (nargin == 1 && ischar(varargin{1}))
		later_ReadPalette = true;
		set(handles.OptionsAutoApply,'checked','off','Enable','off')	% Prevent trying to update an unexisting figure cmap
		set(handles.OptionsApply,'Enable','off')
		set([handles.edit_Zmax handles.edit_Zmin handles.text_MinZ handles.text_MaxZ],'Enable','off');
	end

	if (isempty(handles.home_dir))
		mir_dirs = getappdata(0,'MIRONE_DIRS');
		if (~isempty(mir_dirs))
			handles.home_dir = mir_dirs.home_dir;		% Start in values
			handles.work_dir = mir_dirs.work_dir;
			handles.last_dir = mir_dirs.last_dir;
		else
			handles.home_dir = cd;		handles.work_dir = cd;		handles.last_dir = cd;
		end
	end

	handles.d_path = [handles.home_dir filesep 'data' filesep];

	% Generate lists of available color palettes
	palsML = {'ML -- autumn' 'ML -- bone' 'ML -- colorcube' 'ML -- cool' 'ML -- copper' ...
		'ML -- flag' 'ML -- gray' 'ML -- hot' 'ML -- hsv' 'ML -- jet' 'ML -- jet-improved' 'ML -- lines' ...
		'ML -- pink' 'ML -- prism' 'ML -- summer' 'ML -- winter' 'ML -- vivid'};
	handles.palsGMT = {'GMT -- drywet' 'GMT -- gebco' 'GMT -- globe' 'GMT -- rainbow' ...
		'GMT -- haxby' 'GMT -- no_green' 'GMT -- ocean' 'GMT -- polar' 'GMT -- red2green' ...
		'GMT -- sealand' 'GMT -- seis' 'GMT -- split' 'GMT -- topo' 'GMT -- wysiwyg' ...
		'DEM_screen' 'DEM_print' 'DEM_poster'};
	handles.palsA = {'mag' 'ArcEnCiel' 'circular' 'ChromaDepth' 'Mer' 'MetalChaud' 'Paysage' 'RougeVert' ...
		'Sbm' 'Sismique' 'Terre' 'Terre_Mer' 'Tubulare' 'Tubulare_inv' 'atlas' 'bvr_180' 'bvr_90' ...
		'bvr_clair' 'bvr_sombre' 'pente_90' 'rainbow_hist'};
	handles.palsCAR = {'CAR -- Blue' 'CAR -- Carnation' 'CAR -- Cyan' 'CAR -- Desert' 'CAR -- Earth' 'CAR -- Green' ...
		'CAR -- HotMetal' 'CAR -- Jelly' 'CAR -- Magenta' 'CAR -- MorningGlory' 'CAR -- Mustard' ...
		'CAR -- Ocean' 'CAR -- OceanLight' 'CAR -- Olive' 'CAR -- Oysters' 'CAR -- Pumpkin' 'CAR -- Red' ...
		'CAR -- Rose' 'CAR -- Saturn' 'CAR -- Seafloor' 'CAR -- Space' 'CAR -- SuperNova' ...
		'CAR -- Topographic' 'CAR -- TrackLine' 'CAR -- Yellow' 'CAR -- colors10'};
	handles.palsGIMP = {'Abstract_1' 'Abstract_2' 'Abstract_3' 'Aneurism' 'Blinds' 'Browns' 'Brushed_Aluminium' ...
		'Burning_Paper' 'Burning_Transparency' 'Caribbean_Blues' 'CD' 'CD_Half' 'Cold_Steel' ...
		'Deep_Sea' 'Flare_Glow_Angular_1' 'Flare_Glow_Radial_2' 'Flare_Radial_102' 'Flare_Radial_103' ...
		'Flare_Rays_Size_1' 'Four_bars' 'Full_saturation_spectrum_CCW' 'Full_saturation_spectrum_CW' ...
		'Golden' 'Greens' 'Horizon_1' 'Incandescent' 'Land_1' 'Land_and_Sea' 'Metallic_Something' ...
		'Nauseating_Headache' 'Pastels' 'Pastel_Rainbow' 'Purples' 'Radial_Eyeball_Blue' ...
		'Radial_Eyeball_Brown' 'Radial_Eyeball_Green' 'Radial_Rainbow_Hoop' 'Rounded_edge' ...
		'Shadows_1' 'Shadows_2' 'Shadows_3' 'Skyline' 'Skyline_polluted' 'Sunrise' ...
		'Tropical_Colors' 'Wood_1' 'Wood_2' 'Yellow_Contrast' 'Yellow_Orange'};
	handles.palsT = {'Mag - anomaly' 'SeaLand (m)' 'Bathymetry (m)' 'Topography (m)' ...
		'SST (12-26)' 'SST (0-20)' 'SST (0-35)'};

	handles.palsML = palsML;
	pals = [{'Current'} palsML];
	set(handles.listbox1,'String',pals);

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	% This stupid doesn't allow a frame in the background of an axes. I tried
	% everything but the axes is allways behind the frame. So the trick will be to
	% change it's size here (it was to small in guide) and let frame3D do the job.
	posf_b = get(handles.frame_bot,'Pos');
	posf_t = get(handles.frame_top,'Pos');
	set(handles.frame_top,'Pos',[posf_t(1) posf_t(2) posf_b(3) posf_t(4)]);
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) -----------------------------------------------------

	if (later_ReadPalette)		% When pallete filename was transmited in input
		cmap = FileReadPalette_CB([], handles, [], varargin{1});
	end

	% Show the current colormap in axes
	if (isempty(cmap))
		if strcmp(get(handles.OptionsAutoApply,'checked'),'on')
			cmap = get(handles.hCallingFig,'Colormap');
		else
			cmap = colormap(jet(256));
		end
		if (isempty(cmap)),		cmap = colormap(jet(256));	end		% When handles.hCallingFig is empty
	end

	handles.cmap = cmap;
	handles.cmap_original = cmap;
	colormap(cmap);      I = 1:length(cmap);
	h_img = image(I,'CDataMapping','direct');   set(gca,'YTick',[],'XTick',[]);
	set(h_img,'ButtonDownFcn',{@bdn_pal,handles})

	set(handles.slider_Top,'max',length(cmap),'min',1,'value',length(cmap))
	set(handles.slider_Bottom,'max',length(cmap),'min',1,'value',1)
	handles.no_slider = 1;

	% Choose default command line output for color_palettes_export
	handles.output = hObject;
	guidata(hObject, handles);
	set(hObject,'Visible','on');

	% UIWAIT makes color_palettes wait for user response (see UIRESUME)
	if (nargout)
		set(handles.push_retColorMap,'Visible','on')
		set(handles.edit_Zmax,'Visible','off')      % The other on the row are already hiden by the pushbutton
		uiwait(handles.figure1);
		handles = guidata(handles.figure1);     % Get updated version
		if (handles.killed)				% Don't try to output eventual non-existing variables
			varargout{1} = [];
		else
			varargout{1} = get(handles.figure1,'Colormap');
		end
		delete(handles.figure1);		% The figure can be deleted now
	end

% -----------------------------------------------------------------------------------
function listbox1_CB(hObject, handles)
	contents = get(hObject,'String');
	pal = contents{get(hObject,'Value')};
	handles.thematic = false;
	if (numel(pal) > 8 && strcmp(pal(1:8),'Imported'))
		pal = pal(1:8);
	end

	if (get(handles.radio_T, 'Val'))	% Thematic CPTs are treated in a separate function
		thematic_pal(handles, pal)
		return
	end

switch pal
	case 'Current'
		pal = handles.cmap_original;
	case 'Imported'
		pal = handles.imported_cmap;
	case 'ML -- autumn',		pal = autumn(256);
	case 'ML -- bone',			pal = bone(256);
	case 'ML -- colorcube',		pal = colorcube(256);
	case 'ML -- cool',			pal = cool(256);
	case 'ML -- copper',		pal = copper(256);
	case 'ML -- flag',			pal = flag(256);
	case 'ML -- gray',			pal = gray(256);
	case 'ML -- hot',			pal = hot(256);
	case 'ML -- hsv',			pal = hsv(256);
	case 'ML -- jet',			pal = jet(256);
	case 'ML -- jet-improved',	pal = mkpj(256);
	case 'ML -- lines',			pal = lines(256);
	case 'ML -- pink',			pal = pink(256);
	case 'ML -- prism',			pal = prism(256);
	case 'ML -- summer',		pal = summer(256);
	case 'ML -- white',			pal = white(256);
	case 'ML -- winter',		pal = winter(256);
	case 'ML -- vivid',			pal = vivid(256);
	case 'GMT -- drywet'
		load([handles.d_path 'gmt_other_palettes.mat'],'drywet');		pal = drywet;
	case 'GMT -- gebco'
		load([handles.d_path 'gmt_other_palettes.mat'],'gebco');		pal = gebco;
	case 'GMT -- globe'
		load([handles.d_path 'gmt_other_palettes.mat'],'globo');		pal = globo;
	case 'GMT -- gray'
		load([handles.d_path 'gmt_other_palettes.mat'],'gray');			pal = gray;
	case 'GMT -- haxby'
		load([handles.d_path 'gmt_other_palettes.mat'],'haxby');		pal = haxby;
	case 'GMT -- no_green'
		load([handles.d_path 'gmt_other_palettes.mat'],'no_green');		pal = no_green;
	case 'GMT -- ocean'
		load([handles.d_path 'gmt_other_palettes.mat'],'ocean');		pal = ocean;
	case 'GMT -- polar'
		load([handles.d_path 'gmt_other_palettes.mat'],'polar');		pal = polar;
	case 'GMT -- rainbow'
		load([handles.d_path 'gmt_other_palettes.mat'],'rainbow');		pal = rainbow;
	case 'GMT -- red2green'
		load([handles.d_path 'gmt_other_palettes.mat'],'red2green');	pal = red2green;
	case 'GMT -- sealand'
		load([handles.d_path 'gmt_other_palettes.mat'],'sealand');		pal = sealand;
	case 'GMT -- seis'
		load([handles.d_path 'gmt_other_palettes.mat'],'seis');			pal = seis;
	case 'GMT -- split'
		load([handles.d_path 'gmt_other_palettes.mat'],'split');		pal = split;
	case 'GMT -- topo'
		load([handles.d_path 'gmt_other_palettes.mat'],'topo');			pal = topo;
	case 'GMT -- wysiwyg'
        load([handles.d_path 'gmt_other_palettes.mat'],'wysiwyg');		pal = wysiwyg;
	case 'DEM_screen'
        load([handles.d_path 'gmt_other_palettes.mat'],'DEM_screen');	pal = DEM_screen;
	case 'DEM_print'
        load([handles.d_path 'gmt_other_palettes.mat'],'DEM_print');	pal = DEM_print;
	case 'DEM_poster'
        load([handles.d_path 'gmt_other_palettes.mat'],'DEM_poster');	pal = DEM_poster;
	case 'mag'
        load([handles.d_path 'gmt_other_palettes.mat'],'mag');			pal = mag;
	case 'ArcEnCiel'
		load([handles.d_path 'gmt_other_palettes.mat'],'ArcEnCiel');	pal = ArcEnCiel;
	case 'circular'
		load([handles.d_path 'gmt_other_palettes.mat'],'circular');		pal = circular;
	case 'ChromaDepth'
		load([handles.d_path 'gmt_other_palettes.mat'],'ChromaDepth');	pal = ChromaDepth;
	case 'Mer'
        load([handles.d_path 'gmt_other_palettes.mat'],'Mer');			pal = Mer;
	case 'MetalChaud'
        load([handles.d_path 'gmt_other_palettes.mat'],'MetalChaud');	pal = MetalChaud;
	case 'Paysage'
        load([handles.d_path 'gmt_other_palettes.mat'],'Paysage');		pal = Paysage;
	case 'RougeVert'
        load([handles.d_path 'gmt_other_palettes.mat'],'RougeVert');	pal = RougeVert;
	case 'Sbm'
        load([handles.d_path 'gmt_other_palettes.mat'],'Sbm');			pal = Sbm;
	case 'Sismique'
        load([handles.d_path 'gmt_other_palettes.mat'],'Sismique');		pal = Sismique;
	case 'Terre'
        load([handles.d_path 'gmt_other_palettes.mat'],'Terre');		pal = Terre;
	case 'Terre_Mer'
        load([handles.d_path 'gmt_other_palettes.mat'],'Terre_Mer');	pal = Terre_Mer;
	case 'Tubulare'
        load([handles.d_path 'gmt_other_palettes.mat'],'Tubulare');		pal = Tubulare;
	case 'Tubulare_inv'
        load([handles.d_path 'gmt_other_palettes.mat'],'Tubulare_inv');	pal = Tubulare_inv;
	case 'atlas'
        load([handles.d_path 'gmt_other_palettes.mat'],'atlas');		pal = atlas;
	case 'bvr_180'
        load([handles.d_path 'gmt_other_palettes.mat'],'bvr_180');		pal = bvr_180;
	case 'bvr_90'
        load([handles.d_path 'gmt_other_palettes.mat'],'bvr_90');		pal = bvr_90;
	case 'bvr_clair'
        load([handles.d_path 'gmt_other_palettes.mat'],'bvr_clair');	pal = bvr_clair;
	case 'bvr_sombre'
        load([handles.d_path 'gmt_other_palettes.mat'],'bvr_sombre');	pal = bvr_sombre;
	case 'pente_90'
        load([handles.d_path 'gmt_other_palettes.mat'],'pente_90');		pal = pente_90;
	case 'rainbow_hist'
        load([handles.d_path 'gmt_other_palettes.mat'],'rainbow_hist');	pal = rainbow_hist;
	case 'CAR -- Blue'
        load([handles.d_path 'caris256.mat'],'Blue');					pal = Blue/255;
	case 'CAR -- Carnation'
        load([handles.d_path 'caris256.mat'],'Carnation');				pal = Carnation/255;
	case 'CAR -- Cyan'
        load([handles.d_path 'caris256.mat'],'Cyan');					pal = Cyan/255;
	case 'CAR -- Desert'
        load([handles.d_path 'caris256.mat'],'Desert');					pal = Desert/255;
	case 'CAR -- Earth'
		load([handles.d_path 'caris256.mat'],'Earth');					pal = Earth/255;
	case 'CAR -- Green'
        load([handles.d_path 'caris256.mat'],'Green');					pal = Green/255;
	case 'CAR -- HotMetal'
        load([handles.d_path 'caris256.mat'],'HotMetal');				pal = HotMetal/255;
	case 'CAR -- Jelly'
        load([handles.d_path 'caris256.mat'],'Jelly');					pal = Jelly/255;
	case 'CAR -- Magenta'
        load([handles.d_path 'caris256.mat'],'Magenta');				pal = Magenta/255;
	case 'CAR -- MorningGlory'
        load([handles.d_path 'caris256.mat'],'MorningGlory');			pal = MorningGlory/255;
	case 'CAR -- Mustard'
        load([handles.d_path 'caris256.mat'],'Mustard');				pal = Mustard/255;
	case 'CAR -- Ocean'
        load([handles.d_path 'caris256.mat'],'Ocean');					pal = Ocean/255;
	case 'CAR -- OceanLight'
        load([handles.d_path 'caris256.mat'],'OceanLight');				pal = OceanLight/255;
	case 'CAR -- Olive'
        load([handles.d_path 'caris256.mat'],'Olive');					pal = Olive/255;
	case 'CAR -- Oysters'
        load([handles.d_path 'caris256.mat'],'Oysters');				pal = Oysters/255;
	case 'CAR -- Pumpkin'
        load([handles.d_path 'caris256.mat'],'Pumpkin');				pal = Pumpkin/255;
	case 'CAR -- Red'
        load([handles.d_path 'caris256.mat'],'Red');					pal = Red/255;
	case 'CAR -- Rose'
        load([handles.d_path 'caris256.mat'],'Rose');					pal = Rose/255;
	case 'CAR -- Saturn'
        load([handles.d_path 'caris256.mat'],'Saturn');					pal = Saturn/255;
	case 'CAR -- Seafloor'
        load([handles.d_path 'caris256.mat'],'Seafloor');				pal = Seafloor/255;
	case 'CAR -- Space'
        load([handles.d_path 'caris256.mat'],'Space');					pal = Space/255;
	case 'CAR -- SuperNova'
        load([handles.d_path 'caris256.mat'],'SuperNova');				pal = SuperNova/255;
	case 'CAR -- Topographic'
        load([handles.d_path 'caris256.mat'],'Topographic');			pal = Topographic/255;
	case 'CAR -- TrackLine'
        load([handles.d_path 'caris256.mat'],'TrackLine');				pal = TrackLine/255;
	case 'CAR -- Yellow'
        load([handles.d_path 'caris256.mat'],'Yellow');					pal = Yellow/255;
	case 'CAR -- colors10'
        load([handles.d_path 'caris256.mat'],'colors10');				pal = colors10/255;
	case 'Abstract_1'
        load([handles.d_path 'gimp256.mat'],'Abstract_1');				pal = Abstract_1/255;
	case 'Abstract_2'
        load([handles.d_path 'gimp256.mat'],'Abstract_2');				pal = Abstract_2/255;
	case 'Abstract_3'
        load([handles.d_path 'gimp256.mat'],'Abstract_3');				pal = Abstract_3/255;
	case 'Aneurism'
        load([handles.d_path 'gimp256.mat'],'Aneurism');				pal = Aneurism/255;
	case 'Blinds'
        load([handles.d_path 'gimp256.mat'],'Blinds');					pal = Blinds/255;
	case 'Browns'
        load([handles.d_path 'gimp256.mat'],'Browns');					pal = Browns/255;
	case 'Brushed_Aluminium'
        load([handles.d_path 'gimp256.mat'],'Brushed_Aluminium');		pal = Brushed_Aluminium/255;
	case 'Burning_Paper'
        load([handles.d_path 'gimp256.mat'],'Burning_Paper');			pal = Burning_Paper/255;
	case 'Burning_Transparency'
        load([handles.d_path 'gimp256.mat'],'Burning_Transparency');	pal = Burning_Transparency/255;
	case 'Caribbean_Blues'
        load([handles.d_path 'gimp256.mat'],'Caribbean_Blues');			pal = Caribbean_Blues/255;
	case 'CD'
        load([handles.d_path 'gimp256.mat'],'CD');						pal = CD/255;
	case 'CD_Half'
        load([handles.d_path 'gimp256.mat'],'CD_Half');					pal = CD_Half/255;
	case 'Cold_Steel'
        load([handles.d_path 'gimp256.mat'],'Cold_Steel');				pal = Cold_Steel/255;
	case 'Deep_Sea'
        load([handles.d_path 'gimp256.mat'],'Deep_Sea');				pal = Deep_Sea/255;
	case 'Flare_Glow_Angular_1'
        load([handles.d_path 'gimp256.mat'],'Flare_Glow_Angular_1');	pal = Flare_Glow_Angular_1/255;
	case 'Flare_Glow_Radial_2'
        load([handles.d_path 'gimp256.mat'],'Flare_Glow_Radial_2');		pal = Flare_Glow_Radial_2/255;
	case 'Flare_Radial_102'
        load([handles.d_path 'gimp256.mat'],'Flare_Radial_102');		pal = Flare_Radial_102/255;
	case 'Flare_Radial_103'
        load([handles.d_path 'gimp256.mat'],'Flare_Radial_103');		pal = Flare_Radial_103/255;
	case 'Flare_Rays_Size_1'
        load([handles.d_path 'gimp256.mat'],'Flare_Rays_Size_1');		pal = Flare_Rays_Size_1/255;
	case 'Four_bars'
        load([handles.d_path 'gimp256.mat'],'Four_bars');				pal = Four_bars/255;
	case 'Full_saturation_spectrum_CCW'
        load([handles.d_path 'gimp256.mat'],'Full_saturation_spectrum_CCW');	pal = Full_saturation_spectrum_CCW/255;
	case 'Full_saturation_spectrum_CW'
        load([handles.d_path 'gimp256.mat'],'Full_saturation_spectrum_CW');		pal = Full_saturation_spectrum_CW/255;
	case 'Golden'
        load([handles.d_path 'gimp256.mat'],'Golden');					pal = Golden/255;
	case 'Greens'
        load([handles.d_path 'gimp256.mat'],'Greens');					pal = Greens/255;
	case 'Horizon_1'
        load([handles.d_path 'gimp256.mat'],'Horizon_1');				pal = Horizon_1/255;
	case 'Incandescent'
        load([handles.d_path 'gimp256.mat'],'Incandescent');			pal = Incandescent/255;
	case 'Land_1'
        load([handles.d_path 'gimp256.mat'],'Land_1');					pal = Land_1/255;
	case 'Land_and_Sea'
        load([handles.d_path 'gimp256.mat'],'Land_and_Sea');			pal = Land_and_Sea/255;
	case 'Metallic_Something'
        load([handles.d_path 'gimp256.mat'],'Metallic_Something');		pal = Metallic_Something/255;
	case 'Nauseating_Headache'
        load([handles.d_path 'gimp256.mat'],'Nauseating_Headache');		pal = Nauseating_Headache/255;
	case 'Pastels'
        load([handles.d_path 'gimp256.mat'],'Pastels');					pal = Pastels/255;
	case 'Pastel_Rainbow'
        load([handles.d_path 'gimp256.mat'],'Pastel_Rainbow');			pal = Pastel_Rainbow/255;
	case 'Purples'
        load([handles.d_path 'gimp256.mat'],'Purples');					pal = Purples/255;
	case 'Radial_Eyeball_Blue'
        load([handles.d_path 'gimp256.mat'],'Radial_Eyeball_Blue');		pal = Radial_Eyeball_Blue/255;
	case 'Radial_Eyeball_Brown'
        load([handles.d_path 'gimp256.mat'],'Radial_Eyeball_Brown');	pal = Radial_Eyeball_Brown/255;
	case 'Radial_Eyeball_Green'
        load([handles.d_path 'gimp256.mat'],'Radial_Eyeball_Green');	pal = Radial_Eyeball_Green/255;
	case 'Radial_Rainbow_Hoop'
        load([handles.d_path 'gimp256.mat'],'Radial_Rainbow_Hoop');		pal = Radial_Rainbow_Hoop/255;
	case 'Rounded_edge'
        load([handles.d_path 'gimp256.mat'],'Rounded_edge');			pal = Rounded_edge/255;
	case 'Shadows_1'
        load([handles.d_path 'gimp256.mat'],'Shadows_1');				pal = Shadows_1/255;
	case 'Shadows_2'
		load([handles.d_path 'gimp256.mat'],'Shadows_2');				pal = Shadows_2/255;
	case 'Shadows_3'
		load([handles.d_path 'gimp256.mat'],'Shadows_3');				pal = Shadows_3/255;
	case 'Skyline'
		load([handles.d_path 'gimp256.mat'],'Skyline');					pal = Skyline/255;
	case 'Skyline_polluted'
		load([handles.d_path 'gimp256.mat'],'Skyline_polluted');		pal = Skyline_polluted/255;
	case 'Sunrise'
		load([handles.d_path 'gimp256.mat'],'Sunrise');					pal = Sunrise/255;
	case 'Tropical_Colors'
		load([handles.d_path 'gimp256.mat'],'Tropical_Colors');			pal = Tropical_Colors/255;
	case 'Wood_1'
		load([handles.d_path 'gimp256.mat'],'Wood_1');					pal = Wood_1/255;
	case 'Wood_2'
		load([handles.d_path 'gimp256.mat'],'Wood_2');					pal = Wood_2/255;
	case 'Yellow_Contrast'
		load([handles.d_path 'gimp256.mat'],'Yellow_Contrast');			pal = Yellow_Contrast/255;
	case 'Yellow_Orange'
		load([handles.d_path 'gimp256.mat'],'Yellow_Orange');			pal = Yellow_Orange/255;
end

	% Search eventual color markers and delete them
	hp = findobj(handles.axes1,'Tag','Picos');
	if (~isempty(hp)),		delete(hp);		end
	handles.no_slider = 1;
	guidata(handles.figure1, handles);
	change_cmap(handles,pal);

% -----------------------------------------------------------------------------------
function thematic_pal(handles, pal)
% Deal with the 'Thematic' color tables that may include some imported via OPTcontrol.txt

	handles.thematic = true;
	handles.hinge = false;
	if (numel(pal) > 8 && strcmp(pal(1:8),'Imported'))
		pal = pal(1:8);
	end
	
	if ( strncmp(pal, 'Current', 3) )
		pal = handles.cmap_original;
	elseif ( strncmp(pal, 'Imported', 3) )
		pal = handles.imported_cmap;
		handles.thematic = false;
	elseif ( strncmp(pal, 'Mag - anomaly', 3) )
		load([handles.d_path 'gmt_other_palettes.mat'],'mag');			pal = mag;
		handles.z_intervals = [-800 -600 -500 -400 -300 -200 -100 -50 50 100 200 300 400 500 600;
								-600 -500 -400 -300 -200 -100 -50 50 100 200 300 400 500 600 800]';
	elseif ( strncmp(pal, 'Chlorophyll', 3) )			% Not visible (results of this are lousy and I need to find why)
		load([handles.d_path 'gmt_other_palettes.mat'],'rainbow_hist');	pal = rainbow_hist;
		y = (10).^ [-2 + (0:269-2)*(2 + 2)/(floor(269)-1), 2];		% The same as y = logspace(-2, 2, 269);
		handles.z_intervals = [0 y(1:254); y(1:255)]';
	elseif ( strncmp(pal, 'SeaLand (m)', 3) )			% The "y" case is not used (needs to finish/find a clever algo)
		load([handles.d_path 'gmt_other_palettes.mat'],'Terre_Mer');	pal = Terre_Mer;
		y = [linspace(-7000,-1,146) linspace(0,5500,108) 6000]';
		handles.z_intervals = [y(1:end-1) y(2:end)];
		handles.hinge = 147;
	elseif ( strncmp(pal, 'Bathymetry (m)', 3) )		% The "y" case is not used (needs to finish/find a clever algo)
		load([handles.d_path 'caris256.mat'],'Earth');					pal = Earth/255;
		y = [(-6000:1000:-3000) (-2500:500:-500) (-400:100:-100) -50 0]';
		handles.z_intervals = [y(1:end-1) y(2:end)];
	elseif ( strncmp(pal, 'Topography (m)', 3) )
		load([handles.d_path 'caris256.mat'],'Topographic');			pal = Topographic/255;
		handles.z_intervals = [0 50 100 200 500 (1000:1000:5000); 50 100 200 500 (1000:1000:6000)]';
		handles.hinge = 1;
	elseif (strcmp(pal, 'SST (12-26)'))
		pal = jet(256);
		handles.z_intervals = [linspace(12,25.9,140); linspace(12.1,26,140)]';		% 141 = (26-12)/0.1 + 1
	elseif (strcmp(pal, 'SST (0-20)'))
		pal = jet(256);
		handles.z_intervals = [linspace(0,19.9,200); linspace(0.1,20,200)]';		% 201 = (20-0)/0.1 + 1
	elseif (strcmp(pal, 'SST (0-35)'))
		pal = jet(256);
		handles.z_intervals = [linspace(0,34.8,175); linspace(0.2,35,175)]';		% 176 = (35-0)/0.2 + 1
	else
		for (k = 1:numel(handles.custom_thematic_name))		% Won't be executed if custom_thematic_name is empty
			switch pal
				case handles.custom_thematic_name{k}
					pal = handles.custom_thematic_pal{k,1};
					handles.z_intervals = handles.custom_thematic_pal{k,2};
					continue
			end
		end
	end

	% Search eventual color markers and delete them
	hp = findobj(handles.axes1,'Tag','Picos');
	if (~isempty(hp)),		delete(hp);		end
	handles.no_slider = 1;
	guidata(handles.figure1, handles);
	change_cmap(handles,pal);

% -----------------------------------------------------------------------------------
% --- Executes on slider movement.
function slider_Bottom_CB(hObject, handles)
	val = round(get(hObject,'Value'));

	if ~isempty(handles.pal_top)        % The other slider has been activated
		val_t = round(get(handles.slider_Top,'Value'));
		cmap = handles.pal_top;         % Make a copy of the full current colormap
	else
		cmap = handles.cmap;            % Make a copy of the full current colormap
		val_t = length(handles.cmap);
	end

	for i = 1:val,		cmap(i,:) = handles.cmap(1,:);		end
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
% --- Executes on slider movement.
function slider_Top_CB(hObject, handles)
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
function change_cmap(handles, pal)
% Change the Image's colormap to 'pal'
	if strcmp(get(handles.OptionsAutoApply,'checked'),'on') % otherwise we are just playing with color_palettes alone
		del = 0;
		h = findobj(handles.hCallingFig,'Type','image');
		if ( ndims(get(h,'CData')) == 3 )
			h = handles.figure1;
			waitfor(msgbox('True color images do not use color palettes. So you cannot change it.','Warning','warn'))
			del = 1;
		end
		if (del),	delete(h);		return,		end
	end

	if (handles.thematic && ~handles.hinge)			% Imported GMT palette with Z levels, OR ...
		if (isempty(handles.z_min_orig)),	return,		end						% No Grid, not possible to continue
		len_Pal = size(pal, 1);    
		z_grd = linspace(handles.z_min_orig, handles.z_max_orig, len_Pal)';		% Z[min max] descretized in len_Pal intervals
		z_pal = [handles.z_intervals(:,1); handles.z_intervals(end,end)];		% Palette z intervals
		% Interpolate the grid levels into the palette Z levels
		y = (len_Pal - 1) * ((z_pal - z_pal(1)) / (z_pal(end) - z_pal(1)));		% Spread z_pal into the [0 len_Pal] interval
		%ind_pal = interp1(z_pal, y, z_grd, 'linear', 'extrap');
		ind_pal = interp1(z_pal, y, z_grd, 'linear', NaN);
		indNaN  = isnan(ind_pal);			% Because interp1 only has one 'extrapval' we must trick it like this 
		if (any(indNaN))
			ib = find(diff(indNaN) == -1);		ie = find(diff(indNaN) == 1);	% Find the indices of extrapolated vals
			ind_pal(1:ib) = ind_pal(ib+1);		ind_pal(ie+1:end) = ind_pal(ie);% Replace by first/last valid value
			if (isempty(ib) && isempty(ie)),	ind_pal = zeros(size(ind_pal));	end
		end
		ind_pal = round(ind_pal + 0.5);
		ind_pal(ind_pal < 1) = 1;		ind_pal(ind_pal > len_Pal) = len_Pal;	% Can happen easily for |z_grd| > |z_pal|
		% Map the old pal indices into the new ones
		pal = pal(ind_pal,:);
		set(handles.check_logIt,'Val',0)			% Let no confusions about this
	elseif (handles.thematic && handles.hinge)
		pal = makeCmapBat(handles.z_min_orig, handles.z_max_orig, handles.hinge, pal);
	end

	z_grd = [];				% Most common case. No particular scaling or logaritmization
	if ( get(handles.check_logIt,'Val') )			% If logaritmize
		if ( ~isempty(handles.z_min) && (handles.z_min ~= handles.z_min_orig || handles.z_max ~= handles.z_max_orig) )
			z_grd = linspace( handles.z_min,handles.z_max, size(pal,1) )';
		else
			z_grd = linspace( handles.z_min_orig,handles.z_max_orig, size(pal,1) )';
		end
	elseif ( ~isempty(handles.z_min) && (handles.z_min ~= handles.z_min_orig || handles.z_max ~= handles.z_max_orig) )
		z_grd = linspace(handles.z_min_orig, handles.z_max_orig, size(pal,1))';
	end

	if (~isempty(z_grd))		% We have a non data-orig min/max
		len_Pal = length(pal);    
		if ( ~get(handles.check_logIt,'Val') )
			z_pal = linspace(handles.z_min,handles.z_max,len_Pal)';
		else				% Calculate a logarithm cmap
			if (handles.z_min <= 0)		% We don't want to take logs of negative numbers
				z_min = -handles.z_min;		z_max = handles.z_max + (z_min-handles.z_min);
			else
				z_min = handles.z_min;		z_max = handles.z_max;
			end
			log_maxmin = log(z_max / z_min) / len_Pal;
			z_pal = z_min * exp((1:len_Pal) * log_maxmin);
			if (handles.z_min <= 0)
				z_pal = z_pal - (z_min-handles.z_min);		% Reset the shift applyied above to avoid taking log(negative)
			end
		end
		% Interpolate the grid levels into the User selected extrema
		ind_pal = interp1(z_pal,linspace(0,1,len_Pal),z_grd,'linear','extrap');
		ind_pal = round(ind_pal * len_Pal);
		% Map the old pal indices into the new ones
		pal = interp1(linspace(0,len_Pal-1,length(pal)),pal,ind_pal,'linear','extrap');
	end

	if (~isempty(handles.hCallingFig))		% It is when drag-N-drop a .cpt file
		handMir = guidata(handles.hCallingFig);
		if (isfield(handMir, 'hImg') && strcmp(get(handMir.hImg, 'CDataMapping'), 'scaled'))	% Scaled images may have much shorter cmaps
			clim = get(handMir.axes1, 'CLim');
			pal = interp1(linspace(0,1,size(pal,1)), pal, linspace(0,1,diff(clim)+1), 'linear','extrap');
		end
	end

	pal(pal > 1) = 1;			% Sometimes interpolation gives values that are out of [0,1] range...
	pal(pal < 0) = 0;
	
	set(handles.figure1,'Colormap',pal);

	if strcmp(get(handles.OptionsAutoApply,'checked'),'on')
		if (handles.have_nans),     pal_bg = [handles.bg_color; pal(2:end,:)];
		else                        pal_bg = pal;
		end
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
function new_cmap = makeCmapBat(z_min, z_max, hinge, cmap)
% Put the cmap discontinuity at the zero of bathymetry (coastline)

	nc = length(cmap);
	z_inc = (z_max - z_min) / (nc - 1);
	ind_c = round(abs(0 - z_min) / z_inc + 1);

	nl = hinge;		nu = ind_c;
	if (hinge < 3)			% Among other reasons, there seams to be a BUGish in interp1 that wont let do the else case
		new_cmap_l = repmat(cmap(1,:),nu,1);
	else
		new_cmap_l = interp1(linspace(0,1,nl), cmap(1:nl,:), linspace(0,1,nu));
	end
	new_cmap_u = interp1(linspace(0,1,nc-hinge), cmap(hinge+1:nc,:), linspace(0,1,nc-ind_c));
	new_cmap = [new_cmap_l; new_cmap_u];

% --------------------------------------------------------------------
function FileSavePalette_CB(hObject, handles, opt)
% OPT == [] writes the current cmap as a descrete GMT palette, but with 256 colors (in fact, a continuous cpt)
% OPT == 'master_disc' writes the current cmap as a descrete GMT palette with 16 colors
% OPT == 'master_cont' writes the current cmap as a continuous GMT palette with 16 colors
	if (nargin == 2),   opt = [];   end

	% Get directory history
	if (~isempty(handles.hCallingFig))
		hand = guidata(handles.hCallingFig);
	else
		hand = handles;
	end

	[FileName,PathName] = put_or_get_file(hand, {'*.cpt', 'GMT color palette (*.cpt)'},'Select CPT File name','put','.cpt');
	if isequal(FileName,0),		return,		end

	pal = get(handles.figure1,'Colormap');
	if (handles.have_nans),   pal = pal(2:end,:);   end		% Remove the bg color
	if (strcmp(opt,'master_disc') || strcmp(opt,'master_cont'))
		pal_len = 16;
		dz = 1 / pal_len;
		x = linspace(0,1,pal_len)';
		pal = interp1(linspace(0,1,length(pal)),pal,x);
		z_min = 0;
	else			% current cmap as a descrete GMT palette, but with 256 colors
		pal_len = size(pal,1);
		dz = (handles.z_max - handles.z_min) / pal_len;
		z_min = handles.z_min;
		if (isempty(z_min))			% This happens when the option was used but no grid is loaded
			z_min = 0;
			dz = 1 / pal_len;
		end
	end

	tmp = cell(1, max(pal_len+4, 19));
	tmp{1} = '# Color palette exported by Mirone';
	tmp{2} = '# COLOR_MODEL = RGB';
	if (isempty(opt) || strcmp(opt,'master_disc'))   % Current cmap or descrete GMT master palette with 16 colors
		for i=1:pal_len
			cor = round(pal(i,:)*255);
			cor_str = sprintf([num2str(cor(1),'%.1d') '\t' num2str(cor(2),'%.1d') '\t' num2str(cor(3),'%.1d')]);
			z1 = num2str(z_min+dz*(i-1),'%.4f');
			z2 = num2str(z_min+dz*i,'%.4f');
			tmp{i+2} = sprintf([z1 '\t' cor_str '\t' z2 '\t' cor_str]);
		end
	else                % Continuous GMT master palette with 16 colors
		for i=1:15
			cor_s = round(pal(i,:)*255);
			cor_e = round(pal(i+1,:)*255);
			cor_s_str = sprintf([num2str(cor_s(1),'%.1d') '\t' num2str(cor_s(2),'%.1d') '\t' num2str(cor_s(3),'%.1d')]);
			cor_e_str = sprintf([num2str(cor_e(1),'%.1d') '\t' num2str(cor_e(2),'%.1d') '\t' num2str(cor_e(3),'%.1d')]);
			z1 = num2str(z_min+dz*(i-1),'%.4f');
			z2 = num2str(z_min+dz*i,'%.4f');
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
function OptionsDiscretizePalette_CB(hObject, handles, opt)
	if (nargin == 2),   opt = '16';     end
	n_color = str2double(opt);
	pal = get(handles.figure1,'Colormap');
	if (handles.have_nans),   pal = pal(2:end,:);   end     % Remove the bg color

	xi = linspace(0,1,n_color)';
	x = linspace(0,1,length(pal));
	pal = interp1(x,pal,xi);
	% Now we have to replicate the n_colors until we get 256 because uint8 has that many colors
	pal256 = [];    n_rep = 256 / n_color;
	try			% Wrap it in a try-catch for the case it fails
		for (i = 1:n_color)
			pal256 = [pal256; repmat(pal(i,:),n_rep,1)];
		end
	catch
		pal256 = get(handles.figure1,'Colormap');
	end
	change_cmap(handles,pal256)
	set(handles.listbox1,'Value',1)     % Now we have a new "Current cmap"

% --------------------------------------------------------------------
function FileExit_CB(hObject, handles)
    figure1_CloseRequestFcn(handles.figure1, [])

% --------------------------------------------------------------------
function OptionsApply_CB(hObject, handles)
	if ~isempty(handles.hCallingFig)
		cmap = get(handles.figure1,'Colormap');
		if strcmp(get(handles.OptionsAutoApply,'checked'),'on')
			change_cmap(handles,cmap)
		else
			set(handles.OptionsAutoApply,'checked','on')	% temporarly set it to on
			change_cmap(handles,cmap)						% in order that change_cmap
			set(handles.OptionsAutoApply,'checked','off')	% may work.
		end
	else
		msgbox('Apply what?','Chico Clever')
	end

% --------------------------------------------------------------------
function OptionsAutoApply_CB(hObject, handles)
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
function cmap = FileReadPalette_CB(hObject, handles, opt, opt2)
	if (nargin == 2),   opt = [];	end
	if (nargin < 4),	opt2 = [];	end
	if (isempty(opt2))
		if (~isempty(handles.hCallingFig)),		hand = guidata(handles.hCallingFig);
		else									hand = handles;
		end
		[FileName,PathName] = put_or_get_file(hand, ...
			{'*.cpt;*.CPT', 'CPT files (*.cpt,*.CPT)';'*.*', 'All Files (*.*)'},'Select CPT file','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	else
		[PathName,FileName] = fileparts(opt2);
		fname = opt2;
	end
	try
		[bin, n_column] = guess_file(fname);
		if (isempty(bin))
			errordlg(['Error reading file ' fname],'Error');		return
		end
		if (~isempty(opt))  % Use the cpt Z levels as well
			[cmap,handles.z_intervals] = c_cpt2cmap(['-C' fname]);
		else                % Use only the cpt colors
			cmap = c_cpt2cmap(['-C' fname]);
			handles.z_intervals = [];
		end
	catch
		if (n_column ~= 4)
			errordlg('There was an error reading the CPT file.','Error')
			return
		end
		% Assume a 4 columns file with Z r g b. Allow one comment line only
		numeric_data = text_read(fname);
		cmap = numeric_data(:,2:4);
		if (min(cmap(:)) < 0 || max(cmap(:)) > 255)
			errordlg('Bad Z r g b file. Colors outside the [0 255] inerval','Error'),	return
		end
		if (max(cmap(:)) > 1),	cmap = cmap / 255;		end
		idx1 = linspace(1,256,size(cmap,1));
		cmap = interp1(idx1,cmap,1:256);		% We need it as a 256 cmap
		if (~isempty(opt))		% Use the cpt Z levels
			handles.z_intervals = [numeric_data(1:end-1,1) numeric_data(2:end,1)];
		else
			handles.z_intervals = [];
		end
	end
	handles.cmap = cmap;        handles.imported_cmap = cmap;
	handles.no_slider = 1;
	list = get(handles.listbox1,'String');
	list{1} = ['Imported (' FileName ')'];
	set(handles.listbox1,'String',list,'Value',1);
	guidata(handles.figure1,handles)
	%if (~strcmp(opt,'z_levels')),	handles.z_intervals = [];	end	% Trick to avoid cmap recalculation inside change_cmap()
	if ( ~isempty(handles.hCallingFig) && ~isempty(handles.z_intervals) )
		handles.thematic = true;
	end
	change_cmap(handles, cmap);

% --------------------------------------------------------------------
function push_retColorMap_CB(hObject, handles)
    handles.killed = false;
    guidata(hObject, handles);    uiresume(handles.figure1);

% --------------------------------------------------------------------
function radio_ML_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_GMT handles.radio_MR handles.radio_CAR handles.radio_GIMP handles.radio_T],'Value', 0)
	pals = get(handles.listbox1,'String');
	set(handles.listbox1,'String',[pals(1) handles.palsML],'Value',1)

% --------------------------------------------------------------------
function radio_GMT_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_ML handles.radio_MR handles.radio_CAR handles.radio_GIMP handles.radio_T],'Value', 0)
	pals = get(handles.listbox1,'String');
	set(handles.listbox1,'String',[pals(1) handles.palsGMT],'Value',1)

% --------------------------------------------------------------------
function radio_MR_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_ML handles.radio_GMT handles.radio_CAR handles.radio_GIMP handles.radio_T],'Value', 0)
	pals = get(handles.listbox1,'String');
	set(handles.listbox1,'String',[pals(1) handles.palsA],'Value',1)

% --------------------------------------------------------------------
function radio_CAR_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_ML handles.radio_GMT handles.radio_MR handles.radio_GIMP handles.radio_T],'Value', 0)
	pals = get(handles.listbox1,'String');
	set(handles.listbox1,'String',[pals(1) handles.palsCAR],'Value',1)

% --------------------------------------------------------------------
function radio_GIMP_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_ML handles.radio_GMT handles.radio_MR handles.radio_CAR handles.radio_T],'Value', 0)
	pals = get(handles.listbox1,'String');
	set(handles.listbox1,'String',[pals(1) handles.palsGIMP],'Value',1)

% --------------------------------------------------------------------
function radio_T_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_ML handles.radio_GMT  handles.radio_MR handles.radio_CAR handles.radio_GIMP],'Value', 0)
	pals = get(handles.listbox1,'String');
	set(handles.listbox1,'String',[pals(1) handles.palsT],'Value',1)
	drawnow,	pause(0.1)				% Do not wait until the rest of the function finish executing

	if (handles.check_custom_pal)		% Do this only the first time this button is checked
		opt_file = [handles.home_dir filesep 'data' filesep 'OPTcontrol.txt'];
		if ( exist(opt_file, 'file') == 2 )
			fid = fopen(opt_file, 'r');
			c = (fread(fid,'*char'))';      fclose(fid);
			lines = strread(c,'%s','delimiter','\n');   clear c fid;
			m = numel(lines);		kk = 0;	kkk = 0;
			for (k = 1:m)
				if (~strncmp(lines{k},'MIR_CPT',7)),	continue,	end
				[t, r] = strtok(lines{k}(9:end));
				try
					[cmap, z_int] = c_cpt2cmap(['-C' ddewhite(r)]);
					handles.custom_thematic_name{kk+1} = ddewhite(t);
					handles.custom_thematic_pal{kk+1,1} = cmap;		handles.custom_thematic_pal{kk+1,2} = z_int;
					kk = kk + 1;
				catch
					kkk = kkk + 1;
				end
			end
		end
		if (kkk)
			warndlg(['Waring: ' sprintf('%d ',kkk) 'custom CPT files in OPTcontrol.txt failed to load'], 'Warning')
		end
		if (kk)			% We got new CPTs. Add them to the list
			contents = [get(handles.listbox1,'Str'); handles.custom_thematic_name];
			set(handles.listbox1,'Str', contents)
			if (~isempty(handles.hCallingFig))	% Save this so we can load it in a next live thus preventing the cpt2cmap call crash
				setappdata(handles.hCallingFig, 'thematic_pal', {handles.custom_thematic_name; handles.custom_thematic_pal})
			end
		end

		handles.check_custom_pal = false;		% Do not try the above again
		guidata(handles.figure1, handles)
	end

% --------------------------------------------------------------------
function check_logIt_CB(hObject, handles)
	if (get(hObject,'Val'))
		change_cmap(handles,handles.cmap)
		handles.cmap_prev = handles.cmap;
	else
		handles.cmap = handles.cmap_prev;
		set(handles.figure1,'Colormap',handles.cmap);
		if strcmp(get(handles.OptionsAutoApply,'checked'),'on')
			set(handles.hCallingFig,'Colormap',handles.cmap)         % Change the image colormap
		end
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function bdn_pal(obj,eventdata,handles)
	stype = get(handles.figure1,'selectiontype');
	if (strcmp(stype,'open'))
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
		ind_c = round(px - x(1) / (x(2)-x(1)) * length(handles.cmap)) + 0;
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
		ind_c = round(px - x(1) / (x(2)-x(1)) * length(handles.cmap)) + 0;
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
	set(handles.h_txt_cZ,'String',['Z = ' num2str(z_cur,'%.3f')],'Pos',handles.txt_cZ_pos)

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
function Help_CB(hObject, handles)
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
function edit_Zmin_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),     set(hObject,'String',num2str(handles.z_min));   return; end
	handles.z_min = xx;
	guidata(hObject,handles)
	change_cmap(handles,get(handles.figure1,'Colormap'))

% --------------------------------------------------------------------
function edit_Zmax_CB(hObject, handles)
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
extra = m - nrgsteps^3;
if (extra == 0) && (nrgsteps > 2)
	nbsteps = nrgsteps - 1;
	%extra = m - (nrgsteps^2 * nbsteps);
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
map = map(summap ~= 0,:);

% 4: remove pure colors (ones with two elements zero):

summap = [sum(map(:,[1 2]),2) sum(map(:,[2 3]),2) sum(map(:,[1 3]),2)];
map = map(min(summap,[],2) ~= 0,:);

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

% --------------------------------------------------------------------
function cmap = vivid(varargin)
% VIVID Creates a Personalized Vivid Colormap
%  VIVID(...) Creates a vivid colormap with custom settings
%
%   Inputs:
%       M - (optional) an integer between 1 and 256 specifying the number
%           of colors in the colormap. Default is 128.
%       MINMAX - (optional) is a 1x2 vector with values between 0 and 1
%           representing the intensity range for the colors, which correspond
%           to black and white, respectively. Default is [0.15 0.85].
%       CLRS - (optional) either a Nx3 matrix of values between 0 and 1
%           representing the desired colors in the colormap
%               -or-
%           a string of characters that includes any combination of the
%           following letters:
%               'r' = red, 'g' = green, 'b' = blue
%               'y' = yellow, 'c' = cyan, 'm' = magenta
%               'o' = orange, 'l' = lime green, 'a' = aquamarine
%               's' = sky blue, 'v' = violet, 'p' = pink
%               'k' or 'w' = black/white/grayscale
%
%   Outputs:
%       CMAP - an Mx3 colormap matrix
%
%   Example:
%       % Default Colormap
%       imagesc(sort(rand(200),'descend'));
%       colormap(vivid); colorbar
%
%   Example:
%       % Mapping With 256 Colors
%       imagesc(peaks(500))
%       colormap(vivid(256)); colorbar
%
%   Example:
%       % Mapping With Full Intensity Range
%       imagesc(peaks(500))
%       colormap(vivid([0 1])); colorbar
%
%   Example:
%       % Mapping With Light Colors
%       imagesc(peaks(500))
%       colormap(vivid([.5 1])); colorbar
%
%   Example:
%       % Mapping With Dark Colors
%       imagesc(peaks(500))
%       colormap(vivid([0 .5])); colorbar
%
%   Example:
%       % Mapping With Custom Color Matrix
%       imagesc(peaks(500))
%       clrs = [1 0 .5; 1 .5 0; .5 1 0; 0 1 .5; 0 .5 1; .5 0 1];
%       colormap(vivid(clrs)); colorbar
%
%   Example:
%       % Mapping With Color String
%       imagesc(peaks(500))
%       colormap(vivid('roylgacsbvmp')); colorbar
%
%   Example:
%       % Colormap With Multiple Custom Settings
%       imagesc(sort(rand(300,100),'descend'));
%       colormap(vivid(64,[.1 .9],'rwb')); colorbar
%
% See also: jet, hsv, gray, hot, copper, bone
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 1.0
% Date: 07/25/08


	% Default Color Spectrum
	clrs = [1 0 0;1 .5 0;1 1 0;0 1 0    % Red, Orange, Yellow, Green    
		0 1 1;0 0 1;.5 0 1;1 0 1];      % Cyan, Blue, Violet, Magenta

	% Default Min/Max Intensity Range
	minmax = [0.15 0.85];

	% Default Colormap Size
	m = 128;

	% Process Inputs
	for var = varargin
		input = var{1};
		if ischar(input)
			num_clrs = length(input);
			clr_mat = zeros(num_clrs,3);
			c = 0;
			for k = 1:num_clrs
				c = c + 1;
				switch lower(input(k))
					case 'r', clr_mat(c,:) = [1 0 0];  % red
					case 'g', clr_mat(c,:) = [0 1 0];  % green
					case 'b', clr_mat(c,:) = [0 0 1];  % blue
					case 'y', clr_mat(c,:) = [1 1 0];  % yellow
					case 'c', clr_mat(c,:) = [0 1 1];  % cyan
					case 'm', clr_mat(c,:) = [1 0 1];  % magenta
					case 'p', clr_mat(c,:) = [1 0 .5]; % pink
					case 'o', clr_mat(c,:) = [1 .5 0]; % orange
					case 'l', clr_mat(c,:) = [.5 1 0]; % lime green
					case 'a', clr_mat(c,:) = [0 1 .5]; % aquamarine
					case 's', clr_mat(c,:) = [0 .5 1]; % sky blue
					case 'v', clr_mat(c,:) = [.5 0 1]; % violet
					case {'k','w'}, clr_mat(c,:) = [.5 .5 .5]; % grayscale
					otherwise, c = c - 1;
				end
			end
			clr_mat = clr_mat(1:c,:);
			if ~isempty(clr_mat)
				clrs = clr_mat;
			end
		elseif numel(input) == 1
			m = max(1,min(256,round(real(input))));
		elseif size(input,2) == 3
			clrs = input;
		elseif length(input) == 2
			minmax = max(0,min(1,real(input)));
		end
	end

	% Calculate Parameters
	nc = size(clrs,1);  % number of spectrum colors
	ns = ceil(m/nc);    % number of shades per color
	n = nc*ns;
	d = n - m;

	% Scale Intensity
	sup = 2*minmax;
	sub = 2*minmax - 1;
	high = repmat(min(1,linspace(sup(1),sup(2),ns))',[1 nc 3]);
	low = repmat(max(0,linspace(sub(1),sub(2),ns))',[1 nc 3]);

	% Determine Color Spectrum
	rgb = repmat(reshape(flipud(clrs),1,nc,3),ns,1);
	map = rgb.*high + (1-rgb).*low;

	% Obtain Color Map
	cmap = reshape(map,n,3,1);
	cmap(1:ns:d*ns,:) = [];

% ------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
    handles = guidata(hObject);
	if (~handles.IAmOctave)
		do_uiresume = strcmp(get(hObject, 'waitstatus'), 'waiting');
	else
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	end
	if (do_uiresume)	% The GUI is still in UIWAIT, us UIRESUME
		handles.killed = true;      % User gave up, return nothing
		guidata(hObject, handles);    uiresume(hObject);
	else				% The GUI is no longer waiting, just close it
		delete(hObject)
	end

% --------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function color_palettes_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'CloseRequestFcn',{@figure1_CloseRequestFcn},...
'MenuBar','none',...
'Toolbar','none',...
'Name','Color Palettes',...
'NumberTitle','off',...
'Pos',[520 500 300 388],...
'Resize','off',...
'DoubleBuffer','on',...
'Tag','figure1');

uicontrol('Parent',h1, 'Pos',[5 221 21 128], 'Style','frame', 'Tag','frame_top');
uicontrol('Parent',h1, 'Pos',[5 8 285 207], 'Style','frame', 'Tag','frame_bot');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@color_palettes_uiCB,...
'HorizontalAlignment','left',...
'Pos',[12 14 271 161],...
'String','Color Tables',...
'Style','listbox',...
'Value',1,...
'Tag','listbox1');

uicontrol('Parent',h1, 'Pos',[12 274 271 16],...
'BackgroundColor',[0.9 0.9 0.9],...
'Call',@color_palettes_uiCB,...
'Style','slider',...
'Tag','slider_Bottom');

axes('Parent',h1, 'Units','pixels',...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Pos',[12 294 271 50],...
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

uicontrol('Parent',h1, 'Pos',[12 241 271 16],...
'BackgroundColor',[0.9 0.9 0.9],...
'Call',@color_palettes_uiCB,...
'Style','slider',...
'Tag','slider_Top');

uicontrol('Parent',h1, 'Pos',[16 258 75 15],'String','Stretch Bottom','Style','text');
uicontrol('Parent',h1, 'Pos',[14 226 65 14],'String','Stretch Top','Style','text');

h14 = uimenu('Parent',h1,'Label','File');
h15 = uimenu('Parent',h14,'Label','Read GMT palette');

uimenu('Parent',h15,...
'Call',{@color_palettes_uiCB3,[]},...
'Label','As master','Tag','FileReadPalette');

uimenu('Parent',h15,...
'Call',{@color_palettes_uiCB3,'z_levels'},...
'Label','Use Z levels','Tag','FileReadPalette');

h18 = uimenu('Parent',h14,'Label','Save as GMT palette');

uimenu('Parent',h18, 'Call',{@color_palettes_uiCB3,[]},...
'Label','Grid limits with 256 colors','Tag','FileSavePalette');

uimenu('Parent',h18, 'Call',{@color_palettes_uiCB3,'master_disc'},...
'Label','descrete master cpt (16 colors) ','Tag','FileSavePalette');

uimenu('Parent',h18, 'Call',{@color_palettes_uiCB3,'master_cont'},...
'Label','continuous master cpt (16 colors) ','Tag','FileSavePalette');

uimenu('Parent',h14, 'Call',@color_palettes_uiCB,...
'Label','Exit','Separator','on','Tag','FileExit');

h23 = uimenu('Parent',h1,'Label','Options');

uimenu('Parent',h23, 'Call',@color_palettes_uiCB,...
'Label','Apply','Tag','OptionsApply');

uimenu('Parent',h23, 'Call',@color_palettes_uiCB,...
'Checked','on','Label','Auto Apply','Tag','OptionsAutoApply');

h26 = uimenu('Parent',h23,'Label','Discretize Palette','Separator','on');

uimenu('Parent',h26, 'Call',{@color_palettes_uiCB3,'8'},...
'Label','8 colors', 'Tag','OptionsDiscretizePalette');

uimenu('Parent',h26, 'Call',{@color_palettes_uiCB3,'16'},...
'Label','16 colors', 'Tag','OptionsDiscretizePalette');

uimenu('Parent',h26, 'Call',{@color_palettes_uiCB3,'32'},...
'Label','32 colors', 'Tag','OptionsDiscretizePalette');

uimenu('Parent',h26, 'Call',{@color_palettes_uiCB3,'64'},...
'Label','64 colors', 'Tag','OptionsDiscretizePalette');

uicontrol('Parent',h1, 'Call',@color_palettes_uiCB,...
'Pos',[12 196 45 15], 'String','ML', 'Style','radiobutton', 'Tag','radio_ML','Val',1);

uicontrol('Parent',h1, 'Call',@color_palettes_uiCB,...
'Pos',[70 196 55 15], 'String','GMT', 'Style','radiobutton', 'Tag','radio_GMT','Val',0);

uicontrol('Parent',h1, 'Call',@color_palettes_uiCB,...
'Pos',[133 196 45 15], 'String','MR', 'Style','radiobutton', 'Tag','radio_MR','Val',0);

uicontrol('Parent',h1, 'Call',@color_palettes_uiCB,...
'Pos',[12 176 55 15],'String','CAR', 'Style','radiobutton', 'Tag','radio_CAR','Val',0);

uicontrol('Parent',h1, 'Call',@color_palettes_uiCB,...
'Pos',[70 176 60 15],'String','GIMP', 'Style','radiobutton','Tag','radio_GIMP','Val',0);

uicontrol('Parent',h1, 'Call',@color_palettes_uiCB,...
'Pos',[133 176 80 15], 'String','Thematic', 'Style','radiobutton', 'Tag','radio_T','Val',0);

uicontrol('Parent',h1, 'Call',@color_palettes_uiCB,...
'ToolTip','Take a natural logarithm (ln()) of the current color palette', ...
'Pos',[210 186 85 15],'String','Logaritmize', 'Style','checkbox','Tag','check_logIt','Val',0);

uimenu('Parent',h1, 'Call',@color_palettes_uiCB, 'Label','Help', 'Tag','Help');

uicontrol('Parent',h1, 'Pos',[50 357 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@color_palettes_uiCB,...
'Style','edit',...
'Tooltip','Use a different value to set a fixed color minimum',...
'Tag','edit_Zmin');

uicontrol('Parent',h1, 'Pos',[200 357 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@color_palettes_uiCB,...
'Style','edit',...
'Tooltip','Use a different value to set a fixed color miximum',...
'Tag','edit_Zmax');

uicontrol('Parent',h1,'Pos',[160 359 40 15],'String','Max Z','Style','text','Tag','text_MaxZ');
uicontrol('Parent',h1,'Pos',[10 359 40 15],'String','Min Z','Style','text','Tag','text_MinZ');

uicontrol('Parent',h1, 'Pos',[100 333 85 16], 'Visible','off', ...
'FontSize',10, 'HorizontalAlignment','left',...
'Style','text', 'Tag','h_txt_cZ');

uicontrol('Parent',h1, 'Pos',[10 356 201 23],...
'Call',@color_palettes_uiCB,...
'FontName','Helvetica','FontSize',9,...
'String','Finish and return the Color Map',...
'Visible','off',...
'Tag','push_retColorMap');

function color_palettes_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'], hObject, guidata(hObject));

function color_palettes_uiCB3(hObject, eventdata, opt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'], hObject, guidata(hObject),opt);

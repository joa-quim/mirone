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
	palsML = {'autumn' 'bone' 'colorcube' 'cool' 'copper' 'flag' 'gray' 'hot' 'hsv' 'jet' 'jet-improved' 'lines' ...
		'pink' 'prism' 'summer' 'winter' 'vivid' 'parola' 'viridis' 'magma' 'inferno'};
	handles.palsGMT = {'drywet' 'gebco' 'globe' 'rainbow' 'haxby' 'no_green' 'ocean' 'polar' 'red2green' ...
		'sealand' 'seis' 'split' 'topo' 'wysiwyg' 'DEM_screen' 'DEM_print' 'DEM_poster'};
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
	case 'autumn',		pal = autumn(256);
	case 'bone',		pal = bone(256);
	case 'colorcube',	pal = colorcube(256);
	case 'cool',		pal = cool(256);
	case 'copper',		pal = copper(256);
	case 'flag',		pal = flag(256);
	case 'gray',		pal = gray(256);
	case 'hot',			pal = hot(256);
	case 'hsv',			pal = hsv(256);
	case 'jet',			pal = jet(256);
	case 'jet-improved',	pal = mkpj(256);
	case 'lines',		pal = lines(256);
	case 'pink',		pal = pink(256);
	case 'prism',		pal = prism(256);
	case 'summer',		pal = summer(256);
	case 'white',		pal = white(256);
	case 'winter',		pal = winter(256);
	case 'vivid',		pal = vivid(256);
	case 'parola',		pal = parola(256);
	case 'viridis',		pal = viridis(256);
	case 'magma',		pal = magma(256);
	case 'inferno',		pal = inferno(256);

	case 'drywet'
		load([handles.d_path 'gmt_other_palettes.mat'],'drywet');		pal = drywet;
	case 'gebco'
		load([handles.d_path 'gmt_other_palettes.mat'],'gebco');		pal = gebco;
	case 'globe'
		load([handles.d_path 'gmt_other_palettes.mat'],'globo');		pal = globo;
	case 'gray'
		load([handles.d_path 'gmt_other_palettes.mat'],'gray');			pal = gray;
	case 'haxby'
		load([handles.d_path 'gmt_other_palettes.mat'],'haxby');		pal = haxby;
	case 'no_green'
		load([handles.d_path 'gmt_other_palettes.mat'],'no_green');		pal = no_green;
	case 'ocean'
		load([handles.d_path 'gmt_other_palettes.mat'],'ocean');		pal = ocean;
	case 'polar'
		load([handles.d_path 'gmt_other_palettes.mat'],'polar');		pal = polar;
	case 'rainbow'
		load([handles.d_path 'gmt_other_palettes.mat'],'rainbow');		pal = rainbow;
	case 'red2green'
		load([handles.d_path 'gmt_other_palettes.mat'],'red2green');	pal = red2green;
	case 'sealand'
		load([handles.d_path 'gmt_other_palettes.mat'],'sealand');		pal = sealand;
	case 'seis'
		load([handles.d_path 'gmt_other_palettes.mat'],'seis');			pal = seis;
	case 'split'
		load([handles.d_path 'gmt_other_palettes.mat'],'split');		pal = split;
	case 'topo'
		load([handles.d_path 'gmt_other_palettes.mat'],'topo');			pal = topo;
	case 'wysiwyg'
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
function c = parola(m)
c = [
	0.2081 0.1663 0.5292
	0.2091 0.1721 0.5411
	0.2101 0.1779 0.5530
	0.2109 0.1837 0.5650
	0.2116 0.1895 0.5771
	0.2121 0.1954 0.5892
	0.2124 0.2013 0.6013
	0.2125 0.2072 0.6135
	0.2123 0.2132 0.6258
	0.2118 0.2192 0.6381
	0.2111 0.2253 0.6505
	0.2099 0.2315 0.6629
	0.2084 0.2377 0.6753
	0.2063 0.2440 0.6878
	0.2038 0.2503 0.7003
	0.2006 0.2568 0.7129
	0.1968 0.2632 0.7255
	0.1921 0.2698 0.7381
	0.1867 0.2764 0.7507
	0.1802 0.2832 0.7634
	0.1728 0.2902 0.7762
	0.1641 0.2975 0.7890
	0.1541 0.3052 0.8017
	0.1427 0.3132 0.8145
	0.1295 0.3217 0.8269
	0.1147 0.3306 0.8387
	0.0986 0.3397 0.8495
	0.0816 0.3486 0.8588
	0.0646 0.3572 0.8664
	0.0482 0.3651 0.8722
	0.0329 0.3724 0.8765
	0.0213 0.3792 0.8796
	0.0136 0.3853 0.8815
	0.0086 0.3911 0.8827
	0.0060 0.3965 0.8833
	0.0051 0.4017 0.8834
	0.0054 0.4066 0.8831
	0.0067 0.4113 0.8825
	0.0089 0.4159 0.8816
	0.0116 0.4203 0.8805
	0.0148 0.4246 0.8793
	0.0184 0.4288 0.8779
	0.0223 0.4329 0.8763
	0.0264 0.4370 0.8747
	0.0306 0.4410 0.8729
	0.0349 0.4449 0.8711
	0.0394 0.4488 0.8692
	0.0437 0.4526 0.8672
	0.0477 0.4564 0.8652
	0.0514 0.4602 0.8632
	0.0549 0.4640 0.8611
	0.0582 0.4677 0.8589
	0.0612 0.4714 0.8568
	0.0640 0.4751 0.8546
	0.0666 0.4788 0.8525
	0.0689 0.4825 0.8503
	0.0710 0.4862 0.8481
	0.0729 0.4899 0.8460
	0.0746 0.4937 0.8439
	0.0761 0.4974 0.8418
	0.0773 0.5012 0.8398
	0.0782 0.5051 0.8378
	0.0789 0.5089 0.8359
	0.0794 0.5129 0.8341
	0.0795 0.5169 0.8324
	0.0793 0.5210 0.8308
	0.0788 0.5251 0.8293
	0.0778 0.5295 0.8280
	0.0764 0.5339 0.8270
	0.0746 0.5384 0.8261
	0.0724 0.5431 0.8253
	0.0698 0.5479 0.8247
	0.0668 0.5527 0.8243
	0.0636 0.5577 0.8239
	0.0600 0.5627 0.8237
	0.0562 0.5677 0.8234
	0.0523 0.5727 0.8231
	0.0484 0.5777 0.8228
	0.0445 0.5826 0.8223
	0.0408 0.5874 0.8217
	0.0372 0.5922 0.8209
	0.0342 0.5968 0.8198
	0.0317 0.6012 0.8186
	0.0296 0.6055 0.8171
	0.0279 0.6097 0.8154
	0.0265 0.6137 0.8135
	0.0255 0.6176 0.8114
	0.0248 0.6214 0.8091
	0.0243 0.6250 0.8066
	0.0239 0.6285 0.8039
	0.0237 0.6319 0.8010
	0.0235 0.6352 0.7980
	0.0233 0.6384 0.7948
	0.0231 0.6415 0.7916
	0.0230 0.6445 0.7881
	0.0229 0.6474 0.7846
	0.0227 0.6503 0.7810
	0.0227 0.6531 0.7773
	0.0232 0.6558 0.7735
	0.0238 0.6585 0.7696
	0.0246 0.6611 0.7656
	0.0263 0.6637 0.7615
	0.0282 0.6663 0.7574
	0.0306 0.6688 0.7532
	0.0338 0.6712 0.7490
	0.0373 0.6737 0.7446
	0.0418 0.6761 0.7402
	0.0467 0.6784 0.7358
	0.0516 0.6808 0.7313
	0.0574 0.6831 0.7267
	0.0629 0.6854 0.7221
	0.0692 0.6877 0.7173
	0.0755 0.6899 0.7126
	0.0820 0.6921 0.7078
	0.0889 0.6943 0.7029
	0.0956 0.6965 0.6979
	0.1031 0.6986 0.6929
	0.1104 0.7007 0.6878
	0.1180 0.7028 0.6827
	0.1258 0.7049 0.6775
	0.1335 0.7069 0.6723
	0.1418 0.7089 0.6669
	0.1499 0.7109 0.6616
	0.1585 0.7129 0.6561
	0.1671 0.7148 0.6507
	0.1758 0.7168 0.6451
	0.1849 0.7186 0.6395
	0.1938 0.7205 0.6338
	0.2033 0.7223 0.6281
	0.2128 0.7241 0.6223
	0.2224 0.7259 0.6165
	0.2324 0.7275 0.6107
	0.2423 0.7292 0.6048
	0.2527 0.7308 0.5988
	0.2631 0.7324 0.5929
	0.2735 0.7339 0.5869
	0.2845 0.7354 0.5809
	0.2953 0.7368 0.5749
	0.3064 0.7381 0.5689
	0.3177 0.7394 0.5630
	0.3289 0.7406 0.5570
	0.3405 0.7417 0.5512
	0.3520 0.7428 0.5453
	0.3635 0.7438 0.5396
	0.3753 0.7446 0.5339
	0.3869 0.7454 0.5283
	0.3986 0.7461 0.5229
	0.4103 0.7467 0.5175
	0.4218 0.7473 0.5123
	0.4334 0.7477 0.5072
	0.4447 0.7482 0.5021
	0.4561 0.7485 0.4972
	0.4672 0.7487 0.4924
	0.4783 0.7489 0.4877
	0.4892 0.7491 0.4831
	0.5000 0.7491 0.4786
	0.5106 0.7492 0.4741
	0.5212 0.7492 0.4698
	0.5315 0.7491 0.4655
	0.5418 0.7490 0.4613
	0.5519 0.7489 0.4571
	0.5619 0.7487 0.4531
	0.5718 0.7485 0.4490
	0.5816 0.7482 0.4451
	0.5913 0.7479 0.4412
	0.6009 0.7476 0.4374
	0.6103 0.7473 0.4335
	0.6197 0.7469 0.4298
	0.6290 0.7465 0.4261
	0.6382 0.7460 0.4224
	0.6473 0.7456 0.4188
	0.6564 0.7451 0.4152
	0.6653 0.7446 0.4116
	0.6742 0.7441 0.4081
	0.6830 0.7435 0.4046
	0.6918 0.7430 0.4011
	0.7004 0.7424 0.3976
	0.7091 0.7418 0.3942
	0.7176 0.7412 0.3908
	0.7261 0.7405 0.3874
	0.7346 0.7399 0.3840
	0.7430 0.7392 0.3806
	0.7513 0.7385 0.3773
	0.7596 0.7378 0.3739
	0.7679 0.7372 0.3706
	0.7761 0.7364 0.3673
	0.7843 0.7357 0.3639
	0.7924 0.7350 0.3606
	0.8005 0.7343 0.3573
	0.8085 0.7336 0.3539
	0.8166 0.7329 0.3506
	0.8246 0.7322 0.3472
	0.8325 0.7315 0.3438
	0.8405 0.7308 0.3404
	0.8484 0.7301 0.3370
	0.8563 0.7294 0.3336
	0.8642 0.7288 0.3300
	0.8720 0.7282 0.3265
	0.8798 0.7276 0.3229
	0.8877 0.7271 0.3193
	0.8954 0.7266 0.3156
	0.9032 0.7262 0.3117
	0.9110 0.7259 0.3078
	0.9187 0.7256 0.3038
	0.9264 0.7256 0.2996
	0.9341 0.7256 0.2953
	0.9417 0.7259 0.2907
	0.9493 0.7264 0.2859
	0.9567 0.7273 0.2808
	0.9639 0.7285 0.2754
	0.9708 0.7303 0.2696
	0.9773 0.7326 0.2634
	0.9831 0.7355 0.2570
	0.9882 0.7390 0.2504
	0.9922 0.7431 0.2437
	0.9952 0.7476 0.2373
	0.9973 0.7524 0.2310
	0.9986 0.7573 0.2251
	0.9991 0.7624 0.2195
	0.9990 0.7675 0.2141
	0.9985 0.7726 0.2090
	0.9976 0.7778 0.2042
	0.9964 0.7829 0.1995
	0.9950 0.7880 0.1949
	0.9933 0.7931 0.1905
	0.9914 0.7981 0.1863
	0.9894 0.8032 0.1821
	0.9873 0.8083 0.1780
	0.9851 0.8133 0.1740
	0.9828 0.8184 0.1700
	0.9805 0.8235 0.1661
	0.9782 0.8286 0.1622
	0.9759 0.8337 0.1583
	0.9736 0.8389 0.1544
	0.9713 0.8441 0.1505
	0.9692 0.8494 0.1465
	0.9672 0.8548 0.1425
	0.9654 0.8603 0.1385
	0.9638 0.8659 0.1343
	0.9623 0.8716 0.1301
	0.9611 0.8774 0.1258
	0.9600 0.8834 0.1215
	0.9593 0.8895 0.1171
	0.9588 0.8958 0.1126
	0.9586 0.9022 0.1082
	0.9587 0.9088 0.1036
	0.9591 0.9155 0.0990
	0.9599 0.9225 0.0944
	0.9610 0.9296 0.0897
	0.9624 0.9368 0.0850
	0.9641 0.9443 0.0802
	0.9662 0.9518 0.0753
	0.9685 0.9595 0.0703
	0.9710 0.9673 0.0651
	0.9736 0.9752 0.0597
	0.9763 0.9831 0.0538];

	if (m ~= 256)
		s = size(c,1);
		c = interp1(1:size(c,1), c, linspace(1,s,m), 'linear');
	end
	
% ------------------------------------------------------------------------
function c = magma(m)
c = [
	1.46159096e-03   4.66127766e-04   1.38655200e-02
	2.25764007e-03   1.29495431e-03   1.83311461e-02
	3.27943222e-03   2.30452991e-03   2.37083291e-02
	4.51230222e-03   3.49037666e-03   2.99647059e-02
	5.94976987e-03   4.84285000e-03   3.71296695e-02
	7.58798550e-03   6.35613622e-03   4.49730774e-02
	9.42604390e-03   8.02185006e-03   5.28443561e-02
	1.14654337e-02   9.82831486e-03   6.07496380e-02
	1.37075706e-02   1.17705913e-02   6.86665843e-02
	1.61557566e-02   1.38404966e-02   7.66026660e-02
	1.88153670e-02   1.60262753e-02   8.45844897e-02
	2.16919340e-02   1.83201254e-02   9.26101050e-02
	2.47917814e-02   2.07147875e-02   1.00675555e-01
	2.81228154e-02   2.32009284e-02   1.08786954e-01
	3.16955304e-02   2.57651161e-02   1.16964722e-01
	3.55204468e-02   2.83974570e-02   1.25209396e-01
	3.96084872e-02   3.10895652e-02   1.33515085e-01
	4.38295350e-02   3.38299885e-02   1.41886249e-01
	4.80616391e-02   3.66066101e-02   1.50326989e-01
	5.23204388e-02   3.94066020e-02   1.58841025e-01
	5.66148978e-02   4.21598925e-02   1.67445592e-01
	6.09493930e-02   4.47944924e-02   1.76128834e-01
	6.53301801e-02   4.73177796e-02   1.84891506e-01
	6.97637296e-02   4.97264666e-02   1.93735088e-01
	7.42565152e-02   5.20167766e-02   2.02660374e-01
	7.88150034e-02   5.41844801e-02   2.11667355e-01
	8.34456313e-02   5.62249365e-02   2.20755099e-01
	8.81547730e-02   5.81331465e-02   2.29921611e-01
	9.29486914e-02   5.99038167e-02   2.39163669e-01
	9.78334770e-02   6.15314414e-02   2.48476662e-01
	1.02814972e-01   6.30104053e-02   2.57854400e-01
	1.07898679e-01   6.43351102e-02   2.67288933e-01
	1.13094451e-01   6.54920358e-02   2.76783978e-01
	1.18405035e-01   6.64791593e-02   2.86320656e-01
	1.23832651e-01   6.72946449e-02   2.95879431e-01
	1.29380192e-01   6.79349264e-02   3.05442931e-01
	1.35053322e-01   6.83912798e-02   3.14999890e-01
	1.40857952e-01   6.86540710e-02   3.24537640e-01
	1.46785234e-01   6.87382323e-02   3.34011109e-01
	1.52839217e-01   6.86368599e-02   3.43404450e-01
	1.59017511e-01   6.83540225e-02   3.52688028e-01
	1.65308131e-01   6.79108689e-02   3.61816426e-01
	1.71713033e-01   6.73053260e-02   3.70770827e-01
	1.78211730e-01   6.65758073e-02   3.79497161e-01
	1.84800877e-01   6.57324381e-02   3.87972507e-01
	1.91459745e-01   6.48183312e-02   3.96151969e-01
	1.98176877e-01   6.38624166e-02   4.04008953e-01
	2.04934882e-01   6.29066192e-02   4.11514273e-01
	2.11718061e-01   6.19917876e-02   4.18646741e-01
	2.18511590e-01   6.11584918e-02   4.25391816e-01
	2.25302032e-01   6.04451843e-02   4.31741767e-01
	2.32076515e-01   5.98886855e-02   4.37694665e-01
	2.38825991e-01   5.95170384e-02   4.43255999e-01
	2.45543175e-01   5.93524384e-02   4.48435938e-01
	2.52220252e-01   5.94147119e-02   4.53247729e-01
	2.58857304e-01   5.97055998e-02   4.57709924e-01
	2.65446744e-01   6.02368754e-02   4.61840297e-01
	2.71994089e-01   6.09935552e-02   4.65660375e-01
	2.78493300e-01   6.19778136e-02   4.69190328e-01
	2.84951097e-01   6.31676261e-02   4.72450879e-01
	2.91365817e-01   6.45534486e-02   4.75462193e-01
	2.97740413e-01   6.61170432e-02   4.78243482e-01
	3.04080941e-01   6.78353452e-02   4.80811572e-01
	3.10382027e-01   6.97024767e-02   4.83186340e-01
	3.16654235e-01   7.16895272e-02   4.85380429e-01
	3.22899126e-01   7.37819504e-02   4.87408399e-01
	3.29114038e-01   7.59715081e-02   4.89286796e-01
	3.35307503e-01   7.82361045e-02   4.91024144e-01
	3.41481725e-01   8.05635079e-02   4.92631321e-01
	3.47635742e-01   8.29463512e-02   4.94120923e-01
	3.53773161e-01   8.53726329e-02   4.95501096e-01
	3.59897941e-01   8.78311772e-02   4.96778331e-01
	3.66011928e-01   9.03143031e-02   4.97959963e-01
	3.72116205e-01   9.28159917e-02   4.99053326e-01
	3.78210547e-01   9.53322947e-02   5.00066568e-01
	3.84299445e-01   9.78549106e-02   5.01001964e-01
	3.90384361e-01   1.00379466e-01   5.01864236e-01
	3.96466670e-01   1.02902194e-01   5.02657590e-01
	4.02547663e-01   1.05419865e-01   5.03385761e-01
	4.08628505e-01   1.07929771e-01   5.04052118e-01
	4.14708664e-01   1.10431177e-01   5.04661843e-01
	4.20791157e-01   1.12920210e-01   5.05214935e-01
	4.26876965e-01   1.15395258e-01   5.05713602e-01
	4.32967001e-01   1.17854987e-01   5.06159754e-01
	4.39062114e-01   1.20298314e-01   5.06555026e-01
	4.45163096e-01   1.22724371e-01   5.06900806e-01
	4.51270678e-01   1.25132484e-01   5.07198258e-01
	4.57385535e-01   1.27522145e-01   5.07448336e-01
	4.63508291e-01   1.29892998e-01   5.07651812e-01
	4.69639514e-01   1.32244819e-01   5.07809282e-01
	4.75779723e-01   1.34577500e-01   5.07921193e-01
	4.81928997e-01   1.36891390e-01   5.07988509e-01
	4.88088169e-01   1.39186217e-01   5.08010737e-01
	4.94257673e-01   1.41462106e-01   5.07987836e-01
	5.00437834e-01   1.43719323e-01   5.07919772e-01
	5.06628929e-01   1.45958202e-01   5.07806420e-01
	5.12831195e-01   1.48179144e-01   5.07647570e-01
	5.19044825e-01   1.50382611e-01   5.07442938e-01
	5.25269968e-01   1.52569121e-01   5.07192172e-01
	5.31506735e-01   1.54739247e-01   5.06894860e-01
	5.37755194e-01   1.56893613e-01   5.06550538e-01
	5.44015371e-01   1.59032895e-01   5.06158696e-01
	5.50287252e-01   1.61157816e-01   5.05718782e-01
	5.56570783e-01   1.63269149e-01   5.05230210e-01
	5.62865867e-01   1.65367714e-01   5.04692365e-01
	5.69172368e-01   1.67454379e-01   5.04104606e-01
	5.75490107e-01   1.69530062e-01   5.03466273e-01
	5.81818864e-01   1.71595728e-01   5.02776690e-01
	5.88158375e-01   1.73652392e-01   5.02035167e-01
	5.94508337e-01   1.75701122e-01   5.01241011e-01
	6.00868399e-01   1.77743036e-01   5.00393522e-01
	6.07238169e-01   1.79779309e-01   4.99491999e-01
	6.13617209e-01   1.81811170e-01   4.98535746e-01
	6.20005032e-01   1.83839907e-01   4.97524075e-01
	6.26401108e-01   1.85866869e-01   4.96456304e-01
	6.32804854e-01   1.87893468e-01   4.95331769e-01
	6.39215638e-01   1.89921182e-01   4.94149821e-01
	6.45632778e-01   1.91951556e-01   4.92909832e-01
	6.52055535e-01   1.93986210e-01   4.91611196e-01
	6.58483116e-01   1.96026835e-01   4.90253338e-01
	6.64914668e-01   1.98075202e-01   4.88835712e-01
	6.71349279e-01   2.00133166e-01   4.87357807e-01
	6.77785975e-01   2.02202663e-01   4.85819154e-01
	6.84223712e-01   2.04285721e-01   4.84219325e-01
	6.90661380e-01   2.06384461e-01   4.82557941e-01
	6.97097796e-01   2.08501100e-01   4.80834678e-01
	7.03531700e-01   2.10637956e-01   4.79049270e-01
	7.09961888e-01   2.12797337e-01   4.77201121e-01
	7.16387038e-01   2.14981693e-01   4.75289780e-01
	7.22805451e-01   2.17193831e-01   4.73315708e-01
	7.29215521e-01   2.19436516e-01   4.71278924e-01
	7.35615545e-01   2.21712634e-01   4.69179541e-01
	7.42003713e-01   2.24025196e-01   4.67017774e-01
	7.48378107e-01   2.26377345e-01   4.64793954e-01
	7.54736692e-01   2.28772352e-01   4.62508534e-01
	7.61077312e-01   2.31213625e-01   4.60162106e-01
	7.67397681e-01   2.33704708e-01   4.57755411e-01
	7.73695380e-01   2.36249283e-01   4.55289354e-01
	7.79967847e-01   2.38851170e-01   4.52765022e-01
	7.86212372e-01   2.41514325e-01   4.50183695e-01
	7.92426972e-01   2.44242250e-01   4.47543155e-01
	7.98607760e-01   2.47039798e-01   4.44848441e-01
	8.04751511e-01   2.49911350e-01   4.42101615e-01
	8.10854841e-01   2.52861399e-01   4.39304963e-01
	8.16914186e-01   2.55894550e-01   4.36461074e-01
	8.22925797e-01   2.59015505e-01   4.33572874e-01
	8.28885740e-01   2.62229049e-01   4.30643647e-01
	8.34790818e-01   2.65539703e-01   4.27671352e-01
	8.40635680e-01   2.68952874e-01   4.24665620e-01
	8.46415804e-01   2.72473491e-01   4.21631064e-01
	8.52126490e-01   2.76106469e-01   4.18572767e-01
	8.57762870e-01   2.79856666e-01   4.15496319e-01
	8.63320397e-01   2.83729003e-01   4.12402889e-01
	8.68793368e-01   2.87728205e-01   4.09303002e-01
	8.74176342e-01   2.91858679e-01   4.06205397e-01
	8.79463944e-01   2.96124596e-01   4.03118034e-01
	8.84650824e-01   3.00530090e-01   4.00047060e-01
	8.89731418e-01   3.05078817e-01   3.97001559e-01
	8.94700194e-01   3.09773445e-01   3.93994634e-01
	8.99551884e-01   3.14616425e-01   3.91036674e-01
	9.04281297e-01   3.19609981e-01   3.88136889e-01
	9.08883524e-01   3.24755126e-01   3.85308008e-01
	9.13354091e-01   3.30051947e-01   3.82563414e-01
	9.17688852e-01   3.35500068e-01   3.79915138e-01
	9.21884187e-01   3.41098112e-01   3.77375977e-01
	9.25937102e-01   3.46843685e-01   3.74959077e-01
	9.29845090e-01   3.52733817e-01   3.72676513e-01
	9.33606454e-01   3.58764377e-01   3.70540883e-01
	9.37220874e-01   3.64929312e-01   3.68566525e-01
	9.40687443e-01   3.71224168e-01   3.66761699e-01
	9.44006448e-01   3.77642889e-01   3.65136328e-01
	9.47179528e-01   3.84177874e-01   3.63701130e-01
	9.50210150e-01   3.90819546e-01   3.62467694e-01
	9.53099077e-01   3.97562894e-01   3.61438431e-01
	9.55849237e-01   4.04400213e-01   3.60619076e-01
	9.58464079e-01   4.11323666e-01   3.60014232e-01
	9.60949221e-01   4.18323245e-01   3.59629789e-01
	9.63310281e-01   4.25389724e-01   3.59469020e-01
	9.65549351e-01   4.32518707e-01   3.59529151e-01
	9.67671128e-01   4.39702976e-01   3.59810172e-01
	9.69680441e-01   4.46935635e-01   3.60311120e-01
	9.71582181e-01   4.54210170e-01   3.61030156e-01
	9.73381238e-01   4.61520484e-01   3.61964652e-01
	9.75082439e-01   4.68860936e-01   3.63111292e-01
	9.76690494e-01   4.76226350e-01   3.64466162e-01
	9.78209957e-01   4.83612031e-01   3.66024854e-01
	9.79645181e-01   4.91013764e-01   3.67782559e-01
	9.81000291e-01   4.98427800e-01   3.69734157e-01
	9.82279159e-01   5.05850848e-01   3.71874301e-01
	9.83485387e-01   5.13280054e-01   3.74197501e-01
	9.84622298e-01   5.20712972e-01   3.76698186e-01
	9.85692925e-01   5.28147545e-01   3.79370774e-01
	9.86700017e-01   5.35582070e-01   3.82209724e-01
	9.87646038e-01   5.43015173e-01   3.85209578e-01
	9.88533173e-01   5.50445778e-01   3.88365009e-01
	9.89363341e-01   5.57873075e-01   3.91670846e-01
	9.90138201e-01   5.65296495e-01   3.95122099e-01
	9.90871208e-01   5.72706259e-01   3.98713971e-01
	9.91558165e-01   5.80106828e-01   4.02441058e-01
	9.92195728e-01   5.87501706e-01   4.06298792e-01
	9.92784669e-01   5.94891088e-01   4.10282976e-01
	9.93325561e-01   6.02275297e-01   4.14389658e-01
	9.93834412e-01   6.09643540e-01   4.18613221e-01
	9.94308514e-01   6.16998953e-01   4.22949672e-01
	9.94737698e-01   6.24349657e-01   4.27396771e-01
	9.95121854e-01   6.31696376e-01   4.31951492e-01
	9.95480469e-01   6.39026596e-01   4.36607159e-01
	9.95809924e-01   6.46343897e-01   4.41360951e-01
	9.96095703e-01   6.53658756e-01   4.46213021e-01
	9.96341406e-01   6.60969379e-01   4.51160201e-01
	9.96579803e-01   6.68255621e-01   4.56191814e-01
	9.96774784e-01   6.75541484e-01   4.61314158e-01
	9.96925427e-01   6.82827953e-01   4.66525689e-01
	9.97077185e-01   6.90087897e-01   4.71811461e-01
	9.97186253e-01   6.97348991e-01   4.77181727e-01
	9.97253982e-01   7.04610791e-01   4.82634651e-01
	9.97325180e-01   7.11847714e-01   4.88154375e-01
	9.97350983e-01   7.19089119e-01   4.93754665e-01
	9.97350583e-01   7.26324415e-01   4.99427972e-01
	9.97341259e-01   7.33544671e-01   5.05166839e-01
	9.97284689e-01   7.40771893e-01   5.10983331e-01
	9.97228367e-01   7.47980563e-01   5.16859378e-01
	9.97138480e-01   7.55189852e-01   5.22805996e-01
	9.97019342e-01   7.62397883e-01   5.28820775e-01
	9.96898254e-01   7.69590975e-01   5.34892341e-01
	9.96726862e-01   7.76794860e-01   5.41038571e-01
	9.96570645e-01   7.83976508e-01   5.47232992e-01
	9.96369065e-01   7.91167346e-01   5.53498939e-01
	9.96162309e-01   7.98347709e-01   5.59819643e-01
	9.95932448e-01   8.05527126e-01   5.66201824e-01
	9.95680107e-01   8.12705773e-01   5.72644795e-01
	9.95423973e-01   8.19875302e-01   5.79140130e-01
	9.95131288e-01   8.27051773e-01   5.85701463e-01
	9.94851089e-01   8.34212826e-01   5.92307093e-01
	9.94523666e-01   8.41386618e-01   5.98982818e-01
	9.94221900e-01   8.48540474e-01   6.05695903e-01
	9.93865767e-01   8.55711038e-01   6.12481798e-01
	9.93545285e-01   8.62858846e-01   6.19299300e-01
	9.93169558e-01   8.70024467e-01   6.26189463e-01
	9.92830963e-01   8.77168404e-01   6.33109148e-01
	9.92439881e-01   8.84329694e-01   6.40099465e-01
	9.92089454e-01   8.91469549e-01   6.47116021e-01
	9.91687744e-01   8.98627050e-01   6.54201544e-01
	9.91331929e-01   9.05762748e-01   6.61308839e-01
	9.90929685e-01   9.12915010e-01   6.68481201e-01
	9.90569914e-01   9.20048699e-01   6.75674592e-01
	9.90174637e-01   9.27195612e-01   6.82925602e-01
	9.89814839e-01   9.34328540e-01   6.90198194e-01
	9.89433736e-01   9.41470354e-01   6.97518628e-01
	9.89077438e-01   9.48604077e-01   7.04862519e-01
	9.88717064e-01   9.55741520e-01   7.12242232e-01
	9.88367028e-01   9.62878026e-01   7.19648627e-01
	9.88032885e-01   9.70012413e-01   7.27076773e-01
	9.87690702e-01   9.77154231e-01   7.34536205e-01
	9.87386827e-01   9.84287561e-01   7.42001547e-01
	9.87052509e-01   9.91437853e-01   7.49504188e-01];

	if (m ~= 256)
		s = size(c,1);
		c = interp1(1:size(c,1), c, linspace(1,s,m), 'linear');
	end
	
% ------------------------------------------------------------------------
function c = inferno(m)
c = [
	1.46159096e-03   4.66127766e-04   1.38655200e-02
	2.26726368e-03   1.26992553e-03   1.85703520e-02
	3.29899092e-03   2.24934863e-03   2.42390508e-02
	4.54690615e-03   3.39180156e-03   3.09092475e-02
	6.00552565e-03   4.69194561e-03   3.85578980e-02
	7.67578856e-03   6.13611626e-03   4.68360336e-02
	9.56051094e-03   7.71344131e-03   5.51430756e-02
	1.16634769e-02   9.41675403e-03   6.34598080e-02
	1.39950388e-02   1.12247138e-02   7.18616890e-02
	1.65605595e-02   1.31362262e-02   8.02817951e-02
	1.93732295e-02   1.51325789e-02   8.87668094e-02
	2.24468865e-02   1.71991484e-02   9.73274383e-02
	2.57927373e-02   1.93306298e-02   1.05929835e-01
	2.94324251e-02   2.15030771e-02   1.14621328e-01
	3.33852235e-02   2.37024271e-02   1.23397286e-01
	3.76684211e-02   2.59207864e-02   1.32232108e-01
	4.22525554e-02   2.81385015e-02   1.41140519e-01
	4.69146287e-02   3.03236129e-02   1.50163867e-01
	5.16437624e-02   3.24736172e-02   1.59254277e-01
	5.64491009e-02   3.45691867e-02   1.68413539e-01
	6.13397200e-02   3.65900213e-02   1.77642172e-01
	6.63312620e-02   3.85036268e-02   1.86961588e-01
	7.14289181e-02   4.02939095e-02   1.96353558e-01
	7.66367560e-02   4.19053329e-02   2.05798788e-01
	8.19620773e-02   4.33278666e-02   2.15289113e-01
	8.74113897e-02   4.45561662e-02   2.24813479e-01
	9.29901526e-02   4.55829503e-02   2.34357604e-01
	9.87024972e-02   4.64018731e-02   2.43903700e-01
	1.04550936e-01   4.70080541e-02   2.53430300e-01
	1.10536084e-01   4.73986708e-02   2.62912235e-01
	1.16656423e-01   4.75735920e-02   2.72320803e-01
	1.22908126e-01   4.75360183e-02   2.81624170e-01
	1.29284984e-01   4.72930838e-02   2.90788012e-01
	1.35778450e-01   4.68563678e-02   2.99776404e-01
	1.42377819e-01   4.62422566e-02   3.08552910e-01
	1.49072957e-01   4.54676444e-02   3.17085139e-01
	1.55849711e-01   4.45588056e-02   3.25338414e-01
	1.62688939e-01   4.35542881e-02   3.33276678e-01
	1.69575148e-01   4.24893149e-02   3.40874188e-01
	1.76493202e-01   4.14017089e-02   3.48110606e-01
	1.83428775e-01   4.03288858e-02   3.54971391e-01
	1.90367453e-01   3.93088888e-02   3.61446945e-01
	1.97297425e-01   3.84001825e-02   3.67534629e-01
	2.04209298e-01   3.76322609e-02   3.73237557e-01
	2.11095463e-01   3.70296488e-02   3.78563264e-01
	2.17948648e-01   3.66146049e-02   3.83522415e-01
	2.24762908e-01   3.64049901e-02   3.88128944e-01
	2.31538148e-01   3.64052511e-02   3.92400150e-01
	2.38272961e-01   3.66209949e-02   3.96353388e-01
	2.44966911e-01   3.70545017e-02   4.00006615e-01
	2.51620354e-01   3.77052832e-02   4.03377897e-01
	2.58234265e-01   3.85706153e-02   4.06485031e-01
	2.64809649e-01   3.96468666e-02   4.09345373e-01
	2.71346664e-01   4.09215821e-02   4.11976086e-01
	2.77849829e-01   4.23528741e-02   4.14392106e-01
	2.84321318e-01   4.39325787e-02   4.16607861e-01
	2.90763373e-01   4.56437598e-02   4.18636756e-01
	2.97178251e-01   4.74700293e-02   4.20491164e-01
	3.03568182e-01   4.93958927e-02   4.22182449e-01
	3.09935342e-01   5.14069729e-02   4.23720999e-01
	3.16281835e-01   5.34901321e-02   4.25116277e-01
	3.22609671e-01   5.56335178e-02   4.26376869e-01
	3.28920763e-01   5.78265505e-02   4.27510546e-01
	3.35216916e-01   6.00598734e-02   4.28524320e-01
	3.41499828e-01   6.23252772e-02   4.29424503e-01
	3.47771086e-01   6.46156100e-02   4.30216765e-01
	3.54032169e-01   6.69246832e-02   4.30906186e-01
	3.60284449e-01   6.92471753e-02   4.31497309e-01
	3.66529195e-01   7.15785403e-02   4.31994185e-01
	3.72767575e-01   7.39149211e-02   4.32400419e-01
	3.79000659e-01   7.62530701e-02   4.32719214e-01
	3.85228383e-01   7.85914864e-02   4.32954973e-01
	3.91452659e-01   8.09267058e-02   4.33108763e-01
	3.97674379e-01   8.32568129e-02   4.33182647e-01
	4.03894278e-01   8.55803445e-02   4.33178526e-01
	4.10113015e-01   8.78961593e-02   4.33098056e-01
	4.16331169e-01   9.02033992e-02   4.32942678e-01
	4.22549249e-01   9.25014543e-02   4.32713635e-01
	4.28767696e-01   9.47899342e-02   4.32411996e-01
	4.34986885e-01   9.70686417e-02   4.32038673e-01
	4.41207124e-01   9.93375510e-02   4.31594438e-01
	4.47428382e-01   1.01597079e-01   4.31080497e-01
	4.53650614e-01   1.03847716e-01   4.30497898e-01
	4.59874623e-01   1.06089165e-01   4.29845789e-01
	4.66100494e-01   1.08321923e-01   4.29124507e-01
	4.72328255e-01   1.10546584e-01   4.28334320e-01
	4.78557889e-01   1.12763831e-01   4.27475431e-01
	4.84789325e-01   1.14974430e-01   4.26547991e-01
	4.91022448e-01   1.17179219e-01   4.25552106e-01
	4.97257069e-01   1.19379132e-01   4.24487908e-01
	5.03492698e-01   1.21575414e-01   4.23356110e-01
	5.09729541e-01   1.23768654e-01   4.22155676e-01
	5.15967304e-01   1.25959947e-01   4.20886594e-01
	5.22205646e-01   1.28150439e-01   4.19548848e-01
	5.28444192e-01   1.30341324e-01   4.18142411e-01
	5.34682523e-01   1.32533845e-01   4.16667258e-01
	5.40920186e-01   1.34729286e-01   4.15123366e-01
	5.47156706e-01   1.36928959e-01   4.13510662e-01
	5.53391649e-01   1.39134147e-01   4.11828882e-01
	5.59624442e-01   1.41346265e-01   4.10078028e-01
	5.65854477e-01   1.43566769e-01   4.08258132e-01
	5.72081108e-01   1.45797150e-01   4.06369246e-01
	5.78303656e-01   1.48038934e-01   4.04411444e-01
	5.84521407e-01   1.50293679e-01   4.02384829e-01
	5.90733615e-01   1.52562977e-01   4.00289528e-01
	5.96939751e-01   1.54848232e-01   3.98124897e-01
	6.03138930e-01   1.57151161e-01   3.95891308e-01
	6.09330184e-01   1.59473549e-01   3.93589349e-01
	6.15512627e-01   1.61817111e-01   3.91219295e-01
	6.21685340e-01   1.64183582e-01   3.88781456e-01
	6.27847374e-01   1.66574724e-01   3.86276180e-01
	6.33997746e-01   1.68992314e-01   3.83703854e-01
	6.40135447e-01   1.71438150e-01   3.81064906e-01
	6.46259648e-01   1.73913876e-01   3.78358969e-01
	6.52369348e-01   1.76421271e-01   3.75586209e-01
	6.58463166e-01   1.78962399e-01   3.72748214e-01
	6.64539964e-01   1.81539111e-01   3.69845599e-01
	6.70598572e-01   1.84153268e-01   3.66879025e-01
	6.76637795e-01   1.86806728e-01   3.63849195e-01
	6.82656407e-01   1.89501352e-01   3.60756856e-01
	6.88653158e-01   1.92238994e-01   3.57602797e-01
	6.94626769e-01   1.95021500e-01   3.54387853e-01
	7.00575937e-01   1.97850703e-01   3.51112900e-01
	7.06499709e-01   2.00728196e-01   3.47776863e-01
	7.12396345e-01   2.03656029e-01   3.44382594e-01
	7.18264447e-01   2.06635993e-01   3.40931208e-01
	7.24102613e-01   2.09669834e-01   3.37423766e-01
	7.29909422e-01   2.12759270e-01   3.33861367e-01
	7.35683432e-01   2.15905976e-01   3.30245147e-01
	7.41423185e-01   2.19111589e-01   3.26576275e-01
	7.47127207e-01   2.22377697e-01   3.22855952e-01
	7.52794009e-01   2.25705837e-01   3.19085410e-01
	7.58422090e-01   2.29097492e-01   3.15265910e-01
	7.64009940e-01   2.32554083e-01   3.11398734e-01
	7.69556038e-01   2.36076967e-01   3.07485188e-01
	7.75058888e-01   2.39667435e-01   3.03526312e-01
	7.80517023e-01   2.43326720e-01   2.99522665e-01
	7.85928794e-01   2.47055968e-01   2.95476756e-01
	7.91292674e-01   2.50856232e-01   2.91389943e-01
	7.96607144e-01   2.54728485e-01   2.87263585e-01
	8.01870689e-01   2.58673610e-01   2.83099033e-01
	8.07081807e-01   2.62692401e-01   2.78897629e-01
	8.12239008e-01   2.66785558e-01   2.74660698e-01
	8.17340818e-01   2.70953688e-01   2.70389545e-01
	8.22385784e-01   2.75197300e-01   2.66085445e-01
	8.27372474e-01   2.79516805e-01   2.61749643e-01
	8.32299481e-01   2.83912516e-01   2.57383341e-01
	8.37165425e-01   2.88384647e-01   2.52987700e-01
	8.41968959e-01   2.92933312e-01   2.48563825e-01
	8.46708768e-01   2.97558528e-01   2.44112767e-01
	8.51383572e-01   3.02260213e-01   2.39635512e-01
	8.55992130e-01   3.07038188e-01   2.35132978e-01
	8.60533241e-01   3.11892183e-01   2.30606009e-01
	8.65005747e-01   3.16821833e-01   2.26055368e-01
	8.69408534e-01   3.21826685e-01   2.21481734e-01
	8.73740530e-01   3.26906201e-01   2.16885699e-01
	8.78000715e-01   3.32059760e-01   2.12267762e-01
	8.82188112e-01   3.37286663e-01   2.07628326e-01
	8.86301795e-01   3.42586137e-01   2.02967696e-01
	8.90340885e-01   3.47957340e-01   1.98286080e-01
	8.94304553e-01   3.53399363e-01   1.93583583e-01
	8.98192017e-01   3.58911240e-01   1.88860212e-01
	9.02002544e-01   3.64491949e-01   1.84115876e-01
	9.05735448e-01   3.70140419e-01   1.79350388e-01
	9.09390090e-01   3.75855533e-01   1.74563472e-01
	9.12965874e-01   3.81636138e-01   1.69754764e-01
	9.16462251e-01   3.87481044e-01   1.64923826e-01
	9.19878710e-01   3.93389034e-01   1.60070152e-01
	9.23214783e-01   3.99358867e-01   1.55193185e-01
	9.26470039e-01   4.05389282e-01   1.50292329e-01
	9.29644083e-01   4.11479007e-01   1.45366973e-01
	9.32736555e-01   4.17626756e-01   1.40416519e-01
	9.35747126e-01   4.23831237e-01   1.35440416e-01
	9.38675494e-01   4.30091162e-01   1.30438175e-01
	9.41521384e-01   4.36405243e-01   1.25409440e-01
	9.44284543e-01   4.42772199e-01   1.20354038e-01
	9.46964741e-01   4.49190757e-01   1.15272059e-01
	9.49561766e-01   4.55659658e-01   1.10163947e-01
	9.52075421e-01   4.62177656e-01   1.05030614e-01
	9.54505523e-01   4.68743522e-01   9.98735931e-02
	9.56851903e-01   4.75356048e-01   9.46952268e-02
	9.59114397e-01   4.82014044e-01   8.94989073e-02
	9.61292850e-01   4.88716345e-01   8.42893891e-02
	9.63387110e-01   4.95461806e-01   7.90731907e-02
	9.65397031e-01   5.02249309e-01   7.38591143e-02
	9.67322465e-01   5.09077761e-01   6.86589199e-02
	9.69163264e-01   5.15946092e-01   6.34881971e-02
	9.70919277e-01   5.22853259e-01   5.83674890e-02
	9.72590351e-01   5.29798246e-01   5.33237243e-02
	9.74176327e-01   5.36780059e-01   4.83920090e-02
	9.75677038e-01   5.43797733e-01   4.36177922e-02
	9.77092313e-01   5.50850323e-01   3.90500131e-02
	9.78421971e-01   5.57936911e-01   3.49306227e-02
	9.79665824e-01   5.65056600e-01   3.14091591e-02
	9.80823673e-01   5.72208516e-01   2.85075931e-02
	9.81895311e-01   5.79391803e-01   2.62497353e-02
	9.82880522e-01   5.86605627e-01   2.46613416e-02
	9.83779081e-01   5.93849168e-01   2.37702263e-02
	9.84590755e-01   6.01121626e-01   2.36063833e-02
	9.85315301e-01   6.08422211e-01   2.42021174e-02
	9.85952471e-01   6.15750147e-01   2.55921853e-02
	9.86502013e-01   6.23104667e-01   2.78139496e-02
	9.86963670e-01   6.30485011e-01   3.09075459e-02
	9.87337182e-01   6.37890424e-01   3.49160639e-02
	9.87622296e-01   6.45320152e-01   3.98857472e-02
	9.87818759e-01   6.52773439e-01   4.55808037e-02
	9.87926330e-01   6.60249526e-01   5.17503867e-02
	9.87944783e-01   6.67747641e-01   5.83286889e-02
	9.87873910e-01   6.75267000e-01   6.52570167e-02
	9.87713535e-01   6.82806802e-01   7.24892330e-02
	9.87463516e-01   6.90366218e-01   7.99897176e-02
	9.87123759e-01   6.97944391e-01   8.77314215e-02
	9.86694229e-01   7.05540424e-01   9.56941797e-02
	9.86174970e-01   7.13153375e-01   1.03863324e-01
	9.85565739e-01   7.20782460e-01   1.12228756e-01
	9.84865203e-01   7.28427497e-01   1.20784651e-01
	9.84075129e-01   7.36086521e-01   1.29526579e-01
	9.83195992e-01   7.43758326e-01   1.38453063e-01
	9.82228463e-01   7.51441596e-01   1.47564573e-01
	9.81173457e-01   7.59134892e-01   1.56863224e-01
	9.80032178e-01   7.66836624e-01   1.66352544e-01
	9.78806183e-01   7.74545028e-01   1.76037298e-01
	9.77497453e-01   7.82258138e-01   1.85923357e-01
	9.76108474e-01   7.89973753e-01   1.96017589e-01
	9.74637842e-01   7.97691563e-01   2.06331925e-01
	9.73087939e-01   8.05409333e-01   2.16876839e-01
	9.71467822e-01   8.13121725e-01   2.27658046e-01
	9.69783146e-01   8.20825143e-01   2.38685942e-01
	9.68040817e-01   8.28515491e-01   2.49971582e-01
	9.66242589e-01   8.36190976e-01   2.61533898e-01
	9.64393924e-01   8.43848069e-01   2.73391112e-01
	9.62516656e-01   8.51476340e-01   2.85545675e-01
	9.60625545e-01   8.59068716e-01   2.98010219e-01
	9.58720088e-01   8.66624355e-01   3.10820466e-01
	9.56834075e-01   8.74128569e-01   3.23973947e-01
	9.54997177e-01   8.81568926e-01   3.37475479e-01
	9.53215092e-01   8.88942277e-01   3.51368713e-01
	9.51546225e-01   8.96225909e-01   3.65627005e-01
	9.50018481e-01   9.03409063e-01   3.80271225e-01
	9.48683391e-01   9.10472964e-01   3.95289169e-01
	9.47594362e-01   9.17399053e-01   4.10665194e-01
	9.46809163e-01   9.24168246e-01   4.26373236e-01
	9.46391536e-01   9.30760752e-01   4.42367495e-01
	9.46402951e-01   9.37158971e-01   4.58591507e-01
	9.46902568e-01   9.43347775e-01   4.74969778e-01
	9.47936825e-01   9.49317522e-01   4.91426053e-01
	9.49544830e-01   9.55062900e-01   5.07859649e-01
	9.51740304e-01   9.60586693e-01   5.24203026e-01
	9.54529281e-01   9.65895868e-01   5.40360752e-01
	9.57896053e-01   9.71003330e-01   5.56275090e-01
	9.61812020e-01   9.75924241e-01   5.71925382e-01
	9.66248822e-01   9.80678193e-01   5.87205773e-01
	9.71161622e-01   9.85282161e-01   6.02154330e-01
	9.76510983e-01   9.89753437e-01   6.16760413e-01
	9.82257307e-01   9.94108844e-01   6.31017009e-01
	9.88362068e-01   9.98364143e-01   6.44924005e-01];

	if (m ~= 256)
		s = size(c,1);
		c = interp1(1:size(c,1), c, linspace(1,s,m), 'linear');
	end
	
% ------------------------------------------------------------------------
function c = viridis(m)
c = [
	0.26700401  0.00487433  0.32941519
	0.26851048  0.00960483  0.33542652
	0.26994384  0.01462494  0.34137895
	0.27130489  0.01994186  0.34726862
	0.27259384  0.02556309  0.35309303
	0.27380934  0.03149748  0.35885256
	0.27495242  0.03775181  0.36454323
	0.27602238  0.04416723  0.37016418
	0.2770184   0.05034437  0.37571452
	0.27794143  0.05632444  0.38119074
	0.27879067  0.06214536  0.38659204
	0.2795655   0.06783587  0.39191723
	0.28026658  0.07341724  0.39716349
	0.28089358  0.07890703  0.40232944
	0.28144581  0.0843197   0.40741404
	0.28192358  0.08966622  0.41241521
	0.28232739  0.09495545  0.41733086
	0.28265633  0.10019576  0.42216032
	0.28291049  0.10539345  0.42690202
	0.28309095  0.11055307  0.43155375
	0.28319704  0.11567966  0.43611482
	0.28322882  0.12077701  0.44058404
	0.28318684  0.12584799  0.44496   
	0.283072    0.13089477  0.44924127
	0.28288389  0.13592005  0.45342734
	0.28262297  0.14092556  0.45751726
	0.28229037  0.14591233  0.46150995
	0.28188676  0.15088147  0.46540474
	0.28141228  0.15583425  0.46920128
	0.28086773  0.16077132  0.47289909
	0.28025468  0.16569272  0.47649762
	0.27957399  0.17059884  0.47999675
	0.27882618  0.1754902   0.48339654
	0.27801236  0.18036684  0.48669702
	0.27713437  0.18522836  0.48989831
	0.27619376  0.19007447  0.49300074
	0.27519116  0.1949054   0.49600488
	0.27412802  0.19972086  0.49891131
	0.27300596  0.20452049  0.50172076
	0.27182812  0.20930306  0.50443413
	0.27059473  0.21406899  0.50705243
	0.26930756  0.21881782  0.50957678
	0.26796846  0.22354911  0.5120084 
	0.26657984  0.2282621   0.5143487 
	0.2651445   0.23295593  0.5165993 
	0.2636632   0.23763078  0.51876163
	0.26213801  0.24228619  0.52083736
	0.26057103  0.2469217   0.52282822
	0.25896451  0.25153685  0.52473609
	0.25732244  0.2561304   0.52656332
	0.25564519  0.26070284  0.52831152
	0.25393498  0.26525384  0.52998273
	0.25219404  0.26978306  0.53157905
	0.25042462  0.27429024  0.53310261
	0.24862899  0.27877509  0.53455561
	0.2468114   0.28323662  0.53594093
	0.24497208  0.28767547  0.53726018
	0.24311324  0.29209154  0.53851561
	0.24123708  0.29648471  0.53970946
	0.23934575  0.30085494  0.54084398
	0.23744138  0.30520222  0.5419214 
	0.23552606  0.30952657  0.54294396
	0.23360277  0.31382773  0.54391424
	0.2316735   0.3181058   0.54483444
	0.22973926  0.32236127  0.54570633
	0.22780192  0.32659432  0.546532  
	0.2258633   0.33080515  0.54731353
	0.22392515  0.334994    0.54805291
	0.22198915  0.33916114  0.54875211
	0.22005691  0.34330688  0.54941304
	0.21812995  0.34743154  0.55003755
	0.21620971  0.35153548  0.55062743
	0.21429757  0.35561907  0.5511844 
	0.21239477  0.35968273  0.55171011
	0.2105031   0.36372671  0.55220646
	0.20862342  0.36775151  0.55267486
	0.20675628  0.37175775  0.55311653
	0.20490257  0.37574589  0.55353282
	0.20306309  0.37971644  0.55392505
	0.20123854  0.38366989  0.55429441
	0.1994295   0.38760678  0.55464205
	0.1976365   0.39152762  0.55496905
	0.19585993  0.39543297  0.55527637
	0.19410009  0.39932336  0.55556494
	0.19235719  0.40319934  0.55583559
	0.19063135  0.40706148  0.55608907
	0.18892259  0.41091033  0.55632606
	0.18723083  0.41474645  0.55654717
	0.18555593  0.4185704   0.55675292
	0.18389763  0.42238275  0.55694377
	0.18225561  0.42618405  0.5571201 
	0.18062949  0.42997486  0.55728221
	0.17901879  0.43375572  0.55743035
	0.17742298  0.4375272   0.55756466
	0.17584148  0.44128981  0.55768526
	0.17427363  0.4450441   0.55779216
	0.17271876  0.4487906   0.55788532
	0.17117615  0.4525298   0.55796464
	0.16964573  0.45626209  0.55803034
	0.16812641  0.45998802  0.55808199
	0.1666171   0.46370813  0.55811913
	0.16511703  0.4674229   0.55814141
	0.16362543  0.47113278  0.55814842
	0.16214155  0.47483821  0.55813967
	0.16066467  0.47853961  0.55811466
	0.15919413  0.4822374   0.5580728 
	0.15772933  0.48593197  0.55801347
	0.15626973  0.4896237   0.557936  
	0.15481488  0.49331293  0.55783967
	0.15336445  0.49700003  0.55772371
	0.1519182   0.50068529  0.55758733
	0.15047605  0.50436904  0.55742968
	0.14903918  0.50805136  0.5572505 
	0.14760731  0.51173263  0.55704861
	0.14618026  0.51541316  0.55682271
	0.14475863  0.51909319  0.55657181
	0.14334327  0.52277292  0.55629491
	0.14193527  0.52645254  0.55599097
	0.14053599  0.53013219  0.55565893
	0.13914708  0.53381201  0.55529773
	0.13777048  0.53749213  0.55490625
	0.1364085   0.54117264  0.55448339
	0.13506561  0.54485335  0.55402906
	0.13374299  0.54853458  0.55354108
	0.13244401  0.55221637  0.55301828
	0.13117249  0.55589872  0.55245948
	0.1299327   0.55958162  0.55186354
	0.12872938  0.56326503  0.55122927
	0.12756771  0.56694891  0.55055551
	0.12645338  0.57063316  0.5498411 
	0.12539383  0.57431754  0.54908564
	0.12439474  0.57800205  0.5482874 
	0.12346281  0.58168661  0.54744498
	0.12260562  0.58537105  0.54655722
	0.12183122  0.58905521  0.54562298
	0.12114807  0.59273889  0.54464114
	0.12056501  0.59642187  0.54361058
	0.12009154  0.60010387  0.54253043
	0.11973756  0.60378459  0.54139999
	0.11951163  0.60746388  0.54021751
	0.11942341  0.61114146  0.53898192
	0.11948255  0.61481702  0.53769219
	0.11969858  0.61849025  0.53634733
	0.12008079  0.62216081  0.53494633
	0.12063824  0.62582833  0.53348834
	0.12137972  0.62949242  0.53197275
	0.12231244  0.63315277  0.53039808
	0.12344358  0.63680899  0.52876343
	0.12477953  0.64046069  0.52706792
	0.12632581  0.64410744  0.52531069
	0.12808703  0.64774881  0.52349092
	0.13006688  0.65138436  0.52160791
	0.13226797  0.65501363  0.51966086
	0.13469183  0.65863619  0.5176488 
	0.13733921  0.66225157  0.51557101
	0.14020991  0.66585927  0.5134268 
	0.14330291  0.66945881  0.51121549
	0.1466164   0.67304968  0.50893644
	0.15014782  0.67663139  0.5065889 
	0.15389405  0.68020343  0.50417217
	0.15785146  0.68376525  0.50168574
	0.16201598  0.68731632  0.49912906
	0.1663832   0.69085611  0.49650163
	0.1709484   0.69438405  0.49380294
	0.17570671  0.6978996   0.49103252
	0.18065314  0.70140222  0.48818938
	0.18578266  0.70489133  0.48527326
	0.19109018  0.70836635  0.48228395
	0.19657063  0.71182668  0.47922108
	0.20221902  0.71527175  0.47608431
	0.20803045  0.71870095  0.4728733 
	0.21400015  0.72211371  0.46958774
	0.22012381  0.72550945  0.46622638
	0.2263969   0.72888753  0.46278934
	0.23281498  0.73224735  0.45927675
	0.2393739   0.73558828  0.45568838
	0.24606968  0.73890972  0.45202405
	0.25289851  0.74221104  0.44828355
	0.25985676  0.74549162  0.44446673
	0.26694127  0.74875084  0.44057284
	0.27414922  0.75198807  0.4366009 
	0.28147681  0.75520266  0.43255207
	0.28892102  0.75839399  0.42842626
	0.29647899  0.76156142  0.42422341
	0.30414796  0.76470433  0.41994346
	0.31192534  0.76782207  0.41558638
	0.3198086   0.77091403  0.41115215
	0.3277958   0.77397953  0.40664011
	0.33588539  0.7770179   0.40204917
	0.34407411  0.78002855  0.39738103
	0.35235985  0.78301086  0.39263579
	0.36074053  0.78596419  0.38781353
	0.3692142   0.78888793  0.38291438
	0.37777892  0.79178146  0.3779385 
	0.38643282  0.79464415  0.37288606
	0.39517408  0.79747541  0.36775726
	0.40400101  0.80027461  0.36255223
	0.4129135   0.80304099  0.35726893
	0.42190813  0.80577412  0.35191009
	0.43098317  0.80847343  0.34647607
	0.44013691  0.81113836  0.3409673 
	0.44936763  0.81376835  0.33538426
	0.45867362  0.81636288  0.32972749
	0.46805314  0.81892143  0.32399761
	0.47750446  0.82144351  0.31819529
	0.4870258   0.82392862  0.31232133
	0.49661536  0.82637633  0.30637661
	0.5062713   0.82878621  0.30036211
	0.51599182  0.83115784  0.29427888
	0.52577622  0.83349064  0.2881265 
	0.5356211   0.83578452  0.28190832
	0.5455244   0.83803918  0.27562602
	0.55548397  0.84025437  0.26928147
	0.5654976   0.8424299   0.26287683
	0.57556297  0.84456561  0.25641457
	0.58567772  0.84666139  0.24989748
	0.59583934  0.84871722  0.24332878
	0.60604528  0.8507331   0.23671214
	0.61629283  0.85270912  0.23005179
	0.62657923  0.85464543  0.22335258
	0.63690157  0.85654226  0.21662012
	0.64725685  0.85839991  0.20986086
	0.65764197  0.86021878  0.20308229
	0.66805369  0.86199932  0.19629307
	0.67848868  0.86374211  0.18950326
	0.68894351  0.86544779  0.18272455
	0.69941463  0.86711711  0.17597055
	0.70989842  0.86875092  0.16925712
	0.72039115  0.87035015  0.16260273
	0.73088902  0.87191584  0.15602894
	0.74138803  0.87344918  0.14956101
	0.75188414  0.87495143  0.14322828
	0.76237342  0.87642392  0.13706449
	0.77285183  0.87786808  0.13110864
	0.78331535  0.87928545  0.12540538
	0.79375994  0.88067763  0.12000532
	0.80418159  0.88204632  0.11496505
	0.81457634  0.88339329  0.11034678
	0.82494028  0.88472036  0.10621724
	0.83526959  0.88602943  0.1026459 
	0.84556056  0.88732243  0.09970219
	0.8558096   0.88860134  0.09745186
	0.86601325  0.88986815  0.09595277
	0.87616824  0.89112487  0.09525046
	0.88627146  0.89237353  0.09537439
	0.89632002  0.89361614  0.09633538
	0.90631121  0.89485467  0.09812496
	0.91624212  0.89609127  0.1007168 
	0.92610579  0.89732977  0.10407067
	0.93590444  0.8985704   0.10813094
	0.94563626  0.899815    0.11283773
	0.95529972  0.90106534  0.11812832
	0.96489353  0.90232311  0.12394051
	0.97441665  0.90358991  0.13021494
	0.98386829  0.90486726  0.13689671
	0.99324789  0.90615657  0.1439362];

	if (m ~= 256)
		s = size(c,1);
		c = interp1(1:size(c,1), c, linspace(1,s,m), 'linear');
	end

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

function varargout = mosaicer(varargin)
% Helper window to paste SRTM grid tiles or Wem image tiles from Bing and others

%	Copyright (c) 2004-2013 by J. Luis
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
	mosaicer_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	handles.hRectangle = varargin{1};
	handles.rectOrigSize = [get(handles.hRectangle,'XData') get(handles.hRectangle,'YData')];
	handles.hMirFig = get(get(handles.hRectangle,'Parent'), 'Parent');
	handMir = guidata(handles.hMirFig);
	handles.last_dir = handMir.last_dir;
	handles.path_data = handMir.path_data;
	handles.path_tmp = handMir.path_tmp;

	handles.hPatches30 = [];
	handles.hPatches5 = [];
	handles.hPatches = [];
	handles.hPatchesGMTED = [];
	handles.hPatchesACE = [];
	handles.hPatchImgs = [];
	handles.url_manual = '';

	% Choose what type of patches to plot based on rectangle's area (GRID oriented decision)
	x = diff(get(handles.hRectangle,'XData'));		y = diff(get(handles.hRectangle,'YData'));
	area = max(abs(x)) * max(abs(y));
	if (area < 30)
		handles = draw_srtm_mesh(handles);		% True SRTM tiles
	elseif (area < 150)
		handles = draw_srtm5_mesh(handles);		% The CGIAR 5º blends
		set(handles.radio_srtm5, 'Val', 1),			set(handles.radio_srtm, 'Val', 0)
	else
		handles = draw_srtm30_mesh(handles);		% Sandwell's 30' big ones
		set(handles.radio_srtm30, 'Val', 1),		set(handles.radio_srtm, 'Val', 0)
	end

	% --------------------- Read the directory and cache dir list from mirone_pref ----------------------
	directory_list = [];
	load([handMir.path_data 'mirone_pref.mat']);
	j = false(1,numel(directory_list));						% vector for eventual cleaning non-existing dirs

	if iscell(directory_list)								% When exists a dir list in mirone_pref
		for i = 1:numel(directory_list)
			if ~exist(directory_list{i},'dir'),   j(i) = true;   end
		end
		directory_list(j) = [];								% clean eventual non-existing directories
	end
	if ~isempty(directory_list)								% If there is one left
		if (~strcmp(directory_list{1}, handles.last_dir))
			directory_list = [{handles.last_dir}; directory_list];
		end
	else
		directory_list = {handles.last_dir};
	end
	set(handles.popup_grd_dir,'String',directory_list)
	handles.last_directories = directory_list;
	handles.files_dir = handles.last_directories{1};

	handles = popup_grd_dir_CB(handles.popup_grd_dir, handles, directory_list{1});

	% -------- Try if we already have a cache directory store in prefs
	try			cacheDirs = cacheTilesDir;
	catch,		cacheDirs = [];
	end

	if ( ~isempty(cacheDirs) )			% We have cache dir(s), but make sure it exists
		ind = true(1,numel(cacheDirs));
		for (k = 1:numel(ind))					% dir(s) are stored in a cell array
			if ( exist(cacheDirs{k}, 'dir') == 7)
				ind(k) = false;
			end
		end
		cacheDirs(ind) = [];			% Remove eventual non existent dirs
		if ( isempty(cacheDirs) ),		cacheDirs = {''};		val = 1;	% It can't be an empty var
		else							cacheDirs = [{''}; cacheDirs];		val = 2;
		end
		set(handles.popup_cache_dir,'String',cacheDirs, 'Val', val)
	end
	% ------------------------------------------------------------------------------------------

	% ---------- Read the tilesServers.txt file and swalow its contents ------------------------
	% NOTE: an important part of this section is futuristic code
	handles.serversImageOut = 'http://a0.ortho.tiles.virtualearth.net/tiles/a';		% DEFAULT TO VE
	handles.serversRoadOut = 'http://r0.ortho.tiles.virtualearth.net/tiles/r';
	handles.serversHybridOut = 'http://h0.ortho.tiles.virtualearth.net/tiles/h';
	handles.serversQuadkeyOut = {'0' '1'; '2' '3'};
	handles.serversOrder = ones(3,1);
	
	fid = fopen([handles.path_data 'tilesServers.txt']);
	if (fid > 0)				% Otherwise, non-existent file, revert to defaults
		c = (fread(fid,'*char'))';      fclose(fid);
		servers = strread(c,'%s','delimiter','\n');   clear c fid;
		% Remove comment and also eventual empty lines
		m = numel(servers);		c = false(m,1);
		for (k = 1:m)
			if ( isempty(servers{k}) || servers{k}(1) == '#' ),		c(k) = true;	end
		end
		servers(c) = [];		n_servers = numel(servers);

		handles.servers_image  = cell(1,n_servers);		handles.servers_road = cell(1,n_servers);
		handles.servers_hybrid = cell(1,n_servers);		handles.servers_quadkey = cell(1,n_servers);
		for (k = 1:n_servers)					% Loop over number of servers (not that many, we hope)
			[tok,rem] = strtok(servers{k});
			handles.servers_image{k} = tok;		% Images server
			if (~isempty(rem))					% Maps server
				[tok,rem] = strtok(rem);		handles.servers_road{k} = tok;
			end
			if (~isempty(rem))					% Hybrids server
				[tok,rem] = strtok(rem);		handles.servers_hybrid{k} = tok;
			end
			if (~isempty(rem))					% Quadtree keeword
				tok = strtok(rem);				handles.servers_quadkey{k} = tok;
			end
		end
	end
	% -------------------------------------------------------------------------------------------

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, [handles.text_TitFrame1 handles.text_TitFrame2])
	%------------- END Pro look (3D) ------------------------------

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hMirFig,'dependentFigs',plugedWin);

	set(hObject,'Visible','on');
	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end

% -------------------------------------------------------------------------
function hand = popup_grd_dir_CB(hObject, handles, opt)
% Do whatever in the directory containing the grids of interest
	if (nargin == 2)    opt = [];   end
	if isempty(opt)
		val = get(hObject,'Value');     str = get(hObject, 'String');
		% Put the selected field on top of the String list. This is necessary because the "OK" button will
		tmp = str(val);				str(val) = [];
		new_str = [tmp; str];		set(hObject,'String',new_str); 
		set(hObject,'Value',1)
		if iscell(tmp),				new_dir = tmp{1};
		else						new_dir = tmp;
		end
	else
		new_dir = opt;
	end
	handles.files_dir = new_dir;

	exts = {'hgt' 'zip' 'gz'};		prefix = '';
	if (get(handles.radio_srtm5, 'Val'))		% Here we will only search for zip files
		exts = 'zip';	prefix = 'srtm_';
	elseif (get(handles.radio_srtm30, 'Val'))
		exts{1} = 'srtm';
	end
	if (get(handles.radio_srtm, 'Val'))
		[handles.srtm_files, handles.srtm_compfiles, handles.srtm_ext, handles.srtm_pato, handles.srtm_pato_comp] = ...
			get_fnames_ext(new_dir, exts, prefix);
	elseif (get(handles.radio_srtm30, 'Val'))
		[handles.srtm30_files, handles.srtm30_compfiles, handles.srtm30_ext, handles.srtm30_pato, handles.srtm30_pato_comp] = ...
			get_fnames_ext(new_dir, exts, prefix);
	else
		[handles.srtm5_files, handles.srtm5_compfiles, handles.srtm5_ext, handles.srtm5_pato, handles.srtm5_pato_comp] = ...
			get_fnames_ext(new_dir, exts, prefix);
	end

	if (nargout),	hand = handles;
	else			guidata(handles.figure1,handles)
	end

% -------------------------------------------------------------------------
function push_grd_dir_CB(hObject, handles)
% Change dir to where files of interest are stored. Call popup_grd_dir to finish work
	if (strcmp(computer, 'PCWIN'))
		work_dir = uigetfolder_win32('Select a directory', cd);
	else            % This guy doesn't let to be compiled
		work_dir = uigetdir(cd, 'Select a directory');
	end
	if (isempty(work_dir)),		return,		end
	handles.last_directories = [cellstr(work_dir); handles.last_directories];
	set(handles.popup_grd_dir,'Str',handles.last_directories)
	guidata(hObject, handles);
	popup_grd_dir_CB(hObject, handles, work_dir)

% -------------------------------------------------------------------------
function radio_srtm_CB(hObject, handles)
	radio_types(hObject, handles)
	if (~get(handles.check_web,'Val') && ~isempty(get(handles.edit_url,'Str')))
		set(handles.edit_url,'Str','')		% Hapened if SRTM5 + Web was attemped
	end

% -------------------------------------------------------------------------
function radio_srtm5_CB(hObject, handles)
	radio_types(hObject, handles)

% -------------------------------------------------------------------------
function radio_srtm30_CB(hObject, handles)
	radio_types(hObject, handles)
	if (~get(handles.check_web,'Val') && ~isempty(get(handles.edit_url,'Str')))
		set(handles.edit_url,'Str','')		% Hapened if SRTM5 + Web was attemped
	end

% -------------------------------------------------------------------------
function radio_gmted075_CB(hObject, handles)
	radio_types(hObject, handles)

% -------------------------------------------------------------------------
function radio_gmted150_CB(hObject, handles)
	radio_types(hObject, handles)

% -------------------------------------------------------------------------
function radio_gmted300_CB(hObject, handles)
	radio_types(hObject, handles)

% -------------------------------------------------------------------------
function radio_types(hObject, handles)
% Centralized function that sets/unsets different properties of the 'radios' family
	if (~get(hObject,'Val')),		set(hObject,'Val',1),	return,		end
	nome = get(hObject,'Str');
	if (nome(1) == 'S' && nome(end-1) == '5'),	nome(end) = [];		end		% Don't take the risk to try to recognize the º symbol
	names = {'SRTM' 'radio_srtm' 'hPatches';
		'SRTM30' 'radio_srtm30' 'hPatches30';
		'SRTM5' 'radio_srtm5' 'hPatches5';
		'GMTED075' 'radio_gmted075' 'hPatchesGMTED';
		'GMTED150' 'radio_gmted150' 'hPatchesGMTED';
		'GMTED300' 'radio_gmted300' 'hPatchesGMTED'};
% 		'ACE' 'radio_ace' 'hPatchesACE'};
	ind = strcmp(nome, names(:,1));
	this_radio = names(ind,:);					% Isolate data referring to current name
	names(ind,:) = [];							% Remove current name from the pack list
	for (k = 1:size(names,1))
		set(handles.(names{k,2}), 'Val', 0)		% Set other radios to 0
		if (~isempty(handles.(names{k,3})))
			set(handles.(names{k,3}), 'Vis', 'off')	% Hide already created patches of different type
		end
	end

	if ( ~isempty(handles.(this_radio{3})) )	% Patches already exist. Just set them visible
		set(handles.(this_radio{3}), 'Vis', 'on')
	else
		switch nome
			case 'SRTM',	handles = draw_srtm_mesh(handles);
			case 'SRTM30',	handles = draw_srtm30_mesh(handles);
			case 'SRTM5',	handles = draw_srtm5_mesh(handles);
			case {'GMTED075' 'GMTED150' 'GMTED300'}
							handles = draw_gmted_mesh(handles);
			case 'ACE',		handles = draw_ace_mesh(handles);
		end
		if ( any(strcmp({'SRTM' 'SRTM30' 'SRTM5'}, nome)) )
			handles = popup_grd_dir_CB(handles.popup_grd_dir, handles);		% Find if def dir has files of interest
		end
		guidata(handles.figure1,handles)

		if ( any(strcmp({'GMTED075' 'GMTED150' 'GMTED300' 'ACE'}, nome)) )	% Web only cases. Different treatment
			set(handles.edit_url, 'Str', 'Web only (Via internaly generated VRTs)', 'Enable', 'off')
			set(handles.check_web, 'Val', 1, 'Enable', 'off')
			return
		else
			set([handles.check_web handles.edit_url], 'Enable', 'on')		% Other cases are free to be changed
		end
	end
	if (get(handles.check_web,'val')),		check_web_CB(handles.check_web, handles),	end

% -------------------------------------------------------------------------
function check_web_CB(hObject, handles)
% Set pre-set web addresses of known data sources
	url = '';
	if (get(hObject,'Val'))
		if (get(handles.radio_srtm5,'Val'))		% This one is unbearably slow
			set(handles.edit_url, 'Str', 'This one is unbearably slow')
			set(hObject,'Val',0),	return
		end
		set([handles.popup_grd_dir handles.push_grd_dir],'Enable', 'off')
		if (get(handles.radio_srtm,'Val'))
			url = 'http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Eurasia/';
		elseif (get(handles.radio_srtm30,'Val'))
			url = 'ftp://topex.ucsd.edu/pub/srtm30_plus/srtm30/erm/';
		end
	else
		set([handles.popup_grd_dir handles.push_grd_dir],'Enable', 'on')
	end
	set(handles.edit_url, 'Str', url)

% -------------------------------------------------------------------------
function edit_url_CB(hObject, handles)
% Save entered address
	handles.url_manual = get(hObject,'Str');
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function handles = draw_srtm_mesh(handles)
% Draw 1 degree squares corresponding to the SRTM1|3 tiles
	hAx = get(handles.hRectangle,'Parent');
	x = get(handles.hRectangle,'XData');		y = get(handles.hRectangle,'YData');
	x_min = floor(min(x));	x_max = ceil(max(x));
	y_min = floor(min(y));	y_max = ceil(max(y));
	x = x_min:x_max;		y = y_min:y_max;
	n = numel(x);			m = numel(y);
	hp = zeros(m-1,n-1);
	for (i = 1:m-1)
		yp = [y(i) y(i+1) y(i+1) y(i) y(i)];
		c1 = 'N';			c2 = 'E';			% Default guesses
		if (y(i) < 0),		c1 = 'S';	end
		for (j = 1:n-1)		% col
			xp = [x(j) x(j) x(j+1) x(j+1) x(j)];
			mesh_idx = sprintf('%dx%d', i,j);
			if (x(j) < 0),	c2 = 'W';	end
			tag = sprintf('%s%.2d%s%.3d.hgt', c1, abs(y(i)), c2, abs(x(j)));
			hp(i,j) = patch('Xdata',xp, 'YData',yp, 'Parent',hAx, 'FaceColor','y', 'FaceAlpha',0.5, ...
				'Tag',tag, 'UserData',0, 'ButtonDownFcn',{@bdn_Tile, handles.figure1});
			setappdata(hp(i,j),'MeshIndex',mesh_idx)
		end
	end
	uistack_j(handles.hRectangle,'top')
	handles.hPatches = hp;

% -----------------------------------------------------------------------------------------
function handles = draw_srtm30_mesh(handles)
% Draw patches corresponding to the SRTM30 tiles
	hAx = get(handles.hRectangle,'Parent');
	[x, y, xP, yP] = mosaic_grid(handles.hRectangle, 40, 50);
	n = numel(xP);					m = numel(yP);
	hp = zeros(m-1,n-1);

	for (i = 1:m-1)
		yp = [yP(i) yP(i+1) yP(i+1) yP(i) yP(i)];
		c2 = 'n';			c1 = 'w';			% Default guesses
		if (yp(2) < 0),		c2 = 's';	end
		for (j = 1:n-1)		% col
			xp = [xP(j) xP(j) xP(j+1) xP(j+1) xP(j)];
			mesh_idx = sprintf('%dx%d', i,j);
			if (xp(1) > 0),	c1 = 'e';	end
 			tag = sprintf('%s%.3d%s%.2d.Bathymetry.srtm', c1, abs(xp(1)), c2, abs(yp(2)));
			hp(i,j) = patch('Xdata',xp, 'YData',yp, 'Parent',hAx, 'FaceColor','y', 'FaceAlpha',0.5, ...
				'Tag',tag, 'UserData',0, 'ButtonDownFcn',{@bdn_Tile, handles.figure1});
			setappdata(hp(i,j),'MeshIndex',mesh_idx)
		end
	end
% 	% Southern tile row has a different width
% 	yp = [-90 -60 -60 -90 -90];
% 	for (j = 1:6)
% 		xp = [(j-1)*60 (j-1)*60 j*60 j*60 (j-1)*60] - 180;
% 		c1 = 'w';
% 		if (xp(1) > 0),		c1 = 'e';	end
% 		tag = sprintf('%s%.3d.s60.Bathymetry.srtm', c1, abs(xp(1)));
% 		patch(xp,yp,0,'FaceColor','none','Tag',tag,'UserData',0,'ButtonDownFcn',@bdn_srtm30Tile);
% 	end
	uistack_j(handles.hRectangle,'top')
	handles.hPatches30 = hp;

% -----------------------------------------------------------------------------------------
function handles = draw_srtm5_mesh(handles)
% Draw 5 degree squares corresponding to the SRTM5 tiles
	hAx = get(handles.hRectangle,'Parent');
	x = get(handles.hRectangle,'XData');	y = get(handles.hRectangle,'YData');
	jx = fix(abs(x+180)/5) + 1;		jx = [min(jx) max(jx)];
	x_min = (jx(1)-1) * 5 - 180;	x_max = jx(2) * 5 - 180;
	iy = fix(abs(60-y)/5) + 1;		iy = [min(iy) max(iy)];
	y_max = 60 - (iy(1)-1) * 5;		y_min = 60 - iy(2) * 5;
	xP = x_min:5:x_max;				yP = y_min:5:y_max;
	n = numel(xP);					m = numel(yP);
	hp = zeros(m-1,n-1);
	for (i = 1:m-1)
		yp = [yP(i) yP(i+1) yP(i+1) yP(i) yP(i)];
		for (j = 1:n-1)		% col
			xp = [xP(j) xP(j) xP(j+1) xP(j+1) xP(j)];
			mesh_idx = sprintf('%dx%d', i,j);
			tag = sprintf('srtm_%.2d_%.2d.zip', jx(1)+(j-1), iy(2)-(i-1));
			hp(i,j) = patch('Xdata',xp, 'YData',yp, 'Parent',hAx, 'FaceColor','y', 'FaceAlpha',0.5, ...
				'Tag',tag, 'UserData',0, 'ButtonDownFcn',{@bdn_Tile, handles.figure1});
			setappdata(hp(i,j),'MeshIndex',mesh_idx)
		end
	end
	uistack_j(handles.hRectangle,'top')
	handles.hPatches5 = hp;

% -----------------------------------------------------------------------------------------
function handles = draw_gmted_mesh(handles)
% Draw a single patch of the size of displayed map. We ignore here that GMTED is made of
% 20x30 degrees tiles because of the multiple resolution issue. Res is selected via radio button
	hAx = get(handles.hRectangle,'Parent');
	[x, y, xP, yP, x_min, y_max, jx, iy] = mosaic_grid(handles.hRectangle, 30, 20);
	n = numel(xP);					m = numel(yP);
	hp = zeros(m-1,n-1);
	for (i = 1:m-1)
		yp = [yP(i) yP(i+1) yP(i+1) yP(i) yP(i)];
		for (j = 1:n-1)		% col
			xp = [xP(j) xP(j) xP(j+1) xP(j+1) xP(j)];
			mesh_idx = sprintf('%dx%d', i,j);
			tag = sprintf('gmted_%.2d_%.2d.zip', jx(1)+(j-1), iy(2)-(i-1));
			hp(i,j) = patch('Xdata',xp, 'YData',yp, 'Parent',hAx, 'FaceColor','y', 'FaceAlpha',0.5, ...
							'Tag',tag, 'UserData',0, 'ButtonDownFcn',{@bdn_Tile, handles.figure1});
			setappdata(hp(i,j),'MeshIndex',mesh_idx)
		end
	end
	uistack_j(handles.hRectangle,'top')
	handles.hPatchesGMTED = hp;

% -----------------------------------------------------------------------------------------
function handles = draw_ace_mesh(handles)
% Draw 15 degree squares corresponding to the ACE tiles
	hAx = get(handles.hRectangle,'Parent');
	[x, y, xP, yP, x_min, y_max, jx, iy] = mosaic_grid(handles.hRectangle, 15, 15);
	n = numel(xP);			m = numel(yP);
	hp = zeros(m-1,n-1);
	for (i = 1:m-1)
		yp = [yP(i) yP(i+1) yP(i+1) yP(i) yP(i)];
		for (j = 1:n-1)		% col
			xp = [xP(j) xP(j) xP(j+1) xP(j+1) xP(j)];
			mesh_idx = sprintf('%dx%d', i,j);
			tag = sprintf('srtm_%.2d_%.2d.zip', jx(1)+(j-1), iy(2)-(i-1));
			hp(i,j) = patch('Xdata',xp, 'YData',yp, 'Parent',hAx, 'FaceColor','y', 'FaceAlpha',0.5, ...
							'Tag',tag, 'UserData',0, 'ButtonDownFcn',{@bdn_Tile, handles.figure1});
			setappdata(hp(i,j),'MeshIndex',mesh_idx)
		end
	end
	uistack_j(handles.hRectangle,'top')
	handles.hPatchesACE = hp;

% -----------------------------------------------------------------------------------------
function bdn_Tile(obj, evt, hFig)
	handles = guidata(hFig);
	tag = get(gcbo,'Tag');
	if (get(handles.radio_srtm, 'Val'))
		xx = strcmp(tag, handles.srtm_files);		xz = strcmp(tag, handles.srtm_compfiles);
	elseif (get(handles.radio_srtm30, 'Val'))
		xx = strcmp(tag, handles.srtm30_files);		xz = strcmp(tag, handles.srtm30_compfiles);
	elseif (get(handles.radio_srtm5, 'Val'))
		xx = strcmp(tag, handles.srtm5_files);		xz = strcmp(tag, handles.srtm5_compfiles);
	else
		xx = strcmp(tag, handles.srtm5_files);		xz = strcmp(tag, handles.srtm5_compfiles);
	end

	if (~get(gcbo,'UserData'))			% If not selected    
		set(gcbo,'FaceColor','g','UserData',1)
		if ( (~any(xx) && ~any(xz)) && isempty(get(handles.edit_url,'Str')) )	% If FALSE ==> WEB dl, we are done
			%str = ['The file ' tag ' does not exist in the current directory'];
			set(gcbo,'FaceColor','r')
			pause(1)
			set(gcbo,'FaceColor','y','UserData',0)
		end
	else
		set(gcbo,'FaceColor','y','UserData',0)
	end

% ------------------------------------------------------------------------------------------
function [files, comp_files, comp_ext, patos, patos_comp] = get_fnames_ext(pato, ext, prefix)
% Get the list of all files with extention "EXT" seating in the "PATO" dir
% EXT may be either a char or a cell array. In the first case, only files with extension EXT
% will be returned (that is;  COMP_FILES & COMP_EXT are empty)
% On the second case, extra values of EXT will will be searched as well (that is; files with
% extension *.EXT{1}.EXT{2:numel(EXT)}.
% FILES is a cell array of chars with the names that have extension EXT.
% COMP_FILES is a cell array of chars with the names that had extension EXT{2, or 3, or 4, etc...}.
% NOTE: the last extension is removed. E.G if file was lixo.dat.zip, it will become lixo.dat
% COMP_EXT is a cell array of chars with the extensions corresponding to COMP_FILES.
% An example is the search for files terminating in *.dat or *.dat.zip (EXT = {'dat' 'zip'})

	comp_files = [];		comp_ext = [];		patos_comp = [];	ext2 = [];
	if (nargin == 2),		prefix = '';		end
	if ( (pato(end) ~= '\') || (pato(end) ~= '/') )
		pato(end+1) = '/';
	end
	if (iscell(ext))
		ext1 = ext{1};
		if (numel(ext) > 1),	ext2 = ext(2:end);		end
	else
		ext1 = ext;
	end

	tmp = dir([pato prefix '*.' ext1]);
	files = {tmp(:).name}';
	patos = repmat({pato},numel(files),1);

	if (~isempty(ext2))				% That is, if we have one or more compression types (e.g. 'zip' 'gz')
		for (k = 1:numel(ext2))		% Loop over compression types
			tmp = dir([pato '*.' ext1 '.' ext2{k}]);
			tmp = {tmp(:).name}';
			tmp1 = cell(1, numel(tmp));
			for (m = 1:numel(tmp))	% Loop over compressed files
				[PATH,FNAME,EXT] = fileparts(tmp{m});
				tmp{m} = [PATH,FNAME];
				tmp1{m} = EXT;		% Save File last extension as well
			end
			comp_files = [comp_files; tmp];%#ok
			comp_ext = [comp_ext; tmp1'];%#ok
		end
		patos_comp = repmat({pato},numel(comp_files),1);
	end

% -------------------------------------------------------------------------
function radio_tileImages_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set(handles.slider_zoomFactor, 'Enable', 'on')
	set([handles.radio_srtm handles.radio_srtm30 handles.radio_srtm5], 'Val', 0)
	if (~isempty(handles.hPatches)),		set(handles.hPatches,  'Vis', 'off'),	end
	if (~isempty(handles.hPatches5)),		set(handles.hPatches5, 'Vis', 'off'),	end
	if (~isempty(handles.hPatches30)),		set(handles.hPatches30,'Vis', 'off'),	end
	if (~isempty(handles.hPatchImgs)),		set(handles.hPatchImgs, 'Vis', 'on')
	else									region2tiles(handles)
	end

% -------------------------------------------------------------------------
function slider_zoomFactor_CB(hObject, handles)
% ...
	zoomLevel = round(get(hObject,'Value')) + 1;		% +1 because slider starts at 0
	pixPerKm = 360 / 2^(zoomLevel - 1) / 256 * 111.3;	% Approx
	region2tiles(handles)
	if (pixPerKm < 0.5)		% Report values im meters
		set(hObject,'Tooltip', sprintf('One pixel ~%.1f m', pixPerKm*1000))
	else
		set(hObject,'Tooltip', sprintf('One pixel ~%.3f km', pixPerKm))
	end
	set(handles.text_zoomFactor, 'String', zoomLevel)

% -----------------------------------------------------------------------------------------
function hand = region2tiles(handles,lon,lat,zoomFactor)
% Get tiles BoundingBoxs and plot patches for each tile of the region delimited by LON, LAT
% region2tiles(handles) LON & LAT are fished from rectangle handle and ZOOMFACTOR frm handles
	hAx = get(handles.hRectangle,'Parent');
	if (~isempty(handles.hPatchImgs))
		delete(handles.hPatchImgs),		handles.hPatchImgs = [];
	end

	if (nargin == 1)
		zoomFactor = round(get(handles.slider_zoomFactor,'Value')) + 1;
		lon = get(handles.hRectangle,'XData');		lon = [min(lon) max(lon)]+[1e-6 -1e-6];
		lat = get(handles.hRectangle,'YData');		lat = [min(lat) max(lat)];
		lat(1) = max(lat(1), -85);					lat(2) = min(lat(2), 85);
	end

	url = url2image('tile2url', lon, lat, zoomFactor,'quadonly',1);
	[lims, tiles_bb] = url2image('quadcoord', url);
	[m,n] = size(url);		hp = zeros(m, n);

	k = 1;
	for (j = 1:n)			% col
		xp = [tiles_bb(k,1) tiles_bb(k,1) tiles_bb(k,2) tiles_bb(k,2)];
		for (i = 1:m)
			yp = [tiles_bb(k,3) tiles_bb(k,4) tiles_bb(k,4) tiles_bb(k,3)];
			hp(i,j) = patch('XData',xp,'YData',yp,'Parent',hAx,'FaceColor','y','EdgeColor','k', ...
				'FaceAlpha',0.5, 'UserData',0, 'ButtonDownFcn',@bdn_tile_region, 'HitTest', 'on');
			k = k + 1;
		end
	end

	handles.tiles_bb = tiles_bb;
	handles.hPatchImgs = hp;
	if (nargout),	hand = handles;
	else			guidata(handles.figure1,handles)
	end

% -----------------------------------------------------------------------------------------
function bdn_tile_region(obj,evt)
	stat = get(gcbo,'UserData');
	if (~stat),		set(gcbo,'FaceColor','g','UserData',1)        % If not selected
	else			set(gcbo,'FaceColor','y','UserData',0)
	end

% -------------------------------------------------------------------------
function popup_cache_dir_CB(hObject, handles, opt)
% Choose one dir for caching tile images
% OPT, used by push_cache_dir, is char array

	if (nargin == 2)    opt = [];   end
	if ( ~isempty(opt) )				% Add a new entry to the cache dir list. Otherwise, just normal popup work.
		contents = get(hObject, 'String');
		if (numel(contents) == 1),	rest = [];
		else						rest = contents(2:end);
		end

		cacheTilesDir = [{opt}; rest];			% Also the var that will be saved in 'mirone_pref'
		set(hObject, 'String', [{''}; cacheTilesDir], 'Val', 2)			% Empty, <=> no cache, always on top
		save([handles.path_data 'mirone_pref.mat'],'cacheTilesDir', '-append', '-v6')		% Update the prefs file
	end

% -------------------------------------------------------------------------
function push_cache_dir_CB(hObject, handles)
	if (strcmp(computer, 'PCWIN'))
		cache_dir = uigetfolder_win32('Select a directory', cd);
	else			% This guy doesn't let to be compiled
		cache_dir = uigetdir(cd, 'Select a directory');
	end
	if ~isempty(cache_dir)
		popup_cache_dir_CB(handles.popup_cache_dir, handles, cache_dir)
	end

% -------------------------------------------------------------------------
function toggle_sat_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.toggle_hybrid handles.toggle_map], 'Val', 0)

% -------------------------------------------------------------------------
function toggle_hybrid_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.toggle_sat handles.toggle_map], 'Val', 0)

% -------------------------------------------------------------------------
function toggle_map_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.toggle_sat handles.toggle_hybrid], 'Val', 0)

% -------------------------------------------------------------------------
function push_src_CB(hObject, handles)
% Get the response of the "Tiles servers" window
	order = handles.serversOrder;		% Order of last call
	out = tiles_servers(handles.servers_image, handles.servers_road, handles.servers_hybrid, order);
	% out.order contains the order that the out.whatkind{i} has in the tilesServers file
	% For example [1 2 1] means that Image and Hybrid are to be fetched from the servers specified
	% on the first line of tilesServers.txt whilst Maps server is the second entry of that file
	if ( ~isempty(out) )
		handles.serversImageOut = out.whatkind{1};
		handles.serversRoadOut = out.whatkind{2};
		handles.serversHybridOut = out.whatkind{3};
		handles.serversOrder = out.order;
		guidata(handles.figure1, handles)
	end

% -------------------------------------------------------------------------
function toggle_mesh_CB(hObject, handles)
% Plot the squared mesh appropriate for the selected dataset
	if (~get(hObject,'Val'))			% Make everybody invisible
		set(handles.hPatches,  'Vis', 'off'),		set(handles.hPatches30,'Vis', 'off')
		set(handles.hPatches5, 'Vis', 'off'),		set(handles.hPatchImgs, 'Vis', 'off')
	else
		rectSize = [get(handles.hRectangle,'XData') get(handles.hRectangle,'YData')];
		if (any(handles.rectOrigSize - rectSize))		% Rectangle was edited, must recompute patches
			if (get(handles.radio_srtm,'Val'))
				handles = draw_srtm_mesh(handles);
			elseif (get(handles.radio_srtm30,'Val'))
				handles = draw_srtm30_mesh(handles);
			elseif (get(handles.radio_srtm5,'Val'))
				handles = draw_srtm5_mesh(handles);
			else				% Image tiles
				handles = region2tiles(handles);
			end
		else							% Either make it visible or compute if first time call
			if (get(handles.radio_srtm,'Val'))
				if (~isempty(handles.hPatches)),	set(handles.hPatches, 'Vis', 'on')
				else								handles = draw_srtm_mesh(handles);
				end
			elseif (get(handles.radio_srtm30,'Val'))
				if (~isempty(handles.hPatches30)),	set(handles.hPatches30, 'Vis', 'on')
				else								handles = draw_srtm30_mesh(handles);
				end
			elseif (get(handles.radio_srtm5,'Val'))
				if (~isempty(handles.hPatches5)),	set(handles.hPatches5, 'Vis', 'on')
				else								handles = draw_srtm5_mesh(handles);
				end
			else
				if (~isempty(handles.hPatchImgs)),	set(handles.hPatchImgs, 'Vis', 'on')
				else								handles = region2tiles(handles);
				end
			end				
		end
	end
	guidata(handles.figure1, handles)

% -------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% Pick the right function to do the building work
	if (get(handles.radio_srtm,'Val'))
		n_tiles = mosaic_srtm(handles);
	elseif (get(handles.radio_srtm5,'Val'))
		n_tiles = mosaic_srtm5(handles);
	elseif (get(handles.radio_srtm30,'Val'))
		n_tiles = mosaic_srtm30(handles);
	elseif (get(handles.radio_gmted075,'Val') || get(handles.radio_gmted150,'Val') ||get(handles.radio_gmted300,'Val'))
		n_tiles = mosaic_gmted(handles);
	elseif (get(handles.radio_tileImages,'Val'))
		n_tiles = mosaic_images(handles);
	end
	if (n_tiles == 0)
		warndlg('No tiles were actualy selected, so there is nothing to do.','Warning')
	end

% -------------------------------------------------------------------------
function n_tiles = mosaic_srtm(handles)
% Build a mosaic of SRTM1|3 tiles
	n_tiles = 0;			% Counter to the number of processed tiles
	[fnames,limits] = sort_patches(handles, 3);
	if (isempty(fnames))	return,		end
	m = 1;		n = 1;		% If one tile only
	if iscell(fnames),		[m,n] = size(fnames);	end
	RC = 1201;				% Default to SRTM3 size
	from_web = get(handles.check_web,'Val');

	z_min = 1e100;     z_max = -z_min;
	Z_tot = repmat(single(NaN), m*(RC-1)+1, n*(RC-1)+1);
	aguentabar(0,'title','Reading SRTM files');		k = 1;
	for (i = 1:m)				% Loop over selected tiles (by rows)
		for (j = 1:n)			%           "              (and by columns)
			cur_file = fnames{i,j};
			if (strcmp(cur_file,'void')),	continue,	end        % blank tile

			if (~from_web)
				ii = strcmp(cur_file, handles.srtm_files);
				if (any(ii))		% File is not compressed
					full_name = [handles.srtm_pato{ii} handles.srtm_files{ii}];
				else				% Try with a compressed version
					ii = strcmp(cur_file, handles.srtm_compfiles);
					if any(ii)		% Got a compressed file.
						if (strcmpi(handles.srtm_ext{ii},'.zip'))
							full_name = ['/vsizip/' handles.srtm_pato_comp{ii} handles.srtm_compfiles{ii} handles.srtm_ext{ii}];
						elseif (strcmpi(handles.srtm_ext(ii),'.gz'))
							full_name = ['/vsigzip/' handles.srtm_pato_comp{ii} handles.srtm_compfiles{ii} handles.srtm_ext{ii}];
						end
					end
				end
			else
				aguentabar(-1,'title',['Downloading file: ' cur_file]);
				full_name = ['/vsizip/vsicurl/' get(handles.edit_url,'Str') '/' cur_file '.zip'];
			end

			% -------- Read the file ... or fail graciously --------------------------------------
			if (isempty(full_name))
				warndlg(['Error finding file ' cur_file]),		continue
			end
			try
				Z = gdalread(full_name, '-U', '-s');
			catch
				warndlg(['WARNERROR: GDAL failed to read file: ' full_name],'Error'),	continue
			end
			% ------------------------------------------------------------------------------------

			if (i == 1 && j == 1 && size(Z,1) == 3601)
				RC = 3601;
				Z_tot = repmat(single(NaN), m*RC, n*RC);	% We must resize the whole thing
			elseif (size(Z,1) ~= RC)
				msgbox('Tiles have different sizes. Probably mixing SRTM3 and SRTM1','Error','error'),	return
			end
			Z(Z <= -32768) = NaN;
			[zzz] = grdutils(Z,'-L');		z_min = min(z_min,zzz(1));		z_max = max(z_max,zzz(2));
			i_r = (1+(i-1)*(RC-1)):i*(RC-1)+1;
			i_c = (1+(j-1)*(RC-1)):j*(RC-1)+1;
			Z_tot(i_r, i_c) = Z;
			aguentabar(k/(m*n))
			k = k + 1;
		end
	end
	n_tiles = m * n;

	if (RC == 1201),	inc = 3 / 3600;		% SRTM3c
	else				inc = 1 / 3600;		% SRTM1c
	end
	tmp.head = [limits z_min z_max 0 inc inc];
	tmp.X = linspace(limits(1), limits(2), size(Z_tot,2));
	tmp.Y = linspace(limits(3), limits(4), size(Z_tot,1));
	tmp.name = 'SRTM blend';
	mirone(Z_tot,tmp);

% -------------------------------------------------------------------------
function n_tiles = mosaic_srtm5(handles)
% Build a mosaic of SRTM 5 degs (CGIAR messy) tiles
	n_tiles = 0;			% Counter to the number of processed tiles
	[fnames,limits] = sort_patches(handles, 5);
	if (isempty(fnames))	return,		end
	m = 1;		n = 1;		% If one tile only
	if iscell(fnames),		[m,n] = size(fnames);	end
	RC = 6001;

	z_min = 1e100;     z_max = -z_min;
	Z_tot = repmat(single(NaN), m*(RC-1)+1, n*(RC-1)+1);
	aguentabar(0,'title','Reading CGIAR files');		k = 1;
	for (i = 1:m)				% Loop over selected tiles (by rows)
		for (j = 1:n)			%           "              (and by columns)
			cur_file = fnames{i,j};
			if (strcmp(cur_file,'void')),	continue,	end        % blank tile

			ii = strcmp(cur_file, handles.srtm5_files);		% the 'srtm_files' are actually zip files
			[PATO,FNAME] = fileparts(cur_file);
			if (~any(ii))
				msgbox(['Something screwed on fishing this tile: ' FNAME],'Error','error'),	continue
			end
			full_name  = ['/vsizip/' handles.srtm5_pato{ii} cur_file '/' FNAME '.tif'];
			full_name2 = ['/vsizip/' handles.srtm5_pato{ii} cur_file '/' FNAME '.asc'];

			try		[Z, att] = gdalread(full_name, '-U', '-s');
			catch
				try
					[Z, att] = gdalread(full_name2, '-U', '-s');
				catch
					msgbox([cur_file ' does not contain a valid .asc or .tif grid file'],'Error','error'),	continue
				end
			end
			if (size(Z,1) == 6000)			% Old version of this mess
				msgbox([cur_file ' is from older V4.0 version, which is not suported here (size matters)'],'Error','error')
				return
			end
			Z(Z == att.Band(1).NoDataValue) = NaN;
			z_min = min(z_min,att.GMT_hdr(5));		z_max = max(z_max,att.GMT_hdr(6));
			i_r = (1+(i-1)*(RC-1)):i*(RC-1)+1;
			i_c = (1+(j-1)*(RC-1)):j*(RC-1)+1;
			Z_tot(i_r, i_c) = Z;
			aguentabar(k/(m*n))
			k = k + 1;
		end
	end
	n_tiles = m * n;

	inc =  3/3600;
	tmp.head = [limits z_min z_max 0 inc inc];
	tmp.X = linspace(limits(1), limits(2), size(Z_tot,2));
	tmp.Y = linspace(limits(3), limits(4), size(Z_tot,1));
	tmp.name = 'SRTM_5_blend';
	mirone(Z_tot,tmp);

% -------------------------------------------------------------------------
function n_tiles = mosaic_srtm30(handles)
% Build a mosaic of SRTM30 tiles
	n_tiles = 0;			% Counter to the number of processed tiles
	[fnames,limits] = sort_patches(handles, 30);
	if (isempty(fnames))	return,		end
	m = 1;		n = 1;		% If one tile only
	if iscell(fnames),		[m,n] = size(fnames);	end
	RC = [6000 4800];
 	from_web = get(handles.check_web,'Val');
% 	from_web = false;
	att = '';

	z_min = 1e100;     z_max = -z_min;	%x_min = 1e10;	x_max = -x_min;	y_min = 1e10;	y_max = -y_min;
	if (m * n > 1)
		Z_tot = repmat(single(NaN), m*RC(1), n*RC(2));		% Only if we are really mosaicing
	end
	aguentabar(0,'title','Reading SRTM30 file(s)');		k = 1;
	for (i = 1:m)				% Loop over selected tiles	(by rows)
		for (j = 1:n)			%			"				(and by columns)
			cur_file = fnames{i,j};
			if (strcmp(cur_file,'void')),	continue,	end		% blank tile

			if (~from_web)
				ii = strcmp(cur_file, handles.srtm30_files);
				if (any(ii))		% File is not compressed
					full_name = [handles.srtm30_pato{ii} handles.srtm30_files{ii}];
					name_hdr = write_esri_hdr(full_name,'SRTM30');
				else				% Try with a compressed version ---------- NOT IMPLEMENTED YET --------
					ii = strcmp(cur_file, handles.srtm30_compfiles);
					full_name = '';
					if any(ii)		% Got a compressed file.
						if (strcmpi(handles.srtm30_ext{ii},'.zip'))
							full_name = ['/vsizip/' handles.srtm30_pato_comp{ii} handles.srtm30_compfiles{ii} handles.srtm30_ext{ii}];
						elseif (strcmpi(handles.srtm_ext(ii),'.gz'))
							full_name = ['/vsigzip/' handles.srtm30_pato_comp{ii} handles.srtm30_compfiles{ii} handles.srtm30_ext{ii}];
						end
						name_hdr = write_vrt(full_name,'SRTM30');
					end
				end
			else
				aguentabar(-1,'title',['Downloading file: ' cur_file]);
				full_name = ['/vsicurl/' get(handles.edit_url,'Str') '/' cur_file];
			end

			% -------- Read the file ... or fail graciously --------------------------------------
			if (isempty(full_name))
				warndlg(['Error finding file ' cur_file]),		continue
			end
			try
				Z = gdalread(full_name, '-U', '-s');
			catch
				warndlg(['WARNERROR: GDAL failed to read file: ' full_name],'Error'),	continue
			end
			% ------------------------------------------------------------------------------------

			Z(Z == att.Band(1).NoDataValue) = NaN;
			z_min = min(z_min,att.GMT_hdr(5));		z_max = max(z_max,att.GMT_hdr(6));
			i_r = (1+(i-1)*RC(1)):i*RC(1);
			i_c = (1+(j-1)*RC(2)):j*RC(2);
			if (m * n == 1),	Z_tot = Z;				% One tile only
			else				Z_tot(i_r, i_c) = Z;
			end
			delete(name_hdr);
			aguentabar(k/(m*n)),	k = k + 1;
		end
	end
	n_tiles = m * n;

	if (isempty(att))
		errordlg('Errors & more errors. Nothing was read','Error')
		aguentabar(1),		return
	end
	limits = limits + [att.GMT_hdr(8) -att.GMT_hdr(8) att.GMT_hdr(9) -att.GMT_hdr(9)]/2;	%SRTM30 are pix reg
	tmp.head = [limits z_min z_max 0 att.GMT_hdr(8:9)];
	tmp.X = linspace(limits(1), limits(2), size(Z_tot,2));
	tmp.Y = linspace(limits(3), limits(4), size(Z_tot,1));
	tmp.name = 'SRTM30_blend';
	mirone(Z_tot,tmp);

% -------------------------------------------------------------------------
function n_tiles = mosaic_gmted(handles)
% Get a subregion of a one or a mosaic of GMTED tiles
	n_tiles = 1;			% Used only to NOT trigger a "No tilles selected" warning
	[x, y, xP, yP, x_min, y_max] = mosaic_grid(handles.hRectangle, 30, 20);

	local_pato = handles.path_tmp;		% VRT files will be stored here
	prefix = '/vsicurl/http://igskmncngs506.cr.usgs.gov/gmted/Global_tiles_GMTED/';

	if (get(handles.radio_gmted075,'Val'))
		prefix = [prefix '075darcsec/bln/'];
		inc = 1 / 480;		n_cols = 14400;		n_rows = 9600;		res = '075';
	elseif (get(handles.radio_gmted150,'Val'))
		prefix = [prefix '150darcsec/bln/'];
		inc = 1 / 240;		n_cols = 7200;		n_rows = 4800;		res = '150';
	else
		prefix = [prefix '300darcsec/bln/'];
		inc = 1 / 120;		n_cols = 3600;		n_rows = 2400;		res = '300';
	end

	n = numel(xP);				m = numel(yP);
	names_vrt = cell((m-1)*(n-1),1);
	W = 'W';	S = 'S';
	ii = 0;		k = 1;
	for (i = m:-1:2)			% Row
		if (y(i) >= 0),			S = 'N';	end
		for (j = 1:n-1)			% Col
			if (xP(j) >= 0),	W = 'E';	end
			name_ = sprintf('%.2d%c%.3d%c_20101117_gmted_bln%s.tif', abs(yP(i-1)), S, abs(xP(j)), W, res);
			name = {[local_pato name_]; sprintf('%s/%c%.3d/%s',prefix,W,abs(xP(j)),name_)};
			names_vrt{k,1} = write_vrt(name, [xP(j) yP(i) n_cols n_rows inc inc], [], ...
				'nodata',-32768, 'PixelOffset',2, 'relative',0, 'source','simple');
			names_vrt{k,2} = ii;		% Indices to compute the the '<DstRect xOff yOff later in write_vrt
			names_vrt{k,3} = j - 1;
			k = k + 1;
		end
		ii = ii + 1;
	end

	write_vrt([local_pato 'mater.xxx'], [x_min y_max n_cols n_rows inc inc], names_vrt, 'nodata',-32768, 'PixelOffset',2, 'relative',1);

	% Now do the reading
	opt_R = sprintf('-R%.18g/%.18g/%.18g/%.18g', min(x),max(x),min(y),max(y));
	[Z, X, Y, srsWKT, hand] = read_grid([], [local_pato 'mater.vrt'], 'OVR', opt_R);
	zz = grdutils(Z,'-L');		hand.head(5:6) = zz(:)';
	tmp.head = hand.head;		tmp.X = X;		tmp.Y = Y;
	if (~isempty(srsWKT))		tmp.srsWKT = srsWKT;	end
	tmp.name = 'GMTED_cut';
	tmp.was_int16 = hand.was_int16;
	tmp.Nodata_int16 = hand.Nodata_int16;
	mirone(Z,tmp);

% -------------------------------------------------------------------------
function n_tiles = mosaic_ace(handles)
% Build a mosaic of ACE tiles
% http://tethys.eaprs.cse.dmu.ac.uk/ACE2/links/ACE2_3_SECONDS/ACE2_3S_HEIGHTS/00N000E_3S.ACE2.gz
% http://tethys.eaprs.cse.dmu.ac.uk/ACE2/links/ACE2_9_SECONDS/ACE2_9S_HEIGHTS/00N000E_9S.ACE2.gz
% http://tethys.eaprs.cse.dmu.ac.uk/ACE2/links/ACE2_30_SECONDS/ACE2_30S_BOTH/00N000E_BOTH_30S.ACE2.gz
% http://oceancolor.gsfc.nasa.gov/cgi/l3/A20113612011363.L3m_3D_NSST_4.ico.png?sub=img
% http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20113632011365.L3m_3D_NSST_4.bz2
	n_tiles = 0;			% Counter to the number of processed tiles
	[x, y, xP, yP, x_min, y_max] = mosaic_grid(handles.hRectangle, 15, 15);
	
	prefix = '/vsicurl/http://tethys.eaprs.cse.dmu.ac.uk/ACE2/links/ACE2_30_SECONDS/ACE2_30S_BOTH/';
	local_pato = handles.path_tmp;
	n = numel(xP);				m = numel(yP);
	names_vrt = cell((m-1)*(n-1),1);k = 1;
	ii = 0;			W = 'W';	S = 'S';
	for (i = m:-1:2)			% Row
		if (y(i) >= 0),			S = 'N';	end
		for (j = 1:n-1)			% Col
			if (xP(j) >= 0),	W = 'E';	end
			name_ = sprintf('%.2d%c%.3d%c_BOTH_30S.ACE2', abs(yP(i-1)), S, abs(xP(j)), W);
			name = {[local_pato name_]; [prefix name_ '.gz']; 'simple'};
			names_vrt{k,1} = write_vrt(name, [xP(j) yP(i) 1800 1800 1/120 1/120], [], 'nodata',-12000, 'PixelOffset',2, 'relative',1);
			names_vrt{k,2} = ii;
			names_vrt{k,3} = j - 1;
			k = k + 1;
		end
		ii = ii + 1;
	end

	write_vrt([local_pato 'mater.xxx'], [x_min y_max 1800 1800 1/120 1/120], names_vrt, 'nodata',-12000, 'PixelOffset',2, 'relative',1);
	
	% Now do the reading
	opt_R = sprintf('-R%.18g/%.18g/%.18g/%.18g', min(x),max(x),min(y),max(y));
	[Z, X, Y, srsWKT, hand] = read_grid([], [local_pato 'mater.vrt'], 'OVR', opt_R);
	zz = grdutils(Z,'-L');		hand.head(5:6) = zz(:)';
	tmp.head = hand.head;		tmp.X = X;		tmp.Y = Y;
	if (~isempty(srsWKT))		tmp.srsWKT = srsWKT;	end
	tmp.name = 'ACE_cut';
	tmp.was_int16 = hand.was_int16;
	tmp.Nodata_int16 = hand.Nodata_int16;
	mirone(Z,tmp);

% -----------------------------------------------------------------------------------------
function [x, y, xP, yP, x_min, y_max, jx, iy] = mosaic_grid(hRect, dx, dy)
% ...
	x = get(hRect,'XData');			y = get(hRect,'YData');
	jx = fix(abs(x+180)/dx) + 1;	jx = [min(jx) max(jx)];
	x_min = (jx(1)-1) * dx - 180;	x_max = jx(2) * dx - 180;
	iy = fix(abs(90-y)/dy) + 1;		iy = [min(iy) max(iy)];
	y_max = 90 - (iy(1)-1) * dy;	y_min = 90 - iy(2) * dy;
	xP = x_min:dx:x_max;			yP = y_min:dy:y_max;

% -----------------------------------------------------------------------------------------
function [fnames,limits] = sort_patches(handles, type)
% Sort the tile names (those that were selected) in order that follow a matrix
% with origin at the lower left corner of the enclosing rectangle.
% Also returns the limits of the enclosing rectangle.
% TYPE indicates which kind of srtm files are we dealing with (3, 5 or 30)
	x_min = [];		x_max = [];		y_min = [];		y_max = [];		limits = [];
	if (type == 3 || type == 1)		% Get selected tiles by type
		h = findobj(handles.hPatches,'UserData',1);
		tileW = 1;		tileH = 1;	% Tile size in degrees
	elseif (type == 5)
		h = findobj(handles.hPatches5,'UserData',1);
		tileW = 5;		tileH = 5;
	else
		h = findobj(handles.hPatches30,'UserData',1);
		tileW = 40;		tileH = 50;
	end
	names = get(h,'Tag');							% Get their names
	if (isempty(names)),		fnames = [];		return,	end
	if (~isa(names, 'cell')),	names = {names};	end

	nTiles = numel(names);
	mesh_idx = cell(1,nTiles);
	for (i = 1:nTiles),		mesh_idx{i} = getappdata(h(i),'MeshIndex');		end
	[B,IX] = sort(mesh_idx);
	fnames = names(IX);			% Order tile names according to the sorted mesh_idx

	% Find the map limits of the total collection of tiles
	idx_r = zeros(1,nTiles);		idx_c = zeros(1,nTiles);
	for (i = 1:nTiles)
		% Find the tile coordinates from the file name
		if (iscell(fnames))		[PATH,FNAME] = fileparts(fnames{i});
		else					[PATH,FNAME] = fileparts(fnames);
		end
		if (type == 5)			% CGIAR grid
			lon = (sscanf(FNAME(6:7), '%f') - 1) * 5 - 180;
			lat = 60 - sscanf(FNAME(9:10), '%f') * 5;		% so to get the bottom latitude
		else
			FNAME = upper(FNAME);		% SRTM30 came in lower cases
			x_w = strfind(FNAME,'W');	x_e = strfind(FNAME,'E');
			y_s = strfind(FNAME,'S');	y_n = strfind(FNAME,'N');
			if ~isempty(x_w),			ind_x = x_w(1);		lon_sng = -1;
			elseif  ~isempty(x_e)		ind_x = x_e(1);		lon_sng = 1;
			end

			if ~isempty(y_n),			ind_y = y_n(1);		lat_sng = 1;
			elseif ~isempty(y_s)		ind_y = y_s(1);		lat_sng = -1;
			end
			lon = sscanf(FNAME(ind_x+1:ind_x+3), '%f') * lon_sng;
			if (type <= 3),			lat = sscanf(FNAME(2:ind_x-1), '%f') * lat_sng;
			else					lat = sscanf(FNAME(ind_y+1:ind_y+2), '%f') * lat_sng - tileH;
			end
		end

		if (isempty(x_min))
			x_min = lon;    x_max = lon + tileW;    y_min = lat;    y_max = lat + tileH;
		else
			x_min = min(x_min,lon);		x_max = max(x_max,lon+tileW);
			y_min = min(y_min,lat);		y_max = max(y_max,lat+tileH);
		end
		% Convert the sorted mesh index to row and column vectors
		[t,r] = strtok(B{i},'x');
		idx_r(i) = sscanf(t, '%f');		idx_c(i) = sscanf(r(2:end), '%f');
	end
	limits = [x_min x_max y_min y_max];

	% If we have only one tile (trivial case) there is no need for the following tests
	if (nTiles == 1),	return,		end

	% Build a test mesh index with the final correct order 
	min_r = min(idx_r);		max_r = max(idx_r);
	min_c = min(idx_c);		max_c = max(idx_c);
	k = 1;					n_r = numel(min_r:max_r);		n_c = numel(min_c:max_c);
	nNames = numel(min_r:max_r) * numel(min_c:max_c);
	test_mesh = cell(1, nNames);
	for (i = min_r:max_r)
		for (j = min_c:max_c)
			test_mesh{k} = [num2str(i) 'x' num2str(j)];
			k = k + 1;
		end
	end

	t_names = repmat({'void'},1,nNames);
	BB = B;
	for (i = nTiles+1:nNames),		BB{i} = '0x0';		end

	% Check t_names against fnames to find the matrix correct order
	fail = 0;
	for (i = 1:nNames)
		k = i - fail;
		if (strcmp(BB{k},test_mesh{i}))
			t_names{i} = fnames{k};
		else
			fail = fail + 1;
		end
	end
	fnames = reshape(t_names,n_c,n_r)';

% -----------------------------------------------------------------------------------------
function n_tiles = mosaic_images(handles)
% Test if everything is ok and call mosaicing function with (which outputs to Mirone)

	n_tiles = 0;			% Counter to the number of processed tiles
	ind = get(handles.hPatchImgs, 'UserData');
	if (isempty(ind)),		return,		end			% No patches ploted
	ind = logical(cat(1,ind{:}));
	if (~any(ind)),			return,		end			% No patch selected
	
	% ---------------- Have cache info? -----------------
	val = get(handles.popup_cache_dir,'Value');
	contents = get(handles.popup_cache_dir, 'String');
	str = contents{val};		cacheDir = [];
	if (~isempty(str)),			cacheDir = str;		end

	tiles_bb = handles.tiles_bb(ind,:);
	lon = [min(tiles_bb(:,1)) max(tiles_bb(:,2))] + [1 -1] * 1e-6;		% eps is for not getting neighboring tiles
	lat = [min(tiles_bb(:,3)) max(tiles_bb(:,4))] + [1 -1] * 1e-6;
	zoomLevel = get(handles.slider_zoomFactor, 'Val') + 1;
	if (zoomLevel > 18 && get(handles.toggle_map,'Val'))				% Map have 18 as zoom level limit (hybrid???)
		zoomLevel = 18;
	end
	n_tiles = 1;		% Don't know how many actually processed tiles, but it doesn't really matter.

	[whatkind, source_PN, source_PV] = get_kind(handles);
	url2image('callmir',lon,lat, zoomLevel, 'cache', cacheDir, 'what',whatkind, source_PN, source_PV, 'verbose','yes');

% -----------------------------------------------------------------------------------------
function [whatkind, source_PN, source_PV] = get_kind(handles)
% Get the type of imagery selected by the Satellite Hybrid Map buttons

	source_PN = 'treta';		source_PV = [];		% Dumb value used to default to VE
	if (get(handles.toggle_sat,'Val')),			whatkind = 'aerial';
	elseif (get(handles.toggle_map,'Val'))		whatkind = 'road';
	else										whatkind = 'hybrid';
	end
	if ( whatkind(1) == 'a' && isempty(strfind(handles.serversImageOut, 'virtualearth')) )
		server = handles.serversImageOut;		quadkey = handles.servers_quadkey(handles.serversOrder(1));
		if (~strncmp(server, 'http', 4)),		server = ['http://' server];	end
		source_PN = 'source';					source_PV = {server; quadkey{1}};
	elseif ( whatkind(1) == 'r' && isempty(strfind(handles.serversRoadOut, 'virtualearth')) )
		server = handles.serversRoadOut;		quadkey = handles.servers_quadkey(handles.serversOrder(2));
		if (~strncmp(server, 'http', 4)),		server = ['http://' server];	end
		source_PN = 'source';					source_PV = {server; quadkey{1}};
	elseif ( whatkind(1) == 'h' && isempty(strfind(handles.serversHybridOut, 'virtualearth')) )
		server = handles.serversHybridOut;		quadkey = handles.servers_quadkey(handles.serversOrder(3));
		if (~strncmp(server, 'http', 4)),		server = ['http://' server];	end
		source_PN = 'source';					source_PV = {server; quadkey{1}};
	end

% -----------------------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
% Executes when user closes figure
	handles = guidata(hObject);
	% These leave in the Mirone figure
	delete(handles.figure1),		delete(handles.hPatches)
	delete(handles.hPatches5),		delete(handles.hPatches30)
	delete(handles.hPatchesGMTED),	delete(handles.hPatchImgs)

% ----------------------------------------------------------- 
function mosaicer_LayoutFcn(h1)

set(h1,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Mosaicer',...
'NumberTitle','off',...
'Position',[520 513 376 287],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[5 36 366 111], 'Style','frame');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[10 70 70 17],...
'String','Zoom Level',...
'Style','text');

uicontrol('Parent',h1, 'Position',[5 153 366 125], 'Style','frame');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[10 226 120 14],...
'String','Directory with files',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@mosaicer_uiCB,...
'Position',[10 206 331 21],...
'String',' ',...
'Style','popupmenu',...
'Tooltip','Directory were the SRTM files are stored',...
'Value',1,...
'Tag','popup_grd_dir');

uicontrol('Parent',h1,...
'Call',@mosaicer_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'Position',[341 205 23 23],...
'String','...',...
'Tag','push_grd_dir');

uicontrol('Parent',h1, 'Position',[11 244 65 21],...
'Call',@mosaicer_uiCB,...
'String','SRTM',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_srtm');

uicontrol('Parent',h1, 'Position',[90 244 75 21],...
'Call',@mosaicer_uiCB,...
'String','SRTM30',...
'Style','radiobutton',...
'Tag','radio_srtm30');

uicontrol('Parent',h1, 'Position',[170 244 75 21],...
'Call',@mosaicer_uiCB,...
'String','SRTM5º',...
'Style','radiobutton',...
'Tag','radio_srtm5');

uicontrol('Parent',h1, 'Position',[240 260 75 20],...
'Call',@mosaicer_uiCB,...
'String','GMTED075',...
'Style','radiobutton',...
'Tooltip','7.5 seconds resolution grids',...
'Tag','radio_gmted075');

uicontrol('Parent',h1, 'Position',[240 244 75 19],...
'Call',@mosaicer_uiCB,...
'String','GMTED150',...
'Style','radiobutton',...
'Tooltip','15 seconds resolution grids',...
'Tag','radio_gmted150');

uicontrol('Parent',h1, 'Position',[240 227 75 20],...
'Call',@mosaicer_uiCB,...
'String','GMTED300',...
'Style','radiobutton',...
'Tooltip','30 seconds resolution grids',...
'Tag','radio_gmted300');

% uicontrol('Parent',h1, 'Position',[325 244 75 21],...
% 'Call',@mosaicer_uiCB,...
% 'String','ACE',...
% 'Style','radiobutton',...
% 'Tag','radio_ace');

uicontrol('Parent',h1, 'Position',[10 180 140 20],...
'Call',@mosaicer_uiCB,...
'String','OR Web download',...
'Style','checkbox',...
'Tag','check_web');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@mosaicer_uiCB,...
'HorizontalAlignment','left',...
'Position',[10 160 351 22],...
'String','',...
'Style','edit',...
'Tooltip','Enter the URL of root dir where to download from',...
'Tag','edit_url');

uicontrol('Parent',h1,...
'Call',@mosaicer_uiCB,...
'Position',[11 106 95 23],...
'String','Tile Images',...
'Style','radiobutton',...
'Tag','radio_tileImages');

uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[116 268 122 18],...
'String','SRTM type files',...
'Style','text',...
'Tag','text_TitFrame1');

uicontrol('Parent',h1, 'Position',[81 70 261 17],...
'BackgroundColor',[1 1 1],...
'Call',@mosaicer_uiCB,...
'Min',0, 'Max',23,...
'SliderStep', [1 1] / 23,...
'Style','slider',...
'Tooltip','Select the desired zooming factor.',...
'Enable','off',...
'Tag','slider_zoomFactor');

uicontrol('Parent',h1,...
'FontSize',10,...
'FontWeight','demi',...
'HorizontalAlignment','left',...
'Position',[344 72 20 15],...
'String','1',...
'Style','text',...
'Tag','text_zoomFactor');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@mosaicer_uiCB,...
'Position',[90 45 251 22],...
'String',' ',...
'Style','popupmenu',...
'Tooltip','Select a cache directory where to search/save tiles files',...
'Value',1,...
'Tag','popup_cache_dir');

uicontrol('Parent',h1,...
'Call',@mosaicer_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'Position',[341 45 23 23],...
'String','...',...
'Tooltip','Select a different directory',...
'Tag','push_cache_dir');

uicontrol('Parent',h1,...
'Position',[10 46 80 17],...
'HorizontalAlignment','left',...
'FontSize',9,...
'String','Cache directory',...
'Style','text');

uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[97 135 180 18],...
'String','OR Tile images (e.g. Bing)',...
'Style','text',...
'Tag','text_TitFrame2');

uicontrol('Parent',h1,...
'Call',@mosaicer_uiCB,...
'Position',[95 104 57 21],...
'String','Satellite',...
'Style','togglebutton',...
'Tag','toggle_sat');

uicontrol('Parent',h1,...
'Call',@mosaicer_uiCB,...
'Position',[153 104 57 23],...
'String','Hybrid',...
'Style','togglebutton',...
'Tag','toggle_hybrid');

uicontrol('Parent',h1,...
'Call',@mosaicer_uiCB,...
'Position',[211 104 53 21],...
'String','Map',...
'Style','togglebutton',...
'Tag','toggle_map');

uicontrol('Parent',h1,...
'Call',@mosaicer_uiCB,...
'FontAngle','italic',...
'Position',[310 104 51 21],...
'String','Source',...
'Tooltip','Select data source',...
'Tag','push_src');

uicontrol('Parent',h1,...
'Call',@mosaicer_uiCB,...
'FontSize',9,...
'Position',[6 6 121 21],...
'String','Tile mesh on/off',...
'Style','togglebutton',...
'Value',1,...
'Tooltip','Toggle tile mesh on/off',...
'Tag','toggle_mesh');

uicontrol('Parent',h1,...
'Call',@mosaicer_uiCB,...
'FontWeight','bold',...
'Position',[280 6 91 21],...
'String','OK',...
'Tag','push_OK');

function mosaicer_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

% ----------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------
function varargout = tiles_servers(varargin)
% ... 
	if (numel(varargin) < 3)
		error('tiles_servers: input must have 3 arguments')
	end

	hObject = figure('Vis','off');
	tiles_servers_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject, 'center')

	handles.servers_image  = varargin{1};
	handles.servers_road   = varargin{2};
	handles.servers_hybrid = varargin{3};
	order = varargin{4};

	% ------- Fill edit and popup menus with just recieved info -----------
	set(handles.edit_aerial,  'String', handles.servers_image{order(1)} )
	set(handles.edit_road,    'String', handles.servers_road{order(2)} )
	set(handles.edit_hybrid,  'String', handles.servers_hybrid{order(3)} )
	set(handles.popup_aerial, 'String', handles.servers_image, 'Val', order(1))
	set(handles.popup_road,   'String', handles.servers_road,  'Val', order(2))
	set(handles.popup_hybrid, 'String', handles.servers_hybrid,'Val', order(3))
	handles.aerial = handles.servers_image{order(1)};					% Initializations
	handles.road   = handles.servers_road{order(2)};
	handles.hybrid = handles.servers_hybrid{order(3)};
	handles.aerial_ind = order(1);		handles.road_ind = order(2);	handles.hybrid_ind = order(3);
	% ---------------------------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Vis','on');

	% UIWAIT makes tiles_servers wait for user response
	uiwait(handles.figure1);
	handles = guidata(hObject);
	varargout{1} = handles.output;
	if ( ishandle(handles.figure1) ),	delete(handles.figure1),	end

% -----------------------------------------------------------------------------------------
function edit_aerial_CB(hObject, handles, opt)
	if (nargin == 3),		handles.aerial = opt;		% Called by popupmenu
	else					handles.aerial = get(hObject, 'String');
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function edit_road_CB(hObject, handles, opt)
	if (nargin == 3),		handles.road = opt;			% Called by popupmenu
	else					handles.road = get(hObject, 'String');
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function edit_hybrid_CB(hObject, handles, opt)
	if (nargin == 3),		handles.hybrid = opt;		% Called by popupmenu
	else					handles.hybrid = get(hObject, 'String');
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function popup_aerial_CB(hObject, handles)
% Get selected server and transmit info to corresponding edit box (which will do the rest) 
	val = get(hObject,'Value');
	contents = get(hObject,'String');		server = contents{val};
	set(handles.edit_aerial, 'String', server)
	handles.aerial_ind = val;
	edit_aerial_CB(hObject, handles, server)

% -----------------------------------------------------------------------------------------
function popup_road_CB(hObject, handles)
% Get selected server and transmit info to corresponding edit box (which will do the rest) 
	val = get(hObject,'Value');
	contents = get(hObject,'String');		server = contents{val};
	set(handles.edit_road, 'String', server)
	handles.road_ind = val;
	edit_road_CB(hObject, handles, server)

% -----------------------------------------------------------------------------------------
function popup_hybrid_CB(hObject, handles)
% Get selected server and transmit info to corresponding edit box (which will do the rest) 
	val = get(hObject,'Value');
	contents = get(hObject,'String');		server = contents{val};
	set(handles.edit_hybrid, 'String', server)
	handles.hybrid_ind = val;
	edit_hybrid_CB(hObject, handles, server)

% -----------------------------------------------------------------------------------------
function push__OK__CB(hObject, handles)
	out.whatkind = {handles.aerial; handles.road; handles.hybrid};		% Selected server for "whatkind"
	out.order = [handles.aerial_ind; handles.road_ind; handles.hybrid_ind];		% Indices order to recover correspondent quadkeey
	handles.output = out;		guidata(hObject,handles)
	guidata(handles.figure1, handles)
	uiresume(handles.figure1);

% -----------------------------------------------------------------------------------------
function push_cancel_CB(hObject, handles)
	handles.output = [];
	guidata(hObject,handles)
	uiresume(handles.figure1);

% -----------------------------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure_servers_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	handles.output = [];		% User gave up, return nothing
	guidata(hObject, handles);	uiresume(handles.figure1);

% -----------------------------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure_servers_KeyPressFcn(hObject, eventdata)
	handles = guidata(hObject);
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles.output = [];		% User said no by hitting escape
		guidata(hObject, handles);	uiresume(handles.figure1);
	end

% --- Creates and returns a handle to the GUI figure. 
function tiles_servers_LayoutFcn(h1)

set(h1, 'Position',[520 541 550 175],...
'CloseRequestFcn',@figure_servers_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure_servers_KeyPressFcn,...
'MenuBar','none',...
'Name','Tiles servers',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 133 251 21],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'HorizontalAlignment','left',...
'Style','edit', 'Tag','edit_aerial');

uicontrol('Parent',h1, 'Position',[10 83 251 21],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'HorizontalAlignment','left',...
'Style','edit', 'Tag','edit_road');

uicontrol('Parent',h1, 'Position',[10 33 251 21],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_hybrid');

uicontrol('Parent',h1, 'Position',[10 55 110 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Hybrid request from:',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 106 110 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Map request from:',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 157 110 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Image request from:',...
'Style','text');

uicontrol('Parent',h1, 'Position',[269 133 271 22],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'Style','popupmenu',...
'Tooltip','Read from ''tileServers.txt'' on ''data'' directory',...
'Value',1,...
'Tag','popup_aerial');

uicontrol('Parent',h1, 'Position',[269 83 271 22],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'Style','popupmenu',...
'Tooltip','Read from ''tileServers.txt'' on ''data'' directory',...
'Value',1,...
'Tag','popup_road');

uicontrol('Parent',h1, 'Position',[269 33 271 22],...
'BackgroundColor',[1 1 1],...
'Call',@tiles_servers_uiCB,...
'Style','popupmenu',...
'Tooltip','Read from ''tileServers.txt'' on ''data'' directory',...
'Value',1,...
'Tag','popup_hybrid');

uicontrol('Parent',h1, 'Position',[320 157 150 16],...
'FontName','Helvetica',...
'String','Alernatives imported from file',...
'Tooltip','Read from ''tileServers.txt'' on ''data'' directory',...
'Style','text');

uicontrol('Parent',h1, 'Position',[381 3 70 23],...
'Call',@tiles_servers_uiCB,...
'FontName','Helvetica',...
'FontWeight','bold',...
'String','OK',...
'Tag','push__OK_');

uicontrol('Parent',h1, 'Position',[471 3 70 23],...
'Call',@tiles_servers_uiCB,...
'FontName','Helvetica',...
'FontWeight','bold',...
'String','Cancel',...
'Tag','push_cancel');

function tiles_servers_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

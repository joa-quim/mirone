function varargout = srtm_tool(varargin) 
% Helper window to build a SRTM mosaic

%	Copyright (c) 2004-2012 by J. Luis
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

	hObject = figure('Vis','off');
	srtm_tool_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	mir_dirs = getappdata(0,'MIRONE_DIRS');
	if (~isempty(mir_dirs))
		d_path = [mir_dirs.home_dir filesep 'data' filesep];
		handles.path_data = [mir_dirs.home_dir filesep 'tmp' filesep];
	else
		d_path = [pwd filesep 'data' filesep];
		handles.path_data = ['tmp' filesep];
	end

	% See what kind of SRTMs we are going to deal with
	if ~isempty(varargin)
		if (strcmp(varargin{1},'srtm30'))
			handles.have_srtm30 = 1;    handles.have_srtm3 = 0;    handles.have_srtm1 = 0;
			handles.RC = [];
			set(handles.push_draw_rectangle,'Visible','off')  % We don't need it
		elseif (strcmp(varargin{1},'srtm1'))
			handles.have_srtm30 = 0;    handles.have_srtm3 = 0;    handles.have_srtm1 = 1;
			handles.RC = 3600;          % Number of rows & cols -1
		else    % Default to SRTM 3c
			handles.have_srtm30 = 0;    handles.have_srtm3 = 1;    handles.have_srtm1 = 0;
			handles.RC = 1200;          % Number of rows & cols -1
		end
	else        % Default to SRTM 3c
		handles.have_srtm30 = 0;    handles.have_srtm3 = 1;    handles.have_srtm1 = 0;
		handles.RC = 1200;
	end

	handles.figHandle = hObject;
	handles.patchHandles = [];

	handles.w_map = flipdim(imread([d_path 'etopo2.jpg']),1);
	image([-180 180],[-90 90],handles.w_map);   axis xy;

	load([d_path 'm_coasts.mat']);
	line(ncst(:,1),ncst(:,2),'Linewidth',1,'Color','w','Tag','PoliticalBoundaries');

	% Read the directory list from mirone_pref
	directory_list = [];
	load([d_path 'mirone_pref.mat']);
	j = false(1,length(directory_list));				% vector for eventual cleaning non-existing dirs

	if iscell(directory_list)							% When exists a dir list in mirone_pref
		for i = 1:length(directory_list)
			if ~exist(directory_list{i},'dir'),   j(i) = 1;   end
		end
		directory_list(j) = [];							% clean eventual non-existing directories
		if ~isempty(directory_list)						% If there is one left
			set(handles.popup_dir_list,'String',directory_list)
			handles.last_directories = directory_list;
			handles.files_dir = handles.last_directories{1};
			handles.last_dir = handles.last_directories{1};	% Initialize last_dir to files_dir
		else
			set(handles.popup_dir_list,'String',cd)
			handles.last_directories = cellstr(cd);
			handles.files_dir = cd;
			handles.last_dir = cd;
		end
	else												% mirone_pref had no dir list
		handles.last_directories = cellstr(cd);
		set(handles.popup_dir_list,'String',handles.last_directories)
		handles.files_dir = cd;
		handles.last_dir = cd;
	end

	% Find out what files do we have available here
	if (handles.have_srtm30)
		% Get the list of all SRTM30 files seating in the current work dir
		% A Note here. I added this option latter than the SRTM3 and it uses a different
		% strategy for finding files, if they are compressed and what's the compression.
		[handles.srtm30_files,handles.srtm30_compfiles,handles.srtm30_exts] = ...
			get_fnames_ext(handles.files_dir,{'srtm' 'zip' 'gz'});
		guidata(hObject, handles);
		set_srtm30_mesh(handles);       % Draw the SRTM30 rectangle tiles
	else
		% Get the list of all SRTM files seating in the current work dir
		handles.srtm_files = dir([handles.files_dir filesep '*.hgt']);
		handles.srtm_zipfiles = dir([handles.files_dir filesep '*.hgt.zip']);
		handles.srtm_files = rmfield(handles.srtm_files,{'date','isdir','bytes'});      % We don't need those
		handles.srtm_zipfiles = rmfield(handles.srtm_zipfiles,{'date','isdir','bytes'});
		% Rip the compress type extension from file names (It will be easiear to deal with latter down)
		for i=1:length(handles.srtm_zipfiles)
			[PATH,FNAME] = fileparts(handles.srtm_zipfiles(i).name);
			handles.srtm_zipfiles(i).name = [PATH FNAME];
		end
		% Create a new field to store also the files path (FDS I was not able to dot any other way)
		for i=1:length(handles.srtm_files)
			handles.srtm_files(i).path = [handles.files_dir filesep];
		end
		for i=1:length(handles.srtm_zipfiles)
			handles.srtm_zipfiles(i).path = [handles.files_dir filesep];
		end
	end

	% Get the handles of the warning text and set it to invisible
	handles.warnHandle = findobj(hObject,'Tag','text_warning');
	set(handles.warnHandle,'Visible','off')

	% Choose default command line output for srtm_tool_export
	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end
	set(hObject,'Visible','on');

% -----------------------------------------------------------------------------------------
function toggle_zoom_OnOff_CB(hObject, handles)
	if (get(hObject,'Value')),	zoom_j('on');
	else						zoom_j('off');
	end

% -----------------------------------------------------------------------------------------
function push_draw_rectangle_CB(hObject, handles)
	zoom_j('off');     set(handles.toggle_zoom_OnOff,'Value',0);     % In case it was on
	[p1,p2,h] = rubberbandbox;
	set(h,'Color','red','LineWidth',0.5,'Tag','rectangle')
	cmenuHand = uicontextmenu;
	set(h, 'UIContextMenu', cmenuHand);
	ui_edit_polygon(h)		% Set edition functions
	uimenu(cmenuHand, 'Label', 'Delete rectangle', 'Call', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'SRTM mesh', 'Separator','on', 'Call', {@set_srtm_mesh,h,handles});

% -----------------------------------------------------------------------------------------
function popup_dir_list_CB(hObject, handles, opt)
% OPT is used by push_change_dir (just to save code)
	if (~isempty(handles.patchHandles))
		warndlg('Now is too late to change directory. If you realy want to change dir you have to start again.','Warning')
		set(hObject,'Value',1)
		return
	end
	if (nargin == 2)    opt = [];   end
	if isempty(opt)
		val = get(hObject,'Value');     str = get(hObject, 'String');
		% Put the selected field on top of the String list. This is necessary because the "OK" button will
		tmp = str(val);         str(val) = [];
		new_str = [tmp; str];   set(hObject,'String',new_str); 
		set(hObject,'Value',1)
		if iscell(tmp)			new_dir = tmp{1};
		elseif ischar(tmp)		new_dir = tmp;
		else					return        % ???
		end
	else
		new_dir = opt;
	end
	handles.files_dir = new_dir;

	% If are in the SRTM30 mode, just get the new files for the selected dir and go away.
	if (handles.have_srtm30)
		[handles.srtm30_files,handles.srtm30_compfiles,handles.srtm30_exts] = ...
			get_fnames_ext(handles.files_dir,{'srtm' 'zip' 'gz'});
		h_ones = findobj(handles.figure1,'Type','patch','UserData',1);
		if (~isempty(h_ones))    % De-select this patch because we changed directory
			set(h_ones,'FaceColor','none','UserData',0)
		end
		guidata(handles.figHandle,handles)
		return
	end

	% Get the list of all SRTM_1_or_3 files seating in the new selected dir
	more_srtm_files = dir([new_dir filesep '*.hgt']);
	more_srtm_zipfiles = dir([new_dir filesep '*.hgt.zip']);

	% Rip the .zip from file names (It will be easiear to deal with if they are used)
	for i=1:length(more_srtm_zipfiles)
		[PATH,FNAME] = fileparts(more_srtm_zipfiles(i).name);
		more_srtm_zipfiles(i).name = [PATH FNAME];
	end

	% Update the dir structure
	if ~isempty(more_srtm_files)
		n_old = length(handles.srtm_files);     n_new = length(more_srtm_files);    k = 1;
		for (i = n_old+1:n_old+n_new)
			handles.srtm_files(i).name = [more_srtm_files(k).name];
			handles.srtm_files(i).path = [new_dir filesep];
			k = k + 1;
		end
	end
	if ~isempty(more_srtm_zipfiles)
		n_old = length(handles.srtm_zipfiles);     n_new = length(more_srtm_zipfiles);    k = 1;
		for (i=n_old+1:n_old+n_new)
			handles.srtm_zipfiles(i).name = [more_srtm_zipfiles(k).name];
			handles.srtm_zipfiles(i).path = [new_dir filesep];
			k = k + 1;
		end
	end

	% Clear eventual repeated entries (resulting from playing around with the dir_list)
	[b,m] = unique(cellstr(strvcat(handles.srtm_files.name)));
	i = 1;
	while (i < m(1))
		handles.srtm_files(i).name = [];
		handles.srtm_files(i).path = [];
		i = i + 1;
	end
	[b,m] = unique(cellstr(strvcat(handles.srtm_zipfiles.name)));
	i = 1;
	while (i < m(1))
		handles.srtm_zipfiles(i).name = [];
		handles.srtm_zipfiles(i).path = [];
		i = i + 1;
	end
	guidata(handles.figHandle,handles)

% -----------------------------------------------------------------------------------------
function push_change_dir_CB(hObject, handles)
	if (~isempty(handles.patchHandles))
		warndlg('Now is too late to change directory. If you realy want to change dir you have to start again.','Warning')    
		return
	end
	if (strcmp(computer, 'PCWIN'))
		work_dir = uigetfolder_win32('Select a directory', cd);
	else			% This guy doesn't let to be compiled
		work_dir = uigetdir(cd, 'Select a directory');
	end
	if ~isempty(work_dir)
		handles.last_directories = [cellstr(work_dir); handles.last_directories];
		set(handles.popup_dir_list,'String',handles.last_directories)
		guidata(hObject, handles);
		popup_dir_list_CB(hObject, handles, work_dir)
	end

% -----------------------------------------------------------------------------------------
function push_help_CB(hObject, handles)
	message = {'With this tool you can select several SRTM files and paste them togheter'
		'in a single file. In order to do it first zoom the world map and next hit'
		'the "Draw rectangle" button. Click and drag to draw a rectangle (another'
		'click finish it). Next right click on the rectangle and select "SRTM mesh"'
		' '
		'Selecting (clicking on them) individual squares will turn them redish if'
		'the corresponding file resides in the directory selected by the popup'
		'menu. Note that zip compressed files will also be recognized. A further'
		'click on a selected square will deselect it.'
		' '
		'And that''s all. Hit "OK" and wait.'};
	helpdlg(message,'Help on SRTM tool');

% -----------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)

	if (handles.have_srtm30)
		read_srtm30(handles);    return
	end

	[fnames,limits] = sort_patches(handles);
	if (isempty(fnames))	return,		end
	m = 1;		n = 1;		% If one tile only
	if iscell(fnames),		[m,n] = size(fnames);	end

	z_min = [];     z_max = [];     del_file = 0;     x_inc = 1/handles.RC;     y_inc = 1/handles.RC;
	Z_tot = repmat(single(NaN), m*handles.RC+1, n*handles.RC+1);
	aguentabar(0,'Now wait (reading SRTM files)');     k = 1;
	for (i = 1:m)				% Loop over selected tiles (by rows)
		for (j = 1:n)			%           "              (and by columns)
			if (m*n == 1)   % Trivial case (One tile only)
				cur_file = fnames{1};
				ii = strcmp(cur_file, {handles.srtm_files.name});
				if any(ii)
					full_name = [handles.srtm_files(ii).path cur_file];
				else
					ii = strcmp(cur_file, {handles.srtm_zipfiles.name});
					if any(ii)	% Got a zipped file. Unzip it to the TMP dir
						%str = ['unzip -qq ' handles.srtm_zipfiles(ii).path cur_file '.zip' ' -d ' handles.path_data];
						full_name = [handles.srtm_zipfiles(ii).path cur_file '.zip'];
						if (handles.have_srtm1)		warn = 'warn';
						else						warn = [];
						end
						full_name = decompress(full_name,warn);     % Named with compression ext removed
						if (isempty(full_name))		return,		end	% Error message already issued.
						del_file = 1;
					end
				end
			else
				cur_file = fnames{i,j};
				ii = strcmp(cur_file, {handles.srtm_files.name});
				if any(ii)		% It may be empty when fnames{i,j} == 'void' OR file is compressed
					full_name = [handles.srtm_files(ii).path cur_file];
				else                % Try with a compressed (zipped) version
					ii = strcmp(cur_file, {handles.srtm_zipfiles.name});
					if any(ii)		% Got a zipped file. Unzip it to the TMP dir
						full_name = [handles.srtm_zipfiles(ii).path cur_file '.zip'];
						if (handles.have_srtm1),	warn = 'warn';
						else						warn = [];
						end
						full_name = decompress(full_name,warn);     % Named with compression ext removed
						if (isempty(full_name)),	return,		end % Error message already issued.
						del_file = 1;
					end
				end
			end
			if (strcmp(cur_file,'void'))        % blank tile
				Z = repmat(single(NaN),handles.RC+1,handles.RC+1);
			else
				fid = fopen(full_name,'r','b');
				if (fid == -1)		errordlg([full_name ': not found !!!'],'Error'),	return,	end
				Z = single(rot90(fread(fid,[handles.RC+1 handles.RC+1],'*int16')));
				fclose(fid);
				if (del_file)       % Delete the unziped file
					delete(full_name);		del_file = 0;
				end
				% Test that we are not mixing SRTM 1 & 3 second files
				[nr,nc] = size(Z);
				if (handles.have_srtm1 && (nr ~= 3601 || nc ~= 3601))
					errordlg(['File ' full_name ' is not a SRTM 1 second grid'],'Error');
					return
				elseif  (handles.have_srtm3 && (nr ~= 1201 || nc ~= 1201))
					errordlg(['File ' full_name ' is not a SRTM 3 second grid'],'Error');
					return
				end
				Z(Z <= -32768) = NaN;
			end
			if (isempty(z_min))
				[zzz] = grdutils(Z,'-L');		z_min = zzz(1);		z_max = zzz(2);		clear zzz;
			else
				[zzz] = grdutils(Z,'-L');
				z_min = min(z_min,zzz(1));		z_max = max(z_max,zzz(2));
			end
			i_r = (1+(i-1)*handles.RC):i*handles.RC+1;
			i_c = (1+(j-1)*handles.RC):j*handles.RC+1;
			Z_tot(i_r,i_c) = Z;
			aguentabar(k/m*n)
			k = k + 1;
		end
	end
	aguentabar(1)		% Close it
	tmp.head = [limits z_min z_max 0 x_inc y_inc];
	tmp.X = limits(1):x_inc:limits(2);		tmp.Y = limits(3):y_inc:limits(4);
	tmp.name = 'SRTM blend';
	mirone(Z_tot,tmp);

% -----------------------------------------------------------------------------------------
function [fnames,limits] = sort_patches(handles)
% Sort the tile names (those that were selected) in order that follow a matrix
% with origin at the lower left corner of the enclosing rectangle.
% Also returns the limits of the enclosing rectangle.
	x_min = [];     x_max = [];     y_min = [];     y_max = [];     limits = [];
	h = findobj(handles.figure1,'Type','patch','UserData',1);   % Get selected tiles
	names = get(h,'Tag');                           % Get their names
	if (isempty(names))		fnames = [];	return,	end
	if (~isa(names, 'cell')),	names = {names};	end

	nTiles = numel(names);
	mesh_idx = cell(1,nTiles);
	for (i = 1:nTiles)
		mesh_idx{i} = getappdata(h(i),'MeshIndex');
	end
	[B,IX] = sort(mesh_idx);
	if (numel(IX) > 1)
		fnames = names(IX);			% Order tile names according to the sorted mesh_idx
	else
		fnames = names;				% Only one tile; nothing to sort
	end

	% Find the map limits of the total collection of tiles
	idx_r = zeros(1,nTiles);		idx_c = zeros(1,nTiles);
	for (i = 1:nTiles)
		% Find the tile coordinates from the file name
		if (iscell(fnames))		[PATH,FNAME] = fileparts(fnames{i});
		else					[PATH,FNAME] = fileparts(fnames);
		end
		x_w = strfind(FNAME,'W');	x_e = strfind(FNAME,'E');
		y_s = strfind(FNAME,'S');	y_n = strfind(FNAME,'N');
		if ~isempty(x_w)        ind_x = x_w(1);     lon_sng = -1;
		elseif  ~isempty(x_e)   ind_x = x_e(1);     lon_sng = 1;
		end

		if ~isempty(y_n)		lat_sng = 1;
		elseif ~isempty(y_s)	lat_sng = -1;
		end
		lon = str2double(FNAME(ind_x+1:ind_x+3)) * lon_sng;
		lat = str2double(FNAME(2:ind_x-1)) * lat_sng;

		if (isempty(x_min))
			x_min = lon;    x_max = lon + 1;    y_min = lat;    y_max = lat + 1;
		else
			x_min = min(x_min,lon);		x_max = max(x_max,lon+1);
			y_min = min(y_min,lat);		y_max = max(y_max,lat+1);
		end
		% Convert the sorted mesh index to row and column vectors
		[t,r] = strtok(B{i},'x');
		idx_r(i) = str2double(t);		idx_c(i) = str2double(r(2:end));
	end
	limits = [x_min x_max y_min y_max];

	% If we have only one tile (trivial case) there is no need for the following tests
	if (nTiles == 1),	return,		end

	% Build a test mesh index with the final correct order 
	min_r = min(idx_r);		max_r = max(idx_r);
	min_c = min(idx_c);		max_c = max(idx_c);
	k = 1;					n_r = length(min_r:max_r);		n_c = length(min_c:max_c);
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
	for (i = nTiles+1:nNames)		BB{i} = '0x0';		end

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
function set_srtm_mesh(obj,eventdata,h,handles)
	x = get(h,'XData');		y = get(h,'YData');
	x_min = min(x);			x_max = max(x);
	y_min = min(y);			y_max = max(y);
	x_min = floor(x_min);	x_max = ceil(x_max);
	y_min = floor(y_min);	y_max = ceil(y_max);
	x = x_min:x_max;		y = y_min:y_max;
	n = length(x);			m = length(y);
	hp = zeros(m,n);
	for (i = 1:m-1)
		for (j = 1:n-1)		% col
			xp = [x(j) x(j) x(j+1) x(j+1) x(j)];
			yp = [y(i) y(i+1) y(i+1) y(i) y(i)];
			mesh_idx = [num2str(i) 'x' num2str(j)];
			if (y(i) >= 0)  c1 = 'N';
			else            c1 = 'S';
			end
			if (x(j) >= 0)  c2 = 'E';
			else            c2 = 'W';
			end
			tag = [c1 num2str(abs(y(i)),'%.2d') c2 num2str(abs(x(j)),'%.3d') '.hgt'];
			hp(i,j) = patch(xp,yp,'y','FaceAlpha',0.5,'Tag',tag,'UserData',0,'ButtonDownFcn',{@bdn_srtmTile,handles});
			setappdata(hp(i,j),'MeshIndex',mesh_idx)
		end
	end
	handles.patchHandles = hp;
	guidata(handles.figHandle,handles)

% -----------------------------------------------------------------------------------------
function bdn_srtmTile(obj,eventdata,handles)
	tag = get(gcbo,'Tag');
	handles = guidata(handles.figHandle);       % We may need to have an updated version of handles
	xx = strcmp(tag, {handles.srtm_files.name});
	xz = strcmp(tag, {handles.srtm_zipfiles.name});

	stat = get(gcbo,'UserData');
	if ~stat        % If not selected    
		set(gcbo,'FaceColor','r','UserData',1)
		if ((~any(xx) && ~any(xz)))
			str = ['The file ' tag ' does not exist in the current directory'];
			set(handles.warnHandle,'String',str,'Visible','on')
			pause(1)
			set(handles.warnHandle,'Visible','off')
			set(gcbo,'FaceColor','y','UserData',0)
		end
	else
		set(gcbo,'FaceColor','y','UserData',0)
	end

% -----------------------------------------------------------------------------------------
function set_srtm30_mesh(handles)
% Draw patches corresponding to the SRTM30 tiles

	for (i=1:3)
		yp = [-i*50 -(i-1)*50 -(i-1)*50 -i*50 -i*50] + 90;
		if (yp(2) >= 0)  c2 = 'n';
		else             c2 = 's';
		end
		for (j=1:9)
			xp = [(j-1)*40 (j-1)*40 j*40 j*40 (j-1)*40] - 180;
			if (xp(1) >  0)  c1 = 'e';
			else             c1 = 'w';
			end
			tag = [c1 num2str(abs(xp(1)),'%.3d') c2 num2str(abs(yp(2)),'%.2d') '.Bathymetry.srtm'];
			patch(xp,yp,0,'FaceColor','none','Tag',tag,'UserData',0,'ButtonDownFcn',{@bdn_srtm30Tile,handles});
		end
	end
	% Southern tile row has a different width
	yp = [-90 -60 -60 -90 -90];
	for (j = 1:6)
		xp = [(j-1)*60 (j-1)*60 j*60 j*60 (j-1)*60] - 180;
		if (xp(1) > 0)   c1 = 'e';
		else             c1 = 'w';   end
		tag = [c1 num2str(abs(xp(1)),'%.3d') 's60.Bathymetry.srtm'];
		patch(xp,yp,0,'FaceColor','none','Tag',tag,'UserData',0,'ButtonDownFcn',{@bdn_srtm30Tile,handles});
	end

% -----------------------------------------------------------------------------------------
function bdn_srtm30Tile(obj,eventdata,handles)
	tag = get(gcbo,'Tag');
	handles = guidata(handles.figHandle);		% We may need to have an updated version of handles
	xx = strcmp(tag, handles.srtm30_files);   xz = [];
	if (isempty(xx) || ~any(xx))			% Try to see if we have a gziped version
		xz = strcmp(tag, handles.srtm30_compfiles);
	end
	if (isempty(xz) || ~any(xz))			% Try again to see if we have a ziped version
		xz = strcmp(tag, handles.srtm30_compfiles);
	end

	stat = get(gcbo,'UserData');
	if ~stat		% If not selected
		if ( (isempty(xx) && isempty(xz)) || (~any(xx) && ~any(xz)) )      % File not found
			set(gcbo,'FaceColor','r','FaceAlpha',0.5,'UserData',1)
			str = ['The file ' tag ' does not exist in the current directory'];
			set(handles.warnHandle,'String',str,'Visible','on')
			pause(1)
			set(handles.warnHandle,'Visible','off')        
			set(gcbo,'FaceColor','none','UserData',0)
		else		% Found file
			h_ones = findobj(handles.figure1,'Type','patch','UserData',1);
			h = setxor(gcbo,h_ones);
			if (~isempty(h))		% De-select the other selected patch because we don't do mosaic with SRTM30
				set(h,'FaceColor','none','UserData',0)
			end
			set(gcbo,'FaceColor','g','FaceAlpha',0.7,'UserData',1)
		end
	else
		set(gcbo,'FaceColor','none','UserData',0)
	end

% -----------------------------------------------------------------------------------------
function read_srtm30(handles)
% Build the header string that pretends this is a GTOPO30 file and call gdal to do the job.
% We don't read the file directly here because SRTM30 files are too big for Matlab (a Pro $$$ soft!?!?)
	h = findobj(handles.figure1,'Type','patch','UserData',1);   % Get selected tile
	fname = get(h,'Tag');							% Get it's name
	if (isempty(fname))		return,		end			% No file selected

	name_hdr = write_esri_hdr([handles.files_dir filesep fname],'SRTM30');

	del_file = 0;
	ii = strcmp(fname, handles.srtm30_files);
	if (~isempty(ii) || any(ii))
		full_name = [handles.files_dir filesep fname];
	else                % Compressed file
		ii = strcmp(fname, handles.srtm30_compfiles);
		if (~isempty(ii) || any(ii))	% Got a compressed file. Uncompress it to the current dir
			full_name = [handles.files_dir filesep fname handles.srtm30_exts{ii}];
			full_name = decompress(full_name, 'warn');
			if (isempty(full_name))		return,		end		% Error message already issued
			del_file = 1;
		else
			errordlg('Unknown Error. Should not pass here. Returning.', 'Error');
		end
	end

	aguentabar(0.1,'title',['Reading File ' fname]);   % It might take some time

	[Z,att] =  gdalread(full_name,'-U','-C');	Z = single(Z);
	head = att.GMT_hdr;
	if (~isempty(att.Band.NoDataValue))			Z(Z <= single(att.Band.NoDataValue)) = NaN;		end
	[m,n] = size(Z);
	tmp.head = head;							tmp.name = fname;
	tmp.X = linspace(head(1),head(2),n);		tmp.Y = linspace(head(3),head(4),m);

	% Delete auxiliary files
	delete(name_hdr);
	if (del_file)   delete(full_name);  end     % Original file is compressed
	aguentabar(1,'title','File read')
	mirone(Z,tmp);

% -----------------------------------------------------------------------------------------
function [files,comp_files,comp_ext] = get_fnames_ext(pato, ext)
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

	comp_files = [];		comp_ext = [];		ext2 = [];
	if ( (pato(end) ~= '\') || (pato(end) ~= '/') )
		pato(end+1) = filesep;
	end
	if (iscell(ext))
		ext1 = ext{1};
		if (numel(ext) > 1),	ext2 = ext(2:end);		end
	else
		ext1 = ext;
	end

	tmp = dir([pato '*.' ext1]);
	files = {tmp(:).name}';

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
			comp_files = [comp_files; tmp];
			comp_ext = [comp_ext; tmp1'];
		end
	end

% --- Creates and returns a handle to the GUI figure. 
function srtm_tool_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','SRTM Tool',...
'NumberTitle','off',...
'Position',[520 365 756 435],...
'doublebuffer','on', ...
'RendererMode','auto', ...
'Resize','off',...
'Tag','figure1');

h2 = axes('Parent',h1, 'Units','pixels', 'Position',[30 70 721 361],...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','axes1');

set(get(h2,'xlabel'),'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[0.498613037447989 -0.0650969529085872 1.00005459937205],...
'VerticalAlignment','cap');

set(get(h2,'ylabel'),'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[-0.0395284327323162 0.497229916897507 1.00005459937205],...
'Rotation',90,...
'VerticalAlignment','bottom');

uicontrol('Parent',h1, 'Position',[30 27 76 21],...
'Call',@srtm_tool_uiCB,...
'String','Zoom On/Off',...
'Style','togglebutton',...
'Tag','toggle_zoom_OnOff');

uicontrol('Parent',h1, 'Position',[110 27 90 21],...
'Call',@srtm_tool_uiCB,...
'String','Draw rectangle',...
'Tag','push_draw_rectangle');

uicontrol('Parent',h1, 'Position',[219 28 272 22],...
'BackgroundColor',[1 1 1],...
'Call',@srtm_tool_uiCB,...
'Style','popupmenu',...
'Tooltip','Select the directory where SRTM data resides',...
'Value',1,...
'Tag','popup_dir_list');

uicontrol('Parent',h1, 'Position',[615 27 66 21],...
'Call',@srtm_tool_uiCB,...
'String','Help',...
'Tag','push_help');

uicontrol('Parent',h1, 'Position',[685 27 66 21],...
'Call',@srtm_tool_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1,'FontSize',9, 'Position',[30 4 500 15],...
'FontWeight','demi',...
'HorizontalAlignment','left',...
'String','text',...
'Style','text',...
'Tag','text_warning');

uicontrol('Parent',h1, 'Position',[490 29 18 21],...
'Call',@srtm_tool_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'Tooltip','Select a different directory',...
'Tag','push_change_dir');

function srtm_tool_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

function varargout = cartas_militares(varargin)
% Load GIF files with the Portuguese "Cartas Militares" and display them georeferenced in Mirone
%
% Georeferencing is acomplished in one of 3 ways (by that order):
% 1. A corresponding Ozi Explorer .map file is found on the same directory as the GIF file
% 2. A .gfw or .wrl (with contents like in a .tfw file) file is found on the same directory as the GIF file
% 3. Internal georeferencing based on the known coordinates (in "Coordenadas Militares") of each image
% 
% Thirth method should be prefered (more accurate, I believe)

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

% $Id$

	hObject = figure('Visible','off');
	cartas_militares_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	if (numel(varargin) > 0)
		handMir = varargin{1};
        path_data = handMir.path_data;
	else
	    path_data = [cd filesep 'data' filesep];
	end

    naos = cell(19,1);				% By Columns
    naos{1} = [1:19 22:58];
    naos{2} = [1:19 28:58];
    naos{3} = [1:16 30:58];
    naos{4} = [4:16 35:58];
    naos{5} = [1 43:47 57:58];
    naos{6} = [1 58];
    naos{7} = [1 58];
    naos{8} = 1;
    naos{9} = [56 58];
    naos{10} = 56:58;
    naos{11} = 56:58;
    naos{12} = [1 2 56:58];
    naos{13} = [1:9 29:30 56:58];
    naos{14} = [1:11 16:20 26:30 57:58];
    naos{15} = [1:12 15:22 25:31 36 57:58];
    naos{16} = [1:46 57:58];
    naos{17} = [1:47 57:58];
    naos{18} = [1:49 54:58];
    naos{19} = 1:58;
	
    atlas_file = [path_data 'countries_dp5.bin'];
    paises.ct = country_select(atlas_file,'-Pportugal');
	% Clean up the empty fields in the ct struct (given I could not do it at mex level)
	id = false(numel(paises.ct),1);
	for (k = 1:numel(paises.ct))
        if (isempty(paises.ct(k).Country)),     id(k) = true;    end
	end
	paises.ct(id) = [];

	xMap_min = 72000;		yMap_min = 0;

	set(handles.axes1,'XLim',[0 (19*16000)] + xMap_min + [-1000 1000])
	set(handles.axes1,'YLim',[0 (58*10000)] + yMap_min + [-1000 1000])
	set(handles.axes1,'XTick',[],'YTick',[], 'DataAspectRatio', [1 1 1])

	projStruc.DstProjWKT = ogrproj('+proj=tmerc +lat_0=39.668258333333333 +lon_0=-8.133108611111111 +k=1.0 +x_0=200000 +y_0=300000 +ellps=intl');
	xy_prj = ogrproj([paises.ct.Country(1,:)' paises.ct.Country(2,:)']+0, projStruc);
	line('XData',xy_prj(:,1),'YData',xy_prj(:,2),'Parent',handles.axes1);

	set(hObject,'Vis','on')
	handles.mapa = cell(58,18);
	nr = 1;                     % counter on the number of rectangles
	for (m=58:-1:1)
		yp = [(m-1) (m-1) m m] * 10000 + yMap_min;
		for (n=1:19)            % Loop over columns
 			if (find(naos{n} == m)),    continue,	end
			xp = [(n-1) n n (n-1)] * 16000 + xMap_min; 
			h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',[.7 .7 .7],'LineWidth',1,'Parent',handles.axes1);
			set(h,'ButtonDownFcn',{@bdnTile,hObject},'UserData',[m n])
			handles.mapa{m,n} = sprintf('%.3d',nr);
			text(xp(1)+8000,yp(4)-5000,sprintf('%d',nr),'HorizontalAlignment','center', 'FontSize',5,'HitTest','off','Parent',handles.axes1)
			nr = nr + 1;
		end
		if (~rem(m,5)),			pause(0.01),	end
	end

	% Exceptions to the pure numeric name tile. For example '401A' ---------------------------------------------
    ABs = cell(19,2);
	ABs{1,1} = 22;			ABs{1,2} = {'401A'};
	ABs{2,1} = 19;			ABs{2,2} = {'441B'};
	ABs{3,1} = 30;			ABs{3,2} = {'306B'};
	ABs{4,1} = [4 12 35:37];		ABs{4,2} = {'583A' '515A' '248B' '238A' '227B'};
	ABs{5,1} = 43;			ABs{5,2} = {'162A'};
	ABs{13,1}= [29 30 56];	ABs{13,2}= {'325A' '315A' '9A'};
	ABs{14,1}= [16 20 30];	ABs{14,2}= {'483A' '441A' '315B'};
	ABs{15,1}= [12 31 36];	ABs{15,2}= {'525A' '306A' '248A'};
	ABs{16,1}= [38 40 46];	ABs{16,2}= {'227A' '205A' '142A'};
	ABs{18,1}= 49;			ABs{18,2}= {'108A'};
	ABs{19,1}= 52;			ABs{19,2}= {'67A'};
	
	for (n=1:19)
		xp = [(n-1) n n (n-1)] * 16000 + xMap_min; 
		for (k = 1:numel(ABs{n,1}))
			m = ABs{n,1}(k);
			yp = [(m-1) (m-1) m m] * 10000 + yMap_min;
			h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',[.7 .7 .7],'LineWidth',1,'Parent',handles.axes1);
			set(h,'ButtonDownFcn',{@bdnTile,hObject},'UserData',{ABs{n,2}{k},[m n]})	% Uf, complicated this UD
			text(xp(1)+8000,yp(4)-5000, ABs{n,2}{k}, 'HorizontalAlignment','center', 'FontSize',5,'HitTest','off','Parent',handles.axes1)			
		end
	end
	% -------------------------------------------------------------------------------------------------------------

	% Exceptions TYPE II ------------------------------------------------------------------------------------------
	xp = [76000 83659 83659 76000];		yp = [270000 270000 280125 280125];
	h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',[.7 .7 .7],'LineWidth',1,'Parent',handles.axes1);
	tile.name = '325B';			tile.coords_x = xp;		tile.coords_y = yp;
	set(h,'ButtonDownFcn',{@bdnTile,hObject},'UserData',tile)		% Use a structure in UD to destinguish from the other 2 types
	text(xp(1)+4000,yp(4)-5000, '325B', 'HorizontalAlignment','center', 'FontSize',5,'HitTest','off','Parent',handles.axes1)			

	xp = [96344 104000 104000 96344];		yp = [270000 270000 280000 280000];
	h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',[.7 .7 .7],'LineWidth',1,'Parent',handles.axes1);
	tile.name = '325C';			tile.coords_x = xp;		tile.coords_y = yp;
	set(h,'ButtonDownFcn',{@bdnTile,hObject},'UserData',tile)		% Use a structure in UD to destinguish from the other 2 types

	xp = [9 10 10 9] * 16000 + xMap_min;		yp = [-617 -617 10000 10000];
	h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',[.7 .7 .7],'LineWidth',1,'Parent',handles.axes1);
	tile.name = '611';			tile.coords_x = xp;		tile.coords_y = yp;
	set(h,'ButtonDownFcn',{@bdnTile,hObject},'UserData',tile)		% Use a structure in UD to destinguish from the other 2 types
	% -------------------------------------------------------------------------------------------------------------
	
	% Tooltip
	msg = sprintf(['When you want to get the files directly from Web\n\n' ...
		'As for example at:\nhttp://www.civil.ist.utl.pt/~ruif/HRF/cartas/Cartas%%20Militares%%20GIF_25000/']);
	set(handles.radio_inWeb,'Tooltip',msg)

	% Load the directory list stored in mirone_pref
	load([path_data 'mirone_pref.mat'], 'directory_list');
	set(handles.popup_directory_list,'String', directory_list)
	handles.last_directories = directory_list;
	handles.files_dir = handles.last_directories{1};

	% Initiliaze with this site. This also means that it doesn't recall changes between usage sessions
	set(handles.edit_forWeb,'String', 'http://www.civil.ist.utl.pt/~ruif/HRF/cartas/Cartas%20Militares%20GIF_25000/')

	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end

% -----------------------------------------------------------------------------------------
function bdnTile(obj,event,hFig)
% Do what ever it must to get a referenced image
	handles = guidata(hFig);
	ud = get(gcbo,'UserData');
	proj = [];		knowLimits = false;

	if (isa(ud,'cell'))				% Exception TYPE I (tile names with inbeded letters)
		tile_name = ud{1};
		m = ud{2}(1);		n = ud{2}(2);		% We still need these for computing the corner coordinates
	elseif (isa(ud,'struct'))		% Exception TYPE II (size of tiles not 10x16 km)
		tile_name = ud.name;
		tmp.X = ud.coords_x(1:2);	tmp.Y = ud.coords_y(2:3);
		knowLimits = true;
	else							% Most of the cases
		m = ud(1);			n = ud(2);
		tile_name = handles.mapa{m,n};
	end

	if (get(handles.radio_inloco, 'Value'))
		pato = get(handles.popup_directory_list,'String');
		pato = pato{1};
		if (pato(end) ~= filesep),		pato = [pato filesep];		end
		fname = [pato tile_name '.gif'];
		if (exist(fname,'file') ~= 2)
			fname = [pato tile_name '.sid'];		% Try in Mr. Sid format (MUST ADD A TEST IF MRSID DRIVER EXISTS)
		end
		if (exist(fname,'file') ~= 2)
			h = text(100000,50000, ['FILE  ' tile_name '.GIF  NOT FOUND IN THERE'], 'HorizontalAlignment','left', 'FontSize',25, 'Rotation', 65, 'Color','r');
			pause(2),	delete(h)
			return
		end
		fnameRef = [pato tile_name '.map'];
	else
		url = get(handles.edit_forWeb,'String');
		if (url(end) ~= '/'),			url = [url '/'];		end
		if (~strncmp(url,'http://',7)),	url = ['http://' url];	end
		fname = [url tile_name '.GIF'];
		fnameRef = [];			% We need to find a way to fetch the corresponding .map file
	end
		
	set(handles.figure1,'pointer','watch'),		pause(0.1)
	[img,att] = gdalread(fname, '-U');
	if (isempty(img))
		errordlg('There was an error reading the image. We got nothing.','Error')
		return
	end

	if (~isempty(fnameRef) && exist(fnameRef,'file') == 2)			% Found a Ozi .map reference file
		[gcp, proj] = getProjFromOzi(fnameRef);
		% Compute and apply the affine transformation
		trans = AffineTransform(gcp(:,1:2),gcp(:,3:4));
		x_pt = [1; size(img,2)];    y_pt = [1; size(img,1)];
		X1 = [x_pt y_pt ones(size(x_pt,1),1)];
		U1 = X1 * trans;
		tmp.X = U1(:,1)';		tmp.Y = U1(:,2)';
		if (tmp.Y(1) > tmp.Y(2))			% They just love that stupid origin in the UL corner
			yy = tmp.Y(1);		tmp.Y(1) = tmp.Y(2);	tmp.Y(2) = yy;
		end
	elseif (~isempty(fnameRef) && ~isempty(att.GeoTransform))		% We had a world file
		tmp.X = [att.Corners.LL(1) att.Corners.LR(1)];
		tmp.Y = [att.Corners.LL(2) att.Corners.UL(2)];
	else
		if (~knowLimits)		% Otherwise we already know them
			tmp.X = [(n-1) n] * 16000 + 72000;
			tmp.Y = [(m-1) m] * 10000;
		end
		% Set the proj string for the "Coordenadas militares, datum Lisboa"
		proj = '+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111111 +k=1.0 +x_0=200000 +y_0=300000 +ellps=intl +towgs84=-304.046,-60.576,103.640,0,0,0,0';
	end

	tmp.geog = 0;       tmp.name = fname;
	x_inc = diff(tmp.X) / (size(img,2) - 1);
	y_inc = diff(tmp.Y) / (size(img,1) - 1);
	tmp.head = [tmp.X tmp.Y att.GMT_hdr(5:6) 0 x_inc y_inc];
	
	try
		if (ndims(img) ~= 3)	% The Sided version are rgb
			tmp.cmap = att.Band.ColorMap.CMap(:,1:3);
		end
	catch
		errordlg(['Bad colormap in file ' fname],'ERROR')
		set(handles.figure1,'pointer','arrow')
		return
	end
	if (~isempty(proj))			% If we know the projection, convert it to the WKT form
		tmp.srsWKT = ogrproj(proj);
	end
	mirone(img, tmp);
	set(handles.figure1,'pointer','arrow')

% ----------------------------------------------------------------------------
function popup_directory_list_CB(hObject, handles, opt)
% OPT is used by push_change_dir (just to save code)
	if (nargin == 2),	opt = [];   end
	if isempty(opt)
		val = get(hObject,'Value');     str = get(hObject, 'String');
		% Put the selected field on top of the String list. This is necessary because the "OK" button will
		tmp = str(val);         str(val) = [];
		new_str = [tmp; str];   set(hObject,'String',new_str); 
		set(hObject,'Value',1)
		if iscell(tmp),			new_dir = tmp{1};
		elseif ischar(tmp),		new_dir = tmp;
		else                    return        % ???
		end
	else
		new_dir = opt;
	end
	handles.files_dir = new_dir;
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------
function push_change_dir_CB(hObject, handles)
	if (strcmp(computer, 'PCWIN'))
		work_dir = uigetfolder_win32('Select a directory', cd);
	else            % This guy doesn't let to be compiled
		work_dir = uigetdir(cd, 'Select a directory');
	end
	if (isempty(work_dir)),		return,		end
	handles.last_directories = [cellstr(work_dir); handles.last_directories];
	set(handles.popup_directory_list,'String',handles.last_directories)
	guidata(hObject, handles);
	popup_directory_list_CB(handles.popup_directory_list, handles, work_dir)

% ----------------------------------------------------------------------------
function radio_inloco_CB(hObject, handles)
	if (~get(hObject,'Val')),		set(hObject,'Val',1),	return,		end
	set(handles.radio_inWeb,'Val',0)
	set(handles.edit_forWeb,'Vis','off')
	set([handles.popup_directory_list handles.push_change_dir], 'Vis','on')

% ----------------------------------------------------------------------------
function radio_inWeb_CB(hObject, handles)
	if (~get(hObject,'Val')),		set(hObject,'Val',1),	return,		end
	set(handles.radio_inloco,'Val',0)
	set(handles.edit_forWeb,'Vis','on')
	set([handles.popup_directory_list handles.push_change_dir], 'Vis','off')

% ----------------------------------------------------------------------------
function [gcp, proj] = getProjFromOzi(fname)
	% Decode relevant info in a .map Ozi reference file 
	fid = fopen(fname);
	c = fread(fid,inf,'*char');
	fclose(fid);
	map = strread(c,'%s','delimiter','\n');   clear c fid;
	
	% Now get the GCPs
	ind = strfind(map{10},',');
	% We are decoding strings like this one
	%Point01,xy,  337,  336,in, deg,    ,        ,N,    ,        ,W, grid,   ,     185000,     579000,N
	gcp = zeros(30,4);
	k = 1;
	while ( map{k+9}(ind(3)-1) ~= ' ')		% While we have a pixel<->coord point pair
		gcp(k,1) = str2double( map{k+9}((ind(2)+1) : (ind(3)-1)) ); 
		gcp(k,2) = str2double( map{k+9}((ind(3)+1) : (ind(4)-1)) ); 
		gcp(k,3) = str2double( map{k+9}((ind(14)+1) : (ind(15)-1)) ); 
		gcp(k,4) = str2double( map{k+9}((ind(15)+1) : (ind(16)-1)) ); 
		k = k + 1;
	end

	% Remove non-used pre-allocated gcp
	gcp = gcp(1:(k-1),:);

	% Get mapp projection - THE ONLY TWO CASES THAT I SAW. OBVIOUSLY THIS SHOULD NOT BE APPLYED TO OTHER OZI FILES
	ind = strfind(map{9},',');
	proj = '';
	if ( strncmp(map{9}( ind(1)+1:ind(2)-1 ),'Transverse',10) )		% Transverse Mercator
		proj = '+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111111 +k=1.0 +x_0=200000 +y_0=300000 +ellps=intl +towgs84=-302.581,-61.360,103.047,0,0,0,0';
	elseif ( strncmp(map{9}( ind(1)+1:ind(2)-1 ),'(UTM)',5) )		% (UTM) Universal Transverse Mercator
		proj = '+proj=utm +zone=29 +k=0.9996 +ellps=intl +towgs84=-85.858,-108.681,-120.361,0,0,0,0';
	end
	
% -----------------------------------------------------------------------------------------
function trans = AffineTransform(uv,xy)
% For an affine transformation:
%                     [ A D 0 ]
% [u v 1] = [x y 1] * [ B E 0 ]
%                     [ C F 1 ]
% There are 6 unknowns: A,B,C,D,E,F
% Another way to write this is:
%                   [ A D ]
% [u v] = [x y 1] * [ B E ]
%                   [ C F ]
% Rewriting the above matrix equation:
% U = X * T, where T = reshape([A B C D E F],3,2)
%
% With 3 or more correspondence points we can solve for T,
% T = X\U which gives us the first 2 columns of T, and
% we know the third column must be [0 0 1]'.

	K = 3;      M = size(xy,1);     X = [xy ones(M,1)];
	U = uv;         % just solve for the first two columns of T

	% We know that X * T = U
	if rank(X) >= K
		Tinv = X \ U;
	else
		msg = 'At least %d non-collinear points needed to infer %s transform.';
		errordlg(sprintf(msg,K,'affine'),'Error');
	end

	Tinv(:,3) = [0 0 1]';       % add third column
	trans = inv(Tinv);
	trans(:,3) = [0 0 1]';


% --- Creates and returns a handle to the GUI figure. 
function cartas_militares_LayoutFcn(h1)
set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'MenuBar','none',...
'Name','Cartas Militares',...
'NumberTitle','off',...
'Pos',[520 113 331 670],...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','Call',...
'Tag','figure1');

axes('Parent',h1,...
'Units','pixels',...
'Box','on',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'Pos',[0 50 331 622],...
'Tag','axes1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Pos',[3 3 325 22],...
'Style','edit',...
'Tooltip','Enter web site address where we can find the image files',...
'Tag','edit_forWeb',...
'Visible','off');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'Pos',[10 3 300 22],...
'Style','popupmenu',...
'String', {' '},...
'Tooltip','Select the directory where the image files reside',...
'Value',1,...
'Tag','popup_directory_list');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'Pos',[309 3 21 23],...
'String','...',...
'Tooltip','Select a different directory',...
'Tag','push_change_dir');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'FontName','Helvetica',...
'Pos',[10 30 61 15],...
'String','In loco',...
'Style','radiobutton',...
'Tooltip','When you have the files on disk ',...
'Value',1,...
'Tag','radio_inloco');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'FontName','Helvetica',...
'Pos',[130 30 71 15],...
'String','In Web',...
'Style','radiobutton',...
'Tooltip','When you have get the files from Web',...
'Tag','radio_inWeb');

function main_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

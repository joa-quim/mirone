function  varargout = aux_funs(opt,varargin)
% This contains Mirone's auxiliay functions that are called by the several
% of the Mirone's callback functions. I puted them here to release somehow
% the burden of the non-stop groing length of the Mirone code.

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

switch opt(1:4)
	case 'Stor'		% 'StoreZ'
		StoreZ(varargin{:})
	case 'msg_'		% 'msg_dlg'
		varargout = msg_dlg(varargin{:});
	case 'in_m'		% 'in_map_region'
		[varargout{1:nargout}] = in_map_region(varargin{:});
	case 'clea'		% 'cleanGRDappdata'
		clean_GRDappdata(varargin{:})
	case 'gues'		% 'guessGeog'
		varargout{1} = guessGeog(varargin{:});
	case 'colo'		% 'colormap_bg'
		colormap_bg(varargin{:})
	case 'find'		% 'findFileType'
		varargout{1} = findFileType(varargin{:});
		if (nargout == 2)       % Check if file exists
			ind = (exist(varargin{1},'file') == 2);
			if (~ind)
				[PATH,FNAME,EXT] = fileparts(varargin{1});
				newFilePato = [cd filesep FNAME EXT];
				ind = (exist(newFilePato,'file') == 2);
				if (ind),		ind = newFilePato;		end		% File exists on the directory from where Mirone was called
			end
			varargout{2} = ind;		% Logical if file exists|not exists on Mirone's root directory. Otherwise -> actual full file path
		end
	case 'appP'		% 'appProjectionRef'
		appProjectionRef(varargin{:})
	case 'isPr'		% 'isProj'
		varargout{1} = isProj(varargin{:});
	case 'toPr'		% 'toProjPT'
		toProjPT(varargin{:});
	case 'togC'		% 'togCheck'			% Toggle the check state of a uimenu
		togCheck(varargin{:})
	case 'poly'		% 'polysplit'
		[varargout{1} varargout{2}] = localPolysplit(varargin{:});
	case 'adju'		% 'adjust_lims'
		[varargout{1} varargout{2}] = adjust_lims(varargin{:});
	case 'insi'		% 'insideRect'
		varargout{1} = insideRect(varargin{:});
	case 'stri'		% 'strip_bg_color'
		varargout{1} = strip_bg_color(varargin{:});
	case 'help'
		str = [varargin{1}.home_dir filesep 'doc' filesep 'MironeMan.chm &'];
		if (ispc),      dos(str);
		else            unix(str);
		end
	case 'getF'		% Get projection info (if it's there)
 		varargout{1} = getFigProjInfo(varargin{:});
	case 'min_'		% 'min_max_single'
		[varargout{1:nargout}] = min_max_single(varargin{:});
	otherwise
		if (nargout)
			[varargout{1:nargout}] = ...
				feval(opt, varargin{:});	% NEW. All calls should evolve to use this form
		else
			feval(opt, varargin{:});
		end
end

% --------------------------------------------------------------------
function addUI(handles)
% Add a uimenu. For the time beeing this is used only in one situation so no need for further options
	h = findobj(handles.figure1, 'type','uimenu','Label','Save GMT script');
	hh = findobj(h, 'Label', '.def symbol');		% Check if we already have it or not
	if (isempty(hh))
		uimenu('Parent',h,'Call','write_gmt_symb(guidata(gcbo))','Label','.def symbol','Sep','on');
	end

% --------------------------------------------------------------------
function StoreZ(handles,X,Y,Z)
% If grid size is not to big I'll store it
	if (numel(Z)*4 > handles.grdMaxSize),		return,		end
	if (~handles.IamCompiled && ~isa(Z,'single')),	setappdata(handles.figure1,'dem_z',single(Z));	% TMP
	else					setappdata(handles.figure1,'dem_z',Z);
	end
	setappdata(handles.figure1,'dem_x',X);  setappdata(handles.figure1,'dem_y',Y);

% --------------------------------------------------------------------
function [x,y,indx,indy,hdr_str] = in_map_region(handles, x, y, tol, map_lims, hdr_str)
%   Given X & Y vectors retain only the elements that are inside the current map region
%   OPTIONS:
%   TOL is used normally when ploting lines and serves to extend the map
%       limits so that lines are allowed to be drawn until the image borders
%   MAP_LIMS a 1x4 vector with [x_min x_max y_min y_max]. If not given it will be fetch here
%   NOTE THAT ALL ARGUMENTS MUST BE PROVIDED, EVEN IF THEY ARE EMPTY

	if (isempty(tol)),	tol = 0;		end
	if (nargin == 5),	hdr_str = [];	end
	if (isempty(map_lims))
        x_lim = get(handles.axes1,'Xlim');      y_lim = get(handles.axes1,'Ylim');
	else
        x_lim(1:2) = map_lims(1:2);    y_lim(1:2) = map_lims(3:4);
	end
	if (handles.geog == 2 && x_lim(2) > 180)			% If basemap is in the [0 360] range
		indx = (x < 0);									% and we may need the wrapping. Do it.
		if (isa(x,'double'))
			x(indx) = x(indx) + 360;
		else
			x(indx) = single(double(x(indx)) + 360);	% Very bad luck if they are integers
		end
	end

	if (tol >= 0)
		indx = find((x < x_lim(1)-tol) | (x > x_lim(2)+tol));
		x(indx) = [];           y(indx) = [];
		if (nargout == 2),      clear indx;     end         % Save memory
		indy = find((y < y_lim(1)-tol) | (y > y_lim(2)+tol));
		x(indy) = [];           y(indy) = [];
		set(handles.figure1,'CurrentAxes',handles.axes1)	% This is for the GCP mode be able to plot on the Master Image
	else
		indx = [];		indy = [];				% Set this while we don't do any cipping
		if (nargin == 6)						% Very special case of GMT_DB polygons
			ind_Gi = strfind(hdr_str, 'G = ');		ind_Ge = strfind(hdr_str, ' L =');	
			%ind_Ei = strfind(hdr_str, 'E = ');		ind_Ee = strfind(hdr_str, ' S =');	
			ind_Ri = strfind(hdr_str, 'R = ');		ind_Re = strfind(hdr_str, ' A =');
			pol_lims = sscanf(hdr_str(ind_Ri+4:ind_Re-1), '%f/%f/%f/%f')';	% This polygon limits
			green = sscanf(hdr_str(ind_Gi+4:ind_Ge-1), '%d');				% Greenwich and/or Dateline crossing
			%datelon = sscanf(hdr_str(ind_Ei+4:ind_Ee-1), '%d');
			reco = rectangle_and(map_lims, pol_lims);
			if (isempty(reco) && ~green)		% Try again by checking if there is a [-180 180] [0 360] ambiguity
				if (handles.geog == 1)			% lon [-180 180]
					reco = rectangle_and(map_lims, pol_lims - [360 360 0 0]);
					shift = -360;
				else							% DON'T LIKE THIS CASE. IT SHOULD NEVER OCCUR
					reco = rectangle_and(map_lims, pol_lims + [360 360 0 0]);
					shift = 360;
				end
				if (~isempty(reco))				% Ah ah, got one. Wrap it to the visible side of the world
					x = x + shift;
					hdr_str = sprintf('%s%f/%f/%f/%f%s', ...
						hdr_str(1:ind_Ri+3), pol_lims(1:2)+shift,pol_lims(3:4), hdr_str(ind_Re:end));
				else
					x = [];		y = [];			% No interception
				end
			end
		elseif ~( min(x) >= x_lim(1) || max(x) <= x_lim(2) || min(y) >= y_lim(1) || max(y) <= y_lim(2) )
			x = [];		y = [];		% No interception
			% there is no 'else' case while the rectangle cipping is not implemented
		end
	end

% --------------------------------------------------------------------
function [rect,rect_crop] = rectangle_and(head1, head2)
%   Given two handles.head vectors, compute the intersection rectangle of the two regions
%	RECT id a 5x2 matrix with rectangle coords X in first col and Y in second
%	RECT_CROP is a row vector with [Xll Yll DX DY] (Xll -> lower left point)
	P1.x = [head1(1) head1(1) head1(2) head1(2) head1(1)];	P1.hole = 0;
	P1.y = [head1(3) head1(4) head1(4) head1(3) head1(3)];
	P2.x = [head2(1) head2(1) head2(2) head2(2) head2(1)];	P2.hole = 0;
	P2.y = [head2(3) head2(4) head2(4) head2(3) head2(3)];
	P3 = PolygonClip(P1, P2, 1);				% Intersection of the two rectangles
	if (~isempty(P3))
		x_min = min(P3.x);		x_max = max(P3.x);		y_min = min(P3.y);		y_max = max(P3.y);
		rect = [x_min y_min; x_min y_max; x_max y_max; x_max y_min; x_min y_min];
		if (nargout == 2)
			rect_crop = [x_min y_min x_max-x_min y_max-y_min];
		end
	else
		rect = [];		rect_crop = [];
	end

% --------------------------------------------------------------------
function colormap_bg(handles,Z,pal)
% Insert the background color in the palette for arrays that have NaNs
% [m,n] = size(Z);    dbl_size = m*n*8;
% if (dbl_size > handles.grdMaxSize)  % I'm almost shure that isnan uses doubles
%     n_stripes = round(dbl_size / handles.grdMaxSize) * 2;
%     d_stripe = fix(m / n_stripes);
%     is_it = 0;      tmp = 0;
%     for (i=1:n_stripes-1)
%         is = (i-1)*d_stripe + 1;    ie = is + d_stripe;
%         tmp = any(isnan(Z(is:ie,:)));
%         tmp = any(tmp);     % Make it a single value
%         is_it = is_it + tmp;
%     end
%     tmp = isnan(Z(ie+1:end,:));     tmp = any(tmp);
%     is_it = is_it + tmp;
%     if (is_it)              pal = [handles.bg_color; pal];   end
% else
%     if any(isnan(Z(:)))     pal = [handles.bg_color; pal];   end
% end
% tic;    is_it = any(Z(:)~=Z(:)); toc
% tic;    is_it = any(isnan(Z(:))); toc
% return

%if any(Z(:)~=Z(:))     pal = [handles.bg_color; pal];   end
if ( handles.have_nans && ~isequal(pal(1,:),handles.bg_color) )
    if (size(pal,1) == 256),    pal = [handles.bg_color; pal(2:end,:)];     % Jump firts color to not have
    else                        pal = [handles.bg_color; pal];              % a CMAP with more than 256 colors
    end
end
set(handles.figure1,'Colormap',pal)

% ----------------------------------------------------------------------------------
function out = findFileType(fname)
% From the extension guess what function should be called to open this file
	out = [];	
	if (isempty(fname))		return,		end
	[PATH,FNAME,EXT] = fileparts(fname);
	if ( strcmpi(EXT,'.grd') )
		out = 'gmt';
	elseif ( strcmpi(EXT,'.nc') )		% .nc files can have grids, mgd77 files or any other thing
		s = nc_funs('info',fname);
		try
			if     ( any(strcmp({s.Dimension.Name}, 'id_dim')) ),			out = 'mgg_gmt';
			elseif ( any(strcmp({s.Attribute.Name}, 'SHAPENC_type')) ),		out = 'ncshape';
			else	out = 'gmt';
			end
		catch
			out = 'dono';
		end
		if (out(1) == 'g')				% But the story is not over. We need to know if multi-2D_datasets
			count = 0;
			for (k = 1:numel(s.Dataset))
				if (numel(s.Dataset(k).Size) > 1)
					count = count + 1;
					if (count > 1),		out = 'dono';	break,	end		% Yes, so we experimentaly send it to gdalread
				end
			end
		end
	elseif ( any(strcmpi(EXT,{'.jpg' '.png' '.bmp' '.gif' '.pcx' '.ras' '.ppm' '.pgm' '.xwd' '.shade' '.raw' '.bin'})) )
		out = 'generic';
	elseif ( any(strcmpi(EXT,{'.tif' '.tiff' '.sid' '.kap' '.nos'})) )
		out = 'geotif';
	elseif ( any(strcmpi(EXT,{'.ecw' '.jp2'})) )	% This is a special case (they cause memory fragmentation)
		out = 'ecw';
	elseif ( strcmpi(EXT,'.mat') )
		warning('off','MATLAB:load:variableNotFound')
		s = load(fname,'grd_name');
		if (isfield(s, 'grd_name'))
			out = 'mat';
		else
			s = load(fname,'FitLine');	% Try if it is a Ecran session
			if (isfield(s, 'FitLine')),	out = 'mat';	end
		end
	elseif ( any(strcmpi(EXT,{'.n1' '.n14' '.n15' '.n16' '.n17'})) )
		out = 'multiband';
	elseif ( strcmpi(EXT,'.img') )
		nome = [PATH filesep FNAME '.lbl'];
		if (exist(nome, 'file'))
			out = 'mola';
		else
			out = 'envherd';
		end
	elseif (strcmpi(EXT,'.cpt'))
		out = 'cpt';
	elseif ( any(strcmpi(EXT,{'.dat' '.xy' '.b' '.txt'})) )
		out = 'dat';
	elseif (strcmpi(EXT,'.shp'))
		out = 'shp';
	elseif ( any(strcmpi(EXT,{'.las' '.laz'})) )
		out = 'las';
	elseif (strcmpi(EXT,'.gmt'))
		fid = fopen(fname,'rt');
		ID = fread(fid,7,'*char');      ID = ID';      fclose(fid);
		if (strncmp(ID,'# @VGMT', 7)),	out = 'ogr';
		else							out = 'mgg_gmt';
		end
	elseif ( any(strcmpi(EXT,{'.kml' '.gml' '.dxf' '.gpx' '.dgn' '.csv' '.s57' '.svg'})) )
		out = 'ogr';
	elseif (strcmpi(EXT,'.srtm'))	% While we don't use GMT5, create a header write away & send to GDAL
		try,	write_esri_hdr(fname,'SRTM30');	end
		out = 'dono';				% aka GDAL
	elseif (strcmpi(EXT,'.sww'))	% Might be an ANUGA netCDF. Confirm it
		fid = fopen(fname,'r');
		ID = fread(fid,3,'*char');      ID = ID';      fclose(fid);
		out = 'dono';
		if (strcmpi(ID,'CDF') || strcmpi(ID,'HDF'))
			s = nc_funs('info',fname);
			ind = strcmp({s.Dimension.Name},'number_of_volumes');
			if (any(ind)),		out = 'sww';	end		% Yes, it's ANUGA file.
		end
	else
		% OK, here we'll send ASCII files to load_xyz (after one more test) and binary to GDAL ... and see what happens.
		bin = guess_file(fname);
		if (isempty(bin))			% Should be a "ERROR READING FILE"
			out = 'boom';			% Unknow from caller, but it will trigger the error message
			try						% Do a silent test trying GDAL first and if fails send to BOOM
				a.att = gdalread(fname,'-M');
				out = 'dono';		% It will go try luck with GDAL
			end
		elseif (isa(bin,'struct') || bin == 0)
			try						% Do a silent test trying GDAL first and if fails send to load_xyz
				a.att = gdalread(fname,'-M');
				out = 'dono';		% It will go try luck with GDAL
			catch
				% Here in either case (BIN or ASCII) we send it to try luck with load_xyz
				% but in future I should try what OGRREAD has to say on this matter.
				out = 'dat';
			end
		else
			out = 'dono';			% It will go try luck with GDAL
		end
	end

% --------------------------------------------------------------------
function out = msg_dlg(in,handles)
% Cast a window with a pre-deffined message
% Alternatively, if IN == 50 and handles.no_file = true, create a Global background image and return no error

	out = {0};    msg = [];     h = [];		maybe_create = false;
	if (in == 50)
		maybe_create = true;
		in = 5;					% So that the old case '5' below work as before
	end
	if (handles.no_file)
		if (maybe_create)		% Create a Global bg image
			img = gdalread([handles.home_dir filesep 'data' filesep 'etopo4.jpg'], '-U');
			x_inc = 360 / (size(img,2)-1);			y_inc = 180 / (size(img,1)-1);
			handles.geog = 1;			handles.image_type = 3;
			handles.head = [-180 180 -90 90 0 255 0 x_inc y_inc];
			mirone('show_image', handles,'Base image',[-180 180],[-90 90],img,0,'xy',1,1);	% Is also saves handles
			drawnow
			return
		else
			msg = 'This option requires that you load a file first to serve as background map.';    out = {1};
			msgbox(msg,'Error')
			return
		end
	end
	switch in
		case 1
			if (~handles.geog)
				msg = 'This operation is currently possible only for geographic type data';
				out = {1};
			end
		case 14
			if (~handles.validGrid)
				msg = 'This operation is deffined only for images derived from float (eg DEMs) grids';
				out = {1};
			end
		case 5
			if (~handles.is_projected && ~handles.geog)
				msg = 'This operation is only possible for geographic data OR when the Map Projection is known';
				out = {1};            
			end
	end
	if (~isempty(msg))
		h = msgbox(msg,'Error');
	end
	axes(handles.axes1)     % This is for the GCP mode be able to plot on the Master Image
	if (~isempty(h))        % Bring the figure forward because it was hiden by the previous command
		figure(h)
	end

% ----------------------------------------------------------------------------------
function handles = isProj(handles, opt)
% Se if we have projected grid/images and act acordingly

	if (~ishandle(handles.Projections)),	return,		end		% True when GCPtool
	
	% Fish eventual proj strings
	prjInfoStruc = getFigProjInfo(handles);
    
	if (~handles.geog)
		if (~isempty(prjInfoStruc.projWKT))              % We have a GDAL projected file
			handles.is_projected = 1;		prjStr = false;
		elseif (~isempty(prjInfoStruc.projGMT) || ~isempty(prjInfoStruc.proj4))  % We have a GMT or Proj4 projection selection
			handles.is_projected = 1;		prjStr = false;
		else                                % We know nothing about these coords (e.g. an image)
			handles.is_projected = 0;		prjStr = true;
		end
	else
		handles.is_projected = 0;			prjStr = false;
	end
	if (handles.is_projected),      set(handles.hAxMenuLF, 'Vis', 'on', 'Separator','on')	% To load/not in prj coords
	else                            set(handles.hAxMenuLF, 'Vis', 'off', 'Separator','off')
	end

	child = get(handles.Projections,'Children');
	h1 = findobj(handles.Projections,'-depth',1,'Label','GDAL project');
	h2 = findobj(handles.Projections,'-depth',1,'Label','GMT project');
	h3 = findobj(handles.Projections,'-depth',1,'Label','Point projections');
	if (~prjStr)
		set(child,'Enable','off')
		set([h1 h2 h3],'Enable','on')		% These are always on
	end

	if (handles.is_projected)       % When read a gdal referenced file we need to set this right away
		toProjPT(handles)
		setappdata(handles.figure1,'DispInGeogs',0)     % We need to set it now (this is the default)
	end

	if (nargin == 2),       guidata(handles.figure1, handles);      end

% ----------------------------------------------------------------------------------
function togCheck(varargin)
% Toggles the check state of a uimenu. If nargin == 2, second arg is assumed to contain
% a handles list of others uimenus belonging to a group. In this case, if one ui is 
% checked on, the others will be checked off. The contrary does not apply.
	h = varargin{1};	st = {'on', 'off'};
	id = 0;
	if (strcmp(get(h,'check'), 'on')),	id = 1;		end
	if (nargin == 1)				% Pure toggle on/off
		set(h,'check', st{id+1})
	else							% "Radiobutton mode"
		if (~id)					% IF NOT CHECKED, check arg1 and un-check the others in arg2
			set(varargin{2}, 'check', 'off')
			set(h, 'check', 'on')
		end
	end

% ----------------------------------------------------------------------------------
function toProjPT(handles)
% The following is to deal with eventual coords display in geogs and is used in PIXVAL_STSBAR 

	displayBar = findobj(handles.figure1, 'Tag', 'pixValStsBar');
	if (isempty(displayBar)),	return,		end		% No file, no projection
	dbud = get(displayBar, 'UserData');
	dbud.toProjPT = 0;
	if (handles.is_projected)
		% Fish eventual proj strings
		prjInfoStruc = getFigProjInfo(handles);
		if (~isempty(prjInfoStruc.projWKT) || ~isempty(prjInfoStruc.proj4))
			% In case we have both 'projWKT' takes precedence because it came from file metadata
			if (~isempty(prjInfoStruc.projWKT)),	dbud.projStruc.SrcProjWKT = prjInfoStruc.projWKT;
			else									dbud.projStruc.SrcProjWKT = ogrproj(prjInfoStruc.proj4);
			end
			dbud.toProjPT = 1;
		elseif (~isempty(prjInfoStruc.projGMT))
			lims = getappdata(handles.axes1,'ThisImageLims');
			out = mapproject_m([lims(1) lims(3); lims(2) lims(4)],'-R-180/180/0/80','-I','-F',prjInfoStruc.projGMT{:});    % Convert lims back to geogs
			x_min = min(out(:,1));        x_max = max(out(:,1));
			y_min = min(out(:,2));        y_max = max(out(:,2));
			dbud.opt_R = sprintf('-R%f/%f/%f/%f',x_min, x_max, y_min, y_max);
			dbud.projGMT = prjInfoStruc.projGMT;
			dbud.toProjPT = 2;
		end
		set(displayBar, 'UserData', dbud)
	end

% ----------------------------------------------------------------------------------
function toBandsList(hFig, I, array_name, fname, n_bands, bands_inMemory, reader)
% Create a structure with data to be retrieved by the BANDS_LIST() GUI
% This still needs to be improved. Not much of error testing
	if (nargin == 3)
		fname = [];		n_bands = size(I,3);	bands_inMemory = 1:n_bands;		reader = [];
	end
	if (isempty(n_bands)),			n_bands = size(I,3);		end
	if (isempty(bands_inMemory)),	bands_inMemory = 1:n_bands;		end
	tmp1 = cell(n_bands+1,2);		tmp2 = cell(n_bands+1,2);
	tmp1{1,1} = array_name;			tmp1{1,2} = array_name;
	for (i = 1:n_bands)
		tmp1{i+1,1} = sprintf('band%d',i);
		tmp1{i+1,2} = sprintf('banda%d',i);			% TEMP
		tmp2{i+1,1} = [sprintf('%d',i) 'x1 bip'];	% TEMP
		tmp2{i+1,2} = i;
	end
	tmp = {['+ ' array_name]; I; tmp1; tmp2; fname; bands_inMemory; [size(I,1) size(I,2) n_bands]; reader};
	setappdata(hFig,'BandList',tmp)

% ----------------------------------------------------------------------------------
function prjInfoStruc = getFigProjInfo(handles)
% Se if we have projection info stored in Figure's appdata. NOTE, often they are empty
	prjInfoStruc.projGMT = getappdata(handles.figure1,'ProjGMT');
	prjInfoStruc.projWKT = getappdata(handles.figure1,'ProjWKT');
	prjInfoStruc.proj4 = getappdata(handles.figure1,'Proj4');

% --------------------------------------------------------------------
function [z_min,z_max] = min_max_single(Z)
% Compute the min/max of single precision Z arrays. I need this due to (another) Matlab
% bug that gives wrong results when the Z (single) array has NaNs. Ouput are doubles.
	z_min = double(min(Z(~isnan(Z(:)))));   z_max = double(max(Z(~isnan(Z(:)))));

% --------------------------------------------------------------------
function img = strip_bg_color(handles,img)
% Strip eventual row/columns with color equal to the Figure's background color
	bg_color = uint8(get(handles.figure1,'color')*255);
	c1 = bg_color(1);    c2 = bg_color(2);    c3 = bg_color(3);
	center_row = round(size(img,1) / 2);
	center_col = round(size(img,2) / 2);
	% Strip north
	i = 1;
	while (img(i,center_col,1) == c1 && img(i,center_col,2) == c2 && img(i,center_col,3) == c3)
		i = i + 1;
	end
	i = i - 1;
	img(1:i,:,:) = [];      % Strip the eventual gray band
	% Strip west
	j = 0;
	while (img(center_row,end-j,1) == c1 && img(center_row,end-j,2) == c2 && img(center_row,end-j,3) == c3)
		j = j + 1;
	end
	j = j - 1;
	img(:,end-j:end,:) = [];      % Strip the eventual gray band
	% Strip south
	i = 0;
	while (img(end-i,center_col,1) == c1 && img(end-i,center_col,2) == c2 && img(end-i,center_col,3) == c3)
		i = i + 1;
	end
	i = i - 1;
	img(end-i:end,:,:) = [];      % Strip the eventual gray band
	% Strip east
	j = 1;
	while (img(center_row,j,1) == c1 && img(center_row,j,2) == c2 && img(center_row,j,3) == c3)
		j = j + 1;
	end
	j = j - 1;
	img(:,1:j,:) = [];      % Strip the eventual gray band

% --------------------------------------------------------------------
function [X,Y] = adjust_lims(X,Y,m,n)
% Convert the image limits from pixel reg to grid reg
	dx = (X(2) - X(1)) / n;			dy = (Y(2) - Y(1)) / m;
	X(1) = X(1) + dx/2;				X(2) = X(2) - dx/2;
	Y(1) = Y(1) + dy/2;				Y(2) = Y(2) - dy/2;

% --------------------------------------------------------------------
function geog = guessGeog(lims)
% Make a good guess if LIMS are geographic
	geog = double( ( (lims(1) >= -180 && lims(2) <= 180) || (lims(1) >= 0 && lims(2) <= 360) )...
		&& (lims(3) >= -90 && lims(4) <= 90) );
	if (geog && lims(2) > 180)
		geog = 2;			% We have a [0 360] range
	end

% --------------------------------------------------------------------
function res = insideRect(rect,pt)
% Check which elements of the  [x y] (Mx2) PT array are inside the rectangle RECT
% RECT = [x_min x_max y_min y_max]
% RES is a logical column vector with length = size(PT,1)
% NO ERROR TESTING
    res = ( pt(:,1) >= rect(1) & pt(:,1) <= rect(2) & pt(:,2) >= rect(3) & pt(:,2) <= rect(4) );

% --------------------------------------------------------------------
function clean_GRDappdata(handles)
% If reistering an image against a grid, those cannot be empty (shity solution)
	try
		rmappdata(handles.figure1,'dem_x');     rmappdata(handles.figure1,'dem_y');
		rmappdata(handles.figure1,'dem_z');
	end

% --------------------------------------------------------------------------------
function [track, names, names_ui, vars, x_min, x_max, y_min, y_max] = get_mgg(names, PathName, varargin)
% Get tracks from either the old style .gmt format or new MGD77+ netCDF format

	what2plot = seek_OPTcontrol('MIR_MGG');		% See if OPTcontrol has a request other than the default (MAG)
	n_column = 1;
	[t,r] = strtok(names{1});
	if (~isempty(r)),	n_column = 2;	end

	m = numel(names);		c = false(m,1);
	if (n_column > 1)							% When 2nd column holds the variable to plot info
		lixo = cell(m,1);
		for (k = 1:m)
			[t,r] = strtok(names{k});
			if (t(1) == '#' || isempty(t))		% Jump empty and comment lines
				c(k) = true;	continue
			end
			names{k} = t;
			r = ddewhite(r);
			lixo{k} = r;
		end
		if (any(c)),	names(c) = [];		lixo(c) = [];	end		% Remove eventual comment lines
	else								% Only one column with fnames
		lixo = [];
		for (k = 1:m)
			if ( isempty(names{k}) || names{k}(1) == '#')
				c(k) = true;		continue
			end
		end
		if (any(c)),	names(c) = [];	end		% Remove eventual comment lines
	end

	m = numel(names);		c = false(m,1);		% Count remaining ones
	vars = cell(m,3);		names_ui = names;

	for (k = 1:numel(names))			% Rip the .??? extension
		[PATH, FNAME, EXT] = fileparts(names{k});
		names{k} = FNAME;		names_ui{k} = [FNAME EXT];
		if (~isempty(PATH))				% File names in the list have a path
			names{k} = [PATH filesep names{k}];
		else							% They do not have a path, but we need it. So prepend the list-file path
			names{k} = [PathName names{k}];
		end
		if (exist([names{k} EXT],'file') ~= 2),			c(k) = true;		continue,	end		% File does not exist
		if (~isempty(lixo))					% Scan second column for the names of the variables to plot
			ind = strfind(lixo{k}, ',');	% See how many required fields
			if (isempty(ind))				% Only one. It will be ploted in the 'faa' axes
				vars{k,1} = lixo{k};	vars(k,2:3) = {'mtf1' 'depth'};
			elseif (numel(ind) == 1)
				vars{k,1} = lixo{k}(1:ind(1)-1);		vars{k,2} = lixo{k}(ind(1)+1:end);		vars{k,3} = 'depth';
			elseif (numel(ind) == 2)
				vars{k,1} = lixo{k}(1:ind(1)-1);		vars{k,2} = lixo{k}(ind(1)+1:ind(2)-1);	vars{k,3} = lixo{k}(ind(2)+1:end);
			else
				error('AUX_FUNS:GET_MGG', 'Wrong formating strig for MGD77 variabe selection')
			end
		end
	end
	if (any(c)),	names(c) = [];		names_ui(c) = [];	vars(c,:) = [];		end		% Remove eventual non-existing files

	if (isempty(names))					% Empty names list
		track = [];		x_min = [];		x_max = [];		y_min = [];		y_max = [];
		return
	end

	% EXTRACT NAVIGATION
	if (strcmpi(EXT, '.gmt'))		% Old style .gmt files (many of the above was useless)
		track = gmtlist_m(names{:}, varargin{:});
		for (k = 1:numel(names)),	names{k} = [names{k} EXT];		end		% Restore the extension
	else							% mgd77 netCDF files
		if (numel(names) == 1)		% Extrac the entire navigation
			names{1} = [names{1} EXT];
 			track.longitude = double(nc_funs('varget', names{1}, 'lon'));
 			track.latitude  = double(nc_funs('varget', names{1}, 'lat'));
			if (what2plot.mag),		track.magnetics  = nc_funs('varget', names{1}, 'mtf1');		end
			if (what2plot.grav),	track.gravity    = nc_funs('varget', names{1}, 'faa');		end
			if (what2plot.topo),	track.topography = nc_funs('varget', names{1}, 'depth');	end
			track.info = names{1};
		else						% Since the poor lousy HG gets stupidly slow, plot only one every other 5th point
			track(numel(names)).longitude = [];		% To shut up MLint
			track(numel(names)).latitude = [];
			for (k = 1:numel(names))
				names{k} = [names{k} EXT];
				x = double(nc_funs('varget', names{k}, 'lon'));		track(k).longitude = x(1:5:end);
				x = double(nc_funs('varget', names{k}, 'lat'));		track(k).latitude = x(1:5:end);
				if (what2plot.mag)
					x = nc_funs('varget', names{k}, 'mtf1');		track(k).magnetics  = x(1:5:end);
				end
				if (what2plot.grav)
					x = nc_funs('varget', names{k}, 'faa');			track(k).gravity  = x(1:5:end);
				end
				if (what2plot.topo)
					x = nc_funs('varget', names{k}, 'depth');		track(k).topography  = x(1:5:end);
				end
				track(k).info = names{k};
			end
		end
	end

	for (k = 1:length(track))
		ind = false( size((track(k).longitude)) );
		if (what2plot.mag && ~isempty(track(k).magnetics)),		ind = (~isnan(track(k).magnetics) | ind);	end
		if (what2plot.grav && ~isempty(track(k).gravity)),		ind = (~isnan(track(k).gravity) | ind);		end
		if (what2plot.topo && ~isempty(track(k).topography)),	ind = (~isnan(track(k).topography) | ind);	end

		if (~any(ind))						% When there isn't any data to plot
			warndlg('The selected field is empty in this file, but I will plot the track anyway.','Warning')
		else
			ind = ~ind;
			track(k).longitude(ind) = NaN;
			track(k).latitude(ind)  = NaN;
		end
	end

	if (nargout == 8)
		len_t = numel(track);  x_min = zeros(1,len_t); x_max = x_min;  y_min = x_min;  y_max = x_min;
		for (k = 1:len_t)
			x_min(k) = min(track(k).longitude);		x_max(k) = max(track(k).longitude);
			y_min(k) = min(track(k).latitude);		y_max(k) = max(track(k).latitude);
		end
		x_min = min(x_min);		x_max = max(x_max);
		y_min = min(y_min);		y_max = max(y_max);
		if (isnan(x_min))
			warndlg('This file had all requested records set to NaN. An error further down the road will likely occur','Warning')
		end
	end

% --------------------------------------------------------------------------------
function [data, agency] = mgd77info(fname)
% Get info from a MGD77+ netCDF file. The 'data' array is the way it is to be compatible with that of old .gmt files

	N_recs = numel(nc_funs('varget', fname, 'time'));
	x = nc_funs('varget', fname, 'mtf1');		x = ~isnan(x);		N_mag = numel(x(x));
	x = nc_funs('varget', fname, 'faa');		x = ~isnan(x);		N_grav = numel(x(x));
	x = nc_funs('varget', fname, 'depth');		x = ~isnan(x);		N_topo = numel(x(x));
	
	data = [N_recs N_grav N_mag N_topo];
	x = nc_funs('varget', fname, 'lon');		data(5:6) = double([min(x) max(x)]);
	x = nc_funs('varget', fname, 'lat');		data(7:8) = double([min(x) max(x)]);	clear x

	s = nc_funs('info',fname);
	ind = strcmp({s.Attribute.Name}, 'Survey_Departure_Year');		data(9) = str2double(s.Attribute(ind).Value);
	ind = strcmp({s.Attribute.Name}, 'Survey_Departure_Month');		data(10) = str2double(s.Attribute(ind).Value);
	ind = strcmp({s.Attribute.Name}, 'Survey_Departure_Day');		data(11) = str2double(s.Attribute(ind).Value);
	ind = strcmp({s.Attribute.Name}, 'Survey_Arrival_Year');		data(12) = str2double(s.Attribute(ind).Value);
	ind = strcmp({s.Attribute.Name}, 'Survey_Arrival_Month');		data(13) = str2double(s.Attribute(ind).Value);
	ind = strcmp({s.Attribute.Name}, 'Survey_Arrival_Day');			data(14) = str2double(s.Attribute(ind).Value);
	ind = strcmp({s.Attribute.Name}, 'Source_Institution');			agency = s.Attribute(ind).Value;

% --------------------------------------------------------------------------------
function out = seek_OPTcontrol(KEY)
% Check if there are rquirements about what to plot in the OPTcontrol.txt file.
% If not, default to 'Magnetics'

	out.mag = true;			% This is the default

	mir_dirs = getappdata(0,'MIRONE_DIRS');
	if (~isempty(mir_dirs))
		opt_file = [mir_dirs.home_dir '/data/OPTcontrol.txt'];
	else
		return				% Since we cannot find the OPTcontrol file
	end
	if (~(exist(opt_file, 'file') == 2)),		return,		end		% Nickles
	
	fid = fopen(opt_file, 'r');
	c = (fread(fid,'*char'))';      fclose(fid);
	lines = strread(c,'%s','delimiter','\n');   clear c fid;
	m = numel(lines);
	for (k = 1:m)
		if (~strncmp(lines{k}, KEY, numel(KEY))),	continue,	end
		% If we reach here that's becasuse the KEY line wsa found
		[t,opt] = strtok(lines{k});
		if (isempty(opt))
			opt = 'M';			% Default to good old Magnetics
		end
		opt = ddewhite(opt);
		out.grav = strcmp(opt,'G');
		out.mag  = strcmp(opt,'M');
		out.topo = strcmp(opt,'T');
		break
	end

% --------------------------------------------------------------------------------
function [latcells,loncells,Zcells] = localPolysplit(lat,lon, Z)
%POLYSPLIT Extract segments of NaN-delimited polygon vectors to cell arrays
%
%   [LATCELLS,LONCELLS] = POLYSPLIT(LAT,LON) returns the NaN-delimited
%   segments of the vectors LAT and LON as N-by-1 cell arrays with one
%   polygon segment per cell.  LAT and LON must be the same size and have
%   identically-placed NaNs.  The polygon segments are column vectors if
%   LAT and LON are column vectors, and row vectors otherwise.

% Copyright 1996-2006 The MathWorks, Inc.
% $Revision: 1.4.4.5 $    $Date: 2006/05/24 03:35:26 $

	n_arg = nargin;
	if (n_arg == 2)
		[lat, lon] = localRemoveExtraNanSeps(lat, lon);
	else
		[lat, lon, Z] = localRemoveExtraNanSeps(lat, lon, Z);
	end
	indx = find(isnan(lat(:)));         % Find NaN locations.

	% Simulate the trailing NaN if it's missing.
	if ~isempty(lat) && ~isnan(lat(end))
		indx(end+1,1) = numel(lat) + 1;
	end

	%  Extract each segment into pre-allocated N-by-1 cell arrays, where N is
	%  the number of polygon segments.  (Add a leading zero to the indx array
	%  to make indexing work for the first segment.)
	N = numel(indx);
	latcells = cell(N,1);       loncells = cell(N,1);
	if (n_arg == 3),			Zcells = cell(N,1);
	else						Zcells = [];
	end
	indx = [0; indx];
	for k = 1:N
		iStart = indx(k)   + 1;
		iEnd   = indx(k+1) - 1;
		latcells{k} = lat(iStart:iEnd);
		loncells{k} = lon(iStart:iEnd);
		if (n_arg == 3),	Zcells{k} = Z(iStart:iEnd);		end
	end

% --------------------------------------------------------------------------------
function [xdata, ydata, zdata] = localRemoveExtraNanSeps(xdata, ydata, zdata)
%removeExtraNanSeps  Clean up NaN separators in polygons and lines

	p = find(isnan(xdata(:)'));     % Determing the positions of each NaN.
	
	% Determine the position of each NaN that is not the final element in a sequence of contiguous NaNs.
	q = p(diff(p) == 1);
	
	% If there's a leading sequence of NaNs (a sequence starting with a NaN in
	% position 1), determine the position of each NaN in this sequence.
	if isempty(p),      r = [];
	else                r = find((p - (1:numel(p))) == 0);
	end
	
	% Determine the position of each excess NaN.
	if isempty(r),      s = q;
	else                s = [r q(q > r(end))];
	end
	
	% Remove the excess NaNs.
	xdata(s) = [];      ydata(s) = [];
	if (nargin >= 3),   zdata(s) = [];  end

% -*-*-*-*-*-*-$-$-$-$-$-$-#-#-#-#-#-#-%-%-%-%-%-%-@-@-@-@-@-@-(-)-(-)-(-)-&-&-&-&-&-&-{-}-{-}-{-}-
function appProjectionRef(handles, strWKT)
% If we have a WKT proj store it, otherwise clean eventual predecessors
	if (~isempty(strWKT))					% If we have a WKT projection, store it
		if (strWKT(1) == '+')				% Well, it is a proj4 string, convert it to WKT
			strWKT = ogrproj(strWKT);
		end
		setappdata(handles.figure1,'ProjWKT',strWKT)
		out = decodeProjectionRef(strWKT);				% Decode Proj reference string
		if ( ~isempty(out.datum) || ~isempty(out.ellipsoid) || ~isempty(out.projection) )
			setappdata(handles.axes1,'DatumProjInfo',out)
		end
	else							% Otherwise remove eventual previous one
		if (isappdata(handles.figure1,'ProjWKT')),    rmappdata(handles.figure1,'ProjWKT'),		end
	end

% -*-*-*-*-*-*-$-$-$-$-$-$-#-#-#-#-#-#-%-%-%-%-%-%-@-@-@-@-@-@-(-)-(-)-(-)-&-&-&-&-&-&-{-}-{-}-{-}-
function out = decodeProjectionRef(strProj)
	ind = strfind(strProj,char(10));
	out.datum = [];     out.ellipsoid = [];     out.projection = [];
	if (numel(ind) <= 1),   return;    end

	ind = [0 ind length(strProj)-1];
	for (i=1:numel(ind)-1)
		str = strProj(ind(i)+1:ind(i+1)-1);
		if ~isempty(strfind(str,'GEOGCS'))      % Get datum
			xx = strfind(str,'"');
			if (numel(xx) < 2),     continue;   end
			out.datum = str(xx(1)+1:xx(2)-1);
		end
		if ~isempty(strfind(str,'SPHEROID'))      % Get ellipsoid
			if (numel(xx) < 2),     continue;   end
			xx = strfind(str,'"');
			out.ellipsoid = str(xx(1)+1:xx(2)-1);
		end
		if ~isempty(strfind(str,'PROJECTION'))      % Get ellipsoid
			if (numel(xx) < 2),     continue;   end
			xx = strfind(str,'"');
			out.projection = str(xx(1)+1:xx(2)-1);
		end
	end

% ----------------------------------------------------------------
function [img, pal] = semaforo_green()

ind = cell(16,1);

ind{1} = [53	82	96	97	98	99	104	105	111	112	126	141	142	143	144	152	156	157	158	170	171	172	187	188	189	201	202	203	216	217	232	233	234	235	245	246	247	248	249	262	278	279	280	281	291	292	293	294	308	324	325	326	327	334	336	337	338	339	340	354	370	371	372	373	379	381	382	384	385	386	387	397	400	401	416	417	418	419	420	424	425	426	428	429	431	432	433	442	446	447	463	464	465	466	468	470	472	473	474	478	479	480	487	488	492	493	494	512	516	517	518	519	525	530	532	541	546];
ind{2} = [129	130	131	132	173	179	226	264	310	356	402	449	496	501];
ind{3} = [70	114	124	211	255	260	299	301	306	344	345	388	389	430	445	477	489	491	520	538	562	563	576	577	583	584];
ind{4} = [77	80	86	115	116	117	123	134	160	161	162	181	205	206	207	208	251	252	297	320	366	503	523	548	565	569	571	578	592	595];
ind{5} = [103	145	197	198	199	236	241	242	243	244	282	285	286	288	289	290	328	330	332	333	335	376	377	378	380	383	421	422	423	427	469	471];
ind{6} = [51	55	56	64	66	69	81	84	91	100	118	137	153	165	183	204	212	213	229	256	257	258	261	300	302	303	304	321	346	348	349	352	390	391	392	393	413	435	436	437	457	475	482	539	547	551	570	572	575	587	589	590];
ind{7} = [5	6	7	8	9	10	11	12	13	14	15	16	17	19	20	21	22	23	24	25	26	27	28	29	30	33	34	35	36	37	38	39	40	41	42	43	44];
ind{8} = [18	32	45	57	65	72	120	155	169	200	557	564	603	604	605	606	607	608	609	610	611	612	613	614	615	618	619	620	621	622	623	624	625	626	627	628	630	631	632	633	634	635	637	639	640	641	642];
ind{9} = [174	175	176	177	178	219	220	221	222	223	224	225	265	266	267	268	269	270	271	272	311	312	313	314	315	316	317	318	357	358	359	360	361	362	363	364	403	404	405	406	407	408	409	410	450	451	452	453	454	455	497	498	499	500];
ind{10} = [52	54	67	68	83	110	113	119	127	133	151	166	186	190	259	263	275	295	305	307	309	341	350	351	353	355	367	374	394	395	396	411	434	439	440	441	443	448	459	467	481	483	484	485	486	502	505	510	511	513	514	515	524	526	527	528	529	531	533	540	542	543	558	559	560	574	588	597];
ind{11} = [59	60	61	62	63	73	74	75	76	79	87	88	89	90	107	108	109	121	122	135	136	154	168	182	214	228	274	412	444	458	476	490	504	521	522	535	536	549	550	566	567	568	579	580	581	582	585	593	594	596];
ind{12} = [	92	106	125	138	184	227	230	276	277	322	323	368	369	399	414	460	506	552	561	598	616	629	643];
ind{13} = [	101	102	146	147	148	149	150	191	192	193	194	195	196	237	238	239	240	283	284	287	329	331	375];
ind{14} = [58	71	78	85	159	163	164	167	209	210	215	250	253	254	296	298	342	343	398	534	537	586	591];
ind{15} = [1	2	3	4	31	46	47	48	49	50	93	94	95	139	140	185	231	415	461	462	507	508	509	553	554	555	556	599	600	601	602	617	636	638	644];
ind{16} = [128	180	218	273	319	347	365	438	456	495	544	545	573];

pal = [ ...
		9	10	0
		1	140	0
		98	94	17
		149	137	23
		76	2	0
		65	65	10
		150	146	103
		101	98	62
		0	247	0
		24	35	4
		192	163	53
		58	65	56
		132	2	0
		119	115	13
		99	104	100
		4	75	0];

	img = uint8(zeros(46*14,1));
	for (k = 0:15)
		img(ind{k+1}) = k;
	end
	img = reshape(img,46,14);
	pal = pal / 255;

% --------------------------------------------------------------------------------------
function [img, pal] = semaforo_red()
ind = cell(16,1);

ind{1} = [52	53	68	96	97	98	110	111	112	113	126	141	142	143	156	157	158	171	187	188	201	202	203	216	217	232	233	234	246	247	248	249	262	278	279	280	292	293	294	295	308	324	325	326	338	339	340	341	354	370	371	372	384	385	386	387	397	400	416	417	418	431	432	433	442	446	447	463	464	465	478	479	480	487	488	492	493	494	511	512	525	526	527	531	532	540	541];
ind{2} = [99	144	152	199	235	245	373	383	466	474	514	518];
ind{3} = [78	159	163	164	169	208	209	210	250	253	254	300	342	344	389	584];
ind{4} = [130	131	174	175	176	177	178	179	219	220	221	222	223	224	225	226	265	266	267	268	269	270	271	272	310	311	312	313	314	315	316	317	318	319	356	357	358	359	360	361	362	363	364	365	403	404	405	406	407	408	409	410	449	450	451	452	453	454	455	456	497	498	499	500	501];
ind{5} = [115	116	117	160	161	162	205	206	207	251	252	296	297];
ind{6} = [64	69	91	118	165	204	211	212	213	255	256	257	258	261	299	301	302	303	306	345	346	347	352	388	390	391	399	435	436	538	539	561	572	575	576	587	589	590	597];
ind{7} = [5	6	7	8	9	10	11	12	13	14	15	16	17	19	20	21	22	23	24	25	26	27	28	29	30	33	34	35	36	37	38	39	40	41	42	43	44	50	94];
ind{8} = [1	2	3	4	31	32	46	47	48	49	93	95	123	139	140	185	231	415	461	462	507	508	509	544	553	554	555	556	599	600	601	602	603	604	605	607	609	610	611	613	617	618	619	620	621	622	623	625	626	628	631	633	634	635	637	639	640	641	644];
ind{9} = [54	67	82	119	127	153	166	170	172	259	275	304	305	307	348	349	350	351	353	392	393	395	396	401	413	419	434	438	439	440	441	443	475	481	482	483	484	485	486	513	524	528	529	530	533	551	559	570	573	574	588];
ind{10} = [100	101	102	103	104	145	146	147	148	149	150	151	190	191	192	193	194	195	196	197	198	236	237	238	239	240	241	242	243	244	281	282	283	284	285	286	287	288	289	290	291	327	328	329	330	331	332	333	334	335	336	337	374	375	376	377	378	379	380	381	382	420	421	422	423	424	425	426	427	428	467	468	469	470	471	472	473	515	516	517];
ind{11} = [59	60	61	62	63	73	74	75	76	79	87	88	89	90	107	108	109	121	122	135	136	154	168	182	214	228	444	458	476	490	522	535	536	550	566	567	568	580	582	594];
ind{12} = [18	45	51	92	106	129	132	138	167	180	184	186	200	230	264	273	276	277	322	323	368	369	402	411	414	460	496	506	510	520	543	545	546	552	557	560	562	598	606	614	616	629	630	632	643];
ind{13} = [56	57	430	495	547];
ind{14} = [55	58	70	71	72	84	114	120	124	128	133	155	173	215	218	227	260	298	343	398	445	448	457	477	489	491	502	534	537	563	564	577	583	591	608	612	615	624	627	636	638	642];
ind{15} = [66	81	83	105	125	137	183	189	229	263	309	321	355	367	394	429	437	459	505	519	542	558];
ind{16} = [65	77	80	85	86	134	181	274	320	366	412	503	504	521	523	548	549	565	569	571	578	579	581	585	586	592	593	595	596];

pal = [ ...
		5	14	1
		148	4	0
		117	116	1
		1	101	50
		148	148	0
		75	74	11
		150	146	103
		102	101	78
		32	37	5
		251	0	0
		197	168	55
		67	76	57
		152	73	26
		112	90	28
		60	41	7
		169	135	43];

	img = uint8(zeros(46*14,1));
	for (k = 0:15)
		img(ind{k+1}) = k;
	end
	img = reshape(img,46,14);
	pal = pal / 255;

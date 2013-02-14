function varargout = load_xyz(handles, opt, opt2)
% Read a generic ascii file that can be, or not, a multiseg file
%
%	Multi-segs files accept -G, -S & -W GMT type options.
%		-S<symb>[size] accepts these GMT type symbol codes <a|c|d|h|i|n|p|s|x|+>
%
%	It does also deal with the case of ploting the isochrons.dat
%
%	HANDLES	->	Should be the Mirone handles. However, when this function is used with output
%				HANDLES can be empty([]) or just a simple structure with the field 'no_file' = false
%
%	Optional
%		OPT can be either [] in which case the fiename will be asked here or contain the filename
%		OPT2 can take several values 'arrows' to read a x,y,u,v file
%			'AsLine'		plots a "regular" line
%			'AsPoint'		plots the line vertex only using small filled circles 
%			'AsMaregraph'	plots yellow dots used in the Tsunami modeling tools
%			'FaultTrace'	plots lines/polylines used by the elastic deformation tools
%			'Isochron'		plots a isochrons polyline from the internal db.
%			'FZ'			plots a Fracture Zones polyline from the internal db.
%							Attention, this option needs a non-empty OPT argument
%			'ncshape'		Input file is a netCDF ncshape
%			If not given defaults to 'AsLine'
%
% If first line in file is of the form
%		'>U_N_I_K'	plot a single line NaN separated
%		'>ARROW'	plot an arrow field
%		'>VIMAGE'	tell Fleder to plot a scene with a VIMAGE
%		'>-:'		swap 1st and 2nd columns (assumed as (y,x) -> (x,y))
%		'>CLOSE'	plot patches instead of lines (idependently of pline being closed or not)
%		'>HAVE_INCLUDES'	Signal that this file has multi-segment headers with the form:
%					'> INCLUDE=FULL_PATH_TO_FILE'
%					that will result in inclusion of an external file at the end of the main one.
%					These included files may be multi-segment files as well but need to have the
%					same number of columns as the master file (no test for this though).
%					When used from an empty figure these included files won't be taken into account
%					to determine data BB because the parsing of this option is done afterwards.
%					If the indicated file(s) do not exist this option is silently ignored.
%		'>NESTING'	File contains 'Nesting grids' rectangles (multisegment) to help with nested grids creation 
%					These files have >NESTING DX DY REG in first multisegment header and only the
%					> DX DY REG for the interior rectangles.
%					REG is either 0 (grid registration) or 1 (pixel reg). Default is 0
%					If neither DY and REG is provided, DY = DX and REG = 0 are assumed

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

%	EXAMPLE CODE OF HOW TO CREATE A TEMPLATE FOR UICTX WHEN THESE ARE TOO MANY
% 	cmenuHand = get(h, 'UIContextMenu');
% 	setappdata(handles.figure1, 'cmenuHand', cmenuHand)
% 	cmenuHand = uicontextmenu('Parent',handles.figure1);
% 	set(h, 'UIContextMenu', cmenuHand);
% 	%uimenu(cmenuHand, 'Label', 'Set all UIcontexts', 'Call', {@resetUIctx,h,handles.axes1});
% 	uimenu(cmenuHand, 'Label', 'Set all UIcontexts', 'Call', 'hand=guidata(gco); set(gco, ''UIContextMenu'', getappdata(hand.axes1, ''cmenuHand''))' );
% 
% function resetUIctx(obj,evt,h,hAxes)
% 	cmenuHand = getappdata(hAxes, 'cmenuHand');
% 	set(h, 'UIContextMenu', cmenuHand)

	% ------------------- Some defaults ---------------------------------------------
	tol = 0.5;
	do_project = false;         % We'll estimate below if this holds true
	got_arrow = false;
	got_internal_file = false;
	got_nc = false;				% Flag to signal if a shapenc type file comes in
	orig_no_mseg = false;		% Flag to know if original file had multiseg strings to store
	line_type = 'AsLine';
	tag = 'polyline';
	struc_vimage = [];
	is_bin = false;				% To flag binary files
	BB = [];					% To eventually hold a BoundingBox
	is_GMT_DB = false;			% To flag the special cases when we are dealing with the GMT database polygons
	do_nesting = false;			% To flag when the imported file is 'Nesting squares'
	% -------------------------------------------------------------------------------

	% ------------------- PARSE INPUTS ----------------------------------------------
	if (nargout && isempty(handles))		% When this function is used to just read a file and return its contents
		handles.no_file = false;
	end
	if (nargin >= 2 && isempty(opt))            % Read a ascii file
		[FileName, PathName, handles] = put_or_get_file(handles, ...
			{'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select File','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
		[lix,lix,EXT] = fileparts(FileName);
	elseif (nargin >= 2)		% Read a ascii file of which we already know the name (drag N'drop)
		fname = opt;
		[PathName,lix,EXT] = fileparts(fname);			% We need the 'PathName' below
	end
	if (nargin == 3)
		if (strcmp(opt2, 'AsLine') && strcmpi(EXT, '.nc'))			% A shapenc loaded 'AsLine', but confirm
			drv = aux_funs('findFileType',fname);
			if (strcmp(drv, 'ncshape')),	opt2 = 'ncshape';	end
		end
		if (strcmp(opt2, 'AsArrow'))		got_arrow = true;		% This case does not care about 'line_type'
		elseif (strcmp(opt2, 'ncshape'))	got_nc = true;			%
		else								line_type = opt2;
		end
		if (strncmpi(line_type,'isochron',4) || strcmpi(line_type,'FZ'))
			if (strncmpi(line_type,'isochron',4))
				tag = 'isochron';		fname = [handles.path_data 'isochrons.dat'];
			else
				tag = 'FZ';				fname = [handles.path_data 'fracture_zones.dat'];
			end
			got_internal_file = true;	PathName = handles.path_data;
			line_type = 'i_file';
		end
	end
	% ------------------- END PARSE INPUTS------------------------------------------------

	if (~got_nc)			% Most common cases

		[bin, n_column, multi_seg, n_headers] = guess_file(fname);
		if (isempty(bin))
			errordlg(['Error reading file (probably empty)' fname],'Error'),	return
		end
		if (isa(bin,'struct') || bin ~= 0)				% ---****** BINARY FILE *******---
			if (isa(bin,'struct'))
				bin = guess_bin(bin.nCols, bin.type);	% Ask user to confirm/modify guessing
			else
				bin = guess_bin(false);					% Ask user what's in file
			end
			if (isempty(bin))		% User quit
				varargout = {};		return
			end
			n_column = bin.nCols;
			multi_seg = 0;		n_headers = 0;		is_bin = true;
			if (nargin < 3),	line_type = 'AsPoint';	end		% If not explicitly set line_type
		end

	else
		out_nc = read_shapenc(fname);
		n_column = 2;		% ???
		bin = false;		multi_seg = 0;		n_headers = 0;
		if ( (out_nc.n_PolyOUT + out_nc.n_PolyIN) == 0 )
			warndlg('Warning, no polygons to plot in this shapenc file. Bye','Warning')
			if (nargout),	[varargout{1:nargout}] = {};		end
			return
		end
		BB = out_nc.BB;
		numeric_data = cell(out_nc.n_PolyOUT + out_nc.n_PolyIN, 1);
		kk = 1;
		for (k = 1:out_nc.n_PolyOUT)		% Order will be for each swarm OUT + its INs
			numeric_data{kk} = [out_nc.Poly(k).OUT.lon out_nc.Poly(k).OUT.lat];
			for ( n = 1:numel(out_nc.Poly(k).IN) )
				if (isempty(out_nc.Poly(k).IN(n).lon)),		continue,	end
				kk = kk + 1;
				numeric_data{kk} = [out_nc.Poly(k).IN(n).lon out_nc.Poly(k).IN(n).lat];
			end
			kk = kk + 1;
		end
		
		if (~isempty(out_nc.description))
			desc = ['> ' out_nc.description];
		else
			desc = '> Nikles ';
		end
		if (numel(numeric_data) > 1)
			multi_seg = true;
			multi_segs_str = repmat({desc}, numel(numeric_data),1);
		else
			multi_segs_str = {desc};
		end

		if (numel(out_nc.Poly) == numel(multi_segs_str))% UGLY PATCH FOR WHEN INNER POLYGONS WOULD SCREW ALL
			for (k = 1:numel(multi_segs_str))			% Append the individual description to the Dataset's one
				multi_segs_str{k} = sprintf('%s\n\n    %s', multi_segs_str{k}, out_nc.Poly(k).OUT.desc);
			end
		end
		out_nc.fname = fname;		% Store the names as we will need it in the set_extra_uicb_options()

		% Add the referencing info to multi_segs_str{1} so that later the recogn algo is the same as for ascii files
		if ( ~isempty(out_nc.SRS) && out_nc.SRS(1) == '+' )		% Only proj4 are accepted
			multi_segs_str{1} = [multi_segs_str{1} ' ' out_nc.SRS];
		end
	end

	if (n_column == 1 && multi_seg == 0)			% Take it as a file names list
		fid = fopen(fname);
		c = fread(fid,'*char')';	fclose(fid);
		names = strread(c,'%s','delimiter','\n');   clear c fid;
	else
		names = {fname};
	end

	% ------------------ Section to deal with START NEW FIG or PLOT ON EXISTING FIG ----------
	if (handles.no_file)		% Start empty but below we'll find the true data region
		if (ischar(handles.DefLineColor) && handles.DefLineColor(1) == 'w')
			handles.DefLineColor = 'k';		% To not plot a white line over a white background
		end
		XMin = 1e50;			XMax = -1e50;    YMin = 1e50;            YMax = -1e50;

		for (k = 1:numel(names))
			fname = names{k};
			j = strfind(fname,filesep);
			if (isempty(j)),    fname = sprintf('%s%s',PathName, fname);   end		% Need to add path as well 
			if (isempty(n_headers)),    n_headers = NaN;    end
			if (~got_nc)				% Otherwise data was read already
				if (multi_seg)
					[numeric_data, multi_segs_str] = text_read(fname,NaN,n_headers,'>');
				elseif (~is_bin)
					numeric_data = text_read(fname,NaN,n_headers);
				else					% Try luck with a binary file
					fid = fopen(fname);		numeric_data = fread(fid,['*' bin.type]);	fclose(fid);
					numeric_data = reshape(numeric_data,bin.nCols,numel(numeric_data)/bin.nCols)';
					if (bin.twoD && (bin.nCols > 2)),	numeric_data(:,3:end) = [];		end % Retain only X,Y 
				end
			end

			if (~isa(numeric_data,'cell'))			% File was not multi-segment.
				numeric_data = {numeric_data};
				multi_segs_str = {'> Nikles '};		% Need something (>= 8 chars) to not error further down
				orig_no_mseg = true;
			end
			if (isempty(BB))
				for (i = 1:length(numeric_data))
					XMin = min(XMin,double(min(numeric_data{i}(:,1))));		XMax = max(XMax,double(max(numeric_data{i}(:,1))));
					YMin = min(YMin,double(min(numeric_data{i}(:,2))));		YMax = max(YMax,double(max(numeric_data{i}(:,2))));
				end
			else
				XMin = BB(1);		XMax = BB(2);		YMin = BB(3);		YMax = BB(4);
			end
		end

		if ( multi_seg && strncmp(multi_segs_str{1},'>-:',3) )		% See if we need to swap x<->y
			tmp = XMin;		XMin = YMin;	YMin = tmp;				% Need to swap min/max
			tmp = XMax;		XMax = YMax;	YMax = tmp;
		end

		dx = XMax - XMin;			dy = YMax - YMin;
		if (dx == 0 || dy == 0)
			errordlg('File is has only one point or all XXs are equal or all YYs are equal','Error')
			return
		end
		
		if (XMin > -179.5 && XMax < 359.5 && YMin > -89.5 && YMax < 89.5)
			XMin = XMin - dx / 200;		XMax = XMax + dx / 200;		% Give an extra 0.5% padding margin
			YMin = YMin - dy / 200;		YMax = YMax + dy / 200;
		end
		xx = [XMin XMax];			yy = [YMin YMax];
		region = [xx yy];
		handles.geog = aux_funs('guessGeog',region);

		if (got_internal_file)				% We know it's geog (Global Isochrons or FZs)
			xx = [-180 180];		yy = [-90 90];
			if (~handles.geog)				handles.geog = 1;
			elseif (handles.geog == 2)		xx = [0 360];
			end
			region = [xx yy];
		end
		mirone('FileNewBgFrame_CB', handles, [region handles.geog])	% Create a background
		hMirFig = handles.figure1;
		drawnow						% Otherwise it takes much longer to plot and other shits
	elseif (~nargout)				% Reading over an established region
		XYlim = getappdata(handles.axes1,'ThisImageLims');
		xx = XYlim(1:2);			yy = XYlim(3:4);
		if (handles.is_projected && (got_internal_file || handles.defCoordsIn > 0) )
			do_project = true;
		end
		XMin = XYlim(1);			XMax = XYlim(2);		% In case we need this names below for line trimming
	end
	% --------------------- End of CREATE NEW FIG / PLOT on EXISTING ONE -----------------

	% --------------------------- Main loop over data files ------------------------------
	for (k = 1:numel(names))
		fname = names{k};
		if (handles.no_file && k == 1)			% Rename figure with dragged file name
			[pato,barName] = fileparts(fname);
			old_name = get(hMirFig,'Name');		ind = strfind(old_name, '@');
			set(hMirFig,'Name',[barName old_name(ind-1:end)])
		end

		if (~handles.no_file)					% Otherwise we already read it
			j = strfind(fname,filesep);
			if (isempty(j)),    fname = sprintf('%s%s',PathName, fname);   end		% Need to add path as well 
			if (isempty(n_headers)),    n_headers = NaN;    end
			if (~got_nc)			% Otherwise data was read already
				if (multi_seg)
					[numeric_data, multi_segs_str] = text_read(fname,NaN,n_headers,'>');
				elseif (~is_bin)
					numeric_data = text_read(fname,NaN,n_headers);
				else				% Try luck with a binary file
					fid = fopen(fname);		numeric_data = fread(fid,['*' bin.type]);		fclose(fid);
					numeric_data = reshape(numeric_data,bin.nCols,numel(numeric_data)/bin.nCols)';
					if (bin.twoD && (bin.nCols > 2)),	numeric_data(:,3:end) = [];		end % Retain only X,Y
				end
			end
		end

		if (~isa(numeric_data,'cell'))			% File was not multi-segment. Now pretend it was but with no info
			numeric_data = {numeric_data};
			multi_segs_str = {'> Nikles '};		% Need something in it to not error below
			orig_no_mseg = true;
		end
		n_isoc = 0;     n_segments = length(numeric_data);
		hLine = ones(n_segments,1)*NaN;			% This is the maximum we can have
		n_clear = false(n_segments,1);
		do_patch = false;						% Default to line object

		% Test if conversion into a single, NaN separated, line its wanted 
		if (strncmp(multi_segs_str{1}, '>U_N_I_K', 8))
			for (i = 1:n_segments-1)
				numeric_data{i} = [numeric_data{i}(:,1:2); nan nan];
			end
			numeric_data{1} = cat(1,numeric_data{:});
			% Rip the U_N_I_K identifier
			if (numel(multi_segs_str{1}) > 8)					% We may have line type specifications
				multi_segs_str{1} = ['> ' multi_segs_str{1}(9:end)];
			else
				multi_segs_str{1} = '> ';
			end
			n_segments = 1;				% Pretend we have only one segment

		elseif (strncmpi(multi_segs_str{1}, '>ARROW', 6) || got_arrow)		% ARROW field (the got_arrow can came via varargin)
			if (~got_arrow),	multi_segs_str{1}(2:6) = [];	end			% Rip the ARROW identifier
			got_arrow = true;
			if (n_column < 4)
				errordlg('Files for arrow plot need 4 columns with the traditional (x,y,u,v)','ERROR'),	return
			end
			UV = cell(n_segments,1);
			for (i = 1:n_segments)		% Split the XY & UV columns to be compatible with the other options
				UV{i} = numeric_data{i}(:,3:4);
				numeric_data{i}(:,3:end) = [];
			end
			struc_arrow = struct('spacingChanged',[], 'hQuiver', [], 'hAx', handles.axes1);

		elseif (strncmp(multi_segs_str{1}, '>VIMAGE', 7))
			[z_Vmin, r] = strtok(multi_segs_str{k}(8:end));		z_Vmin = str2double(z_Vmin);
			[z_Vmax, r] = strtok(r);							z_Vmax = str2double(z_Vmax);
			vimage = ddewhite(r);
			if (isnan(z_Vmin) || isnan(z_Vmax))
				errordlg('Load VIMAGE error. First 2 fields must contain Z_START & Z_END info.','Error'),	return
			end
			if (~ischar(vimage) || ~exist(vimage,'file'))
				errordlg('Load VIMAGE error. Third field must contain an existing picture file name.','Error'),	return
			end
			struc_vimage = struct('z_min', z_Vmin, 'z_max', z_Vmax, 'vimage', vimage);

		elseif (strncmp(multi_segs_str{1}, '>-:', 3))				% File has y,x instead of x,y
			multi_segs_str{1}(2:3) = [];							% Rip the swap -: identifier
			for (i = 1:n_segments)				% Swapp 1st and 2th columns. Do differently would be very complex
				tmp = numeric_data{i}(:,1);
				numeric_data{i}(:,1) = numeric_data{i}(:,2);
				numeric_data{i}(:,2) = tmp;
			end
			clear tmp

		elseif (strncmp(multi_segs_str{1}, '>CLOSE', 6))			% Closed or not, plot a patch
			multi_segs_str{1}(2:6) = [];							% Rip the CLOSE identifier
			do_patch = true;

		elseif (strncmpi(multi_segs_str{1}, '>HAVE_INCLUDES', 7))	% This file includs other files
			if (numel(multi_segs_str) - numel(numeric_data) >= 2)	% Test if first segment has a true header
				multi_segs_str(1) = [];								% (Yes, it has). This header was now in excess.
			else
				multi_segs_str{1} = '> Nikles ';					% We need a first header that didn't exist
			end
			heads_to_del = false(numel(multi_segs_str),1);
			for (i = 1:numel(multi_segs_str))
				ind = strfind(multi_segs_str{i},'INCLUDE');
				if (~isempty(ind))
					tok = multi_segs_str{i}(ind(1):end);
					inc_fname = tok(min(numel(tok),9):end);			% Securely try to get the include file name
					heads_to_del(i) = true;							% This is not a true header, so flag it to deletion
					if (~exist(inc_fname,'file')),		continue,	end	% File does not exist or bad file name
					[out_data, out_str] = load_xyz([], inc_fname);
					imp_n_segments = 1;
					if (iscell(out_data)),	imp_n_segments = numel(out_data);	end
					for (imp_i = 1:imp_n_segments)					% Loop over number of segments of this imported file					
						if (iscell(out_data)),	numeric_data{end+1} = out_data{imp_i};
						else					numeric_data{end+1} = out_data;
						end
						n_segments = n_segments + 1;
						multi_segs_str{end+1} = out_str{imp_i};% Import also eventual included headers
					end
				end
			end
			multi_segs_str(heads_to_del) = [];

		elseif (strncmp(multi_segs_str{1}, '>NESTING', 5))			% 
			multi_segs_str{1}(2:8) = [];							% Rip the swap NESTING identifier
			do_nesting = true;
			orig_no_mseg = true;		% SHIT, I should not have to do this, but need it to go to 'line_uicontext'

		elseif (line_type(3) ~= 'P' && ~isempty(strfind(multi_segs_str{1},'-G')) && isempty(strfind(multi_segs_str{1},'-S')) )
			% -G (paint) alone is enough to make it a patch (if ~point)
			do_patch = true;
		end

		% ---------- If OUT is requested there is nothing left to be done here --------------
		if (nargout)
			if (orig_no_mseg),		numeric_data = numeric_data{1};		end
			varargout{1} = numeric_data;
			if (nargout == 2),	varargout{2} = multi_segs_str;		end
			return
		end
		% -----------------------------------------------------------------------------------

		% ------------------ Check if it is a GSHHS or WDBII file ---------------------------
		if (numel(multi_segs_str{1}) > 120 && strcmp(multi_segs_str{1}(2:7),' Id = ') && ...
				numel(strfind(multi_segs_str{1},'=')) > 10)
			tag = 'GMT_DBpolyline';		is_GMT_DB = true;
		end
		% -----------------------------------------------------------------------------------

		% ------------------ Check if first header line has a PROJ4 string ------------------
		[projStr, multi_segs_str{1}] = parseProj(multi_segs_str{1});
		if (~isempty(projStr)),		do_project = true;		end
		% -----------------------------------------------------------------------------------

		drawnow
		for (i = 1:n_segments)		% Loop over number of segments of current file (external loop)
			tmpz = [];
			if (do_project)         % We need to project
				try
					if (~isempty(projStr))
						numeric_data{i} = proj2proj_pts(handles, numeric_data{i}, 'srcProj4', projStr);

						if (k == 1 && i == 1)	% First time. Store proj info if Figure's appdata
							if ( isempty(getappdata(handles.figure1,'ProjWKT')) && ...
								 isempty(getappdata(handles.figure1,'Proj4')) && ...
								 isempty(getappdata(handles.figure1,'ProjGMT')) )
								setappdata(handles.figure1,'Proj4', projStr)
								handles.is_projected = true;
							end
						end
					else
						numeric_data{i} = proj2proj_pts(handles, numeric_data{i});
					end
					if (any( isinf(numeric_data{1}(1:min(20,size(numeric_data{1},1)))) ))
						warndlg('Your data was probably already projected. Right-click on the axes frame and uncheck the ''Load files in Geogs'' ','Warning')
					end
				catch
					errordlg(lasterr,'ERROR');    return
				end
			end

			indx = false;	indy = false;			% Default to no need for map clipping
			if (handles.no_file)
				tmpx = numeric_data{i}(:,1);		tmpy = numeric_data{i}(:,2);
			else
	 			difes = [( double(numeric_data{i}(1,1)) - double(numeric_data{i}(end,1)) ) ...	% Remember R13
					( double(numeric_data{i}(1,2)) - double(numeric_data{i}(end,2)) )];
				if (any(abs(difes) > 1e-5))			% Assume a not closed polygon
					[tmpx, tmpy, indx, indy] = ...	% Get rid of points that are outside the map limits
						aux_funs('in_map_region',handles,numeric_data{i}(:,1),numeric_data{i}(:,2),tol,[xx yy]);
				elseif (~is_GMT_DB)
					% TEMPORARY. For now if the polygon is partially inside we plot it all (no clipping)
					[tmpx, tmpy, indx, indy] = ...
						aux_funs('in_map_region',handles,numeric_data{i}(:,1),numeric_data{i}(:,2),-1,[xx yy]);
				else
					% SPECIAL case of GMT database polygons. Lon my be wrapped around the sphere and header string updated
					[tmpx, tmpy, indx, indy, multi_segs_str{i}] = ...
						aux_funs('in_map_region',handles,numeric_data{i}(:,1),numeric_data{i}(:,2),-1,[xx yy],multi_segs_str{i});
				end
			end
			if (isempty(tmpx)),     n_clear(i) = true;     continue,		end     % Store indexes for clearing vanished segments info
			if ( numel(numeric_data{i}(1,:)) >= 3 )		% If we have a Z column
				tmpz = numeric_data{i}(:,3);
				if (~isempty(indx) || ~isempty(indy))	tmpz(indx) = [];	tmpz(indy) = [];	end	% If needed, clip outside map data			
			end

			[lThick, cor, multi_segs_str{i}] = parseW(multi_segs_str{i}(min(2,numel(multi_segs_str{i})):end)); % First time, we can chop the '>' char
			if (isempty(lThick)),	lThick = handles.DefLineThick;	end		% IF not provided, use default
			if (isempty(cor)),		cor = handles.DefLineColor;		end		%           "

			% ---------- Check if we have a symbols request. If yes turn the 'AsPoint' option on --------
			[marker, markerSize, multi_segs_str{i}] = parseS(multi_segs_str{i});
			if (~isempty(marker))
				line_type = 'AsPoint';
			else
				marker = 'o';	markerSize = 2;		% The old defaults
			end

			if (do_patch)
				Fcor = parseG(multi_segs_str{i});
				if (isempty(Fcor)),		Fcor = 'none';		end
				hPat = patch('XData',tmpx,'YData',tmpy,'Parent',handles.axes1,'Linewidth',lThick,'EdgeColor',cor,'FaceColor',Fcor);
				if (~isempty(tmpz))
					set(hPat,'UserData',tmpz');
				end	
				draw_funs(hPat,'line_uicontext')
				n_clear(i) = true;			% Must delete this header info because it only applyies to lines, not patches

			else					% Line plottings
				% See if we need to wrap arround the earth roundness discontinuity. Using 0.5 degrees from border. 
				if (handles.geog == 1 && ~do_project && (XMin < -179.5 || XMax > 179.5) )
					[tmpy, tmpx, tmpz] = map_funs('trimwrap', tmpy, tmpx, [-90 90], [XMin XMax], tmpz, 'wrap');
				elseif (handles.geog == 2 && ~do_project && (XMin < 0.5 || XMax > 359.5) )
					[tmpy, tmpx, tmpz] = map_funs('trimwrap', tmpy, tmpx, [-90 90], [XMin XMax], tmpz, 'wrap');
				end

				n_isoc = n_isoc + 1;
				if (~got_arrow)
					switch line_type
						case {'AsLine' 'i_file'}		% 'i_file' means internal file (Isochrons or FZs)
							hLine(i) = line('XData',tmpx,'YData',tmpy,'Parent',handles.axes1,'Linewidth',lThick,...
									'Color',cor,'Tag',tag,'Userdata',n_isoc);
							setappdata(hLine(i),'LineInfo',multi_segs_str{i});
						case 'AsPoint'
							Fcor = parseG(multi_segs_str{i});			% See if user wants colored pts
							if (isempty(Fcor)),		Fcor = 'k';		end
							hLine(i) = line('XData',tmpx,'YData',tmpy,'Parent',handles.axes1, 'LineStyle','none', 'Marker',marker,...
								'MarkerEdgeColor','k','MarkerFaceColor',Fcor, 'MarkerSize',markerSize,'Tag','Pointpolyline');
							draw_funs(hLine(i),'DrawSymbol')			% Set marker's uicontextmenu (tag is very important)
						case 'AsMaregraph'
							hLine(i) = line('XData',tmpx,'YData',tmpy,'Parent',handles.axes1, 'LineStyle','none', 'Marker','o',...
								'MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',10,'Tag','Maregraph');
							draw_funs(hLine(i),'DrawSymbol')			% Set marker's uicontextmenu					
						case 'FaultTrace'
							hLine(i) = line('XData',tmpx,'YData',tmpy,'Parent',handles.axes1,'Color',cor,'LineWidth',lThick,'Tag','FaultTrace');
							draw_funs(hLine(i),'line_uicontext')		% Set lines's uicontextmenu
							% Create empty patches that will contain the surface projection of the fault plane
							hp = zeros(1, numel(tmpx)-1);
							for (j = 1:numel(tmpx)-1),	hp(j) = patch('XData', [], 'YData',[]);    end
							setappdata(hLine(i),'PatchHand',hp);
					end
				else
					struc_arrow.color = cor;
					hQuiver = draw_funs([], 'loc_quiver', struc_arrow, tmpx, tmpy, UV{i}(:,1), UV{i}(:,2));
					set(hQuiver,'Tag','Seta','Userdata',n_isoc)
					setappdata(hQuiver(1),'MyHead',hQuiver(2))		% Store the arrows heads handle
					hLine(i) = hQuiver(1);
				end

				if (~isempty(tmpz))		setappdata(hLine(i),'ZData',tmpz');	end
				if (~orig_no_mseg)
					setappdata(hLine(i),'LineInfo',multi_segs_str{i})  % To work with the sessions and will likely replace old mechansim
				end

				% Finish the Vertical image section (if it exists obviously)
				if (~isempty(struc_vimage))
					vimage = getappdata(handles.axes1,'VIMAGE');
					if (isempty(vimage))			% First one
						struc_vimage.hLine = hLine(i);
						setappdata(handles.axes1, 'VIMAGE', struc_vimage)
					else
						struc_vimage.hLine = hLine(i);
						vimage(end+1) = struc_vimage;
						setappdata(handles.axes1, 'VIMAGE', vimage)
					end
				end

			end
		end			% Loop over number of segments
		multi_segs_str(n_clear) = [];		% Clear the unused info

		% In case of Lines (and Isocs) uicontexts have not been set yet. Do it now.
		ind = isnan(hLine);    hLine(ind) = [];      % Clear unused rows in hLine (due to over-dimensioning)
		if ( ~isempty(hLine) && (strcmp(line_type, 'AsLine') || strcmp(line_type, 'i_file') || got_arrow) )
			if (orig_no_mseg)
				draw_funs(hLine,'line_uicontext')		% Here hLine is actually only a scalar
			else
				draw_funs(hLine,'isochron',multi_segs_str)
			end
		end
	end
	% --------------------- End main loop over files -----------------------------------------

	set(handles.figure1,'pointer','arrow')

	if (handles.no_file)		% Be very carefull, do not trust on the 'geog' estimate donne in show_image (too soon)
		geog = handles.geog;	 is_projected = handles.is_projected;
		handles = guidata(handles.figure1);		handles.geog = geog;	 handles.is_projected = is_projected;
	end
	if (got_nc)
		set_extra_uicb_options(handles, hLine, out_nc)	% Reset two Callbacks in UIContextMenu to offer plot/save
	end
	if (do_nesting)
		set_nest_uicb_options(handles, hLine)			% Do some cleaning of unwanted CBs and set the working one.
	end
	guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function set_nest_uicb_options(handles, hLine)
% Set options particular to the NESTING nature of these rectangles.
% Cleaner solution than to create another set of conditionals in draw_funs
	for (k = 1:numel(hLine))
		h = get(get(hLine(k),'UIContextMenu'),'Children');
		h2 = findobj(h,'-depth',0, 'Label','Delete me');		% To be reused
		set(h2,'Label','Adjust nesting dimensions', 'Call','nesting_sizes(gcbo)')
		hC = findobj(h,'-depth',0, 'Label','Delete inside rect');
		if (k > 1)
			set(hC, 'Label', 'Show nesting info','Call','nesting_sizes(gcbo,''Info'')')
		else
			delete(hC),			h = h(ishandle(h));
		end
		delete(findobj(h,'-depth',0, 'Label','Area under polygon')),		h = h(ishandle(h));
		delete(findobj(h,'-depth',0, 'Label','Register Image')),			h = h(ishandle(h));
		delete(findobj(h,'-depth',0, 'Label','Transplant Image here')),		h = h(ishandle(h));
		delete(findobj(h,'-depth',0, 'Label','Get image from Web Map Server')),		h = h(ishandle(h));
		delete(findobj(h,'-depth',0, 'Label','CMT catalog (Web download)')),		h = h(ishandle(h));
		delete(findobj(h,'-depth',0, 'Label','Custom')),		h = h(ishandle(h));
		delete(findobj(h,'-depth',0, 'Label','Create Mask')),	h = h(ishandle(h));
		hC = findobj(h,'-depth',0, 'Label','Copy');
		if (k < numel(hLine))
			delete(hC)
		else
			set(hC, 'Label', 'New nested grid','Call', 'nesting_sizes(gcbo,''New'')');
		end
		set(hLine(k),'Tag','NEST','UserData',k)
	end

% --------------------------------------------------------------------
function set_extra_uicb_options(handles, hLine, out_nc)
% Reuse two meaningless entries (in this context) of the polygon's UIContextMenu to allow
% either ploting the interior points to a shape_nc polygon or directly save it on file.
% Note that the saving, because is done by the draw_funs 'doSave_formated' ... does what it says.
% OUT_NC is the output of the read_shapenc(fname) call + file name

	h1 = get(get(hLine(1),'UIContextMenu'),'Children');
	h2 = findobj(h1,'-depth',0, 'Label','Line azimuths');	% Resuse this entry that has no sense in this context
	set(h2,'Label','Plot interior points')
	set(h2,'Call',{@nc_plotSave_pts, out_nc.fname, out_nc.PolyIndex, 'plot'})
	set(h2,'Pos',2)			% Move it up on the UIContextMenu order (and hope it doesn't screw)
	h2 = findobj(h1,'-depth',0, 'Label','Point interpolation');	% Reuse this one too
	set(h2,'Label','Save interior points to file')
	set(h2,'Call',{@nc_plotSave_pts, out_nc.fname, out_nc.PolyIndex, 'save'})
	set(h2,'Pos',3)

function nc_plotSave_pts(obj, evt, fname, PolyIndex, opt)
% Plot or save the points in the interior of the selected Outer polygon (of an shape_nc file)

	hLine = gco;
	ud = get(hLine,'UserData');
	s = nc_funs('info',fname);
	x = nc_funs('varget', fname, s.Dataset(PolyIndex(ud) - 4).Name);
	y = nc_funs('varget', fname, s.Dataset(PolyIndex(ud) - 3).Name);
	if (opt(1) == 'p')
		hAx = get(hLine, 'Parent');
		h = line('XData',x,'YData',y,'Parent',hAx, 'LineStyle','none', 'Marker','.',...
		'MarkerSize',2, 'Tag','Pointpolyline');
		draw_funs(h,'DrawSymbol')			% Set marker's uicontextmenu
	else
		z = nc_funs('varget', fname, s.Dataset(PolyIndex(ud) - 2).Name);
		draw_funs([], 'doSave_formated', x, y, z)		% The first arg (empty) means no data fishing in handles
	end
% --------------------------------------------------------------------

% --------------------------------------------------------------------
function [cor, str2] = parseG(str)
% Parse the STR string in search of color. If not found or error COR = [].
% STR2 is the STR string less the -Gr/g/b part
	cor = [];   str2 = str;
	ind = strfind(str,' -G');
	if (isempty(ind)),      return;     end		% No -G option
	try									% There are so many ways to have it wrong that I won't bother testing
		[strG, rem] = strtok(str(ind+1:end));
		str2 = [str(1:ind(1)) rem];		% Remove the -G<str> from STR

		strG(1:2) = [];					% Remove the '-G' part from strG
		% OK, now 'strG' must contain the color in the r/g/b form
		ind = strfind(strG,'/');
		if (isempty(ind))				% E.G. -G100 form
			cor = eval(['[' strG ']']);
			cor = [cor cor cor] / 255;
		else
			% This the relevant part in num2str. I think it is enough here
			cor = [eval(['[' strG(1:ind(1)-1) ']']) eval(['[' strG(ind(1)+1:ind(2)-1) ']']) eval(['[' strG(ind(2)+1:end) ']'])];
			cor = cor / 255;
		end
		if (any(isnan(cor))),   cor = [];   end
	end

% --------------------------------------------------------------------
function [thick, cor, str2] = parseW(str)
% Parse the STR string in search for a -Wpen. Valid options are -W1,38/130/255 -W3 or -W100/255/255
% If not found or error THICK = [] &/or COR = [].
% STR2 is the STR string less the -W[thick,][r/g/b] part
	thick = [];     cor = [];   str2 = str;
	ind = strfind(str,' -W');
	if (isempty(ind)),		return,		end		% No -W option
	try                                 % There are so many ways to have it wrong that I won't bother testing
		[strW, rem] = strtok(str(ind+1:end));
		str2 = [str(1:ind(1)) rem];		% Remove the -W<str> from STR

		strW(1:2) = [];					% Remove the '-W' part from strW
		% OK, now 'strW' must contain the pen in the thick,r/g/b form
		ind = strfind(strW,',');
		if (~isempty(ind))				% First thing before the comma must be the line thickness
			thick = eval(['[' strW(1:ind(1)-1) ']']);
			strW = strW(ind(1)+1:end);	% Remove the line thickness part
		else							% OK, no comma. So we have either a thickness XOR a color
			ind = strfind(strW,'/');
			if (isempty(ind))			% No color. Take it as a thickness
				thick = eval(['[' strW ']']);
			else						% A color
				cor = [eval(['[' strW(1:ind(1)-1) ']']) eval(['[' strW(ind(1)+1:ind(2)-1) ']']) eval(['[' strW(ind(2)+1:end) ']'])];
				cor = cor / 255;
				if (any(isnan(cor))),   cor = [];   end
				% We are done here. RETURN
				return
			end
		end
		% Come here when -Wt,r/g/b and '-Wt,' have already been riped
		ind = strfind(strW,'/');
		if (~isempty(ind))
			% This the relevant part in num2str. I think it is enough here
			cor = [eval(['[' strW(1:ind(1)-1) ']']) eval(['[' strW(ind(1)+1:ind(2)-1) ']']) eval(['[' strW(ind(2)+1:end) ']'])];
			cor = cor / 255;
		end
		% Notice that we cannot have -W100 represent a color because it would have been interpret above as a line thickness
		if (any(isnan(cor))),   cor = [];   end
	end

% --------------------------------------------------------------------
function [proj, str2] = parseProj(str)
% Parse the STR string in search for a +proj4 string indicating data SRS
% The proj4 sting should be one single word e.g. +proj4=latlong or enclose in "+proj4=latlong +datum=..."
% In ater case the "" are sriped from the return PROJ variable
% STR2 is the STR string less the ["]+proj4... part
	proj = [];		str2 = str;
	ind = strfind(str, '+proj');
	if (isempty(ind)),		return,		end			% No proj. Go away.
	
	if (str(ind(1)-1) == '"')				% a "+proj4=... ... ..."
		ind2 = strfind(str, '"');
		if (numel(ind2) < 2)
			disp('Error in proj4 string of this file. Ignoring projection request')
			return
		end
		proj = str(ind(1):ind2(2)-1);
		str(ind(1)-1:ind2(2)) = [];			% Strip the proj string
	else
		proj = strtok(str(ind(1):end));		% For example +proj4=longlat
		str(ind(1):ind(1)+numel(proj)-1) = [];	% Strip the proj string
	end
	str2 = str;

% --------------------------------------------------------------------
function [symbol, symbSize, str2] = parseS(str)
% Parse the STR string in search for a -S<symbol>[size]. Valid options are -S...
% If not found or error SYMBOL = [] &/or SYMBSIZE = [].
% STR2 is the STR string less the -S<symb>[size] part
	symbol = [];	symbSize = 2;	str2 = str;
	ind = strfind(str,' -S');
	if (isempty(ind)),		return,		end		% No -S option
	try
		[strS, rem] = strtok(str(ind+1:end));
		str2 = [str(1:ind(1)) rem];		% Remove the -S<str> from STR

		strS(1:2) = [];					% Remove the '-S' part from strS
		% OK, now 'strS' must contain the symbol and optionally its size
		switch strS(1)
			case 'a',		symbol = '*';
			case 'c',		symbol = 'o';
			case 'd',		symbol = 'd';
			case 'h',		symbol = 'h';
			case 'i',		symbol = 'v';
			case 'n',		symbol = 'p';
			case 'p',		symbol = '.';
			case 's',		symbol = 's';
			case 'x',		symbol = 'x';
			case '+',		symbol = '+';
			otherwise,		symbol = 'o';
		end
		if (numel(strS) > 1),	symbSize = str2double(strS(2:end));		end
		if (isnan(symbSize)),	symbSize = 3;	end
	end

% --------------------------------------------------------------------------------
function out = read_shapenc(fname)
% Read from a shapenc type file FNAME
% OUT is a structure with fields:
%	n_swarms			number of data groups (so to speak)
%	type				Geometry Type like in shapefiles e.g. 'PointZ'
%	BB					The 2D BoundingBox
%	SRS					The referencing system. Usualy a proj4 string (empty if no 'spatial_ref')
%	description			Dataset description (or empty if none)
%	n_PolyOUT			Total number of Outer polygons
%	n_PolyIN			Total number of Inner polygons (irrespective of their parents)
%	Poly(n).OUT.lon		N struct array where N is the number of outer polygons (MAX = n_swarms)
%	Poly(n).OUT.lat				"
%	Poly(n).OUT.desc	Description of the individual polygon (may be different from the dataset desc)
%	PolyIndex			Index of the "lonPolyOUT" variables. With this we can later fetch the
%						interior points bacause we know that Z, Lat, Lon are -2,-3,-4 index behind
%						e.g, this will read the Z's of points inside the PolyOUT_1
%						z = nc_funs('varget', fname, s.Dataset(PolyIndex(1)-2).Name);
%
% If N above > 0 for each Outer polygon we can have M Inner polygons, as in
%	out.Poly(n).IN(i).lon		Where i = 1:M
%	out.Poly(n).IN(i).lat
%
% If there are no OUTer polygons than the code search for polylines stored as "lonPolygon_?" and
% returns the findings as if they were "lonPolyOUT", etc...

	s = nc_funs('info',fname);
	ind = strcmp({s.Attribute.Name},'Number_of_main_ensembles');
	out.n_swarms = s.Attribute(ind).Value;

	ind = strcmp({s.Attribute.Name},'SHAPENC_type');
	out.type = s.Attribute(ind).Value;

	ind = strcmp({s.Attribute.Name},'BoundingBox');
	out.BB = s.Attribute(ind).Value;

	ind = strcmp({s.Attribute.Name},'spatial_ref');
	if (~isempty(ind))
		out.SRS = s.Attribute(ind).Value;
	else
		out.SRS = [];
	end

	ind = strcmpi({s.Attribute.Name},'description');
	if (~isempty(ind))
		out.description = s.Attribute(ind).Value;
	else
		out.description = [];
	end

	out.n_PolyOUT = 0;		out.n_PolyIN = 0;	out.PolyIndex = 0;		% They will be updated later

	% Find the outer polygons
	ind = find(strncmp({s.Dataset.Name},'lonPolyOUT',10));
	n_PolyOUT = numel(ind);
	if (n_PolyOUT)
		out.Poly(n_PolyOUT).OUT.lon = [];		out.Poly(n_PolyOUT).OUT.lat = [];
		out.Poly(n_PolyOUT).OUT.desc = [];
		for (k = 1:n_PolyOUT)
			out.Poly(k).OUT.lon = nc_funs('varget', fname, s.Dataset(ind(k)).Name);
			out.Poly(k).OUT.lat = nc_funs('varget', fname, s.Dataset(ind(k)+1).Name);
			ii = find(strcmpi({s.Dataset(ind(k)-1).Attribute.Name},'name'));
			if (~isempty(ii))
				out.Poly(k).OUT.desc = s.Dataset(ind(k)-1).Attribute(ii(1)).Value;		% Description of this polygon
			else
				out.Poly(k).OUT.desc = 'Nickles';
			end
		end
		out.PolyIndex = ind;		% With this we know that Z = ind - 2; lat = ind - 3 and lon = ind - 4
	else
		out.Poly.OUT.desc = '';		% We need this because it's accessed in main function (but not a perfect idea) 
	end

	% Find the inner polygons (if any)
	ind = find(strncmp({s.Dataset.Name},'lonPolyIN',9));
	n_PolyIN = numel(ind);		% Total number of Inner polygons
	if (n_PolyIN)
		% We need to find which of the swarms those internal polygons belong to.
		tmp = zeros(2,n_PolyIN);
		for (k = 1:n_PolyIN)					% Do a first round to find out what is where
			str = s.Dataset(ind(k)).Name(11:end);
			tmp(:,k) = sscanf(str,'%d_%d');		% First row holds the swarm number and second row the polyg number
		end
		[b, m_first, m_last] = local_unique(tmp(1,:));
		n_PolyIN_groups = numel(b);
		nPolys_in_group = tmp(2,m_last);
		for (k = 1:n_PolyIN_groups)				% Loop over number of groups that have Inner polygons
			n = find( b(k) == 1:n_PolyOUT );
			if (~isempty(n))				% This group has Inner polygs
				out.Poly(n).IN(nPolys_in_group(k)).lon = [];	% Pre-allocate
				out.Poly(n).IN(nPolys_in_group(k)).lat = [];
				start = ind(m_first(k))-2;		% minus 2 because of the 2*i below
				for (i = 1:nPolys_in_group(k))
					out.Poly(n).IN(i).lon = nc_funs('varget', fname, s.Dataset(2*i+start).Name);
					out.Poly(n).IN(i).lat = nc_funs('varget', fname, s.Dataset(2*i+1+start).Name);
				end
			else
				out.Poly(n).IN.lon = [];	% This OUTer polygon has no holes (INners)
				out.Poly(n).IN.lat = [];	% Signal that so it will be easier to parse
			end
		end
	else
		for (k = 1:n_PolyOUT)				% So that we don't have any Poly.OUT without at
			out.Poly(k).IN.lon = [];		% least one corresponding Poly.IN, even if empty
			out.Poly(k).IN.lat = [];
		end
	end
	
	if ((n_PolyOUT + n_PolyIN) == 0)
		% Find the outer polygons
		ind = find(strncmp({s.Dataset.Name},'lonPolygon',10));
		n_PolyOUT = numel(ind);
		if (n_PolyOUT)
			out.Poly(n_PolyOUT).OUT.lon = [];		out.Poly(n_PolyOUT).OUT.lat = [];
			% Due to above cases we need the next to exist as well
			out.Poly(n_PolyOUT).IN.lon = [];		out.Poly(n_PolyOUT).IN.lat = [];
			for (k = 1:n_PolyOUT)
				out.Poly(k).OUT.lon = nc_funs('varget', fname, s.Dataset(ind(k)).Name);
				out.Poly(k).OUT.lat = nc_funs('varget', fname, s.Dataset(ind(k)+1).Name);
			end
		end
	end

	out.n_PolyOUT = n_PolyOUT;
	out.n_PolyIN  = n_PolyIN;

% ------------------------------------------------------------------
function [b, ndx_first, ndx_last] = local_unique(a)
% Striped version of unique that outputs 'first' and 'last'
	numelA = numel(a);
	a = a(:);
	[b, ndx] = sort(a);
	db = diff(b);
	d_last = (db ~= 0);
	d_first = d_last;

	d_last(numelA,1) = true;		% Final element is always a member of unique list.
	d_first = [true; d_first];		% First element is always a member of unique list.
	b = b(d_last);					% Create unique list by indexing into sorted list.
	ndx_last  = ndx(d_last);
	ndx_first = ndx(d_first);
function varargout = load_xyz(handles, opt, opt2)
% Read a generic ascii file that can be, or not, a multi-segment file
%
%	Multi-segs files accept -G, -S -W & -: GMT type options plus a proj4 string for referencing.
%		The proj4 string should be one single word e.g. +proj=latlong or enclosed in "+proj=longlat +datum=..."
%		-S<symb>[size] accepts these GMT type symbol codes <a|c|d|h|i|n|p|s|x|+>
%		-S<symb>[size][+s<scale>][+f][+c<cor>[+c<cor>]]		==> Full syntax
%		 Use +s<scale> with files with 3 columns where 3rd column will be used to determine the symbol color
%		 NOTE that to use this option <scale> must be provided, even if == 1.
%		      A bonus of this option is that GE will plot cylinders with height determined by Z * scale
%		 However, if file has 4 columns we can use fourth column to determine the color. To do that
%		 we use the +f flag. When plotted in GE, cylinders height are set from 3rth column and color from 4rth.
%		 If no +c<cor> is provided, use the current cmap or 'jet(64)' if no fig exists yet.
%		 If only one +c<cor> is provided, use it as constant color for all symbols.
%		 Alternatively use the -G<color> if that option is used
%		 If two +c<cor> are provided, interpolate between them and assign them from z_min to z_max.
%
%	It does also deal with the case of ploting the isochrons.dat
%
%	HANDLES	->	Should be the Mirone handles. However, when this function is used with output
%				HANDLES can be empty([]) or just a simple structure with the field 'no_file' = false
%
%	Optional output:
%
%	out = load_xyz([],fname);               Read and return data in file FNAME
%	[out, multi_str] = load_xyz([],fname);  Read file FNAME and return data in OUT and multisegment strings in MULTI_STR
%	... = load_xyz();                       Like the above but asks the file name here
%
%	Optional
%		OPT can be either [] in which case the filename will be asked here or contain the filename
%		OPT2 can take several values
%			'arrows'		to read a x,y,u,v file
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
%		'>-:'		swap 1st and 2nd columns (assumed as (y,x) -> (x,y)) (The -: can now also be anywhere in the string)
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
%		'>HEAVES'	See ecran.m/saveHeaves_CB()
%		'>POLYMESH'	File contains 'Nesting polygons' (multisegment) to help with unstructured grids creation 
%					Example: >POLYMESH -pol=L-1_G-1_P-1.dat -inc=1000 -interp=1 -data= -grid=0 -binary=0 -single=1
%		'>XY'       Send the data read here to the XYtool (Ecran). File can be single or multi-column & multi-segment

%	Copyright (c) 2004-2016 by J. Luis
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
	tol = 0.5;					% Tolerance for line ploting. Will be set to 0 for Points (symbols) 
	do_project = false;         % We'll estimate below if this holds true
	got_arrow = false;
	got_internal_file = false;	% Internal files are those shipped with Mirone (e.g. isochrons)
	got_nc = false;				% Flag to signal if a shapenc type file comes in
	orig_no_mseg = false;		% Flag to know if original file had multiseg strings to store
	line_type = 'AsLine';
	tag  = 'polyline';
	tagP = '';					% Tag for patches. Some specific tagged files may set it
	struc_vimage = [];
	is_bin = false;				% To flag binary files
	BB = [];					% To eventually hold a BoundingBox
	do_nesting = false;			% To flag when the imported file is a 'Nesting squares'
	do_polymesh = false;		% To flag when the imported file is a 'Nesting polygons' for unstructured grids
	isGSHHS = false;			% To flg when we are dealing with a GSHHS polygon
	goto_XY = false;			% To flag when data will be send to the XYtool
	% -------------------------------------------------------------------------------

	% ------------------- PARSE INPUTS ----------------------------------------------
	n_argin = nargin;
	if (nargout)				% A bit convoluted this test but it's necessary for backward compat reasons
		if (n_argin == 0)
			handles = [];	fname = [];		opt = [];
			n_argin = 2;
		end
		if (isempty(handles))	% When this function is used to just read a file and return its contents
			handles.no_file = false;
			handles.last_dir = cd;			% That's the best we can do (this is needed in put_or_get_file)
		end
		varargout{1} = [];					% So that we always have something to return
	end
	if (n_argin >= 2 && isempty(opt))		% Read a ascii file
		[FileName, PathName, handles] = put_or_get_file(handles, ...
			{'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select File','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
		[lix,lix,EXT] = fileparts(FileName);
	elseif (n_argin >= 2)		% Read a ascii file of which we already know the name (drag N'drop)
		fname = opt;
		[PathName,lix,EXT] = fileparts(fname);			% We need the 'PathName' below
	end
	if (n_argin == 3)
		if (strcmp(opt2, 'AsLine') && strcmpi(EXT, '.nc'))			% A shapenc loaded 'AsLine', but confirm
			drv = aux_funs('findFileType',fname);
			if (strcmp(drv, 'ncshape')),	opt2 = 'ncshape';	end
		end
		if (strcmp(opt2, 'AsArrow')),		got_arrow = true;		% This case does not care about 'line_type'
		elseif (strcmp(opt2, 'ncshape')),	got_nc = true;			%
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
		[bin, n_column, multi_seg, n_headers, isGSHHS, GSHHS_str] = guess_file(fname);
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
			multi_seg = 0;		n_headers = 0;		is_bin = true;
			n_column = bin.nCols;
			if (n_argin < 3),	line_type = 'AsPoint';	end		% If not explicitly set line_type
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
		if (~isempty(out_nc.SRS) && out_nc.SRS(1) == '+')		% Only proj4 are accepted
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
		XMin = 1e50;		XMax = -1e50;    YMin = 1e50;	YMax = -1e50;

		for (k = 1:numel(names))
			fname = names{k};
			j = strfind(fname,filesep);
			if (isempty(j)),    fname = sprintf('%s/%s',PathName, fname);   end		% Need to add path as well 
			if (isempty(n_headers)),    n_headers = NaN;    end
			if (~got_nc)				% Otherwise data was read already
				if (multi_seg)
					[numeric_data, multi_segs_str] = text_read(fname,NaN,n_headers,'>');
				elseif (~is_bin)
					numeric_data = text_read(fname,NaN,n_headers);
				else					% Try luck with a binary file
					[numeric_data, multi_segs_str, multi_seg, BB, XMin, XMax, YMin, YMax] = swallow_bin(handles, fname, bin);
				end
			end

			if (~isa(numeric_data,'cell'))			% File was not multi-segment.
				numeric_data = {numeric_data};
				multi_segs_str = {'> Nikles '};		% Need something (>= 8 chars) to not error further down
				orig_no_mseg = true;
			elseif (~got_nc && strncmpi(multi_segs_str{1}, '>XY', 3))	% A request to call the XYtool to display this data
				goto_XY = true;
				multi_segs_str{1}(2:3) = [];		% Rip the XY identifier
			end

			if (isempty(BB) && ~goto_XY)
				for (i = 1:length(numeric_data))
					XMin = min(XMin,double(min(numeric_data{i}(:,1))));		XMax = max(XMax,double(max(numeric_data{i}(:,1))));
					YMin = min(YMin,double(min(numeric_data{i}(:,2))));		YMax = max(YMax,double(max(numeric_data{i}(:,2))));
				end
			elseif (~isempty(BB))		% Also means ~goto_XY
				XMin = BB(1);		XMax = BB(2);		YMin = BB(3);		YMax = BB(4);
			end
		end

		if (~goto_XY)					% goto_XY means we will call Ecran and not Mirone
			if (multi_seg && ~isempty(strfind(multi_segs_str{1},'-:')))	% See if we need to swap x<->y
				tmp = XMin;		XMin = YMin;	YMin = tmp;				% Need to swap min/max
				tmp = XMax;		XMax = YMax;	YMax = tmp;
			end

			% ----- Check for the special case of only one pt or pure Vertical/Horizontal lines ------
			[XMin, XMax, YMin, YMax] = check_smallness(handles, XMin, XMax, YMin, YMax, numeric_data);
			% ----------------------------------------------------------------------------------------

			xx = [XMin XMax];			yy = [YMin YMax];
			region = [xx yy];
			handles.geog = aux_funs('guessGeog',region);

			if (got_internal_file)				% We know it's geog (Global Isochrons or FZs)
				xx = [-180 180];		yy = [-90 90];
				if (~handles.geog),				handles.geog = 1;
				elseif (handles.geog == 2),		xx = [0 360];
				end
				region = [xx yy];
			end
			mirone('FileNewBgFrame_CB', handles, [region handles.geog])		% Create a background
			hMirFig = handles.figure1;
			drawnow						% Otherwise it takes much longer to plot and other shits
		end

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
		if (handles.no_file && ~goto_XY && k == 1)		% Rename figure with dragged file name
			[pato, barName] = fileparts(fname);
			old_name = get(hMirFig,'Name');		ind = strfind(old_name, '@');
			set(hMirFig, 'Name', [barName old_name(ind-1:end)])
		end

		if (~handles.no_file)					% Otherwise we already read it
			j = strfind(fname,filesep);
			if (isempty(j)),    fname = sprintf('%s%s',PathName, fname);   end		% Need to add path as well 
			if (isempty(n_headers)),    n_headers = NaN;    end
			if (~got_nc)			% Otherwise data was read already
				if (isGSHHS)		% Special case of a "GSHHS Master File" that must be dealt first
					[numeric_data, multi_segs_str] = swallow_GSHHS(handles, fname);
					if (isempty(numeric_data))
						warndlg('There is no GSHHG data inside this region.', 'Warning')
						if (nargout),	varargout{1} = [];	end
						return
					end
				else
					if (multi_seg)
						[numeric_data, multi_segs_str] = text_read(fname,NaN,n_headers,'>');
					elseif (~is_bin)
						numeric_data = text_read(fname,NaN,n_headers);
					else				% Try luck with a binary file
						[numeric_data, multi_segs_str, multi_seg] = swallow_bin(handles, fname, bin);
					end
				end
			end
		end

		if (~isa(numeric_data,'cell'))			% File was not multi-segment. Now pretend it was but with no info
			numeric_data = {numeric_data};
			multi_segs_str = {'> Nikles '};		% Need something in it to not error below
			orig_no_mseg = true;
		end
		n_isoc  = 0;     n_segments = length(numeric_data);
		hLine   = zeros(n_segments,1) * NaN;	% This is the maximum we can have
		hPat    = zeros(n_segments,1) * NaN;	% Or this
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

		elseif (strncmp(multi_segs_str{1}, '>CLOSE', 6))			% Closed or not, plot a patch
			multi_segs_str{1}(2:6) = [];							% Rip the CLOSE identifier
			do_patch = true;

		elseif (strncmpi(multi_segs_str{1}, '>HAVE_INCLUDES', 7))	% This file includes other files
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

		elseif (strncmpi(multi_segs_str{1}, '>HEAVES', 5))			% A Tectonic 'Heaves' file to reconstruct
			[t,r] = strtok(multi_segs_str{1});
			n_segments = size(numeric_data{:},1);
			tmp_cell   = numeric_data;
			numeric_data = cell(n_segments,1);	multi_segs_str = cell(n_segments,1);
			for (kh = 1:n_segments)
				numeric_data{kh} = [tmp_cell{1}(kh,1:2); tmp_cell{1}(kh,3:4)];
				multi_segs_str{kh} = '>';
			end
			multi_segs_str{1} = ['> ' r];							% Rip the HEAVES identifier and leave whathever remains
			hLine = ones(n_segments,1)*NaN;			% We need to update these two too
			n_clear = false(n_segments,1);

		elseif (strncmp(multi_segs_str{1}, '>POLYMESH', 9))			% 
			multi_segs_str{1}(2:9) = [];							% Rip the swap POLYMESH identifier
			do_polymesh = true;
			do_patch = true;		tagP = 'polymesh';
			polymesh_family = cell(n_segments,1);	% pre-allocations
			polymesh_conf(n_segments,1) = struct('inc','', 'interp',0, 'fname','', 'is_grid',0, ...
				'is_binary',0, 'single',1, 'pai_grp',1, 'pai_row',1);
			for (kh = 1:n_segments)
				[multi_segs_str{kh}, polymesh_conf(kh), polymesh_family{kh}, msg] = parse_polymesh(multi_segs_str{kh});
				if (~isempty(msg)),		errordlg(msg, 'Error'),		return,		end
			end

		elseif (strncmp(multi_segs_str{1}, '>-', 2) || strncmp(multi_segs_str{1}, '>+', 2) || strncmp(multi_segs_str{1}, '>"+', 3))
			multi_segs_str{1} = ['> ' multi_segs_str{1}(2:end)];	% Open a space btween '>' and next char

		elseif (line_type(3) ~= 'P' && ~isempty(strfind(multi_segs_str{1},'-G')) && isempty(strfind(multi_segs_str{1},'-S')) )
			% -G (paint) alone is enough to make it a patch (if ~point)
			do_patch = true;
		end

		% ------------------ Check if File has y,x instead of x,y --------------------------
		[do, multi_segs_str{1}] = parseSwap(multi_segs_str{1});
		if (do)
			for (i = 1:n_segments)				% Swapp 1st and 2th columns.
				tmp = numeric_data{i}(:,1);
				numeric_data{i}(:,1) = numeric_data{i}(:,2);
				numeric_data{i}(:,2) = tmp;
			end
			clear tmp do
		end
		% -----------------------------------------------------------------------------------

		% -----------------------------------------------------------------------------------
		% ---------- If OUT is requested there is nothing left to be done here --------------
		% -----------------------------------------------------------------------------------
		if (nargout)
			if (orig_no_mseg),		numeric_data = numeric_data{1};		end
			varargout{1} = numeric_data;
			if (nargout == 2),	varargout{2} = multi_segs_str;		end
			return			% Means this only works for one single file
		end
		% -----------------------------------------------------------------------------------

		% -----------------------------------------------------------------------------------
		% ----- If a call to plot XY is requested there is nothing left to be done here -----
		% -----------------------------------------------------------------------------------
		if (goto_XY)
			n_cols = size(numeric_data{1},2);
			if (n_cols == 1)		% If only one column, xx will be 1:npts
				y = numeric_data{1};	x = 1:numel(y);
				ecran(x, y, fname)
				for (i = 2:n_segments)
					ecran('add', 1:numel(numeric_data{i}), numeric_data{i})
				end
			elseif (n_cols == 2)
				ecran(numeric_data{1}(:,1), numeric_data{1}(:,2), fname)
				for (i = 2:n_segments)
					ecran('add', numeric_data{i}(:,1), numeric_data{i}(:,2))
				end
			else
				out = select_cols(numeric_data{1}, 'xy', fname, 1000);
				if (isempty(out)),		return,		end
				h = [];
				if (numel(out) == 4)				% x,y,z and distance request but we don't use it here
					out(2) = out(3);
					h = warndlg('Sorry, computing distance is not supported here. Must open the file directly in XYtool.','Warning');
				end
				if (n_segments > 1)
					h = warndlg('With multiple coluns and multi-segments, only first segment is processed.', 'Warning');
				end
				ecran(numeric_data{1}(:,out(1)), numeric_data{1}(:,out(2)), fname)
				if (~isempty(h)),	figure(h),		end		% Bring warning message to top
			end
			continue				% Continue the loop over input files (main loop)
		end
		% -----------------------------------------------------------------------------------

		% ------------------ Check if it is a GSHHS or WDBII file ---------------------------
		if (isGSHHS)
			tag = 'GMT_DBpolyline';
		end
		% -----------------------------------------------------------------------------------

		% ------------------ Check if first header line has a PROJ4 string ------------------
		[projStr, multi_segs_str{1}] = parseProj(multi_segs_str{1});
		if (~isempty(projStr)),		do_project = true;		end
		% -----------------------------------------------------------------------------------

		if (do_project && handles.no_file)		% If new image set the projection info
			aux_funs('appProjectionRef', handles, projStr)
			handles.is_projected = true; 
		end

		drawnow
		for (i = 1:n_segments)		% Loop over number of segments of current file (external loop)
			tmpz = [];
			if (do_project)         % We need to project
				try
					if (~isempty(projStr))
						[numeric_data{i}, msg] = proj2proj_pts(handles, numeric_data{i}, 'srcProj4', projStr);
						if (~isempty(msg) && strncmp(msg, 'unknown', 7))
							warndlg('Cannot project because destination reference system is unknown.','Warning')
							do_project = false;
						end

						if (k == 1 && i == 1 && do_project)		% First time. Store proj info in Figure's appdata. IF ...
							if (~handles.geog && isempty(getappdata(handles.figure1,'ProjWKT')) && ...
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

			% Check if we have a symbols request. If yes turn the 'AsPoint' option on, but do the potential trimming first
			[marker, markerSize, markerScale, color_by4, cor1_sc, cor2_sc, multi_segs_str{i}] = parseS(multi_segs_str{i});
			if (~isempty(marker) || strcmp(line_type, 'AsPoint') || strcmp(line_type, 'AsMaregraph'))
				tol = 0;
			end
			% ----- Result consequences of this option is resumed further down, but before we have to check trimming -----

			indx = false;	indy = false;			% Default to no need for map clipping
			if (handles.no_file)
				tmpx = numeric_data{i}(:,1);		tmpy = numeric_data{i}(:,2);
			else
	 			difes = [(double(numeric_data{i}(1,1)) - double(numeric_data{i}(end,1)) ) ...	% Remember R13
					(double(numeric_data{i}(1,2)) - double(numeric_data{i}(end,2)))];
				if (any(abs(difes) > 1e-5))			% Assume a not closed polygon
					[tmpx, tmpy, indx, indy] = ...	% Get rid of points that are outside the map limits
						aux_funs('in_map_region',handles,numeric_data{i}(:,1),numeric_data{i}(:,2),tol,[xx yy]);
				elseif (~isGSHHS)
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
			if (numel(numeric_data{i}(1,:)) >= 3)		% If we have a Z column
				tmpz = numeric_data{i}(:,3);
				if (~isempty(indx) || ~isempty(indy)),	tmpz(indx) = [];	tmpz(indy) = [];	end	% If needed, clip outside map data
			end

			multi_segs_str{i} = deblank(multi_segs_str{i});
			[lThick, cor, multi_segs_str{i}] = parseW(multi_segs_str{i}(min(2,numel(multi_segs_str{i})):end)); % First time, we can chop the '>' char
			if (isempty(lThick)),	lThick = handles.DefLineThick;	end		% IF not provided, use default
			if (isempty(cor)),		cor = handles.DefLineColor;		end		%           "

			% ---- Resume case analysis of -S option in multi-seg that had to be intrrrupted above because of trimming ----
			if (color_by4 && numel(numeric_data{i}(1,:)) >= 4)
				tmpz4 = numeric_data{i}(:,4);		% Will be used to color symbols
				if (~isempty(indx) || ~isempty(indy)), tmpz4(indx) = [];	tmpz4(indy) = [];	end	% If needed, clip outside map data
			else
				color_by4 = false;					% bad +f setting. Just ignore it
			end
			if (~isempty(marker))
				if (~isempty(markerScale) && ~isempty(tmpz))
					line_type = 'scaled';
				else
					line_type = 'AsPoint';
				end
			else
				marker = 'o';	markerSize = 2;		% The old defaults
			end
			% --------------------------------------------------------------------------------------------------------------

			if (do_patch)
				Fcor = parseG(multi_segs_str{i});
				if (isempty(Fcor)),		Fcor = 'none';		end
				hPat(i) = patch('XData',tmpx, 'YData',tmpy, 'Parent',handles.axes1, 'Linewidth',lThick, ...
						'EdgeColor',cor, 'FaceColor',Fcor, 'Tag', tagP);
				if (~isempty(tmpz))
					set(hPat(i),'UserData',tmpz');
				end	
				draw_funs(hPat(i),'line_uicontext')
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
							hLine(i) = line('XData',tmpx, 'YData',tmpy, 'Parent',handles.axes1, 'Linewidth',lThick,...
									'Color',cor, 'Tag',tag, 'Userdata',n_isoc);
							if ~((numel(multi_segs_str{i}) <= 2) && strfind(multi_segs_str{i}, '>')) % Sometimes we still have only '>'
								setappdata(hLine(i),'LineInfo',multi_segs_str{i});
								if (isGSHHS)	% Need to store this for writing edited GSHHG polygons.
									setappdata(hLine(i),'GSHHS_str', GSHHS_str);
								end
							end
							setappdata(hLine(i),'was_binary',is_bin);	% To offer option to save as binary too
						case 'AsPoint'
							Fcor = parseG(multi_segs_str{i});			% See if user wants colored pts
							if (isempty(Fcor)),		Fcor = 'k';		end
							hLine(i) = line('XData',tmpx,'YData',tmpy,'Parent',handles.axes1, 'LineStyle','none', 'Marker',marker,...
								'MarkerEdgeColor','k','MarkerFaceColor',Fcor, 'MarkerSize',2,'Tag','Pointpolyline');
							if (~isempty(strtok(multi_segs_str{i})))
								setappdata(hLine(i),'LineInfo',multi_segs_str{i});
							end
							draw_funs(hLine(i),'DrawSymbol')			% Set marker's uicontextmenu (tag is very important)
							setappdata(hLine(i),'was_binary',is_bin);	% To offer option to save as binary too
						case 'AsMaregraph'
							hLine(i) = line('XData',tmpx,'YData',tmpy,'Parent',handles.axes1, 'LineStyle','none', 'Marker','o',...
								'MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',10,'Tag','Maregraph');
							draw_funs(hLine(i),'DrawSymbol')			% Set marker's uicontextmenu					
						case 'FaultTrace'
							hLine(i) = line('XData',tmpx,'YData',tmpy,'Parent',handles.axes1,'Color',cor,'LineWidth',lThick, ...
								'Tag','FaultTrace');
							draw_funs(hLine(i),'line_uicontext')		% Set lines's uicontextmenu
							% Create empty patches that will contain the surface projection of the fault plane
							hp = zeros(1, numel(tmpx)-1);
							for (j = 1:numel(tmpx)-1),	hp(j) = patch('XData', [], 'YData',[]);    end
							setappdata(hLine(i),'PatchHand',hp);
						case 'scaled'
							ind_NaN = isnan(tmpz);
							if (any(ind_NaN))
								tmpx(ind_NaN) = [];	tmpy(ind_NaN) = [];		tmpz(ind_NaN) = [];
								if (color_by4),		tmpz4(ind_NaN) = [];	end
							end
							nPts = numel(tmpx);
							symbSIZES = repmat(markerSize,nPts,1);
							if (isempty(cor1_sc))	% Make another attempt to see if a color was set by -G
								cor1_sc = parseG(multi_segs_str{i});
							end
							z_colCor = tmpz;
							if (color_by4),		z_colCor = tmpz4;	end
							Zmin = min(z_colCor);	Zmax = max(z_colCor);
							dZ = Zmax - Zmin;

							if (~isempty(cor1_sc))			% Case of colors set in -S option (not finished)
								if (isempty(cor2_sc))
									zC = repmat(cor1_sc, nPts, 1);
								else
									rn = (z_colCor - Zmin) / dZ;	% range normalized to [0 1]
									rn(isnan(rn)) = 1;		% Happens when Zmax == Zmin
									rn = rn(:)';
									dc = cor2_sc - cor1_sc;
									zC = [(cor1_sc(1) + dc(1) * rn)' (cor1_sc(2) + dc(2) * rn)' (cor1_sc(3) + dc(3) * rn)'];							
								end
							else
								if (handles.no_file)		% In this case the fig cmap is all whites
									cmap = jet(64);
								else
									cmap = get(handles.figure1,'ColorMap');
								end
								if (dZ == 0)				% Cte color
									zC = repmat(cmap(round(size(cmap,1)/2),:),nPts,1);      % Middle color
								else
									zC = round(((z_colCor - Zmin) / dZ) * (size(cmap,1)-1) + 1);
									zC = cmap(zC,:);
								end
							end
							tmpz = abs(tmpz);				% Currently the Z is only used to make cylinders in GE
							if (markerScale ~= 1),	tmpz = tmpz * markerScale;		end
							if (~isempty(cor1_sc) && isempty(cor2_sc))	% Unique color, we can plot them all in one single line
								hLine(i) = line('XData',tmpx,'YData',tmpy, 'Parent',handles.axes1, ...
									'LineStyle','none', 'Tag','scatter_symbs', 'Marker',marker,'Color','k', ...
									'MarkerFaceColor',cor1_sc, 'MarkerSize',symbSIZES(1));					
								setappdata(hLine(i),'ZData',tmpz)
								draw_funs(hLine(i),'DrawSymbol')			% Set marker's uicontextmenu					
							else
								h = zeros(1,nPts);
								for (ks = 1:nPts)
									h(ks) = line('XData',tmpx(ks),'YData',tmpy(ks),'Parent',handles.axes1, ...
										'LineStyle','none', 'Tag','scatter_symbs', 'Marker',marker,'Color','k', ...
										'MarkerFaceColor',zC(ks,:), 'MarkerSize',symbSIZES(ks));
									setappdata(h(ks),'ZData',tmpz(ks))
									draw_funs(h(ks),'DrawSymbol')
								end
							end
							tmpz = [];				% Because later we test this for other purposes
							orig_no_mseg = true;	% Also cheat here for the same reason
					end
				else
					struc_arrow.color = cor;
					hQuiver = draw_funs([], 'loc_quiver', struc_arrow, tmpx, tmpy, UV{i}(:,1), UV{i}(:,2));
					set(hQuiver,'Tag','Seta','Userdata',n_isoc)
					setappdata(hQuiver(1),'MyHead',hQuiver(2))		% Store the arrows heads handle
					hLine(i) = hQuiver(1);
				end

				if (~isempty(tmpz)),	setappdata(hLine(i),'ZData',tmpz');	end
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

			end		% END do_patch or line
		end			% Loop over number of segments
		multi_segs_str(n_clear) = [];		% Clear the unused info

		% In case of Lines (and Isocs) uicontexts have not been set yet. Do it now.
		hLine(isnan(hLine)) = [];			% Clear unused rows in hLine (due to over-dimensioning)
		hLine(hLine == 0) = [];				% The recursive call is very complicated to manage and may end up with zeros here
		if (~isempty(hLine) && (strcmp(line_type, 'AsLine') || strcmp(line_type, 'i_file') || got_arrow))
			if (orig_no_mseg)
				draw_funs(hLine,'line_uicontext')		% Here hLine is actually only a scalar
			else
				% Sometimes we still have only '>'. If only one segment test for that case and jump if true
				if (numel(multi_segs_str) == 1) && (numel(multi_segs_str{i}) <= 2) && strfind(multi_segs_str{i}, '>')
					draw_funs(hLine,'isochron', '')
				else
					draw_funs(hLine,'isochron', multi_segs_str)
				end
			end
		end

	end
	% --------------------- End main loop over files -----------------------------------------
	if (goto_XY),	return,		end		% Time now to abandon this function (it had a 'continue' when its job was finish)

	set(handles.figure1,'pointer','arrow')

	if (handles.no_file)		% Be very carefull, do not trust on the 'geog' estimate donne in show_image (too soon)
		geog = handles.geog;	 is_projected = handles.is_projected;
		handles = guidata(handles.figure1);		handles.geog = geog;	 handles.is_projected = is_projected;
	end
	if (got_nc)
		set_extra_uicb_options(handles, hLine, out_nc)	% Reset two Callbacks in UIContextMenu to offer plot/save
	end
	if (do_nesting)
		draw_funs([],'set_recTsu_uicontext', hLine)		% Set uicontextmenus appropriate for grid nesting
	elseif (do_polymesh)
		mesher_helper('set_props_from_outside', hPat, polymesh_family, polymesh_conf);
	end
	guidata(handles.figure1,handles)

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
% Parse the STR string in search for a +proj string indicating data SRS
% The proj4 string should be one single word e.g. +proj=latlong or enclosed in "+proj=latlong +datum=..."
% In later case the "" are sripped away from the return PROJ variable
% STR2 is the STR string stripped the ["]+proj... part
	proj = [];		str2 = str;
	ind = strfind(str, '+proj');
	if (isempty(ind)),		return,		end			% No proj. Go away.
	
	if (str(ind(1)-1) == '"')				% a "+proj=... ... ..."
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
	if (proj(6) == '4')						% A wrong +proj4 spelling. Remove the extra '4'
		proj(6) = [];
	end
	str2 = str;

% ---------------------------------------------------------------------------
function [symbol, symbSize, scale, color_by4, cor1, cor2, str2] = parseS(str)
% Parse the STR string in search for a -S<symbol>[size][+s<scale>][+f][+c<cor>[+c<cor>]]
% If not found or error SYMBOL = [] &/or SYMBSIZE = [].
% STR2 is the STR string less the -S... part
	symbol = [];	symbSize = 7;	scale = [];	cor1 = [];	cor2 = [];	color_by4 = false;
	symb_pos = 4;	% default symbol start position in string
	str2 = str;
	ind = strfind(str,'-S');
	if (isempty(ind)),		return,		end		% No -S option
	symbol = 'o';		% OK, now we have a new default
	try                                 % There are so many ways to have it wrong that I won't bother testing
		[strS, rem] = strtok(str(ind:end));
		str2 = [str(1:ind(1)-1) rem];   % Remove the -S<str> from STR

		if (numel(strS) > 2)            % Get the symbol
			symb_g = 'acdhinpsx+';
			symb_m = '*odhvp.sx+';
			ind = strfind(symb_g,strS(3));						
			if (~isempty(ind))
				symbol = symb_m(ind);
				if (symb_g(ind) == '+'),	strS(3) = '_';		end		% Don't let the '+' symbol be mistaken by an a flag
			elseif ((double(strS(3)) >= 48) && (double(strS(3)) <= 57))
				symb_pos = 3;
			end
		end

		ind_p = strfind(strS, '+');
        if (numel(strS) > 3)            % Get size
			last = numel(strS);
			if (~isempty(ind_p)),	last = ind_p(1) - 1;	end
            symbSize = str2double(strS(symb_pos:last));
            if (isnan(symbSize)),    symbSize = 7;   end
        end
		if (~isempty(ind_p))
			ind_p(end+1) = numel(strS) + 1;		% A fake last one to easy up algo
			for (k = 1:numel(ind_p))
				switch strS(ind_p(k)+1)
					case 's'
						scale = str2double(strS(ind_p(k)+2:ind_p(k+1)-1));
					case 'f'
						color_by4 = true;
					case 'c'
						if (isempty(cor1))
							[cor1, count] = sscanf(strS(ind_p(k)+2:ind_p(k+1)-1),'%d/%d/%d');
							cor1 = cor1' / 255;
						else
							[cor2, count] = sscanf(strS(ind_p(k)+2:ind_p(k+1)-1), '%d/%d/%d');
							cor2 = cor2' / 255;
						end
						if (count ~= 3)
							warndlg('Error in parsing symbol color. Baddly formed color string. Ignoring it','WarnError')
							cor1 = [];	cor2 = [];
						end
				end
			end
		end
	end

% --------------------------------------------------------------------------------
function [do, str] = parseSwap(str)
% Parse the STR string in search for a -: flag
	do = false;
	ind = strfind(str,' -:');
	if (isempty(ind)),		return,		end		% No -: option
	do = true;
	str(ind(1):ind(1)+2) = [];			% Remove the -: from STR

% --------------------------------------------------------------------------------
function [XMin, XMax, YMin, YMax] = check_smallness(handles, XMin, XMax, YMin, YMax, numeric_data)
% Check for the special case of only one pt and pure vertical or horizontal lines.
% If any of such case is found, change the limits to accoomodate a small padding zone.
	only_one_pt = false;
	dx = XMax - XMin;			dy = YMax - YMin;
	if (dx == 0 && dy == 0)
		only_one_pt = true;
	else
		n = size(numeric_data{1}, 1);
		for (i = 2:length(numeric_data))
			n = max(n, size(numeric_data{i}, 1));
		end
		if (n == 1),	only_one_pt = true;		end
	end
	if (only_one_pt)			% OK, give it a small padding zone. Easy on geogs but trickier on others
		if (handles.geog)
			XMin = max(-360, XMin - 0.25);		XMax = min(360, XMax + 0.25);
			YMin = max(-90,  YMin - 0.25);		YMax = min(90,  YMax + 0.25);
		else
			XMin = XMin * (1 - 0.02);			XMax = XMax * (1 + 0.02);	% Just a 2% padding zone
			YMin = YMin * (1 - 0.02);			YMax = YMax * (1 + 0.02);
		end
	elseif (dx == 0)
		if (handles.geog)
			XMin = max(-360, XMin - 0.25);		XMax = min(360, XMax + 0.25);
		else
			XMin = XMin * (1 - 0.02);			XMax = XMax * (1 + 0.02);	% Just a 2% padding zone
		end
	elseif (dy == 0)
		if (handles.geog)
			YMin = max(-90,  YMin - 0.25);		YMax = min(90,  YMax + 0.25);
		else
			YMin = YMin * (1 - 0.02);			YMax = YMax * (1 + 0.02);
		end
	end

	if (XMin > -179.5 && XMax < 359.5 && YMin > -89.5 && YMax < 89.5)
		XMin = XMin - dx / 200;		XMax = XMax + dx / 200;		% Give an extra 0.5% padding margin
		YMin = YMin - dy / 200;		YMax = YMax + dy / 200;
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
	if (~isempty(ind) && any(ind))
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

% ---------------------------------------------------------------------
function [numeric_data, multi_segs_str] = swallow_GSHHS(handles, fname)
% Read GMT GSHHS DB files. Specially the GSHHS_f_Level_1.txt

	if (isfield(handles,'ROI_rect'))	% As for example set by deal_opts->load_GMT_DB()
		XYlim = handles.ROI_rect;
	else
		XYlim = getappdata(handles.axes1,'ThisImageLims');
	end
	if (XYlim(1) < 0),	XYlim(1) = XYlim(1) + 360;	end		% At the Scripts GMT Sumit this failed but I'm not able to reproduce
	if (XYlim(2) < 0),	XYlim(2) = XYlim(2) + 360;	end
	[bin,n_column,multi_seg,n_headers] = guess_file(fname);		% Repeated op but never mind, it's a fast one
	[hdrs, ind] = txt2mat(fname, 'NumHeaderLines',n_headers,'Format','%s','GoodLineString',{'>'},'ReadMode','cell');
	nTot = numel(hdrs);
	numeric_data = cell(nTot,1);
	c = true(nTot,1);
	P1.x = [XYlim(1) XYlim(1) XYlim(2) XYlim(2) XYlim(1)];	P1.hole = 0;
	P1.y = [XYlim(3) XYlim(4) XYlim(4) XYlim(3) XYlim(3)];	P2.hole = 0;
	hBar = aguentabar(0,'title','Scanning a GSHHS DB file');
	for (k = 1:nTot)
		idR = strfind(hdrs{k}, 'R =');		idA = strfind(hdrs{k}, 'A =');
		lim = sscanf(hdrs{k}(idR(1)+4:idA(1)-2), '%f/%f/%f/%f');
		% Compute the intersection of two rectangles
		P2.x = [lim(1) lim(1) lim(2) lim(2) lim(1)];
		P2.y = [lim(3) lim(4) lim(4) lim(3) lim(3)];
		P3 = PolygonClip(P1, P2, 1);				% Intersection of the two rectangles
		if (isempty(P3))
			c(k) = false;
			continue
		end
		idN = strfind(hdrs{k}, 'N =');		idG = strfind(hdrs{k}, 'G =');
		nPts = sscanf(hdrs{k}(idN(1)+4:idG(1)-2), '%d');
		numeric_data{k} = ...
			txt2mat(fname,'RowRange',[1,nPts], 'FilePos',ind(k),'NumHeaderLines',1,'NumColumns',2,'Format','%f %f', 'InfoLevel',1);
		if (rem(k,100) == 0)
			aguentabar(k/nTot)
		end
	end
	multi_segs_str = hdrs(c);
	numeric_data(~c) = [];		% Remove unused cells
	if (ishandle(hBar)),	aguentabar(1),	end

% --------------------------------------------------------------------------------------------------------------
function [numeric_data, multi_segs_str, multi_seg, BB, XMin, XMax, YMin, YMax] = swallow_bin(handles, fname, bin)
% This function deals with the reading of binary files.
% BIN is the struct as issued in output by guess_file. If there was a request for spliting a NaN separated
%     file into multisegments, than there is a bit of work to do here.
% BB  is the Bounding Box, but we only compute it here when no multisegment splitting because in later case
%     the main code will compute, as for the other (ASCII) multisegment cases

	XMin = 1e50;		XMax = -1e50;    YMin = 1e50;	YMax = -1e50;	BB = [];
	multi_segs_str = '';	multi_seg = false;

	fid = fopen(fname);		numeric_data = fread(fid,['*' bin.type]);	fclose(fid);
	numeric_data = reshape(numeric_data,bin.nCols,numel(numeric_data)/bin.nCols)';
	if (bin.twoD && (bin.nCols > 2)),	numeric_data(:,3:end) = [];		end % Retain only X,Y
	if (~bin.multiseg && ~handles.version7 && strcmp(bin.type, 'single') && any(isnan(numeric_data(:,1))))
		fds = double(numeric_data(:,1));		% It means fdsse, ML engeneers after 20 years still think that
		XMin = min(fds);	XMax = max(fds);	% min(single([1 2 NaN])) = NaN
		fds = double(numeric_data(:,2));
		YMin = min(fds);	YMax = max(fds);
		BB = [XMin XMax YMin YMax];
		clear fds
	elseif (bin.multiseg)			% A request to split a multiseg file into individual lines
		indNaN = find(isnan(numeric_data(:,1)));
		if (any(indNaN))
			indNaN = indNaN(:);
			difa = diff(indNaN);
			indNaN(difa == 1) = [];	% Remove contiguos NaNs references
			if (indNaN(1) ~= 1),	indNaN = [1; indNaN];	end		% Pretend that first and last were NaNs
			if (indNaN(end) ~= size(numeric_data,1)),	indNaN(end+1) = size(numeric_data,1)+1;	end
			if (~isempty(indNaN))	% It could be if we have NaNs only at first and/or last positions
				num_data = cell(numel(indNaN)-1, 1);
				multi_segs_str = cell(numel(indNaN)-1, 1);
				for (k = 1:numel(indNaN)-1)
					num_data{k} = numeric_data( (indNaN(k)+1:indNaN(k+1)-1), :);
					multi_segs_str{k} = '> Nikles ';		% Not sure if we really need this
				end
				numeric_data = num_data;
				multi_seg = true;	% From now on pretend the file had always been multiseg
			end
		end
	end

% ---------------------------------------------------------------------
function [str, conf, family, msg] = parse_polymesh(str)
% Parse the contents of multi-seg header of a POLYMESH type file for the 'config' params
% ex: -pol=L-1_G-1_P-1.dat -inc=1000 -interp=1 -data= -grid=0 -binary=0 -single=1 -pai_grp=3  pai_row=1

	msg = '';
	conf = struct('inc','', 'interp',0, 'fname', '', 'is_grid',0, 'is_binary',0, 'single',1, 'pai_grp',1, 'pai_row',1);
	ind = strfind(str,'-pol=');
	if (isempty(ind))
		family = 0;
		msg = 'Fatal error: meshpolygon group is broken, one of its elements misses the hierarchy nesting info';
		return
	end
	family = strtok(str(ind(1)+5:end));
	str(ind(1):ind(1)+5+numel(family)-1) = [];		% Remove this entry from the input string
	
	[str, param] = parse_pm_one(str, '-inc=');
	if (~isempty(param)),	conf.inc = param;		end
	
	[str, param] = parse_pm_one(str, '-interp=');
	if (~isempty(param)),	conf.interp = str2double(param);	end
	
	[str, param] = parse_pm_one(str, '-data=');
	if (~isempty(param)),	conf.fname = param;		end
	
	[str, param] = parse_pm_one(str, '-grid=');
	if (~isempty(param)),	conf.is_grid = str2double(param);	end
	
	[str, param] = parse_pm_one(str, '-binary=');
	if (~isempty(param)),	conf.is_binary = str2double(param);	end
	
	[str, param] = parse_pm_one(str, '-single=');
	if (~isempty(param)),	conf.single = str2double(param);	end
	
	[str, param] = parse_pm_one(str, '-pai_grp=');
	if (~isempty(param)),	conf.pai_grp = str2double(param);	end
	
	[str, param] = parse_pm_one(str, '-pai_row=');
	if (~isempty(param)),	conf.pai_row = str2double(param);	end

% ---------------------------------------------------------------------
function [str, param] = parse_pm_one(str, tok)
% Parse a token from the polymesh header. TOK has the form (ex:) '-single='

	param = '';
	ind = strfind(str,tok);
	if (~isempty(ind))
		try					% Wrap it with a try to prevent out of bounds access errors
			param = strtok(str(ind(1)+numel(tok):end));
			if (param(1) == '-'),	param = '';		end		% Happens when pram= (nothing) and strtok gets next token
			% Remove this entry from the input string
			str(ind(1):ind(1)+numel(tok)+numel(param)-1) = [];
		end
	end

% ---------------------------------------------------------------------
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

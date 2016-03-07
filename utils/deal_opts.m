function out = deal_opts(opt, opt2, varargin)
% Function to provide fine (manual) control on some Mirone modules.
%
% The mechanism of this function relies on the contents of the "OPTcontrol.txt" control file.
% That file has keewords that trigger the usage of optional (specialized) features in Mirone.
%
% OPT directs to the option of interest. It may be a cell array of options.
% OPT2 can be a uimenu handle for cases where the OPT case needs one
%
% OPT2 optionally may hold the name of an internal sub-function, in which case that function
% is called with optional extra arguments transmited in varargin

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

	if (nargin >= 2 && ischar(opt2))
		if (nargout)
			out = feval(opt2, varargin{:});
		else
			feval(opt2, varargin{:});
		end
		return
	end

	MIRONE_DIRS = getappdata(0,'MIRONE_DIRS');
	opt_file = [MIRONE_DIRS.home_dir filesep 'data' filesep 'OPTcontrol.txt'];
	out = [];
	hCust = [];

	if (~exist(opt_file, 'file') == 2),		return,		end
	fid = fopen(opt_file, 'r');
	c = (fread(fid,'*char'))';      fclose(fid);
	lines = strread(c,'%s','delimiter','\n');   clear c fid;
	m = numel(lines);
	c = true(1, m);

	for (k = 1:m)
		if (strncmp(lines{k},'MIR_',4)),	c(k) = false;	end
	end
	lines(c) = [];			% Delete non-keyword (comment) lines
	if (isempty(lines)),	return,		end

	if (~isa(opt,'cell'))	% Some calls (e.g. MGG & Microleveling) come as cells
		opt = {opt};
	end

	for (k = 1:numel(lines))		% Loop over the keyword lines of the OPTcontrol.txt file
		for (n = 1:numel(opt))		% Somtimes one single call holds multiple choices
			switch lower(opt{n})
				case 'gmtedit'		% To select what to plot in GMTEDIT slots
					if (strncmp(lines{k}(5:end),'GMTEDIT',7))
						t = strtok(lines{k}(13:end));
						if ( strcmp(t(1:2), '-V') )		% Here we only check for a -V... and do not check for errors
							out = t;
						end
						break
					end

				case 'mgg'			% To find tracks inside calling rectangle
					if (strncmp(lines{k}(5:end),'MGG',3))
						if (isempty(hCust))			% Create a uimenu associated to a rectangle
							hCust = uimenu(opt2, 'Label', 'Custom','Sep','on');
						end
						uimenu(hCust, 'Label', 'MGG tracks', 'Call',  @get_MGGtracks);
						break
					end

				case 'mgg_coe'		% To plot Cross Over Errors of MGG data
					if (strncmp(lines{k}(5:end),'COEs',4))
						out = lines{k}(10:end);
						coeVar = 'mag';			% Default
						[t, r] = strtok(out);
						if (~isempty(r)),		coeVar = ddewhite(r);	end
						item = uimenu(opt2, 'Label', 'Show Cross Over Errors', 'Sep','on');
						uimenu(item, 'Label', 'All', 'Call', {@get_COEs, t, coeVar, 'all'});
						uimenu(item, 'Label', 'External only', 'Call', {@get_COEs, t, coeVar, 'ext'});
						uimenu(item, 'Label', 'Internal only', 'Call', {@get_COEs, t, coeVar, 'int'});
						uimenu(item, 'Label', 'Plot histogram (all)', 'Call', {@get_COEs, t, coeVar, 'hstA'},'Sep','on');
						uimenu(item, 'Label', 'Plot histogram (internal)', 'Call', {@get_COEs, t, coeVar, 'hstI'});
						%uimenu(opt2, 'Label', 'Show COEs', 'Call', {@get_COEs, t, coeVar},'Sep','on');
						break
					end

				case 'microlev'		% To do microleveling on a rectangular ROI
					if (strncmp(lines{k}(5:end),'MICRO',5))
						if (isempty(hCust))	% Create a uimenu associated to a rectangle
							hCust = uimenu(opt2, 'Label', 'Custom','Sep','on');
						end
						uimenu(hCust, 'Label', 'ROI microleveling', 'Call', 'microlev(gco)');
						break
					end

				case 'gmt_db_ids'		% To report the IDs of the GMT database polygons inside rectangle
					if (strncmp(lines{k}(5:end),'GMT_DB_IDS',10))
						if (isempty(hCust))	% Create a uimenu associated to a rectangle
							hCust = uimenu(opt2, 'Label', 'Custom','Sep','on');
						end
						uimenu(hCust, 'Label', 'Show GMT db polygon IDs', 'Call', @sow_GMT_DB_IDs);
						uimenu(hCust, 'Label', 'Load GMT db polygon(s)', 'Call',  @load_GMT_DB);
						break
					end

				case 'gmt_symbol'		% To plot a GMT custom symbol ~inside rectangle
					if (strcmp(lines{k}(5:end),'GMT_SYMBOL'))
						if (isempty(hCust))	% Create a uimenu associated to a rectangle
							hCust = uimenu(opt2, 'Label', 'Custom','Sep','on');
						end
						uimenu(hCust, 'Label', 'Associate to a GMT symbol', 'Call', @assoc_gmt_symbol, 'Sep', 'on');
						h = get(opt2, 'UserData');		% Fish the line handle
						fname = getappdata(h, 'cust_symb');	% If not empty than we have a call from FileOpenSession_CB
						if (~isempty(fname)),	assoc_gmt_symbol([], [], h),	end		% and want only to set the BDF CB
						break
					end
			end
		end
	end

% -----------------------------------------------------------------------------------------
function get_COEs(obj, event, coeFile, coeVar, opt)
% Get a list of MGG COEs involving gco and plot them with uicontexts
% 'coeVar'	is name of variable to process
% 'opt'		is a char option informing which COEs to plot. Valid options are 'all', 'ext', 'int'

	h = gco;		handles = guidata(h);

	lims = getappdata(handles.axes1, 'ThisImageLims');
	opt_R = sprintf(' -R%.4f/%.4f/%.4f/%.4f',lims(1:4));
	MIRONE_DIRS = getappdata(0,'MIRONE_DIRS');
	tmp_file = [MIRONE_DIRS.home_dir filesep 'tmp' filesep 'MGGtmp.txt'];
	fname = getappdata(h, 'FullName');
	[pato, fname] = fileparts(fname);

	% Decide what kind of COEs to plot or histogramize
	doHistogram = false;
	if ((nargin == 4) || strncmp(opt,'all',3))
		opt_Q = '';
	elseif (strncmp(opt,'ext',3))
		opt_Q = ' -Qe';
	elseif (strncmp(opt,'int',3))
		opt_Q = ' -Qi';
	elseif (strncmp(opt,'hst',3))
		if (opt(end) == 'A'),	opt_Q = '';			% All COES
		else					opt_Q = ' -Qi';		% Internal COEs only
		end
		opt_R = '';
		doHistogram = true;
	end

	if (ispc)
		dos( ['x2sys_list ' coeFile ' -TMGD77+ -Fnxyc -C' coeVar ' -S' fname opt_Q opt_R ' > ' tmp_file]);
	else
		unix(['x2sys_list ' coeFile ' -TMGD77+ -Fnxyc -C' coeVar ' -S' fname opt_Q opt_R ' > ' tmp_file]);
	end

	fid = fopen(tmp_file);
	fgetl(fid);		fgetl(fid);		fgetl(fid);		% Jum the 3 header lines
	c = fread(fid,inf,'*char');		fclose(fid);
	[names,x,y,COEs] = strread(c,'%s\t%f\t%f\t%f');
	delete(tmp_file);

	if (isempty(COEs))
		warndlg('This cruise doesn''t have the requested COEs in the crossings data base file.')
		return
	end
	if (~isa(names, 'cell')),	names = {names};	end

	Zmin = min(COEs);		Zmax = max(COEs);
	dZ = Zmax - Zmin;

	if (doHistogram)
		% Try to pick a resonable bin size
		if (abs(dZ) <= 20),			N = 20;		% bimW = 1
		elseif (abs(dZ) <= 50),		N = 25;		% binW = 2
		elseif (abs(dZ) <= 100),	N = 20;		% binW = 5
		elseif (abs(dZ) <= 500),	N = 50;		% binW = 10
		else						N = 100;	% Whatever
		end
		figure,histo_m('hist', COEs, N);
		return
	end
	
	cmap = jet(32);
	if (dZ == 0)			% Cte color
		zC = repmat(cmap(round(size(cmap,1)/2),:),numel(x),1);      % Midle color
	else            
		zC = round(((COEs - Zmin) / dZ) * (size(cmap,1)-1) + 1);
		zC = cmap(zC,:);
	end

	for (k = 1:numel(COEs))
		hS = line('XData',x(k), 'YData',y(k), 'parent',handles.axes1, 'Marker','o', ...
			'MarkerFaceColor',zC(k,:), 'MarkerEdgeColor','k', 'MarkerSize',8, 'Tag','COEpt');

        cmenuHand = uicontextmenu('Parent',handles.figure1);
        set(hS, 'UIContextMenu', cmenuHand);
        uimenu(cmenuHand, 'Label', sprintf('COE = %.1f (with %s)',COEs(k), names{k}));
        uimenu(cmenuHand, 'Label', 'Delete this', 'Call', {@uictxCOE,handles.axes1,hS,'del'}, 'Sep','on');
        uimenu(cmenuHand, 'Label', 'Delete all',  'Call', {@uictxCOE,handles.axes1,hS,'delAll'});
	end

function uictxCOE(obj,evt,hAx,h,tipo)
	if (strcmp(tipo,'del'))
		delete(h)
	elseif (strcmp(tipo,'delAll'))
		hAll = findobj(hAx,'Tag','COEpt');
		delete(hAll)
	end
% -----------------------------------------------------------------------------------------
	
	
% -----------------------------------------------------------------------------------------
function tracks = get_MGGtracks(obj, event, x, y)
% Get a list of MGG tracks that cross the rectangular area enclosed by gco
% OR by the rect defined by the optional input X, Y vars
%
% When an output is requested it will hold either a file name containing a list of the tracks
% or a cell array with the fullpath of each track.
% In either case, the tracks are not ploted.

	get_tracks_only = false;
	if (nargin == 2)
		h = gco;
		x = get(h,'XData');			y = get(h,'YData');
	end
	if (nargout)
		get_tracks_only = true;		tracks = [];
	end

	MIRONE_DIRS = getappdata(0,'MIRONE_DIRS');
	lim = sprintf('-R%.4f/%.4f/%.4f/%.4f',x(1),x(3),y(1:2));
	tmp_file = [MIRONE_DIRS.home_dir filesep 'tmp' filesep 'MGGtracks.txt'];
	if (ispc)
		dos(['x2sys_get -TMGD77+ -Fmtf1 -D -E ' lim ' > ' tmp_file]);
	else
		unix(['x2sys_get -TMGD77+ -Fmtf1 -D -E ' lim ' > ' tmp_file]);
	end

	d = dir(tmp_file);
	if (d.bytes < 5)
		warndlg('No MGG tracks found inside this rectangle.')
		delete(tmp_file);		return
	end

	fid = fopen(tmp_file,'r');	one_file = fgetl(fid);		fclose(fid);
	if (exist(one_file, 'file') == 2)	% We are probably in the same dir as the data files
		if (get_tracks_only),	tracks = ['list_' tmp_file];		return,		end
		mirone('GeophysicsImportGmtFile_CB',guidata(gcbo),['list_' tmp_file]),		return
	end

	X2SYS_HOME = getenv('X2SYS_HOME');
	if (isempty(X2SYS_HOME))
		warndlg('Tracks for the selected area exist, but I have no idea where they are. You need to learn about the $X2SYS_HOME env var')
		return
	end

	fid = fopen([X2SYS_HOME filesep 'MGD77+' filesep 'MGD77+_paths.txt'],'r');
	if (fid < 0)
		warndlg('Tracks for the selected area exist, but I have no idea where they are. Tell it to me via the $X2SYS_HOME MGD77+_paths.txt mechanism')
		return
	end

	c = fread(fid,inf,'*char');		fclose(fid);
	patos = strread(c,'%s','delimiter','\n');
	c = false(1, numel(patos));
	for (k = 1:numel(patos))
		if (patos{k}(1) == '#'),	c(k) = true;	end
	end
	patos(c) = [];

	% Do a final check that we are able to locate the track files and prepend file's path
	fid = fopen(tmp_file,'r');
	c = fread(fid,inf,'*char');		fclose(fid);
	tracks = strread(c,'%s','delimiter','\n');
	c = false(1, numel(tracks));
	for (k = 1:numel(tracks))
		for (m = 1:numel(patos))
			str = [patos{m} filesep tracks{k}];
			if (exist(str, 'file') == 2)
				tracks{k} = str;
				c(k) = true;	break
			end
		end
	end

	if (~any(c))
		warndlg('Tracks for the selected area exis, but I couldn''t find them. Not even with the help of the MGD77+_paths.txt file.')
	else
		if (get_tracks_only),	return,		end		% We already have the "tracks" content
		mirone('GeophysicsImportGmtFile_CB',guidata(gcbo), tracks)
	end
	delete(tmp_file);

% -----------------------------------------------------------------------------------------
function sow_GMT_DB_IDs(obj, event)
% See if there are GMT DB polygons inside the rectangle, and if yes display their IDs

	hRect = gco;		handles = guidata(hRect);

	hGMT_DB = findobj(handles.axes1, 'tag','GMT_DBpolyline');
	rx = get(hRect,'XData');		ry = get(hRect,'YData');
	rx = [min(rx) max(rx)];			ry = [min(ry) max(ry)];
	c = false(numel(hGMT_DB),1);	ID = cell(numel(hGMT_DB),1);
	for (m = 1:numel(hGMT_DB))			% Loop over objects to find out which cross the rectangle
		x = get(hGMT_DB(m),'XData');	y = get(hGMT_DB(m),'YData');
		if (any((x >= rx(1) & x <= rx(2)) & (y >= ry(1) & y <= ry(2))))
			str = getappdata(hGMT_DB(m), 'LineInfo');
			ind_IDi = strfind(str, 'Id = ');		ind_IDe = strfind(str, ' N =');
			ID{m} = str(ind_IDi+5:ind_IDe-1);
			c(m) = true;
		end
	end
	if (any(c))
		ID = ID(c);
		listbox_message(ID, 1, 'add')
	else
		msgbox('Nope, nothing arround here')
	end

% -----------------------------------------------------------------------------------------
function load_GMT_DB(obj, evt)
% Load one of the GMT DB ASCII files and trim contents such that only elements that
% totally or partially cross the rectangle area are displayed.

	hRect = gco;		handles = guidata(hRect);
	[FileName, PathName, handles] = put_or_get_file(handles, ...
		{'*.txt;*.TXT;*.dat;*DAT', 'Data files (*.txt,*.TXT,*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select File','get');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	rx = get(hRect,'XData');		ry = get(hRect,'YData');
	handles.ROI_rect = [min(rx) max(rx) min(ry) max(ry)];	% Signal load_xyz() that we want to use this ROI
	load_xyz(handles, fname)

% -----------------------------------------------------------------------------------------
function assoc_gmt_symbol(obj, evt, h)
% Associate a GMT custom symbol file to the line object that called this function

	if (nargin == 2)	% Otherwise this is a call issued in FileOpenSession_CB and we only want to set the BDF CB
		h = gco;
		str = {'*.def;*.eps','GMT custom symbol (*.def,*.eps)'; '*.*', 'All Files (*.*)'};
		[FileName,PathName] = put_or_get_file(guidata(gco),str,'Select GMT custom symbol','get');
		if isequal(FileName,0),		return,		end
		setappdata(h, 'cust_symb', [PathName FileName])	% OK, now it's up to write_gmt_script to use this info
		set(h, 'LineStyle', '--')		% Use this style to indicate the different nature of this element
	end
	setappdata(h, 'BD_runCB', {'deal_opts', [], 'show_symbName',[], [], h})	% First [] is the 'opt' (not used) arg in this main

% -----------------------------------------------------------------------------------------
function show_symbName(obj, evt, h)
% Print the associated symbol file name for a period of one second
	hAx = get(h, 'Parent');
	fname = getappdata(h, 'cust_symb');
	if (isempty(fname)),	return,		end
	pos = get(hAx, 'CurrentPoint');		pos = pos(1,:);
	hTxt = text(pos(1),pos(2), fname, 'Interpreter', 'none');
	pause(1)
	delete(hTxt)
	refresh(get(hAx, 'Parent'))		% Need this because of the bloody compiled version bugs.

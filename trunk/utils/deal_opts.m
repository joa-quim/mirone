function out = deal_opts(opt, opt2)
% Function to provide fine (manual) control on some Mirone modules.
%
% The mechanism of this function relies on the contents of the "OPTcontrol.txt" control file.
% That file has keewords that trigger the usage of optional (specialized) features in Mirone.
%
% OPT directs to the option of interest (growing number)
% OPT2 is a uimenu handle for cases where the OPT case needs one

	MIRONE_DIRS = getappdata(0,'MIRONE_DIRS');
	opt_file = [MIRONE_DIRS.home_dir filesep 'data' filesep 'OPTcontrol.txt'];
	out = [];

	if ( ~exist(opt_file, 'file') == 2 ),	return,		end
	fid = fopen(opt_file, 'r');
	c = (fread(fid,'*char'))';      fclose(fid);
	lines = strread(c,'%s','delimiter','\n');   clear c fid;
	m = numel(lines);
	c = true(1, m);

	for (k = 1:m)
		if (strncmp(lines{k},'MIR_',4)),	c(k) = false;	end
	end
	lines(c) = [];			% Delete non-keyword (comment) lines
	if (isempty(lines))		return,		end

	switch lower(opt)
		case 'gmtedit'		% To select what to plot in GMTEDIT slots
			for (k = 1:numel(lines))
				if (strncmp(lines{k}(5:end),'GMTEDIT',7))
					t = strtok(lines{k}(13:end));
					if ( strcmp(t(1:2), '-V') )		% Here we only check for a -V... and do not check for errors
						out = t;
					end
					break
				end
			end
		case 'mgg'			% To find tracks inside calling rectangle
			for (k = 1:numel(lines))
				if (strncmp(lines{k}(5:end),'MGG',3))
					t = strtok(lines{k}(9:end));
					if (t == '1'),	out = true;		end
					break
				end
			end
			%if (out),	custom_rect(opt2),	end		% Create a uimenu associated to a rectangle
			if (out)
				c = uimenu(opt2, 'Label', 'Custom','Sep','on');		% Create a uimenu associated to a rectangle
				uimenu(c, 'Label', 'MGG tracks', 'Call',  @get_MGGtracks);
			end

		case 'mgg_coe'		% To ...
			for (k = 1:numel(lines))
				if (strncmp(lines{k}(5:end),'COEs',4))
					out = lines{k}(10:end);
					break
				end
			end
			if (out)
				coeVar = 'mag';			% Default
				[t, r] = strtok(out);
				if (~isempty(r))		coeVar = ddewhite(r);	end
				uimenu(opt2, 'Label', 'Show COEs', 'Call', {@get_COEs, t, coeVar},'Sep','on');
			end
	end

% -----------------------------------------------------------------------------------------
function get_COEs(obj, event, coeFile, coeVar)
% Get a list of MGG COEs involving gco and plot them with uicontexts

	h = gco;		handles = guidata(h);

	lims = getappdata(handles.axes1, 'ThisImageLims');
	opt_R = sprintf('-R%.4f/%.4f/%.4f/%.4f',lims(1:4));
	MIRONE_DIRS = getappdata(0,'MIRONE_DIRS');
	tmp_file = [MIRONE_DIRS.home_dir filesep 'tmp' filesep 'MGGtmp.txt'];
	fname = getappdata(h, 'FullName');
	[pato, fname] = fileparts(fname);

	if (ispc)
		dos( ['x2sys_list ' coeFile ' -TMGD77+ -Fnxyc -C' coeVar ' -S' fname ' ' opt_R ' > ' tmp_file]);
	else
		unix(['x2sys_list ' coeFile ' -TMGD77+ -Fnxyc -C' coeVar ' -S' fname ' ' opt_R ' > ' tmp_file]);
	end
	
	d = dir(tmp_file);
	if (d.bytes < 5)
		warndlg('This cruise doesn''t have any COEs in the crossings data base file.')
		delete(tmp_file);		return
	end

	fid = fopen(tmp_file);
	fgetl(fid);		fgetl(fid);		fgetl(fid);		% Jum the 3 header lines
	c = fread(fid,inf,'*char');		fclose(fid);
	[names,x,y,COEs] = strread(c,'%s\t%f\t%f\t%f');
	if (~isa(names, 'cell'))	names = {names};	end
	delete(tmp_file);

	Zmin = min(COEs);		Zmax = max(COEs);
	dZ = Zmax - Zmin; 
	cmap = jet(32);
	if (dZ == 0)			% Cte color
		zC = repmat(cmap(round(size(cmap,1)/2),:),nPts,1);      % Midle color
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
function get_MGGtracks(obj, event)
% Get a list of MGG tracks that cross the rectangular area enclosed by gco
	h = gco;
	x = get(h,'XData');     y = get(h,'YData');

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
		mirone('GeophysicsImportGmtFile_CB',guidata(gcbo), tracks)
	end
	delete(tmp_file);

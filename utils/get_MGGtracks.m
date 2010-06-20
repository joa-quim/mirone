function get_MGGtracks(obj, event, x, y)
% Get a list of MGG tracks that cross the rectangular area enclosed by gco
% OR byr the rect defined by the optional input X, Y vars

	if (nargin == 2)
		h = gco;
		x = get(h,'XData');     y = get(h,'YData');
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

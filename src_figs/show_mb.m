function varargout = show_mb(varargin)
% Helper window to decide what to do when a multi-beam file, or datalist.mb-1, is loaded

%	Copyright (c) 2004-2018 by J. Luis
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

% $Id: show_MB.m 11331 2018-06-22 18:52:39Z j $

	if (nargin > 1 && ischar(varargin{1}))
		gui_CB = str2func(varargin{1});
		[varargout{1:nargout}] = feval(gui_CB,varargin{2:end});
	else
		h = show_mb_OF(varargin{:});
		if (nargout),	varargout{1} = h;   end
	end

% ---------------------------------------------------------------------------------
function hObject = show_mb_OF(varargin)

	hObject = show_mb_LayoutFcn;
	handles = guihandles(hObject);
	move2side(hObject,'right');

	handMir = varargin{1};
	handles.fnameMB = varargin{2};
	handles.whichFleder = handMir.whichFleder;
	handles.TDRver   = handMir.TDRver;
	handles.path_tmp = handMir.path_tmp;
	handles.no_file  = handMir.no_file;
	handles.hMirFig  = handMir.figure1;
	handles.IamCompiled = handMir.IamCompiled;

	handles.opt_N = '';		handles.opt_C = '';		handles.opt_R = '';		handles.opt_Z = '';	% Defaults that might be used later
	handles.first_datalist_scan = true;
	handles.show_datalist_check = false;
	handles.did_autoclean = false;
	handles.list_files = '';		% Will store the contents of datalists files when those are used.
	handles.tol = 0.02;				% For the autocleaner step.
	handles.N_ITER = 3;
	handles.grid_inc = [];			% If not set, it will be estimated from point spread.
	handles.area_cmd = '';			% If the "Params" fig is called this may hold the mbareaclen command

	if (~exist([handles.fnameMB '.esf'], 'file'))		% No .esf's, no option to hide them
		set(handles.check_showFlagged, 'Vis', 'off')
	end

	%------------ Give a Pro look (3D) to the frame boxes  ------------------------------------
	new_frame3D(hObject, [handles.text_PC handles.text_PI])
	%------------- END Pro look (3D) ----------------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');

	if (nargin > 2),	external_drive(handles, 'show_mb', varargin{3:end}),	end

% ----------------------------------------------------------------------
function edit_symbSize_CB(hObject, handles)
% Just check that impossible figures are not acepted
	str = get(hObject,'String');        s = str2double(str);
	if (isnan(s) || s <= 0),	set(hObject,'String','0.01'),	return,		end

% ----------------------------------------------------------------------
function [out, BB, DS] = push_showPC_CB(hObject, handles, opt)
% Show data as Point Cloud. Flagged points go to a second layer
% OPT, is present, must be either '-C' or '-Ca' to fetch SideScan or Amplitude instead of Bathymetry
% If argout, OUT is the data as Mx3 (no flagged vals), BB is the BoundingBox and DS the point spacing

	if (nargin == 2),	opt = [];	end
	contents = cellstr(get(handles.popup_symb,'String'));
	val = get(handles.popup_symb,'Value');
	switch contents{val}
		case 'Point',    symb = 7;
		case 'Cube',     symb = 3;
		case 'Sphere',   symb = 6;
		case 'Dyamond',  symb = 4;
		case 'Cylinder', symb = 5;
	end
	fname = [handles.path_tmp 'lixo.sd'];
	siz = str2double(get(handles.edit_symbSize,'String'));
	PTparams.Symbol = symb;			PTparams.PointRad = siz;
	par.TDRver = handles.TDRver;	par.proj = 'geog';		par.PTparams = PTparams;

	if (exist([handles.fnameMB '.inf'], 'file') ~= 2)		% If .inf file does not exist, create one
		D = gmtmex(sprintf('mbinfo -G -I%s -O', handles.fnameMB));
		fid = fopen([handles.fnameMB '.inf'],'wt');
		for (k = 1:numel(D.text))
			fprintf(fid, '%s\n', D.text{k});
		end
		fclose(fid);
	end

	set(handles.figure1,'Pointer','watch');
	if (get(handles.check_showFlagged, 'Val'))		% Show all flagged beans
		D = gmtmex(sprintf('mbgetdata -I%s -A100000 %s', handles.fnameMB, opt));
	else
		D = gmtmex(sprintf('mbgetdata -I%s -A-100000 %s', handles.fnameMB, opt));
	end

	% Search gor garbage around the 0,0 point. If finds some, just delte it
	central = round(size(D(1).data,2) / 2);		% Find the central beam
	x = D(1).data(:,central);		y = D(2).data(:,central);		z = D(3).data(:,central);
	ind = x < 0.01 & x > -0.01 & y < 0.01 & y > -0.01;
	if (~strncmp(opt, '-C', 2)),	ind = ind & z > -2000;	end		% This is the bathymetry case
	if (any(ind))
		D(1).data(ind,:) = [];	D(2).data(ind,:) = [];	D(3).data(ind,:) = [];
	end
	% Here we only want to check for flagged guys in the bathymetry case
	if (strncmp(opt, '-C', 2)),	ind = true;		% SS or Amplitude. This trick makes it ignore the flagging issue.
	else,						ind = (D(3).data < 50000);	% This gives us the indices of all non-flagged guys
	end

	if (~all(ind(:)))						% If we have some flagged
		xyz1 = [D(1).data(ind) D(2).data(ind) D(3).data(ind)];
		ind  = ~ind;						% Now ind holds indices of the flagged guys
		xyz2 = [];
		if (~get(handles.check_showFlagged, 'Val'))		% Do not show flagged beans
			xyz2 = [D(1).data(ind) D(2).data(ind) D(3).data(ind)-100000];
			xyz2(xyz2(:,3) == 0,:) = [];	% Do not show the z = 0 pts
		end
		t = [min(xyz1) max(xyz1)];		bb1 = [t(1) t(4) t(2) t(5) t(3) t(6)];	% [xmin xmax ymin ymax zmin zmax]
		if (~isempty(xyz2))					% These are the flagged
			t = [min(xyz2) max(xyz2)];		bb2 = [t(1) t(4) t(2) t(5) t(3) t(6)];
			bbg  = [min(bb1(1),bb2(1)) max(bb1(2),bb2(2)) min(bb1(3),bb2(3)) max(bb1(4),bb2(4)) ...
					min(bb1(5),bb2(5)) max(bb1(6),bb2(6))];
		else
			bbg  = bb1(1:6);
		end

		if (nargout)				% Than we are done, just assign output vars and return
			set(handles.figure1,'Pointer','arrow');
			out = xyz1;		BB = bbg;
			[dx, dy] = get_typical_pt_dist(D);
			DS = sqrt(dx*dx + dy*dy);
			return
		end

		fid  = write_flederFiles('scene_pts', fname, [], 'begin', bbg, par);
		write_flederFiles('scene_pts', fid, xyz1, 'sec', [bb1 bbg], par);
		if (~isempty(xyz2))
			PTparams.Symbol = symb;	PTparams.ColorBy = 0;	par.PTparams = PTparams;
			write_flederFiles('scene_pts', fid, xyz2, 'end', [bb2 bbg], par);
		else
			% Here I must cheat. Apparently .scene with only one set makes it screem that can't read something. So give it something
			write_flederFiles('scene_pts', fid, xyz1(1:4,:), 'end', [bbg bbg], par);	% Need to send in at least 4 pts
		end
	else							% No flagged pts
		xyz1 = [D(1).data(:) D(2).data(:) D(3).data(:)];
		if (strcmp(opt, '-C'))		% Don't know why but SideScan may have many padding NaNs. Remove them
			ind = isnan(xyz1(:,3));
			if (any(ind)),	xyz1(ind, :) = [];	end
		end
		t = [min(xyz1) max(xyz1)];		bbg = [t(1) t(4) t(2) t(5) t(3) t(6)];	% [xmin xmax ymin ymax zmin zmax]

		if (nargout)				% Than we are done, just assign output vars and return
			set(handles.figure1,'Pointer','arrow');
			out = xyz1;		BB = bbg;
			[dx, dy] = get_typical_pt_dist(D);
			DS = sqrt(dx*dx + dy*dy);
			return
		end

		fid  = write_flederFiles('scene_pts', fname, [], 'begin', bbg, par);
		write_flederFiles('scene_pts', fid, xyz1, 'sec', [bbg bbg], par);
		write_flederFiles('scene_pts', fid, xyz1(1:4,:), 'end', [bbg bbg], par);	% Cheat second dataset. Need at least 4 pts
	end
	fname = [fname '.scene'];			% Because name was changed in write_flederFiles()
	comm  = [' -scene ' fname ' &'];	% A SCENE file
	set(handles.figure1,'Pointer','arrow');
	show_fleder(handles, comm)

% ----------------------------------------------------------------------
function push_autoclean_CB(hObject, handles)
% Do an automatic cleaning and optionally show a point cloud with cleaned and flagged
	tol = handles.tol;			N_ITER = handles.N_ITER;

	contents = cellstr(get(handles.popup_symb,'String'));
	val = get(handles.popup_symb,'Value');
	switch contents{val}
		case 'Point',    symb = 7;
		case 'Cube',     symb = 3;
		case 'Sphere',   symb = 6;
		case 'Dyamond',  symb = 4;
		case 'Cylinder', symb = 5;
	end
	fname = [handles.path_tmp 'lixo.sd'];
	siz = str2double(get(handles.edit_symbSize,'String'));
	PTparams.Symbol = symb;			PTparams.PointRad = siz;
	par.TDRver = handles.TDRver;	par.proj = 'geog';		par.PTparams = PTparams;

	% Check if we are dealing with a single file or a datalist file
	datalist = '';
	[pato, fila, ext] = fileparts(handles.fnameMB);
	if (strcmpi(ext, '.mb-1'))
		datalist = handles.fnameMB;
	end

	if (~get(handles.check_showCleaneds, 'Val'))	% No Viz, only compute and exit
		autocleaner(handles, tol, N_ITER, datalist);
		set(handles.push_applyClean, 'Enable', 'on')
		jumping_done
	else
		[xyz1,xyz2] = autocleaner(handles, tol, N_ITER, datalist);
		if(isempty(xyz1) && isempty(xyz2)),		return,		end			% Some error occured
		xyz2(xyz2(:,3) == 0,:) = [];		% Do not show the z = 0 pts
		bb1 = [min(xyz1)' max(xyz1)']';		bb1 = bb1(:)';
		if (~isempty(xyz2))
			bb2 = [min(xyz2)' max(xyz2)']';		bb2 = bb2(:)';
			bbg = [min(bb1(1),bb2(1)) max(bb1(2),bb2(2)) min(bb1(3),bb2(3)) max(bb1(4),bb2(4)) ...
				   min(bb1(5),bb2(5)) max(bb1(6),bb2(6))];
		else
			bbg  = bb1(1:6);
		end

		fid = write_flederFiles('scene_pts', fname, [], 'begin', bbg, par);
		write_flederFiles('scene_pts', fid, xyz1, 'sec', [bb1 bbg], par);
		if (~isempty(xyz2))
			PTparams.Symbol = symb;	PTparams.ColorBy = 0;	par.PTparams = PTparams;
			write_flederFiles('scene_pts', fid, xyz2, 'end', [bb2 bbg], par);
		end
		fname = [fname '.scene'];			% Because name was changed in write_flederFiles()
		comm = [' -scene ' fname ' &'];		% A SCENE file

		handles.did_autoclean = true;		% So we know that we can update the corresponding .esf
		set(handles.push_applyClean, 'Enable', 'on')
		show_fleder(handles, comm)
	end

% ----------------------------------------------------------------------
function show_fleder(handles, comm)
% Show nn Fleder. Free viewer or the other.
	if (handles.whichFleder),	fcomm = ['iview4d' comm];		% Free viewer
	else,						fcomm = ['fledermaus' comm];	% The real thing
	end
	try
		if (isunix)				% Stupid linux doesn't react to a non-existant iview4d
			resp = unix(fcomm);
			if (resp == 0)
				errordlg('I could not find Fledermaus. Hmmm, do you have it?','Error')
			end
		elseif (ispc)
			if (strcmp(handles.TDRver, '2.0')),	fcomm(6) = '3';		end		% Call iview3d
			if (handles.no_file && ishandle(handles.hMirFig)),	delete(handles.hMirFig),	end
			if (handles.IamCompiled),	fcomm = ['start /b ' fcomm];	end
			dos(fcomm);		% s is always 0 (success) as long as iview4d is accessible, even when the comm fails
		else
			errordlg('Unknown platform.','Error'),	return
		end
	catch
		errordlg('I could not find Fledermaus. Hmmm, do you have it?','Error')
	end

%-------------------------------------------------------------------------------------
function push_applyClean_CB(hObject, handles)
% Apply the cleanings stored in the .mask file created previously.
	try
		% Have to use mbset to tell the .par files where esf and processed files should be written
		if (isempty(handles.list_files))		% Single file
			[pato,outfile,ext] = fileparts(handles.fnameMB);
			outfile = [pato filesep outfile 'p' ext];
			cmd = ['mbset -I' handles.fnameMB ' -POUTFILE:' outfile ' -PEDITSAVEFILE:' handles.fnameMB '.esf'];
			gmtmex(cmd)
			gmtmex(['mbflags -I' handles.fnameMB ' -E' handles.fnameMB '.mask'])
		else									% datalist file
			for (k = 1:numel(handles.list_files))
				this_file = handles.list_files{k};
				[pato,outfile_,ext] = fileparts(this_file);
				outfile = [pato filesep outfile_ 'p' ext]; 
				cmd = ['mbset -I' this_file ' -POUTFILE:' outfile ' -PEDITSAVEFILE:' this_file '.esf'];
				gmtmex(cmd)
				gmtmex(['mbflags -I' this_file ' -E' this_file '.mask'])
			end
		end
		jumping_done
	catch
		errordlg(lasterr, 'Error')
	end

% ----------------------------------------------------------------------
function push_params_CB(hObject, handles)
% Will call a new window to help with the autocleaning choices
	% Check if we are dealing with a single file or a datalist file
	[pato, fila, ext] = fileparts(handles.fnameMB);
	if (strcmpi(ext, '.mb-1'))
		list_files = sanitize_datalist(handles);	% Read and sanitize a datalist file. Returns cell array with full fnames
		if (numel(list_files) == 0)
			errordlg('This datalist file is empty. Bye Bye.','Error')
			return
		end
	else
		list_files = {handles.fnameMB};
	end

	info = find_BBs(list_files);		% Get the BBs of each file
	if (size(info.BB,1) == 1)
		global_BB = info.BB;
	else
		mins = min(info.BB);	maxs = max(info.BB);
		global_BB = [mins(1) maxs(2) mins(3) maxs(4) mins(5) maxs(6)];
	end
	
	% ------------------------------ OK, ready to call the Params Fig -------------------------------
	out = mb_cleaning_params(global_BB);
	if (isempty(out) || (~out.gridclean && ~out.areaclean)),	return,		end		% Neither method was selected

	if (out.gridclean)
		handles.grid_inc = out.grid_inc;
		handles.N_ITER = out.grid_n_iter;
		handles.tol = out.grid_tol / (0.7^out.grid_n_iter);		% 0.7 because that's the scaling fact we use in the iterations
		handles.area_cmd = '';		% Make sure we wont run both
	else
		if (isnan(out.opt_S)),		return,		end		% An errors in binsize
		binsize = out.opt_S * 111000;		% Convert back to aproximate meters because that's what MB wants
		cmd = sprintf('mbareaclean %s -S%f', out.opt_R, binsize);
		if (out.opt_D_active)
			cmd = sprintf('%s -D%.6g/%d', cmd, out.opt_D_threshold, out.opt_D_nmin);
		end
		if (out.opt_M_active)
			cmd = sprintf('%s -M%.6g/%d', cmd, out.opt_M_threshold, out.opt_M_nmin);
		end
		if (out.opt_N_active)
			if (out.opt_N_Rev),		out.opt_N_min = -out.opt_N_min;		end
			cmd = sprintf('%s -N%d', cmd, out.opt_N_min);
		end
		if (out.opt_G_active)
			cmd = [cmd ' -G'];
		end
		cmd = sprintf('%s -I%s', cmd, handles.fnameMB);
		handles.area_cmd = cmd;
	end
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------
function push_help_CB(hObject, handles)
	t = sprintf(['Display contents of a MB- readable data file or, depending on the option choosed, ' ...
		'of a datalist.mb-1 file. A first and VERY IMPORTANT point to be aware off is that the ' ...
		'companion .inf file MUST exist for each of the files beeing displayed.\n\n' ...
		'The first group of options let you see the data as a point cloud in a 3D display. Here, you ' ...
		'can use a datalist.mb-1 to show a set of files or an individual file and ask the program to ' ...
		'do an automatic cleaning based on the distance from the points to a grid computed with all ' ...
		'points of file. The result of this exercise is shown in the point cloud display and saved ' ...
		'in a file with the original file name appended with the ''.mask'' extension. This file can be ' ...
		'used with module ''mbflags'' to update the flags file ''.esf''.\n\n' ...
		'The second group let us have a quick view as an image of either a single file or the contents ' ...
		'of a datalist.mb-1 file. This file may additionally have a first commnet line (a line starting ' ...
		'with the # character) with the options -Z, -C and -N of the ''mbswath'' program to control what ' ...
		'is displayed. In such cases the ''Use commands in datalist'' must be checked.']);
	helpdlg(t, 'Help on show MB data')

%-------------------------------------------------------------------------------------
function [xyz, xyzK] = autocleaner(handles, tol, N_ITER, datalist)
% ...
	% If we have a hanging mbareaclean command, run it and return
	if (~isempty(handles.area_cmd))
		gmtmex(handles.area_cmd)			% Run the mbareaclean command
		handles.area_cmd = '';				% Reset to not run it again without uses further request
		set(handles.push_applyClean, 'Enable', 'off')
		guidata(handles.figure1, handles)
		if (get(handles.check_showCleaneds, 'Val'))
			push_showPC_CB([], handles)		% Show the result
		end
		return
	end

	ds = inf;								% Initialize point sapcing
	if (isempty(datalist))					% Single file, simplest case
		D = gmtmex(sprintf('mbgetdata -I%s -A100000', handles.fnameMB));	% -A is to allow fishing the flagged pts
		x = D(1).data(:);		y = D(2).data(:);		z = D(3).data(:);
		old_flags = (z > 50000);			% Need this to be able to restore the true Z's of old flagged PTs
		x(old_flags) = NaN;		y(old_flags) = NaN;
		z(old_flags) = NaN;					% We don't want to use the old flagged
		ind = x < 0.1 & x > -0.1 & y < 0.1 & y > -0.1 & z > -2000;
		if (any(ind))
			x(ind) = NaN;	y(ind) = NaN;	z(ind) = NaN;
		end
		% Compute increment as 3 times the typical point spacing fetch from data.
		[dx, dy] = get_typical_pt_dist(D);
% 		mid_row = round(size(D(2).data, 1) / 2);
% 		dy = abs(median(diff(D(2).data(mid_row,:))));
% 		dx = abs(median(diff(D(1).data(mid_row,:)))) * cos(D(2).data(mid_row,1) * pi/180);
% 		if (dx == 0 || dy == 0)		% Try at another location infile
% 			mid_row = round(size(D(2).data, 1) / 3);
% 			dy = abs(median(diff(D(2).data(mid_row,:))));
% 			dx = abs(median(diff(D(1).data(mid_row,:)))) * cos(D(2).data(mid_row,1) * pi/180);
% 		end
		% If it's still 0, give up
		if (dx == 0 || dy == 0)
			txt = ['This file ' handles.fnameMB ' is really screwed. Can''t guess a decent point spread distance ' ...
					'You will have to use the "Params" option and set a your guess of the bin size, or otherwise ' ...
					'this option can not be used with this file'];
			errordlg(txt, 'Error')
			xyz = [];	xyzK = [];
			return
		end
		x_min = min(x);		x_max = max(x);		y_min = min(y);		y_max = max(y);
		ds = sqrt(dx*dx + dy*dy) * 3;
		nx = round((x_max - x_min) / ds);		ny = round((y_max - y_min) / ds);
		if (nx > 5000 || ny > 5000)
			t = sprintf(['The estimated grid size is unresonable (%d x %d). This is normally due to too bad ' ...
				'coodinates in file\n\nQUITING'],nx,ny);
			errordlg(t,'Error');
			xyz = [];	xyzK = [];
			return
		end
		opt_I = sprintf('-I%f', ds);
		opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g', x_min, x_max, y_min, y_max);
		for (k = 1:N_ITER)
			[x,y,z, ind] = iteration_cleaner(x,y,z, opt_R, opt_I, tol, ind);	% return 'cleaned'
			tol = tol * 0.7;
		end
		save_flags(handles, handles.fnameMB, ind, size(D(1).data))	% Save the new flags in a .mask file

		if (nargout > 0 || true)
			% Now it's time to remove the NaNs that were inserted in place of the cleaned pts
			indNaN = isnan(z);
			x(indNaN) = [];		y(indNaN) = [];		z(indNaN) = [];
			xyz = [x y z];
			D(3).data(old_flags) = D(3).data(old_flags) - 100000;	% Restore the true value of the old flagged pts
			ind = (ind | old_flags);								% Join all flagged
			xyzK = [D(1).data(ind) D(2).data(ind) D(3).data(ind)];	% The Killed points
		end

	else
		if (nargout > 0 || true)
			[xc,yc,zc, xk,yk,zk, ds] = autocleaner_list(handles, tol, N_ITER);
			xyz  = [xc yc zc];	clear xc yc zc		% AWFULL MEMORY WASTE
			xyzK = [xk yk zk];	clear xk yk zk
		else
			autocleaner_list(handles, tol, N_ITER);
		end
	end

	patoMother = fileparts(handles.fnameMB);		% Path of the datlist file
	fid = fopen([patoMother filesep 'data_cleaned.bin'],'wb');
	if (ds < 1)					% Save in doubles
		fwrite(fid, xyz', 'real*8');
	else						% Save in singles
		xyz = single(xyz)';
		fwrite(fid, xyz', 'real*4');
	end
	fclose(fid);

% 	[Z, head] = gmtmbgrid_m(xyz(:,:,1), xyz(:,:,2), xyz(:,:,3), '-I0.0005', opt_R, '-Mz', '-C2');
% 	opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g', min(xc),max(xc),min(yc),max(yc));
% 	[Z, head] = gmtmbgrid_m(xc, yc, zc, '-I0.0005', opt_R, '-Mz', '-C2');
% 	aux_funs('showgrd', struct('geog',1), Z, head, 'Autocleaned')

%-------------------------------------------------------------------------------------
function [dx, dy] = get_typical_pt_dist(D)
% Estimate the typical point spacing fro data mid size

	mid_row = round(size(D(2).data, 1) / 2);
	dy = abs(median(diff(D(2).data(mid_row,:))));
	dx = abs(median(diff(D(1).data(mid_row,:)))) * cos(D(2).data(mid_row,1) * pi/180);
	if (dx == 0 || dy == 0)		% Try at another locations in file
		for (k = 1:3)
			mid_row = round(size(D(2).data, 1) * k / 4);
			dy = abs(median(diff(D(2).data(mid_row,:))));
			dx = abs(median(diff(D(1).data(mid_row,:)))) * cos(D(2).data(mid_row,1) * pi/180);
			if (dx ~= 0 && dy ~= 0),	break,	end
		end
	end

%-------------------------------------------------------------------------------------
function [x,y,z, ind] = iteration_cleaner(x,y,z, opt_R, opt_I, tol, ind_old)
% Do one iteration round and return only the points that are within the condition.
% TOL is a percentage of the whater depth.
	[Z, head] = gmtmbgrid_m(x, y, z, opt_I, opt_R, '-Mz', '-C1');
	zz = grdtrack_m(Z,head,[x y],'-Z')';
	tol = abs(z * tol);				% From here on, 'tol' is a vector with size = numel(z)
	ind = abs(z - zz(:)) > tol;		% Indices of those points that are no further than TOL from surface
	clear tol zz
	x(ind) = NaN;		y(ind) = NaN;	z(ind) = NaN;
	if (~isempty(ind_old))
		ind = ind | ind_old;		% Add this ieration result to the result of previous iterations
	end

%-------------------------------------------------------------------------------------
function todos = sanitize_datalist(handles)
% Read a datalist file and verify that:
% 1. NO BLANKS
% 2. If format number is present, remove it
% 3. If files have no path, prepend the path of the datalist itself
% 4. TODO. Maybe one day try to accept blanks in names when honestly embeded in ""

	fid = fopen(handles.fnameMB,'rt');
	todos = fread(fid,'*char');		fclose(fid);
	todos = strread(todos','%s','delimiter','\n');
	n_files = numel(todos);

	patoMother = fileparts(handles.fnameMB);			% Path of the datlist file
	c = false(n_files,1);
	msg = '';
	for (k = 1:n_files)
		if (todos{k}(1) == '#' || todos{k}(1) == '>')	% Remove comment lines
			c(k) = true;
		else
			todos{k} = strtok(todos{k});		% Remove eventual format number at the end
			[pato, fname, ext] = fileparts(todos{k});
			if (~strncmpi(ext, '.mb', 3))
				msg = ['Your files are probably violating the golden rule: NO SPACES IN NAMES ' ...
					   'or have not a .mbXX extension'];
				c(k) = true;
			elseif (isempty(pato))
				todos{k} = [patoMother filesep fname ext];
			end
		end
	end
	todos(c) = [];

	if (numel(todos) == 0)
		errordlg('The dalist file is empty or full of errors','Error')
	elseif (numel(todos) == 1)
		errordlg('Cannot process datalist files that have only one file.','Error')
	end
	if (~isempty(msg))
		errordlg(msg, 'Error')
	end

%-------------------------------------------------------------------------------------
function [xc,yc,zc, xk,yk,zk, ds] = autocleaner_list(handles, tol, N_ITER)
% Scan all files in a datalist file and clean each file individually against the data of
% all files. This is awfully inefficient since we should work only over sub-regions. Future work.

	list_files = sanitize_datalist(handles);	% Read and sanitize a datalist file. Returns cell array with full fnames
	handles.list_files = list_files;
	guidata(handles.figure1, handles)			% Save this for use in "Apply cleanings"

	ds = inf;
	n_files = numel(list_files);
	if (n_files <= 1)
		xc=[];	yc=[];	zc=[];	xk=[];	yk=[];	zk=[];
		return
	end

	aguentabar('title',sprintf('cleaning %d files',n_files))

	info = find_BBs(list_files);
	fname_tmp = [handles.path_tmp 'datalist.mb-1'];
	for (k = 1:n_files)							% Loop over number of files in datalist
		overlaps = find_overlaps(list_files, k, info);	% Get indices of files that overlap with current one
		if (isempty(overlaps))					% Shit, this doesn't overlap with any other file
			% Call the single file case
			if (nargout > 0)
				[xyz1,xyzK] = autocleaner(handles, tol, N_ITER, '');
				if (isempty(xyz1)),	continue,	end
				if (k == 1)								% First time, create the output arrays
					xc = xyz1(:,1);		yc = xyz1(:,2);		zc = xyz1(:,3);		clear xyz1
					xk = xyzK(:,1);		yk = xyzK(:,2);		zk = xyzK(:,3);		clear xyzK
				else
					xc = [xc; xyz1(:,1)];	yc = [yc; xyz1(:,2)];	zc = [zc; xyz1(:,3)];	clear xyz1 %#ok<AGROW>
					xk = [xk; xyzK(:,1)];	yk = [yk; xyzK(:,2)];	zk = [zk; xyzK(:,3)];	clear xyzK %#ok<AGROW>
				end
			else
				autocleaner(handles, tol, N_ITER, '');	% Do the cleaning silently
			end
			aguentabar(k/n_files)
			continue
		end

		lista_t = list_files(overlaps);
		fid = fopen(fname_tmp, 'wt');			% Open a TMP file to store the names of the 'other' files
		for (n = 1:numel(lista_t))
			fprintf(fid, '%s\n', lista_t{n});	% Write the TMP datalist file with all but current (loop's) file
		end
		fclose(fid);

		Dcurr   = gmtmex(['mbgetdata -I' list_files{k} ' -A100000']);	% -A1e5 is to allow recovering true values
		Dothers = gmtmex(['mbgetdata -I' fname_tmp ' -A']);		% These ones I know that, flagged => NaN
		if (k == 1)								% Compute working inc as 3 times typical pt spacings
			mid_row = round(size(Dcurr(2).data, 1) / 2);
			dy = abs(median(diff(Dcurr(2).data(mid_row,:))));
			dx = abs(median(diff(Dcurr(1).data(mid_row,:)))) * cos(Dcurr(2).data(mid_row,1) * pi/180);
			% If it's still 0, give up
			if (dx == 0 || dy == 0)
				txt = ['This file ' handles.fnameMB ' is really screwed. Can''t guess a decent point spread distance ' ...
					'You will have to use the "Params" option and set a your guess of the bin size, or otherwise ' ...
					'this option can not be used with this file'];
				errordlg(txt, 'Error')
				xc=[];	yc=[];	zc=[];	xk=[];	yk=[];	zk=[];
				return
			end
			ds = sqrt(dx*dx + dy*dy) * 3;
			opt_I = sprintf('-I%f', ds);
		end
		[x,y,z, ind, old_flags] = iteration_cleaner_list(Dcurr(1).data(:),Dcurr(2).data(:),Dcurr(3).data(:), ...
			Dothers(1).data(:),Dothers(2).data(:),Dothers(3).data(:), opt_I, tol, N_ITER);	% X,Y,Z are the cleaned PTs

		save_flags(handles, list_files{k}, ind, size(Dcurr(1).data))	% Time to save flags

		if (nargout > 0)
			Dcurr(3).data(old_flags) = Dcurr(3).data(old_flags) - 100000;	% Restore the true value of the old flagged PTs
			ind = (ind | old_flags);				% Add the old and new flags
			if (k == 1)								% First time, create the output arrays
				xc = x;		yc = y;		zc = z;
				xk = Dcurr(1).data(ind);	yk = Dcurr(2).data(ind);	zk = Dcurr(3).data(ind);
			else
				xc = [xc; x];	yc = [yc; y];		zc = [zc; z]; %#ok<AGROW>
				xk = [xk; Dcurr(1).data(ind)];		yk = [yk; Dcurr(2).data(ind)];	zk = [zk; Dcurr(3).data(ind)]; %#ok<AGROW>
			end
		end
		aguentabar(k/n_files)
	end

%-------------------------------------------------------------------------------------
function [xc,yc,zc, the_ind, old_flags] = iteration_cleaner_list(xc,yc,zc, xo,yo,zo, opt_I, tol, N_ITER)
% XC,... coordinates of current file beeing cleaned
% XO,... coordinates of all other points from the remaining files in the datalist file (all but current)
% TOL is a scalar with the satrting percentage of the water depth. Points further than this from the
% surface are flagged. Furthermore, this process is run N_ITER times where, in eache ietr, TOL is shrink by 0.7
%
% Also note that the flagged PTs of ZC were added 100000. This allows us to find and restore their original values

	old_flags = (zc > 50000);
	zc(old_flags) = NaN;					% Since they are flagged we don't want to use them
	x = [xc; xo];	y = [yc; yo];	z = [zc; zo];
	opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g', min(x),max(x),min(y),max(y));
	tol_0 = tol;
	for (k = 1:N_ITER)
		[Z, head] = gmtmbgrid_m(x, y, z, opt_I, opt_R, '-Mz', '-C1');
		zz = grdtrack_m(Z,head,[xc yc],'-Z')';
		tol = abs(zc * tol_0);				% From here on, 'tol' is a vector with size = numel(z)
		ind = abs(zc - zz(:)) > tol;		% Indices of those points that are no further than TOL from surface
		tol_0 = tol_0 * 0.7;				% Prepare already for next loop
		clear zz
		xc(ind) = NaN;		yc(ind) = NaN;	zc(ind) = NaN;
		if (k < N_ITER)
			x = [xc; xo];	y = [yc; yo];	z = [zc; zo];	% Update the base dataset (killed removed)
		end
		if (k == 1)
			the_ind = ind;
		else
			the_ind = the_ind | ind;		% Add this ieration result to the result of previous iterations
		end
	end

	% Now it's time to remove the NaNs that were inserted in place of the cleaned pts
	indNaN = isnan(zc);
	xc(indNaN) = [];		yc(indNaN) = [];		zc(indNaN) = [];

% ----------------------------------------------------------------------------------------
function out = find_BBs(list)
% Find the BoundingBoxs of all files in the cell array LIST
% The result is stored in the OUT struct, which has BB (numeric) and REGIONS (-R strings) members
	BB = zeros(numel(list),6);
	regions = cell(numel(list),1);
	for (n = 1:numel(list))
		info = gmtmex(['mbinfo -I' list{n}]);
		for (m = numel(info.text):-1:1)
			ind = strfind(info.text{m}, 'Minimum Latitude:');
			if (~isempty(ind)),		break,	end
		end
		if (m == 1)
			errordlg(['Info for file ' list{n} ' is broken. Couldn''t find Latitude bounds'], 'Error')
			error(['Can''t find Latitude bounds for file ' list{n}])
		end
		%Minimum Longitude:     -11.635798657   Maximum Longitude:     -11.391532952
		ind = strfind(info.text{m}, ':');
		t = strtok(info.text{m}(ind(1)+3:end));		BB(n,3) = str2double(t);
		t = strtok(info.text{m}(ind(2)+3:end));		BB(n,4) = str2double(t);
		m = m - 1;		% Do Long now
		ind = strfind(info.text{m}, ':');
		t = strtok(info.text{m}(ind(1)+3:end));		BB(n,1) = str2double(t);
		t = strtok(info.text{m}(ind(2)+3:end));		BB(n,2) = str2double(t);
		m = m + 4;		% Do Depth now
		ind = strfind(info.text{m}, ':');
		t = strtok(info.text{m}(ind(1)+3:end));		BB(n,5) = str2double(t);
		t = strtok(info.text{m}(ind(2)+3:end));		BB(n,6) = str2double(t);
		regions{n} = sprintf('-R%.12g/%.12g/%.12g/%.12g', BB(n,1:4));
	end
	out.BB = BB;		out.regions = regions;

% ----------------------------------------------------------------------------------------
function out = find_overlaps(list, K, info)
% Find the files whose BB partially intercept the file whose index is K
% The result is stored in the OUT vector which will contain the LIST indices that fullfill the condition

	count = false(numel(list), 1);
	for (n = 1:numel(list))
		if (n == K),	continue,	end		% This is the K file
		D = gmtmex(['gmtselect ' info.regions{K}], [info.BB(n,1) info.BB(n,3); info.BB(n,2) info.BB(n,4)]);
		if (~isempty(D) && ~isempty(D.data))
			count(n) = true;
		end
	end
	out = find(count);

% ----------------------------------------------------------------------------------------
function save_flags(handles, fname, ind, dims)
% DIMS = [n_pings n_beams]
	flags_file = [fname '.mask'];
	ind = reshape(ind, dims)';
	fid = fopen(flags_file, 'wb');
	fwrite(fid, dims, 'integer*4');
	fwrite(fid, ind(:), 'uchar');
	fclose(fid);

% ----------------------------------------------------------------------
function radio_imgSimple_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_imgShaded handles.radio_imgAmp handles.radio_imgSS],'Value',0)

% ----------------------------------------------------------------------
function radio_imgShaded_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_imgSimple handles.radio_imgAmp handles.radio_imgSS],'Value',0)

% ----------------------------------------------------------------------
function radio_imgAmp_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_imgSimple handles.radio_imgShaded handles.radio_imgSS],'Value',0)

% ----------------------------------------------------------------------
function radio_imgSS_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_imgSimple handles.radio_imgShaded handles.radio_imgAmp],'Value',0)

% ----------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% Show a quick view of the dataset as an image

	% If we haven't donne this yet (checking if 1st line in datalist has MB options)
	[pato,fname,ext] = fileparts(handles.fnameMB);
	if (strcmpi(ext, '.mb-1') && handles.first_datalist_scan)
		cmd = '';
		fid = fopen(handles.fnameMB,'rt');
		todos = fread(fid,'*char');		fclose(fid);
		if (todos(1) == '#' || todos(1) == '>')
			ind = find(todos == 10);		% Find new lines
			cmd = todos(1:ind(1)-1)';
		end
		if (~isempty(cmd))
			ind = strfind(cmd, ' -Z');
			if (~isempty(ind) && numel(cmd) >= ind(1)+3 && (double(cmd(ind(1)+3)) >= 49) && (double(cmd(ind(1)+3)) <= 53))
				handles.opt_Z = cmd(ind(1):ind(1)+3);
			end
			ind = strfind(cmd, ' -C');
			if (~isempty(ind))
				handles.opt_C = [' ' strtok(cmd(ind(1)+1:end))];
			end
			ind = strfind(cmd, ' -N');
			if (~isempty(ind))
				handles.opt_N = [' ' strtok(cmd(ind(1)+1:end))];
			end
			ind = strfind(cmd, ' -R');
			if (~isempty(ind))
				handles.opt_R = [' ' strtok(cmd(ind(1)+1:end))];
			end
			if (~isempty(handles.opt_C) || ~isempty(handles.opt_N) || ~isempty(handles.opt_R) || ~isempty(handles.opt_Z))
				set(handles.check_datalist, 'Enable','on')
			end
			handles.first_datalist_scan = false;
			handles.show_datalist_check = true;
			guidata(handles.figure1, handles)
		end
	end		

	opt_C = '';		opt_N = '';		opt_R = '';
	opt_Z = ' -Z5';
	if (get(handles.radio_imgSimple, 'Val')),		opt_Z = ' -Z1';
	elseif (get(handles.radio_imgShaded, 'Val')),	opt_Z = ' -Z2';
	elseif (get(handles.radio_imgAmp, 'Val')),		opt_Z = ' -Z4';
	end

	if (get(handles.check_datalist,'Val'))
		opt_C = handles.opt_C;			opt_N = handles.opt_N;	% For these we have buttons here
		if (~isempty(handles.opt_R)),	opt_R = handles.opt_R;		end
		if (~isempty(handles.opt_Z)),	opt_Z = handles.opt_Z;		end
	end

	if (~get(handles.check_trueGrid, 'Val'))
		I = gmtmex(['mbimport -I' handles.fnameMB opt_Z opt_N opt_C opt_R]);
		mirone(I)
	else
		opt_C = '';
		if (get(handles.radio_imgSS, 'Val')),		opt_C = '-C';
		elseif (get(handles.radio_imgAmp, 'Val')),	opt_C = '-Ca';
		end
		[out, BB, DS] = push_showPC_CB(handles.push_showPC, handles, opt_C);
		h = griding_mir([], 'nearneighbor', BB(1:4), out);
		if (DS ~= 0)
			if (get(handles.radio_imgSS, 'Val')),	DS = 2 * DS;	end		% Otherwise the default is too small
			external_drive(guidata(h), 'griding_mir', {sprintf('edit_x_inc,+%.12g',DS), sprintf('edit_y_inc,+%.12g',DS), ...
				sprintf('edit_S1_Neighbor,+%.12g',2.5*DS), 'check_Option_V,1'})	% No -X... because calling it directly
		end

	end
	if (handles.no_file && ishandle(handles.hMirFig)),	delete(handles.hMirFig),	end

%-------------------------------------------------------------------------------------
function jumping_done(N)
% ...
	if (nargin == 0),	N = 3;		end
	sz = get(0,'screensize');
	sz(3) = sz(3) / 2;		sz(4) = sz(4) / 2;
	h = msgbox('DONE');		set(h, 'units', 'pixels');		pos = get(h, 'pos');
	for (k = 1:N)
		pos(1:2) = round(rand(1,2) .* sz(3:4)) + sz(3:4) / 4;
		set(h, 'pos', pos)
		pause(0.8),
	end
	delete(h)
	

%-------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(hObject);
	end

% ----------------------------------------------------------------------
function h1 = show_mb_LayoutFcn()

	h1 = figure('Position',[520 441 250 330],...
		'Color',get(0,'factoryUicontrolBackgroundColor'),...
		'KeyPressFcn',@figure1_KeyPressFcn,...
		'MenuBar','none',...
		'Name','Show MB',...
		'NumberTitle','off',...
		'DoubleBuffer','on',...
		'Resize','off',...
		'HandleVisibility','Callback',...
		'Vis','off',...
		'Tag','figure1');

	uicontrol('Parent',h1, 'Position',[11 307 171 16],...
		'FontAngle','oblique',...
		'FontSize',10,...
		'FontWeight','demi',...
		'String','How to show this data?',...
		'Style','text');

	uicontrol('Parent',h1, 'Position',[201 303 40 23],...
		'Callback',@showMB_uiCB,...
		'FontSize',10,...
		'FontWeight','bold',...
		'String','?',...
		'TooltipString','Some help',...
		'Tag','push_help');

	uicontrol('Parent',h1,'Position',[10 154 231 136],'Style','frame');
	uicontrol('Parent',h1,'Position',[10 9 231 131],'Style','frame');

	uicontrol('Parent',h1, 'Position',[83 281 80 17],...
		'FontSize',10,...
		'String','Point Cloud',...
		'Style','text',...
		'Tag','text_PC');

	uicontrol('Parent',h1, 'Position',[83 130 80 17],...
		'FontSize',10,...
		'String','Plot Image',...
		'Style','text',...
		'Tag','text_PI');

	uicontrol('Parent',h1, 'Position',[50 253 70 19],...
		'BackgroundColor',[1 1 1],...
		'String',{'Point'; 'Cube'; 'Sphere'; 'Dyamond'; 'Cylinder'},...
		'Style','popupmenu',...
		'TooltipString','Symbol used in the point cloud',...
		'Value',1,...
		'Tag','popup_symb');

	uicontrol('Parent',h1, 'Position',[130 253 47 19],...
		'BackgroundColor',[1 1 1],...
		'Callback',@showMB_uiCB,...
		'String','0.005',...
		'Style','edit',...
		'TooltipString','Symbol size in unknown unites',...
		'Tag','edit_symbSize');

	uicontrol('Parent',h1, 'Position',[16 255 31 15],...
		'String','Symb', 'Style','text',...
		'Tag','text_Symb');

	uicontrol('Parent',h1, 'Position',[179 255 30 15],...
		'String','Size', 'Style','text',...
		'Tag','text_size');

	uicontrol('Parent',h1, 'Position',[31 223 100 20],...
		'Callback',@showMB_uiCB,...
		'String','Show point cloud',...
		'TooltipString','Show the point cloud.',...
		'Value',1,...
		'Tag','push_showPC');

	uicontrol('Parent',h1, 'Position',[134 225 90 17],...
		'String','Hide flaggeds',...
		'Style','checkbox',...
		'TooltipString','Do not show the beams flagged in the .esf file (previouly flagged points).',...
		'Tag','check_showFlagged');

	uicontrol('Parent',h1, 'Position',[31 188 118 26],...
		'Callback',@showMB_uiCB,...
		'String','Approximate cleaning',...
		'TooltipString','Do an automatic cleaning, show it and save results in file.',...
		'Tag','push_autoclean');

	uicontrol('Parent',h1, 'Position',[155 202 70 17],...
		'String','and show',...
		'Style','checkbox',...
		'TooltipString','Show the point cloud including the result of this cleanings (Memory consuming)',...
		'Tag','check_showCleaneds');

	uicontrol('Parent',h1, 'Position',[155 180 50 20],...
		'Callback',@showMB_uiCB,...
		'String','Params',...
		'TooltipString','Set parameters for the automatic cleaning',...
		'Tag','push_params');

	uicontrol('Parent',h1, 'Position',[31 160 100 20],...
		'Callback',@showMB_uiCB,...
		'Enable','off',...
		'String','Apply cleanings',...
		'TooltipString','Apply the automatic cleanings to the MB file',...
		'Tag','push_applyClean');

	uicontrol('Parent',h1, 'Position',[30 113 111 15],...
		'Callback',@showMB_uiCB,...
		'String','Simple bathymetry',...
		'Style','radiobutton',...
		'TooltipString','Plain color image',...
		'Tag','radio_imgSimple');

	uicontrol('Parent',h1, 'Position',[30 90 121 15],...
		'Callback',@showMB_uiCB,...
		'String','Shaded bathymetry',...
		'Style','radiobutton',...
		'Value',1,...
		'TooltipString','Shaded illuminated color image',...
		'Tag','radio_imgShaded');

	uicontrol('Parent',h1, 'Position',[30 67 71 15],...
		'Callback',@showMB_uiCB,...
		'String','Amplitude',...
		'Style','radiobutton',...
		'TooltipString','Gray scale amplitude plot',...
		'Tag','radio_imgAmp');

	uicontrol('Parent',h1, 'Position',[30 44 71 15],...
		'Callback',@showMB_uiCB,...
		'String','Side Scan',...
		'Style','radiobutton',...
		'TooltipString','Gray scale Side Scan plot',...
		'Tag','radio_imgSS');

	uicontrol('Parent',h1, 'Position',[190 80 51 15],...
		'String','grd?',...
		'Style','checkbox',...
		'Tooltip','Compute a true grid. But be warned that it may take some time.',...
		'Tag','check_trueGrid');

	uicontrol('Parent',h1, 'Position',[31 21 151 15],...
		'String','Use commands in datalist',...
		'Style','checkbox',...
		'TooltipString','The datalist may have it''s own commands. If checked, use those commands',...
		'Tag','check_datalist');

	uicontrol( 'Parent',h1, 'Position',[190 17 40 24],...
		'Callback',@showMB_uiCB,...
		'FontSize',10,...
		'FontWeight','bold',...
		'String','OK',...
		'Tag','push_OK');

% ----------------------------------------------------------------------------------
function showMB_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
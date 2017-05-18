function varargout = show_MB(varargin)
% Helper window to decide what to do when a multi-beam file, or datalist.mb-1, is loaded

%	Copyright (c) 2004-2017 by J. Luis
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
	show_MB_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right');

	handMir = varargin{1};
	handles.fnameMB = varargin{2};
	handles.whichFleder = handMir.whichFleder;
	handles.TDRver   = handMir.TDRver;
	handles.path_tmp = handMir.path_tmp;
	handles.no_file  = handMir.no_file;
	handles.hMirFig  = handMir.figure1;

	handles.opt_N = '';		handles.opt_C = '';		handles.opt_Z = '';	% Defaults that might be used later
	handles.first_datalist_scan = true;
	handles.show_datalist_check = false;
	handles.did_autoclean = false;
	handles.list_files = '';		% Will store the contents of datalists files when those are used.

	if (strcmp(varargin{3}, 'MBfbt') || strcmp(varargin{3}, 'MBall'))	% Can't make images of these
		set(handles.radio_image, 'Enable', 'off')
	end

	%------------ Give a Pro look (3D) to the frame boxes  ---------------------------------------------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) -------------------------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% ----------------------------------------------------------------------
function radio_pointCloud_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_image,'Value',0)
	set([handles.popup_symb handles.edit_symbSize handles.radio_justSee ...
		handles.radio_autoclean],'Enable','on')
	set([handles.radio_imgSimple handles.radio_imgShaded handles.radio_imgAmp ...
		handles.radio_imgSS handles.check_datalist],'Enable','off')

% ----------------------------------------------------------------------
function edit_symbSize_CB(hObject, handles)
% Just check that impossible figures are not acepted
	str = get(hObject,'String');        s = str2double(str);
	if (isnan(s) || s <= 0),	set(hObject,'String','0.01'),	return,		end

% ----------------------------------------------------------------------
function radio_justSee_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_autoclean,'Value',0)

% ----------------------------------------------------------------------
function radio_autoclean_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_justSee,'Value',0)

% ----------------------------------------------------------------------
function push_params_CB(hObject, handles)
% Will call a new window to help with the autocleaning choices
	msgbox('Sorry, not yet implemented. So far a 3 iterations procedure is applyied.')

% ----------------------------------------------------------------------
function push_applyClean_CB(hObject, handles)
% Apply the cleanings stored in the .mask file created previously.
	if (isempty(handles.list_files))		% Single file
		gmtmex(['mbflags -I' handles.fnameMB ' -E' handles.fnameMB '.mask'])
	else									% datalist file
		for (k = 1:numel(handles.list_files))
			gmtmex(['mbflags -I' handles.list_files{k} ' -E' handles.list_files{k} '.mask'])
		end
	end
	h = msgbox('DONE');
	pause(3)
	delete(h)

% ----------------------------------------------------------------------
function radio_image_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_pointCloud,'Value',0)
	set([handles.radio_imgSimple handles.radio_imgShaded handles.radio_imgAmp ...
		handles.radio_imgSS],'Enable','on', 'Val',0)
	if (handles.show_datalist_check)
		 set(handles.check_datalist,'Enable','on')
	end
	set(handles.radio_imgSimple, 'Val', 1)
	set([handles.popup_symb handles.edit_symbSize handles.radio_justSee ...
		handles.radio_autoclean],'Enable','off')

	if (handles.first_datalist_scan)		% If we haven't donne this yet (checking if 1st line in datalist has MB options)
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
			if (~isempty(handles.opt_C) || ~isempty(handles.opt_N) || ~isempty(handles.opt_Z))
				set(handles.check_datalist, 'Enable','on')
			end
			handles.first_datalist_scan = false;
			handles.show_datalist_check = true;
			guidata(handles.figure1, handles)
		end
	end

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
	msgbox(t, 'Help on show MB data')

% ----------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...
	tol = 0.03;			N_ITER = 3;		% For the aultocleaner cases
	if (get(handles.radio_image, 'Val'))
		show_as_img(handles),	return		% Done, bye bye
	end

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
	if (strcmpi(ext, '.mb-1') && get(handles.radio_autoclean, 'Val'))
		datalist = handles.fnameMB;
	end

	if (get(handles.radio_justSee, 'Val'))	% See the raw points
		D = gmtmex(sprintf('mbgetdata -I%s', handles.fnameMB));
		xyz = [D(1).data(:) D(2).data(:) D(3).data(:)];
		bb = [min(xyz)' max(xyz)']';		bb = bb(:)';
		write_flederFiles('points', fname, xyz, 'first', bb, par);
		comm = [' -data ' fname ' &'];		% A SD file
	else
		[xyz1,xyz2] = autocleaner(handles, tol, N_ITER, datalist);
		if(isempty(xyz1) && isempty(xyz2)),		return,		end			% Some error occured
		bb1 = [min(xyz1)' max(xyz1)']';		bb1 = bb1(:)';
		bb2 = [min(xyz2)' max(xyz2)']';		bb2 = bb2(:)';
		bbg = [min(bb1(1),bb2(1)) max(bb1(2),bb2(2)) min(bb1(3),bb2(3)) max(bb1(4),bb2(4)) ...
			   min(bb1(5),bb2(5)) max(bb1(6),bb2(6))];

		fid = write_flederFiles('scene_pts', fname, [], 'begin', bbg, par);
		write_flederFiles('scene_pts', fid, xyz1, 'sec', [bb1 bbg], par);
		PTparams.Symbol = 3;	PTparams.ColorBy = 0;	par.PTparams = PTparams;
		write_flederFiles('scene_pts', fid, xyz2, 'end', [bb2 bbg], par);
		fname = [fname '.scene'];			% Because name was changed in write_flederFiles()
		comm = [' -scene ' fname ' &'];		% A SCENE file

		handles.did_autoclean = true;		% So we know that we can update the corresponding .esf
		set(handles.push_applyClean, 'Enable', 'on')
	end

	if (handles.whichFleder),	fcomm = ['iview4d' comm];		% Free viewer
	else						fcomm = ['fledermaus' comm];	% The real thing
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
			dos(fcomm);		% s is always 0 (success) as long as iview4d is accessible, even when the comm fails
		else
			errordlg('Unknown platform.','Error'),	return
		end
	catch
		errordlg('I could not find Fledermaus. Hmmm, do you have it?','Error')
	end

%-------------------------------------------------------------------------------------
function [xyz, xyzK] = autocleaner(handles, tol, N_ITER, datalist)
% ...
	if (isempty(datalist))			% Single file, simplest case
		D = gmtmex(sprintf('mbgetdata -I%s', handles.fnameMB));
		x = D(1).data(:);		y = D(2).data(:);		z = D(3).data(:);
		% Compute increment as 3 times the typical point spacing fetch from data's first row.
		dy = abs(median(diff(D(2).data(1,:))));
		dx = abs(median(diff(D(1).data(1,:)))) * cos(y(1));
		opt_I = sprintf('-I%f', sqrt(dx*dx + dy*dy) * 3);
		opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g', min(x),max(x),min(y),max(y));
		ind = [];
		for (k = 1:N_ITER)
			[x,y,z, ind] = iteration_cleaner(x,y,z, opt_R, opt_I, tol, ind);	% return 'cleaned'
			tol = tol * 0.7;
		end

		% Now it's time to remove the NaNs that were inserted in place of the cleaned pts
		indNaN = isnan(x);
		x(indNaN) = [];		y(indNaN) = [];		z(indNaN) = [];
		xyz = [x y z];
		xyzK = [D(1).data(ind) D(2).data(ind) D(3).data(ind)];		% The Killed points

		save_flags(handles, handles.fnameMB, ind, size(D(1).data))	% Save the flags in a .mask file
	else
		[xc,yc,zc, xk,yk,zk] = autocleaner_list(handles, tol, N_ITER);
		xyz  = [xc yc zc];	clear xc yc zc		% AWFULL MEMORY WASTE
		xyzK = [xk yk zk];
	end

% 	[Z, head] = gmtmbgrid_m(x, y, z, '-I0.001', opt_R, '-Mz', '-C3');
% 	aux_funs('showgrd', struct('geog',1), Z, head, 'Autocleaned')
	
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
	n_lines = numel(todos);

	patoMother = fileparts(handles.fnameMB);			% Path of the datlist file
	c = false(n_lines,1);
	msg = '';
	for (k = 1:n_lines)			% Passar isto a uma funcao chamda sanitize_datalist()
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
function [xc,yc,zc, xk,yk,zk] = autocleaner_list(handles, tol, N_ITER)
% Scan all files in a datalist file and clean each file individually against the data of
% all files. This is awfully inefficient since we should work only over sub-regions. Future work.

	list_files = sanitize_datalist(handles);	% Read and sanitize a datalist file. Returns cell array with full fnames
	handles.list_files = list_files;
	guidata(handles.figure1, handles)			% Save this for use in "Apply cleanings"

	n_lines = numel(list_files);
	if (n_lines <= 1)
		xc=[];yc=[];zc=[];xk=[];yk=[];zk=[];
		return
	end

	aguentabar('title',sprintf('cleaning %d files',n_lines))
	fname_tmp = [handles.path_tmp 'datalist.mb-1'];
	for (k = 1:n_lines)
		fid = fopen(fname_tmp, 'wt');			% Open a TMP file to store the names of the 'other' files
		v = 1:n_lines;		v(k) = [];
		lista_t = list_files(v);
		for (n = 1:numel(lista_t))
			fprintf(fid, '%s\n', lista_t{n});	% Write the TMP datalist file with all but current (loop's) file
		end
		fclose(fid);

		Dcurr   = gmtmex(['mbgetdata -I', list_files{k}]);
		Dothers = gmtmex(['mbgetdata -I' fname_tmp]);
		if (k == 1)								% Compute working inc as 3 times typical pt spacings
			dy = abs(median(diff(Dcurr(2).data(1,:))));
			dx = abs(median(diff(Dcurr(1).data(1,:)))) * cos(Dcurr(2).data(1));
			opt_I = sprintf('-I%f', sqrt(dx*dx + dy*dy) * 3);
		end
		[x,y,z, ind] = iteration_cleaner_list(Dcurr(1).data(:),Dcurr(2).data(:),Dcurr(3).data(:), ...
			Dothers(1).data(:),Dothers(2).data(:),Dothers(3).data(:), opt_I, tol, N_ITER);

		save_flags(handles, list_files{k}, ind, size(Dcurr(1).data))		% Time to save flags
		if (k == 1)
			xc = x;		yc = y;		zc = z;
			xk = Dcurr(1).data(ind);	yk = Dcurr(2).data(ind);	zk = Dcurr(3).data(ind);
		else
			xc = [xc; x];	yc = [yc; y];		zc = [zc; z]; %#ok<AGROW>
			xk = [xk; Dcurr(1).data(ind)];		yk = [yk; Dcurr(2).data(ind)];	zk = [zk; Dcurr(3).data(ind)]; %#ok<AGROW>
		end
		aguentabar(k/n_lines)
	end
	
%-------------------------------------------------------------------------------------
function [xc,yc,zc, the_ind] = iteration_cleaner_list(xc,yc,zc, xo,yo,zo, opt_I, tol, N_ITER)
% XC,... coordinates of current file beeing cleaned
% XO,... coordinates of all other points from the remaining files in the datalist file (all but current)
% TOL is a scalar with the satrting percentage of the water depth. Points further than this from the
% surface are flagged. Furthermore, this process is run N_ITER times where, in eache ietr, TOL is shrink by 0.7
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
	indNaN = isnan(xc);
	xc(indNaN) = [];		yc(indNaN) = [];		zc(indNaN) = [];

% ----------------------------------------------------------------------------------------
function save_flags(handles, fname, ind, dims)
% DIMS = [n_pings n_beams]
	flags_file = [fname '.mask'];
	ind = reshape(ind, dims)';
	fid = fopen(flags_file, 'wb');
	fwrite(fid, dims, 'integer*4');
	fwrite(fid, ind(:), 'uchar');
	fclose(fid);

% -------------------------------------------------------------------------------------
function show_as_img(handles)
% Show a quick view of the dataset as an image
	opt_C = '';		opt_N = '';
	opt_Z = ' -Z5';
	if (get(handles.radio_imgSimple, 'Val')),		opt_Z = ' -Z1';
	elseif (get(handles.radio_imgShaded, 'Val')),	opt_Z = ' -Z2';
	elseif (get(handles.radio_imgAmp, 'Val')),		opt_Z = ' -Z4';
	end

	if (get(handles.check_datalist,'Val'))
		opt_C = handles.opt_C;		opt_N = handles.opt_N;	% For these we have buttons here
		if (~isempty(handles.opt_Z))
			opt_Z = handles.opt_Z;
		end
	end
	I = gmtmex(['mbimport -I' handles.fnameMB opt_Z opt_N opt_C]);
	I.image = flipdim(I.image,1);
	mirone(I)
	if (handles.no_file && ishandle(handles.hMirFig)),	delete(handles.hMirFig),	end

%-------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(hObject);
	end

% ----------------------------------------------------------------------
function show_MB_LayoutFcn(h1)

	set(h1, 'Position',[520 441 251 359],...
		'Color',get(0,'factoryUicontrolBackgroundColor'),...
		'KeyPressFcn',@figure1_KeyPressFcn,...
		'MenuBar','none',...
		'Name','Show MB',...
		'NumberTitle','off',...
		'DoubleBuffer','on',...
		'Resize','off',...
		'HandleVisibility','Call',...
		'Tag','figure1');

	uicontrol('Parent',h1,'Position',[10 198 231 131],'Style','frame');
	uicontrol('Parent',h1,'Position',[10 38 231 151],'Style','frame');

	uicontrol('Parent',h1, 'Position',[21 306 79 15],...
		'Call',@showMB_uiCB,...
		'String','Point cloud',...
		'Style','radiobutton',...
		'Tooltip','Show the beams as a point cloud in FlederMaus',...
		'Value',1,...
		'Tag','radio_pointCloud');

	uicontrol('Parent',h1, 'Position',[51 275 70 19],...
		'BackgroundColor',[1 1 1],...
		'String',{'Point'; 'Cube'; 'Sphere'; 'Dyamond'; 'Cylinder'},...
		'Style','popupmenu',...
		'Tooltip','Symbol used in the point cloud',...
		'Value',1,...
		'Tag','popup_symb');

	uicontrol('Parent',h1, 'Position',[131 275 47 19],...
		'BackgroundColor',[1 1 1],...
		'Call',@showMB_uiCB,...
		'String','0.005',...
		'Style','edit',...
		'Tooltip','Symbol size in unknown unites',...
		'Tag','edit_symbSize');

	uicontrol('Parent',h1, 'Position',[17 277 31 15],...
		'String','Symb',...
		'Style','text',...
		'Tag','text_Symb');
	
	uicontrol('Parent',h1, 'Position',[179 277 31 15],...
		'String','Size',...
		'Style','text',...
		'Tag','text_size');

	uicontrol('Parent',h1, 'Position',[51 252 91 15],...
		'Call',@showMB_uiCB,...
		'String','Just see them',...
		'Style','radiobutton',...
		'Tooltip','Just visualize the intire point cloud',...
		'Value',1,...
		'Tag','radio_justSee');

	uicontrol('Parent',h1, 'Position',[51 231 118 15],...
		'Call',@showMB_uiCB,...
		'String','Automatic cleaning',...
		'Style','radiobutton',...
		'Tooltip','Do an automatic cleaning, show it and save results in file.',...
		'Tag','radio_autoclean');

	uicontrol('Parent',h1, 'Position',[71 206 50 20],...
		'Call',@showMB_uiCB,...
		'String','Params',...
		'Tooltip','Set parameters for the automatic cleaning',...
		'Enable','off',...
		'Tag','push_params');

	uicontrol('Parent',h1, 'Position',[138 206 95 20],...
		'Call',@showMB_uiCB,...
		'String','Apply cleanings',...
		'Tooltip','Apply the automatic cleanings to the MB file',...
		'Enable','off',...
		'Tag','push_applyClean');

	uicontrol('Parent',h1, 'Position',[21 165 79 15],...
		'Call',@showMB_uiCB,...
		'String','Plot image',...
		'Style','radiobutton',...
		'Tooltip','Plot a quick view image',...
		'Tag','radio_image');

	uicontrol('Parent',h1, 'Position',[51 142 111 15],...
		'Call',@showMB_uiCB,...
		'String','Simple bathymetry',...
		'Style','radiobutton',...
		'Tooltip','Plain color image',...
		'Enable','off',...
		'Tag','radio_imgSimple');

	uicontrol('Parent',h1, 'Position',[51 119 121 15],...
		'Call',@showMB_uiCB,...
		'String','Shaded bathymetry',...
		'Style','radiobutton',...
		'Tooltip','Shaded illuminated color image',...
		'Enable','off',...
		'Tag','radio_imgShaded');

	uicontrol('Parent',h1, 'Position',[51 96 71 15],...
		'Call',@showMB_uiCB,...
		'String','Amplitude',...
		'Style','radiobutton',...
		'Tooltip','Gray scale amplitude plot',...
		'Enable','off',...
		'Tag','radio_imgAmp');

	uicontrol('Parent',h1, 'Position',[51 73 71 15],...
		'Call',@showMB_uiCB,...
		'String','Side Scan',...
		'Style','radiobutton',...
		'Tooltip','Gray scale Side Scan plot',...
		'Enable','off',...
		'Tag','radio_imgSS');

	uicontrol('Parent',h1, 'Position',[21 50 151 15],...
		'String','Use commands in datalist',...
		'Style','checkbox',...
		'Tooltip','The datalist has it''s own commands. If checked, use those commands',...
		'Enable','off',...
		'Tag','check_datalist');

	uicontrol('Parent',h1, 'Position',[30 337 181 16],...
		'FontAngle','oblique',...
		'FontSize',10,...
		'FontWeight','demi',...
		'String','How to show this data?',...
		'Style','text');

	uicontrol('Parent',h1, 'Position',[11 7 50 23],...
		'String','?',...
		'Call',@showMB_uiCB,...
		'Tooltip','Help/Info',...
		'FontSize',9,...
		'FontWeight','bold',...
		'Tag','push_help');

	uicontrol('Parent',h1, 'Position',[150 7 90 24],...
		'String','RUN',...
		'Call',@showMB_uiCB,...
		'FontSize',10,...
		'FontWeight','bold',...
		'Tag','push_OK');

% ----------------------------------------------------------------------------------
function showMB_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
function varargout = grid_calculator(varargin)
% An array calculator 
%
% The big trick here is to store the arrays in structure members of a 'grid' struct
% container and at same time maintain an updated map between those struct members and
% variables char strings with the same name. This, and lots of parsings, allow not 
% only to run the whole computation as an eval('command') from within Matlab but also
% to break that 'command' in its tokens and, now with lots of parsings, execute them
% in a way that also works in the stand-alone (compiled) version.

%	Copyright (c) 2004-2020 by J. Luis
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

% $Id: grid_calculator.m 10410 2018-05-18 15:32:41Z j $

	old_vis = get(0,'ShowHiddenHandles');
	set(0,'ShowHiddenHandles','on')
	h_figs = findobj('Type','figure');		% Do this before creating the calc figure
	set(0,'ShowHiddenHandles',old_vis)

	hObject = figure('Vis','off');
	grid_calculator_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	handles.grid_patos = [];    % To hold paths of eventual future loaded grids
	handles.loaded_grid = [];   % To hold names of eventual future loaded grids
	handles.name_str = [];
	handles.BL = [];
	handles.reader = '';

	if (~isempty(varargin) && ishandle(varargin{1}))	% Not tested if varargin{1} is a Fig handle
		handMir = guidata(varargin{1});
		handles.home_dir = handMir.home_dir;		handles.last_dir = handMir.last_dir;
		handles.work_dir = handMir.work_dir;		handles.path_tmp = handMir.path_tmp;
		handles.hMirFig  = handMir.figure1;
		handles.IamCompiled = handMir.IamCompiled;
		handles.version7 = handMir.version7;
		BL = getappdata(varargin{1},'BandList');
		if (~isempty(BL))
			handles.name_str = BL{3}(2:end,2);
			h_figs = [];			% Tear off the net for the 'h_figs' fishing
			handles.BL = BL{2};
			handles.reader = BL{end};
			set(hObject, 'Name', 'Bands calculator')
			if (isappdata(handMir.axes1, 'LandSAT8_MTL'))
				set(handles.push_Trad, 'Vis', 'on')		% They would be set visible because h_figs == []
			end
		end
	elseif (~isempty(varargin) && isa(varargin{1}, 'char'))
		do_tests;
		delete(hObject)
		return
	else
		handles.home_dir = cd;		handles.last_dir = cd;		handles.work_dir = cd;
		handles.path_tmp = [handles.home_dir '/tmp'];		handles.IamCompiled = false;
		handles.hMirFig  = [];
	end

	% Fish whatever arrays are in memory now (hopefully)
	if (~isempty(h_figs))
		n = 1;
		for (i = 1:numel(h_figs))
			hand_fig = guidata(h_figs(i));
			% Use a try->catch to easily fish only the "filled" Mir figs
			try
				Z = getappdata(hand_fig.figure1,'dem_z');
				handles.home_dir = hand_fig.home_dir;		% Since I don't know which is the last good one
				handles.last_dir = hand_fig.last_dir;
				handles.work_dir = hand_fig.work_dir;
			catch
				continue;
			end
			if (isempty(Z)),	continue,	end
			name = get(h_figs(i),'Name');
			ind = strfind(name,' @ ');
			if (~isempty(ind))
				name = ddewhite(name(1:ind-1));
			end
			[pathstr,name,ext] = fileparts(name);
			name = strrep(name, ' ', '_');		% No bloody blanks in names
			handles.name_str{n} = [name ext];
			handles.h_figs(n) = h_figs(i);      % Save the figure handles
			if (n == 1 && ~isempty(getappdata(hand_fig.axes1, 'LandSAT8')))
				set(handles.push_Trad, 'Vis', 'on')
			end
			n = n + 1;
		end
	else
		handles.h_figs = handMir.figure1;
	end

	% ------------ Create a Semaforo -----------------------------------------------------
	semaforo(handles, 'green')

	% Fill the listbox with the names of the in-memory arrays
	if (~isempty(handles.name_str))
        set(handles.listbox_inArrays,'String',handles.name_str)
	end

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% ------------------------------------------------------------------------
function do_tests
% Run some tests that compare the ML eval(comm) solution against the local
% parsing one used in the stand-alone version that is not able to use the eval
	grid.a = ones(5,5)*2;		grid.b = ones(5,5)*3;		grid.c = ones(5,5)*4;
	comms = {'grid.a * (grid.b / grid.c) + log10(999)', ...
			'grid.a + grid.b - grid.c + 100', 'grid.a * sin(grid.b)', ...
			'grid.a * grid.b + grid.c / cos(grid.b)', 'grid.a * grid.b + grid.c / grid.b', ...
			'grid.b - grid.a * (grid.b / grid.a)', 'atan(grid.a) * 5 + 4 * 4', ...
			'(grid.a + grid.b) * (grid.c + grid.a)', 'grid.a * -1', 'grid.c + -1'};
	for (k = numel(comms):-1:1)
		comm = move_operator(comms{k});
		commM = strrep(comm,'*','.*');	commM = strrep(commM,'/','./');	commM = strrep(commM,'^','.^');
		resp_ML = eval(commM);
		[resp_JL, msg] = stalone(comm, grid);
		if (~isempty(msg))
			disp([sprintf('grid_calculator:Test %d',k), ' ERROR: ' msg]),	continue
		elseif (isempty(resp_JL))
			disp([sprintf('grid_calculator:Test %d',k), 'Programming error, result is empty']), continue
		end
		difa = resp_ML - resp_JL;
		if (abs(max(difa(:))) < eps*10)
			disp(['Test ' num2str(k) '  PASSED'])
		else
			disp(['Test ' comms{k} '  -->FAILED'])
		end
		%disp(['Test ' num2str(k) '  Max difa = ' num2str(max(difa(:)))])
	end

% ------------------------------------------------------------------------
function listbox_inArrays_CB(hObject, handles, manual)
% if this is a doubleclick, copy the array's name to the edit box
% Using 3 argins has lso the same effect
	if (nargin == 2 && ~strcmp(get(gcbf,'SelectionType'),'open')),	return,		end
	str = get(hObject,'String');
	sel = get(hObject,'Value');
	com = get(handles.edit_command,'String');
	[t,r] = strtok(com);
	if (isempty(t) || (t(1) == '&' && isempty(r)))		% Only one grid name, replace it by the new selection
		set(handles.edit_command,'String', ['&' str{sel}])
	else
		set(handles.edit_command,'String', [com ' &' str{sel}])
	end
	handMir = guidata(handles.h_figs(1));
	if (isappdata(handMir.axes1, 'LandSAT8_MTL'))
		ind = strfind(str{sel}, '_B');						% Find the band number
		if (isempty(ind)),		return,		end				% For sure not a Lansat8 band
		band = str2double(str{sel}(ind(end)+2:end));
		aux_funs('set_LandSat8_band_pars', handMir, band)	% Set this as the active band
	end

% ------------------------------------------------------------------------
function push_digit_CB(hObject, handles)
% Takes care of the 0-9 digits
	set(handles.edit_command,'String', [get(handles.edit_command,'String') get(hObject, 'String')])

% ------------------------------------------------------------------------
function push_equal_CB(hObject, handles)
	%set(handles.edit_command,'String', [get(handles.edit_command,'String') ' = '])
	push_compute_CB(hObject, handles)

% ------------------------------------------------------------------------
function push_op_CB(hObject, handles)
% Takes care of +-*/
	t = [' ' get(hObject,'String'), ' '];
	set(handles.edit_command,'String', [get(handles.edit_command,'String') t])

% ------------------------------------------------------------------------
function push_funs_CB(hObject, handles)
% Takes care of sin,cos,tan.log,.../
	t = [' ' get(hObject,'String'), '( '];
	set(handles.edit_command,'String', [get(handles.edit_command,'String') t])

% ------------------------------------------------------------------------
function push_loadGrid_CB(hObject, handles)
% This function doesn't realy loads the grid. It only stores the grid name.
% True loading is donne in "Compute"
	str1 = {'*.grd;*.GRD;*.nc;*.NC', 'Grid files (*.grd,*.GRD,*.nc,*.NC)';'*.*', 'All Files (*.*)'};
	[FileName,PathName,handles] = put_or_get_file(handles, str1,'Select GMT grid','get');
	if isequal(FileName,0),		return,		end
	str = get(handles.listbox_inArrays,'String');
	str{end+1} = FileName;
	set(handles.listbox_inArrays,'String',str);
	handles.grid_patos{end+1} = PathName;
	handles.loaded_grid{end+1} = FileName;
	guidata(hObject,handles)

% ------------------------------------------------------------------------
function push_leftPar_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '( '])

% ------------------------------------------------------------------------
function push_rightPar_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' )'])

% ------------------------------------------------------------------------
function push_Trad_CB(hObject, handles)
% Compute Brightness temperature of a LandSat8 termal band, OR Radiance/Reflectance
	handMir = guidata(handles.h_figs(1));
	pars = getappdata(handMir.axes1, 'LandSAT8');
	if (isempty(pars)),	warndlg('Sorry, this is not a Landsat8 band.','Warning'),	return,	end

	semaforo(handles, 'red'),	pause(0.001)
	out = push_compute_CB([], handles);			% See if we have an array loaded in the edit box
	semaforo(handles, 'red'),	pause(0.001)	% Red again because it was set green above
	if (~isempty(out))
		T = out.a;
		if (~isa(T, 'uint16'))
			semaforo(handles, 'green')
			errordlg('This is clearly NOT an Landsat 8 band data.', 'Error'),	return
		end
		[p,s] = fileparts(strtok(get(handles.edit_command, 'String')));
		suff = s(2:end);		% Suffix
		%ind = strfind(suff, '_B');
		%band = str2double(suff(ind(end)+2:end));
		%aux_funs('set_LandSat8_band_pars', handMir, band);	% If this is a VRT, update the MTL parameters for this band
	else
		if (numel(get(handles.listbox_inArrays, 'Str')) > 1)
			%errordlg('Have more then one band, so you must select the one you wish to process.', 'Error')
			listbox_inArrays_CB(handles.listbox_inArrays, handles, true)	% Simulate the double click
			push_Trad_CB(hObject, handles)		% Recursive call
			return
		end
		T = getappdata(handMir.figure1,'dem_z');
		[p, suff] = fileparts(get(handMir.figure1,'Name'));
	end
	suff = [' [' suff ']'];

	indNaN = (T == 0);		% We need these guys in all the cases

	s = get(hObject, 'String');
	if (strcmp(s, 'Trad'))
		if (pars.band < 10)
			warndlg('Nope. Brightness temperature is only for bands 10 or 11. Not yours', 'Error')
			semaforo(handles, 'green')
			return
		end
		opt_T = sprintf('-T%f/%f/%f/%f', pars.rad_mul, pars.rad_add, pars.K1, pars.K2);
		T = grdutils(T, opt_T);
		tmp.name = ['Brightness temperture' suff];
	elseif (strcmp(s, 'R(toa)'))
		opt_M = sprintf('-M%f', pars.rad_mul);	opt_A = sprintf('-A%f', pars.rad_add);
		T = grdutils(T, opt_M, opt_A);
		tmp.name = ['Radiance TOA' suff];
	elseif (strcmp(s, 'Rho(toa)'))
		if (pars.reflect_mul == 1)
			warndlg('Computing Reflectance for Thermal bands is not defined.','Error')
			semaforo(handles, 'green')
			return
		end
		s_elev = sin(pars.sun_elev * pi/180);
		opt_M = sprintf('-M%f', pars.reflect_mul/s_elev);	opt_A = sprintf('-A%f', pars.reflect_add/s_elev);
		T = grdutils(T, opt_M, opt_A);
		tmp.name = ['Reflectance TOA' suff];
	elseif (strcmp(s, 'Rho(surf)'))
		if (pars.reflect_max == 0)
			warndlg('Computing "At Surface Reflectance" for Thermal bands is not defined.','Error')
			semaforo(handles, 'green')
			return
		end
		s_elev = sin(pars.sun_elev * pi/180);
		Esun = (pi * pars.sun_dist ^2) * pars.rad_max / pars.reflect_max;
		TAUv = 1.0;		TAUz = s_elev;		Esky = 0.0;		sun_prct=1;
		if (pars.band == 6 || pars.band == 7 || pars.band == 9)
			TAUz = 1.0;
		end
		Sun_Radiance = TAUv * (Esun * s_elev * TAUz + Esky) / (pi * pars.sun_dist ^2);
		darkDN = getDarkDN(T, 0.0001, indNaN);
		opt_M = sprintf('-M%f', pars.rad_mul);		opt_A = sprintf('-A%f', pars.rad_add);
		T = grdutils(T, opt_M, opt_A);
		radiance_dark = pars.rad_mul * darkDN + pars.rad_add;	% 0.01%
		LHaze = radiance_dark - sun_prct * Sun_Radiance / 100;
		grdutils(T, sprintf('-A%f', -LHaze));
		grdutils(T, sprintf('-M%f', 1/Sun_Radiance));
		grdutils(T, '-F</0/0')
		grdutils(T, '-F>/1/1')
		tmp.name = ['Surface Reflectance [COST]' suff];
	end
	tmp.head = handMir.head;
	zMinMax = grdutils(T,'-L');
	tmp.head(5) = zMinMax(1);	tmp.head(6) = zMinMax(2);
	T(indNaN) = NaN;	clear indNaN
	semaforo(handles, 'green')
	mirone(T, tmp, handMir.figure1)
	% ---------------------------------------
	function darkDN = getDarkDN(DN, pct, indNaN)
		DN = DN(DN > 0);
		DN = nth_element(DN(:), round(numel(DN) * pct));
		darkDN = double(DN(round(numel(DN) * pct)));

% ------------------------------------------------------------------------
function out = push_compute_CB(hObject, handles)
% ...
	out = [];
	com = get(handles.edit_command,'String');
	if (isempty(com)),	return,		end

	semaforo(handles, 'red')

	if (~isempty(handles.BL))		% We are dealing with Bands arithmetics
		if (nargout),	out.a = bandArithm(handles, com);		% To be consistent. Other paths also return a struct
		else,			bandArithm(handles, com);
		end
		semaforo(handles, 'green')
		return
	end

	% Those start at 2 because they are meant to be used only when grid names apear repeatedly
	in_g_count = 2;     out_g_count = 2;	grid = [];

	com = move_operator(com);		% Make sure operators are not "glued" to operands (it currently screws names with '-' characters)
	k = strfind(com,'&');

	try								% Wrap it here to report any possible error
		if (~isempty(k))			% We have grids
			for (i = 1:numel(k))	% Loop over grids
				tok = strtok(com(k(i)+1:end));
				if (isempty(handles.name_str))
					n = [];						% Here we know that we don't have any pre-loaded grid
				else
					n = find(strcmp(tok,handles.name_str));     % n != []  when grid is already in memory
				end

				load_it = false;
				if (isempty(n))					% Grid needs to be loaded
					load_it = true;				% Flag to signal that grid must be loaded
					n_load = find(strcmp(tok,handles.loaded_grid));
					if (numel(n_load) > 1)		% Grid name comes out more than once
						n_load = n_load(out_g_count);
						out_g_count = out_g_count + 1;
					end
				elseif (numel(n) > 1)			% Grid name comes out more than once
					n = n(1);
					in_g_count = in_g_count + 1;
				end

				if (~load_it),	hand_fig = guidata(handles.h_figs(n));	end

				if (i == 1)
					if (load_it)
						handtmp = struct('grdMaxSize',1e30, 'ForceInsitu',false, ...
							'IamCompiled', handles.IamCompiled, 'path_tmp',handles.path_tmp);
						[grid_t, tmp.X, tmp.Y, srsWKT, handtmp] = read_grid(handtmp, [handles.grid_patos{n_load} tok], 'GMT');
						if (isempty(grid_t))			% Shit happened
							if (isempty(findobj(0,'type','figure','Name','ERROR')))		% If no error message, print one
								errordlg(['Unknown error while reading grid ' tok],'ERROR')
							end
							return
						end
						tmp.head = handtmp.head;
						%grid_t = double(grid_t);      % grid reading don't outputs doubles
					else
						grid_t = (getappdata(hand_fig.figure1,'dem_z'));
						tmp.X = getappdata(hand_fig.figure1,'dem_x');
						tmp.Y = getappdata(hand_fig.figure1,'dem_y');
						tmp.head = hand_fig.head;
					end
					grid.(char(i+96)) = grid_t;
				else
					if (load_it)
						handtmp = struct('grdMaxSize',1e30, 'ForceInsitu',false, ...
							'IamCompiled', handles.IamCompiled, 'path_tmp',handles.path_tmp);
						[grid_t, tmp.X, tmp.Y, srsWKT, handtmp] = read_grid(handtmp, [handles.grid_patos{n_load} tok], 'GMT');
						tmp.head = handtmp.head;
						%grid_t = double(grid_t);      % grid reading don't outpus doubles
					else
						grid_t = (getappdata(hand_fig.figure1,'dem_z'));
					end
					grid.(char(i+96)) = grid_t;
				end
			end         % Loop over grids

			for (i = 1:numel(k))
				k = strfind(com,'&');				% We need to recompute '&' positions because they change bellow
				[tok, r] = strtok(com(k(i)+1:end));	% Here we get the grid's name
				kf = k(i) + numel(tok);				% Find the position of last char of grid's name
				if (strncmp(r,' grid ', 6))			% Patch against those names like "Cropped grid"
					com(kf+1:kf+6) = [];
				end
				com = [com(1:k(i)) 'grid.' char(i+96) ' ' com(kf+1:end)];
			end

			if (nargout),	out = grid;		return,		end		% Landsat functions process this separatelly.

		end			% IF We have grids

		if (isempty(grid))
			errordlg('This is a grid calculator, but none of your operands is a matrix. Bye.', 'Error')
			semaforo(handles, 'green')
			return
		end

		com = strrep(com,'&','');					% Remove the '&' characters
		com = strrep(com,'e - ','e-');				% base 10 numbers cannot have those spaces
		com = strrep(com,'e + ','e+');

		if (~handles.IamCompiled)
			com_s = strrep(com,   '*', '.*');	com_s = strrep(com_s, '..*', '.*');		% Play safe
			com_s = strrep(com_s, '/', './');	com_s = strrep(com_s, '../', './');
			com_s = strrep(com_s, '^', '.^');	com_s = strrep(com_s, '..^', '.^');
			fnames = fields(grid);
			if (~handles.version7)				% For debugging porposes only
				for (k = 1:numel(fnames)),	grid.(fnames{k}) = double(grid.(fnames{k}));	end
			else
				for (k = 1:numel(fnames)),	grid.(fnames{k}) = single(grid.(fnames{k}));	end
			end
			resp = eval(com_s);
		else
			com_s = strrep(com,   '.*', '*');	com_s = strrep(com_s, './', '/');
			com_s = strrep(com_s, '.^', '^');
			[resp, msg] = stalone(com_s, grid);
			if (isempty(resp) || ~isempty(msg)),	semaforo(handles, 'green'),		end
			if     (~isempty(msg)),		errordlg(['ERROR: ' msg]),	return
			elseif (isempty(resp)),		errordlg('Programming error. Result is empty', 'Error'),	return
			end
		end
		resp = single(resp);
		semaforo(handles, 'green')

		[zzz] = grdutils(resp,'-L');  z_min = zzz(1);     z_max = zzz(2);
		if (~isempty(k))					% We had grids in input
			tmp.head(5:6) = [z_min z_max];
			tmp.name = 'Calculated_grid';
			if (~load_it),	mirone(resp, tmp, hand_fig.figure1);	% For the case we have proj info
			else,			mirone(resp, tmp);
			end
		elseif (numel(resp) > 1)			% 'resp' is a array. Construct a fake grid
			[m,n] = size(resp);
			%resp = single(resp);
			tmp.X = 1:n;        tmp.Y = 1:m;
			%[zzz] = grdutils(resp,'-L');  z_min = zzz(1);     z_max = zzz(2);
			tmp.head = [1 n 1 m z_min z_max 0 1 1];
			tmp.name = 'Calculated_array';
			mirone(resp,tmp);
		else								% Computations that do not involve grids
			txt = sprintf('%.10f',resp);
			while (txt(end) == '0')			% Remove trailing zeros
				txt(end) = [];
			end
			set(handles.edit_command,'String',txt)
		end
	catch
		semaforo(handles, 'green')
		errordlg(lasterr,'Error')
	end

% ------------------------------------------------------------------------
function out = bandArithm(handles, com)
% ...
	com = move_operator(com);			% Make sure operators are not "glued" to operands
	k = strfind(com,'&');
	semaforo(handles, 'red'),	pause(0.001)

	try									% Wrap it here to report any possible error
		if (~isempty(k))				% We have grids
			N = zeros(1, numel(k));
			for (i = 1:numel(k))		% Loop over bands
				tok = strtok(com(k(i)+1:end),' )+-*/^');
				n = find(strcmp(tok,handles.name_str));		% n ~= [] when it is in memory
				%grid.(char(n+96)) = double(handles.BL(:,:,n)) / double(intmax_(class(handles.BL)));
				if (strcmp(handles.reader{1}, 'GDAL'))
					Z = aux_funs('get_layer_n',handles.hMirFig, n);
					if (isempty(Z)),	Z = handles.BL(:,:,1);	end		% Happens when data is uint8
				else
					Z = handles.BL(:,:,n);	% For now datasets read with multibandread must be all in mem
				end

				if (nargout),	out = Z;	return,		end

				%grid.(char(n+96)) = double(Z) / double(intmax_(class(Z)));
				if (numel(k) > 1)			% Promote and Normalize
					if (isa(Z, 'uint16'))
						grid.(char(n+96)) = grdutils(Z, sprintf('-M%f',1 / double(intmax_(class(Z)))) );
					else
						grid.(char(n+96)) = single(Z);
						grdutils(grid.(char(n+96)), sprintf('-M%f',1 / double(intmax_(class(Z)))) );	% Here OP is insitu
					end
				else
					grid.(char(n+96)) = single(Z);
				end
				N(i) = n;
			end
			clear Z

			for (i = 1:numel(k))
				tok = strtok(com(k(i)+1:end),' )+-*/^');		% Here we get the band's name
				kf = k(i) + length(tok);			% Find the position of last char of band's name
				com = [com(1:k(i)) 'grid.' char(N(i)+96) ' ' com(kf+1:end)];
				k = strfind(com,'&');	% Need because the remaining '&' may have changed position
			end
		end

		com = strrep(com,'&','');					% Remove the '&' characters

		if (~handles.IamCompiled)
			com_s = strrep(com,   '*', '.*');	com_s = strrep(com_s, '..*', '.*');		% Play safe
			com_s = strrep(com_s, '/', './');	com_s = strrep(com_s, '../', './');
			com_s = strrep(com_s, '^', '.^');	com_s = strrep(com_s, '..^', '.^');
			resp = eval(com_s);
		else
			com_s = strrep(com,   '.*', '*');	com_s = strrep(com_s, './', '/');
			com_s = strrep(com_s, '.^', '^');
			[resp, msg] = stalone(com_s, grid);
			if (~isempty(msg) || isempty(resp)),	semaforo(handles, 'green'),		end
			if     (~isempty(msg)),		errordlg(['ERROR: ' msg]),	return
			elseif (isempty(resp)),		errordlg('Programming error. Result is empty', 'Error'),	return
			end
		end
		semaforo(handles, 'green')

		if (numel(resp) > 1)			% 'resp' is a array. Construct a fake grid
			resp = single(resp);
			handMir = guidata(handles.hMirFig);
			if (handMir.image_type == 2)	% Non-referenced figs
				resp = flipud(resp);		% Need because grid is YDir direct
			end
			[zzz] = grdutils(resp,'-L');  z_min = zzz(1);     z_max = zzz(2);
			tmp.name = 'Computed_band';
			if (~isempty(handles.hMirFig))
				tmp.head = handMir.head;
				tmp.head(5:6) = [z_min z_max];
				mirone(resp, tmp, handles.hMirFig);		% Can fish cpt and proj stuff from parent
			else
				[m,n] = size(resp);
				tmp.X = 1:n;        tmp.Y = 1:m;
				tmp.head = [1 n 1 m z_min z_max 0 1 1];
				mirone(resp, tmp);
			end
		else							% Computations that do not involve grids
			txt = sprintf('%.10f',resp);
			while (txt(end) == '0')     % Remove trailing zeros
				txt(end) = [];
			end
			set(handles.edit_command,'String',txt)
		end
	catch
		semaforo(handles, 'green')
		errordlg(lasterr,'Error')
	end

% ------------------------------------------------------------------------
function str = move_operator(str)
% Make sure that operators are not "glued" to operands.
% But allow grid names with the '-' character.
% Also ensure that parentsis are not glued

	[str, did] = let_gridnames_have(str, 'forward');	% Temporarily replace +-()= by invisible chars
	k = strfind(str,')');
	if (k),		str = strrep(str,')',' ) ');    end
	k = strfind(str,'(');
	if (k),		str = strrep(str,'(',' ( ');    end
	k = strfind(str,'+');
	if (k),		str = strrep(str,'+',' + ');    end
	k = strfind(str,'-');
	if (k),		str = strrep(str,'-',' - ');    end
	k = strfind(str,'*');
	if (k),		str = strrep(str,'*',' * ');    end
	k = strfind(str,'/');
	if (k),		str = strrep(str,'/',' / ');    end
	k = strfind(str,'\');
	if (k),		str = strrep(str,'\',' \ ');    end
	if (did)	% Reset the invisible chars back to the originals +,-,(,),=
		str = let_gridnames_have(str, 'inv');
	end
	k = strfind(str,'. ');
	if (k)					% Here we want to have things like '.*' and not '. *'
        while (~isempty(k))
            str = [str(1:k(1)) str(min(k(1)+2,length(str)):end)];
            k = strfind(str,'. ');
        end
	end
	k = strfind(str,' '''); % Hard case of the transpose operator. It may fail often
	if (k)                  % Here we don't want to have things like ")'" turned into ") '"
        while (~isempty(k))
            str = [str(1:k(1)-1) str(min(k(1)+1,length(str)):end)];
            k = strfind(str,' ''');
        end    
	end
	k = strfind(str,'(');
	if (k),		str = strrep(str,'(',' ( ');	end
	k = strfind(str,')');
	if (k),		str = strrep(str,')',' ) ');	end
	% Now remove the extra blanks
	k = strfind(str,'  ');
	while (~isempty(k))
		str = strrep(str,'  ',' ');
		k = strfind(str,'  ');
	end

% ------------------------------------------------------------------------
function [str, did] = let_gridnames_have(str, dir)
% Let grid names have the  +,-,(,),= characters in the name. We do this by temporarily replacing
% them with char(1), char(2), ... char(5). 
% When DIR ~= 'forward' this function does the inverse operations

	did = false;
	ind = strfind(str,'&');
	if (isempty(ind)),	return,		end
	did = true;
	if (strncmpi(dir, 'for',3))
		for (k = 1:numel(ind))
			t = strtok(str(ind(k):end));
			t = strrep(t,'+',char(1));		t = strrep(t,'-',char(2));
			t = strrep(t,'(',char(3));		t = strrep(t,')',char(4));
			t = strrep(t,'=',char(5));
			str(ind(k) : ind(k)+numel(t)-1) = t;
		end
	else			% Do the inverse replacement. That is , restore original name
		str = strrep(str,char(1), '+');		str = strrep(str,char(2), '-');
		str = strrep(str,char(3), '(');		str = strrep(str,char(4), ')');
		str = strrep(str,char(5), '=');
	end

% ------------------------------------------------------------------------
function [out, msg] = stalone(comm, grid)
% Main dispatch function to be executed by the stand alone version that knows not about eval(????)
	comm = ddewhite(comm);
	while(strfind(comm, '  '))		% Replace any eventual consecutive blanks by a single one
		comm = strrep(comm,'  ', ' ');
	end

	ind = strfind(comm, 'grid.');
	s_names = sort(comm(ind+5));	% Store the original members names
	[out, msg] = run_inner(comm, grid, s_names);

% ------------------------------------------------------------------------
function [out, msg, grid, s_names] = run_inner(comm, grid, s_names)
% Execute the inner (or single) command and send back the result
% This function may be called several times (recursively)
	out = [];
	ind_l = strfind(comm, '(');		ind_r = strfind(comm, ')');
	if (numel(ind_l) ~= numel(ind_r))
		msg = 'Number of opening and closing parenthesis is not equal';		return
	end
	if (isempty(ind_l) && isempty(ind_r))			% Simplest case. A chain of basic operations
		[ops, msg] = sanitize_negs(comm);			% Split 'comm' in tokens and also try to address the 'MULL by -1' problem
		if (~isempty(msg)),		return,		end
		s = s_names;	% Shorter variable name
		[out, ops, grid, msg, s_names] = exec_n_consume(ops, '^', grid, s);		% Consume all highest priority operations (POWER)
		if (~isempty(msg) || numel(ops) == 1),		return,		end				% Either finished or error

		[out, ops, grid, msg, s_names] = exec_n_consume(ops, '*', grid, s);		% Consume second highest priority operations (MULL)
		if (~isempty(msg) || numel(ops) == 1),		return,		end

		[out, ops, grid, msg, s_names] = exec_n_consume(ops, '/', grid, s);		% Consume second highest priority operations (DIV)
		if (~isempty(msg) || numel(ops) == 1),		return,		end

		[out, ops, grid, msg, s_names] = exec_n_consume(ops, '+', grid, s);		% Consume lowest priority operations (SUM)
		if (~isempty(msg) || numel(ops) == 1),		return,		end

		[out, ops, grid, msg, s_names] = exec_n_consume(ops, '-', grid, s);		% Consume lowest priority operations (SUB)
		if (~isempty(msg) || numel(ops) == 1),		return,		end

		comm = move_operator(cat(2, ops{:}));		% If we reach here it means we had more than one ADD or SUB that must be
		[out, msg, grid, s_names] = run_inner(comm,grid,s_names);	% computed one by one because are not commutative. So restart.

	else											% We have '( ... )', consume their contents recursively
		ind_r = ind_r(ind_r > ind_l(end));			% The only one we care now are those to the right of last '('
		arg = ddewhite(comm(ind_l(end)+1:ind_r(1)-1));
		if (isempty(strfind(arg, ' ')))				% A single arg, must be a function arg. Find and run that function
			k = ind_l(end) - 2;						% Remember, it takes the form e.g. "sqrt ( arg )"
			while(k > 0 && comm(k) ~= ' '),	k = k - 1;	end
			fun = comm(k+1:ind_l(end) - 2);
			arg_n = str2double(arg);
			try
				trail = '';		msg = '';
				if (ind_r(1) < numel(comm)),	trail = comm(ind_r(1)+1:end);	end		% If last char is NOT a ')'
				if (~isnan(arg_n))					% A scalar, for example log( 5 )
					out = feval(fun, arg_n);
					comm = sprintf('%s%.20g%s',comm(1:k), out, trail);		% Put the result in the place of the function call
				else								% A matrix
					this_member = arg(6);
					if (~isa(grid.(this_member), 'double'))
					end
					FUN = upper(fun);
					ind = strcmp(FUN, {'EXP' 'LOG10' 'LOG' 'SIN' 'COS' 'TAN' 'SQRT'});
					if (any(ind))
						grid.(this_member) = grdutils(grid.(this_member), ['-O' FUN]);	% Reuse container
						out = grid.(this_member);		% Because this is this fun's output
					else
						out = single(feval(fun, double(grid.(this_member))));
						this_member = is_this_member_used(comm, this_member);	% Deal with tricky cases where we may need a new name
						grid.(this_member) = out;		% Reuse this container to hold the result just obtained
					end
					comm = sprintf('%sgrid.%s%s',comm(1:k), this_member, trail);	% Put the result in the place of the function call
				end
			catch
				msg = lasterr;		return
			end
			if (~isempty(strfind(comm, '(')) || ~isempty(strfind(comm, '*')) || ~isempty(strfind(comm, '/')) || ...
				~isempty(strfind(comm, '+')) || ~isempty(strfind(comm, '-')))		% Iterate
				[out, msg, grid, s_names] = run_inner(comm,grid,s_names);
			end

		else										% Multiple args. Exec the ( ... ) content
			[out, msg, grid, s_names] = run_inner(arg, grid, s_names);
			if (~isempty(msg)),		return,		end
			trail = '';
			if (ind_r(1) < numel(comm)),	trail = comm(ind_r(1)+1:end);	end		% If last char is NOT a ')'
			if (numel(out) == 1)					% A scalar
				comm = sprintf('%s%.20g%s',comm(1:ind_l(end)-1), out, trail);		% Put the result in the place of ( ... )
			else									% A matrix. A rather more complicated case
				ops = split_ops(arg);
				for (k = 1:numel(ops))
					if (strncmp(ops{k},'grid',3))	% Find the first grid.X operand, which holds the result
						this_member = ops{k}(6);
						this_member = is_this_member_used(comm, this_member);
						s_names = update_gridNames(s_names, comm, ind_l(end)-1, trail);		% Update the s_names unique names
						comm = sprintf('%sgrid.%s%s',comm(1:ind_l(end)-1), this_member, trail);
						break
					end
				end
			end
			if (numel(split_ops(comm)) > 1)			% Keep running until we have only one
				[out, msg, grid, s_names] = run_inner(comm, grid, s_names);
			end
		end
	end

% ------------------------------------------------------------------------
function [out, ops, grid, msg, s_names] = exec_n_consume(ops, op, grid, s_names)
% OPS is the cell string of tokens issued by SPLIT_OPS()
% OP is the operator. e.g. '*', '/', '+', '-' or '^'
% GRID is the structure holding all the array (when they are)
% At the end return the consumed cell OPS that will hold compute result when it was a scalar
% or update the content of the GRID.() member that was the first operand

	k = 2;		msg = '';	out = [];
	while (k < numel(ops))
		if (ops{k} == op)
			[out, msg] = do_mull_add(ops{k-1}, ops{k+1}, ops{k}, grid);
			if (~isempty(msg)),		return,		end
			if (numel(out) == 1)			% We had a scalar operation
				ops{k-1} = sprintf('%.20g', out);
			else							% A matrix op
				new_name = ops{k-1}(6);
				if (~isnan(str2double(new_name)) && ops{k+1}(1) == 'g')	% Scalar OP grid
					new_name = ops{k+1}(6);
				end
				ind = strfind(s_names, new_name);
				if (numel(ind) > 1)			% Struct name is taken because it occurs more than once. Need a new one
					new_name = char(double(s_names(end))+1);		% Pick the next letter in the alphabet
					s_names = sprintf('%s%s', s_names, new_name);	% Update the list of used names
				end
				grid.(new_name) = out;		% Reuse this container to hold the result just obtained
			end
			ops(k:k+1) = [];				% These ones are consumed, remove them
		else
			k = k + 2;
		end
		if (op == '+' || op == '-')			% Because they are not commutative, we have to compute one at at time
			break
		end
	end

% ------------------------------------------------------------------------
function ops = split_ops(comm)
% Split the COMM string in its tokens and return them in a cell array
	ind = [0 strfind(comm, ' ') numel(comm)+1];
	ops = cell(numel(ind)-1, 1);
	for (k = 1:numel(ind)-1)
		ops{k} = comm(ind(k)+1:ind(k+1)-1);
	end

% ------------------------------------------------------------------------
function [r, msg] = do_mull_add(a, b, op, grid)
% Do A +-/*^ B, where A and B can be scalar or matrices (on the constrain that it makes sense)
	msg = '';	r = [];
	a_n = str2double(a);	b_n = str2double(b);
	if (~isnan(a_n) && ~isnan(b_n))							% Trivial case of scalar operation
		if (op == '+'),		r = a_n + b_n;
		elseif (op == '-'),	r = a_n - b_n;
		elseif (op == '*'),	r = a_n * b_n;
		elseif (op == '/'),	r = a_n / b_n;
		else,				r = a_n ^ b_n;
		end
	elseif (strncmp(a,'grid',3) && ~isnan(b_n))				% add|sub|mull|div|pow(Matrix, b)
		if (~isa(grid.(a(6)), 'single')),	grid.(a(6)) = single(grid.(a(6)));	end
		if (op == '+' || op == '-' || op == '*' || op == '/')
			r = grdutils(grid.(a(6)), '-D');	% Deep copy because MULL & ADD are insitu only
		end
		if (op == '+'),		grdutils(r, sprintf('-A%f',b_n));		%r = grid.(a(6)) + b_n;
		elseif (op == '-'),	grdutils(r, sprintf('-A%f',-b_n));		%r = grid.(a(6)) - b_n;
		elseif (op == '*'),	grdutils(r, sprintf('-M%f',b_n));		%r = grid.(a(6)) * b_n;
		elseif (op == '/'),	grdutils(r, sprintf('-M%f',1/b_n));		%r = grid.(a(6)) / b_n;
		else,				r = double(grid.(a(6))) ^ b_n;	%r = grid.(a(6)) ^ b_n;
		end
	elseif (strncmp(a,'grid',3) && strncmp(b,'grid',3))		% add|sub|mull|div|pow(Matrix, Matrix)
		if (~isequal(size(grid.(a(6))), size(grid.(b(6))) ))
			msg = 'matrix dimension are not equal in Add, Sub or Power operator';	return
		end
		try
			if (~isa(grid.(a(6)), 'single')),	grid.(a(6)) = single(grid.(a(6)));	end
			if (~isa(grid.(b(6)), 'single')),	grid.(b(6)) = single(grid.(b(6)));	end
			if (op == '+'),		r = grdutils(grid.(a(6)), grid.(b(6)), '-OADD');	%r = grid.(a(6)) +  grid.(b(6));
			elseif (op == '-'),	r = grdutils(grid.(a(6)), grid.(b(6)), '-OSUB');	%r = grid.(a(6)) -  grid.(b(6));
			elseif (op == '*'),	r = grdutils(grid.(a(6)), grid.(b(6)), '-OMUL');	%r = grid.(a(6)) .* grid.(b(6));
			elseif (op == '/'),	r = grdutils(grid.(a(6)), grid.(b(6)), '-ODIV');	%r = grid.(a(6)) ./ grid.(b(6));
			else,				r = grdutils(grid.(a(6)), grid.(b(6)), '-OEXP');	%r = grid.(a(6)) .^ grid.(b(6));
			end
		catch,	msg = lasterr;			% Who knows what error
		end
	elseif (~isnan(a_n) && strncmp(b,'grid',3))				% add|sub|mull|div|pow(a, Matrix), some are odd but possible
		if (~isa(grid.(b(6)), 'single')),	grid.(b(6)) = single(grid.(b(6)));	end
		if (op == '+' || op == '-' || op == '*')
			r = grdutils(grid.(b(6)), '-D');	% Deep copy because MULL & ADD are insitu only
		end
		if (op == '-'),		grdutils(r, '-M-1');	op = '+';	end		% First by -1 then ADD 
		if (op == '+'),		grdutils(r, sprintf('-A%f',a_n));		%r = a_n + grid.(b(6));
		elseif (op == '*'),	grdutils(r, sprintf('-M%f',a_n));		%r = a_n * grid.(b(6));
		elseif (op == '/'),	r = a_n ./ double(grid.(b(6)));			%r = a_n ./ grid.(b(6));
		else,				r = a_n .^ double(grid.(b(6)));			%r = a_n .^ grid.(b(6));
		end
	end

% ------------------------------------------------------------------------
function [ops, msg] = sanitize_negs(comm)
% Try to deal with the problem of MULL/DIV by negative numbers
	msg = '';
	ops = split_ops(comm);
	if (rem(numel(ops), 2) == 1),	return,		end		% Done here

	ind = strfind(ops, '-');
	for (k = numel(ops):-1:1)
		if (~isempty(ind{k}))		% Found one '-'
			f = NaN;
			if (ind{k} == 1)		% First character in the command, check only if next is a number
				f = str2double(ops{k+1});
			elseif (ops{k-1} == '*' || ops{k-1} == '\' || ops{k-1} == '^' || ops{k-1} == '+' || ops{k-1} == '-')
				f = str2double(ops{k+1});
			end
			if (~isnan(f))			% Glue the '-' and the number that follows it and clear one 'ops' 
				ops{k} = ['-' ops{k+1}];
				ops(k+1) = [];
			end
		end
	end

	if (rem(numel(ops), 2) == 0),	msg = 'Wrong number of operators + operands';	end

% ------------------------------------------------------------------------
function this_member = is_this_member_used(comm, this_member)
% Check that a grid.X member name is free to hold an updated result or not
% It is not free when the grid.X appears more than once in the expression.
% In those cases we simply create a new, unique, grid.Y name
	ind = strfind(comm, ['grid.' this_member]);
	if (numel(ind) > 1)		% Yes, the member name is taken. Than we need to create a new one
		ind = strfind(comm, 'grid.');
		members = sort(comm(ind+5));
		this_member = char(double(members(end))+1);		% Pick the next letter in the alphabet
	end

% ------------------------------------------------------------------------
function s_names = update_gridNames(s_names, comm, n_prev, trail)
% Analise the chunk of comm that is going to be 'consumed' (between N_PREV and beginning of TRAIL)
% identify the grid.X names in it. Remove an occurrence of which from the names list S_NAMES
	ind = strfind(comm, trail);			% Find the end of the to-be-consumed chunk in COMM
	if (isempty(ind)),	ind = numel(comm);	end
	n_prev = max(1, n_prev);			% Due to indices gymnastic it can be zero sometimes
	chunk = comm(n_prev:ind);
	ind = strfind(chunk, 'grid.');
	if (~isempty(ind))					% It should never be but ...
		members = chunk(ind+5);
		for (k = numel(members):-1:1)
			ind = strfind(s_names, members(k));
			s_names(ind(end)) = '';
		end
	end

% ---
function semaforo(handles, cor)
% Change semaforo color. COR must be either 'red' or 'green'
	[img, pal] = aux_funs(['semaforo_' cor]);
	set(handles.push_semaforo, 'CData', ind2rgb8(img, pal))

% ------------------------------------------------------------------------
function push_help_CB(hObject, handles)
	str = sprintf(['This is mainly a grid calculator tool that operates on grids but\n'...
        'it can also be used in simple scalar calculations.\n\n'...
        'Grids are made available in two ways: (1) all grids loaded in different\n'...
        'Mirone figures show up in the list uppon start of this tool. (2) other\n'...
        'grids may be loaded using the "Load Grid" button.\n\n'...
        'A double click on a grid name moves it to the command zone, where its is\n'...
        'used as any other operand. Complicated expressions can be built, including\n'...
        'the use of Matlab commands. Note, however, that there is no error checking\n'...
        'so any comited error will show up as a Matlab error message, which may be\n'...
        'somewhat criptic.\n\n'...
        'You can also generate a grid, for example by entering the rand(200) command.\n'...
        'In those cases the generated grid will have coordinates 1:M by 1:N and a grid step of 1\n\n'...
		'WARNING: The above is holdsfor the MATLAB version but not all of it works with the\n'...
		'compiled (stand-alone) version.']);
	helpdlg(str,'Help')

% ------------------------------------------------------------------------
function imax = intmax_(varargin)
% Local copy because the one in ML doesn't let be compiled

	classname = varargin{1};
	switch (classname)
		case 'int8',			imax = int8(127);
		case 'uint8',			imax = uint8(255);
		case 'int16',			imax = int16(32767);
		case 'uint16',			imax = uint16(65535);
		case 'int32',			imax = int32(2147483647);
		case 'uint32',			imax = uint32(4294967295);
		case 'int64'
			% Turn off the overflow warning because the constant value below
			% is interpreted as the closest double value, namely 2^63.
			% This is out of range and warns of overflow upon conversion to int64.
			s = warning('query','MATLAB:intConvertOverflow');
			warning('off','MATLAB:intConvertOverflow');
			imax = int32(9223372036854775807);      % R13 has no int64
			warning(s);
		case 'uint64'
			% Turn off the overflow warning because the constant value below
			% is interpreted as the closest double value, namely 2^64.
			% This is out of range and warns of overflow upon conversion to uint64.
			s = warning('query','MATLAB:intConvertOverflow');
			warning('off','MATLAB:intConvertOverflow');
			imax = uint32(9223372036854775807);      % R13 has no uint64
			warning(s);
		otherwise
			error('MATLAB:intmax_:invalidClassName','Invalid class name.')
	end

% --- Creates and returns a handle to the GUI figure. 
function grid_calculator_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'), 'Position',[520 602 669 260],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Grid calculator',...
'NumberTitle','off',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 190 651 61],...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Max',3,...
'Style','edit',...
'Tooltip','Enter here a Matlab valid command',...
'Tag','edit_command');

uicontrol('Parent',h1, 'Position',[10 37 311 135],...
'BackgroundColor',[1 1 1],...
'Callback',@grid_calculator_uiCB,...
'Style','listbox',...
'Tooltip','Names of currently available arrays',...
'Value',1,'Tag','listbox_inArrays');

uicontrol('Parent',h1, 'Position',[341 150 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','1','Tag','push_digit');

uicontrol('Parent',h1, 'Position',[374 150 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','2','Tag','push_digit');

uicontrol('Parent',h1, 'Position',[407 150 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','3','Tag','push_digit');

uicontrol('Parent',h1, 'Position',[341 117 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','4','Tag','push_digit');

uicontrol('Parent',h1, 'Position',[374 117 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','5','Tag','push_digit');

uicontrol('Parent',h1, 'Position',[407 117 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','6','Tag','push_digit');

uicontrol('Parent',h1, 'Position',[341 84 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','7','Tag','push_digit');

uicontrol('Parent',h1, 'Position',[374 84 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','8','Tag','push_digit');

uicontrol('Parent',h1, 'Position',[407 84 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','9','Tag','push_digit');

uicontrol('Parent',h1, 'Position',[341 51 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','0','Tag','push_digit');

uicontrol('Parent',h1, 'Position',[374 51 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','.','Tag','push_digit');

uicontrol('Parent',h1, 'Position',[407 51 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','=','Tag','push_equal');

uicontrol('Parent',h1, 'Position',[440 150 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','/','Tag','push_op');

uicontrol('Parent',h1, 'Position',[440 117 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','*','Tag','push_op');

uicontrol('Parent',h1, 'Position',[440 84 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','-','Tag','push_op');

uicontrol('Parent',h1, 'Position',[440 51 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','+','Tag','push_op');

uicontrol('Parent',h1, 'Position',[11 8 111 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','Load Grid',...
'TooltipString','Load ONLY GMT or Surfer grids',...
'Tag','push_loadGrid');

uicontrol('Parent',h1, 'Position',[491 150 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','sin','Tag','push_funs');

uicontrol('Parent',h1, 'Position',[551 150 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','cos','Tag','push_funs');

uicontrol('Parent',h1, 'Position',[611 150 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','tan','Tag','push_funs');

uicontrol('Parent',h1, 'Position',[491 117 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','log10','Tag','push_funs');

uicontrol('Parent',h1, 'Position',[551 117 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','log e','Tag','push_funs');

uicontrol('Parent',h1, 'Position',[611 117 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','exp','Tag','push_funs');

uicontrol('Parent',h1, 'Position',[491 84 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','sqrt','Tag','push_funs');

uicontrol('Parent',h1, 'Position',[550 84 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','abs','Tag','push_funs');

uicontrol('Parent',h1, 'Position',[610 84 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String','(','Tag','push_leftPar');

uicontrol('Parent',h1, 'Position',[637 84 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'String',')','Tag','push_rightPar');

uicontrol('Parent',h1, 'Position',[491 51 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'Tooltip', 'Compute Bright temperature in Celsius', ...
'Vis', 'off', ...
'String','Trad','Tag','push_Trad');

uicontrol('Parent',h1, 'Position',[551 51 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'Tooltip', 'Compute Radiance at Top of Atmosphere', ...
'Vis', 'off', ...
'String','R(toa)','Tag','push_Trad');

uicontrol('Parent',h1, 'Position',[611 51 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',9,...
'Tooltip', 'Compute Reflectance at Top of Atmosphere', ...
'Vis', 'off', ...
'String','Rho(toa)','Tag','push_Trad');

uicontrol('Parent',h1, 'Position',[491 23 50 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',8,...
'Tooltip', 'Compute "At Surface Reflectance" using DOS2 (COST) method', ...
'Vis', 'off', ...
'String','Rho(surf)','Tag','push_Trad');

uicontrol('Parent',h1, 'Position',[469 5 17 49], 'Enable','inactive', 'Tag','push_semaforo');

uicontrol('Parent',h1, 'Position',[221 7 23 23],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'String','?','Tag','push_help');

uicontrol('Parent',h1, 'Position',[569 6 91 21],...
'Callback',@grid_calculator_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','Compute','Tag','push_compute');

function grid_calculator_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

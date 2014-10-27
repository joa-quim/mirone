function varargout = grid_calculator(varargin)
% An array calculator 
%
% The big trick here is to store the arrays in structure members of a 'grid' struct
% container and at same time maintain an updated map between those struct members and
% variables char strings with the same name. This, and lots of parsings, allow not 
% only to run the whole computation as an eval('command') from within Matlab but also
% to break that 'command' in its tokens and, now with lots of parsings, execute them
% in a way that also works in the stand-alone (compiled) version.
%
% This stil does all computations in doubles

%	Copyright (c) 2004-2014 by J. Luis
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
	grid_calculator_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	old_vis = get(0,'ShowHiddenHandles');
	set(0,'ShowHiddenHandles','on')
	h_figs = findobj('Type','figure');
	set(0,'ShowHiddenHandles',old_vis)

	handles.grid_patos = [];    % To hold paths of eventual future loaded grids
	handles.loaded_grid = [];   % To hold names of eventual future loaded grids
	handles.name_str = [];
	handles.BL = [];

	if (~isempty(varargin) && ishandle(varargin{1}))	% Not tested if varargin{1} is a Fig handle
		handMir = guidata(varargin{1});
		handles.home_dir = handMir.home_dir;		handles.last_dir = handMir.last_dir;
		handles.work_dir = handMir.work_dir;		handles.path_tmp = handMir.path_tmp;
		handles.IamCompiled = handMir.IamCompiled;
		BL = getappdata(varargin{1},'BandList');
		if (~isempty(BL))
            for (i = 1:ndims(BL{2}))
                handles.name_str{i} = ['Band_' num2str(i)];
            end
            h_figs = [];			% Tear off the net for the 'h_figs' fishing
            handles.BL = BL{2};
		end
	elseif (~isempty(varargin) && isa(varargin{1}, 'char'))
		do_tests;
		delete(hObject)
		return
	else
		handles.home_dir = cd;		handles.last_dir = cd;		handles.work_dir = cd;
		handles.path_tmp = [handles.home_dir '/tmp'];		handles.IamCompiled = false;
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
				Z = [];        
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
            n = n + 1;
        end
	end

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
function listbox_inArrays_CB(hObject, handles)
% if this is a doubleclick, copy the array's name to the edit box
	if strcmp(get(gcbf,'SelectionType'),'open')
        str = get(hObject,'String');
        sel = get(hObject,'Value');
        com = get(handles.edit_command,'String');
        set(handles.edit_command,'String', [com ' &' str{sel}])
	end

% ------------------------------------------------------------------------
function push_1_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '1'])

% ------------------------------------------------------------------------
function push_2_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '2'])

% ------------------------------------------------------------------------
function push_3_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '3'])

% ------------------------------------------------------------------------
function push_4_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '4'])

% ------------------------------------------------------------------------
function push_5_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '5'])

% ------------------------------------------------------------------------
function push_6_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '6'])

% ------------------------------------------------------------------------
function push_7_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '7'])

% ------------------------------------------------------------------------
function push_8_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '8'])

% ------------------------------------------------------------------------
function push_9_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '9'])

% ------------------------------------------------------------------------
function push_0_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '0'])

% ------------------------------------------------------------------------
function push_dot_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '.'])

% ------------------------------------------------------------------------
function push_equal_CB(hObject, handles)
	%set(handles.edit_command,'String', [get(handles.edit_command,'String') ' = '])
	push_compute_CB(hObject, handles)

% ------------------------------------------------------------------------
function push_devide_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' / '])

% ------------------------------------------------------------------------
function push_mull_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' * '])

% ------------------------------------------------------------------------
function push_minus_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' - '])

% ------------------------------------------------------------------------
function push_plus_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' + '])

% ------------------------------------------------------------------------
function push_loadGrid_CB(hObject, handles)
% This function doesn't realy loads the grid. It only stores the grid name.
% True loading is donne in "Compute"
	str1 = {'*.grd;*.GRD;*.nc;*.NC', 'Grid files (*.grd,*.GRD,*.nc,*.NC)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles, str1,'Select GMT grid','get');
	if isequal(FileName,0),		return,		end
	str = get(handles.listbox_inArrays,'String');
	str{end+1} = FileName;
	set(handles.listbox_inArrays,'String',str);
	handles.grid_patos{end+1} = PathName;
	handles.loaded_grid{end+1} = FileName;
	guidata(hObject,handles)

% ------------------------------------------------------------------------
function push_sin_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' sin( '])

% ------------------------------------------------------------------------
function push_cos_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' cos( '])

% ------------------------------------------------------------------------
function push_tan_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' tan( '])

% ------------------------------------------------------------------------
function push_log10_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' log10( '])

% ------------------------------------------------------------------------
function push_log_e_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' log( '])

% ------------------------------------------------------------------------
function push_exp_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' exp( '])

% ------------------------------------------------------------------------
function push_sqrt_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' sqrt( '])

% ------------------------------------------------------------------------
function push_abs_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' abs( '])

% ------------------------------------------------------------------------
function push_leftPar_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '( '])

% ------------------------------------------------------------------------
function push_rightPar_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' )'])

% ------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
	com = get(handles.edit_command,'String');
	if (isempty(com)),	return,		end

	if (~isempty(handles.BL))		% We are dealing with Bands arithmetics
		bandArithm(handles, com);	return
	end

	% Those start at 2 because they are meant to be used only when grid names apear repeatedly
	in_g_count = 2;     out_g_count = 2;

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
						tmp.head = handtmp.head;
						grid_t = double(grid_t);      % grdread_m allways outpus singles
					else
						grid_t = double(getappdata(hand_fig.figure1,'dem_z'));
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
						grid_t = double(grid_t);      % grdread_m allways outpus singles
					else
						grid_t = double(getappdata(hand_fig.figure1,'dem_z'));
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
		end

		com = strrep(com,'&','');					% Remove the '&' characters
		com = strrep(com,'e - ','e-');				% base 10 numbers cannot have those spaces
		com = strrep(com,'e + ','e+');

		if (~handles.IamCompiled)
			com_s = strrep(com,   '*', '.*');	com_s = strrep(com_s, '..*', '.*');		% Play safe
			com_s = strrep(com_s, '/', './');	com_s = strrep(com_s, '../', './');
			com_s = strrep(com_s, '^', '.^');	com_s = strrep(com_s, '..^', '.^');
			resp = eval(com_s);
		else
			com_s = strrep(com,   '.*', '*');	com_s = strrep(com_s, './', '/');
			com_s = strrep(com_s, '.^', '^');
			[resp, msg] = stalone(com_s, grid);
			if (isempty(resp)),			errordlg('Programming error. Result is empty', 'Error'),	return
			elseif (~isempty(msg))		errordlg(['ERROR: ' msg]),	return
			end
		end

		if (~isempty(k))					% We had grids in input
			tmp.head(5:6) = [min(resp(:)) max(resp(:))];
			tmp.name = 'Computed_grid';
			mirone(single(resp),tmp);
		elseif (numel(resp) > 1)			% 'resp' is a array. Construct a fake grid
			[m,n] = size(resp);
			resp = single(resp);
			tmp.X = 1:n;        tmp.Y = 1:m;
			[zzz] = grdutils(resp,'-L');  z_min = zzz(1);     z_max = zzz(2);
			tmp.head = [1 n 1 m z_min z_max 0 1 1];
			tmp.name = 'Computed_array';
			mirone(resp,tmp);
		else								% Computations that do not involve grids
			txt = sprintf('%.10f',resp);
			while (txt(end) == '0')			% Remove trailing zeros
				txt(end) = [];
			end
			set(handles.edit_command,'String',txt)
		end
	catch
		errordlg(lasterr,'Error')
	end

% ------------------------------------------------------------------------
function bandArithm(handles, com)
% ...
	com = move_operator(com);			% Make sure operators are not "glued" to operands
	k = strfind(com,'&');

	try									% Wrap it here to report any possible error
		if (~isempty(k))				% We have grids
			N = zeros(1, numel(k));
			for (i = 1:numel(k))		% Loop over bands
				tok = strtok(com(k(i)+1:end));
				n = find(strcmp(tok,handles.name_str));		% n ~= [] when it is in memory                
				grid.(char(n+96)) = double(handles.BL(:,:,n));
				N(i) = n;
			end

			for (i = 1:numel(k))
				tok = strtok(com(k(i)+1:end));      % Here we get the band's name
				kf = k(i) + length(tok);            % Find the position of last char of band's name
				com = [com(1:k(i)) 'grid.' char(N(i)+96) ' ' com(kf+1:end)];
			end
		end

		com = strrep(com,'&','');                   % Remove the '&' characters

		if (~handles.IamCompiled)
			com_s = strrep(com,   '*', '.*');	com_s = strrep(com_s, '..*', '.*');		% Play safe
			com_s = strrep(com_s, '/', './');	com_s = strrep(com_s, '../', './');
			com_s = strrep(com_s, '^', '.^');	com_s = strrep(com_s, '..^', '.^');
			resp = eval(com_s);
		else
			com_s = strrep(com,   '.*', '*');	com_s = strrep(com_s, './', '/');
			com_s = strrep(com_s, '.^', '^');
			[resp, msg] = stalone(com_s, grid);
			if (isempty(resp)),			errordlg('Programming error. Result is empty', 'Error'),	return
			elseif (~isempty(msg))		errordlg(['ERROR: ' msg]),	return
			end
		end

		if (numel(resp) > 1)			% 'resp' is a array. Construct a fake grid
			if (max(resp(:)) <= 1),      resp = single(resp * 255);      end
			[m,n] = size(resp);
			tmp.X = 1:n;        tmp.Y = 1:m;
			[zzz] = grdutils(resp,'-L');  z_min = zzz(1);     z_max = zzz(2);
			tmp.head = [1 n 1 m z_min z_max 0 1 1];
			tmp.name = 'Computed_band';
			resp = uint8(resp);
			mirone(resp,tmp);
		else							% Computations that do not involve grids
			txt = sprintf('%.10f',resp);
			while (txt(end) == '0')     % Remove trailing zeros
				txt(end) = [];
			end
			set(handles.edit_command,'String',txt)
		end
	catch
		errordlg(lasterr,'Error')
	end

% ------------------------------------------------------------------------
function str = move_operator(str)
	% Make sure that operators are not "glued" to operands.
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
				trail = '';
				if (ind_r(1) < numel(comm)),	trail = comm(ind_r(1)+1:end);	end		% If last char is NOT a ')'
				if (~isnan(arg_n))					% A scalar
					out = feval(fun, arg_n);
					comm = sprintf('%s%.20g%s',comm(1:k), out, trail);		% Put the result in the place of the function call
				else								% A matrix
					this_member = arg(6);
					out = feval(fun, grid.(this_member));
					this_member = is_this_member_used(comm, this_member);	% Deal with tricky cases where we may need a new name
					grid.(this_member) = out;		% Reuse this container to hold the result just obtained
					s_names = update_gridNames(s_names, comm, k, trail);			% Update the s_names unique names
					comm = sprintf('%sgrid.%s%s',comm(1:k), this_member, trail);	% Put the result in the place of the function call
				end
			catch
				msg = lasterr;		return
			end
			[out, msg, grid, s_names] = run_inner(comm,grid,s_names);	% Now recursively call this fun untill all were consumed

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
		elseif (op == '-')	r = a_n - b_n;
		elseif (op == '*')	r = a_n * b_n;
		elseif (op == '/')	r = a_n / b_n;
		else				r = a_n ^ b_n;
		end
	elseif (strncmp(a,'grid',3) && ~isnan(b_n))				% add|sub|mull|div|pow(Matrix, b)
		if (op == '+'),		r = grid.(a(6)) + b_n;
		elseif (op == '-')	r = grid.(a(6)) - b_n;
		elseif (op == '*')	r = grid.(a(6)) * b_n;
		elseif (op == '/')	r = grid.(a(6)) / b_n;
		else				r = grid.(a(6)) ^ b_n;
		end
	elseif (strncmp(a,'grid',3) && strncmp(b,'grid',3))		% add|sub|mull|div|pow(Matrix, Matrix)
		if (~isequal(size(grid.(a(6))), size(grid.(a(6))) ))
			msg = 'matrix dimension are not equal in Add, Sub or Power operator';	return
		end
		try
			if (op == '+'),		r = grid.(a(6)) +  grid.(b(6));
			elseif (op == '-')	r = grid.(a(6)) -  grid.(b(6));
			elseif (op == '*')	r = grid.(a(6)) .* grid.(b(6));
			elseif (op == '/')	r = grid.(a(6)) ./ grid.(b(6));
			else				r = grid.(a(6)) .^ grid.(b(6));
			end
		catch,	msg = lasterr;			% Quite likely a data type operation error
		end
	elseif (~isnan(a_n) && strncmp(b,'grid',3))				% add|sub|mull|div|pow(a, Matrix), some are odd but possible
		if (op == '+'),		r = a_n + grid.(b(6));			% A bit idiot, but allowed
		elseif (op == '-')	r = a_n - grid.(b(6));			%		""
		elseif (op == '*')	r = a_n * grid.(b(6));
		elseif (op == '/')	r = a_n ./ grid.(b(6));
		else				r = a_n .^ grid.(b(6));
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

% --- Creates and returns a handle to the GUI figure. 
function grid_calculator_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Grid calculator',...
'NumberTitle','off',...
'Position',[520 602 669 198],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Max',3,...
'Position',[10 147 651 41],...
'Style','edit',...
'Tooltip','Enter here any Matlab valid command',...
'Tag','edit_command');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@grid_calculator_uiCB,...
'Position',[10 37 311 91],...
'Style','listbox',...
'Tooltip','Names of currently available arrays',...
'Value',1,'Tag','listbox_inArrays');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[341 105 23 23],...
'String','1','Tag','push_1');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[374 105 23 23],...
'String','2','Tag','push_2');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[407 105 23 23],...
'String','3','Tag','push_3');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[341 72 23 23],...
'String','4','Tag','push_4');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[374 72 23 23],...
'String','5','Tag','push_5');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[407 72 23 23],...
'String','6','Tag','push_6');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[341 39 23 23],...
'String','7','Tag','push_7');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[374 39 23 23],...
'String','8','Tag','push_8');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[407 39 23 23],...
'String','9','Tag','push_9');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[341 6 23 23],...
'String','0','Tag','push_0');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[374 6 23 23],...
'String','.','Tag','push_dot');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[407 6 23 23],...
'String','=','Tag','push_equal');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[440 105 23 23],...
'String','/','Tag','push_devide');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[440 72 23 23],...
'String','*','Tag','push_mull');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[440 39 23 23],...
'String','-','Tag','push_minus');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[440 6 23 23],...
'String','+','Tag','push_plus');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[11 8 111 21],...
'String','Load Grid',...
'TooltipString','Load ONLY  GMT or Surfer grids',...
'Tag','push_loadGrid');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[491 105 50 21],...
'String','sin','Tag','push_sin');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[551 105 50 21],...
'String','cos','Tag','push_cos');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[611 105 50 21],...
'String','tan','Tag','push_tan');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[491 72 50 21],...
'String','log10','Tag','push_log10');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[551 72 50 21],...
'String','log e','Tag','push_log_e');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[611 72 50 21],...
'String','exp','Tag','push_exp');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[491 39 50 21],...
'String','sqrt','Tag','push_sqrt');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[550 39 50 21],...
'String','abs','Tag','push_abs');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[610 39 23 23],...
'String','(','Tag','push_leftPar');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[637 39 23 23],...
'String',')','Tag','push_rightPar');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[589 6 71 21],...
'String','Compute','Tag','push_compute');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[221 7 23 23],...
'String','?','Tag','push_help');

function grid_calculator_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

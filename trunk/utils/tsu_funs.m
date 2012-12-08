function tsu_funs(opt,varargin)
% Helper function to do Tsunami related computations

%	Copyright (c) 2004-2012 by J. Luis
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

	switch opt
		case 'SwanCompute'
			SwanCompute(varargin{:})
		case 'SwanGridBorder'
			SwanGridBorderStations(varargin{:})
		case 'Tsun2'
			Tsun2Compute(varargin{:})
		case 'TTT'
			TTT(varargin{:})
	end
	
% --------------------------------------------------------------------
function TTT(handles,opt)
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	if (nargin == 1),   opt = [];   end
	if isempty(opt)							% Plot a point source
		pt = click_e_point(1,'crosshair');
		h = line('XData',pt(1,1),'YData',pt(1,2), 'Parent',handles.axes1,'Marker','o', ...
			'MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10,'Tag','TTT');
		draw_funs(h,'DrawSymbol')			% Set symbol's uicontextmenu
	elseif (strcmp(opt,'line'))				% Draw a 'Fault source'. Each vertex will be a punctual source.
		[xp,yp] = getline_j(handles.figure1);
		h = line('XData', xp, 'YData', yp,'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,'Tag','TTT');
		draw_funs(h,'line_uicontext')
	elseif (strcmp(opt,'load'))
		str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)'};
		[FileName,PathName] = put_or_get_file(handles,str1,'Select input xy_time file name','get');
		if isequal(FileName,0),		return,		end
		out = load_xyz(handles,[PathName FileName]);
		if (size(out,2) ~= 3)
			errordlg('Wrong choice. For using this option the file MUST have 3 columns (position and time).','Error'); return
		end
		h = line('XData',out(:,1),'YData',out(:,2), 'Parent',handles.axes1, 'Marker','o', ...
			'MarkerFaceColor','y','linestyle','none', 'MarkerEdgeColor','k','MarkerSize',10,'Tag','TTT','UserData',out(:,3));
		setappdata(h,'TTTimes',out)
		draw_funs(h,'DrawSymbol')			% Set symbol's uicontextmenu    
	else					% Compute
		h_src = findobj(handles.axes1,'Type','line','Tag','TTT');
		if (isempty(h_src)),    errordlg('Yes I compute, but ... WHAT?','Error'),	return,		end
		if (numel(h_src) > 1)
			errordlg('More than one source found. This is not allowed neither in single source or Ray tracing modes. Be modest.','Error');     return
		end
		xy_t = getappdata(h_src,'TTTimes');				% See if we have anything in appdata
		if (isempty(xy_t)),		single_src = 1;			% Single source mode
		else					single_src = 0;			% Ray tracing mode
		end

		[X,Y,Z,head] = load_grd(handles);
		if isempty(Z),   return,	end					% An error message was already issued
		if (handles.have_nans),     errordlg('Bathymetry grid cannot have NaNs','Error');   return;    end
		h_info = [head(1:4) head(8:9)];
		set(handles.figure1,'pointer','watch');
		aguentabar('title','Hold on your camels: computing solution')
		if (single_src)			% Single source mode
			xx = get(h_src,'XData');    yy = get(h_src,'YData');
			tt = wave_travel_time(Z,h_info,[xx(1) yy(1)],handles.geog);
			for (k = 2:numel(xx))
				aguentabar((k-1)/numel(xx))
				tt_ = wave_travel_time(Z,h_info,[xx(k) yy(k)],handles.geog);
				tt  = min(tt, tt_);
			end
			aguentabar(1)
			tit = 'Tsunami Travel Times';
		else					% Find ray tracing solution
			xx = get(h_src,'XData');
			yy = get(h_src,'YData');
			tempo = get(h_src,'UserData');
			if (numel(xx) ~= numel(tempo))				% Some(s) station(s) has been killed
				[c,ia,ib] = setxor(xx(:),xy_t(:,1));	% Find which
				tempo(ib) = [];							% Remove the corresponding times
			end
			tmp = 0;
			for (i = 1:numel(xx))
				tt = double(wave_travel_time(Z, h_info, [xx(i) yy(i)], handles.geog));	% Ghrrr
				tmp = tmp + (tt - tempo(i)) .^2;
				aguentabar(i/numel(xx))
			end
			tt = sqrt(tmp/size(xy_t,1));	tit = 'Ray tracing solution';	clear tmp;
		end
		zz = grdutils(tt,'-L');
		head(5:6) = double(zz(1:2));
		tmp.X = X;		tmp.Y = Y;		tmp.head = head;	tmp.name = tit;
		set(handles.figure1,'pointer','arrow');
		mirone(tt, tmp);
	end

% --------------------------------------------------------------------
function SwanCompute(handles)
% The following are flags to signal how much information is already known before
% calling "swan_options". They will be set to 1 by the guessing code. 
% I don't try to guess a "deform_only".
	small = 1e-6;   % Used in postion comparation. It places the accuracy at sub-meter level
	Z_bat = [];     head_bat = [];      Z_src = [];     head_src = [];
	haveBat = 0;	haveSource = 0;
	nothing_in_memory = 0;      bat_and_deform_with_maregs = 0;
	bat_and_deform = 0;         bat_with_maregs = 0;	bat_only = 0;

	if (handles.no_file || ~handles.validGrid)
		nothing_in_memory = 1;
	elseif (strncmp(get(handles.figure1,'Name'), 'Okada deformation', 17) )
		% Simplest guessing case. We have them both here
		haveBat = 1;    haveSource = 1;
		[X,Y,Z_src,head_src] = load_grd(handles);
		if isempty(Z_src),   return,	end			% An error message was already issued
		hFigParent = getappdata(handles.figure1,'hFigParent');
		Z_bat  = getappdata(hFigParent,'dem_z');
		handBat = guidata(hFigParent);		head_bat = handBat.head;
	else
		% We have to assume that this is a bat grid
		[X,Y,Z_bat,head_bat] = load_grd(handles);
		if isempty(Z_bat),   return,	end
		haveBat = 1;
		hFigs = findobj('Type','figure');		h_src = [];
		for (k = 1:numel(hFigs))
			if (strncmp(get(hFigs(k),'Name'), 'Okada deformation', 17) )
				h_src = hFigs(k);
				break
			end
		end
		if (~isempty(h_src))
			Z_src = getappdata(h_src,'dem_z');		handSrc = guidata(h_src);
			head_src = handSrc.head;
			haveSource = 1;
		else
			bat_only = 1;
		end
	end

	% If have both test that they cover the same region.
	if (haveBat && haveSource)		% If bat & source figures test that they are compatible.
		difes = head_src(1:4) - head_bat(1:4);
		if (any(abs(difes) > small))
			msg{1} = 'Bathymetry & Source grids do not cover the same region';
			msg{2} = ['x_min diff = ' num2str(difes(1))];		msg{3} = ['x_max diff = ' num2str(difes(2))];
			msg{4} = ['y_min diff = ' num2str(difes(3))];		msg{5} = ['y_max diff = ' num2str(difes(4))];
			errordlg(msg','Error'),     return
		end
		if ( numel(Z_bat) ~= numel(Z_src) )
			errordlg('Bathymetry and deformation grids have not the same size.','Error'),   return;
		end
		bat_and_deform = 1;
	end

	% See if we have maregraphs (if they exist, that is interpreted as a computation request)
	h_mareg = findobj('Type','line','Tag','Maregraph');
	if (~isempty(h_mareg))
		mareg_pos = [];
		side = zeros(1,numel(h_mareg));
		for (i = 1:numel(h_mareg)),		side(i) = isappdata(h_mareg(i),'Side');    end
		if (any(side) && ~all(side))
			errordlg('ERROR: you cannot mix individual stations with stations on grid borders.','ERROR');   return
		end
		if (all(side))		% We are in the stations on grid borders case
			side = cell(1,numel(h_mareg));
			for (i=1:numel(h_mareg)),		side{i} = getappdata(h_mareg(i),'Side');    end
			side = strrep(side,'W','1');	side = strrep(side,'S','2');
			side = strrep(side,'E','3');	side = strrep(side,'N','4');
			[i,j] = sort(side);				% Now we can sort them to find the correct WSEN order
			h_mareg = h_mareg(j);			% And reorder the handles to follow the expected order
		end
		for (i = 1:numel(h_mareg))
			mareg_pos = [mareg_pos; [get(h_mareg(i),'XData')' get(h_mareg(i),'YData')']];
		end
		if (bat_and_deform)
			bat_and_deform_with_maregs = 1;
		else					% OK, this also implies that we have a possibly valid bat file
			bat_with_maregs = 1;
		end
	end

	% Check if there is rectangle on the bat image. If yes, output grids will be inside that region only
	hLine = findobj('Type','line');
	if (~isempty(hLine) && ~isempty(h_mareg))
		hLine = setxor(hLine, h_mareg);
	end
	got_it = 0;
	if (~isempty(hLine))
		for (k = 1:numel(hLine))
			if ( check_IsRectangle(hLine(k)) )
				got_it = k;			% Found one. That's do one we gonna use
				break
			end
		end
	end
	opt_R = ' ';
	if (got_it)
		x = get(hLine(got_it),'XData');   y = get(hLine(got_it),'YData');
		opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g',min(x), max(x), min(y), max(y));
	end
	% -----------------------------------------------------------------------------
	
	if (bat_and_deform_with_maregs)
        out = swan_options(handles,'bat_and_deform_with_maregs',mareg_pos);
	elseif (bat_with_maregs)
        out = swan_options(handles,'bat_with_maregs',Z_bat,head_bat,mareg_pos);
	elseif (bat_and_deform)
        out = swan_options(handles,'bat_and_deform');
	elseif (bat_only)
        out = swan_options(handles,'bat_only',Z_bat,head_bat);
	else
        out = swan_options(handles);
        nothing_in_memory = 1;
	end
	pause(0.05);			% Give time to swan_options window to die

	if (isempty(out)),		return,		end

	if (nothing_in_memory)
        Z_bat = out.grid_Z_bat;     head_bat = out.grid_head_bat;
        Z_src = out.grid_Z_src;     head_src = out.grid_head_src;
	end
	if (bat_with_maregs || bat_only)
        Z_src = out.grid_Z_src;     head_src = out.grid_head_src;
	end

	opt_O = ' ';   opt_M = ' ';   opt_N = ' ';   opt_m = [];	opt_G = ' ';   opt_S = ' ';		opt_s = ' ';	opt_J = ' ';
	if (isfield(out,'maregraph_xy') && isfield(out,'maregraph_data_name'))
        opt_O = ['-O' out.maregraph_data_name];
	end
	if (isfield(out,'opt_M')),    opt_M = out.opt_M;   end
	if (isfield(out,'opt_N')),    opt_N = out.opt_N;   end
	if (isfield(out,'opt_m')),    opt_m = out.opt_m;   end
	if (isfield(out,'opt_G')),    opt_G = out.opt_G;   end
	if (isfield(out,'opt_S')),    opt_S = out.opt_S;   end
	if (isfield(out,'opt_s')),    opt_s = out.opt_s;   end
	if (isfield(out,'opt_J')),    opt_J = out.opt_J;   end

	if (isempty(Z_bat) || isempty(head_bat) || isempty(Z_src) || isempty(head_src))
        errordlg('ERROR: one or more of the bat/source variables are empty where they souldn''t be.','Error');
        return
	end

	Z_bat = double(Z_bat);      Z_src = double(Z_src);      % make sure they are both doubles
	if (handles.IamCompiled)	opt_e = '-e';
	else						opt_e = '';
	end

	% Make sure we start with zero water on land 
	Z_src(Z_bat > 0) = 0;

	if (isfield(out,'maregraph_xy'))	% Ask for computation of maregraphs
        if (~isempty(opt_m))			% Movie option
            tmovie = swan(Z_bat, head_bat, Z_src, head_src, out.params, out.maregraph_xy, opt_O, ...
                    opt_M, opt_N, opt_G, '-f', opt_R, opt_J, opt_e);
        else
            swan(Z_bat, head_bat, Z_src, head_src, out.params, out.maregraph_xy, opt_O, ...
                    opt_M, opt_N, opt_G, opt_S, opt_s, opt_R, opt_J, opt_e);
        end    
	else								% Compute grids or movie
        if (~isempty(opt_m))			% Movie option
            tmovie = swan(Z_bat, head_bat, Z_src, head_src, out.params, opt_M, opt_N, opt_G, ...
					'-f', opt_S, opt_s, opt_R, opt_J, opt_e);
        else
            swan(Z_bat, head_bat, Z_src, head_src, out.params, opt_M, opt_N, opt_G, opt_S, ...
					opt_s, opt_R, opt_J, opt_e);
        end
	end
	if (isfield(out,'opt_m')),   do_movie(handles,tmovie,'swan');   end

% --------------------------------------------------------------------
function SwanGridBorderStations(handles)
% Get the limits of smaller grid and plot them as a rectangle (individual lines) on current fig
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	small = 1e-5;       % Used in relative origin comparation.
	str_R = [];         str_I =[];
	adjust_w = 0;       adjust_s = 0;       % Flags used when the finer grid has to be adjusted
	adjust_x_inc = 0;   adjust_y_inc = 0;   % in order to fit correctly with coarser grid
	in_w_to_e = 0;      in_s_to_n = 0;

	[X,Y,Z,head] = load_grd(handles);
	if isempty(Z),		return,		end;    % An error message was already issued
	str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select GMT grid','get');
	if isequal(FileName,0),		return,		end
	D = grdinfo_m([PathName FileName],'silent');

	% Verify if, on the overlapping zone, the nodes of the larger grid cuincide with nodes of the smaler
	xoff_w = abs(head(1) - D(1));   %xoff_e = abs(head(2) - D(2));
	yoff_s = abs(head(3) - D(3));   %yoff_n = abs(head(4) - D(4));
	dx_w = xoff_w / D(8);           %dx_e = xoff_e / D(8);
	dy_s = yoff_s / D(9);           %dy_n = yoff_n / D(9);
	if ((dx_w - fix(dx_w)) > small)     % Need to adjust at the west border
		adjust_w = 1;
		if ((dx_w - fix(dx_w)) > head(8)/2) % The adjustment is to the east
			in_w_to_e = 1;
		end
	end
	if ((dy_s - fix(dy_s)) > small)     % Need to adjust at the south border
		adjust_s = 1;
		if ((dy_s - fix(dy_s)) > head(9)/2) % The adjustment is to the north
			in_s_to_n = 1;
		end
	end
	% We also have to test if the grid increments of both grids are multiples
	inc_x_ratio = head(8) / D(8);    inc_y_ratio = head(9) / D(9);
	if ((inc_x_ratio - fix(inc_x_ratio)) > 1e-3)
		adjust_x_inc = 1;
		new_x_inc = (D(2) - D(1)) / (fix(inc_x_ratio) - 1);
	end
	if ((inc_y_ratio - fix(inc_y_ratio)) > 1e-3)
		adjust_y_inc = 1;
		%new_y_inc = (D(4) - D(3)) / (fix(inc_y_ratio) - 1);
	end

	if (adjust_w || adjust_s || adjust_x_inc || adjust_y_inc)
		warndlg('The fine grid doesn''t fit well within the larger grid. Trying to fix it.','Warning')
        if (adjust_w && adjust_s) % We have to adjust all borders
			if (in_w_to_e)		% Move x origin to the east
				new_w = fix(dx_w) * (head(8) + 1);
				new_e = new_w + D(2) - D(1);
			else				% Move x origin to the west
				new_w = fix(dx_w) * head(8);
				new_e = new_w + D(2) - D(1);
			end
			if (in_s_to_n)		% Move y origin to the north
				new_s = fix(dy_s) * (head(9) + 1);
				new_n = new_s + D(4) - D(3);
			else				% Move y origin to the south
				new_s = fix(dy_s) * head(9);
				new_n = new_s + D(4) - D(3);
			end
			str_R = sprintf(' -R%.10f/%.10f/%.10f/%.10f', new_w, new_s, new_e, new_n);
        elseif (adjust_w)        % Need to adjust only est and west borders
			if (in_w_to_e)      % Move x origin to the east
				new_w = fix(dx_w) * (head(8) + 1);
				new_e = new_w + D(2) - D(1);
			else                % Move x origin to the west
				new_w = fix(dx_w) * head(8);
				new_e = new_w + D(2) - D(1);
			end
			str_R = sprintf(' -R%.10f/%.10f/%.10f/%.10f',new_w, D(3), new_e, D(4));
        elseif (adjust_s)        % Need to adjust only south and north borders
			if (in_s_to_n)      % Move y origin to the north
				new_s = fix(dy_s) * (head(9) + 1);
				new_n = new_s + D(4) - D(3);
			else                % Move y origin to the south
				new_s = fix(dy_s) * head(9);
				new_n = new_s + D(4) - D(3);
			end
			str_R = sprintf(' -R%.10f/%.10f/%.10f/%.10f',D(1), new_s, D(2), new_n);
        else
			errordlg('Asneira desconhecida (Unknown error).','Error')
        end
	end

	if (adjust_x_inc)       % Nao vou testar/usar o y_inc. A malha tem (tera?) de ser quadrada
		str_I = sprintf(' -I%.10f',new_x_inc);
	end

	if (~isempty(str_R))        % Finer grid needs adjustment. Do it and return.
		str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
		[fName,pName] = put_or_get_file(handles,str1,'Select GMT grid','get');
		if (isempty(fName)),     return;     end
		fName = [fName '=6'];       % I want them in Surfer format
		str = ['grdsample ' [PathName FileName] str_R ' -G' [pName fName]];
		if (~isempty(str_I)),   str = [str str_I];    end
		if isunix,		s = unix(str);
		elseif ispc,	s = dos(str);
		else			errordlg('Unknown platform.','Error');
		end
		if ~(isequal(s,0))                  % An error as occured
			errordlg('Error running grdsample. Finer grid was not adjusted.','Error')
		end
		return
	end

	% The following Tag is very important to destinguish from MB tracks, which have Tags = MBtrack#
	% West border
	ny = fix((D(4) - D(3)) / head(9)) + 1;
	x1 = repmat(D(1),1,ny);
	y1 = linspace(D(3),D(3)+(ny-1)*head(9),ny);     % y_min -> y_max
	hold on;
	lineHand_w = plot(x1,y1,'-o','Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,...
		'MarkerEdgeColor','w','MarkerFaceColor','k', 'MarkerSize',4, 'Tag','Maregraph');
	draw_funs(lineHand_w,'line_uicontext')        % Set lines's uicontextmenu
	% South border
	nx = fix((D(2) - D(1)) / head(8)) + 1;
	x2 = linspace(D(1),D(1)+(nx-1)*head(8),nx);     % x_min -> x_max
	y2 = repmat(D(3),1,nx);
	lineHand_s = plot(x2,y2,'-o','Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,...
		'MarkerEdgeColor','w','MarkerFaceColor','k', 'MarkerSize',4, 'Tag','Maregraph');
	draw_funs(lineHand_s,'line_uicontext')        % Set lines's uicontextmenu
	% East border
	x3 = repmat(D(2),1,ny);
	lineHand_e = plot(x3,y1,'-o','Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,...
		'MarkerEdgeColor','w','MarkerFaceColor','k', 'MarkerSize',4, 'Tag','Maregraph');
	draw_funs(lineHand_e,'line_uicontext')        % Set lines's uicontextmenu
	% North border
	y3 = repmat(D(4),1,nx);
	lineHand_n = plot(x2,y3,'-o','Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,...
		'MarkerEdgeColor','w','MarkerFaceColor','k', 'MarkerSize',4, 'Tag','Maregraph');
	draw_funs(lineHand_n,'line_uicontext')        % Set lines's uicontextmenu
	hold off;

	% Now we have to compute the index of the maregraphs positions on the finer grid borders
	ind.x = inc_x_ratio * (0:length(x2)-1);
	ind.y = inc_y_ratio * (0:length(y1)-1);
	% Save the index in the corresponding line handles. This way, the user may delete the
	% line maregraphs that he is not interested in. The remaining maregraphs will be fished
	% out in swan_options, where a tsun2.dat file will be created
	set(lineHand_w,'UserData',ind);     setappdata(lineHand_w,'Side','W')
	set(lineHand_s,'UserData',ind);     setappdata(lineHand_s,'Side','S')
	set(lineHand_e,'UserData',ind);     setappdata(lineHand_e,'Side','E')
	set(lineHand_n,'UserData',ind);     setappdata(lineHand_n,'Side','N')

% --------------------------------------------------------------------
function Tsun2Compute(handles, opt)
	if (aux_funs('msg_dlg',14,handles)),	return,		end
	if (nargin == 1), opt = [];     end

	if (strcmp(opt,'write_params'))     % Just do what it says and return
		% See if we have lines of maregraphs corresponding to the finer grid edges
		h_mareg = findobj(handles.axes1,'Type','line','Tag','Maregraph');
		side = [];
		if (~isempty(h_mareg))
			side = zeros(1, numel(h_mareg));	ind = side;
			for (i = 1:numel(h_mareg))
				side(i) = getappdata(h_mareg(i),'Side');
				ind(i) = get(h_mareg(i),'UserData');
			end
		else
			errordlg('You don''t have any maregraphs. Tsun2 needs to be feed with water height trhough maregraphs.','Error')
			return
		end

		if (length(side) > 2)
			errordlg('Tsunamis cannot arrive at more than two edges of the finer grid. (If you are not convinced, think a bit more)','Error')
			return
		elseif (length(side) == 2 && strcmp(side(1),'W') && strcmp(side(1),'E'))
			errordlg('This is completly idiot. The wave arrives on the West & East borders?','Chico Clever')
			return    
		elseif (length(side) == 2 && strcmp(side(1),'S') && strcmp(side(1),'N'))
			errordlg('This is completly idiot. The wave arrives on the North & South borders?','Chico Clever')
			return    
		elseif (isempty(side))
			errordlg('These maregraphs are not of the correct type to use with the tsun2 code.','Error')
			return
		end

		str1 = {'*.par', 'params file (*.par)';'*.*', 'All Files (*.*)'};
		[FileName,PathName] = put_or_get_file(handles,str1,'Select tsun2 params file','put');
		if isequal(FileName,0),		return,		end
		fid = fopen([PathName FileName],'wt');
		fprintf(fid,'%s\n%s\n','# Mirone generated tsun2 parameter file','#');

		% OK, we have only one or two edges. Lets find out whitch are they
		if (length(side) == 1 && strcmp(opt,'write_params'))     % Single edge
			switch side
				case 'W',       fprintf(fid,'%s\n','W');     fprintf(fid,'%d ',ind.y);
				case 'S',       fprintf(fid,'%s\n','S');     fprintf(fid,'%d ',ind.x);
				case 'E',       fprintf(fid,'%s\n','E');     fprintf(fid,'%d ',ind.y);
				case 'N',       fprintf(fid,'%s\n','N');     fprintf(fid,'%d ',ind.x);
			end
		elseif (length(side) == 2 && strcmp(opt,'write_params'))
			if (strcmp(char(side(1)),'S') && strcmp(char(side(2)),'W'))
				fprintf(fid,'%s %s\n','W','S');
				fprintf(fid,'%d ',ind(2).y);   fprintf(fid,'\n');      fprintf(fid,'%d ',ind(1).x);
			elseif (strcmp(char(side(1)),'E') && strcmp(char(side(2)),'S'))
				fprintf(fid,'%s %s\n','S','E');
				fprintf(fid,'%d ',ind(2).x);   fprintf(fid,'\n');      fprintf(fid,'%d ',ind(1).y);
			elseif (strcmp(char(side(1)),'N') && strcmp(char(side(2)),'E'))
				fprintf(fid,'%s %s\n','E','N');
				fprintf(fid,'%d ',ind(2).y);   fprintf(fid,'\n');      fprintf(fid,'%d ',ind(1).x);
			elseif (strcmp(char(side(1)),'N') && strcmp(char(side(2)),'W'))
				fprintf(fid,'%s %s\n','W','N');
				fprintf(fid,'%d ',ind(2).y);   fprintf(fid,'\n');      fprintf(fid,'%d ',ind(1).x);
			end
		end
		fclose(fid);
		return
	end         % END of write_params

	out = swan_options(handles,'Tsun2');
	if isempty(out),	return,		end
	pause(0.05);        % Give time to swan_options window to die

	opt_P = ' ';		extra_args2 = ' ';		opt_N = ' ';	opt_G = ' ';
	opt_I = ' ';		opt_J = ' ';			opt_F = ' ';
	if (isfield(out,'maregraph_xy') && isfield(out,'params_file_name'))
        opt_P = ['-P' out.params_file_name];
	end
	if (isfield(out,'opt_M')),		extra_args2 = '-M';
	elseif (isfield(out,'opt_D')),	extra_args2 = '-D';
	end
	if (isfield(out,'opt_N')),  opt_N = out.opt_N;   end
	if (isfield(out,'opt_G')),  opt_G = out.opt_G;   end
	if (isfield(out,'opt_I')),  opt_I = out.opt_I;   end
	if (isfield(out,'opt_F')),  opt_F = out.opt_F;   end
	if (isfield(out,'opt_O'))			% New (10/2007) output maregraphs option
		opt_O_xy = out.opt_O.xy;		% The maregraphs locations
		fname = out.opt_O.name;			% The maregraphs file name which will be used to compose a new name
		[pato,name] = fileparts(fname);
		fname = [pato filesep name '_maregHeights.dat'];	% The new name
		opt_O_name = ['-O' fname];
	else
		opt_O_xy = ' ';
		opt_O_name = ' ';
	end

	[X,Y,Z_bat,head_bat] = load_grd(handles);   Z_bat = double(Z_bat);
	if (isempty(Z_bat) || isempty(head_bat))
		errordlg('ERROR: one or more of the bat variables are empty where they souldn''t be.','Error');
		return
	end
	if ( abs(head_bat(8) - head_bat(9)) > 1e-3 )
		warndlg('Grid cells are not square. I don''t know the effect of this.','Warning')
	end
	if (~handles.IamCompiled),	tsun2_hand = @tsun2;          % To stop once for all with the bloody version mixing
	else						tsun2_hand = @tsun2_sem_wbar;
	end

	if (isfield(out,'maregraph_xy'))    % Ask for computation of maregraphs (WRONG - This is not an option, but will be in future)
		dt1 = diff(out.maregraph_xy(:,1));    dt2 = diff(dt1);      t0 = out.maregraph_xy(1,1);
		if any(abs(dt2) > 100*eps)      % First column doesn't have the time (eps due to very small dts rounding errors)
			warndlg('The maregraph file does not have a time increment. I will assume it is 1 second.','SEVERE WARNING')
			opt_I = sprintf('-I%f/1',head_bat(8));
		else                    % First column of maregs file has the time. We don't want it
			dt = dt1(1);        % This is the time increment to be used as option to tsun2
			out.maregraph_xy(:,1) = [];
			cfl = head_bat(8) / sqrt(abs(head_bat(5))*9.8);
			if (cfl <= dt)
				msg{1} = 'Your wave propagates faster than the maregraph time increment. Expect divergent results.';
				msg{2} = '';
				msg{3} = ['To solve this your time increment in the maregraph file must be at least < ' num2str(cfl*.9)];
				warndlg(msg,'SEVERE WARNING')
			end
			opt_I = sprintf('-I%f/%f',head_bat(8),dt);
			if (t0 > 1),    opt_J = ['-J' num2str(t0)];   end     % If we start computing at a latter time
		end
		if (isfield(out,'opt_m'))       % Movie option
			tmovie = feval(tsun2_hand,Z_bat, head_bat, out.maregraph_xy, opt_O_xy, opt_P, ...
				extra_args2, opt_N, opt_G, opt_I, opt_J, opt_F, opt_O_name, '-f');
		else
			feval(tsun2_hand,Z_bat, head_bat, out.maregraph_xy, opt_O_xy, opt_P, ...
				extra_args2, opt_N, opt_G, opt_I, opt_J, opt_F, opt_O_name);
		end    
	else                                % Compute grids or movie. Even if asked ignores opt_O 
        if (isfield(out,'opt_m'))       % Movie option
			tmovie = feval(tsun2_hand,Z_bat, head_bat, extra_args2,opt_N,opt_G,opt_I,opt_F,'-f');
        else
			feval(tsun2_hand,Z_bat, head_bat, extra_args2,opt_N,opt_G,opt_I,opt_F);
        end
	end

	if (isfield(out,'opt_m')),   do_movie(handles,tmovie,'tsun2');   end

%------------------------------------------------------------------------------------------
function do_movie(handles,tmovie,opt)
	if (strcmp(opt,'swan')),	is_swan = 1;
	else						is_swan = 0;
	end
	[m,n,k] = size(tmovie);
	opt_E = '-E60/30/0.55/0.6/0.4/10';      % Should be "controlable"
	n_crop = 0;    % Should be a variabe. It's used to crop borders and therefore hide grid borders reflections
	[X,Y,Z0,head] = load_grd(handles);
	if isempty(Z0),   return;     end;    % An error message was already issued
	if (n_crop),    Z0 = double(Z0(n_crop+1:m-n_crop, n_crop+1:n-n_crop));
	else            Z0 = double(Z0);
	end
	h = aguentabar('title','Wait again (computing movie)');

	cmap = flipud(hot(256));
	cmap(1,:) = [0 0 0.8];      cmap(2,:) = [0 0 0.8];
	idx0 = (Z0 <= 0);           Z0(idx0) = 0;
	if (is_swan)
		head(5) = 0;
		R0 = grdgradient_m(Z0,head,opt_E);
		R0 = flipud(R0);    idx0 = flipud(idx0);
	else
		clear idx0;
	end
	Z0 = flipud(Z0);

	i = 1;
	while (i <= k)
		Z = tmovie(:,:,1);
		tmovie(:,:,1) = [];      % Free this page (we don't need it anymore)
		if (n_crop),    Z = double(Z(n_crop+1:m-n_crop, n_crop+1:n-n_crop));
		else            Z = double(Z);
		end
		z_max = max(max(Z));  z_min = min(min(Z));
		head(5) = z_min;      head(6) = z_max;
		R = grdgradient_m(Z,head,opt_E);
		if (is_swan)
            Z(~idx0) = Z0(~idx0);        R(~idx0) = R0(~idx0);
            z_max = head(6);
		end

		dif_z = abs(Z - Z0);    idx = (dif_z > 0.05);
		if (n_crop),    img = repmat(uint8(0),m-2*n_crop,n-2*n_crop);
		else            img = repmat(uint8(0),m,n);
		end
		img(~idx) = uint8(round( (Z(~idx) / z_max)*254 + 1 ));
		img(idx) = 1;
		idx = (Z == 0);    img(idx) = 1;

		img = ind2rgb8(img,cmap);   img = shading_mat(img,R);
		M(i) = im2frame(img);
		aguentabar(i/k);
		i = i + 1;
	end
	if (ishandle(h))	close(h),	end		% In case it was not close inside the loop

	str1 = {'*.avi;*.AVI', 'avi files (*.avi,*.AVI)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select movie name','put');
	if isequal(FileName,0);     return;     end
	set(handles.figure1,'pointer','watch')
	movie2avi_j(M,[PathName FileName],'compression','none','fps',5)
	set(handles.figure1,'pointer','arrow');
	
% -----------------------------------------------------------------------------------------
function res = check_IsRectangle(h)
% Check if h is a handle to a rectangle.
	x = get(h,'XData');   y = get(h,'YData');
	if ~( (x(1) == x(end)) && (y(1) == y(end)) && length(x) == 5 && ...
			(x(1) == x(2)) && (x(3) == x(4)) && (y(1) == y(4)) && (y(2) == y(3)) )
		res = false;
	else
		res = true;
	end
	

function varargout = mesher_helper(varargin)
% Helper window to create an unstructured grid. It creates a GMT script that does the triangulations step

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

	if (isempty(varargin)),		return,		end

	hObject = figure('Vis','off');
	mesher_helper_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right');

	handles.hMirAx     = varargin{1};
	handles.hMirFig    = get(handles.hMirAx, 'Parent');
	[handles.tree, handles.hPolys] = tree_fill(handles.hMirAx);		% hPolys are sorted in the 'tree' order as well 
	handles.tree_index = tree_index(handles);
	n_polys = numel(handles.tree_index);
	handles.tri_sizes  = cell(n_polys, 1);		% Will be empty if nesting polygs do not exist yet
	handles.data_fname = cell(n_polys, 1);		% To hold the data file name associated by each polygon
	handles.interp_first = false(n_polys,1);	% Set to true  to force grid interp step in script
	handles.is_grid    = false(n_polys,1);		% Set to true  to indicate that input data is a grid
	handles.is_binary  = false(n_polys,1);		% Set to true  to indicate that input data is a binary file
	handles.single     = true(n_polys,1);		% Set to false to indicate that input binary data is double
	handles.out_name   = '';					% To hold the output file name
	fill_nodes_popup(handles)
	handles = check_conf(handles);				% Check if polygs have previous settings and update uicontrols

	plugedWin = getappdata(handles.hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];		% Add this figure handle to the carra?as list
	setappdata(handles.hMirFig,'dependentFigs',plugedWin);

	set(hObject,'Visible','on');
	
	% Choose default command line output for atlas_export
	if (nargout),       varargout{1} = hObject;    end
	guidata(hObject, handles);

% --------------------------------------------------------------------------------------
function push_pickBB_CB(hObject, handles)
% Select the root polygon (rectangle) which is supposed to be just a rectangle and turn it into
% a root rectangle. For that, erase its previous and generalist line properties.

	set(handles.hMirFig,'pointer','crosshair')
	hLine = get_polygon(handles.hMirFig);		% Get the line handle
	set(handles.hMirFig,'pointer','arrow');
	if (~isempty(hLine))
		handMir = guidata(hLine);
		cmenuHand = uicontextmenu('Parent', handMir.figure1);
		set(hLine, 'UIContextMenu', cmenuHand)
		uimenu(cmenuHand, 'Label', 'Rectangle limits (edit)', 'Call', 'draw_funs([],''rectangle_limits'')');
		uimenu(cmenuHand, 'Label', 'New nested polygon', 'Call', {@new_nestPolyg, hLine});
		draw_funs([],'set_common_lineProps', hLine, cmenuHand, false)
		set(hLine,'Tag','polymesh')
		setappdata(hLine,'family',[0 1])		% In this case the parent handle does not exist, so the 0
		conf = struct('inc','', 'interp',false, 'fname', '', 'is_grid', false, 'is_binary',false, 'single',true);
		setappdata(hLine, 'config', conf)		% Store (blank) configuration struct
		[handles.tree, handles.hPolys] = tree_fill(handles.hMirAx);
		guidata(handles.figure1, handles)
	end

% --------------------------------------------------------------------------------------
function push_drawBB_CB(hObject, handles)
% Draw the BoundingBox rectangle that will be the root polygon

	[p1,p2,hl] = rubberbandbox(handles.hMirAx);
	handMir = guidata(hl);

	difa = abs(p2 - p1);
	if ((difa(1) < handMir.head(7)/4) || (difa(2) < handMir.head(8)/4))
		delete(hl),		return			% Don't draw ultra small rectangles
	end
	set(hl,'Color',handMir.DefLineColor,'LineWidth',handMir.DefLineThick)	% Use defaults LineThick and DefLineColor

	cmenuHand = uicontextmenu('Parent',handMir.figure1);
	set(hl, 'UIContextMenu', cmenuHand)
	uimenu(cmenuHand, 'Label', 'Rectangle limits (edit)', 'Call', 'draw_funs([],''rectangle_limits'')');
	uimenu(cmenuHand, 'Label', 'New nested polygon', 'Call', {@new_nestPolyg, hl});
	draw_funs([],'set_common_lineProps', hl, cmenuHand, false)
	set(hl,'Tag','polymesh')
	setappdata(hl,'family',[0 1])		% In this case the parent handle does not exist, so the 0
	ui_edit_polygon(hl)
	[handles.tree, handles.hPolys] = tree_fill(handles.hMirAx);
	guidata(handles.figure1, handles)

% --------------------------------------------------------------------------------------
function popup_nodes_CB(hObject, handles)
% Read the selected row, decode it, and update the corresponding triangle size edit box

	str = get(hObject, 'Str');		val = get(hObject, 'Val');
	[lev, grp, pol] = get_address(str{val});
	set(handles.edit_triangSize, 'Str', handles.tri_sizes{val})
	set(handles.edit_dataFile,   'Str', handles.data_fname{val})
	set(handles.check_force,     'Val', handles.interp_first(val))
	hL = handles.tree{lev}{grp}(pol,1);
	h  = copyobj(hL, handles.hMirAx);
	rmappdata(h,'polygon_data')			% Remove the parent's ui_edit_polygon appdata
	ui_edit_polygon(h)					% And set a new one
	delete(findobj(handles.hMirAx, '-depth',1,'Tag','nestEnhancedLine'))	% Delete old one (if exists)
	prop = 'EdgeColor';					% Default color keyword for the most common case (patches)
	if (strcmpi(get(hL, 'type'), 'line')),	prop = 'Color';		end
	lt = get(hL, 'LineWidth');		lc = get(hL, prop);
	set(h,'LineWidth',lt+2, prop,1-lc, 'Tag','nestEnhancedLine')
	uistack_j(h,'bottom')	

% --------------------------------------------------------------------------------------
function edit_triangSize_CB(hObject, handles)
% Get and save the triangle size for this polygon

	str = get(hObject, 'Str');
	pos = get(handles.popup_nodes, 'Val');		s = get(handles.popup_nodes, 'Str');
	if (s{pos}(end) == '*')				% Remove the asterisk to indicate this polygon is processed
		s{pos} = s{pos}(1:end-2);
		set(handles.popup_nodes, 'Str', s)
	end
	handles.tri_sizes{pos} = str;		% For now, we just save what we got (no test/parsing)
	guidata(handles.figure1, handles)
	update_conf(handles, 'inc', str, s{pos})	% Update polygon's config appdata structure

% --------------------------------------------------------------------------------------
function check_force_CB(hObject, handles)
% Save info that this polygon wants GMT to first interpolate into a grid and
% latter convert to xyz and use that data to feed the triangulate step.
	state = get(hObject, 'Val');
	pos = get(handles.popup_nodes, 'Val');
	handles.interp_first(pos) = state;
	update_conf(handles, 'interp', state)	% Update polygon's config appdata structure
	guidata(handles.figure1, handles)

% --------------------------------------------------------------------------------------
function edit_dataFile_CB(hObject, handles)
% Get the name of the data file associated with the current polygon

	fname = get(hObject,'String');
	if (isempty(fname)),	return,		end		% Must also reset the corresponding node to non-processed state	
	if (strcmp(fname, 'NaN'))		% Particular case to indicate that the polygon this applies to is a data hole.
		pos = get(handles.popup_nodes, 'Val');
		handles.data_fname{pos} = fname;
		handles.is_grid(pos) = false;			set(handles.check_isGrid,  'Val', 0)
		handles.is_binary(pos) = false;			set(handles.check_isBinary,'Val', 0)
		handles.check_force(pos) = false;		set(handles.check_force,   'Val', 0)
		update_conf(handles, 'fname', fname)
		update_conf(handles, 'is_grid', false)
		update_conf(handles, 'is_binary', false)
		update_conf(handles, 'interp', false)
		update_conf(handles, 'single', true)
		guidata(handles.figure1, handles)
		return
	end
	push_dataFile_CB(hObject,guidata(gcbo),fname)	% Let this do all the work

% --------------------------------------------------------------------------------------
function push_dataFile_CB(hObject, handles, opt)
% Browse for a data filename to be associated with the currently selected polygon

	if (nargin == 3),	fname = opt;
	else				opt = [];
	end

	if (isempty(opt))           % Otherwise we already know fname from the 4th input argument
		str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
		handMir = guidata(handles.hMirFig);
		[FileName,PathName] = put_or_get_file(handMir,str1,'Select data file','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName,FileName];
	end

	pos = get(handles.popup_nodes, 'Val');
	handles.data_fname{pos} = fname;			% => Needs lots of testing on what this is
	set(handles.edit_dataFile, 'Str', fname)
	update_conf(handles, 'fname', fname)		% Update polygon's config appdata structure

	drv = aux_funs('findFileType',fname);
	if ( any(strcmp(drv, {'gmt' 'geotif' 'ecw' 'envherd' 'dono'})) )
		set(handles.check_isGrid, 'Val',1, 'Vis','on')
		handles.is_grid(pos) = true;
		update_conf(handles, 'is_grid', 1)
	elseif (strcmp(drv, 'dat'))
		% But we also need to know if it's binary or not
		[bin, n_column] = guess_file(fname);
		if (isempty(bin))
			set(handles.edit_dataFile, 'Str', '')
			errordlg(['Error reading file (probably empty)' fname],'Error'),	return
		end
		if (isa(bin,'struct') || bin ~= 0)		% A binary file. For now I wont test any further
			set(handles.check_isBinary, 'Val',1, 'Vis','on')
			set([handle.radio_single handles.radio_double], 'Vis','on')			
			handles.is_binary(pos) = true;
			update_conf(handles, 'is_binary', true)
		elseif (n_column <= 2)
			set(handles.edit_dataFile, 'Str', '')
			errordlg(['This file ' fname ' doesn''t have at least 3 columns (xyz). Bye Bye'],'Error'),	return
		end
		set(handles.check_isGrid, 'Val',0)
		update_conf(handles, 'is_grid', false)
	end
	guidata(handles.figure1, handles)

% --------------------------------------------------------------------------------------
function check_isGrid_CB(hObject, handles)
% ...
	pos = get(handles.popup_nodes, 'Val');
	handles.is_grid(pos) = get(hObject, 'Val');
	update_conf(handles, 'is_grid', get(hObject, 'Val'))

% --------------------------------------------------------------------------------------
function check_isBinary_CB(hObject, handles)
% ...
	pos = get(handles.popup_nodes, 'Val');
	handles.is_binary(pos) = get(hObject, 'Val');
	update_conf(handles, 'is_binary', get(hObject, 'Val'))

% --------------------------------------------------------------------------------------
function radio_single_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val', 1),		return,		end
	set(handles.radio_double, 'Val', 0)
	update_conf(handles, 'single', true)

% --------------------------------------------------------------------------------------
function radio_double_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val', 1),		return,		end
	set(handles.radio_single, 'Val', 0)
	update_conf(handles, 'single', false)

% --------------------------------------------------------------------------------------
function push_doScript_CB(hObject, handles)
% Check that all polygons have the triangle size set and generate the GMT script

	if (isempty(handles.out_name))
		errordlg('Need to provide the output file name', 'Error'),	return
	elseif (any(strcmp(get(handles.popup_nodes, 'Str'), '*')))
		errordlg('One or more polygon are not configured (the ones with a * in Nodes). Please revise and come back', 'Error'),	return
	end

	polys = findobj(handles.hMirAx,'-depth',1,'Tag','polymesh');
	if (numel(polys) ~= numel(handles.hPolys))
		warndlg('You add more polygons since you opened this window, which is now outdated. Need to restart.', 'WarnError')
		delete(handles.figure1)
		mesher_helper(get(handles.hPolys(1), 'Parent'));	% The parent of hPolys is the Mirone axes
		return
	end

	tree_depth = numel(handles.tree);
	sc = cell(2*tree_depth - 1 + 8, 1);		l = 8;
	if (get(handles.radio_batch,'Val'))
		comm = 'REM ';	fdest = '%all_pts%';	sc{1} = '@echo off';	sc{6} = 'SET all_pts=lixo.dat';
	else
		comm = '# ';	fdest = '$all_pts';		sc{1} = '#!/usr/bin';	sc{6} = 'all_pts=lixo.dat';
	end
	sc{3} = [comm 'Coffeeright Mirone Tec'];
	sc{4} = [comm 'Automatic script to generate a Delaunay triangulation of an unstructured grid'];
	bo = ' -bo3f';

	names_polys = cell(numel(handles.hPolys), 1);	nNam = 1;	% To save the names of polygons and later save themselves
	for (lev = 1:tree_depth - 1)			% Loop over all levels - 1 (last level is dealt separately)
		r = handles.tree{lev};
		for (m = 1:numel(r))				% Loop over all groups at this level
			for (n = 1:size(r{m},1))		% Loop over all polygons of this group
% 				P = poly_diff(handles.tree, lev, grp(n,1));
				children = handles.tree{lev+1}{m};	% Mx2 matrix with the M children handles in first column

				n_rows = size(children,1);
				[fiche, polyFname, opt_I, opt_R, is_grid, is_binary, do_interp] = get_GMTopts(handles, [lev m n]);
				names_polys{nNam} = polyFname;	nNam = nNam + 1;	% Save the polygs names according to our logic
				bi = '';
				if (is_binary)				% Not used in all bellow cases
					bi = ' -bi3f';
					if (get(handles.radio_double, 'Val')),	bi(end) = 'd';	end		% WRONG. MUST BE ON A PER POLYG BASIS
				end

				if (n_rows == 0)			% CASE A)	(One polygon with no children)
					if (is_grid)
						sc{l} = ['grd2xyz ' fiche opt_R ' | gmtselect -F' polyFname ' | blockmedian' opt_I opt_R bo ' >> ' fdest];
					else
						if (do_interp)		% Must first interpolate into a grid
							sc{l} = ['surface ' fiche bi opt_R opt_I ' -Glixo.grd -V '];	l = l + 1;
							sc{l} = ['grd2xyz lixo.grd | gmtselect -F' polyFname bo ' >> ' fdest];
						else
							% Consider the possibility of an on purpose empty polygon (a data hole)
							if (~(strcmpi(polyFname, 'void') || strcmpi(polyFname, 'empty')))	% If not empty
								sc{l} = ['gmtselect ' fiche bi ' -F' polyFname ' | blockmedian' opt_I opt_R bo ' >> ' fdest];
							end
						end
					end
				else						% CASE B) (One polygon with one or more children)
					t = '';
					for (k = 1:n_rows)
						tmpName = sprintf('L-%d_G-%d_P-%d.dat', [lev+1 m k]);	% Children's polygons
						t = [t ' | gmtselect -F' tmpName ' -If'];
					end
					if (is_grid)
						sc{l} = ['grd2xyz ' fiche opt_R t ' | blockmedian' opt_I opt_R bo ' >> ' fdest];
					else
						if (do_interp)		% Must first interpolate into a grid
							sc{l} = ['surface ' fiche bi opt_R opt_I ' -Glixo.grd -V '];	l = l + 1;
							sc{l} = ['grd2xyz lixo.grd | gmtselect -F' polyFname t bo ' >> ' fdest];
						else
							sc{l} = ['gmtselect ' fiche bi ' -F' polyFname t ' | blockmedian' opt_I opt_R bo ' >> ' fdest];
						end
					end
				end
				l = l + 1;		sc{l} = '';		l = l + 1;	% Add empty lines for readibility

% 				for (i = 1:numel(P))
% 					sc{l} = [P(i).x P(i).y];	l = l + 1;
% 				end
			end
		end

		if (lev == tree_depth - 1)			% Special case of last level
			r = handles.tree{tree_depth};
			for (m = 1:numel(r))			% Loop over all groups of last level
				for (n = 1:size(r{m},1))	% Loop over all polygons of this group
					[fiche, polyFname, opt_I, opt_R, is_grid, is_binary, do_interp] = get_GMTopts(handles, [lev+1 m n]);
					names_polys{nNam} = polyFname;	nNam = nNam + 1;	% Save the polygs names according to our logic
					if (is_grid)
						sc{l} = ['grd2xyz ' fiche opt_R ' | gmtselect -F' polyFname ' | blockmedian' opt_I opt_R bo ' >> ' fdest];
					else
						bi = '';
						if (is_binary)
							bi = ' -bi3f';
							if (get(handles.radio_double, 'Val')),	bi(end) = 'd';	end		% WRONG
						end
						if (do_interp)		% Must first interpolate into a grid
							sc{l} = ['surface ' fiche bi opt_R opt_I ' -Glixo.grd -V '];	l = l + 1;
							sc{l} = ['grd2xyz lixo.grd | gmtselect -F' polyFname bo ' >> ' fdest];							
						else
							% Consider the possibility of an on purpose empty polygon (a data hole)
							if (~(strcmpi(polyFname, 'void') || strcmpi(polyFname, 'empty')))	% If not empty
								sc{l} = ['gmtselect ' fiche bi ' -F' polyFname ' | blockmedian' opt_I opt_R bo ' >> ' fdest];
							end
						end
					end
					l = l + 1;		sc{l} = '';		l = l + 1;		% Add empty lines for readibility
				end
			end
		end
	end
	sc{l} = ['triangulate -bi3f ' fdest ' > ' handles.out_name];

	% ---------------------------- WRITE the files section --------------------------------------------
	% ====== First the script itself ==========
	handMir = guidata(handles.hMirAx);		% To use in put_or_get_file()
	if (get(handles.radio_batch,'Val'))
		str1 = {'*.bat;*.BAT', 'Batch file (*.bat,*.BAT)'; '*.*', 'All Files (*.*)'};
		[FileName,PathName] = put_or_get_file(handMir,str1,'Select File name','put','.bat');
	else
		str1 = {'*.sh', 'Bash file (*.sh)'; '*.*', 'All Files (*.*)'};
		[FileName,PathName] = put_or_get_file(handMir,str1,'Select File name','put','.sh');
	end
	if isequal(FileName,0),		return,		end
	f_name = [PathName FileName];
	%double2ascii(f_name,sc,'%f','maybeMultis');
	double2ascii(f_name,sc,'%s');

	% ========== Now the polygons  =============
	for (k = 1:numel(names_polys))
		if (isempty(names_polys{k})),	continue,	end				% We might have booked in excess
		[lev, grp, pol] = get_address(names_polys{k}(1:end-4), 1);	% Remember that we need to strip the '.dat'
		x = get(handles.tree{lev}{grp}(pol,1), 'XData');
		y = get(handles.tree{lev}{grp}(pol,1), 'YData');
		f_name = [PathName names_polys{k}];
		double2ascii(f_name,[x(:) y(:)],'%f');
	end

% --------------------------------------------------------------------------------------
function radio_batch_CB(hObject, handles)
	if (~get(hObject,'Val')),		set(hObject,'Val', 1),		return,		end
	set(handles.radio_bash, 'Val', 0)

% --------------------------------------------------------------------------------------
function radio_bash_CB(hObject, handles)
	if (~get(hObject,'Val')),		set(hObject,'Val', 1),		return,		end
	set(handles.radio_batch, 'Val', 0)

% --------------------------------------------------------------------------------------
function edit_outFile_CB(hObject, handles)
% ...
	fname = get(hObject,'String');
	if (isempty(fname)),	return,		end		% Must also reset the corresponding node to non-processed state	
	push_outFile_CB(hObject,guidata(gcbo),fname)	% Let this do all the work

% --------------------------------------------------------------------------------------
function push_outFile_CB(hObject, handles, opt)
% ...
	if (nargin == 3),	fname = opt;
	else				opt = [];
	end

	if (isempty(opt))           % Otherwise we already know fname from the 4th input argument
		str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
		handMir = guidata(handles.hMirFig);
		[FileName,PathName] = put_or_get_file(handMir,str1,'Select data file','put','.dat');
		if isequal(FileName,0),		return,		end
		fname = [PathName,FileName];
	end
	set(handles.edit_outFile, 'Str', fname)
	hPol = handles.tree{1}{1}(1);
	setappdata(hPol, 'fname_out', fname)	% Save this here for eventual future Fig rebuild
	handles.out_name = fname;
	guidata(handles.figure1, handles)

% --------------------------------------------------------------------------------------
function [fname, polyg, opt_I, opt_R, is_grid, is_binary, do_interp] = get_GMTopts(handles, ind)
% Get the options needed to run the gmtselect and blockmedian steps in final script
% IND(1) Level
% IND(2) Group
% IND(3) Mx2 matrix with polygon handles in first column

	hPolyg = handles.tree{ind(1)}{ind(2)}(ind(3));		% Get the polygon referenced by IND
	x = get(hPolyg, 'XData');	y = get(hPolyg, 'yData');
	opt_R = sprintf(' -R%.10g/%.10g/%.10g/%.10g', min(x), max(x), min(y), max(y));
	for (k = 1:numel(handles.tree_index))				% Scan tree_index to find which element is equal to
		if (isequal(ind, handles.tree_index{k}))		% IND. That gives us the adress needed for tri sizes
			val = k;	% Found. But this is kind of ugly. I should be able to compute it directly from IND
			break
		end
	end
	opt_I = sprintf(' -I%s', handles.tri_sizes{val});
	fname = handles.data_fname{val};					% File name of the data associated with this polygon
	polyg = sprintf('L-%d_G-%d_P-%d.dat', ind);			% Name of this polygon
	is_grid   = handles.is_grid(val);
	is_binary = handles.is_binary(val);
	do_interp = handles.interp_first(val);

% --------------------------------------------------------------------------------------
function new_nestPolyg(hObj, evt, hPar)
% Draw a new polygon

	handMir = guidata(hPar);
	[xp,yp] = getline_j(handMir.figure1,'closed');
	if (numel(xp) < 4),		return,		end		% 4 because a straight line has 3 vertex (last one repeats)
	xp = xp(:)';	yp = yp(:)';
	h = patch('XData',xp,'YData',yp,'FaceColor','none','EdgeColor',handMir.DefLineColor,...
		'LineWidth',handMir.DefLineThick,'Tag','polymesh');

	cmenuHand = uicontextmenu('Parent',handMir.figure1);
	set(h, 'UIContextMenu', cmenuHand)
	uimenu(cmenuHand, 'Label', 'New nested polygon', 'Call', {@new_nestPolyg, h});
	uimenu(cmenuHand, 'Label', 'Delete', 'Call', {@del_nestPolyg, h});
	draw_funs([],'set_common_lineProps', h, cmenuHand, false)
	ui_edit_polygon(h)
	ap = getappdata(hPar, 'family');
	lev = ap(2) + 1;
	setappdata(h,'family',[hPar lev])		% Store parent handle and this level of nesting

	conf = struct('inc','', 'interp',false, 'fname', '', 'is_grid', false, 'is_binary',false, 'single',true);
	setappdata(h, 'config', conf)			% Store (blank) configuration struct

% --------------------------------------------------------------------------------------
function del_nestPolyg(hObj, evt, hPol)
% Delete the polygon with handle HPOL and ALL its descendents

	handMir = guidata(hPol);
	tree = tree_fill(handMir.axes1);
	family = getappdata(hPol, 'family');	% family = [hParent, lev];

	for (lev = family(2)+1:numel(tree))		% Loop over remaining levels
		r = tree{lev};
		for (grp = 1:(numel(r)))			% Loop number of groups
			p = r{grp}(1,2);				% parent handle of currently looped polygon
			if (hPol == p)					% If current polyg ascendent handle is equal to the to-be-killed
				p = r{grp}(:,1);			% polygon (argin) those are descendents and must die too.
				delete(p(ishandle(p)));
			end
		end
	end
	delete(hPol)

% 	plugedWin = getappdata(handMir.hMirFig,'dependentFigs');
% 	for (k = 1:numel(plugedWin))
% 		if (strcmp(get(plugedWin(k), 'Name')), 'Mesher helper')
% 			handles = guidata(plugedWin(k));
% 			break
% 		end
% 	end

	hFig = findobj(0, '-depth',1, 'type','figure', 'Name','Mesher helper');		% Find myself
	% If hFig exists it now needs to be updated because we just killed some guys
	if (~isempty(hFig))
		delete(hFig)
		mesher_helper(handMir.axes1);
	end
	
% --------------------------------------------------------------------------------------
function P = poly_diff(tree, lev, hPoly)
% Calculate the union of all LEV+1 polygons that are descendent of HPOLY (a handle)
% and than calculate the difference of both

	P1 = poly_union(tree, lev+1, 1);
	for (k = 2:numel(tree{lev+1}))			% Loop over the remaining number of groups
		P2 = poly_union(tree, lev+1, k);
		P1 = PolygonClip(P1, P2, 3);		% Union of the two polygons
	end
	
	P0.x = get(hPoly, 'XData');
	P0.y = get(hPoly, 'YData');		P0.hole = 0;
	P = PolygonClip(P0, P1, 0);				% Difference of the two polygons

% --------------------------------------------------------------------------------------
function P1 = poly_union(tree, lev, grp)
% Make a union of all polygons at a level LEV that belong to the group GRP

	r = tree{lev};		polys = r{grp};
	P1.x = get(polys(1,1), 'XData');	P2.hole = 0;
	P1.y = get(polys(1,1), 'YData');	P1.hole = 0;

	for (k = 2:size(polys,1))
		P2.x = get(polys(k,1), 'XData');
		P2.y = get(polys(k,1), 'YData');
		P1 = PolygonClip(P1, P2, 3);		% Union of the two polygons
	end

% --------------------------------------------------------------------------------------
function [tree, polys] = tree_fill(hAx)
% Find all "meshy" polygons and construct an hierarchic tree
% HAX is the Mirone axes1 handle that is the parent of all mesh polygons.
% TREE is a cell vector with N_LEVELS elements. Each element is in itself another cell
%      vector with N elements. One for each group (leaves of a parent node). Each of this
%      elements has a Mx2 matrix where first column holds the polygon handle and the second
%      the handle of its parent. So second column should be constant since by deffinition
%      all polygons in one group share the same ancestor polygon.

	polys = findobj(hAx,'-depth',1,'Tag','polymesh');
	if (isempty(polys)),	tree = [];		return,		end

	appd  = zeros(numel(polys), 2);
	for (k = 1:numel(polys))
		appd(k, :) = getappdata(polys(k),'family');		% family = [hParent, lev];
	end
	levs  = appd(:,2);					% All levels
	hPar  = appd(:,1);					% Handles of the Parent of each polygon
	[levs, ind] = sort(levs);
	polys = polys(ind);					% Now polygons are sorted in growing level order as well
	hPar  = hPar(ind);

	tree = cell(levs(end),1);
	for (k = levs(end):-1:2)
		indL = find(levs == k);			% Find all that have the same level
		[this_par,i] = sort(hPar(indL));% By sorting the parent handles we get them already grouped
		this_h = polys(indL);
		this_h = this_h(i);				% The handles of the polygons themselves, reordered to their parents order 
		g = find(diff(this_par) ~= 0);
		g = g(:)';						% Ensure its a row vector
		g = [1 g+1 numel(this_par)+1];	% Because the diff operator lowers 1 size and when one group only g was = []
		ng = numel(g) - 1;				% Number of groups. -1 because g has one more element to simplify the algo 
		t = cell(1, ng);
		for (m = 1:ng)					% Loop over number of groups at this level
			v = g(m):g(m+1)-1;
			t{m} = [this_h(v) this_par(v)];
		end
		tree{k} = t;
	end

	% Now the root polyg (rectangle)
	tree{1} = {[polys(1) NaN]};

% --------------------------------------------------------------------------------------
function ti = tree_index(handles)
% Travel along the cell array TREE and compute a cell vector with as many elements as the
% number of polygons. Each element contains a 1x3 row vector where its elements contain:
%	1-element	level indice (1:n_levels)
%	2-element	group indice (1:n_groups_for_this_level)
%	3-element	indice of the Nth element in this group
%				(1:n_rows of the handles matrix of polygons that make this group)
%
% We use this array to easily navigate into the TREE cell array (which contains the line handles)

	ti = cell(numel(handles.tree), 1);		% Probably not enough but it won't grow too much anyway
	j = 1;
	for (k = 1:numel(handles.tree))			% Loop over number of levels
		ng = numel(handles.tree{k});
		for (m = 1:ng)						% Loop over the number of groups at this level
			for (n = 1:size(handles.tree{k}{m}, 1))	% Loop over number of polygons of this group
				ti{j} = [k m n];
				j = j + 1;
			end
		end
	end

	if (nargout == 0)						% Save result in handles
		handles.tree_index = ti;
		guihandles(handles.figure1, handles);
	end

% --------------------------------------------------------------------------------------
function fill_nodes_popup(handles)
% Fill the Nodes popup with the list of polygons

	s = cell(numel(handles.tree_index), 1);
	for (k = 1:numel(handles.tree_index))
		s{k} = sprintf('L-%d G-%d P-%d *', handles.tree_index{k});
	end
	if (~isempty(s))	% Happens when main function is called from a Mir fig still with no nest polygs
		set(handles.popup_nodes, 'Str', s)
	end

% --------------------------------------------------------------------------------------
function handles = check_conf(handles)
% Check the contents of the 'conf' appdata of all nested poygons and use that info to fill
% the corresponding uicontrols. This is used to not have to enter the same information
% again and again when we close the figure, edit polygons and start this figure again.

	for (k = 1:numel(handles.hPolys))
		conf = getappdata(handles.hPolys(k), 'config');
		if (isempty(conf)),		continue,	end
		set(handles.edit_triangSize, 'Str', conf.inc),	handles.tri_sizes{k}    = conf.inc;
		set(handles.edit_dataFile, 'Str', conf.fname),	handles.data_fname{k}   = conf.fname;
		set(handles.check_force,   'Val', conf.interp),	handles.interp_first(k) = conf.interp;
		set(handles.radio_single,  'Val', conf.single),	handles.single(k)       = conf.single;
		set(handles.radio_double,  'Val', ~conf.single)
		str = get(handles.popup_nodes, 'Str');	s = str{k};
		s(end-1:end) = [];
		str{k} = s;
		set(handles.popup_nodes, 'Str', str)
	end
	if (~isempty(handles.hPolys))		% It's empty when we still do not have any mesh polygon
		handles.out_name = getappdata(handles.hPolys(1), 'fname_out');	% Since the polygons are sorted, this should work
		set(handles.edit_outFile, 'Str', handles.out_name)
		guidata(handles.figure1, handles)
	end

% --------------------------------------------------------------------------------------
function update_conf(handles, field, val, str)
% Update the 'conf' appdata configuration of currently selected polygon (via popup_nodes())
% There is no testing on arguments validity, so it's supposed that FIELD is the right name of
% a field of the struct 'conf' and VAL an appropriate value. For reference, here is the 'conf' deffinition
% conf = struct('inc',str, 'interp',false, 'fname', '', 'is_grid', false, 'is_binary',false, 'single',true);
% STR is an optional argument, holding the currently selected line in popup_nodes(). If not
% provided we'll fish it from within this function.

	if (nargin == 3)
		str = get(handles.popup_nodes, 'Str');	str = str{get(handles.popup_nodes, 'Val')};
	end

	[lev, grp, pol] = get_address(str);			% Get this polygon address in 'tree'
	hPol = handles.tree{lev}{grp}(pol);
	conf = getappdata(hPol, 'config');
	conf.(field) = val;
	setappdata(hPol, 'config', conf);

% --------------------------------------------------------------------------------------
function [lev, grp, pol] = get_address(str, opt)
% Decode a string of the form 'L-2 G-1 P-1' (from the nodes popup) as meaning:
% LEV = 2; GRP = 1; POL = 1
% That is, the coordinates of that polygon in the 'tree' cell array.
% OPT is used when STR == 'L-2_G-1_P-1'

	if (nargin == 2)
		str = strrep(str, '_', ' ');
	end
	[t, r] = strtok(str);	lev = str2double(t(3:end));
	[t, r] = strtok(r);		grp = str2double(t(3:end));
	t = strtok(r);			pol = str2double(t(3:end));

% --------------------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, evt)
% ...
	handles = guidata(hObject);
	try
		delete(findobj(handles.hMirAx, '-depth',1,'Tag','nestEnhancedLine')) % Delete this guy (if exists)
	end
	delete(hObject)

% --------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject)
% ...
	if isequal(get(hObject,'CurrentKey'),'escape')
		figure1_CloseRequestFcn(hObject, []);
	end


% --- Creates and returns a handle to the GUI figure. 
function mesher_helper_LayoutFcn(h1)

set(h1, 'Position',[520 620 390 241],...
	'Color',get(0,'factoryUicontrolBackgroundColor'),...
	'CloseRequestFcn',@figure1_CloseRequestFcn,...
	'KeyPressFcn',@figure1_KeyPressFcn,...
	'MenuBar','none',...
	'Name','Mesher helper',...
	'NumberTitle','off',...
	'Resize','off',...
	'HandleVisibility','callback',...
	'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 209 101 22],...
'Call',@mesh_helper_uiCB,...
'String','Draw BB rectangle',...
'TooltipString','Draw the BoundingBox rectangle',...
'Tag','push_drawBB');

uicontrol('Parent',h1, 'Position',[140 209 101 22],...
'Call',@mesh_helper_uiCB,...
'String','Pick BB rectangle',...
'TooltipString','Pick a generic ractangle and turn it into a BoundinBox root rectangle.',...
'Tag','push_pickBB');

uicontrol('Parent',h1, 'Position',[11 173 130 15],...
'HorizontalAlignment','left',...
'String','Approximate triang side',...
'Style','text');

uicontrol('Parent',h1, 'Position',[290 211 91 20],...
'BackgroundColor',[1 1 1],...
'Call',@mesh_helper_uiCB,...
'String',' ',...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_nodes');

uicontrol('Parent',h1, 'Position',[140 170 71 22],...
'BackgroundColor',[1 1 1],...
'Call',@mesh_helper_uiCB,...
'String','',...
'Style','edit',...
'TooltipString','Approximate size for the triangle''s side. Only respected if data allows the selected size.',...
'Tag','edit_triangSize');

uicontrol('Parent',h1, 'Position',[255 214 35 14],...
'HorizontalAlignment','left',...
'String','Nodes',...
'Style','text');

uicontrol('Parent',h1, 'Position',[229 170 135 23],...
'Call',@mesh_helper_uiCB,...
'String','Force by interpolation',...
'Style','checkbox',...
'TooltipString','Force trangle size by interpolating first into a regular grid. Needed when data is too sparse',...
'Tag','check_force');

uicontrol('Parent',h1, 'Position',[112 213 25 14],...
'FontSize',9,...
'FontWeight','bold',...
'String','OR',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 114 351 22],...
'BackgroundColor',[1 1 1],...
'Call',@mesh_helper_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'TooltipString','If empty, this polygon will use the data file of its parent. Enter ''Void'' or ''Empty'' (no quotes) if this polygon has no data inside it at all',...
'Tag','edit_dataFile');

uicontrol('Parent',h1, 'Position',[360 114 23 23],...
'Call',@mesh_helper_uiCB,...
'FontWeight','bold',...
'String','...',...
'Tag','push_dataFile');

uicontrol('Parent',h1, 'Position',[10 136 240 15],...
'HorizontalAlignment','left',...
'String','Data file name for this polygon (it may be a grid)',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 93 57 21],...
'Call',@mesh_helper_uiCB,...
'String','Is grid?',...
'Style','checkbox',...
'TooltipString','Input data for this polygon is a grid. Automatic guess, feel free to correct',...
'Tag','check_isGrid',...
'Visible','off');

uicontrol('Parent',h1, 'Position',[198 93 71 21],...
'Call',@mesh_helper_uiCB,...
'String','Is binary?',...
'Style','checkbox',...
'TooltipString','Input data for this polygon is binary. Automatic guess, feel free to correct',...
'Tag','check_isBinary',...
'Visible','off');

uicontrol('Parent',h1, 'Position',[271 93 55 21],...
'Call',@mesh_helper_uiCB,...
'String','Single',...
'Style','radiobutton',...
'TooltipString','Binary file is single precision',...
'Value',1,...
'Tag','radio_single',...
'Visible','off');

uicontrol('Parent',h1, 'Position',[327 93 55 21],...
'Call',@mesh_helper_uiCB,...
'String','Double',...
'Style','radiobutton',...
'TooltipString','Binary file is double precision',...
'Tag','radio_double',...
'Visible','off');

uicontrol('Parent',h1, 'Position',[10 44 351 22],...
'BackgroundColor',[1 1 1],...
'Call',@mesh_helper_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tag','edit_outFile');

uicontrol('Parent',h1, 'Position',[360 44 23 23],...
'Call',@mesh_helper_uiCB,...
'FontWeight','bold',...
'String','...',...
'Tag','push_outFile');

uicontrol('Parent',h1, 'Position',[10 10 231 22],...
'Call',@mesh_helper_uiCB,...
'FontSize',9,...
'String','Generate GMT script',...
'TooltipString','Generate a GMT script that will do the job.',...
'Tag','push_doScript');

uicontrol('Parent',h1, 'Position',[271 10 50 23],...
'Call',@mesh_helper_uiCB,...
'String','Batch',...
'Style','radiobutton',...
'TooltipString','Create a DOS batch file',...
'Value',1,...
'Tag','radio_batch');

uicontrol('Parent',h1, 'Position',[329 10 50 23],...
'Call',@mesh_helper_uiCB,...
'String','Bash',...
'Style','radiobutton',...
'TooltipString','Create a bash shell script',...
'Tag','radio_bash');

uicontrol('Parent',h1, 'Position',[10 66 91 15],...
'HorizontalAlignment','left',...
'String','Output file name',...
'Style','text');

function mesh_helper_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

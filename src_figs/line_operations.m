function varargout = line_operations(varargin)
% Wraper figure to perform vectorial operations on line/patch objects

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

	if (isempty(varargin)),		return,		end
	if (~isfield(varargin{1}, 'head')),		return,		end		% Call from an empty fig

	isPC = strncmp(computer,'PC',2);
	floating = isPC;
 
	hMirFig = varargin{1}.figure1;
	hMirAxes = varargin{1}.axes1;

	if (~floating)
		hLineOP = getappdata(hMirFig, 'hLineOP');		% Get the uicontrol handles (if they already exist)
		if (strcmp(get(varargin{1}.lineOP,'checked'), 'off') && isempty(hLineOP))		% First time use
			old_unit = get(hMirAxes,'units');	set(hMirAxes,'units','pixels')
			pos = get(hMirAxes,'pos');			set(hMirAxes,'units',old_unit)
			hObject = figure('Tag','figure1','MenuBar','none','Visible','off');
			h = line_operations_LayoutFcn(hObject,hMirFig, pos,-2);
			handles = guihandles(hObject);
			handles.edit_cmd		= h(1);
			handles.push_pickLine	= h(2);
			handles.popup_cmds		= h(3);
			handles.push_apply		= h(4);
			handles.push_semaforo	= h(5);
			handles.figure1			= hObject;
			set(varargin{1}.lineOP,'checked', 'on')
			setappdata(hMirFig, 'hLineOP', [hObject h])			% Save for an eventual future use

		elseif (strcmp(get(varargin{1}.lineOP,'checked'), 'off') && ~isempty(hLineOP))	% Reuse
			handles.figure1			= hLineOP(1);
			handles.edit_cmd		= hLineOP(2);
			handles.push_pickLine	= hLineOP(3);
			handles.popup_cmds		= hLineOP(4);
			handles.push_apply		= hLineOP(5);
			handles.push_semaforo	= hLineOP(6);
			hObject = hLineOP(1);
			set(handles.push_pickLine,'Tooltip','Have 0 lines to play with')	% Restart clean because we didn't store
			set(handles.push_semaforo,'BackgroundColor',[1 0 0])				% any eventual picked line handles
			set(hLineOP(2:end), 'Vis', 'on')
			set(varargin{1}.lineOP,'checked', 'on')

		else				% Hide, clean and leave
			set(hLineOP, 'Vis', 'off')
			set(varargin{1}.lineOP,'checked', 'off')
			return
		end
	else
		hObject = figure('Tag','figure1','MenuBar','none','Vis','off');
		line_operations_LayoutFcn(hObject);
		handles = guihandles(hObject);
		% Shrink the figure a bit
		pos = get(handles.popup_cmds, 'Pos');		set(handles.popup_cmds, 'Pos', pos - [10 0 40 0])
		pos = get(handles.push_apply, 'Pos');		set(handles.push_apply, 'Pos', pos - [70 0 10 0])
		pos = get(handles.figure1, 'Pos');			set(handles.figure1, 'Pos', pos - [0 0 80 0])
	end

	handles.hMirFig = hMirFig;
	handles.hMirAxes = hMirAxes;
	handles.lt = varargin{1}.DefLineThick;
	handles.lc = varargin{1}.DefLineColor;
	handles.geog = varargin{1}.geog;
	handles.head = varargin{1}.head;
	handles.path_tmp = varargin{1}.path_tmp;
	handles.isPC = isPC;
	IamCompiled = varargin{1}.IamCompiled;
	handles.version7 = varargin{1}.version7;

	handles.known_ops = {'bezier'; 'buffer'; 'bspline'; 'closing'; 'cspline'; 'delete'; 'group'; ...
			 'line2patch'; 'polysimplify'; 'polyunion'; 'polyintersect'; 'polyxor'; 'polyminus'; 'pline'; ...
			'scale'; 'stitch'; 'thicken'; 'toRidge'; 'self-crossings'};
	handles.hLine = [];
	popup_cmds = {'Choose command'; 'bezier N'; 'buffer DIST'; 'bspline'; 'closing DIST'; 'cspline N RES'; ...
			'delete DUP|SMALL N|SPUR'; 'group lines'; 'line2patch'; 'polysimplify TOL'; 'polyunion'; ...
			'polyintersect'; 'polyxor'; 'polyminus'; 'pline [x1 ..xn; y1 .. yn]'; 'scale to [-0.5 0.5]'; ...
			'stitch TOL'; 'thicken N'; 'toRidge 5'; 'self-crossings (find)'};

	if (~IamCompiled)
		handles.known_ops{end+1} = 'hand2Workspace';
		popup_cmds{end+1} = 'hand2Workspace';
	end

	% --- See if we have one GMT database file loaded, as they trigger a dedicated function -----
	hGMT_DB = findobj(handles.hMirAxes, 'tag','GMT_DBpolyline');
	if (~isempty(hGMT_DB))
		handles.known_ops{end+1} = 'GMT_DB';
		popup_cmds{end+1} = 'GMT_DB';
	end
	% -------------------------------------------------------------------------------------------

	set(handles.popup_cmds,'Tooltip', 'Select one operation from this list')
	set(handles.popup_cmds,'String', popup_cmds);

	if (floating)
		move2side(handles.hMirFig, hObject, 'north')
		set(hObject,'Visible','on');
		if (handles.version7 < 8.4),	WindowAPI(hObject, 'TopMost'),	end
	end

	handles.ttips = cell(numel(handles.known_ops)+1,1);
	handles.ttips{1} = 'Select one operation from this list';
	handles.ttips{2} = sprintf(['Fit a Bezier curve to a polyline.\n' ...
								'Replace N by the number of nodes of the Bezier curve [default 100].']);
	handles.ttips{3} = sprintf(['Compute buffer zones arround polylines.\n' ...
								'Replace DIST by the width of the buffer zone.\n' ...
								'Optionaly if map is in geogs apend M, K or N to DIST\n' ...
								'to indicate Meters, Kilometers or NMiles (e.g. 10M)\n\n' ...
								'Optionally you can append NPTS, DIR or GEOD options. Where:\n' ...
								'NPTS -> Number of points used to contruct the circles\n' ...
								'around each polygon vertex. If omitted, default is 13.\n\n' ...
								'DIR -> ''in'', ''out'' or ''both'' selects whether to plot\n' ...
								'inside, outside or both delineations of the buffer zone.\n\n' ...
								'GEOD -> ''geod'', or [A B|F] uses an ellipsoidal model.\n' ...
								'First of this cases uses WGS-84 Ellipsoid. Second form selects\n' ...
								'an ellipsoid with A = semi-major axis and B = semi-minor axis\n' ...
								'or B = ellipsoid Flatenning.']);
	handles.ttips{4} = sprintf(['Smooth line with a B-form spline\n' ...
								'It will open a help control window to\n' ...
								'help with selection of nice parameters']);
	handles.ttips{5} = sprintf(['Compute the equivalent of the image processing ''Closing'' operation\n' ...
								'tha consist in buffer out followed by buffer in by the same amount.\n' ...
								'Replace DIST by the width of the buffer zone.\n' ...
								'For further options to provide finer control, see help of the buffer operation.']);							
	handles.ttips{6} = sprintf(['Smooth line with a cardinal spline.\n' ...
								'Replace N by the downsampling rate. For example a N of 10\n' ...
								'will take one every other 10 vertices of the polyline.\n' ...
								'The optional RES argument represents the number subdivisions\n' ...
								'between consecutive data points. As an example RES = 10 will\n' ...
								'split the downsampled interval in 10 sub-intervald, thus reseting\n' ...
								'the original number of point, excetp at the end of the line.\n' ...
								'If omited RES defaults to 10.']);
	handles.ttips{7} = sprintf(['Delete duplicate lines OR delete lines with with less than N\n' ...
								'vertices, OR remove spurs from lines. For the first case use the form\n' ...
								'delete DUP\n' ...
								'For the second case (where N is the line''s minimum number of vertices)\n' ...
								'delete SMALL N (e.g. delete SMALL 10\n' ...
								'For the third case use the form (here you can also select the line to work on)\n' ...
								'delete SPUR']);
	handles.ttips{8} = sprintf(['Group lines that have exactly the same characteristics.\n' ...
								'Unfortunately Matlab has a very bad memory management and it\n' ...
								'very slow when there are many different lines plotted.\n' ...
								'This option groups same type lines into one single multi-segment\n' ...
								'line. Visually it will look the same but performance jumps.']);
	handles.ttips{9} = 'Convert line objects into patch. Patches, for example, accept fill color.';
	handles.ttips{10} = sprintf(['Approximates polygonal curve with desired precision\n' ...
								'using the Douglas-Peucker algorithm.\n' ...
								'Replace TOL by the desired approximation accuracy.\n' ...
								'When data is in geogs, TOL is the tolerance in km.']);
	handles.ttips{11}  = 'Performs the boolean operation of Union to the selected polygons.';
	handles.ttips{12} = 'Performs the boolean operation of Intersection to the selected polygons.';
	handles.ttips{13} = 'Performs the boolean operation of exclusive OR to the selected polygons.';
	handles.ttips{14} = 'Performs the boolean operation of subtraction to the selected polygons.';
	handles.ttips{15} = sprintf(['Dray a polyline with vertices defined by coords [x1 xn; y1 yn].\n' ...
								'Note: you must use the brackets and semi-comma notation as above.\n' ...
								'Example vector: [1 1.5 3.1; 2 4 8.4]']);
	handles.ttips{16} = 'Scale to the [-0.5 0.5] interval.';
	handles.ttips{17} = sprintf(['Stitch in cascade the lines that are closer than TOL to selected line.\n' ...
								'Replace TOL by the desired maximum distance for lines still be stitched.\n' ...
								'If removed or left as the string "TOL" (no quotes) it defaults to Inf.\n' ...
								'Optionaly if map is in geogs apend M, K or N to TOL to indicate\n' ...
								'Meters, Kilometers or NMiles (e.g. 10M)\n\n' ...
								'Apend ALL to the command to stich all lines without needing to select one.']);
	handles.ttips{18} = sprintf(['Thicken line object to a thickness corresponding to N grid cells.\n' ...
								'The interest of this comes when used trough the "Extract profile"\n' ...
								'option. Since the thickned line stored in its pocked N + 1 parallel\n' ...
								'lines, roughly separate by 1 grid cell size, the profile interpolation\n' ...
								'is carried on those N + 1 lines, which are averaged (stacked) in the end.']);
	handles.ttips{19} = sprintf(['Calculate a new line with vertex siting on top of nearby ridges.\n' ...
								'The parameter N is used to search for ridges only inside a sub-region\n' ...
								'2Nx2N centered on current vertex. Default is 5, but you can change it.']);
	handles.ttips{20}  = 'Detect self-crossings of all plotted lines/patches OR of the selected line.';

	n = 20;			% NEED TO BE EQUAL TO LAST EXPLICITLY ttips{n} above
	if (~IamCompiled)
		n = n + 1;
		handles.ttips{n} = sprintf(['Send the selected object handles to the Matlab workspace.\n' ...
								'Use this if you want to gain access to all handle properties.']);
	end
	if (~isempty(hGMT_DB))
		n = n + 1;
		handles.ttips{n} = sprintf(['Update the GMT coastlines, rivers or boundaries DB.\n' ...
								'Just load the to-be-updated lines and hit "Apply".']);
	end

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hMirFig,'dependentFigs',plugedWin);

	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end

% --------------------------------------------------------------------------------------------------
function push_pickLine_CB(hObject, handles)    
    set(handles.hMirFig,'pointer','crosshair')
    hLine = get_polygon(handles.hMirFig,'multi');        % Get the line handle
	if (numel(hLine) >= 1)
		handles.hLine = unique(hLine);
		set(handles.push_semaforo,'BackgroundColor',[0 1 0])
	else
		set(handles.push_semaforo,'BackgroundColor',[1 0 0])
		handles.hLine = [];
	end
	set(hObject,'Tooltip', ['Have ' sprintf('%d',numel(hLine)) ' lines to play with'])

    set(handles.hMirFig,'pointer','arrow')
	guidata(handles.figure1, handles);

% --------------------------------------------------------------------------------------------------
function popup_cmds_CB(hObject, handles)
	str = get(hObject,'String');	val = get(hObject,'Val');
	set(hObject, 'Tooltip', handles.ttips{val})
	if (val == 1)
		set(handles.edit_cmd,'String', '', 'Tooltip', '' )
	else
		set(handles.edit_cmd,'String', str{val}, 'Tooltip', handles.ttips{val} )
	end

% --------------------------------------------------------------------------------------------------
function push_apply_CB(hObject, handles)
	cmd = get(handles.edit_cmd,'String');
	if (isempty(cmd))
		h = errordlg('Fiu Fiu!! Apply WHAT????','ERROR');
		if (handles.version7 < 8.4 && handles.isPC),	WindowAPI(h, 'TopMost'),	end
		return
	end

	[t, r] = strtok(cmd);
	r = ddewhite(r);	% Remove leading spaces on 'r' that idiot strtok didn't care to do

	jump = false;		% If true, the missing line selection tests below are not executed
	stitch_all = false;
	if (strcmp(t, 'stitch'))	% Test for the case 'stitch TOL ALL' that will stitch all lines
		[tt, rr] = strtok(r);
		rr = ddewhite(rr);		% Ghrr, stupid strtok
		if (strcmpi(rr, 'all'))
			stitch_all = true;
			jump = true;		% Do not execute the missing line selection tests below.
			r = tt;				% This is the TOL val that wi'll need later down
		end
	end
	if (strcmp(t, 'delete'))
		jump = true;
	end

	ind = find(strcmp(t, handles.known_ops));
	if (isempty(ind))
		% Here (will come) a function call which tries to adress the general command issue
		return
	end
	if (isempty(handles.hLine) && ~strcmp(handles.known_ops{ind}, 'pline') && ...
	    ~strcmp(handles.known_ops{ind},'scale') && ~strcmp(handles.known_ops{ind},'GMT_DB') && ...
	    ~strcmp(handles.known_ops{ind},'self-crossings') && ~jump)
		h = errordlg('Fiu Fiu!! Apply WHERE????','ERROR');
		if (handles.version7 < 8.4 && handles.isPC),	WindowAPI(h, 'TopMost'),	end
		return
	end
	if (~strcmp(handles.known_ops{ind},'pline') && ~strcmp(handles.known_ops{ind},'scale') && ...
	    ~strcmp(handles.known_ops{ind},'GMT_DB') && ~strcmp(handles.known_ops{ind},'self-crossings') && ~jump)
		handles.hLine = handles.hLine(ishandle(handles.hLine));
		if (isempty(handles.hLine))
			h = errordlg('Invalid handle. You probably killed the line(s)','ERROR');
			if (handles.version7 < 8.4 && handles.isPC),	WindowAPI(h, 'TopMost'),	end
			return
		end
	end

	switch handles.known_ops{ind}

		case 'bezier'
			n_nodes = validate_args(handles.known_ops{ind}, r);

			vt = linspace(0,1,n_nodes);			bez = zeros(numel(vt), 2);
			x = get(handles.hLine, 'XData');	y = get(handles.hLine, 'YData');
			xy = [x(:) y(:)];

			% Compute each point in the Bezier curve. From Jesus Lucio bezier_.m file in FEX
			np = 1;		n = numel(x) - 1;
			for t = vt
				som = [0 0];
				for (i = 0:n)			% Add the next point multiplied by the corresponding Bernstein polynomial
					som = som + xy(i+1, :) * local_nchoosek(n, i) * (t^i) * ((1 - t)^(n - i));  
				end
				bez(np, :) = som;
				np = np + 1;
			end
			h = line('XData', bez(:,1), 'YData', bez(:,2), 'Parent',handles.hMirAxes, 'Color','r', 'LineWidth',handles.lt,'Tag','polyline');
			draw_funs(h,'line_uicontext')

		case 'buffer'
			[out, msg] = validate_args(handles.known_ops{ind}, r, handles.geog);
			if (~isempty(msg)),		errordlg(msg,'ERROR'),		return,		end

			direction = 'both';		npts = 13;		% Defaults (+ next line)
			geodetic = handles.geog;		% 1 uses spherical aproximation -- OR -- 0 do cartesian calculation
			dist = out.dist;
			if (out.dir),		direction = out.dir;	end
			if (out.npts),		npts = out.npts;		end
			if (out.ab),		geodetic = out.ab;		% Custom ellipsoid
			elseif (out.geod),	geodetic = out.geod;	% WGS84
			end

			if (out.toDegFac)		% Convert to degrees using the simple s = R*theta relation
				dist = (dist * out.toDegFac) / 6371005.076 * 180 / pi;		% Use the Authalic radius
			end
			% Now we will test if user got confused and sent in wrong DIST unites (for geog only)
			if (geodetic)
				lims = getappdata(handles.hMirAxes,'ThisImageLims');
				if ( dist > sqrt( diff(lims(1:2))^2 + diff(lims(3:4))^2)/2 )	% Almost sure a mistake
					msg = sprintf('You probably gave a bad buffer width (DIST) as the buffer line will be out of map.\nShell I stop?');
					r = yes_or_no('string',msg,'title','Warning');
					if (r(1) == 'N'),	return,		end
				end
			end

			for (k = 1:numel(handles.hLine))
				x = get(handles.hLine(k), 'xdata');		y = get(handles.hLine(k), 'ydata');
 				[y, x] = buffer_j(y, x, dist, direction, npts, geodetic);
				if (isempty(x)),	return,		end
				ind = find(isnan(x));
				if (isempty(ind))				% One only, so no doubts
					h = patch('XData',x, 'YData',y, 'Parent',handles.hMirAxes, 'EdgeColor',handles.lc, ...
						'FaceColor','none', 'LineWidth',handles.lt*2, 'Tag','polybuffer');
					draw_funs(h,'line_uicontext')
				else
					% Get the length of the segments by ascending order and plot them in that order.
					% Since the uistack will send the last ploted one to the bottom, that will be the larger polygon
					ind = [0 ind(:)' numel(x)];
					[s, i] = sort(diff(ind));		h = zeros(numel(i),1);
					for (m = 1:numel(i))
						h(m) = patch('XData',x(ind(i(m))+1:ind(i(m)+1)-1), 'YData',y(ind(i(m))+1:ind(i(m)+1)-1), 'Parent',handles.hMirAxes, ...
							 'EdgeColor',handles.lc, 'FaceColor','none', 'LineWidth',handles.lt*2, 'Tag','polybuffer');
					end
					% Do this after so that when it takes time (uistack may be slow) the user won't notice it
					for (m = 1:numel(i))
						uistack_j(h(m),'bottom'),		draw_funs(h(m),'line_uicontext')
					end
					uistack_j(handles.hLine(k),'bottom')	% But keep the processing line at the very bottom
				end
			end

		case 'bspline'
			x = get(handles.hLine(1), 'xdata');		y = get(handles.hLine(1), 'ydata');
			[pp,p] = spl_fun('csaps',x,y);			% This is just to get csaps's p estimate
			y = spl_fun('csaps',x,y,p,x);
			h = line('XData', x, 'YData', y, 'Parent',handles.hMirAxes, 'Color',handles.lc, 'LineWidth',handles.lt, 'Tag','polyline');
			draw_funs(h,'line_uicontext')
			smoothing_param(p, [x(1) x(2)-x(1) x(end)], handles.hMirFig, handles.hMirAxes, handles.hLine, h);

		case 'closing'
			% This shares A LOT of code with 'buffer' and should be merged or something of that kind
			[out, msg] = validate_args(handles.known_ops{ind}, r, handles.geog);
			if (~isempty(msg)),		errordlg(msg,'ERROR'),		return,		end

			direction = 'out';		npts = 13;		% Defaults (+ next line)
			geodetic = handles.geog;		% 1 uses spherical aproximation -- OR -- 0 do cartesian calculation
			dist = out.dist;
			if (out.npts),		npts = out.npts;		end
			if (out.ab),		geodetic = out.ab;		% Custom ellipsoid
			elseif (out.geod),	geodetic = out.geod;	% WGS84
			end

			if (out.toDegFac)		% Convert to degrees using the simple s = R*theta relation
				dist = (dist * out.toDegFac) / 6371005.076 * 180 / pi;		% Use the Authalic radius
			end
			% Now we will test if user got confused and sent in wrong DIST unites (for geog only)
			if (geodetic)
				lims = getappdata(handles.hMirAxes,'ThisImageLims');
				if ( dist > sqrt( diff(lims(1:2))^2 + diff(lims(3:4))^2)/2 )	% Almost sure a mistake
					msg = sprintf('You probably gave a bad buffer width (DIST) as the buffer line will be out of map.\nShell I stop?');
					r = yes_or_no('string',msg,'title','Warning');
					if (r(1) == 'N'),	return,		end
				end
			end

			for (k = 1:numel(handles.hLine))
				x = get(handles.hLine(k), 'xdata');		y = get(handles.hLine(k), 'ydata');
 				[y, x] = buffer_j(y, x, dist, direction, npts, geodetic);
				direction = 'in';				% Now revert the sense of the buffering op
 				[y, x] = buffer_j(y, x, dist, direction, npts, geodetic);
				if (isempty(x)),	return,		end
				ind = find(isnan(x));
				if (isempty(ind))				% One only, so no doubts
					h = patch('XData',x, 'YData',y, 'Parent',handles.hMirAxes, 'EdgeColor',handles.lc, ...
						'FaceColor','none', 'LineWidth',handles.lt*2, 'Tag','polybuffer');
					draw_funs(h,'line_uicontext')
				else
					% Get the length of the segments by ascending order and plot them in that order.
					% Since the uistack will send the last ploted one to the bottom, that will be the larger polygon
					ind = [0 ind(:)' numel(x)];
					[s, i] = sort(diff(ind));		h = zeros(numel(i),1);
					for (m = 1:numel(i))
						h(m) = patch('XData',x(ind(i(m))+1:ind(i(m)+1)-1), 'YData',y(ind(i(m))+1:ind(i(m)+1)-1), 'Parent',handles.hMirAxes, ...
							 'EdgeColor',handles.lc, 'FaceColor','none', 'LineWidth',handles.lt*2, 'Tag','polybuffer');
					end
					% Do this after so that when it takes time (uistack may be slow) the user won't notice it
					for (m = 1:numel(i))
						uistack_j(h(m),'bottom'),		draw_funs(h(m),'line_uicontext')
					end
					uistack_j(handles.hLine(k),'bottom')	% But keep the processing line at the very bottom
				end
			end

		case 'cspline'
			for (k = 1:numel(handles.hLine))
				x = get(handles.hLine(k), 'xdata');		y = get(handles.hLine(k), 'ydata');
				[N, msg] = validate_args(handles.known_ops{ind}, r, numel(x));
				if (~isempty(msg)),		errordlg(msg,'ERROR'),		continue,		end
				np = numel(x);			ind = (1:N(1):np);
				if (x(ind(end)) ~= x(end)),		ind = [ind np];		end		% We want also the last point
				if (numel(N) == 1),			[x,y] = spline_interp(x(ind),y(ind));				% With 10 pts per segment
				else						[x,y] = spline_interp(x(ind),y(ind),'n',N(2));		% With a resolution request
				end
				h = line('XData', x, 'YData', y, 'Parent',handles.hMirAxes, 'Color',handles.lc, 'LineWidth',handles.lt,'Tag','polyline');
				draw_funs(h,'line_uicontext')
			end

		case 'delete'
			% Ok, here already know 't' & 'r' from the strtok at the beguining of this fun
			do_DUP = false;		do_SMALL = false;	do_SPUR = false;
			if (strcmpi(r, 'DUP'))
				do_DUP = true;
			elseif (strncmpi(r, 'SMALL', 5))
				[tt, rr] = strtok(r);		% We are searching here for the 'SMALL N'
				if (~strcmpi(tt, 'SMALL'))
					errordlg('Bad usage. Correct syntax is ''delete SMALL N''','ERROR');	return
				end
				n_vert = round(str2double(rr));	% Minimize chances of stupidity.
				if (isnan(n_vert) || n_vert <= 0)
					errordlg('Bad usage. N must be a valid positive integer','ERROR');	return
				end
				do_SMALL = true;
			elseif (strcmpi(r, 'SPUR'))
				do_SPUR = true;
			elseif (strcmp(r, 'DUP|SMALL'))
				errordlg('Bad usage: DUP|SMALL is not to be used literaly. You must choose one only of DUP or SMALL.', 'Error')
				return
			end
			if (~any([do_DUP do_SMALL do_SPUR]))
				errordlg('Unknown invented option.','ERROR');	return
			end

			hLines = findobj(handles.hMirAxes, 'Type', 'line');
			nLines = numel(hLines);
			if (do_DUP)
				hAguenta = aguentabar(0,'title',sprintf('Removing duplicates from %d lines', nLines));
				n = 0;
				perc = max(round(nLines / 100 * 2), 2);		% Advance every 2% in aguentabar
				t0 = cputime;
				for (k = 1:nLines)
					hCurrLine = hLines(k);
					if (~ishandle(hCurrLine)),	continue,	end
					x = get(hCurrLine, 'XData');	y = get(hCurrLine, 'YData');
					for (m = k+1:nLines)
						if (~ishandle(hLines(m))),	continue,	end
						% Search first for repetition with lines plotted backward one to the other (more likely)
						if (isequal(x(end:-1:1), get(hLines(m), 'XData')) && isequal(y(end:-1:1), get(hLines(m), 'YData')))
							delete(hLines(m)),	n = n + 1;		continue
						end
						if (isequal(x, get(hLines(m), 'XData')) && isequal(y, get(hLines(m), 'YData')))
							delete(hLines(m)),	n = n + 1;		continue
						end
					end
					if (~rem(k, perc))
						hAguenta = aguentabar(k/nLines);
					end
				end
				if (ishandle(hAguenta)),	delete(hAguenta),	end
				t1 = cputime;
				helpdlg(sprintf('Removed %d DUPlicated lines in %.1f sec', n, (t1-t0)),'Info')
			elseif (do_SMALL)
				n = 0;
				for (k = 1:nLines)
					if (numel(get(hLines(k), 'XData')) < n_vert)
						delete(hLines(k)),	n = n + 1;
					end
				end
				helpdlg(sprintf('Removed %d SMALL lines', n),'Info')
			else		% The SPURs case
				if (~isempty(handles.hLine))	% If we have a line selected, use that one only.
					hLines = handles.hLine;		nLines = 1;
				end
				
				for (k = 1:nLines)			% Loop over all lines
					x = get(hLines(k), 'XData');	y = get(hLines(k), 'YData');
					dif_x = x(1:end-2) - x(3:end);
					if (all(dif_x ~= 0))	% No spurs here for shure
						continue
					else
						dif_y = y(1:end-2) - y(3:end);
					end
					ind = find((dif_x == 0) & (dif_y == 0));
					for (m = 1:numel(ind))
						n1 = ind(m) - 1;	n2 = ind(m) + 1;	c = 0;
						while ((n1 > 0 && n2 <= numel(x)) && (dif_x(n1) == 0) && (dif_x(n2) == 0) && ...
						       (dif_y(n1) == 0) && (dif_y(n2) == 0))
						   n1 = n1 + 1;	   n2 = n2 + 1;
						   c  = c + 1;
						end
						if (c)			% Means we have a SPUR here
							x((ind(m) - c):(ind(m) + c)) = [];		y((ind(m) - c):(ind(m) + c)) = [];
							if (isempty(x))
								delete(hLines(k))		% If it is empty, better delete it
							else
								set(hLines(k), 'XData', x);		set(hLines(k), 'YData', y);
							end
						end
					end
				end
				resetSemaf(handles, hObject)	% Make it clear that it is neccessary to explicitly pick another line
				
			end

		case 'group'
			if ( ~strcmp(get(handles.hLine(1),'Type'),'line') )
				h = warndlg('The selected object is not of type LINE. Objects of type PATCH are not currently "groupable"','Warning');
				if (handles.version7 < 8.4 && handles.isPC),	WindowAPI(h, 'TopMost'),	end
				return
			end
			lt = get(handles.hLine, 'LineWidth');		lc = get(handles.hLine, 'Color');
			ls = get(handles.hLine, 'LineStyle');		tag = get(handles.hLine, 'Tag');
			hLines = findobj(handles.hMirAxes, 'Type', 'line');
			hLines = setxor(hLines, handles.hLine(1));
			did_one = false;		% To tell us if at least 2 lines were grouped
			x = get(handles.hLine(1), 'xdata');			y = get(handles.hLine(1), 'ydata');
			for (k = 2:numel(hLines))
				if ( strcmp(get(hLines(k),'Tag'), tag)  && strcmp(get(hLines(k),'LineStyle'), ls) && ...
					isequal(get(hLines(k),'Color'), lc) && isequal(get(hLines(k),'LineWidth'), lt) )
					xi = get(hLines(k), 'xdata');		yi = get(hLines(k), 'ydata');
					x = [x NaN xi];						y = [y NaN yi];
					delete(hLines(k))
					did_one = true;
				end
			end
			if (did_one),	set(handles.hLine(1), 'XData',x, 'YData',y),	end

		case 'hand2Workspace'
			assignin('base','lineHandles',handles.hLine);

		case 'line2patch'
			for (k = 1:numel(handles.hLine))
				if (~strcmp(get(handles.hLine(k),'type'), 'line')),	continue,	end
				x = get(handles.hLine(k), 'xdata');		y = get(handles.hLine(k), 'ydata');
				h = patch('XData',x, 'YData',y, 'Parent',handles.hMirAxes, 'EdgeColor',get(handles.hLine(k),'color'), ...
					'FaceColor','none', 'LineWidth',get(handles.hLine(k),'LineWidth'), 'Tag','line2patch');
				draw_funs(h,'line_uicontext')
				delete(handles.hLine(k))
			end

		case 'polysimplify'
			[tol, msg] = validate_args(handles.known_ops{ind}, r);
			if (~isempty(msg)),		errordlg(msg,'ERROR'),		return,		end
			for (k = 1:numel(handles.hLine))
				x = get(handles.hLine(k), 'xdata');		y = get(handles.hLine(k), 'ydata');
				if (handles.geog),		B = cvlib_mex('dp', [x(:) y(:)], tol, 'GEOG');
				else					B = cvlib_mex('dp', [x(:) y(:)], tol);
				end
				h = line('XData',B(:,1), 'YData',B(:,2), 'Parent',handles.hMirAxes, 'Color',handles.lc, 'LineWidth',handles.lt, 'Tag','polyline');
				draw_funs(h,'line_uicontext')
			end

		case {'polyunion' 'polyintersect' 'polyxor' 'polyminus'}
			if		(handles.known_ops{ind}(5) == 'u'),		ID = 3;
			elseif	(handles.known_ops{ind}(5) == 'i'),		ID = 1;
			elseif	(handles.known_ops{ind}(5) == 'x'),		ID = 2;
			else											ID = 0;
			end
			[lixo, msg] = validate_args(handles.known_ops{ind}, numel(handles.hLine));
			if (~isempty(msg)),		errordlg(msg,'ERROR'),		return,		end
			P1.x = get(handles.hLine(1), 'xdata');			P1.y = get(handles.hLine(1), 'ydata');
			for (k = 2:numel(handles.hLine))
				P2.x = get(handles.hLine(k), 'xdata');		P2.y = get(handles.hLine(k), 'ydata');
				P1 = PolygonClip(P1, P2, ID);
			end
			for (k = 1:length(P1))
				h = patch('XData', P1(k).x, 'YData', P1(k).y, 'Parent',handles.hMirAxes, 'EdgeColor',handles.lc, ...
					'FaceColor','none', 'LineWidth',handles.lt, 'Tag','polyclip');
				draw_funs(h,'line_uicontext')
			end

		case 'pline'
			[out, msg] = validate_args(handles.known_ops{ind}, r);
			if (~isempty(msg)),		errordlg(msg,'ERROR'),		return,		end
			try
				[X{1:out{3}}] = strread(out{1}, repmat('%f',1,out{3}));		% Parse the X coords
				[Y{1:out{3}}] = strread(out{2}, repmat('%f',1,out{3}));		% Parse the Y coords
				x = cat(1,X{:});			y = cat(1,Y{:});
				h = line('XData',x, 'YData',y, 'Parent',handles.hMirAxes, 'Color',handles.lc, 'LineWidth',handles.lt, 'Tag','polyline');
				draw_funs(h,'line_uicontext')
			catch
				h = errordlg(['Something screwd up. Error message is ' lasterr]);
				if (handles.version7 < 8.4 && handles.isPC),	WindowAPI(h, 'TopMost'),	end
			end

		case 'scale'
			xMin = 1e50;	yMin = 1e50;	xMax = -xMin;	yMax = -yMin;
			hLines = findobj(handles.hMirAxes, 'Type', 'line');
			for (k = 1:numel(hLines))
				x = get(hLines(k), 'XData');	y = get(hLines(k), 'YData');
				xMin = min([xMin min(x)]);		xMax = max([xMax max(x)]);
				yMin = min([yMin min(y)]);		yMax = max([yMax max(y)]);
			end
			hPatch = findobj(handles.hMirAxes, 'Type', 'patch');
			for (k = 1:numel(hPatch))
				x = get(hPatch(k), 'XData');	y = get(hPatch(k), 'YData');
				xMin = min([xMin min(x)]);		xMax = max([xMax max(x)]);
				yMin = min([yMin min(y)]);		yMax = max([yMax max(y)]);
			end

			h = mirone;		hNewMirHand = guidata(h);
			hNewMirHand = mirone('FileNewBgFrame_CB', hNewMirHand, [-0.5 0.5 -0.5 0.5 0 1], [], 'Whatever');
			set(hNewMirHand.figure1, 'Name', 'GMT custom symbol')
			scale_x = 1 / (xMax - xMin);	off_x = xMin;
			scale_y = 1 / (yMax - yMin);	off_y = yMin;
			scale = scale_y;
			if (scale_x < scale_y),		scale = scale_x;	end

			% Now a second run to efectively scale the data
			for (k = 1:numel(hLines))
				x = (get(hLines(k), 'XData') - off_x) * scale - 0.5;
				y = (get(hLines(k), 'YData') - off_y) * scale - 0.5;
				h = line('XData',x, 'YData',y, 'Parent', hNewMirHand.axes1, 'LineWidth', get(hLines(k),'LineWidth'), ...
					'Color', get(hLines(k),'Color'), 'LineStyle', get(hLines(k),'LineStyle'), 'Tag', get(hLines(k),'Tag'));
				draw_funs(h,'line_uicontext')
			end
			for (k = 1:numel(hPatch))
				x = (get(hPatch(k), 'XData') - off_x) * scale - 0.5;
				y = (get(hPatch(k), 'YData') - off_y) * scale - 0.5;
				h = patch('XData',x, 'YData',y, 'Parent', hNewMirHand.axes1, 'LineWidth', get(hPatch(k),'LineWidth'), ...
					'EdgeColor', get(hPatch(k),'EdgeColor'), 'FaceColor', get(hPatch(k),'FaceColor'), ...
					'LineStyle', get(hPatch(k),'LineStyle'), 'Tag', get(hPatch(k),'Tag'));
				draw_funs(h,'line_uicontext')
			end

		case 'self-crossings'
			if (isempty(handles.hLine))		% If no line selected, use them all
				hLines = findobj(handles.hMirAxes, 'Type', 'line');
				hPatches = findobj(handles.hMirAxes, 'Type', 'patch');
				hLines = [hLines(:); hPatches(:)];
			else
				hLines = handles.hLine;
			end
			for (k = 1:numel(hLines))
				x = get(hLines(k), 'XData');	y = get(hLines(k), 'YData');
				n_pts = numel(x);	robust = 1;
				if (n_pts < 2),			continue
				elseif (n_pts > 1000),	robust = 1000;		% Do chunking otherwise MEM Kaboom
				end
				[xc, yc, iout] = intersections(x, y, robust);
				if (~isempty(iout) && (iout(1) - 1) < 1e-5),	xc(1) = [];		yc(1) = [];		end		% Avoid annoying false +
				if (~isempty(yc))
					h = line('XData',xc, 'YData',yc, 'Parent',handles.hMirAxes, 'LineStyle','none', 'Marker','o', ...
						'MarkerEdgeColor','k','MarkerFaceColor','y', 'MarkerSize',7, 'Tag','self-crossings');
					draw_funs(h,'line_uicontext')
				end
			end
			resetSemaf(handles, hObject)	% reset line selection to empty and and red color in semaforo icon

		case 'stitch'
			out = validate_args(handles.known_ops{ind}, r);
			tol = out.val;
			if (out.toDegFac)		% Convert to degrees using the simple s = R*theta relation
				tol = (out.val * out.toDegFac) / 6371005.076 * 180 / pi;		% Use the Authalic radius
			end

			hCurrLine = handles.hLine;		xcell = [];
			hLines = findobj(handles.hMirAxes, 'Type', 'line');
			if (numel(hLines) == 1)				% Only one line in the whole plot. Check if it's a multi-segment
				[xcell, ycell] = polysplit(get(hCurrLine,'XData'), get(hCurrLine,'YData'));
				if (isempty(ycell))
					warndlg('There is only one line in Town and with no NaNs breaking. So stitch where?','Warning')
					return
				end
			end

			if (stitch_all)
				for (k = 1:numel(hLines))
					hCurrLine = hLines(k);
					if (~ishandle(hCurrLine)),	continue,	end		% Means this line chunk was already stitched
					hLines_t = findobj(handles.hMirAxes, 'Type', 'line');
					hLines_t = setxor(hLines_t, hCurrLine);
					do_stitching(hLines_t, hCurrLine, tol, xcell);
				end
			else
				hLines = setxor(hLines, hCurrLine);
				do_stitching(hLines, hCurrLine, tol, xcell);
			end

			resetSemaf(handles, hObject)	% Make it clear that it is neccessary to explicitly pick another line

		case 'thicken'
			[out, msg] = validate_args(handles.known_ops{ind}, r, handles.hMirAxes);
			if (~isempty(msg)),		errordlg(msg,'ERROR'),		return,		end
			N = out(1);		hscale = out(2);	vscale = out(3);
			x = get(handles.hLine(1), 'xdata');		y = get(handles.hLine(1), 'ydata');
			
			if (handles.geog)
				[dumb, az] = vdist(y(1:end-1),x(1:end-1), y(2:end),x(2:end));
				az = 90 - az;		% Make it trigonometric
				co = cos(az * pi / 180);	si = sin(az * pi / 180);
			else
				dx = diff(x);				dy = diff(y);
				% calculate the cosine and sine
				hy = (dx.^2 + dy.^2)^.5;
				co =  dx ./ hy;				si =  dy ./ hy;
			end

			thick = N * (handles.head(8) + handles.head(9)) / 2;		% desenrasque, mas foleiro
			% rotate a control line, "th" cells long, based on the slope of the line
% 			foo = [co -si; si  co] * [0	 0; thick/2 -thick/2];
			% Add rotated points to line vertices
			n_pts = numel(x);
			x = x(:);			y = y(:);		% We need them as column vectrs
			x_copy = x;			y_copy = y;
% 			x = [x; x(end:-1:1)];		y = [y; y(end:-1:1)];		% Wrap arround
% 			x(1:n_pts) = x(1:n_pts) + repmat(foo(1,1),n_pts,1);
% 			x(n_pts+1:end) = x(n_pts+1:end) + repmat(foo(1,2),n_pts,1);
% 			y(1:n_pts) = y(1:n_pts) + repmat(foo(2,1),n_pts,1);
% 			y(n_pts+1:end) = y(n_pts+1:end) + repmat(foo(2,2),n_pts,1);
% 			x(end+1) = x(1);			y(end+1) = y(1);			% Close line

			scale = (hscale + vscale) / 2;		% Approximation that works relatively well
 			set(handles.hLine(1), 'LineWidth', thick / scale, 'Tag', 'cellthick')
			refresh
			hui = findobj(get(handles.hLine(1),'UIContextMenu'),'Label', 'Extract profile');
			set(hui, 'Call', 'setappdata(gcf,''StackTrack'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo))')
			try		rmappdata(handles.hMirFig,'TrackThisLine'),		end		% Clear it so that ExtractProfile_CB() in mirone.m knows the way

% 			h = line('XData',x', 'YData',y', 'Parent',handles.hMirAxes, 'Color','w', 'LineWidth',handles.lt, 'Tag','thickned');
% 			draw_funs(h,'line_uicontext')

			% Compute the N+1 lines that will be attached to this thickned line
			dl = thick / N;
			xL = zeros(n_pts, N+1);		yL = zeros(n_pts, N+1);
% 			co = [co co(end)];			si = [si si(end)];
% 			foo1 = zeros(n_pts, 1);		foo2 = zeros(n_pts, 1);
			for (k = 1:N+1)
				th = thick/2 - (k-1) * dl;
				foo = [co(1) -si(1); si(1) co(1)] * [0; th];
				xL(:,k) = x_copy + repmat(foo(1), n_pts, 1);		% One line per column
				yL(:,k) = y_copy + repmat(foo(2), n_pts, 1);
				
				% Tenho de voltar a remoer isto. O problema (em baixo) é que depois não vai ter
				% o mesmo número de elementos quando interpolada no grid_profiler (para stakar)
% 				for (m = 1:n_pts)		% This loop has meaning only for polylines, as angles vray between segments
% 					foo = [co(m) -si(m); si(m) co(m)] * [0; th];
% 					foo1(m) = foo(1);	foo2(m) = foo(2);
% 				end
% 				xL(:,k) = x_copy + foo1;		% One line per column
% 				yL(:,k) = y_copy + foo2;
			end
			set(handles.hLine(1), 'UserData', {xL; yL; thick; [x_copy y_copy]; handles.geog; ...	% We'll next info if line is edited
					'MxN X and Y with M = number_vertex and N = number_lines; THICK = thickness in map units; Mx2 = original line; is geog?'})
			
		case 'toRidge'
			hanMir = guidata(handles.hMirFig);
			[out, msg] = validate_args(handles.known_ops{ind}, r, hanMir);
			if (~isempty(msg)),		errordlg(msg,'ERROR'),		return,		end
			N = out(1);
			[X,Y,Z,head] = load_grd(hanMir);
			if (isempty(Z)),	return,		end
			x = get(handles.hLine(1), 'xdata');		y = get(handles.hLine(1), 'ydata');
			dx = N * head(8);			dy = N * head(9);
			f_name = [handles.path_tmp 'lixo.dat'];		% Temp file
			xR = x * NaN;				yR = y * NaN;
			hR = line('XData',xR, 'YData',yR, 'Parent',handles.hMirAxes,'Linewidth', ...
				get(handles.hLine(1),'Linewidth')+1, 'Color',(1 - get(handles.hLine(1),'Color')));
			draw_funs(hR,'line_uicontext')	% Set edition functions
			for (k = 1:numel(x))		% Loop over vertex to find its closest position on the Ridge
				xRec = [max(head(1), x(k)-dx) min(head(2), x(k)+dx)];
				yRec = [max(head(3), y(k)-dy) min(head(4), y(k)+dy)];
				rect = [xRec(1) xRec(1) xRec(2) xRec(2) xRec(1); yRec(1) yRec(2) yRec(2) yRec(1) yRec(1)];
				% Get a small grid arround the current point
				[X,Y,Z,hdr] = mirone('ImageCrop_CB', hanMir, rect, 'CropaGrid_pure');
				out = grdppa_m(Z, hdr);
				if ( ~all(isnan(out(1,:))) )
					double2ascii(f_name,out','%f\t%f');					% Save as file so we can use it mapproject
					pt = c_mapproject([x(k) y(k)], ['-L' f_name]);		% Project and get points along the line
				else
					pt = [0 0 Inf];		% Make it fall into the next 'else' case
				end
				if ( ~isinf(pt(3)) )
					xR(k) = pt(4);		yR(k) = pt(5);
				else
					xR(k) = x(k);		yR(k) = y(k);
				end
				set(hR, 'XData', xR, 'YData', yR)
			end
			%clear mapproject_m			% Because of the memory leaks
			set(hanMir.figure1,'pointer','arrow')

		case 'GMT_DB'
			update_GMT_DB(handles)
	end

% -------------------------------------------------------------------------------------------------
function resetSemaf(handles, hObject)
% Make it clear that it is neccessary to explicitly pick another line
	if (~isempty(handles.hLine))
		set(handles.push_semaforo,'BackgroundColor',[1 0 0])
		set(hObject,'Tooltip', 'Have 0 lines to play with')
		handles.hLine = [];		guidata(handles.figure1, handles)
	end

% -------------------------------------------------------------------------------------------------
function update_GMT_DB(handles, TOL)
% Take the GMT database polygons displayed in the figure and update them with the
% data from other ordinary polylines also present in the figure. The updating is done
% after a guessing, via distance calculations, of what is closest to which.
% End points of ordinary plines further away than DS (below) are discarded as potential
% updaters to the GMT_DBpolyline line types.
%
% It would be nice to search for the closest GMT_DB / Updater pair and not just the first
% one that is found, but I'm not sure it worth the effort. Some day perhaps.

	ds = 0.02;			% Endpoints must be closer than this to the GMD_DB polygones
	Is = [];			% Will be used to know if nothing found/done
	hGMT_DB = findobj(handles.hMirAxes, 'tag','GMT_DBpolyline');
	hLines = findobj(handles.hMirAxes, 'type', 'line');
	hLines = setxor(hGMT_DB, hLines);	% The hGMT_DB was repeated as they are also of type 'line'
	for (k = 1:numel(hLines))
		[x, y] = check_bombordo(hLines(k));
		%is_UPDATERclosed = (x(1) == x(end) && y(1) == y(end));
		for (n = 1:numel(hGMT_DB))
			xG = get(hGMT_DB(n), 'XData');		yG = get(hGMT_DB(n), 'YData');
			%is_DBclosed = (xG(1) == xG(end) && yG(1) == yG(end));		% YES for the coastlines, but not forcedly for others
			x0 = x(1) - ds;		x1 = x(1) + ds;
			y0 = y(1) - ds;		y1 = y(1) + ds;
			ind = ( xG > x0 & xG < x1 & yG > y0 & yG < y1);
			% 'ind' is not empty when the updating line is close to old one
			if (any(ind))
				id = find(ind);			% Index of points inside the searching region
				dist = vdist(y(1), x(1), yG(id), xG(id), [6378137, 0, 1/298.2572235630]);
				[minS, Is] = min(dist);
				Is = Is + id(1) - 1;

				% Now find the index of the last point of the chunk to be replaced
				x0 = x(end) - ds;		x1 = x(end) + ds;
				y0 = y(end) - ds;		y1 = y(end) + ds;
				ind = ( xG > x0 & xG < x1 & yG > y0 & yG < y1);
				id = find(ind);			% Index of points inside the searching region
				if (isempty(id))
					x0 = xG(end) - ds;		x1 = xG(end) + ds;
					y0 = yG(end) - ds;		y1 = yG(end) + ds;					
					ind = ( x > x0 & x < x1 & y > y0 & y < y1);
					id = find(ind);
					if (isempty(id)),	Is = [];	continue;	end
				end
				dist = vdist(y(end), x(end), yG(id), xG(id), [6378137, 0, 1/298.2572235630]);
				[minE, Ie] = min(dist);
				Ie = Ie + id(1) - 1;
				if (Ie < Is)
					h = warndlg(['The updating polygon had a wrong orientation non detected by the test,' ...
						' or the test failled miserably. In either case, polygon was ignored.'],'Warning');
					if (handles.version7 < 8.4 && handles.isPC),	WindowAPI(h, 'TopMost'),	end
					continue
				end
				xG = [xG(1:Is-1) x xG(Ie+1:end)];
				yG = [yG(1:Is-1) y yG(Ie+1:end)];
				set(hGMT_DB(n), 'XData', xG, 'YData', yG)
				setappdata(hGMT_DB(n),'edited',true)	% Kind of flag to guide during the saving step
				delete(hLines(k))					% Delete this updatter line since it was used already
				break
			end
		end
	end
	if (isempty(Is))
		warndlg(['Nada, Nothing, Nickles. No polygons updated. Possibly the ends points of the updating ' ...
				'polyline was not within 0.02 degrees of the to be updated polygon'],'Warning')
	end

function [x, y] = check_bombordo(hLine)
% Get the data points from the line handle and check that they are "bombordo" (CCW) oriented
% as required by the GMT database coastlines. If not, reverse them.

	x = get(hLine, 'XData');	y = get(hLine, 'YData');
	xt = x(1:2:end-1);			yt = y(1:2:end-1);	% Ensure testing polygon is not closed and save resources
	xt = xt - mean(xt);
	n = numel(xt);
	i = [2:n 1];
	j = [3:n 1 2];
	k = (1:n);
	a = sum(xt(i) .* (yt(j) - yt(k)));
	if (a < 0)				% If area is negative data is CW oriented (leave land on estibordo)
		x = x(end:-1:1);	% and we have to revert it
		y = y(end:-1:1);
	end	
% -------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------------------------------
function do_stitching(hLines, hCurrLine, tol, xcell)
% Core function to do the line stitching work

	nLines = numel(hLines);				doAguenta = false;	hAguenta = [];	x = [];		y = [];
	if (nLines == 0)					% When one single hLine holding a multi-segment pline
		nLines = numel(xcell);
	end

	if (nLines > 100 )					% greater than 100?
		hAguenta = aguentabar(0,'title',sprintf('Stitching %d lines', nLines));
		perc = max(round(nLines / 100 * 2), 2);		% Advance every 2% in aguentabar
		doAguenta = true;
	end

	% If we have a multi-segment situation, deal with it right away and at the end short-circuit the line handles case
	if (~isempty(xcell))
		nCycles = numel(xcell)-1;		% -1 because the first is used in first arg to find_closestline()
		for (k = 1:nCycles)
			hLines = {xcell(2:end) ycell(2:end)};
			[hLineClosest, endType, indOfFound] = find_closestline({xcell(1) ycell(1)}, hLines, tol);
			if (~indOfFound)
				continue
			else
				indOfFound = indOfFound + 1;	% Because we started at second segment
			end
			x1 = xcell{1};				y1 = ycell{1};
			x2 = xcell{indOfFound};		y2 = ycell{indOfFound};
			if (endType == 1)			% Lines grow in oposite directions from a "mid point"
				x = [x2(end:-1:1) x1];	y = [y2(end:-1:1) y1];
			elseif (endType == 2)		% Line 2 ends near the begining of line 1 
				x = [x2 x1];			y = [y2 y1];
			elseif (endType == 3)		% Line 1 ends near the begining of line 2
				x = [x1 x2];			y = [y1 y2];
			else						% Lines grow from the extremeties twards the "mid point"
				x = [x1 x2(end:-1:1)];	y = [y1 y2(end:-1:1)];
			end
			xcell{1} = x;				ycell{1} = y;	% Equivallent of the below set(hCurrLine, 'XData',x, 'YData',y)
			xcell{indOfFound} = [];		ycell{indOfFound} = [];		% Do not process the same segment again.
			if (doAguenta && ~rem(k, perc)),
				hAguenta = aguentabar(k/nLines);
			end
		end
		set(hCurrLine, 'XData',x, 'YData',y)
		nLines = 0;		% TO PREVENT EXECUTING ALSO THE NEXT LOOP
	end

	for (k = 1:nLines)					% Loop over all others but selected line. FOR NON MULTI-SEGMENT CASES
		[hLineClosest, endType, indOfFound] = find_closestline(hCurrLine, hLines, tol);
		if (isempty(hLineClosest))		% Either we found and finished or found nothing
			break
		end
		x1 = get(hCurrLine,'XData');	y1 = get(hCurrLine,'YData');
		x2 = get(hLineClosest,'XData');	y2 = get(hLineClosest,'YData');
		if (endType == 1)				% Lines grow in oposite directions from a "mid point"
			x = [x2(end:-1:1) x1];		y = [y2(end:-1:1) y1];
		elseif (endType == 2)			% Line 2 ends near the begining of line 1 
			x = [x2 x1];				y = [y2 y1];
		elseif (endType == 3)			% Line 1 ends near the begining of line 2
			x = [x1 x2];				y = [y1 y2];
		else							% Lines grow from the extremeties twards the "mid point"
			x = [x1 x2(end:-1:1)];		y = [y1 y2(end:-1:1)];
		end
		set(hCurrLine, 'XData',x, 'YData',y)
		delete(hLineClosest)			% It was assimilated, now delete old one.
		hLines(indOfFound) = [];		% Do not process the same line again.
		if (doAguenta && rem(k,perc)),	hAguenta = aguentabar(k/nLines);		end
	end

	% Now search for repeated points along the stitched line (nothing uncommon)
	ind_x = (diff(x) ~= 0);		ind_y = (diff(y) ~= 0);
	unicos = (ind_x | ind_y);
	if (~isempty(x))					% Prevent case of an initial x = []
		unicos(end+1) = true;
		x = x(unicos);			y = y(unicos);
		set(hCurrLine, 'XData',x, 'YData',y)
	end

	% Files (imported with ogr) can have each line segment repeated which causes strange effects. 
	% This corrects one type of those effects, which is hard to describe with words.
	knees = zeros(2,1);			nKnees = 0;
	for (k = 1:numel(x)-2)		% Do it scalar because is equaly fast (or more) and waste no memory
		if ( (x(k) == x(k+2)) && (y(k) == y(k+2)) )	% Detect points where line returns along same path.
			nKnees = nKnees + 1;
			knees(nKnees) = k;	
		end
	end
	if (nKnees == 2)
		x(knees(1)+1:knees(2)+1) = [];		y(knees(1)+1:knees(2)+1) = [];
		set(hCurrLine, 'XData',x, 'YData',y)
	end

	% We still need to check if first and last pts are whithin TOL, case in which line is closed.
	if (~isempty(x) && (x(1) ~= x(end)) && (y(1) ~= y(end)) && (sqrt((x(1) - x(end))^2 + (y(1) - y(end))^2) <= tol))
		x(end+1) = x(1);		y(end+1) = y(1);
		set(hCurrLine, 'XData',x, 'YData',y)
	end
	if (ishandle(hAguenta)),	delete(hAguenta),	end

% -------------------------------------------------------------------------------------------------
function [hLineClosest, endType, indOfFound] = find_closestline(hMe, hLines, TOL)
% Find among the HLINES vector of line handles which one is closest the line with handle HME.
% TOL is the max distance that the two lines can be apart and still be considered close.
%	If not provided, defaults to Inf.
% HLINECLOSEST	Is the handle of the closest line or [] if any passes the TOL condition.
% ENDTYPE		Is a integer ranging between 1 and 4, or 0 if no closest line is found. Where
%				1	Lines grow in oposite directions from a "mid point"
%				2	Line 2 (HLINECLOSEST) ends near the begining of line 1 (HME)
%				3	Line 1 (HME) ends near the begining of line 2 (HLINECLOSEST)
%				4	Lines grow from the extremeties twards the "mid point"
% INDOFFOUND	Index of HLINECLOSEST in original HLINES vector. That index of closest line.
%
% Alternativelly, HME and HLINES can be a 1x2 cell array of cells, where each cell of the
% outer cell container holds the xx and yy coordinates of each segment. Example:
% hMe={{1 1} {2 1}};	hLines = {{get(line(1),'XData'), get(line(2),'XData')} { "Same for Y here" }}
% The contents of hLines in cell form is normally created from the output of polysplit()

	if (nargin == 2),	TOL = inf;		end
	if (ishandle(hMe))
		x1 = get(hMe,'XData');		y1 = get(hMe,'YData');
		nLines = numel(hLines);
	else			% Not tested but it must be a cell array ... of cells
		x1 = hMe{1}{1};				y1 = hMe{2}{1};
		nLines = numel(hLines{1});
	end

	hLineClosest = [];		endType = 0;	indOfFound = 0;
	minDist = TOL * TOL;	% We do the test with square distances to save calls to sqrt 

	for (k = 1:nLines)
		if (ishandle(hLines(1)))
			x2 = get(hLines(k),'XData');		% These calls are horribly expensive, so ... (what?)
		else
			x2 = hLines{1}{k};
			if (isempty(x2)),	continue,	end	% It's ok, we accept empty cells
		end
		dif_x2 = ([(x1(1) - x2(1)); (x1(1) - x2(end)); (x1(end) - x2(1)); (x1(end) - x2(end))]).^2;
		if (min(dif_x2) > minDist),		continue,	end		% we know enough to drop this line.
		if (ishandle(hLines(1))),		y2 = get(hLines(k),'YData');
		else							y2 = hLines{2}{k};
		end
		if ( (x2(1) == x2(end)) && (y2(1) == y2(end)) )		% Ignore closed polygons
			continue
		end

		dif_y = [(y1(1) - y2(1)); (y1(1) - y2(end)); (y1(end) - y2(1)); (y1(end) - y2(end))];
		dist = sum([dif_x2 dif_y.^2], 2);	% Square of distances between the 4 extremities
		[mimi, I] = min(dist);
		if ( mimi <= minDist )
			indOfFound = k;
			endType = I;
			if (ishandle(hLines(1))),		hLineClosest = hLines(k);
			else							hLineClosest = k;
			end
			if (mimi < 1e-8),	break,	end		% This is so small that we can assume to have found the closest
		end
	end

% --------------------------------------------------------------------------------------------
function res = local_nchoosek(n ,k)   % A faster stripped alternative of nchoosek(n, k)
	res = 1;
	if (k == 0),	return,		end
	if (n == 0),	res = 0;	return,	end
	n_k = n - k;
	for (i = 1:k)
		res = res * ((n_k + i)/i);
	end

% ------------------------------------------------------------------------------------------
function [xcell, ycell] = polysplit(x, y)
% Return the NaN-delimited segments of the vectors X and Y as cell arrays.
% Each element of the cell array contains one segment.

	%  Make sure no NaNs at begining or end.
	% Hopefully well behaved data won't have them and so no memory waist (no local copy)
	while ( isnan(x(1)) || isnan(y(1)) )
		x(1) = [];		y(1) = [];
	end
	while ( isnan(x(end)) || isnan(y(end)) )
		x(end) = [];	y(end) = [];
	end

	indx = find(isnan(x));			% Find were are the NaNs 
	if ( ~isequal(indx, find(isnan(y))) )
		error('X and Y vectors do not have the NaNs in same places')
	end

	ind = false(1, numel(indx)+1);
	xcell = cell(1,numel(indx)+1);		ycell = cell(1,numel(indx)+1);
	xcell{1} = x(1:indx(1)-1);			ycell{1} = y(1:indx(1)-1);
	for (i = 2:numel(indx))				% Do midle chunks (if any)
		xcell{i} = x(indx(i-1)+1:indx(i)-1);
		ycell{i} = y(indx(i-1)+1:indx(i)-1);
		if (isempty(xcell{i})),		ind(i) = true;		end		% They may be empty if consecutive NaNs
	end
	if (numel(indx) > 1)				% Now, if there is one, do the last chunk
		xcell{end} = x(indx(end)+1:end);		ycell{end} = y(indx(end)+1:end);
	end
	
	if (any(ind))						% Remove eventual empty ones
		ycell(ind) = [];		xcell(ind) = [];
	end

% --------------------------------------------------------------------------------------------------
function [out, msg] = validate_args(qual, str, np)
% Check for errors on arguments to operation QUAL and return required args, or error
% WARN: the NP arg has not always the same meaning. For 'buffer' it holds the handles.geog
	msg = [];	out = [];
	switch qual
		case {'buffer' 'closing'}
			out.dir = false;	out.npts = false;		out.ab = false;		out.geod = false;
			convFrom = false;	out.toDegFac = false;	% If we are to convert meters, kilometers or nmiles to degrees
			[t, str] = strtok(str);
			if (strcmp(t, 'DIST'))
				msg = 'The argument "DIST" must be replaced by a numeric value representing the width of the buffer zone';	return
			end
			if (np)		% Here 'np' is actually the handles.geog
				if ( strcmpi(t(end), 'm') || strcmpi(t(end), 'k') || strcmpi(t(end), 'n') )	% Convert any of those to degrees
					convFrom = lower(t(end));		t(end) = [];
				end
			end
			out.dist = abs( str2double(t) );
			if (isnan(out.dist)),		msg = 'BUFFER argument is nonsense';	return,		end
			if (convFrom)		% Right. We must convert the DIST distance given in cartesians to degrees
				% But the problem is that here we don't know the lat so we can't finish the conversion. So, send back the collected info
				if (convFrom == 'm'),		out.toDegFac = 1;
				elseif (convFrom == 'k'),	out.toDegFac = 1000;
				else						out.toDegFac = 1852;
				end
			end

			% OPTIONS
			str = strrep(str,'''','');		% no ' ' around the char variables
			ind = strfind(str, 'in');
			if (~isempty(ind))
				out.dir = str(ind:ind+1);		str(ind:ind+1) = [];
			end
			ind = strfind(str, 'out');
			if (~isempty(ind))
				out.dir = str(ind:ind+2);		str(ind:ind+2) = [];
			end

			% Ok, here we are over with the 'in' or 'out' options. Proceed
			ind = strfind(str, 'geod');
			if (~isempty(ind))
				out.geod = true;		str(ind:ind+3) = [];		% WGS-84
			end

			% See if we have an ellipsoid
			ind1 = strfind(str,'[');
			if (~isempty(ind1) && numel(ind1) ~= 1),	msg = 'Wrong syntax. ''['' symbol must appear one and only one time';	return,	end
			ind2 = strfind(str,']');
			if (~isempty(ind1) && numel(ind2) ~= 1),	msg = 'Wrong syntax. '']'' symbol must appear one and only one time';	return,	end

			% test the geodetic option
			if (ind1)
				[t1, r] = strtok(str(ind1+1:ind2-1));	t2 = strtok(r);
				a = str2double(t1);		b = str2double(t2);		% b actually may be f (flatening)
				if (isnan(a) || isnan(b))
					msg = 'ellipsoid vector is screwed up. Please revise or pay more attention';	return
				end
				out.ab = [a b];
				str(ind1:ind2) = [];		% Remove it from the opt string
			end
			
			% We still may have the circle NPTS option 
			np = round(abs( str2double(strtok(str)) ));
			if (~isnan(np)),	out.npts = np;		end			

		case 'polysimplify'
			d = strtok(str);
			if (strcmp(d, 'TOL'))
				msg = 'The argument "TOL" must be replaced by a numeric value representing the tolerance';	return
			end
			out = abs( str2double(d) );
			if (isnan(out)),		msg = 'POLYSIMPLIFY argument is nonsense';		end

		case 'cspline'
			[N, r] = strtok(str);
			if (isempty(N))
				msg = 'Must provide N (decimation interval)';		return
			end
			out = round(abs( str2double(N) ));
			if (isnan(out))
				if (N(1) == 'N')
					msg = 'The argument "N" is no to be taken literaly. It must contain the decimation interval';
				else
					msg = 'cspline argument is nonsense';
				end
				return
			elseif (out > np)
				msg = 'N is larger than the number of polyline vertices. Smoothed line would desapear';
			end
			N = round(abs( str2double(strtok(r)) ));		% See if we have resolution request
			if (~isnan(N)),		out(2) = N;		end

		case 'bezier'
			if ( isempty(str) || strcmpi(str, 'N') )
				out = 100;
			else
				out = round(abs(str2double(str)));
			end

		case {'polyunion' 'polyintersect' 'polyxor' 'polyminus'}
			if (str < 2)		% str is in fact the number of handles
				msg = 'you want to make the union of a single line???. Wierd!';
			end

		case 'thicken'
			N = strtok(str);
			if (isempty(N)),	N = '10';		end		% Default value
			out = round(abs( str2double(N) ));
			if (isnan(out))
				if (N(1) == 'N')
					msg = 'The argument "N" is no to be taken literaly. It must contain the number of grid cell to thicken line';
				else
					msg = 'thicken argument is nonsense';
				end
				return
			end
			hMirAxes = np;
         	axLims = getappdata(hMirAxes,'ThisImageLims');
			% create a conversion from data to points for the current axis
			oldUnits = get(hMirAxes,'Units');		set(hMirAxes,'Units','points');
			Pos = get(hMirAxes,'Position');			set(hMirAxes,'Units',oldUnits);
			vscale = 1/Pos(4) * diff(axLims(1:2));	hscale = 1/Pos(3) * diff(axLims(3:4));
			%vscale = (vscale + hscale) / 2;			hscale = vscale;	% For not having a X|Y direction dependency
			DAR = get(hMirAxes, 'DataAspectRatio');
			if (DAR(1) == 1 && DAR(1) ~= DAR(2))	% To account for the "Scale geog images at mean lat" effect
				vscale = vscale * DAR(2);		hscale = hscale * DAR(1);
			end
			out = [out hscale vscale];

		case 'pline'
			ind1 = strfind(str,'[');
			if (numel(ind1) ~= 1),	msg = 'Wrong syntax. ''['' symbol must appear one and only one time';	return,	end
			ind2 = strfind(str,']');
			if (numel(ind2) ~= 1),	msg = 'Wrong syntax. '']'' symbol must appear one and only one time';	return,	end
			ind3 = strfind(str,';');
			if (numel(ind3) ~= 1),	msg = 'Wrong syntax. '';'' symbol must appear one and only one time';	return,	end
			[t, r] = strtok(str(ind1+1:ind3-1));	% Count the number of X coords
			if (isempty(t)),		msg = 'There are no X pts inside the brackets. You fill clever?';		return,	end 
			n = 1;
			while (~isempty(r))
				[t, r] = strtok(r);		n = n + 1;
			end
			out = {str(ind1+1:ind3-1); str(ind3+1:ind2-1); n};

		case 'toRidge'
			if (~np.validGrid)
				msg = 'This operation is only possible with grids. Not images';		return
			end
			N = strtok(str);
			if (isempty(N)),	N = '5';		end		% Default value
			out = round(abs( str2double(N) ));
			if (isnan(out))
				if (N(1) == 'N')
					msg = 'The argument "N" is no to be taken literaly. It must contain the number of grid cell arround vertex';
				else
					msg = 'the N argument is nonsense';
				end
			end

		case 'stitch'
			out.toDegFac = 0;		% Means no conversion from cartesian to degrees
			if ( isempty(str) || strcmpi(str, 'TOL') )
				out.val = Inf;
			else
				if ( strcmpi(str(end), 'm') || strcmpi(str(end), 'k') || strcmpi(str(end), 'n') )	% Convert any of those to degrees
					convFrom = lower(str(end));		str(end) = [];
					% Right. We must convert the DIST distance given in cartesians to degrees
					% But the problem is that here we don't know the lat so we can't finish the conversion. So, send back the collected info
					if (convFrom == 'm'),		out.toDegFac = 1;
					elseif (convFrom == 'k'),	out.toDegFac = 1000;
					else						out.toDegFac = 1852;
					end
				end
				out.val = abs(str2double(str));
			end
	end

% --------------------------------------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function h = line_operations_LayoutFcn(h1, hMirFig, axPos, y_off)
	if (nargin == 4)
		ofset = [axPos(1) y_off 0 0];
		h2 = hMirFig;
		floating = false;
	else
		ofset = [0 0 0 0];
		h2 = h1;
		floating = true;
	end
	
set(h1, 'Position',[520 729 470 48],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'Name','Line operations',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','Call',...
'Tag','figure1');

W = 465;	xOff = 3;
if (~floating)
	Wf = min(axPos(3),W+60);		% Frame width
	xOff = (axPos(3) - Wf) / 2;
	h = zeros(1,6);
	h(6) = uicontrol('Parent',h2, 'Position',[xOff 2 Wf axPos(2)] + ofset,'Style', 'frame');
	xOff = xOff + (Wf - W) / 2;		% X offset for the buttons so they are centered in the frame
end

h(1) = uicontrol('Parent',h2, 'Position',[xOff 2 W 24] + ofset,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_cmd');

h(2) = uicontrol('Parent',h2, 'Position',[xOff 26 91 21] + ofset,...
'Call',{@lop_CB, h1},...
'FontName','Helvetica', 'FontSize',9,...
'String','Get line(s)',...
'TooltipString','Have 0 lines to play with',...
'Tag','push_pickLine');

h(3) = uicontrol('Parent',h2, 'Position',[140+xOff 26 201 21] + ofset,...
'BackgroundColor',[1 1 1],...
'Call',{@lop_CB, h1},...
'FontName','Helvetica',...
'Style','popupmenu',...
'TooltipString','See list of possible operations and slect template if wished',...
'Value',1,...
'Tag','popup_cmds');

h(4) = uicontrol('Parent',h2, 'Position',[375+xOff 26 91 21] + ofset,...
'Call',{@lop_CB, h1},...
'FontName','Helvetica', 'FontSize',9,...
'String','Apply',...
'Tag','push_apply');

h(5) = uicontrol('Parent',h2, 'Position',[101+xOff 30 14 14] + ofset,...
'BackgroundColor',[1 0 0],...
'Enable','inactive',...
'FontName','Helvetica', 'FontSize',9,...
'TooltipString','Color informs if we have line(s) to work with.',...
'Tag','push_semaforo');

function lop_CB(hObject, eventdata, h1)
% This function is executed by the Call and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(h1));

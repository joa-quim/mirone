function varargout = line_operations(varargin)
% Wraper figure to perform vectorial operations on line/patch objects
%
%	Copyright (c) 2004-2008 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------
 
	hObject = figure('Tag','figure1','Visible','off');
	line_operations_LayoutFcn(hObject);
	handles = guihandles(hObject);
	
	if (~isempty(varargin))
 	   handles.hMirFig = varargin{1}.figure1;
	   handles.hMirAxes = varargin{1}.axes1;
	   handles.lt = varargin{1}.DefLineThick;
	   handles.lc = varargin{1}.DefLineColor;
	   IamCompiled = varargin{1}.IamCompiled;
	   move2side(handles.hMirFig, hObject, 'bottom')
   end

	handles.known_ops = {'buffer'; 'polysimplify'; 'bspline'; 'cspline'; 'polyunion'; 'polyintersect'; ...
			'polyxor'; 'polyminus'; 'line2patch'; 'pline'; 'hand2Workspace'};
	handles.hLine = [];
	set(handles.popup_cmds,'Tooltip', 'Select one operation from this list')

	handles.ttips = cell(numel(handles.known_ops));
	handles.ttips{1} = 'Select one operation from this list';
	handles.ttips{2} = sprintf(['Compute buffer zones arround polylines.\n' ...
								'Replace DIST by the width of the buffer zone.']);
	handles.ttips{3} = sprintf(['Approximates polygonal curve with desired precision\n' ...
								'using the Douglas-Peucker algorithm.\n' ...
								'Replace TOL by the desired approximation accuracy.\n' ...
								'When data is in geogs, TOL is the tolerance in km.']);
	handles.ttips{4} = sprintf(['Smooth line with a B-form spline\n' ...
								'It will open a help control window to\n' ...
								'help with selection of nice parameters']);
	handles.ttips{5} = sprintf(['Smooth line with a cardinal spline\n' ...
								'Replace N by the downsampling rate. For example a N of 10\n' ...
								'will take one every other 10 vertices of the polyline.\n' ...
								'The second argument, RES represents the number subdivisions\n' ...
								'between consecutive data points. As an example RES = 10 will\n' ...
								'split the downsampled interval in 10 sub-intervald, thus reseting\n' ...
								'the original number of point, excetp at the end of the line.']);
	handles.ttips{6} = 'Performs the boolean operation of Union to the selected polygons.';
	handles.ttips{7} = 'Performs the boolean operation of Intersection to the selected polygons.';
	handles.ttips{8} = 'Performs the boolean operation of exclusive OR to the selected polygons.';
	handles.ttips{9} = 'Performs the boolean operation of subtraction to the selected polygons.';
	handles.ttips{10} = 'Convert line objects into patch. Patches, for example, accept fill color.';
	handles.ttips{11} = sprintf(['Dray a polyline with vertices defined by coords [x1 xn; y1 yn].\n' ...
								'Note: you must use the brackets and semi-comma notation as above.\n' ...
								'Example vector: [1 1.5 3.1; 2 4 8.4]']);
	handles.ttips{12} = sprintf(['Send the selected object handles to the Matlab workspace.\n' ...
								'Use this if you want to gain access to all handle properties.']);

	if (IamCompiled),	handles.known_ops(11) = [];	handles.ttips(12);	end		% regretably

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hMirFig,'dependentFigs',plugedWin);

	set(hObject,'Visible','on');
	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end

% --------------------------------------------------------------------------------------------------
function edit_cmd_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------------------------------------
function push_pickLine_Callback(hObject, eventdata, handles)    
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
    figure(handles.figure1)                 % Bring this figure to front again
	guidata(hObject, handles);

% --------------------------------------------------------------------------------------------------
function popup_cmds_Callback(hObject, eventdata, handles)
	str = get(hObject,'String');	val = get(hObject,'Val');
	set(hObject, 'Tooltip', handles.ttips{val})
	if (val == 1)
		set(handles.edit_cmd,'String', '')
	else
		set(handles.edit_cmd,'String', str{val})
	end

% --------------------------------------------------------------------------------------------------
function push_apply_Callback(hObject, eventdata, handles)
	cmd = get(handles.edit_cmd,'String');
	if (isempty(cmd))
		errordlg('Fiu Fiu!! Apply WHAT????','ERROR'),	return
	end

	[t, r] = strtok(cmd);
	ind = strmatch(t, handles.known_ops);
	if (isempty(ind))
		% Here comes a function call which trys to adress the general command issue
		return
	end
	if (isempty(handles.hLine) && ~strcmp(handles.known_ops{ind}, 'pline'))
		errordlg('Fiu Fiu!! Apply WHERE????','ERROR'),	return
	end
	
	switch handles.known_ops{ind}
		case 'buffer'
			[dist, msg] = validate_args(handles.known_ops{ind}, r);
			if (~isempty(msg)),		errordlg(msg,'ERROR'),		return,		end
			for (k = 1:numel(handles.hLine))
				x = get(handles.hLine(k), 'xdata');		y = get(handles.hLine(k), 'ydata');
 				[y, x] = buffer_j(y, x, dist, 'out');
				h = line('XData',x, 'YData',y, 'Parent',handles.hMirAxes, 'Color',handles.lc, 'LineWidth',handles.lt, 'Tag','polyline');
				draw_funs(h,'line_uicontext')
			end
		case 'polysimplify'
			[tol, msg] = validate_args(handles.known_ops{ind}, r);
			if (~isempty(msg)),		errordlg(msg,'ERROR'),		return,		end
			for (k = 1:numel(handles.hLine))
				x = get(handles.hLine(k), 'xdata');		y = get(handles.hLine(k), 'ydata');
 				B = cvlib_mex('dp', [x(:) y(:)], tol, 'GEOG');
				h = line('XData',B(:,1), 'YData',B(:,2), 'Parent',handles.hMirAxes, 'Color',handles.lc, 'LineWidth',handles.lt, 'Tag','polyline');
				draw_funs(h,'line_uicontext')
			end
		case 'bspline'
			x = get(handles.hLine(1), 'xdata');		y = get(handles.hLine(1), 'ydata');
			[pp,p] = spl_fun('csaps',x,y);			% This is just to get csaps's p estimate
			y = spl_fun('csaps',x,y,p,x);
			h = line('XData', x, 'YData', y, 'Parent',handles.hMirAxes, 'Color',handles.lc, 'LineWidth',handles.lt, 'Tag','polyline');
			draw_funs(h,'line_uicontext')
			smoothing_param(p, [x(1) x(2)-x(1) x(end)], handles.hMirFig, handles.hMirAxes, handles.hLine, h);
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
		case 'line2patch'
			for (k = 1:numel(handles.hLine))
				if (~strcmp(get(handles.hLine(k),'type'), 'line')),	continue,	end
				x = get(handles.hLine(k), 'xdata');		y = get(handles.hLine(k), 'ydata');
				h = patch('XData',x, 'YData',y, 'Parent',handles.hMirAxes, 'EdgeColor',get(handles.hLine(k),'color'), ...
					'FaceColor','none', 'LineWidth',get(handles.hLine(k),'LineWidth'), 'Tag','line2patch');
				draw_funs(h,'line_uicontext')
				delete(handles.hLine(k))
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
				errordlg(['Something screwd up. Error message is ' lasterr])
			end
		case 'hand2Workspace'
			assignin('base','lineHandles',handles.hLine);
	end
	
% --------------------------------------------------------------------------------------------------
function [out, msg] = validate_args(qual, str, np)
% Check for errors on arguments to operation QUAL and return required args, or error
	msg = [];
	switch qual
		case 'buffer'
			d = strtok(str);
			if (strcmp(d, 'DIST'))
				msg = 'The argument "DIST" must be replaced by a numeric value representing the width of the buffer zone';	return
			end
			out = abs( str2double(d) );
			if (isnan(out))
				msg = 'BUFFER argument is nonsense';
			end
		case 'polysimplify'
			d = strtok(str);
			if (strcmp(d, 'TOL'))
				msg = 'The argument "TOL" must be replaced by a numeric value representing the tolerance';	return
			end
			out = abs( str2double(d) );
			if (isnan(out))
				msg = 'POLYSIMPLIFY argument is nonsense';
			end
		case 'cspline'
			[N, r] = strtok(str);
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
		case {'polyunion' 'polyintersect' 'polyxor' 'polyminus'}
			out = [];
			if (str < 2)		% str is in fact the number of handles
				msg = 'you want to make the union of a single line???. Wierd!';
			end
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
	end

% --------------------------------------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function line_operations_LayoutFcn(h1)
set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Line operations',...
'NumberTitle','off',...
'Position',[520 729 470 60],...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@line_operations_uicallback,h1,'edit_cmd_Callback'},...
'HorizontalAlignment','left',...
'Position',[5 5 461 25],...
'Style','edit',...
'Tag','edit_cmd');

uicontrol('Parent',h1,...
'Callback',{@line_operations_uicallback,h1,'push_pickLine_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[5 33 91 23],...
'String','Get line(s)',...
'TooltipString','Have 0 lines to play with',...
'Tag','push_pickLine');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@line_operations_uicallback,h1,'popup_cmds_Callback'},...
'FontName','Helvetica',...
'Position',[142 34 201 22],...
'String',{'Possible commands'; 'buffer DIST'; 'polysimplify TOL'; 'bspline'; 'cspline N RES'; 'polyunion'; 'polyintersect'; 'polyxor'; 'polyminus'; 'line2patch'; 'pline [x1 ..xn; y1 .. yn]'; 'hand2Workspace' },...
'Style','popupmenu',...
'TooltipString','See list of possible operations and slect template if wished',...
'Value',1,...
'Tag','popup_cmds');

uicontrol('Parent',h1,...
'Callback',{@line_operations_uicallback,h1,'push_apply_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[375 33 91 23],...
'String','Apply',...
'Tag','push_apply');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0 0],...
'Enable','inactive',...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[103 38 12 12],...
'TooltipString','Color informs if we have line(s) to work with.',...
'Tag','push_semaforo');

function line_operations_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

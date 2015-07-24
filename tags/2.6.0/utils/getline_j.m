function varargout = getline_j(varargin)
%GETLINE Select polyline with mouse.
%   [X,Y] = GETLINE_J(FIG) lets you select a polyline in the current axes of figure FIG
%   using the mouse.  Coordinates of the polyline are returned in X and Y.  Use normal
%   button clicks to add points to the polyline.  A shift-, or double-click adds a
%   final point and ends the polyline selection. Pressing RETURN or ENTER ends the polyline
%   selection without adding a final point. Pressing DELETE removes the previously
%   selected point from the polyline.
%
%   [X,Y] = GETLINE_J(AX) lets you select a polyline in the axes specified by the handle AX.
%   [X,Y] = GETLINE_J is the same as [X,Y] = GETLINE_J(GCF).
%   [X,Y] = GETLINE_J(HLINE) extends (continue adding points) the editing of the line whose handle is HLINE.
%   [X,Y] = GETLINE_J(...,'closed') animates and returns a closed polygon.
%   [X,Y] = GETLINE_J(...,'freehand') draw a line following the mouse movements.
%   [X,Y] = GETLINE_J(...,'dynamic') calls grdtrack_m to do dynamic profiling.
%   [X,Y] = GETLINE_J(...,'spline') draw an spline interpolated line.
%

%   Grandfathered syntaxes:
%   XY = GETLINE_J(...) returns output as M-by-2 array; first column is X; second column is Y.

%   Havily hacked version of getline that, contrary to the original, let be compiled.
%   Also the right-click button was reprogramed. Now it removes the previously selected point.

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

% $Id: $

	ud.GETLINE_ISCLOSED = false;	ud.GETLINE_FREEHAND = false;
	ud.GETLINE_DYNAMIC  = false;	ud.GETLINE_SPLINE   = false;
	extend_line = false;
	if ((nargin >= 1) && (ischar(varargin{end})))
		str = varargin{end};
		if (str(1) == 'c')					% getline_j(..., 'closed')
			ud.GETLINE_ISCLOSED = true;		varargin = varargin(1:end-1);
		elseif (str(1) == 'f')				% getline_j(..., 'freehand')
			ud.GETLINE_FREEHAND = true;		varargin = varargin(1:end-1);
		elseif (str(1) == 'd')				% getline_j(..., 'dynamic')
			ud.GETLINE_DYNAMIC = true;		varargin = varargin(1:end-1);
		elseif (str(1) == 's')				% getline_j(..., 'spline')
			ud.GETLINE_SPLINE = true;		varargin = varargin(1:end-1);
		end
	end

	ud.GETLINE_X = [];     ud.GETLINE_Y = [];

	if (length(varargin) < 1)
		ud.GETLINE_AX = gca;
		ud.GETLINE_FIG = get(ud.GETLINE_AX, 'Parent');
		while ( ~strcmp('figure', get(ud.GETLINE_FIG,'type')) )	% In case the axes is a uipanel descendent
			ud.GETLINE_FIG = get(ud.GETLINE_FIG,'parent');
		end
	else
		if (~ishandle(varargin{1}))
			error('First argument is not a valid handle');
		end
		switch get(varargin{1}, 'Type')
			case 'figure'
				ud.GETLINE_FIG = varargin{1};
				ud.GETLINE_AX = get(ud.GETLINE_FIG, 'CurrentAxes');
				if (isempty(ud.GETLINE_AX))
					ud.GETLINE_AX = axes('Parent', ud.GETLINE_FIG);
				end
			case 'axes'
				ud.GETLINE_AX = varargin{1};
				ud.GETLINE_FIG = get(ud.GETLINE_AX, 'Parent');
				while ( ~strcmp('figure', get(ud.GETLINE_FIG,'type')) )	% In case the axes is a uipanel descendent
					ud.GETLINE_FIG = get(ud.GETLINE_FIG,'parent');
				end
			case 'line'
				ud.GETLINE_AX = get(varargin{1}, 'Parent');
				ud.GETLINE_FIG = get(ud.GETLINE_AX, 'Parent');
				while ( ~strcmp('figure', get(ud.GETLINE_FIG,'type')) )	% In case the axes is a uipanel descendent
					ud.GETLINE_FIG = get(ud.GETLINE_FIG,'parent');
				end
				ud.GETLINE_X = get(varargin{1}, 'XData');
				ud.GETLINE_Y = get(varargin{1}, 'YData');
				set(varargin{1}, 'XData',[], 'YData',[])		% Don't need those duplicated
				extend_line = true;
			otherwise
				error('First argument should be a figure, axes or line handle');
		end
	end

	xlimorigmode = get(ud.GETLINE_AX,'xlimmode');     ylimorigmode = get(ud.GETLINE_AX,'ylimmode');
	set(ud.GETLINE_AX,'xlimmode','manual');           set(ud.GETLINE_AX,'ylimmode','manual');

	% Remember initial figure state
	old_db = get(ud.GETLINE_FIG, 'DoubleBuffer');
	state = uisuspend_j(ud.GETLINE_FIG);

	% Set up initial callbacks for initial stage
	set(ud.GETLINE_FIG, 'Pointer', 'crosshair', ...
		'WindowButtonDownFcn', {@FirstButtonDown,ud.GETLINE_FIG},...
		'KeyPressFcn', {@KeyPress,ud.GETLINE_FIG}, 'DoubleBuffer', 'on');

	figure(ud.GETLINE_FIG);		% Bring target figure forward

	% Initialize the lines to be used for the drag
	ud.GETLINE_H1 = line('Parent', ud.GETLINE_AX, ...
						'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y, ...
						'Visible', 'off', 'Clipping', 'off', ...
						'Color', 'k', 'LineStyle', '-');

	ud.GETLINE_H2 = line('Parent', ud.GETLINE_AX, ...
						'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y, ...
						'Visible', 'off', 'Clipping', 'off', ...
						'Color', 'w', 'LineStyle', ':');

	setappdata(ud.GETLINE_FIG, 'fromGL', ud);

	if (extend_line)
		pt = get(ud.GETLINE_AX, 'CurrentPoint');		% Add a new point to show right away that we're in extending mode
		newx = pt(1,1);    newy = pt(1,2);
		x = [ud.GETLINE_X newx];    y = [ud.GETLINE_Y newy];
		set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', x, 'YData', y, 'Vis', 'on')
		set(ud.GETLINE_FIG, 'WindowButtonMotionFcn',{@ButtonMotion,ud.GETLINE_FIG})
	end

	% We're ready; wait for the user to do the drag. Wrap the call to waitfor
	% in try-catch so we'll have a chance to clean up after ourselves.
	errCatch = 0;
	try
		waitfor(ud.GETLINE_H1, 'UserData', 'Completed');
	catch
		errCatch = 1;
	end

	% After the waitfor, if GETLINE_H1 is still valid and its UserData is 'Completed', then the user
	% completed the drag.  If not, the user interrupted the action somehow, perhaps by a Ctrl-C in the
	% command window or by closing the figure.

	if (errCatch == 1)
		errStatus = 'trap';
	elseif (~ishandle(ud.GETLINE_H1) || ~strcmp(get(ud.GETLINE_H1, 'UserData'), 'Completed'))
		errStatus = 'unknown';
	else
		errStatus = 'ok';
		x = get(ud.GETLINE_H1,'XData');
		y = get(ud.GETLINE_H1,'YData');
		% If no points were selected, return rectangular empties.
		% This makes it easier to handle degenerate cases in functions that call getline_j.
		if (isempty(x));        x = zeros(0,1);     y = zeros(0,1);   end
	end

	% Delete the animation objects
	if (ishandle(ud.GETLINE_H1));    delete(ud.GETLINE_H1);     end
	if (ishandle(ud.GETLINE_H2));    delete(ud.GETLINE_H2);     end

	% Restore the figure's initial state
	if (ishandle(ud.GETLINE_FIG))
	   uirestore_j(state, 'nochildren');
	   set(ud.GETLINE_FIG, 'DoubleBuffer', old_db);
	end

	set(ud.GETLINE_AX,'xlimmode',xlimorigmode);   set(ud.GETLINE_AX,'ylimmode',ylimorigmode);

	% Depending on the error status, return the answer or generate an error message.
	switch errStatus
	case 'ok'      % Return the answer
		if (nargout >= 2)
			varargout{1} = x;        varargout{2} = y;
		else        % Grandfathered output syntax
			varargout{1} = [x(:) y(:)];
		end
	case 'trap'    % An error was trapped during the waitfor
		%error('Interruption during mouse selection.');
		varargout{1} = [];        varargout{2} = [];
	case 'unknown'
		% User did something to cause the polyline selection to
		% terminate abnormally.  For example, we would get here
		% if the user closed the figure in the middle of the selection.
		%error('Interruption during mouse selection.');
		varargout{1} = [];        varargout{2} = [];
	end

	try  rmappdata(ud.GETLINE_FIG,'fromGL');     end

%-------------------------------------------------------------------------------
function KeyPress(obj,eventdata,hfig)
	ud = getappdata(hfig, 'fromGL');

	if (ud.GETLINE_FREEHAND)        % NextButtonDown is the one who update those
		ud.GETLINE_X = get(ud.GETLINE_H1,'XData');  % but for the freehand case
		ud.GETLINE_Y = get(ud.GETLINE_H1,'YData');  % it didn't do it.
	end
	key = get(hfig, 'CurrentCharacter');
	switch key
	case {char(8), char(127)}  % delete and backspace keys
		% remove the previously selected point
		switch length(ud.GETLINE_X)
		case 0        % nothing to do
		case 1
			ud.GETLINE_X = [];        ud.GETLINE_Y = [];
			% remove point and start over
			set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
			set(hfig, 'WindowButtonDownFcn', {@FirstButtonDown,hfig}, 'WindowButtonMotionFcn', '');
		otherwise        % remove last point
			if (ud.GETLINE_ISCLOSED)
				ud.GETLINE_X(end-1) = [];  ud.GETLINE_Y(end-1) = [];
			else
				ud.GETLINE_X(end) = [];    ud.GETLINE_Y(end) = [];
			end
			set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
		end
	case 'c'              % "c" char key (close line)
		if (length(ud.GETLINE_X) > 2  && ~ud.GETLINE_ISCLOSED)  % don't close a line with less than 2 vertex
			ud.GETLINE_X = [ud.GETLINE_X ud.GETLINE_X(1)];
			ud.GETLINE_Y = [ud.GETLINE_Y ud.GETLINE_Y(1)];
			set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
			set(ud.GETLINE_H1, 'UserData', 'Completed');
		end
	case {char(13), char(3)}   % enter and return keys
		% return control to line after waitfor
		set(ud.GETLINE_H1, 'UserData', 'Completed');     
	end
	setappdata(hfig, 'fromGL', ud);

%----------------------------------------------------------------------------------
function FirstButtonDown(obj,eventdata,hfig)
	ud = getappdata(hfig, 'fromGL');
	pt = get(ud.GETLINE_AX, 'CurrentPoint');
	x = pt(1,1);    y = pt(1,2);

	% check if GETLINE_X,GETLINE_Y is inside of axis
	x_lim = get(ud.GETLINE_AX,'xlim');      y_lim = get(ud.GETLINE_AX,'ylim');
	if (x>=x_lim(1)) && (x<=x_lim(2)) && (y>=y_lim(1)) && (y<=y_lim(2))    % inside axis limits
		ud.GETLINE_X = [ud.GETLINE_X x];    ud.GETLINE_Y = [ud.GETLINE_Y y];
	else    % outside axis limits, ignore this FirstButtonDown
		return
	end

	if (ud.GETLINE_ISCLOSED)
		ud.GETLINE_X = [ud.GETLINE_X ud.GETLINE_X];    ud.GETLINE_Y = [ud.GETLINE_Y ud.GETLINE_Y];
	end

	set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y, 'Visible', 'on');

	if (~strcmp(get(hfig, 'SelectionType'), 'normal'))    % We're done!
		set(ud.GETLINE_H1, 'UserData', 'Completed');
	else    % Let the motion functions take over.
		% We must first reset WindowButtonDownFcn otherwise it will bug when compiled
		set(hfig, 'WindowButtonDownFcn', '', 'WindowButtonMotionFcn',{@ButtonMotion,hfig});
		set(hfig, 'WindowButtonDownFcn', {@NextButtonDown,hfig})
	end
	setappdata(hfig, 'fromGL', ud);

%---------------------------------------------------------------------------------------
function NextButtonDown(obj,eventdata,hfig)
	ud = getappdata(hfig, 'fromGL');

	selectionType = get(ud.GETLINE_FIG, 'SelectionType');
	if (~strcmp(selectionType,'open') && strcmp(selectionType,'normal') && ~ud.GETLINE_FREEHAND)
		% We don't want to add a point on the second click of a double-click
		pt = get(ud.GETLINE_AX, 'CurrentPoint');
		x = pt(1,1);    y = pt(1,2);

		% check if GETLINE_X,GETLINE_Y is inside of axis
		x_lim = get(ud.GETLINE_AX,'xlim');      y_lim = get(ud.GETLINE_AX,'ylim');
		if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2))    % outside axis limits, ignore this ButtonDown
			return
		end

		% Pad mouses are terrible in creating (nearly) duplicate points that are not detected as double-clicks
		% Here we try to detect such events and ignore them. The idea is to detect new vertex that are close to
		% last vertex. The criteria used is to find if they fall inside the red-marker (the red rectangle of the
		% edit mode) of the previous point. Because the computation depends on the zoom level we have to do it here.
		axUnit = get(ud.GETLINE_AX,'Units');	set(ud.GETLINE_AX,'Units','pixels')
		pos = get(ud.GETLINE_AX,'Position');	set(ud.GETLINE_AX,'Units',axUnit)	
		scaleX = get(0,'ScreenPixelsPerInch') * (6/72) / pos(3);
		scaleY = get(0,'ScreenPixelsPerInch') * (6/72) / pos(4);
		dx = 0.6 * diff(x_lim) * scaleX;		dy = 0.6 * diff(y_lim) * scaleY;
		x0 = x - dx;		x1 = x + dx;
		y0 = y - dy;		y1 = y + dy;
		if (ud.GETLINE_X(end) > x0 && ud.GETLINE_X(end) < x1 && ud.GETLINE_Y(end) > y0 && ud.GETLINE_Y(end) < y1)
			return
		end

		if (ud.GETLINE_ISCLOSED)
			ud.GETLINE_X = [ud.GETLINE_X(1:end-1) x ud.GETLINE_X(end)];
			ud.GETLINE_Y = [ud.GETLINE_Y(1:end-1) y ud.GETLINE_Y(end)];
		else
			ud.GETLINE_X = [ud.GETLINE_X x];        ud.GETLINE_Y = [ud.GETLINE_Y y];
			% At this point, if (ud.GETLINE_SPLINE) we'll see the knots polyline until the next mouse mov
		end
		set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
	end

	if (strcmp(get(hfig,'SelectionType'),'alt') && ~ud.GETLINE_FREEHAND)	% Right-click, delete previous point
		pt = get(ud.GETLINE_AX, 'CurrentPoint');
		x = pt(1,1);    y = pt(1,2);
		% check if GETLINE_X,GETLINE_Y is inside of axis
		x_lim = get(ud.GETLINE_AX,'xlim');      y_lim = get(ud.GETLINE_AX,'ylim');
		if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2))    % outside axis limits, ignore this ButtonDown
			return
		end    
		switch numel(ud.GETLINE_X)
			case 0        % nothing to do
			case 1
				ud.GETLINE_X = [];        ud.GETLINE_Y = [];
				% remove point and start over
				set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
				set(hfig, 'WindowButtonDownFcn', {@FirstButtonDown,hfig}, 'WindowButtonMotionFcn', '');
			otherwise        % remove last point
				if (ud.GETLINE_ISCLOSED)
					ud.GETLINE_X(end-1) = [];  ud.GETLINE_Y(end-1) = [];
				else
					ud.GETLINE_X(end) = [];    ud.GETLINE_Y(end) = [];
				end
				set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
		end
	end

	if (strcmp(selectionType, 'extend') || strcmp(selectionType, 'open'))    % We're done (midle button or double-click)
		if (ud.GETLINE_SPLINE)		% Get the final splined curve
			[ud.GETLINE_X, ud.GETLINE_Y] = spline_interp(ud.GETLINE_X, ud.GETLINE_Y);
			set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
		end
		set(ud.GETLINE_H1, 'UserData', 'Completed');
	end
	setappdata(ud.GETLINE_FIG, 'fromGL', ud);

%-----------------------------------------------------------------------------------
function ButtonMotion(obj,eventdata,hFig)
	ud = getappdata(hFig, 'fromGL');

	pt = get(ud.GETLINE_AX, 'CurrentPoint');
	newx = pt(1,1);    newy = pt(1,2);
	if (ud.GETLINE_ISCLOSED && (length(ud.GETLINE_X) >= 3))
		x = [ud.GETLINE_X(1:end-1) newx ud.GETLINE_X(end)];
		y = [ud.GETLINE_Y(1:end-1) newy ud.GETLINE_Y(end)];
	else
		x = [ud.GETLINE_X newx];    y = [ud.GETLINE_Y newy];
	end

	if (~ud.GETLINE_FREEHAND)
		if (ud.GETLINE_SPLINE)
			[x,y] = spline_interp(x,y);
		end
		set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', x, 'YData', y);
		if (ud.GETLINE_DYNAMIC)
			% First logical indicates "not a point interpolation" and second that we are in dynamic mode
			grid_profiler(hFig, x, y, false, true);
		end
	else            % Do a freehand drawing
		set([ud.GETLINE_H1 ud.GETLINE_H2],'xdata',[get(ud.GETLINE_H1,'xdata'),pt(1,1)],...
				'ydata',[get(ud.GETLINE_H1,'ydata'),pt(1,2)]);
	end

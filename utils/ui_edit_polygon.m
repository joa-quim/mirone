function ui_edit_polygon(varargin)
% Interactive edit polylines (closed or open) and patches.
%
% ui_edit_polygon(handle1, handle2, ..., 'move_choice')
% ui_edit_polygon([handle1, handle2, ...],'move_choice')  
%       (..., where handle_i is a handle to a line or patch)
%
% Interactive edit polylines (closed or open) and patches. If obj(handle_i) has allready
% been made editable with ui_edit_polygon, the function will check the edit state
% of obj and turn it off if necessary
%
% This function may also call a callback function registered in 'RunCB' appdata
% that will be executed at a button up. See wbu_EditPolygon() to see how it works.
% Another option is a function, registered in 'BD_runCB' appdata, that will be
% called only when doing a Shift+click on the line element.
%
% 	MOVE_CHOICE		% Controls what mouse selection is used to move the whole polygon.
% 					% If == empty, polygon is moved with a left click. That has often anoying side effects
% 					% If == 'extend', means that we need a Shift-click left mouse button or click both
% 					%		left and right mouse buttons to move it. A bit more cumbersome, but safer.
% 					% If == 'y', OR == 'x' means that Y or X vertices are move while editing
%
%	An alternative way to set the MOVE_CHOICE option is to put it as an appdata of the axes where lines are drawn.
%	As an example of how it works see this snipet showing how it is recovered here
% 			hAx = get(varargin{1},'parent');				% Get parent axes
% 			move_choice = getappdata(hAx, 'MovPolyg');		% Get MOVE_CHOICE option
% 
% -----------------------------------------------------------------------------
% Double clicking on the line or patch displays its vertices as control points.
% To stop edit double click again. 
% To edit vertex, click and drag control point.
% Click and drag outside the line/patch vertices control points moves the object (but see 'MOVE_CHOICE')
%
%  In edit mode, hit the following keybord keys ...
%     "r", or backspace:			remove the active vertex
%     "i":							inserts a new vertex after the active vertex, or in case there is
%									no active vertex, guess the right order where to insert the new
%									vertex in such a way as not to crete zig-zags.
%     "b":							breaks a line in two parts. One from the first point up to the
%									active vertice and the other from the active vertice + 1 till the end.
%									Patches cannot be broken.
%     "c":							close a polyline
%     "e":							put the line in extending mode (add points at the end) by calling getline_j
%     "n":							Move active vertex one step twards end of line
%     "p":							Move active vertex one step twards beguining of line
%     "P":							close a polyline and convert it into a patch
%     escape:						stop edit mode
%     delete:						delete the line/patch object
%
%  When VARARGIN is a handle to a rectangle (line or patch) three things may happen:
%   1. The rectangle is described in counter-clockwise way with origin at lower left point
%       - The editing preserves the rectangular shape
%   2. The rectangle is described in a clockwise way, with origin at lower left point
%       - The editing deforms the rectangular.
%   3. The remain cases are not forseen. Unknown behavior.
%
%  Author(s)
% ----------
%   Joaquim Luis (jluis@ualg.pt) - Original version
%       25-Oct-2005 Updated version with bug corrections and added the rectangle mode descrimination. 
%   Sebastian Hoelz (hoelz@geophysik.tu-berlin.de)
%       10-Nov-2005 Nice improvments and code cleaning
%   JL & SH joint version -- ??-Nov-2005
%	JL	20-Jul-2008 Guess best positions for the insert option
%					Allow choosing what mouse selection is used to move the whole polygon.
%	JL	25-Aug-2008 Added the "e" (extend) keyboard option
%	JL	08-Feb-2010 Added 'y' or 'x' to the move_choice option
%	JL	21-apr-2012 Significant re-write with a large boost in efficiency by not needing to duplicate line

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

% $Id: ui_edit_polygon.m 7945 2016-09-05 22:51:05Z j $

	if (isa(varargin{end},'char'))	% Moving polygon option was transmitted in input
		move_choice = varargin{end};
		varargin(end) = [];
		if ( ~(strncmp(move_choice, 'ex' ,2) || strcmp(move_choice, 'x') || strcmp(move_choice, 'y')) )
			% Besides [], the other possible values are 'extend', 'x' or 'y'
			move_choice = [];
		end
	else							% See if the Moving polygon option is stored in axe's appdata
		try
			hAx = get(varargin{1},'parent');		% The TRY is used only if varargin{1} is not an handle
			move_choice = getappdata(hAx, 'MovPolyg');
		catch
			move_choice = [];
		end
	end
	if ( ~isempty(move_choice) && ~((move_choice(1) == 'x') || (move_choice(1) == 'y')) )
		move_choice = 'a';			% Means 'xy'
	end

	for (i = 1:numel(varargin))				% Argument check
		if (~ishandle(varargin{i}))
			disp(['Warning: Input argument ' sprintf('%d',i) 'is not a valid handle.'])

		elseif ~(strcmp(get(varargin{i},'Type'),'patch') || strcmp(get(varargin{i},'Type'),'line'))
			disp(['Warning: Input argument ' sprintf('%d',i) 'needs to be line or patch.'])

		% If the patch is allready editable by ui_edit_polygon, disable the edit mode (if set on) 
		elseif ~isempty(getappdata(varargin{i},'polygon_data'))
			s = getappdata(varargin{i},'polygon_data');
			if strcmpi(s.controls,'on')
				set(s.h_fig,'selectiontype','open');
				polygonui(s.h_pol)
			end

		else
			if (~strcmp(get(varargin{i},'Marker'),'none'))	% Lines with markers need a different set of functions
				s.duplicate = true;							% Because we can't simply set the markers ON to use for
			else											% edition. Instead we need to create a copy of the line object
				s.duplicate = false;
			end

			% Creating data-structure for polygon
			s.h_pol = varargin{i};
			s.h_vert = [];
			s.hCurrentMarker = [];
			s.h_ax = get(s.h_pol,'parent');
			s.h_fig = get(s.h_ax,'parent');
			s.KeyPress_orig = get(s.h_fig,'KeyPressFcn');

			s.controls = 'off';
			s.vert_index = [];
			s.save_x = [];
			s.save_y = [];

			x = get(s.h_pol,'XData');
			y = get(s.h_pol,'YData');
			s.is_closed_patch = false;
			s.is_patch  = strcmpi(get(s.h_pol,'Type'),'patch');
			s.is_rect = false;			s.keep_rect = false;
			if ( numel(x) == 5 && (x(1) == x(2)) && (x(3) == x(4)) && (y(1) == y(4)) && (y(2) == y(3)) )
				s.is_rect = true;		s.keep_rect = true;
			elseif ( numel(x) == 5 && (x(1) == x(4)) && (x(2) == x(3)) && (y(1) == y(2)) && (y(3) == y(4)) )
				s.is_rect = true;		s.keep_rect = false;
			end
			if ( numel(x) > 1 && (x(1) == x(end)) && (y(1) == y(end)))
				s.is_closed = true;
				if (s.is_patch),		s.is_closed_patch = true;  end
			else
				s.is_closed = false;
			end

			% This is for use in edit_polygon to detect which vertex (if any) was selected
			% The idea is that we are able to detect if the current point click was performed inside
			% the region, in map coordinates, delimited by the square marker. We can thus avoid
			% potentially very expensive operations based on detecting the vertex by finding the
			% closest vertice by minimum distances calculations.
			axUnit = get(s.h_ax,'Units');	set(s.h_ax,'Units','pixels')
			pos = get(s.h_ax,'Position');	set(s.h_ax,'Units',axUnit)	
			s.scaleX = get(0,'ScreenPixelsPerInch') * (6/72) / pos(3);
			s.scaleY = get(0,'ScreenPixelsPerInch') * (6/72) / pos(4);

			s.what_move = move_choice;

			set(s.h_pol,'buttondownfcn',@polygonui);
			setappdata(s.h_pol,'polygon_data',s)
		end
	end

%--------------------------------------------------
function polygonui(varargin)
% Set/unset red squares markers on polyline vertex that are used for edition

	if (~ishandle(varargin{1})),		return,		end
	s = getappdata(varargin{1},'polygon_data');
	stype = get(s.h_fig,'selectiontype');

	% See if we have a registered ButtonDownFcn function. If yes, run it and return
	if (strcmpi(stype, 'extend'))
		RunCB = getappdata(varargin{1}, 'BD_runCB');	% See if we have a callback function to run when just right-click
		if (~isempty(RunCB))
			if (numel(RunCB) == 1),		feval(RunCB{1});
			else						feval(RunCB{1},RunCB{2:end});
			end
		end
		return
	end

	if (~strcmpi(stype,'open')),	return,		end

	switch s.controls
		case 'on'						% Getting out of the edit mode
			s.vert_index = [];			% delete vertice markers
			if (~s.duplicate)
				set(s.h_pol, 'Marker', 'none')
			else
				delete(s.h_vert);		s.h_vert = [];		
			end
			if (~isempty(s.hCurrentMarker) && ishandle(s.hCurrentMarker))
				delete(s.hCurrentMarker);		s.hCurrentMarker = [];
			end
			s.controls = 'off';
			set(s.h_pol,'buttondownfcn',@polygonui);
			set(s.h_fig,'KeyPressFcn',s.KeyPress_orig)
			setappdata(s.h_pol,'polygon_data',s)
			setappdata(s.h_fig,'epActivHand',0)
			refresh(s.h_fig)		% Because of old R13 bugs

		case 'off'					% We are entering in the edit mode
			% Make sure that only one polygon is edited at a time
			h_active = getappdata(s.h_fig,'epActivHand');
			if (h_active),		polygonui(h_active),	end
			setappdata(s.h_fig,'epActivHand',s.h_pol)
			s.controls = 'on';
			if (~s.duplicate)
				set(s.h_pol, 'Marker','s','MarkerEdgeColor','r', 'MarkerFaceColor','none', ...
						'MarkerSize',5,'buttondownfcn',{@edit_polygon,s.h_pol});
			else
				s.h_vert = line('xdata',get(s.h_pol,'XData'),'ydata',get(s.h_pol,'YData'), ...
						'Parent',s.h_ax, 'Marker','s','color','r', 'MarkerFaceColor','none', ...
						'MarkerSize',5,'buttondownfcn',{@edit_polygon,s.h_pol});
				% Now we also need to set the line style to be equal to original so that we can drag it
				set(s.h_vert, 'linestyle',get(s.h_pol,'linestyle'), 'LineWidth',get(s.h_pol,'LineWidth'))
			end
			set(s.h_fig,'KeyPressFcn',{@KeyPress_local, s.h_pol})
			setappdata(s.h_pol,'polygon_data',s)
			refresh(s.h_fig)		% Because of old R13 bugs
	end

%--------------------------------------------------
function edit_polygon(obj,evt,h)
% Edit the polygon whose handle is h
	s = getappdata(h,'polygon_data');
	stype = get(s.h_fig,'selectiontype');
	if strcmp(stype,'open')					% When a line has many vertices, the markers may completely
		polygonui(s.h_pol,evt)				% hide it. So the only way to get out of edition mode is
		return								% to provide an other exit. That's where this call to
	end										% polygonui('markermousedown') comes to hand.        

	x_lim = get(s.h_ax,'xlim');
	y_lim = get(s.h_ax,'ylim');
	current_pt = get(s.h_ax, 'CurrentPoint');

	% Find out which vertice is beeing edited. Do that by finding if 'CurrentPoint' is
	% located inside the square marker region in   >>>> MAP COORDINATES <<<<
	x = get(s.h_pol,'XData');		y = get(s.h_pol,'YData');

	dx = 0.6 *  diff(x_lim) * s.scaleX;
	dy = 0.6 *  diff(y_lim) * s.scaleY;
	x0 = current_pt(1,1) - dx;		x1 = current_pt(1,1) + dx;
	y0 = current_pt(1,2) - dy;		y1 = current_pt(1,2) + dy;
	ind = ( x > x0 & x < x1 & y > y0 & y < y1);
	if (any(ind))
		id = find(ind);
		s.vert_index = id(1);	% Use the index for the unlikely case that more than one was found
	else
		% Click was not on the marker so see if it is a drag request
		move_polygon(s.h_pol)
		return
	end

	s.save_x = x(s.vert_index);			s.save_y = y(s.vert_index);		% needed in the "i"nsert option
	if isempty(s.hCurrentMarker)		% If the black marker doesn't exist, creat it
		s.hCurrentMarker = line('xdata',s.save_x,'ydata', s.save_y, 'parent', s.h_ax,'Marker','s', ...
								'MarkerFaceColor','k','MarkerSize',5,'Tag','current_marker');
		uistack_j(s.hCurrentMarker,'bottom')	% Since it has no ButtonDown and was on top
	else								% The black marker exists, just move it.
		set(s.hCurrentMarker,'XData',s.save_x,'YData',s.save_y)
	end
	setappdata(h,'polygon_data',s);
	state = uisuspend_safe(s.h_fig);		% Remember initial figure state

	set(s.h_fig,'WindowButtonMotionFcn',{@wbm_EditPolygon,h,[x_lim y_lim],s.h_fig},...
		'WindowButtonUpFcn',{@wbu_EditPolygon,h,state});

%--------------------------------------------------
function wbm_EditPolygon(obj,eventdata,h,lim,hFig)
	set(hFig, 'Pointer','fleur')		% I know, but this way fleur pointer shows only when we have a movement
	s = getappdata(h,'polygon_data');
	pt = get(s.h_ax, 'CurrentPoint');
	if (pt(1,1) < lim(1)) || (pt(1,1) > lim(2)) || (pt(1,2) < lim(3)) || (pt(1,2) > lim(4));   return; end
	xx = get(s.h_pol,'XData');      yy = get(s.h_pol,'YData');
	xx = xx(:)';                    yy = yy(:)';    % Make sure they are row vectors
	newx = pt(1,1);                 newy = pt(1,2);

	if (s.is_rect && s.keep_rect)				% We are dealing with a line/patch rectangle
		if (s.vert_index == 1)
			x = [newx newx xx(3) xx(4) newx];   y = [newy yy(2) yy(3) newy newy];
		elseif (s.vert_index == 2)
			x = [newx newx xx(3) xx(4) newx];   y = [yy(1) newy newy yy(4) yy(1)];
		elseif (s.vert_index == 3)
			x = [xx(1) xx(2) newx newx xx(1)];  y = [yy(1) newy newy yy(4) yy(1)];
		elseif (s.vert_index == 4)
			x = [xx(1) xx(2) newx newx xx(1)];  y = [newy yy(2) yy(3) newy newy];
		elseif (s.vert_index == 5)
			x = [newx newx xx(3) xx(4) newx];   y = [newy yy(2) yy(3) newy newy];
		end
	elseif (s.vert_index == 1)					% Selected first vertice
		if (s.is_closed && s.is_closed_patch)
			x = [newx xx(2:end-1) newx];    y = [newy yy(2:end-1) newy];
		elseif (s.is_closed && ~s.is_closed_patch)
			if (s.is_patch)
				x = [newx xx(2:end)];           y = [newy yy(2:end)];
			else								% deformable rectangle (rect given in cw direction)
				x = [newx xx(2:end-1) newx];    y = [newy yy(2:end-1) newy];
			end
		else									% "Open" polyline
			x = [newx xx(2:end)];           y = [newy yy(2:end)];
		end
	elseif (s.vert_index == length(xx) && ~s.is_patch)     % Selected last polyline vertice
		if (s.is_closed)
			x = [newx xx(2:end-1) newx];    y = [newy yy(2:end-1) newy];
		else									% "Open" polyline
			x = [xx(1:end-1) newx];         y = [yy(1:end-1) newy];
		end
	else										% Midle vertices
		x = [xx(1:s.vert_index-1) newx xx(s.vert_index+1:end)];
		y = [yy(1:s.vert_index-1) newy yy(s.vert_index+1:end)];
	end

	set(s.h_pol, 'XData',x, 'YData',y);
	set(s.hCurrentMarker,'XData',newx,'YData',newy)	% Move the current point marker together with the editing point
	if (s.duplicate)
		set(s.h_vert,'XData',x, 'YData',y)			% Move the marker together with the editing point
	end

%--------------------------------------------------
function wbu_EditPolygon(obj,evt,h,state)
	s = getappdata(h,'polygon_data');
	if (s.duplicate)
		uistack_j(s.h_vert,'top')		% Need to do this because the black marker has no ButtonDown and was on top
	end
	setappdata(h,'edited',true)
	uirestore_j(state, 'nochildren');	% Restore the figure's initial state
	RunCB = getappdata(h,'RunCB');		% See if we have a callback function to run after polygon edition
	if (~isempty(RunCB))
		if (numel(RunCB) == 1),		feval(RunCB{1})
		else						feval(RunCB{1},RunCB{2:end})
		end
	end

%--------------------------------------------------
function move_polygon(h)
% Move the polygon whose handle is h.
% Conditionally called by edit_polygon()

	s = getappdata(h,'polygon_data');
	stype = get(s.h_fig,'selectiontype');

	if (strcmp(stype,'open')),		polygonui(s.h_pol,[]),		return,		end
	if (~isempty(s.what_move) && ((stype(1) ~= 'e') && s.what_move(1) ~= 'y') )		% needs "password" to continue
		return
	end

	state = uisuspend_safe(s.h_fig);                 % Remember initial figure state    
	x_lim = get(s.h_ax,'xlim');        y_lim = get(s.h_ax,'ylim');
	current_pt = get(s.h_ax, 'CurrentPoint');
	setappdata(s.h_pol,'old_pt',[current_pt(1,1) current_pt(1,2)])

	set(s.h_fig,'WindowButtonMotionFcn',{@wbm_MovePolygon,h,[x_lim y_lim]},...
		'WindowButtonUpFcn',{@wbu_MovePolygon,h,state}, 'Pointer','fleur');

%--------------------------------------------------
function wbm_MovePolygon(obj,evt,h,lim)
	s = getappdata(h,'polygon_data');
	pt = get(s.h_ax, 'CurrentPoint');
	if (pt(1,1) < lim(1)) || (pt(1,1) > lim(2)) || (pt(1,2) < lim(3)) || (pt(1,2) > lim(4)),	return,	end
	old_pt = getappdata(s.h_pol,'old_pt');
	dx = pt(1,1) - old_pt(1);			dy = pt(1,2) - old_pt(2);
	xx = get(s.h_pol,'XData') + dx;		yy = get(s.h_pol,'YData') + dy;
	setappdata(s.h_pol,'old_pt',[pt(1,1) pt(1,2)])
 	if ( isempty(s.what_move) || s.what_move(1) == 'a' )
		set(s.h_pol, 'XData',xx, 'YData',yy);
 	elseif ( ~isempty(s.what_move) && s.what_move(1) == 'y' )
		set(s.h_pol, 'YData',yy);
	else
		set(s.h_pol, 'YData',xx);
	end

	if (~isempty(s.hCurrentMarker))				% If the black marker exists, move it too
		x = get(s.hCurrentMarker,'XData');		y = get(s.hCurrentMarker,'YData');
		if ( isempty(s.what_move) || s.what_move(1) == 'a' )
			if (isa(xx,'double'))
				x = x + dx;		y = y + dy;
			else
				x = cvlib_mex('addS', x, dx);	% R13 ...
				y = cvlib_mex('addS', y, dy);
			end
	 	elseif ( ~isempty(s.what_move) && s.what_move(1) == 'y' )	% gmtedit, I think
			y = y + dy;
		else
			x = x + dx;
		end
		set(s.hCurrentMarker,'XData',x,'YData',y)
	end

	% When we have a duplicated line (originally a line with Markers), we have to move it too
	if (s.duplicate)
		if ( isempty(s.what_move) || s.what_move(1) == 'a' )
			set(s.h_vert, 'XData',xx, 'YData',yy);
		elseif ( ~isempty(s.what_move) && s.what_move(1) == 'y' )
			set(s.h_vert, 'YData',yy);
		else
			set(s.h_vert, 'YData',xx);
		end
	end

%--------------------------------------------------
function wbu_MovePolygon(obj,evt,h,state)
	s = getappdata(h,'polygon_data');
	if (s.duplicate)
		uistack_j(s.h_vert,'top')			% Need to do this because the black marker has no ButtonDown and was on top
	end
	setappdata(h,'edited',true)
	uirestore_j(state, 'nochildren');		% Restore the figure's initial state
	RunCB = getappdata(h,'RunCB');			% See if we have a callback function to run after polygon edition
	if (~isempty(RunCB))
		if (numel(RunCB) == 1),		feval(RunCB{1})
		else						feval(RunCB{1},RunCB{2:end})
		end
	end

%--------------------------------------------------
function KeyPress_local(obj,evt,h)

if (~ishandle(h)),		return,		end
s = getappdata(h,'polygon_data');
key = get(s.h_fig, 'CurrentCharacter');
z = getappdata(s.h_pol,'ZData');

switch key
    case 'r'								% delete vertex
		if (isempty(s.hCurrentMarker)),	return,		end
		x = get(s.h_pol,'XData');           x(s.vert_index) = [];
		y = get(s.h_pol,'YData');           y(s.vert_index) = [];
		if (~isempty(z)),	z(s.vert_index) = [];	end
		delete(s.hCurrentMarker);         s.hCurrentMarker = [];
		s.vert_index = [];
		if (~s.duplicate),	set(s.h_pol,'XData',x,'YData',y)		% Update data
		else				set([s.h_pol s.h_vert],'XData',x,'YData',y)
		end
		if (~isempty(z)),		setappdata(s.h_pol,'ZData',z),		end
		setappdata(h,'edited',true)

	case 'i'								% insert vertex
		pt = get(s.h_ax, 'CurrentPoint');
		x = get(s.h_pol,'XData');		y = get(s.h_pol,'YData');
		if (size(x,1) > 1),		x=x(:)';	y=y(:)';	z=z(:)';	end
		if (isempty(s.vert_index))			% No current vertex selected
			[x, y, z] = insert_pt(x, y, z, [pt(1) pt(1,2)], s.h_ax);	% Make a good guess of the insertion point
		else								% Add new point after selected vertex (kept for backward compatibility)
			x = [x(1:s.vert_index) pt(1,1) x(s.vert_index+1:end)];
			y = [y(1:s.vert_index) pt(1,2) y(s.vert_index+1:end)];
			if (~isempty(z))
				n = s.vert_index;
				z = [z(1:n) inter_Z(x(n-1:n), y(n-1:n), z(n-1:n), pt) z(n+1:end)];
			end
			s.vert_index = s.vert_index+1;
			set(s.hCurrentMarker,'XData',pt(1,1),'YData',pt(1,2));
		end
		if (~s.duplicate),	set(s.h_pol,'XData',x,'YData',y)		% Update data
		else				set([s.h_pol s.h_vert],'XData',x,'YData',y)
		end
		if (~isempty(z)),		setappdata(s.h_pol,'ZData',z),		end
		setappdata(h,'edited',true)

	case {'b', 'B'}						% break line
		if (isempty(s.vert_index)),		return,		end		% No reference vertice selected
		if (s.is_patch) || isempty(s.hCurrentMarker),		return,		end		% Patches don't break & marker needed
		x = get(s.h_pol,'XData');	y = get(s.h_pol,'YData');                
		x1 = x(1:s.vert_index);     x2 = x(s.vert_index:end);
		y1 = y(1:s.vert_index);     y2 = y(s.vert_index:end);

		if (numel(x1) == 1 || numel(x2) == 1),		return,		end         % We don't break at extremities

		if (~s.duplicate),	set(s.h_pol,'XData',x1,'YData',y1)
		else				set([s.h_pol s.h_vert],'XData',x1,'YData',y1)
		end
		if (~isempty(z))
			z1 = z(1:s.vert_index);     z2 = z(s.vert_index:end);
			setappdata(s.h_pol,'ZData',z1)
		end

		% Now make the a new segment from rest of the original (but without markers)
		lc = get(s.h_pol,'Color');			ls = get(s.h_pol,'LineStyle');
		lw = get(s.h_pol,'LineWidth');		lT = get(s.h_pol,'Tag');
		lI = getappdata(s.h_pol,'LineInfo');
		ud = get(s.h_pol,'UserData');
		% create a new line handle
		tmp = line('XData',x2,'YData',y2,'Parent',s.h_ax,'LineWidth',lw,'Color',lc,'LineStyle',ls, 'UserData',ud);
		if (~isempty(z)),		setappdata(tmp, 'ZData',z2),	end
		if (~isempty(lT)),		set(tmp, 'Tag', lT),			end
		if (~isempty(lI)),		setappdata(tmp, 'LineInfo', lI),end 
		set(tmp,'uicontextmenu',get(s.h_pol,'uicontextmenu'))   % Copy the uicontextmenu
		ui_edit_polygon(tmp)
		s.is_closed = false;		% It's not closed anymore
		uistack_j(tmp,'bottom')		% I'm not yet ready to accept this bloated op
		setappdata(h,'edited',true)
		setappdata(tmp,'edited',true)
		s2 = getappdata(tmp,'polygon_data');		% Here we need to (re)set the correct KeyPress_orig
		s2.KeyPress_orig = s.KeyPress_orig;
		setappdata(s2.h_pol,'polygon_data',s2)

	case 'c'						% close line
		if (s.is_patch || s.is_closed),		return,		end		% Don't close what is already closed
		x = get(s.h_pol,'XData');
		if (length(x) <= 2),	return,		end			% don't close a line with less than 2 vertex 
		y = get(s.h_pol,'YData');
		set(s.h_pol,'XData',[x x(1)],'YData',[y y(1)]);
		if (~isempty(z)),	setappdata(s.h_pol,'ZData',[z z(1)]),		end
		s.is_closed = true;
		setappdata(h,'edited',true)
		cmenuHand = get(h, 'UIContextMenu');
		uimenu(cmenuHand, 'Label', 'Create Mask', 'Call', 'poly2mask_fig(guidata(gcbo),gco)');

	case 'e'						% edit (extend) line with getline_j
		if (s.duplicate),	delete(s.h_vert);		s.h_vert = [];		end
		s.vert_index = [];			% delete vertex markers
		if (~isempty(s.hCurrentMarker) && ishandle(s.hCurrentMarker))
			delete(s.hCurrentMarker);		s.hCurrentMarker = [];
		end
		s.controls = 'off';
		set(s.h_pol,'buttondownfcn',@polygonui);
		set(s.h_fig,'KeyPressFcn',s.KeyPress_orig)
		setappdata(s.h_fig,'epActivHand',0)
		setappdata(s.h_pol,'polygon_data',s)		% Play safe
		if (~isempty(z)),	setappdata(s.h_pol,'ZData',[]),	end		% getline_j doesn't handle 3D lines
		n1 = numel(z);
		[x,y] = getline_j(s.h_pol);
		n2 = numel(x);
		set(s.h_pol, 'XData',x, 'YData',y);
		if (~isempty(z)),	setappdata(s.h_pol,'ZData',[z repmat(z(end),1,n2-n1)]);	end		% Reset the Zs
		setappdata(h,'edited',true)
		return
                       
	case {'n', 'p'}							% Move active vertex one step forward or backward
		if (isempty(s.vert_index))			% No current vertex selected
			return
		end
		x = get(s.h_pol,'XData');		y = get(s.h_pol,'YData');
		if (key == 'n')						% Next vertex
			s.vert_index = min(s.vert_index+1, numel(x));
		else								% Previous vertex
			s.vert_index = max(s.vert_index-1, 1);
		end
		set(s.hCurrentMarker,'XData',x(s.vert_index),'YData',y(s.vert_index));

	case 'P'							% close line -> patch
		if (s.is_patch || numel(get(s.h_pol,'XData')) <= 2);  return;     end	
		p = patch('XData',get(s.h_pol,'XData'), 'YData',get(s.h_pol,'YData'), 'parent',s.h_ax, ...
			'FaceColor','none', 'EdgeColor',get(s.h_pol,'Color'));
		draw_funs(p,'line_uicontext')
		if (~isempty(z)),	setappdata(p,'ZData',z),		end
		s_old = getappdata(s.h_pol,'polygon_data');
		s_old.h_pol = p;							% Need to update for the correct handle
		s_old.is_patch = true;
		if (s.is_closed),	s_old.is_closed_patch = true;	% This not an "auto-closed" patch
		else				s_old.is_closed_patch = false;
		end
		setappdata(p,'polygon_data',s_old);			% Set the corrected appdata
		polygonui(s.h_pol)							% Remove the markers
		delete(s.h_pol)								% Delete the ancestor
		setappdata(h,'edited',true)
		ui_edit_polygon(p)							% Start all over
		return

    case char(27)                   % escape:   stop editing
		set(s.h_fig,'selectiontype','open')
		polygonui(s.h_pol)
		return

    case char(127)                  % delete:   delete polygon
		set(s.h_fig,'selectiontype','open')
		polygonui(s.h_pol)
		delete(s.h_pol)
		return
end

setappdata(s.h_pol,'polygon_data',s)
    
%--------------------------------------------------
function state = uisuspend_safe(h_fig)
% Workaround function to avoid a probable bug in R13 uisuspend
% The point is that we need to set the 'KeyPressFcn' to {@KeyPress_local, s.h_pol}
% but the get(fig, 'KeyPressFcn') command returns a cell array with
% [@KeyPress_local] & [s.h_pol handle] 
% and this turns the uistate structure into a 2x1 with all fields repeated.
% In consequence, uisuspend breaks on the "get(uistate.children,..."
% complaining on "too many inputs"

	KPs = get(h_fig,'KeyPressFcn');
	if (iscell(KPs))
		set(h_fig,'KeyPressFcn',KPs{1})         % Temporarly forget the line/patch handle
	else
		try
			set(h_fig,'KeyPressFcn',KPs(1))     % Temporarly forget the line/patch handle
		end
	end

	state = uisuspend_j(h_fig);         % Remember initial figure state
	state.KeyPressFcn = KPs;            % Set to its true value for use in uirestore

%--------------------------------------------------
function [x, y, z] = insert_pt(x_in, y_in, z_in, pt, hAx)
% Guess where to insert he new point. The idea is to insert the new point
% between the closest and the next closest points with respect to the current point.
% Because of the vectorized programming memory eater, we work on data chunks.

	is_subset = false;		type_cast = false;
	if (numel(x_in) > 1000)		% Ad-hoc
		[x, y, z, idS, idE, is_subset] = get_subset(hAx, x_in, y_in, z_in, pt);
		if (idS == idE)			% A crazy click or some other unknown algo failure. Just return
			x = x_in;	y = y_in;	z = z_in;
			return
		end
	else
		x = x_in;	y = y_in;	z = z_in;		% No worries, this won't consume any extra memory
	end
	
	if (~isa(x,'double'))		% R13 ...
		x = double(x);		y = double(y);		z = double(z);	% Won't be very bad because of the subset
		type_cast = true;
	end

	r = (((pt(1) - x)).^2 + ((pt(2) - y) ).^2);		% Squares of distances
	[r_min,i] = min(r);
	
	if (i == 1)						% Near beginning of line ambiguity
		a = ( (x(2) - x(1))^2 + (y(2) - y(1))^2 );		% square distance from 1rst to 2nth points
		hypot_pitag = ( a + r(i) );
		if ( r(i+1) > hypot_pitag )	% Current point is closer to first point. Add the new one before it
			if (~isempty(z)),	z = [z(1) z];	end			% Do not extrapolate (!?)
			x = [pt(1) x];		y = [pt(2) y];
		else						% Current point lyies between the first two points
			if (~isempty(z)),	z = [z(1) inter_Z(x(1:2), y(1:2), z(1:2), pt) z(2)];	end
			x = [x(1) pt(1) x(2:end)];		y = [y(1) pt(2) y(2:end)];
		end
		
	elseif (i == numel(r))			% Near end of line ambiguity
		a = ( (x(end) - x(end-1))^2 + (y(end) - y(end-1))^2 );		% square distance between the two last points
		hypot_pitag = ( a + r(i) );
		if ( r(i-1) > hypot_pitag )	% Current point is closer to the end point. Add the new one after it
			if (~isempty(z)),	z = [z z(1)];		end			% Do not extrapolate (!?)
			x = [x pt(1)];		y = [y pt(2)];
		else						% New point is between the two last points
			if (~isempty(z)),	z = [z(1:i-1) inter_Z(x(i-1:i), y(i-1:i), z(i-1:i), pt) z(end)];	end
			x = [x(1:i-1) pt(1) x(end)];
			y = [y(1:i-1) pt(2) y(end)];	
		end
	else							% Other cases
		a = ( (x(i) - x(i-1))^2 + (y(i) - y(i-1))^2 );		% square distance between current and before points
		hypot_pitag = ( a + r(i) );
		if ( r(i-1) < hypot_pitag )		% Insert point is in the interval [previous_point closest_point]
			if (~isempty(z)),	z = [z(1:i-1) inter_Z(x(i-1:i), y(i-1:i), z(i-1:i), pt) z(i:end)];	end
			x = [x(1:i-1) pt(1) x(i:end)];
			y = [y(1:i-1) pt(2) y(i:end)];			
		else						% Insert point is in the interval [closest_point next_point]
			if (~isempty(z)),	z = [z(1:i) inter_Z(x(i:i+1), y(i:i+1), z(i:i+1), pt) z(i+1:end)];	end
			x = [x(1:i) pt(1) x(i+1:end)];
			y = [y(1:i) pt(2) y(i+1:end)];
		end
	end

	if (type_cast)		% Cast back to singles (it they were INTs before, we may have a bug here) 
		x = single(x);		y = single(y);		z = single(z);
	end

	if (is_subset)
		x = [x_in(1:idS-1) x x_in(idE+1:end)];
		y = [y_in(1:idS-1) y y_in(idE+1:end)];
		if (~isempty(z)),	z = [z_in(1:idS-1) z z_in(idE+1:end)];	end
	end

%--------------------------------------------------
function [x, y, z, idS, idE, is_subset] = get_subset(hAx, x, y, z, pt)
% Get a subset of the [x y z] vector centered on currrent point PT.
% One problem is to decide the width of the window inside which we will get the data.
% If it is too small we risk to get no points and too big will perhaps find too many
% points and thus we still fall into memory eating monster of the vectorized programming.
% For starting we'll use 1/3 of the size of currently display area wich, don't forget,
% depends on the zooming level. If no points are found in first attempt widen the searching
% window to the size of the displayed data (but still centered on current point).
% If that one still fails, than return the same data as input.

	idS = 0;	idE = 0;
	is_subset = true;		% Will be turn off if no subset is found
	x_lim = get(hAx,'xlim');		y_lim = get(hAx,'ylim');

	dx = diff(x_lim) / 6;	dy = diff(y_lim) / 6;
	x0 = pt(1) - dx;		x1 = pt(1) + dx;
	y0 = pt(2) - dy;		y1 = pt(2) + dy;
	ind = ( x > x0 & x < x1 & y > y0 & y < y1);
	if (any(ind))
		id = find(ind);
		idS = id(1);		idE = id(end);		% Start and end indices of chunk data
		x = x(idS:idE);		y = y(idS:idE);
		if (~isempty(z)),	z = z(idS:idE);		end
	else					% NOTHING??. Make another attempt with a bigger window
		dx = diff(x_lim) / 2;	dy = diff(y_lim) / 2;
		x0 = pt(1) - dx;		x1 = pt(1) + dx;
		y0 = pt(2) - dy;		y1 = pt(2) + dy;
		ind = ( x > x0 & x < x1 & y > y0 & y < y1);
		if (any(ind))		% If this is still empty than return the same data as input
			id = find(ind);
			idS = id(1);		idE = id(end);
			x = x(idS:idE);		y = y(idS:idE);
			if (~isempty(z)),	z = z(idS:idE);		end
		else
			is_subset = false;
		end
	end

%--------------------------------------------------
function z = inter_Z(x, y, z, pt)
% Do a crude estimate of the Z value by linear interpolation both in Z and horizontally
	ri = (pt(1) - x(1))^2 + (pt(2) - y(1))^2;
	r = (x(2) - x(1))^2 + (y(2) - y(1))^2;
	z = z(1) + (z(2)-z(1)) * ri / r;


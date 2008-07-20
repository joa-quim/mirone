function ui_edit_polygon(varargin)
% Interactive edit polylines (closed or open) and patches.
%
% ui_edit_polygon(handle1, handle2, ..., 'move_choice')
% ui_edit_polygon([handle1, handle2, ...],'move_choice')  
%       (..., where handle_i is a handle to a line or patch)
%
% Interactive edit polylines (closed or open) and patches.
% If obj(handle_i) has allready been made editable with
% ui_edit_polygon, the function will check the edit state 
% of obj and turn it off if necessary
%
% 	MOVE_CHOICE			% Controls what mouse selection is used to move the whole polygon.
% 						% If == empty, polygon is moved with a left click. That has often anoying side effects
% 						% If == 'extend', means that we need a Shift-click left mouse button or click both
% 						% left and right mouse buttons to move it. A bit more cumbersome, but safer.
%		An alternative way to set the MOVE_CHOICE option is to put it as an appdata of the axes where lines are drawn.
%		As an example of how it works see this snipet showing how it is recovered here
% 				hAx = get(varargin{1},'parent');				% Get parent axes
% 				move_choice = getappdata(hAx, 'MovPolyg');		% Get MOVE_CHOICE option
% 
% --------------------------------------------------------
% Double clicking on the line or patch displays its vertices as control points.
% To stop edit double click again. 
% To edit vertex, click and drag control point.
% Click and drag outside the line/patch vertices control points moves the object (but see 'MOVE_CHOICE')
%
%  In edit mode, hit the following keybord keys ...
%     "r", "-", or backspace:		remove the active vertex
%     "i" or "+":					inserts a new vertex after the active vertex, or in case there is
%									no active vertex, guess the right order where to insert the new
%									vertex in such a way as not to crete zig-zags.
%     "b":							breaks a line in two parts. One from the first point up to the
%									active vertice and the other from the active vertice + 1 till the end.
%									Patches cannot be broken.
%     "c":							close a polyline
%     "p":							close a polyline and convert it into a patch
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
% -------
%   Joaquim Luis (jluis@ualg.pt) - Original version
%       25-Oct-2005 Updated version with bug corrections and added the rectangle mode descrimination. 
%   Sebastian Hoelz (hoelz@geophysik.tu-berlin.de)
%       10-Nov-2005 Nice improvments and code cleaning
%   JL & SH joint version -- ??-Nov-2005
%	JL	20-Jul-2008 Guess best positions for the insert option
%					Allow choosing what mouse selection is used to move the whole polygon.

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

	if (isa(varargin{end},'char'))	% Moving polygon option was transmitted in input
		move_choice = varargin{end};
		varargin(end) = [];
		if ( ~strcmp(move_choice, 'extend') )		% Besides [], the only other possible value is 'extend'
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
			% Creating data-structure for polygon
			s.h_pol = varargin{i};
			s.h_vert = [];
			s.h_current_marker = [];
			s.h_ax = get(s.h_pol,'parent');
			s.h_fig = get(s.h_ax,'parent');
			s.KeyPress_orig = get(s.h_fig,'KeyPressFcn');
	
			s.controls = 'off';
			s.vert_index = [];
			s.save_x = [];
			s.save_y = [];
	
			x = get(s.h_pol,'XData');
			y = get(s.h_pol,'YData');
			s.is_closed_patch = 0;
			s.is_patch  = strcmpi(get(s.h_pol,'Type'),'patch');
			s.is_rect = 0;      s.keep_rect = 0;
			if ( numel(x) == 5 && (x(1) == x(2)) && (x(3) == x(4)) && (y(1) == y(4)) && (y(2) == y(3)) )
				s.is_rect = 1;      s.keep_rect = 1;
			elseif ( numel(x) == 5 && (x(1) == x(4)) && (x(2) == x(3)) && (y(1) == y(2)) && (y(3) == y(4)) )
				s.is_rect = 1;      s.keep_rect = 0;
			end
			if ( numel(x) > 1 && (x(1) == x(end)) && (y(1) == y(end)))
				s.is_closed = 1;
				if (s.is_patch),    s.is_closed_patch = 1;  end
			else
				s.is_closed = 0;
			end
			
			s.what_move = move_choice;
	
			set(s.h_pol,'buttondownfcn',@polygonui);
			setappdata(s.h_pol,'polygon_data',s)
		end
	end

%--------------------------------------------------
function polygonui(varargin)

	if (~ishandle(varargin{1})),		return,		end
	s = getappdata(varargin{1},'polygon_data');
	stype = get(s.h_fig,'selectiontype');
	if (~strcmpi(stype,'open')),	return,		end

	switch s.controls
		case 'on'
			delete(s.h_vert);        s.h_vert = [];		s.vert_index = [];		% delete vertice markers
			try delete(s.h_current_marker); s.h_current_marker = []; end
			s.controls = 'off';
			set(s.h_pol,'buttondownfcn',@polygonui);
			set(s.h_fig,'KeyPressFcn',s.KeyPress_orig)
			setappdata(s.h_pol,'polygon_data',s)
			setappdata(s.h_fig,'epActivHand',0)

		case 'off'
			% Make sure that only one polygon is edited at a time
			h_active = getappdata(s.h_fig,'epActivHand');
			if (h_active),		polygonui(h_active),	end
			setappdata(s.h_fig,'epActivHand',s.h_pol)

			s.controls = 'on';
			x = get(s.h_pol,'XData');
			y = get(s.h_pol,'YData');
			s.h_vert = line('xdata',x,'ydata',y,'Parent',s.h_ax, 'Marker','s','color','r', 'MarkerFaceColor','none', ...
							'linestyle','none','MarkerSize',5,'buttondownfcn',{@edit_polygon,s.h_pol});
			set(s.h_pol,'buttondownfcn',{@move_polygon,s.h_pol});
			set(s.h_fig,'KeyPressFcn',{@KeyPress_local, s.h_pol})
			setappdata(s.h_pol,'polygon_data',s)
	end

%--------------------------------------------------
function edit_polygon(obj,eventdata,h)
% Edit the polygon whose handle is h
	s = getappdata(h,'polygon_data');
	stype = get(s.h_fig,'selectiontype');
	if strcmp(stype,'open')					% When a line has many vertices, the markers may completely
		polygonui(s.h_pol,eventdata)		% hide it. So the only way to get out of edition mode is
		return								% to provide an other exit. That's where this call to
	end										% polygonui('markermousedown') comes to hand.        
	state = uisuspend_safe(s.h_fig);		% Remember initial figure state
	
	x_lim = get(s.h_ax,'xlim');
	y_lim = get(s.h_ax,'ylim');
	current_pt = get(s.h_ax, 'CurrentPoint');
	% Find out which vertex is beeing edited
	x = get(s.h_pol,'XData');   y = get(s.h_pol,'YData');
	dif_x = x - current_pt(1,1);    dif_y = y - current_pt(1,2);
	dist = sqrt(dif_x.^2 + dif_y.^2);
	[B,IX] = sort(dist);    s.vert_index = IX(1);
	s.save_x = x(s.vert_index);   s.save_y = y(s.vert_index);   % needed in the "i"nsert option
	if isempty(s.h_current_marker)		% If the black marker doesn't exist, creat it
		s.h_current_marker = line('xdata',s.save_x,'ydata', s.save_y, 'parent', s.h_ax,'Marker','s','MarkerEdgeColor','none', ...
								'MarkerFaceColor','k','MarkerSize',5,'Tag','current_marker');
	else								% The black marker exists, just move it.
		set(s.h_current_marker,'XData',s.save_x,'YData',s.save_y)
	end
	setappdata(h,'polygon_data',s);
	
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

	if (s.is_rect && s.keep_rect)                % We are dealing with a line/patch rectangle
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
	elseif (s.vert_index == 1)                  % Selected first vertice
        if (s.is_closed && s.is_closed_patch)
            x = [newx xx(2:end-1) newx];    y = [newy yy(2:end-1) newy];
		elseif (s.is_closed && ~s.is_closed_patch)
			if (s.is_patch)
				x = [newx xx(2:end)];           y = [newy yy(2:end)];
			else                                % deformable rectangle (rect given in cw direction)
				x = [newx xx(2:end-1) newx];    y = [newy yy(2:end-1) newy];
			end
		else								% "Open" polyline
			x = [newx xx(2:end)];           y = [newy yy(2:end)];
		end
	elseif (s.vert_index == length(xx) && ~s.is_patch)     % Selected last polyline vertice
		if (s.is_closed)
			x = [newx xx(2:end-1) newx];    y = [newy yy(2:end-1) newy];
		else								% "Open" polyline
			x = [xx(1:end-1) newx];         y = [yy(1:end-1) newy];
		end
	else                                    % Midle vertices
		x = [xx(1:s.vert_index-1) newx xx(s.vert_index+1:end)];
		y = [yy(1:s.vert_index-1) newy yy(s.vert_index+1:end)];
	end

	set(s.h_pol, 'XData',x, 'YData',y);
	set(s.h_current_marker,'XData',newx,'YData',newy)   % Move the current point marker together with the editing point
	set(s.h_vert,'XData',x, 'YData',y)                  % Move the marker together with the editing point

%--------------------------------------------------
function wbu_EditPolygon(obj,eventdata,h,state)
	s = getappdata(h,'polygon_data');
	uistack(s.h_vert,'top')			% Need to do this because the black marker has no ButtonDown and was on top
	uirestore_fig(state);			% Restore the figure's initial state

%--------------------------------------------------
function move_polygon(obj,eventdata,h)
% Move the polygon whose handle is h
	s = getappdata(h,'polygon_data');
	stype = get(s.h_fig,'selectiontype');

	if (strcmp(stype,'open')),		polygonui(s.h_pol,[]),		return,		end
	if (~isempty(s.what_move) && ~strcmp(stype,s.what_move)),	return,		end		% needs "password" to continue

	state = uisuspend_safe(s.h_fig);                 % Remember initial figure state    
	x_lim = get(s.h_ax,'xlim');        y_lim = get(s.h_ax,'ylim');
	current_pt = get(s.h_ax, 'CurrentPoint');
	setappdata(s.h_pol,'old_pt',[current_pt(1,1) current_pt(1,2)])

	set(s.h_fig,'WindowButtonMotionFcn',{@wbm_MovePolygon,h,[x_lim y_lim]},...
		'WindowButtonUpFcn',{@wbu_MovePolygon,h,state}, 'Pointer','fleur');

%--------------------------------------------------
function wbm_MovePolygon(obj,eventdata,h,lim)
	s = getappdata(h,'polygon_data');
	pt = get(s.h_ax, 'CurrentPoint');
	if (pt(1,1) < lim(1)) || (pt(1,1) > lim(2)) || (pt(1,2) < lim(3)) || (pt(1,2) > lim(4));   return; end
	old_pt = getappdata(s.h_pol,'old_pt');
	xx = get(s.h_pol,'XData');      yy = get(s.h_pol,'YData');
	xv = get(s.h_vert,'XData');     yv = get(s.h_vert,'YData');
	dx = pt(1,1) - old_pt(1);       dy = pt(1,2) - old_pt(2);
	xx = xx + dx;                   yy = yy + dy;
	xv = xv + dx;                   yv = yv + dy;
	setappdata(s.h_pol,'old_pt',[pt(1,1) pt(1,2)])
	set(s.h_pol, 'XData',xx, 'YData',yy);
	set(s.h_vert, 'XData',xv, 'YData',yv);

	if (~isempty(s.h_current_marker))      % If the black marker exists, move it also
		x = get(s.h_current_marker,'XData');      y = get(s.h_current_marker,'YData');
		x = x + dx;     y = y + dy;
		set(s.h_current_marker,'XData',x,'YData',y)
	end

%--------------------------------------------------
function wbu_MovePolygon(obj,eventdata,h,state)
	s = getappdata(h,'polygon_data');
	uistack(s.h_vert,'top')     % Need to do this because the black marker has no ButtonDown and was on top
	uirestore_fig(state);         % Restore the figure's initial state

%--------------------------------------------------
function KeyPress_local(obj,eventdata,h)

if (~ishandle(h)),		return,		end
s = getappdata(h,'polygon_data');
key = get(s.h_fig, 'CurrentCharacter');

switch key
    case {'r', 'R', '-'}				% delete vertex
		if (isempty(s.h_current_marker)),	return,		end
		x = get(s.h_pol,'XData');           x(s.vert_index) = [];
		y = get(s.h_pol,'YData');           y(s.vert_index) = [];
		delete(s.h_current_marker);         s.h_current_marker = [];
		s.vert_index = [];
		set([s.h_pol s.h_vert],'XData',x,'YData',y);	% Update data
                       
	case {'i', 'I', '+'}				% insert vertex
		pt = get(s.h_ax, 'CurrentPoint');
		x = get(s.h_pol,'XData');   x=x(:)';
		y = get(s.h_pol,'YData');   y=y(:)';
		if (isempty(s.vert_index))			% No current vertex selected
			[x, y] = insert_pt(x, y, [pt(1,1) pt(1,2)]);		% Make a good guess of the insertion point
		else								% Add new point after selected vertex (kept for backword compatibility)
			x = [x(1:s.vert_index) pt(1,1) x(s.vert_index+1:end)];
			y = [y(1:s.vert_index) pt(1,2) y(s.vert_index+1:end)];
			s.vert_index = s.vert_index+1;
			set(s.h_current_marker,'XData',pt(1,1),'YData',pt(1,2));
		end
		set([s.h_pol s.h_vert],'XData',x,'YData',y);

	case {'b', 'B'}				% break line
		if (isempty(s.vert_index)),		return,		end		% No reference vertice selected
		if (s.is_patch) || isempty(s.h_current_marker),		return,		end		% Patches don't break & marker needed
		x = get(s.h_pol,'XData');    
		y = get(s.h_pol,'YData');                
		x1 = x(1:s.vert_index);     x2 = x(s.vert_index:end);
		y1 = y(1:s.vert_index);     y2 = y(s.vert_index:end);

		if (numel(x1) == 1 || numel(x2) == 1)         % We don't break at extremities
			return
		end
		set([s.h_pol s.h_vert],'XData',x1,'YData',y1);

		% Now make the a new segment from rest of the original (but without markers)
		lc = get(s.h_pol,'Color');			ls = get(s.h_pol,'LineStyle');
		lw = get(s.h_pol,'LineWidth');		lT = get(s.h_pol,'Tag');
		ud = get(s.h_pol,'UserData');
		% create a new line handle
		tmp = line('XData',x2,'YData',y2,'LineWidth',lw,'Color',lc,'LineStyle',ls, 'UserData',ud);
		if (~isempty(lT)),		set(tmp,'Tag', lT),		end
		set(tmp,'uicontextmenu',get(s.h_pol,'uicontextmenu'))   % Copy the uicontextmenu
		ui_edit_polygon(tmp)
		s.is_closed = 0;            % Its not closed anymore
		%uistack(tmp,'bottom')      % I'm not yet ready to accept this bloated op

	case {'c', 'C'}					% close line
		if (s.is_patch || s.is_closed),		return,		end		% Don't close what is already closed
		x = get(s.h_pol,'XData');
		if (length(x) <= 2); return; end   % don't close a line with less than 2 vertex 
		y = get(s.h_pol,'YData');
		set(s.h_pol,'XData',[x x(1)],'YData',[y y(1)]);
		s.is_closed = 1;

	case {'p', 'P'}					% close line -> patch
		% Don't close what is already closed or line with less than 2 vertices
		if (s.is_patch || s.is_closed || length(get(s.h_pol,'XData')) <= 2);  return;     end	
		p = patch(get(s.h_pol,'xData'),get(s.h_pol,'yData'),1,'parent',s.h_ax);
		%p = patch(get(s.h_pol,'xData'),get(s.h_pol,'yData'),1,'Facecolor','none','parent',s.h_ax);
		s_old = getappdata(s.h_pol,'polygon_data');
		s_old.h_pol = p;							% Need to update for the correct handle
		s.is_patch = 1;
		s.is_closed_patch = 0;						% This not an "auto-closed" patch
		setappdata(p,'polygon_data',s_old);			% Set the corrected appdata
		uictx = get(s.h_pol,'uicontextmenu');		% Get eventual uicontextmenu
		if (~isempty(uictx)),   set(p,'uicontextmenu',uictx);   end
		polygonui(s.h_pol)							% Remove the markers
		delete(s.h_pol)								% Delete the ancestor
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

state = uisuspend_fig(h_fig);         % Remember initial figure state
state.KeyPressFcn = KPs;            % Set to its true value for use in uirestore

%--------------------------------------------------
function [x, y] = insert_pt(x, y, pt)
% Guess where to insert he new point. The idea is to insert the new point
% between the closest and the next closest points with respect to the current point

	r = sqrt(((pt(1) - x)).^2 + ((pt(2) - y) ).^2);
	[temp,i] = min(r);
	
	if (i == 1)						% Near beginning of line ambiguity
		a = sqrt( (x(2) - x(1))^2 + (y(2) - y(1))^2 );		% distance from 1rst to 2nth points
		hypot_pitag = sqrt( a^2 + r(i)^2 );
		if ( r(i+1) > hypot_pitag )		% Current point is closer to first point. Add the new one before it
			x = [pt(1) x];			y = [pt(2) y];
		else						% Current point lyies between the first two points
			x = [x(1) pt(1) x(2:end)];		y = [y(1) pt(2) y(2:end)];
		end
		
	elseif (i == numel(r))			% Near end of line ambiguity
		a = sqrt( (x(end) - x(end-1))^2 + (y(end) - y(end-1))^2 );		% distance between the two last points
		hypot_pitag = sqrt( a^2 + r(i)^2 );
		if ( r(i) > hypot_pitag )		% Current point is closer to the end point. Add the new one after it
			x = [x pt(1)];		y = [y pt(2)];
		else						% New point is between the two last points
			x = [x(1:i-1) pt(1) x(end)];
			y = [y(1:i-1) pt(2) y(end)];	
		end
	else							% Other cases
		a = sqrt( (x(i) - x(i-1))^2 + (y(i) - y(i-1))^2 );		% distance between current and before points
		hypot_pitag = sqrt( a^2 + r(i)^2 );
		if ( r(i-1) < hypot_pitag )		% Insert point is in the interval [previous_point closest_point]
			x = [x(1:i-1) pt(1) x(i:end)];
			y = [y(1:i-1) pt(2) y(i:end)];			
		else						% Insert point is in the interval [closest_point next_point]
			x = [x(1:i) pt(1) x(i+1:end)];
			y = [y(1:i) pt(2) y(i+1:end)];
		end
	end

function ui_edit_polygon(varargin)
% ui_edit_polygon:  Interactive edit polylines (closed or open) and patches.
%
%  ui_edit_polygon: Double clicking on the line or patch displays its vertices as control
%  points. To edit, click and drag the control point. To stop edit double click again.
%
%  Click and drag outside the line/patch vertices control points moves the object
%
%  If the element is a closed polyline or a patch, while in edit mode:
%  Hiting "r" , delete or backspace removes the active vertice
%    "    "i" inserts a new vertice before the active vertice
%    "    "b" breaks the line in two. One from the first point up to the active vertice and the
%  other from the active vertice + 1 till the end. Patches cannot be broken.
%  Hiting "c" closes the polyline.
%  NOTE. Possybly due to the coexistance of this function with others  that
%  also use the WindowButton... functions to do other things, I'm forced to use uisuspend/uirestore
%  A side effect of this is that hiting keys only works while in "click and drag" mode.

% Usage
%      ui_edit_polygon(hand) 
% Where HAND is a line or patch object handle
%
%     Joaquim Luis
%     jluis@ualg.pt

if (nargin == 1 & isstr(varargin{1}))         % ... callbacks
    polygonui(varargin{:});
    return
elseif (nargin == 1 & ishandle(varargin{1}))  % Handle in input
    if ~(strcmp(get(varargin{1},'Type'),'patch') | strcmp(get(varargin{1},'Type'),'line'))
        errordlg('Error in ui_edit_polygon: Invalid handle.','Error')
        return
    end
    h_pol = varargin{1};
else
    errordlg('Error in ui_edit_polygon: Invalid number of inputs.','Error')
    return
end

x = get(h_pol,'XData');     y = get(h_pol,'YData');
%  Set the correct structures
if ( length(x) == 5 & (x(1) == x(2)) & (x(3) == x(4)) & (y(1) == y(4)) & (y(2) == y(3)) )
        s.is_rect = 1;
else    s.is_rect = 0;      end
s.controls = 'off';
s.vert_index = [];
s.save_x = [];
s.save_y = [];
s.h_pol = h_pol;
s.hcenter = [];
s.hvert = [];
s.h_current_marker = [];        % Handle to the (selected) current marker
s.parent = get(s.h_pol,'parent');	
s.uicontext = get(h_pol,'uicontextmenu');   % Save the eventual uicontextmenu
s.old_ud = get(h_pol,'UserData');           % Save the eventual userdata
if (strcmp(get(h_pol,'Type'),'patch'))
        s.is_patch = 1;
else    s.is_patch = 0;     end
if ( (x(1) == x(end)) & (y(1) == y(end)) )
        s.is_closed = 1;
else    s.is_closed = 0;    end
set(h_pol,'userdata',s,'buttondownfcn','ui_edit_polygon(''poligmousedown'')');

% -----------------------------------------------------------------------------------------
function edit_polygon(obj,eventdata,h)
% Edit the polygon whose handle is h
stype = get(gcf,'selectiontype');
if strcmp(stype,'open')                 % When a line has many vertices, the markers may compltely
    ui_edit_polygon('markermousedown')  % hide it. So the only way to get out of edition mode is
    return                              % to provide an other exit. That's where this call to
end                                     % ui_edit_polygon('markermousedown') comes to hand.
state = uisuspend(gcf);                 % Remember initial figure state

s = get(h,'userdata');
x_lim = get(gca,'xlim');        y_lim = get(gca,'ylim');
current_pt = get(gca, 'CurrentPoint');
% Find out which vertice is beeing edited
x = get(s.h_pol,'XData');   y = get(s.h_pol,'YData');
dif_x = x - current_pt(1,1);    dif_y = y - current_pt(1,2);
dist = sqrt(dif_x.^2 + dif_y.^2);
[B,IX] = sort(dist);    s.vert_index = IX(1);
s.save_x = x(s.vert_index);   s.save_y = y(s.vert_index);   % needed in the "i"nsert option
if isempty(s.h_current_marker)      % If the black marker doesn't exist, creat it
    hold on
    s.h_current_marker = plot(s.save_x,s.save_y,'Marker','s','MarkerEdgeColor','none',...
        'MarkerFaceColor','k','MarkerSize',6,'Tag','current_marker');
    hold off
else                                % The black marker exists, just move it.
    set(s.h_current_marker,'XData',s.save_x,'YData',s.save_y)
end
set(h,'userdata',s);

set(gcf,'WindowButtonMotionFcn',{@wbm_EditPolygon,h,[x_lim y_lim]},...
    'WindowButtonUpFcn',{@wbu_EditPolygon,h,state}, 'KeyPressFcn',{@KeyPress,h}, 'Pointer','fleur');
% set(gcf,'WindowButtonMotionFcn',{@wbm_EditPolygon,h,[x_lim y_lim]},...
%     'WindowButtonUpFcn',{@wbu_EditPolygon,h,state}, 'KeyPressFcn','ui_edit_polygon(''KeyPress'',gco)', 'Pointer','crosshair');

% ---------
function wbm_EditPolygon(obj,eventdata,h,lim)
s = get(h,'userdata');
pt = get(gca, 'CurrentPoint');
if (pt(1,1)<lim(1)) | (pt(1,1)>lim(2)) | (pt(1,2)<lim(3)) | (pt(1,2)>lim(4));   return; end
xx = get(s.h_pol,'XData');      yy = get(s.h_pol,'YData');
xx = xx(:)';                    yy = yy(:)';    % Make sure they are row vectors
newx = pt(1,1);                 newy = pt(1,2);

if (s.is_rect)                          % We are dealing with a rectangle
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
elseif (s.vert_index == 1)              % Selected first polyline vertice
    if (s.is_closed | s.is_patch)
        x = [newx xx(2:end-1) newx];    y = [newy yy(2:end-1) newy];
    else                                % "Open" polyline
        x = [newx xx(2:end)];           y = [newy yy(2:end)];
    end
elseif (s.vert_index == length(xx))     % Selected last polyline vertice
    if (s.is_closed | s.is_patch)
        x = [newx xx(2:end-1) newx];    y = [newy yy(2:end-1) newy];
    else                                % "Open" polyline
        x = [xx(1:end-1) newx];         y = [yy(1:end-1) newy];
    end
else                                    % Midle vertices
    x = [xx(1:s.vert_index-1) newx xx(s.vert_index+1:end)];
    y = [yy(1:s.vert_index-1) newy yy(s.vert_index+1:end)];
end

set(s.h_pol, 'XData',x, 'YData',y);
set(s.h_current_marker,'XData',newx,'YData',newy)  % Move the current point marker together with the editing point
set(s.hvert,'XData',x, 'YData',y)             % Move the marker together with the editing point

% ---------
function wbu_EditPolygon(obj,eventdata,h,state)
s = get(h,'userdata');
uistack(s.hvert,'top')      % Need to do this because the black marker has no ButtonDown and was on top
%set(gcf,'WindowButtonMotionFcn','', 'WindowButtonUpFcn','', 'Pointer', 'arrow')
%state.KeyPressFcn = get(gcf,'KeyPressFcn');
uirestore(state);           % Restore the figure's initial state

%---------------------------------------------------------------------------------------------------
function polygonui(action)

switch action
case 'poligmousedown'              % mouse down on line/patch
	stype = get(gcf,'selectiontype');
	s = get(gco,'userdata');
	s.parent = get(gco,'parent');	
	switch stype % shift-click or double click to toggle control points
		case 'open'
				switch s.controls
					case 'on'
                        delete(s.hvert);        s.hvert = [];       % delete vrtice markers
                        delete(s.h_current_marker); s.h_current_marker = [];
						s.controls = 'off';
                        set(s.h_pol,'uicontextmenu',s.uicontext)
						set(gco,'userdata',s)
                        set(gcf,'KeyPressFcn','')
					case 'off'
						s.controls = 'on';
						hand = gco;
                        x = get(hand,'XData');   y = get(hand,'YData');
                        hold on;    
						s.hvert = plot(x,y,'sk','MarkerFaceColor','none','MarkerSize',6,'userdata',s,...
							  'buttondownfcn',{@edit_polygon,hand});
                        hold off;
                        set(s.h_pol,'buttondownfcn',{@move_polygon,hand},'uicontextmenu','');
						set(hand,'userdata',s)
                        set(s.hvert,'userdata',s)
				end
	end
case 'markermousedown'      % User double-clicked on a marker
    s = get(gco,'userdata');
    s.parent = get(gco,'parent');
    delete(s.hvert);        s.hvert = [];       % delete vrtice markers
    try,        % We may, or may not, have a current marker
        sp = get(s.h_pol,'UserData');
        delete(sp.h_current_marker);    sp.h_current_marker = [];
        set(s.h_pol,'UserData',sp);
    end
    s.controls = 'off';
    set(s.h_pol,'userdata',s,'buttondownfcn','ui_edit_polygon(''poligmousedown'')','uicontextmenu',s.uicontext);
    set(gco,'userdata',s)
    set(gcf,'KeyPressFcn','')
    % The following is a patch for an (don't know where) introduced bug. Apparently, for
    % individual symbols (and I hope its only for them) that were moved, XData and YData
    % had repeted values. This makes a mess in other codes that require knowing the point
    % coordinates.
    xx = get(s.h_pol,'XData');  yy = get(s.h_pol,'YData');
    if ((length(xx) == 2) & (xx(1) == xx(2)) & (yy(1) == yy(2)))
        xx(2) = [];     yy(2) = [];
        set(s.h_pol,'XData',xx,'YData',yy)
    end
end

% -----------------------------------------------------------------------------------------
function move_polygon(obj,eventdata,h)
% Move the polygon whose handle is h
s = get(h,'userdata');
stype = get(gcf,'selectiontype');

if (strcmp(stype,'open'))
    delete(s.hvert);        s.hvert = [];       % delete vrtice markers
    delete(s.h_current_marker); s.h_current_marker = [];
	s.controls = 'off';
    set(s.h_pol,'userdata',s,'buttondownfcn','ui_edit_polygon(''poligmousedown'')','uicontextmenu',s.uicontext);
	set(gco,'userdata',s)
    set(gcf,'KeyPressFcn','')
    return
end

state = uisuspend(gcf);                 % Remember initial figure state
x_lim = get(gca,'xlim');        y_lim = get(gca,'ylim');
current_pt = get(gca, 'CurrentPoint');
setappdata(s.h_pol,'old_pt',[current_pt(1,1) current_pt(1,2)])

set(gcf,'WindowButtonMotionFcn',{@wbm_MovePolygon,h,[x_lim y_lim]},...
    'WindowButtonUpFcn',{@wbu_MovePolygon,h,state}, 'Pointer','fleur');

% ---------
function wbm_MovePolygon(obj,eventdata,h,lim)
s = get(h,'userdata');
pt = get(gca, 'CurrentPoint');
if (pt(1,1)<lim(1)) | (pt(1,1)>lim(2)) | (pt(1,2)<lim(3)) | (pt(1,2)>lim(4));   return; end
old_pt = getappdata(s.h_pol,'old_pt');
xx = get(s.h_pol,'XData');      yy = get(s.h_pol,'YData');
xv = get(s.hvert,'XData');      yv = get(s.hvert,'YData');
dx = pt(1,1) - old_pt(1);       dy = pt(1,2) - old_pt(2);
xx = xx + dx;                   yy = yy + dy;
xv = xv + dx;                   yv = yv + dy;
setappdata(s.h_pol,'old_pt',[pt(1,1) pt(1,2)])
set(s.h_pol, 'XData',xx, 'YData',yy);
set(s.hvert, 'XData',xv, 'YData',yv);

if (~isempty(s.h_current_marker))      % If the black marker exists, move it also
    x = get(s.h_current_marker,'XData');      y = get(s.h_current_marker,'YData');
    x = x + dx;     y = y + dy;
    set(s.h_current_marker,'XData',x,'YData',y)
end

% ---------
function wbu_MovePolygon(obj,eventdata,h,state)
s = get(h,'userdata');
uistack(s.hvert,'top')      % Need to do this because the black marker has no ButtonDown and was on top
uirestore(state);           % Restore the figure's initial state

%--------------------------------------------------
% Subfunction KeyPress
%--------------------------------------------------
function KeyPress(obj,eventdata,h)
% Function to deal with key press while in polyline or patch edit mode
% Hiting "r" , delete or backspace removes the active vertice
%   "    "i" inserts a new vertice before the active vertice
%   "    "b" breaks the line in two. One from the first point up to the active vertice and the
%   other from the active vertice + 1 till the end. Patches cannot be broken.
%   "    "c" closes the polyline.
% Note. The code could be a bit more elegant, but I rather repeat some code and save memory.
s = get(h,'UserData');
key = get(get(s.parent,'Parent'), 'CurrentCharacter');
switch key
    case {char(8), char(127), char(114)}  % "r", delete and backspace keys
        x = get(s.h_pol,'XData');           y = get(s.h_pol,'YData');
        len = length(x);
        x(s.vert_index) = [];               y(s.vert_index) = [];
        if (s.is_closed | s.is_patch)       % Test for special cases
            if (s.vert_index == 1)
                x(end) = x(1);  y(end) = y(1);
            elseif (s.vert_index == len)
                x(1) = x(end);  y(1) = y(end);
            end
        end
        set(s.h_pol,'XData',x,'YData',y);   % Update the vertices
        x = get(s.hvert,'XData');           y = get(s.hvert,'YData');
        x(s.vert_index) = [];               y(s.vert_index) = [];
        if (s.is_closed | s.is_patch)       % Test for special cases
            if (s.vert_index == 1)
                x(end) = x(1);  y(end) = y(1);
            elseif (s.vert_index == len)
                x(1) = x(end);  y(1) = y(end);
            end
        end
        set(s.hvert,'XData',x,'YData',y);   % Update the markers
        delete(s.h_current_marker);         s.h_current_marker = [];
    case char(105)              % "i" char key (insert vertice)
        pt = get(gca, 'CurrentPoint');
        newx = pt(1,1);                     newy = pt(1,2);
        xx = get(s.h_pol,'XData');          yy = get(s.h_pol,'YData');
        xx = xx(:)';                        yy = yy(:)';  % From patches they are column vectors
        xx(s.vert_index) = s.save_x;        yy(s.vert_index) = s.save_y;
        if (s.vert_index < length(xx) & s.vert_index > 1)
            x = [xx(1:s.vert_index) newx xx(s.vert_index+1:end)];   % And we would get an error here
            y = [yy(1:s.vert_index) newy yy(s.vert_index+1:end)];
            set(s.h_pol,'XData',x,'YData',y);
            xx = get(s.hvert,'XData');          yy = get(s.hvert,'YData');
            xx = xx(:)';                        yy = yy(:)';
            xx(s.vert_index) = s.save_x;        yy(s.vert_index) = s.save_y;
            x = [xx(1:s.vert_index) newx xx(s.vert_index+1:end)];
            y = [yy(1:s.vert_index) newy yy(s.vert_index+1:end)];
            set(s.hvert,'XData',x,'YData',y);
        elseif (s.vert_index == 1)        % New vertice is also the first vertice
            x = [newx xx];                      y = [newy yy];
            set(s.h_pol,'XData',x,'YData',y);
            xx = get(s.hvert,'XData');          yy = get(s.hvert,'YData');
            xx = xx(:)';                        yy = yy(:)';
            xx(s.vert_index) = s.save_x;        yy(s.vert_index) = s.save_y;
            x = [newx xx];                      y = [newy yy];
            set(s.hvert,'XData',x,'YData',y);
        else        % New vertice is also the last vertice
            x = [xx newx];                      y = [yy newy];
            set(s.h_pol,'XData',x,'YData',y);
            xx = get(s.hvert,'XData');          yy = get(s.hvert,'YData');
            xx = xx(:)';                        yy = yy(:)';
            xx(s.vert_index) = s.save_x;        yy(s.vert_index) = s.save_y;
            x = [xx newx];                      y = [yy newy];
            set(s.hvert,'XData',x,'YData',y);
        end
        s.vert_index = s.vert_index + 1;        
    case char(98)               % "b" char key (break line)
        if (s.is_patch)         % Patches don't break
            return
        end
        xx = get(s.h_pol,'XData');    yy = get(s.h_pol,'YData');
        xx(s.vert_index) = s.save_x;      yy(s.vert_index) = s.save_y;
        if (s.vert_index > 1 & s.vert_index < length(xx))
            x = xx(1:s.vert_index);            y = yy(1:s.vert_index);
            set(s.h_pol,'XData',x,'YData',y,'Tag','polyline');
            % Delete the markers of the second segment
            xx = get(s.hvert,'XData');         yy = get(s.hvert,'YData');
            x = xx(1:s.vert_index);            y = yy(1:s.vert_index);
            set(s.hvert,'XData',x,'YData',y);
            s.is_closed = 0;    % Its not closed anymore
            % Now make the a new segment from rest of the original (but without markers)
            lc = get(s.h_pol,'Color');      ls = get(s.h_pol,'LineStyle');
            lw = get(s.h_pol,'LineWidth');      
            % create a new line handle
            tmp = line('XData',[],'YData',[],'LineWidth',lw,'Color',lc,'LineStyle',ls,'Tag','polyline');
            set(tmp, 'XData',xx(s.vert_index+1:end), 'YData',yy(s.vert_index+1:end),'uicontextmenu',s.uicontext)
            ui_edit_polygon(tmp)
        end
    case char(99)              % "c" char key (close line)
        if (s.is_patch | s.is_closed)         % Don't close what is already closed
            return
        end
        xx = get(s.h_pol,'XData');    yy = get(s.h_pol,'YData');
        if length(xx) > 2    % don't close a line with less than 2 vertex
            x = [xx xx(1)];                     y = [yy yy(1)];
            set(s.h_pol,'XData',x,'YData',y);
        end
end
set(h, 'UserData', s);

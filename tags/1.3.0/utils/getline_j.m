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
%   [X,Y] = GETLINE_J(...,'closed') animates and returns a closed polygon.
%   [X,Y] = GETLINE_J(...,'freehand') draw a line following the mouse movements.
%

%   Grandfathered syntaxes:
%   XY = GETLINE_J(...) returns output as M-by-2 array; first column is X; second column is Y.

%   Havily hacked version of getline that, contrary to the original, let be compiled.
%   Also the right-click button was reprogramed. Now it removes the previously selected point.
%
%   Joaquim Luis

ud.GETLINE_ISCLOSED = 0;
ud.GETLINE_FREEHAND = 0;
if ((nargin >= 1) && (ischar(varargin{end})))
    str = varargin{end};
    if (str(1) == 'c')                  % getline_j(..., 'closed')
        ud.GETLINE_ISCLOSED = 1;        varargin = varargin(1:end-1);
    elseif (str(1) == 'f')              % getline_j(..., 'freehand')
        ud.GETLINE_FREEHAND = 1;        varargin = varargin(1:end-1);
    end
end

ud.GETLINE_X = [];     ud.GETLINE_Y = [];

if (length(varargin) < 1)
    ud.GETLINE_AX = gca;
    ud.GETLINE_FIG = get(ud.GETLINE_AX, 'Parent');
else
    if (~ishandle(varargin{1}))
        varargout{1} = [];      varargout{2} = [];
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
        otherwise
            varargout{1} = [];      varargout{2} = [];
            error('First argument should be a figure or axes handle');
    end
end
    
xlimorigmode = get(ud.GETLINE_AX,'xlimmode');     ylimorigmode = get(ud.GETLINE_AX,'ylimmode');
set(ud.GETLINE_AX,'xlimmode','manual');           set(ud.GETLINE_AX,'ylimmode','manual');

% Remember initial figure state
old_db = get(ud.GETLINE_FIG, 'DoubleBuffer');
state = uisuspend_fig(ud.GETLINE_FIG);

% Set up initial callbacks for initial stage
set(ud.GETLINE_FIG, 'Pointer', 'crosshair', ...
    'WindowButtonDownFcn', {@FirstButtonDown,ud.GETLINE_FIG},...
    'KeyPressFcn', {@KeyPress,ud.GETLINE_FIG}, 'DoubleBuffer', 'on');

% Bring target figure forward
figure(ud.GETLINE_FIG);

% Initialize the lines to be used for the drag
ud.GETLINE_H1 = line('Parent', ud.GETLINE_AX, ...
                  'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y, ...
                  'Visible', 'off', 'Clipping', 'off', ...
                  'Color', 'k', 'LineStyle', '-');

ud.GETLINE_H2 = line('Parent', ud.GETLINE_AX, ...
                  'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y, ...
                  'Visible', 'off', 'Clipping', 'off', ...
                  'Color', 'w', 'LineStyle', ':');

setappdata(ud.GETLINE_FIG, 'FromGetLine_j', ud);

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
   uirestore_fig(state);
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

try  rmappdata(ud.GETLINE_FIG,'FromGetLine_j');     end

%-------------------------------------------------------------------------------
function KeyPress(obj,eventdata,hfig)
ud = getappdata(hfig, 'FromGetLine_j');

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
case char(99)              % "c" char key (close line)
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
setappdata(hfig, 'FromGetLine_j', ud);

%----------------------------------------------------------------------------------
function FirstButtonDown(obj,eventdata,hfig)
ud = getappdata(hfig, 'FromGetLine_j');
pt = get(ud.GETLINE_AX, 'CurrentPoint');
x = pt(1,1);    y = pt(1,2);

% check if GETLINE_X,GETLINE_Y is inside of axis
x_lim = get(ud.GETLINE_AX,'xlim');      y_lim = get(ud.GETLINE_AX,'ylim');
if (x>=x_lim(1)) && (x<=x_lim(2)) && (y>=y_lim(1)) && (y<=y_lim(2))    % inside axis limits
    ud.GETLINE_X = x;    ud.GETLINE_Y = y;
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
setappdata(hfig, 'FromGetLine_j', ud);

%---------------------------------------------------------------------------------------
function NextButtonDown(obj,eventdata,hfig)
ud = getappdata(hfig, 'FromGetLine_j');

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
    
    if (ud.GETLINE_ISCLOSED)
        ud.GETLINE_X = [ud.GETLINE_X(1:end-1) x ud.GETLINE_X(end)];
        ud.GETLINE_Y = [ud.GETLINE_Y(1:end-1) y ud.GETLINE_Y(end)];
    else
        ud.GETLINE_X = [ud.GETLINE_X x];        ud.GETLINE_Y = [ud.GETLINE_Y y];
    end
    set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
end

if (strcmp(get(hfig,'SelectionType'),'alt') && ~ud.GETLINE_FREEHAND)
    % Right-click, delete previous point
    pt = get(ud.GETLINE_AX, 'CurrentPoint');
    x = pt(1,1);    y = pt(1,2);
    % check if GETLINE_X,GETLINE_Y is inside of axis
    x_lim = get(ud.GETLINE_AX,'xlim');      y_lim = get(ud.GETLINE_AX,'ylim');
    if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2))    % outside axis limits, ignore this ButtonDown
        return
    end    
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
end

if (strcmp(selectionType, 'extend') || strcmp(selectionType, 'open'))    % We're done (midle button or double-click)
    set(ud.GETLINE_H1, 'UserData', 'Completed');
end
setappdata(ud.GETLINE_FIG, 'FromGetLine_j', ud);

%-----------------------------------------------------------------------------------
function ButtonMotion(obj,eventdata,hfig)
ud = getappdata(hfig, 'FromGetLine_j');

pt = get(ud.GETLINE_AX, 'CurrentPoint');
newx = pt(1,1);    newy = pt(1,2);
if (ud.GETLINE_ISCLOSED && (length(ud.GETLINE_X) >= 3))
    x = [ud.GETLINE_X(1:end-1) newx ud.GETLINE_X(end)];
    y = [ud.GETLINE_Y(1:end-1) newy ud.GETLINE_Y(end)];
else
    x = [ud.GETLINE_X newx];    y = [ud.GETLINE_Y newy];
end

if (~ud.GETLINE_FREEHAND)
    set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', x, 'YData', y);
else            % Do a freehand drawing
	set([ud.GETLINE_H1 ud.GETLINE_H2],'xdata',[get(ud.GETLINE_H1,'xdata'),pt(1,1)],...
            'ydata',[get(ud.GETLINE_H1,'ydata'),pt(1,2)]);
end
function [x,y,trackHand,barHand] = getline_mb(varargin)
%GETLINE Draw a multibeam track with mouse.
%   [X,Y] = GETLINE(FIG) lets you plan a multibeam track in the current axes of figure FIG
%   using the mouse.  Coordinates of the track are returned in X and Y.  Use normal
%   button clicks to add points to the track.  A shift-, or double-click adds a
%   final point and ends the track selection. Pressing RETURN or ENTER ends the track
%   selection without adding a final point. Pressing BACKSPACEor DELETE removes the previously
%   selected point from the track.
%
%   Callback syntaxes:
%        getline_mb('KeyPress')
%        getline_mb('FirstButtonDown')
%        getline_mb('NextButtonDown')
%        getline_mb('ButtonMotion')

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

if (getappdata(gcbf,'haveUserdata'))    % Seek to a pre-existing figure's Userdata
    ud_old = get(gcbf,'Userdata');
    setappdata(gcbf,'haveUserdata',0)   % temporarly set it to 0
end

xlimorigmode = get(gca,'xlimmode');    ylimorigmode = get(gca,'ylimmode');
set(gca,'xlimmode','manual');   set(gca,'ylimmode','manual');

if ((length(varargin) >= 1) && ischar(varargin{1}))
    % Callback invocation: 'KeyPress', 'FirstButtonDown', 'NextButtonDown', or 'ButtonMotion'.
    feval(varargin{:});    return;
end

ud.GETLINE_X = [];     ud.GETLINE_Y = [];       ud.GETSEGLINE = [];       ud.GETCIRC = [];

if (length(varargin) < 1)
    ud.GETLINE_AX = gca;
    ud.GETLINE_FIG = get(ud.GETLINE_AX, 'Parent');
else
    if (~ishandle(varargin{1}))
        CleanUp(xlimorigmode,ylimorigmode);
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
        CleanUp(xlimorigmode,ylimorigmode);
        error('First argument should be a figure or axes handle');
    end
end

% Remember initial figure state
old_db = get(ud.GETLINE_FIG, 'DoubleBuffer');
state= uisuspend_j(ud.GETLINE_FIG);

% Set up initial callbacks for initial stage
set(ud.GETLINE_FIG, 'Pointer', 'crosshair', ...
    'WindowButtonDownFcn', 'getline_mb(''FirstButtonDown'');',...
    'KeyPressFcn', 'getline_mb(''KeyPress'');', 'DoubleBuffer', 'on');

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

set(gcbf, 'UserData', ud);

% We're ready; wait for the user to do the drag. Wrap the call to waitfor
% in try-catch so we'll have a chance to clean up after ourselves.
errCatch = 0;
try         waitfor(ud.GETLINE_H1, 'UserData', 'Completed');
catch,		errCatch = 1;
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
    % Recover the data that meawhile (and God knows why) was lost in nature
    ud_tmp = getappdata(gcf,'vertices');    try    rmappdata(gcf,'vertices');  end
    if ~isempty(ud_tmp)         % This is basicly an error testing. It happens (for example)
        ud = ud_tmp;            % when the user right-clicks on the first line point.
    end                         % Which stops the line drawing imediatly.
	try		x = ud.GETLINE_X(:);    y = ud.GETLINE_Y(:);
	catch,	x = [];                 y = [];
	end
    % If no points were selected, return rectangular empties.
    % This makes it easier to handle degenerate cases in functions that call getline_mb.
    if (isempty(x));        x = zeros(0,1);    end
    if (isempty(y));        y = zeros(0,1);    end
end

% Delete the animation objects
if (ishandle(ud.GETLINE_H1));    delete(ud.GETLINE_H1);     end
if (ishandle(ud.GETLINE_H2));    delete(ud.GETLINE_H2);     end

% Restore the figure's initial state
if (ishandle(ud.GETLINE_FIG))
   uirestore_j(state, 'nochildren');
   set(ud.GETLINE_FIG, 'DoubleBuffer', old_db);
end

CleanUp(xlimorigmode,ylimorigmode);
% Depending on the error status, return the answer or generate an error message.
switch errStatus
case 'ok'      % Return the answer
    if (nargout == 3 || nargout == 4)
        trackHand = ud.GETSEGLINE;
        barHand = ud.GETCIRC;
    end
case 'trap'    % An error was trapped during the waitfor
    %error('Interruption during mouse selection.');
    x = [];     y = [];    trackHand = [];   barHand = [];
case 'unknown'
    % User did something to cause the polyline selection to
    % terminate abnormally.  For example, we would get here
    % if the user closed the figure in the middle of the selection.
    %error('Interruption during mouse selection.');
    x = [];     y = [];    trackHand = [];   barHand = [];
end

set(gcbf,'Userdata',[])
try     % I have to use "try" because ud_old might not exist, and that make an error below
    if ~isempty(ud_old)         % If a previous Userdata existed, set it back
        set(gcbf,'Userdata',ud_old)
        setappdata(gcbf,'haveUserdata',1)
    end
end

%--------------------------------------------------
% Subfunction KeyPress
%--------------------------------------------------
function KeyPress
ud = get(gcbf, 'UserData');

key = get(ud.GETLINE_FIG, 'CurrentCharacter');
switch key
    case {char(8), char(127)}  % delete and backspace keys
        % remove the previously selected point
        switch length(ud.GETLINE_X)
            case 0        % nothing to do
            case 1
                ud.GETLINE_X = [];        ud.GETLINE_Y = [];
                % remove point and start over
                set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
                delete(ud.GETCIRC(end));        ud.GETCIRC = [];
                set(ud.GETLINE_FIG, 'WindowButtonDownFcn', ...
                        'getline_mb(''FirstButtonDown'');', 'WindowButtonMotionFcn', '');
            otherwise        % remove last point
                ud.GETLINE_X(end) = [];         ud.GETLINE_Y(end) = [];
                delete(ud.GETSEGLINE(end));     ud.GETSEGLINE(end) = [];
                delete(ud.GETCIRC(end));        ud.GETCIRC(end) = [];
                set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
        end
    case {char(13), char(3)}   % enter and return keys
        % return control to line after waitfor
        set(ud.GETLINE_H1, 'UserData', 'Completed');     
        setappdata(gcf,'vertices',ud)      % Save the plolyline vertices
end
set(gcbf, 'UserData', ud);

%--------------------------------------------------
% Subfunction FirstButtonDown
%--------------------------------------------------
function FirstButtonDown
ud = get(gcbf, 'UserData');
[x,y,z] = getcurpt_eu(ud.GETLINE_AX);

% check if GETLINE_X,GETLINE_Y is inside of axis
x_lim = get(ud.GETLINE_AX,'xlim');      y_lim = get(ud.GETLINE_AX,'ylim');
if (x>=x_lim(1)) && (x<=x_lim(2)) && (y>=y_lim(1)) && (y<=y_lim(2))    % inside axis limits
    ud.GETLINE_X = x;    ud.GETLINE_Y = y;
else    % outside axis limits, ignore this FirstButtonDown
    return
end
if isnan(z);    return;  end

set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y, 'Visible', 'on');
sr = getappdata(ud.GETLINE_FIG,'swathRatio');
rad = abs(z) * (sr/2) / 111194.9; % metros -> degrees
[latc,lonc] = circ_geo(y,x,rad);
h_circ = line('XData', lonc, 'YData', latc,'Color','w','LineWidth',.5);
ud.GETCIRC = [ud.GETCIRC h_circ];  % save this circle handle

set(ud.GETLINE_FIG, 'WindowButtonMotionFcn', 'getline_mb(''ButtonMotion'');', ...
      'WindowButtonDownFcn', 'getline_mb(''NextButtonDown'');');
set(gcbf, 'UserData', ud);

%--------------------------------------------------
% Subfunction NextButtonDown
%--------------------------------------------------
function NextButtonDown
ud = get(gcbf, 'UserData');

selectionType = get(ud.GETLINE_FIG, 'SelectionType');
if (~strcmp(selectionType, 'open') && strcmp(selectionType, 'normal'))
    % We don't want to add a point on the second click of a double-click
    [x,y,z] = getcurpt_eu(ud.GETLINE_AX);

    % check if GETLINE_X,GETLINE_Y is inside of axis
    x_lim = get(ud.GETLINE_AX,'xlim');      y_lim = get(ud.GETLINE_AX,'ylim');
    if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2))    % outside axis limits, ignore this ButtonDown
        return
    end    
    if isnan(z);    return;  end
    
    ud.GETLINE_X = [ud.GETLINE_X x];        ud.GETLINE_Y = [ud.GETLINE_Y y];
    set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
    h_line = line('XData', ud.GETLINE_X(end-1:end), 'YData', ud.GETLINE_Y(end-1:end),'LineWidth',1);
    ud.GETSEGLINE = [ud.GETSEGLINE h_line];  % save this segment line handle
    sr = getappdata(ud.GETLINE_FIG,'swathRatio');
    rad = abs(z) * (sr/2) / 111194.9; % metros -> degrees
    az = azimuth_geo(ud.GETLINE_Y(end-1),ud.GETLINE_X(end-1),ud.GETLINE_Y(end),ud.GETLINE_X(end));
    [lat1,lon1] = circ_geo(y,x,rad,[az-90-1 az-90+1],3);
    [lat2,lon2] = circ_geo(y,x,rad,[az+90-1 az+90+1],3);
    h_circ = line('XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)],'Color',[.8 .8 .8],'LineWidth',4);
    %h_circ = line('XData', lonc, 'YData', latc,'Color','w','LineWidth',.5);
    ud.GETCIRC = [ud.GETCIRC h_circ];  % save this circle handle
end

if (strcmp(get(ud.GETLINE_FIG, 'SelectionType'), 'alt'))    % Right-click, delete previous point
    switch length(ud.GETLINE_X)
        case 0        % nothing to do
        case 1
            ud.GETLINE_X = [];        ud.GETLINE_Y = [];
            % remove point and start over
            set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
            delete(ud.GETCIRC(end));        ud.GETCIRC = [];
            set(ud.GETLINE_FIG, 'WindowButtonDownFcn', ...
                    'getline_mb(''FirstButtonDown'');', 'WindowButtonMotionFcn', '');
        otherwise        % remove last point
            ud.GETLINE_X(end) = [];         ud.GETLINE_Y(end) = [];
            delete(ud.GETSEGLINE(end));     ud.GETSEGLINE(end) = [];
            delete(ud.GETCIRC(end));        ud.GETCIRC(end) = [];
            set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', ud.GETLINE_X, 'YData', ud.GETLINE_Y);
    end
end

if (strcmp(selectionType, 'extend') || strcmp(selectionType, 'open'))    % We're done (midle button or double-click)
    set(ud.GETLINE_H1, 'UserData', 'Completed');
    setappdata(gcf,'vertices',ud)      % Save the plolyline vertices
end
set(gcbf, 'UserData', ud);

%-------------------------------------------------
% Subfunction ButtonMotion
%-------------------------------------------------
function ButtonMotion
	ud = get(gcbf, 'UserData');
	[newx, newy] = getcurpt_eu(ud.GETLINE_AX);
	x = [ud.GETLINE_X newx];    y = [ud.GETLINE_Y newy];
	set([ud.GETLINE_H1 ud.GETLINE_H2], 'XData', x, 'YData', y);

%---------------------------------------------------
% Subfunction CleanUp
%--------------------------------------------------
function CleanUp(xlimmode,ylimmode)
	set(gca,'xlimmode',xlimmode);   set(gca,'ylimmode',ylimmode);

%---------------------------------------------------
% Subfunction getcurpt_eu
%--------------------------------------------------
function [x,y,z] = getcurpt_eu(axHandle)
%GETCURPT Get current point.
	pt = get(axHandle, 'CurrentPoint');
	x = pt(1,1);    y = pt(1,2);
	if nargout == 3
		X = getappdata(gcf,'dem_x');    Y = getappdata(gcf,'dem_y');
		Z = getappdata(gcf,'dem_z');    z = abs(bi_linear(X,Y,Z,x,y));
	end

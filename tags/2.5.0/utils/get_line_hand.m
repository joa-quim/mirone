function h_line = get_line_hand(varargin)
%   H_LINE = GET_LINE_HAND(FIG) Selects a line/patch with a mouse double-click and return its handle
%
%   Double-click on the object to select it.
%   Alternatively a left-click selects the object and a right-click confirms its selection
%   A right-click stops the selection operation anyway. Either if an object was selected or not
%   In the latter case H_LINE = []
%
%   Example:
%       x = [1 2 3];
%       y = [1 2 2.5];
%       h = figure;
%       plot(x,y);
%       h_l = get_line_hand(h)
%
%   This function uses ideas from the GETLINE (ITBx) and GraphTools
%
%	Joaquim Luis
%

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

if (nargin ~= 1 || ~ishandle(varargin{1}) || ~strcmp(get(varargin{1},'Type'),'figure'))
    error('get_line_hand error: Input argument must be a figure handle');
end
ud.GETLINE_FIG = varargin{1};       ud.GETLINE_AX = get(ud.GETLINE_FIG,'CurrentAxes');
    
% Remember initial figure state
old_db = get(ud.GETLINE_FIG, 'DoubleBuffer');
state= uisuspend(ud.GETLINE_FIG);

% Set up initial callbacks for initial stage
set(ud.GETLINE_FIG, 'Pointer', 'crosshair', 'WindowButtonDownFcn', {@FirstButtonDown,ud.GETLINE_FIG}, ...
    'DoubleBuffer', 'on');

% Bring target figure forward
figure(ud.GETLINE_FIG);

ud.GETLINE_H1 = line('Parent', ud.GETLINE_AX, 'XData', [], 'YData', [], ...
                  'Tag', 'xxx', 'Visible', 'off');

setappdata(ud.GETLINE_FIG, 'FromGetLineHand', ud);

% We're ready; wait for the user to do the drag. Wrap the call to waitfor
% in try-catch so we'll have a chance to clean up after ourselves.
errCatch = 0;
try         waitfor(ud.GETLINE_H1, 'UserData', 'Completed');
catch       errCatch = 1;   end

% After the waitfor, if GETLINE_H1 is still valid and its UserData is 'Completed', then the user
% completed the drag.  If not, the user interrupted the action somehow, perhaps by a Ctrl-C in the
% command window or by closing the figure.

if (errCatch == 1)
    errStatus = 'trap';
elseif (~ishandle(ud.GETLINE_H1) || ~strcmp(get(ud.GETLINE_H1, 'UserData'), 'Completed'))
    errStatus = 'unknown';
else
    errStatus = 'ok';
    h_line = getappdata(ud.GETLINE_H1,'TrackHandle');
end

% Delete the store object
if (ishandle(ud.GETLINE_H1));    delete(ud.GETLINE_H1);     end

% Restore the figure's initial state
if (ishandle(ud.GETLINE_FIG))
   uirestore(state);
   set(ud.GETLINE_FIG, 'DoubleBuffer', old_db);
   try  rmappdata(ud.GETLINE_FIG,'FromGetLineHand');     end
end

% Depending on the error status, return the answer.
switch errStatus
	case 'ok'                   % Return the answer
	case {'trap' 'unknown'}     % An error was trapped during the waitfor
        h_line = [];
end

%---------------------------------------------------------------------------------------
function FirstButtonDown(obj,eventdata,hfig)
ud = getappdata(hfig, 'FromGetLineHand');
selectionType = get(hfig, 'SelectionType');
% I have to do this test here because there is another mouse click inside get_lineHandle
if (strcmp(selectionType,'alt')) || (strcmp(selectionType,'extend')) || (strcmp(selectionType,'open'))
    % User changed his mind (right click) or ended selection
    set(ud.GETLINE_H1, 'UserData', 'Completed');
    try,    % If we have an (unknown) error in ud, another one would occur here
        if (ishandle(ud.markers)),        delete(ud.markers);     end
    end
else
    current_pt = get(ud.GETLINE_AX, 'CurrentPoint');
    ud.lineHandle = get_lineHandle(current_pt,hfig);
    setappdata(ud.GETLINE_H1,'TrackHandle',ud.lineHandle)
    if isempty(ud.lineHandle)       % User gave up (right or midle click) inside get_lineHandle
        set(ud.GETLINE_H1, 'UserData', 'Completed');
    else
        x = get(ud.lineHandle,'XData');   y = get(ud.lineHandle,'YData');
        hold on;    
        ud.markers = plot(x,y,'kp','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10); 
        hold off
        set(hfig, 'WindowButtonDownFcn', {@NextButtonDown,hfig});
        setappdata(ud.GETLINE_H1,'Xvert',x);    setappdata(ud.GETLINE_H1,'Yvert',y);
    end
end
setappdata(hfig, 'FromGetLineHand', ud);

%---------------------------------------------------------------------------------------
function NextButtonDown(obj,eventdata,hfig)
% Most of what this does is to check if it's finish and if not, call back FirstButtonDown
ud = getappdata(hfig, 'FromGetLineHand');
try     selectionType = get(ud.GETLINE_FIG, 'SelectionType'); 
catch   % Open an exit door, otherwise we are stucked inside whatever caused the error.
    set(ud.GETLINE_H1, 'UserData', 'Completed');
end

if (strcmp(selectionType,'alt')) || (strcmp(selectionType,'extend')) || (strcmp(selectionType,'open'))
    % User changed his mind (right click) or ended selection
    set(ud.GETLINE_H1, 'UserData', 'Completed');
    if (ishandle(ud.markers)),        delete(ud.markers);     end
elseif (strcmp(selectionType, 'normal'))    % left-click
    pt = get(ud.GETLINE_AX, 'CurrentPoint');
    x = pt(1,1);    y = pt(1,2);
    % check if x,y is inside of axis
    x_lim = get(ud.GETLINE_AX,'xlim');      y_lim = get(ud.GETLINE_AX,'ylim');
    if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2))    % outside axis limits, ignore this ButtonDown
        return
    end    
    set(hfig, 'WindowButtonDownFcn', {@FirstButtonDown,hfig}); 
end
setappdata(hfig, 'FromGetLineHand', ud);

%---------------------------------------------------------------------
function lineHandle = get_lineHandle(pt,hfig)
% lineHandle = get_lineHandle(pt,hfig)
% This subfunction and the next are adapted from GraphTools
%---------------------------------------------------------------------
ii = prop_list(hfig,'axes');     h_lines = [];
for (i = 1:length(ii))
    axes(ii(i))
    h_lines = [h_lines; prop_list(hfig,'line',1)];          % try to find lines handles
    if (isempty(h_lines))                                   % Nope, no line
        h_lines = [h_lines; prop_list(hfig,'patch',1)];     % try find patch handles
    end
end

h_lines = h_lines(h_lines ~= 0);

key = 0;    ii = [];    first = 1;  button = 1;    lineHandle = [];
while (key == 0)
    if (first)
        x = pt(1,1);    y = pt(1,2);    first = 0;
    else
        try                                 % We use a try wrap to make it possible to this function
            [x,y,button] = ginput(1);       % to work in R13. Otherwise, if the user had sight troubles
        catch                               % and did not click on a object, uisuspend limitations would
            button = 1;                     % produce an error here if.
            continue
        end
    end
    if (button ~= 1)
        lineHandle = [];    return;
    end
    for (i = 1:length(h_lines)),    ii=[ii,i];      end
    hh = h_lines(ii);
    if (length(hh) > 0)
        if any(hh == gco)
            lineHandle = gco;
            key = 1;
        end
    else
        warndlg('Sorry, there are no lines to be edited!','Warning!')
        lineHandle = [];    return; 
    end   
end

%---------------------------------------------------------------------------
function key = prop_list(hfig,Type,arg1)
%   key = prop_list(hfig,Type,arg1) 
%   arg1 serves only to increase narargin and force to search for axes childs
%   return the handles of the specified property in the current window.
%---------------------------------------------------------------------------
if (nargin == 2)
    y = get(hfig,'Child');
else        % nargin == 3, but not tested
    h_ax = findobj(hfig,'Type','axes');
    y = get(h_ax,'Child');
    store = findobj(y,'Tag','xxx');     % This is the "store" object (to hold the markers)
    y = setxor(y, store);               % It doesn't count, so remove it from the child list
end
ii = [];
for (i = 1:length(y))
   c = get(y(i),'Type');
   if (strcmp(c,Type)), ii = [ii,i];    end
end   
key = y(ii);

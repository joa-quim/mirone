function h_line = get_polygon(varargin)
%   H_LINE = GET_POLYGON(FIG) Selects a line or patch with a mouse left-click and return its handle
%
%   A left-click selects the object. A right-click confirms its selection
%   Alternatively double-click on the object.

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

if (~ishandle(varargin{1}) || ~strcmp(get(varargin{1},'Type'),'figure'))
    errordlg('get_polygon error: First input argument must be a figure handle','Error');    return;
end
ud.GETLINE_FIG = varargin{1};       ud.GETLINE_AX = get(ud.GETLINE_FIG,'CurrentAxes');

if (nargin == 2 && strcmp(varargin{2}, 'multi'))
	multi_selec = true;
else
	multi_selec = false;
end

% Remember initial figure state
old_db = get(ud.GETLINE_FIG, 'DoubleBuffer');
state = uisuspend_j(ud.GETLINE_FIG);

% Set up initial callbacks for initial stage
set(ud.GETLINE_FIG, 'Pointer', 'crosshair', 'WindowButtonDownFcn', {@FirstButtonDown,ud.GETLINE_FIG}, ...
    'DoubleBuffer', 'on');

% Bring target figure forward
figure(ud.GETLINE_FIG);

% This had a totaly different use in the original getline, but had reveled very usefull
% here for it's handle is used to store information that if stored in ud.lineHandle gave
% a very strange bug on a second call of get_polygon (e.g. after clicking the right or
% midle button, it didn't work anymore because ud structure was unknown (???))
ud.GETLINE_H1 = line('Parent', ud.GETLINE_AX, 'XData', [], 'YData', [], 'Tag', 'xxxxxxx', ...
					'Visible', 'off');

ud.multi_selec = multi_selec;
setappdata(ud.GETLINE_FIG, 'FromGetPolygon', ud);

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
    h_line = getappdata(ud.GETLINE_H1,'TrackHandle');
end

% Restore the probable 'buttondownfcn' set by ui_edit_polygon
bdfcn = getappdata(ud.GETLINE_H1,'Eventbdfcn');

% Delete the store object
if (ishandle(ud.GETLINE_H1));    delete(ud.GETLINE_H1);     end

% If by an (user) error we still have the markers, kill them
ud = getappdata(ud.GETLINE_FIG, 'FromGetPolygon');
try
    if (ishandle(ud.markers)),        delete(ud.markers);     end
end
if (multi_selec)
	delete(findobj(ud.GETLINE_AX,'type','line','Tag','StarMarkers'))
end

% Restore the figure's initial state
if (ishandle(ud.GETLINE_FIG))
   uirestore_j(state, 'nochildren');
   set(ud.GETLINE_FIG, 'DoubleBuffer', old_db);
   try, rmappdata(ud.GETLINE_FIG,'FromGetPolygon');     end
end

% Depending on the error status, return the answer.
switch errStatus
	case 'ok'                   % Return the answer
        set(h_line,'buttondownfcn',bdfcn)
	case {'trap' 'unknown'}     % An error was trapped during the waitfor
        h_line = [];
end

%---------------------------------------------------------------------------------------
function FirstButtonDown(obj,eventdata,hfig)
% If MULTI_SELECT select multiple lines. Otherwise, keep only the last selected one
	ud = getappdata(hfig, 'FromGetPolygon');
	selectionType = get(hfig, 'SelectionType');
	% I have to do this test here because there is another mouse click inside get_trackHandle
	if (strcmp(selectionType,'alt')) || (strcmp(selectionType,'extend')) || (strcmp(selectionType,'open'))
        % User changed his mind (right click) or ended selection
        set(ud.GETLINE_H1, 'UserData', 'Completed');
        try				% If we have an (unknown) error in ud, another one would occur here
            if (ishandle(ud.markers)),        delete(ud.markers);     end
        end
	else
        current_pt = get(ud.GETLINE_AX, 'CurrentPoint');
        [ud.lineHandle ud.h_circ] = get_trackHandle(current_pt);
		if (~ud.multi_selec)		% Single line selection
        	setappdata(ud.GETLINE_H1,'TrackHandle',ud.lineHandle)
		else						% Multiple line selection (Ignores the h_circ)
			tmp = getappdata(ud.GETLINE_H1,'TrackHandle');
			setappdata(ud.GETLINE_H1,'TrackHandle',[tmp; ud.lineHandle])
		end
        setappdata(ud.GETLINE_H1,'BarHandle',ud.h_circ)
        if isempty(ud.lineHandle)       % User gave up (right or midle click) inside get_trackHandle
            set(ud.GETLINE_H1, 'UserData', 'Completed');
        else
            x = get(ud.lineHandle,'XData');   y = get(ud.lineHandle,'YData');
            ud.markers = line('XData',x ,'YData',y, 'Parent',ud.GETLINE_AX, 'LineStyle','none','Marker','p', ...
				'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10, 'HitTest','off','Tag','StarMarkers'); 
			if (~ud.multi_selec)		% I wonder if we realy ever need the next line case
				set(hfig, 'WindowButtonDownFcn', {@NextButtonDown,hfig});
			end
            setappdata(ud.GETLINE_H1,'Xvert',x);    setappdata(ud.GETLINE_H1,'Yvert',y);
            % Backup and remove 'buttondownfcn' set by ui_edit_polygon. It will be reset by main function
            bdfcn = get(ud.lineHandle,'buttondownfcn');
            setappdata(ud.GETLINE_H1,'Eventbdfcn',bdfcn)
            set(ud.lineHandle,'buttondownfcn','')
        end
	end
	setappdata(hfig, 'FromGetPolygon', ud);

%---------------------------------------------------------------------------------------
function NextButtonDown(obj,eventdata,hfig)
% Most of what this does is to check if finish and if not call back FirstButtonDown
	ud = getappdata(hfig, 'FromGetPolygon');
	try,    selectionType = get(ud.GETLINE_FIG, 'SelectionType'); 
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
	setappdata(hfig, 'FromGetPolygon', ud);

%---------------------------------------------------------------------
function [lineHandle,barHandles,button] = get_trackHandle(pt)
% get_trackHandle function is used to get the handle of the multi-beam track object.
% [lineHandle,barHandles,button]=get_trackHandle(pt)
% This subfunction and the next are based from GraphTools
%---------------------------------------------------------------------
	ii = prop_list('axes');		h_lines=[];	    lineHandle = [];     barHandles = [];
	for i=1:length(ii)
       axes(ii(i))
       h_lines=[h_lines; prop_list('line',1)];		% find the lines handles
       h_lines=[h_lines; prop_list('patch',1)];		% find the patch handles
	end
	if (numel(h_lines) == 0)
		warndlg('Sorry, but there are no lines to be selected!','Warning!')
		lineHandle = [];		barHandles = [];	return 
	end
	
	% Get rid of the swath width line handles. For that they have to have a 'swath_w#' Tag
	% where # stands for the track number
	tmp1 = h_lines;  tmp2 = zeros(1,numel(h_lines));
	for i = 1:numel(h_lines)
		tag = get(h_lines(i),'Tag');
		if (strncmp(tag, 'swath_w', 7))
			tmp1(i) = 0;	tmp2(i) = h_lines(i);      % get the bar (or circle) handles
		end
	end
	h_lines = h_lines(tmp1 ~= 0);			% get only the tracks handles
	ALLbarHandles = tmp2(tmp2 ~= 0);		% get the handles to ALL bar handles (e.g. all bars of all tracks)
	
	key = 0;	first = 1;		button = 1;
	while (key == 0)
		if (first),		first = 0;
		else            [x,y,button] = click_e_point(1,'crosshair');
		end
		if (button ~= 1),	lineHandle=[];	return,		end
		if (any(h_lines == gco)),	lineHandle = gco;		key=1;		end
	end
	
	% Now find the handles to the bar lines which correspond to the swath width
	% at each vertex of the selected lineHandle
	tag = get(lineHandle,'Tag');
	if (strncmp(tag, 'swath_w', 7))
		nTrack = tag(8:end);        % get the string with the track number
		barHandles = findobj(ALLbarHandles,'Tag',['swath_w' nTrack]);
		barHandles = sort(barHandles);
	end

%---------------------------------------------------------------------------
function key = prop_list(Type,arg1)
%  key=prop_list(Type,arg1) 
%  return the handles of the specified property in the current window.
%---------------------------------------------------------------------------
	if (nargin == 1),	y = get(gcf,'Child');
	else				y = get(gca,'Child');
	end
	ii = false(numel(y),1);
	for i = 1:numel(y)
		c = get(y(i),'Type');
		if strcmp(c,Type),	ii(i) = true;	end
	end   
	key = y(ii);

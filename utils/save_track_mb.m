function [x,y] = save_track_mb(varargin)
%   Selects a multibeam track with mouse and return its vertices or
%   delete it, if varargin = 1.
%   compute it's length in NM, if varargin = 2.
%   A left-click selects the MB track. Any other click confirms either the saving,
%   deleting or compute it's length.
  
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

xlimorigmode = get(gca,'xlimmode');     ylimorigmode = get(gca,'ylimmode');
set(gca,'xlimmode','manual');           set(gca,'ylimmode','manual');

if ((length(varargin) >= 1) && ischar(varargin{1}))
    % Callback invocation: 'FirstButtonDown', 'NextButtonDown', or 'ButtonMotion'.
    feval(varargin{:});    return;
end

ud.GETLINE_X = [];      ud.GETLINE_Y = [];       ud.GETSEGLINE = [];       ud.GETCIRC = [];
ud.GETLINE_AX = gca;    ud.GETLINE_FIG = get(ud.GETLINE_AX, 'Parent');

if (nargin == 0 && nargout ~= 2)
    errordlg('Error in save track. Must give two ouput variables (x,y).','Error');    return;
elseif (nargin == 1 && varargin{1} == 2 && nargout == 0)
    errordlg('Error in computing track length . Must give one ouput variable.','Error');    return;
end
    
% Remember initial figure state
old_db = get(ud.GETLINE_FIG, 'DoubleBuffer');
state = uisuspend_j(ud.GETLINE_FIG);

% Set up initial callbacks for initial stage
set(ud.GETLINE_FIG, 'Pointer', 'crosshair', 'WindowButtonDownFcn', 'save_track_mb(''FirstButtonDown'');', ...
    'KeyPressFcn', 'save_track_mb(''KeyPress'');', 'DoubleBuffer', 'on');

% Bring target figure forward
figure(ud.GETLINE_FIG);

% This had a totaly different use in the original getline, but had reveled very usefull
% here for it's handle is used to store information that if stored in ud.lineHandle gave
% a very strange bug on a second call of save_track_mb (e.g. after clicking the right or
% midle button, it didn't work anymore because ud structure was unknown (???))
ud.GETLINE_H1 = line('Parent', ud.GETLINE_AX, 'XData', [], 'YData', [], 'Tag', 'xxxxxxx', ...
                  'Visible', 'off');

set(gcbf, 'UserData', ud);

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
    if (nargin == 1 && varargin{1} == 1)             % Delete track
        kill_t = getappdata(ud.GETLINE_H1,'TrackHandle');
        if ishandle(kill_t);  delete(kill_t);   end
        kill_b = getappdata(ud.GETLINE_H1,'BarHandle');
        if ishandle(kill_b);  delete(kill_b);   end
        x = [];     y = [];         % Although there are no outputs in this case, but ML is very picky.
    elseif (nargin == 1 && varargin{1} == 2)         % compute track length
        try    lon = getappdata(ud.GETLINE_H1,'Xvert');   lat = getappdata(ud.GETLINE_H1,'Yvert');
            D2R = pi/180;
            earth_rad = 6371;
            %earth_rad = 6366.7;     % NOTE: this radius (in km) is the one that leads to the deffinition of a
                                    % NM = 1.852 km. It is obtained by rad = 1.852*60*360/(2*pi)
            lon = lon * D2R;    lat = lat * D2R;
            lat_i = lat(1:length(lat)-1);   lat_f = lat(2:length(lat));
            lon_i = lon(1:length(lon)-1);   lon_f = lon(2:length(lon));
            tmp = sin(lat_i).*sin(lat_f) + cos(lat_i).*cos(lat_f).*cos(lon_f-lon_i);
            x = sum(acos(tmp) * earth_rad / 1.852);     % Distance in NM
		catch,	x = [];          % When user gave up or an error (?) occured
        end
    else                                            % Save track
        try		x = getappdata(ud.GETLINE_H1,'Xvert');   y = getappdata(ud.GETLINE_H1,'Yvert');
		catch,	x = [];     y = [];     end     % When user gave up or an error (?) occured
    end
end

% Delete the store object
if (ishandle(ud.GETLINE_H1));    delete(ud.GETLINE_H1);     end

% Restore the figure's initial state
if (ishandle(ud.GETLINE_FIG))
   uirestore_j(state, 'nochildren');
   set(ud.GETLINE_FIG, 'DoubleBuffer', old_db);
end

CleanUp(xlimorigmode,ylimorigmode);
% Depending on the error status, return the answer.
switch errStatus
case 'ok'       % Return the answer
case 'trap'     % An error was trapped during the waitfor
    x = [];     y = [];
case 'unknown'  % User did something to cause the polyline selection to terminate abnormally.
    x = [];     y = [];
end

%--------------------------------------------------
% Subfunction KeyPress
%--------------------------------------------------
function KeyPress
ud = get(gcbf, 'UserData');
key = get(ud.GETLINE_FIG, 'CurrentCharacter');
switch key
    case {char(13), char(3)}   % enter and return keys
        % return control to line after waitfor
        set(ud.GETLINE_H1, 'UserData', 'Completed');
        try    % If we have an (unknown) error in ud, another one would occur here
            if (ishandle(ud.markers)),        delete(ud.markers);     end
        end
end
set(gcbf, 'UserData', ud);

%--------------------------------------------------
% Subfunction FirstButtonDown
%--------------------------------------------------
function FirstButtonDown
ud = get(gcbf, 'UserData');
selectionType = get(ud.GETLINE_FIG, 'SelectionType');
% I have to do this test here because there is another mouse click inside get_trackHandle
if (strcmp(selectionType,'alt')) || (strcmp(selectionType,'extend')) || (strcmp(selectionType,'open'))
    % User changed his mind (right click) or ended selection
    set(ud.GETLINE_H1, 'UserData', 'Completed');
    try    % If we have an (unknown) error in ud, another one would occur here
        if (ishandle(ud.markers)),        delete(ud.markers);     end
    end
else
    current_pt = get(ud.GETLINE_AX, 'CurrentPoint');
    [ud.lineHandle ud.h_circ] = get_trackHandle(current_pt);
    setappdata(ud.GETLINE_H1,'TrackHandle',ud.lineHandle)
    setappdata(ud.GETLINE_H1,'BarHandle',ud.h_circ)
    if isempty(ud.lineHandle)       % User gave up (right or midle click) inside get_trackHandle
        set(ud.GETLINE_H1, 'UserData', 'Completed');
    else
        x = get(ud.lineHandle,'XData');   y = get(ud.lineHandle,'YData');
        hold on;    
        ud.markers = plot(x,y,'kp','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10); 
        hold off
        set(ud.GETLINE_FIG, 'WindowButtonDownFcn', 'save_track_mb(''NextButtonDown'');');
        setappdata(ud.GETLINE_H1,'Xvert',x);    setappdata(ud.GETLINE_H1,'Yvert',y);
    end
end

set(gcbf, 'UserData', ud);

%--------------------------------------------------
% Subfunction NextButtonDown
%--------------------------------------------------
function NextButtonDown
% Most of waht this does is to check if finish and if not call back FirstButtonDown
ud = get(gcbf, 'UserData');
try    selectionType = get(ud.GETLINE_FIG, 'SelectionType'); 
catch   % Open an exit door, otherwise we are stucked inside whatever caused the error.
    set(ud.GETLINE_H1, 'UserData', 'Completed');
end

if (strcmp(selectionType,'alt')) || (strcmp(selectionType,'extend')) || (strcmp(selectionType,'open'))
    % User changed his mind (right click) or ended selection
    set(ud.GETLINE_H1, 'UserData', 'Completed');
    if (ishandle(ud.markers)),        delete(ud.markers);     end
elseif (strcmp(selectionType, 'normal'))    % left-click
    [x,y] = getcurpt(ud.GETLINE_AX);
    % check if x,y is inside of axis
    x_lim = get(ud.GETLINE_AX,'xlim');      y_lim = get(ud.GETLINE_AX,'ylim');
    if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2))    % outside axis limits, ignore this ButtonDown
        return
    end    
    set(ud.GETLINE_FIG, 'WindowButtonDownFcn', 'save_track_mb(''FirstButtonDown'');');        
end
set(gcbf, 'UserData', ud);

%---------------------------------------------------
% Subfunction CleanUp
%--------------------------------------------------
function CleanUp(xlimmode,ylimmode)
set(gca,'xlimmode',xlimmode);   set(gca,'ylimmode',ylimmode);

%---------------------------------------------------
% Subfunction getcurpt
%--------------------------------------------------
function [x,y] = getcurpt(axHandle)
%GETCURPT Get current point.
pt = get(axHandle, 'CurrentPoint');
x = pt(1,1);    y = pt(1,2);

%---------------------------------------------------------------------
function [lineHandle,barHandles,button]=get_trackHandle(pt)
%get_trackHandle function is used to get the handle of the multi-beam track object.
%[lineHandle,barHandles,button]=get_trackHandle(pt)
% This subfunction and the next are based from GraphTools
%---------------------------------------------------------------------
ii = prop_list('axes');     h_lines=[];
for i=1:length(ii)
   axes(ii(i))
   h_lines=[h_lines; prop_list('line',1)]; % find the lines handles
end

% Get rid of the swath width line handles. For that they have to have a 'swath_w#' Tag
% where # stands for the track number
tmp1 = h_lines;  tmp2 = zeros(1,length(h_lines));
for i = 1:length(h_lines)
    tag = get(h_lines(i),'Tag');
    if length(tag) < 7 || isempty(tag),     continue;   end         % prevents errors from line elements with tags shorter than 6 char 
    if strcmp(tag(1:7),'swath_w')
        tmp1(i) = 0;     tmp2(i) = h_lines(i);      % get the bar (or circle) handles
    end
end
h_lines = h_lines(find(tmp1 ~= 0));         % get only the tracks handles
ALLbarHandles = tmp2(find(tmp2 ~= 0));         % get the handles to ALL bar handles (e.g. all bars of all tracks)

key=0;      ii=[];     first = 1;  button = 1;    lineHandle = [];     barHandles = [];
while key==0
    if (first),     x = pt(1,1);    y = pt(1,2);    first = 0;
    else           [x,y,button] = click_e_point(1,'crosshair');
    end
    if button~=1, lineHandle=[];    return;     end
    for i=1:length(h_lines)
        ii=[ii,i];
    end
    hh=h_lines(ii);
    if length(hh) > 0
        if any(hh == gco), lineHandle=gco;  key=1;         end
    else
        warndlg('Sorry, there are no tracks to be edited!','Warning!')
        lineHandle=[]; return; 
    end   
end

% Now find the handles to the bar lines which correspond to the swath width
% at each vertex of the selected lineHandle
tag = get(lineHandle,'Tag');
%if ~strcmp(tag(1:7),'MBtrack')
    %errordlg('Error in get_trackHandle: selected track has a wrong Tag','Error');   return
%end
nTrack = tag(8:end);        % get the string with the track number
barHandles = findobj(ALLbarHandles,'Tag',['swath_w' nTrack]);
barHandles = sort(barHandles);

%---------------------------------------------------------------------------
%  key=prop_list(Type,arg1) 
%  return the handles of the specified property in the current window.
%---------------------------------------------------------------------------
function key=prop_list(Type,arg1)
if (nargin == 1),   y = get(gcf,'Child');
else                   y = get(gca,'Child'); end
ii = [];
for (i = 1:length(y))
   c = get(y(i),'Type');
   if strcmp(c,Type), ii = [ii,i]; end
end   
key = y(ii);

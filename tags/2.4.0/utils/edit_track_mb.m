function edit_track_mb(varargin)
%GETLINE Select polyline with mouse.
%   [X,Y] = GETLINE(FIG) lets you select a polyline in the current axes of figure FIG
%   using the mouse.  Coordinates of the polyline are returned in X and Y.  Use normal
%   button clicks to add points to the polyline.  A shift-, right-, or double-click adds a
%   final point and ends the polyline selection. Pressing RETURN or ENTER ends the polyline
%   selection without adding a final point. Pressing BACKSPACEor DELETE removes the previously
%   selected point from the polyline.
%
%   [X,Y] = GETLINE(AX) lets you select a polyline in the axes specified by the handle AX.
%   [X,Y] = GETLINE is the same as [X,Y] = GETLINE(GCF).
%

%   Callback syntaxes:
%        edit_track_mb('KeyPress')
%        edit_track_mb('FirstButtonDown')
%        edit_track_mb('NextButtonDown')
%        edit_track_mb('ButtonMotion')

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
    % Callback invocation: 'KeyPress', 'FirstButtonDown', 'NextButtonDown', or 'ButtonMotion'.
    feval(varargin{:});    return;
end

ud.GETLINE_AX = gca;    ud.GETLINE_FIG = get(ud.GETLINE_AX, 'Parent');
ud.inserted = 0;              % = 1 when inserting vertices
   
% Remember initial figure state
old_db = get(ud.GETLINE_FIG, 'DoubleBuffer');
state= uisuspend_j(ud.GETLINE_FIG);

% Set up initial callbacks for initial stage
set(ud.GETLINE_FIG, 'Pointer', 'crosshair', ...
    'WindowButtonDownFcn', 'edit_track_mb(''FirstButtonDown'');',...
    'KeyPressFcn', 'edit_track_mb(''KeyPress'');', 'DoubleBuffer', 'on');

% Bring target figure forward
figure(ud.GETLINE_FIG);

% This had a totaly different use in the original getline, but had reveled very usefull
% here for it's handle is used to store information that if stored in ud.lineHandle gave
% a very strange bug on a second call of edit_track_mb (e.g. after clicking the right or
% midle button, it didn't work anymore because ud structure was unknown (???))
ud.GETLINE_H1 = line('Parent', ud.GETLINE_AX, 'XData', [], 'YData', [], 'Tag', 'xxxxxxx', ...
                  'Visible', 'off');

set(gcbf, 'UserData', ud);
setappdata(gcbf,'haveUserdata',1);  % Notify other eventual functions that may be called
                                    % while waitfor is active to not overwrite this fucntion's
                                    % userdata. This happens for example when, while in the
                                    % editing mode, the user decides to measure a distance.

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
end

% Delete the animation objects
if (ishandle(ud.GETLINE_H1));    delete(ud.GETLINE_H1);     end

% Restore the figure's initial state
if (ishandle(ud.GETLINE_FIG))
   uirestore_j(state, 'nochildren');
   set(ud.GETLINE_FIG, 'DoubleBuffer', old_db);
end

CleanUp(xlimorigmode,ylimorigmode);
% Depending on the error status, return the answer or generate an error message.
switch errStatus
case 'ok'      % Return the answer
case 'trap'    % An error was trapped during the waitfor
    %x = [];     y = [];
case 'unknown'  % User did something to cause the polyline selection to terminate abnormally.
    %x = [];     y = [];
end

set(gcbf,'Userdata',[])     % Don't need this Userdata anymore, so better remove it

%--------------------------------------------------
% Subfunction KeyPress
%--------------------------------------------------
function KeyPress
ud = get(gcbf, 'UserData');
key = get(ud.GETLINE_FIG, 'CurrentCharacter');
switch key
    case {char(8), char(127), char(114)}  % "r", delete and backspace keys
        % remove the previously selected point (there are alot of particular cases)
        xx = get(ud.lineHandle,'XData');     yy = get(ud.lineHandle,'YData');
        xx(ud.vert_index) = [];              yy(ud.vert_index) = [];
        delete(ud.h_circ(ud.vert_index));    ud.h_circ(ud.vert_index) = [];
        set(ud.lineHandle,'XData',xx,'YData',yy);
        % now actualize the h_circ userdata that contains the vertex order
        for i = ud.vert_index:length(xx),     set(ud.h_circ(i),'Userdata',i);     end
        sr = getappdata(ud.lineHandle,'swathRatio');
        if (ud.vert_index > 1 && ud.vert_index < length(xx))
            az = azimuth_geo(yy(ud.vert_index-1),xx(ud.vert_index-1),yy(ud.vert_index),xx(ud.vert_index));
            [x,y,z] = getcurpt_eu(ud.GETLINE_AX,[xx(ud.vert_index) yy(ud.vert_index)]);
            rad = abs(z) * (sr/2) / 111194.9;
            [lat1,lon1] = circ_geo(yy(ud.vert_index),xx(ud.vert_index),rad,[az-90-1 az-90+1],3);
            [lat2,lon2] = circ_geo(yy(ud.vert_index),xx(ud.vert_index),rad,[az+90-1 az+90+1],3);
            set(ud.h_circ(ud.vert_index),'XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)]);
        elseif (ud.vert_index == 1 && length(xx) > 2)
            az = azimuth_geo(yy(ud.vert_index),xx(ud.vert_index),yy(ud.vert_index+1),xx(ud.vert_index+1));
            [x,y,z] = getcurpt_eu(ud.GETLINE_AX,[xx(ud.vert_index) yy(ud.vert_index)]);
            rad = abs(z) * (sr/2) / 111194.9;
            [lat1,lon1] = circ_geo(yy(ud.vert_index),xx(ud.vert_index),rad,[az-90-1 az-90+1],3);
            [lat2,lon2] = circ_geo(yy(ud.vert_index),xx(ud.vert_index),rad,[az+90-1 az+90+1],3);
            set(ud.h_circ(ud.vert_index),'XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)]);
            % now thirth vertex
            [x,y,z] = getcurpt_eu(ud.GETLINE_AX,[xx(ud.vert_index+1) yy(ud.vert_index+1)]);
            rad = abs(z) * (sr/2) / 111194.9;
            [lat1,lon1] = circ_geo(yy(ud.vert_index+1),xx(ud.vert_index+1),rad,[az-90-1 az-90+1],3);
            [lat2,lon2] = circ_geo(yy(ud.vert_index+1),xx(ud.vert_index+1),rad,[az+90-1 az+90+1],3);
            set(ud.h_circ(ud.vert_index+1),'XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)]);            
        end
        if (ud.vert_index == 2 && length(xx) > 2)     %Also in this case i-1 (that is i=1) bar has to be adjusted
            az = azimuth_geo(yy(1),xx(1),yy(ud.vert_index),xx(ud.vert_index));
            [x,y,z] = getcurpt_eu(ud.GETLINE_AX,[xx(1) yy(1)]);
            rad = abs(z) * (sr/2) / 111194.9;
            [lat1,lon1] = circ_geo(yy(1),xx(1),rad,[az-90-1 az-90+1],3);
            [lat2,lon2] = circ_geo(yy(1),xx(1),rad,[az+90-1 az+90+1],3);
            set(ud.h_circ(1),'XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)]);
        end
        set(ud.GETLINE_FIG, 'WindowButtonDownFcn', 'edit_track_mb(''FirstButtonDown'');', 'WindowButtonMotionFcn', '');
    case char(105)              % "i" char key (insert vertice)
        [newx, newy] = getcurpt_eu(ud.GETLINE_AX);
        xx = get(ud.lineHandle,'XData');    yy = get(ud.lineHandle,'YData');
        xx(ud.vert_index) = ud.save_x;      yy(ud.vert_index) = ud.save_y;
        if (ud.vert_index < length(xx))
            x = [xx(1:ud.vert_index) newx xx(ud.vert_index+1:end)];
            y = [yy(1:ud.vert_index) newy yy(ud.vert_index+1:end)];
        else        % New vertice is also the last vertice
            x = [xx(1:ud.vert_index) newx];
            y = [yy(1:ud.vert_index) newy];
        end
        set(ud.lineHandle,'XData',x,'YData',y);
        ud.inserted =  1;              % inform of an inserted vertex
        set(ud.h_circ(ud.vert_index),'Visible','on');
        ud.vert_index = ud.vert_index + 1;
    case char(98)              % "b" char key (break line)
        xx = get(ud.lineHandle,'XData');    yy = get(ud.lineHandle,'YData');
        xx(ud.vert_index) = ud.save_x;      yy(ud.vert_index) = ud.save_y;
        if (ud.vert_index > 1 && ud.vert_index < length(xx))
            x = xx(1:ud.vert_index);            y = yy(1:ud.vert_index);
            set(ud.lineHandle,'XData',x,'YData',y);
            set(ud.h_circ(ud.vert_index),'Visible','on');
            % Now make the a new track from rest of the original one, and not forget the tags.
            tag = get(ud.lineHandle,'Tag');
            tagL = [tag(1:7) num2str(str2double(tag(8:end)) + 5000)];      % adding 5000 should be enough
            tagB = ['swath_w' num2str(str2double(tag(8:end)) + 5000)];
            uictxm = get(ud.lineHandle,'uicontextmenu');
            tmp = line('XData',[], 'YData',[],'LineWidth',1,'Tag',tagL);  % create a new line handle
            set(tmp, 'XData',xx(ud.vert_index+1:end), 'YData',yy(ud.vert_index+1:end),'uicontextmenu',uictxm)
            set(ud.h_circ(ud.vert_index+1:end),'Tag',tagB)
            xxx = get(tmp,'XData');
            for i=1:length(xxx),    set(ud.h_circ(ud.vert_index+i),'UserData',i);     end
            set(ud.GETLINE_H1, 'UserData', 'Completed');
        end
    case {char(13), char(3)}   % enter and return keys
        % return control to line after waitfor
        set(ud.GETLINE_H1, 'UserData', 'Completed'); 
end
set(gcbf, 'UserData', ud);

%--------------------------------------------------
% Subfunction FirstButtonDown
%--------------------------------------------------
function FirstButtonDown
ud = get(gcbf, 'UserData');
% I have to do this test here because there is another mouse click inside get_trackHandle
if (~strcmp(get(ud.GETLINE_FIG, 'SelectionType'), 'normal'))    % We're done!
    set(ud.GETLINE_H1, 'UserData', 'Completed');
else
    current_pt = get(ud.GETLINE_AX, 'CurrentPoint');
    [ud.lineHandle ud.h_circ] = get_trackHandle(current_pt);
    if isempty(ud.lineHandle)       % User gave up (right or midle click) inside get_trackHandle
        set(ud.GETLINE_H1, 'UserData', 'Completed');
    else
        x = get(ud.lineHandle,'XData');   y = get(ud.lineHandle,'YData');
        if ~isempty(ud.lineHandle)
            % Find out which polyline vertice is closest to "current_pt".
            dif_x = x - current_pt(1,1);    dif_y = y - current_pt(1,2);
            dist = sqrt(dif_x.^2 + dif_y.^2);
            [B,IX] = sort(dist);    ud.vert_index = IX(1);      % tenho que ser mais esperto (potenc. lento)
            ud.save_x = x(ud.vert_index);   ud.save_y = y(ud.vert_index);   % needed in the "i"nsert option
        end
        % I'll use a try statement to prevent an eventual error in get_trackHandle
        % This happens when a regular line (not a MB track) was clicked
        try		set(ud.h_circ(ud.vert_index),'Visible','off'),		end
        set(ud.GETLINE_FIG, 'WindowButtonMotionFcn', 'edit_track_mb(''ButtonMotion'');', ...
            'WindowButtonDownFcn', 'edit_track_mb(''NextButtonDown'');');
    end
end
set(gcbf, 'UserData', ud);

%--------------------------------------------------
% Subfunction NextButtonDown
%--------------------------------------------------
function NextButtonDown
ud = get(gcbf, 'UserData');

selectionType = get(ud.GETLINE_FIG, 'SelectionType');
if (~strcmp(selectionType, 'open'))    % We don't want to add a point on the second click of a double-click
    [x,y,z] = getcurpt_eu(ud.GETLINE_AX);
    xx = get(ud.lineHandle,'XData');      yy = get(ud.lineHandle,'YData');

    % check if x,y is inside of axis
    x_lim = get(ud.GETLINE_AX,'xlim');      y_lim = get(ud.GETLINE_AX,'ylim');
    if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2))    % outside axis limits, ignore this ButtonDown
        return
    end    
    if isnan(z);    return;  end

    %sr = getappdata(ud.GETLINE_FIG,'swathRatio');
    sr = getappdata(ud.lineHandle,'swathRatio');
    rad = abs(z) * (sr/2) / 111194.9;       % metros -> degrees
    if (ud.vert_index == 1)                 % Selected first polyline vertice
        az = azimuth_geo(y,x,yy(2),xx(2));
    else
        az = azimuth_geo(yy(ud.vert_index-1),xx(ud.vert_index-1),y,x);
    end
    [lat1,lon1] = circ_geo(y,x,rad,[az-90-1 az-90+1],3);
    [lat2,lon2] = circ_geo(y,x,rad,[az+90-1 az+90+1],3);
    if (ud.inserted == 0)   % dealing with an existing vertex
        set(ud.h_circ(ud.vert_index),'XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)]);
        set(ud.h_circ(ud.vert_index),'Visible','on');
    else            % Insert new vertex
        tag = get(ud.h_circ(1),'Tag');                  % First vertice allways have the correct tag
        bar_uicm = get(ud.h_circ(1),'UIContextMenu');   % and the uicontextmenu
        ud.h_circ(length(xx)) = line('XData',[], 'YData',[],'Color',[.8 .8 .8],'LineWidth',4,...
            'Tag',tag,'UIContextMenu',bar_uicm,'UserData',ud.vert_index);    % create a new vertex
        % and actualize the h_circ userdata that contains the vertex order
        for i = ud.vert_index+1:length(xx),     set(ud.h_circ(i),'Userdata',i);     end
        if (ud.vert_index < length(xx))
            for i = length(xx):-1:ud.vert_index+1   % h_circ(vert_index+1) is wrong but it will be taken care in the next if block
                xt = get(ud.h_circ(i-1),'XData');   yt = get(ud.h_circ(i-1),'YData');
                set(ud.h_circ(i), 'XData',xt,'YData',yt);
            end
            set(ud.h_circ(ud.vert_index),'XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)]);
        else    % New vertice will be the last vertex
            set(ud.h_circ(ud.vert_index),'XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)]);
        end
        ud.inserted = 0;
    end

    % Now we have to take into account that when moving the i-th vertex, the i+1
    % vertex will have a wrong bar orientation (azimuth has changed)
    if (ud.vert_index < length(xx))
        az = azimuth_geo(y,x,yy(ud.vert_index+1),xx(ud.vert_index+1));
        [x,y,z] = getcurpt_eu(ud.GETLINE_AX,[xx(ud.vert_index+1) yy(ud.vert_index+1)]);     % get z(i+1)
        rad = abs(z) * (sr/2) / 111194.9;
        [lat1,lon1] = circ_geo(yy(ud.vert_index+1),xx(ud.vert_index+1),rad,[az-90-1 az-90+1],3);
        [lat2,lon2] = circ_geo(yy(ud.vert_index+1),xx(ud.vert_index+1),rad,[az+90-1 az+90+1],3);
        set(ud.h_circ(ud.vert_index+1),'XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)]);
    end
    if (ud.vert_index == 2)         % Also in this case i-1 (that is i=1) bar has to be adjusted
        az = azimuth_geo(yy(1),xx(1),yy(2),xx(2));
        [x,y,z] = getcurpt_eu(ud.GETLINE_AX,[xx(1) yy(1)]);     % get z(1)
        rad = abs(z) * (sr/2) / 111194.9;
        [lat1,lon1] = circ_geo(yy(1),xx(1),rad,[az-90-1 az-90+1],3);
        [lat2,lon2] = circ_geo(yy(1),xx(1),rad,[az+90-1 az+90+1],3);
        set(ud.h_circ(1),'XData', [lon1(2) lon2(2)], 'YData', [lat1(2) lat2(2)]);
    end
    set(ud.GETLINE_FIG, 'WindowButtonDownFcn', 'edit_track_mb(''FirstButtonDown'');', 'WindowButtonMotionFcn', '');        
end
set(gcbf, 'UserData', ud);

%-------------------------------------------------
% Subfunction ButtonMotion
%-------------------------------------------------
function ButtonMotion
ud = get(gcbf, 'UserData');
[newx, newy] = getcurpt_eu(ud.GETLINE_AX);
xx = get(ud.lineHandle,'XData');      yy = get(ud.lineHandle,'YData');
if (ud.vert_index == 1)                 % Selected first polyline vertice
    x = [newx xx(2:end)];   y = [newy yy(2:end)];
elseif (ud.vert_index == length(xx))    % Selected last polyline vertice
    x = [xx(1:end-1) newx];   y = [yy(1:end-1) newy];
else
    x = [xx(1:ud.vert_index-1) newx xx(ud.vert_index+1:end)];
    y = [yy(1:ud.vert_index-1) newy yy(ud.vert_index+1:end)];
end

set(ud.lineHandle,'XData',x,'YData',y);
set(gcbf, 'UserData', ud);

%---------------------------------------------------
% Subfunction CleanUp
%--------------------------------------------------
function CleanUp(xlimmode,ylimmode)
set(gca,'xlimmode',xlimmode);   set(gca,'ylimmode',ylimmode);

%---------------------------------------------------
% Subfunction getcurpt_eu
%--------------------------------------------------
function [x,y,z] = getcurpt_eu(axHandle,pt)
%GETCURPT Get current point, or use a given pt to interpolate the grid
if (nargin == 1)
    pt = get(axHandle, 'CurrentPoint');
end
x = pt(1,1);    y = pt(1,2);
if nargout == 3
    X = getappdata(gcf,'dem_x');    Y = getappdata(gcf,'dem_y');
    Z = getappdata(gcf,'dem_z');    z = abs(bi_linear(X,Y,Z,x,y));
end

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
    if (length(tag) < 7),   continue;   end
    if strcmp(tag(1:7),'swath_w')
        tmp1(i) = 0;     tmp2(i) = h_lines(i);     % get the bar (or circle) handles
    end
end
h_lines = h_lines(tmp1 ~= 0);			% get only the tracks handles
ALLbarHandles = tmp2(tmp2 ~= 0);		% get the handles to ALL bar handles (e.g. all bars of all tracks)

key=0;      ii=[];     first = 1;  button = 1;    lineHandle = [];     barHandles = [];
while key==0
	if (first),		x = pt(1,1);    y = pt(1,2);    first = 0;
	else			[x,y,button] = click_e_point(1,'crosshair');
	end
    if button~=1, lineHandle=[];    return;     end
    for i=1:length(h_lines)
        ii=[ii,i];
    end
    hh=h_lines(ii);
    if ~isempty(hh)
        if any(hh == gco), lineHandle=gco;  key=1;         end
    else
        warndlg('Sorry, there are no tracks to be edited!','Warning!')
        lineHandle=[]; return; 
    end   
end

% Now find the handles to the bar lines which correspond to the swath width
% at each vertice of the selected lineHandle
tag = get(lineHandle,'Tag');
if ~strcmp(tag(1:7),'MBtrack')      % Get out as graciously as we can
    errordlg('Error in get_trackHandle: selected line is probably not a multibeam track.','Error');
    ud = get(gcbf, 'UserData');
    set(ud.GETLINE_H1, 'UserData', 'Completed');
    set(gcbf, 'UserData', ud);
    return
end
nTrack = tag(8:end);        % get the string with the track number
barHandles = findobj(ALLbarHandles,'Tag',['swath_w' nTrack]);
barHandles = barHandles(cell2mat(get(barHandles,'Userdata')));      % they were out of order, who knows why

%---------------------------------------------------------------------------
%  key=prop_list(Type,arg1) 
%  return the handles of the specified property in the current window.
%---------------------------------------------------------------------------
function key=prop_list(Type,arg1)
if nargin==1, y=get(gcf,'Child');
else y=get(gca,'Child'); end
ii=[];
for i=1:length(y)
   c=get(y(i),'Type');
   if strcmp(c,Type), ii=[ii,i]; end
end   
key=y(ii);

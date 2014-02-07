function edit_line(varargin)
%GETLINE Select polyline with mouse.
%   EDIT_LINE(FIG) lets you select a polyline in the current axes of figure FIG
%   using the mouse.  Coordinates of the polyline are returned in X and Y.  Use normal
%   button clicks to add points to the polyline.  A shift-, right-, or double-click adds a
%   final point and ends the polyline selection. Pressing RETURN or ENTER ends the polyline
%   selection without adding a final point. Pressing BACKSPACEor DELETE removes the previously
%   selected point from the polyline.
%   EDIT_LINE(...,'closed') edits a closed polygon. However, inside 'FirstButtonDown' there
%   is a test to see if the polygon is closed. So in fact, there isn't much need of this option

%   Callback syntaxes:
%        edit_line('KeyPress')
%        edit_line('FirstButtonDown')
%        edit_line('NextButtonDown')
%        edit_line('ButtonMotion')
%
% IT IS NOT USED ANYMORE. IT WAS SURPASSED BY UI_EDIT_POLYGON

%xlimorigmode = xlim('mode');    ylimorigmode = ylim('mode');
xlimorigmode = get(gca,'xlimmode');    ylimorigmode = get(gca,'ylimmode');
%xlim('manual');                 ylim('manual');
set(gca,'xlimmode','manual');   set(gca,'ylimmode','manual');

if ((nargin >= 1) && (ischar(varargin{end})))
    str = varargin{end};
    if (str(1) == 'c')        % getline(..., 'closed')
        ud.GETLINE_ISCLOSED = 1;        varargin = varargin(1:end-1);
    end
else
    ud.GETLINE_ISCLOSED = 0;
end

if ((length(varargin) >= 1) && ischar(varargin{1}))
    % Callback invocation: 'KeyPress', 'FirstButtonDown', 'NextButtonDown', or 'ButtonMotion'.
    feval(varargin{:});    return;
end

ud.GETLINE_AX = gca;    ud.GETLINE_FIG = get(ud.GETLINE_AX, 'Parent');
   
% Remember initial figure state
old_db = get(ud.GETLINE_FIG, 'DoubleBuffer');
state= uisuspend_j(ud.GETLINE_FIG);

% Set up initial callbacks for initial stage
set(ud.GETLINE_FIG, 'Pointer', 'crosshair', ...
    'WindowButtonDownFcn', 'edit_line(''FirstButtonDown'');',...
    'KeyPressFcn', 'edit_line(''KeyPress'');', 'DoubleBuffer', 'on');

% Bring target figure forward
figure(ud.GETLINE_FIG);

% This had a totaly different use in the original getline, but had reveled very usefull
% here for it's handle is used to store information that if stored in ud.lineHandle gave
% a very strange bug on a second call of edit_line (e.g. after clicking the right or
% midle button, it didn't work anymore because ud structure was unknown (???))
ud.GETLINE_H1 = line('Parent', ud.GETLINE_AX, 'XData', [], 'YData', [], ...
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
    x = [];     y = [];
case 'unknown'  % User did something to cause the polyline selection to terminate abnormally.
    x = [];     y = [];
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
        x = get(ud.lineHandle,'XData');     y = get(ud.lineHandle,'YData');
        x(ud.vert_index) = [];              y(ud.vert_index) = [];
        set(ud.lineHandle,'XData',x,'YData',y);
        set(ud.GETLINE_FIG, 'WindowButtonDownFcn', 'edit_line(''FirstButtonDown'');', 'WindowButtonMotionFcn', '');
    case char(105)              % "i" char key (insert vertice)
        [newx, newy] = getcurpt(ud.GETLINE_AX);
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
        ud.vert_index = ud.vert_index + 1;        
    case char(98)              % "b" char key (break line)
        xx = get(ud.lineHandle,'XData');    yy = get(ud.lineHandle,'YData');
        xx(ud.vert_index) = ud.save_x;      yy(ud.vert_index) = ud.save_y;
        if (ud.vert_index > 1 && ud.vert_index < length(xx))
            x = xx(1:ud.vert_index);            y = yy(1:ud.vert_index);
            set(ud.lineHandle,'XData',x,'YData',y);
            % Now make the a new segment from rest of the original one
            uictxm = get(ud.lineHandle,'uicontextmenu');
            tmp = line('XData',[], 'YData',[],'LineWidth',1);  % create a new line handle
            set(tmp, 'XData',xx(ud.vert_index+1:end), 'YData',yy(ud.vert_index+1:end),'uicontextmenu',uictxm)
            set(ud.GETLINE_H1, 'UserData', 'Completed');
        end
    case char(99)              % "c" char key (close line)
        xx = get(ud.lineHandle,'XData');    yy = get(ud.lineHandle,'YData');
        if length(xx) > 2    % don't close a line with less than 2 vertex
            x = [xx xx(1)];                     y = [yy yy(1)];
            set(ud.lineHandle,'XData',x,'YData',y);
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
% I have to do this test here because there is another mouse click inside get_lineHandle
if (~strcmp(get(ud.GETLINE_FIG, 'SelectionType'), 'normal'))    % We're done!
    set(ud.GETLINE_H1, 'UserData', 'Completed');
else
    current_pt = get(ud.GETLINE_AX, 'CurrentPoint');
    ud.lineHandle = get_lineHandle(current_pt);
    if isempty(ud.lineHandle)       % User gave up (right or midle click) inside get_lineHandle
        set(ud.GETLINE_H1, 'UserData', 'Completed');
    else
        x = get(ud.lineHandle,'XData');   y = get(ud.lineHandle,'YData');
        if length(ud.lineHandle) > 0
            % Check to see if we are dealing with a closed polyline
            if ( (x(1) == x(end)) && (y(1) == y(end)) ),  ud.GETLINE_ISCLOSED = 1;    end
            % Find out which polyline vertice is closest to "current_pt".
            dif_x = x - current_pt(1,1);    dif_y = y - current_pt(1,2);
            dist = sqrt(dif_x.^2 + dif_y.^2);
            [B,IX] = sort(dist);    ud.vert_index = IX(1);      % tenho que ser mais esperto (potenc. lento)
            ud.save_x = x(ud.vert_index);   ud.save_y = y(ud.vert_index);   % needed in the "i"nsert option
        end
        set(ud.GETLINE_FIG, 'WindowButtonMotionFcn', 'edit_line(''ButtonMotion'');', ...
            'WindowButtonDownFcn', 'edit_line(''NextButtonDown'');');
        ud.first_bm = 1;
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
    [x,y] = getcurpt(ud.GETLINE_AX);

    % check if GETLINE_X,GETLINE_Y is inside of axis
    x_lim = get(ud.GETLINE_AX,'xlim');      y_lim = get(ud.GETLINE_AX,'ylim');
    if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2))    % outside axis limits, ignore this ButtonDown
        return
    end    
    
    set(ud.GETLINE_FIG, 'WindowButtonDownFcn', 'edit_line(''FirstButtonDown'');', 'WindowButtonMotionFcn', '');        
end
set(gcbf, 'UserData', ud);

%-------------------------------------------------
% Subfunction ButtonMotion
%-------------------------------------------------
function ButtonMotion
ud = get(gcbf, 'UserData');
[newx, newy] = getcurpt(ud.GETLINE_AX);
xx = get(ud.lineHandle,'XData');      yy = get(ud.lineHandle,'YData');
if (ud.first_bm),    ud.first_bm = 0;   end
if (ud.vert_index == 1)                 % Selected first polyline vertice
    if (ud.GETLINE_ISCLOSED)
        x = [newx xx(2:end-1) newx];    y = [newy yy(2:end-1) newy];
    else
        x = [newx xx(2:end)];           y = [newy yy(2:end)];
    end
elseif (ud.vert_index == length(xx))    % Selected last polyline vertice
    if (ud.GETLINE_ISCLOSED)
        x = [newx xx(2:end-1) newx];    y = [newy yy(2:end-1) newy];
    else
        x = [xx(1:end-1) newx];         y = [yy(1:end-1) newy];
    end
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
%xlim(xlimmode);     ylim(ylimmode);
set(gca,'xlimmode',xlimmode);   set(gca,'ylimmode',ylimmode);

%---------------------------------------------------
% Subfunction getcurpt
%--------------------------------------------------
function [x,y] = getcurpt(axHandle,pt)
%GETCURPT Get current point
if (nargin == 1)
    pt = get(axHandle, 'CurrentPoint');
end
x = pt(1,1);    y = pt(1,2);

%---------------------------------------------------------------------
function [lineHandle,button]=get_lineHandle(pt)
%get_lineHandle function is used to get the handle of a line object.
%[lineHandle,button]=get_lineHandle(pt)
% This subfunction and the next are adapted from GraphTools
%---------------------------------------------------------------------
ii = prop_list('axes');     h_lines=[];
for i=1:length(ii)
   axes(ii(i))
   h_lines=[h_lines; prop_list('line',1)]; % find the lines handles
end

key=0;      ii=[];     first = 1;  button = 1;    lineHandle = [];
while key==0
    if (first),     x = pt(1,1);    y = pt(1,2);    first = 0;
    else           [x,y,button] = click_e_point(1,'crosshair');     end
    if button~=1, lineHandle=[];    return;     end
    for i=1:length(h_lines)
        ii=[ii,i];
    end
    hh=h_lines(ii);
    if length(hh) > 0
        if any(hh == gco), lineHandle=gco;  key=1;         end
    else
        warndlg('Sorry, there are no lines to be edited!','Warning!')
        lineHandle=[]; return; 
    end   
end

%---------------------------------------------------------------------------
%  key=prop_list(Type,arg1) 
%  return the handles of the specified property in the current window.
%---------------------------------------------------------------------------
function key = prop_list(Type,arg1)
if (nargin == 1), y = get(gcf,'Child');
else, y = get(gca,'Child'); end
ii = [];
for i=1:length(y)
   c = get(y(i),'Type');
   if strcmp(c,Type), ii = [ii,i];  end
end   
key = y(ii);

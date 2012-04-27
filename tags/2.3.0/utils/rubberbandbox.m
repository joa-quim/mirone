function [p1,p2,lh1] = rubberbandbox(varargin)
% Function to draw a rubberband box and return the start and end points
%
% Usage: [p1,p2]=rubberbox;     uses current axes
%        [p1,p2]=rubberbox(h);  uses axes refered to by handle, h
%        [p1,p2,lh]=rubberbox(...);  returns handle of the rubberbox and dont't erase it
% Based on an idea of Sandra Martinka's Rubberline
% Edited by Bob Hamans (B.C.Hamans@student.tue.nl) 02-04-2003
% Re-written for Mirone use (J. Luis).
% E.G rectangles have origin at LowerLeft corner and handles to rect lines are cw oriented.

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

%Check for optional argument
switch nargin
    case 0
        h = gca;
    case 1
        h = varargin{1};  %axes(h);
    otherwise
        errordlg('Too many input arguments.','Error');
end

hFig = get(h,'Parent');
set(hFig, 'CurrentAxes', h)			% Was not necessarely the case

% Get current user data
state = uisuspend_j(hFig);			% Remember initial figure state
cudata = get(hFig,'UserData'); 
set(hFig,'Pointer','Crosshair')
% Wait for left mouse button to be pressed
waitforbuttonpress;
h = get(hFig,'CurrentAxes');		% O multi-axes previous current axes may not be the same as current one (eg gmtedit)

p1 = get(h,'CurrentPoint');			% get starting point
p1 = p1(1,1:2);						% extract x and y

lh1 = line('XData',p1(1),'YData',p1(2),'Parent',h,'Color', 'k', 'LineStyle', '-');		%plot starting point
lh2 = line('XData',p1(1),'YData',p1(2),'Parent',h,'Color', 'w', 'LineStyle', ':');

% Save current point (p1) data in a structure to be passed
udata.p1 = p1;    udata.h = h;  udata.lh = [lh1 lh2];
% Set gcf object properties 'UserData' and call function 'wbmf' on mouse motion.
db = get(hFig,'DoubleBuffer');       % Get current DoubleBuffer state
set(hFig,'WindowButtonMotionFcn',{@wbmf,udata},'WindowButtonDownFcn',{@wbd,udata},'DoubleBuffer','on');

waitfor(lh1, 'UserData', 'Completed');

% Get data for the end point
p2 = get(h,'Currentpoint');       %get end point
p2 = p2(1,1:2);                   %extract x and y
set(hFig,'UserData',cudata,'WindowButtonMotionFcn','','DoubleBuffer',db); %reset UserData, etc..
if (nargout == 2)   % Remove the rectangle and don't return it's handle (otherwise do both)
    delete(lh1);    delete(lh2);
    % Make shure that p1 is LowLeft and p2 UpperRight
    x_min = min([p1(1) p2(1)]);     x_max = max([p1(1) p2(1)]);
    y_min = min([p1(2) p2(2)]);     y_max = max([p1(2) p2(2)]);
    p1(1) = x_min;  p1(2) = y_min;  p2(1) = x_max;  p2(2) = y_max;
elseif (nargout == 3)
    delete(lh2);           % Delete this extra handle (only lh1 is returned)
    % Make shure that the rectangle is clockwise oriented with origin at the lower left
    x_min = min([p1(1) p2(1)]);    x_max = max([p1(1) p2(1)]);
    y_min = min([p1(2) p2(2)]);    y_max = max([p1(2) p2(2)]);
    p1(1) = x_min;  p1(2) = y_min;  p2(1) = x_max;  p2(2) = y_max;
    set(lh1,'XData',[x_min,x_min,x_max,x_max,x_min],'YData',[y_min,y_max,y_max,y_min,y_min])
end
uirestore_j(state, 'nochildren');         % Restore the figure's initial state

% -------------------------------------------------------------------------
function wbmf(obj,eventdata,utemp)   % window motion callback function
	ptemp = get(utemp.h,'CurrentPoint');
	ptemp = ptemp(1,1:2);
	% Use 5 point to draw a rectangular rubberband box
	set(utemp.lh(1),'XData',[ptemp(1),ptemp(1),utemp.p1(1),utemp.p1(1),ptemp(1)], ...
        'YData',[ptemp(2),utemp.p1(2),utemp.p1(2),ptemp(2),ptemp(2)]);
	set(utemp.lh(2),'XData',[ptemp(1),ptemp(1),utemp.p1(1),utemp.p1(1),ptemp(1)], ...
        'YData',[ptemp(2),utemp.p1(2),utemp.p1(2),ptemp(2),ptemp(2)]);

function wbd(obj,eventdata,udata)
    set(udata.lh(1),'UserData','Completed')

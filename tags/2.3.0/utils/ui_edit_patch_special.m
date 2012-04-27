function ui_edit_patch_special(varargin)
% ui_edit_patch_special: Interactive move some type of patches.
%                        Currently only used to move focal mechanims
%  ui_edit_patch_special. Double clicking on the patch displays its vertices limits.
%  To move, click and drag on the patch. To stop "moving mode" double click again.

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

if (nargin == 1 && ischar(varargin{1}))         % ... callbacks
    patchui(varargin{:});
    return
elseif (nargin == 1 && ishandle(varargin{1}))  % Handle in input
    if (~strcmp(get(varargin{1},'Type'),'patch'))
        errordlg('Error in ui_edit_patch_special: Invalid handle.','Error')
        return
    end
    h_patch = varargin{1};
else
    errordlg('Error in ui_edit_patch_special: Invalid number of inputs.','Error')
    return
end

%  Set the correct structures
s.controls = 'off';
s.h_axes = get(h_patch,'parent');
s.h_fig = get(s.h_axes,'parent');
s.h_patch = h_patch;
s.markers = [];
s.other_h = getappdata(h_patch,'other_hand');
if (length(s.other_h) == 3)
    s.have_title = 1;
else
    s.have_title = 0;
end
s.limits = getappdata(h_patch,'Limits');
s.parent = get(s.h_patch,'parent');	
s.uicontext = get(h_patch,'uicontextmenu');   % Save the eventual uicontextmenu
s.old_ud = get(h_patch,'UserData');           % Save the eventual userdata
set(h_patch,'userdata',s,'buttondownfcn','ui_edit_patch_special(''poligmousedown'')');

%---------------------------------------------------------------------------------------------------
function patchui(action)

switch action
case 'poligmousedown'              % mouse down on line/patch
	stype = get(gcf,'selectiontype');
    curr_obj = gco;
	s = get(curr_obj,'userdata');
	s.parent = get(curr_obj,'parent');	
	switch stype % shift-click or double click to toggle control points
		case 'open'
				switch s.controls
					case 'on'
                        delete(s.markers);        s.markers = [];       % delete vrtice markers
						s.controls = 'off';
                        set(s.h_patch,'uicontextmenu',s.uicontext)
						set(curr_obj,'userdata',s)
					case 'off'
						s.controls = 'on';
                        hold on;    
						s.markers = plot(s.limits(:,1),s.limits(:,2),'sk','MarkerFaceColor','w', ...
                            'MarkerSize',5,'userdata',s);
                        hold off;
                        set(s.h_patch,'buttondownfcn',{@move_meca,curr_obj},'uicontextmenu','');
						set(curr_obj,'userdata',s)
				end
	end
end

% -----------------------------------------------------------------------------------------
function move_meca(obj,eventdata,h)
% Move the patch whose handle is h
s = get(h,'userdata');
stype = get(s.h_fig,'selectiontype');

if (strcmp(stype,'open'))
    delete(s.markers);        s.markers = [];       % delete vrtice markers
	s.controls = 'off';
    set(s.h_patch,'userdata',s,'buttondownfcn','ui_edit_patch_special(''poligmousedown'')','uicontextmenu',s.uicontext);
	set(gco,'userdata',s)
    return
end

state = uisuspend_j(s.h_fig);					% Remember initial figure state
x_lim = get(s.h_axes,'xlim');        y_lim = get(s.h_axes,'ylim');
current_pt = get(s.h_axes, 'CurrentPoint');
setappdata(s.h_patch,'old_pt',[current_pt(1,1) current_pt(1,2)])

set(s.h_fig,'WindowButtonMotionFcn',{@wbm_MoveMeca,h,[x_lim y_lim]},...
    'WindowButtonUpFcn',{@wbu_MoveMeca,h,state}, 'Pointer','fleur');

% ---------
function wbm_MoveMeca(obj,eventdata,h,lim)
s = get(h,'userdata');
pt = get(s.h_axes, 'CurrentPoint');
if (pt(1,1)<lim(1)) || (pt(1,1)>lim(2)) || (pt(1,2)<lim(3)) || (pt(1,2)>lim(4));   return; end
old_pt = getappdata(s.h_patch,'old_pt');
dx = pt(1,1) - old_pt(1);       dy = pt(1,2) - old_pt(2);

x1 = get(s.h_patch,'XData');    y1 = get(s.h_patch,'YData');        % Get the selected "MecaPart"
x1 = x1 + dx;                   y1 = y1 + dy;

x2 = get(s.other_h(1),'XData'); y2 = get(s.other_h(1),'YData');     % Get the other "MecaPart"
x2 = x2 + dx;                   y2 = y2 + dy;

x3 = get(s.other_h(2),'XData'); y3 = get(s.other_h(2),'YData');     % Get the anchor points
x3(2) = x3(2) + dx;             y3(2) = y3(2) + dy;

xv = get(s.markers,'XData');    yv = get(s.markers,'YData');        % Get markers
xv = xv + dx;                   yv = yv + dy;

set(s.markers, 'XData',xv, 'YData',yv);
set(s.h_patch, 'XData',x1, 'YData',y1);
set(s.other_h(1), 'XData',x2, 'YData',y2);
set(s.other_h(2), 'XData',x3, 'YData',y3);
if (s.have_title)
	pos = get(s.other_h(3),'Position');                     % Get title position
	pos(1) = pos(1) + dx;       pos(2) = pos(2) + dy;
    set(s.other_h(3),'Position',pos);
end
setappdata(s.h_patch,'old_pt',[pt(1,1) pt(1,2)])

% ---------
function wbu_MoveMeca(obj,eventdata,h,state)
s = get(h,'userdata');
x = get(s.markers,'XData');     y = get(s.markers,'YData');
s.limits = [x(:) y(:)];
set(h,'userdata',s)
uirestore_j(state, 'nochildren');         % Restore the figure's initial state

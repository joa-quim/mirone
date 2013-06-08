function varargout = arrow_shape(varargin)
% Helper window to set up arrow shape parameters

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

% $Id$

	if (isempty(varargin))
		errordlg('ARROW_SHAPE: wrong number of input arguments.','Error'),	return
	end

	hObject = figure('Tag','figure1','Visible','off');
	arrow_shape_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center');

	handles.hArrow = varargin{3};
	ud = get(varargin{3}, 'UserData');

	handles.headLength = ud.headLength;
	set(handles.edit_len, 'Str', sprintf('%d',handles.headLength))
	edit_len_CB(handles.edit_len, handles, '')		% It also takes care of eventual SliderStep reseting

	handles.aspectRatio = ud.aspectRatio;
	set(handles.edit_ratio, 'Str', sprintf('%.2f',handles.aspectRatio))
	edit_ratio_CB(handles.edit_ratio, handles, '')

	handles.vFac = ud.vFac;
	set(handles.edit_shape, 'Str', sprintf('%.2f',handles.vFac))
	edit_shape_CB(handles.edit_shape, handles, '')

	handles.lineThickness = get(varargin{3},'LineWidth');
	set(handles.edit_lineThickness, 'Str', sprintf('%.1f',handles.lineThickness))

	% Scale factor for make_arrow
	axLims = [get(handles.axes1,'XLim') get(handles.axes1,'YLim')];
	oldUnits = get(handles.axes1,'Units');			set(handles.axes1,'Units','points');
	Pos = get(handles.axes1,'Position');			set(handles.axes1,'Units',oldUnits);
	vscale = 1/Pos(4) * diff(axLims(3:4));	hscale = 1/Pos(3) * diff(axLims(1:2));
	handles.hscale = (vscale + hscale) / 2;

	% Plot the demo arrow
	[xt, yt] = make_arrow([-1 1; 0 0], handles.hscale, handles.hscale, ud.headLength, ud.vFac, ud.aspectRatio);
	handles.hDemoArrow = patch('XData',xt, 'YData', yt, 'LineWidth',handles.lineThickness);

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% -----------------------------------------------------------------------------
function slider_len_CB(hObject, handles)
	val = round(get(hObject,'Value'));
	set(handles.edit_len,'String',sprintf('%d',val))
	update_arrow(handles)

% -----------------------------------------------------------------------------
function slider_ratio_CB(hObject, handles)
	val = get(hObject,'Value');
	set(handles.edit_ratio,'String',sprintf('%.2f',val))
	update_arrow(handles)

% -----------------------------------------------------------------------------
function slider_shape_CB(hObject, handles)
	val = get(hObject,'Value');
	set(handles.edit_shape,'String',sprintf('%.2f',val))
	update_arrow(handles)

% -----------------------------------------------------------------------------
function edit_len_CB(hObject, handles, opt)
% The OPT (whatever) is just used to NOT call update_arrow()
	x = round(str2double(get(hObject,'String')));
	if (isnan(x))
		set(hObject,'String',handles.headLength),	return
	end
	if (x > get(handles.slider_len,'Max'))		% If value greater than max slide val, augment this
		set(handles.slider_len,'Max',x)
		set(handles.slider_shape,'SliderStep',[1/x 5/x])
	end
	set(handles.slider_len,'Value',x)
	if (nargin == 2),	update_arrow(handles),	end

% -----------------------------------------------------------------------------
function edit_ratio_CB(hObject, handles, opt)
% The OPT (whatever) is just used to NOT call update_arrow()
	x = str2double(get(hObject,'String'));
	if (isnan(x))
		set(hObject,'String',handles.aspectRatio),	return
	end
	if (x > get(handles.slider_ratio,'Max'))		% If value greater than max slide val, augment this
		set(handles.slider_ratio,'Max',x)
		set(handles.slider_ratio,'SliderStep',[1/x 10/x]*1e-1)
	end
	set(handles.slider_ratio,'Value',x)
	if (nargin == 2),	update_arrow(handles),	end

% -----------------------------------------------------------------------------
function edit_shape_CB(hObject, handles, opt)
% The OPT (whatever) is just used to NOT call update_arrow()
	x = str2double(get(hObject,'String'));
	if (isnan(x))
		set(hObject,'String',handles.vFac),	return
	end
	if (x > get(handles.slider_shape,'Max'))		% If value greater than max slide val, augment this
		set(handles.slider_shape,'Max',x)
		set(handles.slider_shape,'SliderStep',[1/x 4/x]*1e-1)
	end
	set(handles.slider_shape,'Value',x)
	if (nargin == 2),	update_arrow(handles),	end

% -----------------------------------------------------------------------------
function edit_lineThickness_CB(hObject, handles)
	x = str2double(get(hObject,'String'));
	if (isnan(x))
		set(hObject,'String',handles.lineThickness)
	end
	update_arrow(handles)

% -----------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% Update the original arrow that was used to call this tool
	len = round(str2double(get(handles.edit_len,'Str')));
	vFac = str2double(get(handles.edit_shape,'Str'));
	rat = str2double(get(handles.edit_ratio,'Str'));
	lth = str2double(get(handles.edit_lineThickness,'Str'));
	ud = get(handles.hArrow, 'UserData');
	[xt, yt] = make_arrow(ud.anchors', ud.hscale, ud.hscale, len, vFac, rat);
	ud.arrow_xy = [xt(:) yt(:)];
	ud.vFac = vFac;			ud.headLength = len;		ud.aspectRatio = rat;
	set(handles.hArrow, 'XData', xt, 'YData', yt, 'LineWidth', lth, 'UserData',ud)

% -----------------------------------------------------------------------------
function update_arrow(handles)
% Update the demo arrow
	len = round(str2double(get(handles.edit_len,'Str')));
	vFac = str2double(get(handles.edit_shape,'Str'));
	rat = str2double(get(handles.edit_ratio,'Str'));
	lth = str2double(get(handles.edit_lineThickness,'Str'));
	[xt, yt] = make_arrow([-1 1; 0 0], handles.hscale, handles.hscale, len, vFac, rat);
	set(handles.hDemoArrow, 'XData', xt, 'YData', yt, 'LineWidth', lth)

% -----------------------------------------------------------------------------
function arrow_shape_LayoutFcn(h1)

set(h1,...
'Position',[520 630 325 165],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Arrow shape',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[136 143 151 15],...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@arrow_shape_uiCB,...
'Style','slider',...
'Max',40,...
'Min',0,...
'SliderStep',[1/40 5/40],...
'TooltipString',sprintf('Arrow head lenght in points\nYou may increase the maximum manually'),...
'Tag','slider_len');

uicontrol('Parent',h1, 'Position',[136 115 151 15],...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@arrow_shape_uiCB,...
'Style','slider',...
'Max',5,...
'Min',0,...
'SliderStep',[0.02 0.2],...
'TooltipString',sprintf('Arrow head length/width aspect ratio\nYou may increase the maximum manually'),...
'Tag','slider_ratio');

uicontrol('Parent',h1, 'Position',[136 87 151 15],...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@arrow_shape_uiCB,...
'Style','slider',...
'Max',2,...
'Min',0,...
'SliderStep',[0.05 0.2],...
'TooltipString',sprintf('Arrow head shape factor\nYou may increase the maximum manually'),...
'Tag','slider_shape');

uicontrol('Parent',h1, 'Position',[286 141 35 20],...
'BackgroundColor',[1 1 1],...
'Callback',@arrow_shape_uiCB,...
'Style','edit',...
'TooltipString',sprintf('Arrow head lenght in points\nYou may increase the maximum manually'),...
'Tag','edit_len');

uicontrol('Parent',h1, 'Position',[286 113 35 20],...
'BackgroundColor',[1 1 1],...
'Callback',@arrow_shape_uiCB,...
'Style','edit',...
'TooltipString',sprintf('Arrow head length/width aspect ratio\nYou may increase the maximum manually'),...
'Tag','edit_ratio');

uicontrol('Parent',h1, 'Position',[286 85 35 20],...
'BackgroundColor',[1 1 1],...
'Callback',@arrow_shape_uiCB,...
'Style','edit',...
'TooltipString',sprintf('Arrow head shape factor\nYou may increase the maximum manually'),...
'Tag','edit_shape');

uicontrol('Parent',h1,'Position',[286 50 35 21],...
'BackgroundColor',[1 1 1],...
'Callback',@arrow_shape_uiCB,...
'Style','edit',...
'TooltipString','Outer line thickness in points',...
'Tag','edit_lineThickness');

uicontrol('Parent',h1, 'Position',[33 143 100 14],...
'HorizontalAlignment','right',...
'String','Arrow head length',...
'Style','text');

uicontrol('Parent',h1, 'Position',[2 116 130 14],...
'HorizontalAlignment','right',...
'String','Arrow head aspect ratio',...
'Style','text');

uicontrol('Parent',h1, 'Position',[4 88 130 14],...
'HorizontalAlignment','right',...
'String','Arrow head shape factor',...
'Style','text');

uicontrol('Parent',h1, 'Position',[223 54 60 14],...
'HorizontalAlignment','right',...
'String','Line Width',...
'Style','text');

axes('Parent',h1,...
'Units','pixels',...
'Position',[5 8 216 70],...
'CameraPosition',[0 0 9.16025403784439],...
'Visible','off',...
'XLim',[-1.01 1.01],...
'XLimMode','manual',...
'YLim',[-0.33333 0.33333],...
'YLimMode','manual',...
'Tag','axes1');

uicontrol('Parent',h1, 'Position',[246 4 76 21],...
'Callback',@arrow_shape_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','Apply',...
'TooltipString','Apply changes to edited object',...
'Tag','push_OK');

function arrow_shape_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

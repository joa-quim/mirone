function varargout = world_is_not_round_enough(varargin)
% Helper tool to convert between-to [-180 180] <-> [0 360] longitude ranges

%	Copyright (c) 2004-2016 by J. Luis
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

% $Id: world_is_not_round_enough.m 7813 2016-02-26 21:48:41Z j $

	if (numel(varargin) == 0),		return,		end
 
	hObject = figure('Tag','figure1','Visible','off');
	world_is_not_round_enough_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	handMir = varargin{1};

	if (~handMir.geog)
		warndlg('This tool makes sense only when using geographical coordinates. Quiting','WarnError')
		delete(hObject),	return
	elseif (handMir.no_file)
		delete(hObject),	return
	end

	handles.hMirFig = handMir.figure1;
	handles.hMirAxes = handMir.axes1;
	handles.hMirImg = handMir.hImg;
	handles.head = handMir.head;
	handles.validGrid = handMir.validGrid;
	handles.isGblobalPixReg = false;
	d_path = handMir.path_data;
	handles.dx = diff(handles.head(1:2));
	handles.square = [handles.head(1) handles.head(3);	handles.head(1) handles.head(4); ...
					  handles.head(2) handles.head(4);	handles.head(2) handles.head(3); ...
					  handles.head(1) handles.head(3)];

	% First check if it's a global pix-reg image
	if ( abs(handles.dx - (360 - handles.head(8))) < 1e-5 )		% Yes, it is.
		handles.isGblobalPixReg = true;
		handles.square(1,1) = handles.square(1,1) - handles.head(8)/2;
		handles.square(2,1) = handles.square(1,1);
		handles.square(3,1) = handles.square(3,1) + handles.head(8)/2;
		handles.square(4,1) = handles.square(3,1);
		handles.square(5,1) = handles.square(1,1);
	end

	if ( (handles.dx < (360 - handles.head(8)/2)) && ~handles.isGblobalPixReg )
		%set([handles.push_to360 handles.push_to180],'Enable','inactive')
		str = sprintf(['This button is blocked because the map is not global\n' ...
						'That is, it does not cover 360 degrees of longitude\n' ...
						'You can still click and drag the rectangle and try luck.']);
		set([handles.push_to360 handles.push_to180],'Enable','off','Tooltip', str)
	else
		set(handles.push_to360,'Tooltip', 'Convert map longitudes from [-180 180] to [0 360].')
		set(handles.push_to180,'Tooltip', 'Convert map longitudes from [0 360] to [-180 180].')
	end
	set(handles.edit_xMin,'String', sprintf('%.10g',handles.square(1,1)) )
	set(handles.edit_xMax,'String', sprintf('%.10g',handles.square(3,1)) )

	% Import the world map
	handles.w_map = flipdim(imread([d_path 'logo_etopo2.jpg']),1);
	siza = size(handles.w_map);
	handles.w_map = [handles.w_map handles.w_map(:,1:fix(siza(2)/2)+1,:) ];
	
	image([-180 360],[-90 90],handles.w_map,'Parent',handles.axes1);
	set(handles.axes1,'YDir','normal');
	handles.hRectLims = patch('Parent',handles.axes1, 'xdata', handles.square(:,1), ...
						'ydata', handles.square(:,2), 'LineWidth',3, 'FaceColor','none');
	set(handles.hRectLims,'buttondownfcn',{@moveRect, handles, handles.hRectLims});

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hMirFig,'dependentFigs',plugedWin);

	% Choose default command line output for world_is_not_round_enough
	set(hObject,'Visible','on')
	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end

% -------------------------------------------------------------------------------------
function edit_xMin_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < -180 || xx >= 360)
		set(hObject,'String', sprintf('%.6f',handles.square(1,1)) )
		set(handles.edit_xMax,'String', sprintf('%.6f',handles.square(3,1)) )
		set(handles.hRectLims, 'xdata', handles.square(:,1), 'ydata', handles.square(:,2))
		x(3) = handles.square(3,1);
	elseif (xx + handles.dx > 360 )	% Do not let it get out
		x = handles.square(:,1);
		x(1) = 360 - handles.dx;
		x(2) = x(1);	x(5) = x(1);
		x(3) = 360;		x(4) = 360;
		set(hObject,'String', sprintf('%.6f',x(1)) )
		set(handles.hRectLims, 'xdata', x, 'ydata', handles.square(:,2))
	else
		x = handles.square(:,1);
		x(1) = xx;		x(2) = xx;		x(5) = xx;
		x(3) = xx + handles.dx;			x(4) = x(3);
		set(handles.hRectLims, 'xdata', x, 'ydata', handles.square(:,2))
	end
	set(handles.edit_xMax,'String', sprintf('%.6f',x(3)) )

% -------------------------------------------------------------------------------------
function push_to360_CB(hObject, handles)
% [-180 180] -> [0 360]
	rect_x = get(handles.hRectLims, 'xdata');
	rect_x = rect_x + 180;
	if (rect_x(4) > 360),	return,		end
	set(handles.hRectLims, 'XData',rect_x);
	push_apply_CB(handles.push_apply, handles)

% -------------------------------------------------------------------------------------
function push_to180_CB(hObject, handles)
% [0 360] -> [-180 180]
	rect_x = get(handles.hRectLims, 'xdata');
	rect_x = rect_x - 180;
	if (rect_x(1) < -180),	return,		end
	set(handles.hRectLims, 'XData',rect_x);
	push_apply_CB(handles.push_apply, handles)

% -------------------------------------------------------------------------------------
function push_apply_CB(hObject, handles)

	rect_x = get(handles.hRectLims, 'xdata');
	eps_x = handles.head(8);
	XYlim = getappdata(handles.hMirAxes,'ThisImageLims');
	if ((abs(rect_x(1) - XYlim(1)) < eps_x) && (abs(rect_x(3) - XYlim(2)) < eps_x))		% Nothing changed, so by by
		return
	end

	dy = diff(handles.square(1:2,2));
	x_min = rect_x(1);		x_max = rect_x(3);
	handMir = guidata(handles.hMirFig);			% Get the Mirone handles
	if (handles.validGrid)
		[X,Y,Z] = load_grd(handMir);
		if (isempty(Z)),	return,		end
	else
		Z = get(handles.hMirImg,'CData');
	end

	[Z, handMir] = to_from_180(handles, handMir, Z, x_min, x_max, dy, eps_x);

	% Update the Mirone window
	if (handles.validGrid)
		X = linspace(x_min, x_max, size(Z,2));
		setappdata(handMir.figure1,'dem_x',X);
		setappdata(handMir.figure1,'dem_z',Z);
		Z = get(handles.hMirImg,'CData');
		Z = to_from_180(handles, handMir, Z, x_min, x_max, dy, eps_x);
		set(handles.hMirAxes,'XLim',[X(1) X(end)])
	else
		X(1) = x_min;		X(2) = x_max;
		set(handles.hMirAxes,'XLim',X)
	end
	handles.head(1:2) = [X(1) X(end)];		guidata(handles.figure1, handles)
	handMir.head(1:2) = [X(1) X(end)];
	if (handles.isGblobalPixReg)			% Need to reset the grid-reg coords (they were extended)
		handMir.head(1:2) = handMir.head(1:2) + [handMir.head(8) -handMir.head(8)]/2;
	end
	guidata(handles.hMirFig, handMir)

	set(handles.hMirImg,'CData', Z, 'XData', X);
	set(handles.edit_xMin,'Str', X(1))
	set(handles.edit_xMax,'Str', X(end))
	setappdata(handles.hMirAxes,'ThisImageLims',[get(handles.hMirAxes,'XLim') get(handles.hMirAxes,'YLim')])

% -----------------------------------------------------------------------------------------
function [Z, handMir] = to_from_180(handles, handMir, Z, x_min, x_max, dy, eps_x)
% Note that the x_min,x_max here it's already the destination because they were already updated
	n_col = size(Z,2);		ncM = fix(n_col / 2);
	if (x_min >= 0 && (abs(x_max - 360) < eps_x))			% [-180 180] -> [0 360]
		Z = [Z(:,ncM+1:end,:) Z(:,1:ncM,:)];
		handMir.geog = 2;
		wrap_lines(handMir, handMir.geog)
	elseif (x_min < 180 && (abs(x_max - 180) < eps_x))		% [0 360] -> [-180 180]
		Z = [Z(:,ncM+1:end,:) Z(:,1:ncM,:)];
		handMir.geog = 1;
		wrap_lines(handMir, handMir.geog)
	end

% -----------------------------------------------------------------------------------------
function wrap_lines(handMir, geog)
% When wrapping around from [-180 180] -> [0 360] or inverse line vertices may be lost. Restore them.
	hLines = findobj(handMir.axes1, 'Type', 'line');
	for (k = 1: numel(hLines))
		x = get(hLines(k), 'XData');
		if (geog == 2)		% Means, it used to be 1 ([-180 180])
			ind = (x < 0);
			x(ind) = x(ind) + 360;
			set(hLines(k), 'XData', x);
		else				% Was [0 360], will be [-180 180]
			ind = (x > 180);
			x(ind) = x(ind) - 360;
			set(hLines(k), 'XData', x);
		end
	end

% -----------------------------------------------------------------------------------------
function moveRect(obj, evt, handles, hRect)
	hFig = handles.figure1;			hAxes = handles.axes1;
	hXmin = handles.edit_xMin;		hXmax = handles.edit_xMax;
	state = uisuspend_j(hFig);            % Remember initial figure state
	x_lim = get(hAxes,'xlim');		y_lim = get(hAxes,'ylim');
	current_pt = get(hAxes, 'CurrentPoint');
	set(hRect,'UserData',[current_pt(1,1) current_pt(1,2)])
	
	set(hFig,'WindowButtonMotionFcn',{@wbm_MoveRect,hRect,[x_lim y_lim], handles.dx, hAxes, hXmin, hXmax},...
        'WindowButtonUpFcn',{@wbu_MoveRect,hRect,state}, 'Pointer','fleur');

% ---------
function wbm_MoveRect(obj,eventdata, h, lim, width, hAxes, hXmin, hXmax)
	pt = get(hAxes, 'CurrentPoint');
	if (pt(1,1)<lim(1)) || (pt(1,1)>lim(2)) || (pt(1,2)<lim(3)) || (pt(1,2)>lim(4)),	return,		end
	old_pt = get(h,'UserData');
	xx = get(h,'XData');
	dx = pt(1,1) - old_pt(1);
	xx = xx + dx;
	if (xx(1) < -180)
		xx(1) = -180;		xx(3) = -180 + width;
		xx(2) = xx(1);		xx(5) = xx(1);		xx(4) = xx(3);
	elseif (xx(3) > 360)
		xx(1) = 360 - width;	xx(3) = 360;
		xx(2) = xx(1);			xx(5) = xx(1);	xx(4) = xx(3);
	end
	set(h,'UserData',[pt(1,1) pt(1,2)])
	set(h, 'XData',xx);
	set(hXmin,'String', sprintf('%.6f',xx(1))),		set(hXmax,'String', sprintf('%.6f',xx(3)))

% ---------
function wbu_MoveRect(obj,eventdata,h,state)
	uirestore_j(state, 'nochildren');           % Restore the figure's initial state
% -----------------------------------------------------------------------------------------


% --- Creates and returns a handle to the GUI figure. 
function world_is_not_round_enough_LayoutFcn(h1)
	IAmAMac = getappdata(0, 'IAmAMac');

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'MenuBar','none',...
'Name','The World is not Round Enough',...
'NumberTitle','off',...
'Position',[520 430 811 373],...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

h2 = axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'Position',[5 49 801 321],...
'Tag','axes1');

h4 = get(h2,'xlabel');

set(h4,'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[0.498751560549313 -0.0732087227414331 1.00005459937205],...
'VerticalAlignment','cap',...
'HandleVisibility','off');

h5 = get(h2,'ylabel');

set(h5,'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[-0.0355805243445693 0.496884735202492 1.00005459937205],...
'Rotation',90,...
'VerticalAlignment','bottom',...
'HandleVisibility','off');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@world_NRE_uiCB,h1,'edit_xMin_CB'},...
'Position',[389 4 91 23],...
'Style','edit',...
'Tooltip','You may try entering a West longitude here, but result is indeterminate.',...
'Tag','edit_xMin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[550 4 91 23],...
'Style','edit',...
'Tag','edit_xMax');

uicontrol('Parent',h1,...
'Position',[349 6 35 15],...
'String','X min',...
'Style','text');

uicontrol('Parent',h1,...
'Position',[510 8 35 15],...
'String','X max',...
'Style','text');

uicontrol('Parent',h1,...
'Call',{@world_NRE_uiCB,h1,'push_apply_CB'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[715 3 90 21],...
'String','Apply',...
'Tag','push_apply');

if (IAmAMac),	pos = [10 3 125 30];
else			pos = [10 3 120 26];
end
uicontrol('Parent',h1,...
'Call',{@world_NRE_uiCB,h1,'push_to360_CB'},...
'FontName','Helvetica',...
'Position',pos,...
'String','[-180 180] -> [0 360]',...
'Tag','push_to360');

if (IAmAMac),	pos = [161 3 130 30];
else			pos = [161 3 120 26];
end
uicontrol('Parent',h1,...
'Call',{@world_NRE_uiCB,h1,'push_to180_CB'},...
'FontName','Helvetica',...
'Position',pos,...
'String','[0 360] -> [-180 180]',...
'Tag','push_to180');


function world_NRE_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

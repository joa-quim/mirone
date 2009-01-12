function varargout = world_is_not_round_enough(varargin)
% Helper tool to convert between-to [-180 180] <-> [0 360] longitude ranges

	if (numel(varargin) == 0),		return,		end
 
	hObject = figure('Tag','figure1','Visible','off');
	world_is_not_round_enough_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'center')

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
	d_path = handMir.path_data;
	handles.dx = diff(handles.head(1:2));
	handles.square = [handles.head(1) handles.head(3);  handles.head(1) handles.head(4); ...
					handles.head(2) handles.head(4); handles.head(2) handles.head(3); handles.head(1) handles.head(3)];
	
	if (handles.dx < (360 - handles.head(8)/2))
		set([handles.push_to360 handles.push_to180],'Enable','inactive')
		str = sprintf(['This button is blocked because the map is not global\n' ...
						'That is, it does not cover 360 degrees of longitude\n' ...
						'You can still click and drag the rectangle and try luck.']);
		set([handles.push_to360 handles.push_to180],'Enable','inactive','Tooltip', str)
	else
		set(handles.push_to360,'Tooltip', 'Convert map longitudes from [-180 180] to [0 360].')
		set(handles.push_to180,'Tooltip', 'Convert map longitudes from [0 360] to [-180 180].')
	end
	set(handles.edit_xMin,'String', sprintf('%.6f',handles.square(1,1)) )
	set(handles.edit_xMax,'String', sprintf('%.6f',handles.square(3,1)) )

	% Import the world map
	handles.w_map = flipdim(imread([d_path 'logo_etopo2.jpg']),1);
	siza = size(handles.w_map);
	handles.w_map = [handles.w_map handles.w_map(:,1:fix(siza(2)/2)+1,:) ];
	
	image([-180 360],[-90 90],handles.w_map);   set(handles.axes1,'YDir','normal');
	handles.hRectLims = patch('Parent',handles.axes1, 'xdata', handles.square(:,1), 'ydata', handles.square(:,2), ...
		'LineWidth',3, 'FaceColor','none');
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
function edit_xMin_Callback(hObject, eventdata, handles)
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
function push_to360_Callback(hObject, eventdata, handles)
	rect_x = get(handles.hRectLims, 'xdata');
	rect_x = rect_x + 180;
	set(handles.hRectLims, 'XData',rect_x);
	push_apply_Callback(handles.push_apply, eventdata, handles)

% -------------------------------------------------------------------------------------
function push_to180_Callback(hObject, eventdata, handles)
	rect_x = get(handles.hRectLims, 'xdata');
	rect_x = rect_x - 180;
	set(handles.hRectLims, 'XData',rect_x);
	push_apply_Callback(handles.push_apply, eventdata, handles)

% -------------------------------------------------------------------------------------
function push_apply_Callback(hObject, eventdata, handles)

	rect_x = get(handles.hRectLims, 'xdata');
	eps_x = handles.head(8);
	XYlim = getappdata(handles.hMirAxes,'ThisImageLims');
	if ( (abs(rect_x(1) - XYlim(1)) < eps_x) && (abs(rect_x(2) - XYlim(2)) < eps_x) )		% Nothing changed, so by by
		return
	end

	dy = diff(handles.square(1:2,2));
	x_min = rect_x(1);		x_max = rect_x(3);
	handMir = guidata(handles.hMirFig);			% Get the Mirone handles
	if (handles.validGrid)
		[X,Y,Z,head] = load_grd(handMir);
		if (isempty(Z)),	return,		end
	else
		Z = get(handles.hMirImg,'CData');
	end

	Z = to_from_180(handles, handMir, Z, x_min, x_max, dy, eps_x);

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
	handMir.head(1:2) = [X(1) X(end)];
	guidata(handles.hMirFig, handMir)

	set(handles.hMirImg,'CData', Z, 'XData', X);
	setappdata(handles.hMirAxes,'ThisImageLims',[get(handles.hMirAxes,'XLim') get(handles.hMirAxes,'YLim')])

% -----------------------------------------------------------------------------------------
function Z = to_from_180(handles, handMir, Z, x_min, x_max, dy, eps_x)

	if ( x_max > 0 && (abs(x_min - 180) < eps_x) )			% 360 -> 180
		[Z_rect,r_c] = cropimg(handles.head(1:2),handles.head(3:4),Z,[x_min handles.square(1,2) 180 dy],'out_grid');
		Z = [Z_rect Z(:, 1:(size(Z,2)-size(Z_rect,2)), :)];
		handMir.geog = 1;		guidata(handMir.figure1,handMir);
	elseif ( x_min < 180 && (abs(x_max - 360) < eps_x) )	% 180 -> 360
		[Z_rect,r_c] = cropimg(handles.head(1:2),handles.head(3:4),Z,[x_min handles.square(1,2) 180-x_min dy],'out_grid');
		Z = [Z_rect Z(:, 1:(size(Z,2)-size(Z_rect,2)), :)];
		handMir.geog = 2;		guidata(handMir.figure1,handMir);
	end

% -----------------------------------------------------------------------------------------
function moveRect(obj,eventdata, handles, hRect)
	hFig = handles.figure1;			hAxes = handles.axes1;
	hXmin = handles.edit_xMin;		hXmax = handles.edit_xMax;
	state = uisuspend_fig(hFig);            % Remember initial figure state
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
	uirestore_fig(state);           % Restore the figure's initial state
% -----------------------------------------------------------------------------------------


% --- Creates and returns a handle to the GUI figure. 
function world_is_not_round_enough_LayoutFcn(h1)

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'MenuBar','none',...
'Name','The World is not Round Enough',...
'NumberTitle','off',...
'Position',[520 430 811 370],...
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
'Callback',{@world_is_not_round_enough_uicallback,h1,'edit_xMin_Callback'},...
'Position',[389 4 91 21],...
'Style','edit',...
'TooltipString','You may try entering a West longitude here, but result is indeterminate.',...
'Tag','edit_xMin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[550 4 91 21],...
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
'Callback',{@world_is_not_round_enough_uicallback,h1,'push_apply_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[715 3 90 23],...
'String','Apply',...
'Tag','push_apply');

uicontrol('Parent',h1,...
'Callback',{@world_is_not_round_enough_uicallback,h1,'push_to360_Callback'},...
'FontName','Helvetica',...
'Position',[10 3 120 24],...
'String','[-180 180] -> [0 360]',...
'Tag','push_to360');

uicontrol('Parent',h1,...
'Callback',{@world_is_not_round_enough_uicallback,h1,'push_to180_Callback'},...
'FontName','Helvetica',...
'Position',[161 3 120 24],...
'String','[0 360] -> [-180 180]',...
'Tag','push_to180');


function world_is_not_round_enough_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

function varargout = shading_params(varargin)
% Helper window to select shading illumination parmeters 

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

	if (isempty(varargin))
		errordlg('Unknown Illumination option','Error'),	return
	end

	hObject = figure('Vis','off');
	shading_params_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center');

	mir_dirs = getappdata(0,'MIRONE_DIRS');
	if (~isempty(mir_dirs))
		f_path = [mir_dirs.home_dir filesep 'data' filesep];
	else
		f_path = [cd filesep 'data' filesep];
	end

	handles.mercedes = 0;       % Flag for False Color option  (when set to 1)
	handles.dirDerivative = 0;  % Flag for Directional Derivative option

	if (~strcmp(varargin{1},'dirDerivative'))
		load([f_path 'mirone_icons.mat'],'um_ico','dois_ico','tres_ico','quatro_ico','cinco_ico','seis_ico','sete_ico');
		h_toolbar = uitoolbar('parent',hObject,'Clipping', 'on', 'BusyAction','queue','HandleVisibility','on',...
			'Interruptible','on','Tag','FigureToolBar','Vis','on');
		handles.ui_grdgrad_A = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'grdgradient_A'},...
			'Tooltip','GMT grdgradient classic', 'CData',um_ico);
		handles.ui_grdgrad_E1 = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'grdgradient_E1'},...
			'Tooltip','GMT grdgradient Lambertian', 'CData',dois_ico);
		handles.ui_grdgrad_E2 = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'grdgradient_E2'},...
			'Tooltip','GMT grdgradient Peucker', 'CData',tres_ico);
		handles.ui_lambert = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'lambertian'},...
			'Tooltip','Lambertian with lighting', 'CData',quatro_ico);
		handles.ui_color = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'color'},...
			'Tooltip','Color (Manip-Raster)', 'CData',cinco_ico);
		handles.ui_gray = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'gray'},...
			'Tooltip','Gray (Manip-Raster)', 'CData',seis_ico);
		handles.ui_falseColor = uitoggletool('parent',h_toolbar,'Click',{@show_needed,'mercedes'},...
			'Tooltip','False color', 'CData',sete_ico);

		handles.ampFactor = 125;		% For the Mercedes illumination with Tierrie's algo

		% The following is for use in toggle_uis(...). It's easier there to deal with
		% numbers, but in other places its preferable to have names. So we have a
		% duplicate. Attention, the order in .ui_tools must reproduce its declarations
		handles.ui_tools = [handles.ui_grdgrad_A handles.ui_grdgrad_E1 handles.ui_grdgrad_E2 ...
			handles.ui_lambert handles.ui_color handles.ui_gray handles.ui_falseColor];
	else
		% With this option the better is realy not to show the rest of the window elements
		pos_f = get(hObject,'Position');			% Original Fig size
		pos_a = get(handles.axes1,'Position');		pos_a(1) = 0;
		set(handles.axes2,'Vis','off')
		pos_textAzim = get(handles.text_azim,'Position');		% Get and change text_azim position
		set(handles.text_azim,'Position',[2 pos_textAzim(2) pos_textAzim(3)-2 pos_textAzim(4)])
		pos_azim = get(handles.edit_azim,'Position');			% Get and change edit_azim position
		set(handles.edit_azim,'Position',pos_azim+[-14 0 0 0])
		pos_ok = get(handles.push_OK,'Position');     pos_ok(1) = 90; pos_ok(3) = 40;    pos_ok(4) = 20;
		set(handles.push_OK,'Position',pos_ok)
		set(hObject,'Position',[pos_f(1) pos_f(2) pos_a(1)+pos_a(3)+14 pos_f(4)],'Name','')   % New figure's size
		handles.dirDerivative = 1;
	end

	% Import background image
	astrolabio = imread([f_path 'astrolabio.jpg']);
	image(astrolabio,'parent',handles.axes1);

	pos = get(handles.axes1,'Position');
	set(handles.axes1,'Vis','off')

	% Draw everything that may be needed for all options. Later, depending on the
	% option selected, only the allowed features will be made visible
	x0 = pos(3)/2;      y0 = pos(4)/2;      radius = pos(3)/2;
	h_line(1) = line('parent',handles.axes1,'XData',[x0 x0],'YData',[y0 0],'Color','r','Tag','red','LineWidth',3,'Userdata',radius);
	if (handles.dirDerivative == 0)         % Otherwise there is no point in creating those
		x1 = x0 + radius * cos(30*pi/180);      y1 = y0 + radius * sin(30*pi/180);
		h_line(2) = line('parent',handles.axes1,'XData',[x0 x1],'YData',[y0 y1],'Color','g','Tag','green','LineWidth',3,'Vis','off');
		x1 = x0 + radius * cos(150*pi/180);     y1 = y0 + radius * sin(150*pi/180);
		h_line(3) = line('parent',handles.axes1,'XData',[x0 x1],'YData',[y0 y1],'Color','b','Tag','blue','LineWidth',3,'Vis','off');
		set(h_line,'Userdata',radius)        % save radius of circumscribed circle (image is square)
		% Now draw, on axes2, a quarter of circle and a line
		t = 0:0.02:pi/2;    x = [0 cos(t) 0];     y = [0 sin(t) 0];
		line('parent',handles.axes2,'XData',x,'YData',y,'HitTest','off','Color','k','LineWidth',1);
		h_line(4) = line('parent',handles.axes2,'XData',[0 cos(30*pi/180)],'YData',[0 sin(30*pi/180)],'Color','k','LineWidth',3,'Vis','off');
		set(h_line(4),'Tag','Elev','Userdata',1)        % save radius of circumscribed circle
	end

	switch varargin{1}
		case 'grdgradient_A'
			set(handles.ui_grdgrad_A,'State','on');
		case 'dirDerivative'
			% Do nothing here. Just to account for this option            
		otherwise
			errordlg('Unknown Illumination option','Error')
	end

	handles.h_line = h_line;
	guidata(hObject, handles);
	show_needed(hObject,[],varargin{1})
	set(hObject,'WindowButtonDownFcn',{@ButtonDown,h_line,handles});

	% Choose default command line output for shading_params_export
	handles.output = hObject;
	guidata(hObject, handles);

	set(hObject,'Vis','on');
	% UIWAIT makes shading_params_export wait for user response (see UIRESUME)
	uiwait(handles.figure1);

	handles = guidata(hObject);
	varargout{1} = handles.output;
	delete(handles.figure1);
	drawnow		% Force killing to happen before Mirone takes over and Matlab decides to do what it pleases.

% -----------------------------------------------------------------------------------------
function show_needed(obj,eventdata,opt)
	handles = guidata(obj);         % Get handles
	h_all = handles.h_line;
	handles.mercedes = 0;
	if (strncmp(opt,'grdgradient',11))
		set(handles.edit_elev,'Enable','off');		set(handles.edit_azim,'Vis','on')
		set(handles.edit_azimR,'Vis','off');		set(handles.edit_azimG,'Vis','off')
		set(handles.edit_azimB,'Vis','off');		set(handles.text_elev,'Enable','on');
		if (strcmp(opt(12:end),'_A'))
			set(handles.edit_azim,'Enable','on');	set(handles.text_azim,'Enable','on');
			set(h_all(1),'Vis','on');				set(h_all(2:4),'Vis','off')
			toggle_uis(handles,1);					set(handles.figure1,'Name','GMT grdgradient')
		elseif (strcmp(opt(12:end),'_E1'))
			set(h_all([1 4]),'Vis','on');			set(h_all(2:3),'Vis','off')
			set(handles.edit_azim,'Enable','on');	set(handles.text_azim,'Enable','on');
			set(handles.edit_elev,'Enable','on');	toggle_uis(handles,2);
			set(handles.figure1,'Name','GMT grdgradient - Lambertian')
		else        % _E2
			set(handles.edit_azim,'Enable','off');	set(handles.text_azim,'Enable','off');
			set(handles.text_elev,'Enable','off');
			set(h_all(1:4),'Vis','off');			toggle_uis(handles,3);
			set(handles.figure1,'Name','GMT grdgradient - Peucker')
		end
	elseif (strcmp(opt,'color') || strcmp(opt,'gray'))
		set(handles.edit_elev,'Enable','on');		set(handles.edit_azim,'Vis','on')
		set(handles.edit_azimR,'Vis','off');		set(handles.edit_azimG,'Vis','off')
		set(handles.edit_azimB,'Vis','off');
		set(handles.text_elev,'Enable','on');
		set(handles.edit_azim,'Enable','on');		set(handles.text_azim,'Enable','on');
		set(h_all(1),'Vis','on');					set(h_all(4),'Vis','on')
		set(h_all(2:3),'Vis','off')
		if (strcmp(opt,'color'))
			toggle_uis(handles,5);					set(handles.figure1,'Name','Color')
		else
			toggle_uis(handles,6);					set(handles.figure1,'Name','Gray')
		end
	elseif (strcmp(opt,'lambertian'))
		set(handles.edit_elev,'Enable','on');		set(handles.edit_azim,'Vis','on')
		set(handles.edit_azimR,'Vis','off');		set(handles.edit_azimG,'Vis','off')
		set(handles.edit_azimB,'Vis','off');
		set(handles.text_elev,'Enable','on');
		set(handles.edit_azim,'Enable','on');		set(handles.text_azim,'Enable','on');
		set(h_all(1),'Vis','on');					set(h_all(4),'Vis','on')
		set(h_all(2:3),'Vis','off')
		toggle_uis(handles,4)
		set(handles.figure1,'Name','Lambertian lighting')
	elseif (strcmp(opt,'mercedes'))
		set(handles.edit_elev,'Enable','on');		set(handles.edit_azim,'Vis','off')
		set(handles.edit_azimR,'Vis','on');			set(handles.edit_azimG,'Vis','on')
		set(handles.edit_azimB,'Vis','on');
		set(handles.text_elev,'Enable','on');
		set(h_all(1:4),'Vis','on')
		handles.mercedes = 1;
		toggle_uis(handles,7)
		set(handles.figure1,'Name','False color')
	elseif (strcmp(opt,'dirDerivative'))			% This is for good because this function won't be called again
		%set(h_all(1),'Vis','on');                   %set(h_all(2:4),'Vis','off')
		set(handles.edit_azimR,'Vis','off');		set(handles.edit_azimG,'Vis','off')
		set(handles.edit_azimB,'Vis','off');		set(handles.text_elev,'Enable','off');
		set(handles.figure1,'Name','Azim')
	end
	guidata(obj,handles)

% -----------------------------------------------------------------------------------------
function toggle_uis(handles,ui)
% Do not let more the one uitoggletool be on the state of pushed
	n = 1:numel(handles.ui_tools);
	n(n == ui) = [];        % Remove current ui index
	set(handles.ui_tools(n),'State','off');
	
	if (ui == 4),	vis = 'on';
	else			vis = 'off';
	end
	set([handles.edit_ambient handles.edit_diffuse handles.edit_specular handles.edit_shine],'Vis',vis);
	set([handles.text_ambient handles.text_diffuse handles.text_reflection handles.text_shine],'Vis',vis);

	if (ui == 7)			% Mercedes
		set([handles.radio_oldAlgo handles.radio_grdgrad], 'Vis', 'on')
		set(handles.text_diffuse, 'Vis','on', 'String','Amp factor')
		set(handles.edit_diffuse, 'Vis','on', 'String',sprintf('%d',handles.ampFactor))
	else
		set([handles.radio_oldAlgo handles.radio_grdgrad], 'Vis', 'off')
		set(handles.text_diffuse, 'String','Diffuse reflection')
		set(handles.edit_diffuse, 'String','0.6')
	end
	handles.current_method = ui;
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function ButtonDown(obj,eventdata,h_all,handles)
% It could be cleverer.
	pt = get(gca, 'CurrentPoint');
	x_lim = get(gca,'xlim');      y_lim = get(gca,'ylim');
	% check if x,y is inside of axis
	if ~((pt(1,1)>=x_lim(1)) && (pt(1,1)<=x_lim(2)) && (pt(1,2)>=y_lim(1)) && (pt(1,2)<=y_lim(2)))    % outside axis limits
		return
	end
	if any(h_all == gco)
		h = h_all(h_all == gco);    % When more than one line handle exists, find only the selected one
		set(gcf,'WindowButtonMotionFcn',{@ButtonMotion,h,handles},'WindowButtonUpFcn',{@ButtonUp,h_all,handles},...
			'Pointer', 'crosshair');
	else
		return
	end

% -----------------------------------------------------------------------------------------
function ButtonMotion(obj,eventdata,h,handles)
	selectionType = get(gcf, 'SelectionType');
	pt = get(gca, 'CurrentPoint');
	if strcmp(selectionType, 'normal')      % right-cick
		xx = get(h,'XData');    yy = get(h,'YData');
		theta = cart2pol(pt(1,1)-xx(1),pt(1,2)-yy(1));
		radius = get(h,'Userdata');
		x2 = xx(1) + radius * cos(theta);      y2 = yy(1) + radius * sin(theta);
		if strcmp(get(h,'Tag'),'Elev') && (theta >= 0 && theta <= pi/2)   % Elevation line
			set(h,'XData',[xx(1) x2],'YData',[yy(1) y2]);
			set(handles.edit_elev,'String',num2str(fix(theta *180/pi)) )
		elseif ~strcmp(get(h,'Tag'),'Elev')     % Azimuth line(s)
			set(h,'XData',[xx(1) x2],'YData',[yy(1) y2]);

			% NOTE to if I ever want to reuse this code. Normally ang_2pi should be = pi/2 - (pi*.....)
			% for the normal y origin at bottm left corner. However, due to the stupid habit of using y=0
			% at top left corner when dealing with images, to get an azimuth angle we have to do like following. 

			% truncate angles into [-pi pi] range
			ang_2pi = pi/2 + ( pi*((abs(theta)/pi) - 2*ceil(((abs(theta)/pi)-1)/2)) * sign(theta) );
			epsilon = -1e-7;        %  Allow points near zero to remain there
			indx = find(ang_2pi < epsilon);
			%  Shift the points in the [-pi 0] range to [pi 2pi] range
			if ~isempty(indx);  ang_2pi(indx) = ang_2pi(indx) + 2*pi;  end;
			if strcmp(get(h,'Tag'),'red')
				if (~handles.mercedes)
					set(handles.edit_azim,'String',num2str(fix(ang_2pi *180/pi)) )
				else
					set(handles.edit_azimR,'String',num2str(fix(ang_2pi *180/pi)) )
				end
			elseif strcmp(get(h,'Tag'),'green')
				set(handles.edit_azimG,'String',num2str(fix(ang_2pi *180/pi)) )
			elseif strcmp(get(h,'Tag'),'blue')
				set(handles.edit_azimB,'String',num2str(fix(ang_2pi *180/pi)) )
			end
		end
	end

% -----------------------------------------------------------------------------------------
function ButtonUp(obj,eventdata,h,handles)
	set(handles.figure1,'WindowButtonMotionFcn','','WindowButtonDownFcn',{@ButtonDown,h,handles},'WindowButtonUpFcn','');
	set(handles.figure1,'Pointer', 'arrow')

% ---------------------------------------------------------------------
function radio_oldAlgo_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_grdgrad, 'Val', 0)
	set(handles.h_line(4),'Vis','on')
	set(handles.edit_elev,'Enable','on');

% ---------------------------------------------------------------------
function radio_grdgrad_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_oldAlgo, 'Val', 0)
	set(handles.h_line(4),'Vis','off')
	set(handles.edit_elev,'Enable','off');

% -----------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
	if (~handles.mercedes)
		out.azim = str2double(get(handles.edit_azim,'String'));
	else
		out.azim(1) = str2double(get(handles.edit_azimR,'String'));
		out.azim(2) = str2double(get(handles.edit_azimG,'String'));
		out.azim(3) = str2double(get(handles.edit_azimB,'String'));
	end
	out.elev = str2double(get(handles.edit_elev,'String'));
	out.ambient = str2double(get(handles.edit_ambient,'String'));
	out.diffuse = str2double(get(handles.edit_diffuse,'String'));
	out.specular = str2double(get(handles.edit_specular,'String'));
	out.shine = str2double(get(handles.edit_shine,'String'));

	if (~handles.dirDerivative)
		% Find out which illumination model has been choosen. 
		% This is needed in Mirone to know what to do with the out vars
		if (strcmp(get(handles.ui_grdgrad_A,'State'),'on'))
			out.illum_model = 1;
		elseif (strcmp(get(handles.ui_grdgrad_E1,'State'),'on'))
			out.illum_model = 2;
		elseif (strcmp(get(handles.ui_grdgrad_E2,'State'),'on'))
			out.illum_model = 3;
		elseif (strcmp(get(handles.ui_lambert,'State'),'on'))
			out.illum_model = 4;
		elseif (strcmp(get(handles.ui_color,'State'),'on'))
			out.illum_model = 5;
		elseif (strcmp(get(handles.ui_gray,'State'),'on'))
			out.illum_model = 6;
		elseif (strcmp(get(handles.ui_falseColor,'State'),'on'))
			out.illum_model = 7;
			out.ampFactor = out.diffuse;
			if (get(handles.radio_oldAlgo, 'Val'))
				out.mercedes_type = 0;			% Original Manip-raster algo
			else
				out.mercedes_type = 1;			% grdgradient type
			end
		else
			errordlg('Uknown illumination model.','Error')
			out = [];
		end
	end

	handles.output = out;
	guidata(hObject,handles);
	uiresume(handles.figure1);

% --------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	if (exist('OCTAVE_VERSION','builtin'))		% To know if we are running under Octave
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	else
		do_uiresume = strcmp(get(handles.figure1, 'waitstatus'), 'waiting');
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		handles.output = [];		% User gave up, return nothing
		guidata(handles.figure1, handles);	uiresume(handles.figure1);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

% --------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata)
% Check for "escape"
	handles = guidata(hObject);
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles.output = [];    % User said no by hitting escape
		guidata(hObject, handles);
		uiresume(handles.figure1);
	end


% --- Creates and returns a handle to the GUI figure. 
function shading_params_LayoutFcn(h1)

set(h1,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','shading_params',...
'NumberTitle','off',...
'Position',[520 400 320 150],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[56 6 34 21],...
'String','0',...
'Style','edit',...
'Tooltip','Azimuth direction',...
'Tag','edit_azim');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[131 6 30 21],...
'String','30',...
'Style','edit',...
'Tooltip','Elevation light direction',...
'Tag','edit_elev');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[186 102 30 21],...
'String','.55',...
'Style','edit',...
'Tag','edit_ambient');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[186 79 30 21],...
'String','.6',...
'Style','edit',...
'Tag','edit_diffuse');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[186 56 30 21],...
'String','.4',...
'Style','edit',...
'Tag','edit_specular');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[186 33 30 21],...
'String','10',...
'Style','edit',...
'Tag','edit_shine');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[219 102 70 16],...
'String','Ambient light',...
'Style','text',...
'Tag','text_ambient');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[219 78 83 17],...
'String','Diffuse reflection',...
'Style','text',...
'Tag','text_diffuse');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[219 57 91 15],...
'String','Specular reflection',...
'Style','text',...
'Tag','text_reflection');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[219 35 77 15],...
'String','Specular shine',...
'Style','text',...
'Tag','text_shine');

uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[13 7 42 16],...
'String','Azimuth',...
'Style','text',...
'Tag','text_azim');

uicontrol('Parent',h1,...
'FontSize',9,...
'Position',[123 104 50 16],...
'String','Elevation',...
'Style','text',...
'Tag','text_elev');

axes('Parent',h1,'Units','pixels','Position',[16 29 91 91],'Tag','axes1','Vis','off');

axes('Parent',h1,...
'Units','pixels',...
'Position',[126 49 51 51],...
'Tag','axes2',...
'Vis','off');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0 0],...
'Position',[10 6 34 21],...
'String','0',...
'Style','edit',...
'Tooltip','Red component azimuth',...
'Tag','edit_azimR');

uicontrol('Parent',h1,...
'BackgroundColor',[0 1 0],...
'Position',[44 6 34 21],...
'String','120',...
'Style','edit',...
'Tooltip','Green component azimuth',...
'Tag','edit_azimG');

uicontrol('Parent',h1,...
'BackgroundColor',[0 0 1],...
'Position',[79 6 34 21],...
'String','240',...
'Style','edit',...
'Tooltip','Blue component azimuth',...
'Tag','edit_azimB');

uicontrol('Parent',h1,...
'Call',@shading_params_uiCB,...
'Position',[187 101 87 23],...
'String','Old algorithm',...
'Style','radiobutton',...
'TooltipString','Use the older algorithm ',...
'Value',1,...
'Tag','radio_oldAlgo');

uicontrol('Parent',h1,...
'Call',@shading_params_uiCB,...
'Position',[187 45 87 23],...
'String','grdgradient',...
'Style','radiobutton',...
'TooltipString','Compute illumination with grdgradient ',...
'Tag','radio_grdgrad');

uicontrol('Parent',h1,...
'Call',@shading_params_uiCB,...
'Position',[229 6 66 21],...
'String','OK',...
'Tag','push_OK');

function shading_params_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));


function varargout = manual_pole_adjust(varargin)
%	Do Euler rotations interactively

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

% $Id: $

	if (isempty(varargin))
		errordlg('MANUAL POLE ADJUST: wrong number of arguments.','Error'),	return
	end

	hObject = figure('Vis','off');
	manual_pole_adjust_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject, 'right')

	handles.hCallingFig = varargin{1};
	handMir = guidata(handles.hCallingFig);
	handles.geog = handMir.geog;
	if (handMir.no_file)
		errordlg('You didn''t even load a file. What are you expecting then?','ERROR')
		delete(hObject);    return
	end
	if (~handMir.geog)
		errordlg('This operation is currently possible only for geographic type data','ERROR')
		delete(hObject);    return
	end

	handles.hLine = [];
	handles.hLine_g = [];
	handles.line_x = [];
	handles.line_y = [];
	handles.have_pole = 0;
	handles.p_lon = 0;      % Default to this non-sense
	handles.p_lat = 0;
	handles.p_omega = 0;
	handles.h_active_line_str = findobj(handles.figure1,'Tag','text_activeLine');      % Get this handle
	handles.path_continent = [pwd filesep 'continents' filesep];

	set(handles.edit_lon, 'String', num2str(handles.p_lon))
	set(handles.edit_lat, 'String', num2str(handles.p_lat))
	set(handles.edit_omega, 'String', num2str(handles.p_omega))

	set(handles.slider_lon, 'Value', handles.p_lon)
	set(handles.slider_lat, 'Value', handles.p_lat)
	set(handles.slider_omega, 'Value', handles.p_omega)

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% --------------------------------------------------------------------------
function edit_lon_CB(hObject, handles)
	x = str2double(get(hObject,'String'));
	set(handles.slider_lon,'Value',x)
	handles.p_lon = x;
	if (isempty(handles.p_lat) || isempty(handles.p_omega))
		handles.have_pole = 0;
	else                            % OK, we have all the pole parameters
		handles.have_pole = 1;
	end
	guidata(hObject, handles);
	if (~handles.have_pole),		return,		end     % We don't still have a pole
	if (isempty(handles.hLine)),	return,		end     % There is nothing to do yet
	apply_rot(handles)

% --------------------------------------------------------------------------
function slider_lon_CB(hObject, handles)
	if (isempty(handles.hLine))
		set(hObject,'Value',0),		return
	end
	if (~handles.have_pole),		return,		end     % We don't still have a pole
	val = get(hObject,'Value');
	set(handles.edit_lon,'String',num2str(val))
	handles.p_lon = val;
	apply_rot(handles)
	guidata(hObject, handles);

% --------------------------------------------------------------------------
function edit_lat_CB(hObject, handles)
	x = str2double(get(hObject,'String'));
	set(handles.slider_lat,'Value',x)
	handles.p_lat = x;
	if (isempty(handles.p_lon) || isempty(handles.p_omega))
		handles.have_pole = 0;
	else                            % OK, we have all the pole parameters
		handles.have_pole = 1;
	end
	guidata(hObject, handles);
	if (~handles.have_pole),		return,		end		% We don't still have a pole
	if (isempty(handles.hLine)),	return,		end		% There is nothing to do yet
	apply_rot(handles)

% --------------------------------------------------------------------------
function slider_lat_CB(hObject, handles)
	if (isempty(handles.hLine))
		set(hObject,'Value',0),		return
	end
	if (~handles.have_pole),	return,		end		% We still do not have a pole
	val = get(hObject,'Value');
	set(handles.edit_lat,'String',num2str(val))
	handles.p_lat = val;
	apply_rot(handles)
	guidata(hObject, handles);

% --------------------------------------------------------------------------
function edit_omega_CB(hObject, handles)
	x = str2double(get(hObject,'String'));
	set(handles.slider_omega,'Value',x)
	handles.p_omega = x;
	if (isempty(handles.p_lon) || isempty(handles.p_lat))
		handles.have_pole = 0;
	else                            % OK, we have all the pole parameters
		handles.have_pole = 1;
	end
	guidata(hObject, handles);
	if (~handles.have_pole),		return,		end		% We don't still have a pole
	if (isempty(handles.hLine)),	return,		end		% There is nothing to do yet
	apply_rot(handles)

% --------------------------------------------------------------------------
function slider_omega_CB(hObject, handles)
	if (isempty(handles.hLine))
		set(hObject,'Value',0),		return
	end
	if (~handles.have_pole),	return,		end     % We don't still have a pole
	val = get(hObject,'Value');
	set(handles.edit_omega,'String',num2str(val))
	handles.p_omega = val;
	apply_rot(handles)
	guidata(hObject, handles);

% --------------------------------------------------------------------------
function push_polesList_CB(hObject, handles)
	fid = fopen([handles.path_continent 'lista_polos.dat'],'rt');
	c = fread(fid,'*char').';
	fclose(fid);
	s = strread(c,'%s','delimiter','\n');

	[s,v] = choosebox('Name','One Euler list',...
						'PromptString','List of poles:',...
						'SelectString','Selected poles:',...
						'ListSize',[450 300],...
						'ListString',s);

	if (v == 1)         % Finite pole
		handles.p_lon = s(1);
		handles.p_lat = s(2);
		handles.p_omega = s(3);
		set(handles.edit_lon, 'String', num2str(s(1)))
		set(handles.edit_lat, 'String', num2str(s(2)))
		set(handles.edit_omega, 'String', num2str(s(3)))
		handles.have_pole = 1;
		guidata(hObject,handles)
	else                % Stage poles or cancel
		handles.have_pole = 0;
		return;
	end

	set(handles.slider_lon, 'Value', s(1))
	set(handles.slider_lat, 'Value', s(2))
	set(handles.slider_omega, 'Value', s(3))

% --------------------------------------------------------------------------
function toggle_pickLine_CB(hObject, handles)
if (get(hObject,'Value'))
    % Test if we have potential target lines and their type
    h_mir_lines = findobj(handles.hCallingFig,'Type','line');     % Fish all objects of type line in Mirone figure
    if (isempty(h_mir_lines))                                       % We don't have any lines
        str = ['If you hited this button on purpose, than you deserve the following insult.',...
                'You #!|"*!%!?~^)--$&.',... 
                'THERE ARE NO LINES IN THAT FIGURE.'];
        errordlg(str,'Chico Clever');   set(hObject,'Value',0);     return;
    end
    % The above test is not enough. For exemple, coastlines are not eligible neither,
    % but is very cumbersome to test all the possibilities of pure non-eligible lines.
    set(handles.hCallingFig,'pointer','crosshair')
    hLine = get_polygon(handles.hCallingFig);          % Get the line handle
    if (~isempty(hLine))
        x = get(hLine,'XData');			y = get(hLine,'YData');
        handles.line_x = x(:);			handles.line_y = y(:);
        % Create a empty line handle that will hold the rotated line
        handles.hLine = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',[],'YData',[], ...
            'LineStyle','-.','LineWidth',2);
		if (get(handles.check_geodetic, 'Val'))
	        handles.hLine_g = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',[],'YData',[], ...
		        'LineStyle',':','LineWidth',2, 'Color','r');
		end
    else
        handles.line_x = [];			handles.line_y = [];
        set(hObject,'Value',0)
    end
    set(handles.hCallingFig,'pointer','arrow')
    set(hObject,'Value',0)
    set(handles.h_active_line_str,'String','GOT A LINE TO WORK WITH')
    figure(handles.figure1)                 % Bring this figure to front again
end
guidata(handles.figure1, handles);

% --------------------------------------------------------------------------
function check_geodetic_CB(hObject, handles)
	if (get(hObject, 'Val'))
		if (isempty(handles.hLine_g) || ~ishandle(handles.hLine_g))		% First time usage or line was killed
	        handles.hLine_g = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',[],'YData',[], ...
		        'LineStyle',':','LineWidth',2, 'Color','r');
			guidata(handles.figure1, handles);
		else
			set(handles.hLine_g, 'Vis', 'on');
		end
	else
		set(handles.hLine_g, 'Vis', 'off');
	end

% --------------------------------------------------------------------------
function apply_rot(handles)
	[rlon,rlat] = rot_euler(handles.line_x, handles.line_y, handles.p_lon, handles.p_lat, handles.p_omega, -1);
	if (handles.geog == 2)
		ind = (rlon < 0);
		rlon(ind) = rlon(ind) + 360;
	end
	if (ishandle(handles.hLine))	% Caution because line may have been deleted
		set(handles.hLine,'XData',rlon,'YData',rlat)
	end
	if (ishandle(handles.hLine_g))	% Compare with the geodetic latitude solution
		if (get(handles.check_geodetic, 'Val'))
			[rlon_g,rlat_g] = rot_euler(handles.line_x,handles.line_y,handles.p_lon,handles.p_lat,handles.p_omega);
			if (handles.geog == 2)
				rlon_g(ind) = rlon_g(ind) + 360;
			end
			set(handles.hLine_g,'XData',rlon_g,'YData',rlat_g)
		end
	end

% -----------------------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles = guidata(hObject);
		draw_funs(handles.hLine, 'line_uicontext')
		if (ishandle(handles.hLine_g))
			draw_funs(handles.hLine_g, 'line_uicontext')
		end
		delete(handles.figure1);
	end

% -----------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	draw_funs(handles.hLine, 'line_uicontext')
	if (ishandle(handles.hLine_g))
		draw_funs(handles.hLine_g, 'line_uicontext')
	end
	delete(handles.figure1);

% --- Creates and returns a handle to the GUI figure. 
function manual_pole_adjust_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'MenuBar','none',...
'Name','Manual pole adjust',...
'NumberTitle','off',...
'Position',[520 660 650 140],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[558 73 81 21],...
'BackgroundColor',[1 1 1],...
'Call',@manual_pole_adjust_uiCB,...
'Style','edit',...
'Tag','edit_lon');

uicontrol('Parent',h1, 'Position',[78 75 471 18],...
'BackgroundColor',[0.9 0.9 0.9],...
'Call',@manual_pole_adjust_uiCB,...
'Max',360,...
'Min',-180,...
'Style','slider',...
'SliderStep',[0.001 0.01],...
'Tag','slider_lon');

uicontrol('Parent',h1, 'Position',[558 42 81 21],...
'BackgroundColor',[1 1 1],...
'Call',@manual_pole_adjust_uiCB,...
'Style','edit',...
'Tag','edit_lat');

uicontrol('Parent',h1, 'Position',[78 44 471 18],...
'BackgroundColor',[0.9 0.9 0.9],...
'Call',@manual_pole_adjust_uiCB,...
'Max',90,...
'Min',-90,...
'Style','slider',...
'SliderStep',[0.001 0.01],...
'Tag','slider_lat');

uicontrol('Parent',h1, 'Position',[558 13 81 21],...
'BackgroundColor',[1 1 1],...
'Call',@manual_pole_adjust_uiCB,...
'Style','edit',...
'Tag','edit_omega');

uicontrol('Parent',h1, 'Position',[78 15 471 18],...
'BackgroundColor',[0.9 0.9 0.9],...
'Call',@manual_pole_adjust_uiCB,...
'Max',90,...
'Min',-90,...
'Style','slider',...
'SliderStep',[0.0005 0.005],...
'Tag','slider_omega');

uicontrol('Parent',h1, 'Position',[10 75 55 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Longitude',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 45 55 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Latitude',...
'Style','text');

uicontrol('Parent',h1, 'Position',[11 17 55 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Angle',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 107 121 21],...
'Call',@manual_pole_adjust_uiCB,...
'String','Pick line from Figure',...
'Tooltip','Allows you to mouse select one line from a Mirone figure',...
'Tag','toggle_pickLine');

uicontrol('Parent',h1, 'Position',[150 108 120 21],...
'Call',@manual_pole_adjust_uiCB,...
'String','Geodetic',...
'Tooltip',sprintf('Plot a second line (red dotted)\nwith the geodetic latitudes solution.'),...
'Style','checkbox',...
'Tag','check_geodetic');

uicontrol('Parent',h1, 'Position',[520 107 121 21],...
'Call',@manual_pole_adjust_uiCB,...
'String','Poles selector',...
'Tooltip','Select a pole from the default list',...
'Tag','push_polesList');

uicontrol('Parent',h1, 'Position',[240 109 231 17],...
'FontSize',10,...
'String','NO ACTIVE LINE',...
'Style','text',...
'Tag','text_activeLine');

function manual_pole_adjust_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

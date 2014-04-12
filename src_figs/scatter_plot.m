function varargout = scatter_plot(varargin)
% Plot scaled symbols. Scaling info (and optional color info) are provided in file's col 4 (5-7)
%
% Optional column 4 hold the individual symbol size
% Optional columns 5-7 hold  individual symbol color
% This function is still somewhat inefficient as it always plot each symbol as a different line

%	Copyright (c) 2004-2014 by J. Luis
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
		errordlg('SCATTER PLOT: wrong number of arguments.','Error'),	return
	end
 
	hObject = figure('Vis','off');
	scatter_plot_LayoutFcn(hObject);
	handles = guihandles(hObject);
 
	handMir = varargin{1};
	handles.handMir = handMir;
	handles.hCallingFig = handMir.figure1;        % This is the Mirone's fig handle
	handles.hCallingAxes = handMir.axes1;
	fname = varargin{2};

	% Position this figure glued to the right of calling figure
	posThis = get(hObject,'Pos');
	posParent = get(handles.hCallingFig,'Pos');
	ecran = get(0,'ScreenSize');
	xLL = posParent(1) + posParent(3) + 6;
	xLR = xLL + posThis(3);
	if (xLR > ecran(3))         % If figure is partially out, bring totally into screen
		xLL = ecran(3) - posThis(3);
	end
	yLL = (posParent(2) + posParent(4)/2) - posThis(4) / 2;
	set(hObject,'Pos',[xLL yLL posThis(3:4)])

	% Load one ico
	load (['data' filesep 'mirone_icons.mat'],'colorwheel_ico');
	set(handles.push_color,'CData',colorwheel_ico)

	handles.hCallingAxes = handMir.axes1;
	handles.no_file = handMir.no_file;
	handles.symbSIZES = [];
	handles.symbCOR = [];
	handles.hSymbs = [];
	handles.symbSYMB = 'o';     % Default symbol
	handles.lastSymbSize = 7;   % Default symbol size
	handles.edit_symbScale = 1;	% To eventualy scale the Z's and so have some control on GE cylinder heights

	numeric_data = text_read(fname,NaN,NaN);
	nCol = size(numeric_data,2);

	if (nCol < 2)
		delete(hObject)
		errordlg('SCATTER PLOT: Input file does not have even 2 columns.','Error')
		return
	elseif (nCol == 2)		% A file with 2 cols should never had landed here. Pass it to load_xyz 
		delete(hObject)
		load_xyz(handMir, fname, 'AsPoint')
		return
	end
	
	if (nCol >= 4)			% 4th column must contain size
		handles.symbSIZES = numeric_data(:,4);
		set(handles.popup_symbSize,'Enable','off')
		set(handles.edit_symbSize,'Enable','off')
	end

	if (nCol == 7)			% 5-7 columns must contain RGB color [0 1]
		cmin = min(min(numeric_data(:,5:7)));
		cmax = max(max(numeric_data(:,5:7)));
		if (cmin >= 0 && cmax <= 255)       % Color given in the [0 255] range
			handles.symbCOR = numeric_data(:,5:7)/255;
			set(handles.push_color,'Enable','off')
		elseif (cmin >= 0 && cmax <= 1)     % Color given in the [0 1] range
			handles.symbCOR = numeric_data(:,5:7);
			set(handles.push_color,'Enable','off')
		else
			warndlg('Color information given is erroneous and will therefore be ignored. It must be either in the [0 1] or [0 255] intervals','Warnerror')
		end
	end

	handles.symbXYZ = numeric_data(:,1:3);

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hCallingFig,'dependentFigs',plugedWin);

	set(hObject,'Vis','on')
	guidata(hObject, handles);
	if (nargout),	varargout{1} = 	hObject;	end

% ------------------------------------------------------------------------------------
function popup_symbol_CB(hObject, handles)
% Select a symbol
%'plus sign' 'circle' 'asterisk' 'cross' 'square' 'diamond' 'upward triangle'
%'downward triangle' 'right triangle' 'left triangle' 'five-pointed star' 'six-pointed star'
	symbs = '+o*xsd^v><ph';
	handles.symbSYMB = symbs(get(hObject,'Value'));
	guidata(handles.figure1,handles)
    
% ------------------------------------------------------------------------------------
function popup_symbSize_CB(hObject, handles)
	contents = get(hObject,'String');
	s = contents{get(hObject,'Value')};
	ss = str2double(s);         % 2 cases may occur here
	if (~isnan(ss))             % its NaN when 'other' was selected
		set(handles.edit_symbSize,'String',ss)
	end

% ------------------------------------------------------------------------------------
function edit_symbSize_CB(hObject, handles)
	ss = fix(str2double(get(hObject,'String')));
	if (isnan(ss))      % User stupidity
		set(hObject,'String',handles.lastSymbSize)
		return
	end

	% Now see if this is an already avilable size in which case it updates the popup
	avail = [4 5 6 7 8 9 10 11 12 14 18 24];
	ind = find(avail == ss);
	if (~isempty(ind))
		set(handles.popup_symbSize,'Val',ind)
	else
		avail = sort([avail ss]);
		avail = {avail};
		avail{end+1} = 'other';
		set(handles.popup_symbSize,'String',avail)
	end
	handles.lastSymbSize = ss;
	guidata(handles.figure1,handles)
    
% ------------------------------------------------------------------------------------
function edit_symbScale_CB(hObject, handles)
% Get a new scale factor that will be applied to the Zs and take effect only on GE plots
	ss = str2double(get(hObject,'String'));
	if (isnan(ss))      % User stupidity
		set(hObject,'String',handles.edit_symbScale)
		return
	end
	handles.edit_symbScale = ss;
	guidata(handles.figure1,handles)

% ------------------------------------------------------------------------------------
function push_plot_CB(hObject, handles)
% 
	if (handles.no_file)			% Start empty but below we'll find the true data region
		geog = 1;					% Not important. It will be confirmed later
		XMin = min(handles.symbXYZ(:,1));       XMax = max(handles.symbXYZ(:,1));        
		YMin = min(handles.symbXYZ(:,2));       YMax = max(handles.symbXYZ(:,2));
		xx = [XMin XMax];           yy = [YMin YMax];
		region = [xx yy];			% 1 stands for geog but that will be confirmed later
		mirone('FileNewBgFrame_CB', handles, [region geog])   % Create a background
	else							% Reading over an established region
		XYlim = getappdata(handles.hCallingAxes,'ThisImageLims');
		xx = XYlim(1:2);            yy = XYlim(3:4);
	end

	% Get rid of points that are outside the map limits
	[x,y,indx,indy] = aux_funs('in_map_region',handles.handMir,handles.symbXYZ(:,1),handles.symbXYZ(:,2),0,[xx yy]);
	if (isempty(x))
		warndlg('There are no points inside the current Window limits.','Warning');
		return
	end
	z = handles.symbXYZ(:,3);
	z(indx) = [];       z(indy) = [];
	if (~isempty(handles.symbCOR))          % Don't forget the color
		zC = handles.symbCOR;
		zC(indx) = [];      zC(indy) = [];
	end

	ind_NaN = isnan(z);
	if (any(ind_NaN))
		x(ind_NaN) = [];	y(ind_NaN) = [];	z(ind_NaN) = [];
		if (~isempty(handles.symbCOR)),		zC(ind_NaN) = [];	end
	end

	nPts = numel(x);

	if (handles.no_file)        % We need to compute the data extent in order to set the correct axes limits
		XMin = min(XMin,min(x));     XMax = max(XMax,max(x));
		YMin = min(YMin,min(y));     YMax = max(YMax,max(y));
		region = [XMin XMax YMin YMax];
		set(handles.hCallingAxes,'XLim',[XMin XMax],'YLim',[YMin YMax])
		setappdata(handles.hCallingAxes,'ThisImageLims',region)
		handles.geog = aux_funs('guessGeog',region);
		guidata(handles.hCallingFig,handles)
	end

	if (isempty(handles.symbSIZES))     % Symbol sizes was not provided in input
		ss = str2double(get(handles.edit_symbSize,'String'));
		symbSIZES = repmat(ss,nPts,1);
	else
		symbSIZES = handles.symbSIZES;
	end

	if (isempty(handles.symbCOR))
		cmap = get(handles.figure1,'ColorMap');
		Zmin = min(z);        Zmax = max(z);
		dZ = Zmax - Zmin;
		if (dZ == 0)        % Cte color
			zC = repmat(cmap(round(size(cmap,1)/2),:),nPts,1);      % Midle color
		else
			zC = round(((z - Zmin) / dZ) * (size(cmap,1)-1) + 1);
			zC = cmap(zC,:);
		end
	end

	z = abs(z);			% Currently the Z is only used (and many times badly) to make cylinders in GE
	if (handles.edit_symbScale ~= 1),	z = z * handles.edit_symbScale;		end
	h = zeros(1,nPts);
	for (k = 1:nPts)
		h(k) = line('XData',x(k),'YData',y(k),'Parent',handles.hCallingAxes, 'LineStyle','none', 'Tag','scatter_symbs',...
			'Marker',handles.symbSYMB,'Color','k','MarkerFaceColor',zC(k,:),'MarkerSize',symbSIZES(k));
		setappdata(h(k),'ZData',z(k)) 
	end
	setUIs(handles,h)

	handles.hSymbs = h;
	guidata(handles.figure1, handles)

% ------------------------------------------------------------------------------------
function push_delete_CB(hObject, handles)
% ...
	if (~isempty(handles.hSymbs))
	    delete(handles.hSymbs)
		handles.hSymbs = [];
		guidata(handles.figure1, handles)
	end

% ------------------------------------------------------------------------------------
function push_color_CB(hObject, handles)
	cmap = color_palettes;
	if (~isempty(cmap))
		set(handles.figure1,'Colormap',cmap)
	end
    
% ------------------------------------------------------------------------------------
function setUIs(handles,h)
	for (k = 1:numel(h))
		cmenuHand = uicontextmenu('Parent', handles.hCallingFig);
		set(h(k), 'UIContextMenu', cmenuHand);
		uimenu(cmenuHand, 'Label', 'Delete this', 'Callback', {@del_line,h(k)});
		uimenu(cmenuHand, 'Label', 'Delete all', 'Callback', {@del_all,handles,h(k)});
		ui_edit_polygon(h(k))
	end
    
% -----------------------------------------------------------------------------------------
function push_help_CB(hObject, handles)
msg = sprintf(['You landed here because the input file had 3 or more columns.\n' ...
		'Third column is used to attribute symbol color. That is done by\n' ...
		'scaling Z values to the current colormap. If ... to be continued\n']);
msgbox(msg,'Help')

% -----------------------------------------------------------------------------------------
function del_all(obj,eventdata,handles,h)
% Delete all objects that share the same tag of h
	tag = get(h,'Tag');
	hAll = findobj(handles.hCallingAxes,'Tag',tag);
	del_line(obj,eventdata,hAll)
    
% -----------------------------------------------------------------------------------------
function del_line(obj,eventdata,h)
% Delete symbols but before check if they are in edit mode
	for (k=1:numel(h))
		if (~isempty(getappdata(h(k),'polygon_data')))
			s = getappdata(h(k),'polygon_data');
			if strcmpi(s.controls,'on')         % Object is in edit mode, so this
				ui_edit_polygon(h(k))           % call will force out of edit mode
			end
		end
		delete(h(k));
	end

% --- Executes on key press over figure1 with no controls selected.%
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(hObject);
	end

% --- Creates and returns a handle to the GUI figure. 
function scatter_plot_LayoutFcn(h1)

set(h1, 'Position',[520 698 204 131],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Scatter plot',...
'NumberTitle','off',...
'PaperSize',[20.98404194812 29.67743169791],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[46 101 121 22],...
'BackgroundColor',[1 1 1],...
'Call',@scatter_plot_uiCB,...
'String',{'plus sign'; 'circle'; 'asterisk'; 'cross'; 'square'; 'diamond'; 'upward triangle'; 'downward triangle'; 'right triangle'; 'left triangle'; 'five-pointed star'; 'six-pointed star' },...
'Style','popupmenu',...
'Tooltip','Select symbol type',...
'Value',2,...
'Tag','popup_symbol');

uicontrol('Parent',h1, 'Position',[174 100 23 23],...
'Call',@scatter_plot_uiCB,...
'FontName','Helvetica',...
'Tooltip','Select a new colormap',...
'Tag','push_color');

uicontrol('Parent',h1, 'Position',[46 69 61 22],...
'BackgroundColor',[1 1 1],...
'Call',@scatter_plot_uiCB,...
'String',{'4'; '5'; '6'; '7'; '8'; '9'; '10'; '11'; '12'; '14'; '18'; '24'; 'other'},...
'Style','popupmenu',...
'Tooltip','Select symbols size in points',...
'Value',4,...
'Tag','popup_symbSize');

uicontrol('Parent',h1, 'Position',[107 70 41 21],...
'BackgroundColor',[1 1 1],...
'Call',@scatter_plot_uiCB,...
'String','7',...
'Style','edit',...
'Tooltip','Enter new size here if not in list',...
'Tag','edit_symbSize');

uicontrol('Parent',h1, 'Position',[7 104 38 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Symbol',...
'Style','text');

uicontrol('Parent',h1, 'Position',[8 72 38 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Size',...
'Style','text');

uicontrol('Parent',h1, 'Position',[174 69 23 23],...
'Call',@scatter_plot_uiCB,...
'FontName','Helvetica',...
'FontSize',12,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'String','?',...
'Tooltip','Help',...
'Tag','push_help');

uicontrol('Parent',h1, 'Position',[8 39 38 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Scale',...
'Style','text');

uicontrol('Parent',h1, 'Position',[46 37 41 21],...
'BackgroundColor',[1 1 1],...
'Call',@scatter_plot_uiCB,...
'String','1',...
'Style','edit',...
'Tooltip','Enter a Z scale factor but attention, this will affect only the GE plot',...
'Tag','edit_symbScale');

uicontrol('Parent',h1, 'Position',[7 8 81 21],...
'Call',@scatter_plot_uiCB,...
'FontName','Helvetica',...
'String','Plot',...
'Tag','push_plot');

uicontrol('Parent',h1, 'Position',[113 8 81 21],...
'Call',@scatter_plot_uiCB,...
'FontName','Helvetica',...
'String','Delete',...
'Tooltip','Delete all symbols previously (in this session) ploted',...
'Tag','push_delete');

function scatter_plot_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

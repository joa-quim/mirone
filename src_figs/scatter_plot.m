function varargout = scatter_plot(varargin)
% M-File changed by desGUIDE 
    
	if (isempty(varargin))
		errordlg('SCATTER PLOT: wrong number of arguments.','Error'),	return
	end
 
	hObject = figure('Tag','figure1','Visible','off');
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

    [numeric_data,multi_segs_str] = text_read(fname,NaN,NaN);
    if (~isempty(multi_segs_str))
        % Do something
    end
    nCol = size(numeric_data,2);
    
    if (nCol >= 4)      % 4th column must contain size
        handles.symbSIZES = numeric_data(:,4);
        set(handles.popup_symbSize,'Enable','off')
        set(handles.edit_symbSize,'Enable','off')
    end
    
    if (nCol == 7)      % 5-7 columns must contain RGB color [0 1]
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

	set(hObject,'Visible','on')
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function popup_symbol_Callback(hObject, eventdata, handles)
% Select a symbol
%'plus sign' 'circle' 'asterisk' 'cross' 'square' 'diamond' 'upward triangle'
%'downward triangle' 'right triangle' 'left triangle' 'five-pointed star' 'six-pointed star'
	symbs = '+o*xsd^v><ph';
	handles.symbSYMB = symbs(get(hObject,'Value'));
	guidata(handles.figure1,handles)
    
% ------------------------------------------------------------------------------------
function popup_symbSize_Callback(hObject, eventdata, handles)
	contents = get(hObject,'String');
	s = contents{get(hObject,'Value')};
	ss = str2double(s);         % 2 cases may occur here
	if (~isnan(ss))             % its NaN when 'other' was selected
		set(handles.edit_symbSize,'String',ss)
	end

% ------------------------------------------------------------------------------------
function edit_symbSize_Callback(hObject, eventdata, handles)
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
function push_plot_Callback(hObject, eventdata, handles)
% 
	if (handles.no_file)        % Start empty but below we'll find the true data region
        XMin = 1e50;            XMax = -1e50;    YMin = 1e50;            YMax = -1e50;
        geog = 1;               % Not important. It will be confirmed later
        XMin = min(handles.symbXYZ(:,1));       XMax = max(handles.symbXYZ(:,1));
        YMin = min(handles.symbXYZ(:,2));       YMax = max(handles.symbXYZ(:,2));
        xx = [XMin XMax];           yy = [YMin YMax];
        region = [xx yy];           % 1 stands for geog but that will be confirmed later
        mirone('FileNewBgFrame_CB', handles, [region geog])   % Create a background
	else                        % Reading over an established region
		XYlim = getappdata(handles.hCallingAxes,'ThisImageLims');
		xx = XYlim(1:2);            yy = XYlim(3:4);
	end

	% Get rid of points that are outside the map limits
	[x,y,indx,indy] = aux_funs('in_map_region',handles.handMir,handles.symbXYZ(:,1),handles.symbXYZ(:,2),0.5,[xx yy]);
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
	h = zeros(1,nPts);
	for (k=1:nPts)
		h(k) = line('XData',x(k),'YData',y(k),'ZData',z(k),'Parent',handles.hCallingAxes,'Tag','scatter_symbs',...
			'Marker',handles.symbSYMB,'Color','k','MarkerFaceColor',zC(k,:),'MarkerSize',symbSIZES(k));
	end
	setUIs(handles,h)

	handles.hSymbs = h;
	guidata(handles.figure1, handles)

% ------------------------------------------------------------------------------------
function push_delete_Callback(hObject, eventdata, handles)
    delete(handles.hSymbs)

% ------------------------------------------------------------------------------------
function push_color_Callback(hObject, eventdata, handles)
    cmap = color_palettes;
    if (~isempty(cmap))
        set(handles.figure1,'Colormap',cmap)
    end
    
% ------------------------------------------------------------------------------------
function setUIs(handles,h)
    for k=1:numel(h)
        cmenuHand = uicontextmenu;
        set(h(k), 'UIContextMenu', cmenuHand);
        uimenu(cmenuHand, 'Label', 'Delete this', 'Callback', {@del_line,h(k)});
        uimenu(cmenuHand, 'Label', 'Delete all', 'Callback', {@del_all,handles,h(k)});
        ui_edit_polygon(h(k))
    end
    
% -----------------------------------------------------------------------------------------
function push_help_Callback(hObject, eventdata, handles)
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
function scatter_plot_LayoutFcn(h1);

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Scatter plot',...
'NumberTitle','off',...
'PaperSize',[20.98404194812 29.67743169791],...
'Position',[520 698 204 101],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@scatter_plot_uicallback,h1,'popup_symbol_Callback'},...
'Position',[46 71 121 22],...
'String',{  'plus sign'; 'circle'; 'asterisk'; 'cross'; 'square'; 'diamond'; 'upward triangle'; 'downward triangle'; 'right triangle'; 'left triangle'; 'five-pointed star'; 'six-pointed star' },...
'Style','popupmenu',...
'TooltipString','Select symbol type',...
'Value',2,...
'Tag','popup_symbol');

uicontrol('Parent',h1,...
'Callback',{@scatter_plot_uicallback,h1,'push_color_Callback'},...
'FontName','Helvetica',...
'Position',[174 70 23 23],...
'TooltipString','Select a new colormap',...
'Tag','push_color');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@scatter_plot_uicallback,h1,'popup_symbSize_Callback'},...
'Position',[46 39 61 22],...
'String',{  '4'; '5'; '6'; '7'; '8'; '9'; '10'; '11'; '12'; '14'; '18'; '24'; 'other' },...
'Style','popupmenu',...
'TooltipString','Select symbols size in points',...
'Value',4,...
'Tag','popup_symbSize');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@scatter_plot_uicallback,h1,'edit_symbSize_Callback'},...
'Position',[106 40 41 21],...
'String','7',...
'Style','edit',...
'TooltipString','Enter new size here if not in list',...
'Tag','edit_symbSize');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Position',[7 74 38 15],...
'String','Symbol',...
'Style','text',...
'Tag','text1');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Position',[8 42 38 15],...
'String','Size',...
'Style','text',...
'Tag','text2');

uicontrol('Parent',h1,...
'Callback',{@scatter_plot_uicallback,h1,'push_plot_Callback'},...
'FontName','Helvetica',...
'Position',[7 8 81 21],...
'String','Plot',...
'Tag','push_plot');

uicontrol('Parent',h1,...
'Callback',{@scatter_plot_uicallback,h1,'push_delete_Callback'},...
'FontName','Helvetica',...
'Position',[113 8 81 21],...
'String','Delete',...
'TooltipString','Delete all symbols previously (in this session) ploted',...
'Tag','push_delete');

uicontrol('Parent',h1,...
'Callback',{@scatter_plot_uicallback,h1,'push_help_Callback'},...
'FontName','Helvetica',...
'FontSize',12,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[174 39 23 23],...
'String','?',...
'TooltipString','Help',...
'Tag','push_help');

function scatter_plot_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

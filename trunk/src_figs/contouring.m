function varargout = contouring(varargin)
% Front end to contour selection

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

	if isempty(varargin)		return,		end

	hObject = figure('Vis','off');
	contouring_LayoutFcn(hObject);
	handles = guihandles(hObject);

	handles.Zmin = [];
	handles.Zmax = [];
	handles.Zinc = [];
	handles.Zstart = [];
	handles.Zsingle = [];
	handles.selected_val = [];

	handles.hCallingFig = varargin{1};
	handles.head = varargin{2};
	old_cont = varargin{3};
	handles.Xmin = handles.head(1);    handles.Xmax = handles.head(2);
	handles.Ymin = handles.head(3);    handles.Ymax = handles.head(4);
	handles.Zmin = handles.head(5);    handles.Zmax = handles.head(6);
	set(handles.edit_Zmin,'String',num2str(handles.Zmin))
	set(handles.edit_Zmax,'String',num2str(handles.Zmax))
	move2side(handles.hCallingFig, hObject, 'right')

	if (~isempty(old_cont))
		list = num2cell(old_cont(:));
		set(handles.listbox_ElevValues,'String',list)
	end

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, handles.txt_CE)
	%------------- END Pro look (3D) ------------------------------

	% Add this figure handle to the carra?as list
	plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hCallingFig,'dependentFigs',plugedWin);

	% Update handles structure
	guidata(hObject, handles);

	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

%-------------------------------------------------------------------------------
function push_GuessIntervals_CB(hObject, handles)
% I don't know how contourc guesses the contours but I can make it work for me
% The trick isto generate a small grid with z ranging from min to max of the
% original grid and ask for the contours.
	X = linspace(handles.Xmin,handles.Xmax,50);
	Y = linspace(handles.Ymin,handles.Ymax,50);
	zz = repmat(linspace(handles.Zmin,handles.Zmax,50), 50,1);
	c = contourc(X,Y,zz);
	limit = size(c,2);
	i = 1;  k = 1;
	while(i < limit)
		z_level(k) = c(1,i);    npoints = c(2,i);
		i = i + npoints + 1;
		k = k + 1;
	end
	list = num2cell(z_level');
	set(handles.listbox_ElevValues,'String',list)

%-------------------------------------------------------------------------------
function edit_ElevStep_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if ( isempty(xx) || isnan(xx))
		set(hObject,'String','');
		handles.Zinc = [];
	else
		handles.Zinc = xx;
	end
	guidata(hObject, handles);

%-------------------------------------------------------------------------------
function edit_StartElev_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if ( isempty(xx) || isnan(xx))
		set(hObject,'String','');
		handles.Zstart = [];
	else
		handles.Zstart = xx;
	end
	guidata(hObject, handles);

%-------------------------------------------------------------------------------
function edit_SingleElev_CB(hObject, handles)
xx = str2double(get(hObject,'String'));
if ( isempty(xx) || isnan(xx))
    set(hObject,'String','');
    handles.Zsingle = [];
else
    handles.Zsingle = xx;
end
guidata(hObject, handles);

%-------------------------------------------------------------------------------
function push_Add_CB(hObject, handles)
	if (isempty(handles.Zsingle)),   return;    end
	list = get(handles.listbox_ElevValues,'String');
	list{end+1} = get(handles.edit_SingleElev,'String');
	[B,ind] = sort(str2num(char(list)));
	list = list(ind);
	set(handles.listbox_ElevValues,'String',list)
	set(handles.edit_SingleElev,'String','')
	handles.Zsingle = [];
	guidata(hObject, handles);

%-------------------------------------------------------------------------------
function push_GenerateIntervals_CB(hObject, handles)
	if (~isempty(handles.Zstart) && ~isempty(handles.Zinc))
		list = num2cell((handles.Zstart:handles.Zinc:handles.Zmax)');
		set(handles.listbox_ElevValues,'String',list)
	else
		errordlg('Generate how? From the empty space?','Chico Clever')
	end

%-------------------------------------------------------------------------------
function listbox_ElevValues_CB(hObject, handles)
	handles.selected_val = get(hObject,'Value');
	guidata(handles.figure1, handles);

%-------------------------------------------------------------------------------
function push_DeleteSelected_CB(hObject, handles)
	if (isempty(handles.selected_val)),		return,		end
	list = get(handles.listbox_ElevValues,'String');
	if (isempty(list)),		return,		end
	this_cont = str2double(list{handles.selected_val});
	list(handles.selected_val) = [];
	set(handles.listbox_ElevValues,'String',list,'Value',1)
	handles.selected_val = [];
	guidata(handles.figure1, handles)

	% Remove this contout from the Mirone figure and update the handMir.which_cont list
	handMir = guidata(handles.hCallingFig);
	h = findobj(handMir.axes1, 'type','line','tag','contour', 'userdata', this_cont);
	if (~isempty(h))
		labHand = getappdata(h,'LabelHands');
		try		delete(labHand),	end
		delete(h)
		id = find(handMir.which_cont == h);
		if (~isempty(id))		% It should not be empty
			handMir.which_cont(id) = [];
			guidata(handMir.figure1, handMir)
		end
	end

%-------------------------------------------------------------------------------
function push_DeleteAll_CB(hObject, handles)
	set(handles.listbox_ElevValues,'String','')
	handles.selected_val = [];
	guidata(handles.figure1, handles)

	handMir = guidata(handles.hCallingFig);
	delete(findobj(handMir.axes1,'tag','contour'))		% Remove all contours from the Mirone figure
	handMir.which_cont = [];			% Reset so that next time Mirone knows it has to rebuild
	guidata(handMir.figure1, handMir)

%-------------------------------------------------------------------------------
function push_Apply_CB(hObject, handles)
	list = get(handles.listbox_ElevValues,'String');
	if (isempty(list)),  return,	end
	list = str2num(char(list{:}));
	handMir = guidata(handles.hCallingFig);
	if (get(handles.check_plotLabels, 'Val'))
		handMir.plotContourLabels = true;
	else
		handMir.plotContourLabels = false;
	end
	mirone('DrawContours_CB',handMir,list);

%-------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(hObject);
	end

	
% --- Creates and returns a handle to the GUI figure. 
function contouring_LayoutFcn(h1)

set(h1, 'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','contouring',...
'NumberTitle','off',...
'Position',[520 473 453 327],...
'Resize','off',...
'HandleVisibility', 'off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[282 38 161 281],...
'Enable','off',...
'Tag','pushbutton9');

uicontrol('Parent',h1, 'Position',[10 38 261 281],...
'Enable','off',...
'Tag','pushbutton8');

uicontrol('Parent',h1, 'Position',[40 268 91 21],...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Style','edit',...
'Tag','edit_Zmin');

uicontrol('Parent',h1, 'Position',[168 269 91 21],...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Style','edit',...
'Tag','edit_Zmax');

uicontrol('Parent',h1, 'Position',[19 272 21 15], 'String','Min:', 'Style','text');
uicontrol('Parent',h1, 'Position',[142 272 23 15], 'String','Max:', 'Style','text');

uicontrol('Parent',h1, 'Position',[21 236 90 17],...
'String','Plot Labels',...
'Style','checkbox',...
'Tooltip','Plot master contour lines',...
'Value',1,...
'Tag','check_plotLabels');

uicontrol('Parent',h1, 'Position',[20 206 241 21],...
'Call',@contouring_uiCB,...
'String','Add Common Charting Intervals',...
'Tag','push_GuessIntervals');

uicontrol('Parent',h1, 'Units','pixels', 'Position',[20 48 241 141], 'Style','frame');

uicontrol('Parent',h1, 'Position',[30 145 75 15],...
'HorizontalAlignment','left',...
'String','Single Elevation',...
'Style','text');

uicontrol('Parent',h1, 'Position',[30 119 83 15],...
'HorizontalAlignment','left',...
'String','Starting Elevation',...
'Style','text');

uicontrol('Parent',h1, 'Position',[30 89 80 15],...
'HorizontalAlignment','left',...
'String','Elevation Step',...
'Style','text');

uicontrol('Parent',h1, 'Position',[120 85 131 21],...
'BackgroundColor',[1 1 1],...
'Call',@contouring_uiCB,...
'Style','edit',...
'Tag','edit_ElevStep');

uicontrol('Parent',h1, 'Position',[120 114 131 21],...
'BackgroundColor',[1 1 1],...
'Call',@contouring_uiCB,...
'Style','edit',...
'Tag','edit_StartElev');

uicontrol('Parent',h1, 'Position',[120 143 71 21],...
'BackgroundColor',[1 1 1],...
'Call',@contouring_uiCB,...
'Style','edit',...
'Tag','edit_SingleElev');

uicontrol('Parent',h1, 'Position',[30 55 221 21],...
'Call',@contouring_uiCB,...
'String','Generate Elevation Intervals',...
'Tag','push_GenerateIntervals');

uicontrol('Parent',h1, 'Position',[200 142 51 21],...
'Call',@contouring_uiCB,...
'String','Add',...
'Tag','push_Add');

uicontrol('Parent',h1, 'Position',[31 177 105 18],...
'String','Custom Elevations',...
'Style','text',...
'Tag','txt_CE');

uicontrol('Parent',h1, 'Position',[290 108 145 180],...
'BackgroundColor',[1 1 1],...
'Call',@contouring_uiCB,...
'Max',2,...
'Style','listbox',...
'Value',1,...
'Tag','listbox_ElevValues');

uicontrol('Parent',h1, 'Position',[290 76 140 21],...
'Call',@contouring_uiCB,...
'String','Delete Selected Values',...
'Tag','push_DeleteSelected');

uicontrol('Parent',h1, 'Position',[290 47 140 21],...
'Call',@contouring_uiCB,...
'String','Delete All Values',...
'Tag','push_DeleteAll');

uicontrol('Parent',h1, 'Position',[295 293 144 15],...
'HorizontalAlignment','left',...
'String','Elevation Values in Range',...
'Style','text');

uicontrol('Parent',h1, 'Position',[363 6 80 21],...
'Call',@contouring_uiCB,...
'String','Apply',...
'Tag','push_Apply');

uicontrol('Parent',h1, 'Position',[20 294 90 15],...
'HorizontalAlignment','left',...
'String','Elevation Range:',...
'Style','text');

function contouring_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

function varargout = contouring(varargin)
% M-File changed by desGUIDE  
% varargin   command line arguments to contouring (see VARARGIN)

%	Copyright (c) 2004-2006 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

	hObject = figure('Tag','figure1','Visible','off','HandleVisibility', 'off');
	contouring_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'east')

	handles.Zmin = [];
	handles.Zmax = [];
	handles.Zinc = [];
	handles.Zstart = [];
	handles.Zsingle = [];
	handles.selected_val = [];

	if ~isempty(varargin)
		handles.hCallingFig = varargin{1};
		handles.head = varargin{2};
		old_cont = varargin{3};
		handles.Xmin = handles.head(1);    handles.Xmax = handles.head(2);
		handles.Ymin = handles.head(3);    handles.Ymax = handles.head(4);
		handles.Zmin = handles.head(5);    handles.Zmax = handles.head(6);
		set(handles.edit_Zmin,'String',num2str(handles.Zmin))
		set(handles.edit_Zmax,'String',num2str(handles.Zmax))
	else
		delete(hObject);    return
	end

	if (~isempty(old_cont))
		list = num2cell(old_cont(:));
		set(handles.listbox_ElevValues,'String',list)
	end

%------------------ Give a Pro look (3D) to the frame boxes --------------------------------
bgcolor = get(0,'DefaultUicontrolBackgroundColor');
framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
h_f = findobj(hObject,'Style','Frame');
for i=1:length(h_f)
    frame_size = get(h_f(i),'Position');
    f_bgc = get(h_f(i),'BackgroundColor');
    usr_d = get(h_f(i),'UserData');
    if abs(f_bgc(1)-bgcolor(1)) > 0.01           % When the frame's background color is not the default's
        frame3D(hObject,frame_size,framecolor,f_bgc,usr_d)
    else
        frame3D(hObject,frame_size,framecolor,'',usr_d)
        delete(h_f(i))
    end
end

% Recopy the text fields on top of previously created frames (uistack is to slow)
h_t = findobj(hObject,'Style','Text');
for i=1:length(h_t)
    usr_d = get(h_t(i),'UserData');
    t_size = get(h_t(i),'Position');   t_str = get(h_t(i),'String');    fw = get(h_t(i),'FontWeight');
    bgc = get (h_t(i),'BackgroundColor');   fgc = get (h_t(i),'ForegroundColor');
    uicontrol('Parent',hObject, 'Style','text', 'Position',t_size,'String',t_str, ...
        'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw,'UserData',usr_d);
end
delete(h_t)
%---------------------------------------------------------------------------------------

	% Add this figure handle to the carra?as list
	plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hCallingFig,'dependentFigs',plugedWin);

	% Update handles structure
	guidata(hObject, handles);

	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

%-------------------------------------------------------------------------------
function pushbutton_GuessIntervals_Callback(hObject, eventdata, handles)
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
		i = i+npoints+1;
		k = k + 1;
	end
	list = num2cell(z_level');
	set(handles.listbox_ElevValues,'String',list)

%-------------------------------------------------------------------------------
function edit_ElevStep_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if ( isempty(xx) || isnan(xx))
		set(hObject,'String','');
		handles.Zinc = [];
	else
		handles.Zinc = xx;
	end
	guidata(hObject, handles);

%-------------------------------------------------------------------------------
function edit_StartElev_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if ( isempty(xx) || isnan(xx))
		set(hObject,'String','');
		handles.Zstart = [];
	else
		handles.Zstart = xx;
	end
	guidata(hObject, handles);

%-------------------------------------------------------------------------------
function edit_SingleElev_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if ( isempty(xx) || isnan(xx))
    set(hObject,'String','');
    handles.Zsingle = [];
else
    handles.Zsingle = xx;
end
guidata(hObject, handles);

%-------------------------------------------------------------------------------
function pushbutton_Add_Callback(hObject, eventdata, handles)
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
function pushbutton_GenerateIntervals_Callback(hObject, eventdata, handles)
if (~isempty(handles.Zstart) && ~isempty(handles.Zinc))
    list = num2cell((handles.Zstart:handles.Zinc:handles.Zmax)');
    set(handles.listbox_ElevValues,'String',list)
else
    errordlg('Generate how? From the empty space?','Chico Clever')
end

%-------------------------------------------------------------------------------
function listbox_ElevValues_Callback(hObject, eventdata, handles)
	handles.selected_val = get(hObject,'Value');
	guidata(handles.figure1, handles);

%-------------------------------------------------------------------------------
function pushbutton_DeleteSelected_Callback(hObject, eventdata, handles)
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
function pushbutton_DeleteAll_Callback(hObject, eventdata, handles)
	set(handles.listbox_ElevValues,'String','')
	handles.selected_val = [];
	guidata(handles.figure1, handles)

	handMir = guidata(handles.hCallingFig);
	delete(findobj(handMir.axes1,'tag','contour'))		% Remove all contours from the Mirone figure
	handMir.which_cont = [];			% Reset so that next time Mirone knows it has to rebuild
	guidata(handMir.figure1, handMir)

%-------------------------------------------------------------------------------
function pushbutton_Apply_Callback(hObject, eventdata, handles)
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

uicontrol('Parent',h1, 'Position',[21 236 80 17],...
'String','Plot Labels',...
'Style','checkbox',...
'TooltipString','Plot master contour lines',...
'Value',1,...
'Tag','check_plotLabels');

uicontrol('Parent',h1, 'Position',[20 206 241 21],...
'Callback',{@contouring_uicallback,h1,'pushbutton_GuessIntervals_Callback'},...
'String','Add Common Charting Intervals',...
'Tag','pushbutton_GuessIntervals');

uicontrol('Parent',h1, 'Units','pixels', 'Position',[20 48 241 141],...
'Style','frame',...
'Tag','frame1');

uicontrol('Parent',h1, 'Position',[30 145 75 15],...
'HorizontalAlignment','left',...
'String','Single Elevation',...
'Style','text');

uicontrol('Parent',h1, 'Position',[30 119 83 15],...
'HorizontalAlignment','left',...
'String','Starting Elevation',...
'Style','text');

uicontrol('Parent',h1, 'Position',[30 89 75 15],...
'HorizontalAlignment','left',...
'String','Elevation Step',...
'Style','text');

uicontrol('Parent',h1, 'Position',[120 85 131 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@contouring_uicallback,h1,'edit_ElevStep_Callback'},...
'Style','edit',...
'Tag','edit_ElevStep');

uicontrol('Parent',h1, 'Position',[120 114 131 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@contouring_uicallback,h1,'edit_StartElev_Callback'},...
'Style','edit',...
'Tag','edit_StartElev');

uicontrol('Parent',h1, 'Position',[120 143 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',{@contouring_uicallback,h1,'edit_SingleElev_Callback'},...
'Style','edit',...
'Tag','edit_SingleElev');

uicontrol('Parent',h1, 'Position',[30 55 221 21],...
'Callback',{@contouring_uicallback,h1,'pushbutton_GenerateIntervals_Callback'},...
'String','Generate Elevation Intervals',...
'Tag','pushbutton_GenerateIntervals');
uicontrol('Parent',h1, 'Position',[200 142 51 21],...
'Callback',{@contouring_uicallback,h1,'pushbutton_Add_Callback'},...
'String','Add',...
'Tag','pushbutton_Add');

uicontrol('Parent',h1, 'Position',[31 177 95 18],...
'String','Custom Elevations',...
'Style','text',...
'Tag','text6');

uicontrol('Parent',h1, 'Position',[290 108 145 180],...
'BackgroundColor',[1 1 1],...
'Callback',{@contouring_uicallback,h1,'listbox_ElevValues_Callback'},...
'Max',2,...
'Style','listbox',...
'Value',1,...
'Tag','listbox_ElevValues');

uicontrol('Parent',h1, 'Position',[290 76 140 21],...
'Callback',{@contouring_uicallback,h1,'pushbutton_DeleteSelected_Callback'},...
'String','Delete Selected Values',...
'Tag','pushbutton_DeleteSelected');

uicontrol('Parent',h1, 'Position',[290 47 140 21],...
'Callback',{@contouring_uicallback,h1,'pushbutton_DeleteAll_Callback'},...
'String','Delete All Values',...
'Tag','pushbutton_DeleteAll');

uicontrol('Parent',h1, 'Position',[308 293 130 15],...
'HorizontalAlignment','left',...
'String','Elevation Values in Range',...
'Style','text');

uicontrol('Parent',h1, 'Position',[363 6 80 21],...
'Callback',{@contouring_uicallback,h1,'pushbutton_Apply_Callback'},...
'String','Apply',...
'Tag','pushbutton_Apply');

uicontrol('Parent',h1, 'Position',[20 294 90 15],...
'HorizontalAlignment','left',...
'String','Elevation Range:',...
'Style','text');

function contouring_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

function varargout = time_stamp(varargin)
% command line arguments to time_stamp

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
 
	hObject = figure('Tag','figure1','Visible','off');
	time_stamp_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	handles.command = cell(15,1);
	handles.command{3} = '-U';
	set(handles.edit_ShowCommand, 'String', [handles.command{3}]);

	% Choose default command line output for time_stamp_export
	handles.output = hObject;
	guidata(hObject, handles);

	set(hObject,'Visible','on');
	% UIWAIT makes time_stamp_export wait for user response (see UIRESUME)
	uiwait(handles.figure1);

	handles = guidata(hObject);
	varargout{1} = handles.output;
	delete(handles.figure1);

function edit_Xpos_CB(hObject, handles)
	xx = get(hObject,'String');
	if isempty(xx)
		set(hObject,'String','0.0');    handles.command{5} = '0.0';
	else
		handles.command{4} = '/'; handles.command{5} = xx;      handles.command{7} = '/';
		handles.command{8} = get(handles.edit_Ypos,'String');       handles.command{10} = '/';
	end
	set(handles.edit_ShowCommand, 'String', [handles.command{3:end}]);
	guidata(hObject, handles);

function edit_Ypos_CB(hObject, handles)
xx = get(hObject,'String');
if isempty(xx)
    set(hObject,'String','-2.5');  handles.command{8} = '-2.5';
else
    handles.command{8} = xx;      handles.command{7} = '/';       handles.command{10} = '/';
    handles.command{5} = get(handles.edit_Xpos,'String');
end
set(handles.edit_ShowCommand, 'String', [handles.command{3:end}]);
guidata(hObject, handles);

% --- Executes on selection change in popup_Unities.
function popup_Unities_CB(hObject, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case 'centimeter'
        if ~isempty(get(handles.edit_Xpos,'String')) && ~isempty(get(handles.edit_Ypos,'String'))
            handles.command{5} = get(handles.edit_Xpos,'String');     handles.command{4} = '/';
            handles.command{6} = 'c';     handles.command{7} = '/';
            handles.command{8} = get(handles.edit_Ypos,'String');
            handles.command{9} = 'c';     handles.command{10} = '/';
        end
    case 'inch'
        if ~isempty(get(handles.edit_Xpos,'String')) && ~isempty(get(handles.edit_Ypos,'String'))
            handles.command{5} = get(handles.edit_Xpos,'String');     handles.command{4} = '/';
            handles.command{6} = 'i';     handles.command{7} = '/';
            handles.command{8} = get(handles.edit_Ypos,'String');
            handles.command{9} = 'i';     handles.command{10} = '/';
        end
    case 'meter'
        if ~isempty(get(handles.edit_Xpos,'String')) && ~isempty(get(handles.edit_Ypos,'String'))
            handles.command{5} = get(handles.edit_Xpos,'String');     handles.command{4} = '/';
            handles.command{6} = 'm';     handles.command{7} = '/';
            handles.command{8} = get(handles.edit_Ypos,'String');
            handles.command{9} = 'm';     handles.command{10} = '/';
        end
    case ''
        if ~isempty(get(handles.edit_Xpos,'String')) && ~isempty(get(handles.edit_Ypos,'String'))
            handles.command{5} = get(handles.edit_Xpos,'String');     handles.command{4} = '/';
            handles.command{6} = '';     handles.command{7} = '/';
            handles.command{8} = get(handles.edit_Ypos,'String');
            handles.command{9} = '';     handles.command{10} = '/';
        end
end
set(handles.edit_ShowCommand, 'String', [handles.command{3:end}]);
guidata(hObject, handles);

function edit_Label_CB(hObject, handles)
xx = get(hObject,'String');
if ~isempty(xx)
    handles.command{11} = ['"' xx '"'];     handles.command{10} = '/';
    set(handles.checkbox_PlotCommand,'Value',0)
else
    handles.command{11} = '';
end
if isempty(handles.command{5}) && isempty(handles.command{8})
    handles.command{10} = '';
end
set(handles.edit_ShowCommand, 'String', [handles.command{3:end}]);
guidata(hObject, handles);

% --- Executes on button press in checkbox_PlotCommand.
function checkbox_PlotCommand_CB(hObject, handles)
if get(hObject,'Value')
    handles.command{11} = 'c';    handles.command{10} = '/';
    set(handles.edit_Label,'String','')
else
    handles.command{11} = '';
end
if isempty(handles.command{5}) && isempty(handles.command{8})
    handles.command{10} = '';
end
set(handles.edit_ShowCommand, 'String', [handles.command{3:end}]);
guidata(hObject, handles);

function edit_ShowCommand_CB(hObject, handles)

% --- Executes on button press in push_Cancel.
function push_Cancel_CB(hObject, handles)
handles.output = '';        % User gave up, return nothing
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in push_OK.
function push_OK_CB(hObject, handles)
handles.output = get(handles.edit_ShowCommand, 'String');
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes on button press in push_Help.
function push_Help_CB(hObject, handles)
message = {'Draw a time stamp on plot. User may specify where the lower left corner of'
    'the stamp should fall on the page relative to lower left corner of plot.'
    'Optionally, append a label (e.g. your sinature), or select to plot the'
    'command string. The GMT parameters UNIX_TIME and UNIX_TIME_POS'
    'can affect the appearance; see the gmtdefaults man page for details.'};
helpdlg(message,'Help on Time Stamp');

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, evt)
	handles = guidata(hObject);
	if (exist('OCTAVE_VERSION','builtin'))		% To know if we are running under Octave
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	else
		do_uiresume = strcmp(get(handles.figure1, 'waitstatus'), 'waiting');
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		handles.output = [];		% User gave up, return nothing
		guidata(hObject, handles);	uiresume(hObject);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, evt)
% Check for "enter" or "escape"
handles = guidata(hObject);
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = '';    % User said no by hitting escape
    guidata(hObject, handles);
    uiresume(handles.figure1);
end   
if isequal(get(hObject,'CurrentKey'),'return')
    handles.output = '';        % User gave up, return nothing
    guidata(hObject, handles);
    uiresume(handles.figure1);
end


% --- Creates and returns a handle to the GUI figure. 
function time_stamp_LayoutFcn(h1)

set(h1,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','time_stamp',...
'NumberTitle','off',...
'Position',[520 638 292 162],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_Xpos_CB'},...
'HorizontalAlignment','left',...
'Position',[51 135 47 21],...
'String','0.0',...
'Style','edit',...
'Tooltip','X Time Stamp position',...
'Tag','edit_Xpos');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_Ypos_CB'},...
'HorizontalAlignment','left',...
'Position',[171 134 47 21],...
'String','-2.5',...
'Style','edit',...
'Tooltip','Y Time Stamp position',...
'Tag','edit_Ypos');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'popup_Unities_CB'},...
'Position',[51 106 82 22],...
'String',{''; 'centimeter'; 'inch'; 'meter'},...
'Style','popupmenu',...
'Tooltip','Measure unities used in Time Stamp positioning',...
'Value',1,...
'Tag','popup_Unities');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_Label_CB'},...
'HorizontalAlignment','left',...
'Position',[51 82 231 21],...
'Style','edit',...
'Tooltip','Plot this label (e.g. your signature)',...
'Tag','edit_Label');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'checkbox_PlotCommand_CB'},...
'Position',[51 60 116 15],...
'String','Plot command string',...
'Style','checkbox',...
'Tooltip','Plot the full GMT command string',...
'Tag','checkbox_PlotCommand');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_ShowCommand_CB'},...
'HorizontalAlignment','left',...
'Position',[51 35 231 21],...
'Style','edit',...
'Tooltip','Show the correspoding GMT command',...
'Tag','edit_ShowCommand');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_Cancel_CB'},...
'Position',[142 6 66 23],...
'String','Cancel',...
'Tag','push_Cancel');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_OK_CB'},...
'Position',[216 6 66 23],...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[11 84 33 18],...
'String','Label',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[10 139 40 15],...
'String','X origin',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[130 137 40 15],...
'String','Y origin',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[11 110 32 15],...
'String','Unities',...
'Style','text');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_Help_CB'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[261 115 21 41],...
'String','?',...
'Tag','push_Help');

function main_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

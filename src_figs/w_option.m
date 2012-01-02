function varargout = w_option(varargin)
% command line arguments to w_option

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
	w_option_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	handles.command = cell(20,1);
	handles.command{4}  = '-W';
	if (isempty(varargin))      % Default Line thickness and color
		handles.command{5}  = '0.5p';
		handles.command{7}  = '0';    handles.command{9}  = '0';    handles.command{11}  = '0';
		handles.command{6}  = '/';    handles.command{8}  = '/';    handles.command{10}  = '/';
		handles.command{12} = '';
	else
		handles.command{5}  = [varargin{1} 'p'];      cor = varargin{2};
		handles.command{7}  = cor{1};   handles.command{9}  = cor{2};   handles.command{11}  = cor{3};
		handles.command{6}  = '/';    handles.command{8}  = '/';    handles.command{10}  = '/';
		handles.command{12} = varargin{3};
		k = strmatch(varargin{1},{'0.5' '1' '2' '3' '5' '7' '10'},'exact');
		if (~isempty(k))                % Check line thickness
			set(handles.popup_LineThickness, 'Value', k);
		end
		if (~isempty(varargin{3}))      % Check line style
			k = strmatch(varargin{3},{'ta' 'to' 't10_2_2_5:5'});
			set(handles.popup_LineType, 'Value', k+1);
		end
		if (~strcmp(cor{1},'0') || ~strcmp(cor{2},'0') || ~strcmp(cor{3},'0'))
			if (strcmp(cor{1},'200') && strcmp(cor{2},'200') && strcmp(cor{3},'200'))
				set(handles.popup_LineColor,'Value',2)
			elseif (strcmp(cor{1},'255') && strcmp(cor{2},'255') && strcmp(cor{3},'255'))
				set(handles.popup_LineColor,'Value',3)
			elseif (strcmp(cor{1},'255') && strcmp(cor{2},'0') && strcmp(cor{3},'0'))
				set(handles.popup_LineColor,'Value',4)
			elseif (strcmp(cor{1},'0') && strcmp(cor{2},'255') && strcmp(cor{3},'0'))
				set(handles.popup_LineColor,'Value',5)
			elseif (strcmp(cor{1},'0') && strcmp(cor{2},'0') && strcmp(cor{3},'255'))
				set(handles.popup_LineColor,'Value',6)
			elseif (strcmp(cor{1},'0') && strcmp(cor{2},'255') && strcmp(cor{3},'255'))
				set(handles.popup_LineColor,'Value',7)
			elseif (strcmp(cor{1},'255') && strcmp(cor{2},'255') && strcmp(cor{3},'0'))
				set(handles.popup_LineColor,'Value',8)
			elseif (strcmp(cor{1},'255') && strcmp(cor{2},'0') && strcmp(cor{3},'255'))
				set(handles.popup_LineColor,'Value',9)
			end
		end
	end

	set(handles.edit_LineThickness_pt, 'String', handles.command{5}(1:end-1));
	set(handles.edit_LineColor_R, 'String', [handles.command{7}]);
	set(handles.edit_LineColor_G, 'String', handles.command{9});
	set(handles.edit_LineColor_B, 'String', handles.command{11});
	set(handles.edit_ShowCommand, 'String', [handles.command{4:end}]);
	pushbutton_Example_CB(hObject, handles)

	a(:,:,1) = color_wheel(10,0);     % 21x21
	a(:,:,2) = color_wheel(10,1);
	a(:,:,3) = color_wheel(10,2);
	set(handles.pushbutton_CustomColor,'CData',a)

	% Choose default command line output for w_option_export
	handles.output = hObject;
	guidata(hObject, handles);

	set(hObject,'Visible','on');
	% UIWAIT makes w_option_export wait for user response (see UIRESUME)
	uiwait(handles.figure1);

	handles = guidata(hObject);
	varargout{1} = handles.output;
	delete(handles.figure1);

%----------------------------------------------------------------
function C = color_wheel(n,i)     % It will eventualy become that
	r = 0:n;  r = [r (n-1):-1:0]'/n;
	theta = (pi+i*2*pi/3)*(-n:n)/n;
	X = r*cos(theta);
	C = (X + 1) / 2;

%-----------------------------------------------------------------
function popup_LineType_CB(hObject, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case 'Solid line'
        handles.command{12} = '';
    case 'Dashed line'
        handles.command{12} = 'ta';
    case 'Doted line'
        handles.command{12} = 'to';
    case 'Dash-Dot line'
        handles.command{12} = 't10_2_2_5:5';
end
set(handles.edit_ShowCommand, 'String', [handles.command{4:end}]);
pushbutton_Example_CB(hObject, handles)
guidata(hObject, handles);

%----------------------------------------------------------------------------------------------
function popup_LineThickness_CB(hObject, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case '0.5  pt'
        set(handles.edit_LineThickness_pt, 'String', '0.5');     handles.command{5} = '0.5p';
    case '1  pt'
        set(handles.edit_LineThickness_pt, 'String', '1');     handles.command{5} = '1p';
    case '2  pt'
        set(handles.edit_LineThickness_pt, 'String', '2');     handles.command{5} = '2p';
    case '3  pt'
        set(handles.edit_LineThickness_pt, 'String', '3');     handles.command{5} = '3p';
    case '5  pt'
        set(handles.edit_LineThickness_pt, 'String', '5');     handles.command{5} = '5p';
    case '7  pt'
        set(handles.edit_LineThickness_pt, 'String', '7');     handles.command{5} = '7p';
    case '10 pt'
        set(handles.edit_LineThickness_pt, 'String', '10');     handles.command{5} = '10p';
end
set(handles.edit_ShowCommand, 'String', [handles.command{4:end}]);
pushbutton_Example_CB(hObject, handles)
guidata(hObject, handles);

%----------------------------------------------------------------------------------------------
function popup_LineColor_CB(hObject, handles)
% Hints: contents = get(hObject,'String') returns popup_LineColor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_LineColor
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case 'Black'
        set(handles.edit_LineColor_R, 'String', '0');     handles.command{7} = '0';
        set(handles.edit_LineColor_G, 'String', '0');     handles.command{9} = '0';
        set(handles.edit_LineColor_B, 'String', '0');     handles.command{11} = '0';
    case 'Light gray'
        set(handles.edit_LineColor_R, 'String', '200');   handles.command{7} = '200';
        set(handles.edit_LineColor_G, 'String', '200');   handles.command{9} = '200';
        set(handles.edit_LineColor_B, 'String', '200');   handles.command{11} = '200';
    case 'White'
        set(handles.edit_LineColor_R, 'String', '255');   handles.command{7} = '255';
        set(handles.edit_LineColor_G, 'String', '255');   handles.command{9} = '255';
        set(handles.edit_LineColor_B, 'String', '255');   handles.command{11} = '255';
    case 'Red'
        set(handles.edit_LineColor_R, 'String', '255');   handles.command{7} = '255';
        set(handles.edit_LineColor_G, 'String', '0');     handles.command{9} = '0';
        set(handles.edit_LineColor_B, 'String', '0');     handles.command{11} = '0';
    case 'Green'
        set(handles.edit_LineColor_R, 'String', '0');     handles.command{7} = '0';
        set(handles.edit_LineColor_G, 'String', '255');   handles.command{9} = '255';
        set(handles.edit_LineColor_B, 'String', '0');     handles.command{11} = '0';
    case 'Blue'
        set(handles.edit_LineColor_R, 'String', '0');     handles.command{7} = '0';
        set(handles.edit_LineColor_G, 'String', '0');     handles.command{9} = '0';
        set(handles.edit_LineColor_B, 'String', '255');   handles.command{11} = '255';
    case 'Cyan'
        set(handles.edit_LineColor_R, 'String', '0');     handles.command{7} = '0';
        set(handles.edit_LineColor_G, 'String', '255');   handles.command{9} = '255';
        set(handles.edit_LineColor_B, 'String', '255');   handles.command{11} = '255';
    case 'Yellow'
        set(handles.edit_LineColor_R, 'String', '255');   handles.command{7} = '255';
        set(handles.edit_LineColor_G, 'String', '255');   handles.command{9} = '255';
        set(handles.edit_LineColor_B, 'String', '0');     handles.command{11} = '0';
    case 'Magenta'
        set(handles.edit_LineColor_R, 'String', '255');   handles.command{7} = '255';
        set(handles.edit_LineColor_G, 'String', '0');     handles.command{9} = '0';
        set(handles.edit_LineColor_B, 'String', '255');   handles.command{11} = '255';
end        
set(handles.edit_ShowCommand, 'String', [handles.command{4:end}]);
pushbutton_Example_CB(hObject, handles)
guidata(hObject, handles);

%----------------------------------------------------------------------------------------------
function pushbutton_CustomColor_CB(hObject, handles)
c = uisetcolor;
if length(c) > 1            % That is, if a color was selected
    c(1) = round(c(1)*255);     c(2) = round(c(2)*255);     c(3) = round(c(3)*255);
    handles.command{7} = num2str(c(1)); handles.command{9} = num2str(c(2)); handles.command{11} = num2str(c(3));
    set(handles.edit_LineColor_R, 'String', num2str(c(1)));
    set(handles.edit_LineColor_G, 'String', num2str(c(2)));
    set(handles.edit_LineColor_B, 'String', num2str(c(3)));
    set(handles.edit_ShowCommand, 'String', [handles.command{4:end}]);
    pushbutton_Example_CB(hObject, handles)
    guidata(hObject, handles);
end

%----------------------------------------------------------------------------------------------
function edit_LineThickness_pt_CB(hObject, handles)
handles.command{5} = [get(hObject,'String') 'p'];
set(handles.edit_ShowCommand, 'String', [handles.command{4:end}]);
pushbutton_Example_CB(hObject, handles)
guidata(hObject, handles);

%----------------------------------------------------------------------------------------------
function edit_LineColor_R_CB(hObject, handles)
	xx = get(hObject,'String');
	if (str2double(xx) > 255 || str2double(xx) < 0)
		set(handles.edit_LineColor_R, 'String', '0');
		xx = '0';
	end
	handles.command{7} = num2str(fix(str2double(xx)));
	set(handles.edit_ShowCommand, 'String', [handles.command{4:end}]);
	pushbutton_Example_CB(hObject, handles)
	guidata(hObject, handles);

%----------------------------------------------------------------------------------------------
function edit_LineColor_G_CB(hObject, handles)
	xx = get(hObject,'String');
	if (str2double(xx) > 255 || str2double(xx) < 0)
		set(handles.edit_LineColor_G, 'String', '0');
		xx = '0';
	end
	handles.command{9} = num2str(fix(str2double(xx)));
	set(handles.edit_ShowCommand, 'String', [handles.command{4:end}]);
	pushbutton_Example_CB(hObject, handles)
	guidata(hObject, handles);

%----------------------------------------------------------------------------------------------
function edit_LineColor_B_CB(hObject, handles)
	xx = get(hObject,'String');
	if (str2double(xx) > 255 || str2double(xx) < 0)
		set(handles.edit_LineColor_G, 'String', '0');
		xx = '0';
	end
	handles.command{11} = num2str(fix(str2double(xx)));
	set(handles.edit_ShowCommand, 'String', [handles.command{4:end}]);
	pushbutton_Example_CB(hObject, handles)
	guidata(hObject, handles);

%----------------------------------------------------------------------------------------------
function pushbutton_Example_CB(hObject, handles)
	x = linspace(0,1,10);   y = ones(1,10);
	if ischar(handles.command{5}(end))    thick = str2double(handles.command{5}(1:end-1));
	else    thick = str2double(handles.command{5});     end
	r = str2double(handles.command{7}) / 255;
	g = str2double(handles.command{9}) / 255;
	b = str2double(handles.command{11}) / 255;
	if strcmp(handles.command{12},'ta')      % dashed
		ltype = '--';
	elseif strcmp(handles.command{12},'to')  % doted
		ltype = ':';
	elseif strcmp(handles.command{12},'')  % solid
		ltype = '-';
	else                                    % dash-dot
		ltype = '-.';
	end
	plot(x,y,ltype,'LineWidth',thick,'Color',[r g b])
	set(gca,'Visible','off')

%----------------------------------------------------------------------------------------------
function pushbutton_OK_CB(hObject, handles)
	handles.output = get(handles.edit_ShowCommand, 'String');
	guidata(hObject,handles)
	uiresume(handles.figure1);

%----------------------------------------------------------------------------------------------
function pushbutton_Cancel_CB(hObject, handles)
	handles.output = '';        % User gave up, return nothing
	guidata(hObject, handles);
	uiresume(handles.figure1);

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
		guidata(handles.figure1, handles);	uiresume(handles.figure1);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, evt)
% Check for "escape"
handles = guidata(hObject);
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = '';    % User said no by hitting escape
    guidata(hObject, handles);
    uiresume(handles.figure1);
end   

% --- Creates and returns a handle to the GUI figure. 
function w_option_LayoutFcn(h1)

set(h1,...
'PaperUnits','centimeters',...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','w_option',...
'NumberTitle','off',...
'Position',[265.768111202607 265.768111202607 392 204],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'popup_LineThickness_CB'},...
'Position',[83 172 91 22],...
'String',{'0.5  pt'; '1  pt'; '2  pt'; '3  pt'; '5  pt'; '7  pt'; '10 pt';},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_LineThickness');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'popup_LineColor_CB'},...
'Position',[83 145 91 22],...
'String',{  'Black'; 'Light gray'; 'White'; 'Red'; 'Green'; 'Blue'; 'Cyan'; 'Yellow'; 'Magenta' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_LineColor');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'popup_LineType_CB'},...
'Position',[83 118 91 22],...
'String',{  'Solid line'; 'Dashed line'; 'Doted line'; 'Dash-Dot line' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_LineType');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_LineThickness_pt_CB'},...
'CData',[],...
'HorizontalAlignment','left',...
'Position',[193 173 35 21],...
'Style','edit',...
'TooltipString','Select a different thickness (in points)',...
'Tag','edit_LineThickness_pt');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0 0],...
'Call',{@main_uiCB,h1,'edit_LineColor_R_CB'},...
'HorizontalAlignment','left',...
'Position',[193 148 31 21],...
'Style','edit',...
'TooltipString','Red component (0-255)',...
'Tag','edit_LineColor_R');

uicontrol('Parent',h1,...
'BackgroundColor',[0 1 0],...
'Call',{@main_uiCB,h1,'edit_LineColor_G_CB'},...
'HorizontalAlignment','left',...
'Position',[224 148 31 21],...
'Style','edit',...
'TooltipString','Green component (0-255)',...
'Tag','edit_LineColor_G');

uicontrol('Parent',h1,...
'BackgroundColor',[0 0 1],...
'Call',{@main_uiCB,h1,'edit_LineColor_B_CB'},...
'HorizontalAlignment','left',...
'Position',[255 148 31 21],...
'Style','edit',...
'TooltipString','Blue component (0-255)',...
'Tag','edit_LineColor_B');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'pushbutton_CustomColor_CB'},...
'Position',[292 148 22 21],...
'TooltipString','Interactive color selection',...
'Tag','pushbutton_CustomColor');

uicontrol('Parent',h1,...
'Enable','off',...
'Position',[193 119 71 21],...
'String','Custom type',...
'Tag','pushbutton_CustomLineType');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[13 26 291 21],...
'Style','edit',...
'Tag','edit_ShowCommand');

axes('Parent',h1,...
'Units','pixels',...
'Position',[13 75 295 26],...
'Tag','axes_example');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[11 174 69 19],...
'String','Line thickness',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[30 150 51 15],...
'String','Line color',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[30 122 51 15],...
'String','Line type',...
'Style','text');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'pushbutton_Example_CB'},...
'Position',[323 77 46 23],...
'String','Show',...
'Tag','pushbutton_Example',...
'Visible','off');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'pushbutton_OK_CB'},...
'Position',[318 10 66 23],...
'String','OK',...
'Tag','pushbutton_OK');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'pushbutton_Cancel_CB'},...
'Position',[318 40 66 23],...
'String','Cancel',...
'Tag','pushbutton_Cancel');

function main_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

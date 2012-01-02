function varargout = pscoast_options_Mir(varargin)
%	command line arguments to pscoast_options_Mir

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
 
	hObject = figure('Vis','off');
	pscoast_options_Mir_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right');

	handles.gray = [.706,.706,.706];
	handles.red = [1 0 0];
	handles.green = [0 1 0];
	handles.blue = [0 0 1];
	handles.ColorDes = [.764,.603,.603];    % Color for desabled options

	handles.command = cell(50,1);
	handles.repeat_I = 0;   handles.repeat_N = 0;

	if (~isempty(varargin))
		handMir = varargin{1};				handles.psc_res = varargin{2};
		handles.psc_opt_W = varargin{3};	handles.psc_type_p = varargin{4};
		handles.psc_type_r = varargin{5};
	end

	handles.path_data = handMir.path_data;

	% See what coastlines resolution are allowed
	str_res = {'crude'; 'low'; 'intermediate'};
	if (strcmp(get(handMir.CoastLineHigh,'Enable'),'on'))
		str_res{end+1} = 'high';
	end
	if (strcmp(get(handMir.CoastLineFull,'Enable'),'on'))
		str_res{end+1} = 'full';
	end
	set(handles.popup_Resolution,'String',str_res)

	% ------------ See if a db resolution was transmited
	if (~isempty(handles.psc_res))
		switch handles.psc_res(3)
			case 'c',   set(handles.popup_Resolution,'Value',1);    handles.command{1} = ' -Dc';
			case 'l',   set(handles.popup_Resolution,'Value',2);    handles.command{1} = ' -Dl';
			case 'i',   set(handles.popup_Resolution,'Value',3);    handles.command{1} = ' -Di';
			case 'h',   set(handles.popup_Resolution,'Value',4);    handles.command{1} = ' -Dh';
			case 'f',   set(handles.popup_Resolution,'Value',5);    handles.command{1} = ' -Df';
		end
	else
		set(handles.popup_Resolution,'Value',2);
	end

	% ------------------------------------------------
	% GMT4 -W color syntax changed from the time I first wrote this function.
	% It was than (e.g.) -W1p/0/0/0 which now writes -W1p,0/0/0. To minimize head-hakes
	% just replace the ',' by an '/' and use old algo.
	handles.coast_ls = '';      handles.coast_lt = [];      handles.coast_cor = cell(3,1);
	if (~isempty(handles.psc_opt_W))
		psc_opt_W = handles.psc_opt_W;			% Need to make a copy
		handles.command{7} = [' ' psc_opt_W];
		psc_opt_W = strrep(psc_opt_W,',','/');
		k = strfind(psc_opt_W,'/');
		handles.coast_lt = psc_opt_W(3:k(1)-2);
		handles.coast_cor = {psc_opt_W(k(1)+1:k(2)-1); psc_opt_W(k(2)+1:k(3)-1); psc_opt_W(k(3)+1:end)};
		k = strfind(handles.coast_cor{3},'t');
		if (~isempty(k))             % We have a line style other than solid
			handles.coast_ls = handles.coast_cor{3}(k:end);
			handles.coast_cor{3} = handles.coast_cor{3}(1:k-1);
		end
	end 
	handles.political_ls = '';  handles.political_lt = [];  handles.political_cor = cell(3,1);
	if (~isempty(handles.psc_type_p))
		psc_type_p = handles.psc_type_p;
		handles.command{27} = [' ' psc_type_p];
		psc_type_p = strrep(psc_type_p,',','/');
		k = strfind(handles.psc_type_p,'/');
		handles.political_lt = psc_type_p(k(1)+1:k(2)-2);
		handles.political_cor = {psc_type_p(k(2)+1:k(3)-1); psc_type_p(k(3)+1:k(4)-1); psc_type_p(k(4)+1:end)};
		k = strfind(handles.political_cor{3},'t');
		if (~isempty(k))             % We have a line style other than solid
			handles.political_ls = handles.political_cor{3}(k:end);
			handles.political_cor{3} = handles.political_cor{3}(1:k-1);
		end
		switch handles.psc_type_p(3)
			case '1',   set(handles.popup_PoliticalBound,'Value',1);
			case '2',   set(handles.popup_PoliticalBound,'Value',2);
			case '3',   set(handles.popup_PoliticalBound,'Value',3);
			case 'a',   set(handles.popup_PoliticalBound,'Value',4);
		end
	end
	handles.rivers_ls = '';  handles.rivers_lt = [];  handles.rivers_cor = cell(3,1);
	if (~isempty(handles.psc_type_r))
		psc_type_r = handles.psc_type_r;
		handles.command{25} = [' ' psc_type_r];
		psc_type_r = strrep(psc_type_r,',','/');
		k = strfind(psc_type_r,'/');
		handles.rivers_lt = handles.psc_type_r(k(1)+1:k(2)-2);
		handles.rivers_cor = {psc_type_r(k(2)+1:k(3)-1); psc_type_r(k(3)+1:k(4)-1); psc_type_r(k(4)+1:end)};
		k = strfind(handles.rivers_cor{3},'t');
		if (~isempty(k))             % We have a line style other than solid
			handles.rivers_ls = handles.rivers_cor{3}(k:end);
			handles.rivers_cor{3} = handles.rivers_cor{3}(1:k-1);
		end
		switch psc_type_r(3:4)
			case '1/',   set(handles.popup_Rivers,'Value',1);
			case '2/',   set(handles.popup_Rivers,'Value',2);
			case '3/',   set(handles.popup_Rivers,'Value',3);
			case '4/',   set(handles.popup_Rivers,'Value',4);
			case '5/',   set(handles.popup_Rivers,'Value',5);
			case '6/',   set(handles.popup_Rivers,'Value',6);
			case '7/',   set(handles.popup_Rivers,'Value',7);
			case '8/',   set(handles.popup_Rivers,'Value',8);
			case '9/',   set(handles.popup_Rivers,'Value',9);
			case '10',   set(handles.popup_Rivers,'Value',10);
			case 'a/',   set(handles.popup_Rivers,'Value',11);
			case 'r/',   set(handles.popup_Rivers,'Value',12);
			case 'i/',   set(handles.popup_Rivers,'Value',13);
			case 'c/',   set(handles.popup_Rivers,'Value',14);
		end
	end
	% --------------------------

	set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);

	% Load background image
	fundo = imread([handles.path_data 'caravela.jpg']);
	image(fundo)

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, findobj(hObject,'Style','Text'))
	%------------- END Pro look (3D) ------------------------------

	% Choose default command line output for pscoast_options_Mir_export
	handles.output = hObject;
	guidata(hObject, handles);

	set(hObject,'Vis','on');
	% UIWAIT makes pscoast_options_Mir_export wait for user response (see UIRESUME)
	uiwait(handles.figure1);

	handles = guidata(hObject);
	varargout{1} = handles.output;
	delete(handles.figure1);

% -----------------------------------------------------------------------------------
function popup_Resolution_CB(hObject, handles)
val = get(hObject,'Value');
switch val;
    case 1,        handles.command{1} = ' -Dc';
    case 2,        handles.command{1} = ' -Dl';
    case 3,        handles.command{1} = ' -Di';
    case 4,        handles.command{1} = ' -Dh';
    case 5,        handles.command{1} = ' -Df';
end
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_HelpOption_D_CB(hObject, handles)
message = {'Selects the resolution of the data set to use: full, high,'
    'intermediate, low, and crude.  The  resolution drops off'
    'by 80% between data sets.'};
helpdlg(message,'Help on Coast lines resolution');

% -----------------------------------------------------------------------------------
function push_Option_W_CB(hObject, handles)
tmp = w_option(handles.coast_lt,handles.coast_cor,handles.coast_ls);
if (isempty(tmp))       return;     end
handles.command{7} = [' ' tmp];     k = strfind(tmp,'/');
handles.coast_lt = tmp(3:k(1)-2);
handles.coast_cor = {tmp(k(1)+1:k(2)-1); tmp(k(2)+1:k(3)-1); tmp(k(3)+1:end)};
k = strfind(handles.coast_cor{3},'t');
if (~isempty(k))             % We have a line style other than solid
    handles.coast_ls = handles.coast_cor{3}(k:end);
    handles.coast_cor{3} = handles.coast_cor{3}(1:k-1);
end
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_HelpOption_W_CB(hObject, handles)
message = {'Draw coastlines. There is not much else to say, the selecting window '
    'should be self explanatory. Press Clear to remove this option'};
helpdlg(message,'Help on Draw coastlines');

% -----------------------------------------------------------------------------------
function push_Option_W_clear_CB(hObject, handles)
handles.command{7} = '';
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function edit_Option_A_MinArea_CB(hObject, handles)
xx = get(hObject,'String');
if isnan(str2double(xx)) 
    errordlg('Not a valid number','Error')
    set(hObject,'String','')
    return
end
handles.command{9} = [' -A' xx];
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function listbox_Option_A_HierarchicalLevel_CB(hObject, handles)
% To be programed when I understand what it means

% -----------------------------------------------------------------------------------
function push_Option_A_clear_CB(hObject, handles)
for i = 9:14     handles.command{i} = '';    end
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_HelpOption_A_CB(hObject, handles)
message = {'Features with an area smaller than "Min Area" in km^2 or of hierarchical'
    'level that is lower than min_level or higher than max_level will not be'
    'plotted [Default is 0/0/4 (all features)]. See DATABASE INFORMATION'
    'in pscoast man page for more details.'};
helpdlg(message,'Help on Min Area');

% -----------------------------------------------------------------------------------
function push_Option_G_CB(hObject, handles)
xx = paint_option([handles.path_data 'GMT_patterns.jpg']);
if ~isempty(xx)     handles.command{16} = ' -G';   handles.command{17} = xx; end
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_Option_G_clear_CB(hObject, handles)
for i = 16:17     handles.command{i} = '';    end
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_HelpOption_G_CB(hObject, handles)
message = {'Select  painting or clipping of "dry" areas. A lengthy help is'
    'available in the option window.'};
helpdlg(message,'Help on Land color fill');

% -----------------------------------------------------------------------------------
function push_Option_S_CB(hObject, handles)
xx = paint_option([handles.path_data 'GMT_patterns.jpg']);
if ~isempty(xx)     handles.command{19} = ' -S';   handles.command{20} = xx; end
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_Option_S_clear_CB(hObject, handles)
for i = 19:20     handles.command{i} = '';    end
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_HelpOption_S_CB(hObject, handles)
message = {'Select  painting or clipping of "wet" areas. A lengthy help is'
    'available in the option window.'};
helpdlg(message,'Help on Water color fill');

% -----------------------------------------------------------------------------------
function push_Option_C_CB(hObject, handles)
xx = paint_option([handles.path_data 'GMT_patterns.jpg']);
if ~isempty(xx)     handles.command{22} = ' -C';   handles.command{23} = xx; end
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_Option_C_clear_CB(hObject, handles)
for i = 22:23     handles.command{i} = '';    end
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_HelpOption_C_CB(hObject, handles)
	message = {'Set  the  shade, color, or pattern for lakes [Default is the fill chosen'
			'for "wet" areas]. A lengthy help is available in the option window.'};
	helpdlg(message,'Help on Lake color fill');

% -----------------------------------------------------------------------------------
function popup_Rivers_CB(hObject, handles)
	val = get(hObject,'Value');
	switch val;
		case 1,     handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -I1');
		case 2,     handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -I2');
		case 3,     handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -I3');
		case 4,     handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -I4');
		case 5,     handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -I5');
		case 6,     handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -I6');
		case 7,     handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -I7');
		case 8,     handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -I8');
		case 9,     handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -I9');
		case 10,    handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -I10');
		case 11,    handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -Ia');
		case 12,    handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -Ir');
		case 13,    handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -Ii');
		case 14,    handles.command{25} = set_PolitRiver(handles, 25, handles.repeat_I, ' -Ic');
	end
	set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
	guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_RiversPen_CB(hObject, handles)
xx = w_option(handles.rivers_lt,handles.rivers_cor,handles.rivers_ls);
if (isempty(xx))       return;     end
if isempty(handles.command{25})
    handles.command{25} = [' -I1/' xx(3:end)];
else
    if (handles.repeat_I)
        k = strfind(handles.command{25},' ');
        handles.command{25} = [handles.command{25}(1:k(end)+3) '/' xx(3:end)];    
    else
        if (~strcmp(handles.command{25}(1:5),' -I10'))
            handles.command{25} = [handles.command{25}(1:5) xx(3:end)];
        else
            handles.command{25} = [handles.command{25}(1:6) xx(3:end)];
        end
    end
end
k = strfind(xx,'/');
handles.rivers_lt = xx(3:k(1)-2);
handles.rivers_cor = {xx(k(1)+1:k(2)-1); xx(k(2)+1:k(3)-1); xx(k(3)+1:end)};
k = strfind(handles.rivers_cor{3},'t');
if (~isempty(k))             % We have a line style other than solid
    handles.rivers_ls = handles.rivers_cor{3}(k:end);
    handles.rivers_cor{3} = handles.rivers_cor{3}(1:k-1);
end
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function checkbox_Option_I_Repeat_CB(hObject, handles)
	if get(hObject,'Value')		handles.repeat_I = 1;
	else						handles.repeat_I = 0;
	end
	guidata(hObject, handles)
    
% -----------------------------------------------------------------------------------
function push_ClearRivers_CB(hObject, handles)
	handles.command{25} = '';
	set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
	guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function popup_PoliticalBound_CB(hObject, handles)
	val = get(hObject,'Value');
	switch val;
		case 1,		handles.command{27} = set_PolitRiver(handles, 27, handles.repeat_N, ' -N1');
		case 2,		handles.command{27} = set_PolitRiver(handles, 27, handles.repeat_N, ' -N2');
		case 3,		handles.command{27} = set_PolitRiver(handles, 27, handles.repeat_N, ' -N3');
		case 4,		handles.command{27} = set_PolitRiver(handles, 27, handles.repeat_N, ' -Na');
	end
	set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
	guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function out = set_PolitRiver(handles,n,repeat,str)
if isempty(handles.command{n})
    out = str;
else
    if (repeat)     out = [handles.command{n} str];
    else            out = [str handles.command{n}(5:end)];    end
end

% -----------------------------------------------------------------------------------
function push_BoundariesPen_CB(hObject, handles)
xx = w_option(handles.political_lt,handles.political_cor,handles.political_ls);
if (isempty(xx))       return;     end
if isempty(handles.command{27})
    handles.command{27} = [' -N1/' xx(3:end)];
else
    if (handles.repeat_N)
        k = strfind(handles.command{27},' ');
        handles.command{27} = [handles.command{27}(1:k(end)+3) '/' xx(3:end)];    
    else
        handles.command{27} = [handles.command{27}(1:5) xx(3:end)];    
    end
end
k = strfind(xx,'/');
handles.political_lt = xx(3:k(1)-2);
handles.political_cor = {xx(k(1)+1:k(2)-1); xx(k(2)+1:k(3)-1); xx(k(3)+1:end)};
k = strfind(handles.political_cor{3},'t');
if (~isempty(k))             % We have a line style other than solid
    handles.political_ls = handles.political_cor{3}(k:end);
    handles.political_cor{3} = handles.political_cor{3}(1:k-1);
end
set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function checkbox_Option_N_Repeat_CB(hObject, handles)
	if get(hObject,'Value')		handles.repeat_N = 1;
	else						handles.repeat_N = 0;
	end
	guidata(hObject, handles)

% -----------------------------------------------------------------------------------
function push_ClearPoliticalBound_CB(hObject, handles)
	handles.command{27} = '';
	set(handles.edit_ShowCommand, 'String', [handles.command{1:end}]);
	guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_Cancel_CB(hObject, handles)
	handles.output = '';        % User gave up, return nothing
	guidata(hObject, handles);
	uiresume(handles.figure1);

% -----------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
	handles.output = get(handles.edit_ShowCommand, 'String');
	guidata(hObject,handles);
	uiresume(handles.figure1);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, evt)
	handles = guidata(hObject);
	if (exist('OCTAVE_VERSION','builtin'))		% To know if we are running under Octave
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	else
		do_uiresume = strcmp(get(hObject, 'waitstatus'), 'waiting');
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		handles.output = [];		% User gave up, return nothing
		guidata(hObject, handles);	uiresume(hObject);
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
function pscoast_options_Mir_LayoutFcn(h1)

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','pscoast options',...
'NumberTitle','off',...
'Position',[266 184 421 360],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[220 247 191 40],'Style','frame');
uicontrol('Parent',h1,'Position',[10 309 131 42],'Style','frame');
uicontrol('Parent',h1,'Position',[10 248 131 40],'Style','frame');
uicontrol('Parent',h1,'Position',[10 76 401 70],'Style','frame');
uicontrol('Parent',h1,'Position',[10 174 131 40],'Style','frame');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[20 319 85 22],...
'String',{'crude'; 'low'; 'intermediate'; 'high'; 'full' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_Resolution');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[110 319 21 22],...
'String','?',...
'Tag','push_HelpOption_D');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[20 255 45 23],...
'String','How-to',...
'Tag','push_Option_W');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[70 255 33 23],...
'String','Clear',...
'Tag','push_Option_W_clear');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[108 255 21 23],...
'String','?',...
'Tag','push_HelpOption_W');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@pscoast_options_Mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[230 254 47 21],...
'Style','edit',...
'Tag','edit_Option_A_MinArea');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@pscoast_options_Mir_uiCB,...
'Enable','off',...
'Position',[280 253 51 21],...
'String',{  '0'; '1'; '2'; '3'; '4' },...
'Style','listbox',...
'Value',1,...
'Tag','listbox_Option_A_HierarchicalLevel');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[336 253 33 22],...
'String','Clear',...
'Tag','push_Option_A_clear');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[374 253 21 22],...
'String','?',...
'Tag','push_HelpOption_A');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[20 182 45 22],...
'String','How-to',...
'Tag','push_Option_G');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[21 208 105 15],...
'String','Paint or clip land (-G)',...
'Style','text');

uicontrol('Parent',h1,'Position',[145 174 131 40],'Style','frame');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[70 182 33 22],...
'String','Clear',...
'Tag','push_Option_G_clear');

uicontrol(...
'Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[108 182 21 22],...
'String','?',...
'Tag','push_HelpOption_G');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[155 182 45 22],...
'String','How-to',...
'Tag','push_Option_S');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[151 207 120 15],...
'String','Paint or clip oceans (-S)',...
'Style','text');

uicontrol('Parent',h1,'Position',[280 174 131 40],'Style','frame');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[205 182 33 22],...
'String','Clear',...
'Tag','push_Option_S_clear');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[243 182 21 22],...
'String','?',...
'Tag','push_HelpOption_S');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[290 182 45 22],...
'String','How-to',...
'Tag','push_Option_C');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[288 207 110 17],...
'String','Paint or clip Lakes (-C)',...
'Style','text');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[340 182 33 22],...
'String','Clear',...
'Tag','push_Option_C_clear');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[378 183 21 22],...
'String','?',...
'Tag','push_HelpOption_C');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[32 115 161 22],...
'String',{'Permanent major rivers'; 'Additional major rivers'; 'Additional rivers'; 'Minor rivers'; 'Intermittent rivers - major'; 'Intermittent rivers - additional'; 'Intermittent rivers - minor'; 'Major canals'; 'Minor canals'; 'Irrigation canals'; 'All rivers and canals'; 'All permanent rivers'; 'All intermittent rivers'; 'All canals' },...
'Style','popupmenu',...
'Tooltip','Draw rivers',...
'Value',1,...
'Tag','popup_Rivers');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[202 114 61 23],...
'String','Pen Attrib',...
'Tooltip','Specify line thickness, color,etc ..',...
'Tag','push_RiversPen');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[270 115 55 22],...
'String','Repeat',...
'Style','checkbox',...
'Tooltip','Repeat this option as often as necessary',...
'Tag','checkbox_Option_I_Repeat');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[342 114 51 23],...
'String','Clear',...
'Tooltip','Remove this option',...
'Tag','push_ClearRivers');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[32 85 161 22],...
'String',{  'National boundaries'; 'State boundaries (NAm)'; 'Marine boundaries'; 'All boundaries' },...
'Style','popupmenu',...
'Tooltip','Draw political boundaries',...
'Value',1,...
'Tag','popup_PoliticalBound');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[202 84 61 23],...
'String','Pen Attrib',...
'Tooltip','Specify line thickness, color,etc ..',...
'Tag','push_BoundariesPen');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[272 84 55 22],...
'String','Repeat',...
'Style','checkbox',...
'Tooltip','Repeat this option as often as necessary',...
'Tag','checkbox_Option_N_Repeat');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[342 84 51 23],...
'String','Clear',...
'Tooltip','Remove this option',...
'Tag','push_ClearPoliticalBound');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[24 279 105 15],...
'String','Draw Coastline (-W)',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[10 6 401 21],...
'Style','edit',...
'Tooltip','Display of GMT command',...
'Tag','edit_ShowCommand');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[330 38 81 29],...
'String','Cancel',...
'Tag','push_Cancel');

uicontrol('Parent',h1,...
'Call',@pscoast_options_Mir_uiCB,...
'Position',[240 38 81 29],...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[80 139 220 15],...
'String','Plot rivers and/or national boundaries (-I, -N)',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[18 343 116 15],...
'String','Coastline resolution',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[270 277 75 15],...
'String','Min area (-A)',...
'Style','text');

axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Position',[0 29 421 332],...
'Tag','axes1');

function pscoast_options_Mir_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

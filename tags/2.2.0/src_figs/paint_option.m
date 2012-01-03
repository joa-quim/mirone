function varargout = paint_option(varargin)
% command line arguments to paint_option

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
	paint_option_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right');

	handles.patternFile = varargin{1};

	handles.gray = [.706,.706,.706];
	handles.red = [1 0 0];  handles.green = [0 1 0];  handles.blue = [0 0 1];
	handles.ColorDes = [.764,.603,.603];    % Color for desabled options
	set(handles.edit_OptionG_gray, 'Backgroundcolor',handles.gray)
	set(handles.edit_OptionG_R, 'Backgroundcolor',handles.red)
	set(handles.edit_OptionG_G, 'Backgroundcolor',handles.green)
	set(handles.edit_OptionG_B, 'Backgroundcolor',handles.blue)
	set(handles.edit_OptionG_Foreground_R, 'Backgroundcolor',handles.red)
	set(handles.edit_OptionG_Foreground_G, 'Backgroundcolor',handles.green)
	set(handles.edit_OptionG_Foreground_B, 'Backgroundcolor',handles.blue)
	set(handles.edit_OptionG_Background_R, 'Backgroundcolor',handles.red)
	set(handles.edit_OptionG_Background_G, 'Backgroundcolor',handles.green)
	set(handles.edit_OptionG_Background_B, 'Backgroundcolor',handles.blue)

handles.output = [];
	handles.command = cell(45,1);
	% There still remains a lot of free cells (1-19) and (43-45).
	% Originally, this file bellonged to the GMT_pscoast_options file, where
	% cell 20 contained '-G'. I'll keep it that way but will only export cells 21
	% to 45. This will allow to use this same program to work with pscoast
	% -G -S -C options (and maybe some other GMT programs like psxy[z]).

	%a(:,:,1) = rand(20,20);    a(:,:,2) = rand(20,20);     a(:,:,3) = rand(20,20);
	a(:,:,1) = color_wheel(10,0);     % 21x21
	a(:,:,2) = color_wheel(10,1);
	a(:,:,3) = color_wheel(10,2);
	set(handles.push_OptionG_CustomColor,'CData',a)
	set(handles.push_OptionG_FgCustomColor,'CData',a)
	set(handles.push_OptionG_BgCustomColor,'CData',a)

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, findobj(hObject,'Style','Text'))
	%------------- END Pro look (3D) -----------------------------------------------------

	% Choose default command line output for paint_option_export
	handles.output = hObject;
	guidata(hObject, handles);

	set(hObject,'Vis','on');
	% UIWAIT makes paint_option_export wait for user response (see UIRESUME)
	uiwait(handles.figure1);

	handles = guidata(hObject);
	varargout{1} = handles.output;
	delete(handles.figure1);

function C = color_wheel(n,i)     % It will eventualy become that
	r = 0:n;		r = [r (n-1):-1:0]'/n;
	theta = (pi+i*2*pi/3)*(-n:n)/n;
	X = r*cos(theta);
	C = (X + 1) / 2;

% --------------------------------------------------------------------------------------------
function popup_Option_G_CB(hObject, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case 'Paint'
		handles.command{20} = '-G';   handles.command{21} = '';
		set(handles.edit_OptionG_gray,'Enable','on','Backgroundcolor',handles.gray)
		set(handles.edit_OptionG_R,'Enable','on','Backgroundcolor',handles.red)
		set(handles.edit_OptionG_G,'Enable','on','Backgroundcolor',handles.green)
		set(handles.edit_OptionG_B,'Enable','on','Backgroundcolor',handles.blue)
		set(handles.edit_OptionG_Foreground_R,'Enable','on','Backgroundcolor',handles.red)
		set(handles.edit_OptionG_Foreground_G,'Enable','on','Backgroundcolor',handles.green)
		set(handles.edit_OptionG_Foreground_B,'Enable','on','Backgroundcolor',handles.blue)
		set(handles.edit_OptionG_Background_R,'Enable','on','Backgroundcolor',handles.red)
		set(handles.edit_OptionG_Background_G,'Enable','on','Backgroundcolor',handles.green)
		set(handles.edit_OptionG_Background_B,'Enable','on','Backgroundcolor',handles.blue)
		set(handles.check_OptionG_BitReverse,'Enable','on')
		set(handles.edit_OptionG_dpi,'Enable','on','Backgroundcolor',[1 1 1])
		set(handles.listbox_OptionG_Patterns,'Enable','on')
		set(handles.push_OptionG_CustomColor,'Enable','on')
		set(handles.push_OptionG_ViewPatterns,'Enable','on')
		set(handles.push_OptionG_LoadRasFile,'Enable','on')
		set(handles.push_OptionG_FgCustomColor,'Enable','on')
		set(handles.push_OptionG_BgCustomColor,'Enable','on')
		set(handles.check_TransparentFg,'Enable','on')
		set(handles.check_TransparentBg,'Enable','on')    
    case 'Clip'
		handles.command{20} = '-G';   handles.command{21} = 'c';
		set(handles.edit_OptionG_gray,'Enable','off','Backgroundcolor',handles.ColorDes)
		set(handles.edit_OptionG_R,'Enable','off','Backgroundcolor',handles.ColorDes)
		set(handles.edit_OptionG_G,'Enable','off','Backgroundcolor',handles.ColorDes)
		set(handles.edit_OptionG_B,'Enable','off','Backgroundcolor',handles.ColorDes)
		set(handles.edit_OptionG_Foreground_R,'Enable','off','Backgroundcolor',handles.ColorDes)
		set(handles.edit_OptionG_Foreground_G,'Enable','off','Backgroundcolor',handles.ColorDes)
		set(handles.edit_OptionG_Foreground_B,'Enable','off','Backgroundcolor',handles.ColorDes)
		set(handles.edit_OptionG_Background_R,'Enable','off','Backgroundcolor',handles.ColorDes)
		set(handles.edit_OptionG_Background_G,'Enable','off','Backgroundcolor',handles.ColorDes)
		set(handles.edit_OptionG_Background_B,'Enable','off','Backgroundcolor',handles.ColorDes)
		set(handles.check_OptionG_BitReverse,'Enable','off')
		set(handles.edit_OptionG_dpi,'Enable','off','Backgroundcolor',handles.ColorDes)
		set(handles.listbox_OptionG_Patterns,'Enable','off')
		set(handles.push_OptionG_CustomColor,'Enable','off')
		set(handles.push_OptionG_ViewPatterns,'Enable','off')
		set(handles.push_OptionG_LoadRasFile,'Enable','off')
		set(handles.push_OptionG_FgCustomColor,'Enable','off')
		set(handles.push_OptionG_BgCustomColor,'Enable','off')
		set(handles.check_TransparentFg,'Value',0);
		set(handles.check_TransparentBg,'Value',0)
		set(handles.check_OptionG_BitReverse,'Value',0)
		set(handles.check_TransparentFg,'Enable','off')
		set(handles.check_TransparentBg,'Enable','off')
		clear_editBox(handles.edit_OptionG_gray)
		clear_editBox(handles.edit_OptionG_R)
		clear_editBox(handles.edit_OptionG_G)
		clear_editBox(handles.edit_OptionG_B)
		clear_editBox(handles.edit_OptionG_dpi)
		clear_editBox(handles.edit_OptionG_Background_R)
		clear_editBox(handles.edit_OptionG_Background_G)
		clear_editBox(handles.edit_OptionG_Background_B)
		clear_editBox(handles.edit_OptionG_Foreground_R)
		clear_editBox(handles.edit_OptionG_Foreground_G)
		clear_editBox(handles.edit_OptionG_Foreground_B)
		for i = 22:42   handles.command{i} = '';    end
end
set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
guidata(hObject, handles);

function clear_editBox(clean)
% Just clean what might be inside an edit box
set(clean, 'String', '')

% --------------------------------------------------------------------------------------------
function edit_OptionG_gray_CB(hObject, handles)
	xx = get(hObject,'String');
	if str2double(xx) > 255 errordlg('Gray value must be in the range 0-255','Error');
		set(handles.edit_OptionG_gray, 'String', '?');     xx = '';   end
	handles.command{22} = xx;     handles.command{21} = '';
	handles.command{20} = '-G';    set(handles.popup_Option_G, 'Value',2)
	for i = 23:42   handles.command{i} = '';  end
	toggle_ClearFg_CB(hObject, handles)
	toggle_ClearBg_CB(hObject, handles)
	set(handles.edit_OptionG_R, 'String', '');  set(handles.edit_OptionG_G, 'String', '');
	set(handles.edit_OptionG_B, 'String', '');
	set(handles.edit_OptionG_dpi, 'String', '');    set(handles.listbox_OptionG_Patterns, 'Value',1)
	set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
	guidata(hObject, handles);

% --------------------------------------------------------------------------------------------
function edit_OptionG_R_CB(hObject, handles)
	xx = get(hObject,'String');
	if str2double(xx) > 255 errordlg('Red value must be in the range 0-255','Error');
		set(handles.edit_OptionG_R, 'String', '?');     xx = '';   end
	handles.command{22} = xx;
	handles.command{20} = '-G';    set(handles.popup_Option_G, 'Value',2)
	handles.command{23} = '/';    handles.command{25} = '/';
	if isempty(get(handles.edit_OptionG_G,'String')) set(handles.edit_OptionG_G, 'String', '?'); end
	if isempty(get(handles.edit_OptionG_B,'String')) set(handles.edit_OptionG_B, 'String', '?'); end
	set(handles.edit_OptionG_gray, 'String', '');
	if strfind(xx, '?')
		xx= strrep(xx,'?','');    handles.command{22} = xx;    set(handles.edit_OptionG_R, 'String', xx);
	end
	set(handles.edit_OptionG_dpi, 'String', '');    set(handles.listbox_OptionG_Patterns, 'Value',1)
	set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
	guidata(hObject, handles);

% --------------------------------------------------------------------------------------------
function edit_OptionG_G_CB(hObject, handles)
	xx = get(hObject,'String');
	if str2double(xx) > 255 errordlg('Green value must be in the range 0-255','Error');
		set(handles.edit_OptionG_G, 'String', '?');     xx = '';   end
	handles.command{24} = xx;     set(handles.popup_Option_G, 'Value',2)
	handles.command{23} = '/';    handles.command{25} = '/';
	if isempty(get(handles.edit_OptionG_R,'String')) set(handles.edit_OptionG_R, 'String', '?'); end
	if isempty(get(handles.edit_OptionG_B,'String')) set(handles.edit_OptionG_B, 'String', '?'); end
	set(handles.edit_OptionG_gray, 'String', '');
	if strfind(xx, '?')
		xx= strrep(xx,'?','');    handles.command{24} = xx;    set(handles.edit_OptionG_G, 'String', xx);
	end
	set(handles.edit_OptionG_dpi, 'String', '');    set(handles.listbox_OptionG_Patterns, 'Value',1)
	set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
	guidata(hObject, handles);

% --------------------------------------------------------------------------------------------
function edit_OptionG_B_CB(hObject, handles)
	xx = get(hObject,'String');
	if str2double(xx) > 255 errordlg('Blue value must be in the range 0-255','Error');
		set(handles.edit_OptionG_B, 'String', '?');     xx = '';   end
	handles.command{26} = xx;     set(handles.popup_Option_G, 'Value',2)
	handles.command{23} = '/';    handles.command{25} = '/';
	if isempty(get(handles.edit_OptionG_R,'String')) set(handles.edit_OptionG_R, 'String', '?'); end
	if isempty(get(handles.edit_OptionG_G,'String')) set(handles.edit_OptionG_G, 'String', '?'); end
	set(handles.edit_OptionG_gray, 'String', '');
	if strfind(xx, '?')
		xx= strrep(xx,'?','');    handles.command{26} = xx;    set(handles.edit_OptionG_B, 'String', xx);
	end
	set(handles.edit_OptionG_dpi, 'String', '');    set(handles.listbox_OptionG_Patterns, 'Value',1)
	set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
	guidata(hObject, handles);

function push_OptionG_CustomColor_CB(hObject, handles)
handles.command{20} = '-G';    handles.command{21} = '';
set(handles.popup_Option_G, 'Value',2)
c = uisetcolor;
if length(c) > 1            % That is, if a color was selected
    c(1) = round(c(1)*255);     c(2) = round(c(2)*255);     c(3) = round(c(3)*255);
    handles.command{22} = num2str(c(1)); handles.command{24} = num2str(c(2)); handles.command{26} = num2str(c(3));
    handles.command{23} = '/';    handles.command{25} = '/';
    set(handles.edit_OptionG_R, 'String', num2str(c(1)));
    set(handles.edit_OptionG_G, 'String', num2str(c(2)));
    set(handles.edit_OptionG_B, 'String', num2str(c(3)));
    set(handles.listbox_OptionG_Patterns, 'Value',1)
    clear_editBox(handles.edit_OptionG_gray);   clear_editBox(handles.edit_OptionG_dpi)
    clear_editBox(handles.edit_OptionG_Background_R)
    clear_editBox(handles.edit_OptionG_Background_G)
    clear_editBox(handles.edit_OptionG_Background_B)
    clear_editBox(handles.edit_OptionG_Foreground_R)
    clear_editBox(handles.edit_OptionG_Foreground_G)
    clear_editBox(handles.edit_OptionG_Foreground_B)
    set(handles.check_TransparentFg,'Value',0);
    set(handles.check_TransparentBg,'Value',0)
    set(handles.check_OptionG_BitReverse,'Value',0)
    for i = 27:42   handles.command{i} = '';    end
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    guidata(hObject, handles);
end

% --------------------------------------------------------------------------------------------
function edit_OptionG_dpi_CB(hObject, handles)
xx = get(hObject,'String');
if isnan(str2double(xx))
    errordlg('DPI resolution must be a valid number','Error');
    clear_editBox(handles.edit_OptionG_dpi)
    return
end

if ~isempty(xx)
    handles.command{20} = '-G';    handles.command{21} = 'p';    set(handles.popup_Option_G, 'Value',2)
    set(handles.edit_OptionG_R, 'String', '');  set(handles.edit_OptionG_G, 'String', '');
    set(handles.edit_OptionG_B, 'String', '');  set(handles.edit_OptionG_gray, 'String', '');
    for i = 22:26   handles.command{i} = '';  end
    handles.command{27} = xx;     handles.command{28} = '/';
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    if get(handles.listbox_OptionG_Patterns, 'Value') ~= 1  % Keep the selected pattern
        handles.command{29} = num2str(get(hObject,'Value')-1);
    else
        set(handles.listbox_OptionG_Patterns, 'Value',2)    % Set the first as the default pattern
        handles.command{29} = '1';
    end
    listbox_OptionG_Patterns_CB(handles.listbox_OptionG_Patterns, handles)
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    guidata(hObject, handles);
else        % In case a previous value in the dpi box is removed
    clear_editBox(handles.edit_OptionG_Background_R)
    clear_editBox(handles.edit_OptionG_Background_G)
    clear_editBox(handles.edit_OptionG_Background_B)
    clear_editBox(handles.edit_OptionG_Foreground_R)
    clear_editBox(handles.edit_OptionG_Foreground_G)
    clear_editBox(handles.edit_OptionG_Foreground_B)
    set(handles.check_TransparentFg,'Value',0);
    set(handles.check_TransparentBg,'Value',0)
    set(handles.check_OptionG_BitReverse,'Value',0)
    set(handles.listbox_OptionG_Patterns, 'Value',1)
    set(handles.popup_Option_G, 'Value',1)
    for i = 20:42   handles.command{i} = '';    end
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}])
    guidata(hObject, handles)
end

% --------------------------------------------------------------------------------------------
function listbox_OptionG_Patterns_CB(hObject, handles)
val = get(hObject,'Value');     handles.command{29} = num2str(val-1);
set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
guidata(hObject, handles);
    
function push_OptionG_LoadRasFile_CB(hObject, handles)
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    [FileName,PathName] = uigetfile('*.ras','Select a Sun raster file');
    if ~isequal(FileName,0)
        handles.command{29} = [PathName FileName];
        set(handles.listbox_OptionG_Patterns, 'Value',1)
        set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
        guidata(hObject, handles);
    end
else
    msgbox('Must select pattern resolution dpi first.','Warning')
end

function check_OptionG_BitReverse_CB(hObject, handles)
% Hint: get(hObject,'Value') returns toggle state of check_OptionG_BitReverse
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    if get(hObject,'Value') handles.command{21} = 'P'; else handles.command{21} = 'p'; end
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    guidata(hObject, handles);
else
      set(hObject,'Value',0)
end

% --------------------------------------------------------------------------------------------
function edit_OptionG_Foreground_R_CB(hObject, handles)
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    xx = get(hObject,'String');
    if str2double(xx) > 255 errordlg('Red value must be in the range 0-255','Error');
        set(handles.edit_OptionG_Foreground_R, 'String', '?');     xx = '';   end
    handles.command{30} = ':';    handles.command{31} = 'F';   handles.command{32} = xx;
    handles.command{33} = '/';    handles.command{35} = '/';
    if isempty(get(handles.edit_OptionG_Foreground_G,'String')) set(handles.edit_OptionG_Foreground_G, 'String', '?'); end
    if isempty(get(handles.edit_OptionG_Foreground_B,'String')) set(handles.edit_OptionG_Foreground_B, 'String', '?'); end
    if strfind(xx, '?')
        xx= strrep(xx,'?','');    handles.command{32} = xx;    set(handles.edit_OptionG_Foreground_R, 'String', xx);
    end
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    guidata(hObject, handles);
else
    clear_editBox(handles.edit_OptionG_Foreground_R)
end

% --------------------------------------------------------------------------------------------
function edit_OptionG_Foreground_G_CB(hObject, handles)
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    xx = get(hObject,'String');
    if str2double(xx) > 255 errordlg('Green value must be in the range 0-255','Error');
        set(handles.edit_OptionG_Foreground_G, 'String', '?');     xx = '';   end
    handles.command{30} = ':';    handles.command{31} = 'F';   handles.command{34} = xx;
    handles.command{33} = '/';    handles.command{35} = '/';
    if isempty(get(handles.edit_OptionG_Foreground_R,'String')) set(handles.edit_OptionG_Foreground_R, 'String', '?'); end
    if isempty(get(handles.edit_OptionG_Foreground_B,'String')) set(handles.edit_OptionG_Foreground_B, 'String', '?'); end
    if strfind(xx, '?')
        xx= strrep(xx,'?','');    handles.command{34} = xx;    set(handles.edit_OptionG_Foreground_G, 'String', xx);
    end
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    guidata(hObject, handles);
else
    clear_editBox(handles.edit_OptionG_Foreground_G)
end

% --------------------------------------------------------------------------------------------
function edit_OptionG_Foreground_B_CB(hObject, handles)
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    xx = get(hObject,'String');
    if str2double(xx) > 255 errordlg('Green value must be in the range 0-255','Error');
        set(handles.edit_OptionG_Foreground_B, 'String', '?');     xx = '';   end
    handles.command{30} = ':';    handles.command{31} = 'F';   handles.command{36} = xx;
    handles.command{33} = '/';    handles.command{35} = '/';
    if isempty(get(handles.edit_OptionG_Foreground_R,'String')) set(handles.edit_OptionG_Foreground_R, 'String', '?'); end
    if isempty(get(handles.edit_OptionG_Foreground_G,'String')) set(handles.edit_OptionG_Foreground_G, 'String', '?'); end
    if strfind(xx, '?')
        xx= strrep(xx,'?','');    handles.command{36} = xx;    set(handles.edit_OptionG_Foreground_B, 'String', xx);
    end
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    guidata(hObject, handles);
else
    clear_editBox(handles.edit_OptionG_Foreground_B)
end

function push_OptionG_FgCustomColor_CB(hObject, handles)
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    c = uisetcolor;
    if length(c) > 1
        handles.command{30} = ':';    handles.command{31} = 'F';
        c(1) = round(c(1)*255);     c(2) = round(c(2)*255);     c(3) = round(c(3)*255);
        handles.command{32} = num2str(c(1)); handles.command{34} = num2str(c(2));
        handles.command{36} = num2str(c(3));
        handles.command{33} = '/';    handles.command{35} = '/';
        set(handles.edit_OptionG_Foreground_R, 'String', num2str(c(1)));
        set(handles.edit_OptionG_Foreground_G, 'String', num2str(c(2)));
        set(handles.edit_OptionG_Foreground_B, 'String', num2str(c(3)));
        set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
        guidata(hObject, handles);
    end
end

% --------------------------------------------------------------------------------------------
function edit_OptionG_Background_R_CB(hObject, handles)
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    xx = get(hObject,'String');
    if str2double(xx) > 255 errordlg('Red value must be in the range 0-255','Error');
        set(handles.edit_OptionG_Background_R, 'String', '?');     xx = '';   end
    handles.command{37} = 'B';   handles.command{38} = xx;
    handles.command{39} = '/';    handles.command{41} = '/';
    if isempty(get(handles.edit_OptionG_Background_G,'String')) set(handles.edit_OptionG_Background_G, 'String', '?'); end
    if isempty(get(handles.edit_OptionG_Background_B,'String')) set(handles.edit_OptionG_Background_B, 'String', '?'); end
    if strfind(xx, '?')
    xx= strrep(xx,'?','');    handles.command{38} = xx;    set(handles.edit_OptionG_Background_R, 'String', xx);
    end
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    guidata(hObject, handles);
else
    clear_editBox(handles.edit_OptionG_Background_R)
end

% --------------------------------------------------------------------------------------------
function edit_OptionG_Background_G_CB(hObject, handles)
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    xx = get(hObject,'String');
    if str2double(xx) > 255 errordlg('Green value must be in the range 0-255','Error');
        set(handles.edit_OptionG_Background_G, 'String', '?');     xx = '';   end
    handles.command{37} = 'B';   handles.command{40} = xx;
    handles.command{39} = '/';    handles.command{41} = '/';
    if isempty(get(handles.edit_OptionG_Background_R,'String')) set(handles.edit_OptionG_Background_R, 'String', '?'); end
    if isempty(get(handles.edit_OptionG_Background_B,'String')) set(handles.edit_OptionG_Background_B, 'String', '?'); end
    if strfind(xx, '?')
        xx= strrep(xx,'?','');    handles.command{40} = xx;    set(handles.edit_OptionG_Background_G, 'String', xx);
    end
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    guidata(hObject, handles);
else
    clear_editBox(handles.edit_OptionG_Background_G)
end

% --------------------------------------------------------------------------------------------
function edit_OptionG_Background_B_CB(hObject, handles)
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    xx = get(hObject,'String');
    if str2double(xx) > 255 errordlg('Red value must be in the range 0-255','Error');
        set(handles.edit_OptionG_Background_R, 'String', '?');     xx = '';   end
    handles.command{37} = 'B';   handles.command{42} = xx;
    handles.command{39} = '/';    handles.command{41} = '/';
    if isempty(get(handles.edit_OptionG_Background_R,'String')) set(handles.edit_OptionG_Background_R, 'String', '?'); end
    if isempty(get(handles.edit_OptionG_Background_G,'String')) set(handles.edit_OptionG_Background_G, 'String', '?'); end
    if strfind(xx, '?')
        xx= strrep(xx,'?','');    handles.command{42} = xx;    set(handles.edit_OptionG_Background_B, 'String', xx);
    end
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    guidata(hObject, handles);
else
    clear_editBox(handles.edit_OptionG_Background_B)
end

function push_OptionG_BgCustomColor_CB(hObject, handles)
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    c = uisetcolor;
    if length(c) > 1
        handles.command{30} = ':';    handles.command{37} = 'B';
        c(1) = round(c(1)*255);     c(2) = round(c(2)*255);     c(3) = round(c(3)*255);
        handles.command{38} = num2str(c(1)); handles.command{40} = num2str(c(2)); handles.command{42} = num2str(c(3));
        handles.command{39} = '/';    handles.command{41} = '/';
        set(handles.edit_OptionG_Background_R, 'String', num2str(c(1)));
        set(handles.edit_OptionG_Background_G, 'String', num2str(c(2)));
        set(handles.edit_OptionG_Background_B, 'String', num2str(c(3)));
        set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
        guidata(hObject, handles);
    end
end

% --- Executes on button press in push_OptionG_ViewPatterns.
function push_OptionG_ViewPatterns_CB(hObject, handles)
	%display_patterns
	mirone(handles.patternFile)

function toggle_ClearFg_CB(hObject, handles)
	for i = 31:36   handles.command{i} = '';    end
	if isempty(handles.command{37}) handles.command{30} = '';    end
	set(handles.edit_OptionG_Foreground_R, 'String', '');
	set(handles.edit_OptionG_Foreground_G, 'String', '');
	set(handles.edit_OptionG_Foreground_B, 'String', '');
	set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);  set(hObject,'Value',0)
	guidata(hObject, handles);

function check_TransparentFg_CB(hObject, handles)
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    if get(hObject,'Value')
        set(handles.edit_OptionG_Foreground_R, 'String', '');
        set(handles.edit_OptionG_Foreground_G, 'String', '');
        set(handles.edit_OptionG_Foreground_B, 'String', '');
        for i = 32:36   handles.command{i} = '';    end
        handles.command{30} = ':';    handles.command{31} = 'F';    handles.command{32} = '-';
    end
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    guidata(hObject, handles);
else
    set(hObject,'Value',0)
end

function check_TransparentBg_CB(hObject, handles)
if ~isempty(get(handles.edit_OptionG_dpi,'String'))
    if get(hObject,'Value')
        set(handles.edit_OptionG_Background_R, 'String', '');
        set(handles.edit_OptionG_Background_G, 'String', '');
        set(handles.edit_OptionG_Background_B, 'String', '');
        for i = 38:42   handles.command{i} = '';    end
        handles.command{30} = ':';    handles.command{37} = 'B';    handles.command{38} = '-';
    end
    set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);
    guidata(hObject, handles);
else
    set(hObject,'Value',0)
end

function toggle_ClearBg_CB(hObject, handles)
	for i = 37:42   handles.command{i} = '';    end
	if isempty(handles.command{31}) handles.command{30} = '';    end
	set(handles.edit_OptionG_Background_R, 'String', '');
	set(handles.edit_OptionG_Background_G, 'String', '');
	set(handles.edit_OptionG_Background_B, 'String', '');
	set(handles.edit_ShowCommand, 'String', [handles.command{21:end}]);  set(hObject,'Value',0)
	guidata(hObject, handles);

% --------------------------------------------------------------------------------------------
function push_HelpFillColor_CB(hObject, handles)
message = {'Specify a gray shade (0-255) or a color (r/g/b in the 0-255 range). The'
    'second form also allows us to use a predefined bit-image pattern (but see'
    'the help box bellow). Instead of color selecting, you can also chose'
    '"clipping" in first the popup menu. In this case all other options will be'
    'desabled. Notice that tooltips are available for all choices and should be'
    'almost enough to explain the meaning of each option.'};
helpdlg(message,'Help on area color fill option');

function push_HelpFillPattern_CB(hObject, handles)
message = {'Instead of a color you can select here to fill the area with a pattern.'
    'In that case you must first select the pattern resolution in dpi (in first'
    'box), than choose a pattern number (default is first one on the list)'
    'The "Pattern" can either be a number in the range 1-90 or the name of'
    'a 1-, 8-, or 24-bit Sun raster file. The dpi parameter sets the resolution'
    'of this image on the page; the area fill is thus made up of a series of'
    'these "tiles". Specifying dpi as 0 will result in highest resolution'
    'obtainable given the present dpi setting in .gmtdefaults. By selecting'
    '"Bit reverse" the image will be bit-reversed, i.e., white and black areas'
    'will be interchanged (only applies to 1-bit images or predefined bit-image'
    'patterns). For these patterns and other 1-bit images one may specify'
    'alternative background and foreground colors by means of the three'
    'RedGreenBlue colored boxes. If applyied, that will replace the default'
    'white and black pixels, respectively. Checking "Transparent" yields a'
    'transparent image where only the back- or foreground pixels will be'
    'painted. Due to PostScript implementation limitations the rasterimages'
    'used here must be less than 146 x 146 pixels in size; for larger images'
    'see psimage. The format of Sun raster files is outlined in GMT manual'
    'Appendix B. Note that under PostScript Level 1 the patterns are filled by'
    'using the polygon as a clip path. Complex clip paths may require more'
    'memory than the PostScript interpreter has been assigned. There is'
    'therefore the possibility that some PostScript interpreters (especially'
    'those supplied with older laserwriters) will run out of memory and abort.'
    'Should that occur we recommend that you use a regular grayshade fill'
    'instead of the patterns. Installing more memory in your printer may or may'
    'not solve the problem!'};
helpdlg(message,'Help on area pattern fill option');

function push_Cancel_CB(hObject, handles)
	handles.output = '';        % User gave up, return nothing
	guidata(hObject, handles);
	uiresume(handles.figure1);

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
function paint_option_LayoutFcn(h1)
set(h1,...
'PaperUnits','centimeters',...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','paint_option',...
'NumberTitle','off',...
'Position',[265.768111202607 265.768111202607 437 230],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 44 345 101],...
'Style','frame',...
'Tag','frame2');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 159 345 51],...
'Style','frame',...
'Tag','frame1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'Position',[20 175 72 19],...
'String',{''; 'Paint'; 'Clip' },...
'Style','popupmenu',...
'Tooltip','Select painting or clipping',...
'Value',1,...
'Tag','popup_Option_G');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'HorizontalAlignment','left',...
'Position',[110 174 30 19],...
'Style','edit',...
'Tooltip','Gray shade component (0-255)',...
'Tag','edit_OptionG_gray');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'HorizontalAlignment','left',...
'Position',[160 174 30 19],...
'Style','edit',...
'Tooltip','Red component (0-255)',...
'Tag','edit_OptionG_R');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'HorizontalAlignment','left',...
'Position',[190 174 30 19],...
'Style','edit',...
'Tooltip','Green component (0-255)',...
'Tag','edit_OptionG_G');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'HorizontalAlignment','left',...
'Position',[220 174 30 19],...
'Style','edit',...
'Tooltip','Blue component (0-255)',...
'Tag','edit_OptionG_B');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'Position',[261 174 20 19],...
'Tooltip','Interactive color selection',...
'Tag','push_OptionG_CustomColor');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[321 172 21 23],...
'String','?',...
'Tooltip','Help on collor painting or clipping',...
'Tag','push_HelpFillColor');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'HorizontalAlignment','left',...
'Position',[20 106 30 19],...
'Style','edit',...
'Tooltip','Set pattern resolution in dpi',...
'Tag','edit_OptionG_dpi');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'Position',[53 106 41 19],...
'String',{  ''; '1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; '11'; '12'; '13'; '14'; '15'; '16'; '17'; '18'; '19'; '20'; '21'; '22'; '23'; '24'; '25'; '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '36'; '37'; '38'; '39'; '40'; '41'; '42'; '43'; '44'; '45'; '46'; '47'; '48'; '49'; '50'; '51'; '52'; '53'; '54'; '55'; '56'; '57'; '58'; '59'; '60'; '61'; '62'; '63'; '64'; '65'; '66'; '67'; '68'; '69'; '70'; '71'; '72'; '73'; '74'; '75'; '76'; '77'; '78'; '79'; '80'; '81'; '82'; '83'; '84'; '85'; '86'; '87'; '88'; '89'; '90' },...
'Style','listbox',...
'Tooltip','Select one of the pre-deffined 90 patterns',...
'Value',1,...
'Tag','listbox_OptionG_Patterns');

uicontrol('Parent',h1,...
'Call',{@paint_uiCB,h1,'push_OptionG_ViewPatterns_CB'},...
'Position',[105 106 71 20],...
'String','View Patterns',...
'Tag','push_OptionG_ViewPatterns');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'Position',[180 106 66 20],...
'String','Sun Ras file',...
'Tooltip','Load a Sun raster file that will be used as a pattern',...
'Tag','push_OptionG_LoadRasFile');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'Position',[258 109 71 15],...
'String','Bit reverse',...
'Style','checkbox',...
'Tooltip','Toggle black pixels to white and vice-versa',...
'Tag','check_OptionG_BitReverse');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'HorizontalAlignment','left',...
'Position',[20 79 30 19],...
'Style','edit',...
'Tooltip','Red component (0-255)',...
'Tag','edit_OptionG_Foreground_R');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'HorizontalAlignment','left',...
'Position',[50 79 30 19],...
'Style','edit',...
'Tooltip','Green component (0-255)',...
'Tag','edit_OptionG_Foreground_G');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'HorizontalAlignment','left',...
'Position',[80 79 30 19],...
'Style','edit',...
'Tooltip','Blue component (0-255)',...
'Tag','edit_OptionG_Foreground_B');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'Position',[120 79 20 19],...
'Tooltip','Foreground interactive color selection',...
'Tag','push_OptionG_FgCustomColor');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'HorizontalAlignment','left',...
'Position',[181 79 30 19],...
'Style','edit',...
'Tooltip','Red component (0-255)',...
'Tag','edit_OptionG_Background_R');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'HorizontalAlignment','left',...
'Position',[211 79 30 19],...
'Style','edit',...
'Tooltip','Green component (0-255)',...
'Tag','edit_OptionG_Background_G');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@paint_uiCB,...
'HorizontalAlignment','left',...
'Position',[241 79 30 19],...
'Style','edit',...
'Tooltip','Blue component (0-255)',...
'Tag','edit_OptionG_Background_B');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'Position',[281 79 20 19],...
'Tooltip','Background interactive color selection',...
'Tag','push_OptionG_BgCustomColor');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'Position',[20 53 79 15],...
'String','Transparent',...
'Style','checkbox',...
'Tooltip','Make Foreground pattern transparent',...
'Tag','check_TransparentFg');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'Position',[111 52 30 20],...
'String','Clear',...
'Style','togglebutton',...
'Tooltip','Remove color in Foreground pattern',...
'Tag','toggle_ClearFg');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'Position',[181 53 79 15],...
'String','Transparent',...
'Style','checkbox',...
'Tooltip','Make Background pattern transparent',...
'Tag','check_TransparentBg');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'Position',[272 52 30 20],...
'String','Clear',...
'Style','togglebutton',...
'Tooltip','Remove color in Background pattern',...
'Tag','toggle_ClearBg');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[323 64 21 37],...
'String','?',...
'Tooltip','Help on pattern fill',...
'Tag','push_HelpFillPattern');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[110 202 121 15],...
'String','Set Color for Painting',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[112 137 131 15],...
'String','OR Fill area with a pattern',...
'Style','text');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'Position',[364 45 66 23],...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1,...
'Call',@paint_uiCB,...
'Position',[364 75 66 23],...
'String','Cancel',...
'Tag','push_Cancel');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[10 9 421 23.4],...
'Style','edit',...
'Tooltip','Show GMT command corresponding to selected options',...
'Tag','edit_ShowCommand');

function paint_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

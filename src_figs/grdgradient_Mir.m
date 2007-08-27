function varargout = grdgradient_Mir(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to grdgradient_Mir (see VARARGIN)
%   The output is a structure with one or more of the following fields (or empty):
%   out.opt_A   -> contains the azimuth(s). It has the form '-Aaz1[/az2]'
%   out.opt_D   -> Direction of gradient '-D<ang>[c][o][n]'
%   out.opt_N   -> Normalization '-N[e][t][amp][/sigma][/offset]'
%   out.opt_L   -> if boundary conditions were selected

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
 
	hObject = figure('Tag','figure1','Visible','off');
	grdgradient_Mir_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'center')

	if ~isempty(varargin)
		handMir  = varargin{1};
		handles.Z = getappdata(handMir.figure1,'dem_z');
	else
        errordlg('GRDGRADIENT: wrong number of arguments.','Error')
        delete(hObject);    return
	end
    
	if (handMir.no_file)
		errordlg('GRDGRADIENT: You didn''t even load a file. What are you expecting then?','ERROR')
        delete(hObject);    return
	end
	if (~handMir.validGrid)
        errordlg('GRDGRADIENT: This operation is deffined only for images derived from DEM grids.','ERROR')
        delete(hObject);    return
	end
	if (isempty(handles.Z))
        errordlg('GRDGRADIENT: Grid was not saved in memory. Increase "Grid max size" and start over.','ERROR')
        delete(hObject);    return
	end

    handles.hMirFig = handMir.figure1;
	handles.head = handMir.head;
	handles.geog = handMir.geog;

    handles.bd_cond = [];

	% Generate a list of 360 lines for the listboxes containing azim & azim2
	az = cell(361,1);       azz = (-1:359)';  az = num2cell(azz);    az{1} = ' ';
	set(handles.listbox_azim1,'String',az);
	set(handles.listbox_azim2,'String',az);

	% Give a Pro look (3D) to the frame boxes 
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
	h_t = handles.text2;
	for i=1:length(h_t)
        usr_d = get(h_t(i),'UserData');
        t_size = get(h_t(i),'Position');   t_str = get(h_t(i),'String');    fw = get(h_t(i),'FontWeight');
        bgc = get (h_t(i),'BackgroundColor');   fgc = get (h_t(i),'ForegroundColor');
        uicontrol('Parent',hObject, 'Style','text', 'Position',t_size,'String',t_str, ...
            'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw,'UserData',usr_d);
	end
	delete(h_t)

	guidata(hObject, handles);

	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% ------------------------------------------------------------------------------------
function listbox_azim1_Callback(hObject, handles)
val = get(hObject,'Value');
if val == 1
    set(handles.checkbox_direction_grad,'Enable','on')
    set(handles.checkbox_magnitude_grad,'Value',0)
    set(handles.listbox_azim2,'Value',1)
else
    set(handles.checkbox_direction_grad,'Value',0);
    set(handles.checkbox_direction_grad,'Enable','off')
    set(handles.popup_direction_grad,'Value',1)
    set(handles.popup_direction_grad,'Enable','off')
end
% The rest will be done by the OK button

% ------------------------------------------------------------------------------------
function listbox_azim2_Callback(hObject, handles)
val = get(hObject,'Value');
if (val > 1 & get(handles.listbox_azim1,'Value') == 1)
    set(hObject,'Value',1)
end
% The rest will be done by the OK button

% ------------------------------------------------------------------------------------
function pushbutton_Help_Azim_Callback(hObject, handles)
message = {'Azimuthal direction for a directional derivative; "Azim1" is the angle in'
    'the x,y plane measured in degrees positive clockwise from north (the +y'
    'direction) toward east (the +x direction). The negative of the directional'
    'derivative, -[dz/dx*sin(azim1) + dz/dy*cos(azim1)], is found; negation'
    'yields positive values when the slope of z(x,y) is downhill in the azim1'
    'direction, the correct sense for shading the illumination of an image (see'
    'grdimage and grdview) by a light source above the x,y plane shining from'
    'the azim1 direction. Optionally, supply two azimuths, -Aazim1/azim2, in'
    'which case the gradients in each of these directions are calculated and'
    'the one larger in magnitude is retained; this is useful for illuminating'
    'data with two directions of lineated structures, e.g. -A0/270 illuminates'
    'from the north (top) and west (left).'};
helpdlg(message,'Help -A option');

% ------------------------------------------------------------------------------------
function checkbox_direction_grad_Callback(hObject, handles)
if get(hObject,'Value')
    set(handles.listbox_azim1,'Value',1);       set(handles.listbox_azim2,'Value',1);
    set(handles.listbox_azim1,'Enable','on');   set(handles.listbox_azim2,'Enable','on');
    set(handles.popup_direction_grad,'Enable','on');
else
    set(handles.listbox_azim1,'Value',1);       set(handles.listbox_azim2,'Value',1);
    set(handles.listbox_azim1,'Enable','off');  set(handles.listbox_azim2,'Enable','off');
    set(handles.checkbox_magnitude_grad,'Value',0);
    set(handles.popup_direction_grad,'Value',1);
    set(handles.popup_direction_grad,'Enable','off');
end

% ------------------------------------------------------------------------------------
function pushbutton_Help_direction_grad_Callback(hObject, handles)
message = {'Find the direction of the gradient of the data. By default, the directions'
    'are measured clockwise from north. Select "trignometric angles" to use'
    'conventional Cartesian angles measured counterclockwise from the'
    'positive x (east) direction. Select "0-180 degrees orientations" to report'
    'orientations (0-180) rather than directions (0-360). Select "add 90'
    'degrees to all angles" to add 90 degrees to all angles (e.g., to give'
    'orientation of lineated features).'};
helpdlg(message,'Help -D option');

% ------------------------------------------------------------------------------------
function popup_BoundaryCondition_Callback(hObject, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	switch str{val};
		case ' ',        handles.bd_cond = [];
		case '',         handles.bd_cond = [];
		case 'x',        handles.bd_cond = ['-Lx'];
		case 'y',        handles.bd_cond = ['-Ly'];
		case 'xy',       handles.bd_cond = ['-Lxy']; 
		case 'g',        handles.bd_cond = ['-Lg'];
	end
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function pushbutton_Help_L_Callback(hObject, handles)
message = {'Boundary condition flag may be x or y or xy indicating data is periodic'
           'in range of x or y or both, or flag may be g indicating geographical'
           'conditions (x and y may be lon and lat). [Default is no boundary conditions].'};
helpdlg(message,'Help boundary condition');

% ------------------------------------------------------------------------------------
function popup_normalization_Callback(hObject, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case ''
        set(handles.edit_normalization_amp,'String','')
        set(handles.edit_normalization_sigma,'String','')
        set(handles.edit_normalization_offset,'String','')
end
% The rest will be done by the OK button

% ------------------------------------------------------------------------------------
function edit_normalization_amp_Callback(hObject, handles)
xx = get(hObject,'String');
if ~isempty(xx) & (get(handles.popup_normalization,'Value') ~= 1)
    if ~(str2double(xx) > 0 & ~isnan(str2double(xx)))
        set(hObject,'String','')
    end
else
    set(hObject,'String','')
    set(handles.edit_normalization_sigma,'String','');
    set(handles.edit_normalization_offset,'String','')
end

% ------------------------------------------------------------------------------------
function edit_normalization_sigma_Callback(hObject, handles)
xx = get(hObject,'String');
if ~isempty(xx) & ~isempty(get(handles.edit_normalization_amp,'String'))
    if ~(str2double(xx) > 0 & ~isnan(str2double(xx)))
        set(hObject,'String','')
    end
else
    set(hObject,'String','')
    set(handles.edit_normalization_offset,'String','')
end

% ------------------------------------------------------------------------------------
function edit_normalization_offset_Callback(hObject, handles)
xx = get(hObject,'String');
if ~isempty(xx) & ~isempty(get(handles.edit_normalization_sigma,'String'))
    if ~(str2double(xx) > 0 & ~isnan(str2double(xx)))
        set(hObject,'String','')
    end
else
    set(hObject,'String','')
end

% ------------------------------------------------------------------------------------
function pushbutton_Help_N_grad_Callback(hObject, handles)
message = {'Normalization. [Default: no normalization.] The actual gradients g are'
    'offset and scaled to produce normalized gradients gn with a maximum output'
    'magnitude of amp. If amp is not given, default amp = 1. If offset is not'
    'given, it is set to the average of g. "Simple" yields gn = amp * (g'
    '- offset) / max(abs(g - offset)). "Laplace" normalizes using a cumulative'
    'Laplace distribution yielding gn = amp * (1.0 - exp(sqrt(2) * (g - offset)'
    '/ sigma)) where sigma is estimated using the L1 norm of (g - offset) if it'
    'is not given. "Cauchy" normalizes using a cumulative Cauchy distribution'
    'yielding gn = (2 * amp / PI) * atan( (g - offset)/sigma) where sigma is'
    'estimated using the L2 norm of (g - offset) if it is not given.'};
helpdlg(message,'Help normalization');

% ------------------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, handles)
	delete(handles.figure1);

% ------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, handles)
	a_set = 0;  d_set = 0;  n_set = 0;  s_set = 0;
	opt_A = ' ';     opt_D = ' ';     opt_L = ' ';     opt_M = ' ';     opt_N = ' ';

	% Check for Azim selection(s)
	val1 = get(handles.listbox_azim1,'Value');
	if (val1 > 1)
        a_set = 1;
        opt_A = ['-A' num2str(val1-2)];
        val2 = get(handles.listbox_azim2,'Value');
        if (val2 > 1)       % Second azimuth was selected also
            opt_A = [opt_A '/' num2str(val2-2)];
        end
	end

	% Check for the Direction Gradient option
	if (get(handles.checkbox_direction_grad,'Value'))
        opt_D = '-D';
        d_set = 1;
        % Check also for sub-options
        val = get(handles.popup_direction_grad,'Value');
        str = get(handles.popup_direction_grad, 'String');
        if (val > 1)
			switch str{val};
                case 'trignometric angles'
                    opt_D = [opt_D 'c'];
                case '0-180 degrees orientations'
                    opt_D = [opt_D 'o'];
                case 'add 90 degrees to all angles'
                    opt_D = [opt_D 'n'];
			end
        end
	end

	if (get(handles.checkbox_magnitude_grad,'Value'))
        s_set = 1;
	end

	% Check for the boundary condition option
	if (~isempty(handles.bd_cond)),    opt_L = handles.bd_cond;     end

	% Check for the Normalized option
	val = get(handles.popup_normalization,'Value');
	str = get(handles.popup_normalization, 'String');
	if (val > 1)
        n_set = 1;    opt_N = '-N';
        switch str{val};
            case 'Simple'
                opt_N = '-N';
            case 'Laplace'
                opt_N = '-Ne';
            case 'Cauchy'
                opt_N = '-Nt';
        end
	end

	% Check for normalization parameters (Relyes on no error in these params selection)
	if (n_set)
		xx = get(handles.edit_normalization_amp,'String');
		if (~isempty(xx))   opt_N = [opt_N xx];    end
		xx = get(handles.edit_normalization_sigma,'String');
		if (~isempty(xx))   opt_N = [opt_N '/' xx];    end
		xx = get(handles.edit_normalization_offset,'String');
		if (~isempty(xx))   opt_N = [opt_N '/' xx];    end    
	end

	% Consistency check
	if ~(a_set || d_set)
        errordlg('Must select one of "Horizontal Light" or "Gradient Direction"','Error');  return
	end
	if (s_set && ~d_set)
        errordlg('Computing scalar magnitudes implyies selecting also "Gradient Direction"','Error');  return
	end
	if ~(a_set || d_set || n_set || s_set)
        errordlg('You haven''t select anything usefull to do.','Chico Clever');   return
	end
    
	if (handles.geog),           opt_M = '-M';        end

	set(handles.figure1,'pointer','watch');     set(handles.hMirFig,'pointer','watch')
	newZ = grdgradient_m(single(-double(handles.Z)),handles.head,opt_A,opt_D,opt_L,opt_M,opt_N);     % should be clever
	zz = grdutils(newZ,'-L');       handles.head(5:6) = zz(1:2);
	set(handles.figure1,'pointer','arrow');     set(handles.hMirFig,'pointer','arrow')
    tmp.X = linspace(handles.head(1),handles.head(2),size(handles.Z,2));
    tmp.Y = linspace(handles.head(3),handles.head(4),size(handles.Z,1));
    tmp.head = handles.head;
    tmp.name = 'Gradient grid';
    mirone(newZ,tmp);
    figure(handles.figure1)         % Don't let this figure forgotten behind the newly created one

% --------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, event)
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(hObject);
	end

% --- Creates and returns a handle to the GUI figure. 
function grdgradient_Mir_LayoutFcn(h1);
set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','grdgradient',...
'NumberTitle','off',...
'Position',[520 599 490 201],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[270 110 210 81],...
'Style','frame',...
'Tag','frame1');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 110 231 81],...
'Style','frame',...
'Tag','frame2');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdgradient_Mir_uicallback,h1,'listbox_azim1_Callback'},...
'Position',[55 118 50 61],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_azim1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdgradient_Mir_uicallback,h1,'listbox_azim2_Callback'},...
'Position',[150 118 50 61],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_azim2');

uicontrol('Parent',h1,...
'Callback',{@grdgradient_Mir_uicallback,h1,'pushbutton_Help_Azim_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[212 122 18 33],...
'String','?',...
'Tag','pushbutton_Help_Azim');

uicontrol('Parent',h1,...
'Callback',{@grdgradient_Mir_uicallback,h1,'checkbox_direction_grad_Callback'},...
'Position',[279 146 103 17],...
'String','Gradient direction',...
'Style','checkbox',...
'TooltipString','Find  the direction of the gradient of the data',...
'Tag','checkbox_direction_grad');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[279 121 168 22],...
'String',{'azimuthal direction'; 'trignometric angles'; '0-180 degrees orientations'; 'add 90 degrees to all angles'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_direction_grad');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[281 82 95 16],...
'String','Boundary condition',...
'Style','text',...
'Tag','text1');

uicontrol('Parent',h1,...
'Callback',{@grdgradient_Mir_uicallback,h1,'pushbutton_Help_direction_grad_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[455 121 18 33],...
'String','?',...
'Tag','pushbutton_Help_direction_grad');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdgradient_Mir_uicallback,h1,'popup_BoundaryCondition_Callback'},...
'HorizontalAlignment','right',...
'Position',[378 79 47 21],...
'String',{  ''; 'x'; 'y'; 'xy'; 'g' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_BoundaryCondition');

uicontrol('Parent',h1,...
'Callback',{@grdgradient_Mir_uicallback,h1,'pushbutton_Help_L_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[442 79 21 23],...
'String','?',...
'Tag','pushbutton_Help_L');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[21 183 120 15],...
'String','Horizontal Light Angles',...
'Style','text',...
'Tag','text2');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[19 143 33 13],...
'String','Azim1',...
'Style','text',...
'Tag','text3');

uicontrol('Parent',h1,...
'Callback',{@grdgradient_Mir_uicallback,h1,'pushbutton_cancel_Callback'},...
'Position',[341 9 66 23],...
'String','Cancel',...
'Tag','pushbutton_cancel');

uicontrol('Parent',h1,...
'Callback',{@grdgradient_Mir_uicallback,h1,'pushbutton_OK_Callback'},...
'Position',[416 9 66 23],...
'String','OK',...
'Tag','pushbutton_OK');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdgradient_Mir_uicallback,h1,'popup_normalization_Callback'},...
'Position',[78 48 72 22],...
'String',{  ''; 'Simple'; 'Laplace'; 'Cauchy' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_normalization');

 uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdgradient_Mir_uicallback,h1,'edit_normalization_amp_Callback'},...
'Position',[157 48 47 21],...
'Style','edit',...
'TooltipString','Amp',...
'Tag','edit_normalization_amp');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdgradient_Mir_uicallback,h1,'edit_normalization_sigma_Callback'},...
'Position',[209 48 47 21],...
'Style','edit',...
'TooltipString','Sigma',...
'Tag','edit_normalization_sigma');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdgradient_Mir_uicallback,h1,'edit_normalization_offset_Callback'},...
'Position',[261 48 47 21],...
'Style','edit',...
'TooltipString','offset',...
'Tag','edit_normalization_offset');

uicontrol('Parent',h1,...
'Callback',{@grdgradient_Mir_uicallback,h1,'pushbutton_Help_N_grad_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[317 47 21 23],...
'String','?',...
'Tag','pushbutton_Help_N_grad');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[12 50 64 15],...
'String','Normalization',...
'Style','text',...
'Tag','text5');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[114 142 33 15],...
'String','Azim2',...
'Style','text',...
'Tag','text4');

uicontrol('Parent',h1,...
'Position',[279 168 151 15],...
'String','Compute scalar magnitudes',...
'Style','checkbox',...
'TooltipString','Compute the scalar magnitudes of  gradient  vectors instead of the directional derivative',...
'Tag','checkbox_magnitude_grad');

function grdgradient_Mir_uicallback(hObject, event, h1, callback_name)
	% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

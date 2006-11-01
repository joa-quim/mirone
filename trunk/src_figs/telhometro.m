function varargout = telhometro(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to telhometro (see VARARGIN)

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
handles = guihandles(hObject);
guidata(hObject, handles);
telhometro_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

%#function telha_m choosebox

movegui(hObject,'center');
global home_dir;    home_dir = pwd;
handles.h_active_line_str = findobj(handles.figure1,'Tag','text_activeLine');      % Get this handle
handles.path_data = [home_dir filesep 'data' filesep];
handles.path_tmp = [home_dir filesep 'tmp' filesep];
handles.path_continent = [home_dir filesep 'continents' filesep];
handles.h_line_orig = [];

if (~isempty(varargin))
    handles.h_calling_fig = varargin{1};
    handles.mirone_axes = get(varargin{1},'CurrentAxes');
    if (length(varargin) == 2)          % Called with the line handle in argument
        handles.h_line_orig = varargin{2};
        set(handles.h_active_line_str,'String','GOT A LINE TO WORK WITH')
    end
else
    errordlg('TELHOMETRO: wrong number of arguments.','Error')
    delete(hObject)
    return
end

% Choose default command line output for telhometro_export
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes earthquakes wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = telhometro_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = telhometro_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------------------
function edit_polesFile_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if isempty(fname)    return;    end
% Let the pushbutton_readPolesFile_Callback do all the work
pushbutton_readPolesFile_Callback(hObject,[],guidata(gcbo),fname)

% -------------------------------------------------------------------------------------
function pushbutton_readPolesFile_Callback(hObject, eventdata, handles, opt)
% Get poles file name
if (nargin == 4)    fname = opt;
else                opt = [];
end

if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
    cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
    last_dir = cfig_handles.last_dir;
    home = cfig_handles.home_dir;
else
    last_dir = [];
end

if (isempty(opt))           % Otherwise we already know fname from the 4th input argument
    if (~isempty(last_dir)),    cd(last_dir);   end
	str1 = {'*.stg;*.dat;*.DAT', 'Data files (*.stg,*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = uigetfile(str1,'Select poles file');  pause(0.05)
    if (~isempty(last_dir)),    cd(home);   end
	if isequal(FileName,0)      return;    end
    fname = [PathName FileName];
end
set(handles.edit_polesFile,'String',fname)

% --------------------------------------------------------------------
function checkbox_revertRot_Callback(hObject, eventdata, handles)
% Nothing to do here. The compute callback will take care

% -------------------------------------------------------------------------------------
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1)

% -----------------------------------------------------------------------------------
function edit_timeStart_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < 0)     set(hObject, 'String','0');     end

% -----------------------------------------------------------------------------------
function edit_timeEnd_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < 0)     set(hObject, 'String','');      end

% --------------------------------------------------------------------
function checkbox_ridgeStartTime_Callback(hObject, eventdata, handles)
% Nothing to do here. The compute callback will take care

% -------------------------------------------------------------------------------------
function pushbutton_compute_Callback(hObject, eventdata, handles)

if (isempty(handles.h_line_orig))
    errordlg('Will you be so kind to let me know what line should I rotate?','Unknown target')
    return
end
lt = getappdata(handles.h_calling_fig,'DefLineThick');
lc = getappdata(handles.h_calling_fig,'DefLineColor');

poles_name = get(handles.edit_polesFile,'String');
if (isempty(poles_name))
    errordlg('No stage poles provided','Error');    return
end

axes(handles.mirone_axes)       % Make the Mirone axes the CurrentAxes
x = get(handles.h_line_orig,'XData');       y = get(handles.h_line_orig,'YData');
linha = [x(:) y(:)];
opt_E = ['-E' poles_name];
if (get(handles.checkbox_revertRot,'Value'))
    opt_I = '-I';
else
    opt_I = ' ';
end

t0 = str2double(get(handles.edit_timeStart,'String'));
if (t0 > 0)     opt_T = ['-T' num2str(t0)];
else            opt_T = ' ';
end

t1 = str2double(get(handles.edit_timeEnd,'String'));
if (t1 > 0)     opt_N = ['-N' num2str(t1)];
else            opt_N = ' ';
end

if (t1 < t0)
    errordlg('You are a bit confused with ages. Time End is < Time Start ','Error')
    return
end

[out_x,out_y,first_mag,n_flow] = telha_m(linha, opt_E, opt_I, '-B', opt_T, opt_N);
clear mex;

% if (get(handles.checkbox_ridgeStartTime,'Value'))
% 	dx = out_x(1,1) - linha(1,1);
% 	dy = out_y(1,1) - linha(1,2);
% 	out_x = out_x - dx;
% 	out_y = out_y - dy;
% end

if (get(handles.checkbox_ridgeStartTime,'Value') & t0 > 0)
    x1 = [linha(1,1) linha(2,1)];       y1 = [linha(1,2) linha(2,2)];
    x2 = [out_x(1,1) out_x(2,1)];       y2 = [out_y(1,1) out_y(2,1)];
    [p_lon,p_lat,omega] = calcBoninEulerPole(x1,y1,x2,y2);
    [rlon,rlat] = rot_euler(out_x(:),out_y(:),p_lon,p_lat,-omega);
    % Now we have to rebuild the "telhas" matrix format
    n_col = size(out_x,2); 
    out_x = reshape(rlon,4,n_col);
    out_y = reshape(rlat,4,n_col);
    clear rlon rlat;
end

% Save this for evental use in write_script
to_save.line = linha;
to_save.opt_E = opt_E;      to_save.opt_I = opt_I;
to_save.opt_T = opt_T;      to_save.opt_N = opt_N;

% Remove last column, it has only zeros
out_x(:,end) = [];              out_y(:,end) = [];
% Split between direct and inverse telhas
out_x1 = out_x(:,1:2:end);      out_x2 = out_x(:,2:2:end);
out_y1 = out_y(:,1:2:end);      out_y2 = out_y(:,2:2:end);
if (first_mag < 0)      % Negative magnetizations correspond to direct polarities in telha.h
    cor_d = 'r';    cor_r = 'b';
else
    cor_d = 'b';    cor_r = 'r';
end

hd = patch(out_x1,out_y1,cor_d,'EdgeColor','none','Tag','tapete');
hr = patch(out_x2,out_y2,cor_r,'EdgeColor','none','Tag','tapete_R');
set(hd,'UserData',to_save)

%set_telhas_uis(hd);     set_telhas_uis(hr)
draw_funs(hd,'telhas_patch')
draw_funs(hr,'telhas_patch')

% --------------------------------------------------------------------
function set_telhas_uis(h)
cmenuHand = uicontextmenu;
set(h, 'UIContextMenu', cmenuHand);
uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');

% --------------------------------------------------------------------
function pushbutton_callMagBarCode_Callback(hObject, eventdata, handles)
MagBarCode([handles.path_data 'Cande_Kent_95.dat'])

% --------------------------------------------------------------------
function pushbutton_polesList_Callback(hObject, eventdata, handles)
fid = fopen([handles.path_continent 'lista_polos.dat'],'rt');
c = fread(fid,'*char').';
fclose(fid);
s = strread(c,'%s','delimiter','\n');

[s,v] = choosebox('Name','One Euler list',...
                    'PromptString','List of poles:',...
                    'SelectString','Selected poles:',...
                    'ListSize',[380 300],...
                    'ListString',s);

if (v == 1)         % Finite pole
    return          % This case is not used here
elseif (v == 2)     % Stage poles
    set(handles.edit_polesFile,'String',s)
end

% -----------------------------------------------------------------------------------
function pushbutton_pickLine_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    % Test if we have potential target lines and their type
    h_mir_lines = findobj(handles.h_calling_fig,'Type','line');     % Fish all objects of type line in Mirone figure
    if (isempty(h_mir_lines))                                       % We don't have any lines
        str = ['If you hited this button on purpose, than you deserve the following insult.',...
                'You #!|"*!%!?~^)--$&.',... 
                'THERE ARE NO LINES IN THAT FIGURE.'];
        errordlg(str,'Chico Clever');   set(hObject,'Value',0);     return;
    end
    % The above test is not enough. For exemple, coastlines are not eligible neither,
    % but is very cumbersome to test all the possibilities of pure non-eligible lines.
    set(handles.h_calling_fig,'pointer','crosshair')
    h_line = get_polygon(handles.h_calling_fig);          % Get the line handle
    handles.h_line_orig = h_line;
    if (~isempty(h_line))
        % Create a empty line handle that will hold the rotated line
        handles.h_line = line('parent',get(handles.h_calling_fig,'CurrentAxes'),'XData',[],'YData',[], ...
            'LineStyle','-.','LineWidth',2);
        set(handles.h_active_line_str,'String','GOT A LINE TO WORK WITH','ForegroundColor',[0 1 0])
    else
        handles.h_line_orig = [];
        set(handles.h_active_line_str,'String','NO ACTIVE LINE','ForegroundColor',[1 0 0])
        set(hObject,'Value',0)
    end
    set(handles.h_calling_fig,'pointer','arrow')
    set(hObject,'Value',0)
    figure(handles.figure1)                 % Bring this figure to front again
else        % What should I do?
end
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')
    delete(handles.figure1);
end


% --- Creates and returns a handle to the GUI figure. 
function telhometro_LayoutFcn(h1,handles);

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','Telhometro',...
'NumberTitle','off',...
'Position',[520 600 463 200],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',[]);

h2 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@telhometro_uicallback,h1,'edit_polesFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[10 105 211 21],...
'Style','edit',...
'Tag','edit_polesFile',...
'UserData','DoRotations');

h3 = uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'pushbutton_readPolesFile_Callback'},...
'FontSize',10,...
'FontWeight','bold',...
'Position',[220 104 21 21],...
'String','...',...
'Tag','pushbutton_readPolesFile',...
'UserData','DoRotations');

h4 = uicontrol('Parent',h1,...
'Position',[10 129 51 15],...
'String','Poles file',...
'Style','text',...
'Tag','text2');

h5 = uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'pushbutton_Cancel_Callback'},...
'Position',[300 13 66 23],...
'String','Cancel',...
'Tag','pushbutton_Cancel');

h6 = uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'pushbutton_compute_Callback'},...
'Position',[389 13 66 23],...
'String','Compute',...
'Tag','pushbutton_compute',...
'UserData','DoRotations');

h7 = uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'checkbox_revertRot_Callback'},...
'Position',[10 81 141 15],...
'String','Revert sense of rotation',...
'Style','checkbox',...
'TooltipString','Revert the sense of rotation defined by the stages poles',...
'Tag','checkbox_revertRot',...
'UserData','DoRotations');

h8 = uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'pushbutton_callMagBarCode_Callback'},...
'Position',[300 63 131 23],...
'String','Magnetic Bar Code',...
'TooltipString','Open the magnetic bar code window',...
'Tag','pushbutton_callMagBarCode');

h9 = uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'pushbutton_polesList_Callback'},...
'Position',[300 103 131 23],...
'String','Poles selector',...
'TooltipString','Construct a stage poles file',...
'Tag','pushbutton_polesList');

h10 = uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'pushbutton_pickLine_Callback'},...
'Position',[10 158 161 23],...
'String','Pick line from Figure',...
'TooltipString','Allows you to mouse select one line from a Mirone figure',...
'Tag','pushbutton_pickLine',...
'UserData','DoRotations');

h11 = uicontrol('Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[1 0 0],...
'Position',[210 161 231 16],...
'String','NO ACTIVE LINE',...
'Style','text',...
'Tag','text_activeLine');

h12 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@telhometro_uicallback,h1,'edit_timeStart_Callback'},...
'Position',[20 25 41 21],...
'String','0',...
'Style','edit',...
'Tag','edit_timeStart',...
'UserData','DoRotations');

h13 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@telhometro_uicallback,h1,'edit_timeEnd_Callback'},...
'Position',[90 25 41 21],...
'Style','edit',...
'Tag','edit_timeEnd',...
'UserData','DoRotations');

h14 = uicontrol('Parent',h1,...
'Position',[17 48 51 15],...
'String','Start time',...
'Style','text',...
'Tag','text19',...
'UserData','DoRotations');

h15 = uicontrol('Parent',h1,...
'Position',[86 48 51 15],...
'String','Stop time',...
'Style','text',...
'Tag','text20',...
'UserData','DoRotations');

h16 = uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'checkbox_ridgeStartTime_Callback'},...
'Position',[140 25 111 15],...
'String','Ridge at Start time',...
'Style','checkbox',...
'TooltipString','Use this when you want to and old (Start time) ridge position',...
'Tag','checkbox_ridgeStartTime');

function telhometro_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

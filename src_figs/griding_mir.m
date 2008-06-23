function varargout = griding_mir(varargin)
% M-File changed by desGUIDE
% varargin   command line arguments to griding_mir (see VARARGIN) 

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
griding_mir_LayoutFcn(hObject);
handles = guihandles(hObject);

movegui(hObject,'center');                      % Reposition the window on screen
global home_dir

if (isempty(home_dir)),	handles.home_dir = pwd;        % Case when this function was called directly
else					handles.home_dir = home_dir;
end

if isempty(home_dir),	f_path = [pwd filesep 'data' filesep];
else					f_path = [home_dir filesep 'data' filesep];
end

% Import icons
load([f_path 'mirone_icons.mat'],'Mfopen_ico');
set(handles.pushbutton_InputFile,'CData',Mfopen_ico)
clear Mfopen_ico;

handles.command = cell(50,1);
handles.x_min = [];             handles.x_max = [];
handles.y_min = [];             handles.y_max = [];
handles.x_inc = [];             handles.y_inc = [];
handles.dms_xinc = 0;           handles.dms_yinc = 0;
handles.one_or_zero = 1;        % For Grid Registration grids, which are the most common cases
handles.h_calling_fig = [];     % Update this bellow when integrated in Mirone

% Inactivate the headers parameters. They will be activated by the header checkbox
set(handles.edit_nHeaders,'Enable','inactive')
set(handles.popup_binInput,'Visible','off')
set(handles.edit_binary_ncolumnIn,'Visible','off')  % Those two have to wait until we
set(handles.pushbutton_Help_H,'Visible','off')      % know how to read binary files

% When called by Mirone varargin must contain: mirone fig handle, "type"
if ~isempty(varargin)
    if ( length(varargin) == 2 && ishandle(varargin{1}) && ischar(varargin{2}) )
        handles.h_calling_fig = varargin{1};
        type = varargin{2};
    else
        type = 'surface';   % Default to surface
    end
else
    type = 'surface';   % Default to surface
end
handles.type = type;

% Choose the default griding_mir_export method
% In Mirone the 'Delauny Triangulation' method is not yet implemented
if strcmp(type,'surface')
    set(hObject,'Name','Surface')
    %set(handles.popup_GridMethod, 'String', {'Minimum Curvature';'Delauny Triangulation';'Near Neighbor'});
    set(handles.popup_GridMethod, 'String', {'Minimum Curvature';'Near Neighbor'});
    handles.command{1} = 'surface ';
    set(handles.edit_S1_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
    set(handles.popup_S2_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
    set(handles.pushbutton_Help_S,'Enable', 'off')
    set(handles.checkbox_Option_F,'Enable', 'off')
elseif strcmp(type,'triangulate')
    set(hObject,'Name','Triangulate')
    set(handles.popup_GridMethod, 'String', {'Delauny Triangulation';'Minimum Curvature';'Near Neighbor'});
    handles.command{1} = 'triangulate ';
    set(handles.edit_S1_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
    set(handles.popup_S2_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
    set(handles.pushbutton_Help_S,'Enable', 'off')
elseif strcmp(type,'nearneighbor')
    set(hObject,'Name','Nearneighbor')
    %set(handles.popup_GridMethod, 'String', {'Near Neighbor';'Delauny Triangulation';'Minimum Curvature'});
    set(handles.popup_GridMethod, 'String', {'Near Neighbor';'Minimum Curvature'});
    set(handles.checkbox_Option_V,'Enable', 'off')
    handles.command{1} = 'nearneighbor ';
else        % Defaults to surface
    set(hObject,'Name','Surface')
    %set(handles.popup_GridMethod, 'String', {'Minimum Curvature';'Delauny Triangulation';'Near Neighbor'});
    set(handles.popup_GridMethod, 'String', {'Minimum Curvature';'Near Neighbor'});
    handles.command{1} = 'surface ';
    set(handles.edit_S1_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
    set(handles.popup_S2_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
    set(handles.checkbox_Option_F,'Enable', 'off')
    set(handles.pushbutton_Help_S,'Enable', 'off')
end

if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
    cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
    handles.last_dir = cfig_handles.last_dir;
    handles.home_dir = cfig_handles.home_dir;
else
    handles.last_dir = pwd;
    handles.home_dir = pwd;
end

%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
bgcolor = get(0,'DefaultUicontrolBackgroundColor');
framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
set(0,'Units','pixels');    set(hObject,'Units','pixels')    % Pixels are easier to reason with
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
    bgc = get (h_t(i),'BackgroundColor');   fgc = get (h_t(i),'ForegroundColor');tag=get(h_t(i),'Tag');
    uicontrol('Parent',hObject, 'Style','text', 'Position',t_size,'String',t_str, ...
        'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw,'UserData',usr_d,'Tag',tag);
end
delete(h_t)
%------------- END Pro look (3D) -------------------------------------------------------

% Choose default command line output for griding_mir_export
handles.output = hObject;
guidata(hObject, handles);

set(hObject,'Visible','on');
% UIWAIT makes griding_mir_export wait for user response (see UIRESUME)
uiwait(handles.figure1);

handles = guidata(hObject);
out = griding_mir_OutputFcn(hObject, [], handles);
varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = griding_mir_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
try,    varargout{1} = handles.output;  end
% The figure can be deleted now
delete(handles.figure1);

% -----------------------------------------------------------------------------------
function checkbox_Option_H_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    if isempty(get(handles.edit_InputFile,'String'))
        handles.command{41} = ' -H';
        set(handles.edit_nHeaders,'Enable','on', 'Backgroundcolor','white')
        %set(handles.popup_binInput,'Enable','on')
    else
        errordlg(['Now it''s to late to chose this option. If you realy want to do this, clean ' ...
                'the "Input Data File" box, hit the "Return" key, and start again but ' ...
                'check this option before load the data file.'],'Error')
        handles.command{41} = '';
        set(hObject,'Value',0)
        set(handles.edit_nHeaders,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
        %set(handles.popup_binInput,'Enable','off')
    end
else
    handles.command{41} = '';     handles.command{42} = '';
    set(handles.edit_nHeaders,'String','1','Enable','off', 'Backgroundcolor',[.764,.603,.603])
    %set(handles.popup_binInput,'Enable','off')
end
guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function edit_nHeaders_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if isnan(str2double(xx))
    errordlg([xx ' is not a valid number'],'Error');    return
end
if ~isempty(xx) && str2double(xx) ~= 0
    handles.command{42} = xx;
elseif str2double(xx) == 0
    handles.command{42} = '';
    set(hObject,'String','1')
else
    handles.command{42} = '';
    set(hObject,'String','1')
end
guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function popup_binInput_Callback(hObject, eventdata, handles)
if ~isempty(get(handles.edit_InputFile,'String'))
    errordlg(['Now it''s to late to choose any of these options. If you realy want to do this, clean ' ...
              'the "Input Data File" box, hit the "Return" key, and start again but first ' ...
              'choose here the input file format before load the data file.'],'Error')
    return
end
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case 'ascii'
        handles.command{38} = '';     handles.command{39} = '';
        set(handles.edit_binary_ncolumnIn,'String','3')
        set(handles.edit_binary_ncolumnIn,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
        set(handles.checkbox_Option_H,'Enable','on')
        set(handles.edit_nHeaders,'Enable','on', 'Backgroundcolor','white')
        set(handles.checkbox_Option_H,'Value',0);   set(handles.edit_nHeaders,'String','1');
        set(handles.edit_binary_ncolumnIn,'Visible','off')
        if ~get(handles.checkbox_Option_H,'Value')
            set(handles.edit_nHeaders,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
            set(hObject,'Enable','off')
        end
    case 'binary double'
        set(handles.edit_binary_ncolumnIn,'Visible','on')
        handles.command{38} = ' -bi';  handles.command{41} = '';    handles.command{42} = '';
        set(handles.edit_binary_ncolumnIn,'Enable','on', 'Backgroundcolor','white')
        set(handles.checkbox_Option_H,'Value',0);   set(handles.edit_nHeaders,'String','1');
        set(handles.checkbox_Option_H,'Enable','off')
        set(handles.edit_nHeaders,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
    case 'binary single'
        set(handles.edit_binary_ncolumnIn,'Visible','on')
        handles.command{38} = ' -bis';  handles.command{41} = '';    handles.command{42} = '';
        %set(handles.edit_binary_ncolumnIn,'Enable','on', 'Backgroundcolor','white')
        set(handles.checkbox_Option_H,'Value',0);   set(handles.edit_nHeaders,'String','1');
        set(handles.checkbox_Option_H,'Enable','off')
        set(handles.edit_nHeaders,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
end
guidata(hObject,handles)

%----------------------------------------------------------------------------------------------
function edit_binary_ncolumnIn_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if ~isempty(xx) && str2double(xx) > 3
    handles.command{39} = xx;
    set(hObject,'String',num2str(fix(str2double(xx))))  % in case user tries to joke
else
    handles.command{39} = '';
    set(hObject,'String','3')
end
guidata(hObject,handles)

%----------------------------------------------------------------------------------------------
function pushbutton_Help_H_Callback(hObject, eventdata, handles)
message = {'Input file has Header record(s). Number of header records can be changed'
           'on the "N? of headers" box. Not used with binary data. This, like all GMT'
           'programs, accept table data input in ASCII or binary data. When using'
           'binary data you must be aware of the fact that GMT has no way of'
           'determining the actual number of columns in the file. You should therefore'
           'pass that information to GMT via the ascii/binary popup menu and its side'
           'box, where the actual number of data columns can be introduced.'};
helpdlg(message,'Help on input data');

%----------------------------------------------------------------------------------------------
function edit_InputFile_Callback(hObject, eventdata, handles, opt)
if (nargin == 4 && ischar(opt))   % OPT is a file name transmited by pushbutton_InputFile_Callback
    xx = opt;
    hObject = handles.edit_InputFile;    % hObject contained the handle to pushbutton_InputFile_Callback
else
    xx = get(hObject,'String');
end
if ~isempty(xx)
    handles.command{3} = [xx];
    str = ['minmax -C ' xx];
    [s,w] = mat_lyies(str);
    if ~(isequal(s,0))                  % An error as occured
        errordlg(w,'GMT Error');        return
    else
        val = string_token(w);      % Decompose the string output from minmax
        if length(val) < 6
            errordlg('File error. Your file doesn''t have at least 3 columns','Error')
            handles.command{3} = '';        set(hObject,'String','')
            guidata(hObject, handles);
            return
        end
        handles.command{6} = val{1};    handles.command{7} = '/';
        handles.command{8} = val{2};    handles.command{9} = '/';
        handles.command{10} = val{3};   handles.command{11} = '/';
        handles.command{12} = val{4};   handles.command{5} = ' -R';
        % compute also handles.x_min, handles.x_max, ...
        xx = test_dms(val{1});      x_min = 0;
        if str2double(xx{1}) > 0
            for i = 1:length(xx)   x_min = x_min + str2double(xx{i}) / (60^(i-1));    end
        else
            for i = 1:length(xx)   x_min = x_min - abs(str2double(xx{i})) / (60^(i-1));   end
        end
        xx = test_dms(val{2});      x_max = 0;
        if str2double(xx{1}) > 0
            for i = 1:length(xx)   x_max = x_max + str2double(xx{i}) / (60^(i-1));    end
        else
            for i = 1:length(xx)   x_max = x_max - abs(str2double(xx{i})) / (60^(i-1));   end
        end
        xx = test_dms(val{3});      y_min = 0;
        if str2double(xx{1}) > 0
            for i = 1:length(xx)   y_min = y_min + str2double(xx{i}) / (60^(i-1));    end
        else
            for i = 1:length(xx)   y_min = y_min - abs(str2double(xx{i})) / (60^(i-1));   end
        end
        xx = test_dms(val{4});      y_max = 0;
        if str2double(xx{1}) > 0
            for i = 1:length(xx)   y_max = y_max + str2double(xx{i}) / (60^(i-1));    end
        else
            for i = 1:length(xx)   y_max = y_max - abs(str2double(xx{i})) / (60^(i-1));   end
        end
        handles.x_min = x_min;  handles.x_max = x_max;  handles.y_min = y_min;  handles.y_max = y_max;
        set(handles.edit_Xmin,'String',val{1});     set(handles.edit_Xmax,'String',val{2});
        set(handles.edit_Ymin,'String',val{3});     set(handles.edit_Ymax,'String',val{4});
        % Until something more inteligent is devised (like using some kind of estatistics to estimate
        % default's Nrow & Ncol) the default value of Nrow = Ncol = 100 will be used.
        x_inc = ivan_the_terrible((x_max - x_min),100,1);
        y_inc = ivan_the_terrible((y_max - y_min),100,1);
        handles.command{14} = ' -I';    handles.command{16} = '/';
        handles.command{15} = num2str(x_inc,8);   handles.command{17} = num2str(y_inc,8);
        set(handles.edit_Xinc,'String',num2str(x_inc,8));     set(handles.edit_Yinc,'String',num2str(y_inc,8));
        set(handles.edit_Ncols,'String','100');     set(handles.edit_Nrows,'String','100');
        guidata(hObject, handles);
    end
else
    % Reset everything to initial state (falta a parte do nearneigh)
    set(handles.edit_Xmin,'String','');     set(handles.edit_Xmax,'String','');
    set(handles.edit_Ymin,'String','');     set(handles.edit_Ymax,'String','');
    set(handles.edit_Xinc,'String','');     set(handles.edit_Yinc,'String','');
    set(handles.edit_Ncols,'String','');    set(handles.edit_Nrows,'String','');
    set(handles.checkbox_ToggleXY,'Value',0);
    %set(handles.popup_binInput,'Value',1);
    %set(handles.edit_binary_ncolumnIn,'String','3');
    for i = 3:length(handles.command)   handles.command{i} = '';  end
    guidata(hObject,handles)
end

%----------------------------------------------------------------------------------------------
function pushbutton_InputFile_Callback(hObject, eventdata, handles)
cd(handles.last_dir)
[FileName,PathName] = uigetfile({'*.dat;*.DAT;*.xyz;*.XYX', 'XYZ files (*.dat,*.DAT,*.xyz,*.XYZ)';'*.*', 'All Files (*.*)'},'Select input data');
pause(0.01);
cd(handles.home_dir);       % allways go home to avoid troubles
if isequal(FileName,0);     return;     end
handles.last_dir = PathName;
handles.command{3} = [PathName FileName];
set(handles.edit_InputFile, 'String',handles.command{3})
guidata(hObject, handles);
edit_InputFile_Callback(hObject, eventdata, handles,[PathName FileName]);

% -----------------------------------------------------------------------------------
function edit_Xmin_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)            % when dd:mm or dd:mm:ss was given
    x_min = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   x_min = x_min + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   x_min = x_min - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.x_min = x_min;
    if ~isempty(handles.x_max) && x_min >= handles.x_max
        errordlg('West Longitude >= East Longitude ','Error in Longitude limits')
        handles.command{6} = '';
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        handles.command{6} = xx;
        handles.command{5} = ' -R';handles.command{7} = '/';handles.command{9} = '/';
        handles.command{11} = '/';
        guidata(hObject, handles);
    end
    nc = get(handles.edit_Ncols,'String');
    if ~isempty(handles.x_max) && ~isempty(nc)       % x_max and ncols boxes are filled
        % Compute Ncols, but first must recompute x_inc
        x_inc = ivan_the_terrible((handles.x_max - x_min),round(abs(str2double(nc))),1);
        xx = floor((handles.x_max - str2double(xx)) / (str2double(get(handles.edit_Xinc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_Xinc,'String',num2str(x_inc,8))
        handles.command{14} = [' -I'];      handles.command{15} = [num2str(x_inc,8)];
        guidata(hObject, handles);
    elseif ~isempty(handles.x_max)      % x_max box is filled but ncol is not, so put to the default (100)
        x_inc = ivan_the_terrible((handles.x_max - x_min),100,1);
        set(handles.edit_Xinc,'String',num2str(x_inc,8))
        set(handles.edit_Ncols,'String','100')
        handles.command{14} = [' -I'];      handles.command{15} = [num2str(x_inc,8)];
        guidata(hObject, handles);
    end
else                % box is empty, so clear also x_inc and ncols
    handles.command{6} = '';      handles.command{15} = '';
    set(handles.edit_Xinc,'String','');     set(handles.edit_Ncols,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% -----------------------------------------------------------------------------------
function edit_Xmax_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    x_max = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   x_max = x_max + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   x_max = x_max - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.x_max = x_max;
    if ~isempty(handles.x_min) && x_max <= handles.x_min 
        errordlg('East Longitude <= West Longitude','Error in Longitude limits')
        handles.command{8} = '';
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        handles.command{8} = xx;
        handles.command{5} = ' -R';	handles.command{7} = '/';handles.command{9} = '/';
        handles.command{11} = '/';
        guidata(hObject, handles);
    end
    nc = get(handles.edit_Ncols,'String');
    if ~isempty(handles.x_min) && ~isempty(nc)       % x_max and ncols boxes are filled
        % Compute Ncols, but first must recompute x_inc
        x_inc = ivan_the_terrible((x_max - handles.x_min),round(abs(str2double(nc))),1);
        xx = floor((handles.x_min - str2double(xx)) / (str2double(get(handles.edit_Xinc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_Xinc,'String',num2str(x_inc,8))
        handles.command{14} = ' -I';      handles.command{15} = [num2str(x_inc,8)];
        guidata(hObject, handles);    
    elseif ~isempty(handles.x_min)      % x_min box is filled but ncol is not, so put to the default (100)
        x_inc = ivan_the_terrible((x_max - handles.x_min),100,1);
        set(handles.edit_Xinc,'String',num2str(x_inc,8))
        set(handles.edit_Ncols,'String','100')
        handles.command{14} = ' -I';      handles.command{15} = num2str(x_inc,8);
        guidata(hObject, handles);
    end
else                % box is empty, so clear also x_inc and ncols
    handles.command{8} = '';      handles.command{15} = '';
    set(handles.edit_Xinc,'String','');     set(handles.edit_Ncols,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% -----------------------------------------------------------------------------------
function edit_Xinc_Callback(hObject, eventdata, handles)
dms = 0;
xx = get(hObject,'String');     val = test_dms(xx);
if isempty(val)
    handles.command{15} = '';    set(hObject, 'String', '');
    guidata(hObject, handles);    return
end
% If it survived then ...
if length(val) > 1    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
x_inc = 0;
for i = 1:length(val)   x_inc = x_inc + str2double(val{i}) / (60^(i-1));    end
if ~isempty(handles.x_min) & ~isempty(handles.x_max)
    % Make whatever x_inc given compatible with GMT_grd_RI_verify
    x_inc = ivan_the_terrible((handles.x_max - handles.x_min), x_inc,2);
    if ~dms         % case of decimal unities
        set(hObject,'String',num2str(x_inc,8))
        ncol = floor((handles.x_max - handles.x_min) / x_inc + 0.5) + handles.one_or_zero;
    else            % inc was in dd:mm or dd:mm:ss format
        ncol = floor((handles.x_max - handles.x_min) / x_inc + 0.5) + handles.one_or_zero;
        ddmm = dec2deg(x_inc);
        set(hObject,'String',ddmm)
    end
    set(handles.edit_Ncols,'String',num2str(ncol))
end
handles.dms_xinc = dms;
handles.command{14} = [' -I'];      handles.command{15} = [num2str(x_inc,8)];
guidata(hObject, handles);
if isempty(get(handles.edit_Yinc,'String'))     set(handles.edit_Yinc,'String',xx);    end

% -----------------------------------------------------------------------------------
function edit_Ncols_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if ~isempty(get(handles.edit_Xmin,'String')) & ~isempty(get(handles.edit_Xmax,'String')) & ...
        ~isempty(get(handles.edit_Xinc,'String')) & ~isempty(xx)
    x_inc = ivan_the_terrible((handles.x_max - handles.x_min),round(abs(str2double(xx))),1);
    if handles.dms_xinc        % x_inc was given in dd:mm:ss format
        ddmm = dec2deg(x_inc);
        set(handles.edit_Xinc,'String',ddmm)
    else                    % x_inc was given in decimal format
        set(handles.edit_Xinc,'String',num2str(x_inc,8));
    end
    handles.command{15} = num2str(x_inc,8);
    guidata(hObject, handles);
end

% -----------------------------------------------------------------------------------
function edit_Ymin_Callback(hObject, eventdata, handles)
% Read value either in decimal or in the dd:mm or dd_mm:ss formats and do some tests
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    y_min = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   y_min = y_min + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   y_min = y_min - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.y_min = y_min;
    if ~isempty(handles.y_max) && y_min >= handles.y_max
        errordlg('South Latitude >= North Latitude','Error in Latitude limits')
        handles.command{10} = '';
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        handles.command{10} = xx;
        handles.command{5} = ' -R'	;handles.command{7} = '/';	handles.command{9} = '/';
        handles.command{11} = '/';
        guidata(hObject, handles);
    end
    nr = get(handles.edit_Nrows,'String');
    if ~isempty(handles.y_max) && ~isempty(nr)       % y_max and nrows boxes are filled
        % Compute Nrowss, but first must recompute y_inc
        y_inc = ivan_the_terrible((handles.y_max - y_min),round(abs(str2double(nr))),1);
        xx = floor((handles.y_max - str2double(xx)) / (str2double(get(handles.edit_Yinc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_Yinc,'String',num2str(y_inc,8))
        handles.command{14} = ' -I';      handles.command{16} = '/';
        handles.command{17} = num2str(y_inc,8);
        guidata(hObject, handles);
    elseif ~isempty(handles.y_max)      % y_max box is filled but nrows is not, so put to the default (100)
        y_inc = ivan_the_terrible((handles.y_max - y_min),100,1);
        set(handles.edit_Yinc,'String',num2str(y_inc,8))
        set(handles.edit_Nrows,'String','100')
        handles.command{14} = ' -I';      handles.command{16} = '/';
        handles.command{17} = num2str(y_inc,8);
        guidata(hObject, handles);
    end
else                % box is empty, so clear also y_inc and nrows
    handles.command{10} = '';      handles.command{16} = '';      handles.command{17} = '';
    set(handles.edit_Yinc,'String','');     set(handles.edit_Nrows,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% -----------------------------------------------------------------------------------
function edit_Ymax_Callback(hObject, eventdata, handles)
% Read value either in decimal or in the dd:mm or dd_mm:ss formats and do some tests
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    y_max = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val),	y_max = y_max + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val),	y_max = y_max - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.y_max = y_max;
    if ~isempty(handles.y_min) && y_max <= handles.y_min 
        errordlg('North Latitude <= South Latitude','Error in Latitude limits')
        handles.command{12} = '';
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        handles.command{12} = xx;
        handles.command{5} = [' -R'];handles.command{7} = ['/'];handles.command{9} = ['/'];
        handles.command{11} = '/';
        guidata(hObject, handles);
    end
    nr = get(handles.edit_Nrows,'String');
    if ~isempty(handles.y_min) & ~isempty(nr)       % y_min and nrows boxes are filled
        % Compute Nrows, but first must recompute y_inc
        y_inc = ivan_the_terrible((y_max - handles.y_min),round(abs(str2double(nr))),1);
        xx = floor((handles.y_min - str2double(xx)) / (str2double(get(handles.edit_Yinc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_Yinc,'String',num2str(y_inc,8))
        handles.command{14} = [' -I'];      handles.command{16} = '/';
        handles.command{17} = [num2str(y_inc,8)];
        guidata(hObject, handles);
    elseif ~isempty(handles.y_min)      % y_min box is filled but nrows is not, so put to the default (100)
        y_inc = ivan_the_terrible((y_max - handles.y_min),100,1);
        set(handles.edit_Yinc,'String',num2str(y_inc,8))
        set(handles.edit_Nrows,'String','100')
        handles.command{14} = [' -I'];      handles.command{16} = '/';
        handles.command{17} = [num2str(y_inc,8)];
        guidata(hObject, handles);
    end
else                % This box is empty, so clear also y_inc and nrows
    handles.command{10} = '';      handles.command{16} = '';      handles.command{17} = '';
    set(handles.edit_Yinc,'String','');     set(handles.edit_Nrows,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% -----------------------------------------------------------------------------------
function edit_Yinc_Callback(hObject, eventdata, handles)
dms = 0;
xx = get(hObject,'String');     val = test_dms(xx);
if isempty(val)
    handles.command{16} = '';  handles.command{17} = '';    set(hObject, 'String', '');
    guidata(hObject, handles);    return
end
% If it survived then ...
if length(val) > 1    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
y_inc = 0;
for i = 1:length(val)   y_inc = y_inc + str2double(val{i}) / (60^(i-1));    end
if ~isempty(handles.y_min) & ~isempty(handles.y_max)
    % Make whatever y_inc given compatible with GMT_grd_RI_verify
    y_inc = ivan_the_terrible((handles.y_max - handles.y_min), y_inc,2);
    if ~dms         % case of decimal unities
        set(hObject,'String',num2str(y_inc,8))
        nrow = floor((handles.y_max - handles.y_min) / y_inc + 0.5) + handles.one_or_zero;
    else            % inc was in dd:mm or dd:mm:ss format
        nrow = floor((handles.y_max - handles.y_min) / y_inc + 0.5) + handles.one_or_zero;
        ddmm = dec2deg(y_inc);
        set(hObject,'String',ddmm)
    end
    set(handles.edit_Nrows,'String',num2str(nrow))
end
handles.dms_yinc = dms;
handles.command{14} = [' -I'];  handles.command{16} = '/';   handles.command{17} = [num2str(y_inc,8)];
guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function edit_Nrows_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if ~isempty(get(handles.edit_Ymin,'String')) & ~isempty(get(handles.edit_Ymax,'String')) & ...
        ~isempty(get(handles.edit_Yinc,'String')) & ~isempty(xx)
    y_inc = ivan_the_terrible((handles.y_max - handles.y_min),round(abs(str2double(xx))),1);
    if handles.dms_yinc        % y_inc was given in dd:mm:ss format
        ddmm = dec2deg(y_inc);
        set(handles.edit_Yinc,'String',ddmm)
    else                    % y_inc was given in decimal format
        set(handles.edit_Yinc,'String',num2str(y_inc,8));
    end
    handles.command{14} = ' -I';     handles.command{16} = '/';  handles.command{17} = num2str(y_inc,8);
    guidata(hObject, handles);
end

% -----------------------------------------------------------------------------------
function checkbox_ToggleXY_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    t_xmin = get(handles.edit_Xmin,'String');   t_xmax = get(handles.edit_Xmax,'String');
    t_ymin = get(handles.edit_Ymin,'String');   t_ymax = get(handles.edit_Ymax,'String');
    if ~isempty(t_xmin) && ~isempty(t_xmax) && ~isempty(t_ymin) && ~isempty(t_ymax)
        set(handles.edit_Xmin,'String',t_ymin);        set(handles.edit_Xmax,'String',t_ymax)
        set(handles.edit_Ymin,'String',t_xmin);        set(handles.edit_Ymax,'String',t_xmax)
        handles.command{6} = t_ymin;       handles.command{8} = t_ymax;
        handles.command{10} = t_xmin;      handles.command{12} = t_xmax;
        handles.command{36} = ' -:';
    else        % There is nothing yet to toggle
        set(hObject,'Value',0)
    end
    t_xinc = get(handles.edit_Xinc,'String');    t_yinc = get(handles.edit_Yinc,'String');
    if ~isempty(t_xinc) && ~isempty(t_yinc)
        set(handles.edit_Xinc,'String',t_yinc);        set(handles.edit_Yinc,'String',t_xinc)
        handles.command{15} = t_yinc;       handles.command{17} = t_xinc;
    end
    t_ncol = get(handles.edit_Ncols,'String');    t_nrow = get(handles.edit_Nrows,'String');
    if ~isempty(t_ncol) && ~isempty(t_nrow)
        set(handles.edit_Ncols,'String',t_nrow);        set(handles.edit_Nrows,'String',t_ncol)
    end
    guidata(hObject,handles)
else
    t_xmin = get(handles.edit_Xmin,'String');   t_xmax = get(handles.edit_Xmax,'String');
    t_ymin = get(handles.edit_Ymin,'String');   t_ymax = get(handles.edit_Ymax,'String');
    if ~isempty(t_xmin) && ~isempty(t_xmax) && ~isempty(t_ymin) && ~isempty(t_ymax)
        set(handles.edit_Xmin,'String',t_ymin);        set(handles.edit_Xmax,'String',t_ymax)
        set(handles.edit_Ymin,'String',t_xmin);        set(handles.edit_Ymax,'String',t_xmax)
        handles.command{6} = t_ymin;       handles.command{8} = t_ymax;
        handles.command{10} = t_xmin;      handles.command{12} = t_xmax;
        handles.command{36} = '';
    end
    t_xinc = get(handles.edit_Xinc,'String');    t_yinc = get(handles.edit_Yinc,'String');
    if ~isempty(t_xinc) && ~isempty(t_yinc)
        set(handles.edit_Xinc,'String',t_yinc);        set(handles.edit_Yinc,'String',t_xinc)
        handles.command{15} = t_yinc;       handles.command{17} = t_xinc;
    end
    t_ncol = get(handles.edit_Ncols,'String');    t_nrow = get(handles.edit_Nrows,'String');
    if ~isempty(t_ncol) & ~isempty(t_nrow)
        set(handles.edit_Ncols,'String',t_nrow);        set(handles.edit_Nrows,'String',t_ncol)
    end
    guidata(hObject,handles)
end

% -----------------------------------------------------------------------------------
function checkbox_Option_F_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.one_or_zero = 0;    handles.command{44} = [' -F'];
else
    handles.one_or_zero = 1;    handles.command{44} = '';
end
guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function pushbutton_Help_R_F_toggle_Callback(hObject, eventdata, handles)
message = {'Min and Max, of "X Direction" and "Y Direction" specify the Region of'
    'interest. To specify boundaries in degrees and minutes [and seconds],'
    'use the dd:mm[:ss.xx] format.'
    '"Spacing" sets the grid size for grid output. You may choose different'
    'spacings for X and Y. Also here you can use the dd:mm[:ss.xx] format.'
    'In "#of lines" it is offered the easyeast way of controling the grid'
    'dimensions (lines & columns).'
    '"Toggle x,y" toggles between (longitude,latitude) and (latitude,longitude)'
    'input/output. Default is (longitude,latitude).'
    '"Pixel registration" forces the grid to be in pixel registration. Default'
    'is grid registration. Read Appendix B of GMT Cookbook to learn the'
    'differences between grid and pixel registrations. This option is only'
    'available with "Nearneighbor" and "Delauny Triangulation" interpolators.'};
helpdlg(message,'Help on Grid Line Geometry');

% -----------------------------------------------------------------------------------
function edit_S1_Neighbor_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if isnan(str2double(xx)) && ~isempty(xx)
    set(hObject, 'String', '');
    errordlg('"Search radius" must be a number','Error');
    return
end
if ~isempty(xx)
    handles.command{19} = ' -S';      handles.command{20} = [num2str(abs(str2double(xx)))];
    guidata(hObject, handles);
else
    handles.command{19} = '';      handles.command{20} = [xx];
    guidata(hObject, handles);
end

% -----------------------------------------------------------------------------------
function popup_S2_Neighbor_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case 'minutes'
        handles.command{21} = 'm';
    case 'seconds'
        handles.command{21} = 'c';
    case 'kilometers'
        handles.command{21} = 'k';
    case 'Kilometers'
        handles.command{21} = 'K';
    otherwise
        handles.command{21} = '';
end
guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function pushbutton_Help_S_Callback(hObject, eventdata, handles)
message = {'Sets the search radius in same units as the grid spacing. Use the'
           'optional popup to select in which unities the distance is given.'
           'Selecting kilometers implies that grid limits and spacing are in degrees.'
           'Use uppercase Kilometers if distances should be calculated using great'
           'circles [kilometrs uses flat Earth].'};
helpdlg(message,'Help -S option');

% -----------------------------------------------------------------------------------
function popup_GridMethod_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case 'Minimum Curvature'
        set(handles.edit_S1_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
        set(handles.popup_S2_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
        set(handles.pushbutton_Help_S,'Enable', 'off')
        set(handles.checkbox_Option_F,'Enable', 'off')
        handles.command{19} = '';     handles.command{20} = '';
        handles.command{21} = '';     handles.command{30} = '';
        handles.command{1} = ['surface '];
        handles.type = 'surface';
    case 'Delauny Triangulation'
        set(handles.edit_S1_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
        set(handles.popup_S2_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
        set(handles.pushbutton_Help_S,'Enable', 'off')
        set(handles.checkbox_Option_F,'Enable', 'on')
        handles.command{19} = '';     handles.command{20} = '';
        handles.command{21} = '';     handles.command{30} = '';
        handles.command{1} = ['triangulate '];
        handles.type = 'triangulate';
    case 'Near Neighbor'
        set(handles.edit_S1_Neighbor,'Enable', 'on', 'Backgroundcolor','white')
        set(handles.popup_S2_Neighbor,'Enable', 'on', 'Backgroundcolor','white')
        set(handles.pushbutton_Help_S,'Enable', 'on')
        set(handles.checkbox_Option_F,'Enable', 'on')
        set(handles.edit_S1_Neighbor,'String','')
        set(handles.popup_S2_Neighbor,'Value',1)
        handles.command{1} = ['nearneighbor '];     handles.command{30} = '';
        handles.type = 'nearneighbor';
end
guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function pushbutton_Grid_Options_Callback(hObject, eventdata, handles)
switch handles.command{1}
    case 'surface '
        handles.command{30} = [' ' surface_options(handles.command{30})];
    case 'triangulate '
        msgbox('Not yet programed')
    case 'nearneighbor '
        handles.command{30} = [' ' nearneighbor_options(handles.command{30})];
end
guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function checkbox_Option_V_Callback(hObject, eventdata, handles)
% The OK button will check this button state

% -----------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
error = VerifyCommand(handles);
if error == 0
    tmp = horzcat(handles.command{1:end});
    [tok,rem] = strtok(tmp);
    out{1} = tok;
    i = 2;
    while (rem)
        [tok,rem] = strtok(rem);
        out{i} = tok;        i = i + 1;
    end
    
    if (get(handles.checkbox_Option_V,'Value'))
        out{end+1} = '-V';
    end
    
    % Call Mirone that will either draw or save the newly computed grid
    set(handles.figure1,'Name','COMPUTING')
    switch handles.type
        case 'surface'
            [Z,head] = surface_m(out{2:end});     % We don't want the last ','
            tit = 'surface interpolation';          name = tit;
            set(handles.figure1,'Name','Surface')
        case 'nearneighbor'
            [Z,head] = nearneighbor_m(out{2:end}); % We don't want the last ','
            tit = 'nearneighbor interpolation';     name = tit;
            set(handles.figure1,'Name','Nearneighbor')
    end
    [ny,nx] = size(Z);
    X = linspace(head(1),head(2),nx);       Y = linspace(head(3),head(4),ny);
    tmp.head = [head(1) head(2) head(3) head(4) head(5) head(6) 0 head(8) head(9)];
    tmp.X = X;    tmp.Y = Y;    tmp.name = 'surface interpolation';
    new_window = mirone(Z,tmp);
end

% -----------------------------------------------------------------------------------
function error = VerifyCommand(handles)
% ERROR TESTING
error = 0;
if isempty(handles.command{6})
    errordlg('Lon Min box is empty','Error');   error = error + 1;    
end
if isempty(handles.command{8})
    errordlg('Lon Max box is empty','Error');   error = error + 1;    
end
if isempty(handles.command{10})
    errordlg('Lat Min box is empty','Error');   error = error + 1;    
end
if isempty(handles.command{12})
    errordlg('Lat Max box is empty','Error');   error = error + 1;    
end
if isempty(handles.command{15})
    errordlg('X_INC box is empty','Error');     error = error + 1;    
end
if isempty(handles.command{17})
    errordlg('Y_INC box is empty','Error');     error = error + 1;    
end

if (strcmp(handles.type,'nearneighbor'))
    if (isempty(get(handles.edit_S1_Neighbor,'String')))
        errordlg('Must give a value for "Search radius".','Error');
        error = error + 1;
    end
end

% --------------------------------------------------------------------
function Menu_Help_Callback(hObject, eventdata, handles)
switch handles.command{1};    
    case 'surface '
        message = {'surface reads randomly-spaced (x,y,z) triples from a file and'
        'produces a binary grdfile of gridded values z(x,y) by solving:'
        '               (1 - T) * L (L (z)) + T * L (z) = 0'
        'where T is a tension factor between 0 and 1, and L indicates the'
        'Laplacian operator. T = 0 gives the "minimum curvature" solution'
        'which is equivalent to SuperMISP and the ISM packages. Minimum'
        'curvature can cause undesired oscillations and false local maxima or'
        'minima (See Smith and Wessel, 1990), and you may wish to use T > 0'
        'to suppress these effects. Experience suggests T ~ 0.25 usually'
        'looks good for potential field data and T should be larger'
        '(T ~ 0.35) for steep topography data. T = 1 gives a harmonic'
        'surface (no maxima or minima are possible except at control data'
        'points). It is recommended that'' the user pre-process the data'
        'with blockmean, blockmedian, or blockmode to avoid spatial aliasing'
        'and eliminate redundant data. You may impose lower and/or upper'
        'bounds on the solution. These may be entered in the form of a fixed'
        'value, a grdfile with values, or simply be the minimum/maximum'
        'input data values.'};
        helpdlg(message,'Help on surface');        
    case 'triangulate '
        message = {'triangulate reads one or more ASCII [or binary] files'
        'containing x,y[,z] and performs Delauney triangulation, i.e., it find how'
        'the points should be connected to give the most equilateral triangulation'
        'possible. If a map projection is chosen then it is applied before the'
        'triangulation is calculated. By default, the output is triplets of point'
        'id numbers that make up each triangle and is written to standard output.'
        'The id numbers refer to the points position in the input file. As an'
        'option, you may choose to create a multiple segment file that can be piped'
        'through psxy to draw the triangulation network. If -G -I -R are set a grid'
        'will be calculated based on the surface defined by the planar triangles.'
        'The actual algorithm used in the triangulations is either that of Watson'
        '[1982] [Default] or Shewchuk [1996] (if installed). This choice is made'
        'during the GMT installation.'};
        helpdlg(message,'Help on triangulate');        
    case 'nearneighbor '
        message = {'nearneighbor reads arbitrarily located (x,y,z[,w]) triples [quadruplets]'
        'from standard input [or xyzfile(s)] and uses a nearest neighbor algorithm'
        'to assign an average value to each node that have one or more points'
        'within a radius centered on the node. The average value is computed'
        'as a weighted mean of the nearest point from each sector inside the'
        'search radius. The weighting function used is w(r) = 1.0 / (1 + d ^ 2),'
        'where d = 3 * r / search_radius and r is distance from the node. This'
        'weight is modulated by the observation points'' weights [if supplied].'};
        helpdlg(message,'Help on nearneighbor');
end

% -----------------------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata)
% Hint: delete(hObject) closes the figure
handles = guidata(hObject);
handles.output = '';        % User gave up, return nothing
guidata(hObject, handles);
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.figure1);
else    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end

% -----------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
% Check for "escape"
handles = guidata(hObject);
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = '';    % User said no by hitting escape
    guidata(hObject, handles);
    if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
        uiresume(handles.figure1);
    else    % The GUI is no longer waiting, just close it
        delete(handles.figure1);
    end
end

% --------------------------------------------------------------------
function about_window_Callback(hObject, eventdata, handles)
if (strcmp(handles.command{1},'surface '))
    about_box('surface_Last_modified__01_June_2008');
elseif (strcmp(handles.command{1},'triangulate '))
    about_box('triangulate_Last_modified_at_01_June_2008');
else
    about_box('nearneighbor_Last_modified_at_01_June_2008');
end

% --- Creates and returns a handle to the GUI figure. 
function griding_mir_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','griding_mir',...
'NumberTitle','off',...
'Position',[265.768111202607 206.47204765826 391 353],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 103 370 48],...
'Style','frame',...
'Tag','frame5');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 43 370 41],...
'Style','frame',...
'Tag','frame4');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 279 370 67],...
'Style','frame',...
'Tag','frame3');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 166 370 93],...
'Style','frame',...
'Tag','frame1');

uicontrol('Parent',h1,...
'Callback',{@griding_mir_uicallback,h1,'checkbox_Option_H_Callback'},...
'CData',[],...
'Position',[21 318 70 15],...
'String','Headers?',...
'Style','checkbox',...
'TooltipString','Are there any header lines in the input file?',...
'Tag','checkbox_Option_H',...
'UserData',[]);

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_nHeaders_Callback'},...
'HorizontalAlignment','left',...
'Position',[170 314 31 20],...
'String','1',...
'Style','edit',...
'TooltipString','How many?',...
'Tag','edit_nHeaders');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'popup_binInput_Callback'},...
'Position',[214 313 91 22],...
'String',{  'ascii'; 'binary double'; 'binary single' },...
'Style','popupmenu',...
'TooltipString','Input data type',...
'Value',1,...
'Tag','popup_binInput');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_binary_ncolumnIn_Callback'},...
'HorizontalAlignment','left',...
'Position',[304 314 25 20],...
'String','3',...
'Style','edit',...
'TooltipString','Number of columns in binary data',...
'Tag','edit_binary_ncolumnIn');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Callback',{@griding_mir_uicallback,h1,'pushbutton_Help_H_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[351 314 22 22],...
'String','?',...
'Tag','pushbutton_Help_H');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_InputFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[20 286 330 22],...
'Style','edit',...
'Tag','edit_InputFile');

uicontrol('Parent',h1,...
'Callback',{@griding_mir_uicallback,h1,'pushbutton_InputFile_Callback'},...
'Position',[350 285 23 23],...
'Tag','pushbutton_InputFile');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_Xmin_Callback'},...
'HorizontalAlignment','left',...
'Position',[76 219 71 21],...
'Style','edit',...
'TooltipString','X min value',...
'Tag','edit_Xmin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_Xmax_Callback'},...
'HorizontalAlignment','left',...
'Position',[152 219 71 21],...
'Style','edit',...
'TooltipString','X max value',...
'Tag','edit_Xmax');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_Xinc_Callback'},...
'HorizontalAlignment','left',...
'Position',[228 219 71 21],...
'Style','edit',...
'TooltipString','DX grid spacing',...
'Tag','edit_Xinc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_Ncols_Callback'},...
'HorizontalAlignment','left',...
'Position',[304 219 65 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Ncols');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_Ymin_Callback'},...
'HorizontalAlignment','left',...
'Position',[76 193 71 21],...
'Style','edit',...
'TooltipString','Y min value',...
'Tag','edit_Ymin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_Ymax_Callback'},...
'HorizontalAlignment','left',...
'Position',[152 193 71 21],...
'Style','edit',...
'TooltipString','Y max value',...
'Tag','edit_Ymax');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_Yinc_Callback'},...
'HorizontalAlignment','left',...
'Position',[228 193 71 21],...
'Style','edit',...
'TooltipString','DY grid spacing',...
'Tag','edit_Yinc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_Nrows_Callback'},...
'HorizontalAlignment','left',...
'Position',[304 193 65 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Nrows');

uicontrol('Parent',h1,...
'Callback',{@griding_mir_uicallback,h1,'checkbox_ToggleXY_Callback'},...
'Position',[76 172 75 19],...
'String','Toggle x,y',...
'Style','checkbox',...
'TooltipString','Toggle x and y columns',...
'Tag','checkbox_ToggleXY');

uicontrol('Parent',h1,...
'Callback',{@griding_mir_uicallback,h1,'checkbox_Option_F_Callback'},...
'Position',[169 173 101 15],...
'String','Pixel registration',...
'Style','checkbox',...
'TooltipString','Make image type grid',...
'Tag','checkbox_Option_F');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Callback',{@griding_mir_uicallback,h1,'pushbutton_Help_R_F_toggle_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[309 171 61 18],...
'String','?',...
'Tag','pushbutton_Help_R_F_toggle');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'edit_S1_Neighbor_Callback'},...
'HorizontalAlignment','left',...
'Position',[60 117 51 21],...
'Style','edit',...
'Tag','edit_S1_Neighbor');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'popup_S2_Neighbor_Callback'},...
'Position',[193 116 91 22],...
'String',{  ''; 'minutes'; 'seconds'; 'kilometers'; 'Kilometers' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_S2_Neighbor');

uicontrol('Parent',h1,...
'Callback',{@griding_mir_uicallback,h1,'pushbutton_Help_S_Callback'},...
'Position',[305 115 66 23],...
'String','Help',...
'Tag','pushbutton_Help_S');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@griding_mir_uicallback,h1,'popup_GridMethod_Callback'},...
'Position',[20 50 181 22],...
'String',{  'Minimum Curvature'; 'Delauny Triangulation'; 'Near Neighbor' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_GridMethod');

uicontrol('Parent',h1,...
'Callback',{@griding_mir_uicallback,h1,'pushbutton_Grid_Options_Callback'},...
'Position',[231 50 141 23],...
'String','Options',...
'Tag','pushbutton_Grid_Options');

uicontrol('Parent',h1,...
'Callback',{@griding_mir_uicallback,h1,'checkbox_Option_V_Callback'},...
'Position',[20 10 66 15],...
'String','Verbose',...
'Style','checkbox',...
'Tag','checkbox_Option_V');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[18 224 55 15],...
'String','X Direction',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[17 198 55 15],...
'String','Y Direction',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[169 241 41 13],...
'String','Max',...
'Style','text');

uimenu('Parent',h1,...
'Callback',{@griding_mir_uicallback,h1,'Menu_Help_Callback'},...
'Label','Help',...
'Tag','Menu_Help');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[91 242 41 13],...
'String','Min',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[246 242 41 13],...
'String','Spacing',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[314 242 51 13],...
'String','# of lines',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[30 75 90 15],...
'String','Griding Method',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[31 338 91 15],...
'String','Input Data File',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[30 252 121 15],...
'String','Griding Line Geometry',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'FontAngle','italic',...
'FontSize',9,...
'Position',[38 142 133 17],...
'String','For Nearneighbor only',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[17 113 40 30],...
'String',{  'Search'; 'Radius' },...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[149 120 41 15],...
'String','Unities',...
'Style','text');

uimenu('Parent',h1,...
'Callback',{@griding_mir_uicallback,h1,'about_window_Callback'},...
'Label','About',...
'Tag','about_window');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[101 317 67 15],...
'String','N? of headers',...
'Style','text',...
'TooltipString','How many?');

uicontrol('Parent',h1,...
'Callback',{@griding_mir_uicallback,h1,'pushbutton_OK_Callback'},...
'FontWeight','bold',...
'Position',[270 7 111 21],...
'String','Compute',...
'Tag','pushbutton_OK');

function griding_mir_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

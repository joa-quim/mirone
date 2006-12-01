function varargout = euler_stuff(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to euler_stuff (see VARARGIN) 

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
euler_stuff_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

%#function telha_m choosebox

movegui(hObject,'center');
global home_dir;    home_dir = pwd;
handles.path_data = [home_dir filesep 'data' filesep];
handles.path_continent = [home_dir filesep 'continents' filesep];
handles.h_line_orig = [];
handles.hLineSelected = [];
handles.p_lon = [];
handles.p_lat = [];
handles.p_omega = [];
handles.edit_pole1Lon = [];
handles.edit_pole1Lat = [];
handles.edit_pole1Ang = [];
handles.edit_pole2Lon = [];
handles.edit_pole2Lat = [];
handles.edit_pole2Ang = [];
handles.ages = [];
handles.do_interp = 0;          % Used to decide which 'compute' function to use
handles.finite_poles = [];      % Used to store a collection of finite poles (issued by choosebox)
handles.noKill = 0;

if (~isempty(varargin))
    handles.h_calling_fig = varargin{1};
    handles.mironeAxes = get(varargin{1},'CurrentAxes');
    if (length(varargin) == 2)          % Called with the line handle in argument
        handles.h_line_orig = varargin{2};
        handles.noKill = 1;
        set(handles.text_activeLine,'String','GOT A LINE TO WORK WITH','ForegroundColor',[0 0.8 0])
    end
else
    errordlg('EULER_STUFF: wrong number of arguments.','Error')
    delete(hObject);    return
end

% This is the tag that all tab push buttons share.  If you have multiple
% sets of tab push buttons, each group should have unique tag.
group_name = 'tab_group';

% This is a list of the UserData values used to link tab push buttons and
% the components on their linked panels.  To add a new tab panel to the group
%  Add the button using GUIDE
%  Assign the Tag based on the group name - in this case tab_group
%  Give the UserData a unique name - e.g. another_tab_panel
%  Add components to GUIDE for the new panel
%  Give the new components the same UserData as teh tab button
%  Add the new UserData name to the below cell array
panel_names = {'DoRotations','AddPoles','InterpPoles'};

% tabpanelfcn('makegroups',...) adds new fields to the handles structure,
% one for each panel name and another called 'group_name_all'.  These fields
% are used by the tabpanefcn when tab_group_handler is called.
handles = tabpanelfcn('make_groups',group_name, panel_names, handles, 1);

% Choose default command line output for euler_stuff_export
handles.output = hObject;
guidata(hObject, handles);

set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = euler_stuff_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = euler_stuff_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------------------
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over tab_group.
function tab_group_ButtonDownFcn(hObject, eventdata, handles)
% Call the tab_group_handler.  This updates visiblity of components as needed to
% hide the components from the previous tab and show components on this tab.
% This also updates the last_tab field in the handles structure to keep track
% of which panel was hidden.
handles = tabpanelfcn('tab_group_handler',hObject, handles, get(hObject, 'Tag'));
% Since this tab uses mostly existing uis, just make the visible here
if (strcmp(get(hObject,'UserData'),'InterpPoles'))
    set(handles.h_Stg_txt,'Visible','on','String','Finite rotation poles file')
    set(handles.edit_polesFile,'Visible','on')
    set(handles.pushbutton_readPolesFile,'Visible','on')
    set(handles.text1,'Visible','on')
    set(handles.edit_agesFile,'Visible','on')
    set(handles.pushbutton_ReadAgesFile,'Visible','on')
    set(handles.listbox_ages,'Visible','on')
    set(handles.pushbutton_polesList,'Visible','on')
    set(handles.pushbutton_compute,'Visible','on')
    handles.do_interp = 1;          % Redirect the 'compute' function
else
    handles.do_interp = 0;
    set(handles.h_Stg_txt,'String','Stage poles file')
end
guidata(hObject, handles);

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
function edit_agesFile_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if isempty(fname)    return;    end
id = strfind(fname,':');
if (~isempty(id))
    handles.ages = eval(fname);
    set(handles.listbox_ages,'String',mat2cell(handles.ages',length(handles.ages),1))
    guidata(hObject, handles);
else
    % Let the pushbutton_ReadAgesFile_Callback do all the work
    pushbutton_ReadAgesFile_Callback(hObject,[],guidata(gcbo),fname)
end

% -------------------------------------------------------------------------------------
function pushbutton_ReadAgesFile_Callback(hObject, eventdata, handles, opt)
% Read a file with ages where to compute the rotations
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
	str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = uigetfile(str1,'Select ages file');   pause(0.05)
    if (~isempty(last_dir)),    cd(home);   end
	if isequal(FileName,0)      return;    end
	fname = [PathName,FileName];
end

hFig = gcf;
[bin,n_column,multi_seg,n_headers] = guess_file(fname);
% If msgbox exist we have to move it from behind the main window. So get it's handle
hMsgFig = gcf;
if (hFig ~= hMsgFig)        uistack(hMsgFig,'top');   end   % If msgbox exists, bring it forward
% If error in reading file
if isempty(bin) & isempty(n_column) & isempty(multi_seg) & isempty(n_headers)
    errordlg(['Error reading file ' fname],'Error');    return
end
if (multi_seg ~= 0)   % multisegments are not spported
    errordlg('Multisegment files are not supported here.','Error');   return
end
if (bin == 0)   % ASCII
    fid = fopen(fname);
    todos = fread(fid,'*char');
    if (n_column == 1)
        handles.ages = strread(todos,'%f');
        handles.age_label = '';
    elseif (n_column == 2)
        [handles.ages handles.age_label] = strread(todos,'%f %s');
    else
        errordlg('Ages file can only have one OR two columns (IF 2 column: first column contains chron name)','Error')
        return
    end
    fclose(fid);
else        % BINARY
    errordlg('Sorry, binary files is not yet suported','Error');   return
end

if (~isempty(handles.age_label))    % We have labeled ages
    s1 = num2str(handles.ages);
    str = cell(size(handles.ages,1),1);
    for (k=1:size(handles.ages,1))
        str{k} = [handles.age_label{k} '    ' s1(k,1:end)];
    end
else
    str = {num2str(handles.ages)};
end
set(handles.listbox_ages,'String',str)

set(handles.edit_agesFile,'String',fname)
guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function listbox_ages_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns listbox_ages contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_ages

% -------------------------------------------------------------------------------------
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
    if (~handles.noKill),   delete(handles.h_line_orig);    end
    delete(handles.figure1)

% -------------------------------------------------------------------------------------
function pushbutton_compute_Callback(hObject, eventdata, handles)

if (handles.do_interp == 1)     % Compute interpolated poles instead
    cumpute_interp(handles)
    return
end

if (isempty(handles.h_line_orig))
    errordlg('Will you be so kind to let me know what line/point should I rotate?','Unknown target')
    return
end

if (get(handles.checkbox_singleRotation,'Value'))
    % Do the rotation using the pole parameters entered in the GUI and return
    if (isempty(handles.p_lon) || isempty(handles.p_lat) || isempty(handles.p_omega))
        return
    end
    for (i=1:numel(handles.h_line_orig))
        lon = get(handles.h_line_orig(i),'XData');
        lat = get(handles.h_line_orig(i),'YData');
        [rlon,rlat] = rot_euler(lon,lat,handles.p_lon,handles.p_lat,handles.p_omega);
        axes(handles.mironeAxes)       % Make the Mirone axes the CurrentAxes
        if (length(rlon) == 1)          % Single point rotation
            smb = get(handles.hLineSelected(i),'Marker');
            smb_fc = get(handles.hLineSelected(i),'MarkerFaceColor');
            smb_ec = get(handles.hLineSelected(i),'MarkerEdgeColor');
            smb_s = get(handles.hLineSelected(i),'MarkerSize');
            smb_t = get(handles.hLineSelected(i),'Linewidth');
            set(handles.h_line_orig(i),'XData',rlon,'YData',rlat,'Marker',smb,'MarkerFaceColor',smb_fc,...
                'MarkerEdgeColor',smb_ec,'MarkerSize',smb_s,'Linewidth',smb_t,'Tag','Rotated Line','Userdata',1);
        else
            lt = get(handles.hLineSelected(i),'LineWidth');
            lc = get(handles.hLineSelected(i),'Color');
            set(handles.h_line_orig(i),'XData',rlon,'YData',rlat,'Linewidth',lt,'Color',lc,'Tag','Rotated Line','Userdata',1);
        end
        line_info = {['Ang = ' num2str(handles.p_omega)]};
        draw_funs(handles.h_line_orig(i),'isochron',line_info)
    end
    handles.h_line_orig = [];       handles.hLineSelected = [];
    guidata(handles.figure1,handles)
    set(handles.text_activeLine,'String','NO ACTIVE LINE','ForegroundColor',[1 0 0])
    return
end

if (isempty(handles.ages))
    errordlg('I need to know the ages at which to compute the rotations','Error');    return
end

poles_name = get(handles.edit_polesFile,'String');
if (isempty(poles_name))
    errordlg('No stage poles provided','Error');    return
end

axes(handles.mironeAxes)       % Make the Mirone axes the CurrentAxes
opt_E = ['-E' poles_name];
opt_I = ' ';
if (get(handles.checkbox_revertRot,'Value'))
    opt_I = '-I';
end

for (i=1:numel(handles.h_line_orig))
	x = get(handles.h_line_orig(i),'XData');       y = get(handles.h_line_orig(i),'YData');
	linha = [x(:) y(:)];
	[out,n_data,n_seg,n_flow] = telha_m(linha, handles.ages, '-P', opt_E, opt_I);
    lt = get(handles.hLineSelected(i),'LineWidth');
    lc = get(handles.hLineSelected(i),'Color');
	if (length(linha) == 2)     % Only one point. Flow line mode
        aa = isnan(out(:,1));
        out = out(~aa,1:2);
        h_line = line('XData',out(:,1),'YData',out(:,2),'Linewidth',lt,'Color',lc,'Tag','Flow Line','Userdata',1);
        stg = get(handles.edit_polesFile,'String');
        [PATH,FNAME,EXT] = fileparts(stg);
        line_info = {['Stage file: ' FNAME EXT]};
	else
		h_line = zeros(n_flow,1);
		for (k=1:n_flow)              % For each time increment
            [x,y] = get_time_slice(out,n_data,n_seg,k);
            h_line(k) = line('XData',x,'YData',y,'Linewidth',lt,'Color',lc,'Tag','Rotated Line','Userdata',k);
		end
        line_info = get(handles.listbox_ages,'String');
	end
	draw_funs(h_line,'isochron',line_info)
end

set(handles.text_activeLine,'String','NO ACTIVE LINE','ForegroundColor',[1 0 0])
if (~handles.noKill),   delete(handles.h_line_orig);    end
handles.h_line_orig = [];       handles.hLineSelected = [];
guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function cumpute_interp(handles)
% Compute interpolated poles.
% This function is far from optimized, but it's not supposed to be applyied to long series either

if (isempty(handles.ages))
    errordlg('I need to know the ages at which to interpolate the poles','Error');    return
end

if ( strcmp(get(handles.edit_polesFile,'String'),'In memory poles') )
    poles = handles.finite_poles;
else
	poles_name = get(handles.edit_polesFile,'String');
	if (isempty(poles_name))
        errordlg('No poles file provided','Error');    return
	end
    poles = read_poles(poles_name);
    if (isempty(poles))      return;     end         % Bad poles file
end

if (size(poles,2) ~= 4)
    errordlg('The poles matrix MUST have 4 columns.','Error');    return
end
poles = sortrows(poles,4);          % To make sure they go from youngest to oldest

ages = handles.ages;
n_ages = length(ages);
pol = zeros(n_ages,4);
n_poles_in = size(poles,1);
n_new_finite = 0;
id = find(ages <= poles(1,4));      % Find ages <= first pole (different case)
if (~isempty(id))
    n_new_finite = length(id);      % Count number of interpolations of the first finite pole
    for (i = 1:n_new_finite)
        pol(i,:) = [poles(1,1:2) poles(1,3)*ages(i)/poles(1,4) ages(i)];
    end
    clear id;
end

id = (ages > poles(end,4));         % Find ages > last pole (they cannot be computed - no extrapolation)
if (~isempty(id))
    ages(id) = [];                  % This will only apply if we have extrapolation ages
    pol(id,:) = [];
    n_ages = length(ages);
end

for (i = n_new_finite+1:n_ages)
    id = find(ages(i) > poles(:,4));
    id = id(end);                   % We can only have one value and it's the last one that counts
    if (~isempty(id))
        t0 = poles(id,4);        t1 = poles(id+1,4);
        stg = finite2stages([poles(id,:); poles(id+1,:)], 1, 0);    % Compute the stage pole between t0 & t1
        frac = (poles(id+1,4) - ages(i)) / (poles(id+1,4) - poles(id,4));
        [pol(i,1) pol(i,2) pol(i,3)] = add_poles(poles(id+1,1),poles(id+1,2),poles(id+1,3),stg(1,1),stg(1,2), frac*stg(1,5));
        pol(i,4) = ages(i);                 % Give it its age
    else
        errordlg(['Error: age = ' num2str(ages(i)) ' does not fit inside poles ages interval'],'Error')
        return
    end
end

% Now the interpolated poles are on the antipodes (don't know why), so revert that
pol(n_new_finite+1:end,2:3) = -pol(n_new_finite+1:end,2:3);     % Change latitude & angle sign
pol(n_new_finite+1:end,1) = pol(n_new_finite+1:end,1) + 180;
id = pol(:,1) > 360;
pol(id,1) = pol(id,1) - 360;

str1 = {'*.dat;*.stg', 'Data file (*.dat,*.stg)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = uiputfile(str1,'Interp poles file');
if isequal(FileName,0);     return;     end
% Open and write to ASCII file
if ispc;        fid = fopen([PathName FileName],'wt');
elseif isunix;  fid = fopen([PathName FileName],'w');
else    errordlg('Unknown platform.','Error');      end
fprintf(fid,'#longitude\tlatitude\tangle(deg)\tage(Ma)\n');
fprintf(fid,'%9.5f\t%9.5f\t%7.4f\t%8.4f\n', pol');
fclose(fid);

% --------------------------------------------------------------------
function [x,y] = get_time_slice(data,n_data,n_seg,n,first)
i1 = (n-1)*(n_data + n_seg) + 2;
i2 = i1 + n_data - 1 + n_seg - 1;
x = data(i1:i2,1);    y = data(i1:i2,2);

% --------------------------------------------------------------------
function pushbutton_callMagBarCode_Callback(hObject, eventdata, handles)
MagBarCode([handles.path_data 'Cande_Kent_95.dat'])

% --------------------------------------------------------------------
function edit_poleLon_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx))
    set(hObject,'String','')
    handles.p_lon = [];
    return
end
handles.p_lon = xx;
guidata(hObject,handles)

% --------------------------------------------------------------------
function edit_poleLat_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx))
    set(hObject,'String','')
    handles.p_lat = [];
    return
end
handles.p_lat = xx;
guidata(hObject,handles)

% --------------------------------------------------------------------
function edit_poleAngle_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx))
    set(hObject,'String','')
    handles.p_omega = [];
    return
end
handles.p_omega = xx;
guidata(hObject,handles)

% --------------------------------------------------------------------
function checkbox_singleRotation_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.edit_poleLon,'Enable','on')
    set(handles.edit_poleLat,'Enable','on')
    set(handles.edit_poleAngle,'Enable','on')
else
    set(handles.edit_poleLon,'Enable','off')
    set(handles.edit_poleLat,'Enable','off')
    set(handles.edit_poleAngle,'Enable','off')
end

% --------------------------------------------------------------------
function pushbutton_polesList_Callback(hObject, eventdata, handles)
fid = fopen([handles.path_continent 'lista_polos.dat'],'rt');
c = fread(fid,'*char').';
fclose(fid);
s = strread(c,'%s','delimiter','\n');

multiple_str = 'multiple_finite';
if (handles.do_interp)      multiple_val = 1;
else                        multiple_val = 0;
end

[s,v] = choosebox('Name','One Euler list',...
                    'PromptString','List of poles:',...
                    'SelectString','Selected poles:',...
                    'ListSize',[380 300],...
                    multiple_str,multiple_val,...
                    'ListString',s);

if (v == 1)         % Finite pole (one only)
    handles.p_lon = s(1);
    handles.p_lat = s(2);
    handles.p_omega = s(3);
    set(handles.edit_poleLon, 'String', num2str(s(1)))
    set(handles.edit_poleLat, 'String', num2str(s(2)))
    set(handles.edit_poleAngle, 'String', num2str(s(3)))
    guidata(hObject,handles)
elseif (v == 2)     % Stage poles
    set(handles.edit_polesFile,'String',s)
elseif (v == 3)     % Multiple finite poles (with ages)
    set(handles.edit_polesFile,'String','In memory poles')
    handles.finite_poles = s;
    guidata(hObject,handles)
end

% -----------------------------------------------------------------------------------
function tab_group_Callback(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------
function pushbutton_tab_bg_Callback(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------
function push_pickLine_Callback(hObject, eventdata, handles)
    % Test if we have potential target lines and their type
    h_mir_lines = findobj(handles.h_calling_fig,'Type','line');     % Fish all objects of type line in Mirone figure
    if (isempty(h_mir_lines))                                       % We don't have any lines
        str = ['If you hited this button on purpose, than you deserve the following insult.',...
                'You #!|"*!%!?~^)--$&.',... 
                'THERE ARE NO LINES IN THAT FIGURE.'];
        errordlg(str,'Chico Clever');     return;
    end
    
    set(handles.h_calling_fig,'pointer','crosshair')
    h_line = get_polygon(handles.h_calling_fig);        % Get the line handle
    tf = ismember(h_line,handles.hLineSelected);        % Check that the line was not already selected
    if (tf)     % Repeated line
        set(handles.h_calling_fig,'pointer','arrow');   figure(handles.figure1);   return;
    end
    if (~isempty(h_line))
        c = get(h_line,'Color');
        t = get(h_line,'LineWidth');
        h = copyobj(h_line,handles.mironeAxes);
        set(h,'LineWidth',t+1,'Color',1-c)
        uistack_j(h,'bottom')
        handles.h_line_orig = [handles.h_line_orig; h];
        % Make a copy of the selected handles to be used in props recovering
        handles.hLineSelected = [handles.hLineSelected; h_line];
    end
    set(handles.h_calling_fig,'pointer','arrow')
    figure(handles.figure1)                 % Bring this figure to front again

	nl = numel(handles.h_line_orig);
	if (nl)
        set(handles.text_activeLine,'String',['GOT ' num2str(nl) ' LINE(S) TO WORK WITH'],'ForegroundColor',[0 0.8 0])
	else
        set(handles.text_activeLine,'String','NO ACTIVE LINE','ForegroundColor',[1 0 0])
	end
	guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function push_rectSelect_Callback(hObject, eventdata, handles)
    % Test if we have potential target lines and their type
    h_mir_lines = findobj(handles.h_calling_fig,'Type','line');     % Fish all objects of type line in Mirone figure
    if (isempty(h_mir_lines))                                       % We don't have any lines
        return
    end
    figure(handles.h_calling_fig)
    [p1,p2,hl] = rubberbandbox;
    delete(hl)
    figure(handles.figure1)         % Bring this figure fowrward again
    h = zeros(numel(h_mir_lines),1);
    hc = h;
    for (i=1:numel(h_mir_lines))    % Loop over lines to find out which cross the rectangle
        x = get(h_mir_lines(i),'XData');
        y = get(h_mir_lines(i),'YData');
        if ( (any(x >= p1(1)) || any(x <= p2(1))) && (any(y >= p1(2)) || any(y <= p2(2))) )
            tf = ismember(h_mir_lines(i),handles.hLineSelected);    % Check that the line was not already selected
            if (tf),    continue;     end                           % Repeated line
            c = get(h_mir_lines(i),'Color');
            t = get(h_mir_lines(i),'LineWidth');
            h(i) = copyobj(h_mir_lines(i),handles.mironeAxes);
            set(h(i),'LineWidth',t+1,'Color',1-c)
            uistack_j(h(i),'bottom')
            hc(i) = h_mir_lines(i);         % Make a copy of the selected handles to be used in props recovering
        end
    end
    h(h == 0) = [];     hc(hc == 0) = [];
    if (~isempty(h))
        handles.h_line_orig = [handles.h_line_orig; h];        % This is a bad name
        handles.hLineSelected = [handles.hLineSelected; hc];
        guidata(handles.figure1,handles)
    end
    set(handles.text_activeLine,'String',['GOT ' num2str(numel(h)) ' LINE(S) TO WORK WITH'],'ForegroundColor',[0 0.8 0])

% -----------------------------------------------------------------------------------
function edit_pole1Lon_Callback(hObject, eventdata, handles)
handles.edit_pole1Lon = str2double(get(hObject,'String'));
if (isnan(handles.edit_pole1Lon))   set(hObject,'String','');   return;     end
guidata(hObject, handles);
if (~got_them_all(handles))     return;     end     % Not yet all parameters of the 2 poles
[lon_s,lat_s,ang_s] = add_poles(handles.edit_pole1Lon,handles.edit_pole1Lat,handles.edit_pole1Ang,...
    handles.edit_pole2Lon,handles.edit_pole2Lat,handles.edit_pole2Ang);
set(handles.edit_pole3Lon,'String',num2str(lon_s,'%.4f'))
set(handles.edit_pole3Lat,'String',num2str(lat_s,'%.4f'))
set(handles.edit_pole3Ang,'String',num2str(ang_s,'%.4f'))

% -----------------------------------------------------------------------------------
function edit_pole1Lat_Callback(hObject, eventdata, handles)
handles.edit_pole1Lat = str2double(get(hObject,'String'));
if (isnan(handles.edit_pole1Lat))   set(hObject,'String','');   return;     end
guidata(hObject, handles);
if (~got_them_all(handles))     return;     end     % Not yet all parameters of the 2 poles
[lon_s,lat_s,ang_s] = add_poles(handles.edit_pole1Lon,handles.edit_pole1Lat,handles.edit_pole1Ang,...
    handles.edit_pole2Lon,handles.edit_pole2Lat,handles.edit_pole2Ang);
set(handles.edit_pole3Lon,'String',num2str(lon_s,'%.4f'))
set(handles.edit_pole3Lat,'String',num2str(lat_s,'%.4f'))
set(handles.edit_pole3Ang,'String',num2str(ang_s,'%.4f'))

% -----------------------------------------------------------------------------------
function edit_pole1Ang_Callback(hObject, eventdata, handles)
handles.edit_pole1Ang = str2double(get(hObject,'String'));
if (isnan(handles.edit_pole1Ang))   set(hObject,'String','');   return;     end
guidata(hObject, handles);
if (~got_them_all(handles))     return;     end     % Not yet all parameters of the 2 poles
[lon_s,lat_s,ang_s] = add_poles(handles.edit_pole1Lon,handles.edit_pole1Lat,handles.edit_pole1Ang,...
    handles.edit_pole2Lon,handles.edit_pole2Lat,handles.edit_pole2Ang);
set(handles.edit_pole3Lon,'String',num2str(lon_s,'%.4f'))
set(handles.edit_pole3Lat,'String',num2str(lat_s,'%.4f'))
set(handles.edit_pole3Ang,'String',num2str(ang_s,'%.4f'))

% -----------------------------------------------------------------------------------
function edit_pole2Lon_Callback(hObject, eventdata, handles)
handles.edit_pole2Lon = str2double(get(hObject,'String'));
if (isnan(handles.edit_pole2Lon))   set(hObject,'String','');   return;     end
guidata(hObject, handles);
if (~got_them_all(handles))     return;     end     % Not yet all parameters of the 2 poles
[lon_s,lat_s,ang_s] = add_poles(handles.edit_pole1Lon,handles.edit_pole1Lat,handles.edit_pole1Ang,...
    handles.edit_pole2Lon,handles.edit_pole2Lat,handles.edit_pole2Ang);
set(handles.edit_pole3Lon,'String',num2str(lon_s,'%.4f'))
set(handles.edit_pole3Lat,'String',num2str(lat_s,'%.4f'))
set(handles.edit_pole3Ang,'String',num2str(ang_s,'%.4f'))

% -----------------------------------------------------------------------------------
function edit_pole2Lat_Callback(hObject, eventdata, handles)
handles.edit_pole2Lat = str2double(get(hObject,'String'));
if (isnan(handles.edit_pole2Lat))   set(hObject,'String','');   return;     end
guidata(hObject, handles);
if (~got_them_all(handles))     return;     end     % Not yet all parameters of the 2 poles
[lon_s,lat_s,ang_s] = add_poles(handles.edit_pole1Lon,handles.edit_pole1Lat,handles.edit_pole1Ang,...
    handles.edit_pole2Lon,handles.edit_pole2Lat,handles.edit_pole2Ang);
set(handles.edit_pole3Lon,'String',num2str(lon_s,'%.4f'))
set(handles.edit_pole3Lat,'String',num2str(lat_s,'%.4f'))
set(handles.edit_pole3Ang,'String',num2str(ang_s,'%.4f'))

% -----------------------------------------------------------------------------------
function edit_pole2Ang_Callback(hObject, eventdata, handles)
handles.edit_pole2Ang = str2double(get(hObject,'String'));
if (isnan(handles.edit_pole2Ang))   set(hObject,'String','');   return;     end
guidata(hObject, handles);
if (~got_them_all(handles))     return;     end     % Not yet all parameters of the 2 poles
[lon_s,lat_s,ang_s] = add_poles(handles.edit_pole1Lon,handles.edit_pole1Lat,handles.edit_pole1Ang,...
    handles.edit_pole2Lon,handles.edit_pole2Lat,handles.edit_pole2Ang);
set(handles.edit_pole3Lon,'String',num2str(lon_s,'%.4f'))
set(handles.edit_pole3Lat,'String',num2str(lat_s,'%.4f'))
set(handles.edit_pole3Ang,'String',num2str(ang_s,'%.4f'))

% -----------------------------------------------------------------------------------
function yeap = got_them_all(handles)
% Check if we have all the 6 parameters (2 poles x 3 params each)
% If at least one of them is empty returns YEAP = 0;

yeap = 1;
if ( isempty(handles.edit_pole1Lon) | isempty(handles.edit_pole1Lat) | isempty(handles.edit_pole1Ang) | ...
        isempty(handles.edit_pole2Lon) | isempty(handles.edit_pole2Lat) | isempty(handles.edit_pole2Ang) )
    yeap = 0;
end

% -----------------------------------------------------------------------------------
function edit_pole3Lon_Callback(hObject, eventdata, handles)
% Nothing to do. Just wait for the result

% -----------------------------------------------------------------------------------
function edit_pole3Lat_Callback(hObject, eventdata, handles)
% Nothing to do. Just wait for the result

% -----------------------------------------------------------------------------------
function edit_pole3Ang_Callback(hObject, eventdata, handles)
% Nothing to do. Just wait for the result

% -----------------------------------------------------------------------------------------
function poles = read_poles(poles_file)
% Read a poles file (with ages also) and store it in a cell array

fid = fopen(poles_file,'r');
c = fread(fid,'*char').';
fclose(fid);
s = strread(c,'%s','delimiter','\n');
ix = strmatch('#',s);

hdr = s(ix);
n_hdr = length(hdr);
n_poles = length(s)-n_hdr;
poles = zeros(n_poles,4);
try
	for (i = 1:n_poles)
         tmp = sscanf(s{i+n_hdr}','%f',4);
         poles(i,1:4) = tmp';
	end
catch
    errordlg(['The file ' poles_file 'is not a properly formated Stage poles file.'],'Error');
    poles = [];
end

% -----------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata, handles)
    if (~handles.noKill),   delete(handles.h_line_orig);    end
    delete(handles.figure1);

% -----------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')
    if (~handles.noKill),   delete(handles.h_line_orig);    end
    delete(handles.figure1);
end

% --- Creates and returns a handle to the GUI figure. 
function euler_stuff_LayoutFcn(h1,handles);

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@euler_stuff_uicallback,h1,'figure1_KeyPressFcn'},...
'CloseRequestFcn',{@euler_stuff_uicallback,h1,'figure1_CloseRequestFcn'},...
'MenuBar','none',...
'Name','Euler stuff',...
'NumberTitle','off',...
'Position',[520 464 472 336],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

h2 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'tab_group_Callback'},...
'Enable','inactive',...
'Position',[102 310 91 23],...
'String','Add poles',...
'ButtonDownFcn',{@euler_stuff_uicallback,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','AddPoles');

h3 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'tab_group_Callback'},...
'Enable','inactive',...
'Position',[10 310 91 23],...
'String','Do Rotations',...
'ButtonDownFcn',{@euler_stuff_uicallback,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','DoRotations');

uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'tab_group_Callback'},...
'Enable','inactive',...
'Position',[194 310 91 23],...
'String','Interpolate poles',...
'ButtonDownFcn',{@euler_stuff_uicallback,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','InterpPoles');

h4 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'pushbutton_tab_bg_Callback'},...
'Enable','inactive',...
'Position',[10 11 451 301],...
'Tag','pushbutton_tab_bg');

h5 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_polesFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[20 226 211 21],...
'Style','edit',...
'Tag','edit_polesFile',...
'UserData','DoRotations');

h6 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback4,h1,[],'pushbutton_readPolesFile_Callback'},...
'FontSize',10,...
'FontWeight','bold',...
'Position',[230 226 21 21],...
'String','...',...
'Tag','pushbutton_readPolesFile',...
'UserData','DoRotations');

h7 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_agesFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[20 145 211 21],...
'Style','edit',...
'Tag','edit_agesFile',...
'ToolTipString','Enter either a filename with ages OR a ML command like: [1:5:30]',...
'UserData','DoRotations');

h8 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback4,h1,[],'pushbutton_ReadAgesFile_Callback'},...
'FontSize',10,...
'FontWeight','bold',...
'Position',[230 145 21 21],...
'String','...',...
'Tag','pushbutton_ReadAgesFile',...
'UserData','DoRotations');

h9 = uicontrol('Parent',h1,...
'Position',[20 171 51 15],...
'String','Age file',...
'Style','text',...
'Tag','text1',...
'UserData','DoRotations');

h10 = uicontrol('Parent',h1,...
'Position',[20 250 131 15],...
'String','Stage poles file',...
'Style','text',...
'Tag','h_Stg_txt',...
'UserData','DoRotations');

h11 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'listbox_ages_Callback'},...
'Position',[20 25 211 101],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_ages',...
'UserData','DoRotations');

h12 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'pushbutton_Cancel_Callback'},...
'Position',[296 27 66 23],...
'String','Cancel',...
'Tag','pushbutton_Cancel');

h13 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'pushbutton_compute_Callback'},...
'Position',[385 27 66 23],...
'String','Compute',...
'Tag','pushbutton_compute',...
'UserData','DoRotations');

h14 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'checkbox_revertRot_Callback'},...
'Position',[20 202 141 15],...
'String','Revert sense of rotation',...
'Style','checkbox',...
'TooltipString','Revert the sense of rotation defined by the stages poles',...
'Tag','checkbox_revertRot',...
'UserData','DoRotations');

h15 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'pushbutton_callMagBarCode_Callback'},...
'Position',[280 109 131 23],...
'String','Magnetic Bar Code',...
'TooltipString','Open the magnetic bar code window',...
'Tag','pushbutton_callMagBarCode',...
'UserData','DoRotations');

h16 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_poleLon_Callback'},...
'Enable','off',...
'Position',[280 226 51 21],...
'Style','edit',...
'TooltipString','Longitude of the Euler pole',...
'Tag','edit_poleLon',...
'UserData','DoRotations');

h17 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_poleLat_Callback'},...
'Enable','off',...
'Position',[340 226 51 21],...
'Style','edit',...
'TooltipString','Latitude of the Euler pole',...
'Tag','edit_poleLat',...
'UserData','DoRotations');

h18 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_poleAngle_Callback'},...
'Enable','off',...
'Position',[400 226 51 21],...
'Style','edit',...
'TooltipString','Angle of rotation',...
'Tag','edit_poleAngle',...
'UserData','DoRotations');

h19 = uicontrol('Parent',h1,...
'Position',[286 251 41 15],...
'String','Lon',...
'Style','text',...
'Tag','text3',...
'UserData','DoRotations');

h20 = uicontrol('Parent',h1,...
'Position',[344 251 41 15],...
'String','Lat',...
'Style','text',...
'Tag','text4',...
'UserData','DoRotations');

h21 = uicontrol('Parent',h1,...
'Position',[404 251 41 15],...
'String','Angle',...
'Style','text',...
'Tag','text5',...
'UserData','DoRotations');

h22 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'checkbox_singleRotation_Callback'},...
'Position',[280 202 91 15],...
'String','Use this Pole',...
'Style','checkbox',...
'Tag','checkbox_singleRotation',...
'UserData','DoRotations');

h23 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'pushbutton_polesList_Callback'},...
'Position',[280 159 131 23],...
'String','Poles selector',...
'Tag','pushbutton_polesList',...
'UserData','DoRotations');

h24 = uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'push_pickLine_Callback'},...
'Position',[20 279 161 23],...
'String','Pick line from Figure',...
'TooltipString','Allows you to mouse select one line from a Mirone figure',...
'Tag','togglebutton_pickLine',...
'UserData','DoRotations');

r=zeros(19,19,3)*NaN;   % Make a crude rectangle icon
r(4:17,3,1:3) = 0;      r(4:17,19,1:3) = 0;     % Verical lines
r(4,3:19,1:3) = 0;      r(17,3:19,1:3) = 0;
uicontrol('Parent',h1,...
'Callback',{@euler_stuff_uicallback,h1,'push_rectSelect_Callback'},...
'Position',[190 279 25 23],...
'CData',r,...
'TooltipString','Select objects inside a rectangular region',...
'Tag','push_rectSelect',...
'UserData','DoRotations');

h25 = uicontrol('Parent',h1,...
'FontSize',10,...
'FontWeight','Bold',...
'Position',[220 283 235 16],...
'String','NO ACTIVE LINE',...
'ForegroundColor',[1 0 0],...
'Style','text',...
'Tag','text_activeLine',...
'UserData','DoRotations');

h26 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_pole1Lon_Callback'},...
'Position',[50 208 51 21],...
'Style','edit',...
'TooltipString','Longitude of the first Euler pole',...
'Tag','edit_pole1Lon',...
'UserData','AddPoles');

h27 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_pole1Lat_Callback'},...
'Position',[110 208 51 21],...
'Style','edit',...
'TooltipString','Latitude of the first Euler pole',...
'Tag','edit_pole1Lat',...
'UserData','AddPoles');

h28 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_pole1Ang_Callback'},...
'Position',[170 208 51 21],...
'Style','edit',...
'TooltipString','Angle of rotation of first pole',...
'Tag','edit_pole1Ang',...
'UserData','AddPoles');

h29 = uicontrol('Parent',h1,...
'Position',[56 233 41 15],...
'String','Lon',...
'Style','text',...
'Tag','text7',...
'UserData','AddPoles');

h30 = uicontrol('Parent',h1,...
'Position',[114 233 41 15],...
'String','Lat',...
'Style','text',...
'Tag','text8',...
'UserData','AddPoles');

h31 = uicontrol('Parent',h1,...
'Position',[174 233 41 15],...
'String','Angle',...
'Style','text',...
'Tag','text9',...
'UserData','AddPoles');

h32 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_pole2Lon_Callback'},...
'Position',[260 209 51 21],...
'Style','edit',...
'TooltipString','Longitude of the second Euler pole',...
'Tag','edit_pole2Lon',...
'UserData','AddPoles');

h33 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_pole2Lat_Callback'},...
'Position',[320 209 51 21],...
'Style','edit',...
'TooltipString','Latitude of the second Euler pole',...
'Tag','edit_pole2Lat',...
'UserData','AddPoles');

h34 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_pole2Ang_Callback'},...
'Position',[380 209 51 21],...
'Style','edit',...
'TooltipString','Angle of rotation of the second pole',...
'Tag','edit_pole2Ang',...
'UserData','AddPoles');

h35 = uicontrol('Parent',h1,...
'Position',[266 234 41 15],...
'String','Lon','Style','text',...
'Tag','text10',...
'UserData','AddPoles');

h36 = uicontrol('Parent',h1,...
'Position',[324 234 41 15],...
'String','Lat',...
'Style','text',...
'Tag','text11',...
'UserData','AddPoles');

h37 = uicontrol('Parent',h1,...
'Position',[384 234 41 15],...
'String','Angle',...
'Style','text',...
'Tag','text12',...
'UserData','AddPoles');

h38 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[60 262 151 17],...
'String','First pole',...
'Style','text',...
'Tag','text13',...
'UserData','AddPoles');

h39 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[263 262 151 17],...
'String','Second pole',...
'Style','text',...
'Tag','text14',...
'UserData','AddPoles');

h40 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[195 176 101 16],...
'String','Result',...
'Style','text',...
'Tag','text15',...
'UserData','AddPoles');

h41 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_pole3Lon_Callback'},...
'Position',[123 121 71 21],...
'Style','edit',...
'TooltipString','Longitude of the resulting Euler pole',...
'Tag','edit_pole3Lon',...
'UserData','AddPoles');

h42 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_pole3Lat_Callback'},...
'Position',[209 121 71 21],...
'Style','edit',...
'TooltipString','Latitude of the resulting Euler pole',...
'Tag','edit_pole3Lat',...
'UserData','AddPoles');

h43 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@euler_stuff_uicallback,h1,'edit_pole3Ang_Callback'},...
'Position',[295 121 71 21],...
'Style','edit',...
'TooltipString','Angle of rotation',...
'Tag','edit_pole3Ang',...
'UserData','AddPoles');

h44 = uicontrol('Parent',h1,...
'Position',[149 146 41 15],...
'String','Lon',...
'Style','text',...
'Tag','text16',...
'UserData','AddPoles');

h45 = uicontrol('Parent',h1,...
'Position',[224 146 41 15],...
'String','Lat',...
'Style','text',...
'Tag','text17',...
'UserData','AddPoles');

h46 = uicontrol('Parent',h1,...
'Position',[299 146 41 15],...
'String','Angle',...
'Style','text',...
'Tag','text18',...
'UserData','AddPoles');

function euler_stuff_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

function euler_stuff_uicallback4(hObject, eventdata, h1, opt, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1),opt);

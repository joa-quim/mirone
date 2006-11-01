function varargout = swan_options(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to swan_options (see VARARGIN) 
%
% varargout will have a structure with the following fields:
%   --- If the grids were not already in memory
%   output.grid_Z_bat
%   output.grid_head_bat
%   output.grid_Z_src
%   output.grid_head_src
%   --- Allways
%   output.params
%   --- Optional
%   output.maregraph_xy
%   output.maregraph_data_name
%   output.opt_D
%   output.opt_G = '-G<namestem> (write grids)
%   output.opt_M = '-M';
%   output.opt_N
%   output.opt_m = '-m' (movie)
%   --- Or empty if user gave up

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
swan_options_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

global home_dir
if isempty(home_dir)    f_data = [pwd filesep 'data' filesep];
else                    f_data = [home_dir filesep 'data' filesep];   end

movegui(hObject,'east')

handles.head_bat = [];
handles.head_src = [];
handles.maregraph_xy = [];
handles.BothGridsInMemory = 0;  % Flag to signal that bat & deform arrays are already in memory
handles.BatGridInMemory = 0;    % Flag to signal that bat array is already in memory
handles.MaregraphInMemory = 1;  % Flag to signal that the maregraphs locations are already in memory
handles.got_params = 0;         % Flag that signals a (possibly correct) params file
handles.h_calling_fig = [];     % Handles to the calling figure
handles.is_tsun2 = 0;           % Flag to signal that options are for use in tsun2 code
handles.n_jump = 0;             % TSUN2 only. If > 0, the maregraphs array will start at n_jump + 1
last_dir = [];

% Import icons
load([f_data 'mirone_icons.mat'],'Mfopen_ico');
set(handles.pushbutton_BatGrid,'CData',Mfopen_ico)
set(handles.pushbutton_SourceGrid,'CData',Mfopen_ico)
set(handles.pushbutton_ParamsFile,'CData',Mfopen_ico)
set(handles.pushbutton_MaregraphPosFile,'CData',Mfopen_ico)
set(handles.pushbutton_MaregraphDataFile,'CData',Mfopen_ico)
clear Mfopen_ico;

if (~isempty(varargin))
    if (strcmp(varargin{1},'bat_and_deform_with_maregs') | strcmp(varargin{1},'bat_and_deform'))
        set(handles.edit_BatGrid,'String','In memory array','Enable','off')
        set(handles.pushbutton_BatGrid,'Enable','off')
        set(handles.edit_SourceGrid,'String','In memory array','Enable','off')
        set(handles.pushbutton_SourceGrid,'Enable','off')
        set(handles.edit_Jump_initialTime,'Enable','off')
        handles.BothGridsInMemory = 1;
        if (strcmp(varargin{1},'bat_and_deform_with_maregs'))
            set(handles.checkbox_WantMaregraphs,'Value',1)
            set(handles.edit_MaregraphPosFile,'String','In memory array','Enable','off')
            set(handles.pushbutton_MaregraphPosFile,'Enable','off')
            handles.MaregraphInMemory = 1;
            handles.maregraph_xy = varargin{2};     % need this for the error test logic
            handles.h_calling_fig = varargin{3};
            cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
            last_dir = cfig_handles.last_dir;
            set(handles.edit_MaregraphDataFile,'String',[last_dir 'maregs.dat'],'Enable','on')
        else
            handles.h_calling_fig = varargin{2};
            set(handles.edit_MaregraphPosFile,'String','','Enable','off')
            set(handles.edit_MaregraphDataFile,'String','','Enable','off')
            set(handles.pushbutton_MaregraphPosFile,'Enable','off')
            set(handles.pushbutton_MaregraphDataFile,'Enable','off')
        end
    elseif (strcmp(varargin{1},'bat_with_maregs') | strcmp(varargin{1},'bat_only'))
        handles.Z_bat = varargin{2};
        handles.head_bat = varargin{3};
        set(handles.edit_BatGrid,'String','In memory array','Enable','off')
        set(handles.pushbutton_BatGrid,'Enable','off')
        set(handles.edit_Jump_initialTime,'Enable','off')
        handles.BatGridInMemory = 1;
        if (strcmp(varargin{1},'bat_with_maregs'))
            set(handles.checkbox_WantMaregraphs,'Value',1)
            set(handles.edit_MaregraphPosFile,'String','In memory array','Enable','off')
            set(handles.pushbutton_MaregraphPosFile,'Enable','off')
            handles.MaregraphInMemory = 1;
            handles.maregraph_xy = varargin{4};     % need this for the error test logic
            handles.h_calling_fig = varargin{5};
            cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
            last_dir = cfig_handles.last_dir;
            set(handles.edit_MaregraphDataFile,'String',[last_dir 'maregs.dat'],'Enable','on')
        else
            handles.h_calling_fig = varargin{4};
            set(handles.edit_MaregraphPosFile,'String','','Enable','off')
            set(handles.edit_MaregraphDataFile,'String','','Enable','off')
            set(handles.pushbutton_MaregraphPosFile,'Enable','off')
            set(handles.pushbutton_MaregraphDataFile,'Enable','off')
        end
    elseif (strcmp(varargin{1},'Tsun2'))
        handles.is_tsun2 = 1;
        set(handles.edit_BatGrid,'String','In memory array','Enable','off')
        set(handles.pushbutton_BatGrid,'Enable','off')
        set(handles.edit_SourceGrid,'String','Not used in tsun2','Enable','off')
        set(handles.pushbutton_SourceGrid,'Enable','off')
        handles.BatGridInMemory = 1;
        set(handles.checkbox_WantMaregraphs,'Value',1)
        set(handles.edit_MaregraphPosFile,'Enable','on')
        set(handles.pushbutton_MaregraphPosFile,'Enable','on')
        set(handles.edit_MaregraphDataFile,'String','Not used in tsun2','Enable','off')
        set(handles.pushbutton_MaregraphDataFile,'Enable','off')
        set(gcf,'Name','Tsun2 options')
        handles.h_calling_fig = varargin{2};
    else
        handles.h_calling_fig = varargin{1};
    end
end

if (~isempty(handles.h_calling_fig))
    if (isempty(last_dir))      % This happens mostly when called in the Tsun2 mode
        cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
        last_dir = cfig_handles.last_dir;
    end
    set(handles.edit_gridNameStem,'String',[last_dir filesep 'tsu_time_'],'Enable','off')
end

% Choose default command line output for swan_options_export
handles.output = hObject;
guidata(hObject, handles);

set(hObject,'Visible','on');
% UIWAIT makes swan_options_export wait for user response (see UIRESUME)
uiwait(handles.figure1);

handles = guidata(hObject);
out = swan_options_OutputFcn(hObject, [], handles);
varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = swan_options_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;
% The figure can be deleted now
delete(handles.figure1);

% -----------------------------------------------------------------------------------------
function edit_BatGrid_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if ~isempty(fname)
    % Because GMT and Surfer share the .grd extension, find out which kind grid we are dealing with
    [fid, msg] = fopen(fname, 'r');
    if (fid < 0)    errordlg([fname ': ' msg],'ERROR');     return;    end
    ID = fread(fid,4,'*char');
    ID = strread(ID,'%s');
    if strcmp(ID,'DSBB') | strcmp(ID,'DSRB') 
        fname = [fname '=6'];
    elseif strcmp(ID,'DSAA')
        warndlg('I don''t know and do not intend to learn how to read ASCII Surfer grids.','Warning')
        return
    end
    [X,Y,handles.Z_bat,handles.head_bat] = grdread_m(fname);
    guidata(hObject,handles)
else
    set(hObject,'String','')
end

% -----------------------------------------------------------------------------------------
function pushbutton_BatGrid_Callback(hObject, eventdata, handles)
[FileName,PathName] = uigetfile({'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'},'Select GMT grid');
pause(0.01);
if isequal(FileName,0);     return;     end
fname = [PathName FileName];

% Because GMT and Surfer share the .grd extension, find out which kind grid we are dealing with
[fid, msg] = fopen(fname, 'r');
if fid < 0
    errordlg([PathName FileName ': ' msg],'ERROR'); return
end
ID = fread(fid,4,'*char');
ID = strread(ID,'%s');
if strcmp(ID,'DSBB') | strcmp(ID,'DSRB') 
    fname = [fname '=6'];
elseif strcmp(ID,'DSAA')
    warndlg('I don''t know and do not intend to learn how to read ASCII Surfer grids.','Warning')
    return
end
[X,Y,handles.Z_bat,handles.head_bat] = grdread_m(fname);
set(handles.edit_BatGrid,'String',fname)
guidata(hObject,handles)

% -----------------------------------------------------------------------------------------
function edit_SourceGrid_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if ~isempty(fname)
    % Because GMT and Surfer share the .grd extension, find out which kind grid we are dealing with
    [fid, msg] = fopen(fname, 'r');
    if fid < 0
        errordlg([fname ': ' msg],'ERROR'); return
    end
    ID = fread(fid,4,'*char');
    ID = strread(ID,'%s');
    if strcmp(ID,'DSBB') | strcmp(ID,'DSRB') 
        fname = [fname '=6'];
    elseif strcmp(ID,'DSAA')
        warndlg('I don''t know and do not intend to learn how to read ASCII Surfer grids.','Warning')
        return
    end
    [X,Y,handles.Z_src,handles.head_src] = grdread_m(fname);
    guidata(hObject,handles)
else
    set(hObject,'String','')
end

% -----------------------------------------------------------------------------------------
function pushbutton_SourceGrid_Callback(hObject, eventdata, handles)

if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
    cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
    last_dir = cfig_handles.last_dir;
    home = cfig_handles.home_dir;
else
    last_dir = [];
end

if (~isempty(last_dir)),    cd(last_dir);   end
[FileName,PathName] = uigetfile({'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'},'Select GMT grid');
pause(0.01);
if (~isempty(last_dir)),    cd(home);   end
if isequal(FileName,0);     return;     end
fname = [PathName FileName];

% Because GMT and Surfer share the .grd extension, find out which kind grid we are dealing with
[fid, msg] = fopen(fname, 'r');
if fid < 0
    errordlg([PathName FileName ': ' msg],'ERROR'); return
end
ID = fread(fid,4,'*char');
ID = strread(ID,'%s');
if strcmp(ID,'DSBB') | strcmp(ID,'DSRB') 
    fname = [fname '=6'];
elseif strcmp(ID,'DSAA')
    warndlg('I don''t know and do not intend to learn how to read ASCII Surfer grids.','Warning')
    return
end
[X,Y,handles.Z_src,handles.head_src] = grdread_m(fname);
set(handles.edit_SourceGrid,'String',fname)
guidata(hObject,handles)

% -----------------------------------------------------------------------------------------
function edit_SwanParams_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if ~isempty(fname)
    fid = fopen(fname,'r');
    if (fid < 0)
        errordlg('Error opening params file (it probably doesn''t exist)','Error');
        set(hObject,'String','')
        handles.got_params = 0;
        return
    end
    while 1
        tline = fgetl(fid);
        if (~strcmp(tline(1),'#'))      % Jump comment lines
            token = strtok(tline);
            handles.params(1) = str2double(token);
        end
    end
    handles.got_params = 1;     % Lets just hope that the info in it is correct
    fclose(fid);
else
    handles.got_params = 0;
end
guidata(hObject,handles)

% -----------------------------------------------------------------------------------------
function pushbutton_ParamsFile_Callback(hObject, eventdata, handles)

if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
    cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
    last_dir = cfig_handles.last_dir;
    home = cfig_handles.home_dir;
else
    last_dir = [];
end

if (~isempty(last_dir)),    cd(last_dir);   end
[FileName,PathName] = uigetfile({'*.par', 'Params file (*.par)';'*.*', 'All Files (*.*)'},'Select Swan parameter file');
pause(0.01);
if (~isempty(last_dir)),    cd(home);   end
if isequal(FileName,0);     return;     end
fname = [PathName FileName];
fid = fopen(fname,'r');
if (fid < 0)
    errordlg('Error opening params file (it probably doesn''t exist)','Error');
    set(hObject,'String','')
    handles.got_params = 0;
    return
end
i = 1;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if (~strcmp(tline(1),'#'))      % Jump comment lines
        token = strtok(tline);
        handles.params(i) = str2double(token);
        i = i + 1;
    end
end
handles.got_params = 1;     % Lets just hope that the info in it is correct
fclose(fid);
set(handles.edit_SwanParams,'String',fname)
guidata(hObject,handles)

% -----------------------------------------------------------------------------------------
function radiobutton_SurfaceLevel_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_MaxWater,'Value',0)
    set(handles.radiobutton_TotalWater,'Value',0)
else
    set(hObject,'Value',1)
    set(handles.radiobutton_MaxWater,'Value',0)
    set(handles.radiobutton_TotalWater,'Value',0)
end

% -----------------------------------------------------------------------------------------
function radiobutton_MaxWater_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_SurfaceLevel,'Value',0)
    set(handles.radiobutton_TotalWater,'Value',0)
else
    set(hObject,'Value',1)
    set(handles.radiobutton_SurfaceLevel,'Value',0)
    set(handles.radiobutton_TotalWater,'Value',0)
end

% -----------------------------------------------------------------------------------------
function radiobutton_TotalWater_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_SurfaceLevel,'Value',0)
    set(handles.radiobutton_MaxWater,'Value',0)
else
    set(hObject,'Value',1)
    set(handles.radiobutton_SurfaceLevel,'Value',0)
    set(handles.radiobutton_MaxWater,'Value',0)
end

% -----------------------------------------------------------------------------------------
function edit_Number_of_cycles_Callback(hObject, eventdata, handles)
% The OK button will get the number of cycles from here. So make sure it gives a possible number
xx = get(hObject,'String');
if isnan(str2double(xx))
    set(hObject,'String','1010')
elseif (str2double(xx) <= 0)
    set(hObject,'String','1010')
else
    set(hObject,'String',num2str(fix(str2double(xx))) ) % Insure that it is an integer
end

% -----------------------------------------------------------------------------------------
function edit_Jump_initialTime_Callback(hObject, eventdata, handles)
% Remove this number of seconds from the beguining of the maregraphs (if they were loaded)
jmp = str2double(get(hObject,'String'));
if (isnan(jmp) | jmp < 0 | isempty(handles.maregraph_xy))
    set(hObject,'String','0');
    handles.n_jump = 0;
else
    [m,n] = size(handles.maregraph_xy);
    dt1 = diff(handles.maregraph_xy(:,1));    dt2 = diff(dt1);
    if any(dt2 ~= 0)            % First column doesn't have the time
    else                        % First column has the time.
        dt = dt1(1);            % This is the time increment
        n_jmp = fix(jmp / dt);  % Number of records to jump
        if ((m - n_jmp) < 10)  % Stupid jump time choice
            set(hObject,'String','0')
            handles.n_jump = 0;
            return
        end
        handles.n_jump = n_jmp;
        set(handles.edit_Number_of_cycles,'String',num2str(m-n_jmp))    % Uppdate max pssible
        guidata(hObject,handles)
    end
end

% -----------------------------------------------------------------------------------------
function checkbox_WantMaregraphs_Callback(hObject, eventdata, handles)
% In tsun2 maregraphs are mandatory (for reading), so never let the user uncheck
if (handles.is_tsun2 == 1),     set(hObject,'Value',1);     return;     end

if (get(hObject,'Value'))
    set(handles.edit_MaregraphPosFile,'Enable','on')
    set(handles.edit_MaregraphDataFile,'Enable','on')
    set(handles.pushbutton_MaregraphPosFile,'Enable','on')
    set(handles.pushbutton_MaregraphDataFile,'Enable','on')
else
    set(handles.edit_MaregraphPosFile,'String','','Enable','off')
    set(handles.edit_MaregraphDataFile,'String','','Enable','off')
    set(handles.pushbutton_MaregraphPosFile,'Enable','off')
    set(handles.pushbutton_MaregraphDataFile,'Enable','off')
end

% -----------------------------------------------------------------------------------------
function edit_MaregraphPosFile_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if isempty(fname),  return;     end
hFig = gcf;
[bin,n_column,multi_seg,n_headers] = guess_file(fname);
% If msgbox exist we have to move it from behind the main window. So get it's handle
hMsgFig = gcf;
if (hFig ~= hMsgFig)        uistack(hMsgFig,'top');   end   % If msgbox exists, bring it forward
% If error in reading file
if isempty(bin) & isempty(n_column) & isempty(multi_seg) & isempty(n_headers)
    errordlg(['Error reading file ' fname],'Error');    return
end
if multi_seg ~= 0   % multisegments are not spported
    errordlg('Multisegment files are yet not supported.','Error');   return
end
if (bin == 0)   % ASCII
    if n_column < 2
        errordlg('File error. Your file doesn''t have at least 2 columns','Error'); return
    end
    handles.maregraph_xy = read_xy(fname,n_column,n_headers);
    if (hFig ~= hMsgFig);       figure(hFig);   end     % gain access to the drawing figure
    [nr,nc] = size(handles.maregraph_xy);
    if (nr == 0)
        errordlg('Your file is empty.','Chico Clever');   return
    end
else        % BINARY
    errordlg('Sorry, reading binary files is not yet programed','Error');   return
end
guidata(hObject,handles)

% -----------------------------------------------------------------------------------------
function pushbutton_MaregraphPosFile_Callback(hObject, eventdata, handles)

if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
    cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
    last_dir = cfig_handles.last_dir;
    home = cfig_handles.home_dir;
else
    last_dir = [];
end

if (~isempty(last_dir)),    cd(last_dir);   end
[FileName,PathName] = uigetfile({'*.dat;*.DAT;*.xy', 'Maregraph location (*.dat,*.DAT,*.xy)';'*.*', 'All Files (*.*)'},'Select Maregraphs position');
pause(0.01);
if (~isempty(last_dir)),    cd(home);   end
if isequal(FileName,0);     return;     end
fname = [PathName FileName];
hFig = gcf;
set(handles.figure1,'pointer','watch')
[bin,n_column,multi_seg,n_headers] = guess_file(fname);
% If msgbox exist we have to move it from behind the main window. So get it's handle
hMsgFig = gcf;
if (hFig ~= hMsgFig)        uistack(hMsgFig,'top');   end   % If msgbox exists, bring it forward
% If error in reading file
if isempty(bin) & isempty(n_column) & isempty(multi_seg) & isempty(n_headers)
    set(handles.figure1,'pointer','arrow');
    errordlg(['Error reading file ' fname],'Error');    return
end
if multi_seg ~= 0   % multisegments are not spported
    set(handles.figure1,'pointer','arrow');
    errordlg('Multisegment files are yet not supported.','Error');   return
end
if (bin == 0)   % ASCII
    if (n_column < 2)
        set(handles.figure1,'pointer','arrow');
        errordlg('File error. Your file doesn''t have at least 2 columns','Error'); return
    end
    set(handles.figure1,'pointer','watch')
    handles.maregraph_xy = read_xy(fname,n_column,n_headers);
    set(handles.figure1,'pointer','arrow');
    set(handles.edit_Jump_initialTime,'String','0')     % Reset this anyway
    if (hFig ~= hMsgFig);       figure(hFig);   end     % gain access to the drawing figure
    [nr,nc] = size(handles.maregraph_xy);
    if (nr == 0)
        errordlg('Your file is empty.','Chico Clever');   return
    end
else        % BINARY
    errordlg('Sorry, reading binary files is not yet programed','Error');   return
end
set(handles.edit_Number_of_cycles,'String',num2str(nr))
set(handles.edit_MaregraphPosFile,'String',fname)
guidata(hObject,handles)

% -----------------------------------------------------------------------------------------
function edit_MaregraphDataFile_Callback(hObject, eventdata, handles)
% Nothing to test or program here. The OK button will read the contents of this box

% -----------------------------------------------------------------------------------------
function pushbutton_MaregraphDataFile_Callback(hObject, eventdata, handles)

if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
    cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
    last_dir = cfig_handles.last_dir;
    home = cfig_handles.home_dir;
else
    last_dir = [];
end

if (~isempty(last_dir)),    cd(last_dir);   end
[FileName,PathName] = uigetfile({'*.dat;*.DAT;*.xy', 'Maregraph data file (*.dat,*.DAT,*.xy)';'*.*', 'All Files (*.*)'},'Select Maregraph');
pause(0.01);
if (~isempty(last_dir)),    cd(home);   end
if isequal(FileName,0);     return;     end
fname = [PathName FileName];
set(handles.edit_MaregraphDataFile,'String',fname)

% -----------------------------------------------------------------------------------------
function checkbox_OutputGrids_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.edit_gridNameStem,'Enable','on')
else
    set(handles.edit_gridNameStem,'String','tsu_time_')
    set(handles.edit_gridNameStem,'Enable','off')
end

% -----------------------------------------------------------------------------------------
function edit_gridNameStem_Callback(hObject, eventdata, handles)
% Nothing to test or program here. The OK button will read the contents of this box

% -----------------------------------------------------------------------------------------
function checkbox_MakeMovie_Callback(hObject, eventdata, handles)
% Nothing to test or program here. The OK button will read the status of this box

% -----------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
% Do some error checking

if (handles.is_tsun2 == 1)    % Tests for the case of tsun2 (and return)
    error = check_errors_tsun2(handles);
    if (error)
        handles.output = [];
        guidata(hObject,handles);
        return
    end
    handles.output.params_file_name = get(handles.edit_SwanParams,'String');  % desenrasque
    % Reading maregraphs was set?
    if (~isempty(handles.maregraph_xy))
        if (handles.n_jump)     % We have a time jump
            handles.maregraph_xy = handles.maregraph_xy(handles.n_jump+1:end,:);
        end
        handles.output.maregraph_xy = handles.maregraph_xy;
    end
    if (get(handles.radiobutton_MaxWater,'Value'))
        handles.output.opt_M = '-M';
    end
    if (get(handles.radiobutton_TotalWater,'Value'))
        handles.output.opt_D = '-D';
    end
    % Do we want a collection of intermediary tsunami grid steps?
    if (get(handles.checkbox_OutputGrids,'Value'))
        handles.output.opt_G = ['-G' get(handles.edit_gridNameStem,'String')];
    end
    % Want a movie? 
    if (get(handles.checkbox_MakeMovie,'Value'))
        handles.output.opt_m = '-m';
    end
    % Get number of cycles
    handles.output.opt_N = ['-N' get(handles.edit_Number_of_cycles,'String')];
    guidata(hObject,handles);
    uiresume(handles.figure1);
    return
end
        
error = check_errors_swan(handles);
if (error)      % fdeu-se
    handles.output = [];
    guidata(hObject,handles);
    return
end

if (handles.BothGridsInMemory == 0)         % Otherwise Mirone already knows about it
    handles.output.grid_Z_bat = handles.Z_bat;
    handles.output.grid_head_bat = handles.head_bat;
    handles.output.grid_Z_src = handles.Z_src;
    handles.output.grid_head_src = handles.head_src;
end

if (handles.BatGridInMemory == 1)
    handles.output.grid_Z_src = handles.Z_src;
    handles.output.grid_head_src = handles.head_src;
end

% This one is allways read in this program
handles.output.params = handles.params;

% See if computation of maregraphs was required
if (~isempty(handles.maregraph_xy))
    handles.output.maregraph_xy = handles.maregraph_xy;
    handles.output.maregraph_data_name = get(handles.edit_MaregraphDataFile,'String');
end

% Check what is going to be computed (surface level OR Max water OR Total water depths)
% Remember, Surface level is the default
if (get(handles.radiobutton_MaxWater,'Value'))
    handles.output.opt_M = '-M';
end
if (get(handles.radiobutton_TotalWater,'Value'))
    handles.output.opt_D = '-D';
end    

% Do we want a collection of intermediary tsunami grid steps?
if (get(handles.checkbox_OutputGrids,'Value'))
    handles.output.opt_G = ['-G' get(handles.edit_gridNameStem,'String')];
end

% Want a movie? 
if (get(handles.checkbox_MakeMovie,'Value'))
    handles.output.opt_m = '-m';
end    

% Get number of cycles
handles.output.opt_N = ['-N' get(handles.edit_Number_of_cycles,'String')];
guidata(hObject,handles);
uiresume(handles.figure1);

% -----------------------------------------------------------------------------------------
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
handles.output = [];        % User gave up, return nothing
guidata(hObject, handles);  uiresume(handles.figure1);

%---------------------------------------------------------------------------------------------------
function error = check_errors_swan(handles)
error = 0;
small = 1e-6;   % Used in postion comparation. It places the accuracy at sub-meter level
if (handles.BothGridsInMemory == 0)     % If arrays are already in Mirone's memory there is no need to test this
    if (isempty(handles.head_bat) & handles.BatGridInMemory == 0)
        errordlg('No bathymetry grid provided (Ha! was it realy necessary?)','Error')
        error = 1;
    end
    if (isempty(handles.head_src))
        errordlg('No tsunami source grid provided (Yes, tsunamis need to be induced)','Error')
        error = 1;
    end
end
if (handles.BatGridInMemory == 1)     %
    if (isempty(handles.head_src))
        errordlg('No tsunami source grid provided (Yes, tsunamis need to be induced)','Error')
        error = 1;
    end
end    

if (handles.got_params == 0)        % Error in params file
    errordlg('Tsunami propagation parameter file is not valid or inexistent','Error')
    error = 1;
end
    
if (get(handles.checkbox_WantMaregraphs,'Value'))       % That is, requested maregraphs computation
    if (isempty(handles.maregraph_xy) & handles.MaregraphInMemory == 0)
        errordlg('Where are your maregraphs? On Moon? (In case you forgot, it has no oceans)','Error')
        error = 1;
    end
    if (isempty(get(handles.edit_MaregraphDataFile,'String')))
        errordlg('You need to tell me the file name where I''ll write the maregraphs water hight','Error')
        error = 1;
    end
end

% Cheeck that at lest one operation was selected
if (~get(handles.checkbox_WantMaregraphs,'Value') & ~get(handles.checkbox_MakeMovie,'Value') & ...
        ~get(handles.checkbox_OutputGrids,'Value'))
    errordlg('You need to select at least one operation','Error')
    error = 1;
end

if (error),     return;     end

% If we reach here, we can now check if grids are compatible. However, if grids where already in
% Mirone's memory, it is there that this test should have been made
if (handles.BothGridsInMemory == 0 | handles.BatGridInMemory == 1)
    if ( abs(handles.head_src(1) - handles.head_bat(1)) > small | ...
           abs(handles.head_src(2) - handles.head_bat(2)) > small | ...
           abs(handles.head_src(3) - handles.head_bat(3)) > small | ...
           abs(handles.head_src(4) - handles.head_bat(4)) > small )
        errordlg('Bathymetry & Source grids do not cover the same region','Error');
        error = 1;
    elseif ( abs(handles.head_src(8) - handles.head_bat(8)) > small | ...
           abs(handles.head_src(9) - handles.head_bat(9)) > small )
        errordlg('Bathymetry & Source grids have different nx and/or ny','Error');
        error = 1;
    end
end

%---------------------------------------------------------------------------------------------------
function error = check_errors_tsun2(handles)
% Check errors when used for tsun2 options
error = 0;

if (handles.got_params == 0)        % Error in params file
    errordlg('Tsunami propagation parameter file is not valid or inexistent','Error')
    error = 1;
end

if (isempty(handles.maregraph_xy))
    errordlg('You need to read a maregraph file','Error')
    error = 1;
end

% Cheeck that at lest one operation was selected
if (~get(handles.checkbox_MakeMovie,'Value') & ~get(handles.checkbox_OutputGrids,'Value'))
    errordlg('You need to select at least one the two operation: grids or movie','Error')
    error = 1;
end

% -----------------------------------------------------------------------------------------
function xy = read_xy(file,n_col,n_head)
% build the format string to read the data n_columns
xy = [];    format = [];    fid = fopen(file,'r');
for (i=1:n_col),    format = [format '%f '];    end
% Jump header lines
for (i = 1:n_head),    tline = fgetl(fid);  end

todos = fread(fid,'*char');
xy = sscanf(todos,format,[n_col inf])';    % After hours strugling agains this FILHO DA PUTA, I may have found
fclose(fid);

% -----------------------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    handles.output = [];        % User gave up, return nothing
    guidata(hObject, handles);    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    handles.output = [];        % User gave up, return nothing
    guidata(hObject, handles);    delete(handles.figure1);
end

% -----------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = [];    % User said no by hitting escape
    guidata(hObject, handles);    uiresume(handles.figure1);
end

% -----------------------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function swan_options_LayoutFcn(h1,handles);

set(h1, 'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','swan_options',...
'NumberTitle','off',...
'Position',[520 401 351 370],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1');

h2 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_BatGrid_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 341 271 21],...
'Style','edit','Tag','edit_BatGrid');

h3 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_BatGrid_Callback'},...
'Position',[321 339 23 23],...
'Tag','pushbutton_BatGrid');

h4 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_SourceGrid_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 311 271 21],...
'Style','edit','Tag','edit_SourceGrid');

h5 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'checkbox_OutputGrids_Callback'},...
'Position',[10 112 77 15],...
'String','Output grids','Style','checkbox',...
'Tag','checkbox_OutputGrids');

h6 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_gridNameStem_Callback'},...
'HorizontalAlignment','left',...
'Position',[95 108 225 21],...
'Style','edit',...
'TooltipString','Grids are numbered after this name stem',...
'Tag','edit_gridNameStem');

h7 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'radiobutton_SurfaceLevel_Callback'},...
'Position',[10 235 82 15],...
'String','Surface level','Style','radiobutton',...
'Value',1,'Tag','radiobutton_SurfaceLevel');

h8 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'radiobutton_MaxWater_Callback'},...
'Position',[131 235 82 15],...
'String','Max water','Style','radiobutton',...
'Tag','radiobutton_MaxWater');

h9 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'radiobutton_TotalWater_Callback'},...
'Position',[249 235 82 15],...
'String','Total water','Style','radiobutton',...
'Tag','radiobutton_TotalWater');

h10 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_Number_of_cycles_Callback'},...
'Position',[270 78 51 21],...
'String','1010','Style','edit',...
'TooltipString','Use this number of cycles from the Maregraph file',...
'Tag','edit_Number_of_cycles');

h11 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_MaregraphPosFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 179 271 21],...
'Style','edit',...
'TooltipString','Name of the file with maregraph locations',...
'Tag','edit_MaregraphPosFile');

h12 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_MaregraphPosFile_Callback'},...
'Position',[321 178 23 23],...
'Tag','pushbutton_MaregraphPosFile');

h13 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_MaregraphDataFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 149 271 21],...
'Style','edit',...
'TooltipString','Name of the file that will contain the maregraphs water height',...
'Tag','edit_MaregraphDataFile');

h14 = uicontrol('Parent',h1,...
'Position',[209 82 60 15],...
'String','Nº of cycles','Style','text','Tag','text1');

h15 = uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[10 345 40 15],...
'String','Bat','Style','text','Tag','text2');

h16 = uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[10 314 40 15],...
'String','Source','Style','text','Tag','text3');

h17 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_MaregraphDataFile_Callback'},...
'Position',[321 148 23 23],...
'Tag','pushbutton_MaregraphDataFile');

h18 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_SwanParams_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 270 271 21],...
'Style','edit',...
'TooltipString','Either Swan or Tsun2 need a parameter file.',...
'Tag','edit_SwanParams');

h19 = uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[10 273 40 15],...
'String','Params','Style','text','Tag','text4');

h20 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_OK_Callback'},...
'FontWeight','bold',...
'Position',[195 8 66 23],...
'String','OK','Tag','pushbutton_OK');

h21 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_Cancel_Callback'},...
'FontWeight','bold',...
'Position',[275 8 66 23],...
'String','Cancel',...
'Tag','pushbutton_Cancel');

h22 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'checkbox_WantMaregraphs_Callback'},...
'Position',[10 202 76 15],...
'String','Maregraphs','Style','checkbox',...
'TooltipString','Check this if you want to compute water height at maregraphs',...
'Tag','checkbox_WantMaregraphs');

h23 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_SourceGrid_Callback'},...
'Position',[321 310 23 23],...
'Tag','pushbutton_SourceGrid');

h24 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_ParamsFile_Callback'},...
'Position',[321 269 23 23],...
'Tag','pushbutton_ParamsFile');

h25 = uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'checkbox_MakeMovie_Callback'},...
'Position',[10 81 85 15],...
'String','Make a movie','Style','checkbox',...
'Tag','checkbox_MakeMovie');

h26 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_Jump_initialTime_Callback'},...
'Position',[270 48 51 21],...
'String','0','Style','edit',...
'TooltipString','Jump these number of seconds from the beguining of Maregraphs file',...
'Tag','edit_Jump_initialTime');

h27 = uicontrol('Parent',h1,'HorizontalAlignment','right',...
'Position',[212 52 55 15],...
'String','Jump initial','Style','text','Tag','text5');

h28 = uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[13 181 30 16],...
'String','In file','Style','text','Tag','text6');

h29 = uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[12 152 35 16],...
'String','Out file','Style','text','Tag','text7');

function swan_options_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

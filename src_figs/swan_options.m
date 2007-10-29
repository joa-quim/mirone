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
	%   --- Optional (not so much)
	%   output.maregraph_xy
	%   output.maregraph_data_name
	%   output.opt_D
	%   output.opt_F	(tsun2 only - friction)
	%   output.opt_O	(tsun2 only - out maregs locations)
	%   output.opt_G = '-G<namestem> (write grids)
	%   output.opt_M = '-M';
	%   output.opt_N
	%   output.opt_m = '-m' (movie)
	%   output.opt_S = '-S' (swan only - write momentums)
	%   output.opt_s = '-s' (swan only - write velocities)
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

	if (isempty(varargin))
		errordlg('Bad call to swan_option. Please warn me (the author)','Error')
		return
	end
	
	hObject = figure('Tag','figure1','Visible','off');
	swan_options_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'east')

	handles.head_bat = [];
	handles.head_src = [];
	handles.maregraph_xy = [];
	handles.BothGridsInMemory = 0;  % Flag to signal that bat & deform arrays are already in memory
	handles.BatGridInMemory = 0;    % Flag to signal that bat array is already in memory
	handles.MaregraphInMemory = 1;  % Flag to signal that the maregraphs locations are already in memory
	handles.got_params = 0;         % Flag that signals a (possibly correct) params file
	handles.is_tsun2 = 0;           % Flag to signal that options are for use in tsun2 code
	handles.n_jump = 0;             % TSUN2 only. If > 0, the maregraphs array will start at n_jump + 1
	last_dir = [];
	
	handMir = varargin{1};
	handles.hCallingFig = handMir.figure1;
	handles.home_dir = handMir.home_dir;
	handles.work_dir = handMir.work_dir;
	handles.last_dir = handMir.last_dir;

	% Import icons
	f_data = [handles.home_dir filesep 'data' filesep];
	load([f_data 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.pushbutton_BatGrid,'CData',Mfopen_ico)
	set(handles.pushbutton_SourceGrid,'CData',Mfopen_ico)
	set(handles.pushbutton_ParamsFile,'CData',Mfopen_ico)
	set(handles.pushbutton_MaregraphPosFile,'CData',Mfopen_ico)
	set(handles.pushbutton_MaregraphDataFile,'CData',Mfopen_ico)
	clear Mfopen_ico;

	if (numel(varargin) > 1)
		if (strcmp(varargin{2},'bat_and_deform_with_maregs') || strcmp(varargin{2},'bat_and_deform'))
			set(handles.edit_BatGrid,'String','In memory array','Enable','off')
			set(handles.pushbutton_BatGrid,'Enable','off')
			set(handles.edit_SourceGrid,'String','In memory array','Enable','off')
			set(handles.pushbutton_SourceGrid,'Enable','off')
			set(handles.edit_Jump_initialTime,'Enable','off')
			handles.BothGridsInMemory = 1;
			if (strcmp(varargin{2},'bat_and_deform_with_maregs'))
				set(handles.checkbox_WantMaregraphs,'Value',1)
				set(handles.edit_MaregraphPosFile,'String','In memory array','Enable','off')
				set(handles.pushbutton_MaregraphPosFile,'Enable','off')
				handles.MaregraphInMemory = 1;
				handles.maregraph_xy = varargin{3};     % need this for the error test logic
				set(handles.edit_MaregraphDataFile,'String',[handMir.last_dir 'maregs.dat'],'Enable','on')
			else
				set(handles.edit_MaregraphPosFile,'String','','Enable','off')
				set(handles.edit_MaregraphDataFile,'String','','Enable','off')
				set(handles.pushbutton_MaregraphPosFile,'Enable','off')
				set(handles.pushbutton_MaregraphDataFile,'Enable','off')
			end
		elseif (strcmp(varargin{2},'bat_with_maregs') || strcmp(varargin{2},'bat_only'))
			handles.Z_bat = varargin{3};
			handles.head_bat = varargin{4};
			set(handles.edit_BatGrid,'String','In memory array','Enable','off')
			set(handles.pushbutton_BatGrid,'Enable','off')
			set(handles.edit_Jump_initialTime,'Enable','off')
			handles.BatGridInMemory = 1;
			if (strcmp(varargin{2},'bat_with_maregs'))
				set(handles.checkbox_WantMaregraphs,'Value',1)
				set(handles.edit_MaregraphPosFile,'String','In memory array','Enable','off')
				set(handles.pushbutton_MaregraphPosFile,'Enable','off')
				handles.MaregraphInMemory = 1;
				handles.maregraph_xy = varargin{5};     % need this for the error test logic
				set(handles.edit_MaregraphDataFile,'String',[handMir.last_dir 'maregs.dat'],'Enable','on')
			else
				set(handles.edit_MaregraphPosFile,'String','','Enable','off')
				set(handles.edit_MaregraphDataFile,'String','','Enable','off')
				set(handles.pushbutton_MaregraphPosFile,'Enable','off')
				set(handles.pushbutton_MaregraphDataFile,'Enable','off')
			end
		elseif (strcmp(varargin{2},'Tsun2'))
			handles.is_tsun2 = 1;
			handles.outTsu_maregs = [];		% Optional output maregraphs file
			set(handles.edit_BatGrid,'String','In memory array','Enable','off')
			set(handles.pushbutton_BatGrid,'Enable','off')
			set(handles.edit_SourceGrid,'String','Not used in tsun2','Enable','off')
			set(handles.pushbutton_SourceGrid,'Enable','off')
			set(handles.checkbox_momentum,'Vis','off')
			set(handles.checkbox_velocity,'Vis','off')
			handles.BatGridInMemory = 1;
			set(handles.checkbox_WantMaregraphs,'Value',1)
			set(handles.edit_MaregraphPosFile,'Enable','on')
			set(handles.pushbutton_MaregraphPosFile,'Enable','on')
			set(handles.txtOutMaregs,'String','Option')
			str = sprintf(['Optional file with output maregraphs location.\n'...
        					'The results will be saved in a file whose name is\n' ...
							'obtained by appending _maregHeights.dat to this one']);
			set(handles.edit_MaregraphDataFile,'Tooltip',str)
			set(handles.edit_friction,'Vis','on')
			set(handles.textFriction,'Vis','on')
			set(handles.figure1,'Name','Tsun2 options')
		end
	end

	if (handMir.last_dir(end) == filesep)
		set(handles.edit_gridNameStem,'String',[handMir.last_dir 'tsu_time_'],'Enable','off')
	else
		set(handles.edit_gridNameStem,'String',[handMir.last_dir filesep 'tsu_time_'],'Enable','off')
	end

	% Choose default command line output for swan_options
	handles.output = hObject;
	guidata(hObject, handles);

	set(hObject,'Visible','on');
	% UIWAIT makes swan_options wait for user response (see UIRESUME)
	uiwait(handles.figure1);

	handles = guidata(hObject);
	varargout{1} = handles.output;
	delete(handles.figure1);

% -----------------------------------------------------------------------------------------
function edit_BatGrid_Callback(hObject, eventdata, handles)
	fname = get(hObject,'String');
	if ~isempty(fname)
		[handles,X,Y,handles.Z_bat,handles.head_bat] = read_gmt_type_grids(handles,fname);
		if (isempty(X)),    set(hObject,'String',''),	return,		end
		guidata(hObject,handles)
	end

% -----------------------------------------------------------------------------------------
function pushbutton_BatGrid_Callback(hObject, eventdata, handles)
	str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select grid','get');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	[handles,X,Y,handles.Z_bat,handles.head_bat] = read_gmt_type_grids(handles,fname);
	if (isempty(X)),    return,		end
	set(handles.edit_BatGrid,'String',fname)
	guidata(hObject,handles)

% -----------------------------------------------------------------------------------------
function edit_SourceGrid_Callback(hObject, eventdata, handles)
	fname = get(hObject,'String');
	if ~isempty(fname)
		[handles,X,Y,handles.Z_src,handles.head_src] = read_gmt_type_grids(handles,fname);
		if (isempty(X)),    set(hObject,'String',''),	return,		end
		guidata(hObject,handles)
	end

% -----------------------------------------------------------------------------------------
function pushbutton_SourceGrid_Callback(hObject, eventdata, handles)
	str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select grid','get');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	[handles,X,Y,handles.Z_src,handles.head_src] = read_gmt_type_grids(handles,fname);
	if (isempty(X)),    return,		end
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

	cfig_handles = guidata(handles.hCallingFig);      % get handles of the calling fig
	last_dir = cfig_handles.last_dir;
	home = cfig_handles.home_dir;

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
		if (~isempty(tline) && tline(1) ~= '#')      % Jump comment lines
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
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 2)
		set(hObject,'String','1010')
	else
		set(hObject,'String',fix(xx) ) % Insure that it is an integer
	end

% -----------------------------------------------------------------------------------------
function edit_Jump_initialTime_Callback(hObject, eventdata, handles)
% Remove this number of seconds from the beguining of the maregraphs (if they were loaded)
jmp = str2double(get(hObject,'String'));
if (isnan(jmp) || jmp < 0 || isempty(handles.maregraph_xy))
	set(hObject,'String','0');
	handles.n_jump = 0;
else
    [m,n] = size(handles.maregraph_xy);
    dt1 = diff(handles.maregraph_xy(:,1));    dt2 = diff(dt1);
    if any(dt2 ~= 0)				% First column doesn't have the time
    else							% First column has the time.
        dt = dt1(1);				% This is the time increment
        n_jmp = fix(jmp / dt);		% Number of records to jump
        if ((m - n_jmp) < 10)		% Stupid choice for jump time
            set(hObject,'String','0')
			n_jmp = 0;
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
	if isempty(fname),  return,		end
	pushbutton_MaregraphPosFile_Callback(handles.pushbutton_MaregraphPosFile, [], handles, fname)
	
% -----------------------------------------------------------------------------------------
function pushbutton_MaregraphPosFile_Callback(hObject, eventdata, handles, opt)

	if (nargin == 3)
		cfig_handles = guidata(handles.hCallingFig);      % get handles of the calling fig
		last_dir = cfig_handles.last_dir;
		home = cfig_handles.home_dir;
	
		if (~isempty(last_dir)),    cd(last_dir);   end
		[FileName,PathName] = uigetfile({'*.dat;*.DAT;*.xy', 'Maregraph location (*.dat,*.DAT,*.xy)';'*.*', 'All Files (*.*)'},'Select Maregraphs position');
		pause(0.01);
		if (~isempty(last_dir)),    cd(home);   end
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	else
		fname = opt;
	end

    set(handles.figure1,'pointer','watch');
	[handles, msg] = getMaregsPos(handles, fname);
	if (~isempty(msg))
		errordlg(msg,'Error')
	    set(handles.figure1,'pointer','arrow');
		return
	end
    set(handles.figure1,'pointer','arrow');
	
	set(handles.edit_Number_of_cycles,'String',size(handles.maregraph_xy,1))
	set(handles.edit_MaregraphPosFile,'String',fname)
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function [handles, msg] = getMaregsPos(handles, fname)
	msg = [];
	hFig = gcf;
	[bin,n_column,multi_seg,n_headers] = guess_file(fname);
	% If msgbox exist we have to move it from behind the main window. So get it's handle
	hMsgFig = gcf;
	if (hFig ~= hMsgFig),		figure(hMsgFig),	end   % If msgbox exists, bring it forward
	% If error in reading file
	if (isempty(bin) && isempty(n_column) && isempty(multi_seg))
		msg = ['Error reading file ' fname];		return
	end
	if multi_seg ~= 0   % multisegments are not spported
		msg = 'Multisegment files are yet not supported.';   return
	end
	if (bin ~= 0)
		msg = 'Sorry, reading binary files is not programed';   return
	end
    if (n_column < 2)
		msg = 'File error. Your file doesn''t have at least 2 columns';		return
    end
    handles.maregraph_xy = read_xy(fname,n_column,n_headers);
    set(handles.edit_Jump_initialTime,'String','0')     % Reset this anyway
    if (hFig ~= hMsgFig),		figure(hFig);   end     % gain access to the drawing figure
    [nr,nc] = size(handles.maregraph_xy);
    if (nr == 0)
		msg = 'Your file is empty.';   return
    end

% -----------------------------------------------------------------------------------------
function edit_MaregraphDataFile_Callback(hObject, eventdata, handles)
	if (handles.is_tsun2)
		fname = get(hObject,'String');
		if isempty(fname),  return,		end
		pushbutton_MaregraphDataFile_Callback(handles.pushbutton_MaregraphDataFile, [], handles, fname)
	end

% -----------------------------------------------------------------------------------------
function pushbutton_MaregraphDataFile_Callback(hObject, eventdata, handles, opt)

	if (nargin == 3)
		cfig_handles = guidata(handles.hCallingFig);	% get handles of the calling fig
		last_dir = cfig_handles.last_dir;
		home = cfig_handles.home_dir;
	
		if (~isempty(last_dir)),    cd(last_dir);   end
		[FileName,PathName] = uigetfile({'*.dat;*.DAT;*.xy', 'Maregraph data file (*.dat,*.DAT,*.xy)';'*.*', 'All Files (*.*)'},'Select Maregraph');
		pause(0.01);
		if (~isempty(last_dir)),    cd(home);   end
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
		
		if (~handles.is_tsun2)		% We are done with this case
			set(handles.edit_MaregraphDataFile,'String',fname)
			return
		end
	else
		fname = opt;
	end

	% If we get here we are in the tsun2 oputput maregraphs option. We must read that file
	bak = handles.maregraph_xy;		% To reuse the code from getMaregsPos() we must backup this
	[handles, msg] = getMaregsPos(handles, fname);
	if (~isempty(msg))
		errordlg(msg,'Error')
		handles.outTsu_maregs = [];		handles.maregraph_xy = bak;			% Get back the original data
		guidata(handles.figure1,handles)
		return
	end
	% OK, we are still playing
	handles.outTsu_maregs = handles.maregraph_xy(:,1:2);	% We only want the locations
	handles.maregraph_xy = bak;			% Get back the original data
	set(handles.edit_MaregraphDataFile,'String',fname)
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_friction_Callback(hObject, eventdata, handles)
	xx = str2double(get(handles.edit_friction,'String'));
	if (isnan(xx) || xx < 0)
		set(handles.edit_friction,'String','0.025')
	end

% -----------------------------------------------------------------------------------------
function checkbox_speeds_Callback(hObject, eventdata, handles)
	% If velocities or momentums, activate the output grids option
	if (get(hObject,'Value'))
		set(handles.checkbox_OutputGrids,'Value',1)
		set(handles.edit_gridNameStem,'Enable','on')
	end
	
% -----------------------------------------------------------------------------------------
function checkbox_OutputGrids_Callback(hObject, eventdata, handles)
	if (get(hObject,'Value'))
		set(handles.edit_gridNameStem,'Enable','on')
	else
		set(handles.edit_gridNameStem,'String','tsu_time_')
		set(handles.edit_gridNameStem,'Enable','off')
	end

% -----------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
% Do some error checking

if (handles.is_tsun2 == 1)    % Tests for the case of tsun2 (and return)
    error = check_errors_tsun2(handles);
    if (error)
        handles.output = [];
        guidata(handles.figure1,handles);
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
    % Get the friction coefficient
	handles.output.opt_F = ['-F' get(handles.edit_friction,'String')];
	
	if (~isempty(handles.outTsu_maregs))
		handles.output.opt_O.xy = handles.outTsu_maregs;
		handles.output.opt_O.name = get(handles.edit_MaregraphDataFile,'String');
	end
	
    guidata(handles.figure1,handles);
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
	
	% Want velocity/momentum?
	if (get(handles.checkbox_velocity, 'Value'))
		handles.output.opt_s = '-s';
	end
	if (get(handles.checkbox_momentum, 'Value'))
		handles.output.opt_S = '-S';
	end
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
        errordlg('Where are your maregraphs? On the Moon? (In case you forgot, it has no oceans)','Error')
        error = 1;
    end
    if (isempty(get(handles.edit_MaregraphDataFile,'String')))
        errordlg('You need to tell me the file name where I''ll write the maregraphs water hight','Error')
        error = 1;
    end
end

% Cheeck that at lest one operation was selected
if (~get(handles.checkbox_WantMaregraphs,'Value') && ~get(handles.checkbox_MakeMovie,'Value') && ...
        ~get(handles.checkbox_OutputGrids,'Value'))
    errordlg('You need to select at least one operation','Error')
    error = 1;
end

if (error),		return,		end

% If we reach here, we can now check if grids are compatible. However, if grids where already in
% Mirone's memory, it is there that this test should have been made
if (handles.BothGridsInMemory == 0 || handles.BatGridInMemory == 1)
	difes = handles.head_src(1:4) - handles.head_bat(1:4);
	if (any(abs(difes) > small))
		msg{1} = 'Bathymetry & Source grids do not cover the same region';
		msg{2} = ['x_min diff = ' num2str(difes(1))];		msg{3} = ['x_max diff = ' num2str(difes(2))];
		msg{4} = ['y_min diff = ' num2str(difes(3))];		msg{5} = ['y_max diff = ' num2str(difes(4))];
		errordlg(msg','Error')
		error = 1;
	end
	if ( abs(handles.head_src(8) - handles.head_bat(8)) > small || ...
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
	if (~get(handles.checkbox_MakeMovie,'Value') && ~get(handles.checkbox_OutputGrids,'Value') ...
			&& isempty(handles.outTsu_maregs))
		errordlg('You need to select at least one the three operation: grids, maregs out or movie','Error')
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
function figure1_CloseRequestFcn(hObject, eventdata)
    handles = guidata(hObject);
	if isequal(get(hObject, 'waitstatus'), 'waiting')
		% The GUI is still in UIWAIT, us UIRESUME
		handles.output = [];        % User gave up, return nothing
		guidata(hObject, handles);    uiresume(hObject);
	else
		% The GUI is no longer waiting, just close it
		handles.output = [];        % User gave up, return nothing
		guidata(hObject, handles);	delete(hObject)
	end

% -----------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
        handles = guidata(hObject);
		handles.output = [];    % User said no by hitting escape
		guidata(hObject, handles);    uiresume(hObject);
	end

% -----------------------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function swan_options_LayoutFcn(h1);

set(h1, 'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Swan options',...
'NumberTitle','off',...
'Position',[520 401 351 370],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_BatGrid_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 341 271 21],...
'Style','edit','Tag','edit_BatGrid');

uicontrol('Parent',h1, 'Position',[321 339 23 23],...
'Callback',{@swan_options_uicallback,h1,'pushbutton_BatGrid_Callback'},...
'Tag','pushbutton_BatGrid');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_SourceGrid_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 311 271 21],...
'Style','edit','Tag','edit_SourceGrid');

uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'checkbox_OutputGrids_Callback'},...
'Position',[10 112 77 15],...
'String','Output grids','Style','checkbox',...
'Tag','checkbox_OutputGrids');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[95 108 225 21],...
'Style','edit',...
'TooltipString','Grids are numbered after this name stem',...
'Tag','edit_gridNameStem');

uicontrol('Parent',h1, 'Position',[10 235 82 15],...
'Callback',{@swan_options_uicallback,h1,'radiobutton_SurfaceLevel_Callback'},...
'String','Surface level','Style','radiobutton',...
'Value',1,'Tag','radiobutton_SurfaceLevel');

uicontrol('Parent',h1, 'Position',[131 235 82 15],...
'Callback',{@swan_options_uicallback,h1,'radiobutton_MaxWater_Callback'},...
'String','Max water','Style','radiobutton',...
'Tag','radiobutton_MaxWater');

uicontrol('Parent',h1, 'Position',[249 235 82 15],...
'Callback',{@swan_options_uicallback,h1,'radiobutton_TotalWater_Callback'},...
'String','Total water','Style','radiobutton',...
'Tag','radiobutton_TotalWater');

uicontrol('Parent',h1, 'Position',[131 215 82 15],...
'Callback',{@swan_options_uicallback,h1,'checkbox_speeds_Callback'},...
'String','Velocity','Style','checkbox',...
'TooltipString','Write velocity grids (u and v with sufixes _U, _V) ',...
'Tag','checkbox_velocity');

uicontrol('Parent',h1, 'Position',[249 215 82 15],...
'Callback',{@swan_options_uicallback,h1,'checkbox_speeds_Callback'},...
'String','Momentum','Style','checkbox',...
'TooltipString','Write momentum grids (with sufixes _Uh, _Vh) ',...
'Tag','checkbox_momentum');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_Number_of_cycles_Callback'},...
'Position',[270 78 51 21],...
'String','1010','Style','edit',...
'TooltipString','Use this number of cycles from the Maregraph file',...
'Tag','edit_Number_of_cycles');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_MaregraphPosFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 179 271 21],...
'Style','edit',...
'TooltipString','Name of the file with maregraph locations',...
'Tag','edit_MaregraphPosFile');

uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_MaregraphPosFile_Callback'},...
'Position',[321 178 23 23],...
'Tag','pushbutton_MaregraphPosFile');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_MaregraphDataFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 149 271 21],...
'Style','edit',...
'TooltipString','Name of the file that will contain the maregraphs water height',...
'Tag','edit_MaregraphDataFile');

uicontrol('Parent',h1, 'Position',[209 82 60 15],...
'String','Nº of cycles','Style','text');

uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[10 345 40 15],...
'String','Bat','Style','text','Tag','text2');

uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[10 314 40 15],...
'String','Source','Style','text');

uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_MaregraphDataFile_Callback'},...
'Position',[321 148 23 23],...
'Tag','pushbutton_MaregraphDataFile');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_SwanParams_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 270 271 21],...
'Style','edit',...
'TooltipString','Either Swan or Tsun2 need a parameter file.',...
'Tag','edit_SwanParams');

uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[10 273 40 15],...
'String','Params','Style','text');

uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_OK_Callback'},...
'FontWeight','bold',...
'Position',[195 8 66 23],...
'String','OK','Tag','pushbutton_OK');

uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_Cancel_Callback'},...
'FontWeight','bold',...
'Position',[275 8 66 23],...
'String','Cancel',...
'Tag','pushbutton_Cancel');

uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'checkbox_WantMaregraphs_Callback'},...
'Position',[10 202 76 15],...
'String','Maregraphs','Style','checkbox',...
'TooltipString','Check this if you want to compute water height at maregraphs',...
'Tag','checkbox_WantMaregraphs');

uicontrol('Parent',h1,...
'Callback',{@swan_options_uicallback,h1,'pushbutton_SourceGrid_Callback'},...
'Position',[321 310 23 23],...
'Tag','pushbutton_SourceGrid');

uicontrol('Parent',h1, 'Position',[321 269 23 23],...
'Callback',{@swan_options_uicallback,h1,'pushbutton_ParamsFile_Callback'},...
'Tag','pushbutton_ParamsFile');

uicontrol('Parent',h1, 'Position',[10 81 85 15],...
'String','Make a movie','Style','checkbox',...
'Tag','checkbox_MakeMovie');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_Jump_initialTime_Callback'},...
'Position',[270 48 51 21],...
'String','0','Style','edit',...
'TooltipString','Jump these number of seconds from the beguining of Maregraphs file',...
'Tag','edit_Jump_initialTime');

uicontrol('Parent',h1,'HorizontalAlignment','right',...
'Position',[212 52 55 15],...
'String','Jump initial','Style','text');

uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[13 181 30 16],...
'String','In file','Style','text');

uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[12 152 35 16],...
'String','Out file','Style','text','Tag','txtOutMaregs');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@swan_options_uicallback,h1,'edit_friction_Callback'},...
'Position',[10 48 51 21],...
'String','0.025','Style','edit',...
'TooltipString','Manning''s Friction coefficient',...
'Visible','off',...
'Tag','edit_friction');

uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[65 50 80 16],'Visible','off',...
'String','Friction','Style','text','Tag','textFriction');

function swan_options_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

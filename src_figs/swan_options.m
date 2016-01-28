function varargout = swan_options(varargin)
% Create an options structure to be used later on swan or tsun2 mexs
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
%   output.opt_J = '-J<time_jump>';
%   output.opt_M = '-M';
%   output.opt_N
%   output.opt_S = '-S' (swan only - write momentums)
%   output.opt_s = '-s' (swan only - write velocities)
%   --- Or empty if user gave up

%	Copyright (c) 2004-2013 by J. Luis
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

% $Id: swan_options.m 3891 2013-03-01 14:25:29Z j $

	if (isempty(varargin))
		errordlg('Bad call to swan_option. Please warn me (the author)','Error')
		return
	end

	hObject = figure('Tag','figure1','Visible','off');
	swan_options_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')
	
	handles.head_bat = [];
	handles.head_src = [];
	handles.maregraph_xy = [];
	handles.BothGridsInMemory = 0;  % Flag to signal that bat & deform arrays are already in memory
	handles.BatGridInMemory = 0;    % Flag to signal that bat array is already in memory
	handles.MaregraphInMemory = 1;  % Flag to signal that the maregraphs locations are already in memory
	handles.got_params = 0;         % Flag that signals a (possibly correct) params file
	handles.is_tsun2 = 0;           % Flag to signal that options are for use in tsun2 code
	handles.n_jump = 0;             % If > 0, the maregraphs array will start at n_jump + 1. OR netCDF files start
	handles.polar = 0;
	handles.dt = [];
	handles.grn = 10;
	handles.cumint = 1;
	
	handMir = varargin{1};
	handles.hCallingFig = handMir.figure1;
	handles.home_dir = handMir.home_dir;
	handles.work_dir = handMir.work_dir;
	handles.last_dir = handMir.last_dir;

	% Import icons
	f_data = [handles.home_dir filesep 'data' filesep];
	load([f_data 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_BatGrid,'CData',Mfopen_ico)
	set(handles.push_SourceGrid,'CData',Mfopen_ico)
	set(handles.push_paramsFile,'CData',Mfopen_ico)
	set(handles.push_MaregraphPosFile,'CData',Mfopen_ico)
	set(handles.push_MaregraphDataFile,'CData',Mfopen_ico)
	clear Mfopen_ico;

	if (numel(varargin) > 1)
		if (strcmp(varargin{2},'bat_and_deform_with_maregs') || strcmp(varargin{2},'bat_and_deform'))
			set(handles.edit_BatGrid,'String',' In memory array','Enable','off')
			set(handles.push_BatGrid,'Enable','off')
			set(handles.edit_SourceGrid,'String',' In memory array','Enable','off')
			set(handles.push_SourceGrid,'Enable','off')
			handles.BothGridsInMemory = 1;
			if (strcmp(varargin{2},'bat_and_deform_with_maregs'))
				set(handles.check_wantMaregs,'Value',1)
				set(handles.edit_MaregraphPosFile,'String',' In memory array','Enable','off')
				set(handles.push_MaregraphPosFile,'Enable','off')
				handles.MaregraphInMemory = 1;
				handles.maregraph_xy = varargin{3};     % need this for the error test logic
				set(handles.edit_MaregraphDataFile,'String',[handMir.last_dir filesep 'maregs.dat'],'Enable','on')
			else
				set(handles.edit_MaregraphPosFile,'String','','Enable','off')
				set(handles.edit_MaregraphDataFile,'String','','Enable','off')
				set(handles.push_MaregraphPosFile,'Enable','off')
				set(handles.push_MaregraphDataFile,'Enable','off')
			end
		elseif (strcmp(varargin{2},'bat_with_maregs') || strcmp(varargin{2},'bat_only'))
			handles.Z_bat = varargin{3};
			handles.head_bat = varargin{4};
			set(handles.edit_BatGrid,'String',' In memory array','Enable','off')
			set(handles.push_BatGrid,'Enable','off')
			handles.BatGridInMemory = 1;
			if (strcmp(varargin{2},'bat_with_maregs'))
				set(handles.check_wantMaregs,'Value',1)
				set(handles.edit_MaregraphPosFile,'String','In memory array','Enable','off')
				set(handles.push_MaregraphPosFile,'Enable','off')
				handles.MaregraphInMemory = 1;
				handles.maregraph_xy = varargin{5};     % need this for the error test logic
				set(handles.edit_MaregraphDataFile,'String',[handMir.last_dir filesep 'maregs.dat'],'Enable','on')
			else
				set(handles.edit_MaregraphPosFile,'String','','Enable','off')
				set(handles.edit_MaregraphDataFile,'String','','Enable','off')
				set(handles.push_MaregraphPosFile,'Enable','off')
				set(handles.push_MaregraphDataFile,'Enable','off')
			end
		elseif (strcmp(varargin{2},'Tsun2'))
			handles.is_tsun2 = 1;
			handles.outTsu_maregs = [];		% Optional output maregraphs file
			set(handles.edit_BatGrid,'String',' In memory array','Enable','off')
			set(handles.push_BatGrid,'Enable','off')
			set(handles.edit_SourceGrid,'String','Not used in tsun2','Enable','off')
			set(handles.push_SourceGrid,'Enable','off')
			set(handles.check_momentum,'Vis','off')
			set(handles.check_velocity,'Vis','off')
			set(handles.radio_outGrids,'Val',1)
			set([handles.radio_anuga handles.radio_most],'Enable','off')
			set([handles.text_minimalist handles.text_dt handles.text_grn handles.text_cumint handles.edit_dt ...
					handles.edit_grn handles.edit_cumint handles.radio_cartesian handles.radio_geog],'Enable','off')
			set(handles.text_parFile,'String','Read parameters from a tsun2.par file')
			set(handles.edit_SwanParams,'Tooltip','Enter a Tsun2 parameter file name')
			handles.BatGridInMemory = 1;
			set(handles.check_wantMaregs,'Value',1)
			set(handles.edit_MaregraphPosFile,'Enable','on')
			set(handles.push_MaregraphPosFile,'Enable','on')
			set(handles.txtOutMaregs,'String','Option')
			str = sprintf(['Optional file with output maregraphs location.\n'...
        					'The results will be saved in a file whose name is\n' ...
							'obtained by appending _maregHeights.dat to this one']);
			set(handles.edit_MaregraphDataFile,'Tooltip',str)
			set(handles.edit_jumpInitial,'Tooltip','Jump these number of seconds from the begining of Maregraphs file')
			set(handles.edit_friction,'Vis','on')
			set(handles.textFriction,'Vis','on')
			set(handles.figure1,'Name','Tsun2 options')
		end
	end

	if ( ~isempty(handles.head_bat) && ~handles.is_tsun2 )
		geog = guessGeog(handles.head_bat(1:4));
		set(handles.radio_cartesian,'Value', ~geog)
		set(handles.radio_geog,'Value', geog)
		% Guess a valid dt
		dx = min(handles.head_bat(8), handles.head_bat(9));		% Grid cells may not be squares
		if (geog)
			dx = dx * 111000;
			handles.polar = 1;			% Southern hemisphere -1
		end
		dt = dx / sqrt(abs(handles.head_bat(5)) * 9.8) * 0.5;
		if (dt > 0)
			set(handles.edit_dt, 'String',sprintf('%.2f',dt))
			handles.dt = dt;
		end
	end

	handles.outPato = handMir.last_dir;
	if (handMir.last_dir(end) ~= filesep),		handles.outPato = [handMir.last_dir filesep];	end

	% --------- Move the uicontrols that are still 'somewhere' to their correct positions
	hhs = findobj(hObject, 'UserData', 'parameter');
	hTab = findobj(handles.tab_group,'UserData','parameter'); 		% remove the "Params" tab push button from the hhs list
    hhs = setdiff(hhs, hTab);
	for (k = 1:numel(hhs))
		hhsPos = get(hhs(k), 'Pos');		hhsPos = hhsPos - [350 0 0 0];
		set(hhs(k), 'Pos', hhsPos)
	end

	%----------- Give a Pro look (3D) to the frame box  ---------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) -------------------------------------------------------

	% ------------------ TABPANEL SECTION ----------------------------------------
	% This is the tag that all tab push buttons share.  If you have multiple
	% sets of tab push buttons, each group should have unique tag.
	group_name = 'tab_group';
	
	% This is a list of the UserData values used to link tab push buttons and
	% the components on their linked panels.  To add a new tab panel to the group
	%  Add the button using GUIDE
	%  Assign the Tag based on the group name - in this case tab_group
	%  Give the UserData a unique name - e.g. another_tab_panel
	%  Add components to GUIDE for the new panel
	%  Give the new components the same UserData as the tab button
	%  Add the new UserData name to the below cell array
	panel_names = {'main','parameter'};
	
	% tabpanelfcn('makegroups',...) adds new fields to the handles structure,
	% one for each panel name and another called 'group_name_all'.  These fields
	% are used by the tabpanefcn when tab_group_handler is called.
	handles = tabpanelfcn('make_groups',group_name, panel_names, handles, 1);
	% ------------------------------------------------------------------------------

	guidata(hObject, handles);

	set(hObject,'Visible','on');
	% UIWAIT makes swan_options wait for user response (see UIRESUME)
	uiwait(handles.figure1);

	handles = guidata(hObject);
	varargout{1} = handles.output;
	delete(handles.figure1);

% -------------------------------------------------------------------------------------
function tab_group_ButtonDownFcn(hObject, handles)
% Call the tab_group_handler.  This updates visiblity of components as needed to
% hide the components from the previous tab and show components on this tab.
% This also updates the last_tab field in the handles structure to keep track
% of which panel was hidden.
    handles = tabpanelfcn('tab_group_handler',hObject, handles, get(hObject, 'Tag'));
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_BatGrid_CB(hObject, handles)
	fname = get(hObject,'String');
	if ~isempty(fname)
		push_BatGrid_CB(handles.push_BatGrid, handles, fname)
	end

% -----------------------------------------------------------------------------------------
function push_BatGrid_CB(hObject, handles, opt)
	if (nargin == 3)
		str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
		[FileName,PathName] = put_or_get_file(handles,str1,'Select grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	else
		fname = opt;
	end

	[handles,X,Y,handles.Z_bat,handles.head_bat] = read_gmt_type_grids(handles,fname);
	if (isempty(X)),    return,		end
	set(handles.edit_BatGrid,'String',fname)
	geog = guessGeog(handles.head_bat(1:4));
	set(handles.radio_cartesian,'Value', ~double(geog))
	set(handles.radio_geog,'Value', double(geog))
	% Guess a valid dt
	dx = handles.head_bat(9);
	if (geog)
		dx = dx * 111000;
		handles.polar = 1;
	end
	dt = dx / sqrt(abs(handles.head_bat(5)) * 9.8) * 0.5;
	if (dt > 0)
		set(handles.edit_dt, 'String',sprintf('%.2f',dt))
		handles.dt = dt;
	end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_SourceGrid_CB(hObject, handles)
	fname = get(hObject,'String');
	if (isempty(fname)),	return,		end
	push_SourceGrid_CB(handles.push_SourceGrid, handles, fname)

% -----------------------------------------------------------------------------------------
function push_SourceGrid_CB(hObject, handles)
	str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select grid','get');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	[handles,X,Y,handles.Z_src,handles.head_src] = read_gmt_type_grids(handles,fname);
	if (isempty(X)),    return,		end
	set(handles.edit_SourceGrid,'String',fname)
	geog = guessGeog(handles.head_src(1:4));
	set(handles.radio_cartesian,'Value', ~double(geog))
	set(handles.radio_geog,'Value', double(geog))
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_SwanParams_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),  return,		end
	push_paramsFile_CB(handles.push_paramsFile, handles, fname)

% -----------------------------------------------------------------------------------------
function push_paramsFile_CB(hObject, handles, opt)

	if ( nargin == 2 )
		[FileName,PathName] = put_or_get_file(handles,{'*.par', 'Params file (*.par)';'*.*', 'All Files (*.*)'},'Select parameter file','get');
		if ( isequal(FileName,0) ),		return,		end
		fname = [PathName FileName];
	else
		fname = opt;
	end

	fid = fopen(fname,'r');
	if (fid < 0)
		errordlg('Error opening params file (it probably doesn''t exist)','Error');
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
	fclose(fid);
	handles.got_params = 1;     % Lets just hope that the info in it is correct
	set(handles.edit_SwanParams,'String',fname)
	if ( ~handles.is_tsun2 )
		set(handles.edit_dt, 'String',handles.params(1))
		set(handles.edit_grn,'String',handles.params(2))
		set(handles.edit_cumint,'String',handles.params(9))
		if ( abs(handles.params(7)) == 1 )
			set(handles.radio_geog, 'Val', 1);			set(handles.radio_cartesian, 'Val', 0)
		else
			set(handles.radio_geog, 'Val', 0);			set(handles.radio_cartesian, 'Val', 1)
		end
	end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function radio_surfLevel_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.radio_maxWater,'Value',0)
		set(handles.radio_totalWater,'Value',0)
	else
		set(hObject,'Value',1)
		set(handles.radio_maxWater,'Value',0)
		set(handles.radio_totalWater,'Value',0)
	end

% -----------------------------------------------------------------------------------------
function radio_maxWater_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.radio_surfLevel,'Value',0)
		set(handles.radio_totalWater,'Value',0)
	else
		set(hObject,'Value',1)
		set(handles.radio_surfLevel,'Value',0)
		set(handles.radio_totalWater,'Value',0)
	end

% -----------------------------------------------------------------------------------------
function radio_totalWater_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.radio_surfLevel,'Value',0)
		set(handles.radio_maxWater,'Value',0)
	else
		set(hObject,'Value',1)
		set(handles.radio_surfLevel,'Value',0)
		set(handles.radio_maxWater,'Value',0)
	end

% -----------------------------------------------------------------------------------------
function edit_Number_of_cycles_CB(hObject, handles)
	% The OK button will get the number of cycles from here. So make sure it gives a possible number
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 2)
		set(hObject,'String','1010')
	else
		set(hObject,'String',fix(xx) ) % Insure that it is an integer
	end

% -----------------------------------------------------------------------------------------
function edit_jumpInitial_CB(hObject, handles)
% Remove this number of seconds from the beguining of the maregraphs (if they were loaded)
	jmp = str2double(get(hObject,'String'));
	if (isnan(jmp) || jmp < 0)
		set(hObject,'String','0');		handles.n_jump = 0;
		guidata(handles.figure1,handles)
		return
	end

	if ( ~isempty(handles.maregraph_xy) )
		m = size(handles.maregraph_xy,1);
     	dt1 = diff(handles.maregraph_xy(:,1));    dt2 = diff(dt1);
		if any(dt2 ~= 0)				% First column doesn't have the time
		else							% First column has the time.
			dt = dt1(1);				% This is the time increment
			n_jmp = fix(jmp / dt);		% Number of records to jump
			if ((m - n_jmp) < 10)		% Stupid choice for jump time
				set(hObject,'String','0')
				n_jmp = 0;
			end
			set(handles.edit_Number_of_cycles,'String',num2str(m-n_jmp))    % Uppdate max pssible
			handles.n_jump = n_jmp;
		end
	else
		handles.n_jump = jmp;			% To be used in ANUGA or MOST outputs. It's a jump time, not cycles
	end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function check_wantMaregs_CB(hObject, handles)
	% In tsun2 maregraphs are mandatory (for reading), so never let the user uncheck
	if (handles.is_tsun2 == 1),     set(hObject,'Value',1);     return;     end

	if (get(hObject,'Value'))
		set(handles.edit_MaregraphPosFile,'Enable','on')
		set(handles.edit_MaregraphDataFile,'Enable','on')
		set(handles.push_MaregraphPosFile,'Enable','on')
		set(handles.push_MaregraphDataFile,'Enable','on')
	else
		set(handles.edit_MaregraphPosFile,'String','','Enable','off')
		set(handles.edit_MaregraphDataFile,'String','','Enable','off')
		set(handles.push_MaregraphPosFile,'Enable','off')
		set(handles.push_MaregraphDataFile,'Enable','off')
	end

% -----------------------------------------------------------------------------------------
function edit_MaregraphPosFile_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),  return,		end
	push_MaregraphPosFile_CB(handles.push_MaregraphPosFile, handles, fname)
	
% -----------------------------------------------------------------------------------------
function push_MaregraphPosFile_CB(hObject, handles, opt)

	if (nargin == 3)
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.dat;*.DAT;*.xy', 'Maregraph location (*.dat,*.DAT,*.xy)';'*.*', 'All Files (*.*)'},'Select Maregraphs position','get');
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
	if (isa(bin,'struct') || bin ~= 0)
		msg = 'Sorry, reading binary files is not programed';   return
	end
    if (n_column < 2)
		msg = 'File error. Your file doesn''t have at least 2 columns';		return
    end
    handles.maregraph_xy = read_xy(fname,n_column,n_headers);
    set(handles.edit_jumpInitial,'String','0')     % Reset this anyway
    if (hFig ~= hMsgFig),		figure(hFig);   end     % gain access to the drawing figure
    nr = size(handles.maregraph_xy,1);
    if (nr == 0)
		msg = 'Your file is empty.';   return
    end

% -----------------------------------------------------------------------------------------
function edit_MaregraphDataFile_CB(hObject, handles)
	if (handles.is_tsun2)
		fname = get(hObject,'String');
		if isempty(fname),  return,		end
		push_MaregraphDataFile_CB(handles.push_MaregraphDataFile, handles, fname)
	end

% -----------------------------------------------------------------------------------------
function push_MaregraphDataFile_CB(hObject, handles, opt)

	if (nargin == 3)
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.dat;*.DAT;*.xy', 'Maregraph data file (*.dat,*.DAT,*.xy)';'*.*', 'All Files (*.*)'},'Select Maregraph','get');
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
function edit_friction_CB(hObject, handles)
	xx = str2double(get(handles.edit_friction,'String'));
	if (isnan(xx) || xx < 0)
		set(handles.edit_friction,'String','0.025')
	end

% -----------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% Do some error checking

if (handles.is_tsun2 == 1)    % Tests for the case of tsun2 (and return)
    error = check_errors_tsun2(handles);
    if (error)
        handles.output = [];
        guidata(handles.figure1,handles);
        return
    end
    handles.output.params_file_name = get(handles.edit_SwanParams,'String');	% desenrasque
    % Reading maregraphs was set?
    if (~isempty(handles.maregraph_xy))
        if (handles.n_jump)     % We have a time jump
            handles.maregraph_xy = handles.maregraph_xy(handles.n_jump+1:end,:);
        end
        handles.output.maregraph_xy = handles.maregraph_xy;
    end
    if (get(handles.radio_maxWater,'Value')),		handles.output.opt_M = '-M';	end
    if (get(handles.radio_totalWater,'Value')),		handles.output.opt_D = '-D';	end
	
    % Do we want a collection of intermediary tsunami grid steps?
    if (get(handles.radio_outGrids,'Value')),	handles.output.opt_G = ['-G' get(handles.edit_gridNameStem,'String')];	end
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
    guidata(handles.figure1,handles);
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

% This one is either read in this program or filled with default values
if (handles.got_params)
	handles.output.params = handles.params;
else		% Fill the rest os the swan.par parameters with default values
	handles.output.params = zeros(1,22);
	handles.output.params(1) = handles.dt;
	handles.output.params(2) = handles.grn;
	handles.output.params(7) = handles.polar;
	handles.output.params(8) = 1;
	handles.output.params(9) = handles.cumint;
	handles.output.params(18:21) = 2;
	if (handles.polar && handles.output.grid_head_src(4) < 0)	% We still must find if = 1 (north) or = -1 (south lats)
		handles.output.params(7) = -1;
	end
end

% See if computation of maregraphs was required
if (~isempty(handles.maregraph_xy))
    handles.output.maregraph_xy = handles.maregraph_xy;
    handles.output.maregraph_data_name = get(handles.edit_MaregraphDataFile,'String');
end

% Check what is going to be computed (surface level OR Max water OR Total water depths - OR ANUGA/MOST)
% Remember, Surface level is the default
if (get(handles.radio_maxWater,'Value')),	handles.output.opt_M = '-M';	end
if (get(handles.radio_totalWater,'Value')),	handles.output.opt_D = '-D';	end    

% Do we want a collection of intermediary tsunami grid steps?
if (get(handles.radio_outGrids,'Value'))
	handles.output.opt_G = ['-G' get(handles.edit_gridNameStem,'String')];
	% Want velocity/momentum?
	if (get(handles.check_velocity, 'Value')),		handles.output.opt_s = '-s';	end
	if (get(handles.check_momentum, 'Value')),		handles.output.opt_S = '-S';	end

elseif (get(handles.radio_anuga,'Value'))
	handles.output.opt_G = ['-A' get(handles.edit_gridNameStem,'String')];

elseif (get(handles.radio_most,'Value'))
	handles.output.opt_G = ['-n' get(handles.edit_gridNameStem,'String')];
end

% Want to Jump initial time on output grids/netCDFs?
if ( handles.n_jump ),	handles.output.opt_J = sprintf('-J%.3f', handles.n_jump);	end

% Get number of cycles
handles.output.opt_N = ['-N' get(handles.edit_Number_of_cycles,'String')];
guidata(handles.figure1,handles)
uiresume(handles.figure1);

% -----------------------------------------------------------------------------------------
function push_Cancel_CB(hObject, handles)
	handles.output = [];			% User gave up, return nothing
	guidata(handles.figure1, handles);		uiresume(handles.figure1);

%---------------------------------------------------------------------------------------------------
function error = check_errors_swan(handles)
error = 0;
small = 1e-6;   % Used in postion comparation. It places the accuracy at sub-meter level
if (handles.BothGridsInMemory == 0)     % If arrays are already in Mirone's memory there is no need to test this
    if (isempty(handles.head_bat) && handles.BatGridInMemory == 0)
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

% Cheeck that at lest one operation was selected
if ( ~(get(handles.check_wantMaregs,'Val') || get(handles.radio_outGrids,'Val') || ...
        get(handles.radio_anuga,'Val') || get(handles.radio_most,'Val')) )
    errordlg('You need to select at least one thing to do','Error')
    error = 1;
end
    
if (get(handles.check_wantMaregs,'Value'))       % That is, requested maregraphs computation
    if (isempty(handles.maregraph_xy) && handles.MaregraphInMemory == 0)
        errordlg('Where are your maregraphs? On the Moon? (In case you forgot, it has no oceans)','Error')
        error = 1;
    end
    if (isempty(get(handles.edit_MaregraphDataFile,'String')))
        errordlg('You need to tell me the file name where I''ll write the maregraphs water hight','Error')
        error = 1;
    end
end

if ( ~handles.got_params && isempty(handles.dt) )        % Error in params file
	errordlg('Must read a parameter file or fill all fields of "Minimalist parameters setting"','Error')
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
	if (~get(handles.radio_outGrids,'Value') && isempty(handles.outTsu_maregs))
		errordlg('You need to select at least one the operation: grids or maregs','Error')
		error = 1;
	end

% -----------------------------------------------------------------------------------------
function radio_outGrids_CB(hObject, handles)
	if ( ~get(hObject,'Value') )
		set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
				handles.check_velocity handles.check_momentum], 'Enable','off')
		set(handles.edit_gridNameStem, 'String','', 'Tooltip','')
		return
	end
	set([handles.radio_anuga handles.radio_most], 'Value',0)
	set(handles.edit_gridNameStem,'String',[handles.outPato 'tsu_time_'])
	set(handles.edit_gridNameStem,'Tooltip','Grids are numbered using this name stem')
	set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
			handles.check_velocity handles.check_momentum], 'Enable','on')
	set(handles.radio_surfLevel, 'Val', 1)

% -----------------------------------------------------------------------------------------
function radio_anuga_CB(hObject, handles)
	if ( ~get(hObject,'Value') )
		set(handles.edit_gridNameStem, 'String','', 'Tooltip','')
		return
	end
	set([handles.radio_outGrids handles.radio_most], 'Value',0)
	set(handles.edit_gridNameStem,'String',[handles.outPato 'swan_tsu.sww'])
	set(handles.edit_gridNameStem,'Tooltip','Output file name')
	set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
			handles.check_velocity handles.check_momentum], 'Val',0)
	set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
			handles.check_velocity handles.check_momentum], 'Enable','off')

% -----------------------------------------------------------------------------------------
function radio_most_CB(hObject, handles)
	if ( ~get(hObject,'Value') )
		set(handles.edit_gridNameStem, 'String','', 'Tooltip','')
		return
	end
	set([handles.radio_outGrids handles.radio_anuga], 'Value',0)
	set(handles.edit_gridNameStem,'String',[handles.outPato 'swan_tsu'])
	set(handles.edit_gridNameStem,'Tooltip','Basename for MOST files (do not give extension)')
	set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
			handles.check_velocity handles.check_momentum], 'Val',0)
	set([handles.radio_surfLevel handles.radio_maxWater handles.radio_totalWater ...
			handles.check_velocity handles.check_momentum], 'Enable','off')
	
% -----------------------------------------------------------------------------------------
function edit_dt_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)
		set(hObject,'String',handles.dt),	return
	end
	% Check that xx is a valid (CFL condition) dt
	dx = handles.head_bat(9);
	geog = get(handles.radio_geog,'Value');
	if (geog),		dx = dx * 111000;		end
	dt = dx / sqrt(abs(handles.head_bat(5)) * 9.8);
	if ( xx > dt )
		warndlg(['The value entered here doesn''t conform to the CFL condition (' ...
				sprintf('%.3f',dt) ') and will probably make the simulation diverge.'],'Warning')
	end
	handles.dt = xx;
	if (handles.got_params),	handles.params(1) = xx;		end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_grn_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)
		set(hObject,'String',handles.grn),	return
	end
	handles.grn = xx;
	if (handles.got_params),	handles.params(2) = xx;		end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_cumint_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)
		set(hObject,'String',handles.cumint),	return
	end
	handles.cumint = xx;
	if (handles.got_params),	handles.params(9) = xx;		end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function radio_cartesian_CB(hObject, handles)
	if ( ~get(hObject,'Value') )
		set(hObject,'Value',1),		return
	end
	set(handles.radio_geog,'Value',0)
	handles.polar = 0;
	if (handles.got_params),	handles.params(7) = 0;		end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function radio_geog_CB(hObject, handles)
	if ( ~get(hObject,'Value') )
		set(hObject,'Value',1),		return
	end
	set(handles.radio_cartesian,'Value',0)
	handles.polar = 1;					% Later on we must decide if = 1 (north hemisphere) or -1 (south)
	if (handles.got_params),	handles.params(7) = 1;		end		% Changed after importing params. Accept it
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
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

% -----------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
        handles = guidata(hObject);
		handles.output = [];    % User said no by hitting escape
		guidata(hObject, handles);    uiresume(hObject);
	end

% -----------------------------------------------------------------------------------------
function xy = read_xy(file,n_col,n_head)
	% build the format string to read the data n_columns
	format = [];    fid = fopen(file,'r');
	for (i=1:n_col),    format = [format '%f '];    end
	% Jump header lines
	for (i = 1:n_head),    tline = fgetl(fid);  end

	todos = fread(fid,'*char');
	xy = sscanf(todos,format,[n_col inf])';
	fclose(fid);

% --------------------------------------------------------------------
function geog = guessGeog(lims)
    % Make a good guess if LIMS are geographic
    geog = double( ( (lims(1) >= -180 && lims(2) <= 180) || (lims(1) >= 0 && lims(2) <= 360) )...
        && (lims(3) >= -90 && lims(4) <= 90) );


% --- Creates and returns a handle to the GUI figure. 
function swan_options_LayoutFcn(h1)

set(h1,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Swan options',...
'NumberTitle','off',...
'Position',[520 359 371 365],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 333 66 21],...
'Enable','inactive',...
'String','Main',...
'ButtonDownFcn',{@swan_options_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','main');

uicontrol('Parent',h1, 'Position',[80 333 80 21],...
'Enable','inactive',...
'String','Parameters',...
'ButtonDownFcn',{@swan_options_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','parameter');

uicontrol('Parent',h1,'Enable','inactive','Position',[10 74 351 261],'Tag','push_bg');

uicontrol('Parent',h1, 'Position',[20 194 331 61],...
'Style','frame',...
'UserData','main');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@swan_options_uiCB,h1,'edit_BatGrid_CB'},...
'HorizontalAlignment','left',...
'Position',[60 304 271 21],...
'Style','edit',...
'Tag','edit_BatGrid',...
'UserData','main');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'push_BatGrid_CB'},...
'Position',[331 302 23 23],...
'Tag','push_BatGrid',...
'UserData','main');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@swan_options_uiCB,h1,'edit_SourceGrid_CB'},...
'HorizontalAlignment','left',...
'Position',[60 274 271 21],...
'Style','edit',...
'Tag','edit_SourceGrid',...
'UserData','main');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'push_SourceGrid_CB'},...
'Position',[331 273 23 23],...
'Tag','push_SourceGrid',...
'UserData','main');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'radio_outGrids_CB'},...
'Position',[31 230 100 15],...
'String','Output grids',...
'Style','radiobutton',...
'Tooltip','Save output  as a series of Surfer format grids',...
'Tag','radio_outGrids',...
'UserData','main');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'radio_anuga_CB'},...
'Position',[141 230 100 15],...
'String','ANUGA .sww',...
'Style','radiobutton',...
'Tooltip','Save output as a single netCDF file using the ANUGA''s .sww format',...
'Tag','radio_anuga',...
'UserData','main');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'radio_most_CB'},...
'Position',[261 230 80 15],...
'String','MOST .nc',...
'Style','radiobutton',...
'Tooltip','Save output as 3 netCDF files as produced by MOST',...
'Tag','radio_most',...
'UserData','main');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[60 204 271 21],...
'Style','edit',...
'Tooltip','Grids are created using this name stem',...
'Tag','edit_gridNameStem',...
'UserData','main');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'radio_surfLevel_CB'},...
'Enable','off',...
'Position',[19 169 100 15],...
'String','Surface level',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_surfLevel',...
'UserData','main');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'radio_maxWater_CB'},...
'Enable','off',...
'Position',[140 169 82 15],...
'String','Max water',...
'Style','radiobutton',...
'Tag','radio_maxWater',...
'UserData','main');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'radio_totalWater_CB'},...
'Enable','off',...
'Position',[258 169 90 15],...
'String','Total water',...
'Style','radiobutton',...
'Tag','radio_totalWater',...
'UserData','main');

uicontrol('Parent',h1, 'Position',[140 149 82 15],...
'Enable','off',...
'String','Velocity',...
'Style','checkbox',...
'Tooltip','Write velocity grids (u and v with sufixes _U, _V)',...
'Tag','check_velocity',...
'UserData','main');

uicontrol('Parent',h1, 'Position',[258 149 90 15],...
'Enable','off',...
'String','Momentum',...
'Style','checkbox',...
'Tooltip','Write momentum grids (with sufixes _Uh, _Vh)',...
'Tag','check_momentum',...
'UserData','main');

uicontrol('Parent',h1, 'Position',[20 134 90 15],...
'Call',{@swan_options_uiCB,h1,'check_wantMaregs_CB'},...
'String','Maregraphs',...
'Style','checkbox',...
'Tooltip','Check this if you want to compute water height at maregraphs',...
'Tag','check_wantMaregs',...
'UserData','main');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@swan_options_uiCB,h1,'edit_MaregraphPosFile_CB'},...
'HorizontalAlignment','left',...
'Position',[60 110 271 21],...
'Style','edit',...
'Tooltip','Name of the file with maregraph locations',...
'Tag','edit_MaregraphPosFile',...
'UserData','main');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'push_MaregraphPosFile_CB'},...
'Position',[331 109 23 23],...
'Tag','push_MaregraphPosFile',...
'UserData','main');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@swan_options_uiCB,h1,'edit_MaregraphDataFile_CB'},...
'HorizontalAlignment','left',...
'Position',[60 80 271 21],...
'Style','edit',...
'Tooltip','Name of the file that will contain the maregraphs water height',...
'Tag','edit_MaregraphDataFile',...
'UserData','main');

uicontrol('Parent',h1, 'Position',[331 80 23 23],...
'Call',{@swan_options_uiCB,h1,'push_MaregraphDataFile_CB'},...
'Tag','push_MaregraphDataFile',...
'UserData','main');

uicontrol('Parent',h1,...
'Position',[120 49 75 15],...
'String','N? of cycles',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@swan_options_uiCB,h1,'edit_Number_of_cycles_CB'},...
'Position',[189 46 51 21],...
'String','1010',...
'Style','edit',...
'Tooltip','Use this number of cycles from the Maregraph file',...
'Tag','edit_Number_of_cycles');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[20 308 40 15],...
'String','Bat',...
'Style','text',...
'UserData','main');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[20 277 40 15],...
'String','Source',...
'Style','text',...
'UserData','main');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@swan_options_uiCB,h1,'edit_SwanParams_CB'},...
'HorizontalAlignment','left',...
'Position',[380 104 271 21],...
'Style','edit',...
'Tooltip','Enter a Swan parameter file name.',...
'Tag','edit_SwanParams',...
'UserData','parameter');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'push_paramsFile_CB'},...
'Position',[651 103 23 23],...
'Tag','push_paramsFile',...
'UserData','parameter');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@swan_options_uiCB,h1,'edit_jumpInitial_CB'},...
'CData',[],...
'Position',[308 46 51 21],...
'String','0',...
'Style','edit',...
'Tooltip','Output is produced only after this modeling time (seconds) has elapsed',...
'Tag','edit_jumpInitial');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[247 49 60 15],...
'String','Jump initial ',...
'Style','text');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[23 113 30 16],...
'String','In file',...
'Style','text',...
'UserData','main');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[22 84 35 16],...
'String','Out file',...
'Style','text',...
'Tag','txtOutMaregs',...
'UserData','main');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@swan_options_uiCB,h1,'edit_friction_CB'},...
'CData',[],...
'Position',[10 21 51 21],...
'String','0.025',...
'Style','edit',...
'Tooltip','Manning''''s Friction coefficient',...
'Tag','edit_friction',...
'Visible','off');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[63 20 40 19],...
'String','Friction',...
'Style','text',...
'Tag','textFriction',...
'Visible','off');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[30 206 30 16],...
'String','Name',...
'Style','text',...
'UserData','main');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@swan_options_uiCB,h1,'edit_dt_CB'},...
'Position',[433 265 41 23],...
'Style','edit',...
'Tooltip','Time step for the modeling run',...
'Tag','edit_dt',...
'UserData','parameter');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@swan_options_uiCB,h1,'edit_grn_CB'},...
'Position',[620 266 31 21],...
'String','10',...
'Style','edit',...
'Tooltip','save grids at this cycle interval (it depends on dt)',...
'Tag','edit_grn',...
'UserData','parameter');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'radio_cartesian_CB'},...
'Position',[380 191 82 15],...
'String','Cartesian',...
'Style','radiobutton',...
'Tooltip','Select this if your grids are in cartesian (eg UTM) coordinates',...
'Value',1,...
'Tag','radio_cartesian',...
'UserData','parameter');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'radio_geog_CB'},...
'Position',[501 191 90 15],...
'String','Geographic',...
'Style','radiobutton',...
'Tooltip','Select this if your grids are in geographical coordinates',...
'Tag','radio_geog',...
'UserData','parameter');

uicontrol('Parent',h1, 'Position',[375 261 55 30],...
'HorizontalAlignment','right',...
'String',{'Time step'; '(seconds)' },...
'Style','text',...
'Tag','text_dt',...
'UserData','parameter');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[551 261 65 30],...
'String',{'Saving step'; '(cycle units)' },...
'Style','text',...
'Tag','text_grn',...
'UserData','parameter');

uicontrol('Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'Position',[410 304 195 16],...
'String','Minimalist parameters setting',...
'Style','text',...
'Tag','text_minimalist',...
'UserData','parameter');

uicontrol('Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'Position',[380 133 270 32],...
'String',{'OR'; 'Read parametrs from swan.par file'},...
'Style','text',...
'Tag','text_parFile',...
'UserData','parameter');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@swan_options_uiCB,h1,'edit_cumint_CB'},...
'Position',[620 230 31 21],...
'String','1',...
'Style','edit',...
'Tooltip','cumint * dt = time at which the virtual tide gauges are writen',...
'Tag','edit_cumint',...
'UserData','parameter');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[485 222 130 31],...
'String',{'Maregraph saving step'; '(time = Time step * this)'},...
'Style','text',...
'Tag','text_cumint',...
'UserData','parameter');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'push_OK_CB'},...
'FontWeight','bold',...
'Position',[214 8 66 21],...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1,...
'Call',{@swan_options_uiCB,h1,'push_Cancel_CB'},...
'FontWeight','bold',...
'Position',[294 8 66 21],...
'String','Cancel',...
'Tag','push_Cancel');

function swan_options_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

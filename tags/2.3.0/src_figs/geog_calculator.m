function varargout = geog_calculator(varargin)
% Window wraper to grdproject
%
% NOTE: I use the GMT4.0 datum definitions. In it WGS84 has an ID=220, but due to the
% fact that indices in Matlab start at 1, I sometimes have to refer at it as ID=221
% In case of datum aditions to GMT, those numbers have to be reviewed. 

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
	geog_calculator_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right');

	handles.x_left = [];        handles.y_left = [];    handles.z_left = [];
	handles.x_right = [];       handles.y_right = [];   handles.z_right = [];
	handles.fileDataLeft = [];  % Will eventually host data read from a file
	handles.fileDataRight = []; % Will eventually host data read from a file
	handles.which_conv = 1;     % Default to "Interactive Conversions"
	handles.is_ellipsoidHeight = 0; % Used to compute ellipsoidal height conv as well
	handles.dms_xinc = 0;
	handles.dms_yinc = 0;
	handles.MeasureUnit = {'meters' 'kilometers' 'nautic miles' 'miles'};
	handles.DegreeFormat1 = {'DD.xxxxxx' 'DD:MM' 'DD:MM.xxxx' 'DD:MM:SS' 'DD:MM:SS.xx'};
	handles.DegreeFormat2 = {'DD.xxxxxx' 'DD MM' 'DD MM.xxxx' 'DD MM SS' 'DD MM SS.xx'};

	handles.by_mirone2grid = 0;     handles.by_mirone2all = 0;
	if ~isempty(varargin)
		if (numel(varargin) == 2 && isstruct(varargin{1}) && ischar(varargin{2}))
			handMir = varargin{1};
			if (strcmp(varargin{2},'onlyGrid'))
				if (handMir.no_file)
					errordlg('GEOG CALCULATOR: You didn''t even load a file. What are you expecting then?','ERROR')
                    delete(hObject);    return
				end
				handles.gridLeft = getappdata(handMir.figure1,'dem_z');
				if ( isempty(handles.gridLeft) )
					errordlg('GEOG CALCULATOR: Grid is not in memory. Check your "Grid max size" and restart.','ERROR')
                    delete(hObject);    return
				end
				handles.gridLeftHead = handMir.head;
				handles.x_min = handles.gridLeftHead(1);
				handles.x_max = handles.gridLeftHead(2);
				handles.y_min = handles.gridLeftHead(3);
				handles.y_max = handles.gridLeftHead(4);
				handles.one_or_zero = ~handles.gridLeftHead(7); 
				handles.by_mirone2grid = 1;
			elseif (strcmp(varargin{2},'allTypes'))
				handles.by_mirone2all = 1;
			end
		elseif (numel(varargin) == 1 && isstruct(varargin{1}))
			handMir = varargin{1};
		end
	else
		errordlg('GEOG CALCULATOR: wrong number of arguments.','Error')
		delete(hObject);    return
	end

	% Add this Fig to the carra?as list
	plugedWin = getappdata(handMir.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handMir.figure1,'dependentFigs',plugedWin);

	handles.last_dir = handMir.last_dir;
	handles.home_dir = handMir.home_dir;
	handles.work_dir = handMir.work_dir;
	handles.version7 = handMir.version7;
	handles.path_data = handMir.path_data;

	% Import icons
	load([handles.path_data 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_gridLeft,'CData',Mfopen_ico)
	set(handles.push_gridRight,'CData',Mfopen_ico)
	set(handles.push_fileLeft,'CData',Mfopen_ico)
	set(handles.push_fileRight,'CData',Mfopen_ico)
	clear Mfopen_ico;

	coord_system_left = [];
	coord_system_right = [];
	load([handles.path_data 'mirone_pref.mat']);

	try
		if (iscell(directory_list)),    handles.work_dir = directory_list{1};
		else                            handles.work_dir = directory_list;
		end
	catch
		handles.work_dir = pwd;
	end

	% Check that the coord_system structure has no errors. If it has, load the default value.
	% The result is used to update the handles structure.
	% Firs the left side structure
	[handles,coord_system_left] = check_coord_system(handles,coord_system_left,'_left');
	handles.coord_system_left = coord_system_left;
	% And now the right side structure
	[handles,coord_system_right] = check_coord_system(handles,coord_system_right,'_right');
	handles.coord_system_right = coord_system_right;

	handles.all_ellipsoides = DefineEllipsoide;     % This is already in mirone_prefs
	handles.proj_groups = {'Longitude / Latitude';
			'Universal Transversal Mercator';...
			'World Projections';...
			'All GMT Projections';...
			'Portuguese Mess'};

	handles.all_datums = load_datums;

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	%h_t = [handles.text_CSleft handles.hTxt_in handles.hTxt_out handles.text_CSright];
	h_t = new_frame3D(hObject, [handles.text_CSleft handles.hTxt_in handles.hTxt_out handles.text_CSright]);
	handles.hTxt_in = h_t(2);		% We need to ressurect these two
	handles.hTxt_out = h_t(3);
	%------------- END Pro look (3D) -----------------------------------------------------

	%----------- Recall previous settings stored in mirone_pref -------------------
	% First the left side
	handles.txt_info_l_pos = get(handles.h_txt_info_l,'Position');
	set(handles.h_txt_info_l,'String',handles.proj_info_txt_left,'Position',handles.txt_info_l_pos)
	if (handles.is_geog_left)
		set(handles.popup_UnitesLeft,'String',handles.DegreeFormat1,'Value',handles.DegreeFormat1_val_left)
		set(handles.h_text_xLeft,'String','Longitude')
		set(handles.h_text_yLeft,'String','Latitude')
	else
		set(handles.popup_UnitesLeft,'String',handles.MeasureUnit,'Value',handles.MeasureUnit_val_left)
		set(handles.h_text_xLeft,'String','East/West')
		set(handles.h_text_yLeft,'String','North/South')
	end

	% And now the right side
	handles.txt_info_r_pos = get(handles.h_txt_info_r,'Position');
	set(handles.h_txt_info_r,'String',handles.proj_info_txt_right,'Position',handles.txt_info_r_pos)
	if (handles.is_geog_right)
		set(handles.popup_UnitesRight,'String',handles.DegreeFormat1,'Value',handles.DegreeFormat1_val_right)
		set(handles.h_text_xRight,'String','Longitude')
		set(handles.h_text_yRight,'String','Latitude')
	else
		set(handles.popup_UnitesRight,'String',handles.MeasureUnit,'Value',handles.MeasureUnit_val_right)
		set(handles.h_text_xRight,'String','East/West')
		set(handles.h_text_yRight,'String','North/South')
	end

	% This is the tag that all tab push buttons share.  If you have multiple
	% sets of tab push buttons, each group should have unique tag.
	group_name = 'tab_group';

	% This is a list of the UserData values used to link tab push buttons and the
	% components on their linked panels. To add a new tab panel to the group
	%  Add the button using GUIDE
	%  Assign the Tag based on the group name - in this case tab_group
	%  Give the UserData a unique name - e.g. another_tab_panel
	%  Add components to GUIDE for the new panel
	%  Give the new components the same UserData as the tab button
	%  Add the new UserData name to the below cell array

	if (handles.by_mirone2grid)
		panel_names = {'GridConv','interactive','FileConv'};
		h = findobj(hObject,'Style','pushbutton','String','File Conversions');
		set(h,'Visible','off')
		h = findobj(hObject,'Style','pushbutton','String','Interactive Conversions');
		set(h,'Visible','off')
		set(handles.edit_gridLeft,'Enable','off','String','In memory array')
		set(handles.push_gridLeft,'Enable','off')
		set(handles.edit_xIncLeft,'String',sprintf('%.12g',handles.gridLeftHead(8)))
		set(handles.edit_yIncLeft,'String',sprintf('%.12g',handles.gridLeftHead(9)))
		[m,n] = size(handles.gridLeft);
		set(handles.edit_nColsLeft,'String',sprintf('%d',n))
		set(handles.edit_nRowsLeft,'String',sprintf('%d',m))
		% Now the right side
		set(handles.edit_nColsRight,'String',sprintf('%d',n))
		set(handles.edit_nRowsRight,'String',sprintf('%d',m))
		set(handles.push_right2left,'Visible','off')
		bgcolor = [0.8314 0.81569 0.7843];
		set(handles.edit_gridRight,'Enable','off','BackgroundColor',bgcolor)
		set(handles.push_gridRight,'Enable','off')
		handles.which_conv = 3;
	else
		panel_names = {'interactive','FileConv','GridConv'};
		set(handles.hTxt_in,'Visible','off')
		set(handles.hTxt_out,'Visible','off')
	end

	% tabpanelfcn('makegroups',...) adds new fields to the handles structure,
	% one for each panel name and another called 'group_name_all'.  These fields
	% are used by the tabpanefcn when tab_group_handler is called.
	handles = tabpanelfcn('make_groups',group_name, panel_names, handles, 1);

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

%-------------------------------------------------------------------------------------
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over tab_group.
function tab_group_ButtonDownFcn(hObject, handles)
% Call the tab_group_handler.  This updates visiblity of components as needed to
% hide the components from the previous tab and show components on this tab.
% This also updates the last_tab field in the handles structure to keep track
% of which panel was hidden.
handles = tabpanelfcn('tab_group_handler',hObject, handles, get(hObject, 'Tag'));
str = get(hObject,'String');
switch str
    case 'Interactive Conversions'
        set(handles.push_right2left,'Visible','on')
        handles.which_conv = 1;
        if (handles.is_geog_left)
            set(handles.popup_UnitesLeft,'String',handles.DegreeFormat1,...
                'Value',handles.DegreeFormat1_val_left)
        end
        if (handles.is_geog_right)
            set(handles.popup_UnitesRight,'String',handles.DegreeFormat1,...
                'Value',handles.DegreeFormat1_val_right)
        end
	    set(handles.hTxt_in,'Visible','off')
	    set(handles.hTxt_out,'Visible','off')
    case 'File Conversions'
        set(handles.push_right2left,'Visible','off')
        handles.which_conv = 2;
        if (handles.is_geog_left)
            set(handles.popup_UnitesLeft,'String',handles.DegreeFormat2,...
                'Value',handles.DegreeFormat2_val_left)
        end
        if (handles.is_geog_right)
            set(handles.popup_UnitesRight,'String',handles.DegreeFormat2,...
                'Value',handles.DegreeFormat2_val_right)
        end
	    set(handles.hTxt_in,'Visible','on')
	    set(handles.hTxt_out,'Visible','on')
    case 'Grid Conversions'
        set(handles.push_right2left,'Visible','off')
	    set(handles.hTxt_in,'Visible','on')
	    set(handles.hTxt_out,'Visible','on')
        handles.which_conv = 3;
end
guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function x = get_editValue(hObject)
xx = get(hObject,'String');     val = test_dms(xx);     x = 0;
if isempty(xx),      x = [];     return;     end     % in case of z it may be on porpuse
if (str2double(val{1}) > 0)
    for i=1:length(val),   x = x + str2double(val{i}) / (60^(i-1));    end
else
    for i=1:length(val),   x = x - abs(str2double(val{i})) / (60^(i-1));   end
end
if (isnan(x)),   x = [];     set(hObject, 'String', '');     end

%-------------------------------------------------------------------------------------
function edit_zLeft_CB(hObject, handles)
% xx = str2double(get(hObject,'String'));
% if isnan(xx) | isempty(xx)  set(hObject, 'String', '');  return;    end    % Just a stupid user error
% handles.z_left = xx;    guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function popup_UnitesLeft_CB(hObject, handles)
fname = [handles.path_data 'mirone_pref.mat'];
load([handles.path_data 'mirone_pref.mat']);
val = get(hObject,'Value');
if (handles.is_geog_left)                       % That is, we are in geogs
    if (handles.which_conv == 1)
        handles.DegreeFormat1_val_left = val;
        coord_system_left.DegreeFormat1_val = val;
    elseif (handles.which_conv == 2)
        handles.DegreeFormat2_val_left = val;
        coord_system_left.DegreeFormat2_val = val;
    end
else
    handles.MeasureUnit_val_left = val;
    coord_system_left.MeasureUnit_val = val;
end
if (~handles.version7)                  % Update mirone_pref.  R<=13
    save(fname,'coord_system_left','-append');
else
    save(fname,'coord_system_left','-append','-v6');
end
guidata(hObject,handles)

%-------------------------------------------------------------------------------------
function popup_UnitesRight_CB(hObject, handles)
fname = [handles.path_data 'mirone_pref.mat'];
load([handles.path_data 'mirone_pref.mat']);
val = get(hObject,'Value');
if (handles.is_geog_right)                       % That is, we are in geogs
    if (handles.which_conv == 1)
        handles.DegreeFormat1_val_right = val;
        coord_system_right.DegreeFormat1_val = val;
    elseif (handles.which_conv == 2)
        handles.DegreeFormat2_val_right = val;
        coord_system_right.DegreeFormat2_val = val;
    end
else
    handles.MeasureUnit_val_right = val;
    coord_system_right.MeasureUnit_val = val;
end
if (~handles.version7)                  % Update mirone_pref.  R<=13
    save(fname,'coord_system_right','-append');
else
    save(fname,'coord_system_right','-append','-v6');
end
guidata(hObject,handles)

%-------------------------------------------------------------------------------------
function push_DefCoordLeft_CB(hObject, handles)
	fname = [handles.path_data 'mirone_pref.mat'];
	coord_system_left = coordinate_system(handles.coord_system_left, handles.all_datums);
	if (isempty(coord_system_left)),   return,		end
	handles.coord_system_left = coord_system_left;
	handles.projection_left = coord_system_left.projection;
	handles.datum_val_left = coord_system_left.datum_val;
	handles.map_scale_factor_left = coord_system_left.map_scale_factor;
	handles.system_FE_FN_left = coord_system_left.system_FE_FN;
	pos = handles.txt_info_l_pos;
	%string = {['Projection -> ' coord_system_left.ProjName];...
	string = {['System -> ' coord_system_left.SysName];...
			['Projection -> ' coord_system_left.ProjName];...
			['Datum ->  ' handles.all_datums{coord_system_left.datum_val,1}];...
			['Ellipsoide ->  ' handles.all_datums{coord_system_left.datum_val,2}];...
			['J<options> ->  ' coord_system_left.projection]};
	[outstring,newpos] = textwrap(handles.h_txt_info_l,string);
	pos(4) = newpos(4);
	set(handles.h_txt_info_l,'String',outstring,'Position',[pos(1),pos(2),pos(3),pos(4)])
	coord_system_left.proj_info_txt = outstring;
	coord_system_left.proj_info_pos = pos;
	if (isempty(handles.projection_left))           % That is, we are in geogs
		if (handles.which_conv == 1)
			set(handles.h_text_xLeft,'String','Longitude')
			set(handles.h_text_yLeft,'String','Latitude')
		end
		set(handles.popup_UnitesLeft,'String',handles.DegreeFormat1,'Value',handles.DegreeFormat1_val_left)
		handles.is_geog_left = 1;    coord_system_left.is_geog = 1;
	else
		if (handles.which_conv == 1)
			set(handles.h_text_xLeft,'String','East/West')
			set(handles.h_text_yLeft,'String','North/South')
		end
		set(handles.popup_UnitesLeft,'String',handles.MeasureUnit,'Value',handles.MeasureUnit_val_left)
		handles.is_geog_left = 0;    coord_system_left.is_geog = 0;
	end
	if (~handles.version7)                  % Update mirone_pref.  R<=13
		save(fname,'coord_system_left','-append');
	else
		save(fname,'coord_system_left','-append','-v6');
	end
	guidata(hObject,handles)

%-------------------------------------------------------------------------------------
function push_DefCoordRight_CB(hObject, handles)
	fname = [handles.path_data 'mirone_pref.mat'];
	coord_system_right = coordinate_system(handles.coord_system_right,handles.all_datums);
	if (isempty(coord_system_right)),	return,		end
	handles.coord_system_right = coord_system_right;
	handles.projection_right = coord_system_right.projection;
	handles.datum_val_right = coord_system_right.datum_val;
	handles.map_scale_factor_right = coord_system_right.map_scale_factor;
	handles.system_FE_FN_right = coord_system_right.system_FE_FN;
	pos = handles.txt_info_r_pos;
	%string = {['Projection -> ' coord_system_right.ProjName];...
	string = {['System -> ' coord_system_right.SysName];...
			['Projection -> ' coord_system_right.ProjName];...
			['Datum ->  ' handles.all_datums{coord_system_right.datum_val,1}];...
			['Ellipsoide ->  ' handles.all_datums{coord_system_right.datum_val,2}];...
			['J<options> ->  ' coord_system_right.projection]};
	[outstring,newpos] = textwrap(handles.h_txt_info_r,string);
	pos(4) = newpos(4);
	set(handles.h_txt_info_r,'String',outstring,'Position',[pos(1),pos(2),pos(3),pos(4)])
	coord_system_right.proj_info_txt = outstring;
	coord_system_right.proj_info_pos = pos;
	if (isempty(handles.projection_right))           % That is, we are in geogs
        if (handles.which_conv == 1)
			set(handles.h_text_xRight,'String','Longitude')
			set(handles.h_text_yRight,'String','Latitude')
        end
        set(handles.popup_UnitesRight,'String',handles.DegreeFormat1,'Value',handles.DegreeFormat1_val_right)
        handles.is_geog_right = 1;    coord_system_right.is_geog = 1;
	else
		if (handles.which_conv == 1)
			set(handles.h_text_xRight,'String','East/West')
			set(handles.h_text_yRight,'String','North/South')
		end
		set(handles.popup_UnitesRight,'String',handles.MeasureUnit,'Value',handles.MeasureUnit_val_right)
		handles.is_geog_right = 0;    coord_system_right.is_geog = 0;
	end
	if (~handles.version7)                  % Update mirone_pref.  R<=13
		save(fname,'coord_system_right','-append');
	else
		save(fname,'coord_system_right','-append','-v6');
	end
	guidata(hObject,handles)

% -------------------------------------------------------------------------------------
function push_left2right_CB(hObject, handles)
	msg_err = [];   zl = [];
	if (handles.which_conv == 1)        % Interactive Conversions
		xl = get_editValue(handles.edit_xLeft);
		yl = get_editValue(handles.edit_yLeft);
		zl = get_editValue(handles.edit_zLeft);
		x_c = xl;   y_c = yl;
		if (~isempty(zl)),  in = [xl yl zl];
		else                in = [xl yl];   end
		if (isempty(xl))
			msg_err = 'You have to agree that I cannot convert a NOTHING (see the left input X box).';
		end
		if (isempty(yl))
			msg_err = 'You have to agree that a NOTHING cannot be converted (see the left input Y box).';
		end
	elseif (handles.which_conv == 2)    % File Conversions
		if (isempty(handles.fileDataLeft))
			msg_err = 'For converting a data file I need to know ... the data!';
		else
			if (get(handles.checkbox_ToggleXYLeft,'Value'))
				in = [handles.fileDataLeft(:,2) handles.fileDataLeft(:,1) handles.fileDataLeft(:,3:end)];
			else
				in = handles.fileDataLeft;
			end
			x_max = max(in(:,1));        x_min = min(in(:,1));
			y_max = max(in(:,2));        y_min = min(in(:,2));
			% Find means used to estimate the "good" opt_R
			x_c = (x_min + x_max) / 2;   y_c = (y_min + y_max) / 2;
			[n_row,n_col] = size(handles.fileDataLeft);
		end
	elseif (handles.which_conv == 3)
		msg_err = transform_grid(handles);		% Do the work there and return.
		if (~isempty(msg_err)),		errordlg(msg_err,'Error'),		end
		return
	end

	if (strcmp(get(handles.h_txt_info_l,'String'),'Nikles') | strcmp(get(handles.h_txt_info_r,'String'),'Nikles'))
		msg_err = 'You need to select BOTH Input and Output coordinate systems.';
	end

	if (~isempty(msg_err)),   errordlg(msg_err,'Error'),	return,		end

	small = 1;
	opt_R = sprintf('-R%.10g/%.10g/%.10g/%.10g', x_c-small, x_c+small, y_c-small, y_c+small);

	% Check if we need a datum conversion as well
	if (handles.datum_val_left ~= handles.datum_val_right)
		if (handles.is_ellipsoidHeight),     opt_T = '-Th';
		elseif ((handles.which_conv == 1) && ~isempty(zl)),  opt_T = '-Th';
		else    opt_T = '-T';
		end
		opt_T = sprintf('%s%d/%d',opt_T,handles.datum_val_left-1,handles.datum_val_right-1);
	else
		opt_T = ' ';
	end

	% Need this only in the projection-to-projection case
	if (handles.datum_val_left == 221 && handles.datum_val_right == 221)
		is_wgs84 = 1;
	else	
		is_wgs84 = 0;
	end

	% Check if we have a MAP_SCALE_FACTOR
	opt_SF = ' ';
	if (~isempty(handles.map_scale_factor_right) && handles.map_scale_factor_right ~= 1)
		opt_SF = sprintf('--MAP_SCALE_FACTOR=%.4f',handles.map_scale_factor_right);
	end

	% Check if we have false eastings/northings
	opt_C = '-C';
	if (~isempty(handles.system_FE_FN_right))
		opt_C = sprintf( '-C%.3f/%.3f',handles.system_FE_FN_right(1), handles.system_FE_FN_right(2) );
	elseif (~isempty(handles.system_FE_FN_left))		% Look at here as well
		opt_C = sprintf( '-C%.3f/%.3f',handles.system_FE_FN_left(1),  handles.system_FE_FN_left(2) );
	end

	if (~isempty(handles.projection_right))
		% Check if output is in other unites than meters
		str = get(handles.popup_UnitesRight, 'String');
		switch lower(str{handles.MeasureUnit_val_right});
			case 'meters',      opt_F = '-F';
			case 'kilometers',  opt_F = '-Fk';
			case 'nautic miles',opt_F = '-Fn';
			case 'miles',       opt_F = '-Fm';
			otherwise,          opt_F = '-F';
		end
	else
		opt_F = '-F';
	end

	opt_Rg = '-R-180/180/0/80';
	if (isempty(handles.projection_left) && isempty(handles.projection_right & ~strcmp(opt_T,' ')))
		% Do only a datum conversion, but first... test
		if (strcmp(opt_T,' '))
			errordlg('Start and ending datums are equal. Nothing to do.','Error')
			return
		end
		out = mapproject_m(in, opt_T);
		set(hObject,'TooltipString',sprintf('%s\n%s','Last command:', opt_T))
	elseif (~isempty(handles.projection_left) && ~isempty(handles.projection_right))
		% More complicated. To go from one projection to other we have to pass by geogs.
		% First do the inverse conversion (get result in geogs)
		opt_J = handles.projection_left;
		% First apply the trick to get a good estimation of -R
		if (~is_wgs84)      % Otherwise, no need to waste time with datum conversions
			opt_T = sprintf('-T%d/220',handles.datum_val_left-1);
		else    opt_T = ' ';
		end
		if (y_c < 0),   opt_Rg = '-R-180/180/-80/0';    end         % Patches over inventions, not good
		opt_SF = ' ';
		if (~isempty(handles.map_scale_factor_left) && handles.map_scale_factor_left ~= 1)
			opt_SF = sprintf('--MAP_SCALE_FACTOR=%.4f',handles.map_scale_factor_left);
		end
		tmp = mapproject_m([x_c y_c], opt_T, opt_J, opt_Rg, opt_SF, opt_C, opt_F, '-I');
		opt_R = sprintf('-R%.3f/%.3f/%.3f/%.3f', tmp(1)-small, tmp(1)+small, tmp(2)-small, tmp(2)+small);
		% Here we have to re-check the need for ellipsoidal height transformations
		if (~strcmp(opt_T,' '))
			if (handles.is_ellipsoidHeight),    opt_T = ['-Th' opt_T(3:end)];
			else                                opt_T = ['-T' opt_T(3:end)];
			end
		end
		out = mapproject_m(in, opt_T, opt_J, opt_R, opt_SF, opt_C, opt_F, '-I');
		comm1 = [opt_R ' ' opt_J ' ' opt_C ' ' opt_F ' ' opt_T ' ' opt_SF ' -I'];
		% And now do a direct conversion to the final destination
		opt_J = handles.projection_right;		opt_T = ' ';
		if (~is_wgs84)
			if (handles.is_ellipsoidHeight)
				opt_T = sprintf('-Th220/%d',handles.datum_val_right-1);
			else
				opt_T = sprintf('-T220/%d',handles.datum_val_right-1);
			end
		end
		try     % If right value exists, we need it. Otherwise, don't use shifts
			opt_C = sprintf( '-C%.3f/%.3f',handles.system_FE_FN_right(1),  handles.system_FE_FN_right(2) );
		catch
			opt_C = '-C';
		end
		opt_SF = ' ';
		if (~isempty(handles.map_scale_factor_right) && handles.map_scale_factor_right ~= 1)
			opt_SF = sprintf('--MAP_SCALE_FACTOR=%.4f',handles.map_scale_factor_right);
		end
		out = mapproject_m(out, opt_T, opt_J, opt_R, opt_SF, opt_C, opt_F);
		comm = sprintf('%s\n%s\n%s','Last command:', comm1, [opt_R ' ' opt_J ' ' opt_C ' ' opt_F ' ' opt_T ' ' opt_SF]);
		set(hObject,'TooltipString',comm)
	elseif (isempty(handles.projection_right) && ~isempty(handles.projection_left))
		%       GEOG at right                            NON-GEOG at left.  Do a inverse transformation
		opt_J = handles.projection_left;
		% First apply the trick to get a good estimation of -R
		if (y_c < 0),   opt_Rg = '-R-180/180/-80/0';    end         % Patches over inventions, not good
		tmp = mapproject_m([x_c y_c], opt_T, opt_J, opt_Rg, opt_SF, opt_C, opt_F, '-I');
		opt_R = ['-R' sprintf('%.10g/%.10g/%.10g/%.10g', tmp(1)-small, tmp(1)+small, tmp(2)-small, tmp(2)+small)];
		out = mapproject_m(in, opt_T, opt_J, opt_R, opt_SF, opt_C, opt_F, '-I');
		comm = sprintf('%s\n%s','Last command:', [opt_R ' ' opt_J ' ' opt_C ' ' opt_F ' -I ' opt_T ' ' opt_SF]);
		set(hObject,'TooltipString',comm)
	elseif (~isempty(handles.projection_right) && isempty(handles.projection_left))
		%       NON-GEOG at right                           GEOG at left.  Do a direct transformation
		opt_J = handles.projection_right;
		out = mapproject_m(in, opt_T, opt_J, opt_R, opt_SF, opt_C, opt_F);
		comm = sprintf('%s\n%s','Last command:', [opt_R ' ' opt_J ' ' opt_C ' ' opt_F ' ' opt_T ' ' opt_SF]);
		set(hObject,'TooltipString',comm)
	else
		return
	end                 % Otherwise -> BOOM

% ------ Update Box or write result file ---------------------------------------------
	if (handles.which_conv == 1)        % Interactive Conversions
        output_format(handles,out(1),out(2),'right');
        if (~isnan(zl));    set(handles.edit_zRight,'String',sprintf('%.6g',out(3)));  end
	elseif (handles.which_conv == 2)    % File Conversions
        fname = get(handles.edit_fileRight,'String');
        if (isempty(fname))
			[FileName,PathName] = put_or_get_file(handles, ...
				{'*.dat;*.DAT', 'Coord file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select data file','put','.dat');
			if isequal(FileName,0);     return;     end
			fname = [PathName FileName];
        end
        % Now write the results in the file fname
        format = repmat('%f\t',1,n_col);
        format = format(1:end-2);               % We don't want the last '\t'
        try
			double2ascii(fname,out,format);
			msgbox('File successefully writen on disk','Sorte');
        catch
			errordlg('There was an error (origin unknown) while writing file on disk','Error');
        end
	end

%-------------------------------------------------------------------------------------
function push_right2left_CB(hObject, handles)
	% Given the current mapproject limitation in invertion transformations I use a trick that
	% seams to work. The trick consists in making a first transformation with a global -R
	% and than use the output from it to find a new -R that allows a correct (?) transformation
	msg_err = [];
	if (handles.which_conv == 1)        % Interactive Conversions
		xr = get_editValue(handles.edit_xRight);
		yr = get_editValue(handles.edit_yRight);
		zr = get_editValue(handles.edit_zRight);
		x_c = xr;		y_c = yr;
		if (~isempty(zr)),	in = [xr yr zr];
		else				in = [xr yr];
		end
        if (isempty(xr))
			msg_err = 'You have to agree that I cannot convert a NOTHING (see the right input X box).';
        end
        if (isempty(yr))
			msg_err = 'You have to agree that a NOTHING cannot be converted (see the right input Y box).';
        end
	elseif (handles.which_conv == 2)    % File Conversions
		if (isempty(handles.fileDataRight))
			msg_err = 'For converting a data file I need to know ... the data!';
		else
			if (get(handles.checkbox_ToggleXYRight,'Value'))
				in = [handles.fileDataRight(:,2) handles.fileDataRight(:,1) handles.fileDataRight(:,3:end)];
			else
				in = handles.fileDataRight;
			end
			x_max = max(in(:,1));        x_min = min(in(:,1));
			y_max = max(in(:,2));        y_min = min(in(:,2));
			% Find means used to estimate the "good" opt_R
			x_c = (x_min + x_max) / 2;   y_c = (y_min + y_max) / 2;
			[n_row,n_col] = size(handles.fileDataRight);
		end
	elseif (handles.which_conv == 3)
		msg_err = transform_grid(handles);        % Do the work there and return.
		if (~isempty(msg_err)),     errordlg(msg_err,'Error'),	end
		return
	end
	
	if (~isempty(msg_err)),     errordlg(msg_err,'Error'),		end

	small = 1;
	opt_R = sprintf('-R%.10g/%.10g/%.10g/%.10g', x_c-small, x_c+small, y_c-small, y_c+small);

	% Check if we need a datum conversion as well
	if (handles.datum_val_left ~= handles.datum_val_right)
		if (handles.is_ellipsoidHeight),    opt_T = '-Th';
		elseif ((handles.which_conv == 1) && ~isempty(zr)),  opt_T = '-Th';
		else    opt_T = '-T';
		end
		opt_T = sprintf('%s%d/%d',opt_T,handles.datum_val_right-1,handles.datum_val_left-1);
	else
		opt_T = ' ';
	end

	% This is need only in the projection-to-projection case
	if (handles.datum_val_left == 221 && handles.datum_val_right == 221)
		is_wgs84 = 1;
	else    
		is_wgs84 = 0;
	end

	% Check if we have a MAP_SCALE_FACTOR
	opt_SF = ' ';
	if (~isempty(handles.map_scale_factor_left) && handles.map_scale_factor_left ~= 1)
		opt_SF = sprintf('--MAP_SCALE_FACTOR=%.4f',handles.map_scale_factor_left);
	end

	% Check if we have false eastings/northings
	opt_C = '-C';
	if (~isempty(handles.system_FE_FN_right))
		opt_C = sprintf( '-C%.4f/%.3f',handles.system_FE_FN_right(1), handles.system_FE_FN_right(2) );
	elseif (~isempty(handles.system_FE_FN_left))		% Look at here as well
		opt_C = sprintf( '-C%.4f/%.3f',handles.system_FE_FN_left(1),  handles.system_FE_FN_left(2) );
	end

	if (~isempty(handles.projection_left))
		% Check if output is in other unites than meters
		str = get(handles.popup_UnitesLeft, 'String');
		switch lower(str{handles.MeasureUnit_val_right});
            case 'meters',      opt_F = '-F';
            case 'kilometers',  opt_F = '-Fk';
            case 'nautic miles',opt_F = '-Fn';
            case 'miles',       opt_F = '-Fm';
            otherwise,          opt_F = '-F';
		end
	else
        opt_F = '-F';
	end

opt_Rg = '-R-180/180/0/80';
if (isempty(handles.projection_left) && isempty(handles.projection_right & ~strcmp(opt_T,' ')))
    % Do only a datum conversion, but first... test
    if (strcmp(opt_T,' '))
        errordlg('Start and ending datums are equal. Nothing to do.','Error')
        return
    end
    out = mapproject_m(in, opt_T);
    set(hObject,'TooltipString',sprintf('%s\n%s','Last command:', opt_T))
elseif (~isempty(handles.projection_left) && ~isempty(handles.projection_right))
	% More complicated. To go from one projection to other we have to pass by geogs.
	% First do the inverse conversion (that is, get geogs).
	opt_J = handles.projection_right;
	% First apply the trick to get a good estimation of -R
	if (~is_wgs84)		% Otherwise, no need to waste time with datum conversions
		opt_T = sprintf('-T%d/220',handles.datum_val_right-1);
	else    opt_T = ' ';
	end
	if (y_c < 0),   opt_Rg = '-R-180/180/-80/0';    end         % Patches over inventions, not good
	opt_SF = ' ';
	if (~isempty(handles.map_scale_factor_right) && handles.map_scale_factor_right ~= 1)
		opt_SF = sprintf('--MAP_SCALE_FACTOR=%.4f',handles.map_scale_factor_right);
	end
	tmp = mapproject_m([x_c y_c], opt_T, opt_J, opt_Rg, opt_SF, opt_C, opt_F, '-I');
	opt_R = sprintf('-R%.3f/%.3f/%.3f/%.3f', tmp(1)-small, tmp(1)+small, tmp(2)-small, tmp(2)+small);
	% Here we have to re-check the need for ellipsoidal height transformations
	if (~strcmp(opt_T,' '))
		if (handles.is_ellipsoidHeight),    opt_T = ['-Th' opt_T(3:end)];
		else                                opt_T = ['-T' opt_T(3:end)];
		end
	end
	out = mapproject_m(in, opt_T, opt_J, opt_R, opt_SF, opt_C, opt_F, '-I');
	comm1 = [opt_R ' ' opt_J ' ' opt_C ' ' opt_F ' ' opt_T ' ' opt_SF ' -I'];
	% And now do a direct conversion to the final destination
	opt_J = handles.projection_left;		opt_T = ' ';
	if ( ~is_wgs84 && (handles.datum_val_left ~= 221) )
		if (handles.is_ellipsoidHeight)
			opt_T = sprintf('-Th220/%d',handles.datum_val_left-1);
		else
			opt_T = sprintf('-T220/%d',handles.datum_val_left-1);
		end
	end
	try     % If left value exists, we need it. Otherwise, don't use shifts
		opt_C = sprintf( '-C%.3f/%.3f',handles.system_FE_FN_left(1),  handles.system_FE_FN_left(2) );
	catch
		opt_C = '-C';
	end
	opt_SF = ' ';
	if (~isempty(handles.map_scale_factor_left) && handles.map_scale_factor_left ~= 1)
		opt_SF = sprintf('--MAP_SCALE_FACTOR=%.4f',handles.map_scale_factor_left);
	end
	out = mapproject_m(out, opt_T, opt_J, opt_R, opt_SF, opt_C, opt_F);
	comm = sprintf('%s\n%s\n%s','Last command:', comm1, [opt_R ' ' opt_J ' ' opt_C ' ' opt_F ' ' opt_T ' ' opt_SF]);
	set(hObject,'TooltipString',comm)
elseif (~isempty(handles.projection_right) && isempty(handles.projection_left))
    %       GEOG at left                            NON-GEOG at right.  Do direc transformation
    opt_J = handles.projection_right;
    % First apply the trick to get a good estimation of -R
    if (y_c < 0),   opt_Rg = '-R-180/180/-80/0';    end         % Patches over inventions, not good
    tmp = mapproject_m([x_c y_c], opt_T, opt_J, opt_Rg, opt_SF, opt_C, opt_F, '-I');
    opt_R = ['-R' sprintf('%.10g/%.10g/%.10g/%.10g', tmp(1)-small, tmp(1)+small, tmp(2)-small, tmp(2)+small)];
    out = mapproject_m(in, opt_T, opt_J, opt_R, opt_SF, opt_C, opt_F, '-I');
    comm = sprintf('%s\n%s','Last command:', [opt_R ' ' opt_J ' ' opt_C ' ' opt_F ' -I ' opt_T ' ' opt_SF]);
    set(hObject,'TooltipString',comm)
elseif (isempty(handles.projection_right) && ~isempty(handles.projection_left))
    %       NON-GEOG at left                           GEOG at right.  Do inverse transformation
    opt_J = handles.projection_left;
    out = mapproject_m(in, opt_T, opt_J, opt_R, opt_SF, opt_C, opt_F);
    comm = sprintf('%s\n%s','Last command:', [opt_R ' ' opt_J ' ' opt_C ' ' opt_F ' ' opt_T ' ' opt_SF]);
    set(hObject,'TooltipString',comm)
else
    return
end                 % Otherwise -> BOOM

if (handles.which_conv == 1)        % Interactive Conversions
    output_format(handles,out(1),out(2),'left');
    if (~isnan(zr));    set(handles.edit_zLeft,'String',sprintf('%.6g',out(3)));  end
elseif (handles.which_conv == 2)    % File Conversions
    fname = get(handles.edit_fileLeft,'String');
    if (isempty(fname))
        [FileName,PathName] = put_or_get_file(handles, ...
			{'*.dat;*.DAT', 'Coord file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select data file', 'put','.dat');
        if isequal(FileName,0);     return;     end
        fname = [PathName FileName];
    end
    % Now write the results in the file fname
    format = repmat('%f\t',1,n_col);
    format = format(1:end-2);               % We don't want the last '\t'
    try
        double2ascii(fname,out,format);
        msgbox('File successefully writen on disk','Sorte');
    catch
        errordlg('There was an error (origin unknown) while writing file on disk','Error');
    end
end

%-------------------------------------------------------------------------------------
function msg_err = transform_grid(handles)
	msg_err = [];   opt_F = ' ';
	if (isempty(handles.gridLeft))
		msg_err = 'For converting a grid I need to know the ... guess what?';
	else
		if (get(handles.checkbox_toggleRegist,'Value'))
			opt_F = '-F';
		end
		x_max = handles.gridLeftHead(2);        x_min = handles.gridLeftHead(1);
		y_max = handles.gridLeftHead(4);        y_min = handles.gridLeftHead(3);
		% Find means used to estimate the "good" opt_R for ...
		x_c = (x_min + x_max) / 2;   y_c = (y_min + y_max) / 2;
	end
	
	if ( ~handles.by_mirone2grid && isempty(get(handles.edit_gridRight,'String')) )
		msg_err = 'Need output grid name.';
	end

	if (strcmp(get(handles.h_txt_info_l,'String'),'Nikles') | strcmp(get(handles.h_txt_info_r,'String'),'Nikles'))
		msg_err = 'You need to select BOTH Input and Output coordinate systems.';
	end

	str_c = get(handles.edit_nColsRight,'String');
	str_r = get(handles.edit_nRowsRight,'String');
	if (isempty(str_c) || isempty(str_r))
		msg_err = 'You managed to make desapear the number of columns or rows. Congratulations';
	else
		opt_N = ['-N' str_c '/' str_r];
	end

	% If we got an error, return here
	if (~isempty(msg_err)),   return;   end

	opt_R = sprintf('-R%.10g/%.10g/%.10g/%.10g',x_min,x_max,y_min,y_max);

	% Check if we have a MAP_SCALE_FACTOR
	if (~isempty(handles.map_scale_factor_right))
		opt_SF = sprintf('--MAP_SCALE_FACTOR=%.6f',handles.map_scale_factor_right);
	else    opt_SF = ' ';
	end

	% Check if we have false eastings/northings
	if (~isempty(handles.system_FE_FN_right))
		opt_C = sprintf('-C%.9g/%.9g',handles.system_FE_FN_right(1), handles.system_FE_FN_right(2));
	else    opt_C = '-C';
	end

	if (~isempty(handles.projection_right))
		% Check if output is in other unites than meters
		str = get(handles.popup_UnitesRight, 'String');
		switch lower(str{handles.MeasureUnit_val_right});
			case 'meters',      opt_A = '-A';
			case 'kilometers',  opt_A = '-Ak';
			case 'nautic miles',opt_A = '-An';
			case 'miles',       opt_A = '-Am';
			otherwise,          opt_A = '-A';
		end
	else
		opt_A = '-A';
	end
	
	if (length(opt_A) > 2),     opt_mpF = ['-F' opt_A(3)];  % mapproject option -F
	else                        opt_mpF = '-F';
	end

	set(handles.figure1,'pointer','watch')
	opt_Rg = '-R-180/180/0/80';       opt_I = '';
	if (~isempty(handles.projection_left) && ~isempty(handles.projection_right))
		% More complicated. To go from one projection to other we have to pass by geogs.
		% First do the inverse conversion (get result in geogs)
		opt_J = handles.projection_left;
		% First apply the trick to get a good estimation of -R
		if (y_c < 0),   				opt_Rg = '-R-180/180/-80/0';		% Patches over inventions, not good
		elseif (x_c == 0 && y_c == 0)	opt_Rg = '-R-180/180/-80/80';		% Antartica cases
		end
		tmp = mapproject_m([x_c y_c], opt_J, opt_Rg, opt_SF, opt_C, opt_mpF, '-I');
		opt_R = sprintf('-R%.10g/%.10g/%.10g/%.10g', tmp(1)-1, tmp(1)+1, tmp(2)-1, tmp(2)+1);

		% Now find the exact lims in geogs
		tmp = mapproject_m([x_min y_min; x_max y_max], opt_J, opt_R, opt_SF, opt_C, opt_mpF, '-I');
		opt_R = sprintf('-R%.10g/%.10g/%.10g/%.10g',tmp(1,1),tmp(2,1),tmp(1,2),tmp(2,2));

		% Convert to Geogs
		Z = grdproject_m(handles.gridLeft, handles.gridLeftHead, opt_J, opt_R,...
			opt_SF, opt_C, opt_F, opt_A, opt_N, '-I');

		% And now do a direct conversion to the final destination
		opt_J = handles.projection_right;
		try     % If right value exists, we need it. Otherwise, don't use shifts
			opt_C = sprintf('-C%.9g/%.9g',handles.system_FE_FN_right(1), handles.system_FE_FN_right(2));
		catch
			opt_C = '-C';
		end
		[Z,head] = grdproject_m(Z, handles.gridLeftHead, opt_J, opt_R, opt_SF, opt_C, opt_F, opt_A, opt_N);
	elseif (isempty(handles.projection_right) && ~isempty(handles.projection_left))
		%       GEOG at right                            NON-GEOG at left.  Do a inverse transformation
		opt_J = handles.projection_left;
		% First apply the trick to get a good estimation of -R
		if (y_c < 0),   opt_Rg = '-R-180/180/-80/0';    end         % Patches over inventions, not good
		tmp = mapproject_m([x_c y_c], opt_J, opt_Rg, opt_SF, opt_C, opt_mpF, '-I');
		opt_R = sprintf('-R%.10g/%.10g/%.10g/%.10g', tmp(1)-1, tmp(1)+1, tmp(2)-1, tmp(2)+1);
	
		% Now find the exact lims in geogs
		tmp = mapproject_m([x_min y_min; x_max y_max], opt_J, opt_R, opt_SF, opt_C, opt_mpF, '-I');
		opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f',tmp(1,1),tmp(2,1),tmp(1,2),tmp(2,2));

		[Z,head] = grdproject_m(handles.gridLeft, handles.gridLeftHead, opt_A, opt_J, opt_R,...
			opt_SF, opt_C, opt_F, opt_N, '-I');
		opt_I = '-I';		% Used to inform Mirone that we came from an inverse transform
	elseif (~isempty(handles.projection_right) && isempty(handles.projection_left))
		%       NON-GEOG at right                           GEOG at left.  Do a direct transformation
		opt_J = handles.projection_right;
		[Z,head] = grdproject_m(handles.gridLeft, handles.gridLeftHead, opt_J, opt_R, opt_SF,...
			opt_C, opt_F, opt_A, opt_N);
	else    
        return
	end                % Otherwise -> BOOM

	comm = ['grdproject ' opt_J ' ' opt_R ' ' opt_SF ' ' opt_C ' ' opt_A ' ' opt_I];
	set(handles.push_left2right,'TooltipString',comm)
	if (numel(Z) < 4)
		warndlg('The projected grid has less than 4 elements. This prj thing f again didn''t it?','Warnerr')
		return
	end
	set(handles.edit_xIncRight,'String',head(8))
	set(handles.edit_yIncRight,'String',head(9))
	setappdata(handles.figure1,'headProj',head)			% Store projected head for use in eventual new nRow|Cols

	set(handles.figure1,'pointer','arrow')
	if (~handles.by_mirone2grid)
		tit = 'Grid converted by grdproject';
		fname = get(handles.edit_gridRight,'String');
		% Defaults and srsWKT fishing are set in nc_io
		misc = struct('x_units',[],'y_units',[],'z_units',[],'z_name',[],'desc',[], ...
			'title',tit,'history',[],'srsWKT',[], 'strPROJ4',[]);
		nc_io(fname, 'w', handles, Z, misc)
	else    % Output the converted grid to Mirone
		[ny,nx] = size(Z);
		zMinMax = grdutils(Z,'-L');
		tmp.head = [head(1:4) zMinMax(1:2)' head(7:9)];
		tmp.X = linspace(head(1),head(2),nx);		tmp.Y = linspace(head(3),head(4),ny);
		tmp.name = 'Projected grid';				tmp.ProjGMT = opt_J;
		mirone(Z,tmp);
	end

%-------------------------------------------------------------------------------------
function edit_fileLeft_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname)               % Reset the file computing section
		set(handles.edit_fileRight,'String','')
		set(handles.checkbox_nHeadersLeft,'Value',0)
		set(handles.edit_nHeadersLeft,'String','','Enable','inactive')
		handles.input_file = [];
		return
	end
	% Let the push_fileLeft_CB do all the work
	push_fileLeft_CB(handles.push_fileLeft, handles, fname)

%-------------------------------------------------------------------------------------
function push_fileLeft_CB(hObject, handles, opt)
	if (nargin == 3),   fname = opt;
	else    opt = [];
	end

	if (isempty(opt))    % Otherwise we already know fname from the 4th input argument
		[FileName,PathName] = put_or_get_file(handles,{'*.dat;*.DAT', 'Mag file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select data file','get');
		if isequal(FileName,0);     return;     end
		fname = [PathName FileName];
		set(handles.edit_fileLeft,'String',fname)
	end

	hFig = gcf;

	% See if number of headers are requested
	n_headers = [];
	if (get(handles.checkbox_nHeadersLeft,'Value'))
		xx = get(handles.edit_nHeadersLeft,'String');
		if (~isempty(xx))
			n_headers = str2double(xx);
		end
	end

	set(hFig,'Name',['Reading file:' fname])
	if (isempty(n_headers))
		[handles.fileDataLeft,date,headerlines,str_col] = text_read(fname,NaN);
	else
		[handles.fileDataLeft,date,headerlines,str_col] = text_read(fname,NaN,n_headers);
	end
	set(hFig,'Name','Geographic Computator')
	if (headerlines > 0)        % If file had headers, update respective fields
		set(handles.checkbox_nHeadersLeft,'Value',1)
		set(handles.edit_nHeadersLeft,'String',headerlines)
	end

	% If msgbox exist we have to move it from behind the main window. So get it's handle
	hMsgFig = gcf;
	if (hFig ~= hMsgFig)
		figure(hMsgFig);    % If error msgbox exists, bring it forward
	end

	if (isempty(handles.fileDataLeft))
		errordlg('There is nothing numerically usefull inside your file. Quiting','Error')
		return
	end
	
	[m,n] = size(handles.fileDataLeft);
	if (isempty(handles.fileDataLeft) || n < 2)   % Even if the rest is not empty we have to quit and return
		errordlg('File doesn''t have any recognized numeric data OR only one column (Quiting).','Error');
		set(handles.edit_fileLeft,'String','')
		set(handles.checkbox_nHeadersLeft,'Value',0)
		set(handles.edit_nHeadersLeft,'String','','Enable','inactive')
		set(handles.checkbox_ToggleXYLeft,'Value',0)
		set(handles.checkbox_ellipsoidHeightsLeft,'Enable','off');
		return
	end

%handles.DegreeFormat2 = {'DD.xxxxxx' 'DD MM' 'DD MM.xxxx' 'DD MM SS' 'DD MM SS.xx'};
if (handles.is_geog_left)
	switch handles.DegreeFormat2_val_left
		case 1,     min_expected_col = 2;
		case 2,     min_expected_col = 4;
		case 3,     min_expected_col = 4;
		case 4,     min_expected_col = 6;
		case 5,     min_expected_col = 6;
	end

	if (n < min_expected_col)
		errordlg('File does not have enough columns to satisfy the inicated degree format','Error')
		set(handles.edit_fileLeft,'String','')
		set(handles.checkbox_nHeadersLeft,'Value',0)
		set(handles.edit_nHeadersLeft,'String','','Enable','inactive')
		set(handles.checkbox_ToggleXYLeft,'Value',0)
		set(handles.checkbox_ellipsoidHeightsLeft,'Enable','off');
		handles.fileDataLeft = [];        guidata(hObject, handles);
		return
    else
		if (min_expected_col == 4)
			x = dms2degree(handles.fileDataLeft(:,1),handles.fileDataLeft(:,2));
			y = dms2degree(handles.fileDataLeft(:,3),handles.fileDataLeft(:,4));
			handles.fileDataLeft(:,1) = x;      handles.fileDataLeft(:,2) = y;
			handles.fileDataLeft(:,3:4) = [];
		elseif (min_expected_col == 6)
			x = dms2degree(handles.fileDataLeft(:,1),handles.fileDataLeft(:,2),handles.fileDataLeft(:,3));
			y = dms2degree(handles.fileDataLeft(:,4),handles.fileDataLeft(:,5),handles.fileDataLeft(:,6));
			handles.fileDataLeft(:,1) = x;      handles.fileDataLeft(:,2) = y;
			handles.fileDataLeft(:,3:6) = [];
		end
	end
else
	min_expected_col = 2;
end

if (n > min_expected_col),     set(handles.checkbox_ellipsoidHeightsLeft,'Enable','on');   end

[pathstr,name] = fileparts(fname);
fname = [pathstr filesep name '_conv.dat'];     % Proposed output filename
set(handles.edit_fileRight,'String',fname)
guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function checkbox_nHeadersLeft_CB(hObject, handles)
if (get(hObject,'Value')),   set(handles.edit_nHeadersLeft,'Enable','on')
else                        set(handles.edit_nHeadersLeft,'Enable','off','String','')
end

%-------------------------------------------------------------------------------------
function edit_nHeadersLeft_CB(hObject, handles)
n = str2double(get(hObject,'String'));
if (isnan(n))   set(hObject,'String','');    end

%-------------------------------------------------------------------------------------
function checkbox_ToggleXYLeft_CB(hObject, handles)
% This case is handled by the push_left2right

%-------------------------------------------------------------------------------------
function checkbox_ellipsoidHeightsLeft_CB(hObject, handles)
handles.is_ellipsoidHeight = get(hObject,'Value');
guidata(hObject,handles)

%-------------------------------------------------------------------------------------
function push_fileRight_CB(hObject, handles)
	[FileName,PathName] = put_or_get_file(handles, ...
		{'*.dat;*.DAT', 'Mag file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select data file', 'put','.dat');
	if isequal(FileName,0),		return,		end
	set(handles.edit_fileRight,'String',[PathName FileName])
	guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function edit_fileRight_CB(hObject, handles)
% File name will be read by the push_left2right

%-------------------------------------------------------------------------------------
function out = output_format(handles,x,y,side)
% Format the output and update the edit boxes with the result.
if (strcmp(side,'left'))
    is_geog = handles.is_geog_left;
    val = handles.DegreeFormat1_val_left;
else
    is_geog = handles.is_geog_right;
    val = handles.DegreeFormat1_val_right;
end

if (is_geog)                       % That is, we are in geogs
    switch handles.DegreeFormat1{val}
        case 'DD:MM'
            ang = degree2dms(x,'DDMM',0,'str');
            str_x = [ang.dd ':' ang.mm];
            ang = degree2dms(y,'DDMM',0,'str');
            str_y = [ang.dd ':' ang.mm];
        case 'DD:MM.xxxx'
            ang = degree2dms(x,'DDMM.x',4,'str');
            str_x = [ang.dd ':' ang.mm];
            ang = degree2dms(y,'DDMM.x',4,'str');
            str_y = [ang.dd ':' ang.mm];
        case 'DD:MM:SS'
            ang = degree2dms(x,'DDMMSS',0,'str');
            str_x = [ang.dd ':' ang.mm ':' ang.ss];
            ang = degree2dms(y,'DDMMSS',0,'str');
            str_y = [ang.dd ':' ang.mm ':' ang.ss];
        case 'DD:MM:SS.xx'
            ang = degree2dms(x,'DDMMSS.x',2,'str');
            str_x = [ang.dd ':' ang.mm ':' ang.ss];
            ang = degree2dms(y,'DDMMSS.x',2,'str');
            str_y = [ang.dd ':' ang.mm ':' ang.ss];
        otherwise       % Just the DD.xxxxxxx
            str_x = sprintf('%.7f',x);
            str_y = sprintf('%.7f',y);
    end
else        % Cartesian unities
    str_x = sprintf('%.6f',x);
    str_y = sprintf('%.6f',y);
end

if (nargout == 0)
	if (strcmp(side,'left'))
        set(handles.edit_xLeft,'String',str_x)
        set(handles.edit_yLeft,'String',str_y)
	else
        set(handles.edit_xRight,'String',str_x)
        set(handles.edit_yRight,'String',str_y)
	end
	out = [];
else
    out = {str_x str_y};
end

%-------------------------------------------------------------------------------------
function tab_group_CB(hObject, handles)

%-------------------------------------------------------------------------------------
function [handles,out] = check_coord_system(handles,coord_system,side)
% Currently the coord_system structure must contain the following fields:
%              group_val
%             system_val
%              datum_val
%        cilindrical_val
%          azimuthal_val
%              conic_val
%       miscelaneous_val
%               ProjName    % Projection Name (used in the text info)
%       map_scale_factor
%        MeasureUnit_val    % Currently selected measure unitie
%      DegreeFormat1_val    % Degree format value for point conversions
%      DegreeFormat2_val    % Degree format value for file conversions
%           system_FE_FN    % False eastings/northings
%                is_geog    % Signals if we have geogs coordinates
%             projection    % -J string
%     ProjParameterValue    % projection parameters
%          proj_info_txt
%          proj_info_pos    % This one has to be dealt inside the opening function (don't know here a def value)
%
% If any of those is missing, assign it a default value

if (isempty(coord_system))   % If it doesn't exist, create an empty one
    coord_system = struct([]);
end

if (~isfield(coord_system,'group_val')),     out.group_val = 1;
else		out.group_val = coord_system.group_val;
end
if (~isfield(coord_system,'system_val')),    out.system_val = 1;
else        out.system_val = coord_system.system_val;
end
if (~isfield(coord_system,'datum_val')),     out.datum_val = 221;   % Default to wgs84
else        out.datum_val = coord_system.datum_val;
end
if (~isfield(coord_system,'cilindrical_val')),   out.cilindrical_val = 1;
else        out.cilindrical_val = coord_system.cilindrical_val;
end
if (~isfield(coord_system,'azimuthal_val')),     out.azimuthal_val = 1;
else        out.azimuthal_val = coord_system.azimuthal_val;
end
if (~isfield(coord_system,'conic_val')),         out.conic_val = 1;
else        out.conic_val = coord_system.conic_val;
end
if (~isfield(coord_system,'miscelaneous_val')),  out.miscelaneous_val = 1;
else        out.miscelaneous_val = coord_system.miscelaneous_val;
end
if (~isfield(coord_system,'ProjName')),          out.ProjName = 'Unknown';
else        out.ProjName = coord_system.ProjName;
end
if (~isfield(coord_system,'map_scale_factor')),  out.map_scale_factor = [];
else        out.map_scale_factor = coord_system.map_scale_factor;
end
if (~isfield(coord_system,'system_FE_FN')),      out.system_FE_FN = [];
else        out.system_FE_FN = coord_system.system_FE_FN;
end
if (~isfield(coord_system,'projection')),        out.projection = [];
else        out.projection = coord_system.projection;
end
if (~isfield(coord_system,'ProjParameterValue')), out.ProjParameterValue = [];
else        out.ProjParameterValue = coord_system.ProjParameterValue;
end
if (~isfield(coord_system,'proj_info_txt')),     out.proj_info_txt = 'Nikles';
else		out.proj_info_txt = coord_system.proj_info_txt;
end
if (~isfield(coord_system,'MeasureUnit_val')),   out.MeasureUnit_val = 1;
else        out.MeasureUnit_val = coord_system.MeasureUnit_val;
end
if (~isfield(coord_system,'DegreeFormat1_val')),   out.DegreeFormat1_val = 1;
else        out.DegreeFormat1_val = coord_system.DegreeFormat1_val;
end
if (~isfield(coord_system,'DegreeFormat2_val')),   out.DegreeFormat2_val = 1;
else        out.DegreeFormat2_val = coord_system.DegreeFormat2_val;
end
if (~isfield(coord_system,'is_geog')),   out.is_geog = 1;
else        out.is_geog = coord_system.is_geog;
end

% This is my solution to cat 2 structures. There must be a clever way.
hand_cell = struct2cell(handles);       % Convert handles struct to cell
out_cell = struct2cell(out);            % Convert coord_system struct to cell
both_cell = [hand_cell; out_cell];      % Cat them
names_hand = fieldnames(handles);       % Get handles field names
names_out = fieldnames(out);            % Get coord_system field names
for (i=1:length(names_out))             % Append the "side" to the coord_system field names
	names_out{i} = [names_out{i} side];
end
both_names = [names_hand; names_out];   % Cat the handles and the coord_system field names
handles = cell2struct(both_cell,both_names,1);  % Finaly, rebuild the handles structure.

%-------------------------------------------------------------------------------------
function edit_gridLeft_CB(hObject, handles)
	fname = get(hObject,'String');
	if ~isempty(fname)			% Let the push_gridLeft do all the work
		push_gridLeft_CB(handles.push_gridLeft, handles, fname)
	else
		set(hObject,'String','')
	end

%-------------------------------------------------------------------------------------
function push_gridLeft_CB(hObject, handles, opt)
	if (nargin == 3),   fname = opt;
	else    opt = [];
	end

	if (isempty(opt))    % Otherwise we already know fname from the 4th input argument
		[FileName,PathName] = put_or_get_file(handles,{'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'},'Select GMT grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName,FileName];
	end

	[handles,X,Y,Z,head] = read_gmt_type_grids(handles,fname);
	handles.gridLeft = Z;
	handles.x_min = head(1);        handles.x_max = head(2);
	handles.y_min = head(3);        handles.y_max = head(4);
	handles.one_or_zero = ~head(7); 
	[m,n] = size(handles.gridLeft);
	handles.nRow_in = m;    handles.nCol_in = n;
	set(handles.edit_gridLeft, 'String',fname)
	set(handles.edit_xIncLeft,'String',sprintf('%.12g',head(8)))
	set(handles.edit_yIncLeft,'String',sprintf('%.12g',head(9)))
	set(handles.edit_nColsLeft,'String',sprintf('%d',n))
	set(handles.edit_nRowsLeft,'String',sprintf('%d',m))
	set(handles.edit_nColsRight,'String',sprintf('%d',n))
	set(handles.edit_nRowsRight,'String',sprintf('%d',m))
	handles.gridLeftHead = head;
	if (~handles.by_mirone2grid)
		[pathstr,name] = fileparts(fname);
		fname = [pathstr filesep name '_conv.grd'];     % Proposed output grid name
		set(handles.edit_gridRight,'String',fname)
	end
	guidata(hObject,handles)

%-------------------------------------------------------------------------------------
function edit_xIncRight_CB(hObject, handles)
	dms = 0;
	xx = get(hObject,'String');     val = test_dms(xx);
	if (isempty(val)),		set(hObject, 'String', '');		return,		end
	% If it survived then ...
	if (numel(val) > 1),		dms = 1;	end		% inc given in dd:mm or dd:mm:ss format
	x_inc = 0;
	for (i = 1:numel(val)),		x_inc = x_inc + str2double(val{i}) / (60^(i-1));    end
	
	% Make whatever y_inc given compatible with GMT_grd_RI_verify
	headProj = getappdata(handles.figure1,'headProj');
	if (~isempty(headProj))			% When grid was already projected once
		x_min = headProj(1);		x_max = headProj(2);
	else
		errordlg('Hmmm, unknon error here. Sorry.','Error'),	return
	end
	x_inc = ivan_the_terrible((x_max - x_min), x_inc,2);
	if ~dms         % case of decimal unities
		set(hObject,'String',sprintf('%.12g',x_inc))
		ncol = floor((x_max - x_min) / x_inc + 0.5) + handles.one_or_zero;
	else            % inc was in dd:mm or dd:mm:ss format
		ncol = floor((x_max - x_min) / x_inc + 0.5) + handles.one_or_zero;
		ddmm = dec2deg(x_inc);
		set(hObject,'String',ddmm)
	end
	set(handles.edit_nColsRight,'String',sprintf('%d',ncol))
	handles.dms_xinc = dms;
	guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function edit_yIncRight_CB(hObject, handles)
	dms = 0;
	xx = get(hObject,'String');		val = test_dms(xx);
	if (isempty(val)),		set(hObject, 'String', '');		return,		end
	% If it survived then ...
	if (numel(val) > 1),		dms = 1;	end		% inc given in dd:mm or dd:mm:ss format
	y_inc = 0;
	for i = 1:numel(val),   y_inc = y_inc + str2double(val{i}) / (60^(i-1));    end
	
	% Make whatever y_inc given compatible with GMT_grd_RI_verify
	headProj = getappdata(handles.figure1,'headProj');
	if (~isempty(headProj))			% When grid was already projected once
		y_min = headProj(3);		y_max = headProj(4);
	else
		errordlg('Hmmm, unknon error here. Sorry.','Error'),	return
	end
	y_inc = ivan_the_terrible((y_max - y_min), y_inc,2);
	if ~dms         % case of decimal unities
		set(hObject,'String',sprintf('%.12g',y_inc))
		ncol = floor((y_max - y_min) / y_inc + 0.5) + handles.one_or_zero;
	else            % inc was in dd:mm or dd:mm:ss format
		ncol = floor((y_max - y_min) / y_inc + 0.5) + handles.one_or_zero;
		ddmm = dec2deg(y_inc);
		set(hObject,'String',ddmm)
	end
	set(handles.edit_nRowsRight,'String',ncol)
	handles.dms_yinc = dms;
	guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function edit_nColsRight_CB(hObject, handles)
	xx = get(hObject,'String');
	if (isnan(str2double(xx)))
		set(hObject,'String',handles.nCol_in);    return
	end
	if ( isnan(str2double(get(handles.edit_xIncRight,'String'))) )
		return		% Happens before first projection on grid proj
	end
	headProj = getappdata(handles.figure1,'headProj');
	if (~isempty(headProj))			% When grid was already projected once
		x_min = headProj(1);		x_max = headProj(2);
	else
		x_min = handles.x_min;		x_max = handles.x_max;
	end

	x_inc = ivan_the_terrible((x_max - x_min),round(abs(str2double(xx))),1);
	if handles.dms_xinc         % x_inc was given in dd:mm:ss format
		ddmm = dec2deg(x_inc);
		set(handles.edit_xIncRight,'String',ddmm)
	else                        % x_inc was given in decimal format
		set(handles.edit_xIncRight,'String',sprintf('%.12g',x_inc));
	end

%-------------------------------------------------------------------------------------
function edit_nRowsRight_CB(hObject, handles)
	xx = get(hObject,'String');
	if (isnan(str2double(xx)))
		set(hObject,'String',handles.nRow_in);		return
	end
	if ( isnan(str2double(get(handles.edit_yIncRight,'String'))) )
		return		% Happens before first projection on grid proj
	end
	headProj = getappdata(handles.figure1,'headProj');
	if (~isempty(headProj))			% When grid was already projected once
		y_min = headProj(3);		y_max = headProj(4);
	else
		y_min = handles.y_min;		y_max = handles.y_max;
	end
	
	y_inc = ivan_the_terrible((y_max - y_min),round(abs(str2double(xx))),1);
	if handles.dms_yinc        % y_inc was given in dd:mm:ss format
		ddmm = dec2deg(y_inc);
		set(handles.edit_yIncRight,'String',ddmm)
	else                    % y_inc was given in decimal format
		set(handles.edit_yIncRight,'String',sprintf('%.12g',y_inc));
	end
	guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function push_gridRight_CB(hObject, handles)
	str1 = {'*.grd;*.GRD','netCDF int2 grid format (*.grd,*.GRD)'; '*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select output GMT grid','put','.grd');
	if isequal(FileName,0),		return,		end
	set(handles.edit_gridRight,'String',[PathName FileName])
	guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function edit_optDPI_CB(hObject, handles)
	xx = get(hObject,'String');     dpi = str2double(xx);
	if (isempty(xx) || isnan(dpi))
		set(hObject,'String','');    return
	end
	n = round(dpi/2.54);
	set(handles.edit_nColsRight,'String',sprintf('%d',n))
	set(handles.edit_nRowsRight,'String',sprintf('%d',n))
	x_inc = ivan_the_terrible((handles.x_max - handles.x_min),n,1);
	set(handles.edit_xIncRight,'String',sprintf('%.12g',x_inc));  % With the DPI x_inc = y_inc
	set(handles.edit_yIncRight,'String',sprintf('%.12g',x_inc));

%-------------------------------------------------------------------------------------
function all_datums = load_datums()
% Just load the datums
all_datums = {'Adindan - Burkina Faso', 'Clarke-1880','',-118, -14, 218;
        'Adindan - Cameroon', 'Clarke-1880', '',-134, -2, 210;
        'Adindan - Ethiopia', 'Clarke-1880', '', -165, -11, 206;
        'Adindan - Mali', 'Clarke-1880', '', -123, -20, 220;
        'Adindan - MEAN', 'Clarke-1880', 'MEAN FOR Ethiopia; Sudan', -166, -15, 204;
        'Adindan - Senegal', 'Clarke-1880', '', -128, -18, 224;
        'Adindan - Sudan', 'Clarke-1880', '', -161, -14, 205;
        'Afgooye - Somalia', 'Krassovsky', '', -43, -163, 45;
        'Ain el Abd 1970 - Bahrain', 'International-1924', '', -150, -250, -1;
        'Ain el Abd 1970 - Saudi Arabia', 'International-1924', '', -143, -236, 7;
        'American Samoa 1962 - American Samoa Islands', 'Clarke-1866', '', -115, 118, 426;
        'Anna 1 Astro 1965 - Cocos Islands', 'Australian', '', -491, -22, 435;
        'Antigua Island Astro 1943 - Antigua', 'Clarke-1880', 'Antigua (Leeward Islands)', -270, 13, 62;
        'Arc 1950 - Botswana', 'Clarke-1880', '', -138, -105, -289;
        'Arc 1950 - Burundi', 'Clarke-1880', '', -153, -5, -292;
        'Arc 1950 - Lesotho', 'Clarke-1880', '', -125, -108, -295;
        'Arc 1950 - Malawi', 'Clarke-1880', '', -161, -73, -317;
        'Arc 1950 - MEAN', 'Clarke-1880', 'MEAN FOR Botswana; Lesotho; Malawi; Swaziland; Zaire; Zambia; Zimbabwe', -143, -90, -294;
        'Arc 1950 - Swaziland', 'Clarke-1880', '', -134, -105, -295;
        'Arc 1950 - Zaire', 'Clarke-1880', '', -169, -19, -278;
        'Arc 1950 - Zambia', 'Clarke-1880', '', -147, -74, -283;
        'Arc 1950 - Zimbabwe', 'Clarke-1880', '', -142, -96, -293;
        'Arc 1960 - MEAN', 'Clarke-1880', 'MEAN FOR Kenya; Tanzania', -160, -6, -302;
        'Arc 1960 - Kenya', 'Clarke-1880', '', -157, -2, -299;
        'Arc 1960 - Taanzania', 'Clarke-1880', '', -175, -23, -303;
        'Ascension Island 1958 - Ascension', 'International-1924', '', -205, 107, 53;
        'Astro Beacon E 1945 - Iwo Jima', 'International-1924', '', 145, 75, -272;
        'Astro DOS 71/4 - St Helena Island', 'International-1924', '', -320, 550, -494;
        'Astro Tern Island (FRIG) 1961 - Tern', 'International-1924', '', 114, -116, -333;
        'Astronomical Station 1952 - Marcus Island', 'International-1924', '', 124, -234, -25;
        'Australian Geodetic 1966 - Australia/Tasmania', 'Australian', '', -133, -48, 148;
        'Australian Geodetic 1984 - Australia/Tasmania', 'Australian', '', -134, -48, 149;
        'Ayabelle Lighthouse - Djibouti', 'Clarke-1880', '', -79, -129, 145;
        'Bellevue (IGN) - Efate & Erromango Islands', 'International-1924', '', -127, -769, 472;
        'Bermuda 1957 - Bermuda', 'Clarke-1866', '', -73, 213, 296;
        'Bissau - Guinea-Bissau', 'International-1924', '', -173, 253, 27;
        'Bogota Observatory - Colombia', 'International-1924', '', 307, 304, -318;
        'Bukit Rimpah - Indonesia (Bangka & Belitung Ids)', 'Bessel', '', -384, 664, -48;
        'Camp Area Astro - Antarctica (McMurdo Camp Area)', 'International-1924', '', -104, -129, 239;
        'Campo Inchauspe - Argentina', 'International-1924', '', -148, 136, 90;
        'Canton Astro 1966 - Phoenix Islands', 'International-1924', '', 298, -304, -375;
        'Cape - South Africa', 'Clarke-1880', '', -136, -108, -292;
        'Cape Canaveral - Bahamas/Florida', 'Clarke-1866', '', -2, 151, 181;
        'Carthage - Tunisia', 'Clarke-1880', '', -263, 6, 431;
        'Chatham Island Astro 1971 - New Zealand (Chatham Island)', 'International-1924', '', 175, -38, 113;
        'Chua Astro - Paraguay', 'International-1924', '', -134, 229, -29;
        'Corrego Alegre - Brazil', 'International-1924', '', -206, 172, -6;
        'Dabola - Guinea', 'Clarke-1880', '', -83, 37, 124;
        'Datum 73 - Portugal', 'International-1924', '', -223.237, 110.193, 36.649;
        'Datum Lisboa - Portugal', 'International-1924', '', -304.046, -60.576, 103.640	;
        'Deception Island - Deception Island/Antarctia', 'Clarke-1880', '', 260, 12, -147;
        'Djakarta (Batavia) - Indonesia (Sumatra)', 'Bessel', '', -377, 681, -50;
        'DOS 1968 - New Georgia Islands (Gizo Island)', 'International-1924', '', 230, -199, -752;
        'Easter Island 1967 - Easter Island', 'International-1924', '', 211, 147, 111;
        'Estonia; Coordinate System 1937 - Estonia', 'Bessel', '', 374, 150, 588;
        'European 1950 - Cyprus', 'International-1924', '', -104, -101, -140;
        'European 1950 - Egypt', 'International-1924', '', -130, -117, -151;
        'European 1950 - UK', 'International-1924', 'England; Channel Islands; Scotland; Shetland Islands', -86, -96, -120;
        'European 1950 - UK', 'International-1924', 'England; Ireland; Scotland; Shetland Islands', -86, -96, -120;
        'European 1950 - Finland/Norway', 'International-1924', '', -87, -95, -120;
        'European 1950 - Greece', 'International-1924', '', -84, -95, -130;
        'European 1950 - Iran', 'International-1924', '', -117, -132, -164;
        'European 1950 - Italy (Sardinia)', 'International-1924', '', -97, -103, -120;
        'European 1950 - Italy (Sicily)', 'International-1924', '', -97, -88, -135;
        'European 1950 - Malta', 'International-1924', '', -107, -88, -149;
        'European 1950 - MEAN', 'International-1924', 'MEAN FOR Austria; Belgium; Denmark; Finland;France; W Germany; Gibraltar; Greece; Italy; Luxembourg; Netherlands;Norway; Portugal; Spain; Sweden; Switzerland', -87, -98, -121;
        'European 1950 - MEAN FOR Aus/Dk/Fr/WGerm/Nd/Sch', 'International-1924', 'MEAN FOR Austria; Denmark; France; W Germany; Netherlands; Switzerland', -87, -96, -120;
        'European 1950 - ', 'International-1924', 'MEAN FOR Iraq; Israel; Jordan; Lebanon; Kuwait; Saudi Arabia; Syria', -103, -106, -141;
        'European 1950 - Portugal/Spain', 'International-1924', '', -86.277, -108.879, -120.181;
        'European 1950 - Tunisia', 'International-1924', '', -112, -77, -145;
        'European 1979 - MEAN', 'International-1924', 'MEAN FOR Austria; Finland; Netherlands; Norway; Spain; Sweden; Switzerland', -86, -98, -119;
        'Fort Thomas 1955 - Nevis/St. Kitts (Leeward Islands)', 'Clarke-1880', '', -7, 215, 225;
        'Gan 1970 - Republic of Maldives', 'International-1924', '', -133, -321, 50;
        'Geodetic Datum 1949 - New Zealand', 'International-1924', '', 84, -22, 209;
        'Graciosa Base SW - ', 'International-1924', 'Azores (Faial; Graciosa; Pico; Sao Jorge; Terceira)', -111.789, 154.869, -44.675;
        'Guam 1963 - Guam', 'Clarke-1866', '', -100, -248, 259;
        'Gunung Segara - Indonesia (Kalimantan)', 'Bessel', '', -403, 684, 41;
        'GUX 1 Astro - Guadalcanal Island', 'International-1924', '', 252, -209, -751;
        'Herat North - Afghanistan', 'International-1924', '', -333, -222, 114;
        'Hermannskogel Datum - Croatia-Serbia/Bosnia-Herzegovina', 'Bessel (Namibia)', '', 653, -212, 449;
        'Hjorsey 1955 - Iceland', 'International-1924', '', -73, 46, -86;
        'Hong Kong 1963 - Hong Kong', 'International-1924', '', -156, -271, -189;
        'Hu-Tzu-Shan - Taiwan', 'International-1924', '', -637, -549, -203;
        'Indian - Bangladesh', 'Everest-1830', '', 282, 726, 254;
        'Indian - India/Nepal', 'Everest-1830-Kalianpur', '', 295, 736, 257;
        'Indian - Pakistan', 'Everest-1830-Pakistan', '', 283, 682, 231;
        'Indian 1954 - Thailand', 'Everest-1830', '', 217, 823, 299;
        'Indian 1960 - Vietnam (Con Son Island)', 'Everest-1830', '', 182, 915, 344;
        'Indian 1960 - Vietnam (Near 16N)', 'Everest-1830', '', 198, 881, 317;
        'Indian 1975 - Thailand', 'Everest-1830', '', 210, 814, 289;
        'Indonesian - Indonesia', 'Indonesian', '', -24, -15, 5;
        'Ireland 1965 - Ireland', 'Modified-Airy', '', 506, -122, 611;
        'ISTS 061 Astro 1968 - South Georgia Islands', 'International-1924', '', -794, 119, -298;
        'ISTS 073 Astro 1969 - Diego Garcia', 'International-1924', '', 208, -435, -229;
        'Johnston Island 1961 - Johnston Island', 'International-1924', '', 189, -79, -202;
        'Kandawala - Sri Lanka', 'Everest-1830', '', -97, 787, 86;
        'Kerguelen Island 1949 - Kerguelen', 'International-1924', '', 145, -187, 103;
        'Kertau 1948 - West Malaysia & Singapore', 'Everest-1830-Kertau', '', -11, 851, 5;
        'Kusaie Astro 1951 - Caroline Islands', 'International-1924', '', 647, 1777, -1124;
        'Korean Geodetic System - South Korea', 'GRS-80', '', 0, 0, 0;
        'L.C. 5 Astro 1961 - Cayman Brac Island', 'Clarke-1866', '', 42, 124, 147;
        'Leigon - Ghana', 'Clarke-1880', '', -130, 29, 364;
        'Liberia 1964 - Liberia', 'Clarke-1880', '', -90, 40, 88;
        'Luzon - Philippines (Excluding Mindanao)', 'Clarke-1866', '', -133, -77, -51;
        'Luzon - Philippines (Mindanao)', 'Clarke-1866', '', -133, -79, -72;
        'M''Poraloko - Gabon', 'Clarke-1880', '', -74, -130, 42;
        'Macau - Macau', 'International-1924', '', -203.792, -303.371, -154.020;
        'Mahe 1971 - Mahe Island', 'Clarke-1880', '', 41, -220, -134;
        'Massawa - Ethiopia (Eritrea)', 'Bessel', '', 639, 405, 60;
        'Merchich - Morocco', 'Clarke-1880', '', 31, 146, 47;
        'Midway Astro 1961 - Midway Islands', 'International-1924', '', 912, -58, 1227;
        'Minna - Cameroon', 'Clarke-1880', '', -81, -84, 115;
        'Minna - Nigeria', 'Clarke-1880', '', -92, -93, 122;
        'Montserrat Island Astro 1958 - Montserrat (Leeward Islands)', 'Clarke-1880', '', 174, 359, 365;
        'Nahrwan - Oman (Masirah Island)', 'Clarke-1880', '', -247, -148, 369;
        'Nahrwan - Saudi Arabia', 'Clarke-1880', '', -243, -192, 477;
        'Nahrwan - United Arab Emirates', 'Clarke-1880', '', -249, -156, 381;
        'Naparima BWI - Trinidad & Tobago', 'International-1924', '', -10, 375, 165;
        'NAD 1927 - Alaska (Excluding Aleutian Ids)', 'Clarke-1866', '', -5, 135, 172;
        'NAD 1927 - Alaska (Aleutian Ids East of 180W)', 'Clarke-1866', '', -2, 152, 149;
        'NAD 1927 - Alaska (Aleutian Ids West of 180W)', 'Clarke-1866', '', 2, 204, 105;
        'NAD 1927 - Bahamas (Except San Salvador Id)', 'Clarke-1866', '', -4, 154, 178;
        'NAD 1927 - Bahamas (San Salvador Island)', 'Clarke-1866', '', 1, 140, 165;
        'NAD 1927 - Canada (Alberta; British Columbia)', 'Clarke-1866', '', -7, 162, 188;
        'NAD 1927 - Canada (Manitoba; Ontario)', 'Clarke-1866', '', -9, 157, 184;
        'NAD 1927 - Canada (New Brunswick/Newfoundland/Nova Scotia/Quebec)', 'Clarke-1866', '', -22, 160, 190;
        'NAD 1927 - Canada (Northwest Territories/Saskatchewan)', 'Clarke-1866', '', 4, 159, 188;
        'NAD 1927 - Canada (Yukon)', 'Clarke-1866', '', -7, 139, 181;
        'NAD 1927 - Canal Zone', 'Clarke-1866', '', 0, 125, 201;
        'NAD 1927 - Cuba', 'Clarke-1866', '', -9, 152, 178;
        'NAD 1927 - Greenland (Hayes Peninsula)', 'Clarke-1866', '', 11, 114, 195;
        'NAD 1927 - MEAN1', 'Clarke-1866', 'MEAN FOR Antigua; Barbados; Barbuda; Caicos Islands; Cuba; Dominican Republic; Grand Cayman; Jamaica; Turks Islands', -3, 142, 183;
        'NAD 1927 - MEAN2', 'Clarke-1866', 'MEAN FOR Belize; Costa Rica; El Salvador; Guatemala; Honduras; Nicaragua', 0, 125, 194;
        'NAD 1927 - MEAN FOR Canada', 'Clarke-1866', '', -10, 158, 187;
        'NAD 1927 - MEAN FOR CONUS', 'Clarke-1866', '', -8, 160, 176;
        'NAD 1927 - ', 'Clarke-1866', 'MEAN FOR CONUS (East of Mississippi; River Including Louisiana; Missouri; Minnesota)', -9, 161, 179;
        'NAD 1927 - ', 'Clarke-1866', 'MEAN FOR CONUS (West of Mississippi; River Excluding Louisiana; Minnesota; Missouri)', -8, 159, 175;
        'NAD 1927 - Mexico', 'Clarke-1866', '', -12, 130, 190;
        'NAD 1983 - Alaska (Excluding Aleutian Ids)', 'GRS-80', '', 0, 0, 0;
        'NAD 1983 - Aleutian Ids', 'GRS-80', '', -2, 0, 4;
        'NAD 1983 - Canada', 'GRS-80', '', 0, 0, 0;
        'NAD 1983 - CONUS', 'GRS-80', '', 0, 0, 0;
        'NAD 1983 - Hawaii', 'GRS-80', '', 1, 1, -1;
        'NAD 1983 - Mexico/Central America', 'GRS-80', '', 0, 0, 0;
        'North Sahara 1959 - Algeria', 'Clarke-1880', '', -186, -93, 310;
        'Observ Met 1939 - Azores (Corvo & Flores)', 'International-1924', '', -425, -169, 81;
        'Old Egyptian 1907 - Egypt', 'Helmert-1906', '', -130, 110, -13;
        'Old Hawaiian - Hawaii', 'Clarke-1866', '', 89, -279, -183;
        'Old Hawaiian - Kauai', 'Clarke-1866', '', 45, -290, -172;
        'Old Hawaiian - Maui', 'Clarke-1866', '', 65, -290, -190;
        'Old Hawaiian - MEAN', 'Clarke-1866', 'MEAN FOR Hawaii; Kauai; Maui; Oahu', 61, -285, -181;
        'Old Hawaiian - Oahu', 'Clarke-1866', '', 58, -283, -182;
        'Oman - Oman', 'Clarke-1880', '', -346, -1, 224;
        'ORD SURV 1936 - England', 'Airy', '', 371, -112, 434;
        'ORD SURV 1936 - England/Isle of Man/Wales', 'Airy', '', 371, -111, 434;
        'ORD SURV 1936 - MEAN', 'Airy', 'MEAN FOR England; Isle of Man; Scotland; Shetland Islands; Wales', 375, -111, 431;
        'ORD SURV 1936 - Scotland/Shetland Islands', 'Airy', '', 384, -111, 425;
        'ORD SURV 1936 - Wales', 'Airy', '', 370, -108, 434;
        'Pico de las Nieves - Canary Islands', 'International-1924', '', -307, -92, 127;
        'Pitcairn Astro 1967 - Pitcairn Island', 'International-1924', '', 185, 165, 42;
        'Point 58 - MEAN FOR Burkina Faso & Niger', 'Clarke-1880', '', -106, -129, 165;
        'Pointe Noire 1948 - Congo', 'Clarke-1880', '', -148, 51, -291;
        'Porto Santo - Porto Santo/Madeira Islands', 'International-1924', '', -503.362, -247.604, 312.651;
        'Prov S.A. 1956 - Bolivia', 'International-1924', '', -270, 188, -388;
        'Prov S.A. 1956 - Chile (Northern; Near 19S)', 'International-1924', '', -270, 183, -390;
        'Prov S.A. 1956 - Chile (Southern; Near 43S)', 'International-1924', '', -305, 243, -442;
        'Prov S.A. 1956 - Colombia', 'International-1924', '', -282, 169, -371;
        'Prov S.A. 1956 - Ecuador', 'International-1924', '', -278, 171, -367;
        'Prov S.A. 1956 - Guyana', 'International-1924', '', -298, 159, -369;
        'Prov S.A. 1956 - MEAN', 'International-1924', 'MEAN FOR Bolivia; Chile; Colombia; Ecuador; Guyana; Peru; Venezuela', -288, 175, -376;
        'Prov S.A. 1956 - Peru', 'International-1924', '', -279, 175, -379;
        'Prov S.A. 1956 - Venezuela', 'International-1924', '', -295, 173, -371;
        'Prov S.Chilean 1963 - Chile (Near 53S) (Hito XVIII)', 'International-1924', '', 16, 196, 93;
        'Puerto Rico - Puerto Rico/Virgin Islands', 'Clarke-1866', '', 11, 72, -101;
        'Pulkovo 1942 - Russia', 'Krassovsky', '', 28, -130, -95;
        'Qatar National - Qatar', 'International-1924', '', -128, -283, 22;
        'Qornoq - Greenland (South)', 'International-1924', '', 164, 138, -189;
        'Reunion - Mascarene Islands', 'International-1924', '', 94, -948, -1262;
        'Rome 1940 - Italy (Sardinia)', 'International-1924', '', -225, -65, 9;
        'S-42 (Pulkovo 1942) - Hungary', 'Krassovsky', '', 28, -121, -77;
        'S-42 (Pulkovo 1942) - Poland', 'Krassovsky', '', 23, -124, -82;
        'S-42 (Pulkovo 1942) - Czechoslavakia', 'Krassovsky', '', 26, -121, -78;
        'S-42 (Pulkovo 1942) - Latvia', 'Krassovsky', '', 24, -124, -82;
        'S-42 (Pulkovo 1942) - Kazakhstan', 'Krassovsky', '', 15, -130, -84;
        'S-42 (Pulkovo 1942) - Albania', 'Krassovsky', '', 24, -130, -92;
        'S-42 (Pulkovo 1942) - Romania', 'Krassovsky', '', 28, -121, -77;
        'S-JTSK - Czechoslavakia (Prior 1 JAN 1993)', 'Bessel', '', 589, 76, 480;
        'Santo (DOS) 1965 - Espirito Santo Island', 'International-1924', '', 170, 42, 84;
        'Sao Braz - Azores (Sao Miguel/Santa Maria)', 'International-1924', '', -203, 141, 53;
        'Sao Tome - Sao Tome Island', 'International-1924', '', -176, -650, -82;
        'Sapper Hill 1943 - East Falkland Island', 'International-1924', '', -355, 21, 72;
        'Schwarzeck - Namibia', 'Bessel (Namibia)', '', 616, 97, -251;
        'Selvagem Grande 1938 - Salvagem Islands', 'International-1924', '', -289, -124, 60;
        'Sierra Leone 1960 - Sierra Leone', 'Clarke-1880', '', -88, 4, 101;
        'South-American - Argentina', 'South-American', '', -62, -1, -37;
        'South-American - Bolivia', 'South-American', '', -61, 2, -48;
        'South-American - Brazil', 'South-American', '', -60, -2, -41;
        'South-American - Chile', 'South-American', '', -75, -1, -44;
        'South-American - Colombia', 'South-American', '', -44, 6, -36;
        'South-American - Ecuador', 'South-American', '', -48, 3, -44;
        'South-American - Ecuador (Baltra/Galapagos)', 'South-American', '', -47, 26, -42;
        'South-American - Guyana', 'South-American', '', -53, 3, -47;
        'South-American - MEAN', 'South-American', 'MEAN FOR Argentina; Bolivia; Brazil; Chile; Colombia; Ecuador; Guyana; Paraguay; Peru; Trinidad & Tobago; Venezuela', -57, 1, -41;
        'South-American - Paraguay', 'South-American', '', -61, 2, -33;
        'South-American - Peru', 'South-American', '', -58, 0, -44;
        'South-American - Trinidad & Tobago', 'South-American', '', -45, 12, -33;
        'South-American - Venezuela', 'South-American', '', -45, 8, -33;
        'South Asia - Singapore', 'Modified-Fischer-1960', '', 7, -10, -26;
        'Tananarive Observatory 1925 - Madagascar', 'International-1924', '', -189, -242, -91;
        'Timbalai 1948 - Brunei/E. Malaysia (Sabah Sarawak)', 'Everest-1830-Timbalai', '', -679, 669, -48;
        'Tokyo - Japan', 'Bessel', '', -148, 507, 685;
        'Tokyo - MEAN', 'Bessel', 'MEAN FOR Japan; South Korea; Okinawa', -148, 507, 685;
        'Tokyo - Okinawa', 'Bessel', '', -158, 507, 676;
        'Tokyo - South Korea', 'Bessel', '', -147, 506, 687;
        'Tristan Astro 1968 - Tristao da Cunha', 'International-1924', '', -632, 438, -609;
        'Viti Levu 1916 - Fiji (Viti Levu Island)', 'Clarke-1880', '', 51, 391, -36;
        'Voirol 1960 - Algeria', 'Clarke-1880', '', -123, -206, 219;
        'Wake Island Astro 1952 - Wake Atoll', 'International-1924', '', 276, -57, 149;
        'Wake-Eniwetok 1960 - Marshall Islands', 'Hough', '', 102, 52, -38;
        'WGS 1972 - Global Definition', 'WGS-72', '', 0, 0, 0;
        'WGS 1984 - Global Definition', 'WGS-84', '', 0, 0, 0;
        'Yacare - Uruguay', 'International-1924', '', -155, 171, 37;
        'Zanderij - Suriname', 'International-1924', '', -265, 120, -358};


%-------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(hObject);
	end


% --- Creates and returns a handle to the GUI figure. 
function geog_calculator_LayoutFcn(h1)

set(h1, 'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Geographic Computator',...
'NumberTitle','off',...
'Position',[520 470 650 330],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 306 151 21],...
'Enable','inactive',...
'String','Interactive Conversions',...
'ButtonDownFcn',{@geog_calculator_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','interactive');

uicontrol('Parent',h1, 'Position',[161 306 120 21],...
'Enable','inactive',...
'String','File Conversions',...
'ButtonDownFcn',{@geog_calculator_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','FileConv');

uicontrol('Parent',h1, 'Position',[281 306 120 21],...
'Enable','inactive',...
'String','Grid Conversions',...
'ButtonDownFcn',{@geog_calculator_uiCB,h1,'tab_group_ButtonDownFcn'},...
'Tag','tab_group',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'Call',{@geog_calculator_uiCB,h1,'push_tab_bg_CB'},...
'Enable','inactive',...
'Position',[10 177 631 131],...
'Tag','push_tab_bg');

uicontrol('Parent',h1,'Position',[340 40 301 121],'Style','frame');
uicontrol('Parent',h1,'Position',[15 189 295 101],'Style','frame');
uicontrol('Parent',h1,'Position',[339 188 295 101],'Style','frame');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[90 256 151 21],...
'Style','edit',...
'Tag','edit_xLeft',...
'UserData','interactive');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[29 259 58 15],...
'String','Longitude',...
'Style','text',...
'Tag','h_text_xLeft',...
'UserData','interactive');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[90 231 151 21],...
'Style','edit',...
'Tag','edit_yLeft',...
'UserData','interactive');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[29 234 58 15],...
'String','Latitude',...
'Style','text',...
'Tag','h_text_yLeft',...
'UserData','interactive');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'edit_zLeft_CB'},...
'HorizontalAlignment','left',...
'Position',[90 206 151 21],...
'Style','edit',...
'Tag','edit_zLeft',...
'UserData','interactive');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[29 209 58 15],...
'String','Height',...
'Style','text',...
'Tag','text_heightLeft',...
'UserData','interactive');

uicontrol('Parent',h1,'Position',[10 40 301 121],'Style','frame');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[419 257 151 21],...
'Style','edit',...
'Tag','edit_xRight',...
'UserData','interactive');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[419 232 151 21],...
'Style','edit',...
'Tag','edit_yRight',...
'UserData','interactive');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[419 207 151 21],...
'Style','edit',...
'Tag','edit_zRight',...
'UserData','interactive');

uicontrol('Parent',h1, 'Position',[20 154 120 14],...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.627450980392157],...
'String','Coordinate System',...
'Style','text',...
'Tag','text_CSleft');

uicontrol('Parent',h1, 'Position',[36 281 51 15],...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.627450980392157],...
'String','Input',...
'Style','text',...
'Tag','hTxt_in');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[16 75 289 75],...
'String','Input type',...
'Style','text',...
'Tag','h_txt_info_l');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[358 260 58 15],...
'String','East/West',...
'Style','text',...
'Tag','h_text_xRight',...
'UserData','interactive');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[358 235 58 15],...
'String','North/South',...
'Style','text',...
'Tag','h_text_yRight',...
'UserData','interactive');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[358 210 58 15],...
'String','Height',...
'Style','text',...
'Tag','text_heightRight',...
'UserData','interactive');

uicontrol('Parent',h1, 'Position',[359 154 120 14],...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.627450980392157],...
'String','Coordinate System',...
'Style','text',...
'Tag','text_CSright');

uicontrol('Parent',h1, 'Position',[360 281 51 15],...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.627450980392157],...
'String','Output',...
'Style','text',...
'Tag','hTxt_out');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[346 76 289 75],...
'String','Output type',...
'Style','text',...
'Tag','h_txt_info_r');

uicontrol('Parent',h1,'Position',[324 189 3 101],'Style','frame');
uicontrol('Parent',h1,'Position',[324 11 3 151],'Style','frame');

uicontrol('Parent',h1,'Position',[270 227 23 23],...
'Call',{@geog_calculator_uiCB,h1,'push_fileLeft_CB'},...
'Tag','push_fileLeft',...
'UserData','FileConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'edit_fileLeft_CB'},...
'HorizontalAlignment','left',...
'Position',[30 228 240 22],...
'Style','edit',...
'Tag','edit_fileLeft',...
'UserData','FileConv');

uicontrol('Parent',h1,...
'Call',{@geog_calculator_uiCB,h1,'checkbox_nHeadersLeft_CB'},...
'Position',[30 260 65 15],...
'String','Headers?',...
'Style','checkbox',...
'TooltipString','Are there any header lines in the input file?',...
'Tag','checkbox_nHeadersLeft',...
'UserData','FileConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'edit_nHeadersLeft_CB'},...
'HorizontalAlignment','left',...
'Position',[179 256 31 20],...
'Style','edit',...
'TooltipString','How many?',...
'Tag','edit_nHeadersLeft',...
'UserData','FileConv');

uicontrol('Parent',h1,...
'Call',{@geog_calculator_uiCB,h1,'checkbox_ToggleXYLeft_CB'},...
'Position',[224 257 75 19],...
'String','Toggle x,y',...
'Style','checkbox',...
'TooltipString','Toggle x and y columns',...
'Tag','checkbox_ToggleXYLeft',...
'UserData','FileConv');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[110 260 67 15],...
'String','N? of headers',...
'Style','text',...
'TooltipString','How many?',...
'Tag','text_nHeadersLeft',...
'UserData','FileConv');

uicontrol('Parent',h1,...
'Call',{@geog_calculator_uiCB,h1,'checkbox_ellipsoidHeightsLeft_CB'},...
'Enable','off',...
'Position',[30 204 180 18],...
'String','Ellipsoidal heights in third column?',...
'Style','checkbox',...
'Tag','checkbox_ellipsoidHeightsLeft',...
'UserData','FileConv');

uicontrol('Parent',h1,...
'Call',{@geog_calculator_uiCB,h1,'push_fileRight_CB'},...
'Position',[595 226 23 23],...
'Tag','push_fileRight',...
'UserData','FileConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'edit_fileRight_CB'},...
'HorizontalAlignment','left',...
'Position',[355 227 240 22],...
'Style','edit',...
'Tag','edit_fileRight',...
'UserData','FileConv');

uicontrol('Parent',h1,...
'Call',{@geog_calculator_uiCB,h1,'push_DefCoordLeft_CB'},...
'FontAngle','oblique',...
'Position',[30 47 141 29],...
'String','Define Coordinate System',...
'Tag','push_DefCoordLeft');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'popup_UnitesLeft_CB'},...
'Position',[177 47 91 22],...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_UnitesLeft');

uicontrol('Parent',h1,...
'Call',{@geog_calculator_uiCB,h1,'push_DefCoordRight_CB'},...
'FontAngle','oblique',...
'Position',[360 49 141 29],...
'String','Define Coordinate System',...
'Tag','push_DefCoordRight');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'popup_UnitesRight_CB'},...
'Position',[507 50 91 22],...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_UnitesRight');

uicontrol('Parent',h1,...
'Call',{@geog_calculator_uiCB,h1,'push_left2right_CB'},...
'FontSize',11,...
'FontWeight','bold',...
'ForegroundColor',[0 0.250980392156863 0.501960784313725],...
'Position',[100 7 66 21],...
'String','>>',...
'Tag','push_left2right');

uicontrol('Parent',h1,...
'Call',{@geog_calculator_uiCB,h1,'push_right2left_CB'},...
'FontSize',11,...
'FontWeight','bold',...
'ForegroundColor',[0 0.250980392156863 0.501960784313725],...
'Position',[461 7 66 21],...
'String','<<',...
'Tag','push_right2left');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[97 250 81 21],...
'Style','edit',...
'TooltipString','DX grid spacing',...
'Tag','edit_xIncLeft',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[97 224 81 21],...
'Style','edit',...
'TooltipString','DY grid spacing',...
'Tag','edit_yIncLeft',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'Position',[183 250 50 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_nColsLeft',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'Position',[183 224 50 21],...
'Style','edit',...
'TooltipString','Number of rows in the grid',...
'Tag','edit_nRowsLeft',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[119 272 41 13],...
'String','Spacing',...
'Style','text',...
'Tag','text34',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[182 272 51 13],...
'String','# of lines',...
'Style','text',...
'Tag','text35',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[37 254 55 15],...
'String','X Direction',...
'Style','text',...
'Tag','text36',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[37 229 55 15],...
'String','Y Direction',...
'Style','text',...
'Tag','text37',...
'UserData','GridConv');

uicontrol('Parent',h1, 'Position',[270 195 23 23],...
'Call',{@geog_calculator_uiCB,h1,'push_gridLeft_CB'},...
'Tag','push_gridLeft',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'edit_gridLeft_CB'},...
'HorizontalAlignment','left',...
'Position',[30 196 240 22],...
'Style','edit',...
'Tag','edit_gridLeft',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'edit_xIncRight_CB'},...
'HorizontalAlignment','left',...
'Position',[412 250 81 21],...
'Style','edit',...
'TooltipString','DX grid spacing',...
'Tag','edit_xIncRight',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'edit_yIncRight_CB'},...
'HorizontalAlignment','left',...
'Position',[412 224 81 21],...
'Style','edit',...
'TooltipString','DY grid spacing',...
'Tag','edit_yIncRight',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'edit_nColsRight_CB'},...
'Position',[499 250 50 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_nColsRight',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'edit_nRowsRight_CB'},...
'Position',[499 224 50 21],...
'Style','edit',...
'TooltipString','Number of rows in the grid',...
'Tag','edit_nRowsRight',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[435 272 41 13],...
'String','Spacing',...
'Style','text',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[498 272 51 13],...
'String','# of lines',...
'Style','text',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[354 254 55 15],...
'String','X Direction',...
'Style','text',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[354 229 55 15],...
'String','Y Direction',...
'Style','text',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'Call',{@geog_calculator_uiCB,h1,'push_gridRight_CB'},...
'Position',[595 195 23 23],...
'Tag','push_gridRight',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[355 196 240 22],...
'Style','edit',...
'Tag','edit_gridRight',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@geog_calculator_uiCB,h1,'edit_optDPI_CB'},...
'Position',[572 251 41 20],...
'Style','edit',...
'TooltipString','Set the resolution for the new grid in dots per inch',...
'Tag','edit_optDPI',...
'UserData','GridConv');

uicontrol('Parent',h1,...
'Position',[572 228 65 15],...
'String','Toggle',...
'Style','checkbox',...
'TooltipString','Toggle between pixel and gridline registration ',...
'Tag','checkbox_toggleRegist',...
'UserData','GridConv');

function geog_calculator_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

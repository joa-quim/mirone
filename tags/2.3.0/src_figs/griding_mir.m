function varargout = griding_mir(varargin)
% Wrapper figure to call apropriate interpolation MEX

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
	griding_mir_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center');				% Reposition the window on screen

	dirs = getappdata(0,'MIRONE_DIRS');
	if (isempty(dirs))
		handles.home_dir = pwd;
		handles.last_dir = handles.home_dir;
		handles.work_dir = handles.home_dir;
		f_path = [pwd filesep 'data' filesep];
	else
		handles.home_dir = dirs.home_dir;
		handles.last_dir = dirs.last_dir;
		handles.work_dir = dirs.work_dir;
		f_path = [handles.home_dir filesep 'data' filesep];
	end

	% Import icons
	load([f_path 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_InputFile,'CData',Mfopen_ico)
	clear Mfopen_ico;

	handles.command = cell(50,1);
	handles.x_min = [];				handles.x_max = [];
	handles.y_min = [];				handles.y_max = [];
	handles.x_min_or = [];			handles.x_max_or = [];
	handles.y_min_or = [];			handles.y_max_or = [];
	handles.x_inc = [];				handles.y_inc = [];
	handles.dms_xinc = 0;			handles.dms_yinc = 0;
	handles.IamCompiled = false;
	handles.one_or_zero = 1;		% For Grid Registration grids, which are the most common cases
	handles.hMirFig = [];			% Update this bellow when integrated in Mirone

	% Inactivate the headers parameters. They will be activated by the header checkbox
	set(handles.edit_nHeaders,'Enable','inactive')
	set(handles.popup_binInput,'Visible','off')
	set(handles.edit_binary_ncolumnIn,'Visible','off')  % Those two have to wait until we
	set(handles.push_Help_H,'Visible','off')      % know how to read binary files

	% When called by Mirone varargin must contain: mirone fig handle, "type"
	if ~isempty(varargin)
		if ( length(varargin) == 2 && ishandle(varargin{1}) && ischar(varargin{2}) )
			handles.hMirFig = varargin{1};
			type = varargin{2};
			%handles.IamCompiled = handMir.IamCompiled;
		else
			type = 'surface';		% Default to surface
		end
	else
		type = 'surface';			% Default to surface
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
		set(handles.push_Help_S,'Enable', 'off')
		set(handles.check_Option_F,'Enable', 'off')
	elseif strcmp(type,'triangulate')
		set(hObject,'Name','Triangulate')
		set(handles.popup_GridMethod, 'String', {'Delauny Triangulation';'Minimum Curvature';'Near Neighbor'});
		handles.command{1} = 'triangulate ';
		set(handles.edit_S1_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
		set(handles.popup_S2_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
		set(handles.push_Help_S,'Enable', 'off')
	elseif strcmp(type,'nearneighbor')
		set(hObject,'Name','Nearneighbor')
		%set(handles.popup_GridMethod, 'String', {'Near Neighbor';'Delauny Triangulation';'Minimum Curvature'});
		set(handles.popup_GridMethod, 'String', {'Near Neighbor';'Minimum Curvature'});
		set(handles.check_Option_V,'Enable', 'off')
		handles.command{1} = 'nearneighbor ';
	else			% Defaults to surface
		set(hObject,'Name','Surface')
		%set(handles.popup_GridMethod, 'String', {'Minimum Curvature';'Delauny Triangulation';'Near Neighbor'});
		set(handles.popup_GridMethod, 'String', {'Minimum Curvature';'Near Neighbor'});
		handles.command{1} = 'surface ';
		set(handles.edit_S1_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
		set(handles.popup_S2_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
		set(handles.check_Option_F,'Enable', 'off')
		set(handles.push_Help_S,'Enable', 'off')
	end

	if (~isempty(handles.hMirFig))						% If we know the handle to the calling fig
		handMir = guidata(handles.hMirFig);				% get handles of the calling fig
		handles.last_dir = handMir.last_dir;
		handles.home_dir = handMir.home_dir;
		handles.work_dir = handMir.work_dir;
		handles.IamCompiled = handMir.IamCompiled;
	end

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, [handles.txt_GM handles.txt_IDF handles.txt_GLG handles.txt_FNo])
	%------------- END Pro look (3D) -----------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% -----------------------------------------------------------------------------------
function check_Option_H_CB(hObject, handles)
if get(hObject,'Value')
    if isempty(get(handles.edit_InputFile,'String'))
        handles.command{41} = ' -H';
        set(handles.edit_nHeaders,'Enable','on', 'Backgroundcolor','white')
    else
        errordlg(['Now it''s to late to chose this option. If you realy want to do this, clean ' ...
                'the "Input Data File" box, hit the "Return" key, and start again but ' ...
                'check this option before load the data file.'],'Error')
        handles.command{41} = '';
        set(hObject,'Value',0)
        set(handles.edit_nHeaders,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
    end
else
    handles.command{41} = '';     handles.command{42} = '';
    set(handles.edit_nHeaders,'String','1','Enable','off', 'Backgroundcolor',[.764,.603,.603])
end
guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function edit_nHeaders_CB(hObject, handles)
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
function popup_binInput_CB(hObject, handles)
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
        set(handles.check_Option_H,'Enable','on')
        set(handles.edit_nHeaders,'Enable','on', 'Backgroundcolor','white')
        set(handles.check_Option_H,'Value',0);   set(handles.edit_nHeaders,'String','1');
        set(handles.edit_binary_ncolumnIn,'Visible','off')
        if ~get(handles.check_Option_H,'Value')
            set(handles.edit_nHeaders,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
            set(hObject,'Enable','off')
        end
    case 'binary double'
        set(handles.edit_binary_ncolumnIn,'Visible','on')
        handles.command{38} = ' -bi';  handles.command{41} = '';    handles.command{42} = '';
        set(handles.edit_binary_ncolumnIn,'Enable','on', 'Backgroundcolor','white')
        set(handles.check_Option_H,'Value',0);   set(handles.edit_nHeaders,'String','1');
        set(handles.check_Option_H,'Enable','off')
        set(handles.edit_nHeaders,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
    case 'binary single'
        set(handles.edit_binary_ncolumnIn,'Visible','on')
        handles.command{38} = ' -bis';  handles.command{41} = '';    handles.command{42} = '';
        %set(handles.edit_binary_ncolumnIn,'Enable','on', 'Backgroundcolor','white')
        set(handles.check_Option_H,'Value',0);   set(handles.edit_nHeaders,'String','1');
        set(handles.check_Option_H,'Enable','off')
        set(handles.edit_nHeaders,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
end
guidata(hObject,handles)

%----------------------------------------------------------------------------------------------
function edit_binary_ncolumnIn_CB(hObject, handles)
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
function push_Help_H_CB(hObject, handles)
message = {'Input file has Header record(s). Number of header records can be changed'
           'on the "N? of headers" box. Not used with binary data. This, like all GMT'
           'programs, accept table data input in ASCII or binary data. When using'
           'binary data you must be aware of the fact that GMT has no way of'
           'determining the actual number of columns in the file. You should therefore'
           'pass that information to GMT via the ascii/binary popup menu and its side'
           'box, where the actual number of data columns can be introduced.'};
message_win('create',message,'figname','Help on input data');

%----------------------------------------------------------------------------------------------
function edit_InputFile_CB(hObject, handles, opt)
if (nargin == 3 && ischar(opt))   % OPT is a file name transmited by push_InputFile_CB
    xx = opt;
    hObject = handles.edit_InputFile;    % hObject contained the handle to push_InputFile_CB
else
    xx = get(hObject,'String');
end
if ~isempty(xx)
    handles.command{3} = xx;
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
        set(handles.edit_x_min,'String',val{1});     set(handles.edit_x_max,'String',val{2});
        set(handles.edit_y_min,'String',val{3});     set(handles.edit_y_max,'String',val{4});
        % Until something more inteligent is devised (like using some kind of estatistics to estimate
        % default's Nrow & Ncol) the default value of Nrow = Ncol = 100 will be used.
        x_inc = ivan_the_terrible((x_max - x_min),100,1);			% This will be recomputed in dim_funs()
        y_inc = ivan_the_terrible((y_max - y_min),100,1);
        handles.command{14} = ' -I';    handles.command{16} = '/';
        handles.command{15} = num2str(x_inc,8);   handles.command{17} = num2str(y_inc,8);
        set(handles.edit_x_inc,'String',num2str(x_inc,8));     set(handles.edit_y_inc,'String',num2str(y_inc,8));
        set(handles.edit_Ncols,'String','100');     set(handles.edit_Nrows,'String','100');
    end
else
    % Reset everything to initial state (falta a parte do nearneigh)
    set(handles.edit_x_min,'String','');     set(handles.edit_x_max,'String','');
    set(handles.edit_y_min,'String','');     set(handles.edit_y_max,'String','');
    set(handles.edit_x_inc,'String','');     set(handles.edit_y_inc,'String','');
    set(handles.edit_Ncols,'String','');    set(handles.edit_Nrows,'String','');
    set(handles.check_ToggleXY,'Value',0);
    %set(handles.popup_binInput,'Value',1);
    %set(handles.edit_binary_ncolumnIn,'String','3');
    for i = 3:length(handles.command)   handles.command{i} = '';  end
end 
guidata(handles.figure1, handles)

%----------------------------------------------------------------------------------------------
function push_InputFile_CB(hObject, handles)

	if (~isempty(handles.hMirFig) && ishandle(handles.hMirFig))			% If we know it and it exists
        hand = guidata(handles.hMirFig);		% get handles of the calling fig
	else
        hand = handles;
	end

    [FileName,PathName] = put_or_get_file(hand,{ ...
			'*.dat;*.DAT;*.xyz;*.XYX', 'XYZ files (*.dat,*.DAT,*.xyz,*.XYZ)';'*.*', 'All Files (*.*)'},'Select input data','get');
    if isequal(FileName,0),		return,		end

	handles.command{3} = [PathName FileName];
	set(handles.edit_InputFile, 'String',handles.command{3})
	guidata(hObject, handles);
	edit_InputFile_CB(hObject, handles,[PathName FileName]);

% -----------------------------------------------------------------------------------
function edit_x_min_CB(hObject, handles)
	dim_funs('xMin', hObject, handles)

% -----------------------------------------------------------------------------------
function edit_x_max_CB(hObject, handles)
	dim_funs('xMax', hObject, handles)

% -----------------------------------------------------------------------------------
function edit_x_inc_CB(hObject, handles)
	dim_funs('xInc', hObject, handles)

% -----------------------------------------------------------------------------------
function edit_Ncols_CB(hObject, handles)
	dim_funs('nCols', hObject, handles)

% -----------------------------------------------------------------------------------
function edit_y_min_CB(hObject, handles)
	dim_funs('yMin', hObject, handles)

% -----------------------------------------------------------------------------------
function edit_y_max_CB(hObject, handles)
	dim_funs('yMax', hObject, handles)

% -----------------------------------------------------------------------------------
function edit_y_inc_CB(hObject, handles)
	dim_funs('yInc', hObject, handles)

% -----------------------------------------------------------------------------------
function edit_Nrows_CB(hObject, handles)
	dim_funs('nRows', hObject, handles)

% -----------------------------------------------------------------------------------
function check_ToggleXY_CB(hObject, handles)
if get(hObject,'Value')
    t_xmin = get(handles.edit_x_min,'String');   t_xmax = get(handles.edit_x_max,'String');
    t_ymin = get(handles.edit_y_min,'String');   t_ymax = get(handles.edit_y_max,'String');
    if ~isempty(t_xmin) && ~isempty(t_xmax) && ~isempty(t_ymin) && ~isempty(t_ymax)
        set(handles.edit_x_min,'String',t_ymin);        set(handles.edit_x_max,'String',t_ymax)
        set(handles.edit_y_min,'String',t_xmin);        set(handles.edit_y_max,'String',t_xmax)
        handles.command{6} = t_ymin;       handles.command{8} = t_ymax;
        handles.command{10} = t_xmin;      handles.command{12} = t_xmax;
        handles.command{36} = ' -:';
    else        % There is nothing yet to toggle
        set(hObject,'Value',0)
    end
    t_xinc = get(handles.edit_x_inc,'String');    t_yinc = get(handles.edit_y_inc,'String');
    if ~isempty(t_xinc) && ~isempty(t_yinc)
        set(handles.edit_x_inc,'String',t_yinc);        set(handles.edit_y_inc,'String',t_xinc)
        handles.command{15} = t_yinc;       handles.command{17} = t_xinc;
    end
    t_ncol = get(handles.edit_Ncols,'String');    t_nrow = get(handles.edit_Nrows,'String');
    if ~isempty(t_ncol) && ~isempty(t_nrow)
        set(handles.edit_Ncols,'String',t_nrow);        set(handles.edit_Nrows,'String',t_ncol)
    end
    guidata(hObject,handles)
else
    t_xmin = get(handles.edit_x_min,'String');   t_xmax = get(handles.edit_x_max,'String');
    t_ymin = get(handles.edit_y_min,'String');   t_ymax = get(handles.edit_y_max,'String');
    if ~isempty(t_xmin) && ~isempty(t_xmax) && ~isempty(t_ymin) && ~isempty(t_ymax)
        set(handles.edit_x_min,'String',t_ymin);        set(handles.edit_x_max,'String',t_ymax)
        set(handles.edit_y_min,'String',t_xmin);        set(handles.edit_y_max,'String',t_xmax)
        handles.command{6} = t_ymin;       handles.command{8} = t_ymax;
        handles.command{10} = t_xmin;      handles.command{12} = t_xmax;
        handles.command{36} = '';
    end
    t_xinc = get(handles.edit_x_inc,'String');    t_yinc = get(handles.edit_y_inc,'String');
    if ~isempty(t_xinc) && ~isempty(t_yinc)
        set(handles.edit_x_inc,'String',t_yinc);        set(handles.edit_y_inc,'String',t_xinc)
        handles.command{15} = t_yinc;       handles.command{17} = t_xinc;
    end
    t_ncol = get(handles.edit_Ncols,'String');    t_nrow = get(handles.edit_Nrows,'String');
    if ~isempty(t_ncol) && ~isempty(t_nrow)
        set(handles.edit_Ncols,'String',t_nrow);        set(handles.edit_Nrows,'String',t_ncol)
    end
    guidata(hObject,handles)
end

% -----------------------------------------------------------------------------------
function check_Option_F_CB(hObject, handles)
	if get(hObject,'Value')
		handles.one_or_zero = 0;    handles.command{44} = ' -F';
	else
		handles.one_or_zero = 1;    handles.command{44} = '';
	end
	guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function push_Help_R_F_toggle_CB(hObject, handles)
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
	message_win('create',message,'figname','Help on Grid Line Geometry');

% -----------------------------------------------------------------------------------
function edit_S1_Neighbor_CB(hObject, handles)
	xx = get(hObject,'String');
	if isnan(str2double(xx)) && ~isempty(xx)
		set(hObject, 'String', '');
		errordlg('"Search radius" must be a number','Error');
		return
	end
	if ~isempty(xx)
		handles.command{19} = ' -S';      handles.command{20} = num2str(abs(str2double(xx)));
		guidata(hObject, handles);
	else
		handles.command{19} = '';      handles.command{20} = xx;
		guidata(hObject, handles);
	end

% -----------------------------------------------------------------------------------
function popup_S2_Neighbor_CB(hObject, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	switch str{val};
		case 'minutes',		handles.command{21} = 'm';
		case 'seconds',		handles.command{21} = 'c';
		case 'kilometers',	handles.command{21} = 'k';
		case 'Kilometers',	handles.command{21} = 'K';
		otherwise,			handles.command{21} = '';
	end
	guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function push_Help_S_CB(hObject, handles)
message = {'Sets the search radius in same units as the grid spacing. Use the'
           'optional popup to select in which unities the distance is given.'
           'Selecting kilometers implies that grid limits and spacing are in degrees.'
           'Use uppercase Kilometers if distances should be calculated using great'
           'circles [kilometrs uses flat Earth].'};
message_win('create',message,'figname','Help -S option');

% -----------------------------------------------------------------------------------
function popup_GridMethod_CB(hObject, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case 'Minimum Curvature'
        set(handles.edit_S1_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
        set(handles.popup_S2_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
        set(handles.push_Help_S,'Enable', 'off')
        set(handles.check_Option_F,'Enable', 'off')
        handles.command{19} = '';     handles.command{20} = '';
        handles.command{21} = '';     handles.command{30} = '';
        handles.command{1} = 'surface ';
        handles.type = 'surface';
    case 'Delauny Triangulation'
        set(handles.edit_S1_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
        set(handles.popup_S2_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
        set(handles.push_Help_S,'Enable', 'off')
        set(handles.check_Option_F,'Enable', 'on')
        handles.command{19} = '';     handles.command{20} = '';
        handles.command{21} = '';     handles.command{30} = '';
        handles.command{1} = 'triangulate ';
        handles.type = 'triangulate';
    case 'Near Neighbor'
        set(handles.edit_S1_Neighbor,'Enable', 'on', 'Backgroundcolor','white')
        set(handles.popup_S2_Neighbor,'Enable', 'on', 'Backgroundcolor','white')
        set(handles.push_Help_S,'Enable', 'on')
        set(handles.check_Option_F,'Enable', 'on')
        set(handles.edit_S1_Neighbor,'String','')
        set(handles.popup_S2_Neighbor,'Value',1)
        handles.command{1} = 'nearneighbor ';     handles.command{30} = '';
        handles.type = 'nearneighbor';
end
guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function push_Grid_Options_CB(hObject, handles)
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
function push_OK_CB(hObject, handles)
% I will still use the old technique

	tmp = horzcat(handles.command{1:end});
	[tok,rem] = strtok(tmp);
	out{1} = tok;
	i = 2;
	while (rem)
		[tok, rem] = strtok(rem);
		out{i} = tok;        i = i + 1;
	end

	if (get(handles.check_Option_V,'Value'))
		out{end+1} = '-V';
	end
	
	x_min = get(handles.edit_x_min,'String');   x_max = get(handles.edit_x_max,'String');
	y_min = get(handles.edit_y_min,'String');   y_max = get(handles.edit_y_max,'String');
	if isempty(x_min) || isempty(x_max) || isempty(y_min) || isempty(y_max)
		errordlg('One or more grid limits are empty. Open your yes.','Error'),		return
	end

	nx = str2double(get(handles.edit_Ncols,'String'));
	ny = str2double(get(handles.edit_Nrows,'String'));
	if (isnan(nx) || isnan(ny))      % I think this was already tested, but ...
		errordlg('One (or two) of the grid dimensions are not valid. Do your best.','Error'),	return
	end

	if ( strcmp(handles.type,'nearneighbor') && isempty(get(handles.edit_S1_Neighbor,'String')) )
		errordlg('Must give a value for "Search radius".','Error'),		return
	end

	if (~get(handles.check_ToggleXY,'Val'))
		opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f',handles.x_min, handles.x_max, handles.y_min, handles.y_max);
	else
		opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f',handles.y_min, handles.y_max, handles.x_min, handles.x_max);
	end
	opt_I = ['-I' get(handles.edit_x_inc,'string') '/' get(handles.edit_y_inc,'string')];
	if (handles.IamCompiled),	opt_e = '-e';
	else						opt_e = '';
	end
	out{3} = opt_R;			out{4} = opt_I;

	set(handles.figure1,'Name','COMPUTING')
	switch handles.type
		case 'surface'
% 			[Z,head] = gmtmbgrid_m(out{2:end}, '-Mz');		% NOT READY because it doesn't read from file neither -:
			[Z,head] = surface_m(out{2:end});
			tit = 'surface interpolation';
			set(handles.figure1,'Name','Surface')
		case 'nearneighbor'
			[Z,head] = nearneighbor_m(out{2:end}, opt_e);	% We don't want the last ','
			tit = 'nearneighbor interpolation';
			set(handles.figure1,'Name','Nearneighbor')
	end
	if (isnan(head(5))),	head(5) = min(min(Z));		end	% It happens and needs to be investigated
	if (isnan(head(6))),	head(6) = max(max(Z));		end
	[ny,nx] = size(Z);		clear tmp
	X = linspace(head(1),head(2),nx);       Y = linspace(head(3),head(4),ny);
	tmp.head = head;	tmp.X = X;		tmp.Y = Y;		tmp.name = tit;
	mirone(Z,tmp);

% --------------------------------------------------------------------
function Menu_Help_CB(hObject, handles)
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
		message_win('create',message,'figname','Help on surface');        
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
		message_win('create',message,'figname','Help on triangulate');        
	case 'nearneighbor '
		message = {'nearneighbor reads arbitrarily located (x,y,z[,w]) triples [quadruplets]'
		'from standard input [or xyzfile(s)] and uses a nearest neighbor algorithm'
		'to assign an average value to each node that have one or more points'
		'within a radius centered on the node. The average value is computed'
		'as a weighted mean of the nearest point from each sector inside the'
		'search radius. The weighting function used is w(r) = 1.0 / (1 + d ^ 2),'
		'where d = 3 * r / search_radius and r is distance from the node. This'
		'weight is modulated by the observation points'' weights [if supplied].'};
		message_win('create',message,'figname','Help on nearneighbor');        
end

% -----------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
% Check for "escape"
	handles = guidata(hObject);
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(handles.figure1)
	end

% --------------------------------------------------------------------
function about_window_CB(hObject, handles)
	if (~isempty(handles.hMirFig))
		handMir = guidata(handles.hMirFig);
		if (strcmp(handles.command{1},'surface '))
			about_box(handMir, 'surface Last modified 27 Oct 2010', '2.0.1');
		elseif (strcmp(handles.command{1},'triangulate '))
			about_box(handMir, 'triangulate Last modified 27 Oct 2010', '2.0.1');
		else
			about_box(handMir, 'nearneighbor Last modified 27 Oct 2010', '2.0.1');
		end
	end

% --- Creates and returns a handle to the GUI figure. 
function griding_mir_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','griding_mir',...
'NumberTitle','off',...
'Position',[265 206 391 353],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 103 370 48], 'Style','frame');
uicontrol('Parent',h1, 'Position',[10 43 370 41], 'Style','frame');
uicontrol('Parent',h1, 'Position',[10 279 370 67], 'Style','frame');
uicontrol('Parent',h1, 'Position',[10 166 370 93], 'Style','frame');

uicontrol('Parent',h1,...
'Call',@griding_mir_uiCB,...
'Position',[21 318 80 15],...
'String','Headers?',...
'Style','checkbox',...
'Tooltip','Are there any header lines in the input file?',...
'Tag','check_Option_H');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[170 314 31 20],...
'String','1',...
'Style','edit',...
'Tooltip','How many?',...
'Tag','edit_nHeaders');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'Position',[214 313 91 22],...
'String',{  'ascii'; 'binary double'; 'binary single' },...
'Style','popupmenu',...
'Tooltip','Input data type',...
'Value',1,...
'Tag','popup_binInput');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[304 314 25 20],...
'String','3',...
'Style','edit',...
'Tooltip','Number of columns in binary data',...
'Tag','edit_binary_ncolumnIn');

uicontrol('Parent',h1,...
'BackgroundColor',[0.8313725591 0.81568629 0.784314],...
'Call',@griding_mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[351 314 22 22],...
'String','?',...
'Tag','push_Help_H');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[20 286 330 22],...
'Style','edit',...
'Tag','edit_InputFile');

uicontrol('Parent',h1,...
'Call',@griding_mir_uiCB,...
'Position',[350 285 23 23],...
'Tag','push_InputFile');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[76 219 71 21],...
'Style','edit',...
'Tooltip','X min value',...
'Tag','edit_x_min');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[152 219 71 21],...
'Style','edit',...
'Tooltip','X max value',...
'Tag','edit_x_max');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[228 219 71 21],...
'Style','edit',...
'Tooltip','DX grid spacing',...
'Tag','edit_x_inc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[304 219 65 21],...
'Style','edit',...
'Tooltip','Number of columns in the grid',...
'Tag','edit_Ncols');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[76 193 71 21],...
'Style','edit',...
'Tooltip','Y min value',...
'Tag','edit_y_min');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[152 193 71 21],...
'Style','edit',...
'Tooltip','Y max value',...
'Tag','edit_y_max');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[228 193 71 21],...
'Style','edit',...
'Tooltip','DY grid spacing',...
'Tag','edit_y_inc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[304 193 65 21],...
'Style','edit',...
'Tooltip','Number of columns in the grid',...
'Tag','edit_Nrows');

uicontrol('Parent',h1,...
'Call',@griding_mir_uiCB,...
'Position',[76 172 85 19],...
'String','Toggle x,y',...
'Style','checkbox',...
'Tooltip','Toggle x and y columns',...
'Tag','check_ToggleXY');

uicontrol('Parent',h1,...
'Call',@griding_mir_uiCB,...
'Position',[169 173 115 15],...
'String','Pixel registration',...
'Style','checkbox',...
'Tooltip','Make image type grid',...
'Tag','check_Option_F');

uicontrol('Parent',h1,...
'Call',@griding_mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[309 171 61 18],...
'String','?',...
'Tag','push_Help_R_F_toggle');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Position',[60 117 51 21],...
'Style','edit',...
'Tag','edit_S1_Neighbor');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'Position',[193 116 91 22],...
'String',{''; 'minutes'; 'seconds'; 'kilometers'; 'Kilometers' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_S2_Neighbor');

uicontrol('Parent',h1,...
'Call',@griding_mir_uiCB,...
'Position',[305 115 66 21],...
'String','Help',...
'Tag','push_Help_S');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@griding_mir_uiCB,...
'Position',[20 50 181 22],...
'String',{'Minimum Curvature'; 'Delauny Triangulation'; 'Near Neighbor'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_GridMethod');

uicontrol('Parent',h1,...
'Call',@griding_mir_uiCB,...
'Position',[231 50 141 21],...
'String','Options',...
'Tag','push_Grid_Options');

uicontrol('Parent',h1,...
'Position',[20 10 75 15],...
'String','Verbose',...
'Style','checkbox',...
'Tag','check_Option_V');

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
'Call',@griding_mir_uiCB,...
'Label','Help',...
'Tag','Menu_Help');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[91 242 41 13],...
'String','Min',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[246 242 45 13],...
'String','Spacing',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[314 242 51 13],...
'String','# of lines',...
'Style','text');

uicontrol('Parent',h1, 'Position',[30 75 90 15],...
'String','Griding Method',...
'Tag','txt_GM',...
'Style','text');

uicontrol('Parent',h1, 'Position',[31 338 91 15],...
'String','Input Data File',...
'Tag','txt_IDF',...
'Style','text');

uicontrol('Parent',h1, 'Position',[30 252 125 15],...
'String','Griding Line Geometry',...
'Tag','txt_GLG',...
'Style','text');

uicontrol('Parent',h1, 'Position',[38 142 133 17],...
'FontAngle','italic',...
'FontSize',9,...
'String','For Nearneighbor only',...
'Tag','txt_FNo',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[17 113 40 30],...
'String',{'Search'; 'Radius'},...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[149 120 41 15],...
'String','Unities',...
'Style','text');

uimenu('Parent',h1,...
'Call',@griding_mir_uiCB,...
'Label','About',...
'Tag','about_window');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[101 317 67 15],...
'String','N? of headers',...
'Style','text',...
'Tooltip','How many?');

uicontrol('Parent',h1,...
'Call',@griding_mir_uiCB,...
'FontWeight','bold',...
'Position',[270 7 111 21],...
'String','Compute',...
'Tag','push_OK');

function griding_mir_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

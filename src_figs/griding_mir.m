function varargout = griding_mir(varargin)
% Wrapper figure to call apropriate interpolation MEX

%	Copyright (c) 2004-2018 by J. Luis
%
%             DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%                     Version 2, December 2004
% 
%  Everyone is permitted to copy and distribute verbatim or modified
%  copies of this license document, and changing it is allowed as long
%  as the name is changed.
% 
%             DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%    TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
% 
%   0. You just DO WHAT THE FUCK YOU WANT TO.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id: griding_mir.m 11426 2019-07-04 19:28:17Z j $

	if (nargin > 1 && ischar(varargin{1}))
		gui_CB = str2func(varargin{1});
		[varargout{1:nargout}] = feval(gui_CB,varargin{2:end});
	else
		h = griding_mir_OF(varargin{:});
		if (nargout),	varargout{1} = h;   end
	end

% ---------------------------------------------------------------------------------
function hObject = griding_mir_OF(varargin)

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
	handles.xyz = [];
	handles.one_or_zero = 1;		% For Grid Registration grids, which are the most common cases
	handles.hMirFig = [];			% Update this bellow when integrated in Mirone
	handles.opt_b = '';				% For binary files

	% Inactivate the headers parameters. They will be activated by the header checkbox
	set(handles.edit_nHeaders,'Enable','inactive')
	set(handles.popup_binInput,'Visible','off')
	set(handles.edit_binary_ncolumnIn,'Visible','off')  % Those two have to wait until we
	set(handles.push_Help_H,'Visible','off')			% know how to read binary files

	% When called by Mirone varargin must contain: mirone fig handle, "type"
	if ~isempty(varargin)
		if (numel(varargin) >= 2 && ischar(varargin{2}))
			handles.hMirFig = varargin{1};
			type = varargin{2};
			if (numel(varargin) >= 4)
				handles.x_min = double(varargin{3}(1));		handles.x_max = double(varargin{3}(2));
				handles.y_min = double(varargin{3}(3));		handles.y_max = double(varargin{3}(4));
				handles.xyz = varargin{4};
				% Until something more inteligent is devised (like using some kind of estatistics to estimate
				% default's Nrow & Ncol) the default value of Nrow = Ncol = 100 will be used.
				handles.x_inc = ivan_the_terrible((handles.x_max - handles.x_min),100,1);
				handles.y_inc = ivan_the_terrible((handles.y_max - handles.y_min),100,1);
				set(handles.edit_x_min,'Str',sprintf('%.12g',handles.x_min));
				set(handles.edit_x_max,'Str',sprintf('%.12g',handles.x_max));
				set(handles.edit_y_min,'Str',sprintf('%.12g',handles.y_min));
				set(handles.edit_y_max,'Str',sprintf('%.12g',handles.y_max));
				set(handles.edit_x_inc,'Str',sprintf('%.12g',handles.x_inc));
				set(handles.edit_y_inc,'Str',sprintf('%.12g',handles.y_inc));
				set(handles.edit_Ncols,'Str','100');	set(handles.edit_Nrows,'Str','100');
			end
		else
			type = 'surface';		% Default to surface
		end
	else
		type = 'surface';			% Default to surface
	end
	handles.type = type;

	% Choose the default griding_mir method
	% In Mirone the 'Delauny Triangulation' method is not yet implemented
	set(handles.popup_GridMethod, 'String', {'Minimum Curvature - surface';'Minimum Curvature - mbgrid';'Delauny Triangulation';'Near Neighbor';'Median';'Mean'});
	if (strcmp(type,'surface') || strcmp(type,'mbgrid'))
		if (strcmp(type,'surface'))
			set(hObject,'Name','Surface'),		handles.command{1} = 'surface ';
		else
			set(hObject,'Name','gmtmbgrid'),	handles.command{1} = 'gmtmbgrid ';
			set(handles.popup_GridMethod, 'Value', 2)
		end
		set(handles.push_Help_S,'Enable', 'off')
		set(handles.check_Option_F,'Enable', 'off')
	elseif strcmp(type,'triangulate')
		set(hObject,'Name','Triangulate')
		set(handles.popup_GridMethod, 'Value', 3)
		handles.command{1} = 'triangulate ';
		set(handles.push_Help_S,'Enable', 'off')
	elseif strcmp(type,'nearneighbor')
		set(hObject,'Name','Nearneighbor')
		set(handles.popup_GridMethod, 'Value', 4)
		handles.command{1} = 'nearneighbor ';
	elseif strcmp(type,'blockmedian')
		set(hObject,'Name','Median')
		set(handles.popup_GridMethod, 'Value', 5)
		handles.command{1} = 'blockmedian ';
	elseif strcmp(type,'blockmean')
		set(hObject,'Name','Mean')
		set(handles.popup_GridMethod, 'Value', 6)
		handles.command{1} = 'blockmean ';
	else			% Defaults to surface
		set(hObject,'Name','Surface')
		handles.command{1} = 'surface ';
		set(handles.check_Option_F,'Enable', 'off')
		set(handles.push_Help_S,'Enable', 'off')
	end

	if (~strcmp(type,'nearneighbor'))
		set(handles.edit_S1_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
		set(handles.popup_S2_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
	end

	if (~isempty(handles.hMirFig))						% If we know the handle to the calling fig
		handMir = guidata(handles.hMirFig);				% get handles of the calling fig
		handles.last_dir = handMir.last_dir;
		handles.home_dir = handMir.home_dir;
		handles.work_dir = handMir.work_dir;
	end

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, [handles.txt_GM handles.txt_IDF handles.txt_GLG handles.txt_FNo])
	%------------- END Pro look (3D) -----------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');

	if (nargin)
		do_it = false;
		for (k = 1:nargin)
			if (isa(varargin{k}, 'cell')),	do_it = true;	break,	end
		end
		if (do_it)
			external_drive(handles, 'griding_mir', varargin{k:end})
		end
	end

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
		set(handles.edit_nHeaders,'String','0','Enable','off', 'Backgroundcolor',[.764,.603,.603])
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
	if (~isempty(get(handles.edit_InputFile,'String')))
		errordlg(['Now it''s to late to choose any of these options. If you realy want to do this, clean ' ...
				  'the "Input Data File" box, hit the "Return" key, and start again but first ' ...
				  'choose here the input file format before load the data file.'],'Error')
		return
	end
	val = get(hObject,'Value');     str = get(hObject, 'String');
	switch str{val}
		case 'ascii'
			handles.command{38} = '';     handles.command{39} = '';
			set(handles.edit_binary_ncolumnIn,'String','3')
			set(handles.edit_binary_ncolumnIn,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
			set(handles.check_Option_H,'Enable','on')
			set(handles.edit_nHeaders,'Enable','on', 'Backgroundcolor','white')
			set(handles.check_Option_H,'Value',0);   set(handles.edit_nHeaders,'String','0');
			set(handles.edit_binary_ncolumnIn,'Visible','off')
			if ~get(handles.check_Option_H,'Value')
				set(handles.edit_nHeaders,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
				set(hObject,'Enable','off')
			end
		case 'binary double'
			set(handles.edit_binary_ncolumnIn,'Visible','on')
			handles.command{38} = ' -bi';  handles.command{41} = '';    handles.command{42} = '';
			set(handles.edit_binary_ncolumnIn,'Enable','on', 'Backgroundcolor','white')
			set(handles.check_Option_H,'Value',0);   set(handles.edit_nHeaders,'String','0');
			set(handles.check_Option_H,'Enable','off')
			set(handles.edit_nHeaders,'Enable','off', 'Backgroundcolor',[.764,.603,.603])
		case 'binary single'
			set(handles.edit_binary_ncolumnIn,'Visible','on')
			handles.command{38} = ' -bis';  handles.command{41} = '';    handles.command{42} = '';
			%set(handles.edit_binary_ncolumnIn,'Enable','on', 'Backgroundcolor','white')
			set(handles.check_Option_H,'Value',0);   set(handles.edit_nHeaders,'String','0');
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

	if (nargin == 3 && ischar(opt))		% OPT is a file name transmited by push_InputFile_CB
		xx = opt;
	else
		xx = get(hObject,'String');
	end

	if (isempty(xx))
		% Reset everything to initial state (falta a parte do nearneigh)
		set(handles.edit_x_min,'String','');	set(handles.edit_x_max,'String','');
		set(handles.edit_y_min,'String','');	set(handles.edit_y_max,'String','');
		set(handles.edit_x_inc,'String','');	set(handles.edit_y_inc,'String','');
		set(handles.edit_Ncols,'String','');	set(handles.edit_Nrows,'String','');
		set(handles.check_ToggleXY,'Value',0);
		%set(handles.popup_binInput,'Value',1);
		%set(handles.edit_binary_ncolumnIn,'String','3');
		for (i = 3:length(handles.command)),	handles.command{i} = '';	end
		guidata(handles.figure1, handles)
		return
	end

	[bin, n_column] = guess_file(xx);
	if (isempty(bin))
		errordlg(['Error reading file (probably empty)' xx],'Error'),	return
	end
	opt_b = '';
	if (isa(bin,'struct') || bin ~= 0)				% =====---****** BINARY FILE *******---=====
		if (isa(bin,'struct'))
			bin = guess_bin(bin.nCols, bin.type);	% Ask user to confirm/modify guessing
		else
			bin = guess_bin(false);					% Ask user what's in file
		end
		if (isempty(bin))		% User quit
			set(hObject,'String',''),	return
		elseif (bin.nCols < 3)
			errordlg('Number of data columns vcannot be less than 3 (without inventions)','Error')
			set(hObject,'String',''),	return
		end
		opt_b = sprintf(' -bi%d%c', bin.nCols, bin.type(1));
	end
	handles.opt_b = opt_b;

	% Try to deal with fck spaces
	if (~isempty(strfind(xx, ' '))),		xx = ['"' xx '"'];		end
	handles.command{3} = xx;
	
	d_info = gmtmex(['gmtinfo -C ' xx opt_b]);
	if (numel(d_info.data) < 6)
		errordlg('Your data file does not have at least 3 columns.', 'Error')
	end
	x_min = d_info.data(1);		x_max = d_info.data(2);
	y_min = d_info.data(3);		y_max = d_info.data(4);
	val{1} = sprintf('%0.12g', x_min);			val{2} = sprintf('%0.12g', x_max);
	val{3} = sprintf('%0.12g', y_min);			val{4} = sprintf('%0.12g', y_max);
	val{5} = sprintf('%0.12g', d_info.data(5));	val{6} = sprintf('%0.12g', d_info.data(6));

	handles.command{6} = val{1};    handles.command{7} = '/';
	handles.command{8} = val{2};    handles.command{9} = '/';
	handles.command{10} = val{3};   handles.command{11} = '/';
	handles.command{12} = val{4};   handles.command{5} = ' -R';
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
	guidata(hObject, handles)
	edit_InputFile_CB(handles.edit_InputFile, handles, [PathName FileName]);

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
		handles.command{19} = ' -S';	handles.command{20} = num2str(abs(str2double(xx)));
		guidata(hObject, handles);
	else
		handles.command{19} = '';		handles.command{20} = xx;
		guidata(hObject, handles);
	end

% -----------------------------------------------------------------------------------
function popup_S2_Neighbor_CB(hObject, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	switch str{val}
		case 'degrees',		handles.command{21} = 'd';
		case 'minutes',		handles.command{21} = 'm';
		case 'seconds',		handles.command{21} = 's';
		case 'meters',		handles.command{21} = 'e';
		case 'kilometers',	handles.command{21} = 'k';
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
	switch str{val}
		case {'Minimum Curvature - surface' 'Minimum Curvature - mbgrid'}
			set(handles.edit_S1_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
			set(handles.popup_S2_Neighbor,'Enable', 'off', 'Backgroundcolor',[.764,.603,.603])
			set(handles.push_Help_S,'Enable', 'off')
			set(handles.check_Option_F,'Enable', 'off')
			handles.command{19} = '';     handles.command{20} = '';
			handles.command{21} = '';     handles.command{30} = '';
			if (strcmp(str{val}, 'Minimum Curvature - surface'))
				handles.command{1} = 'surface ';		handles.type = 'surface';
			else
				handles.command{1} = 'gmtmbgrid ';		handles.type = 'gmtmbgrid';
			end
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
			handles.command{1} = 'nearneighbor ';		handles.command{30} = '';
			handles.type = 'nearneighbor';
		case {'Median' 'Mean'}
			set(handles.edit_S1_Neighbor,'Enable', 'on', 'Backgroundcolor','white')
			set(handles.popup_S2_Neighbor,'Enable', 'on','Backgroundcolor','white')
			set(handles.push_Help_S,'Enable', 'on')
			set(handles.check_Option_F,'Enable', 'on')
			set(handles.edit_S1_Neighbor,'String','')
			set(handles.popup_S2_Neighbor,'Value',1)
			if (strcmp(str{val}, 'Median'))
				handles.command{1} = 'blockmedian ';	handles.command{30} = '';
				handles.type = 'blockmedian';
			else
				handles.command{1} = 'blockmean ';		handles.command{30} = '';
				handles.type = 'blockmean';
			end
	end
	guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function push_Grid_Options_CB(hObject, handles)
	switch handles.command{1}
		case {'surface ', 'gmtmbgrid '}
			out = surface_options(handles.command{30});
		case 'triangulate '
			msgbox('Not yet programed')
			out = [];
		case 'nearneighbor '
			out = nearneighbor_options(handles.command{30});
	end
	if (~isempty(out))
		handles.command{30} = [' ' out];
	else
		handles.command{30} = '';
	end
	guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function check_plotPts_CB(hObject, handles)
% ...
	if (isempty(handles.xyz))
		set(hObject, 'Val', 0)
		warndlg('Sorry, but this option is only available when data was transmitted in input.','Warning')
	end

% -----------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% Still using the old technique

	tmp = horzcat(handles.command{1:end});
	[tok,rem] = strtok(tmp);
	out = cell(3,1);		% Maybe too short but will shut up Lint
	out{1} = tok;
	i = 2;
	while (numel(rem) > 1)
		[tok, rem] = strtok(rem);
		if (tok(1) == '"')				% Fck spaces
			ind = strfind(rem, '"');	% Find closing quote
			tok = [tok rem(1:ind(1))];	% Build the name
			rem(1:ind(1)+1) = [];		% Finally, remove the trailing chunk from rem
		end
		out{i} = tok;		i = i + 1;
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

	if (strcmp(handles.type,'nearneighbor') && isempty(get(handles.edit_S1_Neighbor,'String')))
		errordlg('Must give a value for "Search radius".','Error'),		return
	end

	if (~get(handles.check_ToggleXY,'Val'))
		opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g',handles.x_min, handles.x_max, handles.y_min, handles.y_max);
	else
		opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g',handles.y_min, handles.y_max, handles.x_min, handles.x_max);
	end
	opt_I = ['-I' get(handles.edit_x_inc,'string') '/' get(handles.edit_y_inc,'string')];
	out{3} = opt_R;			out{4} = opt_I;

	if (get(handles.check_Option_V,'Value'))
		out{end+1} = '-Vl';
	end

	method = handles.type;

	% We are going to use the eventual presence of a -c<n_cells> option to secretely swapp to gmtmbgrid
	opt_c = ' ';
	for (k = 1:numel(out))
		ind  = strfind(out{k}, '-c');
		if (~isempty(ind))
			opt_c = out{k};		opt_c(2) = 'C';		% it was 'c'
			n_cells = str2double(opt_c(3:end));
			if (n_cells >= 0)
				if (strcmp(method, 'surface'))
					out{k} = sprintf('-M%dc', n_cells);	% Because the equivalent in surface is -M
				end
				break
			end
		end
	end

	set(handles.figure1,'Name','COMPUTING')
	set(handles.figure1,'Pointer','watch');
	switch method
		case 'surface'
			if (~isempty(out{2}))				% Then out{2} is presumably a file name
				[Z,head] = c_surface(out{2:end});
			else
				out{2} = handles.xyz;
				[Z,head] = c_surface(out{2:end});
			end
			tit = 'surface interpolation';
			set(handles.figure1,'Name','Surface')

		case 'nearneighbor'
			if (~isempty(handles.xyz))
				[Z,head] = c_nearneighbor(handles.xyz, out{2:end}, '-N4');	% The -N option must be parameterizable
			else
				[Z,head] = c_nearneighbor(out{2:end}, '-N4');
			end
			tit = 'nearneighbor interpolation';
			set(handles.figure1,'Name','Nearneighbor')

		case 'gmtmbgrid'
 			% Data points were transmitted in input. The must-be-doubles is horrible here
			if (~isempty(handles.xyz))
				if (size(handles.xyz,1) > 3)	% xyz is a cols array
					[Z, head] = gmtmbgrid_m(double(handles.xyz(:,1)), double(handles.xyz(:,2)), ...
						double(handles.xyz(:,3)), opt_I, opt_R, '-Mz', opt_c);
				else
					[Z, head] = gmtmbgrid_m(double(handles.xyz(1,:)), double(handles.xyz(2,:)), ...
						double(handles.xyz(3,:)), opt_I, opt_R, '-Mz', opt_c);
				end
			else
				D = gmtmex(['read -Td ' out{2} handles.opt_b]);
				if (~isa(D.data,'double')),		D.data = double(D.data);	end
				[Z, head] = gmtmbgrid_m(D.data(:,1), D.data(:,2), D.data(:,3), opt_I, opt_R, '-Mz', opt_c);
			end
			Z = single(Z);
			tit = 'mbgrid interpolation';
			set(handles.figure1,'Name','gmtmbgrid')

		case {'blockmedian' 'blockmean' 'triangulate'}
			extra_opt = '-Az';
			if (strncmp(method, 'tri', 3)),		extra_opt = '-G';	end
			if (~isempty(handles.xyz))
				G = gmtmex(sprintf('%s %s %s %s', method, opt_R, opt_I, extra_opt), handles.xyz);
			else
				G = gmtmex(sprintf('%s %s %s %s %s', method, out{2}, opt_R, opt_I, extra_opt));
			end
			Z = G.z;
			head = [G.range G.registration G.inc];
			tit = [method ' interpolation'];
			set(handles.figure1,'Name', method)
	end
	set(handles.figure1,'Pointer','arrow');

	if (isnan(head(5))),	head(5) = min(min(Z));		end	% It happens and needs to be investigated
	if (isnan(head(6))),	head(6) = max(max(Z));		end
	[ny,nx] = size(Z);		clear tmp
	X = linspace(head(1),head(2),nx);       Y = linspace(head(3),head(4),ny);
	tmp.head = head;	tmp.X = X;		tmp.Y = Y;		tmp.name = tit;
	hMirFig = mirone(Z,tmp);
	
	if (get(handles.check_plotPts, 'Val'))
		handMir = guidata(hMirFig);
		h = line('XData',handles.xyz(:,1), 'YData',handles.xyz(:,2), 'Parent', handMir.axes1, ...
			'LineStyle','none', 'Marker','o', 'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0], ...
			'MarkerSize',3,'Tag','Pointpolyline');
		draw_funs(h,'DrawSymbol')			% Set marker's uicontextmenu (tag is very important)
	end

% --------------------------------------------------------------------
function Menu_Help_CB(hObject, handles)
switch handles.command{1}
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

% --- Creates and returns a handle to the GUI figure. 
function griding_mir_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'), 'Position',[265 206 391 353],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','griding_mir',...
'NumberTitle','off',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 103 370 48], 'Style','frame');
uicontrol('Parent',h1, 'Position',[10 43 370 41], 'Style','frame');
uicontrol('Parent',h1, 'Position',[10 279 370 67], 'Style','frame');
uicontrol('Parent',h1, 'Position',[10 166 370 93], 'Style','frame');

uicontrol('Parent',h1, 'Position',[21 318 80 15],...
'Callback',@griding_mir_uiCB,...
'String','Headers?',...
'Style','checkbox',...
'Tooltip','Are there any header lines in the input file?',...
'Tag','check_Option_H');

uicontrol('Parent',h1, 'Position',[170 315 31 18],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'String','0',...
'Style','edit',...
'Tooltip','How many?',...
'Tag','edit_nHeaders');

uicontrol('Parent',h1, 'Position',[214 310 91 23],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'String',{'ascii'; 'binary double'; 'binary single'},...
'Style','popupmenu',...
'Tooltip','Input data type',...
'Value',1,...
'Tag','popup_binInput');

uicontrol('Parent',h1, 'Position',[304 315 25 18],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'String','3',...
'Style','edit',...
'Tooltip','Number of columns in binary data',...
'Tag','edit_binary_ncolumnIn');

uicontrol('Parent',h1, 'Position',[351 314 22 22],...
'BackgroundColor',[0.8313725591 0.81568629 0.784314],...
'Callback',@griding_mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'String','?',...
'Tag','push_Help_H');

uicontrol('Parent',h1, 'Position',[20 286 330 22],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Enter input file name',...
'Tag','edit_InputFile');

uicontrol('Parent',h1, 'Position',[350 285 23 23],...
'Callback',@griding_mir_uiCB,...
'Tooltip','Browse for th input file',...
'Tag','push_InputFile');

uicontrol('Parent',h1, 'Position',[76 219 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','X min value',...
'Tag','edit_x_min');

uicontrol('Parent',h1, 'Position',[152 219 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','X max value',...
'Tag','edit_x_max');

uicontrol('Parent',h1, 'Position',[228 219 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','DX grid spacing',...
'Tag','edit_x_inc');

uicontrol('Parent',h1, 'Position',[304 219 65 21],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Number of columns in the grid',...
'Tag','edit_Ncols');

uicontrol('Parent',h1, 'Position',[76 193 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Y min value',...
'Tag','edit_y_min');

uicontrol('Parent',h1, 'Position',[152 193 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Y max value',...
'Tag','edit_y_max');

uicontrol('Parent',h1, 'Position',[228 193 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','DY grid spacing',...
'Tag','edit_y_inc');

uicontrol('Parent',h1, 'Position',[304 193 65 21],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Number of columns in the grid',...
'Tag','edit_Nrows');

uicontrol('Parent',h1, 'Position',[76 172 85 19],...
'Callback',@griding_mir_uiCB,...
'String','Toggle x,y',...
'Style','checkbox',...
'Tooltip','Toggle x and y columns',...
'Tag','check_ToggleXY');

uicontrol('Parent',h1, 'Position',[169 173 115 15],...
'Callback',@griding_mir_uiCB,...
'String','Pixel registration',...
'Style','checkbox',...
'Tooltip','Make image type grid',...
'Tag','check_Option_F');

uicontrol('Parent',h1, 'Position',[309 171 61 18],...
'Callback',@griding_mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'String','?',...
'Tag','push_Help_R_F_toggle');

uicontrol('Parent',h1, 'Position',[60 117 51 21],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_S1_Neighbor');

uicontrol('Parent',h1, 'Position',[193 116 91 22],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'String',{''; 'degrees'; 'minutes'; 'seconds'; 'meters'; 'kilometers' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_S2_Neighbor');

uicontrol('Parent',h1, 'Position',[305 115 66 21],...
'Callback',@griding_mir_uiCB,...
'String','Help',...
'Tag','push_Help_S');

uicontrol('Parent',h1, 'Position',[20 50 181 22],...
'BackgroundColor',[1 1 1],...
'Callback',@griding_mir_uiCB,...
'String',{'Minimum Curvature'; 'Delauny Triangulation'; 'Near Neighbor'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_GridMethod');

uicontrol('Parent',h1, 'Position',[231 50 141 21],...
'Callback',@griding_mir_uiCB,...
'String','Options',...
'Tag','push_Grid_Options');

uicontrol('Parent',h1, 'Position',[20 10 75 15],...
'String','Verbose',...
'Style','checkbox',...
'Tag','check_Option_V');

uicontrol('Parent',h1, 'Position',[18 224 55 15],...
'Enable','inactive',...
'String','X Direction',...
'Style','text');

uicontrol('Parent',h1, 'Position',[17 198 55 15],...
'Enable','inactive',...
'String','Y Direction',...
'Style','text');

uicontrol('Parent',h1, 'Position',[169 241 41 13],...
'Enable','inactive',...
'String','Max',...
'Style','text');

uimenu('Parent',h1,...
'Callback',@griding_mir_uiCB,...
'Label','Help',...
'Tag','Menu_Help');

uicontrol('Parent',h1, 'Position',[91 242 41 13],...
'Enable','inactive',...
'String','Min',...
'Style','text');

uicontrol('Parent',h1, 'Position',[246 242 45 13],...
'Enable','inactive',...
'String','Spacing',...
'Style','text');

uicontrol('Parent',h1, 'Position',[314 242 51 13],...
'Enable','inactive',...
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

uicontrol('Parent',h1, 'Position',[17 113 40 30],...
'Enable','inactive',...
'String',{'Search'; 'Radius'},...
'Style','text');

uicontrol('Parent',h1, 'Position',[149 120 41 15],...
'Enable','inactive',...
'String','Unities',...
'Style','text');

uicontrol('Parent',h1, 'Position',[101 317 67 15],...
'HorizontalAlignment','left',...
'String','N? of headers',...
'Style','text',...
'Tooltip','How many?');

uicontrol('Parent',h1, 'Position',[110 10 60 15],...
'Callback',@griding_mir_uiCB,...
'String','Plot pts',...
'Tooltip','Overlay data points in solution',...
'Style','checkbox',...
'Tag','check_plotPts');

uicontrol('Parent',h1, 'Position',[270 7 111 21],...
'Callback',@griding_mir_uiCB,...
'FontWeight','bold',...
'String','Compute',...
'Tag','push_OK');

function griding_mir_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

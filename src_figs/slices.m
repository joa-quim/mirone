function varargout = slices(varargin)
% SLICES show individual layers of 3D netCDF file and apply several processing algos to satellite data 

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

% $Id$

% For compiling one need to include the aqua_suppfuns.m aquaPlugin.m files.

	hObject = figure('Tag','figure1','Visible','off');
	slices_LayoutFcn(hObject);
	handles = guihandles(hObject);

	got_a_file_to_start = [];		run_aquaPlugin = false;		handles.hMirFig = [];
	if ( numel(varargin) > 0 && ~ischar(varargin{1}) )	% Expects Mirone handles as first arg
		handMir = varargin{1};
		if (numel(varargin) == 2)
			if ( exist(varargin{2}, 'file') == 2 )		% An input file name
				got_a_file_to_start = varargin{2};
				handles.hMirFig = handMir.hMirFig;		% This fig already has first layer
			end
		end
		handles.home_dir = handMir.home_dir;
		handles.last_dir = handMir.last_dir;
		handles.work_dir = handMir.work_dir;
		handles.DefineEllipsoide = handMir.DefineEllipsoide;	% Potentially need in aquaPlugin
		handles.DefineMeasureUnit = handMir.DefineMeasureUnit;	%				"
		handles.IamCompiled = handMir.IamCompiled;
		handles.path_tmp = handMir.path_tmp;
        d_path = handMir.path_data;
	else
		if (numel(varargin) >= 1)		% File name in input
			if ( exist(varargin{1}, 'file') == 2 )
				got_a_file_to_start = varargin{1};
			end
			if (numel(varargin) == 2)		% Run aquaPlugin after loading input file
				if (ischar(varargin{2}))
					run_aquaPlugin = varargin{2};	% Name of the control script
				else
					run_aquaPlugin = true;			% Search the control script in OPTcontrol.txt
				end
			end
		end
		handles.home_dir = cd;
		handles.last_dir = handles.home_dir;
		handles.work_dir = handles.home_dir;
		handles.DefineEllipsoide = [6378137, 0, 1/298.2572235630];	% Defaults to WGS-84
 		handles.DefineMeasureUnit = 'u';							% Defaults to 'user' units
		handles.path_tmp = [handles.home_dir filesep 'tmp' filesep];
		d_path = [handles.home_dir filesep 'data' filesep];
		% Need to know if "IamCompiled". Since that info is in Mirone handles, we need to find it out here
		try			s.s = which('mirone');			handles.IamCompiled = false;
		catch,		handles.IamCompiled = true;
		end
	end

	% -------------- Import/set icons --------------------------------------------
	load([d_path 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_inputName,'CData',Mfopen_ico)
	set(handles.push_fname1,'CData',Mfopen_ico)
	set(handles.push_fname2,'CData',Mfopen_ico)
	set(handles.push_fname3,'CData',Mfopen_ico)

	% Start by making all 'floating' uicontrols invisible
	handles.floaters = [handles.radio_slope handles.radio_pValue handles.text_boxSize handles.edit_boxSize ...
		handles.check_integDim handles.text_nCells handles.edit_nCells handles.text_periods ...
		handles.edit_periods handles.check_doTrends handles.text_bounds handles.edit_subSetA ...
		handles.edit_subSetB handles.text_scale handles.edit_scale handles.text_What handles.popup_what];
	set(handles.floaters, 'Vis','off')

	% Make a list of visibles on a per operation case
	handles.cases = cell(8,6);
	handles.cases{1,1} = [handles.text_boxSize handles.edit_boxSize handles.check_integDim handles.check_doTrends ...
		handles.text_bounds handles.edit_subSetA handles.edit_subSetB];
	handles.cases{1,2} = {'' '' [151 244 139 21] [262 200 81 21] '' '' ''};		% New positions
	handles.cases{1,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{1,4} = {'on' 'Restrict the analysis to the interior of this polygon'};
	handles.cases{1,5} = {'on' 'Optional second polygon'};
	handles.cases{1,6} = {'on' 'Check against a quality flags file'};

	handles.cases{2,1} = [handles.radio_slope handles.radio_pValue handles.check_integDim handles.text_scale ...
		handles.edit_scale handles.text_bounds handles.edit_subSetA handles.edit_subSetB];
	handles.cases{2,2} = {'' '' [151 198 137 21] '' '' '' '' ''};
	handles.cases{2,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{2,4} = {'on' 'Output file'};
	handles.cases{2,5} = {'off' 'Not Used yet'};
	handles.cases{2,6} = {'on' 'Check against a quality flags file'};

	handles.cases{3,1} = [handles.text_nCells handles.edit_nCells];
	handles.cases{3,2} = {'' ''};
	handles.cases{3,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{3,4} = {'on' 'Output file'};
	handles.cases{3,5} = {'off' 'Not Used'};
	handles.cases{3,6} = {'on' 'Check against a quality flags file'};

	handles.cases{4,1} = [handles.text_periods handles.edit_periods handles.check_integDim handles.text_nCells ...
		handles.edit_nCells handles.text_bounds handles.edit_subSetA handles.edit_subSetB handles.text_What ...
		handles.popup_what];
	handles.cases{4,2} = {[11 250 55 16] [64 248 81 21] [137 198 115 21] '' '' '' '' '' '' ''};
	handles.cases{4,3} = {handles.text_bounds 'Bounds'};	% Set handle String prop to second element
	handles.cases{4,4} = {'on' 'Output file'};
	handles.cases{4,5} = {'on' 'Control file'};
	handles.cases{4,6} = {'on' 'Check against a quality flags file'};

	handles.cases{5,1} = [handles.text_bounds handles.edit_subSetA handles.edit_subSetB handles.text_What ...
		handles.popup_what];
	handles.cases{5,2} = {'' '' '' '' ''};
	handles.cases{5,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{5,4} = {'on' 'Output file'};
	handles.cases{5,5} = {'on' 'External polygon file'};
	handles.cases{5,6} = {'on' 'Check against a quality flags file'};

	handles.cases{6,1} = [handles.text_periods handles.edit_periods handles.text_bounds handles.edit_subSetA ...
		handles.edit_subSetB handles.text_What handles.popup_what];
	handles.cases{6,2} = {[11 250 55 16] [64 248 171 21] '' '' '' '' ''};
	handles.cases{6,3} = {handles.text_bounds 'Bounds'};	% Set handle String prop to second element
	handles.cases{6,4} = {'on' 'Output file'};
	handles.cases{6,5} = {'off' ''};
	handles.cases{6,6} = {'off' ''};

	handles.cases{7,1} = [handles.text_bounds handles.edit_subSetA handles.edit_subSetB];
	handles.cases{7,2} = {'' '' ''};
	handles.cases{7,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{7,4} = {'on' 'Output file'};
	handles.cases{7,5} = {'on' 'To correlate file'};
	handles.cases{7,6} = {'off' ''};

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, handles.text_opt)
	%------------- END Pro look (3D) -----------------------------------------------------

	%------------- Some ToolTips  --------------------------------------------------------
	set(handles.edit_scale,'Tooltip',sprintf('Scale rate of change by this value\n(useful when input data has monthly data)'))
	set(handles.edit_qualFlag,'Tooltip',sprintf(['Threshold quality value. Only values of\n' ...
		'quality >= this will be taken into account.\n' ...
		'NOTE: For MODIS use a NEGATIVE value.\n' ...
		'Than, values are retained if quality <= abs(value)']));
	%-------------------------------------------------------------------------------------
	
	handles.sliceNumber = 0;
	handles.is_coards = false;		% Updated, if the case, in aqua_suppfuns
	handles.useLandPhoto = false;	% To not error elsewhere
	handles.radio_shade = [];		%	"

	% ----- Resize the Fig into a minimal size. If needed, it will be expanded later -----
	handles.full_size = get(hObject, 'Pos');
	handles.short_size = true;
	h = findobj(handles.figure1);
	for (k = 2:numel(h))
		pos = get(h(k),'Pos');
		set(h(k),'Pos', pos - [0 handles.full_size(4)-180 0 0])
	end
	set(hObject, 'Pos', [handles.full_size(1:3) 180])
	move2side(hObject,'right');
	set(hObject,'Visible','on');
	%-------------------------------------------------------------------------------------

	% ----------------------- If we got a file in input ------------------------
	if (~isempty(got_a_file_to_start))
		push_inputName_CB(handles.push_inputName, handles, got_a_file_to_start)
		handles = guidata(handles.figure1);		% Get updated handles
		% And now, if asked for, run the aquaPlugin
		if (run_aquaPlugin),	aquaPlugin(handles, run_aquaPlugin),	end
	end
	% --------------------------------------------------------------------------

	guidata(hObject, handles);
	if (nargout),   varargout{1} = hObject;     end

% -------------------------------------------------------------------------------------
function edit_inputName_CB(hObject, handles)
    fname = get(hObject,'String');
	if (~isempty(fname)),	push_inputName_CB(handles.push_inputName, handles, fname),	end

% -------------------------------------------------------------------------------------
function push_inputName_CB(hObject, handles, opt)
% ...
	if (nargin == 2)		% Direct call
		[FileName, PathName, handles] = put_or_get_file(handles, ...
			{'*.nc;*.NC', 'Data files (*.nc,*.NC)';'*.*', 'All Files (*.*)'},'file','get');
		if isequal(FileName,0),		return,		end
		
	else					% File name on input
		[PathName,FNAME,EXT] = fileparts(opt);
		if (~isempty(PathName)),	PathName = [PathName filesep];	end		% To be coherent with the 'if' branch
		FileName = [FNAME EXT];
	end
	pause(0.01);	handles.fname = [PathName FileName];

	if (exist(handles.fname, 'file') ~= 2)
		errordlg(['File: ' handles.fname ' does not exist.'],'Error')
		handles.fname = [];
		return
	end
	set(handles.edit_inputName,'String',handles.fname)

	% Check if it's a netCDF file before decide what to do
	fid = fopen(handles.fname, 'r');
	ID = fread(fid,4,'*char');      ID = ID';      fclose(fid);
	if (strncmpi(ID,'CDF',3) || strncmpi(ID,'HDF',3)  || strcmpi(ID(2:4),'HDF'))
		s = nc_funs('info',handles.fname);
	else			% Some other format. Get it opened with gdalread in aqua_suppfuns
		warndlg('This file is not netCDF. Expect troubles.','WarnError')
		aqua_suppfuns('forGDAL_hdr',handles)
		return
	end

	ind = strcmp({s.Dimension.Name},'number_of_volumes');
	if (any(ind))			% An ANUGA file. A job for aquamoto
		delete(handles.figure1)
		aquamoto(handles.fname)
		return
	end

	[X,Y,Z,head,misc] = nc_io(handles.fname,'R');
	if (numel(head) == 9 && isfield(misc,'z_id'))
		if (numel(misc.z_dim) <= 2)
			errordlg('This netCDF file is not 3D. Use Mirone directly to read it.','Error')
			return
		end
		handles.nc_info = s;		% Save the nc file info
		handles = aqua_suppfuns('coards_hdr',handles,X,Y,head,misc);
		st = [1 10] / (handles.number_of_timesteps - 1);
		set(handles.slider_layer,'Min',1,'Max',handles.number_of_timesteps,'Val',1,'SliderStep',st) 
		guidata(handles.figure1, handles)
		set([handles.edit_sliceNumber handles.slider_layer],'Enable','on')
		slider_layer_CB(handles.slider_layer, handles)	% Indirect way of saying to display first layer.
		set(handles.text_Info,'String',sprintf('Time steps = %d',handles.number_of_timesteps))
	end

% -------------------------------------------------------------------------------------
function slider_layer_CB(hObject, handles)
	handles.sliceNumber = round(get(handles.slider_layer,'Value')) - 1;
	set(handles.edit_sliceNumber,'String', handles.sliceNumber+1)		% Update slice n box
	if (handles.is_coards)
		aqua_suppfuns('coards_slice', handles)
	elseif (handles.is_otherMultiband)
		aqua_suppfuns('forGDAL_slice', handles)
	end

% -------------------------------------------------------------------------------------
function edit_sliceNumber_CB(hObject, handles)
% REVER, O HANDLES.SLICENUMBER NAO E GUARDADO
	xx = fix(str2double(get(hObject,'String')));		% Make sure its an int
	if (isnan(xx) || xx < 1 || xx > handles.number_of_timesteps)
		handles.sliceNumber = 0;		set(hObject,'String','1')
	else
		set(hObject,'String',xx);		set(handles.slider_layer,'Val',xx)		% Update slider
		handles.sliceNumber = xx - 1;
	end
	slider_layer_CB(handles.slider_layer, handles)		% Update image (also saves handles)	

% -------------------------------------------------------------------------------------
function popup_cases_CB(hObject, handles)
% Take care of moving uicontrols around an setting apropriate values for selected case.

	val = get(hObject, 'Val') - 1;
	set(handles.floaters, 'Vis','off')
	if (val == 0),		return,		end		% The NO choice

	% =============== When first time use, grow figure to full size ==========
	if (handles.short_size)
		h = findobj(handles.figure1);
		set(handles.figure1, 'Vis','off')
		for (k = 2:numel(h))
			pos = get(h(k),'Pos');
			set(h(k),'Pos', pos + [0 handles.full_size(4)-180 0 0])
		end
		set(handles.figure1, 'Pos', handles.full_size)
		move2side(handles.figure1,'right');
		set(handles.figure1, 'Vis','on')
		handles.short_size = false;
		guidata(handles.figure1, handles)
	end
	% =========================================================================

	for (k = 1:numel(handles.cases{val,2}))		% Loop over floaters uicontrols of case 'val' for:
		if ( ~isempty(handles.cases{val,2}{k}) )
			set(handles.cases{val,1}(k), 'pos', handles.cases{val,2}{k})	% repositioning
		end
		set(handles.cases{val,3}{1}, 'Str', handles.cases{val,3}{2})		% Reset uicontrol String
	end

	for (k = 4:6)		% Loop over the 3 file names options
		if (strcmp(handles.cases{val,k}(1), 'on'))
			set(handles.(sprintf('text_fname%d',k-3)), 'Str', handles.cases{val,k}{2}, 'Enable', 'on')
			set(handles.(sprintf('edit_fname%d',k-3)), 'Enable', 'on')
			set(handles.(sprintf('push_fname%d',k-3)), 'Enable', 'on')
		else
			set(handles.(sprintf('text_fname%d',k-3)), 'Enable', 'off')
			set(handles.(sprintf('edit_fname%d',k-3)), 'Enable', 'off')
			set(handles.(sprintf('push_fname%d',k-3)), 'Enable', 'off')
		end
	end

	% Exceptions
	str1 = 'Integration in longitude';
	if (val ~= 1),		str1 = 'Spline interpolation';		end
	set(handles.check_integDim, 'Str', str1)
	if (val == 4)
		set(handles.edit_periods, 'Str', '1:12')
		set(handles.edit_periods,'Tooltip',sprintf(['A vector with the months uppon which the mean\n' ...
			'is to be computed. Example:\n' ...
			'months = 1:12 ==> Computes yearly mean\n' ...
			'months = 6:8  ==> Computes June-July-August seazonal means']));
	else
		set(handles.edit_periods, 'Str', '8')
		set(handles.edit_periods,'Tooltip',sprintf(['Number of days of the composit period (e.g. 3, 8, 30)\n' ...
			'Alternatively can be a vector with the periods to compute.\n' ...
			'For example this will compute 3 days composits of Jan2012\n' ...
			'2012.0:3/365:(2012+29/365)']));
	end
	str = get(handles.text_bounds, 'Str');
	if (str(1) == 'S')		% The 'Subset' case
		set(handles.edit_subSetA,'Str','0','Tooltip','Jump this number of layers from from start')
		set(handles.edit_subSetB,'Str','0','Tooltip','Exclude this number of layers, counting from end')
	else					% The 'Bounds' case
		set(handles.edit_subSetA,'Str','0','Tooltip',sprintf(['MIN value allowed on the Z function\n' ...
			'during the spline interpolation.\n' ...
			'The defaults 0 0 mean NO constraints are applied']))
		set(handles.edit_subSetB,'Str','0','Tooltip',sprintf(['MAX value allowed on the Z function\n' ...
			'during the spline interpolation.\n' ...
			'The defaults 0 0 mean NO constraints are applied']))
	end

	set(handles.cases{val,1}, 'Vis','on')

% -------------------------------------------------------------------------------------
function edit_scale_CB(hObject, handles)
% Get the size of the integration box in "Per cell rate of change"
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || isinf(xx) || xx <= 0),	set(hObject,'String','1'),	end
	set(hObject,'String',sprintf('%d',xx))		% Make sure it is an int

% -------------------------------------------------------------------------------------
function edit_boxSize_CB(hObject, handles)
% Get the size of the integration box in Zonal Means
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || isinf(xx) || xx <= 0),	set(hObject,'String','1'),	end

% -------------------------------------------------------------------------------------
function edit_subSetA_CB(hObject, handles)
% This guy may either accept only ints (the Subset case) or floats
	xx = str2double(get(hObject,'String'));
	str = get(handles.text_bounds, 'Str');
	if (str(1) == 'S')		% The 'Subset' case
		if (isnan(xx) || isinf(xx) || xx < 0),	set(hObject,'String','0'),	return,		end
		set(hObject,'String',sprintf('%d',xx))		% Make sure it is an int
	else					% The 'Bounds' case
		if (isnan(xx) || isinf(xx)),	set(hObject,'String','0'),	end
	end

% -------------------------------------------------------------------------------------
function edit_subSetB_CB(hObject, handles)
% This guy may either accept only ints (the Subset case) or floats
	xx = str2double(get(hObject,'String'));
	str = get(handles.text_bounds, 'Str');
	if (str(1) == 'S')		% The 'Subset' case
		if (isnan(xx) || isinf(xx) || xx < 0),	set(hObject,'String','0'),	return,		end
		set(hObject,'String',sprintf('%d',xx))		% Make sure it is an int
	else					% The 'Bounds' case
		if (isnan(xx)),	set(hObject,'String','0'),	end
	end

% -------------------------------------------------------------------------------------
function edit_fname1_CB(hObject, handles)
    fname = get(hObject,'String');
	if (~isempty(fname)),	push_fname1_CB(handles.push_fname1, handles, fname),	end

% -------------------------------------------------------------------------------------
function push_fname1_CB(hObject, handles, fname)
% ...
	if (nargin == 2),	fname = [];		end
	fname = helper_getFile(handles, fname, 1);	% Last arg is used by callee to know who called it
	if (isempty(fname)),	set(handles.edit_fname1,'Str',''),	return,		end
	set(handles.edit_fname1, 'Str', fname)

% -------------------------------------------------------------------------------------
function edit_fname2_CB(hObject, handles)
    fname = get(hObject,'String');
	if (~isempty(fname)),	push_fname2_CB(handles.push_fname2, handles, fname),	end

% -------------------------------------------------------------------------------------
function push_fname2_CB(hObject, handles, fname)
% ...
	if (nargin == 2),	fname = [];		end
	fname = helper_getFile(handles, fname, 2);	% Last arg is used by callee to know who called it
	if (isempty(fname)),	set(handles.edit_fname2,'Str',''),	return,		end
	set(handles.edit_fname2, 'Str', fname)

% -------------------------------------------------------------------------------------
function edit_fname3_CB(hObject, handles)
    fname = get(hObject,'String');
	if (~isempty(fname)),	push_fname3_CB(handles.push_fname3, handles, fname),	end

% -------------------------------------------------------------------------------------
function push_fname3_CB(hObject, handles, fname)
% ...
	if (nargin == 2),	fname = [];		end
	fname = helper_getFile(handles, fname, 3);	% Last arg is used by callee to know who called it
	if (isempty(fname)),	set(handles.edit_fname3,'Str',''),	return,		end
	set(handles.edit_fname3, 'Str', fname)

% -------------------------------------------------------------------------------------
function edit_periods_CB(hObject, handles)
% Hard to test if user screws. We will check that the string has only numbers, ':' and '/'
	str = get(hObject, 'Str');
	if (isempty(str))
		errordlg('The Periods string is empty. If it resulted from a previous error, please do not insist.','Error')
		return
	end
	str = ddewhite(str);					% Make sure no blanks at extremities
	ind = strfind(str, ' ');
	if (~isempty(ind))
		errordlg('There is one or more empty spaces in the Periods string. That is intolerable.','Error')
		set(hObject, 'Str', ''),	return
	end
	str = strrep(str,':',' ');		str = strrep(str,'/',' ');	% Replace ':' '/' by spaces
	erro = false;		n = 0;
	[t,r] = strtok(str);
	while (~isempty(t) && ~erro )
		if (isnan(str2double(t))),  erro = true;
		else                        [t,r] = strtok(r);
		end
		n = n + 1;
	end
	if (n == 0 || n > 3),		erro = true;	end
	if (erro)
		errordlg('There is an error in the Periods string. It is up to you to find it.','Error')
		set(hObject, 'Str', '')
	end

% -------------------------------------------------------------------------------------
function popup_what_CB(hObject, handles)


% -------------------------------------------------------------------------------------
function edit_qualFlag_CB(hObject, handles)
% Hmm, don't know if anything worth doing here.

% -------------------------------------------------------------------------------------
function radio_slope_CB(hObject, handles)
	if (~get(hObject,'Val')),		set(hObject,'Val', 1),		return,		end
	set(handles.radio_pValue, 'Val', 0)

% -------------------------------------------------------------------------------------
function radio_pValue_CB(hObject, handles)
	if (~get(hObject,'Val')),		set(hObject,'Val', 1),		return,		end
	set(handles.radio_slope, 'Val', 0)

% -------------------------------------------------------------------------------------
function edit_nCells_CB(hObject, handles)
% Get the size of the integration box in Zonal Means
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || isinf(xx) || xx <= 0),	set(hObject,'String','200'),	end
	set(hObject,'String',sprintf('%d',xx))		% Make sure it is an int

% -------------------------------------------------------------------------------------
function push_help_CB(hObject, handles)


% -------------------------------------------------------------------------------------
function push_runPlugin_CB(hObject, handles)
% THIS IS A SPECIAL CALLBACK THAT CALLS A FUNCTION NAMED 'aquaPlugin' THAT MAY RESIDE
% ANYWHERE IN THE PATH WORLD. IT'S UP TO THE USER TO DIFFINE ITS CONTENTS.
%	OK, because 'aquaPlugin' checks for the state of the 'check_plugFun' checkbox to
%	allow executing via the OPTcontrol.txt pointed script, and we don't have one with
%	that name (it's in 'aquamoto'), we will do the next trick.
	val_back = get(handles.check_integDim,'Val');		% Use this one temporarily
	set(handles.check_integDim,'Val',1)
	handles.check_plugFun = handles.check_integDim;
	aquaPlugin(handles)							% That's all it should be needed
	set(handles.check_integDim,'Val',val_back)

% -------------------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
%	This function checks the contents of the various uicontrols and creates a file in the
%	'tmp' directory with the parametrs necessary to run each case. The execution itself
%	is acomplish by a recursive call using that script name as argument.
%
%	Because the number of cases is not equal between this function and the 'aquaPlugin'
%	function, it's necessary to map between the two sets.
%	ATTENTION: any change in AQUAPLUGIN must be reflected here
%
%	cases_plugin_map = {...		% From AQUAPLUGIN
% 		'zonal' ...         % 1 - Compute zonal means
% 		'tvar' ...          % 2 - Compute the Temp time rate of a file with annual means by fit of a straight line (Load entire file in memory)
%       'applyFlags' ...    % 3 - Check against its 3D flags companion and replace values < FLAG to NaN
% 		'yearMeanFlag' ...  % 4 - Compute yearly averages from monthly data but checked against a quality flag file
% 		'polygAVG' ...      % 5 - Compute averages of whatever inside polygons (if any)
% 		'flagsStats' ...    % 6 - Compute per/pixel annual or month counts of pixel values with a quality >= flag
% 		'pass_by_count' ... % 7 - Check the curently active 3D file against a count file
% 		'do_math' ...       % 8 - Perform some basic agebraic operations with the 3D planes
% 		'conv2vtk' ...      % 9 - Convert a 3D netCDF file into a VTK format
% 		'L2_periods' ...    % 10 - Calculate composites of L2 products over fixed periods
% 		'corrcoef' ...      % 11 - Calculate correlation coefficients between 2 3D arrays
% 		};
% 		
% 	cases_this_map = { ...		% THIS FUNCTION
% 		'Zonal Means' ...
% 		'Per cell rate of change ' ...
% 		'Seasonal Averages' ...
% 		'Polygonal Averages'; ...
% 		'L2 MODIS Averages' ...
% 		'Per cell correlation coefficient' ...
% 		};

	val = get(handles.popup_cases, 'Val') - 1;
	if (val == 0),		return,		end		% The NO choice

	% Map the case numbers from this function to the one in 'aquaPlugin'. THEY MUST BE IN SYNC
	map2map(1) = 1;		map2map(2) = 2;		map2map(3) = 3;		map2map(4) = 4;
	map2map(5) = 5;		map2map(6) = 10;	map2map(7) = 11;
	val_in_plugin = map2map(val);

	script_name = sprintf('%sCASE%d_sc.txt', handles.path_tmp, val);
	fid = fopen(script_name, 'w');
	fprintf(fid,'# Control script to drive operations in the aquaPlugin function.\n');
	fprintf(fid,'# First line must start with ''case N'', where N is as explained in aquaPlugin.\n');
	fprintf(fid,'# Character arguments, like file names, must be prefixed with the ''char'' key\n');
	fprintf(fid,'# Empties, even for strings, must be set as []\n');
	fprintf(fid,'# Args for run of CASE %d\n\ncase %d\n', val_in_plugin, val_in_plugin);

	if (val == 1)		% Zonal means
		fprintf(fid,'# The ''DLAT'' variable\n%s\n',get(handles.edit_boxSize,'Str'));
		fprintf(fid,'# INTEG_LON (Logical). If TRUE, integration is done along longitude\n%d\n', get(handles.check_integDim,'Val'));
		fprintf(fid,'# DO_TRENDS If FALSE compute zonal integrations, otherwise compute trends of the zonal integrations (per DLAT)\n');
		fprintf(fid,'%d\n', get(handles.check_doTrends,'Val'));
		comm = '# The ''Subset'' var (example: [0 19] -> (82 90)); ([9 9] -> (91 00)); ([19 0] -> (01 09))';
		helper_writeFile(handles, fid, comm, 'subset')
		comm = '# Optional polygon file delimiting an area where the analysis will be carried on';
		helper_writeFile(handles, fid, comm, 'fname', 1)
		comm = '# Optional second polygon file. If it points to a valid file, this function is called twice and results are subtracted';
		helper_writeFile(handles, fid, comm, 'fname', 2)
		comm = '# name of a netCDF file with quality flags. If not provided no quality check is done.';
		helper_writeFile(handles, fid, comm, 'flags')

	elseif (val == 2)	% tvar
		fprintf(fid,'# The ''slope'' var (logical) indicates if compute slope of linear fit (TRUE) or the p parameter (FALSE)\n');
		fprintf(fid,'%d\n', get(handles.check_integDim,'Val'));
		comm = '# The ''Subset'' var (example: [0 19] -> (82 90)); ([9 9] -> (91 00)); ([19 0] -> (01 09))';
		helper_writeFile(handles, fid, comm, 'subset')
		comm = '# name of a netCDF file with quality flags. If not provided no quality check is done.';
		helper_writeFile(handles, fid, comm, 'flags')
		fprintf(fid,'# The ''spline'' var (logical). To compute spline interpolate the missing values\n');
		fprintf(fid,'%d\n', get(handles.check_integDim,'Val'));
		fprintf(fid,'# The ''scale'' var. Scale the final rate by this value\n%s\n',get(handles.edit_scale,'Str'));
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)

	elseif (val == 3)	% Apply flags
		comm = '# name of a netCDF file with quality flags. If not provided no quality check is done.';
		helper_writeFile(handles, fid, comm, 'flags')
		fprintf(fid,'# The ''Ncells'' variable\n%s\n',get(handles.edit_nCells,'Str'));
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)

	elseif (val == 4)	% Seazonal Averages
		comm = '# The ''Periods'' variable (normaly, the number of days of each period)';
		helper_writeFile(handles, fid, comm, 'periods')
		comm = '# name of a netCDF file with quality flags. If not provided no quality check is done.';
		helper_writeFile(handles, fid, comm, 'flags')
		fprintf(fid,'# The ''Ncells'' variable\n%s\n',get(handles.edit_nCells,'Str'));
		fprintf(fid,'# Not used here but need to set as empty\n[]\n');
		fprintf(fid,'# The ''spline'' var (logical or 1x2 vector). To compute spline interpolate the missing values\n');
		if (get(handles.check_integDim,'Var'))
			if (strcmp(get(handles.edit_subSetA,'Str'),'0') && strcmp(get(handles.edit_subSetB,'Str'),'0'))
				fprintf(fid,'1\n');
			else
				fprintf(fid,sprintf('[%s %s]\n',get(handles.edit_subSetA,'Str'), get(handles.edit_subSetB,'Str')));
			end
		end
		fprintf(fid,'# The ''What'' variable that controls what statistic to compute\n');
		valW = get(handles.popup_what,'Val');		str = get(handles.popup_what,'Str');
		switch (str{valW})
			case MEAN,		fprintf(fid,'0\n');
			case MIN,		fprintf(fid,'1\n');
			case MAX,		fprintf(fid,'2\n');
			otherwise,		fprintf(fid,'0\n');
		end
		comm = '# Name of a file with Lon,Lat locations where to output the entire time series';
		helper_writeFile(handles, fid, comm, 'fname', 2)
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)

	elseif (val == 5)	% Polygonal Averages
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)
		fprintf(fid,'# The ''What'' variable that controls what statistic to compute\n');
		valW = get(handles.popup_what,'Val');		str = get(handles.popup_what,'Str');
		switch (str{valW})
			case MEAN,		fprintf(fid,'mean\n');
			case MEDIAN,	fprintf(fid,'median\n');
			case MIN,		fprintf(fid,'min\n');
			case MAX,		fprintf(fid,'max\n');
			case STD,		fprintf(fid,'std\n');
			otherwise,		fprintf(fid,'mean\n');
		end
		comm = '# The ''fnamepolys'' var - File name of a x,y polygon or a list of polygons. Empty = fish polys from fig';
		helper_writeFile(handles, fid, comm, 'fname', 2)
		comm = '# The ''Subset'' var (example: [0 19] -> (82 90)); ([9 9] -> (91 00)); ([19 0] -> (01 09))';
		helper_writeFile(handles, fid, comm, 'subset')
		comm = '# name of a netCDF file with quality flags. If not provided no quality check is done.';
		helper_writeFile(handles, fid, comm, 'flags')

	elseif (val == 6)	% L2 MODIS Averages
		comm = '# The ''Periods'' variable (normaly, the number of days of each period)';
		helper_writeFile(handles, fid, comm, 'periods')
		valW = get(handles.popup_what,'Val');		str = get(handles.popup_what,'Str');
		switch (str{valW})
			case MEAN,		fprintf(fid,'0\n');
			case MEDIAN,	fprintf(fid,'0\n');
			case MIN,		fprintf(fid,'1\n');
			case MAX,		fprintf(fid,'2\n');
			case STD,		fprintf(fid,'0\n');
			otherwise,		fprintf(fid,'0\n');
		end
		comm = '# A 2 elements vector with the MIN and MAX values allowed on the Z function (default [0 inf])';
		helper_writeFile(handles, fid, comm, 'subset')
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)

	elseif (val == 7)	% Per cell correlation coefficient
		comm = '# Name of the other netCDF file whose correlation with loaded array will be estimated';
		helper_writeFile(handles, fid, comm, 'fname', 2)
		comm = '# The ''Subset'' var (example: [0 19] -> (82 90)); ([9 9] -> (91 00)); ([19 0] -> (01 09))';
		helper_writeFile(handles, fid, comm, 'subset')
		fprintf(fid,'# The ''spline'' - NOT USED YET.\n');
		fprintf(fid,'0\n');
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)
	end
	
	fclose(fid);
	if (~isempty(handles.fname))
		aquaPlugin(handles, script_name);
	end

% -------------------------------------------------------------------------------------
function helper_writeFile(handles, fid, comm, opt1, opt2)
% Write pieces of the script file holding the parametrization of selected case
	fprintf(fid, [comm '\n']);	
	if (strncmp(opt1, 'fname',3))
		fname = get(handles.(sprintf('edit_fname%d',opt2)), 'Str');
		if (isempty(fname)),    fprintf(fid,'[]\n');
		else                    fprintf(fid,'char %s\n', fname);
		end
	elseif (strncmp(opt1, 'subset',3))
		fprintf(fid, sprintf( '[%s %s]\n', get(handles.edit_subSetA,'Str'), get(handles.edit_subSetB,'Str') ));
	elseif (strncmp(opt1, 'flags',3))
		fname = get(handles.edit_fname3, 'Str');
		if (isempty(fname)),    fprintf(fid,'[]\n');
		else                    fprintf(fid,'char %s\n', fname);
		end
		fprintf( fid,sprintf('# The ''quality'' value. Ignored if fname = []\n%s\n', get(handles.edit_qualFlag,'Str')) );
	elseif (strncmp(opt1, 'periods',3))
		fprintf(fid,sprintf('%s\n',get(handles.edit_periods,'Str')));
	end

% -------------------------------------------------------------------------------------
function fname = helper_getFile(handles, fname, n)
% Get and check existence of a input filename
	if (isempty(fname))		% Direct call
		[FileName, PathName] = put_or_get_file(handles, ...
			{'*.nc;*.NC;*.dat', 'Data files (*.nc,*.NC,*.dat)';'*.*', 'All Files (*.*)'},'file','get');
		if isequal(FileName,0),		return,		end
		
	else					% File name on input
		[PathName,FNAME,EXT] = fileparts(fname);
		if (~isempty(PathName)),	PathName = [PathName filesep];	end
		FileName = [FNAME EXT];
	end
	pause(0.01);	fname = [PathName FileName];

	checa = true;
	if (n == 1)		% This is box is used both for read/write. Only first case needs to be checked for existence
		str = get(handles.text_fname1,'Str');
		if (str(1) == 'O')		% Output file. Do not check if file exists
			checa = false;
		end
	end
	if (checa && exist(fname, 'file') ~= 2)
		errordlg(['File: ' fname ' does not exist.'],'Error')
		fname = [];
	end

% ----- Insert those other guys here --------------

%--------------------------------------------------

% -------------------------------------------------------------------------------------
% ------------------- The GUI figure. 
function slices_LayoutFcn(h1)

set(h1, 'Position',[520 310 381 440],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Slices',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback');

uicontrol('Parent',h1, 'Position',[10 396 340 21],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','Enter the name of 3D netCDF file',...
'Tag','edit_inputName');

uicontrol('Parent',h1, 'Position',[350 396 21 21],...
'Call',@slices_uiCB,...
'TooltipString','Browse for a 3D netCDF file',...
'Tag','push_inputName');

uicontrol('Parent',h1, 'Position',[10 370 160 15],...
'String','Scale color to global min/max',...
'Style','checkbox',...
'TooltipString','If checked, color palette is computed using global min/max',...
'Tag','check_globalMinMax');

uicontrol('Parent',h1, 'Position',[0 332 382 3], 'Style','frame');
uicontrol('Parent',h1, 'Position',[5 40 372 190], 'Style','frame');

uicontrol('Parent',h1, 'Position',[10 201 55 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Bounds',...
'Style','text',...
'Tag','text_bounds');

uicontrol('Parent',h1, 'Position',[328 363 44 15],...
'HorizontalAlignment','left',...
'String','Layer nº',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 344 321 17],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'Style','slider',...
'Enable','inactive',...
'Tag','slider_layer');

uicontrol('Parent',h1, 'Position',[155 418 200 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Info',...
'Style','text',...
'Tag','text_Info');

uicontrol('Parent',h1, 'Position',[330 343 40 20],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'String','1',...
'Style','edit',...
'TooltipString','Slice number (to go to a specific one, enter a valid slice number here)',...
'Enable','inactive',...
'Tag','edit_sliceNumber');

uicontrol('Parent',h1, 'Position',[10 418 125 17],...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'String','Input netCDF 3D file',...
'Style','text');

uicontrol('Parent',h1, 'Position',[320 303 51 22],...
'Call',@slices_uiCB,...
'FontSize',10,...
'ForegroundColor',[0 0 1],...
'String','Help',...
'Tag','push_help');

uicontrol('Parent',h1, 'Position',[10 280 361 20],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'String',{' '; 'Zonal Means'; 'Per cell rate of change '; 'Apply quality flags'; 'Seasonal Averages'; 'Polygonal Averages'; ...
	'L2 MODIS Averages'; 'Per cell correlation coefficient'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_cases');

uicontrol('Parent',h1, 'Position',[10 248 55 21],...
'Call',@slices_uiCB,...
'String','Slope',...
'Value',1, ...
'Style','radiobutton',...
'Tag','radio_slope');

uicontrol('Parent',h1, 'Position',[80 248 55 21],...
'Call',@slices_uiCB,...
'String','p value',...
'Style','radiobutton',...
'Tag','radio_pValue');

uicontrol('Parent',h1, 'Position',[151 244 139 21],...
'String','Integration in longitude',...
'Style','checkbox',...
'Value',1,...
'Tag','check_integDim');

uicontrol('Parent',h1, 'Position',[262 232 81 21],...
'String','Do trends',...
'Style','checkbox',...
'Tag','check_doTrends');

uicontrol('Parent',h1, 'Position',[70 244 31 21],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'String','1',...
'Style','edit',...
'Tag','edit_boxSize');

uicontrol('Parent',h1, 'Position',[13 245 55 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Box size',...
'Style','text',...
'Tag','text_boxSize');

uicontrol('Parent',h1, 'Position',[57 199 31 21],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'String','0',...
'Style','edit',...
'Tag','edit_subSetA');

uicontrol('Parent',h1, 'Position',[87 199 31 21],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'String','0',...
'Style','edit',...
'Tag','edit_subSetB');

uicontrol('Parent',h1, 'Position',[242 73 74 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Quality value',...
'Style','text',...
'Tag','text_qualFlag');

uicontrol('Parent',h1, 'Position',[11 232 52 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Periods',...
'Style','text',...
'Tag','text_periods');

uicontrol('Parent',h1, 'Position',[64 230 171 21],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'Style','edit',...
'Tag','edit_periods');

uicontrol('Parent',h1, 'Position',[300 196 71 23],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'String',{'MEAN'; 'MEDIAN'; 'MIN'; 'MAX'; 'STD'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_what');

uicontrol('Parent',h1, 'Position',[259 199 40 16],...
'FontSize',9,...
'HorizontalAlignment','right',...
'String','What?',...
'Style','text',...
'Tag','text_What');

uicontrol('Parent',h1, 'Position',[10 300 161 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Choose processing method',...
'Style','text');

uicontrol('Parent',h1, 'Position',[11 171 280 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Restrict the analysis to the interior of this polygon',...
'Style','text',...
'Tag','text_fname1');

uicontrol('Parent',h1, 'Position',[11 121 151 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Optional second polygon ',...
'Style','text',...
'Tag','text_fname2');

uicontrol('Parent',h1, 'Position',[11 71 184 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Check against a quality flags file',...
'Style','text',...
'Tag','text_fname3');

uicontrol('Parent',h1, 'Position',[150 221 76 16],...
'FontAngle','italic',...
'FontSize',9,...
'String','Optional',...
'Style','text',...
'Tag','text_opt');

uicontrol('Parent',h1, 'Position',[229 249 91 16],...
'FontSize',9,...
'HorizontalAlignment','right',...
'String','Inpaint N cells',...
'Style','text',...
'Tag','text_nCells');

uicontrol('Parent',h1, 'Position',[320 247 51 21],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'String','200',...
'Style','edit',...
'TooltipString','Data gaps groups smaller than this are filled by interpolation',...
'Tag','edit_nCells');

uicontrol('Parent',h1, 'Position',[298 199 38 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Scale',...
'Style','text',...
'Tag','text_scale');

uicontrol('Parent',h1, 'Position',[336 197 31 21],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'String','1',...
'Style','edit',...
'Tag','edit_scale');

uicontrol('Parent',h1, 'Position',[10 150 341 21],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_fname1');

uicontrol('Parent',h1, 'Position',[351 150 21 21],...
'Call',@slices_uiCB,...
'Tag','push_fname1');

uicontrol('Parent',h1, 'Position',[10 100 341 21],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_fname2');

uicontrol('Parent',h1, 'Position',[351 100 21 21],...
'Call',@slices_uiCB,...
'Tag','push_fname2');

uicontrol('Parent',h1, 'Position',[319 71 31 21],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'Style','edit',...
'String','6',...
'Tag','edit_qualFlag');

uicontrol('Parent',h1, 'Position',[10 50 341 21],...
'BackgroundColor',[1 1 1],...
'Call',@slices_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_fname3');

uicontrol('Parent',h1, 'Position',[351 50 21 21],...
'Call',@slices_uiCB,...
'Tag','push_fname3');

uicontrol('Parent',h1, 'Position',[9 9 131 22],...
'Call',@slices_uiCB,...
'FontSize',9,...
'String','Run Plugin function',...
'Tag','push_runPlugin');

uicontrol('Parent',h1, 'Position',[276 8 101 22],...
'Call',@slices_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','Compute',...
'Tag','push_compute');

function slices_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

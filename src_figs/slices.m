function varargout = slices(varargin)
% SLICES show individual layers of 3D netCDF file and apply several processing algos to satellite data 
%
%	To read a file and call aquaPlugin and direct it to use a "control script"
%		slices file.nc 'file_name_of_control_script'
%	To read a file and tell aquaPlugin to search the control script name in the OPTcontrol.txt file:
%		slices('file.nc', 0)

%	Copyright (c) 2004-2018 by J. Luis
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

% $Id: slices.m 10389 2018-04-27 16:13:46Z j $

% For compiling one need to include the aqua_suppfuns.m aquaPlugin.m files.

% To add more 'cases' see instructions in push_compute_CB()

	hObject = figure('Tag','figure1','Visible','off');
	slices_LayoutFcn(hObject);
	handles = guihandles(hObject);

	got_a_file_to_start = [];		run_aquaPlugin = false;		handles.hMirFig = [];	handles.ncVarName = [];
	if (numel(varargin) > 0 && ~ischar(varargin{1}))	% Expects Mirone handles as first arg
		handMir = varargin{1};
		if (numel(varargin) == 2)
			if (exist(varargin{2}, 'file') == 2)		% An input file name
				got_a_file_to_start = varargin{2};
				handles.hMirFig = handMir.hMirFig;		% This fig already has first layer (or not!)
			end
		end
		handles.home_dir = handMir.home_dir;
		handles.last_dir = handMir.last_dir;
		handles.work_dir = handMir.work_dir;
		handles.DefineEllipsoide = handMir.DefineEllipsoide;	% Potentially need in aquaPlugin
		handles.DefineMeasureUnit = handMir.DefineMeasureUnit;	%				"
		handles.IamCompiled = handMir.IamCompiled;
		handles.path_tmp = handMir.path_tmp;
		try		handles.ncVarName = handMir.ncVarName;	end		% May exist only for files with 1 or multiple 3D datasets
		handles.path_data = handMir.path_data;		% Used also to fish the deflation_level in nc_io()
	else
		if (numel(varargin) >= 1)		% File name in input
			if (exist(varargin{1}, 'file') == 2)
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
		handles.path_data = [handles.home_dir filesep 'data' filesep];
		% Need to know if "IamCompiled". Since that info is in Mirone handles, we need to find it out here
		try			s.s = which('mirone');			handles.IamCompiled = false;
		catch,		handles.IamCompiled = true;
		end
	end

	% -------------- Import/set icons --------------------------------------------
	load([handles.path_data 'mirone_icons.mat'],'Mfopen_ico');
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
	handles.cases = cell(10,6);
	handles.cases{1,1} = [handles.text_boxSize handles.edit_boxSize handles.check_integDim handles.check_doTrends ...
	                      handles.text_bounds handles.edit_subSetA handles.edit_subSetB];
	handles.cases{1,2} = {'' '' [151 244 139 21] [262 200 81 21] '' '' ''};		% New positions
	handles.cases{1,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{1,4} = {'on' 'Restrict the analysis to the interior of this polygon'};
	handles.cases{1,5} = {'on' 'Optional second polygon'};
	handles.cases{1,6} = {'on' 'Filter with quality flags file (opt)'};

	handles.cases{2,1} = [handles.radio_slope handles.radio_pValue handles.check_integDim handles.text_scale ...
	                      handles.edit_scale handles.text_bounds handles.edit_subSetA handles.edit_subSetB];
	handles.cases{2,2} = {'' '' [151 198 137 21] '' '' '' '' ''};
	handles.cases{2,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{2,4} = {'on' 'Output file'};
	handles.cases{2,5} = {'on' 'Apply this Land mask (opt)'};
	handles.cases{2,6} = {'on' 'Filter with quality flags file (opt)'};

	handles.cases{3,1} = [handles.text_nCells handles.edit_nCells];
	handles.cases{3,2} = {'' ''};
	handles.cases{3,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{3,4} = {'on' 'Output file'};
	handles.cases{3,5} = {'on' 'Not Used'};
	handles.cases{3,6} = {'on' 'Filter with quality flags file (opt)'};

	handles.cases{4,1} = [handles.text_periods handles.edit_periods handles.check_integDim handles.text_nCells ...
	                      handles.edit_nCells handles.text_bounds handles.edit_subSetA handles.edit_subSetB ...
	                      handles.text_What handles.popup_what];
	handles.cases{4,2} = {[11 250 55 16] [64 248 81 23] [137 198 115 21] '' '' '' '' '' '' ''};
	handles.cases{4,3} = {handles.text_bounds 'Bounds'};	% Set handle String prop to second element
	handles.cases{4,4} = {'on' 'Output file'};
	handles.cases{4,5} = {'on' 'Apply this Land mask (opt)'};		% 'Control xy file (optional)'
	handles.cases{4,6} = {'on' 'Filter with quality flags file (opt)'};

	handles.cases{5,1} = [handles.text_bounds handles.edit_subSetA handles.edit_subSetB handles.text_What ...
	                      handles.popup_what];
	handles.cases{5,2} = {'' '' '' '' ''};
	handles.cases{5,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{5,4} = {'on' 'Output file'};
	handles.cases{5,5} = {'on' 'External polygon file'};
	handles.cases{5,6} = {'on' 'Filter with quality flags file (opt)'};

	handles.cases{6,1} = [handles.text_periods handles.edit_periods handles.text_What handles.popup_what];
	handles.cases{6,2} = {[11 250 55 16] [64 248 201 23] '' ''};
	handles.cases{6,3} = {handles.text_bounds 'Bounds'};	% Set handle String prop to second element
	handles.cases{6,4} = {'on' 'Output file'};
	handles.cases{6,5} = {'on' 'Apply this Land mask (optional)'};
	handles.cases{6,6} = {'off' ''};

	handles.cases{7,1} = [handles.text_periods handles.edit_periods handles.text_nCells ...
	                      handles.edit_nCells handles.text_What handles.popup_what];
	handles.cases{7,2} = {[11 250 55 16] [64 248 81 23] '' '' '' ''};
	handles.cases{7,3} = {handles.text_bounds 'Bounds'};	% Set handle String prop to second element
	handles.cases{7,4} = {'on' 'Output file'};
	handles.cases{7,5} = {'on' 'Apply this Land mask (opt)'};
	handles.cases{7,6} = {'on' 'Filter with quality flags file (opt)'};

	handles.cases{8,1} = [handles.text_bounds handles.edit_subSetA handles.edit_subSetB handles.radio_pValue];
	handles.cases{8,2} = {'' '' '' ''};
	handles.cases{8,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{8,4} = {'on' 'Output file'};
	handles.cases{8,5} = {'on' 'To correlate file'};
	handles.cases{8,6} = {'off' ''};

	handles.cases{9,1} = [handles.text_bounds handles.edit_subSetA handles.edit_subSetB];
	handles.cases{9,2} = {'' '' ''};
	handles.cases{9,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{9,4} = {'off' ''};
	handles.cases{9,5} = {'off' ''};
	handles.cases{9,6} = {'off' ''};

	handles.cases{10,1} = [handles.text_bounds handles.edit_subSetA handles.edit_subSetB];
	handles.cases{10,2} = {'' '' ''};
	handles.cases{10,3} = {handles.text_bounds 'Subset'};	% Set handle String prop to second element
	handles.cases{10,4} = {'on' 'Output file'};
	handles.cases{10,5} = {'off' ''};
	handles.cases{10,6} = {'off' ''};

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
	if (strncmpi(ID,'CDF',3) || strncmpi(ID,'HDF',3) || strcmpi(ID(2:4),'HDF'))
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

	if (~isempty(handles.ncVarName))	% Here we have a request to load a certain variable from a multi-container
		fname = {handles.fname; handles.ncVarName};
	else
		fname = handles.fname;
	end

	[X,Y,Z,head,misc] = nc_io(fname,'R');
	if (numel(head) == 9 && isfield(misc,'z_id'))
		if (numel(misc.z_dim) <= 2)
			errordlg('This netCDF file is not 3D. Use Mirone directly to read it.','Error')
			return
		end
		handles.nc_info = s;		% Save the nc file info
		if (any(strcmp({s.Attribute.Name},'TSU')))		% An NSWING TSUnami file
			handles.IamTSU = true;
		else
			handles.IamTSU = false;
		end
		handles = aqua_suppfuns('coards_hdr',handles,X,Y,head,misc,false);
		st = [1 10] / (handles.number_of_timesteps - 1);
		set(handles.slider_layer,'Min',1,'Max',handles.number_of_timesteps,'Val',1,'SliderStep',st) 
		guidata(handles.figure1, handles)
		set([handles.edit_sliceNumber handles.slider_layer],'Enable','on')
		slider_layer_CB(handles.slider_layer, handles)	% Indirect way of saying to display first layer.
		set(handles.text_Info,'String',sprintf('Time steps = %d',handles.number_of_timesteps))

		% Store the nc_info and z_id in Mirone handles so we can access it from there as well
		handles = guidata(handles.figure1);
		try			% Wrap it for the case hMirFig does not exist or is not valid anymore
			handMir = guidata(handles.hMirFig);
			handMir.netcdf_z_id = handles.netcdf_z_id;
			handMir.nc_info = handles.nc_info;
			handMir.time_z  = handles.time;
			handMir.nLayers = misc.z_dim(1);
			if (isa(fname, 'cell')),	fname = fname{1};	end		% take care when it's a cell
			handMir.grdname = fname;		% At least grid_profiler() needs this for the 3D interpolations
			guidata(handles.hMirFig, handMir)
		end
	end


% -------------------------------------------------------------------------------------
function push_globalMinMax_CB(hObject, handles)
% ...

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
		if (~isempty(handles.cases{val,2}{k}))
			set(handles.cases{val,1}(k), 'pos', handles.cases{val,2}{k})	% repositioning
		end
		set(handles.cases{val,3}{1}, 'Str', handles.cases{val,3}{2})		% Reset uicontrol String
	end

	if (val == 2)
		set(handles.check_average3x3, 'Vis', 'on')
	else
		set(handles.check_average3x3, 'Vis', 'off')
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
	if (val ~= 1)
		set(handles.check_integDim, 'Val', 0, 'Str', 'Spline interpolation')
	else
		set(handles.check_integDim, 'Val', 1, 'Str', 'Integration in longitude')
	end

	if (val == 4)
		set(handles.edit_periods, 'Str', '1:12')
		set(handles.edit_periods,'Tooltip',sprintf(['A vector with the months uppon which the mean\n' ...
			'is to be computed. Example:\n' ...
			'months = 1:12 ==> Computes yearly mean\n' ...
			'months = 6:8  ==> Computes June-July-August seazonal means']));
	elseif (val == 7)
		set(handles.edit_periods, 'Str', '1:1')
		set(handles.edit_periods,'Tooltip',sprintf(['A scalar or vector with the months to compute\n' ...
			'the climatology. Example:\n' ...
			'months = 1:1 ==> Computes January climatology\n' ...
			'months = 6:8 ==> Computes June-July-August seazonal climatology\n\n' ...
			'You can also tell it to compute all 12 climatologies at once with\n' ...
			'All or AllMonths (case insensitive) in which case a 12 layers nc file is computed.']));
	else
		set(handles.edit_periods, 'Str', '8')
		set(handles.edit_periods,'Tooltip',sprintf(['Number of days of the composit period (e.g. 3, 8, 30)\n' ...
			'Alternatively can be a vector with the periods to compute.\n' ...
			'For example this will compute 3 days composits of Jan2012\n' ...
			'2012.0:3/365:(2012+31/365)\n' ...
			'Another option is to provide a vector with the limits boundaries. For ex:\n' ...
			'[1 183 365]\n' ...
			'Says we want two periods: one from [1 183] and the other [184 365]\n' ...
			' \n' ...
			'Periods can also be given with Month names using the standard\n' ...
			'three letters code and the year. The year is needed to know if leap year.\n' ...
			'E.g Jan2013.\n' ...
			'You can also give a month range as in Jan2013:Mar2013\n' ...
			'Or request all months of a year (ex: 2014) with Months2014']));
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
	[pato,f,ext] = fileparts(fname);
	if (isempty(pato))			% Than output file will be written to same dir as input data
		name = get(handles.edit_inputName,'String');
		pato = fileparts(name);
		if (~isempty(pato))
			if (isempty(ext)),	fname = [fname '.nc'];	end
			fname = [pato filesep fname];
		end
	end
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
	if (isempty(fname)),	return,		end
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

	% Case where user sent in a period with month names. E.g. 'Jan2013:Mar2013'
	if ((numel(str) == 7 || numel(str) == 15) && any(strcmpi(str(1),{'J' 'F' 'M' 'A' 'S' 'O' 'N' 'D'})))
		year = str2double(str(4:7));
		be = get_months_doys(year, str(1:3));
		if (numel(str) == 15)
			if (str(8) ~= ':')
				errordlg('Badly formated period string with month names. Must be (e.g.) Jan2013:Mar2013','Error')
				error('Badly formated period string with month names. Must be (e.g.) Jan2013:Mar2013')
			end
			be2 = get_months_doys(year, str(9:11));
			be(2) = be2(2);			be(4) = be2(4);
		end
		str = sprintf('%d:%d:%d', be(1), (be(2)-be(1)), be(2));
		set(hObject, 'Str', str)		% Make it as if user had entered a numeric period
	elseif (strncmpi(str, 'all', 3))	% Allow 'all', 'allmonths', 'all?:?' where ? is a month number. For CLIMATOL only
		if ((numel(str) == 3) || strcmp(str, 'allmonths'))
			set(hObject, 'Str', '101:112')
			str = '101:112';
		end
	elseif (numel(str) == 10 && strncmpi(str,'Months',6))		% A 'Months2012' type request
		year = str2double(str(7:10));
		mon = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
		b = ones(1,13);
		for (k = 1:12)
			be = get_months_doys(year, mon{k});
			b(k+1) = be(2);				% Keep the end of the period
		end
		str = sprintf('[%d %d %d %d %d %d %d %d %d %d %d %d %d]', b);
		set(hObject, 'Str', str)		% Make it as if user had entered a numeric period
	end

	% Make sure that if brackets come in pairs
	if (str(1) == '[' && str(end) ~= ']'),	set(hObject, 'Str', [str ']']);	end
	if (str(1) ~= '[' && str(end) == ']'),	set(hObject, 'Str', ['[' str]);	end

	str = strrep(str,':',' ');		str = strrep(str,'/',' ');	% Replace ':' '/' by spaces
	str = strrep(str,'[','');		str = strrep(str,']','');	% Do not test these to see if they are numeric
	erro = false;		n = 0;
	[t,r] = strtok(str);
	while (~isempty(t) && ~erro)
		if (isnan(str2double(t))),	erro = true;
		else,						[t,r] = strtok(r);
		end
		n = n + 1;
	end
	if (n == 0),	erro = true;	end
	if (erro)
		errordlg('There is an error in the Periods string. It is up to you to find it.','Error')
		set(hObject, 'Str', '')
	end

% -------------------------------------------------------------------------------------
function [out, isleapyear] = get_months_doys(year, month)
% Get first and last DOYs for the MONTH period
% MONTH is named after the typical three letters names (EN or PT).
% OUT is a [1 x 4] vector with the DOYs in the first two elements and with
% decimal year values in the last 2 elements.

	isleapyear = ((~rem(year, 4) && rem(year, 100)) || ~rem(year, 400));
	switch (lower(month(1:3)))
		case 'jan',				out = [1 31];
		case {'feb' 'fev'},		out = [32 59    + isleapyear];
		case 'mar',				out = [60 90]   + isleapyear;
		case {'apr' 'abr'},		out = [91 120]  + isleapyear;
		case {'may' 'mai'},		out = [121 151] + isleapyear;
		case 'jun',				out = [152 181] + isleapyear;
		case 'jul',				out = [182 212] + isleapyear;
		case {'aug' 'ago'},		out = [213 243] + isleapyear;
		case {'sep' 'set'},		out = [244 273] + isleapyear;
		case {'oct' 'out'},		out = [274 304] + isleapyear;
		case 'nov',				out = [305 334] + isleapyear;
		case {'dec' 'dez'},		out = [335 365] + isleapyear;
		otherwise
			errordlg('Wrong month 3 letters name', 'Error')
			error('Wrong month 3 letters name')
	end

	% Now compute also the decimal year values corresponding to the above period
	out(3) = year + (out(1)-1) / (365 + isleapyear);
	out(4) = year + (out(2)-1) / (365 + isleapyear);

% -------------------------------------------------------------------------------------
function popup_what_CB(hObject, handles)
% ...
	if (get(hObject, 'Val') == 2)
		val = get(handles.popup_cases, 'Val') - 1;
		if (val == 4 || val == 6 || val == 7)
			warndlg(sprintf(['Sorry, option too comlicated.\n\n   (YES, too complicated)\n\n' ...
			'Try to calculate the median WHITHOUT having all\n' ...
			'the data in memory and you''ll see what I mean)'],'Warning'))
			set(hObject, 'Val', 1)
		end
	end

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
	if (isnan(xx) || isinf(xx) || xx < 0),	xx = 0;		end
	set(hObject,'String',sprintf('%d',xx))		% Make sure it is an int

% -------------------------------------------------------------------------------------
function push_help_CB(hObject, handles)


% -------------------------------------------------------------------------------------
function push_runPlugin_CB(hObject, handles)
% THIS IS A SPECIAL CALLBACK THAT CALLS A FUNCTION NAMED 'aquaPlugin' THAT MAY RESIDE
% ANYWHERE IN THE PATH WORLD. IT'S UP TO THE USER TO DEFINE ITS CONTENTS.
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
%	'tmp' directory with the parameters necessary to run each case. The execution itself
%	is accomplished by a recursive call using that script name as argument.
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
% 		'pass_by_count' ... % 7 - Check the currently active 3D file against a count file
% 		'do_math' ...       % 8 - Perform some basic algebraic operations with the 3D planes
% 		'conv2vtk' ...      % 9 - Convert a 3D netCDF file into a VTK format
% 		'L2_periods' ...    % 10 - Calculate composites of L2 products over fixed periods
% 		'corrcoef' ...      % 11 - Calculate correlation coefficients between 2 3D arrays
%		'count_blooms' ...	% 12 - Count chlor_a blooms
% 		};
% 		
% 	cases_this_map = { ...		% THIS FUNCTION
%		'Zonal Means' ...
%		'Per cell rate of change' ...
%		'Apply quality flags' ...
%		'Seasonal Averages' ...
%		'Polygonal Averages' ...
%		'L2 MODIS Averages' ...
%		'Climatologies (monthly)' ...
%		'Per cell correlation coefficient' ...
%		'Count non-NaNs cells' ...
%		'Find Chlor_a blooms' ...
% 		};

	val = get(handles.popup_cases, 'Val') - 1;
	if (val == 0),		return,		end		% The NO choice

	% Map the case numbers from this function to the one in 'aquaPlugin'. THEY MUST BE IN SYNC
	map2map(1) = 1;		map2map(2) = 2;		map2map(3) = 3;		map2map(4) = 4;
	map2map(5) = 5;		map2map(6) = 10;	map2map(7) = 4;		map2map(8) = 11;
	map2map(9) = 8;		map2map(10)= 12;
	val_in_plugin = map2map(val);

	script_name = sprintf('%sCASE%d_sc.txt', handles.path_tmp, val);
	fid = fopen(script_name, 'w');
	fprintf(fid,'# Control script to drive operations in the aquaPlugin function.\n');
	fprintf(fid,'# First line must start with ''case N'', where N is as explained in aquaPlugin.\n');
	fprintf(fid,'# Character arguments, like file names, must be prefixed with the ''char'' key\n');
	fprintf(fid,'# Empties, even for strings, must be set as []\n');
	fprintf(fid,'# Args for run of CASE %d\n\ncase %d\n', val_in_plugin, val_in_plugin);

	if (val == 1)		% Zonal means -- CASE 1
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

	elseif (val == 2)	% tvar -- CASE 2
		fprintf(fid,'# The ''slope'' var (logical) indicates if compute slope of linear fit (TRUE) or the p parameter (FALSE)\n');
		fprintf(fid,'%d\n', get(handles.radio_slope,'Val'));
		comm = '# The ''Subset'' var (example: [0 19] -> (82 90)); ([9 9] -> (91 00)); ([19 0] -> (01 09))';
		helper_writeFile(handles, fid, comm, 'subset')
		comm = '# name of a netCDF file with quality flags. If not provided no quality check is done.';
		helper_writeFile(handles, fid, comm, 'flags')
		fprintf(fid,'# The ''spline'' var (logical). To compute spline interpolate the missing values\n');
		fprintf(fid,'%d\n', get(handles.check_integDim,'Val'));
		fprintf(fid,'# The ''scale'' var. Scale the final rate by this value\n%s\n',get(handles.edit_scale,'Str'));
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)
		fprintf(fid,'# The ''do3x3'' var (logical). To do calculations over 3x3 windows.\n');
		fprintf(fid,'%d\n', get(handles.check_average3x3, 'Val'));
		comm = '# Name of an optional Land mask grid file';
		helper_writeFile(handles, fid, comm, 'fname', 2)

	elseif (val == 3)	% Apply flags -- CASE 3
		comm = '# name of a netCDF file with quality flags. If not provided no quality check is done.';
		helper_writeFile(handles, fid, comm, 'flags')
		fprintf(fid,'# The ''Ncells'' variable\n%s\n',get(handles.edit_nCells,'Str'));
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)

	elseif (val == 4)	% Seazonal Averages -- CASE 4
		comm = '# The ''Periods'' variable (normaly, the number of days of each period)';
		helper_writeFile(handles, fid, comm, 'periods')
		comm = '# name of a netCDF file with quality flags. If not provided no quality check is done.';
		helper_writeFile(handles, fid, comm, 'flags')
		fprintf(fid,'# The ''Ncells'' variable\n%s\n',get(handles.edit_nCells,'Str'));
		fprintf(fid,'# Not used here but need to set as empty\n[]\n');
		fprintf(fid,'# The ''spline'' var (logical or 1x2 vector). To compute spline interpolate the missing values\n');
		if (get(handles.check_integDim,'Val'))
			if (strcmp(get(handles.edit_subSetA,'Str'),'0') && strcmp(get(handles.edit_subSetB,'Str'),'0'))
				fprintf(fid,'1\n');
			else
				fprintf(fid,sprintf('[%s %s]\n',get(handles.edit_subSetA,'Str'), get(handles.edit_subSetB,'Str')));
			end
		else
			fprintf(fid,'0\n');
		end
		fprintf(fid,'# The ''What'' variable that controls what statistic to compute\n');
		valW = get(handles.popup_what,'Val');		str = get(handles.popup_what,'Str');
		switch (str{valW})
			case 'MEAN',	fprintf(fid,'0\n');
			case 'MEDIAN',	fprintf(fid,'1\n');
			case 'MIN',		fprintf(fid,'2\n');
			case 'MAX',		fprintf(fid,'3\n');
			case 'STD',		fprintf(fid,'4\n');
			otherwise,		fprintf(fid,'0\n');
		end
		comm = '# Name of a file with Lon,Lat locations where to output the entire time series';
		helper_writeFile(handles, fid, comm, 'fname', 2)
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)
		comm = '# Name of an optional Land mask grid file';
		helper_writeFile(handles, fid, comm, 'fname', 2)

	elseif (val == 5)	% Polygonal Averages -- CASE 5
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)
		fprintf(fid,'# The ''What'' variable that controls what statistic to compute\n');
		valW = get(handles.popup_what,'Val');		str = get(handles.popup_what,'Str');
		switch (str{valW})
			case 'MEAN',	fprintf(fid,'[]\n');
			case 'MEDIAN',	fprintf(fid,'char median\n');
			case 'MIN',		fprintf(fid,'char min\n');
			case 'MAX',		fprintf(fid,'char max\n');
			case 'STD',		fprintf(fid,'char std\n');
			otherwise,		fprintf(fid,'[]\n');
		end
		comm = '# The ''fnamepolys'' var - File name of a x,y polygon or a list of polygons. Empty = fish polys from fig';
		helper_writeFile(handles, fid, comm, 'fname', 2)
		comm = '# The ''Subset'' var (example: [0 19] -> (82 90)); ([9 9] -> (91 00)); ([19 0] -> (01 09))';
		helper_writeFile(handles, fid, comm, 'subset')
		comm = '# name of a netCDF file with quality flags. If not provided no quality check is done.';
		helper_writeFile(handles, fid, comm, 'flags')

	elseif (val == 6)	% L2 MODIS Averages -- CASE 10
		comm = '# The ''Periods'' variable (normaly, the number of days of each period)';
		helper_writeFile(handles, fid, comm, 'periods')
		valW = get(handles.popup_what,'Val');		str = get(handles.popup_what,'Str');
		switch (str{valW})
			case 'MEAN',	fprintf(fid,'# Which statistics\n0\n');
			case 'MEDIAN',	fprintf(fid,'# Which statistics\n1\n');
			case 'MIN',		fprintf(fid,'# Which statistics\n2\n');
			case 'MAX',		fprintf(fid,'# Which statistics\n3\n');
			case 'STD',		fprintf(fid,'# Which statistics\n4\n');
			otherwise,		fprintf(fid,'# Which statistics\n0\n');
		end
		comm = '# A 2 elements vector with the MIN and MAX values allowed on the Z function (default [0 inf])';
		helper_writeFile(handles, fid, comm, 'subset')
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)
		comm = '# Name of an optional Land mask grid file';
		helper_writeFile(handles, fid, comm, 'fname', 2)

	elseif (val == 7)	% Climatologies -- CASE 4
		comm = '# The ''Periods'' variable (The month(s) that we want the climatology)';
		helper_writeFile(handles, fid, comm, 'periods')
		comm = '# name of a netCDF file with quality flags. If not provided no quality check is done.';
		helper_writeFile(handles, fid, comm, 'flags')
		fprintf(fid,'# The ''Ncells'' variable\n%s\n',get(handles.edit_nCells,'Str'));
		fprintf(fid,'# Not used here but need to set as empty\n[]\n');
		fprintf(fid,'# The ''spline'' var (a character with the keyword CLIMA)\nchar CLIMA\n');
		fprintf(fid,'# The ''What'' variable that controls what statistic to compute\n');
		valW = get(handles.popup_what,'Val');		str = get(handles.popup_what,'Str');
		switch (str{valW})
			case 'MEAN',	fprintf(fid,'0\n');
			case 'MEDIAN',	fprintf(fid,'1\n');
			case 'MIN',		fprintf(fid,'2\n');
			case 'MAX',		fprintf(fid,'3\n');
			case 'STD',		fprintf(fid,'4\n');
			otherwise,		fprintf(fid,'0\n');
		end
		fprintf(fid,'# Not used here but need to set as empty\n[]\n');
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)
		comm = '# Name of an optional Land mask grid file';
		helper_writeFile(handles, fid, comm, 'fname', 2)

	elseif (val == 8)	% Per cell correlation coefficient -- CASE 11
		comm = '# Name of the other netCDF file whose correlation with loaded array will be estimated';
		helper_writeFile(handles, fid, comm, 'fname', 2)
		comm = '# The ''Subset'' var (example: [0 19] -> (82 90)); ([9 9] -> (91 00)); ([19 0] -> (01 09))';
		helper_writeFile(handles, fid, comm, 'subset')
		fprintf(fid,'# The ''p-value''.\n');
		fprintf(fid,'%d\n', get(handles.radio_pValue,'Val'));
		comm = '# Name of the netCDF file where to store the result. If not provided, it will be asked here';
		helper_writeFile(handles, fid, comm, 'fname', 1)

	elseif (val == 9)	% count the acumulated number of non-NaNs -- CASE 8
		fprintf(fid,'# OPT.\n');
		fprintf(fid,'char count\n');
		comm = '# The ''Subset'' var';
		helper_writeFile(handles, fid, comm, 'subset')

	elseif (val == 10)	% count bloom events -- CASE 12
		fprintf(fid,'# Threshold number to be considered a bloom.\n3\n');
		comm = '# The ''Subset'' var';
		helper_writeFile(handles, fid, comm, 'subset')
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
		if (isempty(fname)),	fprintf(fid,'[]\n');
		else,					fprintf(fid,'char %s\n', fname);
		end
	elseif (strncmp(opt1, 'subset',3))
		fprintf(fid, sprintf('[%s %s]\n', get(handles.edit_subSetA,'Str'), get(handles.edit_subSetB,'Str')));
	elseif (strncmp(opt1, 'flags',3))
		fname = get(handles.edit_fname3, 'Str');
		if (isempty(fname)),	fprintf(fid,'[]\n');
		else,					fprintf(fid,'char %s\n', fname);
		end
		fprintf(fid,sprintf('# The ''quality'' value. Ignored if fname = []\n%s\n', get(handles.edit_qualFlag,'Str')));
	elseif (strncmp(opt1, 'periods',3))
		fprintf(fid,sprintf('%s\n',get(handles.edit_periods,'Str')));
	end

% -------------------------------------------------------------------------------------
function fname = helper_getFile(handles, fname, n)
% Get and check existence of a input filename
	checa = true;
	method = 'get';
	if (n == 1)		% This box is used both for read/write. Only first case needs to be checked for existence
		str = get(handles.text_fname1,'Str');
		if (str(1) == 'O')		% Output file. Do not check if file exists
			checa = false;
			method = 'put';
		end
	end

	if (isempty(fname))		% Direct call
		[FileName, PathName] = put_or_get_file(handles, ...
			{'*.nc;*.NC;*.dat', 'Data files (*.nc,*.NC,*.dat)';'*.*', 'All Files (*.*)'},'file',method);
		if isequal(FileName,0),		return,		end

	else					% File name on input
		[PathName,FNAME,EXT] = fileparts(fname);
		if (~isempty(PathName)),	PathName = [PathName filesep];	end
		FileName = [FNAME EXT];
	end
	pause(0.01);	fname = [PathName FileName];

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
'Callback',@slices_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Enter the name of 3D netCDF file',...
'Tag','edit_inputName');

uicontrol('Parent',h1, 'Position',[350 396 21 21],...
'Callback',@slices_uiCB,...
'Tooltip','Browse for a 3D netCDF file',...
'Tag','push_inputName');

uicontrol('Parent',h1, 'Position',[10 370 160 15],...
'String','Scale color to global min/max',...
'Style','checkbox',...
'Tooltip','If checked, color palette is computed using global min/max',...
'Tag','check_globalMinMax');

uicontrol('Parent',h1, 'Position',[180 368 120 19],...
'Callback',@slices_uiCB,...
'String','Recompute Min/Max',...
'Tooltip','Force recomputing of global min/max',...
'Visible','off',...
'Tag','push_globalMinMax');

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
'Callback',@slices_uiCB,...
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
'Callback',@slices_uiCB,...
'String','1',...
'Style','edit',...
'Tooltip','Slice number (to go to a specific one, enter a valid slice number here)',...
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
'Callback',@slices_uiCB,...
'FontSize',10,...
'ForegroundColor',[0 0 1],...
'String','Help',...
'Tag','push_help');

uicontrol('Parent',h1, 'Position',[10 280 361 20],...
'BackgroundColor',[1 1 1],...
'Callback',@slices_uiCB,...
'String',{' '; 'Zonal Means'; 'Per cell rate of change '; 'Apply quality flags'; 'Seasonal Averages'; 'Polygonal Averages'; ...
	'L2 MODIS Averages'; 'Climatologies (months)'; 'Per cell correlation coefficient'; 'Count non-NaNs cells'; ...
	'Find Chlor_a blooms'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_cases');

uicontrol('Parent',h1, 'Position',[10 248 55 21],...
'Callback',@slices_uiCB,...
'String','Slope',...
'Value',1, ...
'Style','radiobutton',...
'Tag','radio_slope');

uicontrol('Parent',h1, 'Position',[80 248 55 21],...
'Callback',@slices_uiCB,...
'String','p value',...
'Style','radiobutton',...
'Tag','radio_pValue');

uicontrol('Parent',h1, 'Position',[151 244 139 21],...
'String','Integration in longitude',...
'Style','checkbox',...
'Value',0,...
'Tag','check_integDim');

uicontrol('Parent',h1, 'Position',[190 250 150 21],...
'String','Average over 3x3 windows',...
'Style','checkbox',...
'Value',0,...
'Tooltip', 'Do the calculations over 3x3 windows. That is, compute averages values (but slower).',...
'Visible', 'off',...
'Tag','check_average3x3');

uicontrol('Parent',h1, 'Position',[262 232 81 21],...
'String','Do trends',...
'Style','checkbox',...
'Tag','check_doTrends');

uicontrol('Parent',h1, 'Position',[70 244 31 21],...
'BackgroundColor',[1 1 1],...
'Callback',@slices_uiCB,...
'String','1',...
'Style','edit',...
'Tag','edit_boxSize');

uicontrol('Parent',h1, 'Position',[13 245 55 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Box size',...
'Style','text',...
'Tag','text_boxSize');

uicontrol('Parent',h1, 'Position',[57 199 41 21],...
'BackgroundColor',[1 1 1],...
'Callback',@slices_uiCB,...
'String','0',...
'Style','edit',...
'Tag','edit_subSetA');

uicontrol('Parent',h1, 'Position',[97 199 41 21],...
'BackgroundColor',[1 1 1],...
'Callback',@slices_uiCB,...
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

uicontrol('Parent',h1, 'Position',[64 230 171 23],...
'BackgroundColor',[1 1 1],...
'Callback',@slices_uiCB,...
'Style','edit',...
'Tag','edit_periods');

uicontrol('Parent',h1, 'Position',[300 196 71 23],...
'BackgroundColor',[1 1 1],...
'Callback',@slices_uiCB,...
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
'Callback',@slices_uiCB,...
'String','0',...
'Style','edit',...
'Tooltip','Data gaps groups smaller than this are filled by interpolation',...
'Tag','edit_nCells');

uicontrol('Parent',h1, 'Position',[298 199 38 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','Scale',...
'Style','text',...
'Tag','text_scale');

uicontrol('Parent',h1, 'Position',[336 197 31 21],...
'BackgroundColor',[1 1 1],...
'Callback',@slices_uiCB,...
'String','1',...
'Style','edit',...
'Tag','edit_scale');

uicontrol('Parent',h1, 'Position',[10 150 341 21],...
'BackgroundColor',[1 1 1],...
'Callback',@slices_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_fname1');

uicontrol('Parent',h1, 'Position',[351 150 21 21],...
'Callback',@slices_uiCB,...
'Tag','push_fname1');

uicontrol('Parent',h1, 'Position',[10 100 341 21],...
'BackgroundColor',[1 1 1],...
'Callback',@slices_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_fname2');

uicontrol('Parent',h1, 'Position',[351 100 21 21],...
'Callback',@slices_uiCB,...
'Tag','push_fname2');

uicontrol('Parent',h1, 'Position',[319 71 31 21],...
'BackgroundColor',[1 1 1],...
'Callback',@slices_uiCB,...
'Style','edit',...
'String','6',...
'Tag','edit_qualFlag');

uicontrol('Parent',h1, 'Position',[10 50 341 21],...
'BackgroundColor',[1 1 1],...
'Callback',@slices_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_fname3');

uicontrol('Parent',h1, 'Position',[351 50 21 21],...
'Callback',@slices_uiCB,...
'Tag','push_fname3');

uicontrol('Parent',h1, 'Position',[9 9 141 21],...
'Callback',@slices_uiCB,...
'FontSize',9,...
'String','Run Plugin fun (opt)',...
'Tag','push_runPlugin');

uicontrol('Parent',h1, 'Position',[276 8 101 22],...
'Callback',@slices_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','Compute',...
'Tag','push_compute');

function slices_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
function varargout = aqua_suppfuns(opt, varargin)
% Supplement functions to allow using Aquamoto with plain netCDF coards grids

%	Copyright (c) 2004-2016 by J. Luis
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

% $Id: slices.m 10389 2018-04-27 16:13:46Z j $

	switch opt
		case 'coards_hdr',		[varargout{1:nargout}] = init_header_params(varargin{:});
		case 'coards_slice',	coards_sliceShow(varargin{:})
		case 'forGDAL_hdr',		[varargout{1:nargout}] = init_header_gdal(varargin{:});
		case 'forGDAL_slice',	gdal_sliceShow(varargin{:})
	end

% --------------------------------------------------------------------------
function out = init_header_params(handles,X,Y,head,misc,getAllMinMax)
% Use the OUT option when using this function to get several usefull info 
% about the netCDF file but NOT using a GUI.
%
% The 'getAllMinMax' when set to TRUE (or not provided) will scan the entire file to
% compute all individual layers 'actual_range'. However, this may slow quite a bit the
% loading of big files and apparently is not very used/useful because most of the time
% is spent loading the layer and the min/max computation is extremely fast.

	if (nargin < 6),		getAllMinMax = true;		end
	handles.x = X;		handles.y = Y;
	handles.time = [];
	handles.number_of_timesteps = misc.z_dim(1);		% ... NEEDS THINKING
	
	handles.x_min = head(1);			handles.x_max = head(2);
	handles.y_min = head(3);			handles.y_max = head(4);

	% ------------- Finish slider configurations -------------
	s = handles.nc_info;
	if (handles.number_of_timesteps > 1)
		st = [1 10] / (handles.number_of_timesteps - 1);
		id = find(strcmpi('time',{s.Dataset.Name}));			% ONLY WHEN 3RTH DIM IS CALLED time
		if (isempty(id))
			id = find(strcmpi('depth',{s.Dataset.Name}));		% ... or 'depth'
		end
		if (~isempty(id))
			handles.time = double(nc_funs('varget', handles.fname, s.Dataset(id).Name));
			if (numel(handles.time) ~= handles.number_of_timesteps)		% We have a 4D singleton array. Get 'time' from second dim
				id = strcmpi(s.Dimension(end-1).Name, {s.Dataset.Name});
				handles.time = double(nc_funs('varget', handles.fname, s.Dataset(id).Name));
			end
		else
			handles.time = 1:handles.number_of_timesteps;
		end
		slMax = handles.number_of_timesteps;
	else
		slMax = 1+eps;	st = [1 1];		handles.time = 1;		% Defaults for no crashing
	end

	% ------ Compute individual and global min/maxs ----------------------------------
	handles.zMinMaxs = zeros(handles.number_of_timesteps,2);
	if (handles.IamTSU)			% This one needs special doing since we want min/max over ocean only
		aguentabar(0,'title','Computing water surface min/max')
		zBat = nc_funs('varget', handles.fname, 'bathymetry');
		ind = (zBat < 0);	clear zBat
		for (k = 1:handles.number_of_timesteps)
			Z = nc_funs('varget', handles.fname, s.Dataset(misc.z_id).Name, [(k-1) 0 0], [1 s.Dataset(misc.z_id).Size(end-1:end)]);
			Z = Z(ind);
			indNaN = isnan(Z);
			if (any(indNaN)),	Z(indNaN) = [];	end
			if (~isempty(Z))
				handles.zMinMaxs(k,:) = [min(Z) max(Z)];
			end
			aguentabar(k/handles.number_of_timesteps);
		end
		handles.zMinMaxsGlobal = [min(handles.zMinMaxs(:,1)) max(handles.zMinMaxs(:,2))];
		k = 1;
		while ((k <= handles.number_of_timesteps) && isequal(handles.zMinMaxs(k,:), [0 0]))
			k = k + 1;
		end
		if (k == handles.number_of_timesteps + 1)
			k = 1;		% To not error below
			warndlg('WARNING: All layers are zero. You better go home.', 'Warning')
		end
		head(5:6) = handles.zMinMaxs(k,:);			% Take the first non zero slice min/max
		set(handles.radio_shade, 'Val', 1);		set(handles.radio_noShade, 'Val', 0)	% Set shading by default
	elseif (getAllMinMax)
		aguentabar(0,'title','Computing global min/max')
		for (k = 1:handles.number_of_timesteps)
			Z = nc_funs('varget', handles.fname, s.Dataset(misc.z_id).Name, [(k-1) 0 0], [1 s.Dataset(misc.z_id).Size(end-1:end)]);
			if (isa(Z, 'single'))			% min/max are bugged when NaNs in singles
				zz = grdutils(Z,'-L');
				handles.zMinMaxs(k,:) = [zz(1) zz(2)];
			else
				handles.zMinMaxs(k,:) = [double(min(Z(:))) double(max(Z(:)))];
			end
			aguentabar(k/handles.number_of_timesteps);
		end
		handles.zMinMaxsGlobal = [min(handles.zMinMaxs(:,1)) max(handles.zMinMaxs(:,2))];
		while (isequal(handles.zMinMaxs(k,:), [0 0]) && (k <= handles.number_of_timesteps))
			k = k + 1;
		end
		if (k > handles.number_of_timesteps)
			warndlg('WARNING: All layers are zero. You better go home.', 'Warning')
		end
		head(5:6) = handles.zMinMaxs(k,:);			% Take the first non zero slice min/max
	else
		ind = strcmp({s.Dataset(misc.z_id).Attribute.Name},'actual_range');
		if (any(ind))
			handles.zMinMaxsGlobal = s.Dataset(misc.z_id).Attribute(ind).Value;
		else
			handles.zMinMaxsGlobal = [0 1];
		end
		head(5:6) = handles.zMinMaxsGlobal;
	end
	handles.minWater = handles.zMinMaxsGlobal(1);
	handles.maxWater = handles.zMinMaxsGlobal(2);
	handles.geog = aux_funs('guessGeog',head(1:4));
	% ---------------------------------------------------------------------------------

	if (~handles.IamTSU)
		handles.cmapLand = jet(256);		% Reset the default colormap (default's Aquamoto is a specific one)
	end
	handles.head = head;
	handles.illumComm = [];					% New file. Reset illum state.
	handles.imgBat = [];
	handles.netcdf_z_id = misc.z_id;
	handles.is_coards = true;

	% -------------------- See if we have a projection ----------------------------------
	if (~isempty(misc.strPROJ4)),	handles.strPROJ4 = misc.strPROJ4;
	else,							handles.strPROJ4 = [];
	end
	if (~isempty(misc.srsWKT)),		handles.srsWKT = misc.srsWKT;
	else,							handles.srsWKT = [];
	end

	if (nargout)
		out = handles;
	else
		set(handles.edit_Ncols,'String',sprintf('%d',misc.z_dim(end)))
		set(handles.edit_Nrows,'String',sprintf('%d',misc.z_dim(end-1)))
		set(handles.slider_layer,'Min',1,'Max',slMax,'Val',1,'SliderStep',st) 	
		set(handles.edit_globalWaterMin,'String',handles.zMinMaxsGlobal(1))
		set(handles.edit_globalWaterMax,'String',handles.zMinMaxsGlobal(2))
		set(handles.hTabAnuga,'String','netCDF')
		set_common(handles, handles.head)
		guidata(handles.figure1,handles)
		% Store the nc_info and z_id in Mirone handles so we can access it from there as well
		if (~isempty(handles.hMirFig))	% This is another patch for the mess of the hMirFig existing or not
			hFig = handles.hMirFig;		% and how the code flow in places relies on that.
		else
			hFig = handles.hMirFig_safe;
		end
		handMir = guidata(hFig);
		handMir.netcdf_z_id = misc.z_id;
		handMir.nc_info = handles.nc_info;
		handMir.time_z  = handles.time;
		handMir.nLayers = misc.z_dim(1);
		guidata(hFig, handMir)
	end

% --------------------------------------------------------------------------
function out = init_header_gdal(handles)
% Read a multiband file with gdal and fill the header parameters

	handles.illumComm = [];					% New file. Reset illum state.
	handles.imgBat = [];
	handles.is_otherMultiband = true;
	att = gdalread(handles.fname,'-M','-C');
	X = linspace(att.GMT_hdr(1), att.GMT_hdr(2), att.RasterXSize);
	Y = linspace(att.GMT_hdr(3), att.GMT_hdr(4), att.RasterYSize);
	handles.number_of_timesteps = att.RasterCount;
	handles.time = 1:att.RasterCount;
	handles.head = att.GMT_hdr;
	handles.x_min = X(1);		handles.x_max = X(2);
	handles.y_min = Y(1);		handles.y_max = Y(2);
	handles.x = X;				handles.y = Y;
	handles.flip_on_read = true;
	if (isempty(att.GeoTransform)),		handles.flip_on_read = false;	end

	st = [1 10] / (att.RasterCount - 1);

	% ------ Compute individual and global min/maxs ------------------------
	handles.zMinMaxs = zeros(att.RasterCount, 2);
	for (k = 1:att.RasterCount)
		if ( isnan(att.Band(1).MinMax(1)) )		% Shit, nothing usable here
			if (~nargout),	set(handles.check_globalMinMax,'Enable', 'off'),	end
			break
		end
		handles.zMinMaxs(k,:) = att.Band(k).MinMax(:);
	end
	handles.zMinMaxsGlobal = [min(handles.zMinMaxs(:,1)) max(handles.zMinMaxs(:,2))];
	handles.minWater = handles.zMinMaxsGlobal(1);
	handles.maxWater = handles.zMinMaxsGlobal(2);
	handles.geog = aux_funs('guessGeog',att.GMT_hdr(1:4));
	% ---------------------------------------------------------------------------------

	handles.cmapLand = jet(256);			% Reset the default colormap (default's Aquamoto is a specific one)

	% -------------------- See if we have a projection ----------------------------------
	if (~isempty(att.ProjectionRef)),	handles.srsWKT = att.ProjectionRef;
	else,								handles.srsWKT = [];
	end
	handles.strPROJ4 = [];

	if (nargout)
		out = handles;
	else
		set(handles.slider_layer,'Min',1,'Max',att.RasterCount,'Val',1,'SliderStep',st)
		set( handles.edit_Ncols,'String',sprintf('%d',att.RasterXSize) )
		set( handles.edit_Nrows,'String',sprintf('%d',att.RasterYSize) )
		set( handles.check_splitDryWet,'Enable', 'off' )
		set( handles.push_runIn,'Enable', 'off' )
		set( handles.slider_transparency,'Enable', 'off' )
		set(handles.edit_globalWaterMin,'String',handles.zMinMaxsGlobal(1))
		set(handles.edit_globalWaterMax,'String',handles.zMinMaxsGlobal(2))
		set(handles.hTabAnuga,'String','GDALish')
		set_common(handles, handles.head)
		guidata(handles.figure1,handles)
	end

% --------------------------------------------------------------------------------------------
function gdal_sliceShow(handles, att)
% Read the slice with gdalread and send in the array to coards_sliceShow() to do the rest

	opt_U = ' ';
	if (handles.flip_on_read),		opt_U = '-U';	end
	opt_B = sprintf('-B%d', handles.sliceNumber + 1);
	Z = gdalread(handles.fname, opt_B, opt_U, '-C');
	if (isa(Z, 'uint16') || isa(Z, 'int16') || isa(Z, 'double'))
		Z = single(Z);
	end
	coards_sliceShow(handles, Z)

% --------------------------------------------------------------------------------------------
function coards_sliceShow(handles, Z)
% ...

	cmap_slice = [];		% Will have a cmap when doing individual frames of a global cmap

	if (nargin == 1)		% Otherwise we suposedly already know Z (from gdalread)
		if (isempty(handles.fname))
			errordlg('Hey Lou. What about a walk on the Wild Side? Maybe you''ll find a little file there that you can use here!','Chico clever')
			return
		end

		z_id = handles.netcdf_z_id;
		s = handles.nc_info;					% Retrieve the .nc info struct 
		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [handles.sliceNumber 0 0], [1 s.Dataset(z_id).Size(end-1:end)]);
		dims = s.Dataset(z_id).Dimension;		% Variable names of dimensions z variable
		y_id = find(strcmp({s.Dataset.Name}, dims{end-1}));
		if (~isempty(y_id))						% Check if we need to flipud(Z)
			Y = double(nc_funs('varget', handles.fname, s.Dataset(y_id).Name));
			if (Y(2) < Y(1)),	Z = flipud(Z);	end
		end
		try		handles.head(5:6) = handles.zMinMaxs(handles.sliceNumber+1,:);	end		%  handles.head was from 1st layer
	end

	have_nans = 0;
	if (isa(Z,'single') && ~handles.IamTSU)		% Actually I don't know if this is used in any circumnstances.
		have_nans = grdutils(Z,'-N');			% No worry, very fast
	end
	if (have_nans && handles.useLandPhoto)
		alphaMask = alloc_mex(size(Z),'uint8');	% Create an image mask of Dry/Wets
		alphaMask(~isnan(Z)) = 255;				% nan pixeis will be transparent
	end

	try			% SHOULD BE CONVERTED INTO A PROPER TEST (AQUAMOTO OR NOT-AQUAMOTO)
		splitDryWet = get(handles.check_splitDryWet,'Val');		% See if we need to build wet and dry images, or only one
	catch
		splitDryWet = false;
		handles.cmapLand = jet(256);			% Was not yet deffined if a TSU file opened in SLICES
	end

	if (splitDryWet && handles.IamTSU)
		zBat = nc_funs('varget', handles.fname, 'bathymetry');
		dife = cvlib_mex('absDiff', zBat, Z);
		indLand = (dife < 1e-2);				% The 1e-2 SHOULD be parameterized
		zero = 0;
		if     (zero > handles.head(6)),	zero = handles.head(6) - eps;
		elseif (zero < handles.head(5)),	zero = handles.head(5) + eps;
		end
		Z(indLand) = zero;						% Smash these too to not steal the color dynamics
		% Because of the inundation the handles.zMinMaxs may not be correct, so must compute them again
		zz = grdutils(Z,'-L');		handles.head(5:6) = zz(:)';		% This may contradict the 'zero' above but 'bad luck'
		if (handles.compute_runIn && isempty(handles.indMaxWater))
			handles.indMaxWater = compute_indMaxWater(handles, zBat);
		end

		if (handles.useLandPhoto)
			alfa = 255;		% Means land will be completely opac and water 100% transparent
			[nr, nc] = size(indLand);
			if (~isempty(handles.hMirFig) && ishandle(handles.handMir.hImg))	% A Mirone figure with image already exists
				alphaMask = get(handles.handMir.hImg,'AlphaData');
				if (numel(alphaMask) == 1)		% No AlphaMask yet, but we'll need one further down
					alphaMask(nr, nc) = uint8(0);
					alphaMask(~indLand) = alfa;
				else							% AlphaMask exists, but we need to update it to reflect this slice water level
					alphaMask(indLand) = 0;		alphaMask(~indLand) = alfa;
				end
			else								% The Mirone figure will be created later. Compute the AlphaMask
				alphaMask(nr, nc) = uint8(0);
				alphaMask(~indLand) = alfa;
			end
		end

		% For global color scales we need to compute this particular frame cmap
		if (size(handles.cmapWater,1) > 2 && get(handles.check_globalMinMax, 'Val'))
			percent_local_extrema = (handles.head(5:6) - handles.zMinMaxsGlobal(1)) / diff(handles.zMinMaxsGlobal);
			ind_colors = round(percent_local_extrema * size(handles.cmapWater,1));
			cmap_slice = handles.cmapWater(ind_colors(1):ind_colors(2),:);
			cmap_slice = interp1(cmap_slice, linspace(1,size(cmap_slice,1), size(handles.cmapWater,1)));
		end
	end

	% ----- Open or update a Mirone window with the slice display ----
	if (isempty(handles.hMirFig) || ~ishandle(handles.hMirFig))			% First run or killed Mirone window
		tmp.X = handles.x;		tmp.Y = handles.y;		tmp.head = handles.head;	tmp.cmap = handles.cmapLand;
		if (~isempty(cmap_slice)),	tmp.cmap = cmap_slice;				% Restoring a particular frame of a global cmap
		elseif (handles.IamTSU),	tmp.cmap = handles.cmapWater;		% A possible first fig
		end
		tmp.name = sprintf('Layer = %g',handles.time(handles.sliceNumber+1));
		if (~isempty(handles.srsWKT)),		tmp.srsWKT = handles.srsWKT;	end
		hFigs = findobj(0,'type','figure');
		nLayers = 0;
		if (numel(hFigs) == 2)	% Often we have an empty Mir fig but that is very difficult to use here. So blow it
			inds = [isempty(getappdata(hFigs(1), 'IAmAMirone')) isempty(getappdata(hFigs(2), 'IAmAMirone'))];
			hFigs = hFigs(~inds);							% Only one of them is a Mirone fig
			if (~isempty(hFigs))							% It happened once in debug but might happen in a real future case
				handThis = guidata(hFigs);
				if (~isempty(handThis))
					nLayers  = handThis.nLayers;
					if (~handThis.validGrid),		delete(hFigs),	clear handThis,	end
				end
			end
		end

		reset = false;
		if (~ishandle(handles.hMirFig))						% Save this before handles.hMirFig is reset
			nLayers = handles.handMir.nLayers;
			reset = true;
		end

		handles.hMirFig = mirone(Z, tmp);
		move2side(handles.figure1,handles.hMirFig,'left')
		handles.handMir = guidata(handles.hMirFig);			% Get the handles of the now existing Mirone fig

		if (reset || nLayers > 0)							% Figure was killed but the folowing info exists and is needed
			handles.handMir.nLayers = nLayers;
			handles.handMir.netcdf_z_id = handles.netcdf_z_id;
			handles.handMir.nc_info = handles.nc_info;
			handles.handMir.time_z  = handles.time;
			handles.handMir.grdname = handles.fname;		% At least grid_profiler() needs this for the 3D interpolations
		end

		handles.firstLandPhoto = true;
		if (handles.useLandPhoto)
			h = image('XData',handles.geoPhotoX,'YData',handles.geoPhotoY, 'CData',handles.geoPhoto, ...
			          'Parent',handles.handMir.axes1, 'Tag', 'LandPhoto');
			uistack_j(h,'bottom')
			handles.firstLandPhoto = false;
			set(handles.handMir.hImg,'AlphaData',alphaMask)	% 'alphaMask' was updated ... maybe somewhere
		end
		handles.imgBat = [];				% Make sure this one is always reset when user kills fig

	else									% We already have a Mirone image. Update it with this new slice
		handles.handMir = guidata(handles.hMirFig);			% Get updated handles to see if illum has changed
		if (~(isa(Z,'uint8') || isa(Z,'int8')))
			setappdata(handles.handMir.figure1,'dem_z',Z);	% Update grid so that coursor display correct values
		end													% Have to do it here because minmax arg to scalet8 CHANGES Z
		if (~isempty(cmap_slice))							% A particular frame of a global cmap
			set(handles.handMir.figure1,'Colormap',cmap_slice)
		end

		if (handles.useLandPhoto)							% External Land image
			if (handles.firstLandPhoto)						% First time, create the background image
				h = image('XData',handles.geoPhotoX,'YData',handles.geoPhotoY, 'CData',handles.geoPhoto, ...
				          'Parent',handles.handMir.axes1, 'Tag', 'LandPhoto');
				uistack_j(h,'bottom')
				handles.firstLandPhoto = false;
			end
			set(handles.handMir.hImg,'AlphaData',alphaMask)	% 'alphaMask' was updated ... somewhere
		end

		if (~get(handles.check_globalMinMax, 'Val')),		minmax = [];		% Use Slice's min/max
		else,							minmax = handles.zMinMaxsGlobal;
		end
		if (isa(Z,'int8') && ~isempty(minmax)),	minmax = [];	end	% We don't want to scale a 1 byte array

		if (isempty(minmax) && ~isempty(handles.handMir.img_with_minmax))	% Happens when using a Thematic palette
			minmax = handles.handMir.img_with_minmax;
		end

		if (~isa(Z,'uint8'))
			if (~isempty(minmax)),		img = scaleto8(Z, 8, minmax);
			else,						img = scaleto8(Z);
			end
		else
			img = Z;
		end

		if (handles.IamTSU && splitDryWet)
			recomp = false;
			if (~isempty(handles.landIllumComm_bak) && ~isequal(handles.landIllumComm_bak, handles.landIllumComm))
				recomp = true;
			end
			if (recomp || isempty(handles.imgBat))		% First time, compute it (not necessarily shaded)
				% Put the cmap discontinuity at the zero of bat
				zz = grdutils(zBat,'-L');	% I'm lazy to fish this inside the s info struct
				head = handles.head;	head(5:6) = [zz(1) zz(2)];
				handles.cmapBat = aquamoto('makeCmapBat', handles, head, handles.cmapLand, 1);
				handles.imgBat  = ind2rgb8(scaleto8(zBat), handles.cmapBat);
				if (get(handles.radio_shade, 'Val'))
					R = aquamoto('illumByType', handles, zBat, head, handles.landIllumComm);
					handles.imgBat = shading_mat(handles.imgBat, R, 'no_scale');
				end
			end
			img = aquamoto('do_imgWater', handles, 1, Z, handles.imgBat, indLand);

		elseif (get(handles.radio_shade, 'Val'))		% handles.radio_shade was set to [] in slices to not error here
			indVar = 1;									% FAR FROM SURE THAT THIS IS CORRECT
			img = ind2rgb8(img, handles.cmapLand);		% img is now RGB
			head = handles.head;
			if (~isempty(handles.ranges{indVar})),		head(5:6) = handles.ranges{indVar};		end
			if (head(5) == 0 && head(6) == 0)
				warndlg('Something screwed up. Don''t know this grid min/max and must comput it now.','Warning')
				zz = grdutils(Z, '-L');		head(5:6) = [zz(1) zz(2)];
			end
			R = aquamoto('illumByType', handles, Z, head, handles.landIllumComm);
			img = shading_mat(img,R,'no_scale');		% and now it is illuminated
		end

		set(handles.handMir.hImg, 'CData', img)
		set(handles.handMir.figure1, 'Name', sprintf('Level = %.10g',handles.time(handles.sliceNumber+1)))
		setappdata(handles.handMir.figure1,'dem_x',handles.x);		% Don't get bad surprises (like load another file)
		setappdata(handles.handMir.figure1,'dem_y',handles.y);

		% See if we have a colorbar to update
		if (strcmp(get(handles.handMir.PalAt,'Check'),'on'))
			show_palette(handles.handMir, 'At', 'update')
		elseif (strcmp(get(handles.handMir.PalIn,'Check'),'on'))
			show_palette(handles.handMir, 'In', 'update')
		end

		% See if image is illuminated
		if (handles.handMir.Illumin_type >= 1 && handles.handMir.Illumin_type <= 4)
			illumComm = getappdata(handles.handMir.figure1,'illumComm');
			if (handles.handMir.Illumin_type == 1)
				if (handles.geog),	R = grdgradient_m(Z,handles.head,'-M',illumComm,'-Nt','-a1');
				else,				R = grdgradient_m(Z,handles.head,illumComm,'-Nt','-a1');
				end
			else
				R = grdgradient_m(Z,handles.head,illumComm, '-a1');
			end
			img = ind2rgb8(img,get(handles.handMir.figure1,'Colormap'));
			mex_illuminate(img,R)		% Operates insitu too
			clear R
			set(handles.handMir.hImg,'CData',img),		refresh(handles.handMir.figure1)	% Crazzy beast does not always update the image !!!!!!!!!
		end
	end

    guidata(handles.figure1,handles)

	% Save also the updated header in Mirone handles
	handles.handMir.head = handles.head;
    guidata(handles.handMir.figure1,handles.handMir)

	% See if we have a displayed Colorbar and if yes update its values. The thing is done with quite some trickery.
	if (~isempty(handles.hMirFig) && ishandle(handles.hMirFig))
		try			% Let this try be here for a while before removing it no shit happens inside this
			if (splitDryWet && handles.IamTSU && size(handles.cmapWater,1) > 2)		% Cannot allow the 'blue' case because it's F bugged by TMW
				where = '';
				set(0,'ShowHiddenHandles','on')
				if (strcmp(get(handles.handMir.PalAt,'Check'),'on'))
					hAxPal  = findobj(handles.handMir.figure1,'Type','axes','Tag','MIR_CBat');
					handPal = handles.handMir.PalAt;
					where   = 'At';				% Color bar location relative to mirone axes1
					ud_old  = get(handPal, 'UserData');
				elseif (strcmp(get(handles.handMir.PalIn,'Check'),'on'))
					hAxPal  = findobj(handles.handMir.figure1,'Type','axes','Tag','MIR_CBin');
					handPal = handles.handMir.PalIn;
					where   = 'In';				% Color bar location relative to mirone axes1
				end
				set(0,'ShowHiddenHandles','off')
				if (~isempty(where))						% Otherwise what we might have is a Floating Fig with the cbar
					delete(hAxPal)							% Delete the old cbar
					set(handPal,'Check','off')				% So that a new color bar is created instead of delting the old one
					show_palette(handles.handMir, where)	% __-*++ ===> Create a new colorbar with correct values <=== ++*-__
					if (strcmp(where, 'At'))				% Only this needs the 'origFigWidth' info preserved
						ud_new = get(handPal, 'UserData');	% 'show_palette' changed it
						set(handPal, 'UserData', [ud_new(1) ud_old(2)]);		% ud_old(2) has the true 'origFigWidth'
					end
					set(handPal,'Check','on')
				end
			end
		catch
			disp(lasterr)
		end
	end

% ---------------------------------------------------------------------------
function indLand = compute_indMaxWater(handles, zBat)
% Compute the mask of max water
	aguentabar(0, 'title','Computing max water height ...', 'CreateCancelBtn');
	z_id = handles.netcdf_z_id;
	s = handles.nc_info;					% Retrieve the .nc info struct
	Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [0 0 0], [1 s.Dataset(z_id).Size(end-1:end)]);
	dife = cvlib_mex('absDiff', zBat, Z);
	indLand = (dife < 1e-2);				% The 1e-2 SHOULD be parameterized
	for (k = 1:numel(handles.time) - 1)
		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [k 0 0], [1 s.Dataset(z_id).Size(end-1:end)]);
		dife = cvlib_mex('absDiff', zBat, Z);
		ind = (dife < 1e-2);
		indLand = indLand & ind;
		h = aguentabar((k+1) / numel(handles.time));
		if (isnan(h)),	indLand = [];	break,	end
	end
	indLand = ~indLand;

% --------------------------------------------------------------------------
function set_common(handles, head)
% Common settingd to both 'init_header' functions
	set( handles.edit_x_min,'String',sprintf('%.8g',head(1)) )
	set( handles.edit_x_max,'String',sprintf('%.8g',head(2)) )
	set( handles.edit_y_min,'String',sprintf('%.8g',head(3)) )
	set( handles.edit_y_max,'String',sprintf('%.8g',head(4)) )
	set( handles.edit_x_inc,'String',sprintf('%.8g',head(8)) )
	set( handles.edit_y_inc,'String',sprintf('%.8g',head(9)) )

	set(handles.slider_layer,'Enable','on')
	set(handles.edit_sliceNumber,'Enable','on')
	set(handles.text_Info,'String',sprintf('Time steps = %d',handles.number_of_timesteps))

	set(handles.radio_multiLayer, 'Val', 1)
	set(handles.edit_multiLayerInc, 'Enable', 'on')
	set(handles.radio_timeGridsList,'Val',0)
	set([handles.textResize handles.popup_resize], 'Enable', 'off')
	set([handles.radio_stage handles.radio_xmoment handles.radio_ymoment handles.check_derivedVar], 'Enable', 'off')
	set([handles.edit_x_min handles.edit_x_max handles.edit_y_min handles.edit_y_max ...
		handles.edit_x_inc handles.edit_y_inc handles.edit_Ncols handles.edit_Nrows], 'Enable', 'inactive')
function aquaPlugin(handles, auto)
% Plugin function that is called by Aquamoto. Use this function to write custom code
% to solve particular problems taking advantage from the fact that a LOT of information
% about the netCDF files is stored in HANDLES.
%
% OPTIONS
%
% AUTO	- If logical and TRUE, search the OPTcontrol.txt file for the MIR_AQUAPLUG key that
%		  should point into a file name of a control script.
%		- If it is a string, than that is interpreted as the name of the control script.
%
%		A "control script" is a file with the EXACT arguments to select and run one of
%		main functions here as pointed by the CASOS cell array below.
%
%		One way of executing this functionality is to check the "Seek OPTcontrol.txt" checkbox
%		In which case the OPTcontrol.txt will be scanned for the name of the control script.
%		This works both for the ML and the standalone version.
%		The other way, restricted to the ML version, is to run in the Matlab command line:
%				aquamoto file.nc 'file_name_of_control_script'
%		OR
%				aquamoto('file.nc', 0)
%		In the later case the control script name is searched in the OPTcontrol.txt file

%	Copyright (c) 2004-2018 by J. Luis
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

% $Id: slices.m 10389 2018-04-27 16:13:46Z j $

	if (isempty(handles.fname))
		errordlg('Fast trigger, you probably killed my previous encarnation. Now you have to start again. Bye.','Error')
		return
	end
	internal_master = true;		% To know if flow control is determined by the contents of an external file (def NO).

	casos = {'zonal' ...			% 1 - Compute zonal means
			'tvar' ...				% 2 - Compute the Temp time rate of a file with annual means by fit of a straight line (Load entire file in memory)
			'applyFlags' ...		% 3 - Check against its 3D flags companion and replace values < FLAG to NaN
			'yearMeanFlag' ...		% 4 - Compute yearly averages from monthly data but checked against a quality flag file
			'polygAVG' ...			% 5 - Compute averages of whatever inside polygons (if any)
			'flagsStats' ...		% 6 - Compute per/pixel annual or month counts of pixel values with a quality >= flag
			'pass_by_count' ...		% 7 - Check the curently active 3D file against a count file
			'do_math' ...			% 8 - Perform some basic algebric operations with the 3D planes
			'conv2vtk' ...			% 9 - Convert a 3D netCDF file into a VTK format
			'L2_periods' ...		% 10 - Calculate composites of L2 products over fixed periods
			'corrcoef' ...			% 11 - Calculate correlation coefficients between 2 3D arrays
			'count_blooms' ...		% 12 - Count chlor_a blooms
			};

	qual = casos{11};			% <== Active by MANUAL selection. May be override by next section

	n_args = nargin;
	if (isfield(handles,'check_plugFun') && get(handles.check_plugFun, 'Val'))	% This way, the stand alone version can work too
		n_args = 2;		auto = true;
	end

	if (n_args == 2)				% Go figure out if we have a controlling script
		out = script_control(handles, auto);
		if (~isempty(out))
			qual = casos{out{1}};	% The rest will be applied blindly. If it screws, screws
			internal_master = false;
		else
			errordlg('You directed AQUAPLUGIN to work on script mode but OPTcontrol.txt misses necessary info','Error')
			return
		end
	end

	switch qual
		case 'zonal'				% CASE 1
			integ_lon = true;
			dlat = 1.0;
			trends = false;			% If true compute the trends (per stripe) of the zonal integration
			sub_set = [0 0];		% [jump_start stop_before_end], make it [] or [0 0] to be ignored
			fnamPoly1 = 'C:\a1\pathfinder\plataforma_poly.dat';		% If this name is uncorrect, another will be asked
			fnamPoly2 = 'poly_largo.dat';	%fnamPoly2= [];	% If it exists, compute difference of zonal integrations
			fnameFlag  = 'C:\a1\pathfinder\qual_82_09.nc';	% If not empty check againts this file (For monthly data)
			quality = 6;			% Retain only values of quality >= this (or <= abs(this) when MODIS). Ingored if fname = []
			if (internal_master)
				zonal(handles, dlat, integ_lon, trends, sub_set, fnamPoly1, fnamPoly2, fnameFlag, quality)
			else
				zonal(handles, out{2:end})
			end
		case 'tvar'					% CASE 2
			slope = true;			% TRUE to compute slope of linear fit, FALSE to compute p-value parameter
			sub_set = [0 0];		% [jump_start stop_before_end], make it [] or [0 0] to be ignored
			fname  = 'C:\a1\pathfinder\qual_82_09.nc';	% If not empty check againts this file (For monthly data)
			quality = 6;			% Retain only values of quality >= this (or <= abs(this) when MODIS). Ingored if fname = []
			splina = false;			% Fill missing monthly data by a spline interpolation. Ignored if fname = [].
			scale = 12;				% Scale rate of change by this value (useful when input data has monthly data).
			if (internal_master)
				calcGrad(handles, slope, sub_set, fname, quality, splina, scale)
			else
				calcGrad(handles, out{2:end})
			end
		case 'applyFlags'			% CASE 3
			fnameFlag  = 'C:\a1\pathfinder\qual_82_09.nc';
			flag = 6;				% Compute yearly means
			if (internal_master),	applyFlags(handles, fnameFlag, flag, 200)
			else,					applyFlags(handles, out{2:end})
			end
		case 'yearMeanFlag'			% CASE 4
			ano = 1:12;				% Compute yearly (ano = 1:12) or seasonal means (ano = start_month:end_month)
			fname  = 'C:\a1\MODIS\algas\Algas_qual_nsst_TERRA_00_10.nc';
			%fname  = 'C:\a1\pathfinder\qual_82_09.nc';
			quality = 0;			% Retain only values of quality >= this (or <= abs(this) when MODIS)
			nCells = 200;			% Holes (bad data) smaller than this are filled by interpolation
			% Where to save track of filled holes. Ignored if nCells = 0 OR fname3 = []
			fname3 = [];	%'C:\a1\pathfinder\qual7_85_07_Interp200_Q6.nc';
			%splina = true;
			splina = [12 30];		% Fill missing monthly data by a spline interpolation taken over two years (out limits set to NaN)
			tipoStat = 0;			% 0, Compute MEAN, 1 -> Median; 2 -> MINimum; 3 -> MAXimum; 4 -> STD of the ANO period
			% If not empty, it must contain the name of a Lon,Lat file with locations where to output time series
			chkPts_file = [];	%chkPts_file = 'C:\a1\pathfinder\chkPts.dat';
			if (internal_master)
				calc_yearMean(handles, ano, fname, quality, nCells, fname3, splina, tipoStat, chkPts_file)
			else
				calc_yearMean(handles, out{2:end})
			end
		case 'polygAVG'				% CASE 5
			fnameOut = [];			% If not empty, file name where to save the result (otherwise, asked at the end)
			op = [];				% Type average (or other). [] means doing average. 'median' will do a median
			fnamePolys = [];		% If not empty, file name of polygon or list of polygons (otherwise, fished from Mirone fig)
			sub_set = [0 0];		% [jump_start stop_before_end], make it [] or [0 0] to be ignored
			fnameFlag  = 'C:\a1\pathfinder\qual_82_09.nc';	% If not empty check againts this file
			quality = 6;			% Retain only values of quality >= this (or <= abs(this) when MODIS). Ingored if fnameFlag = []
			if (internal_master),	calc_polygAVG(handles, fnameOut, op, fnamePolys, sub_set, fnameFlag, quality)
			else,					calc_polygAVG(handles, out{2:end})
			end
		case 'flagsStats'			% CASE 6
			ano = 1:12;				% Compute yearly stats
			%opt = '';				% Make the counting on a per month basis
			opt = 'per_year';		% Make the counting on a per year basis
			if (internal_master),	calc_flagsStats(handles, ano, 7, opt)
			else,					calc_flagsStats(handles, out{2:end})
			end
		case 'pass_by_count'		% CASE 7
			count = 11;
			fname = 'C:\a1\pathfinder\countPerYear_flag7_Interp200.nc';
			if (internal_master),	pass_by_count(handles, count, fname)
			else,					pass_by_count(handles, out{2:end})
			end
		case 'do_math'				% CASE 8
			opt = 'diffstd';		% Sum all layers
			if (internal_master)
				do_math(handles, opt, [19 0], 'C:\a1\MODIS\mediaAnual_TERRA_NSST_Interp200_Q0.nc', [1 1])
			else
				do_math(handles, out{2:end})
			end
		case 'conv2vtk'				% CASE 9
			write_vtk(handles)
		case 'L2_periods'			% CASE 10
			period = 3;				% Number of days in each period
			regMinMax = [0 inf];	% If want to limit the admissible Z values ([min max])
			tipoStat = 0;			% 0, Compute MEAN, 1 compute MINimum and 2 compute MAXimum of the ANO period
			grd_out = 'C:\SVN\mironeWC\tmp\lixoL2.nc';
			if (internal_master)
				calc_L2_periods(handles, period, tipoStat, regMinMax, grd_out)
			else
				calc_L2_periods(handles, out{2:end})
			end
		case 'corrcoef'				% CASE 11
			secondArray = 'C:\a1\pathfinder\lixoCloros.nc';
			grd_out = 'C:\SVN\mironeWC\tmp\lixoR.nc';
			sub_set = [0 0];		% [jump_start stop_before_end], make it [] or [0 0] to be ignored
			if (internal_master),	calc_corrcoef(handles, secondArray, sub_set, false, grd_out)
			else,					calc_corrcoef(handles, out{2:end})
			end
		case 'count_blooms'			% CASE 12
			if (internal_master),	count_blooms(handles)
			else,					count_blooms(handles, out{2:end})
			end
	end

% ----------------------1-------2--------3---------4---------5---------6-----------7----------8---------9----
function out = zonal(handles, dlat, integ_lon, do_trends, sub_set, fnamePoly1, fnamePoly2, fnameFlag, quality)
% Compute zonal means from a multi-layer file
%
% DLAT 			width of the box in the direction orthogonal to INTEG_LON
% INTEG_LON 	(LOGICAL) If true, integration is done along longitude
% DO_TRENDS		(LOGICAL) If false compute zonal integrations, otherwise compute trends of the zonal integrations (per DLAT)
%				Note that X coords represent different things in the above two cases:
%				- Layer number for the zonal integration case and is up tp the user to make it correspond to a date.
%				- Spatial coordinate orthogonal to the integration direction for the DO_TRENDS == TRUE case
%
% OPTIONS: Since the number of options is variable some make mandatory that prev args exist. In that case use [] if needed
%
% SUB_SET	->  A two columns row vec with number of the offset of years where analysis start and stop.
%				For example [3 1] Starts analysis on forth year and stops on the before last year.
%				[0 0] Means using the all dataset.
%
% FNAMEPOLY1	Optional polygon file delimiting an area where the analysis will be carried on.
%
% FNAMEPOLY2	Optional second polygon file. If it points to a valid file. This function is called 
%				twice and results are subtracted
%
% FNAMEFLAG		name of a netCDF file with quality flags. Obviously this file must be of
%				the same size as the series under analysis. If not provided no quality check is done.
%
% QUALITY		Threshold quality value. Only values of quality >= FLAG will be taken into account
%				NOTE: For MODIS use negative FLAG. Than, values are retained if quality <= abs(FLAG)

	if (nargin < 4)
		errordlg('ZONAL: called with less than ninimum number of arguments', 'Error'),	return
	end
	if (nargin == 4)
		sub_set = [0 0];		fnamePoly1 = [];	fnamePoly2 = [];	fnameFlag = [];		quality = 7;
	elseif (nargin == 5)
		fnamePoly1 = [];		fnamePoly2 = [];	fnameFlag = [];		quality = 7;
	elseif (nargin == 6)
		fnamePoly2 = [];		fnameFlag = [];		quality = 7;
	elseif (nargin == 7)
		fnameFlag = [];			quality = 7;	% This 'quality' def is only to not error in one case.
	end

	if (numel(sub_set) == 2)
		jump_start = sub_set(1);		stop_before_end = sub_set(2);
	else
		jump_start = 0;					stop_before_end = 0;
	end

	if (~isempty(fnamePoly1))
		if (exist(fnamePoly1,'file') ~= 2)		% If given name does not exist, give another chance
			[FileName,PathName] = put_or_get_file(handles, ...
				{'*.dat;*.DAT', 'Data files (*.dat)';'*.*', 'All Files (*.*)'},'Enter polygon file','get');
			if (isequal(FileName,0)),		return,		end
			fnamePoly1 = [PathName FileName];
		end
		S = load(fnamePoly1);
		x = S(:,1);		y = S(:,2);
	else
		x = [];			y = [];
	end

	[z_id, s, rows, cols] = get_ncInfos(handles);

	%------------- Check for quality flags request -------------------
	if (~isempty(fnameFlag))
		[s_flags, z_id_flags, msg] = checkFlags_compat(fnameFlag, handles.number_of_timesteps, rows, cols);
		if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
		do_flags = true;
		if (quality > 0 ),	growing_flag = true;		% PATHFINDER flags
		else,				growing_flag = false;		% MODIS flags
		end
	else
		do_flags = false;
	end

	% ------------- Build the vectors to deal with the zonal integration -------------
	if (integ_lon)
		N_spatialSize = rows;		% Number of points in the spatial dim
		integDim = 2;							% Dimension along which we are going to integrate
		ini = fix(handles.head(3) / dlat) * dlat;
		fim = fix(handles.head(4) / dlat) * dlat + dlat;
		vecD = (ini:dlat:fim);
		Y = linspace(handles.head(3),handles.head(4), N_spatialSize);
	else
		N_spatialSize = cols;
		integDim = 1;
		ini = fix(handles.head(1) / dlat) * dlat;
		fim = fix(handles.head(2) / dlat) * dlat + dlat;
		vecD = (ini:dlat:fim);
		Y = linspace(handles.head(1),handles.head(2), N_spatialSize);
	end
	nStripes = numel(vecD) - 1;
	indStripe = ones(numel(vecD),1);

	for (k = 2:nStripes)
		ind = find(Y >= vecD(k));
		if (~isempty(ind)),		indStripe(k) = ind(1);	end
	end
	indStripe(end) = N_spatialSize;

	aguentabar(0,'title','Computing zonal means','CreateCancelBtn')

	nSeries = handles.number_of_timesteps - (jump_start + stop_before_end);	% Number of layers to be used in this run
	series_vec = (jump_start:(nSeries - 1 + jump_start)) + 1;		% Add 1 so it never starts at 0 (no good for indices)
	allSeries = zeros(nStripes, nSeries);
	if (integ_lon),		N_tot = cols + 1e-10;		% Add eps so that we never have divisions by zero
	else,				N_tot = rows + 1e-10;
	end
	mask = [];

	for (k = series_vec)			% Loop over time layers
		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [k-1 0 0], [1 rows cols]);
		this_has_nans = false;

		if (do_flags)
			flags = nc_funs('varget', fnameFlag, s_flags.Dataset(z_id_flags).Name, [k-1 0 0], [1 rows cols]);
			if (growing_flag),		Z(flags < quality) = NaN;	% Pathfinder style (higher the best) quality flag
			else,					Z(flags > quality) = NaN;	% MODIS style (lower the best) quality flag
			end
		end

		% NaNify polygon exterior points?
		if (~isempty(fnamePoly1) && k == series_vec(1))
			mask = img_fun('roipoly_j',handles.head(1:2),handles.head(3:4),double(Z),x,y);
		end
		if (~isempty(mask)),	Z(~mask) = NaN;		end

		ind = isnan(Z);						% This may, or may not, be equal to 'mask'
		if (any(ind(:))),		Z(ind) = 0;		this_has_nans = true;		end

		tmp = sum(Z, integDim);				% Add along integration dim
		if (this_has_nans)
			tmp2 = sum(ind,integDim);		% Get total number of NaNs along interp dim
			tmp = tmp ./ (N_tot - tmp2);	% Now get the number of valid values along interp dim
		else
			tmp = tmp / N_tot;
		end
		% Now add all inside each stripe
		for (m = 1:nStripes)
			tmp2 = tmp(indStripe(m):indStripe(m+1));
			tmp2(tmp2 == 0) = [];			% Not so unlikely
			allSeries(m,k) = sum(tmp2) / numel(tmp2);
		end

		h = aguentabar(k/nSeries);
		if (isnan(h)),	break,	end
	end

	if (isnan(h)),	return,		end	

	allSeries(allSeries == 0) = nan;		% NaN is more reasonable to denote data absence

	if (~isempty(fnamePoly2) && exist(fnamePoly2,'file') == 2)
		out2 = zonal(handles, dlat, integ_lon, false, sub_set, fnamePoly2, [], fnameFlag, quality);
		allSeries = double(out2) - allSeries;
	end
	allSeries = single(allSeries);

	% ------------ If no argout, show result in a Mirone/Ecran window ------------
	if (~nargout)
		zz = grdutils(allSeries,'-L');
		head = [1 nSeries vecD(1) vecD(end) zz(1) zz(2) 0 1 dlat];
		clear tmp;				% To shut up a useless ML warning due to what will happen at next line
		tmp.X = 1:nSeries;		tmp.Y = linspace((vecD(1)+dlat/2), (vecD(end)-dlat/2), nStripes);
		if (~do_trends)			% 2D, Mirone
			tmp.head = [head(1:2) tmp.Y(1) tmp.Y(end) head(5:end)];
			tmp.geo = 0;		tmp.name = 'Zonal integration';
			mirone(allSeries, tmp)
		else					% 1D, Ecran
			trend = zeros(1,nStripes);
			for (k = 1:nStripes)
				p = polyfit(tmp.X, double(allSeries(k,:)), 1);
				trend(k) = p(1);
			end
			ind = find(~isnan(trend));			% Remove any eventual leading or trailing NaNs
			trend = trend(ind(1):ind(end));
			tmp.Y = tmp.Y(ind(1):ind(end));
 			ecran(handles, tmp.Y, trend, 'Slope of line fit')
		end
	else
		out = allSeries;
	end
	
% -------------------1-------2-------3---------4---------5-------6-------7-------8-------9--------10-----
function calcGrad(handles, slope, sub_set, fnameFlag, quality, splina, scale, grd_out, do_3x3, mask_file)
% Compute the rate of change of a file by fitting a LS straight line. The file can be one of the
% already computed yearly means, in which case last three input arguments do not apply,
% OR the full time series. In this case optional checking against quality flags and spline
% interpolation can be done using info transmitted via the last three args.
%
% NOTE: THIS IS A HIGHLY MEMORY CONSUMPTION ROUTINE AS ALL DATA IS LOADED IN MEMORY
%
% SLOPE		Logical indicating if compute slope of linear fit (TRUE) or the p parameter (FALSE)
%
% SUB_SET -> A two elements row vec with number of the offset of years where analysis start and stop.
%			For example [3 1] Starts analysis on forth year and stops on the before last year.
%			[0 0] Means using the all dataset.
%
% OPTIONS:
% FNAMEFLAG	name of a netCDF file with quality flags. Obviously this file must be of
% 			the same size as the series under analysis.
%
% QUALITY	Threshold quality value. Only values of quality >= FLAG will be taken into account
%			NOTE: For MODIS use negative FLAG. Than, values are retained if quality <= abs(FLAG)
%
% SPLINA	Logical that if true instruct to spline interpolate the missing monthly values
%			before computing the rate of change.
%
% SCALE		Scale the final rate by this value. The idea of all this is that input data
%			can have monthly means and we want to compute time rate of change per year.
%			In this case use SCALE = 12.
%			Ignored if SLOPE is FALSE
%
% GRD_OUT	Name of the netCDF file where to store the result. If not provided, open Mirone Fig.
%
% DO_3x3	Logical. If true, compute over 3x3 windows
%
% MASK_FILE	name of Land mask file. If it has NaNs it will be multiplied but if only 1/0's, 0's will become NaNs

	if (nargin <= 6),	scale = 1;		end
	if (nargin < 8),	grd_out = [];	end
	if (nargin < 9),	do_3x3 = false;	end
	if (nargin < 10),	mask_file = '';	end
	if (~slope),		scale = 1;		end		% Make sure to not scale p-values

	do_blockMean = false;	% Must be turned into an input parameter
	do_flags = false;		% Will be set to true if we do a checking against a quality flgas file
	get_profiles_in_polygon = false;			% Save all profiles (along third dim) located inside the polygonal area
	n_anos = handles.number_of_timesteps;
	[z_id, s, rows, cols] = get_ncInfos(handles);

	% Find if we are dealing with a Pathfinder V5.2 daily file
	[is_PFV52, tempos] = PFV52(handles, s);		% if not a PFV5.2, tempos = handles.time(:)

	if (nargin >= 3 && (numel(sub_set) == 2))
		jump_anos = sub_set(1);		stop_before_end_anos = sub_set(2);
	else
		jump_anos = 0;				stop_before_end_anos = 0;
	end

	% For PFV5.2 daily data, accept the start-stop as time in years. e.g. [1986 2011]
	if (is_PFV52 && jump_anos >= 1982 && stop_before_end_anos <= 2014)	% Though we only have them till 2011
		ind = find(tempos > jump_anos);
		jump_anos = ind(1);
		ind = find(tempos > stop_before_end_anos);
		stop_before_end_anos = numel(tempos) - ind(1);
	end

	n_anos = n_anos - (jump_anos + stop_before_end_anos);	% Number of layers to be used in this run

	% When files are too big (not difficult when number of layers is high) must do a patch processing
	one_Gb = 1024*1024*1024;		% Use patches of 1 Gb
	if (rows * cols * n_anos * 4 > one_Gb)
		slice_cols = round(one_Gb / (rows * n_anos * 4));
		slice_cols = 1:slice_cols:cols;
		if (slice_cols(end) < cols),	slice_cols(end+1) = cols;	end
		slicing = true;
	else
		slice_cols = [1 cols];
		slicing = false;
	end
	read_cols = cols;

	if (nargin >= 4 && ~isempty(fnameFlag))
		if (slicing && ~handles.IamCompiled)
			slice_cols = [1 cols];
			slicing = false;
			warndlg('File very big that would be processed by patches, but that is not implemented whid Flags. Trying wirh a single patch.', 'Warning')
		elseif (slicing && handles.IamCompiled)
			errordlg('File too big. Processing it with quality flags is not implemented','Error'),	return
		end
		[s_flags, z_id_flags, msg] = checkFlags_compat(fnameFlag, handles.number_of_timesteps, rows, cols);
		if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
		flags = alloc_mex(rows, cols, n_anos, 'uint8');
		for (m = 1:n_anos)
			flags(:,:,m) = nc_funs('varget', fnameFlag, s_flags.Dataset(z_id_flags).Name, [(m - 1 + jump_anos) 0 0], [1 rows cols]);
		end
		do_flags = true;
		if (nargin == 4)			% Default to Pathfinder max quality
			quality = 7;	splina = false;
		elseif (nargin == 5)
			splina = false;
		end
		if (quality > 0 ),	growing_flag = true;		% PATHFINDER flags
		else,				growing_flag = false;		% MODIS flags
		end
	else
		splina = false;		flags = [];	growing_flag = false;	% Not used
	end

	if (is_PFV52)			% For Pathfinder V5.2 (daily) use the true time coordinates
		x = tempos;		yy = [];
		if (splina),	yy = tempos;			end
		threshold = 0.10;	% AD-HOC threshold percentage of the n points below which do NOT compute slope
	else
		x = (0:n_anos-1)';	yy = [];
		if (splina),	yy = (0:n_anos-1)';		end
		threshold = 0.66;	% 66%
	end

try
	nSlices = numel(slice_cols) - 1;
	for (ns = 1:nSlices)		% Loop over number of slices (which may be only one)
		msg = 'Loading data ...';
		if (slicing)
			read_cols = slice_cols(ns+1) - slice_cols(ns) + 1;
			if (ns == 1)
				Tvar = zeros(rows, cols);	% To gether the final result. Allocate on first usage
			end
			msg = sprintf('Loading (big) data. %d of %d ...', ns, nSlices);
		end
		Tmed = alloc_mex(rows, read_cols, n_anos, 'single');
		aguentabar(0,'title', msg)
		if (~slicing)
			for (m = 1:n_anos)
		 		Tmed(:,:,m) = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [(m - 1 + jump_anos) 0 0], [1 rows cols]);
			end
		else
			for (m = 1:n_anos)
		 		t = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [(m - 1 + jump_anos) 0 0], [1 rows cols]);
				Tmed(:, :, m) = t(:, slice_cols(ns):slice_cols(ns+1));
			end
			clear t
		end

		% ---- save profiles of points, located inside polygon of Mirone fig, as a multi-segment file
		if (get_profiles_in_polygon)
			profiles_in_polygon(handles, Tmed, n_anos)
			return
		end
		% -----------------------------------------------------------------------------------------

		if (do_blockMean)
			aguentabar(0,'title','Compute block means','CreateCancelBtn')
			for (m = 1:n_anos)
				if (do_flags)
					slice = Tmed(:,:,m);
					fslice = flags(:,:,m);
					if (growing_flag),	slice(fslice < quality) = NaN;
					else,				slice(fslice > quality) = NaN;
					end
					Tmed(:,:,m) = mirblock(slice, '-A3');
				else
					Tmed(:,:,m) = mirblock(Tmed(:,:,m), '-A3');
				end
				h = aguentabar(m/n_anos);
				if (isnan(h)),	return,		end
			end
			if (do_flags)
				clear slice fslice
			end
		end

		if (is_PFV52)
			aguentabar(0,'title','Computing and removing Seazonal cycle','CreateCancelBtn')
			Tavg = tideman(handles, Tmed, tempos, 52);
			Tmed = remove_seazon(handles, Tmed, Tavg, tempos);
			clear Tavg
		end

		if (~slicing)			%Do it all in one passage
			if (do_3x3)
				[Tvar, stopit] = get_slopes_conn8(Tmed, flags, x, rows, cols, slope, do_flags, quality, growing_flag, threshold);
			else
				[Tvar, stopit] = get_slopes(Tmed, flags, x, yy, rows, cols, slope, do_flags, quality, growing_flag, splina, threshold);
			end
		else
			if (do_3x3)
				[Tvar_slice, stopit] = get_slopes_conn8(Tmed, flags, x, rows, read_cols, slope, do_flags, quality, growing_flag, threshold);
			else
				[Tvar_slice, stopit] = get_slopes(Tmed, flags, x, yy, rows, read_cols, slope, do_flags, quality, growing_flag, splina, threshold);
			end
			if (stopit),	return,		end
			clear Tmed
			Tvar(:, slice_cols(ns):slice_cols(ns+1)) = Tvar_slice;
		end
	end
catch
	disp(lasterror)
end

	clear Tmed
	if (stopit),	return,		end

	if (scale ~= 1),		cvlib_mex('CvtScale', Tvar, double(scale),0);		end
	Tvar = single(Tvar);

	if (~isempty(mask_file))
		Tvar = apply_mask(mask_file, [], Tvar);		% Apply the Land mask file
	end

	zz = grdutils(Tvar,'-L');  handles.head(5:6) = [zz(1) zz(2)];
	tmp.head = handles.head;
	if (isempty(grd_out))			% Show result in a Mirone figure
		tmp.X = linspace(tmp.head(1),tmp.head(2),cols);
		tmp.Y = linspace(tmp.head(3),tmp.head(4),rows);
		tmp.name = 'Time gradient (deg/year)';
		mirone(Tvar, tmp)
	else							% Got output name from input arg
		handles.was_int16 = false;
		handles.computed_grid = true;
		handles.geog = 1;
		nc_io(grd_out, 'w', handles, Tvar)
	end

% --------------------------------------------------------------------------------------
function [Tvar, stopit] = get_slopes(Tmed, flags, x, yy, rows, cols, slope, do_flags, quality, growing_flag, splina, threshold)
% Compute the best fit slopes (or p-values) on each node
	Tvar = zeros(rows, cols) * NaN;
	stopit = false;
	n_anos = numel(x);
	atLeast = round(n_anos * threshold);
	aguentabar(0,'title','Computing the Time rate','CreateCancelBtn')

	try
	for (n = 1:cols)
		for (m = 1:rows)
			y = double(squeeze(Tmed(m,n,:)));
			if (do_flags)
				this_flag = squeeze(flags(m,n,:));
				if (growing_flag),		y(this_flag < quality) = NaN;	% Pathfinder style (higher the best) quality flag
				else,					y(this_flag > quality) = NaN;	% MODIS style (lower the best) quality flag
				end
			end
			ind = isnan(y);
			y(ind) = [];
			if (numel(y) < atLeast),	continue,	end	% Completely ad-hoc test (it also jumps land cells)

			if (splina)
				if (~all(ind))		% Otherwise we have them all and so nothing to interp
					akimaspline(x(~ind), y, x, yy);
					y = yy;										
					if (ind(1) || ind(end))				% Cases when where we would have extrapolations
						if (ind(1))
							ki = 1;
							while (ind(ki)),	ki = ki + 1;	end
							y(1:ki) = NaN;				% Reset extraped values to NaN
						end
						if (ind(end))
							kf = numel(ind);
							while (ind(kf)),	kf = kf - 1;	end
							y(kf:numel(ind)) = NaN;		% Reset extraped values to NaN
						end
						ind = isnan(y);					% Get new nan indices again
						y(ind) = [];
					else
						ind = false(n_anos,1);			% Pretend no NaNs for the rest of the code below
					end
				end
			end

 			%p = polyfit(x(~ind),y,1);
			%z=[xvalues(1:4);ones(1,4)]'\yvalues';
			if (slope)		% Compute sople of linear fit
				p = trend1d_m([x(~ind) y],'-L','-N2r');
				Tvar(m,n) = p(1);
			else			% Compute p value
	 			p = trend1d_m([x(~ind) y],'-L','-N2r','-R','-P');
				if (p(1) < -0.5 || p(1) > 1),	continue,	end		% Another ad-hoc (CLIPPING)
				Tvar(m,n) = p(4);
			end
		end
		h = aguentabar(n/cols);
		if (isnan(h)),	stopit = true;	break,	end
	end
	aguentabar(1)
	catch
		disp(lasterror)
	end

% --------------------------------------------------------------------------------------
function [Tvar, stopit] = get_slopes_conn8(Tmed, flags, x, rows, cols, slope, do_flags, quality, growing_flag, threshold)
% Compute the best fit slopes (or p-values) centered on each node but using a 8-connection.
% That is, on each node the value will be an average of all the 8 data points arround, plut itself.
	Tvar = zeros(rows, cols) * NaN;
	stopit = false;
	aguentabar(0,'title','Computing the Time rate','CreateCancelBtn')

	% First do the left and right columns with the connection 1 algo. Won't do the Top & Bottom rows
	atLeast = round(numel(x) * threshold);		% Minimum number of points needed to do the fit
	for (n = [1 cols])
		for (m = 2:rows-1)
			y = double(Tmed(m,n,:));
			y = y(:);
			if (do_flags)
				this_flag = flags(m,n,:);
				if (growing_flag),		y(this_flag < quality) = NaN;	% Pathfinder style (higher the best) quality flag
				else,					y(this_flag > quality) = NaN;	% MODIS style (lower the best) quality flag
				end
			end
			ind = isnan(y);
			y(ind) = [];
			if (numel(y) < atLeast),	continue,	end				% Completely ad-hoc test (it also jumps land cells)

			if (slope)		% Compute sople of linear fit
				p = trend1d_m([x(~ind) y],'-L','-N2r');
				Tvar(m,n) = p(1);
			else			% Compute p value
	 			p = trend1d_m([x(~ind) y],'-L','-N2r','-R','-P');
				if (p(1) < -0.5 || p(1) > 1),	continue,	end		% Another ad-hoc (CLIPPING)
				Tvar(m,n) = p(4);
			end
		end
	end
	
	atLeast = round(numel(x) * threshold * 9);		% 9 because we are doing  a 3x3 window
	x = repmat(x, 1, 9)';
	x = x(:);
	for (n = 2:cols-1)
		for (m = 2:rows-1)
			y = double(Tmed([m-1 m m+1],[n-1 n n+1],:));
			y = y(:);
			if (do_flags)
				this_flag =flags([m-1 m m+1],[n-1 n n+1],:);
				if (growing_flag),		y(this_flag < quality) = NaN;	% Pathfinder style (higher the best) quality flag
				else,					y(this_flag > quality) = NaN;	% MODIS style (lower the best) quality flag
				end
			end
			ind = isnan(y);
			y(ind) = [];
			if (numel(y) < atLeast),	continue,	end				% Completely ad-hoc test (it also jumps land cells)

			if (slope)		% Compute sople of linear fit
				p = trend1d_m([x(~ind) y],'-L','-N2r');
				Tvar(m,n) = p(1);
			else			% Compute p value
	 			p = trend1d_m([x(~ind) y],'-L','-N2r','-R','-P');
				if (p(1) < -0.5 || p(1) > 1),	continue,	end		% Another ad-hoc (CLIPPING)
				Tvar(m,n) = p(4);
			end
		end
		h = aguentabar(n/cols);
		if (isnan(h)),	stopit = true;	break,	end
	end
	aguentabar(1)

% --------------------------------------------------------------------------------------
function profiles_in_polygon(handles, Tmed, n_anos)
% Save profiles of points, located inside polygon of Mirone fig, as a multi-segment file

	hFigs = findobj(0,'type','figure');						% Fish all figures
	IAmAMir = zeros(1, numel(hFigs));
	for (k = 1:numel(hFigs))								% Get the first Mirone figure with something in it
		if (~isempty(getappdata(hFigs(k), 'IAmAMirone')))
			handMir = guidata(hFigs(k));
			if (handMir.no_file),	continue,	end			% A virgin Mirone bar figure
			IAmAMir(k) = 1;		break,	
		end
	end
	if (sum(IAmAMir) ~= 1)
		errordlg('Did not find any valid Mirone figure with data displayed.','Error'),	return
	end

	hLine = findobj(handMir.axes1,'Type','line');
	if (isempty(hLine)),	hLine = findobj(handMir.axes1,'Type','patch');	end		% Try once more
	x = get(hLine,'XData');		y = get(hLine,'YData');
	if (isempty(x))
		errordlg('The Mirone figure needs to have at least one polygon loaded.','Error'),	return
	end
	mask = img_fun('roipoly_j',handles.head(1:2),handles.head(3:4),get(handMir.hImg,'CData'),x,y);
	B = img_fun('find_holes',mask);
	col_min = min(B{1}(:,2));		col_max = max(B{1}(:,2));
	row_min = min(B{1}(:,1));		row_max = max(B{1}(:,1));
	row_vec = row_min:row_max;
	col_vec = col_min:col_max;
	x = (0:n_anos-1)';
	k = 1;
	stack = cell(1, 3);		% Obviously not enough but will shut up MLint
	for (n = col_vec)
		IN = inpolygon(repmat(n, numel(row_vec), 1), row_vec, B{1}(:,2), B{1}(:,1));		% See if ...
		this_row = 1;
		for (m = row_vec)
			if (~IN(this_row)),		continue,	end				% This pixel is outside polygon POI
			this_row = this_row + 1;
			y = double(squeeze(Tmed(m,n,:)));
			ind = isnan(y);
			y(ind) = [];
			if (numel(y) < n_anos/2),		continue,	end		% Completely ad-hoc test
			p = trend1d_m([x(~ind) y],'-L','-N2r','-R','-P');
			stack{k,1} = [x(~ind)+1 y];	% x,temp
			stack{k,2} = p(1);			% Slope
			stack{k,3} = p(4);			% p-value
			k = k + 1;
		end
	end

	str1 = {'*.dat;*.DAT', 'Symbol file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select Output File name','put','.dat');
	if isequal(FileName,0),		return,		end
	f_name = [PathName FileName];
	double2ascii(f_name, stack, '%.0f\t%.3f', 'maybeMultis');
	[PATH, FNAME, EXT] = fileparts(f_name);
	f_name = [PATH filesep FNAME '_mp' EXT];		% Write a second file with 2 columns where 1st col is slope and 2nth is p-value
	xy = [cat(1,stack{:,2}) cat(1, stack{:,3})];
	double2ascii(f_name, xy);

% ------------------------------------------------------------------------------
function applyFlags(handles, fname, flag, nCells, grd_out)
% Check a 3D file against its 3D flags companion and replace values < FLAG to NaN.
% Optionaly interpolate to fill gaps smaller than NCELLS.
%
% FNAME 	name of a netCDF file with quality flags. Obviously this file must be of
% 			the same size as the series under analysis.
%
% FLAG		Threshold quality value. Only values of quality >= FLAG will be taken into account
%			NOTE: For MODIS use negative FLAG. Than, values are retained if quality <= abs(FLAG)
%
% OPTIONS:
% NCELLS	Holes (bad data) groups smaller than this are filled by interpolation (default = 0). 
%			For example if NCELL = 200 groups of equal or less than a total of 200 will be filled.
%
% GRD_OUT	Name of the netCDF file where to store the result. If not provided, it will be asked here.

	if (nargin < 3),	error('calc_yearMean:Not enough input arguments'),	end

	[z_id, s, rows, cols] = get_ncInfos(handles);

	[s_flags, z_id_flags, msg] = checkFlags_compat(fname, handles.number_of_timesteps, rows, cols);
	if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
	if (nargin == 3),	nCells = 0;		end			% If not provided, defaults to no gaps fill

	if (nargin < 5)
		txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
		[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
		if isequal(FileName,0),		return,		end
		grd_out = [PathName FileName];
	end

	% -------------------------------------- END PARSING SECTION --------------------------------------------

	if (flag > 0),		growing_flag = true;				% Pathfinder style (higher the best) quality flag
	else,				growing_flag = false;	flag = -flag;	% MODIS style (lower the best) quality flag
	end
	pintAnoes = (nCells > 0);

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;
	n_layers = handles.number_of_timesteps;

	aguentabar(0,'title','Applying flags.','CreateCancelBtn');

	for (m = 1:n_layers)
		if (m == 1)						% First layer
			Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m-1 0 0], [1 rows cols]);
		else
			Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m-1 0 0], [1 rows cols]);
		end

		Z_flags = nc_funs('varget', fname, s_flags.Dataset(z_id_flags).Name, [m-1 0 0], [1 rows cols]);
		if (growing_flag),		Z(Z_flags < flag) = NaN;	% Pathfinder style (higher the best) quality flag
		else,					Z(Z_flags > flag) = NaN;	% MODIS style (lower the best) quality flag
		end

		ind = isnan(Z);
		if (pintAnoes && any(ind(:)))		% If fill spatial holes is requested
			Z = inpaint_nans(handles, Z, ind, nCells);			% Select interp method inside inpaint_nans()
		end

		% Write this layer to file
		if (m == 1),		nc_io(grd_out, sprintf('w%d/time',n_layers), handles, reshape(Z,[1 size(Z)]))
		else,				nc_io(grd_out, sprintf('w%d', m-1), handles, Z)
		end

		h = aguentabar(m/n_layers,'title','Applying flags.');	drawnow
		if (isnan(h)),	break,	end
	end

% ------------------------1-------2-------3------4------5-------6--------7-------8---------9-----------10--------11----
function calc_yearMean(handles, months, fname2, flag, nCells, fname3, splina, tipoStat, chkPts_file, grd_out, mask_file)
% Compute anual means or climatologies from monthly data (well, and daily too).
%
% MONTHS 	a vector with the months uppon which the mean is to be computed
%		example: 	months = 1:12		==> Computes yearly mean
%					months = 6:8		==> Computes June-July-August seazonal means
%			Default (if months = []) 1:12
%
% OPTIONS:
% FNAME2 	name of a netCDF file with quality flags. Obviously this file must be of
% 			the same size as the series under analysis.
%
% FLAG		Threshold quality value. Only values of quality >= FLAG will be taken into account
%			NOTE: For MODIS use negative FLAG. Than, values are retained if quality <= abs(FLAG)
%
% NCELLS	Holes (bad data) groups smaller than this are filled by interpolation. 
%			For example if NCELL = 200 groups of equal or less than a total of 200 will be filled.
%
% FNAME3 	Optional name of a netCDF file where interpolated nodes will be set to FLAG
%			and the others retain their FNAME2 value. This corresponds to the promotion
%			of interpolated nodes to quality FLAG.
%
% SPLINA	Logical that if true instructs to spline interpolate the missing monthly values
%			before computing the yearly mean. Optionaly, it may be a 2 elements vector with
%			the MIN and MAX values allowed on the Z function (default [0 32]).
%			----------------------- OR (to CLIMA) --------------------
%			A string containing 'CLIMA' to instruct to compute climatologies from montly data
%			This is a uggly hack but I don't want to add yet another input argument.
%
% TIPOSTAT	Variable to control what statistic to compute.
%			0 Compute MEAN of MONTHS period. 1 Compute MINimum and 2 compute MAXimum
%
% CHKPTS_FILE	(Optional)
%			Name of a file with Lon,Lat locations where to output the entire time series.
%			Output name file is constructed from the input name and appended '_tseries'.
%			Fitst column has the month number, even columns the original data and odd columns
%			the data with the holes (NaNs) interpolated with an Akima spline function.
%
% GRD_OUT	Name of the netCDF file where to store the result. If not provided, it will be asked here.
%
% MASK_FILE	name of Land mask file. If it has NaNs it will be multiplied but if only 1/0's, 0's will become NaNs

	% Variables that are not always used but need to exist
	Tmed = [];	ZtoSpline = [];	contanoes = [];		total_months = [];	n_pad_months = [];
	z_id_flags = [];	s_flags= [];

	do_flags = false;		track_filled = false;		do_saveSeries = false;	do_climatologies = false;
	[z_id, s, rows, cols] = get_ncInfos(handles);

	% Find if we are dealing with a Pathfinder V5.2 daily file
	is_PFV52 = PFV52(handles, s);

	if (isempty(months)),	months = 1:12;		end

	if (nargin >= 6 && ~isempty(fname3)),		track_filled = true;	end		% Keep track of interpolated nodes

	if (nargin >= 3 && ~isempty(fname2))				% We have a quality-flag ghost file to check
		[s_flags, z_id_flags, msg] = checkFlags_compat(fname2, handles.number_of_timesteps, rows, cols);
		if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
		if (nargin == 3),	nCells = 0;		flag = 7;	end			% If not provided, defaults to best quality
		do_flags = true;
	end

	if (nargin < 7)
		splina = false;		tipoStat = 0;	chkPts_file = [];	grd_out = [];
	elseif (nargin < 8)
		tipoStat = 0;	chkPts_file = [];	grd_out = [];
	end
	if (nargin < 11),	mask_file = '';		end			% Need to test this before the "chkPts_file"
	if (~isempty(mask_file)),	chkPts_file = [];	end

	% -------------- Test if output time series at locations provided in the CHKPTS_FILE --------------------
	if (nargin >= 9 && ~isempty(chkPts_file))
		if (exist(chkPts_file,'file') == 2)
			pts_pos = text_read(chkPts_file);
			indTimeSeries = zeros(size(pts_pos,1), 2);
			for (k = 1:size(pts_pos,1))
				indTimeSeries(k,1) = round((pts_pos(k,1) - handles.head(1)) / handles.head(8)) + 1;
				indTimeSeries(k,2) = round((pts_pos(k,2) - handles.head(3)) / handles.head(9)) + 1;
			end
			%timeSeries = zeros(s.Dataset(3).Size, k);		% Total number of layers
			timeSeries = [(1:s.Dataset(3).Size)' zeros(s.Dataset(3).Size, 2*k)];	% Total N of layers + N of check pts
			indTSCurr_o = 1;		indTSCurr_s = 1;

			if (rem(size(timeSeries,1), 12) ~= 0)
				warndlg('Output time series works only with complete years of monthly data. Ignoring request','Warning')
			else
				do_saveSeries = true;
			end
		end
	end
	if (nargin < 10 || isempty(grd_out))	% Note: old and simple CASE 3 in main cannot send here the output name 
		txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
		[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
		if isequal(FileName,0),		return,		end
		grd_out = [PathName FileName];
	end
	mask = [];		% If needed more than once, this var will hold the amsking array

	% -------------------------------------------------------------------------------------------------------
	% -------------------------------------- END PARSING SECTION --------------------------------------------
	% -------------------------------------------------------------------------------------------------------

	if (flag > 0),		growing_flag = true;				% Pathfinder style (higher the best) quality flag
	else,				growing_flag = false;	flag = -flag;	% MODIS style (lower the best) quality flag
	end
	pintAnoes = (nCells > 0);

	% The following limits are used to clip unreasonable temperatures computed during the spline interpolation
	% When no time (spline) interpolation is used, they are simply ignored
	if (isa(splina, 'char'))
		if (strcmpi(splina,'CLIMA'))
			do_climatologies = true;
		end
		splina = false;					% Make it logic again for use in IF tests
	elseif (numel(splina) == 2)
		regionalMIN = splina(1);		regionalMAX = splina(2);
		splina = true;					% Make it logic again for use in IF tests
	else
		regionalMIN = 0;				regionalMAX = 32;
	end

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;

	if (rem(handles.number_of_timesteps, 12) == 0 && handles.number_of_timesteps < 336)	% SHIT
		n_anos = handles.number_of_timesteps / 12;
	else
		n_anos = 1;		% TEMPORARY	FOR COMPUTING YEARLY MEANS FROM DAILY DATA --- NON SECURED AND NON DOCUMENTED
	end
	if (is_PFV52)		% This relies on the fact that the apropriate 'description' global attribute has been set
		anos = fix(handles.time);
		n_anos = max(anos) - min(anos) + 1;
		if (n_anos == 2 && anos(1) ~= anos(2))		% Crazy NOAAs have year of 2006 start a 2005.99xxx
			n_anos = 1;
		end

		if (months(end) > n_anos && months(end) < 20)	% A crude test to detect a year-request-overflow (it may screw)
			months = months(1):n_anos;
			if (isempty(months))
				errordlg('Nonsense: Requested period is not covered by the input data','Error'),		return
			else
				warndlg('You requested more years than the data actually has. Reseting max to max years in data.','Warning')
			end
		end
		if (months(end) <= n_anos)					% Convert a period given in number of years to true years.
			months = (anos(1) - months(1) + 1):(anos(end) - months(end) + 1);
		end
		if (months(1) >= anos(1) && months(end) <= anos(end))	% In this case 'months' are actually the years
			new_months = [];
			for (k = 1:numel(months))			% So compute the new 'months' vector
				this_year = months(k);
				ind = find(anos == this_year);
				new_months = [new_months(1:end) ind(:)'];
			end
			months = new_months;
		end
	end

	% Take care of the case when output file has no extension
	[p,f,e] = fileparts(grd_out);
	if (isempty(e))
		if (n_anos == 1),	e = '.grd';
		else,				e = '.nc';
		end
		grd_out = [p filesep f e];
	end

	if (splina)
		n_pad_months = 7;		% Example: if months = 7:9 interpolation domain is 7-n_pad_months-1:9+n_pad_months
		if (n_anos == 1)
			n_pad_months = 0;
			warndlg('Time series is too short to do the spline interpolation option.','WARNING')
		end
		total_months = numel(months) + 2*n_pad_months;
		ZtoSpline = alloc_mex(rows, cols, total_months, 'single', NaN);
	end

	aguentabar(0,'title','Computing means.','CreateCancelBtn');
	if (~splina)
		Tmed = zeros(rows, cols);		% Temp media para cada um dos anos
	end
	in_break = false;					% Inner loop cancel option
	last_processed_month = 0;		already_processed = 0;

	if (do_climatologies)
		multi_climat = false;		n_outer_loop = 1;
		if (months(1) > 100)			% Compute a climat for each element in 'months'
			multi_climat = true;	n_outer_loop = numel(months);
		end
	else
		n_outer_loop = n_anos;
	end

	for (m = 1:n_outer_loop)

		if (splina)
			if (m == 1)
				past_months = max(1,months(1)-n_pad_months) - 1;
				this_months = past_months+1 : (months(end)+n_pad_months);
			elseif (m == n_anos)
				past_months = (m - 1)*12 + months(1) - n_pad_months - 1;
				this_months = past_months+1 : past_months+(n_pad_months + numel(months) + min(n_pad_months, 12-months(end)));
			else
				past_months = (m - 1)*12 + months(1) - n_pad_months - 1;
				this_months = past_months+1 : past_months+(numel(months)+2*n_pad_months);
			end
		elseif (do_climatologies)
			% Litle manip to create a row vector with all months of interest across the years
			if (multi_climat)
				this_months = repmat(months(m)-100, n_anos, 1);
			else
				this_months = repmat(months, n_anos, 1);
			end
			v = ((0:n_anos) * 12)';
			for (k = 1:n_anos)
				this_months(k,:) = this_months(k,:) + v(k);
			end
			this_months = this_months';
			this_months = this_months(:)';
			contanoes = zeros(rows, cols);
		else
			this_months = (m - 1) * 12 + months;
			contanoes = zeros(rows, cols);
		end

		if (tipoStat > 0)
			l_inc = this_months(2) - this_months(1);
			s2 = struct('fname',handles.fname, 'info',s, 'n_layers', handles.number_of_timesteps, 'rows',rows, ...
				 'cols',cols, 'layerOI',this_months(1), 'layer_inc',l_inc, 'z_id',z_id);
			tmp = doM_or_M_or_M([], this_months(1), l_inc, handles.number_of_timesteps, Inf, Inf, tipoStat, s2);
		else
			% For averages and for the time being (not break compat), continue to use the old code in form of function
			[Tmed, contanoes, ZtoSpline, already_processed] = ...
				calc_average_old(handles, s, s_flags, Tmed, ZtoSpline, contanoes, this_months, splina, last_processed_month, ...
				nCells, n_anos, already_processed, total_months, months, n_pad_months, do_flags, flag, ...
				fname2, fname3, pintAnoes, track_filled, rows, cols, z_id, z_id_flags, growing_flag, m);
			if (~splina)					% Do not interpolate along time. Compute averages with all non NaNs
				cvlib_mex('div', Tmed, contanoes);		% The mean for current year
				tmp = single(Tmed);
			end
		end

		last_processed_month = this_months(end);

		if (tipoStat == 0 && splina)					% Fill missing month data by interpolation based on non-NaN data
			hh = aguentabar(eps,'title','Splining it.');	drawnow
			if (isnan(hh)),		break,		end			% Over time loop said: break
			n_meses = numel(this_months);

			if (m == 1),	first_wanted_month = months(1);				% First year in the stack
			else,			first_wanted_month = n_pad_months + 1;
			end
			last_wanted_month = first_wanted_month + numel(months) - 1;

			yy = zeros(1, n_meses)';				% A kind of pre-allocation to be used by the akimaspline MEX

			% ----------- Test if we are saving time series in array for later saving ----------------------
			if (do_saveSeries)
				[timeSeries, indTSCurr_o] = getTimeSeries(ZtoSpline, timeSeries, indTimeSeries, indTSCurr_o, ...
											true, first_wanted_month, last_wanted_month);
			end
			% ----------------------------------------------------------------------------------------------

			for (i = 1:cols)
				for (j = 1:rows)
					if (already_processed),		break,	end
					y = double(squeeze(ZtoSpline(j,i,1:n_meses)));
					ind = ~isnan(y);
					% If have NaNs inside months of interest and the overall series has enough points, interp in the missing positions
					if (all(ind(first_wanted_month:last_wanted_month)))		% We have them all, so nothing to interp
						ZtoSpline(j,i,1:n_meses) = single(y);
					elseif (~any(ind))		% They are all NaNs -- Almost sure a land pixel
						continue
					elseif (numel(ind(~ind)) <= round((numel(this_months) - (n_pad_months * (m ~= 1))) / 2))
						% At least > 1/2 number of valid pts not counting first n_pad_months that were already interpolated
						x = this_months(ind);			y0 = y(ind);
 						akimaspline(x, y0, this_months, yy);
						y(first_wanted_month:last_wanted_month) = yy(first_wanted_month:last_wanted_month);
						if (~ind(1) || ~ind(end))			% Cases when we would have extrapolations
							if (~ind(1))
								ki = 1;
								while (~ind(ki)),	ki = ki + 1;	end
								y(1:ki) = NaN;				% Reset extraped values to NaN
							end
							if (~ind(end))
								kf = numel(ind);
								while (~ind(kf)),	kf = kf - 1;	end
								y(kf:numel(ind)) = NaN;		% Reset extraped values to NaN
							end
						end
						ZtoSpline(j,i,1:n_meses) = single(y);
					else									% Less than half valid points. We'll make them be all NaNs
						if (~tipoStat),	 ZtoSpline(j,i,1:n_meses) = y * NaN;	end		% For MIN & MAX we want to retain data
					end
				end
				hh = aguentabar(i/(cols+1));	drawnow
				if (isnan(hh)),		in_break = true;	break,		end		% Over time loop said: break (comment: FCK ML SHIT)
			end				% End loops over this 2D layer

			if (in_break),		break,		end			% Fck no gotos paranoia obliges to this recursive break

			% ----------- Test if we are saving time series in array for later saving ----------------------
			if (do_saveSeries)
				[timeSeries, indTSCurr_s] = getTimeSeries(ZtoSpline, timeSeries, indTimeSeries, indTSCurr_s, ...
											false, first_wanted_month, last_wanted_month);
			end
			% ----------------------------------------------------------------------------------------------

			% Now we can finaly compute the season MEAN or MIN or MAX
			tmp = doM_or_M_or_M(ZtoSpline, first_wanted_month, 1, last_wanted_month, regionalMIN, regionalMAX, tipoStat);
		end							% End interpolate along time (SPLINA)

		tmp(tmp == 0) = NaN;		% Reset the NaNs

		if (~isempty(mask_file))
			[tmp, mask] = apply_mask(mask_file, mask, tmp);		% First time reads from MASK_FILE, second on uses MASK
		end

		if (in_break),		break,		end		% Fckng no gotos paranoia obliges to this recursive break

% 		% Clip obvious bad data based on cheap median statistics
% 		aguentabar(0.5,'title','Filtering obvious bad data based on cheap statistics.');	drawnow
% 		medianas = c_grdfilter(tmp,[1 size(tmp,2) 1 size(tmp,1) 0 50 0 1 1], '-D0', '-Fm11');
% 		difa = cvlib_mex('absDiff', tmp, medianas);
% 		tmp(difa > 0.75) = NaN;					% 0.75 is probably still too permissive

		% Write this layer to file
		if (m == 1)
			if (n_outer_loop == 1)		% If one single layer don't make it 3D
				% Defaults and srsWKT fishing are set in nc_io
				misc = struct('x_units',[],'y_units',[],'z_units',[],'z_name',[],'desc',[], ...
					'title','Yearly or Seazonal average','history',[],'srsWKT',[], 'strPROJ4',[]);
				zz = grdutils(tmp,'-L');
				handles.head(5:6) = [double(zz(1)) double(zz(2))];
				nc_io(grd_out, 'w', handles, tmp, misc)
			else
				nc_io(grd_out, sprintf('w%d/time',n_outer_loop), handles, reshape(tmp,[1 size(tmp)]))
			end
		else
			nc_io(grd_out, sprintf('w%d', m-1), handles, tmp)
		end

		h = aguentabar(m/n_outer_loop,'title','Computing means.');	drawnow
		if (isnan(h)),	break,	end

		if (~splina && m < n_outer_loop),	cvlib_mex('CvtScale', Tmed, 0.0, 0.0);	end			% Reset it to zeros
	end

	if (do_saveSeries)			% Save the time series file. The name is build from that of locations file
		[pato, fname, ext] = fileparts(chkPts_file);
		fname = [fname '_tseries' ext];
		if (~isempty(pato)),	fname = [pato filesep fname];	end
		double2ascii(fname, timeSeries, ['%d' repmat('\t%.4f',[1 size(timeSeries,2)-1])]);
	end

% ------------------------------------------------------------------------------
function [Tmed, contanoes, ZtoSpline, already_processed] = ...
		calc_average_old(handles, s, s_flags, Tmed, ZtoSpline, contanoes, this_months, splina, last_processed_month, ...
		nCells, n_anos, already_processed, total_months, months, n_pad_months, do_flags, flag, ...
		fname2, fname3, pintAnoes, track_filled, rows, cols, z_id, z_id_flags, growing_flag, m)
% Chunk of code that calculates the average in the old way and converted to a function.

	counter = 0;
	for (n = this_months)
		counter = counter + 1;
		if (m == 1)						% First year
			Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [n-1 0 0], [1 rows cols]);
		else
			if (splina && n <= last_processed_month)
				already_processed = already_processed + 1;
				offset = 1;
				if (m == 2),	offset = min(1, months(1)-n_pad_months);	end		% Because for the 1st year the series may be shorter
				offset = total_months - (last_processed_month - n) + offset - 1;
				ZtoSpline(:,:,already_processed) = ZtoSpline(:,:,offset);
			else
				Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [n-1 0 0], [1 rows cols]);
				already_processed = 0;
			end
		end

		if (do_flags && ~already_processed)
			Z_flags = nc_funs('varget', fname2, s_flags.Dataset(z_id_flags).Name, [n-1 0 0], [1 rows cols]);
			if (growing_flag),		Z(Z_flags < flag) = NaN;	% Pathfinder style (higher the best) quality flag
			else,					Z(Z_flags > flag) = NaN;	% MODIS style (lower the best) quality flag
			end
		end

		ind = isnan(Z);

		if (pintAnoes && ~already_processed && any(ind(:)))		% If fill spatial holes is requested
			if (track_filled),		ind0 = ind;		end			% Get this Z level original NaNs mask
			Z = inpaint_nans(handles, Z, ind, nCells);			% Select interp method inside inpaint_nans()
			ind = isnan(Z);
			if (track_filled && counter <= 12)					% Write updated quality file (The 'splina' case has counters >> 12)
				mn = (m - 1)*12 + counter - 1;					% This will work only for entire years (not seasons)
				Z_flags(ind0 & ~ind) = flag;					% Promote interpolated pixels to quality 'flag'
				grdutils(Z_flags,'-c');							% Shift by -128 so it goes well with the uint8 add_off elsewere
				if (mn == 0),		nc_io(fname3, sprintf('w%d/time',n_anos*numel(months)), handles, reshape(Z_flags,[1 size(Z_flags)]))
				else,				nc_io(fname3, sprintf('w%d', mn), handles, Z_flags)
				end
			end
		end

		if (~splina)				% Do not interpolate along time (months)
			Z(ind) = 0;				% Transmutate the Anoes
			contanoes = contanoes + ~ind;
			cvlib_mex('add', Tmed, double(Z));
		else						% Pack this year into a 3D temporary variable, to be processed later.
			if (~already_processed),	ZtoSpline(:,:,counter) = Z;		end
		end

		if (n_anos == 1 && numel(this_months) > 12)			% For the secret daily data case
			aguentabar(n / (numel(this_months) + 1)),		drawnow
		end
	end								% End loop over months

% ------------------------1----------2-------3----------4---------5---------6----
function calc_L2_periods(handles, period, tipoStat, regMinMax, grd_out, mask_file)
% Compute averages for 1, 3, 8, month periods of L2 data processed by empilhador
%
% PERIOD	Number of days of the composit period (e.g. 3, 8, 30). A variable number
%			of layers per day is allowed.
%			Alternatively PERIOD can be a vector with the periods to compute, and they don't
%			need to have all the same interval. But ATTENTION that the elements of this
%			vector need to cover the dates in the 3rth dim of the input nc file.
%			This allows selecting a fix start. E.G this will compute 3 days composits of
%			the first 29 days of January 2012.
%				PERIOD=2012.0:3/365:(2012+29/365);
%			Off course dates in input must fall in this period, otherwise output is NaNs only
%
% TIPOSTAT	Variable to control what statistic to compute.
%			0 Compute MEAN of MONTHS period. 1 Compute MINimum, 2 compute MAXimum, 3 compute STD
%
% REGMINMAX	A 2 elements vector with the MIN and MAX values allowed on the Z function (default [0 inf])
%
% GRD_OUT	Name of the netCDF file where to store the result. If not provided, it will be asked here.
%
% MASK_FILE	name of Land mask file. If it has NaNs it will be multiplied but if only 1/0's, 0's will become NaNs

	if (nargin < 5 || isempty(grd_out))			% Note: old and simple CASE 3 in main cannot send here the output name 
		txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
		[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
		if isequal(FileName,0),		return,		end
		grd_out = [PathName FileName];
	end
	if (nargin < 4)
		regionalMIN = 0;	regionalMAX = inf;
	else
		regionalMIN = regMinMax(1);		regionalMAX = regMinMax(2);
		if (regionalMIN == 0 && regionalMAX == 0),	regionalMIN = Inf;	regionalMAX = Inf;	end
	end
	if (nargin < 6),	mask_file = '';		end
	mask = [];			% If needed more than once, this var will hold the masking array

	[z_id, s, rows, cols] = get_ncInfos(handles);

	% Find if we are dealing with a Pathfinder V5.2 daily file
	[is_PFV52, tempos] = PFV52(handles, s);	% if not a PFV5.2, tempos = handles.time(:)

	if (numel(period) == 1)
		periods = fix(tempos(1)) : period : fix(tempos(end))+period-1;	% Don't risk to loose an incomplete last interval
		half_period = repmat(period/2, 1, numel(periods));		% for naming layers in nc file
	else
		periods = period;
		% Need to add 1 because HISTC counts in the interval edge(k) <= x < edge(k+1) and, e.g., day 3.7 is still day 3
		periods(2:end) = periods(2:end) + 1;
		half_period = [diff(periods(:)')/2 (periods(end) - periods(end-1))/2];	% Repeat last value
	end

	N = histc(fix(tempos), periods);
	N(end) = [];		periods(end) = [];		half_period(end) = [];	% Last N is for >= edge(end) and we don't want it
	if (isempty(periods))
		warndlg('There is nothing inside the period(s) you have requested. Bye.', 'Warning'),	return
	end

	handles.was_int16 = false;		% I have to get rid of the need to set this

	aguentabar(0,'title','Computing period means.','CreateCancelBtn');

	% C is the counter to the current layer number being processed.
	c = find(fix(tempos) < periods(1));		% Find the starting layer number
	if (isempty(c)),	c = 1;				% We start at the begining of file.
	else,				c = c(end) + 1;		% We start somewhere at the middle of file.
	end
	
	n_periods = numel(periods);
	for (m = 1:n_periods)
		if (N(m) ~= 0)
% 			Z = alloc_mex(rows, cols, N(m), 'single', NaN);
% 			for (n = 1:N(m))			% Loop over the days in current period
% 				Z(:,:,n) = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [c-1 0 0], [1 rows cols]);
% 				c = c + 1;
% 			end
% 			tmp = doM_or_M_or_M(Z, 1, size(Z,3), regionalMIN, regionalMAX, tipoStat);
% 			clear Z;					% Free memory (need because at the alloc time above, two of them would exist)
			s2 = struct('fname',handles.fname, 'info',s, 'n_layers', N(m), 'rows',rows, 'cols',cols, ...
				'layerOI',c, 'layer_inc',1, 'z_id',z_id);
			[tmp, s2] = doM_or_M_or_M([], 1, 1, N(m), regionalMIN, regionalMAX, tipoStat, s2);
			c = s2.layerOI;				% This is crutial because it tells us where we are in the layer stack
			tmp(tmp == 0) = NaN;		% Reset the NaNs
			if (~isempty(mask_file))
				[tmp, mask] = apply_mask(mask_file, mask, tmp);		% First time reads from MASK_FILE, second on uses MASK
			end
			
			zzz = grdutils(tmp,'-L');
			handles.head(5) = min(handles.head(5), zzz(1));		handles.head(6) = max(handles.head(6), zzz(2));
		else
			tmp = alloc_mex(rows, cols, 1, 'single', NaN);
		end

		% Compute the mean time of this bin and use it to name the layer
		thisLevel = periods(m) + half_period(m);

		% Write this layer to file, but must treate compiled version differently since
		% it is not able to write UNLIMITED files
		if (true || ~handles.IamCompiled)		% TEMP. Later, if it works, we'll simply delete the other branch
			if (m == 1),	nc_io(grd_out, sprintf('w-%f/time',thisLevel), handles, reshape(tmp,[1 size(tmp)]))
			else,			nc_io(grd_out, sprintf('w%d\\%f', m-1, thisLevel), handles, tmp)
			end
		else
			if (m == 1)
				handles.levelVec = periods + half_period;
     			nc_io(grd_out,sprintf('w%d/time',n_periods), handles, reshape(tmp,[1 size(tmp)]))
			else
				nc_io(grd_out, sprintf('w%d', m-1), handles, tmp)
			end
		end

		h = aguentabar(m/n_periods,'title','Computing period means.');	drawnow
		if (isnan(h)),	break,	end
	end
	
% --------------------------------------------------------------------------------------
function [Z, mask] = apply_mask(mask_file, mask, Z)
% When MASK is empty, read data from MASK_FILE otherwise just use MASK to mask out Z.
% Masking is done either by multiplication by MASK or by logical op when MASK is a logical (or only 1 & 0's)
	if (isempty(mask))
		try
			G = gmtmex(['read -Tg ' mask_file]);
			if ~((size(G.z,1) == size(Z,1)) && (size(G.z,2) == size(Z,2)))
				errordlg('The sizes of the Land mask grid and the input array differ. Ignoring masking request.', 'Error')
				return
			end
		catch
			errordlg(sprintf('Error reading file %s\n%s', mask_file, lasterror),'Error')
			return
		end
		if (G.range(5) == 0 && G.range(6) == 1)		% A 1/0's mask
			if (nargout == 2)			% If it's going to be reused, better convert it to logicals right away
				mask = logical(G.z);
				clear G
			end
			Z(mask) = NaN;
		else
			Z = Z .* G.z;
		end
	else
		if (isa(mask, 'logical'))
			Z(mask) = NaN;
		else
			Z = Z .* mask;
		end
	end

% ------------------------------------------------------------------------------------
function Tavg = tideman(handles, Temps, tempos, nPeriods)
% Apply the method suggested by TideMan in this post
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/292502#887770
% but use weeks instead of months which would introduce a significant phase
% offset because a month is too long period for this method

	if (nargin == 3),	nPeriods = 10;	end
	rPeriod = 1 / nPeriods;
	
	decTime = tempos - fix(tempos);				% The decimal part only
	Tavg = zeros(size(Temps,1), size(Temps,2), nPeriods);
	for (k = 1:nPeriods)
		s = (k - 1) * rPeriod;		e = k * rPeriod;
		ind = (decTime >= s & decTime < e);
		thisStackPeriod = Temps(:,:,ind);		% All layers that fall inside this (stacking) period
		Tavg(:,:,k) = doM_or_M_or_M(thisStackPeriod);
	end

% ------------------------------------------------------------------------------------
function T_measured = remove_seazon(handles, T_measured, Tavg, tempos)
% T_MEASURED is the array with measured temperatures (or something else)
% TAVG is the mean cycle computed by the tideman function
% Subtract both (after reinterpolations) and return result in inplace modified T_MEASURED

	nYears = numel(find(diff(fix(tempos)) ~= 0)) + 1;
	nPeriods = size(Tavg,3);			% Number of intervals in which the seazon model has been computed
	rPeriod = 1 / nPeriods;
	ty = linspace(0+rPeriod/2, 1-rPeriod/2, nPeriods)';	% Decimal time centered in the middle of each interval
	t = zeros(nPeriods * nYears, 1);
	for (k = 1:nYears)
		t((k-1)*nPeriods+1 : k*nPeriods) = fix(tempos(1)) + (k-1) + ty;
	end

	y = zeros(1, numel(tempos));
	for (m = 1:size(T_measured,1))
		for (n = 1:size(T_measured,2))
			T_seazon = double(repmat(squeeze(Tavg(m,n,:)), nYears, 1));		% Replicate seazonal cycle
			tmp = double(squeeze(T_measured(m,n,:)));
			indNan = isnan(tmp);
			%y = interp1(t, T_seazon, tempos, 'linear', 'extrap');
			akimaspline(t, T_seazon, tempos, y);		% we don't need a spline but it's much faster
			y(indNan) = NaN;
			T_measured(m,n,:) = single(tmp - y(:));
		end
	end

% ----------------------------------------------------------------------
function [out, s] = get_layer(Z, layer, s)
% Get the layer in one of the following two instances
%	1 - S is empty and Z is [m n p]
%	2 - S is a structure with info on how to read the layer directly from the netCDF file

	if (isempty(s))
		out = Z(:,:,layer);
	else
		out = nc_funs('varget', s.fname, s.info.Dataset(s.z_id).Name, [s.layerOI-1 0 0], [1 s.rows s.cols]);
		s.layerOI = s.layerOI + s.layer_inc;
	end

% ----------------------------------------------------------------------
function [out, s] = doM_or_M_or_M(Z, first_level, lev_inc, last_level, regionalMIN, regionalMAX, tipo, s)
% Compute either the MEAN (TIPO = 0) or the MIN (TIPO = 2), MAX (3) or STD (4) of the period selected
% by the first_level:lev_inc:last_level vector. Normaly a year but can be a season as well.
% NOTE1: This function was only used when SPLINA (see above in calc_yearMean()) up to Mirone 2.2.0
% NOTE2: It is now (2.5.0dev) used again by the tideman function (and other calls)
%
% Because of the potentially very large memory comsumption, we can do the data file reading from
% within this function. In that case Z can be empty and S must be a struct with:
%	S = struct('fname',handles.fname, 'info',s, 'n_layers', N(m), 'rows',rows, 'cols',cols, ...
%	           'layerOI',c, layer_inc,n, 'z_id',z_id);
% whre N_LAYERS is the numbers of layers to be read. 'layerOI' flags the layer to be read and must be
% incremented after each layer reading (which is done by the GET_LAYER() function).
% LAYER_INC is the increment between layers. It will be = 1 for consecutive layers, or 12 for climatologies
% We also return S because it keeps the trace of the current layer (the layerOI field), which is needed
% when there are multiple calls to this function.

	if (nargin == 1)			% Compute average of all layers without any constraint
		first_level = 1;		lev_inc = 1;			last_level = size(Z,3);
		regionalMIN = Inf;		regionalMAX = Inf;		tipo = 0;
	end
	if (nargin < 8),	s = [];	end

	if (tipo == 0)				% Compute the MEAN of the considered period
		% We don't use nanmean_j here because of the regionalMIN|MAX
		[out, s] = get_layer(Z, first_level, s);
		out(out < regionalMIN | out > regionalMAX) = NaN;
		ind = isnan(out);
		contanoes = alloc_mex(size(ind,1), size(ind,2), 'single');
		cvlib_mex('add', contanoes, single(~ind));
		out(ind) = 0;						% Mutate NaNs to 0 so that they don't screw the adition
		for (n = (first_level+lev_inc):lev_inc:last_level)
			[tmp, s] = get_layer(Z, n, s);
			if (~isinf(regionalMIN)),	tmp(tmp < regionalMIN) = NaN;	end
			if (~isinf(regionalMAX)),	tmp(tmp > regionalMAX) = NaN;	end
			ind = isnan(tmp);
			tmp(ind) = 0;
			cvlib_mex('add', contanoes, single(~ind));
			cvlib_mex('add', out, tmp);
		end
		cvlib_mex('div', out, contanoes);			% The mean
	elseif (tipo == 1)						% Median
		% too complicated if the entire dataset is not on memory 
	elseif (tipo == 2 || tipo == 3)			% ...
		v = version;
		if (str2double(v(1)) > 6)			% With R14 and above we use the built in min and max
			if (tipo == 2),		fh = @min;	% Minimum of the selected period
			else,				fh = @max;	% Maximum of the selected period
			end
		else								% But for R13 and compiled we must avoid the BUGGY NaNs comparisons
			if (tipo == 2),		fh = @min_nan;		test = 1e10;
			else,				fh = @max_nan;		test = -1e10;
			end
		end
		if (isempty(s))
			out = feval(fh, Z(:,:,first_level:lev_inc:last_level),[],3);
		else
			[out, s] = get_layer(Z, first_level, s);
			for (k = first_level+lev_inc:lev_inc:s.n_layers)
				[tmp, s] = get_layer(Z, k, s);
				out = feval(fh, out, tmp);
			end
			if (str2double(v(1)) <= 6)
				out(out == test) = NaN;		% Reset the ests values back to NaNs
			end
		end
	elseif (tipo == 4)			% STD
		out = nanstd_j(Z, first_level, lev_inc, last_level, s);
	end

% ----------------------------------------------------------------------
function out = min_nan(A, B)
% Compute minimum of a 3D array or two 2D arrays that have NaNs.
% This function is used only by R13 and compiled versionn that very bugged in whta concerns NaNs comparisons
% Fot the two arguments case, the caller function is responsible to reset the 1e10 to NaNs again.
% This to allow calls in a loop and avoid intermediate wasting replacements.
	if (nargin == 1)
		A(isnan(A)) = 1e10;
		out = min(A, [], 3);
		out(out == 1e10) = nan;
	else
		A(isnan(A)) = 1e10;		B(isnan(B)) = 1e10;
		out = min(A, B);
	end

% ----------------------------------------------------------------------
function out = max_nan(A, B)
% Compute maximum of a 3D array or two 2D arrays that have NaNs.
% This function is used only by R13 and compiled versionn that very bugged in whta concerns NaNs comparisons
% Fot the two arguments case, the caller function is responsible to reset the -1e10 to NaNs again.
% This to allow calls in a loop and avoid intermediate wasting replacements.
	if (nargin == 1)
		A(isnan(A)) = -1e10;
		out = max(A, [], 3);
		out(out == -1e10) = nan;
	else
		A(isnan(A)) = -1e10;	B(isnan(B)) = -1e10;
		out = max(A, B);
	end

% ----------------------------------------------------------------------
function out = nanmean_j(Z, first_level, lev_inc, last_level, s)
% ...
	if (nargin == 1)
		first_level = 1;	last_level = size(Z,3);		lev_inc = 1;	s = [];
	end
	[out, s] = get_layer(Z, first_level, s);
	ind = isnan(out);
	contanoes = alloc_mex(size(ind,1), size(ind,2), 'single');
	cvlib_mex('add', contanoes, single(~ind));
	out(ind) = 0;						% Mutate NaNs to 0 so that they don't screw the adition
	for (n = (first_level+lev_inc):lev_inc:last_level)
		[tmp, s] = get_layer(Z, n, s);
		ind = isnan(tmp);
		tmp(ind) = 0;
		cvlib_mex('add', contanoes, single(~ind));
		cvlib_mex('add', out, tmp);
	end
	cvlib_mex('div', out, contanoes);			% The mean

% ----------------------------------------------------------------------
function out = nanstd_j(Z, first_level, lev_inc, last_level, s)
% Compute the STD taking into account the presence of NaNs
% This is a bit more convoluted for memory efficiency concearns (somethig TMW does not care)
%
% S is a structure with info to read the layers from file instead of relying in Z (that hould be [] than)
% For further info, see help section of the doM_or_M_or_M() function

	if (nargin == 1)
		first_level = 1;	last_level = size(Z,3);		lev_inc = 1;	s = [];
	end

	this_mean = nanmean_j(Z, first_level, lev_inc, last_level, s);
	if (isa(this_mean, 'single'))		% Need this gimnastic because cvlib_mex screws if types are different
		out = alloc_mex(size(this_mean,1), size(this_mean,2), 'single');
	else
		out = alloc_mex(size(this_mean,1), size(this_mean,2), 'double');
	end
	denom = zeros(size(this_mean));			% Swallow the thing but this one has to be done with doubles
	for (n = first_level:lev_inc:last_level)
		[t, s] = get_layer(Z, n, s);
		denom = denom + ~isnan(t);
		cvlib_mex('sub', t, this_mean);
		cvlib_mex('mul', t, t);				% The squares of (xi - xm)
		t(isnan(t)) = 0;
		cvlib_mex('add', out, t)
	end

	denom = max(denom-1, 1);		% divide by (n-1). But when n == 0 or 1, we'll return ones
	denom(denom == 0) = NaN;		% When all NaNs return NaN, and thus avoid a divide by 0
	if (isa(this_mean, 'single')),	denom = single(denom);	end

	cvlib_mex('div', out, denom)	
	cvlib_mex('pow', out, 0.5)

% ----------------------------------------------------------------------
function calc_corrcoef(handles, secondArray, sub_set, pValue, grd_out)
% Compute the correlation coefficien between loaded array and 'secondArray'
% 
% NOTE: THIS IS A HIGHLY MEMORY CONSUMPTION ROUTINE AS ALL DATA IS LOADED IN MEMORY
%
% SECONDARRAY	name of the other netCDF file whose correlation with loaded array will be estimated.
% 			Obviously this file must be of the same size as the series under analysis.
%
% OPTIONS:
% SUB_SET -> A two columns row vec with number of the offset of years where analysis start and stop.
%			For example [3 1] Starts analysis on forth year and stops on the before last year.
%			[0 0] Means using the all dataset.
%
% PVALUE	Logical that if true instruct to compute also the p-value of the estimate
%
% GRD_OUT	Name of the netCDF file where to store the result. If not provided, open Mirone Fig.

	n_anos = handles.number_of_timesteps;
	[z_idA, sA, rows, cols] = get_ncInfos(handles);
	if (nargin < 5),	grd_out = [];	end

	if (nargin >= 3 && (numel(sub_set) == 2))
		jump_anos = sub_set(1);		stop_before_end_anos = sub_set(2);
	else
		jump_anos = 0;				stop_before_end_anos = 0;
	end

	n_anos = n_anos - (jump_anos + stop_before_end_anos);	% Number of layers to be used in this run

	[sB, z_idB, msg] = checkFlags_compat(secondArray, handles.number_of_timesteps, rows, cols);
	if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
	
	arrayB = alloc_mex(rows, cols, n_anos, 'single');
	arrayA = alloc_mex(rows, cols, n_anos, 'single');
	for (m = 1:n_anos)
 		arrayA(:,:,m) = nc_funs('varget', handles.fname, sA.Dataset(z_idA).Name, [(m - 1 + jump_anos) 0 0], [1 rows cols]);
		arrayB(:,:,m) = nc_funs('varget', secondArray, sB.Dataset(z_idB).Name, [(m - 1 + jump_anos) 0 0], [1 rows cols]);
	end

	h = aguentabar(0,'title','Computing correlation coefficients','CreateCancelBtn');

	try
	Rgrid = alloc_mex(rows, cols, 'single', NaN);	% Works on R13 as well
	if (pValue),	Pgrid = alloc_mex(rows, cols, 'single', NaN);	end
	for (m = 1:rows)
		for (n = 1:cols)
			Var1 = double(squeeze(arrayA(m,n,:)));
			Var2 = double(squeeze(arrayB(m,n,:)));
			
			ind = isnan(Var1);
			Var1(ind) = [];		Var2(ind) = [];
			ind = isnan(Var2);
			Var1(ind) = [];		Var2(ind) = [];
			% Completely ad-hoc test (it also jumps land cells)
			if ((numel(Var1) < n_anos*0.66)),		continue,	end

%  			%[R_, p] = corrcoef([Var1 Var2]);	%
			thisN = numel(Var1);
			VarAB = double([Var1 Var2]);
			VarAB = VarAB - repmat(sum(VarAB,1) / thisN, thisN, 1);	% Remove mean
			R = (VarAB' * VarAB) / (thisN - 1);		% Covariance matrix
			d = sqrt(diag(R));		% sqrt first to avoid under/overflow
			R = R ./ (d*d');		% Correlation matrix

			Rgrid(m,n) = single(R(2));

			if (pValue)
	 			p = trend1d_m(VarAB,'-L','-N2r','-R','-P');
				if (p(1) < -0.5 || p(1) > 1),	continue,	end		% Another ad-hoc (CLIPPING)
				Pgrid(m,n) = p(4);

				% To compute p, the Matlab way (but too recent Matlab)
				% Tstat = +/-Inf and p = 0 if abs(r) == 1, NaN if r == NaN.
% 				Tstat = R(2) .* sqrt((thisN-2) ./ (1 - R(2).^2));
% 				p = zeros(2);
% 				p(2) = 2*tpvalue(-abs(Tstat), thisN-2);
% 				p = p + p' + diag(diag(R)); % Preserve NaNs on diag.
% 				Pgrid(m,n) = single(p(2));
			end

		end
		if (rem(m, 10) == 0),	h = aguentabar(m/rows);		end
		if (isnan(h)),	break,	end
	end
	if (isnan(h)),	return,		end
	if (ishandle(h)),	aguentabar(1),		end		% Make sure it goes away.

catch
	disp(lasterror)
end
	
	clear arrayA arrayB

	zz = grdutils(Rgrid,'-L');  handles.head(5:6) = [zz(1) zz(2)];
	tmp.head = handles.head;
	if (isempty(grd_out))	% Show result in a Mirone figure
		tmp.X = linspace(tmp.head(1),tmp.head(2),cols);
		tmp.Y = linspace(tmp.head(3),tmp.head(4),rows);
		tmp.name = 'Correlation';
		mirone(Rgrid, tmp)
	else					% Got output name from input arg
		handles.was_int16 = false;
		handles.computed_grid = true;
		handles.geog = 1;
		nc_io(grd_out, 'w', handles, Rgrid)
	end

	if (pValue)
		zz = grdutils(Pgrid,'-L');  handles.head(5:6) = [zz(1) zz(2)];
		tmp.head = handles.head;
		tmp.X = linspace(tmp.head(1),tmp.head(2),cols);
		tmp.Y = linspace(tmp.head(3),tmp.head(4),rows);
		tmp.name = 'p-values';
		mirone(Pgrid, tmp)
	end

% ----------------------------------------------------------------------
% function p = tpvalue(x,v)
% %TPVALUE Compute p-value for t statistic.

% normcutoff = 1e7;
% if length(x)~=1 && length(v)==1
%    v = repmat(v,size(x));
% end
% 
% % Initialize P.
% p = NaN(size(x));
% nans = (isnan(x) | ~(0<v)); % v == NaN ==> (0<v) == false
% 
% % First compute F(-|x|).
% %
% % Cauchy distribution.  See Devroye pages 29 and 450.
% cauchy = (v == 1);
% p(cauchy) = .5 + atan(x(cauchy))/pi;
% 
% % Normal Approximation.
% normal = (v > normcutoff);
% p(normal) = 0.5 * erfc(-x(normal) ./ sqrt(2));
% 
% % See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1.
% gen = ~(cauchy | normal | nans);
% p(gen) = betainc(v(gen) ./ (v(gen) + x(gen).^2), v(gen)/2, 0.5)/2;
% 
% % Adjust for x>0.  Right now p<0.5, so this is numerically safe.
% reflect = gen & (x > 0);
% p(reflect) = 1 - p(reflect);
% 
% % Make the result exact for the median.
% p(x == 0 & ~nans) = 0.5;

% ----------------------------------------------------------------------
function [tSeries, indTSCurr] = getTimeSeries(ZtoSpline, tSeries, indTS, indTSCurr, orig, first_month, last_month)
% Fill the TSERIES array with the original as well as the spline interpolated time series.
% This is intended mostly to provide a way to check the goodness (or not) of the interp mechanism
%
% ZtoSpline The 3D array with yearly (+- pad months) time series -- before (orig = true) and after splining
% tSeries	a Mx(2*n+1) array to store the original data (even columns) and the spline interpolated (odd coluns)
% indTS		array indices of where to extract the time series (a Mx2 array)
% indTSCurr Start index to write on the tSeries vector. Needs to be incremented by 12 at the end of this function
% orig		Logical meaning that if true we are writting original series or, otherwise interpolated, where NANs, values

	if (orig),		iStart = 2;			% Remember that first column is always the time
	else,			iStart = 3;
	end
	for (k = 1:size(indTS,1))	% Loop over number of check points
		tSeries(indTSCurr:(indTSCurr+11), iStart + 2*(k-1)) = ZtoSpline(indTS(k,2), indTS(k,1), first_month:last_month);
	end
	indTSCurr = indTSCurr + 12;			% Remember, this works only for yearly means from month data

% ----------------------------------------------------------------------
function pass_by_count(handles, count, fname2)
% Check the curently active 3D file against a count file
%
% COUNT		Threshold count value. Nodes on in-memory file that have a count on the
% 			corresponding node of FNAME2 < COUNT are set to NaN
%
% OPTIONS:
% FNAME2 	name of a netCDF file with the count quality flags. Asked if not provided
%			Obviously this file must be of the same size as the series under analysis.

	txt1 = 'netCDF grid format (*.nc,*.grd)';		txt2 = 'Select output netCDF grid';
	[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
	if isequal(FileName,0),		return,		end
	grd_out = [PathName FileName];
	
	if (nargin == 2)		% No count grid transmitted. Ask for it
		[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},'Select input netCDF file','get');
		if isequal(FileName,0),		return,		end
		fname2 = [PathName FileName];
	else					% Got a name. Check that it exists
		if (exist(fname2,'file') ~= 2)
			errordlg(['Blheak!! ' fname2 ' does not exist (even if you think so). Bye Bye'],'Error'),	return
		end
	end

	z_id = handles.netcdf_z_id;
	s = handles.nc_info;				% Retrieve the .nc info struct
	rows = s.Dataset(z_id).Size(end-1);
	cols = s.Dataset(z_id).Size(end);

	s_flags = nc_funs('info',fname2);
	[X,Y,Z,head,misc] = nc_io(fname2,'R');
	z_id_flags = misc.z_id;
	if ~(numel(head) == 9 && isfield(misc,'z_id'))
		errordlg(['Blheak!! ' fname2 ' is is not a file with presumably with a count of quality flags. By'],'Error'),	return
	end
	if (misc.z_dim(1) < handles.number_of_timesteps)
		errordlg('Buhhuu!! The count flags file has less "planes" than the-to-be-counted-file. By','Error'),	return
	end
	if (~isequal([rows cols], [s_flags.Dataset(z_id_flags).Size(end-1) s_flags.Dataset(z_id_flags).Size(end)]))
		errordlg('Buhhuu!! quality flags and the-to-be-counted-file have not the same size. By','Error'),		return
	end

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;
	
	aguentabar(0,'title',['NaNify countings < ' sprintf('%d',count)],'CreateCancelBtn')

	n_layers = handles.number_of_timesteps;
	for (m = 1:n_layers)
		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m-1 0 0], [1 rows cols]);
		Z_flags = nc_funs('varget', fname2, s_flags.Dataset(z_id_flags).Name, [m-1 0 0], [1 rows cols]);
		Z(Z_flags < count) = NaN;

		if (m == 1),		nc_io(grd_out, sprintf('w%d/time',n_layers), handles, reshape(Z,[1 size(Z)]))
		else,				nc_io(grd_out, sprintf('w%d', m-1), handles, Z)
		end

		h = aguentabar(m/n_layers);
		if (isnan(h)),	break,	end
	end

% ------------------------------------------------------------------------------------------------
function [is_PFV52, times, anos] = PFV52(handles, s)
% Find if we are dealing with one Pathfinder V5.2 daily file.
% If yes and user request it, convert and return also the time vector in the form of decimal years
% In case file is not a PFV5.5, TIMES = HANDLES.TIME

	is_PFV52 = false;
	for (k = 1:numel(s.Attribute))
		if (strcmp(s.Attribute(k).Name, 'description') && ~isempty(strfind(s.Attribute(k).Value,'Pathfinder 5.2 daily')))
			is_PFV52 = true;
			break
		end
	end

	if (nargout > 1)		% This relies on the fact that the apropriate 'description' global attribute has been set
		anos = fix(handles.time(:));
		if (is_PFV52)
			n_days_in_year = datenummx(anos+1,1,1) - datenummx(anos,1,1);
			times = (handles.time - anos) * 1000;		% Decimal day of the year (wrongly aka julian day)
			times = anos + times ./ n_days_in_year;		% Now we have decimal years
		else
			times = handles.time(:);
		end
	end

% ----------------------------------------------------------------------
function Z = inpaint_nans(handles, Z, bw, nCells)
% Interpolate holes in Z that are smaller than NCELLS in size
%
% BW is a logicall array which maps where Z has NaNs
	
	if (nargin == 3),	nCells = 100;	end

	use_surface = true;
	use_bicubic = false;
	pad = 4;				% Number of cells to increase the rectangle that encloses each hole
	head = handles.head;
	[rows, cols] = size(Z);
	
	% Retain only <= handles.nCells sized of connected groups
	bw2 = img_fun('bwareaopen', bw, nCells);
	bw = xor(bw, bw2);
	clear bw2;
	
	B = img_fun('find_holes',bw);

	opt_I = ' ';
	if (use_surface),		opt_I = sprintf('-I%.10f/%.10f',head(8),head(9));	end

	n_buracos = numel(B);
	for (i = 1:n_buracos)
		% Get rectangles arround each hole
		x_min = min(B{i}(:,2));			x_max = max(B{i}(:,2));
		y_min = min(B{i}(:,1));			y_max = max(B{i}(:,1));
		x_min = max(1,x_min-pad);		x_max = min(x_max+pad,cols);
		y_min = max(1,y_min-pad);		y_max = min(y_max+pad,rows);
		x_min = head(1) + (x_min-1)*head(8);    x_max = head(1) + (x_max-1)*head(8);
		y_min = head(3) + (y_min-1)*head(9);    y_max = head(3) + (y_max-1)*head(9);

		rect_crop = [x_min y_min (x_max-x_min) (y_max-y_min)];
		[Z_rect, r_c] = cropimg(head(1:2),head(3:4),Z,rect_crop,'out_grid');
		[bw_rect, l.lixo] = cropimg(head(1:2),head(3:4),bw,rect_crop,'out_grid');
		Z_rect = double(Z_rect);      % It has to be (GHRRRRRRRRRRRRR)

		%X = x_min:head(8):x_max;	Y = y_min:head(9):y_max;
		X = linspace(x_min, x_max, size(Z_rect, 2));		% Safer against round off errors 
		Y = linspace(y_min, y_max, size(Z_rect, 1));
		[XX,YY] = meshgrid(X,Y);
		XX(bw_rect) = [];			YY(bw_rect) = [];		Z_rect(bw_rect) = [];

		if (use_surface)
			opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', X(1), X(end), Y(1), Y(end));
			%Z_rect = c_surface(XX(:), YY(:), Z_rect(:), opt_R, opt_I, '-T.25');
			Z_rect = gmtmbgrid_m(XX(:), YY(:), Z_rect(:), opt_R, opt_I, '-T.25', '-Mz');
		elseif (use_bicubic)
			Z_rect = griddata_j(XX(:), YY(:), Z_rect(:), X, Y', 'cubic');
		else
			Z_rect = griddata_j(XX(:), YY(:), Z_rect(:), X, Y', 'linear');
		end

		% Inprint the processed rectangle back into orig array
		if (isa(Z,'single')),		Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
		elseif (isa(Z,'int16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = int16(Z_rect);
		elseif (isa(Z,'uint16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = uint16(Z_rect);
		else						Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
		end
	end

% ------------------------1--------2-------3------4----------5--------6---------8----
function calc_polygAVG(handles, fnameOut, op, fnamePolys, sub_set, fnameFlag, quality)
% This function search for polygons (patches or closed lines) and computes averages
% of whatever quantity is respresented inside those polygones. The result is saved
% in an ASCII file whose first two columns contain the the polygons (x,y) centroid.
%
% OPTIONS:
%
% FNAMEOUT		Name of the ouput file. If not provided, it will be asked for here.
%
% OP			Operation to apply to the data inside each polygon. If not provided or is [],
%				defaults to 'mean'. That is compute the mean of all values inside polygon.
%				Otherwise it must be the name of a Matlab function that can executed via
%				the function handles mechanism. For example 'median'.
%
% FNAMEPOLYS	A file name of a (x,y) polygon or the name of a list of polygons (one per line).
%				If given no attempt is made to fish the polygons from line handles.
%
% SUB_SET		A two columns row vec with number of the offset of years where analysis start and stop.
%				For example [3 1] Starts analysis on forth year and stops on the before last year.
%				[0 0] Means using the all dataset.
%
% FNAMEFLAG		name of a netCDF file with quality flags. Obviously this file must be of
%				the same size as the series under analysis. If not provided no quality check is done.
%
% QUALITY		Threshold quality value. Only values of quality >= FLAG will be taken into account
%				NOTE: For MODIS use negative FLAG. Than, values are retained if quality <= abs(FLAG)

	% -------------------------------- Options parsing -------------------------------------
	if (nargin < 4)
		fnamePolys = [];		sub_set = [0 0];	fnameFlag = [];
	end
	if (nargin == 1)
		fnameOut = [];			fhandle = @local_avg;
	elseif (nargin == 2)
		fhandle = @local_avg;	fnamePolys = [];
	elseif (nargin >= 3)
		if (ischar(op)),		fhandle = str2func(op);
		else,					fhandle = @local_avg;
		end
		if (nargin == 4)
			sub_set = [0 0];	fnameFlag = [];
		elseif (nargin == 5)
			fnameFlag = [];
		end
	end
	% -------------------------------- END of options parsing -----------------------------------

	if (numel(sub_set) == 2)
		jump_start = sub_set(1);		stop_before_end = sub_set(2);
	else
		jump_start = 0;					stop_before_end = 0;
	end

	[z_id, s, rows, cols] = get_ncInfos(handles);

	%------------- Check for quality flags request -------------------
	if (~isempty(fnameFlag))
		[s_flags, z_id_flags, msg] = checkFlags_compat(fnameFlag, handles.number_of_timesteps, rows, cols);
		if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
		do_flags = true;
		if (quality > 0 ),	growing_flag = true;		% PATHFINDER flags
		else,				growing_flag = false;		% MODIS flags
		end
	else
		do_flags = false;
	end

	% --------------- Fish the polygon handles from figure (if not provided) ------------
	if (isempty(fnamePolys))

		handMir = guidata(handles.hMirFig);
		hLine = findobj(handMir.axes1,'Type','line');
		hLine = [hLine; findobj(handMir.axes1,'Type','patch')];
		if (isempty(hLine))
			errordlg('polygAVG: No polygon file provided and no polygon in fig either. Bye Bye.','Error')
			return
		end

		polys = cell(1, numel(hLine));
		N = 1;
		for (k = 1:numel(hLine))
			x = get(hLine(k),'XData');   y = get(hLine(k),'YData');
			if (numel(x) >= 3 && x(1) == x(end) && y(1) == y(end) )
				polys{N} = [x(:) y(:)];
				N = N + 1;
			end
		end
		N = N - 1;					% There was one too much incement above
		polys(N+1:end) = [];		% Remove unused

	else							% Got a filename or a list of polygons files in input
		[bin, n_column, multi_seg, n_headers] = guess_file(fnamePolys);
		if (isempty(bin))
			errordlg(['Error reading file (probaby empty)' fnamePolys],'Error'),	return
		end
		if (n_column == 1 && multi_seg == 0)			% Take it as a file names list
			fid = fopen(fnamePolys);
			c = fread(fid,'*char')';	fclose(fid);
			names = strread(c,'%s','delimiter','\n');   clear c fid;
		else
			names = {fnamePolys};
		end

		for (k = 1:numel(names))
			fname = names{k};
			if (isempty(n_headers)),    n_headers = NaN;    end
			if (multi_seg)
				polys = text_read(fname,NaN,n_headers,'>');
			else
				polys = {text_read(fname,NaN,n_headers)};
			end
		end

		N = numel(polys);
	end

	if (N == 0)
		errordlg('Fiu Fiu! No closed polygons to compute whaterver average value inside. Bye Bye.','Error')
		return
	end
	% --------------------------- END of polygons fishing section -----------------------------------

	nLayers = handles.number_of_timesteps - (jump_start + stop_before_end);	% Number of layers to be used in this run
	series_vec = (jump_start:(nLayers - 1 + jump_start)) + 1;		% Add 1 so it never starts at 0 (no good for indices)
	avg = zeros(nLayers,N) * NaN;
	THRESH = 0.25;						% Minimum percentage of valid points inside poly

	aguentabar(0,'title','Compute poligonal averages','CreateCancelBtn')

	for (m = series_vec)				% Loop over layers ensemble
		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m-1 0 0], [1 rows cols]);

		if (do_flags)
			flags = nc_funs('varget', fnameFlag, s_flags.Dataset(z_id_flags).Name, [m-1 0 0], [1 rows cols]);
			if (growing_flag),		Z(flags < quality) = NaN;	% Pathfinder style (higher the best) quality flag
			else,					Z(flags > quality) = NaN;	% MODIS style (lower the best) quality flag
			end
		end

		for (k = 1:N)				% Loop over all polygons
			x = polys{k}(:,1);			y = polys{k}(:,2);
			xp(1) = min(x);				xp(2) = max(x);
			yp(1) = min(y);				yp(2) = max(y);
			rect_crop = [xp(1) yp(1) (xp(2) - xp(1)) (yp(2) - yp(1))];
			x_lim = [xp(1) xp(2)];		y_lim = [yp(1) yp(2)];

			% Get a BoundingBox rect to save work in mask computing
			[Z_rect, l.lixo] = cropimg(handles.head(1:2),handles.head(3:4),Z,rect_crop,'out_grid');
			mask = img_fun('roipoly_j',x_lim,y_lim,Z_rect,x,y);

			% Test for a minimum of valid elements inside polygon
			zz = Z_rect(mask);
			zz = zz(:);
			ind = isnan(zz);
			if (~any(ind))
				%avg(m,k) = sum(double(zz)) / numel(zz);
				avg(m,k) = feval(fhandle, zz);
			else			% Accept/Reject based on % of valid numbers
				nAnoes = sum(ind);		nInPoly = numel(zz);
				if (nAnoes / nInPoly < THRESH)
					zz = zz(~ind);
					%avg(m,k) = sum(double(zz)) / numel(zz);
					avg(m,k) = feval(fhandle, zz);
				end
			end
		end
		h = aguentabar(m/nLayers);
		if (isnan(h)),	break,	end

	end

	if (isnan(h)),	return,		end
	
% 	% -------------- We still need to determine polygon's area to compute the areal average
% 	for (k = 1:N)
% 		x = polys{k}(:,1);		y = polys{k}(:,2);
% 		if (handles.geog)
% 			area = area_geo(y,x);    % Area is reported on the unit sphere
% 			area = area * 4 * pi * (6371005^2);
% 		else
% 			area = polyarea(x,y);   % Area is reported in map user unites
% 		end
% 		avg(:,k) = avg(:,k) / (area * 1e-6);	% per km^2
% 	end
	
	% --------------- Now finaly save the result in a file	------------------
	if (isempty(fnameOut))
		[FileName,PathName] = put_or_get_file(...
			handles,{'*.dat;*.DAT','ASCII file'; '*.*', 'All Files (*.*)'},'Output file','put');
		if isequal(FileName,0),		return,		end
		[PATH,FNAME,EXT] = fileparts([PathName FileName]);
		if isempty(EXT),	fname = [PathName FNAME '.dat'];
		else,				fname = [PathName FNAME EXT];
		end
	else
		fname = fnameOut;
	end
	
	% Open and write to ASCII file
	if (ispc),		fid = fopen(fname,'wt');
	elseif (isunix),fid = fopen(fname,'w');
	else,			error('aquamoto: Unknown platform.');
	end

	% Calculate a a rough polygon centroid
	centro = zeros(N,2);
	for (k = 1:N)				% Loop over polygons
		centro(k,:) = mean(polys{k});
	end
	
	fprintf(fid, ['#  \t', repmat('%g(X)\t', [1,N]) '\n'], centro(:,1));
	fprintf(fid, ['# T\t', repmat('%g(Y)\t', [1,N]) '\n'], centro(:,2));

	try			t = handles.time;		t = t(:);		% Layers's times
	catch,		t = (1:size(avg,1))';
	end

	fprintf(fid,['%.2f\t' repmat('%f\t',[1,N]) '\n'], [t avg]');
	fclose(fid);

% ----------------------------------------------------------------------
function out = local_avg(x)
% A simple average function, but that convert X do doubles and will be assigned a function handle
	out = sum(double(x)) / numel(x);

% ----------------------------------------------------------------------
function calc_flagsStats(handles, months, flag, opt)
% Compute per/pixel annual counts of pixel values with a quality >= flag
% Perfect locations will have a count of 12. Completely cloudy => count = 0.
%
% MONTHS 	is a vector with the months uppon which the mean is to be computed
%		example: 	months = 1:12		==> Computes yearly mean
%					months = 6:8		==> Computes June-July-August seazonal means
% FLAG		Threshold quality value. Only values of quality >= FLAG will be taken into account
% OPT		== 'per_year' output a file with N years planes (count per year)
%			otherwise ouputs a file with N months planes (each plane has a mounthly count)

	if (nargin < 3),		flag = 7;	opt = 'per_year';		end			% If not provided, defaults to best quality
	if (nargin < 4),		opt = 'per_year';		end

	txt1 = 'netCDF grid format (*.nc,*.grd)';
	txt2 = 'Select output netCDF grid';
	[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
	if isequal(FileName,0),		return,		end
	grd_out = [PathName FileName];

	s = handles.nc_info;				% Retrieve the .nc info struct
	rows = s.Dataset(handles.netcdf_z_id).Size(end-1);
	cols = s.Dataset(handles.netcdf_z_id).Size(end);

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;

	n_anos = handles.number_of_timesteps / 12;

	aguentabar(0,'title','Compute flag quality counts','CreateCancelBtn')
	goodCount = zeros([rows, cols]);			% 

	if (strcmp(opt, 'per_year'))		% Make the counting on a per year basis
		for (m = 1:n_anos)
			contanoes = zeros(rows, cols);
			for (n = months)
				mn = (m - 1)*12 + n - 1;
				Z = nc_funs('varget', handles.fname, s.Dataset(handles.netcdf_z_id).Name, [mn 0 0], [1 rows cols]);
				ind = Z < flag;
				Z(ind) = 0;
				Z(~ind) = 1;				% Bellow threshold quality are set to zero
				contanoes = contanoes + ~ind;
				goodCount(:,:) = goodCount(:,:) + double(Z);
			end
			tmp = int16(goodCount(:,:));
	
			if (m == 1),	nc_io(grd_out,sprintf('w%d/time',n_anos), handles, reshape(tmp,[1 size(tmp)]))
			else,			nc_io(grd_out, sprintf('w%d', m-1), handles, tmp)
			end
			
			h = aguentabar(m/n_anos);
			if (isnan(h)),	break,	end
			goodCount = goodCount * 0;			% Reset it to zeros
		end

	else			% Make the counting on a per month basis
		for (n = months)
			contanoes = zeros(rows, cols);
			for (m = 1:n_anos)
				mn = (m - 1)*12 + n - 1;
				Z = nc_funs('varget', handles.fname, s.Dataset(handles.netcdf_z_id).Name, [mn 0 0], [1 rows cols]);
				ind = Z < flag;
				Z(ind) = 0;
				Z(~ind) = 1;				% Bellow threshold quality are set to zero
				contanoes = contanoes + ~ind;
				goodCount(:,:) = goodCount(:,:) + double(Z);
			end
			tmp = int16(goodCount(:,:));
	
			if (n == 1),	nc_io(grd_out,sprintf('w%d/time',numel(months)), handles, reshape(tmp,[1 size(tmp)]))
			else,			nc_io(grd_out, sprintf('w%d', n-1), handles, tmp)
			end
			
			h = aguentabar(n/numel(months));
			if (isnan(h)),	break,	end	
			goodCount = goodCount * 0;			% Reset it to zeros
		end
	end

% ----------------------------------------------------------------------
function do_math(handles, opt, subSet1, fname2, subSet2)
% Perform some basic algebraic operations on 3D planes
%
% OPT can be one of:
% 'sum'       => add all layers
% 'count'     => count the acumulated number of non-NaNs of all layers
% 'diffstd'   => Compute mean and sdt of A - B (inmemory - fname2)
%                The array read from FNAME2 does not need to be exactly comptible with A.
%                When it is not a call to grdsample will take care of compatibilization.
%
% SUBSET1 & SUBSET2 are (optional) two columns row vec with number of the offset layers
% where analysis starts and stop.
% For example [3 1] Starts analysis on forth layer and stops on the before last year.
% [0 0] Means using the all dataset (the default)
%
%	NOT FINISHED, NEED INPUT PARSING

	grid1 = [];		grid2 = [];
	s1 = handles.nc_info;				% Retrieve the .nc info struct
	rows1 = s1.Dataset(handles.netcdf_z_id).Size(end-1);
	cols1 = s1.Dataset(handles.netcdf_z_id).Size(end);
	nLayers = handles.number_of_timesteps;

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;

	% Find if we are dealing with a Pathfinder V5.2 daily file
	[is_PFV52, tempos] = PFV52(handles, handles.nc_info);		% if not a PFV5.2, tempos = handles.time(:)

	if (strcmp(opt, 'sum'))
		% WARNING: SUBSET1 not implemented
		grid1 = nc_funs('varget', handles.fname, s1.Dataset(handles.netcdf_z_id).Name, [0 0 0], [1 rows1 cols1]);
		is_int8 = isa(grid1, 'int8');		is_uint8 = isa(grid1, 'uint8');
		if (is_int8 || is_uint8),			grid1 = int16(grid1);		end		% To avoid ovelflows
		for (m = 2:nLayers)
			Z = nc_funs('varget', handles.fname, s1.Dataset(handles.netcdf_z_id).Name, [m-1 0 0], [1 rows1 cols1]);
			if (is_int8 || is_uint8),		Z = int16(Z);			end			% To avoid ovelflows
			cvlib_mex('add',grid1, Z);
		end
		figName1 = 'Stack Sum';

	elseif (strcmpi(opt, 'count'))			% Count the acumulated number of non NaNs along the 3rth dim of the 3D array
		
		% Check (GUESS!) if subsets were given as years instead of layer numbers
		if (is_PFV52 && subSet1(1) >= 1982 && subSet1(2) <= 2014)	% Though we only have them till 2011
			ind = find(tempos > subSet1(1));
			subSet1(1) = ind(1);
			ind = find(tempos > subSet1(2));
			subSet1(2) = nLayers - ind(1);
		end
	
		aguentabar(0,'title','Counting...','CreateCancelBtn')
		grid1 = alloc_mex(rows1, cols1, 'single');
		ini = 1+subSet1(1);		fim = nLayers-subSet1(2);	N = fim - ini + 1;	N10 = fix(N / 10);	c = 0;
		for (m = ini:fim)
			Z = nc_funs('varget', handles.fname, s1.Dataset(handles.netcdf_z_id).Name, [m-1 0 0], [1 rows1 cols1]);
			ind = isnan(Z);
			Z(~ind) = 1;	Z(ind) = 0;
			cvlib_mex('add',grid1, Z);
			p = fix(((m - ini) / N10));
			if (p)
				c = c + 1;				p = p + c - 1;		N10 = N10 + fix(N / 10);
				aguentabar(p/10)
			end
		end
		grid1(grid1 == 0) = NaN;
		aguentabar(1)
		figName1 = 'Stack Count';

	elseif (strcmpi(opt, 'diffstd'))
		if (nargin == 4),	subSet2 = [0 0];		end
		s2 = nc_funs('info', fname2);
		[X1,Y1,Z,head1] = nc_io(handles.fname, 'R');
		[X2,Y2,Z,head2,misc2] = nc_io(fname2,'R');
		rows2 = s2.Dataset(misc2.z_id).Size(end-1);
		cols2 = s2.Dataset(misc2.z_id).Size(end);
		nLayers2 = misc2.z_dim(1);
		n1 = max(1,subSet1(1)+1):(nLayers  - subSet1(2));	% Use these steps of file 1
		n2 = max(1,subSet2(1)+1):(nLayers2 - subSet2(2));	% Use these steps of file 2
		nSteps = min(numel(n1), numel(n2));					% Minimum number of steps that satisfies both files
		grid1 = alloc_mex(rows1, cols1, 'single');
		aguentabar(0,'title','Computing means and STDs')
		for (k = 0:nSteps-1)								% OK, now first compute meam of (A - B)
			tmp = nc_funs('varget', handles.fname, s2.Dataset(handles.netcdf_z_id).Name, [subSet1(1)+k 0 0], [1 rows1 cols1]);
			z2 = nc_funs('varget', fname2, s2.Dataset(misc2.z_id).Name, [subSet2(1)+k 0 0], [1 rows2 cols2]);
			z2 = make_compatible(z2, head2, head1);
			cvlib_mex('sub', tmp, z2);
			cvlib_mex('add', grid1, tmp);
			aguentabar((k+1)/(2*nSteps));
		end
		cvlib_mex('CvtScale', grid1, 1/nSteps, 0);			% Divide by N --- MEAN of (A- B)
		grid2 = alloc_mex(rows1, cols1, 'single');
		for (k = 0:nSteps-1)								% And now the STD
			tmp = nc_funs('varget', handles.fname, s1.Dataset(handles.netcdf_z_id).Name, [subSet1(1)+k 0 0], [1 rows1 cols1]);
			z2 = nc_funs('varget', fname2, s2.Dataset(misc2.z_id).Name, [subSet2(1)+k 0 0], [1 rows2 cols2]);
			z2 = make_compatible(z2, head2, head1);
			cvlib_mex('sub', tmp, z2);
			cvlib_mex('sub', tmp, grid1);					% Subtract mean
			cvlib_mex('mul', tmp, tmp);						% Square
			cvlib_mex('add', grid2, tmp);
			aguentabar(k/nSteps + 0.5);
		end
		cvlib_mex('CvtScale', grid2, 1/nSteps, 0);			% Divide by N
		cvlib_mex('pow', grid2, 0.5);						% sqrt --- STD of (A - B)
		aguentabar(1);
		clear tmp z2
		figName1 = 'Mean of differences';
		figName2 = 'STD of differences';
	end

	tmp.head = handles.head;
	if (isa(grid1,'single'))
		zz = grdutils(grid1,'-L');  tmp.head(5:6) = [zz(1) zz(2)];		% Singles & NaNs = BUGs in R13
	else
		tmp.head(5:6) = [double(min(grid1(:))) double(max(grid1(:)))];
	end
	tmp.X = linspace(tmp.head(1),tmp.head(2),cols1);
	tmp.Y = linspace(tmp.head(3),tmp.head(4),rows1);
	tmp.name = figName1;
	mirone(grid1, tmp)
	if (~isempty(grid2))
		tmp.head(5:6) = [double(min(grid2(:))) double(max(grid2(:)))];
		tmp.name = figName2;		mirone(grid2, tmp)
	end

% -----------------------------------------------------------------------------------------
function Z2 = make_compatible(Z2, head2, head1)
% Reinterpolate Z2 grid to be compatible in terms of -R and resolution with HEAD1 params
	opt_I = sprintf('-I%.8f/%.8f', head1(8:9));
	opt_R = sprintf('-R%.8f/%.8f/%.8f/%.8f', head1(1:4));
	Z2 = c_grdsample(Z2, head2, opt_R, opt_I, '-Q', '-Lg');
	
% ------------------------------------------------------------------------------------
function count_blooms(handles, Ncount, sub_set, grd_out)
% Given a 3D array with chlorophyl concentration, count events where there are at least
% NCOUNT consecutive growing elements along the 3rth dimension of the 3D array.
% For example, when cell (i,j,k) grows for 3 consecutive k, we store the value 3 (or higher
% the condition holds for k > 3), otherwise we write 0 at the cell position.
% This is a basic way of detecting blooms
%
% OPTIONS:
% NCOUNT  -> Minimum number of sucessive growing events so that number is stored in output file
%
% SUB_SET -> A two elements row vec with number of the offset of years where analysis start and stop.
%			For example [3 1] Starts analysis on forth year and stops on the before last year.
%			[0 0] Means using the all dataset.
%
% GRD_OUT	Name of the netCDF file where to store the result. Asked here if not provided.
%
% NOTE: Low memory footprint. Read laye, process it and write result

	if (nargin == 1)
		Ncount = 3;
	end
	if (nargin < 4 || (nargin == 4 && isempty(grd_out)))
		txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
		[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
		if isequal(FileName,0),		return,		end
		grd_out = [PathName FileName];
	end

	if (nargin >= 3 && (numel(sub_set) == 2))
		jump_slices = sub_set(1);		stop_before_end_slices = sub_set(2);
	else
		jump_slices = 0;				stop_before_end_slices = 0;
	end

	n_slices = handles.number_of_timesteps;
	[z_id, s, rows, cols] = get_ncInfos(handles);
	n_slices = n_slices - (jump_slices + stop_before_end_slices);	% Number of layers to be used in this run

	aguentabar(0,'title','Counting blooms')
 	layer1 = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [jump_slices 0 0], [1 rows cols]);
 	layer2 = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [1 + jump_slices 0 0], [1 rows cols]);
	tmp = int16(layer2 > layer1);
	layer1 = layer2;

	for (m = 1:n_slices-2)		% Loop over remaining slices (first two were already consumed above)
 		layer2 = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [(m + 1 + jump_slices) 0 0], [1 rows cols]);
		ind = (layer2 > layer1);
		cvlib_mex('addS', tmp, 1)% Increase counter of consecutive growings
		tmp(~ind) = 0;			% Reset to zero those that interrupt the growing escalade
		t = tmp;
		ind = (tmp < Ncount);	% Find those that have a counting lower than the threshold.
		t(ind) = 0;				% Those that didn't (yet?) reach the threshould are written as zero

		thisLevel = handles.time(m+2);
		if (m == 1),	nc_io(grd_out, sprintf('w-%f/time',thisLevel), handles, reshape(t,[1 size(t)]))
		else,			nc_io(grd_out, sprintf('w%d\\%f', m-1, thisLevel), handles, t)
		end
		layer1 = layer2;

		h = aguentabar(m/(n_slices-2));	drawnow
		if (isnan(h)),	break,	end
	end

% -----------------------------------------------------------------------------------------
function fid = write_vtk(handles)
% Write a 3D netCDF into a VTK format

	txt1 = 'VTK format (*.vtk)';
	txt2 = 'Select output VRT file';
	[FileName,PathName] = put_or_get_file(handles,{'*.vtk',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.vtk');
	if isequal(FileName,0),		return,		end
	fname_out = [PathName FileName];

	s = handles.nc_info;				% Retrieve the .nc info struct
	rows = s.Dataset(handles.netcdf_z_id).Size(end-1);
	cols = s.Dataset(handles.netcdf_z_id).Size(end);
	nLayers = handles.number_of_timesteps;

	fid = fopen(fname_out, 'wb','b');
	fprintf(fid, '# vtk DataFile Version 2.0\n');
	fprintf(fid, 'converted from A B\n');
	fprintf(fid, 'BINARY\n');
	fprintf(fid, 'DATASET RECTILINEAR_GRID\n');
	fprintf(fid, 'DIMENSIONS %d %d %d\n', cols, rows, nLayers);
	fprintf(fid, 'X_COORDINATES %d float\n', cols);
	X = linspace(handles.head(1), handles.head(2), cols);
	fwrite(fid, X, 'real*4');
	fprintf(fid, 'Y_COORDINATES %d float\n', rows);
	X = linspace(handles.head(3), handles.head(4), rows);
	fwrite(fid, X, 'real*4');
	fprintf(fid, 'Z_COORDINATES %d float\n', nLayers);
	fwrite(fid, 1:nLayers, 'real*4');
	fprintf(fid, 'POINT_DATA %d\n', cols * rows * nLayers);
	fprintf(fid, 'SCALARS dono float 1\n');
	fprintf(fid, 'LOOKUP_TABLE default\n');

	for (m = 1:nLayers)
		Z = nc_funs('varget', handles.fname, s.Dataset(handles.netcdf_z_id).Name, [m-1 0 0], [1 rows cols]);
		Z = double(Z');
		fwrite(fid, Z(:), 'real*4');
	end
		
% ----------------------------------------------------------------------
function [z_id, s, rows, cols] = get_ncInfos(handles)
% Since this is done in nearly all functions, centralize it here
	z_id = handles.netcdf_z_id;
	s = handles.nc_info;			% Retrieve the .nc info struct
	rows = s.Dataset(z_id).Size(end-1);
	cols = s.Dataset(z_id).Size(end);
		
% ----------------------------------------------------------------------
function [s_flags, z_id_flags, msg] = checkFlags_compat(fname, number_of_timesteps, rows, cols)
% Check that the quality flags file FNAME is compatible with the other 3 input params
	msg = [];
	s_flags = nc_funs('info',fname);
	[X,Y,Z,head,misc] = nc_io(fname,'R');
	z_id_flags = misc.z_id;
	if ~(numel(head) == 9 && isfield(misc,'z_id'))
		msg = ['Blheak!! ' fname ' is is not a file with presumably with quality flags. By'];	return
	end
	if (numel(misc.z_dim) <= 2)
		msg = ['Ghrrr!! The ' fname ' is is not a 3D file. By'];		return
	end
	if (misc.z_dim(1) < number_of_timesteps)
		msg = 'Buhhuu!! The quality flags file has less "planes" than the-to-be-flagged-file. By';	return
	end
	if (~isequal([rows cols], [s_flags.Dataset(z_id_flags).Size(end-1) s_flags.Dataset(z_id_flags).Size(end)]))
		msg = 'Buhhuu!! quality flags and the-to-be-flagged-file have not the same size. By';
	end

% -------------------------------------------------------------------------------
function out = script_control(handles, auto)
% See if the OPTcontrol.txt file has an entry pointing to a file with parameters
% to run the aquaPlugin. If it has, read and parse that file.
%
% AUTO	If it's a char, than interpret it as the control script name directly
%		Any other type tells the program to search the name in OPTcontrol.txt

	if (~ischar(auto))
		opt_file = [handles.home_dir filesep 'data' filesep 'OPTcontrol.txt'];
		out = [];
		if (exist(opt_file, 'file') ~= 2),	return,		end

		fid = fopen(opt_file, 'r');
		c = (fread(fid,'*char'))';      fclose(fid);
		lines = strread(c,'%s','delimiter','\n');   clear c fid;
		fname = [];
		for (k = 1:numel(lines))
			if (~strncmp(lines{k},'MIR_AQUAPLUG',7)),	continue,	end
			fname = ddewhite(lines{k}(13:end));
			if (exist(fname,'file') ~= 2)
				errordlg(['Script file for aquaPlugin ' fname ' does not exist. Ignoring request'],'Error')
				return
			end
		end
	else
		fname = auto;
	end

	if (isempty(fname)),	return,		end		% OPTcontrol has not a MIR_AQUAPLUG key. Something will screw 

	try				% Wrap it in a try-catch so we have a chance to figure out the reason of eventual error
		fid = fopen(fname, 'r');
		if (fid < 0)
			errordlg(['File ' fname ' does not exist or can''t be open'])
			out = [];		return
		end
		c = (fread(fid,'*char'))';      fclose(fid);
		out = strread(c,'%s','delimiter','\n');   clear c fid;
		ind = true(1,numel(out));
		for (k = 1:numel(out))
			if (isempty(out{k}) || out{k}(1) == '#'),	continue,	end		% Jump comments
			ind(k) = false;
		end
		out(ind) = [];		% Remove the comment lines (if any)

		[t, r] = strtok(out{1});	% Fish the case selected in the control file.
		out{1} = str2double(r);
		for (k = 2:numel(out))
			if (strncmpi(out{k},'char',4))	% If we have a string argument, use it as it is
				[t, r] = strtok(out{k});
				out{k} = ddewhite(r);
			else
				out{k} = eval(out{k});		% Numeric arguments must be 'evaled'
			end
		end
	catch
		errordlg(lasterr, 'Errror')
		out = [];
	end

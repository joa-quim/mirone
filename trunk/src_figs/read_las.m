function varargout = read_las(varargin)
% Helper window to read, filter and save or display LIDAR data read with libLAS

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

% We are now using the LASlib library instead of the older libLAS
% The mex readers, for the time being, have the same syntax so we
% can swapp the readers by simply toggling (search-replace) the mexs
%	laszreader_mex <-> lasreader_mex

% $Id: $

	hObject = figure('Vis','off');
	read_las_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	handles.hMirFig = [];
	handles.fname   = [];
	handles.bbox    = [];
	handles.classes = [];
	handles.IDs     = [];
	handles.xMinTouched = false;		% Used to know if we must send in -R option to mex
	handles.xMaxTouched = false;
	handles.yMinTouched = false;
	handles.yMaxTouched = false;
	handles.zMinTouched = false;
	handles.zMaxTouched = false;
	handles.classTypes = {'1 Unclassified' '2 Ground' '2 Low Vegetation' '4 Medium Vegetation'...
		'5 High Vegetation' '6 Building' '7 Low Point (noise)' '8 Model Key-point (mass point)' '9 Water'};

	if (~isempty(varargin))         % When called from a Mirone window
		if (ishandle(varargin{1}))
			handles.hMirFig = varargin{1};
		end
		if (numel(varargin) == 2)
			handles.fname = varargin{2};
			set(handles.edit_LASfile,'String',handles.fname)
			try			% Wrap it in a try because we are not sure file exists
				bbox = laszreader_mex(handles.fname,'-B');
				set(handles.edit_x_min,'Str',sprintf('%.12g',bbox(1)))
				set(handles.edit_x_max,'Str',sprintf('%.12g',bbox(2)))
				set(handles.edit_y_min,'Str',sprintf('%.12g',bbox(3)))
				set(handles.edit_y_max,'Str',sprintf('%.12g',bbox(4)))
				set(handles.edit_z_min,'Str',sprintf('%.9g',bbox(5)))
				set(handles.edit_z_max,'Str',sprintf('%.9g',bbox(6)))
				handles.bbox = bbox;
			catch
				errordlg(['File ' handles.fname ' does not exist'],'Error')
				delete(hObject),	return
			end
		end
	end

	if (~isempty(handles.hMirFig))
		handMir = guidata(handles.hMirFig);
		handles.home_dir = handMir.home_dir;
		handles.work_dir = handMir.work_dir;
		handles.last_dir = handMir.last_dir;
		handles.path_data = handMir.path_data;
		handles.whichFleder = handMir.whichFleder;
	else
		handles.home_dir = cd;
		handles.work_dir = cd;		handles.last_dir = cd;	% To not compromize put_or_get_file
		handles.path_data = [cd filesep 'data' filesep];	% Wil fail if called from outside Mir home
		handles.whichFleder = [];
	end

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, [handles.text_color handles.text_clip])
	%------------- END Pro look (3D) -----------------------------------------------------
	
	% ------------ Create a Semaforo -----------------------------------------------------
	[semaforo, pal] = aux_funs('semaforo_green');
	handles.hSemaforo = image(semaforo,'Parent', handles.axes1);
	set(handles.axes1, 'XTick',[], 'YTick', [])
	set(handles.figure1, 'Colormap', pal),		drawnow

	% Add this figure handle to the carraças list
	if (~isempty(handles.hMirFig))
		plugedWin = getappdata(handles.hMirFig,'dependentFigs');
		plugedWin = [plugedWin hObject];
		setappdata(handles.hMirFig,'dependentFigs',plugedWin);
	end

	guidata(hObject, handles);
	set(hObject,'Vis','on');
	if (nargout),	varargout{1} = hObject;		end

% -------------------------------------------------------------------------------------------------
function push_LASfile_CB(hObject, handles)
	[FileName,PathName] = put_or_get_file(handles, ...
		{'*.las;*.LAS;*.laz;*.LAZ', 'LIDAR file (*.las,*.LAS,*.laz,*.LAZ)';'*.*', 'All Files (*.*)'},'Select LAS file','get');
	if isequal(FileName,0),		return,		end
	% Let the edit_LASfile do the rest of the work;
	edit_LASfile_CB(handles.edit_LASfile, handles, [PathName FileName])

% -------------------------------------------------------------------------------------------------
function edit_LASfile_CB(hObject, handles, fname)
	if (nargin == 2),	fname = [];		end

	if (isempty(fname))		fname = get(hObject, 'Str');
	else					set(hObject, 'Str', fname)
	end

	if (exist(fname, 'file') == 2)
		handles.fname = fname;
	else
		errordlg(['File ' fname ' does not exist'],'Error')
		set(hObject, 'Str', '')
		handles.fname = [];
	end
	if (isempty(handles.bbox) && ~isempty(handles.fname))
		bbox = laszreader_mex(handles.fname,'-B');
		set(handles.edit_x_min,'Str',sprintf('%.12g',bbox(1))),	set(handles.edit_x_max,'Str',sprintf('%.12g',bbox(2)))
		set(handles.edit_y_min,'Str',sprintf('%.12g',bbox(3))),	set(handles.edit_y_max,'Str',sprintf('%.12g',bbox(4)))
		set(handles.edit_z_min,'Str',sprintf('%.9g',bbox(5))),	set(handles.edit_z_max,'Str',sprintf('%.9g',bbox(6)))
		handles.bbox = bbox;
	end
	guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------------------
function push_getClass_CB(hObject, handles)
	if (isempty(handles.classes) && ~isempty(handles.fname))
		semaforo_toggle(handles, 'red')
		handles.classes = laszreader_mex(handles.fname,'-S');
		set(handles.listbox1, 'Str', handles.classTypes(handles.classes))		% Info
		set(handles.popup_class, 'Str', handles.classes)		% For the case they will be wanted later
		guidata(handles.figure1, handles);
		semaforo_toggle(handles, 'green')
	end

% -------------------------------------------------------------------------------------------------
function radio_colorDepth_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_colorIntens handles.radio_colorClass handles.radio_colorReturn handles.radio_colorID],'Val',0)

% -------------------------------------------------------------------------------------------------
function radio_colorIntens_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_colorDepth handles.radio_colorClass handles.radio_colorReturn handles.radio_colorID],'Val',0)

% -------------------------------------------------------------------------------------------------
function radio_colorClass_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_colorDepth handles.radio_colorIntens handles.radio_colorReturn handles.radio_colorID],'Val',0)
	if (isempty(handles.classes) && ~isempty(handles.fname))	% We will need this info later
		semaforo_toggle(handles, 'red')
		handles.classes = laszreader_mex(handles.fname,'-S');
		set(handles.listbox1, 'Str', handles.classTypes(handles.classes))		% Info
		set(handles.popup_class, 'Str', handles.classes)		% For the case they will be wanted later
		guidata(handles.figure1, handles);
		semaforo_toggle(handles, 'green')
	end

% -------------------------------------------------------------------------------------------------
function radio_colorReturn_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_colorDepth handles.radio_colorIntens handles.radio_colorClass handles.radio_colorID],'Val',0)

% -------------------------------------------------------------------------------------------------
function radio_colorID_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_colorDepth handles.radio_colorIntens handles.radio_colorClass handles.radio_colorReturn],'Val',0)
	if (isempty(handles.IDs) && ~isempty(handles.fname))	% We will need this info later
		semaforo_toggle(handles, 'red')
		handles.IDs = laszreader_mex(handles.fname,'-D');
		guidata(handles.figure1, handles);
		semaforo_toggle(handles, 'green')
	end

% -------------------------------------------------------------------------------------------------
function check_clipClass_CB(hObject, handles)
	if (get(hObject,'Value'))		% We need to test several cases
		if (isempty(handles.fname))		set(hObject,'Val',0),	return,		end		% No file yet
		if (isempty(handles.classes))	% Don't know them yet, time to do it
			semaforo_toggle(handles, 'red')
			handles.classes = laszreader_mex(handles.fname,'-S');
			set(handles.listbox1, 'Str', handles.classtypes(handles.classes))		% Info
			guidata(handles.figure1, handles);
			semaforo_toggle(handles, 'green')
		end
		set(handles.popup_class, 'Str', handles.classes, 'Enable','on') 
	else
		set(handles.popup_class, 'Val',1, 'Str', ' ', 'Enable','off') 
	end

% -------------------------------------------------------------------------------------------------
function check_clipIntens_CB(hObject, handles)
	if (get(hObject,'Value'))		% We need to test several cases
		if (isempty(handles.fname))		set(hObject,'Val',0),	return,		end		% No file yet
		set(handles.edit_clipIntens, 'Enable','on') 
	else
		set(handles.edit_clipIntens, 'Enable','off') 
	end

% -------------------------------------------------------------------------------------------------
function edit_clipIntens_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1)
		set(hObject,'String',''),		set(handles.check_clipIntens,'Val',0)
		return
	end

% -------------------------------------------------------------------------------------------------
function check_clipReturns_CB(hObject, handles)
	if (get(hObject,'Value'))		% We need to test several cases
		if (isempty(handles.fname))		set(hObject,'Val',0),	return,		end		% No file yet
		set(handles.popup_returns, 'Enable','on') 
	else
		set(handles.popup_returns, 'Enable','off') 
	end

% -------------------------------------------------------------------------------------------------
function check_clipID_CB(hObject, handles)
	if (get(hObject,'Value'))
		warndlg('Not implemented yet.','Warning')
		set(hObject,'Val',0)
	end

% -------------------------------------------------------------------------------------------------
function edit_clipID_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)
		set(hObject,'String',''),		set(handles.check_clipID,'Val',0)
		return
	end

% -------------------------------------------------------------------------------------------------
function edit_x_min_CB(hObject, handles)
	handles.xMinTouched = true;		guidata(handles.figure1, handles)

% -------------------------------------------------------------------------------------------------
function edit_x_max_CB(hObject, handles)
	handles.xMaxTouched = true;		guidata(handles.figure1, handles)

% -------------------------------------------------------------------------------------------------
function edit_y_min_CB(hObject, handles)
	handles.yMinTouched = true;		guidata(handles.figure1, handles)

% -------------------------------------------------------------------------------------------------
function edit_y_max_CB(hObject, handles)
	handles.yMaxTouched = true;		guidata(handles.figure1, handles)

% -------------------------------------------------------------------------------------------------
function edit_z_min_CB(hObject, handles)
	handles.zMinTouched = true;		guidata(handles.figure1, handles)

% -------------------------------------------------------------------------------------------------
function edit_z_max_CB(hObject, handles)
	handles.zMaxTouched = true;		guidata(handles.figure1, handles)

% -------------------------------------------------------------------------------------------------
function push_BINfile_CB(hObject, handles)
	[FileName,PathName] = put_or_get_file(handles, ...
		{'*.bin;*.dat', 'Binary file (*.bin,*.dat)';'*.*', 'All Files (*.*)'},'Select Binary file','get');
	if isequal(FileName,0),		return,		end
	% Let the edit_BINfile do the rest of the work;
	edit_BINfile_CB(handles.edit_BINfile, handles, [PathName FileName])

% -------------------------------------------------------------------------------------------------
function edit_BINfile_CB(hObject, handles, fname)
	if (nargin == 2),	fname = [];		end

	if (isempty(fname))		fname = get(hObject, 'Str');
	else					set(hObject, 'Str', fname)
	end
% 	fname = [fname '.vtk'];

	semaforo_toggle(handles, 'red')
	xyz = get_data(handles);
	if (isempty(xyz))
		warndlg('No data points with current selection','Warning'),		return
	end

	frmt = 'real*4';
	if (get(handles.check_double, 'Val'))		frmt = 'real*8';	end
	fid = fopen(fname, 'wb');
	if (~isa(xyz,'cell'))
% 		fprintf(fid, '# vtk DataFile Version 2.0\n');
% 		fprintf(fid, 'converted from A B\n');
% 		fprintf(fid, 'BINARY\n');
% 		fprintf(fid,'DATASET POLYDATA\n');
% 		fprintf(fid,sprintf('POINTS %d float\n', length(xyz)));
		fwrite(fid,xyz,frmt);
	else							% Write a multisegment file
		for (k = 1:numel(xyz))
			fwrite(fid,[NaN NaN NaN],frmt);
			fwrite(fid,xyz{k},frmt);
		end
	end
	fclose(fid);
	semaforo_toggle(handles, 'green')

% -------------------------------------------------------------------------------------------------
function push_grdTool_CB(hObject, handles)
% Nikles - invisible

% -------------------------------------------------------------------------------------------------
function push_FLEDERfile_CB(hObject, handles)
	[FileName,PathName] = put_or_get_file(handles, ...
		{'*.sd;*.scene', 'Binary file (*.sd,*.scene)';'*.*', 'All Files (*.*)'},'Select Fleder file','get');
	if isequal(FileName,0),		return,		end
	% Let the edit_FLEDERfile do the rest of the work;
	edit_FLEDERfile_CB(handles.edit_FLEDERfile, handles, [PathName FileName])

% -------------------------------------------------------------------------------------------------
function edit_FLEDERfile_CB(hObject, handles, fname)
	if (nargin == 2),	fname = [];		end

	semaforo_toggle(handles, 'red')
	xyz = get_data(handles);
	if (isempty(xyz))
		warndlg('No data points with current selection','Warning'),		return
	end
	
	if (isempty(fname))		fname = get(hObject, 'Str');
	else					set(hObject, 'Str', fname)
	end
	write_flederFiles('points', fname, xyz, 'first', handles.bbox);
	semaforo_toggle(handles, 'green')

% -------------------------------------------------------------------------------------------------
function push_goFleder_CB(hObject, handles)
% Create one Fleder file in ...$mirone/tmp and call Fleder on it.

	fname = [handles.home_dir filesep 'tmp' filesep 'lixo.sd'];
	semaforo_toggle(handles, 'red')
	xyz = get_data(handles);
	if (isempty(xyz))
		warndlg('No data points with current selection','Warning'),		return
	end
	write_flederFiles('points', fname, xyz, 'first', handles.bbox);		pause(0.1)
	semaforo_toggle(handles, 'green')
	
	if (isempty(handles.whichFleder))	handles.whichFleder = 1;	end
	comm = [' -data ' fname ' &'];
	if (handles.whichFleder),	fcomm = ['iview4d' comm];			% Free viewer
	else						fcomm = ['fledermaus' comm];		% The real thing
	end
	try
		if (isunix)				% Stupid linux doesn't react to a non-existant iview4d
			resp = unix(fcomm);
			if (resp == 0)
				errordlg('I could not find Fledermaus. Hmmm, do you have it?','Error')
			end
		elseif (ispc)	dos(fcomm);
		else			errordlg('Unknown platform.','Error'),	return
		end
	catch
		errordlg('I could not find Fledermaus. Hmmm, do you have it?','Error')
	end

	pause(1)
	builtin('delete',fname);

% -------------------------------------------------------------------------------------------------
function xyz = get_data(handles)
% Get the data points from file taking into account possible filter settings

	[out, colorBy] = parse_before_go(handles);
	if (~colorBy.do)			% Simpler case. No data spliting
		if (~isempty(out))		xyz = laszreader_mex(handles.fname, out{:});
		else					xyz = laszreader_mex(handles.fname);
		end
	elseif (colorBy.do && colorBy.split)
		% Here we have to deal with a convoluted logic as we'll allow Clip options as well
		if (colorBy.class)
			xyz = cell(1, numel(handles.classes));
			for (k = 1:numel(handles.classes))
				opt_C = sprintf('-C%d', handles.classes(k));
				if (~isempty(out) && out{1}(2) == 'C')		% Ghrr, already had one -C request.
					out(1) = [];
				end
				if (~isempty(out))							% Case with Clippings too
					xyz{k} = laszreader_mex(handles.fname, opt_C, out{:});
				else
					xyz{k} = laszreader_mex(handles.fname, opt_C);
				end
			end
		elseif (colorBy.ID)
			xyz = cell(1, numel(handles.IDs));
			for (k = 1:numel(handles.IDs))
				opt_D = sprintf('-D%d', handles.IDs(k));
				if (~isempty(out) && out{end}(2) == 'D')	% Ghrr, already had one -D request.
					out(end) = [];
				end
				if (~isempty(out))							% Case with Clippings too
					xyz{k} = laszreader_mex(handles.fname, opt_D, out{:});
				else
					xyz{k} = laszreader_mex(handles.fname, opt_D);
				end
			end
		end
	end

% ---------------------------------------------------------------------------------------
function [out, colorBy] = parse_before_go(handles)
% Parse all possible constraints 

	out = [];	colorBy.do = false;		colorBy.split = false;
	if (get(handles.check_clipClass, 'Val'))
		str = get(handles.popup_class, 'Str');		val = get(handles.popup_class,'Val');
		out{1} = ['-C' str(val,:)];			%Option -C
	end
	if (get(handles.check_clipIntens, 'Val'))
		xx = round(abs(str2double(get(handles.edit_clipIntens,'String'))));
		if (isnan(xx))
			warndlg('Trash or empty in the Intensity edit box. Ignoring this request','Warning')
		else
			out{end+1} = sprintf('-I%d', xx);	%Option -I
		end
	end
	if (get(handles.check_clipReturns, 'Val'))
		n = get(handles.popup_returns,'Val');
		if (n == 2)		n = 10;		end			% Last return is coded 10
		out{end+1} = sprintf('-C%d',n);			%Option -N
	end

	if (any([handles.xMinTouched handles.xMaxTouched handles.yMinTouched		% Before we go for -R
			handles.yMaxTouched handles.zMinTouched handles.zMaxTouched]))		% make sure it worth
		x1 = str2double(get(handles.edit_x_min,'Str'));		x2 = str2double(get(handles.edit_x_max,'Str'));
		y1 = str2double(get(handles.edit_y_min,'Str'));		y2 = str2double(get(handles.edit_y_max,'Str'));
		z1 = str2double(get(handles.edit_z_min,'Str'));		z2 = str2double(get(handles.edit_z_max,'Str'));
		if ( any(abs( [x1 x2 y1 y2] - handles.bbox(1:4) ) > 1e-6) &&  all(abs([z1 z2] - handles.bbox(5:6)) < 1e-4) )
			out{end+1} = sprintf('-R%.8f/%.8f/%.8f/%.8f',x1, x2, y1, y2);		% A 2D -R option
		elseif ( any(abs( [x1 x2 y1 y2 z1 z2 ] - handles.bbox ) > 1e-6) )
			out{end+1} = sprintf('-R%.8f/%.8f/%.8f/%.8f/%.8f/%.8f',x1, x2, y1, y2, z1, z2);	% The 3D -R option
		end
	end

	if (get(handles.check_clipID, 'Val'))
		xx = round(abs(str2double(get(handles.edit_clipID,'String'))));
		if (isnan(xx))
			warndlg('Trash or empty in the Source ID edit box. Ignoring this request','Warning')
		else
			out{end+1} = sprintf('-D%d', xx);	%Option -D
		end
	end

	if (get(handles.radio_colorClass, 'Val'))
		colorBy.do = true;		colorBy.split = true;
		colorBy.class = true;	colorBy.ID = false;
	elseif (get(handles.radio_colorID, 'Val'))
		colorBy.do = true;		colorBy.split = true;
		colorBy.ID = true;		colorBy.class = false;
	end

% ---------------------------------------------------------------------------------------
function semaforo_toggle(handles, cor)
% Toggle semaforo light between red and green

	if (cor(1) == 'r')
		[semaforo, pal] = aux_funs('semaforo_red');
	else
		[semaforo, pal] = aux_funs('semaforo_green');
	end
	set(handles.hSemaforo,'CData',semaforo)
	set(handles.figure1, 'Colormap', pal),		drawnow

% ---------------------------------------------------------------------------------------
function read_las_LayoutFcn(h1)

set(h1, 'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Read LAS',...
'NumberTitle','off',...
'Position',[520 396 421 405],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'Visible','on');

uicontrol('Parent',h1, 'Position',[10 124 401 122],'Style','frame');
uicontrol('Parent',h1, 'Position',[270 258 141 111],'Style','frame');

axes('Parent',h1,'Units','pixels','Position',[375 65 14 46],...
'XTick', [], 'YTick', [], 'Visible', 'off', 'Tag','axes1');

uicontrol('Parent',h1, 'Position',[10 382 31 15],...
'HorizontalAlignment','left',...
'String','File',...
'Style','text');

uicontrol('Parent',h1, 'Position',[40 380 351 21],...
'BackgroundColor',[1 1 1],...
'Call',@read_las_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tooltip','Enter LS file name',...
'Tag','edit_LASfile');

uicontrol('Parent',h1, 'Position',[387 379 23 23],...
'Call',@read_las_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'Tag','push_LASfile');

uicontrol('Parent',h1, 'Position',[10 348 231 21],...
'Call',@read_las_uiCB,...
'String','Get Classifications (opt)',...
'Tooltip','Scan LAS file to find and retrieve existing classifications list',...
'Tag','push_getClass');

uicontrol('Parent',h1, 'Position',[10 258 231 82],...
'BackgroundColor',[1 1 1],...
'String',' ',...
'Style','listbox',...
'Tag','listbox1');

uicontrol('Parent',h1, 'Position',[281 340 70 21],...
'Call',@read_las_uiCB,...
'String','Depth',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_colorDepth');

uicontrol('Parent',h1, 'Position',[281 320 70 23],...
'Call',@read_las_uiCB,...
'String','Intensity',...
'Enable','off',...
'Style','radiobutton',...
'Tag','radio_colorIntens');

uicontrol('Parent',h1, 'Position',[281 300 90 23],...
'Call',@read_las_uiCB,...
'String','Classification',...
'Style','radiobutton',...
'Tag','radio_colorClass');

uicontrol('Parent',h1, 'Position',[281 280 80 23],...
'Call',@read_las_uiCB,...
'String','Returns',...
'Enable','off',...
'Style','radiobutton',...
'Tag','radio_colorReturn');

uicontrol('Parent',h1, 'Position',[281 261 80 21],...
'Call',@read_las_uiCB,...
'String','Source ID',...
'Style','radiobutton',...
'Tag','radio_colorID');

uicontrol('Parent',h1, 'Position',[20 215 90 21],...
'Call',@read_las_uiCB,...
'String','Classification',...
'Style','checkbox',...
'Tooltip','Select points by classification',...
'Tag','check_clipClass');

uicontrol('Parent',h1, 'Position',[110 215 51 20],...
'BackgroundColor',[1 1 1],...
'String',' ',...
'Style','popupmenu',...
'Tooltip','Select classification',...
'Value',1,...
'Enable','off',...
'Tag','popup_class');

uicontrol('Parent',h1, 'Position',[20 186 70 21],...
'Call',@read_las_uiCB,...
'String','Intensity',...
'Style','checkbox',...
'Tag','check_clipIntens');

uicontrol('Parent',h1,'Position',[110 187 51 21],...
'BackgroundColor',[1 1 1],...
'Call',@read_las_uiCB,...
'Style','edit',...
'Tooltip','Intensities lower than this will be clipped',...
'Tag','edit_clipIntens');

uicontrol('Parent',h1, 'Position',[20 159 70 21],...
'Call',@read_las_uiCB,...
'String','Returns',...
'Style','checkbox',...
'Tooltip','Choose First or Last return',...
'Tag','check_clipReturns');

uicontrol('Parent',h1, 'Position',[110 160 51 20],...
'BackgroundColor',[1 1 1],...
'String',{'First'; 'Last'},...
'Style','popupmenu',...
'Tooltip','Choose First or Last return',...
'Enable','off',...
'Value',1,...
'Tag','popup_returns');

uicontrol('Parent',h1, 'Position',[20 132 80 21],...
'Call',@read_las_uiCB,...
'String','Source ID',...
'Style','checkbox',...
'Tag','check_clipID');

uicontrol('Parent',h1, 'Position',[110 131 51 21],...
'BackgroundColor',[1 1 1],...
'Call',@read_las_uiCB,...
'Style','edit',...
'Tooltip','Source IDs different that this are removed',...
'Tag','edit_clipID');

uicontrol('Parent',h1, 'Position',[236 191 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@read_las_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tag','edit_x_min');

uicontrol('Parent',h1, 'Position',[322 191 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@read_las_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tag','edit_x_max');

 uicontrol('Parent',h1,'Position',[236 165 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@read_las_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tag','edit_y_min');

uicontrol('Parent',h1, 'Position',[322 165 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@read_las_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tag','edit_y_max');

 uicontrol('Parent',h1,'Position',[181 195 55 15],...
'Enable','inactive',...
'HorizontalAlignment','left',...
'String','X Direction',...
'Style','text');

uicontrol('Parent',h1, 'Position',[180 169 55 15],...
'Enable','inactive',...
'HorizontalAlignment','left',...
'String','Y Direction',...
'Style','text');

uicontrol('Parent',h1, 'Position',[342 212 41 13],...
'Enable','inactive',...
'String','Max',...
'Style','text');

uicontrol('Parent',h1, 'Position',[258 212 41 13],...
'Enable','inactive',...
'String','Min',...
'Style','text');

uicontrol('Parent',h1, 'Position',[236 136 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@read_las_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tag','edit_z_min');

uicontrol('Parent',h1, 'Position',[322 136 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@read_las_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tag','edit_z_max');

uicontrol('Parent',h1, 'Position',[180 140 55 15],...
'HorizontalAlignment','left',...
'String','Z Direction',...
'Style','text');

uicontrol('Parent',h1, 'Position',[316 358 60 18],...
'FontSize',9,...
'String','Color by:',...
'Style','text',...
'Tag','text_color');

uicontrol('Parent',h1, 'Position',[190 235 60 18],...
'FontSize',9,...
'String','Clip by:',...
'Style','text',...
'Tag','text_clip');

uicontrol('Parent',h1, 'Position',[10 78 311 21],...
'BackgroundColor',[1 1 1],...
'Call',@read_las_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tooltip','Enter output file (binary) name',...
'Tag','edit_FLEDERfile');

uicontrol('Parent',h1, 'Position',[320 77 23 23],...
'Call',@read_las_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'Tag','push_FLEDERfile');

uicontrol('Parent',h1, 'Position',[10 34 311 21],...
'BackgroundColor',[1 1 1],...
'Call',{@read_las_uiCB,h1,'edit_BINfile_CB'},...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tooltip','Enter output file (binary) name',...
'Tag','edit_BINfile');

uicontrol('Parent',h1, 'Position',[320 33 23 23],...
'Call',@read_las_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'Tag','push_BINfile');

uicontrol('Parent',h1, 'Position',[10 56 110 15],...
'HorizontalAlignment','left',...
'String','Save to Binary File',...
'Style','text');

uicontrol('Parent',h1, 'Position',[350 33 61 23],...
'String','Double',...
'Style','checkbox',...
'Tooltip','Check to write in doube precision (def is single)',...
'Tag','check_double');

uicontrol('Parent',h1, 'Position',[10 4 110 21],...
'Call',@read_las_uiCB,...
'String','Call Gridding Tool',...
'Tag','push_grdTool',...
'Visible','off');

uicontrol('Parent',h1, 'Position',[230 4 181 21],...
'Call',@read_las_uiCB,...
'String','View in Fledermaus',...
'Tooltip','Visualize LAS data with Fledermaus',...
'Tag','push_goFleder');

uicontrol('Parent',h1, 'Position',[10 100 110 15],...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'String','Write Fledermaus File ',...
'Style','text');

function read_las_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

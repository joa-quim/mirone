function varargout = empilhador(varargin)
% Stacks a bunch of grids into a single 3D file
%
% WARNING: FOR COMPILING THIS WE NEED TO INCLUDE THE HDF_FUNS.M SRC
%
% NOTE: The gotFromMETA and getZ functions are callable directly by mirone

%	Copyright (c) 2004-2014 by J. Luis
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

	if (nargin > 1 && ischar(varargin{1}))
		gui_CB = str2func(varargin{1});
		[varargout{1:nargout}] = feval(gui_CB,varargin{2:end});
	else
		h = empilhador_OF(varargin{:});
		if (nargout),	varargout{1} = h;   end
	end

% ---------------------------------------------------------------------------------
function hObject = empilhador_OF(varargin)
	hObject = figure('Vis','off');
	empilhador_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	if (~isempty(varargin) && isa(varargin{1},'struct'))
		handMir = varargin{1};
		handles.home_dir = handMir.home_dir;
		handles.last_dir = handMir.last_dir;
		handles.work_dir = handMir.work_dir;
		handles.deflation_level = handMir.deflation_level;
		handles.IamCompiled = handMir.IamCompiled;		% Need to know due to crazy issue of nc_funs
		handles.path_tmp = handMir.path_tmp;
		handles.path_data = handMir.path_data;
	else
		mir_dirs = getappdata(0,'MIRONE_DIRS');
		if (~isempty(mir_dirs))
			handles.home_dir = mir_dirs.home_dir;		% Start in values
			handles.work_dir = mir_dirs.work_dir;
			handles.last_dir = mir_dirs.last_dir;
			handles.path_tmp = [handles.home_dir filesep 'tmp' filesep];
			handles.path_data = [handles.home_dir filesep 'data' filesep];
		else		
			handles.home_dir = cd;
			handles.last_dir = handles.home_dir;
			handles.work_dir = handles.home_dir;
			handles.path_tmp = [pwd filesep 'tmp' filesep];
			handles.path_data = [pwd filesep 'data' filesep];
		end
		handles.IamCompiled = false;
	end
	handles.nameList = [];
	handles.OneByOneNameList = [];	% For when files are loaded one by one (risky)
	handles.OneByOneFirst = true;	% Safety valve to deal with the load one by one case
	handles.testedDS = false;		% To test if a Sub-Dataset request is idiot
	handles.automatic = false;		% Set when called from command line and create a Fig
	handles.Interactive = false;	% Used when need to reinterpolate L2 files. FALSE means do not
									% call the helper window that asks questions (use default ans)

	% -------------- Import/set icons --------------------------------------------
	load([handles.path_data 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_namesList, 'CData',Mfopen_ico)

	%------------ Give a Pro look (3D) to the frame box -----------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) ------------------------------

	set(hObject,'Visible','on');
	guidata(hObject, handles);

	% When calling empilhador by command line (that is, in a non-interactive way) but one that creates a Fig
	if (~isempty(varargin) && ischar(varargin{1}))
		handles.automatic = true;
		push_namesList_CB(handles.push_namesList, handles, varargin{:})
	end

% -----------------------------------------------------------------------------------------
function edit_namesList_CB(hObject, handles)
	fname = get(hObject,'String');
	push_namesList_CB(handles.push_namesList, handles, fname)

% -----------------------------------------------------------------------------------------
function push_namesList_CB(hObject, handles, opt)
% Read an ascii file with a list of names to process. The file list may have 1,2 or 3 columns
%
% OPT	contains either the filename with the fields described bellow or is empty, in which
%		case that name will be asked here.
%		Alternatively it can contain a wildcard string that will be expanded to create the
%		filename automatically in the same dir of the data. The contents of the wildcard should
%		follow what is explained below, but three possible examples are provided right now (no quotes)
%		'C:\a1\*.L2_LAC_SST4 ? -sst4'
%		'C:\a1\*.hdf * -qual'
%		'C:\a1\pathfinder\2010\201*.nc ? sds1'
%
% First column holds the filename that can be absolute, relative or have no path info.
%		If only one column in file the "Time" info below is computed as (1:number_of_files)
%
% Second column is normally the "Time" info. That is, a number that will be used as the 'time'
%		in the 3D netCDF file. In case of L2 scene products one can use '?' to instruct the
%		program to get the time from the HDF file name, which is assume to be YYYYDDDHHMMSS....
%		The '?' mechanism has been extended to the GHRSST PATHFINDER 5.2 netCDF files that have
%		names of the form YYYYMMDDHHMMSS-...
%
% Third column is used when the HDF file has sub-datasets and we want to select one in
%		particular. In that case use the (clumsy) construct: 'sdsN' as in sds3 where 'sds' is
%		fix and N is the number of the sub-dataset as it appears when one do gdalinfo or in the
%		table window that pops up when one attempts to open the HDF file directly in Mirone.
%		NEWs: It is now possible to provide the SDS name instead of the 'sdsN' described above.
%		But to distinguish both mechanisms the SDS name must be pre-pended with a '-', as in -sst4
%
% The list file may still have one header line (a line that starts with a '#') with any of these:
%	-F<flagSDS_minQuality>
%		Which is a composition of two pieces of information. The first, 'flagSDS' is the
%		sub-dataset name of a quality flags array to be applied to the main dataset.
%		Currently this option applies only to PATHFINDER 5.2 netCDF files.
%	-R<region>
%		A GMT style region selector: e.g. -R-20/-5/33/45
%		Using this option is equivalent to fill in the "Use sub-region" boxes
%	-G<output>
%		Give the output file name. If no path is provided the output is written in the same
%		directory as the data files
%	-D<description> or -D"desc ription"
%		A description string that will go into the global attribute 'description'.
%		Use the quoted form if the description has more than one word.
%
% The above options may be combined in a single command that can be executed from the command line.
%	Example:
%	empilhador('empilhador_OF','C:\a1\pathfinder\2010\201*.nc ? sds1 -R-15/-5/35/45 -Fquality_level_4 -Glixo.nc');
%
% Note that for this case to work we had to call it through the opening function, 'empilhador_OF',
% because we need to create one Figure and use the handles structure.
% When -G is used as above, the program executes automatically till the end and deletes the Fig
%
% Besides '#' header line, the list file may include other files by using the notation (one per line)
%	>include full_path_to_the_include_file
%
%	In this case no testing is done in the included file and it is assumed that it has EXACTLY the same
%	construction of the main list file or, in case there are only include directives, that they ALL
%	share the same format in terms of number of columns. Any eventual lines starting with '#' in the
%	included files are simply ingnored.


	if (nargin == 2)		% Direct call
		str1 = {'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.DAT,*.txt,*.TXT)';'*.*', 'All Files (*.*)'};
		[FileName,PathName,handles] = put_or_get_file(handles, str1,'File with grids list','get');
		if isequal(FileName,0),		return,		end
	else					% File name on input
		opt = check_wildcard_fname(opt);	% Check if a lazy 'path/*.xxx X X' request
		[PathName,FNAME,EXT] = fileparts(opt);
		PathName = [PathName filesep];      % To be coherent with the 'if' branch
		FileName = [FNAME EXT];
	end
	fname = [PathName FileName];

	[bin, n_column] = guess_file(fname);

	if isempty(bin)					% If error in reading file
		errordlg(['Error reading file ' fname],'Error'),	return
	elseif (n_column == 0)			% A hole in guess_file logic (error message already issued)
		return
	elseif (bin)					% Binary file. Assume it's a target file, not a name list
		[PATH,FNAME,EXT] = fileparts(fname);
		handles.OneByOneNameList{end+1} = fname;		% Save the full name
		str = get(handles.listbox_list, 'Str');
		str{end+1} = [FNAME EXT];
		set(handles.listbox_list, 'Str', str);
		handles.OneByOneFirst = true;					% Repetitive but ensures that things are always updated
		guidata(handles.figure1, handles)
		return
	end

	[handles, names, n_column] = parse_list_file(handles, fname, n_column);

	m = length(names);
	handles.strTimes = cell(m,1);		% To hold time steps as strings
	SDSinfo = cell(m,1);				% To hold Sub Datasets info
	handles.SDSinfo = [];				% If above exists, it will be copied here
	c = false(m,1);
	caracol = false(m,1);				% Case name list has '@' to pause for a CD change

	n_msg = 1;							% Will hold the "change CD messages" counter
	if (n_column > 1)					% When 2nd column holds the 3D numbering
		for (k = 1:m)
			[t,r] = strtok(names{k});
			if (numel(t) < 2 || t(1) == '#')	% Jump empty or comment lines
				c(k) = true;
				continue
			elseif (t(1) == '@')
				caracol(k) = true;
				if (~isempty(t(2:end))),	handles.changeCD_msg{n_msg} = t(2:end);		% The '@' was glued with the message
				else						handles.changeCD_msg{n_msg} = r;
				end
				n_msg = n_msg + 1;
				continue
			end

			names{k} = ddewhite(t);
			if (n_column == 2)			% Names & numeric label format OR names & sdsname
				r = ddewhite(r);
				if (r(1) == '?')
					r = squeeze_time_from_name(names{k});
				end
				handles.strTimes{k} = r;
			else						% Names, numeric label & SDS info format
				[t,r] = strtok(r);
				t = ddewhite(t);
				if (t(1) == '?')		% Means get the numeric label as time extracted from file name (OceanColor products)
					t = squeeze_time_from_name(names{k});
				elseif (t(1) == '*')	% 
					t = sprintf('%d', k);
				end
				handles.strTimes{k} = t;
				SDSinfo{k} = ddewhite(r);
			end
		end
	else								% Only one column with fnames
		for (k = 1:m)
			if ( isempty(names{k}) ),	continue,			end		% Jump empty lines
			if ( names{k}(1) == '#'),	c(k) = true;		continue,	end
			if ( names{k}(1) == '@')
				caracol(k) = true;
				handles.changeCD_msg{n_msg} = names{k}(2:end);
				n_msg = n_msg + 1;
				continue
			end
			handles.strTimes{k} = sprintf('%d',k);
		end
	end

	if (any(c))					% Remove eventual comment lines
		names(c) = [];			handles.strTimes(c) = [];		caracol(c) = [];	SDSinfo(c) = [];
	end
	m = numel(names);			% Count remaining ones
	
	if (m == 0)
		errordlg('Beautiful service. Your list file had nothing but a bunch of trash.','Error'),	return
	end

	% -------------- Check if we have a Sub-Datasets request ------------------
	if ( n_column == 3 && (strncmpi(SDSinfo{1}, 'sds', 3) || SDSinfo{1}(1) == '-') )
		handles.SDSinfo = SDSinfo;
	elseif (n_column == 2 && (strncmpi(handles.strTimes{1}, 'sds', 3) || handles.strTimes{1}(1) == '-') )
		% Two cols with SDS info in the second
		handles.SDSinfo = handles.strTimes;
		for (k = 1:m),		handles.strTimes{k} = sprintf('%d',k);		end
	end
	if (~isempty(handles.SDSinfo))
		if (strncmpi(handles.SDSinfo{1}, 'sds', 3))				% Old way. Kept for backward compatibility
			handles.SDSthis = str2double(handles.SDSinfo{1}(4:end));
			handles.SDSinfo = handles.SDSthis;					% Make it numeric to help test inside att_get()
		else
			handles.SDSthis = handles.SDSinfo{1}(2:end);		% Should have the SDSname prepended with a '-' sign
		end
	end
	% --------------------------------------------------------------------------

	handles.shortNameList = cell(m,1);      % To hold grid names with path striped
	for (k = 1:m)
		[PATH,FNAME,EXT] = fileparts(names{k});
		if (isempty(PATH))
			handles.shortNameList{k} = names{k};
			names{k} = [PathName names{k}];
		else
			handles.shortNameList{k} = [FNAME EXT];
		end
	end

	% Check that the files provided in list do exist
	if (~any(caracol))			% It makes no sense to test existance of files in other DVDs
		c = false(m,1);
		for (k = 1:m)
			c(k) = (exist(names{k},'file') ~= 2);		% Flag to kill all non-existant files
		end
		names(c) = [];		handles.shortNameList(c) = [];	handles.strTimes(c) = [];
	end

	handles.nameList = names;
	handles.caracol = caracol;
	set(handles.edit_namesList, 'String', fname)
	set(handles.listbox_list,'String',handles.shortNameList)
	guidata(handles.figure1,handles)

	if (~isempty(handles.outname) && handles.automatic)	% Automatic + output name means start the job immediately
		push_compute_CB(handles.push_compute, handles)
		delete(handles.figure1)
	end

% -----------------------------------------------------------------------------------------
function [handles, names, n_column] = parse_list_file(handles, fname, n_column)
% Read a list file, call parse_header() to check for options in header line
% Search for '>include fname' lines that instruct to include other list files
% We need to send in the N_COLUMNS because the main list file may be made of only
% '>include' directives, in which case we don't know yet the true columns number.

	fid = fopen(fname);
	if (fid < 0)
		error('EMPILHADOR:PARSE_LIST_FILE', 'Error opening list file')
	end
	c = fread(fid,'*char')';	fclose(fid);
	ind = strfind(c,'>');
	names = strread(c,'%s','delimiter','\n');

	handles.SDSflag = [];				% For when we want to apply also a quality flag (GHRSST .nc)
	handles.outname = [];				% For when -G<outname> is used in the header
	handles.desc_attrib = '';			% For when -D"description" is used in the header

	handles = parse_header(handles, names);		% Check if we have a header line with further instructions

	if (~isempty(ind))				% We likely have include instructions
		% We need to do some gymnastic here to include the includes in the right order and remove the '> lines
		names_copy = '';	inc_name = '';
		ini = 1;
		for (k = 1:numel(names))
			off = 0;
			if (strncmp(names{k},'>include',8) || strncmp(names{k},'> include',9))
				if (strncmp(names{k},'> include',9)),	off = 1;	end		% Accept also this variant
				inc_name = ddewhite(names{k}(9+off:end));
				fid = fopen(inc_name);
				if (fid < 0),		continue,	end			% Could not open file
				c = fread(fid,'*char')';		fclose(fid);
				names_ = strread(c,'%s','delimiter','\n');
				names_copy = [names_copy(1:end); names(ini:k-1); names_];		% Update list but without the '>' line
				ini = k + 1;								% Next time it will start on the line after the '>'
			end
		end
		if (ini < numel(names))
			names_copy = [names_copy; names(ini:end)];		% Add the rest after the last '>include' line
		end
		if (~isempty(names_copy)),		names = names_copy;	end

		if (n_column == 1 && ~isempty(inc_name))		% Very likely a situation where the main list file has only
			[bin, n_column] = guess_file(inc_name);		% '>include' lines, so count the columns from one of the included
		end
	end

% -----------------------------------------------------------------------------------------
function handles = parse_header(handles, names)
% Check if we have a header line with any of -F<flags>, -R<region>, -G<output> or -D<description>
% Note: -D"bla bla" is also valid

	if (names{1}(1) ~= '#')			% No header
		return
	end

	ind = strfind(names{1}, '-F');	% Do we have a Flag request?
	if (~isempty(ind))
		t = strtok(names{1}(ind+2:end));
		ind  = strfind(t, '_');
		handles.SDSflag = t(1:ind(end)-1);	% SDS number of the quality flags array (SubDataset)
		handles.minQuality = round(sscanf(t(ind(end)+1:end), '%f'));
		if (isnan(handles.minQuality))
			warndlg('The quality value in -F quality option is obviously wrong or non-existent. Ignoring.','Warning')
			handles.SDSflag = [];
		end
	end

	ind = strfind(names{1}, '-R');	% Do we have a region request?
	if (~isempty(ind))
		t = strtok(names{1}(ind+2:end));
		[A, count, errmsg] = sscanf(t,'%f/%f/%f/%f');
		if (count ~= 4 || ~isempty(errmsg))
			warndlg('Error: the -R string is badly formed. Ignoring it','WarnError')
		else
			if (A(1) >= A(2) || A(3) >= A(4))
				warndlg(['Error: either West or South are > than East or North in -R option ->' t],'WarnError')
			else
				set(handles.edit_west, 'String', A(1)); 
				set(handles.edit_east, 'String', A(2)); 
				set(handles.edit_south,'String', A(3)); 
				set(handles.edit_north,'String', A(4)); 
				set(handles.check_region, 'Val', 1)
				set([handles.edit_north handles.edit_south handles.edit_west handles.edit_east],'Enable','on')
			end
		end
	end

	ind = strfind(names{1}, '-G');	% Do we have an output filename request?
	if (~isempty(ind))
		handles.outname = strtok(names{1}(ind+2:end));
	end

	ind = strfind(names{1}, '-D');	% Do we have description string request?
	if (~isempty(ind))
		if (names{1}(ind+2) == '"')		% See if we have a description enclosed by quotes
			ind2 = strfind(names{1}(ind+3:end), '"') + ind + 2;	% Seek for the other quote (")
			handles.desc_attrib = names{1}(ind+3:ind2-1);
		else
			handles.desc_attrib = strtok(names{1}(ind+2:end));
		end
	end

% -------------------------------------------------------------------------------------------
function fname = check_wildcard_fname(strin)
% Check if user gave a 'path/*.xxx ? -sdsname [OPTs]' on input and if yes, create a resp file

	if ( isempty(strfind(strin, '*')) )		% Normal fname input. Nothing to do here
		fname = strin;		return
	end

	% -------------- First check if any of -D, -F, -R or -G is also present --------------
	opt_F = [];		opt_G = [];		opt_R = [];		opt_D = [];
	ind = strfind(strin, '-F');
	if (~isempty(ind))
		[opt_F, r] = strtok(strin(ind:end));
		strin = [strin(1:ind-2) r];		% Remove the -F<...> string from input
	end
	ind = strfind(strin, '-G');
	if (~isempty(ind))
		[opt_G, r] = strtok(strin(ind:end));
		strin = [strin(1:ind-2) r];		% Remove the -G<...> string from input
	end
	ind = strfind(strin, '-R');
	if (~isempty(ind))
		[opt_R, r] = strtok(strin(ind:end));
		strin = [strin(1:ind-2) r];		% Remove the -R<...> string from input
	end
	ind = strfind(strin, '-D');
	if (~isempty(ind))
		if (strin(ind+2) == '"')		% See if we have a description enclosed by quotes
			ind2 = strfind(strin(ind+3:end), '"') + ind + 2;	% Seek for the other quote (")
			opt_D = strin(ind+3:ind2-1);
			strin(ind:ind2) = [];			% Remove the -D"..." string from input
		else
			[opt_D, r] = strtok(strin(ind:end));
			strin = [strin(1:ind-2) r];		% Remove the -D<...> string from input
		end
	end
	% --------------------------------------------------------------------------------

	[t, r] = strtok(strin);
	dirlist = dir(t);
	if (isempty(dirlist))
		errordlg('Wildcard search return an empty result. Have to abort here','ERROR')
		error('Wildcard search return an empty result. Have to abort here')
	end

	PATO = fileparts(t);
	fname = [PATO '/' 'automatic_list.txt'];
	fid = fopen(fname,'w');
	if (fid < 0)
		error('Fail to create the resp file in data directory. Permissions problem?')
	end
	if ( ~isempty(opt_F) || ~isempty(opt_G) || ~isempty(opt_R) || ~isempty(opt_D) )	% Insert -D, -F, -G or -R options
		fprintf(fid, '# %s\n', [opt_F ' ' opt_G ' ' opt_R ' ' opt_D]);
	end

	% Here we assume that the OS returns a list already lexically sorted
	[t, r] = strtok(r);
	for (k = 1:numel(dirlist))
		fprintf(fid, '%s/%s', PATO, dirlist(k).name);	% The file name
		if (~isempty(t))
			fprintf(fid, '\t%s', t);					% Normally the "Time" info
			if (~isempty(r))
				fprintf(fid, '\t%s\n', r);				% The SDS name or number
			else
				fprintf(fid, '\n');
			end
		else
			fprintf(fid, '\n');
		end
	end
	fclose(fid);

% -----------------------------------------------------------------------------------------
function t = squeeze_time_from_name(name)
% ... Read the name of L2 daily scene product and convert it into a time string.
% The name algo is simple YYYYDDDHHMMSS where DDD is day of the year.

	% Example names: A2012024021000.L2_LAC_SST4 S1998001130607.L2_MLAC_OC.x.hdf
	% 20100109005439-NODC-L3C_GHRSST-SSTskin-AVHRR_Pathfinder-PFV5.2_NOAA18_G_2010009_night-v02.0-fv01.0.nc
	[PATH,FNAME,EXT] = fileparts(name);
	indDot = strfind(FNAME,'.');
	if (~isempty(indDot) && strcmpi(FNAME(16:17), 'L2'))	% Second case type name
		FNAME(indDot(1):end) = [];
	elseif (~isempty(EXT) && strcmpi(EXT(2:3), 'L2'))		% First case type name (nothing to do)
	elseif (strcmpi(EXT,'.nc') && FNAME(15) == '-')			% A GS... 2.0 netCDF PATHFINDER file
		FNAME(15:end) = [];
	else
		errordlg(sprintf('This "%s" is not a MODIS type name',name),'ERROR'),	return
	end
	% Compose name as YYYY.xxxxx where 'xxxxx' is the decimal day of year truncated to hour precision (or minute)
	if (double(FNAME(1)) >= 48 && double(FNAME(1)) <= 57)	% Numbers. See, char(48:57)
		% Example 20100109005439 will became 2010.0090379, 2010 (year) 009 (day of year) 0.0379 (decimal part of DOY) 
		dn = datenummx( sscanf(FNAME(1:4),'%f'), sscanf(FNAME(5:6),'%f'), sscanf(FNAME(7:8),'%f'), ...
			sscanf(FNAME(9:10),'%f'), sscanf(FNAME(11:12),'%f'), sscanf(FNAME(13:14),'%f') ) ...
			- datenummx( sscanf(FNAME(1:4),'%f'), 0, 0);
		td = sprintf('%f',(dn - fix(dn)));		% Decimal part of the DateNum as a string
		t = sprintf('%s.%03d%s', FNAME(1:4), fix(dn), td(3:6));
	else
		%t = sprintf('%s.%f',FNAME(2:5),sscanf(FNAME(6:8),'%f') + sscanf(FNAME(9:10),'%f')/24);
		t = sprintf('%f',sscanf(FNAME(6:8),'%f') + sscanf(FNAME(9:10),'%f')/24);	% Decimal 'Julian day'
	end

% -----------------------------------------------------------------------------------------
function radio_conv2netcdf_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_multiBand handles.radio_conv2vtk],'Val',0)

% -----------------------------------------------------------------------------------------
function radio_conv2vtk_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_multiBand handles.radio_conv2netcdf],'Val',0)

% -----------------------------------------------------------------------------------------
function radio_multiBand_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_conv2netcdf handles.radio_conv2vtk],'Val',0)

% -----------------------------------------------------------------------------------------
function check_L2_CB(hObject, handles)
% To deal with MODIS L2 products
	if (get(hObject,'Val'))
		set(handles.check_L2conf ,'Vis','on')			% Make the config file option visible
		set(handles.check_region,'Val',1)				% Do this for the user 
		check_region_CB(handles.check_region, handles)	% Simulate that it has been checked
	else
		set(handles.check_L2conf ,'Vis','off')
	end

% -----------------------------------------------------------------------------------------
function check_L2conf_CB(hObject, handles)
% Scan the the L2config.txt file for the -R string and fill the region boxes with it
	if (get(hObject,'Val'))
		fake_R = 'nikles';
		opt_R = sniff_in_OPTcontrol(fake_R);
		if (strcmp(opt_R, fake_R))
			errordlg('The data/L2config.txt file is corrupted and non usable.','Error')
		else
			ind = strfind(opt_R,'/');
			try
				w = opt_R(3:ind(1)-1);				e = opt_R(ind(1)+1:ind(2)-1);
				s = opt_R(ind(2)+1:ind(3)-1);		n = opt_R(ind(3)+1:end);
				set(handles.edit_west, 'Str', w);	set(handles.edit_east, 'Str', e);
				set(handles.edit_south, 'Str', s);	set(handles.edit_north, 'Str', n);
			catch
				errordlg('The -R string in data/L2config.txt file is screwed and useless.','Error')
				set(handles.edit_west, 'Str', '');	set(handles.edit_east, 'Str', '');
				set(handles.edit_south, 'Str', '');	set(handles.edit_north, 'Str', '');
			end
		end
	end

% -----------------------------------------------------------------------------------------
function check_region_CB(hObject, handles)
	if (get(hObject,'Val'))
		set([handles.edit_north handles.edit_south handles.edit_west handles.edit_east],'Enable','on')
	else
		set([handles.edit_north handles.edit_south handles.edit_west handles.edit_east],'Enable','off')
	end

% -----------------------------------------------------------------------------------------
function edit_north_CB(hObject, handles)
	x1 = str2double(get(hObject,'String'));
	if (isnan(x1)),		set(hObject,'String',''),	return,		end
	x2 = get(handles.edit_south,'String');
	if (~isempty(x2))
		x2 = str2double(x2);
		if (x2 >= x1) 
			errordlg('North Latitude <= South Latitude','Error in Latitude limits')
			set(hObject,'String','')
		end
	end

% -----------------------------------------------------------------------------------------
function edit_south_CB(hObject, handles)
	x1 = str2double(get(hObject,'String'));
	if (isnan(x1)),		set(hObject,'String',''),	return,		end
	x2 = get(handles.edit_north,'String');
	if (~isempty(x2))
		x2 = str2double(x2);
		if (x2 <= x1) 
			errordlg('South Latitude >= North Latitude','Error in Latitude limits')
			set(hObject,'String','')
		end
	end

% -----------------------------------------------------------------------------------------
function edit_west_CB(hObject, handles)
	x1 = str2double(get(hObject,'String'));
	if (isnan(x1)),		set(hObject,'String',''),	return,		end
	x2 = get(handles.edit_east,'String');
	if (~isempty(x2))
		x2 = str2double(x2);
		if (x2 <= x1) 
			errordlg('East Longitude <= West Longitude','Error in Longitude limits')
			set(hObject,'String','')
		end
	end

% -----------------------------------------------------------------------------------------
function edit_east_CB(hObject, handles)
	x1 = str2double(get(hObject,'String'));
	if (isnan(x1)),		set(hObject,'String',''),	return,		end
	x2 = get(handles.edit_west,'String');
	if (~isempty(x2))
		x2 = str2double(x2);
		if (x2 >= x1) 
			errordlg('East Longitude <= West Longitude','Error in Longitude limits')
			set(hObject,'String','')
		end
	end

% -----------------------------------------------------------------------------------------
function edit_stripeWidth_CB(hObject, handles)
	x1 = str2double(get(hObject,'String'));
	if (isnan(x1)),		set(hObject,'String','0.5'),	end

% -----------------------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
% Test for obvious errors and start computation

	if (~isempty(handles.OneByOneNameList) && handles.OneByOneFirst)	% Files we entered one by one. Must trick to reuse code
		lixoName = [handles.path_tmp 'listName_lixo.txt'];
		fid = fopen(lixoName, 'w');
		for (k = 1:numel(handles.OneByOneNameList))
			fprintf(fid, '%s\n', handles.OneByOneNameList{k});	% Bloody thing doesn't let write all at once
		end
		fclose(fid);
		% Now call push_namesList_CB as if we had a file with the names list. Clever me, no?
		push_namesList_CB([], handles, lixoName);
		builtin('delete', lixoName);
		handles = guidata(handles.figure1);		% Get the updated version
		handles.OneByOneFirst = false;			% To get out of the otherwise dead-end logic
		guidata(handles.figure1, handles)
	end

	if (isempty(handles.nameList))
		errordlg('No files to work on. You either didn''t give them or all names are wrong.','Error')
		return
	end

	got_R = false;		west = [];			east = [];		south = [];		north = [];
	if (get(handles.check_region, 'Val'))
		north = str2double(get(handles.edit_north,'String')); 
		west = str2double(get(handles.edit_west,'String')); 
		east = str2double(get(handles.edit_east,'String')); 
		south = str2double(get(handles.edit_south,'String')); 
		if ( any(isnan([west east south north])) )
			errordlg('One or more of the region limits was not provided','Error'),	return
		end
		got_R = true;
		% See if we have an L2 request and if user wants to use our default settings for that
		if (get(handles.check_L2, 'Val') && ~get(handles.check_L2conf, 'Val'))
			if (get(handles.radio_multiBand, 'Val'))
				errordlg('No, creating multi-band file is not possible with L2 MODIS files','Error'),	return
			end
			fid = fopen([handles.path_tmp 'L2config.txt'], 'w');
			fprintf(fid,'# Config file created by Empilhador with default settings for L2 MODIS files\n\n');
			fprintf(fid,'# Variables used in converting L2 satellite images from sensor to geographical coordinates.\n');
			fprintf(fid,'# First one contains parameters of interpolation (-C is gmtmbgrid only) and second are quality flags\n');
			fprintf(fid,'MIR_EMPILHADOR -I0.01 -C3 -R%.12g/%.12g/%.12g/%.12g\n', west,east,south,north);
			fprintf(fid,['#MIR_EMPILHADOR_F  ATMFAIL,LAND,HIGLINT,HILT,HISATZEN,STRAYLIGHT,CLDICE,' ...
				'COCCOLITH,HISOLZEN,LOWLW,CHLFAIL,NAVWARN,MAXAERITER,CHLWARN,ATMWARN,NAVFAIL,FILTER\n']);
			fclose(fid);
		end
	end

	cut2cdf(handles, got_R, west, east, south, north)

% -----------------------------------------------------------------------------------------
function cut2tif(handles, got_R, west, east, south, north, FileName)
% Save into a multi-band GeoTIFF file

	[pato, fname, EXT] = fileparts(FileName);
	if (isempty(EXT)),		FileName = [FileName '.tiff'];	end
	fname = FileName;

	att = gdalread(handles.nameList{1}, '-M');

	opt_R = ' ';		head = att.GMT_hdr;
	% If user wants a sub-region
	if (got_R)			% We must give the region in pixels since the image is trully not georeferenced (comment for nasa HDF)
		cp = round(([west east] - head(1)) / head(8));
		rp = round(([south north] - head(3)) / head(9));
		if (cp(1) < 0 || cp(2) > att.RasterXSize)		% Almost sure it should be >=
			msg = 'Sub-region West/Est is outside that grid''s limits';
			errordlg(msg, 'ERROR'),		error(msg)
		end
		if (rp(1) < 0 || rp(2) > att.RasterYSize)		% Almost sure it should be >=
			msg = 'Sub-region South/North is outside that grid''s limits';
			errordlg(msg, 'ERROR'),		error(msg)
		end
		head(1) = head(1) + cp(1)*head(8);		head(2) = head(1) + cp(2)*head(8);
		head(3) = head(3) + rp(1)*head(9);		head(4) = head(3) + rp(2)*head(9);
		rows = att.RasterYSize;
		rp = rows - rp -1;		rp = [rp(2) rp(1)];
		opt_R = sprintf('-r%d/%d/%d/%d',cp(1:2),rp(1:2));
	end

    nSlices = numel(handles.nameList);
    img = gdalread(handles.nameList{1}, opt_R);
	n_row = size(img, 1);
	n_col = size(img, 2);
	for (k = 2:nSlices)
		set(handles.listbox_list,'Val',k),		pause(0.01)			% Show advance
    	Z = gdalread(handles.nameList{k}, opt_R);
		ny = size(Z, 1);	nx = size(Z, 2);
		if ((nx ~= n_col) || (ny ~= n_row))
			errordlg('This image has not the same size as precedentes.','ERROR'),	return
		end
		n_col = nx;		n_row = ny;
		img = cat(3, img, Z);
	end
	clear Z
	
	hdr.name = fname;		hdr.driver = 'GTiff';
	hdr.projWKT = att.ProjectionRef;
	hdr.Xinc = head(8);		hdr.Yinc = head(9);
	hdr.ULx = head(1);		hdr.ULy = head(4);
	if (~isempty(att.GCPvalues)),	hdr.gcp = att.GCPvalues;	end
	
	gdalwrite(img,hdr)
	set(handles.listbox_list,'Val',1)

% -----------------------------------------------------------------------------------------
function cut2cdf(handles, got_R, west, east, south, north)
% Saveas as a multi-layer netCDF file

	if (get(handles.radio_conv2netcdf,'Val'))		% netCDF format
		this_ext = '.nc';							txt0 = '*.nc;*.grd';
		txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
	elseif (get(handles.radio_conv2vtk,'Val'))		% VTK
		this_ext = '.vtk';							txt0 = '*.vtk';
		txt1 = 'VTK format (*.vtk)';				txt2 = 'Select output VRT file';
	else											% Multi-band
		txt0 = '*.tiff;*.tif';
		txt1 = '(Geo)Tiff format (*.tiff)';			txt2 = 'Select output Tiff file';
	end
	if (isempty(handles.outname))
		[FileName,PathName] = put_or_get_file(handles,{txt0,txt1; '*.*', 'All Files (*.*)'},txt2,'put');
		if isequal(FileName,0),		return,		end
	else
		[PathName, FileName, EXT] = fileparts(handles.outname);
		if (~isempty(EXT)),		FileName = [FileName EXT];	end		% Don't loose the eventual extension
		if (isempty(PathName))					% If it's empty we use the path of the data files
			PathName = [fileparts(handles.nameList{1}) filesep];
		else
			PathName = [PathName filesep];
		end
	end
	[pato, fname, EXT] = fileparts(FileName);
	if (isempty(EXT)),		FileName = [fname this_ext];	end
	grd_out = [PathName FileName];

	if (get(handles.radio_multiBand,'Val'))		% Multi-band. We now pass the hand to its own function
		cut2tif(handles, got_R, west, east, south, north, grd_out)
		return
	end

	% Read relevant metadata. Attention, if we have a subdataset request ATT holds the attribs of the SDS
	[head, opt_R, slope, intercept, base, is_modis, is_linear, is_log, att, do_SDS] = ...
		get_headerInfo(handles, handles.nameList{1}, got_R, west, east, south, north);
	handles = guidata(handles.figure1);			% The get_att() function may have changed handles

	handles.geog = 1;			handles.head = head;
	handles.was_int16 = 0;		handles.computed_grid = 0;

	if (get(handles.radio_conv2vtk,'Val')),		fid = write_vtk(handles, grd_out, 'hdr');	end

	% The attributes
	misc = struct('x_units',[],'y_units',[],'z_units',[],'z_name',[],'desc',[], ...
		'title',[],'history',[],'srsWKT',[], 'strPROJ4',[]);
	if (~isempty(handles.desc_attrib))
		misc.desc = handles.desc_attrib;
	end

	nSlices = numel(handles.nameList);
	n_rows = 0;		n_cols = 0;
	n_cd = 1;
	empties = cell(nSlices, 1);			% To hold the names of the empty/failures files/SDSs
	have_empties = false;

	nL = 1;			% Counter to the effective number of non-empty layers
	for (k = 1:nSlices)
		set(handles.listbox_list,'Val',k),		pause(0.01)			% Show advance

		if (handles.caracol(k))						% Ai, we need to change CD
			msg = handles.changeCD_msg{n_cd};
			resp = yes_or_no('string',['OK, this one is over. ' msg '  ... and Click "Yes" to continue']);
			n_cd = n_cd + 1;
			if (strcmp(resp, 'No')),	return
			else						continue
			end
		end

		try
			if (do_SDS && k > 1)		% If we have an SDS request, get the attribs of that SDS (needed in getZ)
				[att, do_SDS] = get_headerInfo(handles, handles.nameList{k}, got_R, west, east, south, north);
				handles = guidata(handles.figure1);		% The get_att() function may have changed handles
			end
		catch
			str = sprintf('Error reading header of file: %s\nIgnoring it\n', handles.nameList{k});
			fprintf(str)
			warndlg(str,'WarnError')
			empties{k} = handles.nameList{k};	have_empties = true;
			continue
		end

		curr_fname = handles.nameList{k};
		if (isfield(handles, 'uncomp_name') && ~isempty(handles.uncomp_name))
			curr_fname = handles.uncomp_name;		% File was compressed, so we have to use the uncompressed version
		end

		try
			% In the following, if any of slope, intercept or base changes from file to file ... f
			NoDataValue = att.Band(1).NoDataValue;		% Backup it because it might be changed for other (now unknow) reasons.
			[Z, handles.have_nans, att, was_empty_name] = ...
				getZ(curr_fname, att, is_modis, is_linear, is_log, slope, intercept, base, opt_R, handles);
			att.Band(1).NoDataValue = NoDataValue;
			if (~isempty(was_empty_name))
				empties{k} = was_empty_name;	have_empties = true;
			end
		catch
			str = sprintf('Error reading file: %s\nIgnoring it\n', curr_fname);
			str = strrep(str, '\','/');		% We need to escape these '\' to not bother fprintf()
			str = sprintf('%s\n%s\n\n', str, lasterr);
			fprintf(str)
			warndlg(str,'WarnError')
			Z = alloc_mex(n_rows, n_cols, 'single', NaN);
			empties{k} = curr_fname;	have_empties = true;
		end

		% Check if all grids have the same size
		if (k == 1)
			n_rows = size(Z,1);		n_cols = size(Z,2);
		elseif (~isequal([n_rows n_cols],[size(Z,1) size(Z,2)]))
			warndlg(['The grid ' curr_fname ' has different size than precedents. Jumping it.'],'Warning')
			if (~isempty(handles.uncomp_name)),		try		delete(handles.uncomp_name);	end,	end
			continue
		end

		if (isfield(att, 'hdrModisL2_NEED2READ'))			% Grid was likely reinterpolated. Update header
			handles.head = att.GMT_hdr;
		end

		if (get(handles.radio_conv2vtk,'Val'))				% Write this layer of the VTK file and continue
			if (isempty(empties{k}))
				write_vtk(fid, grd_out, Z);
			end
			if (~isempty(handles.uncomp_name)),		try		delete(handles.uncomp_name);	end,	end
			continue
		end
		
		if (isa(Z,'int8') && (min(Z(:)) >= 0))
			grdutils(Z,'-c');								% Shift by -128 so it goes well with the uint8 add_off elsewere
		end

		if (isa(Z,'single'))
			zz = grdutils(Z,'-L');		handles.head(5:6) = [zz(1) zz(2)];
		else			% min/max is bugged when NaNs in singles
			handles.head(5:6) = [double(min(Z(:))) double(max(Z(:)))];
		end

		% Must treat compiled version differently since, biggest of mysteries, nc_funs
		% than hangs when writing unlimited variables
		if (true || ~handles.IamCompiled)			% Let's first short circuit this before total remove of other branch

			% For now we let this case reach here, but in future we should make this test/decision right after getZ()
			% The problem of deleting the uncompressed file must be solved too
			if (isempty(empties{k}))
				if (~isempty(handles.strTimes)),		t_val = handles.strTimes{k};
				else									t_val = sprintf('%d',k - 1);
				end
				if (nL == 1)
					nc_io(grd_out, ['w-' t_val '/time'], handles, reshape(Z,[1 size(Z)]), misc)
				else
					kk = nL - n_cd;			% = k - 1 - (n_cd - 1)	We need this when we had "@ change CD" messages
					nc_io(grd_out, sprintf('w%d\\%s', kk, t_val), handles, Z)
				end
				nL = nL + 1;		% Counts effective number of non-empty layers
			end
		else
			if (k == 1)
				handles.levelVec = str2double(handles.strTimes);
     			nc_io(grd_out,sprintf('w%d/time',nSlices), handles, reshape(Z,[1 size(Z)]), misc)
			else
				kk = k - n_cd;			% = k - 1 - (n_cd - 1)	We need this when we had "@ change CD" messages
				nc_io(grd_out, sprintf('w%d', kk), handles, Z)
			end
		end

		if (isfield(handles, 'uncomp_name') && ~isempty(handles.uncomp_name))
			try		delete(handles.uncomp_name);	end
		end
	end

	if (get(handles.radio_conv2vtk,'Val')),		fclose(fid);	end
	set(handles.listbox_list,'Val',1)

	if (have_empties)			% Right. Get the names of empty arrays
		c = false(nSlices,1);	% c vector with true for the files where we have data
		for (k = 1:nSlices)
			if (isempty(empties{k})),	c(k) = true;	end
		end
		empties(c) = [];		% The remainings of this are the failures/empties
		message_win('create',empties, 'figname','Names of no-data files', 'edit','yes')
	end

% -----------------------------------------------------------------------------------------
function [head, opt_R, slope, intercept, base, is_modis, is_linear, is_log, att, do_SDS] = ...
			get_headerInfo(handles, name, got_R, west, east, south, north)
% Get several direct and indirect (computed) informations about the file NAME or one of its subdatasets.
% The [att, do_SDS] = get_headerInfo(...) form is also supported and used when processing L2 files.

	[att, do_SDS, uncomp_name] = get_att(handles, name);

	% GDAL wrongly reports the corners as [0 nx] [0 ny] when no SRS
	if ( isequal((att.Corners.LR - att.Corners.UL), [att.RasterXSize att.RasterYSize]) && ~all(att.Corners.UL) )
		att.GMT_hdr(1:4) = [1 att.RasterXSize 1 att.RasterYSize];
	end

	% This case needs it
	if (~isempty(uncomp_name)),		att.fname = uncomp_name;
	else							att.fname = name;
	end
	[head , slope, intercept, base, is_modis, is_linear, is_log, att, opt_R] = ...
		getFromMETA(att, got_R, handles, west, east, south, north);
	
	if (nargout <= 2)			% Short form
		head = att;
		if (nargout == 2),		opt_R = do_SDS;		end
	end

% -----------------------------------------------------------------------------------------
function [att, indSDS, uncomp_name] = get_att(handles, name)
% Get the attributes of the root file or, in case we have one, of the requested subdataset
% This is a 'private' function that is called only by get_headerInfo()

	indSDS = 0;
	[att, uncomp_name] = get_baseNameAttribs(name);		% It calls deal_with_compressed()

	if ( att.RasterCount == 0 && ~isempty(att.Subdatasets) )	
		indSDS = 1;
		if (~isempty(handles.SDSinfo))
			if (~isnumeric(handles.SDSinfo))			% The SDS info is in its name form. Must convert to number
				ind = find_in_subdatasets(att.Subdatasets, handles.SDSinfo{1}(2:end));	% Remember, first char is '-'
				if (~ind)
					errordlg('The provided name of the Subdataset does not exist in file. Bye.','Error')
					error('The provided name of the Subdataset does not exist in file.')
				end
				handles.SDSthis = (ind + 1) / 2;		% Do this calc because we need here the SDS number from top of file
			end
			ind = strfind(att.Subdatasets{handles.SDSthis * 2 - 1}, '=');
			indSDS = handles.SDSthis * 2 - 1;
		elseif (strncmp(att.DriverShortName, 'HDF4', 4))	% Some MODIS files
			ind = strfind(att.Subdatasets{1}, '=');
		else
			errordlg('File has Sub-Datasets but you told me nothing about it.','ERROR')
			error('File has Sub-Datasets but you told me nothing about it.')
		end
		AllSubdatasets = att.Subdatasets;				% Copy this for keeping it as a subdataset field too
		FileName = att.Subdatasets{indSDS}(ind+1:end);	% First "ind" chars are of the form SUBDATASET_1_NAME=
		att = gdalread(FileName,'-M','-C');				% Try again
		att.AllSubdatasets = AllSubdatasets;			% A non-standard that is also in some cases set in Mirone
	end

	if (~isempty(uncomp_name))
		handles.uncomp_name = uncomp_name;				% To know whether to use the uncompressed name
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function [head , slope, intercept, base, is_modis, is_linear, is_log, att, opt_R] = ...
	getFromMETA(att, got_R, handles, west, east, south, north)
% Get complementary data from the att struct. This is mostly a helper function to get_headerInfo()
% but it is detached from it because in this way it can be called by exterior code. Namelly by Mirone.
% The cases addressed here are some of the ones raised by HDF files.
% WARNING: the ATT struct must have and extra field att.fname (to be eventualy used by hdfread)
% NOTE1: When called from outside (e.g Mirone) use only the [...] = getFromMETA(att) form
% NOTE2: For HDF files the ATT struct will be added the field 'hdrInfo' (returned when hdrfinfo)

	if (nargin == 1),	got_R = false;		end

	opt_R = ' ';	is_modis = false;		is_linear = false;		is_log = false;
	slope = 1;		intercept = 0;			base = 1;
	modis_or_seawifs = false;				is_HDFEOS = false;		is_ESA = false;
	att.hdrInfo = [];
	head = att.GMT_hdr;

	if ( ~isempty(att.Metadata) && ~isempty(search_scaleOffset(att.Metadata, 'HDFEOSVersion')) )
		is_HDFEOS = true;
		if ( isnan(search_scaleOffset(att.Metadata, 'ENVISAT')) )	% Poor trick to find ESA (well, ENVISAT) products
			is_ESA = true;
		end
	end

	if ( ~is_HDFEOS && ~isempty(att.Metadata) && ...
			~isempty(search_scaleOffset(att.Metadata, 'MODIS')) || ~isempty(search_scaleOffset(att.Metadata, 'SeaWiFS')) )
		modis_or_seawifs = true;
	end

	if ( modis_or_seawifs && strncmp(att.DriverShortName, 'HDF4', 4) && ~isempty(search_scaleOffset(att.Metadata, 'Level-2')) )
		out = search_scaleOffset(att.Metadata, 'slope');
		if (~isempty(out))		% Otherwise, no need to search for a 'intercept'
			slope = out;
			if (slope ~= 1)
				out = search_scaleOffset(att.Metadata, 'intercept');
				if (~isempty(out)),		intercept = out;	end
				is_linear = true;
			end
		end
		head(1:4) = [1 att.RasterXSize 1 att.RasterYSize];
		att.Band(1).NoDataValue = -32767;	% Shity format doesn't declare this null part (good for SST)
		if (isfield(att, 'subDsName'))		% Known values (found by file inspection)
			NoDataValue = guess_nodataval(att.subDsName);
			if (~isempty(NoDataValue)),		att.Band(1).NoDataValue = NoDataValue;		end
		end
		is_modis = true;					% We'll use this knowledge to 'avoid' Land pixels = -32767  
		% Compressed HDF crash R13 hdf_funs() and we probably only need it in rare ocasions.
		% So signal it to the extreme cases where we still have to call it to do it then.
		% Perhaps it will still crash, but at least it won't on all other cases where we don't need hdf_funs().
		att.hdrModisL2_NEED2READ = att.fname;
		if (got_R)			% For L2 files we cannot find -R here (before sensor to geog coords conversion)
							% So we store the croping info in 'att' to use later after reinterpolation
			att.crop_info.opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g',west,east,south,north);
			att.crop_info.limits = [west east south north];
			got_R = false;	% So that last block in this function won't try to execute.
		end
	elseif ( modis_or_seawifs && strncmp(att.DriverShortName, 'HDF4', 4) )
		x_max = search_scaleOffset(att.Metadata, 'Easternmost');	% Easternmost Latitude=180
		x_min = search_scaleOffset(att.Metadata, 'Westernmost');	% Westernmost Latitude=-180
		y_max = search_scaleOffset(att.Metadata, 'Northernmost');	% Northernmost Latitude=90
		y_min = search_scaleOffset(att.Metadata, 'Southernmost');	% Southernmost Latitude=-90
		dx = (x_max - x_min) / att.RasterXSize;
		dy = dx;
		x_min = x_min + dx/2;		x_max = x_max - dx/2;	% Orig data was pixel registered
		y_min = y_min + dx/2;		y_max = y_max - dx/2;
		head(1:4) = [x_min x_max y_min y_max];
		head(8:9) = dx;
		att.Corners.UL = [x_min y_max];			
		att.Corners.LR = [x_max y_min];			
		if (got_R)			% We must give the region in pixels since the image is trully not georeferenced
			rows = att.RasterYSize;
		end
		head(7) = 0;		% Make sure that grid reg is used

		% Get the the scaling equation and its parameters
		if ( ~isempty(search_scaleOffset(att.Metadata, 'linear')) )
			slope = search_scaleOffset(att.Metadata, 'Slope');	% att.Metadata{47} -> Slope=0.000717185
			intercept = search_scaleOffset(att.Metadata, 'Intercept');
			is_linear = true;
		elseif ( ~isempty(search_scaleOffset(att.Metadata, 'logarithmic')) )
			base = search_scaleOffset(att.Metadata, 'Base');	% att.Metadata{41} -> Base=10
			slope = search_scaleOffset(att.Metadata, 'Slope');	% att.Metadata{41} -> Slope=0.000717185
			intercept = search_scaleOffset(att.Metadata, 'Intercept');
			is_log = true;
		end

		att.Band(1).NoDataValue = 65535;						% Shity format doesn't declare this.
		nv = search_scaleOffset(att.Metadata, 'Fill');			% But sometimes (some SeaWifs) it exists
		if (~isempty(nv)),	att.Band(1).NoDataValue = nv;	end
		if (head(5) == att.Band(1).NoDataValue),	head(5) = NaN;	end		% Force later recomputing of array min/max
		att.GMT_hdr = head;			% We need this updated
		is_modis = true;			% We'll use this knowledge to 'avoid' Land pixels = 65535

	elseif ( ~is_HDFEOS && ~modis_or_seawifs && strncmp(att.DriverShortName, 'HDF4', 4)  )		% TEMP -> SST PATHFINDER
		finfo = hdf_funs('hdfinfo', att.fname);
		if (strcmpi(finfo.SDS.Attributes(11).Name, 'slope'))
			slope = double(finfo.SDS.Attributes(11).Value);		% = 0.075;
			intercept = double(finfo.SDS.Attributes(12).Value);	% = -3.0;
		else
			out = search_scaleOffset(finfo.SDS.Attributes, 'slope');
			if (~isempty(out))		% Otherwise, no need to search for a 'intercept'
				slope = out;
				out = search_scaleOffset(finfo.SDS.Attributes, 'intercept');
				if (~isempty(out)),		intercept = out;	end
			end
			if (slope == 1)			% We may have a netCDF style naming. Check
				out = search_scaleOffset(finfo.SDS.Attributes, 'scale_factor');
				if (~isempty(out))
					slope = out;
					out = search_scaleOffset(finfo.SDS.Attributes, 'add_off', 7);	% Sometimes is ADD_OFF, others ADD_OFFSET and no-one is killed for that
					if (~isempty(out)),		intercept = out;	end
				end
			end
		end
		lat = finfo.SDS.Dims(1).Scale;		% Get the latitudes
		lon = finfo.SDS.Dims(2).Scale;		% Get the longitudes
		if (isnumeric(lat))
			x_min = lon(1);			x_max = lon(end);
			y_min = min(lat(1),lat(end));
			y_max = max(lat(1),lat(end));
			att.GMT_hdr(1:4) = [x_min x_max y_min y_max];	% We need this updated
			att.Corners.UL = [x_min y_max];			
			att.Corners.LR = [x_max y_min];			
			dx = lon(3) - lon(2);	dy = abs(lat(3) - lat(2));
			att.GMT_hdr(7:9) = [0 dx dy];
		else				% If not, use array size as coordinates
			x_min = 1;		x_max = finfo.SDS.Dims(2).Size;
			y_min = 1;		y_max = finfo.SDS.Dims(1).Size;
			dx = 1;			dy = 1;
		end
		head(1:4) = [x_min x_max y_min y_max];
		head(8:9) = [dx dy];
		if (got_R)			% We must give the region in pixels since the image is trully not georeferenced
			rows = finfo.SDS.Dims(1).Size;					% Number of rows
		end
		head(7) = 0;		% Make sure that grid reg is used
		is_linear = true;
		att.hdrInfo = finfo;

	elseif ( ~is_HDFEOS && ~modis_or_seawifs && strcmp(att.DriverShortName, 'netCDF') )		% GHRSST -> PATHFINDER
		if (got_R)
			rows = att.RasterYSize;
			x_min = head(1);	y_min = head(3);
			dx = head(8);		dy = head(9);
		end

	elseif ( is_HDFEOS && ~is_ESA)		% This case might not be complete as yet.
		x_min = search_scaleOffset(att.Metadata, 'WESTBOUNDINGCOORDINATE');
		x_max = search_scaleOffset(att.Metadata, 'EASTBOUNDINGCOORDINATE');
		y_min = search_scaleOffset(att.Metadata, 'SOUTHBOUNDINGCOORDINATE');
		y_max = search_scaleOffset(att.Metadata, 'NORTHBOUNDINGCOORDINATE');
		rows = search_scaleOffset(att.Metadata, 'DATAROWS');
		cols = search_scaleOffset(att.Metadata, 'DATACOLUMNS');
		dx = (x_max - x_min) / cols;
		dy = (y_max - y_min) / rows;
		att.GMT_hdr(1:4) = [x_min x_max y_min y_max];	% We need this updated
		att.GMT_hdr(7:9) = [0 dx dy];
		head = att.GMT_hdr;
		att.Corners.UL = [x_min y_max];			
		att.Corners.LR = [x_max y_min];			
		is_linear = true;
	elseif ( is_HDFEOS && is_ESA)		% One more incredible messy HDF product. Nothing (spatialy) reliable inside.
		ind = strfind(att.fname, '_');
		tmp = att.fname(ind(end)+1:end);
		indH = strfind(tmp, 'H');		indV = strfind(tmp, 'V');		indDot = strfind(tmp, '.');
		if (isempty(indH) || isempty(indV))
			errordlg('Sorry but this is not a 5 degrees tile of the super non-documented ESA product. Don''t know how to proceed.','Error')
			return
		end
		% The following is crazy. ESA actually uses a grid registration schema but calls it pixel reg.
		% Hence the left side of each tile is aligned with a multiple of 5 but lacks the last col/row.
		x_min = sscanf(tmp(indH+1:indV-1), '%d') * 5 - 180;
		y_max = 90 - sscanf(tmp(indV+1:indDot-1), '%d') * 5;
		rows = att.RasterYSize;
		dx = 5 / att.RasterXSize;	dy = 5 / att.RasterYSize;
		x_max = x_min + 5 - dx;		y_min = y_max - 5 + dx;
		att.GMT_hdr(1:4) = [x_min x_max y_min y_max];	% We need this updated
		att.GMT_hdr(7:9) = [0 dx dy];
		head = att.GMT_hdr;
		att.Corners.UL = [x_min y_max];		att.Corners.LL = [x_min y_min];
		att.Corners.LR = [x_max y_min];		att.Corners.UR = [x_max y_max];
		att.DriverShortName = sprintf('ATENTION: Spatial info displayed here may NOT be\n reliable due to ESA awfull lack of info in file\n\n%s\n',att.DriverShortName);
		att.ProjectionRef = ogrproj('+proj=longlat +ellps=wgs84 +nodefs');
	else					% Other types
		if (got_R),		rows = att.RasterYSize;		end
		x_min = head(1);	y_min = head(3);
		dx = head(8);		dy = head(9);
	end

	% If user wants a sub-region
	if (got_R)			% We must give the region in pixels since the image is trully not georeferenced (comment for nasa HDF)
		cp = round(([west east] - x_min) / dx);
		rp = round(([south north] - y_min) / dx);
		if (cp(1) < 0 || cp(2) > att.RasterXSize)		% Almost sure it should be >=
			msg = 'Sub-region West/Est is outside that grid''s limits';
			errordlg(msg, 'ERROR'),		error(msg)
		end
		if (rp(1) < 0 || rp(2) > att.RasterYSize)		% Almost sure it should be >=
			msg = 'Sub-region South/North is outside that grid''s limits';
			errordlg(msg, 'ERROR'),		error(msg)
		end
		head(1) = x_min + cp(1)*dx;		head(2) = x_min + cp(2)*dx;
		head(3) = y_min + rp(1)*dy;		head(4) = y_min + rp(2)*dy;
		rp = rows - rp -1;		rp = [rp(2) rp(1)];
		opt_R = sprintf('-r%d/%d/%d/%d',cp(1:2),rp(1:2));
	end

% ----------------------------------------------------------------------------------------
function out = search_scaleOffset(attributes, what, N)
% Search for the WHAT attribute in ATTRIBUTES. If find return its VALUE.
% Used to search for slope/intercept or scale_factor/add_offset in HDF files
	out = [];
	if (isa(attributes, 'struct'))
		if (nargin == 2)						% Exact search for WHAT
			for (k = numel(attributes):-1:1)					% Start from the bottom because they are likely close to it 
				if ( strcmpi(attributes(k).Name, what) )
					out = double(attributes(k).Value);
					break
				end
			end
		else									% Search with a strncmp. Motivated by the uterly stupid play with ADD_OFF & ADD_OFFSET
			for (k = numel(attributes):-1:1)				% Start from the bottom because they are likely close to it 
				if ( strncmpi(attributes(k).Name, what, N) )
					out = double(attributes(k).Value);
					break
				end
			end
		end
	else					% Not tested but it must be a cell array (the att.Metadata)
		for (k = 1:numel(attributes))
			id = strfind(attributes{k}, what);
			if ( ~isempty(id) )
				id_eq = strfind(attributes{k},'=');			% Find the '=' sign
				if (numel(id_eq) > 1),	continue,	end		% Sometimes comments have also the '=' char
				out = str2double(attributes{k}(id_eq+1:end));
				break
			end
		end
	end

% -----------------------------------------------------------------------------------------
function NoDataValue = guess_nodataval(DsName)
% Make an educated guess of the no-data value of HDF4 (MODIS) files
	NoDataValue = [];
	if (strncmp(DsName, 'sst', 3)),			NoDataValue = -32767;
	elseif (strcmp(DsName, 'l3m_data')),	NoDataValue = 65535;		% L3 SST
	elseif (strcmp(DsName, 'chlor_a')),		NoDataValue = -1;
	elseif (strncmp(DsName, 'nLw_',4)),		NoDataValue = 0;			% ???
	elseif (strcmp(DsName, 'K_490')),		NoDataValue = -5000;
	elseif (strcmp(DsName, 'tau_869')),		NoDataValue = 0;
	elseif (strcmp(DsName, 'eps_78')),		NoDataValue = 0;
	elseif (strcmp(DsName, 'angstrom_531')),NoDataValue = -32767;
	end

% -----------------------------------------------------------------------------------------
function [Z, have_nans, att, was_empty_name] = getZ(fname, att, is_modis, is_linear, is_log, slope, intercept, base, opt_R, handles)
% ATT may be still unknown (empty). In that case it will be returned by read_gdal()
% HANDLES, is transmitted only within internal calls of this function (that is, not from outside calls)

	uncomp_name = [];			IamInteractive = true;
	if (nargin < 9),	opt_R = ' ';	end
	if (nargin == 10 && ~handles.Interactive),	IamInteractive = false;		end		% External calls may be interactive

	if (nargin == 10 && ~isempty(handles.SDSinfo))		% We have a Sub-Dataset request
		clear_att = false;
		if (~isempty(att) && ~isempty(att.Subdatasets))
			ind = strfind(att.Subdatasets{handles.SDSthis * 2 - 1}, '=');
			fname = att.Subdatasets{handles.SDSthis * 2 - 1}(ind+1:end);	% First "ind" chars are of the form SUBDATASET_?_NAME=
		elseif (~isempty(att) && isempty(att.Subdatasets))
			ind = strfind(att.AllSubdatasets{handles.SDSthis * 2 - 1}, '=');
			fname = att.AllSubdatasets{handles.SDSthis * 2 - 1}(ind+1:end);
		else		% isempty(att) = true
			[fname, uncomp_name] = deal_with_compressed(fname);		% MUST GET RID OF THIS (read compressed directly)
			att = gdalread(fname, '-M');	clear_att = true;
			if (~isempty(uncomp_name)),		uncomp_name = fname;	end
		end
		if (clear_att),		att = [];	end
		IamCompiled = handles.IamCompiled;
	else
		% Need to know if "IamCompiled". Since that info is in handles, we need to find it out here
		if (nargin == 10 && isfield(handles, 'IamCompiled'))
			IamCompiled = handles.IamCompiled;
		else
			try			which('mirone');			IamCompiled = false;
			catch,		IamCompiled = true;
			end
		end
	end

	saveNoData = false;
	if (is_modis && ~isempty(att))
		NoDataValue = att.Band(1).NoDataValue;		% Save this value to reset it after read_gdal()
		saveNoData = true;
	end

	GMT_hdr = [];
	if (~isempty(att))					% Make copies of these
		GMT_hdr = att.GMT_hdr;
		Corners.LR = att.Corners.LR;
		Corners.UL = att.Corners.UL;
	end

	[Z, att, known_coords, have_nans, was_empty_name] = read_gdal(fname, att, IamCompiled, IamInteractive, '-C', opt_R, '-U');

	if (nargin == 10 && ~isempty(handles.SDSflag) && ~isempty(was_empty_name))
		ind = strfind(fname, ':');
		flag_fname = [fname(1:ind(end)) handles.SDSflag];
		att_flag = gdalread(flag_fname, '-M');
		F = read_gdal(flag_fname, att_flag, IamCompiled, IamInteractive, '-C', opt_R, '-U');
		Z(F < handles.minQuality) = NaN;
	end

	% See if we knew the image coordinates but that knowledge is lost in new att
	if (~known_coords && ~isempty(GMT_hdr))		% For the moment we only know for sure for L2 georeferenced products
		if ((isequal(att.GMT_hdr(8:9), [1 1]) && ~isequal(GMT_hdr(8:9), [1 1])) || (diff(GMT_hdr(1:2)) > 359.9))
			att.GMT_hdr = GMT_hdr;				% Recover the header info
			att.Corners.LR = Corners.LR;
			att.Corners.UL = Corners.UL;
		end
	end

	if (~isempty(uncomp_name)),		delete(uncomp_name);		end		% Delete uncompressed file.

	if (is_modis && ~isempty(att.Band(1).NoDataValue) && saveNoData)	% att.Band... is isempty when all work has been done before
		att.Band(1).NoDataValue = NoDataValue;		% Shity format doesn't inform on the no-data.
	end

	if (isempty(have_nans))			% It's only non-empty when processing L2 products with reinterpolation
		[Z, have_nans, att] = sanitizeZ(Z, att, is_modis, is_linear, is_log, slope, intercept, base);
	end

% -----------------------------------------------------------------------------------------
function [Z, have_nans, att] = sanitizeZ(Z, att, is_modis, is_linear, is_log, slope, intercept, base)
% Take care of possible scaling/offset transformations.
% Have it separate in a function because when processing L2 products we need to apply this right
% before the interpolation (in read_gdal()). This is too early with regard to the normal work-flow
% where it is applied at the end of function getZ

	if (nargin == 2)	% Go again to getFromMETA but without changing again att
		[head, slope, intercept, base, is_modis, is_linear, is_log] = getFromMETA(att);
	end

	ind = [];
	if (~isempty(att.Band(1).NoDataValue) && (att.Band(1).NoDataValue == -9999))		% TEMP -> PATHFINDER
		if ( ~isempty(att.Metadata) && ~isempty(search_scaleOffset(att.Metadata, 'dsp_SubImageName=QUAL')) )
			% Quality flags files cannot be NaNified. 
			% However, they should NOT have a scaling equation either. If quality is [0 7] why scalling?
			is_linear = false;
			% Furthermore, we also need to recast the flag array into int8 (it was uint8) because netCDF doesn't know UINT8
			Z = int8(Z);
		else
			ind = (Z == 0);
		end
	elseif (~isempty(att.Band(1).NoDataValue) && ~isnan(att.Band(1).NoDataValue))
		ind = (Z == (att.Band(1).NoDataValue));
	elseif (isnan(att.Band(1).NoDataValue))		% The nodata is NaN, replace NaNs in Z by zero
		ind = isnan(Z);
	end

	if ( is_modis && (isa(Z, 'int8') || isa(Z, 'uint8')) )
		is_linear = false;		ind = [];
		if (isa(Z, 'uint8')),	Z = int8(Z);	end		% netCDF doesn't know UINT8
	end

	% See if we must apply a scaling equation
	if (is_linear && (slope ~= 1 || intercept ~= 0))
		if (~isa(Z,'single')),		Z = single(Z);		end
		cvlib_mex('CvtScale',Z, slope, intercept)
	elseif (is_log)
		Z = single(base .^ (double(Z) * slope + intercept));
	end
	have_nans = 0;
	if (~isempty(ind))
		if (~isa(Z,'single') && ~isa(Z,'double'))		% Otherwise NaNs would be converted to 0
			Z = single(Z);
		end
		Z(ind) = NaN;		have_nans = 1;
	elseif (isempty(att.Band(1).NoDataValue) && isa(Z,'single'))		% It comes from the interpolation
		att.Band(1).NoDataValue = NaN;
		have_nans = grdutils(Z, '-N');
	end

% -----------------------------------------------------------------------------------------
function [Z, att, known_coords, have_nans, was_empty_name] = read_gdal(full_name, att, IamCompiled, IamInteractive, varargin)
% Help function to gdalread that deals with cases when file is compressed.
% ATT is the GDALREAD returned attributes. If empty, we'll get it here
% VARARGIN will normally contain one or more of '-C', opt_R, '-U'
% WARNING: If exist(att.hdfInfo) than att.fname should exist as well (both non standard)
%
% KNOWN_COORDS	Is a logical that when true the informs the caller that we already know the coordinates for sure
%				and no attempt should be made to fish them from the matadata info.
%
% WAS_EMPTY_NAME Is a char string with the name the file (or SDS when an SDS is requested) IF the corresponding
%                array is empty. If the array is not empty, WAS_EMPTY_NAME itself is empty.
%                This applies only to L2 files that are interpolated. We use it to check suspicious empty frames.
%
% This function is called only by getZ()

	have_nans = [];			% Will only become ~[] when input is a L2 product to be interpolated and referenced
	NoDataValue = [];		% Some defaults
	known_coords = false;	% If we know for sure the coords (as for georefed L2 products) tell that to caller
	was_empty_name = '';	% Will hold the file/SDS name in case it's empty (L2 only)
	[full_name, uncomp_name] = deal_with_compressed(full_name);

	opt_e = '';
	if (IamCompiled),	opt_e = '-e';	end		% Use aguentabar.dll

	if (isempty(att))
		att = get_baseNameAttribs(full_name);
	end

	if (att.RasterCount == 0 && ~isempty(att.Subdatasets) && strncmp(att.DriverShortName, 'HDF4', 4))	% Some MODIS files
		sds_num = 1;
		handles = guidata(gcf);
		if (~isempty(handles.SDSinfo)),		sds_num = handles.SDSthis * 2 - 1;		end		% We have a SubDataset request
		ind = strfind(att.Subdatasets{sds_num}, '=');
		full_name = att.Subdatasets{sds_num}(ind+1:end);				% First "ind" chars are of the form SUBDATASET_1_NAME=
		ind = strfind(att.Subdatasets{sds_num+1}, ' ');					% Need to guess the no-data value
		subDsName = att.Subdatasets{sds_num+1}(ind(1)+1:ind(2)-1);		% Get dataset name
		NoDataValue = guess_nodataval(subDsName);
		att.subDsName = subDsName;										% A non-standard name
	end

	% att.hdrInfo and att.hdrModisL2_NEED2READ are not default fields of the ATT struct
	try		fname = att.fname;
	catch,	fname = [];
	end
	if (isfield(att, 'hdrInfo') && ~isempty(att.hdrInfo) && (strcmp(att.hdrInfo.SDS.Name,'sst')) && isempty(strfind([varargin{:}], '-C')) )
		% Only particular case dealt now. This is HORRIBLE. Patches, over patches, over patches.
		Z = hdf_funs('hdfread', att.fname, att.hdrInfo.SDS.Name, 'index', {[1 1],[1 1], [att.RasterYSize att.RasterXSize]});
		Z = flipud(Z);
	else
		opt_L = ' ';	GCPvalues = [];
		if (isfield(att, 'hdrModisL2_NEED2READ') && ~isempty(att.hdrModisL2_NEED2READ))
			% First check if lon, lat are of the same size of Z. If yes, than there is no need to reinterpolate
			% the location grids ... to their same positions. This has the further advantage that all readings are
			% done with GDAL and so will work also with the stand-alone version.
			lonID = find_in_subdatasets(att.AllSubdatasets, 'longitude');
			ind = strfind(att.AllSubdatasets{lonID},'=');			% Still must rip the 'SUBDATASET_XX_NAME='
			lon_full = gdalread(att.AllSubdatasets{lonID}(ind+1:end), '-L');
			if (isequal(size(lon_full), [att.RasterYSize att.RasterXSize]))
				latID = find_in_subdatasets(att.AllSubdatasets, 'latitude');
				ind = strfind(att.AllSubdatasets{latID},'=');
				lat_full = gdalread(att.AllSubdatasets{latID}(ind+1:end), '-L');

			else			% Bad luck. We need to go through the hdfread way.
				att.hdrModisL2 = hdf_funs('hdfinfo', att.hdrModisL2_NEED2READ);		% MAY CRASH R13 when hdf files are compressed.
				% These boys think they are very funy. Oh it's so cute to write the file from right-to-left !!!
				Vg_index = numel(att.hdrModisL2.Vgroup);	% The uncomprehensible MESS never ends. Assume last has things of interst
				lon = fliplr( hdf_funs('hdfread', att.fname, att.hdrModisL2.Vgroup(Vg_index).SDS(1).Name, 'index', {[],[], []}) );
				lat = fliplr( hdf_funs('hdfread', att.fname, att.hdrModisL2.Vgroup(Vg_index).SDS(2).Name, 'index', {[],[], []}) );
				cntl_pt_cols = hdf_funs('hdfread', att.fname, att.hdrModisL2.Vgroup(Vg_index).SDS(3).Name, 'index', {[],[], []});

				cntl_pt_cols = double(cntl_pt_cols);
				lat = double(lat);		lon = double(lon);
				lat_full = zeros(size(lon,1), cntl_pt_cols(end));
				lon_full = zeros(size(lon,1), cntl_pt_cols(end));
				cols_vec = 1:cntl_pt_cols(end);
				for (k = 1:size(lon,1))
					lon_full(k,:) = akimaspline(cntl_pt_cols, lon(k,:), cols_vec);
					lat_full(k,:) = akimaspline(cntl_pt_cols, lat(k,:), cols_vec);
				end
				clear lon lat cols_vec
			end
			x_min = double(min([lon_full(1) lon_full(1,end) lon_full(end,1) lon_full(end)]));	% double because R6.5
			x_max = double(max([lon_full(1) lon_full(1,end) lon_full(end,1) lon_full(end)]));
			y_min = double(min([lat_full(1) lat_full(1,end) lat_full(end,1) lat_full(end)]));
			y_max = double(max([lat_full(1) lat_full(1,end) lat_full(end,1) lat_full(end)]));
			if (any([x_min x_max y_min y_max] == -999))		% CAN WE BELIVE THIS???? BUT HAPPENS!!!!!!!!!!!
				% YES, L2 MODIS files can have -999 as coordinates. This is unbelievable but happens.
				% The remedy is to recompute limits and forget the f.. coordinates that will be trimmed by -R below
				indNotFckCoords = (lon_full ~= -999);
				x_min = double( min(min(lon_full(indNotFckCoords))) );
				x_max = double( max(max(lon_full(indNotFckCoords))) );
				y_min = double( min(min(lat_full(indNotFckCoords))) );
				y_max = double( max(max(lat_full(indNotFckCoords))) );
				clear indNotFckCoords
			end
			opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', x_min, x_max, y_min, y_max);

			AllSubdatasets = att.AllSubdatasets;		% Copy of all subdatsets names in this sub-dataset att
			full_name = att.Name;
			NoDataValue = att.Band(1).NoDataValue;		% Save this that has been ESTIMATED before
			if (strcmp(varargin{end}, '-U')),	varargin(end) = [];		end		% We also don't want to UpDown
			opt_L = '-L';
		end

		if (nargout == 2)
			[Z, att] = gdalread(full_name, varargin{:}, opt_L);		% This ATT may be of a subdataset
		else
			Z = gdalread(full_name, varargin{:}, opt_L);
		end
		if (strcmp(att.DriverShortName, 'netCDF'))					% GHRSST 2.0 PATHFINDER ??
			if (all(Z(:) == 0))		% When files are compressed GDAL screws and returns all zeros
				handles_tmp.IamCompiled = IamCompiled;	handles_tmp.grdMaxSize = 1e12;
				handles_tmp.ForceInsitu = false;
				Z = read_grid(handles_tmp, full_name, 'GMT');		% Read the grid with our own functions
			else
				[Z, did_scale, att] = handle_scaling(Z, att);		% See if we need to apply a scale/offset
			end
		end
		if (isempty(att.GCPvalues) && ~isempty(GCPvalues)),		att.GCPvalues = GCPvalues;		end
		if ((att.Band(1).NoDataValue == -1) && (min(Z(:)) == -32767))	% DIRTY PATCH to avoid previous bad nodata guessing
			att.Band(1).NoDataValue = -32767;
			if (~isempty(NoDataValue)),		NoDataValue = -32767;	end
		end
		if (~isempty(NoDataValue)),		att.Band(1).NoDataValue = NoDataValue;	end		% Recover ESTIMATED value

		% For MODIS L2 products, we may still have many things to do
		if (strcmp(opt_L,'-L'))		% We are using this as an indication that the file is MODIS (need clever solution)

			% Here we must sanitize Z in case we must apply a scale/offset transform and/or NaNifying NodataValues
			[Z, have_nans, att] = sanitizeZ(Z, att);
			if (have_nans),		NoDataValue = NaN;		end

			if (IamInteractive)
				what = l2_choices(AllSubdatasets);		% Call secondary GUI to select what to do next
			else
				what = struct('georeference',1,'nearneighbor',0,'mask',0,'coastRes',0,'quality','');	% sensor coords
				ID = find_in_subdatasets(AllSubdatasets, 'qual_sst', 8);	% Check if we have a quality flags array
				if (ID)
					ind = strfind(AllSubdatasets{ID}, '=');			% Yes we have. Use it if not overruled by info in L2config
					what.qualSDS = AllSubdatasets{ID}(ind+1:end);
					what.quality = 0;
				end
				if (isfield(att, 'crop_info'))			% We have a crop request
					opt_R = att.crop_info.opt_R;
				end
			end

			% Go check if -R or quality flags request exists in L2config.txt file
			[opt_R, opt_I, opt_C, bitflags, flagsID, despike] = sniff_in_OPTcontrol(opt_R, att);	% Output opt_R gets preference

			if (isempty(what))							% User killed the window, but it's too late to stop so pretend ...
				what =  struct('georeference',1,'nearneighbor',1,'mask',0,'coastRes',0,'quality','');	% sensor coords
			end
			if (isempty(bitflags) && ~isempty(what.quality) && what.quality < 2)	% We have a GUI quality request
				qual = gdalread(what.qualSDS, opt_L);
				if (isequal(size(qual),size(Z)))		% Because we don't want to apply 'qual' to either lon or lat arrays
					Z(qual > what.quality) = NoDataValue;
				end
				clear qual
			end

			if (despike)					% MODIS SST are horribly spiked every other 10 vertical positions in
				Z = clipMySpikes(Z);		% sensor coordinates. This functions signifficantly reduces that effect.
			end

			if (~isempty(what) && what.georeference)	% OK, let's interpolate it into a regular geog grid

				if (isempty(opt_I)),	opt_I = '-I0.01';	end
				if (isempty(opt_C)),	opt_C = '-C3';		end		% For gmtmbgrid only

				if (~isempty(bitflags))
					ind = strfind(att.AllSubdatasets{flagsID},'=');	% Still must rip the 'SUBDATASET_XX_NAME='
					Zf = gdalread(att.AllSubdatasets{flagsID}(ind+1:end), varargin{:}, opt_L);
					c = false(size(Zf));
					Zf = uint32(Zf);
					for (k = 1:numel(Zf))
						if (any(bitget(Zf(k),bitflags))),	c(k) = true;	end
					end
					clear Zf
					if (~isfield(att, 'subDsName') || ~strcmp(att.subDsName, 'l2_flags'))		% 'Normal' cases 
						Z(c) = [];		lon_full(c) = [];		lat_full(c) = [];	% Remove the flagged values
					else										% Make a grid of the flagged vales themselfs
						Z = c;
						opt_C = '-C1';	% Retain in cell values only
					end
					clear c
				end

				if (isnan(NoDataValue)),		ind = isnan(Z);
				elseif (~isempty(NoDataValue)),	ind = (Z == NoDataValue);
				else
					ind = false;	% Just to not error below
				end
				Z(ind) = [];		lon_full(ind) = [];		lat_full(ind) = [];

				if (isempty(Z))			% Instead of aborting we let it go with fully NaNified array
					was_empty_name = full_name;		% Send back to the caller the file (or sds) name of the empty thing
					[Z, head] = nearneighbor_m(single(1e3), single(1e3), single(0), opt_R, opt_e, '-N1', opt_I, '-S0.02');
				else
					if (what.nearneighbor)
						lon_full = single(lon_full);			lat_full = single(lat_full);	Z = single(Z);
						[Z, head] = nearneighbor_m(lon_full(:), lat_full(:), Z(:), opt_R, opt_e, '-N2', opt_I, '-S0.04');
					else
						if (~isa(lon_full,'double')),	lon_full = double(lon_full);	lat_full = double(lat_full);	end
						%[Z, head] = gmtmbgrid_m(lon_full(:), lat_full(:), double(Z(:)), opt_I, opt_R, '-Mz', opt_C);
						[Z, head, was_empty] = smart_grid(lon_full(:), lat_full(:), double(Z(:)), opt_I, opt_R, opt_C);
						if (was_empty),		was_empty_name = full_name;		end
					end
					if (what.mask)
						opt_D = {'-Dc' '-Dl' '-Di' '-Dh' '-Df'};
						opt_D = opt_D{what.coastRes};
						mask = grdlandmask_m(opt_R, opt_D, opt_e, '-I0.01', '-V');
						Z(mask) = NaN;
					end
				end
				att.GMT_hdr = head;
				known_coords = true;				% Signal that coordinates are known and should not be guessed again
				att.Band(1).NoDataValue = [];		% Don't waist time later trying to NaNify again
				x_min = head(1) - head(8)/2;		x_max = head(2) + head(8)/2;		% Goto pixel registration
				y_min = head(3) - head(9)/2;		y_max = head(4) + head(9)/2;		% But not for att.GMT_hdr(7)
			end
			att.RasterXSize = size(Z,2);		att.RasterYSize = size(Z,1);
			att.Band.XSize = size(Z,2);			att.Band.YSize = size(Z,1);
			att.Corners.LL = [x_min y_min];		att.Corners.UL = [x_min y_max];		% CONFIRMAR
			att.Corners.UR = [x_max y_max];		att.Corners.LR = [x_max y_min];
		end
	end

	if (~isempty(fname)),		att.fname = fname;		end		% Needed in some HDF cases
	if (~isempty(uncomp_name)),	delete(uncomp_name);	end		% Delete uncompressed file

% ----------------------------------------------------------------------------------------
function [Z, head, was_empty] = smart_grid(x, y, z, opt_I, opt_R, opt_C)
% Interpolate the x,y,z points into the smallest grid that is compatible with
% opt_I & opt_R. We do this to accelerate the process of interpolating the L2
% data, which often covers only a small subregion of opt_R.
%
% All x,y,z are expected to be double precision column vectors.
%
% Return Z that is an array of size determined by opt_R and opt_I
% We do this interpolating in a smaller subregion determined by x,y and next
% insert the sub-grid in the array Z.

	was_empty = false;

	x_min = min(x);		x_max = max(x);
	y_min = min(y);		y_max = max(y);
	
	ind = strfind(opt_R, '/');
	inc = str2double(opt_I(3:end));
	Rx_min = str2double(opt_R(3:ind(1)-1));			% Full Region min/max
	Rx_max = str2double(opt_R(ind(1)+1:ind(2)-1));
	Ry_min = str2double(opt_R(ind(2)+1:ind(3)-1));
	Ry_max = str2double(opt_R(ind(3)+1:end));
	
	n_cols = round((Rx_max - Rx_min) / inc + 1);
	n_rows = round((Ry_max - Ry_min) / inc + 1);
	X = linspace(Rx_min, Rx_max, n_cols);
	Y = linspace(Ry_min, Ry_max, n_rows);
	
	% Find rows & columns of Z encompassing all points in x,y
	c1 = 1;		c2 = n_cols;		r1 = 1;		r2 = n_rows;
	t = find((X - x_min) > 0);
	if (~isempty(t) && (t(1) > 1)),	c1 = t(1) - 1;	end
	t = find((X - x_max) > 0);
	if (~isempty(t)),	c2 = t(1);	end
	t = find((Y - y_min) > 0);
	if (~isempty(t) && (t(1) > 1)),	r1 = t(1) - 1;	end
	t = find((Y - y_max) > 0);
	if (~isempty(t)),	r2 = t(1);	end

	if (c2 == 1 || r2 == 1)			% It means all data is outside the opt_R (working area) region
		Z = alloc_mex(n_rows, n_cols, 'single', NaN);
		head = [X(1) X(end) Y(1) Y(end) 0 0 0 inc inc];
		was_empty = true;
		return
	end
	
	if (c1 <= 3 && c2 >= n_cols - 3 && r1 <= 3 && r2 >= n_cols - 3)		% If almost the full Region, just do it all
		[Z, head] = gmtmbgrid_m(x, y, z, opt_I, opt_R, '-Mz', opt_C);
		Z = single(Z);
		return
	end
	
	opt_subR = sprintf('-R%.10f/%.10f/%.10f/%.10f', X(c1), X(c2), Y(r1), Y(r2));
	
	[zz, head] = gmtmbgrid_m(x, y, z, opt_I, opt_subR, '-Mz', opt_C);
	head(1:4) = [Rx_min Rx_max Ry_min Ry_max];
	head(5:6) = [min(zz(:)) max(zz(:))];
	Z = alloc_mex(n_rows, n_cols, 'single', NaN);
	if (~isequal(head(5:6), [0 0]))
		Z(r1:r2, c1:c2) = zz;
	else
		was_empty = true;				% Means no data fall inside the opt_subR sub-region
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, did_scale, att, have_new_nans] = handle_scaling(data, att)
%	If there is a scale factor and/or  add_offset attribute, convert the data
%	to single precision and apply the scaling. 

	have_scale_factor = false;			have_add_offset = false;	did_scale = false;	have_new_nans = false;
	if ( att.Band(1).ScaleOffset(1) ~= 1 )
		have_scale_factor = true;
		scale_factor = att.Band(1).ScaleOffset(1);
	end
	if ( att.Band(1).ScaleOffset(2) ~= 0 )
		have_add_offset = true;
		add_offset = att.Band(1).ScaleOffset(2);
	end

	% Return early if we don't have either one.
	if ~(have_scale_factor || have_add_offset),		return,		end

	if (~isa(data,'single')),		data = single(data);	end
	if (~isnan(att.Band(1).NoDataValue))
		data(data == att.Band(1).NoDataValue) = NaN;
		att.Band(1).NoDataValue = NaN;
		have_new_nans = true;
	end

	did_scale = true;
	if (have_scale_factor && have_add_offset)
		data = cvlib_mex('CvtScale',data, scale_factor, add_offset);
	elseif (have_scale_factor)
		data = cvlib_mex('CvtScale',data, scale_factor);
	else
		data = cvlib_mex('addS',data, add_offset);
	end

	if (have_new_nans)
		att.GMT_hdr(5:6) = grdutils(data,'-L');
	end

% -----------------------------------------------------------------------------------------
function Z = clipMySpikes(Z)
% MODIS L2 SST show an incredible noise level peaking at every other 10 rows of data in
% sensor coordinates. One simple but quite efficient way of reducing (but not eliminating)
% is to replace the screwed value by the average of its neighbors.

	n_rows = size(Z,1);
	indSpikesC = 11:10:n_rows-5;		% Index o spiky values
	indSpikesB = 10:10:n_rows-5;		% Index of the Before spikies
	indSpikesA = 12:10:n_rows-5;		% Index of the After spikies
	for (k = 1:size(Z,2))
		col = Z(:,k);
		col(indSpikesC) = (col(indSpikesA) + col(indSpikesB)) / 2;	% Replace spikies by average of neighbors
		Z(:,k) = col;					% and put it back
	end

% -----------------------------------------------------------------------------------------
function [att, uncomp_name] = get_baseNameAttribs(full_name)
% Get the file's metadata and also tests if an SDS was required but does not exist
	[full_name, uncomp_name] = deal_with_compressed(full_name);
	att = gdalread(full_name, '-M');
	if (att.RasterCount > 0)
		handles = guidata(gcf);
		if (~isempty(handles.SDSinfo) && handles.SDSthis > 1 && ~handles.testedDS)
			handles.testedDS = true;
			guidata(handles.figure1, handles)
			errordlg('Input File has no Sub-Datasets so the silly sub-dataset request forced me to abort.','Error')
			error('empilhador: File has no Sub-Datasets so the silly sub-dataset request forced me to abort')
		end
	end

% -----------------------------------------------------------------------------------------
function [opt_R, opt_I, opt_C, bitflags, flagsID, despike] = sniff_in_OPTcontrol(old_R, att)
% Check the L2config file for particular requests in terms of -R, -I or quality flags
% OPT_R is what the L2config has in
%
% BITFLAGS is a vector with the bit number corresponding to the flgs keys in L2config.txt.
%			Returns [] when no bitflags keywords are found.
% FLAGSID is the subdatset number adress of the flags array (l2_flags)
%			Return 0 when no l2_flags array is found.
% DESPIKE	MODIS L2 SST show an incredible noise level peaking at every other 10 rows of data in
%			sensor coordinates. If the keyword MIR_EMPILHADOR_C exists, DESPIKE is set to true.

	got_flags = false;		bitflags = [];		flagsID = 0;	despike = false;
	opt_I = [];				opt_C = [];
	opt_R = old_R;			% In case we return without finding a new -R

	% OK, this is a transitional code. Before we used to seek the info in OPTcontrol.txt but now
	% the user can choose to use my default values that were written in .../tmp/L2config.txt by
	% push_compute_CB, or alternatively by redirecting the reading to the .../data/L2config.txt
	handles = fish_handles;
	if (~isempty(handles) &&  isfield(handles,'check_L2') && get(handles.check_L2, 'Val'))
		% This should now be the main branch
		if (get(handles.check_L2conf, 'Val'))
			opt_file = [handles.path_data 'L2config.txt'];		% Use User edited config file
			if (~(exist(opt_file, 'file') == 2))
				errordlg(['GHrrrr. You told me read ' opt_file ' but that file does not exist'],'Error')
				error('EMPILHADOR:Sniff_in_Config', 'Oh, were is your head?')
			end
		else
			opt_file = [handles.path_tmp 'L2config.txt'];		% Use our automatically generated file
		end
	else
		mir_dirs = getappdata(0,'MIRONE_DIRS');
		if (~isempty(mir_dirs))
			opt_file = [mir_dirs.home_dir '/data/L2config.txt'];
		else
			return				% Since we cannot find the OPTcontrol file
		end
	end
	if (~(exist(opt_file, 'file') == 2)),		return,		end		% Nickles
	
	fid = fopen(opt_file, 'r');
	c = (fread(fid,'*char'))';      fclose(fid);
	lines = strread(c,'%s','delimiter','\n');   clear c fid;
	m = numel(lines);
	for (k = 1:m)
		if (~strncmp(lines{k},'MIR_EMPILHADOR',9)),	continue,	end
		opt = ddewhite(lines{k}(15:end));
		got_one = false;
		[t,r] = strtok(opt);
		while (t)
			if (strncmp(t,'-R',2)),		opt_R = t;		got_one = true;
			elseif (strncmp(t,'-I',2)),	opt_I = t;		got_one = true;
			elseif (strncmp(t,'-C',2)),	opt_C = t;		got_one = true;
			end
			[t,r] = strtok(ddewhite(r));
		end
		if (got_one),	continue,	end			% Done with this line
		if (strcmp(lines{k}(15:16),'_F'))		% We have a bitflags request
			got_flags = true;					% If it comes here means that we have a flags request
			flaglist = opt;
		elseif (strcmp(lines{k}(15:16),'_C'))	% Despike MODIS SST
			% The key MIR_EMPILHADOR_C has currently one argument only, 'AVG', but this may change
			despike = true;
		end
	end

	if (nargout <= 3),		return,		end		% Used when only MIR_EMPILHADOR was scanned (and no att sent in)

	if (got_flags)
		% Before anything else find the 'fl_flags' array ID. If not found go away right away
		flagsID = find_in_subdatasets(att.AllSubdatasets, 'l2_flags');
		if (~flagsID)
			warndlg('You requested for a FLAGS masking (via L2config.txt) but this file does not have one "l2_flags" array','Warning')
			return
		end

		fmap = {'ATMFAIL' 1;
				'LAND' 2;
				'HIGLINT' 4;
				'HILT' 5;
				'HISATZEN' 6;
				'STRAYLIGHT' 9;
				'CLDICE' 10;
				'COCCOLITH' 11;
				'HISOLZEN' 13;
				'LOWLW' 15;
				'CHLFAIL' 16;
				'NAVWARN' 17;
				'MAXAERITER' 20;
				'CHLWARN' 22;
				'ATMWARN' 23;
				'NAVFAIL' 26;
				'FILTER' 27;
				'SSTWARN' 28;
				'SSTFAIL' 29};
		
		loc = [0 strfind(flaglist, ',') numel(flaglist)+1];		% Find the ',' separator. Add 2 to easy algo
		c = false(1, size(fmap,1));								% Vector with as many elements as input flags
		fmap_names = fmap(:,1);
		for (k = 1:numel(loc)-1)		% Loop over number of input flag keys
			ind = strcmp(flaglist(loc(k)+1 : loc(k+1)-1), fmap_names);
			n = find(ind);
			if (~isempty(n)),	c(n) = true;	end
		end
		bitflags = [fmap{c,2}];
	end

% -----------------------------------------------------------------------------------------
function handles = fish_handles()
% Secure function to get the handles structure. GCF is just too risky
	handles = [];
	[hObj, hFig] = gcbo;
	if (~isempty(hFig))					% It will be empty when calling Empilhador fucntions directly.
		handles = guidata(hFig);		% That is, when the Empilhador figure was never created
	end

% -----------------------------------------------------------------------------------------
function ID = find_in_subdatasets(AllSubdatasets, name, ncmp)
% Find the position in the subdatasets array containing the array called 'name'
% NCMP	is an optional argument containing the N if the strNcmp function.
%		If not provided, defaults to numel(name)

	if (nargin == 2),	ncmp = numel(name);		end
	got_it = false;		ID = 0;
	for (k = 2:2:numel(AllSubdatasets))
		ind = strfind(AllSubdatasets{k}, ' ');
		if ( strncmp(AllSubdatasets{k}(ind(1)+1:end), name, ncmp) )
			got_it = true;		break
		end
	end
	if (got_it),	ID = k - 1;		end		% -1 because the subdataset name is one position before

% -----------------------------------------------------------------------------------------
function [full_name, uncomp_name] = deal_with_compressed(full_name)
% Check if FULL_NAME is a compressed file. If it is, uncompress it and return the uncompressed
% name as FULL_NAME.
% UNCOMP_NAME informs if decompressing was done. If it is not empty, it means file was decompressed.
% Use this info to eventualy remove the FULL_NAME (note that in this case this not the original file name)

	str_d = [];		do_warn = true;		cext = [];
	[PATH,fname,EXT] = fileparts(full_name);
	uncomp_name = [PATH filesep fname];			% Only used if file is compressed
	if (strcmpi(EXT,'.bz2'))
		str_d = ['bzip2 -d -q -f -c ' full_name ' > ' uncomp_name];		cext = EXT;
	elseif (strcmpi(EXT,'.zip') || strcmpi(EXT,'.gz'))
		str_d = ['gunzip -q -N -f -c ' full_name ' > ' uncomp_name];		cext = EXT;
	end

	if (~isempty(str_d))     % File is compressed.
		[pato,fname,EXT] = fileparts(fname);	% Need to remove the true extension
	
		if (do_warn),	aguentabar(0.5,'title',['Uncompressing ' fname EXT cext]);	end
		if (isunix),	s = unix(str_d);
		elseif ispc,	s = dos(str_d);
		else			errordlg('Unknown platform.','Error');	error('Unknown platform.')
		end
		if ~(isequal(s,0))						% An error as occured
			errordlg(['Error decompressing file ' full_name],'Error');
			if (do_warn),   aguentabar(1,'title','By'),		end
			error(['Error decompressing file ' full_name])
		end
		if (do_warn),	aguentabar(1,'title','Donne'),		end
		full_name = uncomp_name;				% The uncompressed file name
	else
		uncomp_name = [];						% So that calling function knows file was not compressed
	end

% -----------------------------------------------------------------------------------------
function fid = write_vtk(handles, grd_out, arg3)
% Help function to write in the VTK format
	if (isa(arg3,'char') && strcmp(arg3,'hdr'))
		nx = (handles.head(2) - handles.head(1)) / handles.head(8) + ~handles.head(7);
		ny = (handles.head(4) - handles.head(3)) / handles.head(9) + ~handles.head(7);
		nz = numel(handles.nameList);
% 		fid = fopen(grd_out, 'wa');
		fid = fopen(grd_out, 'wb','b');
		fprintf(fid, '# vtk DataFile Version 2.0\n');
		fprintf(fid, 'converted from A B\n');
% 		fprintf(fid, 'ASCII\n');
		fprintf(fid, 'BINARY\n');
		fprintf(fid, 'DATASET RECTILINEAR_GRID\n');
		fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
		fprintf(fid, 'X_COORDINATES %d float\n', nx);
		X = linspace(handles.head(1), handles.head(2), nx);
% 		fprintf(fid, '%.6f\n', X);
		fwrite(fid, X, 'real*4');
		fprintf(fid, 'Y_COORDINATES %d float\n', ny);
		X = linspace(handles.head(3), handles.head(4), ny);
% 		fprintf(fid, '%.6f\n', X);
		fwrite(fid, X, 'real*4');
		fprintf(fid, 'Z_COORDINATES %d float\n', nz);
		X = str2double(handles.strTimes);
% 		fprintf(fid, '%.6f\n', X);
		fwrite(fid, X, 'real*4');
		fprintf(fid, 'POINT_DATA %d\n', nx * ny * nz);
		fprintf(fid, 'SCALARS dono float 1\n');
		fprintf(fid, 'LOOKUP_TABLE default\n');
	else
		fid = handles;			% First arg actually contains the file id
		Z = double((arg3)');
% 		ind = isnan(Z);
% 		if (any(ind(:)))
% 			Z(ind) = 0;
% 		end
% 		fprintf(fid, '%.6f\n', Z(:));		
		fwrite(fid, Z(:), 'real*4');		
	end

% ---------------------------------------------------------------------
function empilhador_LayoutFcn(h1)

set(h1,'Position',[520 532 540 280],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Toolbar','none',...
'Name','Empilhador',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[340 39 191 117],'Style','frame');
uicontrol('Parent',h1,'Position',[458 40 20 15],'String','S','Style','text');
uicontrol('Parent',h1,'Position',[350 98 20 15],'String','W','Style','text');
uicontrol('Parent',h1,'Position',[502 99 20 15],'String','E','Style','text');
uicontrol('Parent',h1,'Position',[456 122 20 15],'String','N','Style','text');

uicontrol('Parent',h1,'Position',[6 254 501 21],...
'BackgroundColor',[1 1 1],...
'Call','empilhador(''edit_namesList_CB'',gcbo,guidata(gcbo))',...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Name of an ascii file with the grids list. One grid name per row',...
'Tag','edit_namesList');

uicontrol('Parent',h1,'Position',[507 253 23 23],...
'Call','empilhador(''push_namesList_CB'',gcbo,guidata(gcbo))',...
'Tooltip','Browse for a grids list file',...
'Tag','push_namesList');

uicontrol('Parent',h1, 'Position',[344 233 165 15],...
'Call','empilhador(''radio_conv2netcdf_CB'',gcbo,guidata(gcbo))',...
'String','Convert to 3D netCDF file',...
'Style','radiobutton',...
'Tooltip','Take a list of files and create a single 3D netCDF file',...
'Value',1,...
'Tag','radio_conv2netcdf');

uicontrol('Parent',h1, 'Position',[344 214 150 15],...
'Call','empilhador(''radio_conv2vtk_CB'',gcbo,guidata(gcbo))',...
'String','Convert to 3D VTK file',...
'Style','radiobutton',...
'Tooltip','Take a list of files and create a single 3D VTK file',...
'Value',0,...
'Tag','radio_conv2vtk');

uicontrol('Parent',h1, 'Position',[344 195 140 15],...
'Call','empilhador(''radio_multiBand_CB'',gcbo,guidata(gcbo))',...
'String','Make multi-band image',...
'Style','radiobutton',...
'Tooltip','Take a list of image files and create a single multi-band image',...
'Tag','radio_multiBand');

uicontrol('Parent',h1,'Position',[345 167 85 15],...
'Call','empilhador(''check_L2_CB'',gcbo,guidata(gcbo))',...
'String','L2 magic',...
'Style','checkbox',...
'Tooltip','Do all hard work to process and reference MODIS L2 products',...
'Tag','check_L2');

uicontrol('Parent',h1,'Position',[430 167 155 15],...
'Call','empilhador(''check_L2conf_CB'',gcbo,guidata(gcbo))',...
'String','Use config file',...
'Style','checkbox',...
'Vis','off',...
'Tooltip','Fetch fine control from (MironeRoot)/data/L2config.txt file',...
'Tag','check_L2conf');

uicontrol('Parent',h1,'Position',[350 133 110 15],...
'Call','empilhador(''check_region_CB'',gcbo,guidata(gcbo))',...
'String','Use sub-region?',...
'Style','checkbox',...
'Tooltip','Perform computations inside a data sub-region',...
'Tag','check_region');

uicontrol('Parent',h1,'Position',[400 104 71 21],...
'BackgroundColor',[1 1 1],...
'Call','empilhador(''edit_north_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_north');

uicontrol('Parent',h1,'Position',[350 79 71 21],...
'BackgroundColor',[1 1 1],...
'Call','empilhador(''edit_west_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_west');

uicontrol('Parent',h1,'Position',[450 79 71 21],...
'BackgroundColor',[1 1 1],...
'Call','empilhador(''edit_east_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_east');

uicontrol('Parent',h1,'Position',[400 54 71 21],...
'BackgroundColor',[1 1 1],...
'Call','empilhador(''edit_south_CB'',gcbo,guidata(gcbo))',...
'Enable','off',...
'Style','edit',...
'Tag','edit_south');

uicontrol('Parent',h1,'Position',[6 9 325 236],...
'BackgroundColor',[1 1 1],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_list');

uicontrol('Parent',h1,'Position',[440 5 90 21],...
'Call','empilhador(''push_compute_CB'',gcbo,guidata(gcbo))',...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'String','Compute',...
'Tag','push_compute');

% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
function varargout = l2_choices(varargin)
% ... 
	hObject = figure('Vis','off');
	l2_choices_LayoutFcn(hObject);
	handles = guihandles(hObject);

	handles.out.georeference = 0;
	handles.out.nearneighbor = 1;
	handles.out.mask = 0;
	handles.out.coastRes = 0;
	handles.out.quality = '';

	AllSubdatasets = varargin{1};
	ID = find_in_subdatasets(AllSubdatasets, 'qual_sst', 8);
	if (ID)
		set([handles.popup_quality handles.text_quality],'Enable','on')
		ind = strfind(AllSubdatasets{ID}, '=');
		handles.out.qualSDS = AllSubdatasets{ID}(ind+1:end);
	end

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	% UIWAIT makes l2_choices wait for user response
	uiwait(handles.figure1);
	handles = guidata(handles.figure1);
	varargout{1} = handles.out;
	delete(handles.figure1),	drawnow

% -----------------------------------------------------------------------
function radio_sensor_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_interpMin handles.radio_interpNear handles.check_landMask],'Enable','off')
	set(handles.radio_georef,'Val',0)

% -----------------------------------------------------------------------
function radio_georef_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set([handles.radio_interpMin handles.radio_interpNear handles.check_landMask],'Enable','on')
	set(handles.radio_sensor,'Val',0)

% -----------------------------------------------------------------------
function radio_interpNear_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set(handles.radio_interpMin,'Val',0)

% -----------------------------------------------------------------------
function radio_interpMin_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end
	set(handles.radio_interpNear,'Val',0)

% -----------------------------------------------------------------------
function check_landMask_CB(hObject, handles)
	if ( get(hObject,'Val') )
		info = getappdata(0,'gmt_version');		% It should never be empty
		if (info.crude ~= 'y')
			errordlg('You do not have GMT installed. So, no coastlines no masking.','Error')
			set(hObject,'Val',0),	return
		elseif (info.full == 'y'),			handles.out.coastRes = 5;
		elseif (info.high == 'y'),			handles.out.coastRes = 4;
		elseif (info.intermediate == 'y'),	handles.out.coastRes = 3;
		elseif (info.low == 'y'),			handles.out.coastRes = 2;
		else								handles.out.coastRes = 1;
		end
		guidata(handles.figure1, handles)
	end

% -----------------------------------------------------------------------
function push_OK_CB(hObject, handles)	
	handles.out.georeference = get(handles.radio_georef, 'Val');
	handles.out.nearneighbor = get(handles.radio_interpNear, 'Val');
	handles.out.mask = get(handles.check_landMask, 'Val');
	if ( strcmp(get(handles.popup_quality,'Enable'),'on') )		% Only if we are using it
		ind = get(handles.popup_quality, 'Val');
		if (ind ~= 1)
			handles.out.quality = 3 - ind;
		end
	end
	guidata(handles.figure1, handles)
	uiresume(handles.figure1);

% -----------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
% The GUI is still in UIWAIT, do UIRESUME
	handles = guidata(hObject);
	handles.out = '';		% User gave up, return nothing
	guidata(handles.figure1, handles);
	uiresume(handles.figure1);

% -----------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
% Check for "escape"
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles = guidata(hObject);
		handles.out = '';	% User said no by hitting escape
		guidata(handles.figure1, handles);
		uiresume(handles.figure1);
	end

% --- Creates and returns a handle to the GUI figure. 
function l2_choices_LayoutFcn(h1)

set(h1, 'Position',[520 400 291 148],...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','L2 product choices',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[2 124 265 22],...
'FontAngle','oblique',...
'FontName','Helvetica',...
'FontSize',11,...
'FontWeight','demi',...
'ForegroundColor',[1 0.50196 0],...
'String','L2 Level product - Needs decisions',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 101 161 16],...
'Call',@l2_choices_uiCB,...
'FontName','Helvetica',...
'String','Plot in Sensor coordinates',...
'Style','radiobutton',...
'Tooltip','Plot data as it is in file (faster, but deformed)',...
'Value',1,...
'Tag','radio_sensor');

uicontrol('Parent',h1, 'Position',[10 78 231 16],...
'Call',@l2_choices_uiCB,...
'FontName','Helvetica',...
'String','Compute georeferenced grid (takes time)',...
'Style','radiobutton',...
'Tooltip','Reinterpolate data to get a georeferenced grid. (It may take 1 minute)',...
'Tag','radio_georef');

uicontrol('Parent',h1, 'Position',[31 55 131 16],...
'Call',@l2_choices_uiCB,...
'Enable','off',...
'FontName','Helvetica',...
'String','Minimum curvature',...
'Style','radiobutton',...
'Value',1,...
'Tooltip','Interpolation method',...
'Tag','radio_interpMin');

uicontrol('Parent',h1, 'Position',[31 36 110 16],...
'Call',@l2_choices_uiCB,...
'Enable','off',...
'FontName','Helvetica',...
'String','Nearneighbor',...
'Style','radiobutton',...
'Tooltip','Interpolation method',...
'Value',0,...
'Tag','radio_interpNear');

uicontrol('Parent',h1, 'Position',[10 6 51 22],...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'ListboxTop',0,...
'String',{'2'; '1'; '0'},...
'Style','popupmenu',...
'Tooltip','Select the least quality level. 2 - worst - means all values. 0 - only the best',...
'Value',3,...
'Tag','popup_quality');

uicontrol('Parent',h1, 'Position',[64 10 75 15],...
'Enable','off',...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','Quality factor',...
'Style','text',...
'Tag','text_quality');

uicontrol('Parent',h1, 'Position',[180 48 110 15],...
'Call',@l2_choices_uiCB,...
'Enable','off',...
'String','Apply Land Mask',...
'Style','checkbox',...
'Tooltip','Mask reinterpolated Land pixels (need high definition coastlines installed)',...
'Tag','check_landMask');

uicontrol('Parent',h1, 'Position',[215 7 66 23],...
'Call',@l2_choices_uiCB,...
'FontName','Helvetica',...
'FontSize',9,...
'String','OK',...
'Tag','push_OK');

function l2_choices_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

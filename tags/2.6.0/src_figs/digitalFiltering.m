function varargout = digitalFiltering(varargin)
% digitalFiltering applies methods of digital image analysis to images AND grids.
%
% This includes a broad suite of smoothing (low-pass) filters, as well as contrast
% enhancement filters, edge enhancement filters, edge detection filters, general high-pass filters, etc.
%
%   USAGE:
%       DIGITALFILTERING(H_IMG), where H_IMG is the handle to an image object, filters the
%       image pointed by H_IMG by fishing its CData content. The Apply button will also
%       update the original image.
%       DIGITALFILTERING(IMG), where IMG is an MxN or MxNx3 array containing an image,
%       filters that image.
%       IMG = DIGITALFILTERING(IMG), returns the filtered array but this array can also be
%       saved as an image file with the "Save image" button.
%       DIGITALFILTERING(H_IMG,GRD,CMAP), where GRD is a MxN double or single array and H_IMG is a
%       handle to an image object containing the GRD representation, filters the GRD array and
%       updated the original image as well as the "Filtered image" axes of this GUI. The GRD
%       image's representation is computed with a linear scaling on the [0 255] interval.
%       ( img = uint8(round( ((GRD - GRD_z_min) / (GRD_z_max - GRD_z_min))*255 )) )
%       DIGITALFILTERING(IMG,GRD,CMAP), where GRD is a MxN double or single array; IMG is either
%       a uint8 array containing the GRD's image representation or a [] argument and CMAP is an optional
%       colormap, filters the GRD array and shows the result in the "Filtered image" axes.
%       [GRD, IMG] = DIGITALFILTERING(IMG,GRD), returns the filtered GRD array and its
%       image representation.
%
%       The "Apply" button performs the filtering and updates the "Filtered image" axes as well as
%       the parent image (if it applyies), but does not return anything.
%
%       The "Apply n return" button performs the filtering, updates the parent image (if it applyies),
%       and returns the output if it was requested.
%
%       The "Cancel" button closes this figure and undo any previous change it might have done into
%       a parent figure.
%
%       The left panel contains the available filters. Double click to expand the fields
%
%       Filters imported with the "Load external filter" button will be placed in the
%       "General User-defined (mxn)" fied. Note: you may need to click again on this field
%       to update for the imported filter. The coefficients of this type may be entered
%       directly on the table. Click twice on the table element to enter value. First click
%       change to edit mode and the second one put the focus on table element.
%       The default filters are defined inside the loadDefFilters.m function. You can also add
%       more filters here. See the m file for howto basic instructions.
%
%   EXAMPLES:
%       Operate on parent window that have an image displayed on it.
%           img = imread('peppers.png');
%           h_img = imshow(img);
%           digitalFiltering(h_img)
%
%       Operate on a double array
%           [X,Y,Z] = peaks(128);
%           grd = digitalFiltering([],Z);
%       Do some filtering and see the result
%           surf(X,Y,grd)
%
%       Operate on a double array, see changes on parent figure and return the filtered array
%           [X,Y,Z] = peaks(128);
%           h_img = imagesc(Z);
%           grd = digitalFiltering(h_img,Z,jet(256));
%
%       Demonstration mode
%           digitalFiltering

%   THIS GUI NEEDS THE IMAGE PROCESSING TOOLBOX TO RUN

%   CREDITS
%       This GUI was designed trying to mimic the "Digital Filtering" module of the Golden Software
%       SURFER program. The filter coefficients and their descriptions come from this module.
%
%       The left panel tree expander uses much code and ideas of the StructBrowser_gui from
%       H.Lahdili,(hassan.lahdili@crc.ca)
%
%       The table hosting the coefficients is an adaptation of the mltable from
%       Morris Maynard
%
%   AUTHOR
%       Joaquim Luis    24-Jan-2006
%
%       M-File changed by desGUIDE and manualy edited after

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

	hObject = figure('Vis','off');
	digitalFiltering_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	% Choose default command line output for digitalFiltering
	handles.output.cancel = 0;		% It will be set to one if the cancel button is hit or the fig is killed
	handles.img_orig = [];
	handles.grd_orig = [];
	handles.h_calling_fig = [];		% handle of the calling figure (if it will be the case)
	handles.h_calling_img = [];		% handle of the calling image (if it will be the case)
	handles.tree_indice = 1;		% Index to the select tree item in its full expanded form.
	y_dir = 0;						% Axes YDir property
	cmap = [];

	% Find out what is going out in the first and second arguments (if output was requested)
	n_argout = nargout;     n_argin = nargin;
	if (n_argout == 1 || n_argout == 2 && n_argin >= 2)
		first_out = 'grd';					% Signal that the first output arg (if requested) is the grid
	elseif (n_argout == 1 && n_argin == 1)
		first_out = 'img';					% Signal that the output arg (if requested) is the image
	elseif (n_argout > 0)
		errordlg('Requested output is not possible with the current choice of inputs.','ERROR')
		delete(hObject)
		return
	end

	if (n_argin >= 1)
		if (ishandle(varargin{1}))              % First argument is an image handle, or ... error
			if ( ~strcmp(get(varargin{1},'type'),'image') )
				errordlg('ERROR: only image graphical objects are supported.','ERROR')
				delete(hObject)
				return
			end
			handles.img_orig = get(varargin{1},'CData');
			handles.h_calling_img = varargin{1};
			handles.h_calling_fig = get(get(handles.h_calling_img,'Parent'),'Parent');
			y_dir = get(get(varargin{1},'Parent'),'YDir');
			if (ndims(handles.img_orig) == 2)   % If image is indexed, get the colormap
				cmap = get(handles.h_calling_fig,'ColorMap');
			end
		else
			handles.img_orig = varargin{1};     % First argument is an image array
		end
		if (n_argin >= 2)
			handles.grd_orig = varargin{2};
			if (isempty(handles.img_orig))      % We don't have a grid image representation. Create it
				handles.img_orig = scaleto8(handles.grd_orig);
			end
		end
		if (n_argin == 3)
			cmap = varargin{3};
		end

		set(handles.edit_imageIn,'Enable','off')
		set(handles.push_imageIn,'Enable','off')
		set(handles.text_loadImage,'Enable','off')
	else				% Demo mode
		set(handles.push_apply_and_return,'Visible','off')
		handles.img_orig = imread('peppers.png');
	end

	if (~isempty(handles.h_calling_fig))
		handMir = guidata(handles.h_calling_fig);
		handles.home_dir = handMir.home_dir;
		handles.work_dir = handMir.work_dir;
		handles.last_dir = handMir.last_dir;
	else
		handles.home_dir = cd;
		handles.work_dir = cd;		handles.last_dir = cd;	% To not compromize put_or_get_file
	end

	% Put the image in the right axes
	[m,n,k] = size(handles.img_orig);
	bg_aspect = m / n;
	handles.h_img1 = image(handles.img_orig,'Parent',handles.axes2);
	set(handles.axes2,'PlotBoxAspectRatio',[1 bg_aspect 1],'Visible','off')
	handles.h_img2 = image(handles.img_orig,'Parent',handles.axes3);
	set(handles.axes3,'PlotBoxAspectRatio',[1 bg_aspect 1],'Visible','off')

	if (ischar(y_dir))
		set(handles.axes2,'YDir',y_dir)
		set(handles.axes3,'YDir',y_dir)
	end

	if (~isempty(cmap))				% If we have a colormap, apply it
		set(hObject,'ColorMap',cmap)
	end

	% Load the predefined filters
	filter_struct = load_defFilters;
	sz_argin = size(filter_struct,2);
	struct_names = cell(1, sz_argin);
	struct_values = cell(1, sz_argin);
	for I = 1:sz_argin
		struct_names{I} = filter_struct{1, I};  % Structure descriptive name
		struct_values{I} = filter_struct{2, I}; % 
		if isstruct(struct_values{I})
			struct_names{I} = ['+ ' struct_names{I}];
		end
	end
	handles.all_names = filter_struct(3);       % {mx2} struct field name & name to show up in the tree
	handles.filt_desc = filter_struct(4);       % {mx1} text description of struct field

	% Populate the Structure listbox
	set(handles.listbox1,'string',struct_names)

	handles.h_txt_info = findobj(hObject,'Tag','text_info');
	handles.txt_info_pos = get(handles.h_txt_info,'Position');

	% save structures names and values to handles
	handles.struct_names = struct_names';
	handles.struct_values = struct_values';

	set(handles.listbox_nRows,'String',3:2:35,'Value',1,'Enable','off')
	set(handles.listbox_nCols,'String',3:2:35,'Value',1,'Enable','off')

	% Initialize table calling arguments
	cell_data = {0 0 0; 0 1 0; 0 0 0};
	handles.custom_filt = [0 0 0; 0 1 0; 0 0 0];    % Default (null filter) custom filter
	columninfo.format = '%.3g';
	handles.format = columninfo.format;
	rowHeight = 16;
	gFont.size=9;
	gFont.name='Helvetica';

	mltable_j(hObject, handles.tblParams, 'CreateTable', columninfo, rowHeight, cell_data, gFont);
	% so clicking outside the table will finish edit in progress...
	endfcn = sprintf('mltable_j(%14.13f, %14.13f, ''SetCellValue'');', hObject, handles.tblParams);
	set(hObject,'buttondownfcn',endfcn);

	guidata(hObject, handles);				% Update handles structure
	set(hObject,'Visible','on');

	% Check if we have to do a uiwait and what goes out (if anything)
	%if (~isempty(handles.h_calling_img) || n_argout > 0)	% If called to operate in a parent figure
	if (n_argout > 0)						% If called with output, we must wait
		uiwait(handles.figure1);							% or with output, we must wait
		handles = guidata(hObject);
		if (handles.output.cancel)			% Don't try to output eventual non-existing variables
			if (n_argout == 1),		varargout{1} = [];		end
			if (n_argout == 2),		varargout{1} = [];		varargout{2} = [];  end
			delete(handles.figure1);		% The figure can be deleted now
			return
		end
		if (n_argout == 2)					% Easy case
			varargout{1} = handles.output.grd;
			varargout{2} = handles.output.img;
			delete(handles.figure1);		% The figure can be deleted now
		elseif (n_argout == 1)				% Not so easy case. We must find what is going out: grid or image?
			if (strcmp(first_out,'img'))
				varargout{1} = handles.output.img;
			else
				varargout{1} = handles.output.grd;
			end
			delete(handles.figure1);		% The figure can be deleted now
		end
	end

% ------------------------------------------------------------------------
function listbox1_CB(hObject, handles)
% Poor man simulation of a tree viewer

	index_struct = get(hObject,'Value');
	struct_names = handles.struct_names;
	struct_values = handles.struct_values;
	todos = get(hObject,'String');

	indent = '       ';
	root_1 = struct_names{index_struct};
	is_indent = strfind(root_1, indent);
	if (isempty(is_indent)),    level = 0;
	else                        level = (is_indent(end) - 1)/7 + 1;
	end

	struct_val = struct_values{index_struct};
	all_names = handles.all_names{1};

	if isstruct(struct_val)
		fields =  fieldnames(struct_val);
		names = cell(1,length(fields));
		for i = 1:length(fields)
			idx = strcmp(fields{i},all_names(:,1));
			names{i} = all_names{idx,2};
			if isstruct(struct_val(1).(fields{i}))
				fields{i} = ['+ ' fields{i}];
				names{i} = ['+ ' names{i}];            
			end
		end
	end

	% Display filter info
	name_clean = ddewhite(struct_names{index_struct});
	if (name_clean(1) == '+' || name_clean(1) == '-'),		name_clean = name_clean(3:end);		end
	idx = find(strcmp(name_clean,all_names(:,1)));
	handles.tree_indice = idx;
	pos = handles.txt_info_pos;
	[outstring,newpos] = textwrap(handles.h_txt_info,handles.filt_desc{1}(idx));
	set(handles.h_txt_info,'String',outstring,'Position',[pos(1),pos(2),pos(3),newpos(4)])

	if (isnumeric(struct_val))
		ind = strcmp('              General User-defined (mxn)',todos(index_struct));
		if (any(ind))				% Whatever the custom filter is (it might have been loaded from a file), use it
			struct_val = handles.custom_filt;
		end
		update_table(handles,struct_val)        % Update the table
		set_enable_OnOff(handles,idx)
	end

	% If double-click, and is struct, expand structure, and show fields
	if (strcmp(get(handles.figure1, 'SelectionType'), 'open')) % if double click
		idxP = strfind(struct_names{index_struct}, '+');
		idxM = strfind(struct_names{index_struct}, '-');
		if ~isempty(idxP)
			[struct_names, struct_values] = expand_struct(struct_names, struct_values, ...
				index_struct, fields, level, idxP);
		elseif  ~isempty(idxM)
			[struct_names, struct_values] = shrink_struct(struct_names, struct_values, ...
				index_struct, fields, level, idxM);
		end
		names = cell(length(struct_names),1);
		for (i = 1:length(struct_names))
			name_clean = ddewhite(struct_names{i});
			if (name_clean(1) == '+' || name_clean(1) == '-'),       name_clean = name_clean(3:end);     end
			id1 = strcmp(name_clean,all_names(:,1));				% Find index to pretended name
			id2 = strfind(struct_names{i},name_clean);				% Find index of starting text (after the blanks)
			names{i} = [struct_names{i}(1:id2-1) all_names{id1,2}];        
		end
		set(handles.listbox1,'String',names);
		handles.struct_names = struct_names;
		handles.struct_values = struct_values;
	end

	guidata(hObject, handles);

% ------------------------------------------------------------------------
function cell_array = indent_cell(cell_array, level)

	indent = '       ';				indent_app = [];
	for (k = 1:level+1),			indent_app = [indent_app indent];   end
	for (i=1:length(cell_array)),	cell_array{i} = [indent_app cell_array{i}];     end

% ------------------------------------------------------------------------
function [struct_names, struct_values] = expand_struct(struct_names, struct_values, idx, fields, level, idxP)
% expand structure if '+' is double-clicked and update the structure tree

	size_val = size(struct_values{idx});
	if (size_val(1) ~= 1),  struct_values{idx} = (struct_values{idx})';     end
	N = size_val(2);
	names_be = struct_names(1:idx);
	names_af = struct_names(idx + 1:length(struct_names));
	values_be = struct_values(1:idx);
	values_af = struct_values(idx + 1:length(struct_names));
	if N == 1           % if the structure is of size 1 x 1
		names_app = indent_cell(fields, level);
		values_app = cell(1,length(fields));
		for i = 1:length(fields)
			if (fields{i}(1) == '+' || fields{i}(1) == '-')
				fields{i} = fields{i}(3:end);
			end
			values_app{i} = struct_values{idx}.(fields{i});
		end
		struct_names = [names_be; names_app; names_af];
		struct_values = [values_be; values_app'; values_af];
		struct_names{idx}(idxP) = '-';
	else                % if the structure is of size 1 x N
		names_app = cell(N,1);
		values_app = cell(N,1);
		struct_name = struct_names{idx};
		struct_name = remove_indent(struct_name);
		for (j = 1:N)
			names_app(j) = indent_cell(cellstr(strcat(struct_name,'(', num2str(j),')')), level);
		end
		for (j = 1:N)
			values_app{j} = struct_values{idx}(j) ;
		end
		struct_names = [names_be; names_app; names_af];
		struct_values = [values_be; values_app; values_af];
		struct_names{idx}(idxP) = '-';
	end

% ------------------------------------------------------------------------
function [struct_names, struct_values] = shrink_struct(struct_names, struct_values, idx, fields, level, idxM)
% shrink structure if '- ' is double-clicked
	struct_names{idx}(idxM) = '+';
	indent = '       ';
	if ((idxM-1)/7 - level) == 0
		num_steps = 0;
		is_indent_select = strfind(struct_names{idx}, indent);
		if ~isempty(is_indent_select)
			for (k = idx+1 : length(struct_names))
				is_indent = strfind(struct_names{k}, indent);
				if (isempty(is_indent) || is_indent(end) - is_indent_select(end) <= 0),  break;  end
				num_steps = num_steps + 1;
			end
		else
			for (k = idx+1 : length(struct_names))
				is_indent = strfind(struct_names{k}, indent);
				if (isempty(is_indent)),    break;      end
				num_steps = num_steps + 1;
			end
		end
	else
		num_steps = length(fields);
	end
	names_be = struct_names(1:idx);
	names_app = struct_names(idx+num_steps+1:length(struct_names));
	values_be = struct_values(1:idx);
	values_app = struct_values(idx+num_steps+1:length(struct_names));
	struct_names = [names_be; names_app];
	struct_values = [values_be; values_app];

% ------------------------------------------------------------------------
function str1 = remove_indent(str0)
% remove indent keeping '+ ' and '- '
	c = find(isspace(str0));
	if (~isempty(c)),   str1 = str0(max(c)+1:end);
	else                str1 = str0;
	end

% --------------------------------------------------------------------
function listbox_nRows_CB(hObject, handles)
% Increase/decrease number of rows of matrix filters that allow edition

	info = get(handles.tblParams, 'userdata');
	[mt,nt] = size(info.txtCells);      % Size of current table
	nr = get(hObject,'Value') * 2 + 1;  % Remember that rows start at 3 and are odd
	nc = get(handles.listbox_nCols,'Value') * 2 + 1;

	if (mt > nr)                     % Must decrease number of rows
		mltable_j(handles.figure1, handles.tblParams, 'DelRow',mt-nr)     % Remove mt-nr rows
	elseif (mt < nr)                 % Must increase number of rows
		mltable_j(handles.figure1, handles.tblParams, 'AddRow',nr-mt)     % Add nr-mt rows
	end

	info = get(handles.tblParams, 'userdata');      % Get the updated version
	cur_field_name = handles.all_names{1}{handles.tree_indice};
	switch cur_field_name
		case 'GaussLowPass',    f = filts_coef(cur_field_name,[nr nc],1);   % 1 = std
		case 'InvDist',         f = filts_coef(cur_field_name,[nr nc],1);   % 1 = power of inverse distance
		case 'GeneralUD'
			f = zeros(nr,nc);
			f(median(1:2:nr),median(1:2:nc)) = 1;
	end
	info.data = num2cell(f);
	for (i=1:numel(f))
		set(info.txtCells(i),'String',num2str(f(i),handles.format))
	end
	set(handles.tblParams, 'userdata',info);

% --------------------------------------------------------------------
function listbox_nCols_CB(hObject, handles)
% Increase/decrease number of columns of matrix filters that allow edition

	info = get(handles.tblParams, 'userdata');
	[mt,nt] = size(info.txtCells);      % Size of current table
	nc = get(hObject,'Value') * 2 + 1;  % Remember that cols start at 3 and are odd
	nr = get(handles.listbox_nRows,'Value') * 2 + 1;

	if (nt > nc)                     % Must decrease number of cols
		mltable_j(handles.figure1, handles.tblParams, 'DelCol',nt-nc)     % Remove nt-nc columns
	elseif (nt < nc)                 % Must increase number of cols
		mltable_j(handles.figure1, handles.tblParams, 'AddCol',nc-nt)     % Add nc-nt columns
	end

	info = get(handles.tblParams, 'userdata');      % Get the updated version
	cur_field_name = handles.all_names{1}{handles.tree_indice};
	switch cur_field_name
		case 'GaussLowPass',    f = filts_coef(cur_field_name,[nr nc],1);
		case 'InvDist',         f = filts_coef(cur_field_name,[nr nc],2);   % A pow (=2) nao ta ainda param
		case 'GeneralUD'
			f = zeros(nr,nc);
			f(median(1:2:nr),median(1:2:nc)) = 1;
	end
	info.data = num2cell(f);
	for (i=1:numel(info.txtCells))
		set(info.txtCells(i),'String',num2str(f(i),handles.format))
	end
	set(handles.tblParams, 'userdata',info);

% --------------------------------------------------------------------
function update_table(handles,struct_val)
% Update the table in terms of number of rows/columns and its contents
	info = get(handles.tblParams, 'userdata');
	[mt,nt] = size(info.txtCells);      % Size of current table
	[mf,nf] = size(struct_val);         % Size of current filter
	nc = get(handles.listbox_nCols,'Value') * 2 + 1;        % Remember that row & cols start at 3 and
	nr = get(handles.listbox_nRows,'Value') * 2 + 1;        % increase by two. That is they are odd from 3 on
	if (mf ~= nr)       % Need to update the Rows listbox value
		set(handles.listbox_nRows,'Value',fix((mf-1)/2));
	end
	if (nf ~= nc)       % Need to update the Cols listbox value
		set(handles.listbox_nCols,'Value',fix((nf-1)/2))
	end
	if (mt > mf)                     % Must decrease number of rows
		mltable_j(handles.figure1, handles.tblParams, 'DelRow',mt-mf)     % Remove mt-mf rows
	elseif (mt < mf)                 % Must increase number of rows
		mltable_j(handles.figure1, handles.tblParams, 'AddRow',mf-mt)     % Add mf-mt rows
	end

	if (nt > nf)                     % Must decrease number of cols
		mltable_j(handles.figure1, handles.tblParams, 'DelCol',nt-nf)     % Remove nt-nf columns
	elseif (nt < nf)                 % Must increase number of cols
		mltable_j(handles.figure1, handles.tblParams, 'AddCol',nf-nt)     % Add nf-nt columns
	end

	info = get(handles.tblParams, 'userdata');      % Get the updated version
	info.data = num2cell(struct_val);
	for ( i = 1:numel(info.txtCells) )
		set(info.txtCells(i),'String',num2str(struct_val(i),handles.format))
	end
	set(handles.tblParams, 'userdata',info);

% --------------------------------------------------------------------
function set_enable_OnOff(handles,index)
% Set the listboxs and table on/off dependinf on the filter specifities
	cur_field_name = handles.all_names{1}{index};
	info = get(handles.tblParams, 'userdata');
	info.no_edit = 1;
	switch cur_field_name
		case {'MovingAverage' 'DistWeighting' 'InvDist' 'GaussLowPass' 'GeneralUD'}
			set(handles.listbox_nRows,'Enable','on')
			set(handles.listbox_nCols,'Enable','on')
			if ( strcmp(cur_field_name,'GeneralUD') )
				set(handles.tblParams,'HitTest','on')       % This is the only one editable
				info.no_edit = 0;
			else
				set(handles.tblParams,'HitTest','off')
			end
		otherwise
			set(handles.listbox_nRows,'Enable','off')
			set(handles.listbox_nCols,'Enable','off')
	end
	set(handles.tblParams, 'userdata',info);

% --------------------------------------------------------------------
function edit_imageIn_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),   return;    end
	push_imageIn_CB(gcbo,guidata(gcbo),fname)

% --------------------------------------------------------------------
function push_imageIn_CB(hObject, handles, fname)
	if (nargin == 3),   fname = [];   end
	if (isempty(fname))
		[FileName,PathName] = put_or_get_file(handles,{ ...
			'*.jpg', 'JPEG image (*.jpg)'; ...
			'*.png', 'Portable Network Graphics(*.png)'; ...
			'*.bmp', 'Windows Bitmap (*.bmp)'; ...
			'*.hdf', 'Hieralchical Data Format (*.hdf)'; ...
			'*.gif', 'GIF image (*.gif)'; ...
			'*.pcx', 'Windows Paintbrush (*.pcx)'; ...
			'*.ras', 'SUN rasterfile (*.ras)'; ...
			'*.tif', 'Tagged Image File (*.tif)'; ...
			'*.xwd', 'X Windows Dump (*.xwd)'; ...
			'*.*', 'All Files (*.*)'}, ...
			'Select image format','get');
		if isequal(FileName,0);     return;     end
		fname = [PathName FileName];
		set(handles.edit_imageIn,'String',fname)
	end
	try             % Use a try because if the name was given via edit box it may be wrong
		[handles.img_orig,cmap] = imread(fname);
	catch
		errordlg(['Error: -> ' fname ' does not exist or is not a valid image file'],'Error')
		return
	end

	[m,n,k] = size(handles.img_orig);
	if (k == 1 && isempty(cmap))
		set(handles.figure1,'Colormap', gray(256));
	elseif (~isempty(cmap))
		set(handles.figure1,'Colormap', cmap);
	end

	% Compute image aspect ratio and set axes 'PlotBoxAspectRatio' to it
	aspect = m / n;
	handles.h_img1 = image(handles.img_orig,'Parent',handles.axes2);
	set(handles.axes2,'PlotBoxAspectRatio',[1 aspect 1],'Visible','off')
	handles.h_img2 = image(handles.img_orig,'Parent',handles.axes3);
	set(handles.axes3,'PlotBoxAspectRatio',[1 aspect 1],'Visible','off')

	guidata(handles.figure1,handles)

% ------------------------------------------------------------------------
function [grd,img] = push_apply_CB(hObject, handles)
	grd = [];
	info = get(handles.tblParams, 'userdata');
	f = cell2mat(info.data);				% Fish out the filter coefs

	if (isempty(handles.grd_orig)),		img = handles.img_orig;			% We are wowrking on a image
	else								img = handles.grd_orig;			% We are working on a grid (2D array)
	end

	set(handles.figure1,'pointer','watch')
	flags = 6;			% This corresponds to the 'conv' option to imfilter_mex
	if (~isempty(handles.grd_orig))			% We are working on a grid.
		grd = LocalImfilter(handles.grd_orig,f,'symmetric',flags);   % Filter it
		img = scaleto8(grd);
	else									% We are working only with an image
		if (ndims(handles.img_orig) == 2)
			img = LocalImfilter(handles.img_orig,f,'symmetric',flags);
		else								% Due to a imfilter effeciency bug, this is much faster.
			img(:,:,1) = LocalImfilter(handles.img_orig(:,:,1),f,'symmetric',flags);
			img(:,:,2) = LocalImfilter(handles.img_orig(:,:,2),f,'symmetric',flags);
			img(:,:,3) = LocalImfilter(handles.img_orig(:,:,3),f,'symmetric',flags);
		end
	end

	if (get(handles.radio_enhancedImg,'Val'))	% Compute the enhanced image
		if (f( fix(numel(f)/2) + 1) > 0)
			img = img_fun('imadd', handles.img_orig, img);
		else
			img = img_fun('imsubtract', handles.img_orig, img);
		end
	end
	set(handles.h_img2,'CData',img);			% Put the image in the lower axes

	if (~isempty(handles.h_calling_img) && ishandle(handles.h_calling_img))   
		try
			set(handles.h_calling_img,'CData',img);		% Put the filtered image in the original figure
		catch
			warndlg('Why did you kill the calling figure when you requested to process its image?','Stupid atittude')
		end
	end
	set(handles.figure1,'pointer','arrow')

% ------------------------------------------------------------------------
function b = LocalImfilter(a, h, boundary, flags)
% Minimalist code necessary to prepare things to the mex

	rank_a = ndims(a);		rank_h = ndims(h);

	im_size = size(a);
	%Calculate the number of pad pixels
	filter_center = floor((size(h)+1)/2);
	pad = [size(h)-filter_center ones(1,rank_a-rank_h)];
	im_size = int32(im_size);

	%Starting point in padded image, zero based.
	start = int32(pad);

	a = img_fun('padarray',a,pad,boundary,'both');      % Pad image

	%Create connectivity matrix.  Only use nonzero values of the filter.
	conn_logical = h~=0;
	conn = double( conn_logical );  %input to the mex file must be double
	nonzero_h = h(conn_logical);
	b = imfilter_mex(a,im_size,h,nonzero_h,conn,start,flags);

% ------------------------------------------------------------------------
function push_apply_and_return_CB(hObject, handles)
% Apply filter, and return filtered image and grid (if it exists).
% Notice that this button was made visible only when this option may apply.
	[handles.output.grd, handles.output.img] = push_apply_CB(hObject, handles);

	if (~isempty(handles.h_calling_img))            % Update the caller figure image
		set(handles.h_calling_img,'CData',handles.output.img)
	end

	if (exist('OCTAVE_VERSION','builtin'))		% To know if we are running under Octave
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	else
		do_uiresume = strcmp(get(handles.figure1, 'waitstatus'), 'waiting');
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		guidata(handles.figure1, handles);	uiresume(handles.figure1);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

% ------------------------------------------------------------------------
function h = filts_coef(varargin)
	type = varargin{1};
	siz = (varargin{2}-1)/2;

	switch type
		case 'GaussLowPass'         % Gaussian filter
			std = varargin{3};
			[x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
			h = exp(-(x.*x + y.*y)/(2*std*std));
			h(h < eps*max(h(:))) = 0;
		case 'invDist'              % Inverse distance filter
			pow = varargin{3};
			[x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
			h = sqrt(x.*x + y.*y);
			h(siz(1)+1,siz(2)+1) = 1;   % To avoid devide by zero warnings
			h = 1 ./ (h .^ pow);
			h(siz(1)+1,siz(2)+1) = 2;   % Central weight
	end

% ------------------------------------------------------------------------
function push_cancel_CB(hObject, handles)
	if (exist('OCTAVE_VERSION','builtin'))		% To know if we are running under Octave
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	else
		do_uiresume = strcmp(get(handles.figure1, 'waitstatus'), 'waiting');
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		handles.output.cancel = 1;      % User gave up, return nothing
		set(handles.h_calling_img,'CData',handles.img_orig)     % Reset calling fig original's image
		guidata(handles.figure1, handles);	uiresume(handles.figure1);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

% ------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	push_cancel_CB(hObject, handles)

% ------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles = guidata(hObject);
		push_cancel_CB(hObject, handles)
	end

% ------------------------------------------------------------------------
function push_saveFiltImg_CB(hObject, handles)
	str1 = {'*.bmp', 'Windows Bitmap (*.bmp)'; ...
		'*.hdf', 'Hieralchical Data Format (*.hdf)'; ...
		'*.jpg', 'JPEG image (*.jpg)'; ...
		'*.pcx', 'Windows Paintbrush (*.pcx)'; ...
		'*.png', 'Portable Network Graphics(*.png)'; ...
		'*.ras', 'SUN rasterfile (*.ras)'; ...
		'*.tif', 'Tagged Image File (*.tif)'; ...
		'*.xwd', 'X Windows Dump (*.xwd)'};
	[FileName,PathName] = put_or_get_file(handles, str1,'Select image format','put');
	if isequal(FileName,0),		return,		end

	[PATH,FNAME,EXT] = fileparts([PathName FileName]);
	if (isempty(EXT) && ~isempty(FileName))
		msgbox('Sorry, but you have to give the filename extention','Error'); return
	end
	set(handles.figure1,'pointer','watch')
	img = get(findobj(handles.axes3,'Type','image'),'CData');
	if strcmp(EXT,'.jpg') || strcmp(EXT,'.JPG') || strcmp(EXT,'.jpeg') || strcmp(EXT,'.JPEG')
		if isa(img, 'uint8')
			if (ndims(img) == 2)
				imwrite(img,get(handles.figure1,'Colormap'),[PathName FileName],'Quality',100);
			else    % RGB and colormap is forbiden by jpeg norm
				imwrite(img,[PathName FileName],'Quality',100);
			end
		else
			imwrite(img,[PathName FileName],'Quality',100);
		end
	else        % All other image formats
		if isa(img, 'uint8')
			try
				imwrite(img,get(handles.figure1,'Colormap'),[PathName FileName]);
			catch       % For example RGB images canot be saved as pcx
				msgbox('Format not supported for this image','Warning');    return
			end
		else
			try
				imwrite(img,[PathName FileName]);
			catch
				msgbox('Format not supported for this image','Warning');    return
			end
		end
	end
	set(handles.figure1,'pointer','arrow')

% ------------------------------------------------------------------------
function push_loadFilter_CB(hObject, handles)
% Load a filter from an external file. It must have a .dat extension

str1 = {'*.dat;*.DAT', 'Data file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
[FileName,PathName] = put_or_get_file(handles,str1,'Select filter file','get');
if isequal(FileName,0),		return,		end

try
	filt = text_read([PathName FileName]);
catch
    errordlg(lasterr,'ERROR');    return
end

handles.custom_filt = filt;					% Store the custom filter here.
update_table(handles,handles.custom_filt)		% Update the table
guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function radio_filteredImg_CB(hObject, handles)
	if ( ~get(hObject,'Val')),	set(hObject,'Val', 1),		return,		end
	set(handles.radio_enhancedImg, 'Val', 0)

% -----------------------------------------------------------------------------------------
function radio_enhancedImg_CB(hObject, handles)
	if ( ~get(hObject,'Val')),	set(hObject,'Val', 1),		return,		end
	set(handles.radio_filteredImg, 'Val', 0)

% ------------------------------------------------------------------------
function defFilters_struct = load_defFilters
% Define a structure called LinConvFilt with the filter coefficients.
% To add more fields (to host more filters) see how LinConvFilt is
% created and proceed. After don't forget to update the "nomes" and
% "texto" cell arrays. The "texto" contains the filter descriptive that
% is shown at the bottom of DigitalFiltering_gui figure. If you
% don't have any descriptive text you must still provide a empty entry.
% See comments before the cell arrays initialization.

LinConvFilt.UserDefFilt.LowPass.MovAverage = ones(5,5);
%LinConvFilt.UserDefFilt.LowPass.DistWeight = eye(5,5);
[x,y] = meshgrid(-2:2,-2:2);
h = sqrt(x.*x + y.*y);
h(3,3) = 1;             % To avoid devide by zero warnings
h = 1 ./ h;
h(3,3) = 2;             % Central weight
LinConvFilt.UserDefFilt.LowPass.InvDist = h;
h     = exp(-(x.*x + y.*y)/(2*1*1));
h(h < eps*max(h(:))) = 0;
LinConvFilt.UserDefFilt.LowPass.GaussLowPass = h;
tmp = zeros(5);     tmp(3,3) = 1;
LinConvFilt.UserDefFilt.GeneralUD = tmp;

% ------------------------

LinConvFilt.PreDefFilters.LowPassG.Gauss = [1 2 1; 2 4 2; 1 2 1];
LinConvFilt.PreDefFilters.LowPassG.FiveNodeP = [0 1 0; 1 1 1; 0 1 0];
LinConvFilt.PreDefFilters.LowPassG.FiveNodeT = [1 0 1; 0 1 0; 1 0 1];
LinConvFilt.PreDefFilters.LowPassG.NineNode = ones(3,3);
LinConvFilt.PreDefFilters.LowPassG.LowPass1 = [1 1 1; 1 2 1; 1 1 1];
LinConvFilt.PreDefFilters.LowPassG.LowPass2 = [1 1 1; 1 4 1; 1 1 1];
LinConvFilt.PreDefFilters.LowPassG.LowPass3 = [1 1 1; 1 12 1; 1 1 1];

LinConvFilt.PreDefFilters.HighPass.MeanRemov = [-1 -1 -1; -1 9 -1; -1 -1 -1];
LinConvFilt.PreDefFilters.HighPass.HighPass1 = [0 -1 0; -1 5 -1; 0 -1 0];
LinConvFilt.PreDefFilters.HighPass.HighPass2 = [1 -2 1; -2 5 -2; 1 -2 1];
LinConvFilt.PreDefFilters.HighPass.HighPass3 = [0 -1 0; -1 20 -1; 0 -1 0];

LinConvFilt.PreDefFilters.Order1Der.RobRowDetect = [-1 0 0; 0 1 0; 0 0 0];
LinConvFilt.PreDefFilters.Order1Der.RobColDetect = [0 0 -1; 0 1 0; 0 0 0];
LinConvFilt.PreDefFilters.Order1Der.PreRowDetect = [1 0 -1; 1 0 -1; 1 0 -1];
LinConvFilt.PreDefFilters.Order1Der.PreColDetect = [-1 -1 -1; 0 0 0; 1 1 1];
LinConvFilt.PreDefFilters.Order1Der.SobRowDetect = [1 0 -1; 2 0 -2; 1 0 -1];
LinConvFilt.PreDefFilters.Order1Der.SobColDetect = [-1 -2 -1; 0 0 0; 1 2 1];
LinConvFilt.PreDefFilters.Order1Der.FreiRowDetect = [1 0 -1; sqrt(2) 0 -sqrt(2); 1 0 -1];
LinConvFilt.PreDefFilters.Order1Der.FreiColDetect = [-1 -sqrt(2) -1; 0 0 0; 1 sqrt(2) 1];

LinConvFilt.PreDefFilters.Order2Der.Lapla1 = [0 -1 0; -1 4 -1; 0 -1 0];
LinConvFilt.PreDefFilters.Order2Der.Lapla2 = [-1 -1 -1; -1 8 -1; -1 -1 -1];
LinConvFilt.PreDefFilters.Order2Der.Lapla3 = [1 -2 1; -2 4 -2; 1 -2 1];
LinConvFilt.PreDefFilters.Order2Der.Lapla4 = [-1 0 -1; 0 4 0; -1 0 -1];
LinConvFilt.PreDefFilters.Order2Der.LapDif = [0 -1 0; -1 5 -1; 0 -1 0];
LinConvFilt.PreDefFilters.Order2Der.DifGauss = [0 0 -1 -1 -1 0 0; 0 -2 -3 -3 -3 -2 0; -1 -3 5 5 5 -3 -1;...
     -1 -3 5 16 5 -3 -1; -1 -3 5 5 5 -3 -1; 0 -2 -3 -3 -3 -2 0; 0 0 -1 -1 -1 0 0];

LinConvFilt.PreDefFilters.ShiftDif.Horiz = [0 -1 0; 0 1 0; 0 0 0];
LinConvFilt.PreDefFilters.ShiftDif.Vert = [0 0 0; -1 1 0; 0 0 0];

LinConvFilt.PreDefFilters.GradDir.G_E = [-1 1 1; -1 -2 1; -1 1 1];
LinConvFilt.PreDefFilters.GradDir.G_SE = [-1 -1 1; -1 -2 1; 1 1 1];
LinConvFilt.PreDefFilters.GradDir.G_S = [-1 -1 -1; 1 -2 1; 1 1 1];
LinConvFilt.PreDefFilters.GradDir.G_SW = [1 -1 -1; 1 -2 -1; 1 1 1];
LinConvFilt.PreDefFilters.GradDir.G_W = [1 1 -1; 1 -2 -1; 1 1 -1];
LinConvFilt.PreDefFilters.GradDir.G_NW = [1 1 1; 1 -2 -1; 1 -1 -1];
LinConvFilt.PreDefFilters.GradDir.G_N = [1 1 1; 1 -2 1; -1 -1 -1];
LinConvFilt.PreDefFilters.GradDir.G_NE = [1 1 1; -1 -2 1; -1 -1 1];

LinConvFilt.PreDefFilters.Emboss.E_E = [-1 0 1; -1 1 1; -1 0 1];
LinConvFilt.PreDefFilters.Emboss.E_SE = [-1 -1 0; -1 1 1; 0 1 1];
LinConvFilt.PreDefFilters.Emboss.E_S = [-1 -1 -1; 0 1 0; 1 1 1];
LinConvFilt.PreDefFilters.Emboss.E_SW = [0 -1 -1; 1 1 -1; 1 1 0];
LinConvFilt.PreDefFilters.Emboss.E_W = [1 0 -1; 1 1 -1; 1 0 -1];
LinConvFilt.PreDefFilters.Emboss.E_NW = [1 1 0; 1 1 -1; 0 -1 -1];
LinConvFilt.PreDefFilters.Emboss.E_N = [1 1 1; 0 1 0; -1 -1 -1];
LinConvFilt.PreDefFilters.Emboss.E_NE = [0 1 1; -1 1 1; -1 -1 0];

% NonLinFilt.stdfilt = ones(3);

% First column contains the struct field name, second the name that will
% show up in the tree. Notice how the first column cuincides with the
% LinConvFilt structure fieldnames. This aspect is crutial.
nomes = {'Linear Convolution filters' 'Linear Convolution filters';...  % The first must always repeat
		'UserDefFilt' 'User defined'; ...
		'LowPass' 'Low-pass filters'; ...
		'MovAverage' 'Moving Average (mxn)'; ...
		'InvDist' 'Inverse Distance (mxn)';...
		'GaussLowPass' 'Gaussian Low-pass (mxn)';...
		'GeneralUD' 'General User-defined (mxn)';...
		'PreDefFilters' 'Predefined Filters';...
		'LowPassG' 'Low-pass Filters';...
		'Gauss' 'Gaussian (3x3)';...
		'FiveNodeP' '5-node + Averaging (3x3)';...
		'FiveNodeT' '5-node X Averaging (3x3)';...
		'NineNode' '9-node Averaging (3x3)';...
		'LowPass1' 'Low-pass 1 (3x3)';...
		'LowPass2' 'Low-pass 2 (3x3)';...
		'LowPass3' 'Low-pass 3 (3x3)';...
		'HighPass' 'High-pass Filters';...
		'MeanRemov' 'Mean Removal (3x3)';...
		'HighPass1' 'High-pass 1 (3x3)';...
		'HighPass2' 'High-pass 2 (3x3)';...
		'HighPass3' 'High-pass 3 (3x3)';...
		'Order1Der' '1 Order Derivative Filters';...
		'RobRowDetect' 'Roberts Row Detector 3 (3x3)';...
		'RobColDetect' 'Roberts Col Detector 3 (3x3)';...
		'PreRowDetect' 'Prewitt Row Detector 3 (3x3)';...
		'PreColDetect' 'Prewitt Col Detector 3 (3x3)';...
		'SobRowDetect' 'Sobel Row Detector 3 (3x3)';...
		'SobColDetect' 'Sobel Col Detector 3 (3x3)';...
		'FreiRowDetect' 'FreiChen Row Detector 3 (3x3)';...
		'FreiColDetect' 'FreiChen Col Detector 3 (3x3)';...
		'Order2Der' '2 Order Derivative Filters';...
		'Lapla1' 'Laplacian 1 (3x3)';...
		'Lapla2' 'Laplacian 2 (3x3)';...
		'Lapla3' 'Laplacian 3 (3x3)';...
		'Lapla4' 'Laplacian 4 (3x3)';...
		'LapDif' 'Laplacian difference (3x3)';...
		'DifGauss' 'Difference of Gaussian (7x7)';...
		'ShiftDif' 'Shift and Difference Filters';...
		'Horiz' 'Horizontal (3x3)';...
		'Vert' 'Vertical (3x3)';...
		'GradDir' 'Gradient Directional Filters';...
		'G_E' 'East (3x3)';...
		'G_SE' 'Southeast (3x3)';...
		'G_S' 'South (3x3)';...
		'G_SW' 'Southwest (3x3)';...
		'G_W' 'West (3x3)';...
		'G_NW' 'Northwest (3x3)';...
		'G_N' 'North (3x3)';...
		'G_NE' 'Northeast (3x3)';...
		'Emboss' 'Embossing Filters';...
		'E_E' 'East (3x3)';...
		'E_SE' 'Southeast (3x3)';...
		'E_S' 'South (3x3)';...
		'E_SW' 'Southwest (3x3)';...
		'E_W' 'West (3x3)';...
		'E_NW' 'Northwest (3x3)';...
		'E_N' 'North (3x3)';...
		'E_NE' 'Northeast (3x3)';...
%         'Nonlinear Filters' 'Nonlinear Filters';...
%         'stdfilt' 'Standard Deviation (mxn)';...
		};

% Descriptive text that will show in a textbox. The number of elements must
% be exactly equal to 'nomes'
texto = {'Output nodes are weighted sums of neighboring nodal values.';...
		'User defined, variale size, rectangular neighborhoods.';...
		'Block out higher frequencies and reduce noise.';...
		'Moving average in a rectangular neighborhood.';...
		'Iso-weight contours are concentric ellipses.';...
		'Weights form a 2D bell-shaped curve.';...
		'User defined linear convolution filter.';...
		'Fixed size and predefined weights.';...
		'Block out higher frequencies and reduce noise.';...
		'Small version of the gaussian low-pass filter.';...
		'N-S-E-W moving average.';...
		'X pattern moving average.';...
		'Moving average in small square neighborhood.';...
		'';...
		'';...
		'';...
		'Amplify higher frequencies.';...
		'Remove local mean.';...
		'';...
		'';...
		'';...
		'Common gradient-based filters for finding horizontal and vertical edges.';...
		'';...
		'';...
		'';...
		'';...
		'';...
		'';...
		'';...
		'';...
		'Common second order derivative edge localization filters.';...
		'';...
		'';...
		'';...
		'';...
		'';...
		'';...
		'Subtract a spatially shifted copy of the grid from the original grid.';...
		'Enhance horizontal edges.';...
		'Enhance veertical edges.';...
		'Highlight areas of rapid changes -- i.e. high slopes.';...
		''; ''; ''; ''; ''; ''; ''; '';
		'Directional edge enhancement (''embossing'') filters.';...
		''; ''; ''; ''; ''; ''; ''; '';
%		'General functions of neighboring nodes'; ...
%		'Standard deviation of nodes in neighborhood';
	};

defFilters_struct = {'Linear Convolution filters'; LinConvFilt; nomes; texto};

% ------------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function digitalFiltering_LayoutFcn(h1)
set(h1,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Digital Filtering',...
'NumberTitle','off',...
'Position',[520 370 800 420],...
'DoubleBuffer','on',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[14 9 227 46],'Style','frame','Tag','frame2');

h3 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@digitalFiltering_uiCB,...
'FontSize',9,...
'Position',[10 69 251 335],...
'Style','listbox',...
'Value',1,...
'Tag','listbox1');

uicontrol('Parent',h1,'FontSize',10,'Position',[31 404 199 16],...
'String','Filters','Style','text','FontName','Helvetica');

h6 = uicontextmenu('Parent',h1,...
'Call',@digitalFiltering_uiCB,...
'Tag','plot_select_menu');

set(h3,'uicontextmenu',h6)

uimenu('Parent',h6,...
'Call',@digitalFiltering_uiCB,...
'Label','plot selected',...
'Tag','plot_selected_menu');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[20 13 211 40],...
'String','Nikles',...
'Style','text',...
'Tag','text_info');

axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Position',[270 70 271 141],...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'Tag','tblParams');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[270 240 34 15],...
'String','Rows:',...
'FontName','Helvetica',...
'Style','text',...
'Tag','text_rows');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[430 241 34 15],...
'String','Cols:',...
'FontName','Helvetica',...
'Style','text',...
'Tag','text_cols');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@digitalFiltering_uiCB,...
'Position',[310 219 75 51],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_nRows');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@digitalFiltering_uiCB,...
'Position',[466 219 75 51],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_nCols');

axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'Color',get(0,'defaultaxesColor'),...
'Position',[550 229 241 171],...
'XLim',get(0,'defaultaxesXLim'),...
'XLimMode','manual',...
'XTick',[],...
'XTickMode','manual',...
'YLim',get(0,'defaultaxesYLim'),...
'YLimMode','manual',...
'YTick',[],...
'YTickMode','manual',...
'Tag','axes2');

axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'Color',get(0,'defaultaxesColor'),...
'Position',[550 20 241 171],...
'XLim',get(0,'defaultaxesXLim'),...
'XLimMode','manual',...
'XTick',[],...
'XTickMode','manual',...
'YLim',get(0,'defaultaxesYLim'),...
'YLimMode','manual',...
'YTick',[],...
'YTickMode','manual',...
'Tag','axes3');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@digitalFiltering_uiCB,...
'HorizontalAlignment','left',...
'Position',[270 350 251 21],...
'Style','edit',...
'Tooltip','Load an Image',...
'Tag','edit_imageIn');

uicontrol('Parent',h1,...
'Call',@digitalFiltering_uiCB,...
'FontSize',12,...
'FontWeight','bold',...
'Position',[521 348 20 23],...
'String','...',...
'Tooltip','Browse for an image',...
'Tag','push_imageIn');

uicontrol('Parent',h1,'Position',[634 405 67 15],...
'String','Original Image','Style','text','FontName','Helvetica');

uicontrol('Parent',h1,'Position',[636 196 70 15],...
'String','Filtered Image','Style','text','FontName','Helvetica');

uicontrol('Parent',h1,...
'Call',@digitalFiltering_uiCB,...
'Position',[270 289 125 21],...
'String','Load external filter',...
'FontName','Helvetica',...
'Tooltip','Import an external (MxN) filter',...
'Tag','push_loadFilter');

uicontrol('Parent',h1,...
'Call',@digitalFiltering_uiCB,...
'Position',[415 289 125 21],...
'String','Save filtered image',...
'FontName','Helvetica',...
'Tag','push_saveFiltImg');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[270 375 80 15],...
'String','Load an image',...
'FontName','Helvetica',...
'Style','text',...
'Tag','text_loadImage');

uicontrol('Parent',h1,...
'Call',@digitalFiltering_uiCB,...
'FontName','Helvetica',...
'Position',[266 41 135 19],...
'String','Show filtered image',...
'Style','radiobutton',...
'Tooltip','Display only the result of filtering the image',...
'Value',1,...
'Tag','radio_filteredImg');

uicontrol('Parent',h1,...
'Call',@digitalFiltering_uiCB,...
'FontName','Helvetica',...
'Position',[403 41 145 19],...
'String','Show enhanced image',...
'Style','radiobutton',...
'Tooltip','Enhance image by applying to it the result of filtering',...
'Value',0,...
'Tag','radio_enhancedImg');

uicontrol('Parent',h1,...
'Call',@digitalFiltering_uiCB,...
'Position',[270 10 50 21],...
'String','Apply',...
'FontName','Helvetica',...
'Tooltip','Filter image (or grid)',...
'Tag','push_apply');

uicontrol('Parent',h1,...
'Call',@digitalFiltering_uiCB,...
'Position',[350 10 80 21],...
'String','Apply n return',...
'FontName','Helvetica',...
'Tooltip','Stop and return filtered image and grid (if required)',...
'Tag','push_apply_and_return');

uicontrol('Parent',h1,...
'Call',@digitalFiltering_uiCB,...
'Position',[481 9 60 21],...
'String','Cancel',...
'FontName','Helvetica',...
'Tooltip','Stop and undo any previous changes',...
'Tag','push_cancel');

function digitalFiltering_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

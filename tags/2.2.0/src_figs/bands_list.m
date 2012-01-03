function varargout = bands_list(varargin)
% Compose a false color image by band selection

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

	if (numel(varargin) == 0),		return,		end
 
	hObject = figure('Vis','off');
	bands_list_LayoutFcn(hObject);
	handles = guihandles(hObject);

	handles.Rband = [];		% To hold the band number that will be puted here
	handles.Gband = [];		%               "
	handles.Bband = [];		%               "
	handles.frame_movel_pos = get(handles.frame_movel,'Pos');
	handles.edit_Rband_pos = get(handles.edit_Rband,'Pos');

	handles.hMirFig = varargin{1};
	bandList = getappdata(varargin{1},'BandList');
	if (isempty(bandList))
		errordlg('ERROR: There is no list of bands available. Do you know what you are doing?','ERROR')
		delete(hObject);    return
	end
	move2side(handles.hMirFig,hObject,'left')	% Reposition this figure
	handles.struct_names = bandList(1);
	set(handles.listbox1,'String',bandList(1))
	handles.image_bands = bandList{2};	% A MxNxP uint8 array, where P is the number of bands in memory (~= ntotal bands)
	handles.all_names = bandList{3};	% A Mx2 cell array with struct field names & names to show up in the tree
	handles.band_desc = bandList{4};	% A Mx2 cell array with the data description and band number (per struct field)
	handles.fname = bandList{5};		% File name
	handles.bands_inMemory = bandList{6};    % A vector with the band numbers already in memory
	handles.dims = bandList{7};			% A [n_row n_col nBands_total] vector with the 2D image dimensions and the
										% TOTAL number of bands in the dataset (some of them may not be on memory)
	handles.reader = bandList{8};		% A string with 'GDAL' if the dataset was openend with it OR:
										% a Mx2 cell array with {byte_resolution, header, interleave, endian};
										% and a 1 to 3 elements with the subsets option of multiband read.

	for (i = 1:handles.dims(3))
		tmp.(sprintf('band%d',i)) = i;	
	end   
	handles.struct_values = {tmp};		% TENHO DE MUDAR ESTE NOME

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) ------------------------------

	% Add this figure handle to the carra?as list
	plugedWin = getappdata(handles.hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hMirFig,'dependentFigs',plugedWin);

	handMir = guidata(handles.hMirFig);        % Retrive Mirone handles
	handles.image_type_orig = handMir.image_type;
	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% ------------------------------------------------------------------------
function radio_gray_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.radio_RGB,'Value',0)
		pos = handles.frame_movel_pos;
		pos = [pos(1) pos(2)+pos(4)/2 pos(3) pos(4)/2];
		set(handles.edit_Rband,'Pos',pos+[10 5 -20 -28])
		set(handles.radio_R,'Vis','off')
		set(handles.radio_G,'Vis','off')
		set(handles.radio_B,'Vis','off')
		set(handles.edit_Gband,'Vis','off')
		set(handles.edit_Bband,'Vis','off')
		set(handles.text_toGray,'Vis','on')
	else
		set(hObject,'Value',1)
	end

% ------------------------------------------------------------------------
function radio_RGB_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.radio_gray,'Value',0)
		set(handles.edit_Rband,'Pos',handles.edit_Rband_pos)
		set(handles.radio_R,'Vis','on')
		set(handles.radio_G,'Vis','on')
		set(handles.radio_B,'Vis','on')
		set(handles.edit_Gband,'Vis','on')
		set(handles.edit_Bband,'Vis','on')
		set(handles.text_toGray,'Vis','off')
	else
		set(hObject,'Value',1)
	end

% ------------------------------------------------------------------------
function radio_R_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.radio_G,'Value',0)
		set(handles.radio_B,'Value',0)
	else
		set(hObject,'Value',1)
	end

% ------------------------------------------------------------------------
function radio_G_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.radio_R,'Value',0)
		set(handles.radio_B,'Value',0)
	else
		set(hObject,'Value',1)
	end

% ------------------------------------------------------------------------
function radio_B_CB(hObject, handles)
if (get(hObject,'Value'))
    set(handles.radio_R,'Value',0)
    set(handles.radio_G,'Value',0)
else
    set(hObject,'Value',1)
end

% --------------------------------------------------------------------------
function push_pca_CB(hObject, handles)
	[m,n,k] = size(handles.image_bands);
	q = min(k, 6);					% for memory sake only 6 components are computed
	P = princomp(handles.image_bands, q, true);
	lamb = diag(P.Cy);				% The eigenvalues
	P = reshape(P.Y, m, n, q);
% 	P = reshape( fastica(reshape(handles.image_bands, m * n, k)', 'lastEig',q)', m, n, q );		lamb = ones(q,1);
	P8 = alloc_mex(m,n,q,'uint8');
	for (i = 1:q)
		P8(:,:,i) = scaleto8(P(:,:,i),-8);		% Each component must be scaled independently of the others
	end
	hFig = mirone;		handNewMir = guidata(hFig);		set(hFig,'Name', 'PCA image')
	handMir = guidata(handles.hMirFig);        % Retrive Mirone handles
	handNewMir.head = handMir.head;
	handNewMir.geog = handMir.geog;
	handNewMir.image_type = handMir.image_type;
	mirone('FileOpenGDALmultiBand_CB', handNewMir, 'PCA image', P8)
	% Rename bands according to their varianves (eigenvalues)
	tmp = getappdata(hFig,'BandList');
	for (i = 1:q)
		tmp{3}(i+1, 2)={sprintf('Var(%d)_%.0f',i,lamb(i))};
	end
	setappdata(hFig,'BandList',tmp)

% --------------------------------------------------------------------------
function push_Load_CB(hObject, handles)
% TENHO QUE TESTAR SE CASO FOR PRECISO LOADAR BANDAS ELAS SEJAM GDALICAS

if (get(handles.radio_RGB,'Value') && ...
        (isempty(handles.Rband) || isempty(handles.Gband) || isempty(handles.Bband)))
    errordlg('Error: you must select three bands','ERROR')
    return
end
if ( get(handles.radio_gray,'Value') && (isempty(handles.Rband)) )
    errordlg('Error: Are you blind? You must select a band. Kinda obvious no?','ERROR')
    return
end

handMir = guidata(handles.hMirFig);        % Retrive Mirone handles
img = get(handMir.hImg,'CData');
head = [];              % It will be changed only if we load a composition of non uint8

if (get(handles.radio_RGB,'Value'))       % RGB - pure image for sure (is it??)
    [b,b] = ismember([handles.Rband handles.Gband handles.Bband],handles.bands_inMemory);
    idx = (b == 0);		% ISMEMBER returns zeros for elements of A not in B
    b(idx) = [];		% If they exis, clear them
    if (length(b) == 3)         % All bands are aready in memory
        img = handles.image_bands(:,:,b);
    elseif (length(b) == 2)     % Two bands are in memory. Need to load the third
        id = intersect([handles.Rband handles.Gband handles.Bband], b);
        if (isequal(id,[1 2]) || isequal(id,[2 1]))          % Red & Green in memory
            img(:,:,1:2) = handles.image_bands(:,:,b);
            img(:,:,3) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Bband)]);
        elseif (isequal(id,[1 3]) || isequal(id,[3 1]))      % Red & Blue in memory
            img(:,:,1) = handles.image_bands(:,:,b(1));
            img(:,:,2) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Gband)]);
            img(:,:,3) = handles.image_bands(:,:,b(2));
        elseif (isequal(id,[2 3]) || isequal(id,[3 2]))      % Green & Blue in memory
            img(:,:,1) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Rband)]);
            img(:,:,2:3) = handles.image_bands(:,:,b);
        end
    elseif (length(b) == 1)     % One band in memory. Need to load the other two
        id = find([handles.Rband handles.Gband handles.Bband] == b);
        if (id == 1)            % In memory band is the Red layer
            img(:,:,1) = handles.image_bands(:,:,b);
            if (handles.Gband < handles.Bband)  % OK, band numbers are in ascending order
                img(:,:,2:3) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Gband) ',' num2str(handles.Bband)]);
            else                % gdalread cannot read e.g. -B2,1 so we have to make two calls
                img(:,:,2) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Gband)]);
                img(:,:,3) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Bband)]);
            end
        elseif (id == 2)        % In memory band is the Green layer
            img(:,:,1) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Rband)]);
            img(:,:,2) = handles.image_bands(:,:,b);
            img(:,:,3) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Bband)]);
        else                    % In memory band is the Blue layer
            if (handles.Rband < handles.Gband)  % OK, band numbers are in ascending order
                img(:,:,1:2) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Rband) ',' num2str(handles.Gband)]);
            else
                img(:,:,1) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Rband)]);
                img(:,:,2) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Gband)]);
            end
            img(:,:,3) = handles.image_bands(:,:,b);
        end
    else                        % No bands in memory. Need to load them all
        if (handles.Rband < handles.Gband && handles.Gband < handles.Bband)
            img = gdalread(handles.fname,'-S', ['-B' num2str(handles.Rband) ',' ...
                    num2str(handles.Gband) ',' num2str(handles.Bband)]);
        else                % Just read one band at a time
            img(:,:,1) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Rband)]);
            img(:,:,2) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Gband)]);
            img(:,:,3) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Bband)]);
        end
    end

    set(handMir.hImg,'CData',img)
    try
		rmappdata(handles.hMirFig,'dem_x');		rmappdata(handles.hMirFig,'dem_y');		rmappdata(handles.hMirFig,'dem_z');
    end
    image_type = handles.image_type_orig;		% Reset indicator that this is an image only
    computed_grid = 0;      % Reset this also
    was_int16 = 0;
else                        % GRAY SCALE, which can be an image or a > uint8 image band that needs scaling
    %b = intersect(handles.Rband,handles.bands_inMemory);
    [b,b] = ismember(handles.Rband,handles.bands_inMemory);
    idx = (b == 0);         % ISMEMBER returns zeros for elements of A not in B
    b(idx) = [];            % If they exis, clear them
    if (~isempty(b))        % The band is on memory
        img = handles.image_bands(:,:,handles.Rband);
    else                    % Need to load band
        % Here it will fail if the file is not to be read by GDAL
        img = gdalread(handles.fname,'-S', ['-B' num2str(handles.Rband)]);
%         img = gdalread(handles.fname, ['-B' num2str(handles.Rband)]);
%         img = scaleto8(img,8,0);        % Need to give the noDataValue
    end

    if (isa(img,'uint8') || islogical(img))
        image_type = handles.image_type_orig;		% Reset indicator that this is an image only
        computed_grid = 0;							% Reset this also
        was_int16 = 0;

    else                        % Not uint8, so we need scalings
        % Now we are going to load the band not scaled and treat it as a GMT grid
        if (~iscell(handles.reader))    % reader is GDAL
            Z = gdalread(handles.fname, ['-B' num2str(handles.Rband)]);
        else                            % reader is MULTIBANDREAD - SHIT, WHAT A MESS OF TESTS WE NEED TO DO
			if (length(handles.reader(2)) >= 2)     % Two (or three) subsets have been choosed. Ignore third one
                if (strcmp(handles.reader{2}(1), 'Row') && strcmp(handles.reader{2}(2), 'Column'))      % Row & Column
                    Z = multibandread_j(handles.fname, [handles.dims(1) handles.dims(2) handles.dims(3)],...
                        handles.reader{1}{1},handles.reader{1}{2},handles.reader{1}{3},handles.reader{1}{4},...
                        handles.reader{2}(1),handles.reader{2}(2),{'Band','Direct',handles.Rband});
                elseif (strcmp(handles.reader{2}(1), 'Row') && ~strcmp(handles.reader{2}(2), 'Column')) % Row only
                    Z = multibandread_j(handles.fname, [handles.dims(1) handles.dims(2) handles.dims(3)],...
                        handles.reader{1}{1},handles.reader{1}{2},handles.reader{1}{3},handles.reader{1}{4},...
                        handles.reader{2}(1),{'Band','Direct',handles.Rband});
                elseif (~strcmp(handles.reader{2}(1), 'Row') && strcmp(handles.reader{2}(2), 'Column')) % Column only
                    Z = multibandread_j(handles.fname, [handles.dims(1) handles.dims(2) handles.dims(3)],...
                        handles.reader{1}{1},handles.reader{1}{2},handles.reader{1}{3},handles.reader{1}{4},...
                        handles.reader{2}(2),{'Band','Direct',handles.Rband});
                end
			else                                    % One subsect selectet. Ignore it if it was the 'Band'
                if (strcmp(handles.reader{2}(1), 'Row'))            % Row only
                    Z = multibandread_j(handles.fname, [handles.dims(1) handles.dims(2) handles.dims(3)],...
                        handles.reader{1}{1},handles.reader{1}{2},handles.reader{1}{3},handles.reader{1}{4},...
                        handles.reader{2}(1),{'Band','Direct',handles.Rband});
                elseif (strcmp(handles.reader{2}(1), 'Column'))     % Column only
                    Z = multibandread_j(handles.fname, [handles.dims(1) handles.dims(2) handles.dims(3)],...
                        handles.reader{1}{1},handles.reader{1}{2},handles.reader{1}{3},handles.reader{1}{4},...
                        handles.reader{2}(2),{'Band','Direct',handles.Rband});
                else                                                % Neither Row or Column
                    Z = multibandread_j(handles.fname, [handles.dims(1) handles.dims(2) handles.dims(3)],...
                        handles.reader{1}{1},handles.reader{1}{2},handles.reader{1}{3},handles.reader{1}{4},...
                        {'Band','Direct',handles.Rband});
                end
			end
        end
        X = 1:handles.dims(2);    Y = 1:handles.dims(1);
        head = handMir.head;
        head(5:6) = [double(min(min(Z))) double(max(max(Z)))];
        setappdata(handles.hMirFig,'dem_z',Z);  setappdata(handles.hMirFig,'dem_x',X);
        setappdata(handles.hMirFig,'dem_y',Y);
        image_type = 1;         % Pretend this a GMT grid
        computed_grid = 1;      % But set to computed_grid to avoid attempts to reload it with grdread_m
        if (isa(Z,'uint16') || isa(Z,'int16'))
            was_int16 = 1;
        else
            was_int16 = 0;
        end
    
    end         % end if is uint8
    
	set(handMir.hImg,'CData',img)
	set(handles.hMirFig,'ColorMap',gray(256))
end

if (~isempty(head)),    handMir.head = head;    end
handMir.image_type = image_type;
handMir.computed_grid = computed_grid;
handMir.was_int16 = was_int16;
guidata(handles.hMirFig,handMir)           % Save those in Mirone handles

% --------------------------------------------------------------------------
function listbox1_CB(hObject, handles)
index_struct = get(hObject,'Value');
struct_names = handles.struct_names;
struct_values = handles.struct_values;

indent = '       ';
root_1 = struct_names{index_struct};
is_indent = strfind(root_1, indent);
if (isempty(is_indent)),	level = 0;
else						level = (is_indent(end) - 1)/7 + 1;
end
    
struct_val = struct_values{index_struct};
all_names = handles.all_names;

if isa(struct_val,'struct')
    fields =  fieldnames(struct_val);
    %names = cell(1,length(fields));
    for i = 1:length(fields)
        %idx = strmatch(fields{i},all_names(:,1),'exact');
        %names{i} = all_names{idx,2};
        if isstruct(getfield(struct_val(1), fields{i}))
            fields{i} = ['+ ' fields{i}];
            %names{i} = ['+ ' names{i}];            
        end
    end
end

% Display info
name_clean = ddewhite(struct_names{index_struct});
if (name_clean(1) == '+' || name_clean(1) == '-'),       name_clean = name_clean(3:end);     end
idx = strmatch(name_clean,all_names(:,1),'exact');
%handles.tree_indice = idx;
set(handles.edit_dimsDesc,'String',handles.band_desc{idx,1})

if (isnumeric(struct_val))   
	handles = order_bands(handles,idx);
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
	else
		if (get(handles.radio_gray,'Value'))
			guidata(hObject, handles);
			push_Load_CB(handles.push_Load, handles)
		end
		return
	end
	names = cell(length(struct_names),1);
	for (i = 1:length(struct_names))
		name_clean = ddewhite(struct_names{i});
		if (name_clean(1) == '+' || name_clean(1) == '-'),       name_clean = name_clean(3:end);     end
		id1 = strmatch(name_clean,all_names(:,1),'exact');    % Find index to pretended name
		id2 = findstr(name_clean,struct_names{i});              % Find index of starting text (after the blanks)
		names{i} = [struct_names{i}(1:id2-1) all_names{id1,2}];        
	end
	set(handles.listbox1,'String',names);
	handles.struct_names = struct_names;
	handles.struct_values = struct_values;
end

guidata(hObject, handles);

% ------------------------------------------------------------------------
function cell_array = indent_cell(cell_array, level)
	indent = '       ';             indent_app = [];
	for (k = 1:level+1),            indent_app = [indent_app indent];   end
	for (i=1:length(cell_array)),   cell_array{i} = [indent_app cell_array{i}];     end

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
        values_app{i} = getfield(struct_values{idx}, fields{i});
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
    
% ------------------------------------------------------------------------
function handles = order_bands(handles,idx)
% Put selected band (pointed by idx) on the box corresponding to active radiobutton

	n_band = handles.band_desc{idx,2};
	if (get(handles.radio_RGB,'Value'))
		% put the clicked band on the box that has the active radiobutton
		r = get(handles.radio_R,'Value');
		g = get(handles.radio_G,'Value');
		
		if (r)
			set(handles.edit_Rband,'String',handles.all_names{idx,2})
			set(handles.radio_R,'Value',0)
			set(handles.radio_G,'Value',1)
			handles.Rband = n_band;
		elseif (g)
			set(handles.edit_Gband,'String',handles.all_names{idx,2})
			set(handles.radio_G,'Value',0)
			set(handles.radio_B,'Value',1)
			handles.Gband = n_band;
		else
			set(handles.edit_Bband,'String',handles.all_names{idx,2})
			set(handles.radio_B,'Value',0)
			set(handles.radio_R,'Value',1)
			handles.Bband = n_band;
		end
	else        % Single Band
		set(handles.edit_Rband,'String',handles.all_names{idx,2})
		handles.Rband = n_band;
	end

% ---------------------------------------------------------------------------
function P = princomp(X, q, smaller)
%PRINCOMP Obtain principal-component vectors and related quantities.
%   P = PRINCOMP(X, Q) Computes the principal-component vectors of
%   the vector population contained in the rows of X, a matrix of
%   size K-by-n where K is the number of vectors and n is their
%   dimensionality. Q, with values in the range [0, n], is the number
%   of eigenvectors used in constructing the principal-components
%   transformation matrix. P is a structure with the following fields:
%
%     P.Y      K-by-Q matrix whose columns are the principal-
%              component vectors.
%     P.A      Q-by-n principal components transformation matrix
%              whose rows are the Q eigenvectors of Cx corresponding
%              to the Q largest eigenvalues. 
%     P.X      K-by-n matrix whose rows are the vectors reconstructed
%              from the principal-component vectors. P.X and P.Y are
%              identical if Q = n.
%				NOTE COMPUTED IF NARGIN == 3
%
%     P.ems    The mean square error incurred in using only the Q
%              eigenvectors corresponding to the largest
%              eigenvalues. P.ems is 0 if Q = n.  
%     P.Cx     The n-by-n covariance matrix of the population in X.
%     P.mx     The n-by-1 mean vector of the population in X.
%     P.Cy     The Q-by-Q covariance matrix of the population in
%              Y. The main diagonal contains the eigenvalues (in
%              descending order) corresponding to the Q eigenvectors.

%   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
%   Digital Image Processing Using MATLAB, Prentice-Hall, 2004
%   $Revision: 1.4 $  $Date: 2003/10/26 23:16:16 $


if (nargin == 3),	smaller = true;			% Do not compute P.X
else				smaller = false;
end

X = imstack2vectors(X);

[K, n] = size(X);		was_double = true;
if (~isa(X,'double'))
	X = double(X);		was_double = false;
end

% Obtain the mean vector and covariance matrix of the vectors in X.
[P.Cx, P.mx] = covmatrix(X);
P.mx = P.mx'; % Convert mean vector to a row vector.

% Obtain the eigenvectors and corresponding eigenvalues of Cx.  The
% eigenvectors are the columns of n-by-n matrix V.  D is an n-by-n
% diagonal matrix whose elements along the main diagonal are the
% eigenvalues corresponding to the eigenvectors in V, so that X*V = D*V.
[V, D] = eig(P.Cx);

% Sort the eigenvalues in decreasing order.  Rearrange the eigenvectors to match. 
d = diag(D);
[d, idx] = sort(d);
d = flipud(d);
idx = flipud(idx);
D = diag(d);
V = V(:, idx);

% Now form the q rows of A from first q columns of V.
P.A = V(:, 1:q)';

% Compute the principal component vectors.
Mx = repmat(P.mx, K, 1);	% M-by-n matrix.  Each row = P.mx.
%P.Y = P.A*(X - Mx)';		% q-by-K matrix.
P.Y = P.A*( cvlib_mex('sub', X, Mx) )';		% q-by-K matrix.

if (~smaller)		% Obtain the reconstructed vectors.
	P.X = (P.A'*P.Y)' + Mx;
end

% Save some memory
if (~was_double),	clear X,	end		% The X was a local variable
clear Mx

% Convert P.Y to K-by-q array and P.mx to n-by-1 vector.
P.Y = P.Y';
P.mx = P.mx';

% The mean square error is given by the sum of all the 
% eigenvalues minus the sum of the q largest eigenvalues.
d = diag(D);
P.ems = sum(d(q + 1:end));

% Covariance matrix of the Y's:
P.Cy = P.A * P.Cx * P.A';

% -------------------------------------------------------------------------------
function [C, m] = covmatrix(X)
%COVMATRIX Computes the covariance matrix of a vector population.
%   [C, M] = COVMATRIX(X) computes the covariance matrix C and the
%   mean vector M of a vector population organized as the rows of
%   matrix X. C is of size N-by-N and M is of size N-by-1, where N is
%   the dimension of the vectors (the number of columns of X).

%   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
%   Digital Image Processing Using MATLAB, Prentice-Hall, 2004
%   $Revision: 1.4 $  $Date: 2003/05/19 12:09:06 $

[K, n] = size(X);
X = double(X);
if n == 1		% Handle special case.
	C = 0; 
	m = X;
else			% Compute an unbiased estimate of m.
	m = sum(X, 1)/K;
	% Subtract the mean from each row of X.
	X = X - m(ones(K, 1), :);
	% Compute an unbiased estimate of C. Note that the product is
	% X'*X because the vectors are rows of X.	
	C = (X'*X)/(K - 1);
	m = m'; % Convert to a column vector.	 
end

% -------------------------------------------------------------------------------
function [X, R] = imstack2vectors(S, MASK)
%IMSTACK2VECTORS Extracts vectors from an image stack.
%   [X, R] = imstack2vectors(S, MASK) extracts vectors from S, which
%   is an M-by-N-by-n stack array of n registered images of size
%   M-by-N each (see Fig. 11.24). The extracted vectors are arranged
%   as the rows of array X. Input MASK is an M-by-N logical or
%   numeric image with nonzero values (1s if it is a logical array)
%   in the locations where elements of S are to be used in forming X
%   and 0s in locations to be ignored. The number of row vectors in X
%   is equal to the number of nonzero elements of MASK. If MASK is
%   omitted, all M*N locations are used in forming X.  A simple way
%   to obtain MASK interactively is to use function roipoly. Finally,
%   R is an array whose rows are the 2-D coordinates containing the
%   region locations in MASK from which the vectors in S were
%   extracted to form X.  

%   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
%   Digital Image Processing Using MATLAB, Prentice-Hall, 2004
%   $Revision: 1.6 $  $Date: 2003/11/21 14:37:21 $

% Preliminaries.
[M, N, n] = size(S);
if (nargin == 1),		MASK = true(M, N);
else					MASK = (MASK ~= 0);
end

% Find the set of locations where the vectors will be kept before
% MASK is changed later in the program.
if (nargout == 2)
	[I, J] = find(MASK);
	R = [I, J];
end

% Now find X.

% First reshape S into X by turning each set of n values along the third
% dimension of S so that it becomes a row of X. The order is from top to
% bottom along the first column, the second column, and so on.
Q = M*N;
X = reshape(S, Q, n);

% Now reshape MASK so that it corresponds to the right locations 
% vertically along the elements of X.
MASK = reshape(MASK, Q, 1);

% Keep the rows of X at locations where MASK is not 0.
X = X(MASK, :);

% ---------------------------------------------------------------------------	
% --- Creates and returns a handle to the GUI figure. 
function bands_list_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Bands List',...
'NumberTitle','off',...
'Position',[520 357 271 443],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[5 68 260 102], 'Style','frame', 'Tag','frame_movel');
uicontrol('Parent',h1, 'Position',[5 32 260 32], 'Style','frame');
uicontrol('Parent',h1, 'Position',[5 177 260 33], 'Style','frame');

uicontrol('Parent',h1,...
'Call',@bands_list_uiCB,...
'Position',[12 185 90 15],...
'String','Gray Scale',...
'Style','radiobutton',...
'Tag','radio_gray');

uicontrol('Parent',h1,...
'Call',@bands_list_uiCB,...
'Position',[102 185 79 15],...
'String','RGB Color',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_RGB');

uicontrol('Parent',h1,...
'Call',@bands_list_uiCB,...
'FontName','Helvetica',...
'Position',[181 182 80 21],...
'Tooltip', 'Compute Principal Components', ...
'String','Compute PCA',...
'Tag','push_pca');

uicontrol('Parent',h1,...
'Call',@bands_list_uiCB,...
'Position',[10 142 35 15],...
'String','R',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_R');

uicontrol('Parent',h1,...
'Call',@bands_list_uiCB,...
'Position',[10 111 35 15],...
'String','G',...
'Style','radiobutton',...
'Tag','radio_G');

uicontrol('Parent',h1,...
'Call',@bands_list_uiCB,...
'Position',[10 82 35 15],...
'String','B',...
'Style','radiobutton',...
'Tag','radio_B');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[40 139 219 21],...
'Style','edit',...
'Tag','edit_Rband');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[40 109 219 21],...
'Style','edit',...
'Tag','edit_Gband');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[40 79 221 21],...
'Style','edit',...
'Tag','edit_Bband');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[44 38 217 21],...
'Style','edit',...
'Tag','edit_dimsDesc');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[14 40 25 15],...
'String','Dims',...
'Style','text');

uicontrol('Parent',h1,...
'Call',@bands_list_uiCB,...
'FontName','Helvetica',...
'Position',[90 5 66 21],...
'String','Load',...
'Tag','push_Load');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@bands_list_uiCB,...
'Position',[5 219 260 221],...
'Style','listbox',...
'Value',1,...
'Tag','listbox1');

uicontrol('Parent',h1,'Position',[99 149 80 15],...
'String','Selected Band','Style','text',...
'Tag','text_toGray','Visible','off');


function bands_list_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

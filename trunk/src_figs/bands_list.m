function varargout = bands_list(varargin)
% M-File changed by desGUIDE 
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bands_list (see VARARGIN) 
 
hObject = figure('Tag','figure1','Visible','off');
handles = guihandles(hObject);
guidata(hObject, handles);
bands_list_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

movegui(hObject,'west')
handles.Rband = [];     % To hold the band number that will be puted here
handles.Gband = [];     %               "
handles.Bband = [];     %               "
handles.frame_movel_pos = get(handles.frame_movel,'Pos');
handles.edit_Rband_pos = get(handles.edit_Rband,'Pos');

if (length(varargin) > 0)
    handles.h_mirone_fig = varargin{1};
    bandList = getappdata(varargin{1},'BandList');
    if (isempty(bandList))
        errordlg('ERROR: There is no list of bands available. Do you know what you are doing?','ERROR')
        delete(hObject);    return
    end
    handles.struct_names = bandList(1);
    set(handles.listbox1,'String',bandList(1))
    handles.image_bands = bandList{2};  % A MxNxP uint8 array, where P is the number of bands in memory (~= ntotal bands)
    handles.all_names = bandList{3};    % A Mx2 cell array with struct field names & names to show up in the tree
    handles.band_desc = bandList{4};    % A Mx2 cell array with the data description and band number (per struct field)
    handles.fname = bandList{5};        % File name
    handles.bands_inMemory = bandList{6};    % A vector with the band numbers already in memory
    handles.dims = bandList{7};         % A [n_row n_col nBands_total] vector with the 2D image dimensions and the
                                        % TOTAL number of bands in the dataset (some of them may not be on memory)
    handles.reader = bandList{8};       % A string with 'GDAL' if the dataset was openend with it OR:
                                        % a Mx2 cell array with {byte_resolution, header, interleave, endian};
                                        % and a 1 to 3 elements with the subsets option of multiband read.
    for ( i = 1:handles.dims(3) )
        tmp.(['band' sprintf('%d',i)]) = i;
    end
    handles.struct_values = {tmp};      % TENHO DE MUDAR ESTE NOME
else
    delete(hObject);    return
end

%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
bgcolor = get(0,'DefaultUicontrolBackgroundColor');
framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
set(0,'Units','pixels');    set(hObject,'Units','pixels')    % Pixels are easier to reason with
h_f = [handles.frame2 handles.frame3];      % The third frame (handles.frame_movel) cannot be killed
for i=1:length(h_f)
    frame_size = get(h_f(i),'Position');
    f_bgc = get(h_f(i),'BackgroundColor');
    usr_d = get(h_f(i),'UserData');
    if abs(f_bgc(1)-bgcolor(1)) > 0.01           % When the frame's background color is not the default's
        frame3D(hObject,frame_size,framecolor,f_bgc,usr_d)
    else
        frame3D(hObject,frame_size,framecolor,'',usr_d)
        delete(h_f(i))
    end
end
%------------- END Pro look (3D) -------------------------------------------------------

% Choose default command line output for bands_list_export
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes bands_list_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = bands_list_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = bands_list_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ------------------------------------------------------------------------
function radiobutton_gray_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_RGB,'Value',0)
    pos = handles.frame_movel_pos;
    pos = [pos(1) pos(2)+pos(4)/2 pos(3) pos(4)/2];
    set(handles.frame_movel,'Pos',pos)
    set(handles.edit_Rband,'Pos',pos+[10 5 -20 -28])
    set(handles.radiobutton_R,'Visible','off')
    set(handles.radiobutton_G,'Visible','off')
    set(handles.radiobutton_B,'Visible','off')
    set(handles.edit_Gband,'Visible','off')
    set(handles.edit_Bband,'Visible','off')
    set(handles.text_toGray,'Visible','on')
else
    set(hObject,'Value',1)
end

% ------------------------------------------------------------------------
function radiobutton_RGB_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_gray,'Value',0)
    set(handles.frame_movel,'Pos',handles.frame_movel_pos)
    set(handles.edit_Rband,'Pos',handles.edit_Rband_pos)
    set(handles.radiobutton_R,'Visible','on')
    set(handles.radiobutton_G,'Visible','on')
    set(handles.radiobutton_B,'Visible','on')
    set(handles.edit_Gband,'Visible','on')
    set(handles.edit_Bband,'Visible','on')
    set(handles.text_toGray,'Visible','off')
else
    set(hObject,'Value',1)
end

% ------------------------------------------------------------------------
function radiobutton_R_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_G,'Value',0)
    set(handles.radiobutton_B,'Value',0)
else
    set(hObject,'Value',1)
end

% ------------------------------------------------------------------------
function radiobutton_G_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_R,'Value',0)
    set(handles.radiobutton_B,'Value',0)
else
    set(hObject,'Value',1)
end

% ------------------------------------------------------------------------
function radiobutton_B_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_R,'Value',0)
    set(handles.radiobutton_G,'Value',0)
else
    set(hObject,'Value',1)
end

% ------------------------------------------------------------------------
function edit_Rband_Callback(hObject, eventdata, handles)

% ------------------------------------------------------------------------
function edit_Gband_Callback(hObject, eventdata, handles)

% ------------------------------------------------------------------------
function edit_Bband_Callback(hObject, eventdata, handles)

% ------------------------------------------------------------------------
function edit_dimsDesc_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------------
function pushbutton1_Callback(hObject, eventdata, handles)
% TENHO QUE TESTAR SE CASO FOR PRECISO LOADAR BANDAS ELAS SEJAM GDALICAS

if (get(handles.radiobutton_RGB,'Value') && ...
        (isempty(handles.Rband) || isempty(handles.Gband) || isempty(handles.Bband)))
    errordlg('Error: you must select three bands','ERROR')
    return
end
if ( get(handles.radiobutton_gray,'Value') && (isempty(handles.Rband)) )
    errordlg('Error: Are you blind? You must select a band. Kinda obvious no?','ERROR')
    return
end

h_img = findobj(handles.h_mirone_fig,'Type','image');
img = get(h_img,'CData');
handles_mir = guidata(handles.h_mirone_fig);        % Retrive Mirone handles
head = [];              % It will be changed only if we load a composition of non uint8

if (get(handles.radiobutton_RGB,'Value'))       % RGB - pure image for sure (is it??)
    [b,b] = ismember([handles.Rband handles.Gband handles.Bband],handles.bands_inMemory);
    idx = (b == 0);     % ISMEMBER returns zeros for elements of A not in B
    b(idx) = [];        % If they exis, clear them
    if (length(b) == 3)         % All bands are aready in memory
        img = handles.image_bands(:,:,b);
    elseif (length(b) == 2)     % Two bands are in memory. Need to load the third
        id = intersect([handles.Rband handles.Gband handles.Bband], b);
        if (isequal(id,[1 2]) | isequal(id,[2 1]))          % Red & Green in memory
            img(:,:,1:2) = handles.image_bands(:,:,b);
            img(:,:,3) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Bband)]);
        elseif (isequal(id,[1 3]) | isequal(id,[3 1]))      % Red & Blue in memory
            img(:,:,1) = handles.image_bands(:,:,b(1));
            img(:,:,2) = gdalread(handles.fname,'-S', ['-B' num2str(handles.Gband)]);
            img(:,:,3) = handles.image_bands(:,:,b(2));
        elseif (isequal(id,[2 3]) | isequal(id,[3 2]))      % Green & Blue in memory
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

    set(h_img,'CData',img)
    try rmappdata(handles.h_mirone_fig,'dem_x');    rmappdata(handles.h_mirone_fig,'dem_y');
        rmappdata(handles.h_mirone_fig,'dem_z');    rmappdata(handles.h_mirone_fig,'GMThead');
    end
    image_type = 2;         % Reset indicator that this is an image only
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
    
    if (isa(img,'uint8'))       % If we had GOTOs everything would be simpler and cleaner
        image_type = 2;         % Reset indicator that this is an image only
        computed_grid = 0;      % Reset this also
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
        head = handles_mir.head;
        head(5:6) = [double(min(min(Z))) double(max(max(Z)))];
        setappdata(handles.h_mirone_fig,'dem_z',Z);  setappdata(handles.h_mirone_fig,'dem_x',X);
        setappdata(handles.h_mirone_fig,'dem_y',Y);  setappdata(handles.h_mirone_fig,'GMThead',head);
        image_type = 1;         % Pretend this a GMT grid
        computed_grid = 1;      % But set to computed_grid to avoid attempts to reload it with grdread_m
        if (isa(Z,'uint16') || isa(Z,'int16'))
            was_int16 = 1;
        else
            was_int16 = 0;
        end
    
    end         % end if is uint8
    
    set(h_img,'CData',img)
    set(handles.h_mirone_fig,'ColorMap',gray(256))
end

if (~isempty(head)),    handles_mir.head = head;    end
handles_mir.image_type = image_type;
handles_mir.computed_grid = computed_grid;
handles_mir.was_int16 = was_int16;
guidata(handles.h_mirone_fig,handles_mir)           % Save those in Mirone handles

% --------------------------------------------------------------------------
function listbox1_Callback(hObject, eventdata, handles)
index_struct = get(hObject,'Value');
struct_names = handles.struct_names;
struct_values = handles.struct_values;

indent = '       ';
root_1 = struct_names{index_struct};
is_indent = strfind(root_1, indent);
if (isempty(is_indent)),    level = 0;
else                        level = (is_indent(end) - 1)/7 + 1;
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
function cell_array =  indent_cell(cell_array, level)

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
if (get(handles.radiobutton_RGB,'Value'))
    % put the clicked band on the box that has the active radiobutton
	r = get(handles.radiobutton_R,'Value');
	g = get(handles.radiobutton_G,'Value');
	
	if (r)
        set(handles.edit_Rband,'String',handles.all_names{idx,2})
        set(handles.radiobutton_R,'Value',0)
        set(handles.radiobutton_G,'Value',1)
        handles.Rband = n_band;
	elseif (g)
        set(handles.edit_Gband,'String',handles.all_names{idx,2})
        set(handles.radiobutton_G,'Value',0)
        set(handles.radiobutton_B,'Value',1)
        handles.Gband = n_band;
	else
        set(handles.edit_Bband,'String',handles.all_names{idx,2})
        set(handles.radiobutton_B,'Value',0)
        set(handles.radiobutton_R,'Value',1)
        handles.Bband = n_band;
	end
else        % Single Band
    set(handles.edit_Rband,'String',handles.all_names{idx,2})
    handles.Rband = n_band;
end

% --- Creates and returns a handle to the GUI figure. 
function bands_list_LayoutFcn(h1,handles);

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Bands List',...
'NumberTitle','off',...
'Position',[520 357 271 443],...
'Resize','off',...
'Tag','figure1');

h2 = uicontrol('Parent',h1,...
'Position',[5 68 260 102],...
'Style','frame',...
'Tag','frame_movel');

h3 = uicontrol('Parent',h1,...
'Position',[5 32 260 32],...
'Style','frame',...
'Tag','frame3');

h4 = uicontrol('Parent',h1,...
'Position',[5 177 260 33],...
'Style','frame',...
'Tag','frame2');

h5 = uicontrol('Parent',h1,...
'Callback',{@bands_list_uicallback,h1,'radiobutton_gray_Callback'},...
'Position',[15 185 79 15],...
'String','Gray Scale',...
'Style','radiobutton',...
'Tag','radiobutton_gray');

h6 = uicontrol('Parent',h1,...
'Callback',{@bands_list_uicallback,h1,'radiobutton_RGB_Callback'},...
'Position',[122 185 79 15],...
'String','RGB Color',...
'Style','radiobutton',...
'Value',1,...
'Tag','radiobutton_RGB');

h7 = uicontrol('Parent',h1,...
'Callback',{@bands_list_uicallback,h1,'radiobutton_R_Callback'},...
'Position',[12 142 25 15],...
'String','R',...
'Style','radiobutton',...
'Value',1,...
'Tag','radiobutton_R');

h8 = uicontrol('Parent',h1,...
'Callback',{@bands_list_uicallback,h1,'radiobutton_G_Callback'},...
'Position',[12 111 25 15],...
'String','G',...
'Style','radiobutton',...
'Tag','radiobutton_G');

h9 = uicontrol('Parent',h1,...
'Callback',{@bands_list_uicallback,h1,'radiobutton_B_Callback'},...
'Position',[12 82 25 15],...
'String','B',...
'Style','radiobutton',...
'Tag','radiobutton_B');

h10 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bands_list_uicallback,h1,'edit_Rband_Callback'},...
'HorizontalAlignment','left',...
'Position',[40 139 219 21],...
'Style','edit',...
'Tag','edit_Rband');

h11 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bands_list_uicallback,h1,'edit_Gband_Callback'},...
'HorizontalAlignment','left',...
'Position',[40 109 219 21],...
'Style','edit',...
'Tag','edit_Gband');

h12 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bands_list_uicallback,h1,'edit_Bband_Callback'},...
'HorizontalAlignment','left',...
'Position',[40 79 221 21],...
'Style','edit',...
'Tag','edit_Bband');

h13 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bands_list_uicallback,h1,'edit_dimsDesc_Callback'},...
'HorizontalAlignment','left',...
'Position',[44 38 217 21],...
'Style','edit',...
'Tag','edit_dimsDesc');

h14 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[14 40 25 15],...
'String','Dims',...
'Style','text',...
'Tag','text1');

h15 = uicontrol('Parent',h1,...
'Callback',{@bands_list_uicallback,h1,'pushbutton1_Callback'},...
'Position',[90 5 66 23],...
'String','Load',...
'Tag','pushbutton1');

h16 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@bands_list_uicallback,h1,'listbox1_Callback'},...
'Position',[5 219 260 221],...
'Style','listbox',...
'Value',1,...
'Tag','listbox1');

h17 = uicontrol('Parent',h1,'Position',[99 149 80 15],...
'String','Selected Band','Style','text',...
'Tag','text_toGray','Visible','off');


function bands_list_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

function out = select_cols(varargin)
% SELECT_COLS select column order by which the input array may be ordered
%
% SELECT_COLS() --> demo mode
% SELECT_COLS(array),           where array is MxN matrix, runs this function in the 'xy' mode
% SELECT_COLS(array,'igrf')     runs in the 'igrf' mode
% SELECT_COLS(array,'xy')       runs in the 'xy' mode
% SELECT_COLS(array,'xyz')      runs in the 'xyz' mode
% SELECT_COLS(...,...,fname)    use 'fname' in the title window
% SELECT_COLS(...,...,...,n)    do not copy more than n rows into the listboxes. If n greater than
%                               the actual number of array rows, all of them are copied. 
%                               Use Inf (NOT 'Inf') to use all rows. (Why this? SPEED people, SPEED)
% NOTE about the modes ('igrf' or 'xy' or 'xyz').
%		'xy'    Use this mode for (e.g.) 2D plots. Than X = OUT(1) and Y = OUT(2)
%		'xyz'   Use this mode if you want 3 variables. Than X = OUT(1), Y = OUT(2) and Z = OUT(3)
%				There is a checkbox in this mode that can be used to tell the reader to compute
%				hypot(X,Y). If the box is checked OUT will have a fourth element (with 0).
%		'igrf'  Use this mode for selecting the needed variables to input in the IGRF routine.
%				The bare minimum for this is lon,lat,date & altitude but these last two may be easely
%				set before calling the IGRF routine. The 'field' data is needed if computing the mag anomaly
% OUT   contains a row vector with the order by which the array columns must be arranged in
%       order to have X,Y[,Z] in the 'xy[z]' modes or lon,lat,F,Date,Altitude in the 'igrf' mode

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

demo = 0;
if (nargin == 0)        % Demo mode
    cols = 5;   rows = 1;   demo = 1;   fname = 'demo';
    text_nl = num2cell((1:1000)');
    text = [text_nl text_nl text_nl text_nl text_nl];
    str_nc = {'A','B','C','D','E'};
end

if (nargin > 0 && isnumeric(varargin{1}))
    data_in = varargin{1};
    [rows,cols] = size(data_in);
elseif (nargin == 1 && ~isnumeric(varargin{1}))       % Should we issue an error message?
    out = [];   return
end
    
if (nargin > 1 && ischar(varargin{2}) && length(varargin{2}) < 5)
    if (strcmpi(varargin{2},'igrf'))
        igrf_mode = 1;  generic_mode = 0;   xy_mode = 0;    xyz_mode = 0;
    elseif (strcmpi(varargin{2},'xy'))
        igrf_mode = 0;  generic_mode = 1;   xy_mode = 1;    xyz_mode = 0;
    elseif (strcmpi(varargin{2},'xyz'))
        igrf_mode = 0;  generic_mode = 1;   xyz_mode = 1;   xy_mode = 0;
    end
else        % Default to 'xy' mode
    igrf_mode = 0;  generic_mode = 1;   xy_mode = 1;        xyz_mode = 0;
end

if (nargin > 2 && ischar(varargin{3}))	fname = varargin{3};
else									fname = [];
end

if (nargin > 3)
    if (isinf(varargin{4}))
        use_rows = rows;
    elseif (isnumeric(varargin{4}))
        if (varargin{4} > rows)
            use_rows = rows;
        else
            use_rows = varargin{4};
        end
    end
else
    use_rows = rows;
end

screen = get(0,'ScreenSize');
fig_size = [10 100 screen(3)-200 round(screen(4)/4)];
listbox_num_w = 71;     % Listbox width for the # of lines
list_w = 91;            % Data listbox width
list_slider_w = 20;     % Listbox slider width got by trial & error (in my computer) 
list_h = fig_size(4) - 23;   % 23 = trial & error

if (~demo)      % Otherwise necessary variables were already set
	text_nl = num2cell(1:use_rows);
	text = num2cell(data_in(1:use_rows,:));
	tmp = char(65:65+cols-1);   tmp(end+1) = '#';
	str_nc = cellstr(tmp(:))';  clear tmp;
end

% -------- Create the figure
h_fig = figure('Name',['Columns Selector (' fname ')'], 'NumberTitle','off', 'Resize','off',...
    'Visible','off', 'Units','pixels', 'Position',fig_size, 'Menubar','none');

% -------- Make the listboxes
ListBox(1) = uicontrol(h_fig, 'Style', 'listbox', 'Units','pixels','BackgroundColor',[1 1 1],...
    'Position',[5 5 listbox_num_w list_h], 'String', text_nl, 'HorizontalAlignment', 'left');
pos0 = [listbox_num_w+1 5 list_w list_h];
for i=1:cols
    pos = [(pos0(1)+(i-1)*list_w) pos0(2) (list_w+list_slider_w+1) pos0(4)];
    ListBox(i+1) = uicontrol(h_fig, 'Style', 'listbox', 'Units','pixels','BackgroundColor',[1 1 1],...
        'Position',pos, 'String', text(:,i), 'HorizontalAlignment', 'left', 'Enable', 'inactive');
end
clear text;
set(ListBox(1),'Callback',{@scroll_lines,ListBox})
% -------- Make (inactive) pushbuttons with column names
for i=1:cols
    pos = [(pos0(1)+(i-1)*list_w+1) fig_size(4)-21 (list_w+list_slider_w-list_slider_w) 23];
    push_but(i) = uicontrol(h_fig, 'Style', 'pushbutton', 'Units','pixels', 'FontWeight','Bold',...
        'Position',pos, 'String', str_nc{i}, 'Tag', 'ColNames', 'Enable', 'inactive');
end
% -------- Make the labels for the popups
if (igrf_mode)
    % Here I cannot escape to do more a bunch of tests to try to guess what's inside the matrix columns
    label_pop = {'Lon','Lat','Field','Date','Altitude'};
    min_mat = min(data_in(1:min(size(data_in,1),300),:));		% Limit the min/max search. We don't want to
    max_mat = max(data_in(1:min(size(data_in,1),300),:));		% spend ages with thousands lines files.
    i_lon = find(min_mat >= -180 & max_mat <= 360);				% Search for longitude col (loosy guess)
    i_lat = find(min_mat >= -90 & max_mat <= 90);				% Search for latitude col (loosy guess)
    i_field = find(min_mat > 25000 & max_mat < 80000);			% Search for the Total field col (relatively trustful)
    i_date = find(min_mat > 1900 & max_mat < 2020);				% Search for the date col (quite trustful)
    i_alt = find(min_mat > -10 & max_mat < 1500);				% Search for altitude (km) col (loosy guess)
    if (length(i_alt) == 1 && (i_alt == 1 || i_alt == 2))		% Altitude values may be confused with lon/lat
        i_alt = [];
    elseif (length(i_alt) == 2)
		if ((i_alt(1) == 1 && i_alt(2) == 2) || (i_alt(1) == 2 && i_alt(2) == 1))   i_alt = [];
		elseif ((i_alt(1) == 1 || i_alt(1) == 2))   i_alt = i_alt(2);
		end
    elseif (length(i_alt) == 3)
        if ((i_alt(1) == 1 && i_alt(2) == 2) || (i_alt(1) == 2 && i_alt(2) == 1))   i_alt = i_alt(3); end
    end
    % Now if any of the above is empty, we'll assume that the corresponding info was not provided
    n_rec_cols = 2;     emp_field = 1;  emp_date = 1;   emp_alt = 1;
    if (isempty(i_lon))         % Cannot be
        out = [];
        errordlg('Your file doesn''t have a column with longitudes.','Error'); return
    end
    if (isempty(i_lat))         % Cannot be
        out = [];
        errordlg('Your file doesn''t have a column with latitudes.','Error'); return
    end
    if (~isempty(i_field))
        n_rec_cols = n_rec_cols + 1;        emp_field = 0;
    end
    if (~isempty(i_date))
        n_rec_cols = n_rec_cols + 1;        emp_date = 0;
    end
    if (~isempty(i_alt))
        n_rec_cols = n_rec_cols + 1;        emp_alt = 0;
    end
    pop_def_value = 1:cols+1;   % 
    if (n_rec_cols < 5)
        % Some of the possible fields are missing. We have to find which and ...
        if (cols == 2)
			pop_def_value(3) = length(str_nc);
			pop_def_value(4) = length(str_nc);
			pop_def_value(5) = length(str_nc);
        elseif (cols >= 3)
			if (emp_field)		pop_def_value(3) = length(str_nc);
			else				pop_def_value(3) = i_field(1);
			end
			if (emp_date)		pop_def_value(4) = length(str_nc);
			else				pop_def_value(4) = i_date(1);
			end
			if (emp_alt)		pop_def_value(5) = length(str_nc);
			else				pop_def_value(5) = i_alt(1);
			end
        end
    else        % (n_rec_cols == 5)
        pop_def_value(3) = i_field(1);
        pop_def_value(4) = i_date(1);
        pop_def_value(5) = i_alt(1);
    end
    n_labels = 5;
else            % XY[Z] mode 
    if (cols == 1)
        out = [];
        errordlg('Idiot selection. Your file has only one column.','Chico Clever'); return
    elseif (cols == 2)
        label_pop = {'X','Y'};          n_labels = 2;
    elseif (cols >= 3 & xy_mode)
        label_pop = {'X','Y'};          n_labels = 2;
    elseif (cols >= 3 & xyz_mode)
        label_pop = {'X','Y','Z'};      n_labels = 3;
    end
end
for (i=1:n_labels)     % -------- Well, finally do the labels
    pos_label = [(pos(1)+list_w+21), (fig_size(4)-19-i*30), 41, 23];     % pos(1) = x of last listbox
    label_ui(i) = uicontrol(h_fig, 'Style', 'text', 'Units','pixels',...
        'BackgroundColor',get(h_fig,'Color'),...
        'Position',pos_label, 'String', label_pop{i}, 'HorizontalAlignment', 'right');
end
% -------- Make the popups (for column order selecting)
% We need also to guess the default Value of each popup
if (~igrf_mode)     % Otherwise it has already been deffined above
    pop_def_value = 1:n_labels;
    str_nc(end) = [];       % We don't want the # symb that was at the end
    ind = 1:length(str_nc);
else        % Take into account that cols may be inferior to 5
    ind = 1:length(str_nc)-1;
    ind(length(str_nc):n_labels) = length(str_nc);
end
for (i=1:n_labels)      % Number of labels has to be the same as number of popups
    pos_pop = [(pos_label(1)+pos_label(3)+7), (fig_size(4)-15-i*30), 71, 23];
    popup_but(i) = uicontrol(h_fig, 'Style', 'popupmenu', 'Units','pixels','BackgroundColor','white',...
        'Position',pos_pop, 'String', str_nc, 'Value',pop_def_value(ind(i)), 'HorizontalAlignment','left');
end
if (igrf_mode & cols == 2)     % Only lon and lat cols may be interchanged
    set(popup_but(3:5),'Enable','inactive');
end
% -------- Make one/two checkboxes for selecting what to write in file (only for the igrf or xyz cases)
chk_but = [];
if (igrf_mode)
    tip1 = 'Add a column to the file containg the computed IGRF Total field';
    tip2 = 'Add a column to the file containg the computed anomaly (measured field - IGRF)';
    chk_but(1) = uicontrol(h_fig, 'Style','checkbox', 'Units','pixels', 'Value',1,...
        'Position',[pos_pop(1) 70 101 18], 'String','Write total field', 'TooltipString',tip1);
    chk_but(2) = uicontrol(h_fig, 'Style','checkbox', 'Units','pixels', 'Value',1,...
        'Position',[pos_pop(1) 45 101 18], 'String','Write anomaly', 'TooltipString',tip2);
    if (cols == 2)
        set(chk_but(2),'Value',0,'Enable','inactive')   % No Total Field no anomaly
    end
else
	if (xyz_mode)
	    chk_but = uicontrol(h_fig, 'Style','checkbox', 'Units','pixels', 'Value',0,...
			'Position',[pos_pop(1) 70 101 18], 'String','A & B to dist', 'Tooltip','Tell reader to compute hypot(A,B)');
	end
end
% -------- Make the OK button
OK_but = uicontrol(h_fig, 'Style', 'pushbutton', 'Units','pixels',...
		'Position',[pos_pop(1) 8 91 25],  'String', 'OK', 'FontWeight','Bold',...
		'Callback', {@OK_push, popup_but, ListBox, chk_but, igrf_mode, xyz_mode},'Tag','ok');

% Cut the eventualy remaining figure width to an apropriate size.
% NOTE: A slider should be introduced if the final width is wider than the screen
pos = get(OK_but,'Position');

set(h_fig,'Position',[1 1 (pos(1)+pos(3)+10) fig_size(4)])

if ( (pos(1) + pos(3)) > screen(3) )
	first_listBox_pos = get(ListBox(1),'Position');
	slider_pos = [first_listBox_pos(1) first_listBox_pos(2)-5 (first_listBox_pos(1) + pos(1)+pos(3)+5) 10];
	h_all_uis = findobj(gcf,'Style','listbox');     % Get the Listbox handles
	pos_all_uis = get(h_all_uis,'Position');
	
	h_tmp = findobj(gcf,'Style','pushbutton');      % Get the pushbutton column names handles
	h_all_uis = [h_all_uis; h_tmp];
	pos_all_uis = [pos_all_uis; get(h_tmp,'Position')];
	
	h_tmp = findobj(gcf,'Style','popupmenu');       % Get the popupmenu handles
	h_all_uis = [h_all_uis; h_tmp];
	pos_all_uis = [pos_all_uis; get(h_tmp,'Position')];
	
	h_tmp = findobj(gcf,'Style','checkbox');        % Get the checkbox handles
	h_all_uis = [h_all_uis; h_tmp];
	pos = get(h_tmp,'Position');
	if (~isa(pos,'cell')),		pos_all_uis = [pos_all_uis; {pos}];
	else						pos_all_uis = [pos_all_uis; pos];
	end
	clear h_tmp pos;
	cb_slider = {@move_all_uis,h_all_uis,pos_all_uis};
	uicontrol('style','slider','units','pixels','position',slider_pos,...
		'callback',cb_slider,'min',0,'max',slider_pos(3));
end

move2side(h_fig,'left');
set(h_fig,'Visible','on');
uiwait;
but = gco;
if strcmp(get(but,'tag'),'ok')
    out = get(OK_but,'UserData');       % Output was stored there
    delete(h_fig);
else                                    % Figure was killed
    out = [];
end

% -----------------------------------------------------------------------------------------
function move_all_uis(obj,eventdata,h,orig_pos)
	m = length(h);
	for i=1:m
		new_pos = [orig_pos{i}(1)-get(gcbo,'value') orig_pos{i}(2) orig_pos{i}(3) orig_pos{i}(4)];
		set(h(i),'Position',new_pos);
	end

% -----------------------------------------------------------------------------------------
function OK_push(obj, eventdata, h_pop, h_list, h_chk, igrf_mode, xyz_mode)
	n_pops = numel(h_pop);    % How many popups?
	col = zeros(1, n_pops);
	for (i = 1:n_pops)            % See what was choosen by each popup
		col(i) = get(h_pop(i),'Value');
	end
	
	if (igrf_mode)
		contents = get(h_pop(3),'String');
		if ( strcmp(contents{get(h_pop(3),'Value')},'#') ),		col(3) = 0;     end
		contents = get(h_pop(4),'String');
		if ( strcmp(contents{get(h_pop(4),'Value')},'#') ),		col(4) = 0;     end
		contents = get(h_pop(5),'String');
		if ( strcmp(contents{get(h_pop(5),'Value')},'#') ),		col(5) = 0;     end
		% See if Total field and Anomaly are still selected for writing in file
		write_total_field = 0;			write_anom = 0;
		if (get(h_chk(1),'Value')),		write_total_field = 1;   end
		if (get(h_chk(2),'Value')),		write_anom = 1;   end
		set(obj,'UserData',[col write_total_field write_anom])
	else
		if (xyz_mode && get(h_chk(1),'val'))
			col = [col 0];				% This extra value is to be interpreted by the reader
		end								% to compute hypot(col(1),col(2))
		set(obj,'UserData',col)
	end
	uiresume

% -----------------------------------------------------------------------------------------
function scroll_lines(obj,eventdata,h_list)
% Scrolls the data fields according to value of the # lines listbox
	set(h_list(2:end),'Value',get(gcbo,'value'))

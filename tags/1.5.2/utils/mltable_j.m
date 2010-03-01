function data = mltable_j(fig, hObj, action, columnInfo, rowHeight, cell_data, gFont)
% function data = mltable_j(fig, hObj, action, columnInfo, rowHeight, cell_data, gFont)
%
% Author: Morris Maynard
% Based on code by Gregory Gershanok
%
% Manages a table with editing and scrolling capability
% Features: varying column widths, entry formatting, insert/delete rows,
%           editable or read-only cells, font control, scaled numeric
%           display, multiple tables per figure,
%
% Usage:
% Supply the parent figure, the handle of the axes object to use, the
% 'CreateTable' action, info about columns, and the cell data 
%
% 
% columninfo.format = '%4.6g';
% columninfo.weight =      [ 1, 1, 1, 1, 1];
% rowHeight = 16;
% gFont.size=9;
% gFont.name='Helvetica';
% 
% mltable_j(fig, tbl, 'CreateTable', columninfo, rowHeight, cell_data, gFont);
%
% To use in a GUIDE-created figure:
%
% Create a figure including a blank "axes" object with the tag 'tblParams'
% Put the lines starting with the "cell_data" line above into your figure's
% OpeningFcn, but replace the mltable_j line with:
%
% mltable_j(gcf, handles.tblParams, 'CreateTable', columninfo, rowHeight, cell_data, gFont);
% so clicking outside the table will finish edit in progress...
% endfcn = sprintf('mltable_j(%14.13f, ''SetCellValue'');', hObject, handles.tblParams);
% set(hObject,'buttondownfcn',endfcn);
%
% To access the data edited by the table:
%
% info = get(tbl, 'userdata');
% data = info.data;
%

%-------------------------------------------------------------------------
% NOTE: This is a modified version of MLTABLE to work with the DIGITALFILTERING_GUI
%       Many original features were removed (not needed here) and some vectorization.
%
%   Joaquim Luis    24-Jan-2006
%

global MINROWS;
MINROWS = 3;

if (nargin == 3),   columnInfo = 1;    end  % 'columnInfo' works as here as the number of row/col to add/remove
    
switch(action)
  case 'CreateTable'
    data = createTable(fig, hObj, columnInfo, rowHeight, cell_data, gFont);
  case 'ResizeTable'
    resizeTable(fig, hObj);
  case 'ScrollData'
    scrollData(hObj);
  case 'EditCell'
    editCell(fig, hObj);
  case 'SetCellValue'
    setCellValue(hObj);
  case 'AddRow'
    addRow(fig, hObj,columnInfo);       % 'columnInfo' is here the number of rows to add
  case 'DelRow'
    delRow(fig, hObj,columnInfo);       % 'columnInfo' is here the number of rows to remove
  case 'AddCol'
    addCol(fig, hObj,columnInfo);       % 'columnInfo' is here the number of columns to add
  case 'DelCol'
    delCol(fig, hObj,columnInfo);       % 'columnInfo' is here the number of columns to remove
end

%-------------------------------------------------------------------------
function data = createTable(fig, hObj, columnInfo, rowHeight, cell_data, gFont)
% Initially creates the table
%-------------------------------------------------------------------------
global MINROWS
data.figure = fig;
% get axes position in pixel coordinates
set(hObj, 'units', 'pixels');
set(hObj, 'visible', 'on');
pos_ax = get(hObj, 'position');
% set up grid info structure
ds = size(cell_data);
data.maxRows = ds(1);
data.maxCols = ds(2);
if (data.maxRows < MINROWS)
    blanks = cell(1, ds(2));
    for ii = data.maxRows+1:MINROWS
        cell_data = [cell_data; blanks];
    end
    data.maxRows = MINROWS;
end
data.data = cell_data;
data.axes = hObj;
data.userModified = zeros(ds);
data.rowHeight = rowHeight;
data.columnInfo = columnInfo;
data.n_RowsSlider = ds(1);
data.numRows= ds(1);
data.numCols= ds(2);
data.ltGray = [92 92 92]/255;
data.OffscreenPos = [-1000 -1000 30 20];
data.selectedRow = 0;
data.selectedCol = 0;
data.gFont = gFont;
data.no_edit = 0;       % When == 1 cells are not allowed to edit

% use 0...1 scaling on table x and y positions
set(fig, 'CurrentAxes', data.axes);
set(data.axes, 'box', 'on', 'DrawMode', 'fast');
set(data.axes, 'xlimmode', 'manual', 'xlim', [0 1], 'ylim', [0 1], ...
               'xtick', [], 'ytick', [], 'xticklabelmode', 'manual', 'xticklabel', []);
           
pos_ax(3) = pos_ax(3) - 10; % width of slider
set(data.axes, 'position', pos_ax, 'LineWidth', 2);
% callback for starting editing 
editfcn = sprintf('mltable_j(%14.13f, %14.13f, ''EditCell'');',fig, hObj);
set(data.axes, 'ButtonDownFcn', editfcn);
% callback for scrolling table
scrfcn = sprintf('mltable_j(%14.13f, %14.13f, ''ScrollData'');',fig, hObj);
data.slider = uicontrol('style', 'slider', 'units', 'pixels',...
    'position', [pos_ax(1)+pos_ax(3)+2 pos_ax(2) 16 pos_ax(4)],...
    'Callback', scrfcn);

set(hObj, 'UserData', data);
% so clicking outside the table will finish edit in progress
endfcn = sprintf('mltable_j(%14.13f, %14.13f, ''SetCellValue'');', fig, hObj);
set(fig,'buttondownfcn',endfcn);

resizeTable(fig, hObj);

%-------------------------------------------------------------------------
function data = resizeTable(fig, hObj)
% fit table within boundaries and update scrollbar
% at this time doesn't handle figure resize
%-------------------------------------------------------------------------
data = get(hObj, 'UserData');
if (isempty(data)),   return;   end

cla(hObj);

set(hObj, 'units', 'pixels');
set(fig,'CurrentAxes',hObj);  

pos_ax = get(hObj,'position');
data.n_RowsSlider = floor((pos_ax(4)-(2*data.rowHeight))/data.rowHeight);
if (data.n_RowsSlider > data.maxRows)
    data.n_RowsSlider = data.maxRows;
end

% See if we need a vertical slider
if(data.n_RowsSlider < data.maxRows)
    set(data.slider,'Units','pixels');
	set(data.slider, 'visible', 'on',...
                   'position', [pos_ax(1)+pos_ax(3)+1 pos_ax(2) 16 pos_ax(4)], ...
                   'min', 0, 'max', data.maxRows - data.n_RowsSlider, 'value', data.maxRows - data.n_RowsSlider, ...
                   'sliderstep', [1/(data.maxRows-data.n_RowsSlider) data.n_RowsSlider/(data.maxRows-data.n_RowsSlider)]);
else  
  set(data.slider, 'visible', 'off');
end

% get gui units for rows and columns
% average column width and row height
d_x = 1 / data.numCols;
%d_y = 1 / data.numRows;
% minimum adjust unit
unit_d_h = 1 / pos_ax(3);
%unit_d_v = 1 / (pos_ax(4) - data.rowHeight);

% Horizontal line positions
lx_h = ones(2, data.n_RowsSlider);
lx_h(1, :) = 0;
ly_h = [1:data.n_RowsSlider; 1:data.n_RowsSlider] / data.n_RowsSlider;

% Vertical line positions
ly_v = ones(2, data.numCols);
ly_v(1, :) = 0;
lx_v = [d_x*data.numCols 2:data.numCols; d_x*data.numCols 2:data.numCols] / data.numCols;

% draw initial grid
data.vertLines  = line(lx_v, ly_v);
data.vertLines1  = line(lx_v(:, 1:data.numCols-1), ly_v(:, 1:data.numCols-1));
data.horizLines = line(lx_h, ly_h);
set(data.horizLines, 'color', data.ltGray);
set(data.vertLines, 'color', data.ltGray, 'LineWidth', 2);
set(data.vertLines1, 'color', [1 1 1], 'LineWidth', 0.5);

% now display text in grid     
txt_x = [0:data.numCols-1] / data.numCols + 4*unit_d_h;
%txt_x = [0:data.numCols-1]/data.numCols + d_x/2.3;

data.txt_x = txt_x;
data.txtCells = zeros(data.numRows, data.numCols);
uictx = get(hObj,'UIContextMenu');

n_right = min(data.numRows,data.n_RowsSlider);
for j = 1:data.numRows
    %txt_y = (data.n_RowsSlider-j) / data.n_RowsSlider * ones(1, data.numCols);
    txt_y = (n_right-j) / n_right * ones(1, data.numCols);
	data.txtCells(j, :) = text(txt_x, txt_y, 'a','Clipping','on');
	for i = 1:data.numCols
        nums = data.data{j, i};
        nums = num2str(nums, data.columnInfo.format);
        set(data.txtCells(j, i), 'string', nums, 'Position', [txt_x(i), txt_y(i)], 'UIContextMenu', uictx);
	end
end
set(data.txtCells(:, :), 'FontSize', data.gFont.size, 'FontName', data.gFont.name, 'FontWeight', 'normal', ...
                         'HorizontalAlignment','left','VerticalAlignment', 'bottom');
set(data.txtCells(1:data.numRows, 1:data.numCols), 'buttondownfcn', get(data.axes, 'ButtonDownFcn'));

set(hObj, 'UserData', data);
scrollData(hObj);
set(hObj, 'units', 'normalized');

%-------------------------------------------------------------------------
function data = scrollData(hObj)
% handle scrollbar (slider) callback
%-------------------------------------------------------------------------
data = get(hObj, 'UserData');
if (isempty(data)),   return;   end

set(data.slider,'Units','pixels');
if isfield(data,'editBox') && ishandle(data.editBox)
    delete(data.editBox);
end

% handle non-scroll case in case slider was switched off
if(strcmp(get(data.slider, 'visible'), 'off') == 1)
	ind0 = 0;
	for i = 1:data.numRows
		for j = 1:data.numCols
            nums = data.data{i, j};
            nums = num2str(nums, data.columnInfo.format);
            set(data.txtCells(i, j), 'string', nums);
		end
	end
else	
	val = get(data.slider, 'Value');
	val0 = data.maxRows - data.n_RowsSlider;
	ind0 = round(val0-val);
	% move the text to give illusion of scrolling
    for i = ind0+1:(ind0+data.n_RowsSlider)
        for j = 1:data.numCols
            nums = data.data{i, j};
            nums = num2str(nums, data.columnInfo.format);
            set(data.txtCells(i-ind0, j), 'string', nums);
        end
    end
end

data.ind0 = ind0;       % save scroll position
data.hpatch = remakepatch(data.selectedRow, data);
set(data.slider,'Units','normalized');
set(hObj, 'UserData', data);

% -----------------------------------------------------------------------
function hpatch = remakepatch(row, data)
% -----------------------------------------------------------------------
hpatch = []; % empty return if nothing to do
% see if selected row is visible
if ((row - data.ind0) <= (data.numRows) && (row - data.ind0) > 0)
    hpatch = data.hpatch;       % return previous patch
else                            % if patch is no longer visible, delete it
    if isfield(data,'hpatch') && ~isempty(data.hpatch) && ishandle(data.hpatch)
        delete(data.hpatch);
    end
end

%-------------------------------------------------------------------------
function data = editCell(fig,hObj)
% put an edit control over the selected cell
%-------------------------------------------------------------------------

data = get(hObj, 'UserData');
if (isempty(data)), return;     end
if (data.no_edit),  return;     end

set(hObj, 'units', 'pixels');
pt = get(hObj, 'CurrentPoint');
pt = pt(1,:); % strip out 2nd axis info
pos_ax = get(hObj, 'position');
pt(1) = pos_ax(1) + (pt(1) * pos_ax(3));
pt(2) = pt(2) .* pos_ax(4);

d_x = pos_ax(3) / data.numCols;
d_y = pos_ax(4) / data.numRows;
  
% find column index
col = -1;   p1 = 0;
for i = 1:data.numCols
    p2 = p1 + d_x;
    if ((p1 <= (pt(1)-pos_ax(1))) && (p2 >= (pt(1)-pos_ax(1))))
        col = i;
        break;
    else 
        p1 = p2;
    end
end
if (col == -1)
    set(hObj, 'units', 'normalized');
    return
end  

% find row index
row = data.numRows - (floor(pt(2) / d_y));

if (row < 1) % could be header row
    set(hObj, 'units', 'normalized');
    return
end  

data.selectedCol = col;
data.selectedRow = row + data.ind0;

if (isfield(data, 'editBox') && ishandle(data.editBox))
    delete(data.editBox);
end
    
data.hpatch = remakepatch(data.selectedRow, data);

if (row <= data.numRows)
	unit_d_h = 1/pos_ax(3);
    %unit_d_v = 1/pos_ax(4);
    ebtxt = data.data{row + data.ind0, col};
    ebtxt = num2str(ebtxt, data.columnInfo.format);
	% set the edt control contents and position
	% callback for entering cell data
    endfcn = sprintf('mltable_j(%14.13f, %14.13f, ''SetCellValue'');', fig, hObj);
	data.editBox = uicontrol('style', 'edit', 'units', 'pixels', 'Callback', endfcn);
	set(data.editBox, 'FontSize', data.gFont.size, 'FontName', data.gFont.name, 'FontWeight', 'normal');
	set(data.editBox, 'string', ebtxt, 'UserData', [row col]);
	ext_eb = get(data.editBox, 'extent');
	ext_eb(4) = d_y + 3;
	pos = [(pos_ax(1)-unit_d_h+p1) ,...
            pos_ax(2) + ...                             % start of table 
              ceil((data.numRows - row) * d_y) + ...    % cvt index to row #, get offset from ystart
                 (ceil(d_y) - ext_eb(4))/2, ...          % add half of the ctrl height??
            floor(d_x)-unit_d_h, ext_eb(4)];
	set(data.editBox, 'Position', pos, 'HorizontalAlignment' ,'Center');
	set(fig, 'CurrentObject', data.editBox);
end
set(hObj, 'UserData', data);
set(hObj, 'units', 'normalized');

%-------------------------------------------------------------------------
function data = setCellValue(hObj)
% when edit control calls back, update data in cell
%-------------------------------------------------------------------------
data = get(hObj, 'UserData');
if (isempty(data)),   return;   end

if (~isfield(data,'editBox') || ~ishandle(data.editBox)),   return;     end

ind = get(data.editBox, 'UserData');
if (isempty(ind)),   return;    end

nums = get(data.editBox, 'string');
row = ind(1) + data.ind0;       col = ind(2);
d_old = data.data{row, col};
num = sscanf(nums, '%f');
if(isempty(num))
    errordlg('Please enter a valid number', 'Error', 'modal');
    return;
end     
if (d_old == num),  delete(data.editBox);   return;     end
nums = num2str(num, data.columnInfo.format);
data.data{row, col} = num;

if (ishandle(data.editBox)),    delete(data.editBox);   end
if ishandle(data.txtCells(row - data.ind0, col))
    set(data.txtCells(row - data.ind0, col), 'string', nums);
    set(hObj,'UserData',data);
end

%-------------------------------------------------------------------------
function data = addRow(fig, hObj, n_rows)
% insert N_ROWS (empty) into the table
%-------------------------------------------------------------------------
data = get(hObj, 'UserData');
if (isempty(data)),   return;   end

% increase the number of rows
data.maxRows = data.maxRows + n_rows;
data.numRows = data.maxRows;            % <-----------------
data.data = [data.data; repmat({0},n_rows,size(data.data,2))];
set(hObj, 'UserData', data);
resizeTable(fig, hObj);

%-------------------------------------------------------------------------
function data = delRow(fig, hObj, n_remove)
% delete N_REMOVE rows from table
%-------------------------------------------------------------------------
global MINROWS

data = get(hObj, 'UserData');
if (isempty(data)),     return;   end
if (data.maxRows < 2),  return;   end

% decrease the number of rows
data.maxRows = data.maxRows - n_remove;
data.numRows = data.maxRows;            % <-----------------
data.data = data.data(1:size(data.data,1)-n_remove,:);
set(hObj, 'UserData', data);
resizeTable(fig, hObj);

%-------------------------------------------------------------------------
function data = addCol(fig, hObj, n_cols)
% insert N_COLS into the table
%-------------------------------------------------------------------------
data = get(hObj, 'UserData');
if (isempty(data)),   return;   end

data.maxCols = data.maxCols + n_cols;        % increase the number of cols
data.numCols = data.maxCols;
data.data = [data.data repmat({''},size(data.data,1),n_cols)];
set(hObj, 'UserData', data);
resizeTable(fig, hObj);

%-------------------------------------------------------------------------
function data = delCol(fig, hObj, n_remove)
% delete N_REMOVE columns from table
%-------------------------------------------------------------------------

data = get(hObj, 'UserData');
if (isempty(data)),     return;   end
if (data.maxCols < 2),  return;   end

data.maxCols = data.maxCols - n_remove;             % decrease the number of cols
data.numCols = data.maxCols;
data.data = data.data(:,1:size(data.data,2)-n_remove);
set(hObj, 'UserData', data);
resizeTable(fig, hObj);

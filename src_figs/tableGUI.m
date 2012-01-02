function out = tableGUI(varargin)
%   TABLEGUI - Spreadsheet like display and edition of a generic 2D array. By generic it is
%   mean that the array can be a numeric MxN matrix or a MxN cell array (with mixed
%   number and text strings). This function imitates the table cells with edit boxes
%   which may become slow if the number of elements required is large. However, it
%   works fast for small matrices. If the default number of rows is exceeded a vertical
%   slider is created. The slider works by changing the position of the table elements,
%   which again may become slow if the default number of visible rows is large. Otherwise
%   it works pretty fast.
%
%   USAGE:
%       OUT = TABLEGUI(varargin)
%
%   Inputs are in property/value pairs. All properties are strings but the values are
%   of different types depending on the case.
%
%   PROPERTY                        VALUE                                       TYPE
%
%   'array'         It can be either an numeric matrix or a MxN                 numeric or cell array
%                   cell array. 
%   'NumRows'       Total number of rows to create. Use this when the           integer
%                   'array' option is not used, or when you want to
%                   create extra empty rows.
%   'NumCol'        Total number of columns to create. Use this when            integer
%                   the 'array' option is not used, but ignored if it was.
%   'MAX_ROWS'      Number of visible rows. If NumRows > MAX_ROWS               integer
%                   a vertical slider is created (DEF = 10)
%   'RowHeight'     editboxes height in pixels (DEF = 20)                       integer
%   'ColWidth'      editboxes width in pixels (DEF = 60)                        scalar or row vector
%                   If a vector is transmited, it must contains NumCol
%                   elements which set individual column widths.
%   'bd_size'       space between editboxes width in pixels (DEF = 0)           integer
%   'HorAlin'       editboxes 'HorizontalAlignment' property. It can            string
%                   be either 'left', 'center' or 'right' (DEF = 'center')
%   'HdrButtons'    create a first row of buttons to hold column                string - either '' or 'y'
%                   names. Give an empty string ('') if you don't want
%                   column names (DEF = 'Y').
%   'ColNames'      Cell array of strings for column names. If not              1xN cell array of strings
%                   provided the columns are named 'A', 'B', ...
%   'RowNumbers'    Add a first column with row numbers. Note that this         string - either '' or 'y'
%                   column is set to 'inactive' and its not transmited
%                   on the output (DEF = '').
%   'RowNames'      Add a column with row names. Note that this column          1xN cell array of strings
%                   is set to 'inactive' and its not transmited on the
%                   output (DEF = ''). Warning: do not abuse on the 
%                   Names length. Added by Martin Furlan
%   'RowNamesWidth' Width of buttons containing the RowNames
%                   in pixels (DEF = 60)                                        scalar
%   'checks'        If = 'Y' it creates a vertical line of checkboxes           string - either '' or 'y'
%                   This affects what is send as output. Only rows that
%                   have it's checkbox checked will be returned.
%   'FigName'       Name that appears in the title bar (DEF = 'Table').         string
%   'position'      Screen location to be used in the call to MOVEGUI           string
%                   See doc of that function for valid position
%                   strings (DEF = 'east').
%   'modal'         By default the window works in MODAL mode. Give an          string - either '' or 'y'
%                   empty string ('') if you don't want it to be MODAL.
%                   In this later case the output OUT, if requested, will
%                   contain the figure handle but see more about this below.
%
%   OUT - the output - contains different things depending whether or not the
%       figure works in MODAL mode. If yes, OUT is a MxN cell array with the
%       elements retrived from the contents of the edit boxes. Otherwise, OUT
%       will contain the figure's handle. This handle has the 'UserData' property
%       filled with a structure (called hand) which contains the handles of all
%       uicontrols. Use this option if you want to interact with the TABLEGUI
%       figure inside your own code.
%
%   Examples:
%     - Display a 12x6 numeric matrix with two extra blank rows appended to the end
%       out = tableGUI('array',rand(12,6),'ColNames',{'1' '2' '3' '4' '5' '6'},'NumRows',14);
%
%     - Create a cell array with the first column filled with the row number and use
%       columns with different widths. Create also check boxes.
%       zz=cell(4,5);   zz(:,1) = num2cell(1:4)';
%       out = tableGUI('array',zz,'ColNames',{'N','A','B','C','D'},'ColWidth',[20 60 60 60 60],'checks','y');
%
%     - Create a similar table as in the previous example but with the row numbers option.
%       out = tableGUI('array',cell(4,5),'RowNumbers','y','checks','y');
%
%     - Display the Control Points of the Image Processing Toolbox example
%       "Registering an Aerial Photo to an Orthophoto"
%       load westconcordpoints  % load some points that were already picked
%       gcp = [base_points input_points];
%       out = tableGUI('array',gcp,'RowNumbers','y','ColNames',{'Base Points - X','Base Points - Y',...
%               'Input Points - X','Input Points - Y'},'ColWidth',110,'FigName','GCP Table');
%
%     - Create an empty 12x6 empty table
%       out = tableGUI;
%
%   Acknowledgment
%       This function uses the parse_pv_pairs of John d'Errico
%
%   AUTHOR
%       Joaquim Luis (jluis@ualg.pt)   17-Feb-2006
%
%   Revision
%       17-Sep-2006 - There was an error when only one check box existed.
%
%   Addition
%       6-Dec-2006 - A column with row names added by Martin Furlan
%       (martin.furlan@iskra-ae.com)
%
%       13-Mar-2007 - Added option to control the RowNames width.

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

hand.NumRows = 12;          hand.NumCol = 6;
hand.MAX_ROWS = 10;         hand.left_marg = 10;
hand.RowHeight = 20;        hand.bd_size = 0;
hand.checks = '';           hand.HdrButtons = 'Y';
hand.HorAlin = 'center';    hand.FigName = 'Table';
hand.array = cell(hand.NumRows,hand.NumCol);
hand.modal = 'y';           hand.position = 'east';
hand.RowNumbers = '';       d_col = 0;
hand.RowNames = ''; 

if (nargin == 0)        % Demo
    hand.ColNames = {'A' 'B' 'C' 'D' 'E' 'F'};
    hand.ColWidth = [50 80 80 80 80 50];
else
    hand.ColNames = '';	    hand.ColWidth = [];
    hand.array = [];        def_NumRows = hand.NumRows;
    hand.RowNamesWidth = 60;
    hand = parse_pv_pairs(hand,varargin);
    if (~isempty(hand.array))
        if (numel(hand.array) == 1 && numel(hand.array{1}) > 1)
            error('The "array" argument must be a MxN cell array and not a {MxN} cell')
        end
        [NumRows,hand.NumCol] = size(hand.array);
        if (~iscell(hand.array))    % We need as a cell array to be more general
            hand.array = num2cell(hand.array);
        end
        if (NumRows < hand.NumRows && hand.NumRows ~= def_NumRows)     % Extra rows requested
            hand.array = [hand.array; cell(hand.NumRows-NumRows,hand.NumCol)];
        else
            hand.NumRows = NumRows;
        end
        if (hand.NumRows < hand.MAX_ROWS),    hand.MAX_ROWS = hand.NumRows;     end
        
    else                % 'array' not transmited
        hand.array = cell(hand.NumRows,hand.NumCol);
    end
    
    if (isempty(hand.ColNames) && ~isempty(hand.HdrButtons))      % By default columns are labeled 'A','B',...
        hand.ColNames = cell(1,hand.NumCol);
        for (i = 1:hand.NumCol),     hand.ColNames{1,i} = char(i+64);  end
    end
    if (size(hand.array,2) > size(hand.ColNames,2))
        error('"ColNames" argument has less elements than the number of columns is "array"')
    end
    if (isempty(hand.ColWidth))                    % Use default value for button width
        hand.ColWidth = repmat(60,1,hand.NumCol);
    elseif (numel(hand.ColWidth) == 1)             % 'ColWidth' was a scalar
        hand.ColWidth = repmat(hand.ColWidth,1,hand.NumCol);
    end
    
    if (~isempty(hand.RowNumbers))                 % Row numbering was requested
        hand.ColWidth = [35 hand.ColWidth];
        hand.NumCol = hand.NumCol + 1;
        hand.array = [cell(hand.NumRows,1) hand.array];
        d_col = 1;
    end
end

arr_pos_xi = hand.left_marg + [0 (cumsum(hand.ColWidth+hand.bd_size))];
arr_pos_xi(end) = [];      % We don't want the last element
arr_pos_xw = hand.ColWidth;

% ---------------- Create the figure ----------------------------------
fig_height = min(hand.NumRows,hand.MAX_ROWS) * (hand.RowHeight+hand.bd_size);
if (~isempty(hand.HdrButtons)),    fig_height = fig_height + 22;   end     % Make room for header buttons
if (~isempty(hand.modal)),          fig_height = fig_height + 30;   end     % Make room for OK,Cancel buttons
pos = [5 75 sum(arr_pos_xw)+hand.left_marg+(hand.NumCol-1)*hand.bd_size+15 fig_height];  % The 15 is for the slider
nW = hand.RowNamesWidth;        % Short name for the row names width buttons
if (~isempty(hand.checks)),     pos(3) = pos(3) + 15;       end         % Account for checkboxes size
if (~isempty(hand.RowNames)),   pos(3) = pos(3) + nW;       end         % Account for row names size

hand.hFig = figure('unit','pixels','NumberTitle','off','Menubar','none','resize','on','position', ...
    pos,'Name',hand.FigName,'Resize','off','Visible','off');
move2side(hand.hFig,hand.position)

hand.arr_pos_y = (fig_height-hand.RowHeight-hand.bd_size - (0:hand.NumRows-1)*(hand.RowHeight+hand.bd_size))';
if (~isempty(hand.HdrButtons)),     hand.arr_pos_y = hand.arr_pos_y - 22;   end

if (~isempty(hand.checks) && isempty(hand.RowNames))    % Create the checkboxes uicontrols
    arr_pos_xi = arr_pos_xi + 15;                       % Make room for them
    hand.hChecks = zeros(hand.NumRows,1);
    hand.Checks_pos_orig = [ones(hand.NumRows,1)*7 (hand.arr_pos_y+3) ones(hand.NumRows,1)*15 ones(hand.NumRows,1)*15];
end
if (~isempty(hand.RowNames) && isempty(hand.checks))    % Create the row names
    arr_pos_xi = arr_pos_xi + nW;                       % Make room for them
    hand.RowNames_pos_orig = [ones(hand.NumRows,1)*7 (hand.arr_pos_y) ones(hand.NumRows,1)*nW ones(hand.NumRows,1)*20];
end
if (~isempty(hand.checks) && ~isempty(hand.RowNames))   % Create both the checkboxes uicontrols and the row names
    arr_pos_xi = arr_pos_xi + 15 + nW + 2;              % Make room for them
    hand.hChecks = zeros(hand.NumRows,1);
    hand.Checks_pos_orig = [ones(hand.NumRows,1)*7 (hand.arr_pos_y+3) ones(hand.NumRows,1)*15 ones(hand.NumRows,1)*15];
    hand.RowNames_pos_orig = [ones(hand.NumRows,1)*7 + 15 + 5 (hand.arr_pos_y) ones(hand.NumRows,1)*nW ones(hand.NumRows,1)*20];
end

hand.hEdits = zeros(hand.NumRows,hand.NumCol);
hand.Edits_pos_orig = cell(hand.NumRows,hand.NumCol);

% ---------------- Create the edit uicontrols ---------------------------
for (i = 1:hand.NumRows)
    if (~isempty(hand.checks))
        hand.hChecks(i) = uicontrol('Style','checkbox','unit','pixels','position', ...
            hand.Checks_pos_orig(i,:),'Value',1);
    end
    for (j = 1:hand.NumCol)
        hand.Edits_pos_orig{i,j} = [arr_pos_xi(j) hand.arr_pos_y(i) arr_pos_xw(j) hand.RowHeight];
        hand.hEdits(i,j) = uicontrol('Style','edit','unit','pixels','backgroundcolor','w','position', ...
            [arr_pos_xi(j) hand.arr_pos_y(i) arr_pos_xw(j) hand.RowHeight],'String',hand.array{i,j},...
            'HorizontalAlignment',hand.HorAlin);
    end
    if (~isempty(hand.RowNumbers))
        set(hand.hEdits(i,1),'String',i,'Enable','inactive','Background',[200 200 145]/255,'UserData',i)
    else
        set(hand.hEdits(i,1),'UserData',i)
    end
end
if (~isempty(hand.HdrButtons))         % Create the header pushbutton uicontrols
    for (j = 1:hand.NumCol-d_col)        % The d_col is to account for an eventual 'RowNumbers' option
        uicontrol('Style','pushbutton','unit','pixels','Enable','inactive','position', ...
            [arr_pos_xi(j+d_col) hand.arr_pos_y(1)+hand.RowHeight hand.ColWidth(j+d_col) 20],'String',hand.ColNames{j})
    end
end
if (~isempty(hand.RowNames))           % Create the header pushbutton uicontrols
    for (i = 1:length(hand.RowNames))
        uicontrol('Style','pushbutton','unit','pixels','Enable','inactive','position', ...
            hand.RowNames_pos_orig(i,:),'String',hand.RowNames{i}) 
   end
end

% ---------------- See if we need a slider ---------------------------
pos_t = get(hand.hEdits(1,hand.NumCol),'pos');       % Get top right edit position
pos_b = get(hand.hEdits(hand.MAX_ROWS,1),'pos');    % Get last visible edit position
if (hand.NumRows > hand.MAX_ROWS)
    set(hand.hEdits(hand.MAX_ROWS+1:hand.NumRows,1:hand.NumCol),'Visible','off')    % Hide those who are out of view
    if (~isempty(hand.checks))
        set(hand.hChecks(hand.MAX_ROWS+1:hand.NumRows),'Visible','off')
    end
    pos = [pos_t(1)+pos_t(3) pos_b(2) 15 pos_t(2)+pos_t(4)-pos_b(2)];
    sld_step = 1 / (hand.NumRows-1);
    sld_step(2) = 5 * sld_step(1);
    hand.hSlid = uicontrol('style','slider','units','pixels','position',pos,...
        'min',1,'max',hand.NumRows,'Value',hand.NumRows,'SliderStep',sld_step);
    set(hand.hSlid,'callback',{@slider_Callback,hand})
    set(hand.hSlid,'UserData',hand.NumRows)    % Store current value
end

% ---------------- See if the window is MODAL ---------------------------
if (~isempty(hand.modal))
	uicontrol('Style','pushbutton','unit','pixels','String','OK','position',...
        [pos_t(1)+pos_t(3)-110 5 40 20],'FontName','Helvetica','FontSize',9,...
        'callback','uiresume','tag','OK');
	uicontrol('Style','pushbutton','unit','pixels','String','Cancel','position', ...
        [pos_t(1)+pos_t(3)-60 5 60 20],'FontName','Helvetica','FontSize',9,...
        'callback','uiresume','tag','cancel');
    uiwait(hand.hFig)       % It also sets the Figure's visibility 'on'
    but = gco;
    if strcmp(get(but,'tag'),'OK')
        out = reshape(get(hand.hEdits,'String'),hand.NumRows,hand.NumCol);
        if (~isempty(hand.checks))
            if(length(hand.hChecks)~=1)
                unchecked = (cell2mat(get(hand.hChecks,'Value')) == 0);
                out(unchecked,:) = [];      % Remove unchecked rows
            end
        end
        if (~isempty(hand.RowNumbers)) % Do not output the row numbers
            out = out(:,2:end);
        end
        delete(hand.hFig)
    elseif strcmp(get(but,'tag'),'cancel')
        out = [];   delete(hand.hFig)
    else        % Figure was killed
        out = [];
    end
else
    set(hand.hFig,'Visible','on','UserData',hand)
    if (nargout),   out = hand.hFig;    end
end

% ---------------------------------------------------------------------------
function slider_Callback(obj,event,hand)

val = round(get(hand.hSlid,'Value'));
old_val = get(hand.hSlid,'UserData');
ds = val - old_val;

if (ds < 0)                                         % Slider moved down
    n = hand.NumRows - val + 1;    d_col = hand.NumRows - val;
    if (n+hand.MAX_ROWS-1 > hand.NumRows)             % Case we jumped into the midle zone
        adj = (n+hand.MAX_ROWS-1 - hand.NumRows);
        n = n - adj;    d_col = d_col - adj;
    end
    for (i = n:min(n+hand.MAX_ROWS-1,hand.NumRows))   % Update positions
        for (j = 1:hand.NumCol)
            pos = hand.Edits_pos_orig{i,j};
            set(hand.hEdits(i,j),'pos',[pos(1) hand.arr_pos_y(i-d_col) pos(3:4)],'Visible','on')
        end
        if (~isempty(hand.checks))                  % If we have checkboxes
            pos = hand.Checks_pos_orig(i,:);
            set(hand.hChecks(i),'pos',[pos(1) hand.arr_pos_y(i-d_col)+3 pos(3:4)],'Visible','on')            
        end
    end
    if (i == get(hand.hEdits(hand.NumRows,1),'UserData')) % Bottom reached. Jump to there
        val = 1;    set(hand.hSlid,'Value',val)         % This also avoids useless UIs repositioning
    end
elseif (ds > 0)                                     % Slider moved up
    n = hand.NumRows - val + 1;    k = hand.MAX_ROWS;
    if (n < hand.MAX_ROWS)                          % Case we jumped into the midle zone
        adj = (hand.MAX_ROWS - n - 0);
        n = n + adj;
    end
    for (i = n:-1:max(n-hand.MAX_ROWS+1,1))         % Update positions
        for (j = 1:hand.NumCol)
            pos = hand.Edits_pos_orig{i,j};
            set(hand.hEdits(i,j),'pos',[pos(1) hand.arr_pos_y(k) pos(3:4)],'Visible','on')        
        end
        if (~isempty(hand.checks))                  % If we have checkboxes
            pos = hand.Checks_pos_orig(i,:);
            set(hand.hChecks(i),'pos',[pos(1) hand.arr_pos_y(k)+3 pos(3:4)],'Visible','on')            
        end
        k = k - 1;
    end
    set(hand.hEdits(n+1:end,1:end),'Visible','off')
    if (~isempty(hand.checks)),     set(hand.hChecks(n+1:end),'Visible','off');    end
    if (i == get(hand.hEdits(1,1),'UserData'))      % Reached Top. Jump to there
        set(hand.hSlid,'Value',hand.NumRows)          % This also avoids useless UIs repositioning
        val = hand.NumRows;
    end
end
set(hand.hSlid,'UserData',val)                      % Save old 'Value'

% ----------------------------------------------------------------------------
function params = parse_pv_pairs(params,pv_pairs)
% parse_pv_pairs: parses sets of property value pairs, allows defaults
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.
%
% Example usage:
% First, set default values for the parameters. Assume we have four
% parameters that we wish to use optionally in the function examplefun.
%
%  - 'viscosity', which will have a default value of 1
%  - 'volume', which will default to 1
%  - 'pie' - which will have default value 3.141592653589793
%  - 'description' - a text field, left empty by default
%
% The first argument to examplefun is one which will always be supplied.
%
%   function examplefun(dummyarg1,varargin)
%   params.Viscosity = 1;
%   params.Volume = 1;
%   params.Pie = 3.141592653589793
%
%   params.Description = '';
%   params=parse_pv_pairs(params,varargin);
%   params
%
% Use examplefun, overriding the defaults for 'pie', 'viscosity'
% and 'description'. The 'volume' parameter is left at its default.
%
%   examplefun(rand(10),'vis',10,'pie',3,'Description','Hello world')
%
% params = 
%     Viscosity: 10
%        Volume: 1
%           Pie: 3
%   Description: 'Hello world'
%
% Note that capitalization was ignored, and the property 'viscosity' was truncated
% as supplied. Also note that the order the pairs were supplied was arbitrary.

n = length(pv_pairs) / 2;

if n ~= floor(n)
    error 'Property/value pairs must come in PAIRS.'
end
if (n <= 0),    return;     end     % just return the defaults

if ~isstruct(params)
    error 'No structure for defaults was supplied'
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for i=1:n
	p_i = lower(pv_pairs{2*i-1});
	v_i = pv_pairs{2*i};
	
	ind = strcmp(p_i, lpropnames);
    if (~any(ind))
	    ind = find(strncmp(p_i,lpropnames,length(p_i)));
        if isempty(ind)
            error(['No matching property found for: ',pv_pairs{2*i-1}])
	    elseif (length(ind) > 1)
            error(['Ambiguous property name: ',pv_pairs{2*i-1}])
        end
    end
    p_i = propnames{ind};
	params = setfield(params,p_i,v_i);      % override the corresponding default in params
end


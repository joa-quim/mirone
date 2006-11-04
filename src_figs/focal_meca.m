function varargout = focal_meca(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to focal_meca (see VARARGIN)

%	Copyright (c) 2004-2006 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

hObject = figure('Tag','figure1','Visible','off');
handles = guihandles(hObject);
guidata(hObject, handles);
focal_meca_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
movegui(hObject,'center');

global home_dir
f_data = [home_dir filesep 'data' filesep];

if (~isempty(varargin))
    handles.mirone_fig = varargin{1};    handles.MironeAxes = varargin{2};
    zz = get(varargin{2},'XLim');
    handles.x_min = zz(1);    handles.x_max = zz(2);
    zz = get(varargin{2},'YLim');
    handles.y_min = zz(1);    handles.y_max = zz(2);
end

handles.date = [];

% Import icons
load([f_data 'mirone_icons.mat'],'Mfopen_ico');
set(handles.pushbutton_readFile,'CData',Mfopen_ico)
clear Mfopen_ico;

% Fill the listbox fields with the currently available reading filters
%str = {'lon,lat,dep,strike,dip,rake,mag,[lon0,lat0,title]'; 'ISF formated catalog (ascii)';};
str = {'Aki & Richard''s convention file '; ...
        'Harvards''s CMT convention file '; ...
        'ISF formated catalog (ascii)';};
set(handles.listbox_readFilter,'String',str);
set(handles.checkbox_plotDate,'Enable','off')

% ------------- Give a Pro look (3D) to the frame boxes --------------------
bgcolor = get(0,'DefaultUicontrolBackgroundColor');
framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
set(0,'Units','pixels');    set(hObject,'Units','pixels')    % Pixels are easier to reason with
h_f = findobj(hObject,'Style','Frame');
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

% Choose default command line output for focal_meca_export
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes focal_meca_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = focal_meca_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = focal_meca_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------------------
function listbox_readFilter_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns listbox_readFilter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_readFilter
switch get(hObject,'Value')
    case 1
		str = sprintf(['ASCII file with lon,lat,depth,strike,dip,rake,mag.\n'...
            '8th and 9th columns are optional. If present, they\n'...
            'will determine where the beach ball will be ploted.']);
    case 2
		str = sprintf(['ASCII file with lon,lat,depth,strike1,dip1,rake1,\n'...
            'strike2,dip2,rake2,mantissa and exponent of moment in N-m.\n'...
            '12th and 13th columns are optional. If present, they\n'...
            'will determine where the beach ball will be ploted.']);
    case 3
		str = sprintf(['Read an ISF formated catalog file (like the ones\n'...
            'you can get from www.isc.ac.uk) and extract\n'...
            'the included (if any) focal mechanisms.']);
end
set(hObject,'TooltipString',str)

% -------------------------------------------------------------------------------------
function pushbutton_readFile_Callback(hObject, eventdata, handles)
% OK. Now read the earthquakes_export file and retain only the requested interval
item = get(handles.listbox_readFilter,'Value');     % Get the reading filter number
switch item
    case 1      % Read a file formated with the Aki & Richard convention
        str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
        filter = item;
    case 2      % Read a file formated with the CMT convention
        str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
        filter = item;
    case 3
        str1 = {'*.isf;*.ISF', 'Data files (*.isf,*.ISF)';'*.*', 'All Files (*.*)'};
        filter = item;
end

% Get file name
[FileName,PathName] = uigetfile(str1,'Select focal file');
pause(0.05)
if isequal(FileName,0)      return;    end
fname = [PathName,FileName];

handles.date = [];      % Allways reset

try
    set(gcf,'Pointer','watch')
	if (filter == 1 | filter == 2)      % Aki & Richard or CMT file
        [numeric_data,n_column,error] = read_file(fname);
        if (error)  return;     end
		% Get rid of events that are outside the map limits
		ind = find(numeric_data(:,1) < handles.x_min | numeric_data(:,1) > handles.x_max);
        numeric_data(ind,:) = [];
		ind = find(numeric_data(:,2) < handles.y_min | numeric_data(:,2) > handles.y_max);
        numeric_data(ind,:) = [];
        if (all(isempty(numeric_data)))     % Nothing inside region
            return
        end
        if (filter == 1)                % Aki & Richard
            if (~(n_column == 7 | n_column == 9 | n_column == 10))
                errordlg('Wrong number of columns for an A&R file','Error');    return
            end
            % [lon lat depth str1 dip1 rake1 mag]
            handles.data = numeric_data(:,1:7);
            mag = numeric_data(:,7);
            switch n_column
                case 7
                    handles.plot_pos = numeric_data(:,1:2);
                case 9
                    handles.plot_pos = numeric_data(:,8:9);
                case 10
                    handles.plot_pos = numeric_data(:,8:9);
                    %n = size(numeric_data,2);
                    %handles.date = cell(n,1);
                    %for (k=1:n)
                        %handles.date{k} = ??;
                    %end
            end
        else                % CMT convention
            if (~(n_column == 11 | n_column == 13 | n_column == 14))
                errordlg('Wrong number of columns for an CMT file','Error');    return
            end
            handles.mantiss_exp = numeric_data(:,10:11);
            mag = (log10(numeric_data(:,10)) + numeric_data(:,11) - 9.1) * 2 / 3;    % In fact Mw
            % [lon lat depth str1 dip1 rake1 str2 dip2 rake2 mag]
            handles.data = [numeric_data(:,1:9) mag];
            switch n_column
                case 11
                    handles.plot_pos = numeric_data(:,1:2);
                case 13
                    handles.plot_pos = numeric_data(:,12:13);
                case 14
                    handles.plot_pos = numeric_data(:,12:13);
            end
        end        		
    elseif (filter == 3)        % Read a ISF formated catalog
        opt_R = ['-R' num2str(handles.x_min) '/' num2str(handles.x_max) '/' num2str(handles.y_min) '/' num2str(handles.y_max)];
        [out_d,out_i] = read_isf(fname,opt_R,'-M');
        if (isempty(out_d))     % Nothing inside region
            return
        end
        handles.mantiss_exp = [out_d(10,:)' out_d(11,:)'];
        mag = (log10(out_d(10,:)) + out_d(11,:) - 9.1) * 2 / 3;    % In fact Mw
        handles.data = [out_d(1,:)' out_d(2,:)' out_d(3,:)' out_d(4,:)' out_d(5,:)' out_d(6,:)' out_d(7,:)' ...
                out_d(8,:)' out_d(9,:)' mag'];
        handles.plot_pos = [out_d(1,:)' out_d(2,:)'];
        clear out_d;
        n = size(out_i,2);
        handles.date = cell(n,1);
        for (k=1:n)
            handles.date{k} = [int2str_m(out_i(3,k)) '/' int2str_m(out_i(2,k)) '/' int2str_m(out_i(1,k))];
        end
        clear out_i;
        set(handles.checkbox_plotDate,'Enable','on')
	end
	handles.got_userFile = 1;
	handles.usr_DepthMin   = min(handles.data(:,3));    handles.usr_DepthMax = max(handles.data(:,3));
	handles.usr_MagMin     = min(mag);                  handles.usr_MagMax = max(mag);
	set(handles.edit_MagMin,'String',num2str(floor(handles.usr_MagMin)))
	set(handles.edit_MagMax,'String',num2str(ceil(handles.usr_MagMax)))
	set(handles.edit_DepthMin,'String',num2str(floor(handles.usr_DepthMin)))
	set(handles.edit_DepthMax,'String',num2str(ceil(handles.usr_DepthMax)))
	guidata(hObject,handles)
    set(gcf,'Pointer','arrow')
catch   % In case of error, set the pointer back to "normal" 
    set(gcf,'Pointer','arrow')
    w{1} = 'An error occured while reading file. Check that it has the apropriate format.';
    w{2} = '';
    w{3} = 'Alternatively, if you are sure that the format is correct check that there';
    w{4} = 'are no empty spaces at the end of your data lines. This may cause an';
    w{5} = 'error in decoding the ascii file.';
    warndlg(w,'Warning')
end

% -------------------------------------------------------------------------------------
function edit_MagMin_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < 1 | xx > 10)    set(hObject,'String','1');     end

% -------------------------------------------------------------------------------------
function edit_MagMax_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < 1 | xx > 10)    set(hObject,'String','10');    end

% -------------------------------------------------------------------------------------
function edit_Mag5_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < 0)    set(hObject,'String','1');   end

% -------------------------------------------------------------------------------------
function edit_DepthMin_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < 0)    set(hObject,'String','0');   end

% -------------------------------------------------------------------------------------
function edit_DepthMax_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx > 900)    set(hObject,'String','900');end

% -------------------------------------------------------------------------------------
function checkbox_depSlices_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.popup_dep0_33,'Enable','on');       set(handles.popup_dep33_70,'Enable','on')
    set(handles.popup_dep70_150,'Enable','on');     set(handles.popup_dep150_300,'Enable','on')
    set(handles.popup_dep300,'Enable','on')
else
    set(handles.popup_dep0_33,'Enable','off');      set(handles.popup_dep33_70,'Enable','off')
    set(handles.popup_dep70_150,'Enable','off');    set(handles.popup_dep150_300,'Enable','off')
    set(handles.popup_dep300,'Enable','off')
end

% -------------------------------------------------------------------------------------
function checkbox_plotDate_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox_plotDate

% -------------------------------------------------------------------------------------
function popup_dep0_33_Callback(hObject, eventdata, handles)
% Nothing to do here

% -------------------------------------------------------------------------------------
function popup_dep33_70_Callback(hObject, eventdata, handles)
% Nothing to do here

% -------------------------------------------------------------------------------------
function popup_dep70_150_Callback(hObject, eventdata, handles)
% Nothing to do here

% -------------------------------------------------------------------------------------
function popup_dep150_300_Callback(hObject, eventdata, handles)
% Nothing to do here

% -------------------------------------------------------------------------------------
function popup_dep300_Callback(hObject, eventdata, handles)
% Nothing to do here

% -------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
MagMin = str2double(get(handles.edit_MagMin,'String'));
MagMax = str2double(get(handles.edit_MagMax,'String'));
DepthMin = str2double(get(handles.edit_DepthMin,'String'));
DepthMax = str2double(get(handles.edit_DepthMax,'String'));
item = get(handles.listbox_readFilter,'Value');             % Get the reading filter number

if (isnan(MagMin))          MagMin = 1;         end
if (isnan(MagMax))          MagMax = 10;        end
if (isnan(DepthMin))        DepthMin = 0;       end
if (isnan(DepthMax))        DepthMax = 900;     end

if (handles.got_userFile)    % We have a user seismicity file
	% Retain only the requested interval
	ind1 = find(handles.data(:,3) < DepthMin | handles.data(:,3) > DepthMax);
    handles.data(ind1,:) = [];      handles.plot_pos(ind1,:) = [];
    if (~isempty(handles.date))     handles.date(ind1,:) = [];  end
    if (item == 1)                      % A&R.
		ind2 = find(handles.data(:,7) < MagMin | handles.data(:,7) > MagMax);
	elseif (item == 2 | item == 3)      % CMT file or ISF catalog
        handles.mantiss_exp(ind1,:) = [];   % This risked to heve been left behind
		ind2 = find(handles.data(:,10) < MagMin | handles.data(:,10) > MagMax);
        handles.mantiss_exp(ind2,:) = [];
    end
    handles.data(ind2,:) = [];      handles.plot_pos(ind2,:) = [];
    if (~isempty(handles.date))     handles.date(ind2,:) = [];  end
else
    errordlg('Plot What? Your christmas ballons?','Chico Clever');  return;
end

if (all(isempty(handles.data)))
    warndlg('There were no events left.','Warning');  return;
end

if (get(handles.checkbox_depSlices,'Value'))    % We have a depth slice request
    do_depSlices = 1;
    contents = get(handles.popup_dep0_33,'String');     cor_str{1} = contents{get(handles.popup_dep0_33,'Value')};
    contents = get(handles.popup_dep33_70,'String');    cor_str{2} = contents{get(handles.popup_dep33_70,'Value')};
    contents = get(handles.popup_dep70_150,'String');   cor_str{3} = contents{get(handles.popup_dep70_150,'Value')};
    contents = get(handles.popup_dep150_300,'String');  cor_str{4} = contents{get(handles.popup_dep150_300,'Value')};
    contents = get(handles.popup_dep300,'String');      cor_str{5} = contents{get(handles.popup_dep300,'Value')};
else
    do_depSlices = 0;
end

% ------------ OK, now we are ready to plot the mechanisms
oldunit = get(handles.MironeAxes,'Units');
set(handles.MironeAxes,'Units','centimeters')      % normalized
pos = get(handles.MironeAxes,'Position');
set(handles.MironeAxes,'Units',oldunit)
y_lim = get(handles.MironeAxes,'YLim');
handles.size_fac = (y_lim(2) - y_lim(1)) / (pos(4) - pos(2)) * 0.5;  % Scale facor
Mag5 = get(handles.edit_Mag5,'String');    % Size (cm) of a mag 5 event
handles.Mag5 = str2num(Mag5);
setappdata(handles.mirone_fig,'MecaMag5',Mag5)    % For eventual use in 'write_script'
n_meca = size(handles.data(:,1),1);
axes(handles.MironeAxes)
h_pat = zeros(n_meca,3);
plot_text = get(handles.checkbox_plotDate,'Value');
for (k=1:n_meca)
	if (item == 1)                      % Aki & Richard file
        [c,d] = patch_meca(handles.data(k,4), handles.data(k,5), handles.data(k,6));
        mag = handles.data(k,7);
	elseif (item == 2 | item == 3)      % CMT file or ISF catalog
        [c,d] = patch_meca(handles.data(k,4), handles.data(k,5), handles.data(k,6), ...
            handles.data(k,7), handles.data(k,8), handles.data(k,9));
        mag = handles.data(k,10);
	end
    dim = handles.size_fac * mag / 5 * handles.Mag5;    % Scale the balls against the selected Mag 5 size
    c = c * dim;    d = d * dim;
    cx = c(:,1) + handles.plot_pos(k,1);
    cy = c(:,2) + handles.plot_pos(k,2);
    dx = d(:,1) + handles.plot_pos(k,1);
    dy = d(:,2) + handles.plot_pos(k,2);
    h_pat(k,3) = line('XData',[handles.data(k,1) handles.plot_pos(k,1)], ...
        'YData',[handles.data(k,2) handles.plot_pos(k,2)], 'Linestyle','-', 'Marker','o', ...
        'MarkerSize',6, 'MarkerFaceColor','k', 'Tag','FocalMecaAnchor');
    if (~do_depSlices)      % Paint all compressive quadrants with black
        h_pat(k,1) = patch(cx,cy, [0 0 0],'Tag','FocalMeca');
    else
        cor = find_color(handles.data(k,3), cor_str);
        h_pat(k,1) = patch(cx,cy, cor,'Tag','FocalMeca');
    end
    h_pat(k,2) = patch(dx,dy, [1 1 1],'Tag','FocalMeca');
    if (plot_text)          % Plot event text identifier (normaly its date)
        offset = handles.size_fac * mag / 5 * (handles.Mag5 + 0.2);  % text offset regarding the beach ball (2 mm) 
        ht = text(handles.plot_pos(k,1),handles.plot_pos(k,2)+offset,handles.date{k},'HorizontalAlignment', ...
            'Center','VerticalAlignment','Bottom','FontSize',8,'Tag','TextMeca');
        draw_funs(ht,'DrawText');
        if (item == 1)                      % A&R. For eventual use in 'write_script'
            setappdata(h_pat(k,1),'psmeca_com',[handles.data(k,1:7) handles.plot_pos(k,1:2) ht]);
        elseif (item == 2 | item == 3)      % CMT or ISF catalog
            setappdata(h_pat(k,1),'psmeca_com',[handles.data(k,1:9) handles.mantiss_exp(k,:) handles.plot_pos(k,1:2) ht]);
        end
        setappdata(h_pat(k,1),'other_hand',[h_pat(k,2) h_pat(k,3) ht]); % For using in the uiedit
        setappdata(h_pat(k,2),'other_hand',[h_pat(k,1) h_pat(k,3) ht]); % For using in the uiedit
    else        % NO text info
        if (item == 1)                      % A&R. For eventual use in 'write_script'
            setappdata(h_pat(k,1),'psmeca_com',[handles.data(k,1:7) handles.plot_pos(k,1:2)]);
        elseif (item == 2 | item == 3)      % CMT or ISF catalog
            setappdata(h_pat(k,1),'psmeca_com',[handles.data(k,1:9) handles.mantiss_exp(k,:) handles.plot_pos(k,1:2)]);
        end
        setappdata(h_pat(k,1),'other_hand',[h_pat(k,2) h_pat(k,3)]);    % For using in the uiedit
        setappdata(h_pat(k,2),'other_hand',[h_pat(k,1) h_pat(k,3)]);    % For using in the uiedit
    end
    lim_x = [handles.plot_pos(k,1) handles.plot_pos(k,1) handles.plot_pos(k,1) handles.plot_pos(k,1)] + [-1 -1 1 1]*dim;
    lim_y = [handles.plot_pos(k,2) handles.plot_pos(k,2) handles.plot_pos(k,2) handles.plot_pos(k,2)] + [-1 1 1 -1]*dim;
    setappdata(h_pat(k,1),'Limits',[lim_x(:) lim_y(:)]);            % For using in the uiedit
    setappdata(h_pat(k,2),'Limits',[lim_x(:) lim_y(:)]);            % For using in the uiedit
    set_uicontext(h_pat(k,1));    set_uicontext(h_pat(k,2));
end
hand = guidata(handles.mirone_fig);     % Get the Mirone's handles structure
hand.have_focal = handles.Mag5;         % Signal that we have focal mechanisms and store the Mag5 size symbol
guidata(handles.mirone_fig,hand)        % Save the updated Mirone handles
delete(handles.figure1);

% -------------------------------------------------------------------------------------
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1);

% -------------------------------------------------------------------------------------
function cor = find_color(z, id)
if (z < 33)                 cor = id{1};
elseif (z >= 33 & z < 70)   cor = id{2};
elseif (z >= 70 & z < 150)  cor = id{3};
elseif (z >= 150 & z < 300) cor = id{4};
else                        cor = id{5};
end

% -------------------------------------------------------------------------------------
function [numeric_data,n_column,error] = read_file(fname)
error = 0;
hFig = gcf;
[bin,n_column,multi_seg,n_headers] = guess_file(fname);
% If msgbox exist we have to move it from behind the main window. So get it's handle
hMsgFig = gcf;
if (hFig ~= hMsgFig)        uistack(hMsgFig,'top');   end   % If msgbox exists, bring it forward
% If error in reading file
if isempty(bin) & isempty(n_column) & isempty(multi_seg) & isempty(n_headers)
    errordlg(['Error reading file ' fname],'Error');
    error = 1;  return
end

if (bin == 0)   % ASCII
    if (isempty(n_headers))     n_headers = NaN;    end
    if (multi_seg)
        [numeric_data,multi_segs_str,headerlines] = text_read(fname,NaN,n_headers,'>');
    else
        [numeric_data,multi_segs_str,headerlines] = text_read(fname,NaN,n_headers);
    end
    if (hFig ~= hMsgFig);       figure(hFig);   end     % gain access to the drawing figure
    if (iscell(numeric_data))
        n_segments = length(numeric_data);
    else
        n_segments = 1;
    end
else        % BINARY
    errordlg('Sorry, reading binary files is not yet programed','Error');
    error = 1;  return
end

% -------------------------------------------------------------------------------------
function set_uicontext(h)
% Set uicontexts to the Meca patches

cmenuHand = uicontextmenu;
set(h, 'UIContextMenu', cmenuHand);
uimenu(cmenuHand, 'Label', 'Delete this', 'Callback', {@del_Meca,h,'this'});
uimenu(cmenuHand, 'Label', 'Delete all', 'Callback', {@del_Meca,h,'all'});
%uimenu(cmenuHand, 'Label', 'Resize', 'Callback', {@resize_Meca,h});
ui_edit_patch_special(h)

% -------------------------------------------------------------------------------------
function del_Meca(obj,eventdata,h,opt)
% Delete one or all focal mechanisms
if (strcmp(opt,'this'))
    delete(getappdata(h,'other_hand'))
    delete(h)
else
    delete(findobj('Type','patch','Tag','FocalMeca'));
    delete(findobj('Type','line','Tag','FocalMecaAnchor'));
    delete(findobj('Type','text','Tag','TextMeca'));
end

% -------------------------------------------------------------------------------------
function resize_Meca(obj,eventdata,h)
% % Resize the focal mechanisms
% handles = guidata(gcf);
% h_all = findobj('Type','patch','Tag','FocalMeca');
% n_meca = length(h_all) / 2;     % Each ball has two patches
% mag = (log10(numeric_data(:,10)) + numeric_data(:,11) - 9.1) * 2 / 3;    % In fact Mw
% for (k=1:n_meca)
%     meca_com = getappdata(h_all(k),'psmeca_com');
%     dim = handles.size_fac * mag / 5 * handles.Mag5;    % Scale the balls against the selected Mag 5 size
%     c = c * dim;    d = d * dim;
%     cx = c(:,1) + handles.plot_pos(k,1);
%     cy = c(:,2) + handles.plot_pos(k,2);
%     dx = d(:,1) + handles.plot_pos(k,1);
%     dy = d(:,2) + handles.plot_pos(k,2);
% end

%--------------------------------------------------------------------------
function s = int2str_m(x)
%INT2STR Convert integer to string.
%   S = INT2STR(X) rounds the elements of the matrix X to
%   integers and converts the result into a string matrix.
%   Return NaN and Inf elements as strings 'NaN' and 'Inf', respectively.

%   Copyright 1984-2002 The MathWorks, Inc. 
%
% Hacked to work with TRUE integers as well

if (~isa(x,'double'))   % It's so simple. Now it works with true integers and not only
    x = double(x);      % the "pretend-to-be-integer-but-is-double" ML integers
end

x = round(real(x));
if (length(x) == 1)     % handle special case of single infinite or NaN element
   s = sprintf('%.1f',x);
   if (~strcmp(s, '-Inf') & ~strcmp(s, 'Inf') & ~strcmp(s, 'NaN'))
     s(end-1:end) = [];
   end
else
   s = '';
   [m,n] = size(x);
   % Determine elements of x that are finite.
   xfinite = x(isfinite(x));
   % determine the variable text field width quantity
   d = max(1,max(ceil(log10(abs(xfinite(:))+(xfinite(:)==0)))));
   clear('xfinite')
   % delimit string array with one space between all-NaN or all-Inf columns
   if any(isnan(x(:)))|any(isinf(x(:)))
      d = max([d;3]);
   end
   % walk through numbers array and convert elements to strings
   for i = 1:m
      t = [];
      for j = 1:n
         t = [t sprintf('%*.0f',d+2,x(i,j))];
      end
      s = [s; t];
   end
   % trim leading spaces from string array within constraints of rectangularity.
   if ~isempty(s)
      while all(s(:,1) == ' ')
         s(:,1) = []; 
      end
   end
end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')
    delete(handles.figure1);
end


% --- Creates and returns a handle to the GUI figure. 
function focal_meca_LayoutFcn(h1,handles);

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','Focal mechanisms',...
'NumberTitle','off',...
'Position',[520 445 390 355],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

h2 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'listbox_readFilter_Callback'},...
'Position',[50 290 251 61],...
'String',{  'Listbox' },...
'Style','listbox',...
'Value',1,...
'Tag','listbox_readFilter');

h3 = uicontrol('Parent',h1,...
'Callback',{@focal_meca_uicallback,h1,'pushbutton_readFile_Callback'},...
'FontWeight','bold',...
'Position',[300 310 23 23],...
'TooltipString','Browse for wanted file',...
'Tag','pushbutton_readFile');

h4 = uicontrol('Parent',h1,...
'Position',[10 196 371 80],...
'String',{  '' },...
'Style','frame',...
'Tag','frame1');

h5 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'edit_MagMin_Callback'},...
'Position',[71 243 47 21],...
'Style','edit',...
'TooltipString','Do not plot events weeker than this',...
'Tag','edit_MagMin');

h6 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'edit_MagMax_Callback'},...
'Position',[195 243 47 21],...
'Style','edit',...
'TooltipString','Do not plot events stronger than this',...
'Tag','edit_MagMax');

h7 = uicontrol('Parent',h1,...
'Position',[18 236 51 30],...
'String',{  'Minimum'; 'magnitude' },...
'Style','text',...
'Tag','text1');

h8 = uicontrol('Parent',h1,...
'Position',[143 237 51 30],...
'String',{  'Maximum'; 'magnitude' },...
'Style','text',...
'Tag','text2');

h9 = uicontrol('Parent',h1,...
'Position',[10 44 371 141],...
'String',{  '' },...
'Style','frame',...
'Tag','frame2');

h10 = uicontrol('Parent',h1,...
'Position',[19 140 42 30],...
'String',{  'Minimum'; 'depth' },...
'Style','text',...
'Tag','text3');

h11 = uicontrol('Parent',h1,...
'Position',[135 139 51 30],...
'String',{  'Maximum'; 'depth' },...
'Style','text',...
'Tag','text4');

h12 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'edit_Mag5_Callback'},...
'Position',[323 244 47 21],...
'String','0.8',...
'Style','edit',...
'TooltipString','The beach balls will be scaled to this value',...
'Tag','edit_Mag5');

h13 = uicontrol('Parent',h1,...
'Callback',{@focal_meca_uicallback,h1,'checkbox_plotDate_Callback'},...
'Position',[69 208 101 15],...
'String','Plot event date',...
'Style','checkbox',...
'TooltipString','Plot time information',...
'Tag','checkbox_plotDate');

h14 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'edit_DepthMin_Callback'},...
'Position',[63 145 47 21],...
'String','0',...
'Style','edit',...
'TooltipString','Do not plot events shalower than this',...
'Tag','edit_DepthMin');

h15 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'edit_DepthMax_Callback'},...
'Position',[187 145 47 21],...
'Style','edit',...
'TooltipString','Do not plot events deeper than this',...
'Tag','edit_DepthMax');

h16 = uicontrol('Parent',h1,...
'Callback',{@focal_meca_uicallback,h1,'checkbox_depSlices_Callback'},...
'Position',[19 110 211 15],...
'String','Use different colors for depth intervals',...
'Style','checkbox',...
'TooltipString','Destinguish the epicenter depths by color',...
'Tag','checkbox_depSlices');

h17 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'popup_dep0_33_Callback'},...
'Enable','off',...
'Position',[18 53 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'TooltipString','Symbol color for this depth interval',...
'Value',1,...
'Tag','popup_dep0_33');

h18 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[21 76 47 16],...
'String','0-33 km',...
'Style','text',...
'Tag','text5');

h19 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'popup_dep33_70_Callback'},...
'Enable','off',...
'Position',[91 53 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'TooltipString','Symbol color for this depth interval',...
'Value',2,...
'Tag','popup_dep33_70');

h20 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[94 76 54 16],...
'String','33-70 km',...
'Style','text',...
'Tag','text6');

h21 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'popup_dep70_150_Callback'},...
'Enable','off',...
'Position',[165 53 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'TooltipString','Symbol color for this depth interval',...
'Value',3,...
'Tag','popup_dep70_150');

h22 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[166 76 61 16],...
'String','70-150 km',...
'Style','text',...
'Tag','text7');

h23 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[312 76 55 16],...
'String','> 300 km',...
'Style','text',...
'Tag','text8');

h24 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'popup_dep150_300_Callback'},...
'Enable','off',...
'Position',[238 53 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'TooltipString','Symbol color for this depth interval',...
'Value',4,...
'Tag','popup_dep150_300');

h25 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'popup_dep300_Callback'},...
'Enable','off',...
'Position',[310 53 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'TooltipString','Symbol color for this depth interval',...
'Value',5,...
'Tag','popup_dep300');

h26 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[235 76 68 16],...
'String','150-300 km',...
'Style','text',...
'Tag','text9');

h27 = uicontrol('Parent',h1,...
'Callback',{@focal_meca_uicallback,h1,'pushbutton_Cancel_Callback'},...
'Position',[224 10 66 23],...
'String','Cancel',...
'Tag','pushbutton_Cancel');

h28 = uicontrol('Parent',h1,...
'Callback',{@focal_meca_uicallback,h1,'pushbutton_OK_Callback'},...
'Position',[315 10 66 23],...
'String','OK',...
'Tag','pushbutton_OK');

h29 = uicontrol('Parent',h1,...
'Position',[255 238 68 30],...
'String',{  'Magnitude 5'; 'size (cm)' },...
'Style','text',...
'Tag','text10');

function focal_meca_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

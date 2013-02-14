function out = gmtedit(varargin)
% Revival of the ancient Sunview gmtedit.
% gmtedit lets the user edit a gmtfile by mousing. The use should be
% fairly obvious from the menues and buttons' labels. Clicking on a point
% changes its status from good (green) to bad (red), a second click turns
% it back into a good one again.
%
%   Usage:
%       GMTEDIT with no arguments opens a empty GUI from where files can be loaded.
%       GMTEDIT(varargin)
%           All varargin(s) are optional
%           varargin{i} = 'FILE' opens the .gmt file 'FILE' (FILE must include the .gmt extension)
%           varargin{i} = '-L<width>' sets the <width> windowwidth in km (default is 300 km)
%           varargin{i} = '-G' force geodetical longitudes (0-360) output [Default is -180/+180];
%                           I'm sorry, this is against GMT defaults but I hate the [0-360] range.
%       OUT, if given will contain this figure handle.
%

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

f_name = '';
got_inFile = 0;
def_width_km = 200;         % Default width in km (approximatly 1/2 of a day at 10 knots)
opt_G = '-G';               % Default output to [-180/+180] range
begin = 1;                  % This means that the display starts at the begining of track
center_win = [];

if (nargin > 0)
    f_name = varargin{1};
    [PATH,FNAME,EXT] = fileparts(f_name);
    if (isempty(EXT))       % Remember that here we need the extension
        f_name = [f_name '.gmt'];
    end
    if (exist(f_name,'file') == 2)
        got_inFile = 1;
        varargin(1) = [];
    end
    for (k=1:length(varargin))
        switch (varargin{k}(1:2))
            case '-L'
                def_width_km = str2double(varargin{k}(3:end));
                if (isnan(def_width_km))        % Nonsense use of -L option
                    def_width_km = 200;
                end
            case '-G'
                opt_G = ' ';
            case '-P'
                str = varargin{k}(3:end);
                [tok,rem] = strtok(str,'/');
                lon = str2double(tok);
                lat = str2double(rem(2:end));
                begin = 0;
                center_win = [lon lat];
        end
    end
end

sc_size = get(0,'ScreenSize');
fig_height = fix(sc_size(4)*.9);
fp = [1 sc_size(4)-fig_height sc_size(3) fig_height];
marg_l = 60;                            % Left margin
marg_r = 15;                            % Right margin
marg_tb = 10;                           % Top & Bottom margins
marg_ax = 15;                           % Margin between axes
ax_height = fix(fig_height * .295);
ax_width = sc_size(3) - marg_l - marg_r;

hf = figure('name','gmtedit','resize','off','numbertitle','off', 'visible','off', ...
    'position',fp, 'DoubleBuffer','on', 'Tag','figure1', 'closerequestfcn','delete(gcbf)');
	
% Use system color scheme for figure:
set(hf,'Color',get(0,'defaultUicontrolBackgroundColor'));

% Apply this trick to get the icons. Let's hope that this is not version/OS dependent 
hA = findall(hf);
hh = findobj(hA,'TooltipString','Open File');
openFile_img = get(hh(1),'CData');
hh = findobj(hA,'TooltipString','Save Figure');
saveFile_img = get(hh(1),'CData');
hh = findobj(hA,'TooltipString','Zoom In');
zoomIn_img = get(hh(1),'CData');
hh = findobj(hA,'TooltipString','Zoom Out');
zoomOut_img = get(hh(1),'CData');
set(hf,'menubar','none')            % Set the menubar to none

pos = [marg_l fig_height-ax_height-marg_tb ax_width ax_height];
h_a1 = axes('Parent',hf, 'Units','pixels', 'Position',pos, 'XLim',[0 def_width_km], 'Tag','axes1');

pos = [marg_l fig_height-2*(ax_height+marg_tb)-marg_ax ax_width ax_height];
h_a2 = axes('Parent',hf, 'Units','pixels', 'Position',pos, 'XLim',[0 def_width_km], 'Tag','axes2');

pos = [marg_l fig_height-3*(ax_height+marg_tb)-2*marg_ax ax_width ax_height];
h_a3 = axes('Parent',hf,'Units','pixels', 'Position',pos, 'XLim',[0 def_width_km], 'Tag','axes3');

handles = guihandles(hf);
handles.opt_G = opt_G;
handles.home_dir = pwd;
handles.info = [];
handles.h_broken = [];
handles.begin = begin;
handles.center_win = center_win;

if (got_inFile)
    handles.last_dir = PATH;
else
    handles.last_dir = handles.home_dir;    % This means, last_dir is not saved between sessions
end

% Load some icons from mirone_icons.mat
load([handles.home_dir filesep 'data' filesep 'mirone_icons.mat'],'rectang_ico','info_ico');

%SetAxesNumericType(handles.axes1);        % Set axes uicontextmenus
set(get(h_a1,'YLabel'),'String','Gravity anomaly (mGal)')
set(get(h_a2,'YLabel'),'String','Magnetic anomaly (nT)')
set(get(h_a3,'YLabel'),'String','Bathymetry (m)')

handles.def_width_km = def_width_km;
handles.max_x_data = 10000;
scroll_plots(def_width_km,[0 10000])     % The [0 10000] will reset inside scroll_plots to a more appropriate val
movegui(hf,'north')

h_toolbar = uitoolbar('parent',hf,'Clipping', 'on', 'BusyAction','queue','HandleVisibility','on',...
   'Interruptible','on','Tag','FigureToolBar','Visible','on');
uipushtool('parent',h_toolbar,'Click',{@import_clickedcallback,f_name},'Tag','import',...
   'cdata',openFile_img,'TooltipString','Open gmt file');
uipushtool('parent',h_toolbar,'Click',@save_clickedcallback,'Tag','save',...
   'cdata',saveFile_img,'TooltipString','Save gmt file');
% uitoggletool('parent',h_toolbar,'Click',@zoom_clickedcallback,'Tag','zoom',...
%    'cdata',zoom_img,'TooltipString','Zoom');
uipushtool('parent',h_toolbar,'Click',@info_clickedcallback,'Tag','info','cdata',info_ico,...
   'TooltipString','Cruise Info');
uipushtool('parent',h_toolbar,'Click',@rectang_clickedcallback,'Tag','rectang','cdata',rectang_ico,...
   'TooltipString','Rectangular region','Separator','on');
uipushtool('parent',h_toolbar,'Click',@rectangMove_clickedcallback,'cdata',rectang_ico,'TooltipString','Select for moving');
uipushtool('parent',h_toolbar,'Click',{@changeScale_clickedcallback,'inc'},...
   'cdata',zoomIn_img,'TooltipString','Increase scale','Separator','on');
uipushtool('parent',h_toolbar,'Click',{@changeScale_clickedcallback,'dec'},...
   'cdata',zoomOut_img,'TooltipString','Decrease scale');


% Create empty lines just for the purpose of having their handles
handles.h_gl = line('XData',[],'YData',[],'Color','k','Parent',h_a1);
handles.h_gm = line('XData',[],'YData',[],'LineStyle','none','Marker','s', ...
    'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',4,'Parent',h_a1);
handles.h_ml = line('XData',[],'YData',[],'Color','k','Parent',h_a2);
handles.h_mm = line('XData',[],'YData',[],'LineStyle','none','Marker','s', ...
    'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',4,'Parent',h_a2);
handles.h_tl = line('XData',[],'YData',[],'Color','k','Parent',h_a3);
handles.h_tm = line('XData',[],'YData',[],'LineStyle','none','Marker','s', ...
    'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',4,'Parent',h_a3);

guidata(hf, handles);

% Now that we have a figure and its handles, we can open the .gmt file if it was requested on input
if (got_inFile)
    import_clickedcallback(hf,[],f_name)
end

% Add or remove red Markers
set(hf,'WindowButtonDownFcn',@add_MarkColor)

set(hf,'Visible','on')

% Choose default command line output for gmtedit
if (nargout == 1)   out = hf;   end

% --------------------------------------------------------------------
function import_clickedcallback(hObject, eventdata, opt)
handles = guidata(hObject);     % get handles
if (isempty(opt))
	[FileName,PathName] = put_or_get_file(handles,{'*.gmt;*.GMT', 'gmt files (*.gmt,*.GMT)'},'Select gmt File','get');
	if isequal(FileName,0),		return,		end
    f_name = [PathName FileName];
else
    f_name = opt;
end

[PATH,FNAME,EXT] = fileparts(f_name);
handles.f_name = [FNAME '.gmt'];
f_name = [PATH filesep FNAME];          % Rip the .gmt extension (we can't have here)
set(handles.figure1,'Name',['gmtedit  ' f_name])

set(handles.figure1,'Pointer','watch')
track = gmtlist_m(f_name,'-Fsxygmtd',handles.opt_G);
% Save those for use when saving into a new file
handles.time = track.time;
handles.lon = track.longitude;
handles.lat = track.latitude;
handles.year = track.year;
handles.info = track.info;
if (length(track.agency) ~= 10)         % Ensures that agency is exactly 10 chars
    agency = '          ';              % 10 blanks
    len = min(length(track.agency),10);
    agency(1:len) = track.agency(1:len);
    handles.agency = agency;
else
    handles.agency = track.agency;
end

% Search for de-activated (red marked) points handles - They may exist if another file was already loaded
h_gn = findobj(handles.figure1,'Type','Line','tag','GravNull');
h_mn = findobj(handles.figure1,'Type','Line','tag','MagNull');
h_tn = findobj(handles.figure1,'Type','Line','tag','TopNull');
if (~isempty(h_gn))     set(h_gn,'XData',[],'YData',[]);    end
if (~isempty(h_mn))     set(h_mn,'XData',[],'YData',[]);    end
if (~isempty(h_tn))     set(h_tn,'XData',[],'YData',[]);    end

% See if any "broken line" was left behind - They may exist if another file was already loaded and processed
if (~isempty(handles.h_broken))         % If we have broken lines, delete them
    for (k=1:length(handles.h_broken))
        set(handles.h_broken(k),'Xdata',[],'YData',[])
    end
    rmfield(handles,'h_broken');
    handles.h_broken = [];
end

if (~all(isnan(track.gravity)))
    set(handles.h_gm,'XData',track.distance,'YData',track.gravity, 'Tag','orig_grav')
else
    set(handles.h_gm,'XData',[],'YData',[])    
end
if (~all(isnan(track.magnetics)))
    set(handles.h_mm,'XData',track.distance,'YData',track.magnetics, 'Tag','orig_mag')
else
    set(handles.h_mm,'XData',[],'YData',[])    
end
if (~all(isnan(track.topography)))
    set(handles.h_tm,'XData',track.distance,'YData',track.topography, 'Tag','orig_topo')
else
    set(handles.h_tm,'XData',[],'YData',[])    
end

% Update the slider Max propertie
hs = findobj(handles.figure1,'style','slider');
handles.max_x_data = track.distance(end);
max_s = handles.max_x_data-handles.def_width_km;
if (max_s < 0)      % I already had one case like this. A very short track
    max_s = handles.max_x_data;
end
val = track.distance(1);
%set(hs,'Max',handles.max_x_data-handles.def_width_km,'value',track.distance(1))

if (~handles.begin)         % Start the display at a user selected coordinate
    x = handles.lon - handles.center_win(1);    [zz,id1] = min(abs(x));     clear x;
    y = handles.lat - handles.center_win(2);    [zz,id2] = min(abs(y));     clear y;
    % id1 and id2 are not forcedly equal. Find out the "best"
    r1 = sqrt((handles.lon(id1) - handles.center_win(1))^2 + (handles.lat(id1) - handles.center_win(2))^2);
    r2 = sqrt((handles.lon(id2) - handles.center_win(1))^2 + (handles.lat(id2) - handles.center_win(2))^2);
    id = id1;
    if (r1 ~= r2)
        if (r2 < r1)
            id = id2;
        end
    end
    x_lim = track.distance(id) + [-handles.def_width_km/2 handles.def_width_km/2];
    set(findall(gcf,'Type','axes'),'xlim',x_lim)
    val0 = track.distance(id)-handles.def_width_km;
    if (val0 > track.distance(1))
        val = val0;
    end
    %disp(['X_aqui = ' num2str(track.longitude(id)) '   Y_aqui = ' num2str(track.latitude(id))])
end

% Update the slider properties
set(hs,'Max',max_s,'value',val)
set(handles.figure1,'Pointer','arrow')
guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function save_clickedcallback(hObject, eventdata)
handles = guidata(hObject);     % get handles
NODATA = -32000;
[FileName,PathName] = put_or_get_file(handles, handles.f_name,'Select gmt File', 'put','.gmt');
if isequal(FileName,0),		return,		end
f_name = [PathName FileName];

% Search for de-activated points handles (the reds)
h_gn = findobj(handles.figure1,'Type','Line','tag','GravNull');
h_mn = findobj(handles.figure1,'Type','Line','tag','MagNull');
h_tn = findobj(handles.figure1,'Type','Line','tag','TopNull');

% And the corresponding red marker values values
x_gn = get(h_gn,'XData');       y_gn = get(h_gn,'YData');
x_mn = get(h_mn,'XData');       y_mn = get(h_mn,'YData');
x_tn = get(h_tn,'XData');       y_tn = get(h_tn,'YData');

% Get G,M,T values
x_g = get(handles.h_gm,'XData');         y_g = get(handles.h_gm,'YData');
x_m = get(handles.h_mm,'XData');         y_m = get(handles.h_mm,'YData');
x_t = get(handles.h_tm,'XData');         y_t = get(handles.h_tm,'YData');

if (~isempty(handles.h_broken))         % If we have broken lines we must join them
    x_broken = [];    y_broken = [];
    for (k=1:length(handles.h_broken))
        x_broken = [x_broken get(handles.h_broken(k),'XData')];
        y_broken = [y_broken get(handles.h_broken(k),'YData')];
    end
    [x_m,id] = sort([x_m x_broken]);
    y_m = [y_m y_broken];
    y_m = y_m(id);                      % Otherwise the y's would be out of order
end

set(gcf,'Pointer','watch')
x_lim = get(gca,'XLim');
if (~isempty(x_gn))
    id_x = zeros(length(x_gn),1);
    for (k=1:length(x_gn))
        tmp = find((x_g - x_gn(k)) == 0);   % Find the gravity points that were marked
        id_x(k) = tmp(1);               % Old files often have repeated coords
    end
    y_g(id_x) = NODATA;                 % Remove them
end
if (~isempty(x_mn))
    id_x = zeros(length(x_mn),1);
    for (k=1:length(x_mn))
        tmp = find((x_m - x_mn(k)) == 0);   % Find the magnetic points that were marked
        id_x(k) = tmp(1);               % Old files often have repeated coords
    end
    y_m(id_x) = NODATA;                 % Remove them
end
if (~isempty(x_tn))
    id_x = zeros(length(x_tn),1);
    for (k=1:length(x_tn))
        tmp = find((x_t - x_tn(k)) == 0);   % Find the topo points that were marked
        id_x(k) = tmp(1);               % Old files often have repeated coords
    end
    y_t(id_x) = NODATA;                 % Remove them
end

n_rec = length(handles.lon);
if (isempty(y_g))                       % Original file had no gravity data
    y_g = repmat(int16(NODATA),1,n_rec);
else                                    % Replace eventual NaNs with NODATA
    y_g(isnan(y_g*10)) = NODATA;        % But before convert to GU units
    y_g = int16(y_g);
end
if (isempty(y_m))                       % Original file had no magnetic data
    y_m = repmat(int16(NODATA),1,n_rec);
else                                    % Replace eventual NaNs with NODATA
    y_m(isnan(y_m)) = NODATA;
    y_m = int16(y_m);
end
if (isempty(y_t))                       % Original file had no topography data
    y_t = repmat(int16(NODATA),1,n_rec);
else                                    % Replace eventual NaNs with NODATA
    y_t(isnan(y_t)) = NODATA;
    y_t = int16(y_t);
end

tempo = int32(handles.time);
lat  = int32(handles.lat * 1e6);    % And convert back to millidegrees
lon  = int32(handles.lon * 1e6);

fid = fopen(f_name,'wb');
fwrite(fid,[int32(handles.year) n_rec],'int32');
fwrite(fid,handles.agency,'schar');
% This is STUPIDLY slow but I didn't find any other way to do it.
for (k=1:n_rec)
    fwrite(fid,[tempo(k) lat(k) lon(k)],'int32');
    fwrite(fid,[y_g(k) y_m(k) y_t(k)],'int16');
end
fclose(fid);

set(gcf,'Pointer','arrow')

% --------------------------------------------------------------------
function add_MarkColor(hObject, eventdata)
% Add a red Marker over the closest (well, near closest) clicked point.
handles = guidata(hObject);     % get handles

button = get(handles.figure1, 'SelectionType');
if (~strcmp(button,'normal'))     return;     end     % Accept only left-clicks

in_grav = 0;    in_mag = 0;     in_topo = 0;
ax = gca;
pt = get(ax, 'CurrentPoint');
if (strcmp(get(ax,'Tag'),'axes1'))
    in_grav = 1;        opt = 'GravNull';
elseif (strcmp(get(ax,'Tag'),'axes2'))
    in_mag = 1;         opt = 'MagNull';
elseif (strcmp(get(ax,'Tag'),'axes3'))
    in_topo = 1;        opt = 'TopNull';
end
if ((in_grav + in_mag + in_topo) == 0)      return;     end     % Click was outside axes

if (in_grav)
    hM = findobj(handles.figure1,'Type','Line','tag','GravNull');
    x = get(handles.h_gm,'XData');      y = get(handles.h_gm,'YData');
elseif (in_mag)
    hM = findobj(handles.figure1,'Type','Line','tag','MagNull');
    x = get(handles.h_mm,'XData');      y = get(handles.h_mm,'YData');
else
    hM = findobj(handles.figure1,'Type','Line','tag','TopNull');
    x = get(handles.h_tm,'XData');      y = get(handles.h_tm,'YData');
end

x_lim = get(ax,'XLim');                 y_lim = get(ax,'YLim');
dx = diff(x_lim) / 20;                  % Search only betweem +/- 1/10 of x_lim
id = (x < (pt(1,1)-dx) | x > (pt(1,1)+dx));
x(id) = [];             y(id) = [];     % Clear outside-2*dx points to speed up the search code
XScale = diff(x_lim);       YScale = diff(y_lim);

r = sqrt(((pt(1,1)-x)./XScale).^2+((pt(1,2)-y)./YScale).^2);
[temp,i] = min(r);
pt_x = x(i);                pt_y = y(i);

xr = get(hM,'XData');       yr = get(hM,'YData');
id = find(xr == pt_x);
if (isempty(id))            % New Marker
    if (isempty(hM))        % First red Marker on this axes
        line(pt_x,pt_y,'Marker','s','MarkerFaceColor','r','MarkerSize',4,'LineStyle','none','Tag',opt);
    else
        xr = [xr pt_x];     yr = [yr pt_y];
        set(hM,'XData',xr, 'YData', yr)
    end
else                        % Marker already exists. Kill it
    xr(id) = [];            yr(id) = [];
    set(hM,'XData',xr, 'YData', yr)
end

% --------------------------------------------------------------------------------------------------
function zoom_clickedcallback(obj,eventdata)
if (strcmp(get(obj,'State'),'on'))
    zoom_j xon;
else
    zoom_j off;
end

% --------------------------------------------------------------------------------------------------
function info_clickedcallback(obj,eventdata)
handles = guidata(obj);     % get handles
if (isempty(handles.info))      return;     end
data = handles.info;
str{1} = ['N_recs = ' num2str(data(1)) ', N_grav = ' num2str(data(2)) ', N_mag = ' num2str(data(3)) ...
    ', N_top = ' num2str(data(4))];
str{2} = ['E: = ' num2str(data(5)) '  W: = ' num2str(data(6))];
str{3} = ['S: = ' num2str(data(7)) '  N: = ' num2str(data(8))];
str{4} = ['Start day,month,year: = ' num2str(data(9)) '  ' num2str(data(10)) '  ' num2str(data(11))];
str{5} = ['End   day,month,year: = ' num2str(data(12)) '  ' num2str(data(13)) '  ' num2str(data(14))];
msgbox(str,'Cruise Info')

% --------------------------------------------------------------------------------------------------
function rectang_clickedcallback(obj,eventdata)
handles = guidata(obj);     % get handles
try
    [p1,p2] = rubberbandbox;
catch       % Don't know why but uisuspend sometimes breaks
    set(gcf,'Pointer','arrow');
    return
end

in_grav = 0;    in_mag = 0;     in_topo = 0;
ax = gca;
if (strcmp(get(ax,'Tag'),'axes1'))
    in_grav = 1;        opt = 'GravNull';
elseif (strcmp(get(ax,'Tag'),'axes2'))
    in_mag = 1;         opt = 'MagNull';
elseif (strcmp(get(ax,'Tag'),'axes3'))
    in_topo = 1;        opt = 'TopNull';
end

if (in_grav)
    hM = findobj(handles.figure1,'Type','Line','tag','GravNull');
    x = get(handles.h_gm,'XData');      y = get(handles.h_gm,'YData');
elseif (in_mag)
    hM = findobj(handles.figure1,'Type','Line','tag','MagNull');
    x = get(handles.h_mm,'XData');      y = get(handles.h_mm,'YData');
else
    hM = findobj(handles.figure1,'Type','Line','tag','TopNull');
    x = get(handles.h_tm,'XData');      y = get(handles.h_tm,'YData');
end

% Search for "regular" points inside the rectangle
id = find(x >= p1(1,1) & x <= p2(1,1) & y >= p1(1,2) & y <= p2(1,2));
if (isempty(id))    return;     end     % Nothing inside rect
golo_x = x(id);             golo_y = y(id);

% Now search also for eventual pre-existing red markers inside the rectangle
xr = get(hM,'XData');       yr = get(hM,'YData');
if (~isempty(xr))           % We have red markers. See if any is inside rect
    id = find(xr >= p1(1,1) & xr <= p2(1,1) & yr >= p1(1,2) & yr <= p2(1,2));
    if (~isempty(id))       % Yes we have, so remove them for not creating duplicates
        xr(id) = [];        yr(id) = [];
    end
end

% Finaly plot the points cought inside the rect
if (isempty(hM))        % First red Markers on this axes
    line(golo_x,golo_y,'Marker','s','MarkerFaceColor','r','MarkerSize',4,'LineStyle','none','Tag',opt);
else
    golo_x = [golo_x xr];    golo_y = [golo_y yr];
    set(hM,'XData',golo_x, 'YData', golo_y)
end

% --------------------------------------------------------------------------------------------------
function rectangMove_clickedcallback(obj,eventdata)
handles = guidata(obj);     % get handles
try
    [p1,p2] = rubberbandbox;
catch
    return
end

if (~strcmp(get(gca,'Tag'),'axes2'))     % This option is meant only to magnetic data
    return
end

x = get(handles.h_mm,'XData');      y = get(handles.h_mm,'YData');

% Search for "regular" points inside the rectangle
id = find(x >= p1(1,1) & x <= p2(1,1) & y >= p1(1,2) & y <= p2(1,2));
if (isempty(id))    return;     end     % Nothing inside rect

% create a new line and set it the line_uicontext in order that they may be moved
golo_x = x(id);             golo_y = y(id);

x(id) = [];         y(id) = [];         % Remove these from the main line
set(handles.h_mm,'XData',x,'YData',y)

n = length(handles.h_broken);
handles.h_broken(n+1) = line(golo_x,golo_y,'Marker','s','MarkerFaceColor','b','MarkerSize',4,'LineStyle','-');
ui_edit_polygon(handles.h_broken(n+1))      % Set edition functions

guidata(obj,handles)

% --------------------------------------------------------------------------------------------------
function changeScale_clickedcallback(obj,eventdata,opt)
handles = guidata(obj);     % get handles

if (strcmp(opt,'inc'))
    handles.def_width_km = handles.def_width_km - 50;
    if (handles.def_width_km < 50)
        return
    end
else            % Decrease scale
    handles.def_width_km = handles.def_width_km + 50;
end

x_lim = get(gca,'Xlim');
set(findall(gcf,'Type','axes'),'Xlim',x_lim(1)+[0 handles.def_width_km])

h_slider = findobj(handles.figure1,'style','slider');
cb = get(h_slider,'callback');
% Here we make use of the knowledge that the "cb" is a string of the form:
% set(findall(gcf,'Type','axes'),'xlim',get(gcbo,'value')+[0 ???])
% where '???' is the width of the currently displyed axes. And that's what we need to change
new_cb = [cb(1:59) num2str(handles.def_width_km) '])'];
set(h_slider,'callback',new_cb)

% Now update the slider 'Max' propertie
new_max = handles.max_x_data - handles.def_width_km;
set(h_slider,'Max',new_max)

val = get(h_slider,'Value');
if (val > new_max)
    set(h_slider,'Value',new_max)
end

min_s = get(h_slider,'Min');
max_s = get(h_slider,'Max');
if (min_s >= max_s)
    disp(['DEBUG. Error in slider. Min = ' num2str(min_s) ' Max = ' num2str(max_s)])
    aa=0;
end

guidata(obj,handles)

% --------------------------------------------------------------------
function scroll_plots(width,x);
%   This function uses the idea of the scrollplotdemo from Steven Lord (http://www.mathworks.com/matlabcentral)
%   scroll_plots(width,x);
%   width: window width in x units
%   x: absicssae vector

% Generate constants for use in uicontrol initialization
pos = get(gca,'position');

% This will create a slider which is just underneath the axis
% but still leaves room for the axis labels above the slider
Newpos=[pos(1) 5 pos(3) 20];

S=['set(findall(gcf,''Type'',''axes''),''xlim'',get(gcbo,''value'')+[0 ' num2str(width) '])'];

% Creating Uicontrol with initial value of the minimum of x
uicontrol('style','slider','units','pixels','position',Newpos, ...
    'callback',S,'min',x(1),'max',x(end)-width,'value',x(1));

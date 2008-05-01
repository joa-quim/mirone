function varargout = animate_seismicity(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to animate_seismicity (see VARARGIN)

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
animate_seismicity_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
movegui(hObject,'center')

if (~isempty(varargin))
    handles.mirone_fig = varargin{1};
    handles.h_events = varargin{2};
else
    delete(hObject)
    return
end

handles.events_time = getappdata(handles.h_events,'SeismicityTime');
if (isempty(handles.events_time))
    delete(hObject)
    return
end

min_t = min(handles.events_time);
max_t = max(handles.events_time);
handles.StartYear = fix(min_t);
handles.EndYear = fix(max_t);

% Find start and end Month & Day
tmp = min_t - fix(min_t);
jd0 = fix(tmp * (365 + isleapyear(handles.StartYear))) + 1;
tmp = max_t - fix(max_t);
jd1 = fix(tmp * (365 + isleapyear(handles.EndYear))) + 1;
[handles.StartMonth,handles.StartDay] = jd2monday(jd0,handles.StartYear);
[handles.EndMonth,handles.EndDay] = jd2monday(jd1,handles.EndYear);

% Fill the time editboxes
set(handles.edit_StartYear,'String',num2str(handles.StartYear))
set(handles.edit_StartMonth,'String',num2str(handles.StartMonth))
set(handles.edit_StartDay,'String',num2str(handles.StartDay))
set(handles.edit_EndYear,'String',num2str(handles.EndYear))
set(handles.edit_EndMonth,'String',num2str(handles.EndMonth))
set(handles.edit_EndDay,'String',num2str(handles.EndDay))

% Set the slider vars
sl_step = [1 10]/((max_t - min_t)*365);
if (sl_step > 1)        % It happens when we have a < one day swarm
    sl_step(1) = 1;    sl_step(2) = 1;
end
set(handles.slider_time,'Min',min_t,'Max',max_t,'Value',min_t,'SliderStep',sl_step)

% Make a copy of the event coordinates. We need a copy because the h_events will be changed often
handles.x_bak = get(handles.h_events,'XData');      handles.y_bak = get(handles.h_events,'YData');

handles.fig_name = get(hObject,'Name');

% Choose default command line output for animate_seismicity_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes animate_seismicity_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = animate_seismicity_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = animate_seismicity_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% -----------------------------------------------------------------------------------------
function edit_StartYear_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < handles.StartYear)
    set(hObject,'String',num2str(handles.StartYear))
end

% -----------------------------------------------------------------------------------------
function edit_StartMonth_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < 1)
    set(hObject,'String',num2str(handles.StartMonth))
end

% -----------------------------------------------------------------------------------------
function edit_StartDay_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < 1)
    set(hObject,'String',num2str(handles.StartDay))
end

% -----------------------------------------------------------------------------------------
function edit_EndYear_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx > handles.EndYear)
    set(hObject,'String',num2str(handles.EndYear))
end

% -----------------------------------------------------------------------------------------
function edit_EndMonth_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx > 12)
    set(hObject,'String',num2str(handles.EndMonth))
end

% -----------------------------------------------------------------------------------------
function edit_EndDay_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx > 31)
    set(hObject,'String',num2str(handles.EndDay))
end

% -----------------------------------------------------------------------------------------
function edit_timeStep_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < 0)
    set(hObject,'String','2')
end

% -----------------------------------------------------------------------------------------
function edit_lifeSpan_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < str2double(get(handles.edit_timeStep,'String')))
    set(hObject,'String','10')
end

% -----------------------------------------------------------------------------------------
function edit_frameDelay_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) | xx < 0)
    set(hObject,'String','0.25')
end

% -----------------------------------------------------------------------------------------
function pushbutton_stop_Callback(hObject, eventdata, handles)
set(handles.h_events,'XData',handles.x_bak,'YData',handles.y_bak);
delete(handles.figure1);

% -----------------------------------------------------------------------------------------
function pushbutton_run_Callback(hObject, eventdata, handles)

StartYear = str2double(get(handles.edit_StartYear,'String'));
EndYear = str2double(get(handles.edit_EndYear,'String'));
StartMonth = str2double(get(handles.edit_StartMonth,'String'));
EndMonth = str2double(get(handles.edit_EndMonth,'String'));
StartDay = str2double(get(handles.edit_StartDay,'String'));
EndDay = str2double(get(handles.edit_EndDay,'String'));
dt = str2double(get(handles.edit_timeStep,'String')) / (365 + isleapyear(StartYear));
frame_int = str2double(get(handles.edit_frameDelay,'String'));
vida = str2double(get(handles.edit_lifeSpan,'String'));
killed = 0;

x = handles.x_bak;      y = handles.y_bak;

% Retain only the requested time interval
lower_date = dec_year(StartYear,StartMonth,StartDay);
upper_date = dec_year(EndYear,EndMonth,EndDay+0.999);   % + 0.999 to use the entire current day
ind = (handles.events_time < lower_date | handles.events_time > upper_date);
handles.events_time(ind) = [];      y(ind) = [];      x(ind) = [];

t_up = min(handles.events_time);
jd = fix((t_up - StartYear) * (365 + isleapyear(StartYear)));
t_up = StartYear + jd / (365 + isleapyear(StartYear));
t_end = max(handles.events_time);
vida = vida / (365 + isleapyear(handles.StartYear));
set(handles.h_events,'XData',[],'YData',[])
dbf = get(handles.mirone_fig,'Renderer');
set(handles.mirone_fig,'Renderer','zbuffer')

while (t_up <= t_end)
    t_up = t_up + dt;
    t_low = t_up - vida;
    year = fix(t_up);
    jd = fix((t_up - year) * (365 + isleapyear(handles.StartYear))) + 1;
    [month,day] = jd2monday(jd,year);
    tit = [handles.fig_name '  ' num2str(day,'%02d') '/' num2str(month,'%02d') '/' num2str(year)];
    ind = find(handles.events_time <= t_up & handles.events_time >= t_low);
    set(handles.h_events,'XData',x(ind),'YData',y(ind));
    pause(frame_int)
    try
        set(handles.figure1,'Name',tit)
    catch
        killed = 1;
        break
    end
end
set(handles.h_events,'XData',handles.x_bak,'YData',handles.y_bak);
set(handles.mirone_fig,'Renderer',dbf)
if (~killed)
    set(handles.figure1,'Name',handles.fig_name)
end

%--------------------------------------------------------------------------
function slider_time_Callback(hObject, eventdata, handles)
t = get(hObject,'Value');
year = fix(t);
vida = str2double(get(handles.edit_lifeSpan,'String')) / (365 + isleapyear(year));
dbf = get(handles.mirone_fig,'Renderer');
set(handles.mirone_fig,'Renderer','zbuffer')

jd = fix((t - year) * (365 + isleapyear(year))) + 1;
t = year + jd / (365 + isleapyear(year));   % Count to the end of the day
[month,day] = jd2monday(jd,year);
ind = find( (handles.events_time <= t) & (handles.events_time >= t-vida) );
set(handles.h_events,'XData',handles.x_bak(ind),'YData',handles.y_bak(ind));
tit = [handles.fig_name '  ' num2str(day,'%02d') '/' num2str(month,'%02d') '/' num2str(year)];
set(handles.figure1,'Name',tit)
set(handles.mirone_fig,'Renderer',dbf)

% -----------------------------------------------------------------------------------------
function pushbutton_help_Callback(hObject, eventdata, handles)

% -----------------------------------------------------------------------------
function [month, day] = jd2monday(jday,year)
%JD2MONDAY Julian day number to Julian calendar date.
%
%   [MONTH, DAY] = JD2MONDAY(JDAY,YEAR) returns the
%   Julian calendar date (month, day) corresponding to the Julian day JDAY.

%   Author:      Peter J. Acklam
%   Hacked to work with the "fake" JD that start at 1 at 1fst January of each year.

t = ( ~rem(year, 4) & rem(year, 100) ) | ~rem(year, 400);       % Check for leap-years
tt = (~t & jday > 59);
jday(tt) = jday(tt) + 1;            % Trick to make this algo work also for non leap-years

c = jday + 32081;
d = floor((4 * c + 3) / 1461);
e = c - floor((1461 * d) / 4);
m = floor((5 * e + 2) / 153);

day   = e - floor((153 * m + 2) / 5) + 1;
month = m + 3 - 12 * floor(m / 10);

% -----------------------------------------------------------------------------
function yd = dec_year(varargin)
%   DEC_YEAR(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) returns the ordinal year
%   number plus a fractional part depending on the month, day, and time of day
%
%   Any missing MONTH or DAY will be replaced by 1.  HOUR, MINUTE or SECOND
%   will be replaced by zeros.

%   Adapted from timeutil functions of Peter J. Acklam by Joaquim Luis

argv = { 1 1 1 0 0 0 };
argv(1:nargin) = varargin;
[year, month, day, hour, minute, second] = deal(argv{:});

days_in_prev_months = [0 31 59 90 120 151 181 212 243 273 304 334]';

% Day in given month.
try
    yd = days_in_prev_months(month) ...               % days in prev. months
         + ( isleapyear(year) & ( month > 2 ) ) ...   % leap day
         + day ...                                    % day in month
         + ( second + 60*minute + 3600*hour )/86400;  % part of day
catch
    yd = [];    return
end
yd = year + (yd - 1) ./ (365 + isleapyear(year));

%--------------------------------------------------------------------------
function t = isleapyear(year)
t = ( ~rem(year, 4) & rem(year, 100) ) | ~rem(year, 400);

%--------------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);


% --- Creates and returns a handle to the GUI figure. 
function animate_seismicity_LayoutFcn(h1,handles);

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Seismicity Swing',...
'NumberTitle','off',...
'Position',[520 630 461 170],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

h2 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@animate_seismicity_uicallback,h1,'edit_StartYear_Callback'},...
'Position',[50 137 47 21],...
'Style','edit',...
'Tag','edit_StartYear');

h3 = uicontrol('Parent',h1,'Position',[10 127 36 33],'String',{'Start'; 'year'},'Style','text','Tag','text1');

h4 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@animate_seismicity_uicallback,h1,'edit_StartMonth_Callback'},...
'Position',[174 137 47 21],...
'Style','edit',...
'Tag','edit_StartMonth');

h5 = uicontrol('Parent',h1,'Position',[130 131 41 30],'String',{'Start'; 'month'},'Style','text','Tag','text2');

h6 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@animate_seismicity_uicallback,h1,'edit_StartDay_Callback'},...
'Position',[300 137 47 21],...
'Style','edit',...
'Tag','edit_StartDay');

h7 = uicontrol('Parent',h1,'Position',[264 131 36 30],'String',{'Start'; 'day'},'Style','text','Tag','text3');

h8 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@animate_seismicity_uicallback,h1,'edit_EndYear_Callback'},...
'Position',[50 103 47 21],...
'Style','edit',...
'Tag','edit_EndYear');

h9 = uicontrol('Parent',h1,'Position',[10 93 36 30],'String',{'End'; 'year'},'Style','text','Tag','text4');

h10 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@animate_seismicity_uicallback,h1,'edit_EndMonth_Callback'},...
'Position',[174 103 47 21],...
'Style','edit',...
'Tag','edit_EndMonth');

h11 = uicontrol('Parent',h1,'Position',[130 94 41 30],'String',{'End'; 'month'},'Style','text','Tag','text5');

h12 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@animate_seismicity_uicallback,h1,'edit_EndDay_Callback'},...
'Position',[300 103 47 21],...
'Style','edit',...
'Tag','edit_EndDay');

h13 = uicontrol('Parent',h1,'Position',[270 98 27 30],'String',{'End'; 'day'},'Style','text','Tag','text6');

h14 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@animate_seismicity_uicallback,h1,'edit_timeStep_Callback'},...
'Position',[50 38 47 21],...
'String','2',...
'Style','edit',...
'TooltipString','Pick events at this time rate',...
'Tag','edit_timeStep');

h15 = uicontrol('Parent',h1,'Position',[30 63 91 15],'String','Time step (days)','Style','text','Tag','text7');

h16 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@animate_seismicity_uicallback,h1,'edit_lifeSpan_Callback'},...
'Position',[174 38 47 21],...
'String','10',...
'Style','edit',...
'TooltipString','Events will stay on screen for these number of days',...
'Tag','edit_lifeSpan');

h17 = uicontrol('Parent',h1,'Position',[151 63 91 15],'String','Life span (days)','Style','text','Tag','text8');

h18 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@animate_seismicity_uicallback,h1,'edit_frameDelay_Callback'},...
'Position',[300 38 47 21],...
'String','0.25',...
'Style','edit',...
'TooltipString','Next event(s) will show after this interval (seconds)',...
'Tag','edit_frameDelay');

h19 = uicontrol('Parent',h1,'Position',[282 63 79 15],'String','Frame int (sec)','Style','text','Tag','text9');

h20 = uicontrol('Parent',h1,...
'Callback',{@animate_seismicity_uicallback,h1,'pushbutton_stop_Callback'},...
'Position',[385 8 66 23],...
'String','Stop',...
'Tag','pushbutton_stop');

h21 = uicontrol('Parent',h1,...
'Callback',{@animate_seismicity_uicallback,h1,'pushbutton_run_Callback'},...
'Position',[384 113 66 23],...
'String','Run',...
'Tag','pushbutton_run');

h22 = uicontrol('Parent',h1,...
'Callback',{@animate_seismicity_uicallback,h1,'pushbutton_help_Callback'},...
'FontSize',10,...
'FontWeight','bold',...
'Position',[403 61 23 23],...
'String','?',...
'Tag','pushbutton_help');

h23 = uicontrol('Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',{@animate_seismicity_uicallback,h1,'slider_time_Callback'},...
'Position',[49 8 302 15],...
'String',{''},...
'Style','slider',...
'TooltipString','Show events in the interval [time-Life span; time]',...
'Tag','slider_time');

function animate_seismicity_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

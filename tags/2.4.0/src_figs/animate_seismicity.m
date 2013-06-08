function varargout = animate_seismicity(varargin)
% Make a seismicity movie

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
	if (isempty(varargin))		return,		end

	hObject = figure('Tag','figure1','Visible','off');
	animate_seismicity_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	handles.mirone_fig = varargin{1};
	handles.h_events = varargin{2};

	handles.events_time = getappdata(handles.h_events,'SeismicityTime');
	if (isempty(handles.events_time))
		delete(hObject),	return
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

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.mirone_fig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.mirone_fig,'dependentFigs',plugedWin);

	% Update handles structure
	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = handles.output;  end

% -----------------------------------------------------------------------------------------
function edit_StartYear_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < handles.StartYear),	set(hObject,'String',num2str(handles.StartYear)),	end

% -----------------------------------------------------------------------------------------
function edit_StartMonth_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1),	set(hObject,'String',num2str(handles.StartMonth)),	end

% -----------------------------------------------------------------------------------------
function edit_StartDay_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1),	set(hObject,'String',num2str(handles.StartDay)),	end

% -----------------------------------------------------------------------------------------
function edit_EndYear_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx > handles.EndYear),	set(hObject,'String',num2str(handles.EndYear)),		end

% -----------------------------------------------------------------------------------------
function edit_EndMonth_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx > 12),	set(hObject,'String',num2str(handles.EndMonth)),	end

% -----------------------------------------------------------------------------------------
function edit_EndDay_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx > 31),		set(hObject,'String',num2str(handles.EndDay)),	end

% -----------------------------------------------------------------------------------------
function edit_timeStep_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0),	set(hObject,'String','2'),	end

% -----------------------------------------------------------------------------------------
function edit_lifeSpan_CB(hObject, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx) || xx < str2double(get(handles.edit_timeStep,'String')))
    set(hObject,'String','10')
end

% -----------------------------------------------------------------------------------------
function edit_frameDelay_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0),	set(hObject,'String','0.25'),	end

% -----------------------------------------------------------------------------------------
function push_run_CB(hObject, handles)

	StartYear = str2double(get(handles.edit_StartYear,'String'));
	EndYear = str2double(get(handles.edit_EndYear,'String'));
	StartMonth = str2double(get(handles.edit_StartMonth,'String'));
	EndMonth = str2double(get(handles.edit_EndMonth,'String'));
	StartDay = str2double(get(handles.edit_StartDay,'String'));
	EndDay = str2double(get(handles.edit_EndDay,'String'));
	dt = str2double(get(handles.edit_timeStep,'String')) / (365 + isleapyear(StartYear));
	frame_int = str2double(get(handles.edit_frameDelay,'String'));
	vida = str2double(get(handles.edit_lifeSpan,'String'));
	stoped = 0;

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
        tit = [handles.fig_name '  ' sprintf('%02d/%02d/%d', day, month, year)];
        ind = find(handles.events_time <= t_up & handles.events_time >= t_low);
        set(handles.h_events,'XData',x(ind),'YData',y(ind));
        pause(frame_int)
		if ( ~get(handles.push_stop,'Userdata') )
			set(handles.figure1,'Name',tit)
		else
			stoped = 1;
			set(handles.push_stop,'Userdata',0)		% Reset for the next time
			break
		end
	end
	set(handles.h_events,'XData',handles.x_bak,'YData',handles.y_bak);
	set(handles.mirone_fig,'Renderer',dbf)
	if (~stoped)
		set(handles.figure1,'Name',handles.fig_name)
	end

% -----------------------------------------------------------------------------------------
function push_stop_CB(hObject, handles)
	set(hObject,'UserData', 1)

%--------------------------------------------------------------------------
function slider_time_CB(hObject, handles)
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
function push_help_CB(hObject, handles)

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
% --- Creates and returns a handle to the GUI figure. 
function animate_seismicity_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Seismicity Swing',...
'NumberTitle','off',...
'Position',[520 630 401 140],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_StartYear_CB'},...
'Position',[44 115 47 21],...
'Style','edit',...
'Tag','edit_StartYear');

uicontrol('Parent',h1,'Position',[4 105 36 33],'String',{'Start'; 'year'},'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_StartMonth_CB'},...
'Position',[160 115 47 21],...
'Style','edit',...
'Tag','edit_StartMonth');

uicontrol('Parent',h1,'Position',[116 109 41 30],'String',{'Start'; 'month'},'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_StartDay_CB'},...
'Position',[267 115 47 21],...
'Style','edit',...
'Tag','edit_StartDay');

uicontrol('Parent',h1,'Position',[231 109 36 30],'String',{'Start'; 'day'},'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_EndYear_CB'},...
'Position',[44 81 47 21],...
'Style','edit',...
'Tag','edit_EndYear');

uicontrol('Parent',h1,'Position',[4 71 36 30],'String',{'End'; 'year'},'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_EndMonth_CB'},...
'Position',[160 81 47 21],...
'Style','edit',...
'Tag','edit_EndMonth');

uicontrol('Parent',h1,'Position',[116 72 41 30],'String',{'End'; 'month'},'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_EndDay_CB'},...
'Position',[267 81 47 21],...
'Style','edit',...
'Tag','edit_EndDay');

uicontrol('Parent',h1,'Position',[237 76 27 30],'String',{'End'; 'day'},'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_timeStep_CB'},...
'Position',[44 25 47 21],...
'String','2',...
'Style','edit',...
'TooltipString','Pick events at this time rate',...
'Tag','edit_timeStep');

uicontrol('Parent',h1,'Position',[24 50 91 15],'String','Time step (days)','Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_lifeSpan_CB'},...
'Position',[160 25 47 21],...
'String','10',...
'Style','edit',...
'TooltipString','Events will stay on screen for these number of days',...
'Tag','edit_lifeSpan');

uicontrol('Parent',h1,'Position',[137 50 91 15],'String','Life span (days)','Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_frameDelay_CB'},...
'Position',[267 25 47 21],...
'String','0.25',...
'Style','edit',...
'TooltipString','Next event(s) will show after this interval (seconds)',...
'Tag','edit_frameDelay');

uicontrol('Parent',h1,'Position',[249 50 79 15],'String','Frame int (sec)','Style','text');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_stop_CB'},...
'Position',[331 6 66 23],...
'String','Stop',...
'UserData', 0, ...
'Tag','push_stop');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_run_CB'},...
'Position',[330 94 66 23],...
'String','Run',...
'Tag','push_run');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_help_CB'},...
'FontSize',10,...
'FontWeight','bold',...
'Position',[349 51 23 23],...
'String','?',...
'Tag','push_help');

uicontrol('Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Call',{@main_uiCB,h1,'slider_time_CB'},...
'Position',[10 3 310 15],...
'Style','slider',...
'TooltipString','Show events in the interval [time-Life span; time]',...
'Tag','slider_time');

function main_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

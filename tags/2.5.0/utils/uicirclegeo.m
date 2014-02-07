function [out1,out2] = uicirclegeo(varargin)
% UICIRCLEGEO:  Display of interactive small circle defined via mouse clicks.
%
%  UICIRCLEGEO lets the user click and move the cursor to define the center and
%  perimeter of the small circle. Double clicking on the circle displays two
%  control points (center point and radial) and the associated control window.
%  To translate the circle, click and drag the center control. To resize the
%  circle, click and drag the radial (o) control.  The control buttons can be
%  hidden by double clicking on the circle. The control window allows modification
%  of certain parameters of the small circle. The center point (latitude and longitude),
%  radius, and number of points can be changed using the appropriate edit boxes.
%  The distance units (kilometers, miles, nautical miles, radians) can be modified
%  using the appropriate pop-up menus.
%
%  UICIRCLEGEO(clon, clat, hAxes)	hAxes is optional. If not given defaults to gca
%  h = UICIRCLEGEO(...) returns the handles of the circles drawn.
%
%  [lat,lon] = UICIRCLEGEO(...) returns the latitude and longitude
%  matrices corresponding to the circle drawn.

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

	know_center = 0;

	if nargin == 1 && ischar(varargin{1})         % GUI callbacks
		circleui(varargin{:});
		return
	elseif (nargin >= 2 && isnumeric(varargin{1}) && isnumeric(varargin{2}))  % Circle center in input
		know_center = 1;
	end

	if ( nargin == 3 && ishandle(varargin{3}) && strcmp(get(varargin{3},'type'), 'axes') )
		s.h_axes = varargin{3};
		s.h_fig = get(s.h_axes, 'Parent');
	else
		s.h_fig = gcf;		s.h_axes = gca;
	end

	state = uisuspend_j(s.h_fig);   % Remember initial figure state

	%  Plot the circles
	set(s.h_fig,'Pointer', 'crosshair');
	w = waitforbuttonpress;
	if w == 0			% A mouse click
		if ~know_center
			h_circ = circFirstButtonDown;
		else			% Center was transmited in input
			h_circ = circFirstButtonDown(varargin{1},varargin{2});
		end
	else
		set(s.h_fig,'Pointer', 'arrow');
		out1 = [];    return
	end

	% I have to use a waitfor to give time to operate (that is, select a circle) inside the
	% the functions called by circFirstButtonDown. I don't understand why, bu if I don't
	% do this, we are just kicked out of the program whithout having time to do anything.
	try
		waitfor(h_circ, 'Tag', 'Completed');
		uirestore_j(state, 'nochildren');         % Restore the figure's initial state
	end

	lon_lat_rad = getappdata(h_circ,'LonLatRad');

	%  Get the circle definition.
	lat1 = lon_lat_rad(2);      lon1 = lon_lat_rad(1);
	xx = getappdata(h_circ,'FirstEndPoint');        % Retrieve the circle's first end point
	lon2 = xx(1);       lat2 = xx(2);
	rad = geo2dist([lon1 lon2],[lat1 lat2],'deg');
	az = azimuth_geo(lat1,lon1,lat2,lon2);
	rmappdata(h_circ,'FirstEndPoint')

	%  Set the correct structures
	s.npts = 180;       % Use 180 points by default
	s.clat = lat1;      s.clon = lon1;
	s.rlat = lat2;      s.rlon = lon2;
	s.rad = rad;        s.az = az;
	s.controls = 'off';
	s.hcontrol = [];    s.hcenter = [];
	s.hcirc = h_circ;   s.hend = [];
	s.parent = get(s.hcirc,'parent');	
	s.num = rand;
	s.omega = [];       % Used when the circle is about an Euler pole
	s.plates = [];      %                   ''
	set(h_circ,'userdata',s,'buttondownfcn','uicirclegeo(''circlemousedown'')');
	set(h_circ,'Tag','circleGeo')

	%  Set output arguments
	if nargout == 1
		out1 = h_circ;
	elseif nargout == 2
		out1 = get(h_circ,'XData');
		out2 = get(h_circ,'YData');
	end

%---------------------------------------------------------------------------------------------------
function h = circFirstButtonDown(lon,lat)
% Draw the circle
    hFig = get(0,'CurrentFigure');		handles = guidata(hFig);
    lc = handles.DefLineColor;			lt = handles.DefLineThick;
    
	h = line('XData', [], 'YData', [], 'Parent',handles.axes1, 'Tag','','Color',lc,'LineWidth',lt);
	if (nargin == 2),		pt(1,1) = lon;		pt(1,2) = lat;
	else					pt = get(handles.axes1, 'CurrentPoint');
	end
	x_lim = get(handles.axes1,'xlim');    y_lim = get(handles.axes1,'ylim');
	set(hFig,'WindowButtonMotionFcn',{@wbm_circle,handles.axes1,[pt(1,1) pt(1,2)],h,[x_lim y_lim], handles.geog},...
        'WindowButtonDownFcn',{@wbd_circle,hFig,handles.axes1,h});

% ----------------------------------------------
function wbm_circle(obj, eventdata, hAxes, center, h, lim, lon_type)
	pt = get(hAxes, 'CurrentPoint');
	if (pt(1,1) < lim(1)) || (pt(1,1) > lim(2)) || (pt(1,2) < lim(3)) || (pt(1,2) > lim(4)),	return,	end
	rad = geo2dist([pt(1,1) center(1)],[pt(1,2) center(2)],'deg');
	[latc,lonc] = circ_geo(center(2), center(1), rad, [], 180, lon_type == 2);
	% Find the eventual Date line discontinuity and insert a NaN on it
	ind = find(abs(diff(lonc)) > 100);   % 100 is good enough
	if (~isempty(ind))
        if (length(ind) == 2)
            lonc = [lonc(1:ind(1)) NaN lonc(ind(1)+1:ind(2)) NaN lonc(ind(2)+1:end)];
            latc = [latc(1:ind(1)) NaN latc(ind(1)+1:ind(2)) NaN latc(ind(2)+1:end)];
        elseif (length(ind) == 1)
            lonc = [lonc(1:ind) NaN lonc(ind+1:end)];   latc = [latc(1:ind) NaN latc(ind+1:end)];
        end
	end
	set(h, 'XData', lonc, 'YData', latc);
	setappdata(h,'LonLatRad',[center(1) center(2) rad])

% ----------------------------------------------
function wbd_circle(obj,eventdata,hFig,hAxes,h)
	pt = get(hAxes, 'CurrentPoint');
	setappdata(h,'FirstEndPoint',[pt(1,1) pt(1,2)])		% Save the circle's first end point
	set(hFig,'WindowButtonMotionFcn','', 'WindowButtonDownFcn','', 'Pointer', 'arrow')
	set(h,'Tag','Completed')			% Signal waitfor that we are done

% -----------------------------------------------------------------------------------------
function move_circle(obj,eventdata,h)
% Translate the circle
	s = get(h,'userdata');		handles = guidata(s.h_fig);
	state = uisuspend_j(s.h_fig);   % Remember initial figure state
	center = getappdata(h,'LonLatRad');
	hcenter = s.hcenter;		hend = s.hend;
	x_lim = get(s.h_axes,'xlim');      y_lim = get(s.h_axes,'ylim');
	set(s.h_fig,'WindowButtonMotionFcn',{@wbm_MoveCircle,h,center,hcenter,hend,s,[x_lim y_lim], handles.geog},...
        'WindowButtonUpFcn',{@wbu_MoveCircle,h,hcenter,hend,state}, 'Pointer', 'crosshair');

% ---------
function wbm_MoveCircle(obj,eventdata,h,center,hcenter,hend,s,lim, lon_type)
	pt = get(s.h_axes, 'CurrentPoint');
	if (pt(1,1)<lim(1)) || (pt(1,1)>lim(2)) || (pt(1,2)<lim(3)) || (pt(1,2)>lim(4)),	return,		end
	[latc,lonc] = circ_geo(pt(1,2), pt(1,1), center(3), [], s.npts, lon_type == 2);
	% Find the eventual date line discontinuity and insert a NaN on it
	ind = find(abs(diff(lonc)) > 100);   % 100 is good enough
	if (~isempty(ind))
        if (length(ind) == 2)
            lonc = [lonc(1:ind(1)) NaN lonc(ind(1)+1:ind(2)) NaN lonc(ind(2)+1:end)];
            latc = [latc(1:ind(1)) NaN latc(ind(1)+1:ind(2)) NaN latc(ind(2)+1:end)];
        elseif (length(ind) == 1)
            lonc = [lonc(1:ind) NaN lonc(ind+1:end)];   latc = [latc(1:ind) NaN latc(ind+1:end)];
        end
	end
	set(h, 'XData', lonc, 'YData', latc);
	set(hcenter,'XData',pt(1,1), 'YData', pt(1,2))				% Move the circle's center marker
	[s.rlat,s.rlon] = circ_geo(pt(1,2),pt(1,1),s.rad,s.az,1, lon_type == 2);	% Compute new end point
	set(hend,'XData',s.rlon, 'YData', s.rlat)					% Move the circle's end marker
	
	% update data in controls window if it exists
	if ishandle(s.hcontrol)
        set(findobj(s.hcontrol,'tag','lat'),'string',num2str(pt(1,2)));
        set(findobj(s.hcontrol,'tag','lon'),'string',num2str(pt(1,1)));
	end
	set(h,'userdata',s);

% ---------
function wbu_MoveCircle(obj,eventdata,h,hcenter,hend,state)
	% check if x,y is inside of axis
	s = get(h,'userdata');
	pt = get(s.h_axes, 'CurrentPoint');
	x = pt(1,1);    y = pt(1,2);
	x_lim = get(s.h_axes,'xlim');      y_lim = get(s.h_axes,'ylim');
	if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2));   return;     end
	
	% Update the new circles center and end in userdata
	s = get(h,'userdata');
	s.clon = x;     s.clat = y;
	s.rlon = get(hend,'XData');     s.rlat = get(hend,'YData');
	set(h,'userdata',s);
	setappdata(h,'LonLatRad',[s.clon s.clat s.rad])
	if ishandle(s.hcontrol)     % That is, if user didn't kill the controls window
        set(s.hcontrol,'userdata',s)
	end
	uirestore_j(state, 'nochildren');         % Restore the figure's initial state

% -----------------------------------------------------------------------------------------
function resize_circle(obj,eventdata,h)
% Resize the circle
	s = get(h,'userdata');		handles = guidata(s.h_fig);
	state = uisuspend_j(s.h_fig);		% Remember initial figure state
	x_lim = get(s.h_axes,'xlim');		y_lim = get(s.h_axes,'ylim');
	center = [s.clon s.clat];
	hcenter = s.hcenter;				hend = s.hend;
	set(s.h_fig,'WindowButtonMotionFcn',{@wbm_ResizeCircle,h,center,hend,[x_lim y_lim], handles.geog},...
		'WindowButtonUpFcn',{@wbu_ResizeCircle,h,hcenter,hend,state}, 'Pointer', 'crosshair');

% ---------
function wbm_ResizeCircle(obj,eventdata,h,center,hend,lim, lon_type)
	s = get(h,'userdata');
	pt = get(s.h_axes, 'CurrentPoint');
	if (pt(1,1) < lim(1)) || (pt(1,1) > lim(2)) || (pt(1,2) < lim(3)) || (pt(1,2) > lim(4)),	return,	end
	rad = geo2dist([pt(1,1) center(1)],[pt(1,2) center(2)],'deg');
	s.rad = rad;                                            % The radius is stored in degrees
	s.rlon = pt(1,1);
	s.rlat = pt(1,2);
	[latc,lonc] = circ_geo(center(2),center(1),rad,[],s.npts, lon_type == 2);
	% Find the eventual date line discontinuity and insert a NaN on it
	ind = find(abs(diff(lonc)) > 100);   % 100 is good enough
	if (~isempty(ind))
		if (length(ind) == 2)
			lonc = [lonc(1:ind(1)) NaN lonc(ind(1)+1:ind(2)) NaN lonc(ind(2)+1:end)];
			latc = [latc(1:ind(1)) NaN latc(ind(1)+1:ind(2)) NaN latc(ind(2)+1:end)];
		elseif (length(ind) == 1)
			lonc = [lonc(1:ind) NaN lonc(ind+1:end)];   latc = [latc(1:ind) NaN latc(ind+1:end)];
		end
	end
	set(h, 'XData', lonc, 'YData', latc);
	set(hend,'XData',pt(1,1), 'YData', pt(1,2))             % Move the end marker together with resing circle

	% update data in controls window if it exists
	if ishandle(s.hcontrol)
		objtype = get(s.hcontrol,'tag');
		switch objtype
			case 'sccontrol'
				obja = findobj(s.hcontrol,'style','popupmenu','tag','units');
				str = get(obja, 'String');
				units = str{get(obja, 'Val')};
				switch units
					case 'Kilometers',		rad1 = geo2dist(rad,'deg','km');
					case 'Miles',			rad1 = geo2dist(rad,'deg','mi');
					case 'Nautical Miles',	rad1 = geo2dist(rad,'deg','nm');
					case 'Radians',			rad1 = rad*pi/180;
				end
				set(findobj(s.hcontrol,'tag','rad'),'string',num2str(rad1));
		end
		set(s.hcontrol,'userdata',s)
	end
	set(h,'userdata',s);

% ---------
function wbu_ResizeCircle(obj,eventdata,h,hcenter,hend,state)
% Update the azimuth between center and end poin
	s = get(h,'userdata');
	s.az = azimuth_geo(s.clat,s.clon,s.rlat,s.rlon);
	set(h,'userdata',s);
	setappdata(h,'LonLatRad',[s.clon s.clat s.rad])
	uirestore_j(state, 'nochildren')		% Restore the figure's initial state

%---------------------------------------------------------------------------------------------------
function circleui(action)

this_obj = gcbo;		this_gco = gco;
switch action
	case 'createcontrols'			% create GUI controls
		s = get(this_gco,'userdata');
		h = figure('units','char','pos',[20 5 36.4 16.8],...
			'numbertitle','off','name','Small Circles','tag','sccontrol',...
			'resize','off','HandleVisibility','Call','Menubar','none');
		framecolor = [0.8944 0.8944 0.8944];	
		uicontrol('style','frame','units','char', 'pos',[1.1866 10.9113 34.0 5.2375], 'backgroundcolor',framecolor);
		uicontrol('style','frame','units','char', 'pos',[1.1866 3.5 34.0 7.1], 'backgroundcolor',framecolor);
		uicontrol('style','text','string','Circle center', 'units','char','pos',[2.4 14.6154 18.4 1.3087],...
						'fontweight','bold','horizontalalignment','left',...
						'backgroundcolor',framecolor);
		uicontrol('style','text','string','Size', 'units','char','pos',[2.4 9.1538 18.4 1.3087],...
						'fontweight','bold','horizontalalignment','left',...
						'backgroundcolor',framecolor);
		uicontrol('style','text','string','Lat', 'units','char','pos',[2.4 13.0769 18.4 1.3087],...
						'fontweight','bold','horizontalalignment','left',...
						'backgroundcolor',framecolor);
		uicontrol('style','text','string','Lon', 'units','char','pos',[2.4 11.3846 18.4 1.3087],...
						'fontweight','bold','horizontalalignment','left',...
						'backgroundcolor',framecolor);
		uicontrol('style','text','string','Units', 'units','char','pos',[2.4 7.3 18.4 1.3087],...
						'fontweight','bold','horizontalalignment','left',...
						'backgroundcolor',framecolor);
		uicontrol('style','text','string','Radius', 'units','char','pos',[2.4 5.6154 18.4 1.3087],...
						'fontweight','bold','horizontalalignment','left',...
						'backgroundcolor',framecolor);
		uicontrol('style','text','string','Npts', 'units','char','pos',[2.4 4.0769 18.4 1.3087],...
						'fontweight','bold','horizontalalignment','left',...
						'backgroundcolor',framecolor);
		uicontrol('style','edit','units','char','pos', [11.8 13.0769 21.6016 1.4615],...
						'fontsize',9,'fontweight','bold','tag','lat','string','1',...
						'Call','uicirclegeo(''changelat'')');
		uicontrol('style','edit','units','char','pos', [11.8 11.3846 21.6016 1.4615],...
						'fontsize',9,'fontweight','bold','tag','lon','string','2',...
						'Call','uicirclegeo(''changelon'')');
		uicontrol('style','edit','units','char','pos', [11.8 5.6923 21.6016 1.4615],...
						'fontsize',9,'fontweight','bold','tag','rad','string','3',...
						'Call','uicirclegeo(''changeradius'')');
		uicontrol('style','edit','units','char','pos', [11.8 3.9231 21.6016 1.4615],...
						'fontsize',9,'fontweight','bold','tag','npts','string','4',...
						'Call','uicirclegeo(''npts'')');
		uicontrol('style','push','units','char','pos', [11.4 0.3077 15.2 1.3077],...
						'string','Close','fontweight','bold',...
						'Call','uicirclegeo(''close'')');
		popstr = {'Kilometers','Miles','Nautical Miles','Radians'};				   
		uicontrol('style','popup','units','char','pos', [11.8 7.6 21.6016 1.4615],...
						'string',popstr,'tag','units','fontsize',9,...
						'fontweight','bold',...
						'Call','uicirclegeo(''changeunits'')');
		% update userdata
		s.hcontrol = h;
		s.parent = get(s.hcirc,'parent');	
		set(findobj(h,'tag','lat'),'string',num2str(s.clat));
		set(findobj(h,'tag','lon'),'string',num2str(s.clon));
		set(findobj(h,'tag','rad'),'string',num2str(geo2dist(s.rad,'deg','km')));
		set(findobj(h,'tag','npts'),'string',num2str(s.npts));
		% update the values in the objects
		set(s.hcontrol,'userdata',s)
		set(s.hcirc,'userdata',s)
		set(s.hcenter,'userdata',s)
		set(s.hend,'userdata',s)
	case 'npts'								% change the number of points
		par = get(this_obj,'parent');		s = get(par,'userdata');
		% return original string if user gives stupid data
		npts = str2double(get(this_obj,'string'));
		if (isempty(npts) || npts < 4)
			set(this_obj,'string',num2str(s.npts)),		return
		end	
		s.npts = npts;
		[latc,lonc] = circ_geo(s.clat,s.clon,s.rad,[],npts);    % compute the new radial points
		set(s.hcirc, 'XData', lonc, 'YData', latc);
		% update the userdata
		if ishandle(s.hcontrol),    set(s.hcontrol,'userdata',s);	end
		% save the userdata	
		set(par,'userdata',s)
		set(s.hcirc,'userdata',s)
		set(s.hcenter,'userdata',s)
		set(s.hend,'userdata',s)
	case 'changelat'						% change the latitude
		par = get(this_obj,'parent');		s = get(par,'userdata');
		lat = str2double(get(this_obj,'string'));
		if isempty(lat)     % return original string if user said nonsense
			set(this_obj,'string',num2str(s.clat));     return
		end	
		y_lim = get(gca,'ylim');
		if (lat < y_lim(1) || lat > y_lim(2))       % If out of map, do nothing
			set(this_obj,'string',num2str(s.clat));     return
		end	
		s.clat = lat;
		s.az = azimuth_geo(s.clat,s.clon,s.rlat,s.rlon);
		% compute the new radial points
		[latc,lonc] = circ_geo(s.clat,s.clon,s.rad,[],s.npts);
		set(s.hcirc, 'XData', lonc, 'YData', latc);
		[s.rlat,s.rlon] = circ_geo(s.clat,s.clon,s.rad,s.az,1);     % Compute new end point
		set(s.hcenter,'XData',s.clon, 'YData', s.clat)              % Move the circle's center marker
		set(s.hend,'XData',s.rlon, 'YData', s.rlat)                 % Move the circle's end marker
		% update the userdata
		if ishandle(s.hcontrol),    set(s.hcontrol,'userdata',s);	end
		% save the userdata	
		set(par,'userdata',s)
		set(s.hcirc,'userdata',s)
		set(s.hcenter,'userdata',s)
		set(s.hend,'userdata',s)
		setappdata(s.hcirc,'LonLatRad',[s.clon s.clat s.rad])
	case 'changelon'						% change the longitude
		par = get(this_obj,'parent');		s = get(par,'userdata');
		lon = str2double(get(this_obj,'string'));
		% return original string if user entries wrong data
		if isempty(lon)
			set(this_obj,'string',num2str(s.clon)),			return
		end	
		x_lim = get(gca,'xlim');
		if (lon < x_lim(1) || lon > x_lim(2))			% If out of map, do nothing
			set(this_obj,'string',num2str(s.clat)),		return
		end	
		s.clon = lon;
		s.az = azimuth_geo(s.clat,s.clon,s.rlat,s.rlon);
		[latc,lonc] = circ_geo(s.clat,s.clon,s.rad,[],s.npts);
		set(s.hcirc, 'XData', lonc, 'YData', latc);
		[s.rlat,s.rlon] = circ_geo(s.clat,s.clon,s.rad,s.az,1);     % Compute new end point
		set(s.hcenter,'XData',s.clon, 'YData', s.clat)              % Move the circle's center marker
		set(s.hend,'XData',s.rlon, 'YData', s.rlat)                 % Move the circle's end marker
		% update the userdata
		if ishandle(s.hcontrol),    set(s.hcontrol,'userdata',s);	end
		set(par,'userdata',s)
		set(s.hcirc,'userdata',s)
		set(s.hcenter,'userdata',s)
		set(s.hend,'userdata',s)
		setappdata(s.hcirc,'LonLatRad',[s.clon s.clat s.rad])
	case 'changeradius'						% change the radius
		par = get(this_obj,'parent');		s = get(par,'userdata');
		obja = findobj(par,'tag','units');
		str = get(obja, 'String');
		units = str{get(obja, 'Val')};
		radius = str2double(get(this_obj,'string'));
		% return original string if user gives wrong data
		if (isempty(radius) || radius <= 0)
			switch units
				case 'Kilometers',		rad1 = geo2dist([s.clon s.rlon],[s.clat s.rlat],'km');
				case 'Miles',			rad1 = geo2dist([s.clon s.rlon],[s.clat s.rlat],'mi');
				case 'Nautical Miles',	rad1 = geo2dist([s.clon s.rlon],[s.clat s.rlat],'nm');
				case 'Radians',			rad1 = geo2dist([s.clon s.rlon],[s.clat s.rlat],'km') / 6371;
			end
			set(this_obj,'string',rad1)
			return
		end	
		% update data in control window
		switch units                        % All angles must be in degrees for use in circ_geo
			case 'Kilometers',		rad1 = dist2rad(radius,'km')*180/pi;
			case 'Miles',			rad1 = dist2rad(radius,'mi')*180/pi;
			case 'Nautical Miles',	rad1 = dist2rad(radius,'nm')*180/pi;
			case 'Radians',			rad1 = radius*180/pi;
		end
		s.rad = rad1;

		% compute the new radial points
		[latc,lonc] = circ_geo(s.clat,s.clon,rad1,[],s.npts);
		set(s.hcirc, 'XData', lonc, 'YData', latc);
		set(s.hcenter,'XData',s.clon, 'YData', s.clat)      % Move the circle's center marker
		az = azimuth_geo(s.clat,s.clon,s.rlat,s.rlon);      % Compute the azimuth between center and end
		[s.rlat,s.rlon] = circ_geo(s.clat,s.clon,rad1,az,1);% Compute new end point
		set(s.hend,'XData',s.rlon, 'YData', s.rlat)         % Move the circle's end marker

		% update the userdata
		if ishandle(s.hcontrol),    set(s.hcontrol,'userdata',s);	end
		% save the userdata	
		set(par,'userdata',s)
		set(s.hcirc,'userdata',s)
		set(s.hcenter,'userdata',s)
		set(s.hend,'userdata',s)
		setappdata(s.hcirc,'LonLatRad',[s.clon s.clat s.rad])
	case 'changeunits'              % change the radius units	
		par = get(this_obj,'parent');
		s = get(par,'userdata');
		str = get(this_obj, 'String');
		units = str{get(this_obj, 'Val')};
		switch units
			case 'Kilometers',		rad1 = geo2dist([s.clon s.rlon],[s.clat s.rlat],'km');
			case 'Miles',			rad1 = geo2dist([s.clon s.rlon],[s.clat s.rlat],'mi');
			case 'Nautical Miles',	rad1 = geo2dist([s.clon s.rlon],[s.clat s.rlat],'nm');
			case 'Radians',			rad1 = geo2dist([s.clon s.rlon],[s.clat s.rlat],'km') / 6371;
		end
		set(findobj(par,'tag','rad'),'string',num2str(rad1));
	case 'close'                    % close the control window
		par = get(this_obj,'parent');
		s = get(par,'userdata');
		% delete the objects
		if ishandle(s.hcenter),		delete(s.hcenter),	end
		if ishandle(s.hend),		delete(s.hend),		end
		if ishandle(s.hcontrol)
			objtype = get(s.hcontrol,'tag');
			switch objtype
				case 'sccontrol'	
					delete(s.hcontrol); 
			end		
		end
		s.controls = 'off';
		s.hcontrol = [];
		s.hcenter = [];
		s.hend = [];
		if ishandle(s.hcirc),	set(s.hcirc,'userdata',s),	end	
	case 'circlemousedown'              % mouse down on circle
		stype = get(get(get(this_gco,'Parent'),'Parent'),'selectiontype');	% Same as get(gcf,'selectiontype')
		s = get(this_gco,'userdata');
		s.parent = get(this_gco,'parent');	
		switch stype		% shift-click to toggle control points
			case 'open'
					switch s.controls
						case 'on'
							% kill the control window if its open
							if ishandle(s.hcontrol)
								objtype = get(s.hcontrol,'tag');
								if strcmp(objtype,'sccontrol')
									delete(s.hcontrol);     s.hcontrol = [];
								end
							end
							if (ishandle(s.hcenter)),	delete(s.hcenter);	end		% delete center
							s.hcenter = [];
							if (ishandle(s.hend)),		delete(s.hend);		end		% delete end
							s.hend = [];
							s.controls = 'off';
							set(this_gco,'userdata',s)
						case 'off'
							s.controls = 'on';
							hcirc = this_gco;
							hcent = line(s.clon,s.clat,'Marker','p','markerfacecolor','r','MarkerEdgeColor','k',...
									'MarkerSize',8,'userdata',s,'buttondownfcn',{@move_circle,hcirc});
							hend = line(s.rlon,s.rlat,'Marker','o','markerfacecolor','r','markerfacecolor','r',...
									'MarkerEdgeColor','k','userdata',s,'buttondownfcn',{@resize_circle,hcirc});
							s.hcenter = hcent;
							s.hend = hend;
							set(this_gco,'userdata',s)
							% display the control window
							uicirclegeo('createcontrols')
							% turn off handle visibility
							hcontrol = findobj('tag','sccontrol');
							for i = 1:length(hcontrol)
								r = get(hcontrol(i),'userdata');
								if r.clat == s.clat && r.num == s.num
									s.hcontrol = hcontrol(i);
									set(hcontrol(i),'handlevisibility','off')
								end	
							end	
							set(hcirc,'userdata',s)
							set(s.hcenter,'userdata',s)
							set(s.hend,'userdata',s)
					end
		end
end

%-------------------------------------------------------------------------------------------
function dist = geo2dist(lon,lat,opt)
% This is a quick (no error testing) computation of the spherical distances.
% dist = geo2dist(lon,lat,opt) where LON and LAT contain a two point vector in geog coordinates
% dist = geo2dist(ang,units,opt) where ANG is a angular distance, UNITS is either 'deg' or 'rad'
% In all cases OPT must contain one othe following: NM, KM, MI or DEG

	D2R = pi/180;
	ang = [];   units = [];

	if (nargin == 3 && ischar(lat))		% Input is an angle and not two points
		ang = lon;     units = lat;
	end

	if isempty(ang)						% Two points as input
		lon_i = lon(1)*D2R;		lon_f = lon(2)*D2R;
		lat_i = lat(1)*D2R;		lat_f = lat(2)*D2R;
		c = sin(lat_i)*sin(lat_f) + cos(lat_i)*cos(lat_f)*cos(lon_f-lon_i);
		ang = acos(c);
	else
		if strcmpi(units(1),'d')		% angle was given in degrees
			ang = ang * D2R;
		end
	end
	switch lower(opt)
		case 'nm',		dist = ang * 6371 / 1.852;			% Distance in Nautical Miles
		case 'mi',		dist = ang * 6371 / 1.6093;			% Distance in Miles
		case 'km',		dist = ang * 6371;					% Distance in km
		case 'deg',		dist = ang / D2R;					% Arc distance in degrees
	end

%-------------------------------------------------------------------------------------------
function teta = dist2rad(dist,opt)
% Convert arc distances DIST to angle. Again spherical approximation.
	switch lower(opt)
		case 'nm',		teta = dist * 1.852 / 6371;
		case 'mi',		teta = dist * 1.6093 / 6371;
		case 'km',		teta = dist / 6371;
	end

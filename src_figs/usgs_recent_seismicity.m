function usgs_recent_seismicity(varargin)
% Helper window to select and download a recent seismicity CVS file from USGS
%
% It slices by magnitude (intervals of 1) but not depth.
% If no base map, create one that holds all the events. Otherwise throw away
% the events outside the map limits.
% The events time is not being stored yet.

%	Copyright (c) 2004-2014 by J. Luis
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

% $Id: usgs_recent_seismicity.m 4440 2014-05-01 18:39:30Z j $

	handMir = guidata(varargin{1});		% Not tested but must be a Mirone fig handle

	hFig = figure('Pos',[520 541 401 160], 'MenuBar','none',...
		'Color',get(0,'factoryUicontrolBackgroundColor'),...
		'Name','Plot USGS recent seismicity',...
		'NumberTitle','off',...
		'Resize','off',...
		'HandleVisibility','callback', 'Vis','off');
	move2side(hFig,'right')

	handles.path_data = handMir.path_data;
	handles.path_tmp = handMir.path_tmp;
	handles.no_file  = handMir.no_file;
	handles.geog     = handMir.geog;		% If calling Fig is not empty and not geog we need to create a new Mir fig
	handles.hMirFig  = handMir.figure1;
	handles.hMirAx   = handMir.axes1;	
	handles.figure1  = hFig;

	uicontrol('Parent',hFig, 'Pos',[10  133 90 15],'Style','text','Str','Past Hour','FontSize',10,'HorizontalAlignment','left')
	uicontrol('Parent',hFig, 'Pos',[110 133 90 15],'Style','text','Str','Past Day','FontSize',10,'HorizontalAlignment','left')
	uicontrol('Parent',hFig, 'Pos',[210 133 90 15],'Style','text','Str','Past 7 Days','FontSize',10,'HorizontalAlignment','left')
	uicontrol('Parent',hFig, 'Pos',[310 133 90 15],'Style','text','Str','Past 30 Day','FontSize',10,'HorizontalAlignment','left')

	names = {'All', 'M2.5+', 'M4.5+', 'Significant'};
	handles.hRadios = zeros(4,4);
	for (m = 1:4)
		for (n = 1:4)
			handles.hRadios(m,n) = uicontrol('Parent',hFig, 'Pos',[(10 + (n-1)*100) (40 + (m-1)*22) 80 22], ...
				'Call',@radios_CB,'Str', names{m}, 'Style','radiobutton', 'UserData',[m n]);
		end
	end
	set(handles.hRadios(1,2),'Val',1)		% Make the "Past Day (All)" the default one

% 	uicontrol('Parent',hFig, 'Pos',[210 9 82 23], 'Str','Plot control', ...
% 		'Tooltip','Call a new window to allow fine plot control','Style','checkbox')
	uicontrol('Parent',hFig, 'Pos',[323 9 69 21], 'FontName','Helvetica',...
		'Call',@push_OK_CB, 'FontWeight','bold','String','OK')

	guidata(hFig, handles)
	set(hFig,'Vis','on')

% --------------------------------------------------------------------------------------
function push_OK_CB(obj, evt)
% Find the selected Radio, download the file and plot it
	handles = guidata(obj);
	break_too = false;
	for (m = 1:4)			% Loop over magnitude classification
		for (n = 1:4)		% Loop over period
			if (get(handles.hRadios(m,n),'Val'))	% Found the selected one
				break_too = true;
				break
			end
		end
		if (break_too),		break,	end				% Dirty trick because of stupid paranoia against GOTOs
	end

	% Create the file name
	periods = {'hour', 'day', 'week', 'month'};
	mags = {'all_', '2.5_', '4.5_', 'significant_'};
	name = [mags{m} periods{n}];

	dest_fiche = [handles.path_tmp 'usgs.csv'];
	url = ['http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/' name '.csv'];
	if (strncmp(computer,'PC',2))
		dos(['wget "' url '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
	else
		unix(['wget "' url '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
	end
	r = csv2cell(dest_fiche);
	if (isempty(r))
		errordlg(['Some network problem occured and the ' url ' was not downloaded'], 'Error'),		return
	end
	lon = str2double(r(2:end,3));	lat = str2double(r(2:end,2));
	mag = str2double(r(2:end,5));	dep = str2double(r(2:end,4));

	if (handles.no_file || ~handles.geog)
		w = max(-180, round(min(lon)-5));	e = min(180, round(max(lon)+5));	% With a padding of 5 degs but not Earthflow
		s = max(-90,  round(min(lat)-5));	n = min(90,  round(max(lat)+5));
		pix_x = round(getPixel_coords(5400, [-180 180], [w e]));			% Now convert to the correct indices of big image 
		pix_y = 2700 - round(getPixel_coords(2700, [-90 90], [s n])) + 1;	% Again the Y origin shit
		pix_y = pix_y(2:-1:1);
		opt_r = sprintf('-r%d/%d/%d/%d', pix_x(1:2), pix_y(1:2));
		img = gdalread([handles.path_data 'etopo4.jpg'], opt_r, '-U');

		tmp.X = [w e];      tmp.Y = [s n];		tmp.geog = 1;
		x_inc = diff(tmp.X) / (size(img,2)-1);
		y_inc = diff(tmp.Y) / (size(img,1)-1);
		tmp.head = [tmp.X tmp.Y 0 255 0 x_inc y_inc];
		tmp.name = ['USGS seismicity ' name];

		h = mirone(img, tmp);				% Create a new Mir fig
		if (ishandle(handles.hMirFig) && handles.no_file)		% Delete the old, empty, Mir fig
			delete(handles.hMirFig)
		end
		handMir = guidata(h);
		handles.hMirFig = handMir.figure1;
		handles.hMirAx  = handMir.axes1;
		guidata(handles.figure1, handles)	% Update the handles 
	else
		mirFigLimits = getappdata(handles.hMirAx,'ThisImageLims');
		if (handles.geog == 2),		lon = lon + 180;		end
		[lon, lat, indx, indy] = ...	% Get rid of points that are outside the map limits
			aux_funs('in_map_region',handles,lon,lat,0.01,mirFigLimits);
		if (isempty(lon))
			warndlg('No events inside map region', 'Warning')
			return
		end
		if (~isempty(indx) || ~isempty(indy))		% If needed, clip outside map data
			mag(indx) = [];		mag(indy) = [];
			dep(indx) = [];		dep(indy) = [];
		end
	end

	[mag_s, ind_mag] = slice_events(mag, dep);	% 9 Mag slices and 6 Depth slices
	ss = [2 3 4 5 6 7 8 10 12];				% Symbol sizes

	for (k = 1:numel(ss))
		if (isempty(ind_mag{k})),		continue,	end
		h = line('Parent',handles.hMirAx, 'XData',lon(ind_mag{k}), 'YData',lat(ind_mag{k}), 'LineStyle','none', 'Marker','o', ...
			'MarkerSize',ss(k), 'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'Tag','Earthquakes');
		setappdata(h,'SeismicityDepth',uint16(dep(ind_mag{k}) * 10));	% Save events depth
		setappdata(h,'SeismicityMag',  uint8 (mag(ind_mag{k}) * 10));	% Save events magnitude
		%setappdata(h,'SeismicityTime',year_dec);		% Save events time
		draw_funs(h,'Earthquakes',[])
	end

% --------------------------------------------------------------------------------------
function radios_CB(obj, evt)
% Make sure that only current selection is checked
	handles = guidata(obj);
	if (~get(obj,'Val')),	set(obj,'Val',1),	return,		end
	for (k = 1:numel(handles.hRadios))
		set(handles.hRadios(k),'Val',0)		% Easier to set them all to 0 and next put current one to 1
	end
	set(obj,'Val',1)

% -------------------------------------------------------------------------------------
function pix_coords = getPixel_coords(img_length, XData, axes_coord)
% Convert coordinates from axes (real coords) to image (pixel) coordinates.
% IMG_LENGTH is the image width (n_columns)
% XDATA is the image's [x_min x_max] in axes coordinates
% AXES_COORD is the (x,y) coordinate of the point(s) to be converted

	slope = (img_length - 1) / (XData(end) - XData(1));
	if ((XData(1) == 1) && (slope == 1))
		pix_coords = axes_coord;
	else
		pix_coords = slope * (axes_coord - XData(1)) + 1;
	end

% -------------------------------------------------------------------------------------
function  [mag_s, ind_mag, dep_s, ind_dep] = slice_events(mag, dep)
% Slice the events by magnitude and depth
	% Slice by magnitudes of 1 centered at integer mags
	ind_mag = cell(9,1);	mag_s = cell(9,1);
	m0 = [0   1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5];
	m1 = [1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5];
	for (k = 1:9)
		ind_mag{k} = find(mag > m0(k) & mag <= m1(k));
		if (~isempty(ind_mag{k})),		mag_s{k} = mag(ind_mag{k});		end
	end
	
	% Now, if requested, slice by depths
	if (nargout > 2)
		ind_dep = cell(6,1);	dep_s = cell(6,1);
		d0 = [-1 33  70 150 300  500];
		d1 = [33 70 150 300 500 1000];
		for (k = 1:6)
			ind_dep{k} = find(dep > d0(k) & dep <= d1(k));
			if (~isempty(ind_dep{k})),		dep_s{k} = dep(ind_dep{k});		end
		end
	end
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

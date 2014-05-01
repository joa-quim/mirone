function varargout = usgs_recent_seismicity(varargin)
% Helper window to select and download a recent seismicity CVS file from USGS

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

% $Id$

	handMir = guidata(varargin{1});		% Not tested but must be a Mirone fig handle

	hFig = figure('Pos',[520 541 401 180], 'MenuBar','none',...
		'Color',get(0,'factoryUicontrolBackgroundColor'),...
		'Name','Plot USGS recent seismicity',...
		'NumberTitle','off',...
		'Resize','off',...
		'HandleVisibility','callback', 'Vis','off');
	move2side(hFig,'right')

	handles.hRadios = zeros(5,4);
	handles.figure1 = hFig;
	handles.path_tmp = handMir.path_tmp;
	handles.hMirFig  = handMir.figure1;
	handles.hMirAx   = handMir.axes1;	

	uicontrol('Parent',hFig, 'Pos',[10  155 90 15],'Style','text','Str','Past Hour','FontSize',10,'HorizontalAlignment','left')
	uicontrol('Parent',hFig, 'Pos',[110 155 90 15],'Style','text','Str','Past Day','FontSize',10,'HorizontalAlignment','left')
	uicontrol('Parent',hFig, 'Pos',[210 155 90 15],'Style','text','Str','Past 7 Days','FontSize',10,'HorizontalAlignment','left')
	uicontrol('Parent',hFig, 'Pos',[310 155 90 15],'Style','text','Str','Past 30 Day','FontSize',10,'HorizontalAlignment','left')

	names = {'All', 'M1.0+', 'M2.5+', 'M4.5+', 'Significant'};
	for (m = 1:5)
		for (n = 1:4)
			handles.hRadios(m,n) = uicontrol('Parent',hFig, 'Pos',[(10 + (n-1)*100) (40 + (m-1)*22) 80 22], ...
				'Call',@radios_CB,'Str', names{m}, 'Style','radiobutton', 'UserData',[m n]);
		end
	end

	uicontrol('Parent',hFig, 'Pos',[227 9 82 23], 'Str','Plot control', ...
		'Tooltip','Call a new window to allow fine plot control','Style','checkbox')
	uicontrol('Parent',hFig, 'Pos',[323 9 69 21], 'FontName','Helvetica',...
		'Call',@push_OK_CB, 'FontWeight','bold','String','OK')

	guidata(hFig, handles)
	set(hFig,'Vis','on')

% --------------------------------------------------------------------------------------
function push_OK_CB(obj, evt)
% Find the selected Radio, download the file and plot it
	handles = guidata(obj);
	break_too = false;
	for (m = 1:5)
		for (n = 1:4)
			if (get(handles.hRadios(m,n),'Val'))	% Found the selected one
				break_too = true;
				break
			end
		end
		if (break_too),		break,	end				% Dirty trick because of stupid paranoia against GOTOs
	end

	% Create the file name
	periods = {'hour', 'day', 'week', 'month'};
	mags = {'all_', '1.0_', '2.5_', '4.5_', 'significant_'};
	name = [mags{m} periods{n}]

	dest_fiche = [handles.path_tmp 'usgs.csv'];
	url = ['http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/' name '.csv'];
	system(['wget "' url '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
	r = csv2cell(dest_fiche);
	lon = str2double(r(2:end,3));	lat = str2double(r(2:end,2));
	mag = uint8(str2double(r(2:end,5)) * 10);
	dep = uint16(str2double(r(2:end,4)) * 10);
	h = line('Parent',handles.hMirAx, 'XData',lon, 'YData',lat, 'LineStyle','none', 'Marker','o', ...
		'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'Tag','Earthquakes');
	setappdata(h_quakes,'SeismicityDepth',dep);			% Save events depth
	setappdata(h_quakes,'SeismicityMag',mag);			% Save events magnitude
	%setappdata(h_quakes,'SeismicityTime',year_dec);		% Save events time

	draw_funs(h,'Earthquakes',[])

% --------------------------------------------------------------------------------------
function radios_CB(obj, evt)
% ...
	handles = guidata(obj);
	if (~get(obj,'Val')),	set(obj,'Val',1),	return,		end
	for (k = 1:numel(handles.hRadios))
		set(handles.hRadios(k),'Val',0)		% Easier to set them all to 0 and next put current one to 1
	end
	set(obj,'Val',1)

function varargout = find_clusters(varargin)
% Helper Windows to find seismic clusters

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

	hObject = figure('Vis','off');
	find_clusters_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right');

	handles.mirone_fig = varargin{1};
	handles.h_polyg = varargin{2};

	handles.h_events = findobj(handles.mirone_fig,'Tag','Earthquakes');
	if (isempty(handles.h_events))      % Should issue an error message
        delete(hObject);    return
	end

	handles.clusters = [];

	% Get the region that encloses the polygon
	x = get(handles.h_polyg,'XData');   y = get(handles.h_polyg,'YData');
	x_min = min(x);     x_max = max(x);
	y_min = min(y);     y_max = max(y);
	handles.region = [x_min x_max y_min y_max];

	% Update handles structure
	guidata(hObject, handles);

	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

%--------------------------------------------------------------------------
function edit_dt_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),      set(hObject,'String','3');      end

%--------------------------------------------------------------------------
function edit_nEvents_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),     set(hObject,'String','30');      end

%--------------------------------------------------------------------------
function edit_maxSTD_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),      set(hObject,'String','1.5');      end

%--------------------------------------------------------------------------
function radio_bigEvent_CB(hObject, handles)
	if (get(hObject,'Value'))
        set(handles.radio_swarmCenter,'Value',0)
	else
        set(hObject,'Value',1)
	end

%--------------------------------------------------------------------------
function radio_swarmCenter_CB(hObject, handles)
	if (get(hObject,'Value'))
        set(handles.radio_bigEvent,'Value',0)
	else
        set(hObject,'Value',1)
	end

%--------------------------------------------------------------------------
function push_compute_CB_(hObject, handles)
	events_time = getFromAppdata(handles.h_events,'SeismicityTime');
	events_mag = getFromAppdata(handles.h_events,'SeismicityMag');
	events_dep = getFromAppdata(handles.h_events,'SeismicityDepth');
	if (isempty(events_time))
        errordlg('This seismicity plot somehow lost its time information.','Error');    return;
	end

	gap = str2double(get(handles.edit_dt,'String'));
	nEvents = str2double(get(handles.edit_nEvents,'String'));

	x = get(handles.h_events,'XData');      y = get(handles.h_events,'YData');
	IN = inpolygon(x,y,get(handles.h_polyg,'XData'),get(handles.h_polyg,'YData'));
	x(~IN) = [];    y(~IN) = [];    events_time(~IN) = [];
	events_mag(~IN) = [];           events_dep(~IN) = [];

	[events_time,ind] = sort(events_time);
	dt = diff(events_time);
	jumps = find(dt > gap/365) + 1;         % Find time jumps in the seismicity
	clus1 = (diff(jumps) >= nEvents);       % Find which jumps have at least nEvents (logical vector with 0s & 1s)
	clus2 = find(clus1);                    % Isolate the 1s (jumps(clus) == starting index of selected jumps)
	juju = [jumps(clus2)'; jumps(clus2+1)'];% [2 x n_clusters] array where 1st row = cluster_start & 2nd row = cluster_end+2
	ju = juju(:);
	if ( (length(dt) - jumps(end)) >= nEvents)  % Check that we are not loosing one last cluster
        ju = [ju; jumps(end); length(events_time)+2];
	end
	n_clust = length(ju) / 2;				% Divide by 2 because ju contains the begining AND end of each swarm

	set(handles.text_nFound,'String',['Found ' num2str(n_clust) ' potential clusters'],'Visible','on')
	if (n_clust == 0)						% No clusters found
        set(handles.popup_nSwarm,'String',' ','Value',1);
        handles.clusters = [];
        return
	end

	x = x(ind)';     y = y(ind)';
	events_time = events_time(ind);
	events_mag = events_mag(ind);
	events_dep = events_dep(ind);
	handles.clusters = cell(n_clust,1);
	handles.times = cell(n_clust,1);
	handles.mags = cell(n_clust,1);
	handles.depths = cell(n_clust,1);
	l = 1;
	for (k=1:2:2*n_clust)                   % Jump 2 because each swarm starts a 1, 3, 5, ... and stops at 2, 4, 6, ...
        handles.clusters{l} = [x(ju(k):ju(k+1)-2) y(ju(k):ju(k+1)-2)];  % -2 because diff reduces vector length by 1
        handles.times{l} = events_time(ju(k):ju(k+1)-2);
        handles.mags{l} = events_mag(ju(k):ju(k+1)-2);
        handles.depths{l} = events_dep(ju(k):ju(k+1)-2);
        l = l + 1;
	end

	std_max = str2double(get(handles.edit_maxSTD,'String'));
	str = cell(n_clust+1,1);                % + 1 to account for the first empty line
	for (k=1:n_clust)
        [std_dist,dist] = std_geo(handles.clusters{k}(:,1),handles.clusters{k}(:,2));
        ind = (dist > std_max*std_dist);
        handles.clusters{k}(ind,:) = [];
        handles.times{k}(ind,:) = [];
        handles.mags{k}(ind,:) = [];
        handles.depths{k}(ind,:) = [];
        str{k+1} = [num2str(k) 'th swarm (' num2str(length(handles.clusters{k})) ')'];
	end
	set(handles.popup_nSwarm,'String',str,'Value',1)

	guidata(hObject, handles);

%--------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
	events_time = getFromAppdata(handles.h_events,'SeismicityTime');
	events_mag = getFromAppdata(handles.h_events,'SeismicityMag');
	events_dep = getFromAppdata(handles.h_events,'SeismicityDepth');
	if (isempty(events_time))
        errordlg('This seismicity plot somehow lost its time information.','Error');    return;
	end

	gap = str2double(get(handles.edit_dt,'String'));
	nEvents = str2double(get(handles.edit_nEvents,'String'));
	std_max = str2double(get(handles.edit_maxSTD,'String'));

	x = get(handles.h_events,'XData');      y = get(handles.h_events,'YData');
	if (iscell(x))
        x = cat(2,x{:});        y = cat(2,y{:});
	end
	IN = inpolygon(x,y,get(handles.h_polyg,'XData'),get(handles.h_polyg,'YData'));
	IN = ~IN;       % Negate bacause what we are interested in is to remove the out-polygon-points
	x(IN) = [];     y(IN) = [];    events_time(IN) = [];
	events_mag(IN) = [];           events_dep(IN) = [];

	if (get(handles.radio_bigEvent,'Value'))
        mode = 'bigest';
	elseif (get(handles.radio_swarmCenter,'Value'))
        mode = 'center';
	end

	[handles.clusters,handles.times,handles.mags,handles.depths] = ...
        swarms(handles,x,y,events_time,events_mag,events_dep,gap,nEvents,std_max,mode,1);
	
	n_clust = length(handles.clusters);
	set(handles.text_nFound,'String',['Found ' num2str(n_clust) ' potential clusters'],'Visible','on')
	if (~n_clust)                           % No clusters found
        set(handles.popup_nSwarm,'String',' ','Value',1);
        handles.clusters = [];
        guidata(hObject, handles);
        return
	end

	str = cell(n_clust+1,1);                % + 1 to account for the first empty line
	for (k=1:n_clust)
        str{k+1} = [num2str(k) 'th swarm (' num2str(length(handles.clusters{k})) ')'];
	end
	set(handles.popup_nSwarm,'String',str,'Value',1)

	guidata(hObject, handles);

%--------------------------------------------------------------------------
function [clusters,times,mags,depths] = swarms(handles,x,y,events_time,events_mag,events_dep,gap,nEvents,std_max,mode,first)
% MODE = 'center'  rejects points that are further than STD_MAX from the geographical mean of each cluster
% MODE = 'bigest'  rejects points that are further than STD_MAX from the location of the bigest event of each cluster
% FIRST = 1 Use allways this when calling this function. Internally, FIRST is set to 0 and recursively use 
%           this function to refine the clusters determination.

	if (first)                              % At first time call we are not sure that data is already sorted
        [events_time,ind] = sort(events_time);
	end

	dt = diff(events_time);
	jumps = find(dt > gap/365) + 1;         % Find time jumps in the seismicity
	if (isempty(jumps))                     % No jumps > gap, It means no cluster or all data is a single cluster.
        clusters = [];  times = [];     mags = [];  depths = [];
        return
	end
	clus1 = (diff(jumps) >= nEvents);       % Find which jumps have at least nEvents (logical vector with 0s & 1s)
	clus2 = find(clus1);                    % Isolate the 1s (jumps(clus) == starting index of selected jumps)
	juju = [jumps(clus2)'; jumps(clus2+1)'];% [2 x n_clusters] array where 1st row = cluster_start & 2nd row = cluster_end+2
	ju = juju(:);
	if ( (length(dt) - jumps(end)) >= nEvents)  % Check that we are not loosing one last cluster
        ju = [ju; jumps(end); length(events_time)+2];
	end
	n_clust = length(ju) / 2;               % Divide by 2 because ju contains the begining AND end of each swarm

	if (n_clust == 0)						% No clusters found
        clusters = [];  times = [];     mags = [];  depths = [];
        return
	end

	if (first)								% Need to order also those vars according to index of sorting events_time
        x = x(ind)';     y = y(ind)';		% transpose because x & y come from a get(...,XData) and are row vectors
        events_mag = events_mag(ind);
        events_dep = events_dep(ind);
	end

	clusters = cell(n_clust,1);     times = cell(n_clust,1);
	mags = cell(n_clust,1);         depths = cell(n_clust,1);

	clipa = get(handles.checkbox_mainShock,'Value');

	l = 1;
	for (k=1:2:2*n_clust)                   % Jump 2 because each swarm starts a 1, 3, 5, ... and stops at 2, 4, 6, ...
        clusters{l} = [x(ju(k):ju(k+1)-2) y(ju(k):ju(k+1)-2)];  % -2 because diff reduces vector length by 1
        times{l} = events_time(ju(k):ju(k+1)-2);
        mags{l} = events_mag(ju(k):ju(k+1)-2);
        depths{l} = events_dep(ju(k):ju(k+1)-2);
        if (clipa && first)                  % If it is to clip, do it only on the first round
            [mmag,ind_tm] = max(mags{l});
            clusters{l}(1:ind_tm-1,:) = [];
            times{l}(1:ind_tm-1) = [];
            mags{l}(1:ind_tm-1) = [];
            depths{l}(1:ind_tm-1) = [];
        end
        l = l + 1;
	end

	% Select the center mode to use in std_geo
	if (strcmp(mode,'center'))
        center = cell(n_clust,1);
	elseif (strcmp(mode,'bigest'))
        center = cell(n_clust,1);
        for (k=1:n_clust)
            [mm,ind] = max(mags{k});
            center{k} = clusters{k}(ind,:); 
        end
	end

	for (k=1:n_clust)
        [std_dist,dist] = std_geo(clusters{k}(:,1),clusters{k}(:,2),center{k});
        ind = (dist > std_max*std_dist);
        clusters{k}(ind,:) = [];
        times{k}(ind,:) = [];
        mags{k}(ind,:) = [];
        depths{k}(ind,:) = [];
	end

	% The loop above may have rejected some points. A situation occurs when the rejected points
	% may have created a internal (to each cluster) time gap > "gap". The clever solution would
	% be to search only for cluster that have been changed. However, that would be quite cumbersome
	% and after all, this algo seams fast and not very memory consuming. So we will use a brute
	% force aproach. Cat all clusters together and re-run the swarms function.
	%
	% But there is more. If the clip option was used, many events may have gone in each cluster

	if (first)
        % x = cat(1,clusters{:}(:,1));      % This gives an error. WHY?
        zz = cat(1,clusters{:});
        x = zz(:,1);    y = zz(:,2);    clear zz;
        events_time = cat(1,times{:});
        events_mag = cat(1,mags{:});
        events_dep = cat(1,depths{:});
        [clusters,times,mags,depths] = ...
            swarms(handles,x,y,events_time,events_mag,events_dep,gap,nEvents,std_max,mode,0);
	end

%--------------------------------------------------------------------------
function radio_separate_CB(hObject, handles)
	if (get(hObject,'Value'))
        set(handles.radio_same,'Value',0)
	else
        set(handles.radio_same,'Value',1)
	end

%--------------------------------------------------------------------------
function radio_same_CB(hObject, handles)
	if (get(hObject,'Value'))
        set(handles.radio_separate,'Value',0)
	else
        set(handles.radio_separate,'Value',1)
	end

%--------------------------------------------------------------------------
function push_plot_CB(hObject, handles)

	if (isempty(handles.clusters))
        warndlg('No, I won''t plot. The reason why should be pretty obvious.','Warning')
        return
	end

	if (get(handles.radio_same,'Value'))    
        h_mir = new_mir(handles.mirone_fig, handles.h_polyg, handles);
        set(h_mir,'Name','All swarms')
        hold on;
		for (k=1:length(handles.clusters))
            h_quakes = plot(handles.clusters{k}(:,1),handles.clusters{k}(:,2),'kp','Marker','o', ...
                'MarkerFaceColor',rand(1,3),'MarkerEdgeColor','k','MarkerSize',5,'Tag','Earthquakes');
            setappdata(h_quakes,'SeismicityTime',handles.times{k});     % Save events time
            setappdata(h_quakes,'SeismicityMag',handles.mags{k});       % Save events mags
            setappdata(h_quakes,'SeismicityDepth',handles.depths{k});   % Save events depths
            draw_funs(h_quakes,'Earthquakes',[])
		end
        hold off;
	else            % Plot one swarm per figure
		h_mir = zeros(1, numel(handles.clusters));
        for (k = 1:numel(handles.clusters))
            h_mir(k) = new_mir(handles.mirone_fig, handles.h_polyg, handles);
            set(h_mir(k),'Name',[num2str(k) 'th  swarm (' num2str(length(handles.clusters{k})) ')'])
            hold on;
            h_quakes = plot(handles.clusters{k}(:,1),handles.clusters{k}(:,2),'kp','Marker','o', ...
                'MarkerFaceColor',rand(1,3),'MarkerEdgeColor','k','MarkerSize',5,'Tag','Earthquakes');
            setappdata(h_quakes,'SeismicityTime',handles.times{k});     % Save events time
            setappdata(h_quakes,'SeismicityMag',handles.mags{k});       % Save events mags
            setappdata(h_quakes,'SeismicityDepth',handles.depths{k});   % Save events depths
            draw_funs(h_quakes,'Earthquakes',[])
            hold off;
        end
	end

%--------------------------------------------------------------------------
function h_mir = new_mir(h_mirone_fig, h_polyg, handles)
	% Create a new Mirone figure with the image extents determined by the polygon limits
	% It makes use of the Mirone ImageCrop_CB function, and also of the "register image"
	% solution of draw_funs.

	mirone('ImageCrop_CB',guidata(h_mirone_fig),h_polyg,'CropaWithCoords');
	set(0,'ShowHiddenHandles','on')
	h_mir = findobj(0,'Type','figure','Name','Croped Image');
	set(0,'ShowHiddenHandles','off')

	if (numel(h_mir) > 1)       % If for bad luck we have more than one, pick the
        h_mir = h_mir(1);       % first because it seams that it's the newest
	end
	figure(h_mir)               % Make it the current figure (crutial)

%--------------------------------------------------------------------------
function popup_nSwarm_CB(hObject, handles)
	val = get(hObject,'Value');
	if (val == 1),   return;     end;    % First in list is empty
	contents = get(hObject,'String');
	tit = contents{get(hObject,'Value')};
	tit = [tit(1:end-1) ' events )'];
	val = val - 1;
	h_mir = new_mir(handles.mirone_fig, handles.h_polyg, handles);
	set(h_mir,'Name',tit)

	% hold on
	% h_quakes = plot(handles.clusters{val}(:,1),handles.clusters{val}(:,2),'kp','Marker','o', ...
	%     'MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',5,'Tag','Earthquakes');
	% hold off

	mirAx = get(h_mir,'CurrentAxes');
	h_quakes = line('XData',handles.clusters{val}(:,1),'YData',handles.clusters{val}(:,2),'Parent',mirAx,'LineStyle','none','Marker','o', ...
        'MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',5,'Tag','Earthquakes');

	setappdata(h_quakes,'SeismicityTime',handles.times{val});     % Save events time
	setappdata(h_quakes,'SeismicityMag',handles.mags{val});       % Save events mags
	setappdata(h_quakes,'SeismicityDepth',handles.depths{val});   % Save events depths
	draw_funs(h_quakes,'Earthquakes',[])

%--------------------------------------------------------------------------
function push_save_CB(hObject, handles)

	if (isempty(handles.clusters))
        warndlg('No, I won''t save. The reason why should be pretty obvious.','Warning')
        return
	end

	%save_seismicity(handles.mirone_fig, h_events)

	handMir = guidata(handles.mirone_fig);				% Get the Mirone handles
	[FileName,PathName] = put_or_get_file(handMir, ...
		{'*.dat;*.DAT', 'Seimicity file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'}, 'Select File name','put','.dat');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	fid = fopen(fname, 'w');
	if (fid < 0),    errordlg(['Can''t open file:  ' fname],'Error');    return;     end

	for (k=1:length(handles.clusters))
		year = fix(handles.times{k});       % integer part
		frac = handles.times{k} - year;     % fraction part
		jdYear = date2jd(year);             % True Julian day of the 1st January of YEAR
		frac = frac .* (365 + isleapyear(year));
		[dumb,month, day, hour, minute] = jd2date((jdYear+frac),[]);

		fprintf(fid,'>\n');
		fprintf(fid,'%.3f\t%.3f\t%.1f\t%.1f\t%d\t%02d\t%02d\t%02d\t%02d\n',[handles.clusters{k} double(handles.mags{k})/10,...
			double(handles.depths{k})/10 year month day hour minute]');
	end
	fclose(fid);

%--------------------------------------------------------------------------
function t = isleapyear(year)
	t = ( ~rem(year, 4) & rem(year, 100) ) | ~rem(year, 400);

%--------------------------------------------------------------------------
function [std_dist, dist] = std_geo(lon,lat,method,center)
	%  [std_dist,dist] = STD_GEO(lon,lat) computes the average standard distance for
	%  geographic data using the spherical aproximation. Distances are returned in degrees.
	%  DIST is an array with the same size as LON and contains the distance between LON,LAT
	%  and the mean LON, mean LAT
	%
	%  METHOD = 'linear' computes the average distance; 'quadratic' computes the sqrt of the
	%  average squared distance; Default is 'quadratic'
	%
	%  OPTION
	%  CENTER = [] makes the function ignore this option (created for progamatically reasons)
	%  CENTER = [center_x center_y] computes the distances and std with respect to this point.
	%           That is, it doesn't compute the geographical mean position mean LON and mean LAT 
	
	D2R = pi / 180;

	nins = nargin;
	if (nins == 2)
        method = 'quadratic';
        do_mean = 1;
	elseif (nins == 3)
        if (ischar(method) && isempty(strmatch(method,{'linear' 'quadratic'})) )
            error('STD_GEO: Unrecognized method string');
        elseif (isnumeric(method) && isempty(method))
            method = 'quadratic';
            do_mean = 1;
        elseif (length(method) == 2)
            method = 'quadratic';
            lon_mean = method(1);
            lat_mean = method(2);
            do_mean = 0;
        else
            error('STD_GEO: Uncompreensible 3th input var. It must either be a "mthod" string or a [x,y] point.')
        end
	elseif (nins == 4)
        if (isempty(method))
            do_mean = 1;
        elseif (length(method) == 2)
            lon_mean = method(1);
            lat_mean = method(2);
            do_mean = 0;
        else
            error('STD_GEO: Uncompreensible 4th input var. It must either be "NaN" or [center_x center_y]')
        end
	end

	% Make sure lon & lat are column vectors
	if (size(lat,1)==1),    lat=lat(:);     lon=lon(:);     end

	if (do_mean)        %  Compute the mean location
        [lon_mean,lat_mean] = mean_geog(lon,lat);
	end

	dist = sin(lat_mean*D2R) * sin(lat*D2R) + cos(lat_mean*D2R) * cos(lat*D2R) .* cos((lon - lon_mean)*D2R);
	dist(dist >  1) =  1;     dist(dist < -1) = -1;
	dist = acos(dist) / D2R;

	% Compute the average range using the appropriate method
	switch  method
        case 'linear',       std_dist = mean(dist);
        case 'quadratic',    std_dist = sqrt(mean(dist.^2));
	end

%--------------------------------------------------------------------------
function [lon_mean,lat_mean] = mean_geog(lon,lat)
	% MEAN_GEOG(lon,lat) computes means for geographic data using the spherical aproximation

	D2R = pi / 180;
	lat = lat * D2R;        lon = lon * D2R;    % Convert to radians.

	%  Compute the centroid by summing all cartesian data.
	[x,y,z] = sph2cart(lon,lat,ones(size(lat)));
	[lon_mean,lat_mean] = cart2sph(sum(x),sum(y),sum(z));

	% Set longitude in [-pi; pi] range
	lon_mean = pi*((abs(lon_mean)/pi) - 2*ceil(((abs(lon_mean)/pi)-1)/2)) .* sign(lon_mean);
	lon_mean = lon_mean / D2R;
	lat_mean = lat_mean / D2R;

% -------------------------------------------------------------------------------
function out = getFromAppdata(hand,tag)
% If seismicity has been added from more than one file length(hand) > 1 and
% we cannot jus do out = getappdata(hand,tag)
	out = [];
	for i = 1:length(hand)
        out = [out; getappdata(hand(i),tag)];
	end

% --- Creates and returns a handle to the GUI figure. 
function find_clusters_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Find Clusters',...
'NumberTitle','off',...
'Position',[520 597 399 203],...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','Call',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'Position',[21 154 47 21],...
'String','3',...
'Style','edit',...
'Tooltip','Events separated more than this number of days do not belong to the same cluster',...
'Tag','edit_dt');

uicontrol('Parent',h1,'Position',[15 180 65 15],'String','Max time gap','Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'Position',[131 154 47 21],...
'String','30',...
'Style','edit',...
'Tooltip','A swarm must contain at least this number of events',...
'Tag','edit_nEvents');

uicontrol('Parent',h1,'Position',[125 180 65 15],'String','Min events','Style','text','Tag','text2');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'Position',[241 154 47 21],...
'String','1.5',...
'Style','edit',...
'Tooltip','Events distant more than this standard deviation do not belong to the same cluster',...
'Tag','edit_maxSTD');

uicontrol('Parent',h1,'Position',[235 180 57 15],'String','Max STD','Style','text','Tag','text3');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[10 130 116 15],...
'String','Center at mainshock',...
'Style','radiobutton',...
'Tooltip','Reject events that distant more than Max STD from the main shock position',...
'Tag','radio_bigEvent');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[137 130 87 15],...
'String','Swarm center',...
'Style','radiobutton',...
'Tooltip','Reject events that distant more than Max STD from the swarm center',...
'Value',1,'Tag','radio_swarmCenter');

uicontrol('Parent',h1,'Position',[250 130 111 15],...
'String','Clip at mainshock',...
'Style','checkbox',...
'Tooltip','Reject events occuring before the main shock',...
'Tag','checkbox_mainShock');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[321 152 66 23],'String','Compute','Tag','push_compute');

uicontrol('Parent',h1,'FontSize',10,...
'HorizontalAlignment','left','Position',[12 101 171 16],...
'String','Found 0 potential swarms',...
'Style','text','Tag','text_nFound','Visible','off');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[10 73 184 15],...
'String','Plot swarms in separate windows',...
'Style','radiobutton',...
'Value',1,'Tag','radio_separate');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[10 43 192 15],...
'String','Plot all swarms in the same window',...
'Style','radiobutton',...
'Tooltip','Will use different colors to distinguish swarms',...
'Tag','radio_same');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'Position',[220 46 91 22],...
'Style','popupmenu',...
'String',{' '},...
'Tooltip','Plot one candidate swarm',...
'Value',1,'Tag','popup_nSwarm');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[320 56 66 23],...
'String','Plot',...
'Tooltip','Plot the swarms in separate windows',...
'Tag','push_plot');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[320 15 66 23],...
'String','Save',...
'Tooltip','Save the the swarms as a multisegment file',...
'Tag','push_save');

uicontrol('Parent',h1,'Position',[220 70 81 15],'String','Select swarm','Style','text','Tag','text9');

function main_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

function histos_seis(h_polyg,opt)
% HISTOS_SEIS plot several types of histograms and oher line plots using seismic data
%
% H_POLYG is either a handle to a closed polygon, OR a handle to the seismicity dots
% OPT = 'GR' compute a Guttemberg & Richter relation
% OPT = 'CH' compute a Cumulative Number Histogram
% OPT = 'CM' compute a Cumulative Moment Histogram
% OPT = 'MH' compute a Magnitude Histogram
% OPT = 'TH' compute a Time Histogram
% OPT = 'HH' compute a Hour of the day Histogram
% OPT = 'TM' Plot a Time Magnitude relation
% OPT = 'TD' Plot a Time Depth relation
% OPT = 'BV' Mc and b estimate
% OPT = 'OL' Fit Omori Law
% OPT = 'HT' Display a Table histogram

%	Copyright (c) 2004-2022 by J. Luis
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

	h_mir_fig = gcf;
	h_events = findobj(h_mir_fig,'Tag','Earthquakes');
	if (any(strcmp( opt,{'GR' 'CM' 'MH' 'TM' 'HT'}) ))
		events_mag = double(getFromAppdata(h_events,'SeismicityMag')) / 10;
		if (isempty(events_mag)),		return,		end
	else
		events_time = getFromAppdata(h_events,'SeismicityTime');
		if (isempty(events_time)),		return,		end
	end
	tag = get(h_polyg,'tag');		IN = [];
	if (strcmp(tag,'SeismicPolyg'))			% Case of a closed polygon
		x = get(h_events,'XData');		y = get(h_events,'YData');
		if (iscell(x))
			x = cat(2,x{:});			y = cat(2,y{:});
		end
		IN = inpolygon(x,y,get(h_polyg,'XData'),get(h_polyg,'YData'));
	end
	
	if (strcmp(opt,'GR') || strcmp(opt,'MH'))
		% The other options creates their own window
		h_fig = figure('NumberTitle','off','Visible','off','Color',get(0,'factoryUicontrolBackgroundColor'));
	end

	if (strcmp(opt,'GR'))			% Gutt & Richter
		if (~isempty(IN)),      events_mag(~IN) = [];   end         % INpolygon option was used
		minmag = min(events_mag);	maxmag = ceil(10*max(events_mag))/10;
		if (minmag > 100)			% Hydrofone SL magnitude
			X = minmag:1:maxmag;
		else
			X = minmag:0.1:maxmag;
		end
		N = histo_m('hist',events_mag,X);
		N = cumsum(N(length(N):-1:1));    % N for M >= (counted backwards)
		semilogy(X(end:-1:1),N,'kp','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5)
		tit = 'Guttenberg and Richter';    lab_x = 'Magnitude';    lab_y = 'Cumulative Number';
	elseif (strcmp(opt,'CM'))        % Cumulative Moment histogram
		events_time = getFromAppdata(h_events,'SeismicityTime');
		if (~isempty(IN)),      events_mag(~IN) = [];   events_time(~IN) = [];   end     % INpolygon option was used
		IN = (events_mag <= 0);
		events_mag(IN) = [];    events_time(IN) = [];   % Remove events with 0 magnitude
		[events_time, events_mag] = merge_classes(events_time, events_mag);		% Take case of the case where we have > 1 class.
		c = cumsum( 10.^(1.5*events_mag + 16.1));
		call_ecran(events_time, c, false, 'Cumulative Moment', 'Time', 'Cumulative Moment')
		return
	elseif (strcmp(opt,'MH'))        % Magnitude histogram
		if (~isempty(IN)),      events_mag(~IN) = [];   end         % INpolygon option was used
		events_mag(~(events_mag > 0)) = [];             % Remove events with 0 magnitude
		X = floor(min(events_mag)):0.2:ceil(max(events_mag));
		[N,xout] = histo_m('hist',events_mag,X);
		histo_m('bar',xout,N,'hist');
		tit = 'Magnitude histogram';    lab_x = 'Magnitude';    lab_y = 'Number';
	elseif (strcmp(opt,'HT'))        % Histogram table
		if (~isempty(IN)),      events_mag(~IN) = [];  end			% INpolygon option was used
		X = (.25:.5:9.25);
		[N,xout] = histo_m('hist',events_mag,X);
		txt_mags = {'[0 - 0.5]' ']0.5 - 1.0]' ']1.0 - 1.5]' ']1.5 - 2.0]' ']2.0 - 2.5]' ']2.5 - 3.0]' ']3.0 - 3.5]' ...
			']3.5 - 4.0]' ']4.0 - 4.5]' ']4.5 - 5.0]' ']5.0 - 5.5]' ']5.5 - 6.0]' ']6.0 - 6.5]' ...
			']6.5 - 7.0]' ']7.0 - 7.5]' ']7.5 - 8.0]' ']8.0 - 8.5]' ']8.5 - 9.0]' ']9.0 - 9.5]'};
		tabHist = cell(19,2);
		tabHist(:,1) = txt_mags';
		for (m=1:19),       tabHist{m,2} = N(m);    end
		tableGUI('array',tabHist,'ColWidth',[90 90],'ColNames',{'Magnitudes' 'Number of events'},...
			'FigName','Kind of histogram','RowNumbers','y','MAX_ROWS',19,'modal','');
		return
	elseif (strcmp(opt,'CH'))        % Cumulative histogram
		if (~isempty(IN)),      events_time(~IN) = [];  end         % INpolygon option was used
		events_time = sort(events_time);			% If we have classes of sizes time is not monotonically increasing
		call_ecran(events_time, 1:length(events_time), false, 'Cumulative Number', 'Time', 'Cumulative Number')
		return
	elseif (strcmp(opt,'TH'))        % Time histogram
		if (~isempty(IN)),      events_time(~IN) = [];  end         % INpolygon option was used
		call_ecran(events_time, [], true, 'Time histogram', 'Time', 'Number')
		return
		% elseif (strcmp(opt,'HH'))        % Hour of day histogram
		%     if (~isempty(IN))       events_time(~IN) = [];  end     % INpolygon option was used
		%     X = 0:1:24;
		%     N = cumsum(histo_m('hist',events_time,X));
		%     plot(X,N,'b','LineWidth',2)
		%     tit = 'Cumulative histogram';  lab_x = 'Time in years';    lab_y = 'Cumulative Number';
	elseif (strcmp(opt,'TM'))        % Time Magnitude relation
		events_time = getFromAppdata(h_events,'SeismicityTime');
		if (~isempty(IN)),      events_mag(~IN) = [];   events_time(~IN) = [];   end    % INpolygon option was used
		IN = (events_mag <= 0);
		events_mag(IN) = [];    events_time(IN) = [];   % Remove events with 0 magnitude
		[events_time, events_mag] = merge_classes(events_time, events_mag);		% Take case of the case where we have > 1 class.
		call_ecran(events_time, events_mag, false, 'Magnitude vs Time', 'Time', 'Magnitude', {'LineStyle','none', 'Marker','o', 'MarkerFaceColor','k', 'MarkerSize',2})
		return
	elseif (strcmp(opt,'TD'))        % Time Depth relation
		events_dep = double(getFromAppdata(h_events,'SeismicityDepth')) / 10;
		if (~isempty(IN)),   events_time(~IN) = [];   events_dep(~IN) = [];     end      % INpolygon option was used
		[events_time, events_dep] = merge_classes(events_time, events_dep);		% Take case of the case where we have > 1 class.
		call_ecran(events_time, -events_dep, false, 'Depth vs Time', 'Time', 'Depth', {'LineStyle','none', 'Marker','o', 'MarkerFaceColor','k', 'MarkerSize',2})
		return
	elseif (strcmp(opt,'BV'))        % b-value
		events_mag = double(getFromAppdata(h_events,'SeismicityMag')) / 10;
		if (~isempty(IN)),   events_time(~IN) = [];   events_mag(~IN) = [];     end      % INpolygon option was used
		IN = (events_mag <= 0);
		events_mag(IN) = [];    events_time(IN) = [];   % Remove events with 0 magnitude
		[events_time, events_mag] = merge_classes(events_time, events_mag);		% Take case of the case where we have > 1 class.
		ud.events_time = events_time;    ud.events_mag = events_mag;    ud.h_mir_fig = h_mir_fig;
		mc_e_b_estimate('in',ud)
	elseif (strcmp(opt,'OL'))        % Omori law
		events_mag = double(getFromAppdata(h_events,'SeismicityMag')) / 10;
		if (~isempty(IN)),   events_time(~IN) = [];   events_mag(~IN) = [];     end      % INpolygon option was used
		IN = (events_mag <= 0);
		events_mag(IN) = [];    events_time(IN) = [];   % Remove events with 0 magnitude
		[events_time, events_mag] = merge_classes(events_time, events_mag);		% Take case of the case where we have > 1 class.
		calc_omori(events_time,events_mag,h_mir_fig)
	end

	if (strcmp(opt,'GR') || strcmp(opt,'MH'))
		set(h_fig,'Name',tit,'Visible','on')
		set(get(gca,'XLabel'),'String',lab_x,'FontSize',10)
		set(get(gca,'YLabel'),'String',lab_y,'FontSize',10)
	end


% -------------------------------------------------------------------------------
function out = getFromAppdata(hand,tag)
% If seismicity has been added from more than one file length(hand) > 1 and
% we cannot jus do out = getappdata(hand,tag)
	out = [];
	for i = 1:length(hand)
    	out = [out; getappdata(hand(i),tag)];
	end

% -------------------------------------------------------------------------------
function [evt_t, evt_x] = merge_classes(evt_t, evt_x)
% When events were plotted in size classes the time jumps between classes (each goes from t_start to t_stop)
% This makes a mess in histograms, so we must merge the classes back again and asure the increases monotonically.
	ind = (diff(evt_t) < 0);
	if (any(ind))			% Yes, we have classes and hence the time jumps in 'evt_t'
		[evt_t, ind] = sort(evt_t);
		evt_x = evt_x(ind);
	end

% -------------------------------------------------------------------------------
function call_ecran(evt_t, evt_x, is_hist, title, lab_x, lab_y, PV)
% Call the ecran figure and set the xx axes as time
	Y = fix(evt_t);
	dec = (evt_t - Y);
	days_in_year = (datenum(Y+1,1,1) - datenum(Y,1,1)) .* dec;
	serial_date = datenum(Y,1,1) + days_in_year;

	if (is_hist)					% Make an histogram
		hFig = ecran;				% Create a new Ecran Fig to hold the histogram
		delete(findobj(hFig, 'Tag', 'isocs_but'))	% Others uis should be removed too
		set(hFig, 'Name', title)
		n_bins = sshist(serial_date);		% Automatic bining
		hP = histo_m('hist', serial_date, n_bins, 'hands');
		set(hP, 'Tag', 'Histogram')			% To be used by write_gmt_script
		setappdata(hP, 'xy', serial_date)	% Store it for use in write_gmt_script
	else
		if (nargin == 7)
			hFig = ecran(serial_date, evt_x, title, PV);	% Use also the PV property/value pairs
		else
			hFig = ecran(serial_date, evt_x, title);
		end
	end
	hAxes = findobj(hFig, 'Tag', 'axes1');
	setappdata(hAxes, 'LabelFormatType', 'Date');		% Tell pixval_stsbar to display XX coords as Date-Time
	h = findobj(hFig,'Tag','add_uictx');
	cb = get(h, 'Callback');
	feval(cb, h, guidata(hFig))			% Call the ecran's add_uictx_CB function

	set(get(hAxes,'XLabel'),'String',lab_x,'FontSize',10)
	set(get(hAxes,'YLabel'),'String',lab_y,'FontSize',10)

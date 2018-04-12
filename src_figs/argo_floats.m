function varargout = argo_floats(varargin)
% Helper figure to plot ARGO buoys position and data in function of time

%	Copyright (c) 2004-2018 by J. Luis
%
%             DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%                     Version 2, December 2004
% 
%  Everyone is permitted to copy and distribute verbatim or modified
%  copies of this license document, and changing it is allowed as long
%  as the name is changed.
% 
%             DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%    TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
% 
%   0. You just DO WHAT THE FUCK YOU WANT TO.
%
% BUT NOTE: This license does not apply to the included argofiles() and argodata()
%			functions that have their own author and license (BSD)
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id$

	if (nargin > 1 && ischar(varargin{1}))
		gui_CB = str2func(varargin{1});
		[varargout{1:nargout}] = feval(gui_CB,varargin{2:end});
	else
		h = argo_floats_OF(varargin{:});
		if (nargout),	varargout{1} = h;   end
	end

% ---------------------------------------------------------------------------------
function hObject = argo_floats_OF(varargin)

	hObject = argo_floats_LayoutFcn;
	handles = guihandles(hObject);

 	handles.hMirFig = varargin{1};
	move2side(hObject, 'right')

	% Set default start and end dates. Don't really care about being 100% accurate at end date (end of month)
	mos = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
	[yr, mo] = datevec(now);
	t1 = sprintf('01-%s-%d', mos{mo}, yr);	
	t2 = t1;
	if (mo == 2),	t2(1:2) = '28';
	else,			t2(1:2) = '30';
	end
	set(handles.edit_startDate, 'Str', t1)
	set(handles.edit_endDate,   'Str', t2)
	handles.startDate_back = t1;
	handles.endDate_back = t2;

	set(hObject,'Vis','on');
	guidata(hObject, handles);
	
	if (nargin > 1),	external_drive(handles, 'argo_floats', varargin{2:end}),	end

% -------------------------------------------------------------------------------------
function edit_startDate_CB(hObject, handles)
	t1 = datenum(get(hObject, 'Str'));
	if (isnan(t1) || t1 < datenum('01-Jan-1998'))
		set(hObject, 'Str', handles.startDate_back)
	end

% -------------------------------------------------------------------------------------
function push_startDate_CB(hObject, handles)
	new_date = uisetdate;
	if (isempty(new_date)),		return,		end
	yr = str2double(new_date(8:end));
	if (yr < 1998)
		errordlg('Argo data is not availabe for this date.','Error');     return
	end
	set(handles.edit_startDate, 'Str', new_date)

% -------------------------------------------------------------------------------------
function edit_endDate_CB(hObject, handles)
	t1 = datenum(get(hObject, 'Str'));
	if (isnan(t1) || t1 < datenum('01-Jan-1998'))
		set(hObject, 'Str', handles.endDate_back)
	end

% -------------------------------------------------------------------------------------
function push_endDate_CB(hObject, handles)
	new_date = uisetdate;
	if (isempty(new_date)),		return,		end
	yr = str2double(new_date(8:end));
	if (yr < 1998)
		errordlg('Argo data is not availabe for this date.','Error');     return
	end
	set(handles.edit_endDate, 'Str', new_date)

% -------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...
try
	t1 = datenum(get(handles.edit_startDate, 'Str'));
	t2 = datenum(get(handles.edit_endDate, 'Str'));
	
	if (t1 > t2)
		errordlg('Start date is posterior to end date. Can''t work, right?', 'Error')
		return
	end

	handMir = guidata(handles.hMirFig);
	if (aux_funs('msg_dlg',50,handMir))		% If no_file create one.
		handMir = guidata(handMir.figure1);
	end
	x_lim = get(handMir.axes1,'Xlim');	y_lim = get(handMir.axes1,'Ylim');
	str = get(handles.popup_basins, 'Str');		val = get(handles.popup_basins, 'Val');
	if (val == 4)
		warndlg('The "All" basins is not yet implemented. Reverting to atlantic.', 'Warning')
		val = 1;
	end
	basin = str{val};
	[ncfiles,lat,lon,t] = argofiles(y_lim, x_lim, [t1 t2], basin);
	if (isempty(lat))
		warndlg('Something went wrong. No points inside data region', 'Warning'),	return
	end

	h = line('Parent',handMir.axes1,'XData',lon,'YData',lat,'LineStyle','none','Marker','o','MarkerFaceColor',[1 0.5490196 0],...
        'MarkerEdgeColor','w','MarkerSize',7,'Tag','ARGO');
	setappdata(h, 'time', t),	setappdata(h, 'ncfiles', ncfiles)
	draw_funs(h,'DrawSymbol')		% Set symbol's uicontextmenu
catch
	disp(lasterr)
end

% -------------------------------------------------------------------------------------
function [ncfiles,lat,lon,t] = argofiles(latlim,lonlim,tlim,basin)
% argofiles returns urls of ARGO NetCDF files within specified
% geographic and temporal limits.
% 
% Usage note: Processing time is not affected by geographic limits, but
% a wide temporal range will slow things down. Expect roughly a few seconds
% per month in your specified tlim. 
% 
%% Syntax
% 
% ncfiles = argodata(latlim,lonlim,tlim,basin)
% [ncfiles,lat,lon,t] = argodata(...)
% 
%% Description 
% 
% ncfiles = argodata(latlim,lonlim,tlim,basin) searches the internet for ARGO
% profiles that start within the geographic confines latlim,lonlim and temporal
% range tlim.  latlim and lonlim can be two-element arrays to specify limits or a  
% polygon of three or more points if you only want floats within an arbitrary shape. 
% Times tlim can be a single date or a starting and ending date and can be
% in any format recognized by Matlab's datenum function. Due to the way files are 
% structured on the NODC you must also specify a basin as 'pacific', 'atlantic', 
% or 'indian'.  Output ncfiles is a cell array of urls.
% 
% [ncfiles,lat,lon,t] = argodata(...) also returns arrays of lat, lon, and times for
% each nc file. 
% 
%% Example 1: In the Nino 3.4 box and in a range of days
% 
% latlim = [-5 5]; 
% lonlim = [-170 -120]; 
% t = [datenum('dec 25, 2005') datenum('jan 1, 2006')]; 
% f = argofiles(latlim,lonlim,t,'pacific'); 
% 
%% Example 2: In an arbitrary polygon on a specific day
% Perhaps you're interested in all the ARGO floats within some polygon on a specific day.
% Below we define a polygon by lats,lons, then find all the profiles inside
% the polygon on Valentines Day 2009, then plot all seven profile locations
% inside the polygon on that day: 
% 
% lats = [32;33;35;36;35;33;30;28;25;23;21;23;25;27;30;33;35;36;35;33]; 
% lons = [-47;-50;-53;-56;-61;-63;-63;-59;-54;-50;-48;-45;-39;-36;-33;-32;-33;-37;-41;-45];
% plot(lons,lats,'ro-') 
% 
% [ncfiles,lat,lon] = argofiles(lats,lons,'feb. 14, 2009','atlantic'); 
% plot(lon,lat,'gp')
% 
%% Example 3: Anywhere in the world on a given day
% Where were all the ARGO floats on Flag Day, 2000?  We'll have to query
% each ocean basin separately to answer this question: 
% 
% [nca,lata,lona] = argofiles([-90 90],[-180 180],'June 14, 2010','atlantic'); 
% [nci,lati,loni] = argofiles([-90 90],[-180 180],'June 14, 2010','indian'); 
% [ncp,latp,lonp] = argofiles([-90 90],[-180 180],'June 14, 2010','pacific'); 
% 
% hold on
% plot(lona,lata,'b.','markersize',10)
% plot(loni,lati,'k.','markersize',10)
% plot(lonp,latp,'r.','markersize',10)
%  
%% Author Info: 
% The argofiles function was written by Chad A. Greene of the University of
% Texas at Austin's Institute for Geophysics (UTIG), December 2015. 
% http://www.chadagreene.com.
% 
% See also argodata.
%
% Licenced under the BSD (FEX: 54503)
%
% JL: This function is modified to work in R13 and depends on GMT

	%% Error checks: 
	if (ismember(lower(basin),{'pacific','atlantic','indian'}) == 0)
		error('Unrecognized ocean basin. Must be ''pacific'', ''atlantic'', or ''indian')
	end

	%% Format latlim,lonlim: 

	% Wrap longitudes to 180: 
	lonlim(lonlim > 180) = lonlim(lonlim > 180) - 360; 

	% If two elements are input, assume user wants floats bound by geoquad:  
	if numel(latlim)==2
		latv = [latlim(1) latlim(1) latlim(2) latlim(2)]; 
		lonv = [lonlim(1) lonlim(2) lonlim(2) lonlim(1)]; 
	else
		latv = latlim; 
		lonv = lonlim; 
	end

	%% Format input dates: 

	% Make sure dates are in datenum format (will still work even if dates are already datenum):  
	dn = datenum(tlim); 

	% Adjust date ranges outside data bounds: 
	dn(dn < datenum('jan 1, 1998')) = datenum('jan 1, 1998'); 
	dn(dn > now) = now; 
	dn = dn(1):dn(end); 

	% Get year and months associated with every date in the temporal range: 
	[yr,mo] = datevec(dn);

	%% Define basin and abbreviation: 

	switch lower(basin) 
		case 'pacific'
			bsn = 'pa'; 
		case 'atlantic'
			bsn =  'at'; 
		case 'indian' 
			bsn = 'in'; 
		otherwise
			error('Unrecognized basin.') 
	end

	%% Loop through each possible month: 

	% Initialize arrays that we'll try to populate:  
	ncfiles = {};	latout = [];	lonout = [];	tout = []; 
	mos = {'-Jan-'; '-Feb-'; '-Mar-'; '-Apr-'; '-May-'; '-Jun-'; '-Jul-'; '-Aug-'; '-Sep-'; '-Oct-'; '-Nov-'; '-Dec-'};

	for (k = 1:numel(yr))  
		if (k > 1)				% If this month has already been logged, try the next one
			if (mo(k) == mo(k-1)),	continue,	end
		end

		try
			% Download data file for the month: 
			url = ['https://data.nodc.noaa.gov/argo/inv/basins/',basin,'/',num2str(yr(k),'%04.f'),'/',bsn,num2str(yr(k),'%04.f'),num2str(mo(k),'%02.f'),'_argoinv.txt'];

			s = gmtmex(['gmtconvert -h1 ' url]);
			s = char(s.text)';
			[file_urls,tmp,lat,lon] = strread(s, '%*s %s %s %*s %f %*f %f %*s %*s','delimiter', ',');

			% Remove empty rows:
			ind = cellfun('isempty',tmp);
			file_urls(ind) = [];
			lat(ind) = [];
			lon(ind) = [];
			tmp(ind) = [];

			% Convert date strings to datenums:
			%t = datenum(tmp,'yyyy-mm-dd');
			% R13 datenum is incredibly bugged so must change the times that come as: '2018-04-01T16:56:20Z'
			% into the 'dd-mmm-yyyy' format (e.g. '01-Apr-2018'). It all look awfully inefficient
			tt = char(tmp);
			ind = str2num(tt(:,6:7));
			mes = char(mos(ind));
			t = [tt(:,9:10) mes tt(:,1:4)];
			t = datenum(t,'dd-mmm-yyyy');

			% Get indices of all files with user-specified polygon and timeframe
			[inpoly,onpoly] = inpolygon(lon,lat,lonv,latv);
			inpoly = inpoly+onpoly;
			ind = t >= dn(1) & t <= dn(end) & inpoly;

			% Log any matching filenames:
			ncfiles(length(ncfiles)+1:length(ncfiles)+sum(ind),1) = strcat('https://data.nodc.noaa.gov/argo/',file_urls(ind));

			% Write geo coordinate and time arrays only if user wants them
			if (nargout > 1)
				latout = [latout;lat(ind)];
				lonout = [lonout;lon(ind)];
				tout = [tout;t(ind)];
			end
		end
	end

	if ~iscellstr(ncfiles)
		ncfiles = cellstr(ncfiles); 
	end

	% Give latout, lonout, and tout more intuitive names: 
	lat = latout; 
	lon = lonout; 
	t = tout; 

% -------------------------------------------------------------------------------------
function [lat,lon,t,P,T,S,pn] = argodata(ncfiles,varargin) 
% argodata downloads Argo profile data and imports it into your Matlab workspace.  
% 
%% Syntax 
% 
% [lat,lon,t,P,T,S,pn] = argodata(ncfiles) 
% [lat,lon,t,P,T,S,pn] = argodata(ncfiles,'keep') 
% 
%% Description
% 
% [lat,lon,t,P,T,S,pn] = argodata(ncfiles) downloads data from urls listed in the 
% cell array ncfiles (generated by the argofiles function).  If the data exist locally
% on your machine, ncfiles may be a cell array list of local filenames. Function 
% outputs lat, lon, t, and pn are geo coordinates, times in datenum format, and 
% platform numbers pn.  P, T, and S are cell arrays of pressure, temperature, and salinity.    
%
% [lat,lon,t,P,T,S,pn] = argodata(ncfiles,'keep') prevents automatic deletion of .nc files. 
% Select the 'keep' option if you want to access the raw NetCDF data or keep the files
% for later. 
% 
%% Example: All profiles in the Nino 3.4 box in December 2009   
% To plot all Argo data within the Nino 3.4 box collected in May 2009
% start by using the argofiles function to get urls of the data files: 
% 
%     latlim = [-5 5]; 
%     lonlim = [-170 -120]; 
%     t = [datenum('may 1, 2009') datenum('may 31, 2009')]; 
%     f = argofiles(latlim,lonlim,t,'pacific'); 
% 
% With a list of NetCDF file urls obtained above, use argodata to download the data and 
% import get measurements into the Matlab workspace.  Give it a minute to download all the data: 
% 
%     [lat,lon,t,P,T,S,pn] = argodata(f); 
% 
% Now plot all the profiles as a scattered T-S diagram with pressure as the
% color-scaled variable: 
% 
%     scatter(cell2mat(S(:)),cell2mat(T(:)),30,cell2mat(P(:)),'filled')
%     xlabel 'salinity'
%     ylabel 'temperature'
%     cb = colorbar; 
%     ylabel(cb,' pressure ')
%     set(cb,'ydir','reverse')
% 
%% Author Info: 
% The argodata function was written by Chad A. Greene of the University of
% Texas at Austin's Institute for Geophysics (UTIG), December 2015. 
% http://www.chadagreene.com.
%
% See also argofiles.
%
% Licenced under the BSD (FEX: 54503)
%
% JL: This function is modified to work in R13 and depends on Mirone code to read netCDF

	%% Input checks: 

	if (~isa(ncfiles, 'cell')),		ncfiles = {ncfiles};	end

	% Keep .nc files or automatically delete them? 
	if any(strcmpi(varargin,'keep'))
		keepfiles = true; 
	else
		keepfiles = false; 
	end

	%% Look for nc files: 

	% Are the nc files on the machine locally? 
	if exist(ncfiles{1},'file')==2
		islocal = true;
	else 
		islocal = false; 
	end

	% Preallocate: 
	ln = numel(ncfiles); 
	lat = zeros(ln,1) * NaN; 
	lon = zeros(ln,1) * NaN; 
	t = zeros(ln,1) * NaN; 
	P = cell(1,ln); 
	T = cell(1,ln); 
	S = cell(1,ln); 
	pn = zeros(ln,1) * NaN; 

	for k = 1:ln
		try 
			if islocal 
				filename = ncfiles{k}; 
			else
				if keepfiles
					[pato,fnm,ext] = fileparts(ncfiles{k}); 
					filename = strcat(fnm,ext); 
				else
					filename = 'argotmp.nc'; 
				end
				if (ispc),	dos(['wget "' ncfiles{k} '" -q --tries=2 --connect-timeout=4 --no-check-certificate -O "' filename '"']);
				else,		unix(['wget ''' ncfiles{k} ''' -q --tries=2 --connect-timeout=--no-check-certificate  -O "' filename '"']);
				end
			end
			x = nc_funs('read', filename,'latitude');		% Some files are fcked and have double values
			lat(k) = x(1);
			x = nc_funs('read', filename,'longitude'); 
			lon(k) = x(1);
			x = nc_funs('read', filename,'juld');
			t(k) = datenum(1950,1,x(1));
			P{k} = nc_funs('read', filename,'pres_adjusted');
			if (~all(isnan(P{k}(1,:))))
				T{k} = nc_funs('read', filename,'temp_adjusted');
				S{k} = nc_funs('read', filename,'psal_adjusted');
			else					% Then try these
				P{k} = nc_funs('read', filename,'pres');
				T{k} = nc_funs('read', filename,'temp');
				S{k} = nc_funs('read', filename,'psal');
			end
			pn(k) = str2double(nc_funs('read', filename,'platform_number')');
		catch
			disp(['Error accessing ',ncfiles{k}, '   ', lasterr])
		end
	end

	%% Clean up: 
	if exist('argotmp.nc','file')==2
		try
			delete('argotmp.nc')
		end
	end


% --- Creates and returns a handle to the GUI figure. 
function h1 = argo_floats_LayoutFcn()

	h1 = figure('Position',[748 837 230 150],...
	'Color', get(0,'factoryUicontrolBackgroundColor'),...
	'MenuBar','none',...
	'Name','Argo floats',...
	'NumberTitle','off',...
	'Resize','off',...
	'HandleVisibility','callback',...
	'Vis','off',...
	'Tag','figure1');

	uicontrol('Parent',h1, 'Position',[10 109.5 191 21],...
	'Callback',@argoFloats_uiCB,...
	'Style','edit',...
	'BackgroundColor',[1 1 1],...
	'TooltipString','Begining of period for the data display (rounded to months)',...
	'Tag','edit_startDate');

	uicontrol('Parent',h1, 'Position',[200 109.5 21 21],...
	'Callback',@argoFloats_uiCB,...
	'TooltipString','Call the callendar tool',...
	'Tag','push_startDate');

	uicontrol('Parent',h1, 'Position',[11 131 70 16],...
	'HorizontalAlignment','left',...
	'String','Start date',...
	'Style','text',...
	'FontSize',9);

	uicontrol('Parent',h1, 'Position',[10 60 191 21],...
	'Callback',@argoFloats_uiCB,...
	'Style','edit',...
	'BackgroundColor',[1 1 1],...
	'TooltipString','End of period for the data display (rounded to months)',...
	'Tag','edit_endDate');

	uicontrol('Parent',h1, 'Position',[200 60 21 21],...
	'Callback',@argoFloats_uiCB,...
	'TooltipString','Call the callendar tool',...
	'Tag','push_endDate');

	uicontrol('Parent',h1, 'Position',[11 81 70 16],...
	'HorizontalAlignment','left',...
	'String','End date',...
	'Style','text',...
	'FontSize',9);

	uicontrol('Parent',h1, 'Position',[11 19 75 19],...
	'String',{'atlantic' 'pacific' 'indian' 'all'}, ...
	'Style','popupmenu',...
	'Value',1,...
	'Tag','popup_basins');

	uicontrol('Parent',h1, 'Position',[11 37 40 16],...
	'HorizontalAlignment','left',...
	'String','Basins',...
	'Style','text');

	uicontrol('Parent',h1, 'Position',[140.5 10.5 81 21],...
	'Callback',@argoFloats_uiCB,...
	'String','OK',...
	'FontSize',9,...
	'FontWeight','bold',...
	'Tag','push_OK');

function argoFloats_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

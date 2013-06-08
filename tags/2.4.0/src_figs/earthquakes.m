function varargout = earthquakes(varargin)
% Helper window to plot seimicity

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

	if (isempty(varargin))
		errordlg('EARTHQUAKES: wrong number of input arguments.','Error'),	return
	end
 
	hObject = figure('Vis','off');
	earthquakes_LayoutFcn(hObject);
	handles = guihandles(hObject);
    
    handles.hMirFig = varargin{1};
    handMir = guidata(handles.hMirFig);
	if (handMir.no_file)
        errordlg('This option requires that you load a file first to serve as a background map.','ERROR')
        delete(hObject);    return
	end
	move2side(handMir.figure1, hObject,'right')

	handles.home_dir = handMir.home_dir;
	handles.last_dir = handMir.last_dir;
	handles.work_dir = handMir.work_dir;

    if (~handMir.is_projected && ~handMir.geog)
        errordlg('This operation is only possible for geographic data OR when the Map Projection is known','ERROR')
        delete(hObject);    return
    end
    handles.is_projected = handMir.is_projected;
    handles.defCoordsIn = handMir.defCoordsIn;
    handles.mironeAxes = handMir.axes1;
    zz = get(handles.mironeAxes,'XLim');
    handles.x_min = zz(1);    handles.x_max = zz(2);
    zz = get(handles.mironeAxes,'YLim');
    handles.y_min = zz(1);    handles.y_max = zz(2);

	% Add this figure handle to the carra?as list
	plugedWin = getappdata(handles.hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hMirFig,'dependentFigs',plugedWin);

    handles.path_data = handMir.path_data;

	if (nargin == 1)        % Use the default file shiped with Mirone
        set(handles.listbox_readFilter,'String','Not useful here','Enable','off')
        set(handles.push_externalFile,'Visible','off')
        handles.use_default_file = 1;
	else
		% Import icons
        load([handles.path_data 'mirone_icons.mat'],'Mfopen_ico');
        set(handles.push_externalFile,'CData',Mfopen_ico)
        % Fill the listbox fields with the currently available reading filters
        str = {'ISF formated catalog (ascii)'; 'Posit file'; 'lon,lat,mag,dep,yy,mm,dd,hh,mm,ss'; 'lon,lat,dep,mag,yy,mm,dd'};
        set(handles.listbox_readFilter,'String',str)
        handles.use_default_file = 0;
	end

	handles.got_userFile = 0;
	handles.have_mag_nans = 0;
	handles.have_dep_nans = 0;
    handles_fake.figure1 = handles.hMirFig;					% Create a fake handles only for
    handles_fake.axes1 = handles.mironeAxes;				% proj2proj_pts() satisfaction
    handles_fake.geog = handMir.geog;
    handles.handles_fake = handles_fake;

    if (handles.is_projected && handles.defCoordsIn > 0)        % We need a proj job here
        tmp = [handles.x_min handles.y_min; handles.x_max handles.y_max];
        lims = [handles.x_min handles.x_max handles.y_min handles.y_max 0 ];
        tmp = proj2proj_pts(handles.handles_fake,tmp, 'srcProj4','+proj=longlat', 'lim',lims);
        x_min = tmp(1,1);           x_max = tmp(2,1);
        y_min = tmp(1,2);           y_max = tmp(2,2);
        handles.lims_geogs = [x_min x_max y_min y_max];		% We'll need this if reading an external file
    else
        x_min = handles.x_min;      x_max = handles.x_max;
        y_min = handles.y_min;      y_max = handles.y_max;
    end

	if (handles.use_default_file)
		% Read the Mirone's default earthquakes file
		fid = fopen([handles.path_data 'quakes.dat'],'r');
		todos = fread(fid,'*char');
		[year mo day lat lon depth mag] = strread(todos,'%d %d %d %f %f %f %f');
		fclose(fid);    clear todos
        year_dec = dec_year(year,mo,day);
        
		% Get rid of events that are outside the map limits
		[lon,lat,indx,indy] = aux_funs('in_map_region', handMir, lon, lat, 0, [handles.x_min handles.x_max handles.y_min handles.y_max]);
		year(indx) = [];	mo(indx) = [];		day(indx) = [];
        depth(indx) = [];	mag(indx) = [];		year_dec(indx) = [];
		year(indy) = [];	mo(indy) = [];		day(indy) = [];
        depth(indy) = [];	mag(indy) = [];		year_dec(indy) = [];
		
		handles.def_StartYear = min(year);
		handles.def_EndYear = max(year);
		handles.def_StartMonth = min(mo);
		handles.def_EndMonth = max(mo);
		handles.def_StartDay = min(day);
		handles.def_EndDay = max(day);
		handles.def_MagMin = min(mag);
		handles.def_MagMax = max(mag);
		handles.def_DepthMin = min(depth);
		handles.def_DepthMax = max(depth);
		handles.default_dat = [lon lat depth mag];
		handles.default_date = [day mo year year_dec];
		set_lims(handles,'def')
	end

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) ------------------------------

	% Choose default command line output for earthquakes_export
	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% -------------------------------------------------------------------------------------------------
function edit_StartYear_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if isnan(xx),    set(hObject,'String','1900');   end

% -------------------------------------------------------------------------------------------------
function edit_StartMonth_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1 || xx > 12)
		set(hObject,'String','1')
	end

% -------------------------------------------------------------------------------------------------
function edit_StartDay_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1 || xx > 31)
		set(hObject,'String','1')
	end

% -------------------------------------------------------------------------------------------------
function edit_EndYear_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if isnan(xx),   set(hObject,'String','2010');   end

% -------------------------------------------------------------------------------------------------
function edit_EndMonth_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1 || xx > 12)
		set(hObject,'String','12')
	end

% -------------------------------------------------------------------------------------------------
function edit_EndDay_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1 || xx > 31)
		set(hObject,'String','31')
	end

% -------------------------------------------------------------------------------------------------
function edit_MagMin_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1 || xx > 10)
		set(hObject,'String','1')
	end

% -------------------------------------------------------------------------------------------------
function edit_MagMax_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1 || xx > 10)
		set(hObject,'String','10')
	end

% -------------------------------------------------------------------------------------------------
function edit_DepthMin_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)
		set(hObject,'String','0')
	end

% -------------------------------------------------------------------------------------------------
function edit_DepthMax_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx > 900),    set(hObject,'String','900');       end

% -----------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% Find out the time interval selected
	StartYear = str2double(get(handles.edit_StartYear,'String'));
	EndYear = str2double(get(handles.edit_EndYear,'String'));
	StartMonth = str2double(get(handles.edit_StartMonth,'String'));
	EndMonth = str2double(get(handles.edit_EndMonth,'String'));
	StartDay = str2double(get(handles.edit_StartDay,'String'));
	EndDay = str2double(get(handles.edit_EndDay,'String'));
	MagMin = str2double(get(handles.edit_MagMin,'String'));
	MagMax = str2double(get(handles.edit_MagMax,'String'));
	DepthMin = str2double(get(handles.edit_DepthMin,'String'));
	DepthMax = str2double(get(handles.edit_DepthMax,'String'));

	if (isnan(StartYear)),      StartYear = 1900;   end
	if (isnan(EndYear)),        EndYear = 2015;     end
	if (isnan(StartMonth)),     StartMonth = 1;     end
	if (isnan(EndMonth)),       EndMonth = 12;      end
	if (isnan(StartDay)),       StartDay = 1;       end
	if (isnan(EndDay)),         EndDay = 31;        end
	if (isnan(MagMin)),         MagMin = 1;         end
	if (isnan(MagMax)),         MagMax = 10;        end
	if (isnan(DepthMin)),       DepthMin = 0;       end
	if (isnan(DepthMax)),       DepthMax = 900;     end

	if (handles.use_default_file)
        lon = handles.default_dat(:,1);			lat = handles.default_dat(:,2);
        depth = handles.default_dat(:,3);		mag = handles.default_dat(:,4);
        year_dec = handles.default_date(:,4);
	elseif (handles.got_userFile)				% We have a user seismicity file
        lon = handles.external_dat(:,1);		lat = handles.external_dat(:,2);
        depth = handles.external_dat(:,3);		mag = handles.external_dat(:,4);
        year_dec = handles.external_date(:,2);
	else
        errordlg('Plot What? Your fears?','Chico Clever');  return;
	end
	
	lower_date = dec_year(StartYear,StartMonth,StartDay);
	upper_date = dec_year(EndYear,EndMonth,EndDay+0.999);   % 0.999 to use the entire current day
	ind = (year_dec < lower_date | year_dec > upper_date);
	year_dec(ind) = [];     lat(ind) = [];      lon(ind) = [];  depth(ind) = [];    mag(ind) = [];
	
	ind = (mag < MagMin | mag > MagMax);
	lat(ind) = [];  lon(ind) = [];  depth(ind) = [];    mag(ind) = [];  year_dec(ind) = [];
	ind = (depth < DepthMin | depth > DepthMax);
	lat(ind) = [];  lon(ind) = [];  depth(ind) = [];    mag(ind) = [];  year_dec(ind) = [];

    axes(handles.mironeAxes)       % Make Mirone axes active here
    
    if (handles.is_projected && handles.defCoordsIn > 0)        % We need a proj job here
        lims = [handles.x_min handles.x_max handles.y_min handles.y_max];
        tmp = proj2proj_pts(handles.handles_fake,[lon lat],  'srcProj4','+proj=longlat', 'lim',lims);
        lon = tmp(:,1);           lat = tmp(:,2);
    end

	if (min(mag) > 100)         % That's the case for hydrophone SL magnitudes
        mag_save = mag;         % Not converted to uint8 as well ???
	else
        mag_save = uint8(mag*10);
	end

% See if user only wants equal symbols (simple case)
if (~get(handles.check_magSlices,'Value') && ~get(handles.check_depSlices,'Value'))
	h_quakes = line('XData',lon,'YData',lat,'Parent',handles.mironeAxes,'Marker','o','LineStyle','none',...
          'MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',4,'Tag','Earthquakes');
    setappdata(h_quakes,'SeismicityTime',year_dec);         % Save events time
    setappdata(h_quakes,'SeismicityDepth',int16(depth*10)); % Save events depth
    setappdata(h_quakes,'SeismicityMag',mag_save);          % Save events magnitude
	draw_funs(h_quakes,'Earthquakes',[])
    return      % We are donne. Bye Bye
end

if (get(handles.check_magSlices,'Value'))					% We have a magnitude slice request
    cont = get(handles.popup_mag03,'String');	s(1) = str2double(cont{get(handles.popup_mag03,'Value')});
    cont = get(handles.popup_mag35,'String');	s(2) = str2double(cont{get(handles.popup_mag35,'Value')});
    cont = get(handles.popup_mag56,'String');	s(3) = str2double(cont{get(handles.popup_mag56,'Value')});
    cont = get(handles.popup_mag67,'String');	s(4) = str2double(cont{get(handles.popup_mag67,'Value')});
    cont = get(handles.popup_mag78,'String');	s(5) = str2double(cont{get(handles.popup_mag78,'Value')});
    cont = get(handles.popup_mag8,'String');	s(6) = str2double(cont{get(handles.popup_mag8,'Value')});
    id{1} = find(mag < 3);						id{2} = find(mag >= 3 & mag < 5);
    id{3} = find(mag >= 5 & mag < 6);			id{4} = find(mag >= 6 & mag < 7);
    id{5} = find(mag >= 7 & mag < 8);			id{6} = find(mag >= 8);
    j = 1;
    for (k=1:length(id))
        if (~isempty(id{k}))
            data_s{j} = [lon(id{k}) lat(id{k}) depth(id{k})];
            grand(j) = s(k);
            depth_s{j} = int16( depth(id{k})*10 );
            mag_s{j} = mag_save(id{k});
            year_s{j} = year_dec(id{k});
            j = j + 1;
        end
    end
    if (handles.have_mag_nans && get(handles.check_allMagnitudes,'Value'))
        id = isnan(mag);
        data_s{1} = [data_s{:}; lon(id) lat(id) depth(id)];
    end
    clear id;
end

if (get(handles.check_depSlices,'Value'))    % We have a depth slice request
    val(1) = get(handles.popup_dep0_33,'Value');    val(2) = get(handles.popup_dep33_70,'Value');
    val(3) = get(handles.popup_dep70_150,'Value');  val(4) = get(handles.popup_dep150_300,'Value');
    val(5) = get(handles.popup_dep300,'Value');
    id{1} = find(depth < 33);                   id{2} = find(depth >= 33 & depth < 70);
    id{3} = find(depth >= 70 & depth < 150);    id{4} = find(depth >= 150 & depth < 300);
    id{5} = find(depth >= 300);
    j = 1;      cor_str = get(handles.popup_dep0_33,'String');
    for (k=1:length(id))
        if (~isempty(id{k}))
            data_d{j} = [lon(id{k}) lat(id{k})];
            color(j) = cor_str{val(k)}(1);
            depth_d{j} = int16( depth(id{k})*10 );
            mag_d{j} = mag_save(id{k});
            year_d{j} = year_dec(id{k});
            j = j + 1;
        end
    end    
end

if (get(handles.check_magSlices,'Value') && ~get(handles.check_depSlices,'Value'))		% Mag slices
    for (k = 1:length(data_s))
		h_quakes = line('XData',data_s{k}(:,1),'YData',data_s{k}(:,2),'LineStyle','none','Marker','o',...
            'MarkerFaceColor','r','Parent',handles.mironeAxes,'MarkerEdgeColor','k',...
            'MarkerSize',grand(k),'Tag','Earthquakes');
        setappdata(h_quakes,'SeismicityDepth',depth_s{k});		% Save events depth
        setappdata(h_quakes,'SeismicityMag',mag_s{k});			% Save events magnitude
        setappdata(h_quakes,'SeismicityTime',year_s{k});		% Save events time
        draw_funs(h_quakes,'Earthquakes',[])
    end
elseif (~get(handles.check_magSlices,'Value') && get(handles.check_depSlices,'Value')) % Depth slices
    for (k = 1:length(data_d))
		h_quakes = line('XData',data_d{k}(:,1),'YData',data_d{k}(:,2),'LineStyle','none','Marker','o',...
            'MarkerFaceColor',color(k),'Parent',handles.mironeAxes,'MarkerEdgeColor','k',...
            'MarkerSize',5,'Tag','Earthquakes');
        setappdata(h_quakes,'SeismicityDepth',depth_d{k});		% Save events depth
        setappdata(h_quakes,'SeismicityMag',mag_d{k});			% Save events magnitude
        setappdata(h_quakes,'SeismicityTime',year_d{k});		% Save events time
        draw_funs(h_quakes,'Earthquakes',[])
    end
else        % Both magnitude and depth slices
    for (k = 1:length(data_s))    %
        id{1} = find(data_s{k}(:,3) < 33);
        id{2} = find(data_s{k}(:,3) >= 33 & data_s{k}(:,3) < 70);
        id{3} = find(data_s{k}(:,3) >= 70 & data_s{k}(:,3) < 150);
        id{4} = find(data_s{k}(:,3) >= 150 & data_s{k}(:,3) < 300);
        id{5} = find(data_s{k}(:,3) >= 300);
        for (m=1:5)
            if (isempty(id{m})),     continue;      end
		    h_quakes = line('XData',data_s{k}(id{m},1),'YData',data_s{k}(id{m},2),'LineStyle','none',...
                'Marker','o','MarkerFaceColor',color(m),'Parent',handles.mironeAxes,...
                'MarkerEdgeColor','k','MarkerSize',grand(k),'Tag','Earthquakes');
            %setappdata(h_quakes,'SeismicityDepth',data_s{k}(id{m},3));		% Save events depth
            setappdata(h_quakes,'SeismicityDepth',depth_s{k}(id{m}));		% Save events depth
            setappdata( h_quakes,'SeismicityMag',mag_s{k}(id{m}) );			% Save events magnitude
            setappdata( h_quakes,'SeismicityTime',year_s{k}(id{m}) );		% Save events time
            draw_funs(h_quakes,'Earthquakes',[])
        end
    end
end

% -----------------------------------------------------------------------------
function set_lims(handles,opt)
	if (strcmp(opt,'def'))      % Set the limits corresponding to the default file
		set(handles.edit_StartYear,'String',num2str(handles.def_StartYear))
		set(handles.edit_EndYear,'String',num2str(handles.def_EndYear))
		set(handles.edit_StartMonth,'String',num2str(handles.def_StartMonth))
		set(handles.edit_EndMonth,'String',num2str(handles.def_EndMonth))
		set(handles.edit_StartDay,'String',num2str(handles.def_StartMonth))
		set(handles.edit_EndDay,'String',num2str(handles.def_EndDay))
		set(handles.edit_MagMin,'String',num2str(handles.def_MagMin))
		set(handles.edit_MagMax,'String',num2str(handles.def_MagMax))
		set(handles.edit_DepthMin,'String',num2str(handles.def_DepthMin))
		set(handles.edit_DepthMax,'String',num2str(handles.def_DepthMax))
	elseif (strcmp(opt,'usr'))  % Set the limits corresponding to the user's file
		set(handles.edit_StartYear,'String',num2str(handles.usr_StartYear))
		set(handles.edit_EndYear,'String',num2str(handles.usr_EndYear))
		set(handles.edit_StartMonth,'String',num2str(handles.usr_StartMonth))
		set(handles.edit_EndMonth,'String',num2str(handles.usr_EndMonth))
		set(handles.edit_StartDay,'String',num2str(handles.usr_StartMonth))
		set(handles.edit_EndDay,'String',num2str(handles.usr_EndDay))
		set(handles.edit_MagMin,'String',num2str(handles.usr_MagMin))
		set(handles.edit_MagMax,'String',num2str(handles.usr_MagMax))
		set(handles.edit_DepthMin,'String',num2str(handles.usr_DepthMin))
		set(handles.edit_DepthMax,'String',num2str(handles.usr_DepthMax))
	else						% Set these fields to empty
		set(handles.edit_StartYear,'String','');	set(handles.edit_EndYear,'String','')
		set(handles.edit_StartMonth,'String','');	set(handles.edit_EndMonth,'String','')
		set(handles.edit_StartDay,'String',''); 	set(handles.edit_EndDay,'String','')
		set(handles.edit_MagMin,'String','');   	set(handles.edit_MagMax,'String','')
		set(handles.edit_DepthMin,'String',''); 	set(handles.edit_DepthMax,'String','')
	end

% -----------------------------------------------------------------------------
function push_externalFile_CB(hObject, handles)
% OK. Now read the earthquakes_export file and retain only the requested interval
	item = get(handles.listbox_readFilter,'Value');     % Get the reading filter number
	switch item
		case 1
			str1 = {'*.isf;*.ISF', 'Data files (*.isf,*.ISF)';'*.*', 'All Files (*.*)'};
			filtro = 1;
		case 2
			str1 = {'*.posit;*.POSIT', 'Data files (*.posit,*.POSIT)';'*.*', 'All Files (*.*)'};
			filtro = 2;
		case 3
			str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
			filtro = 3;
		case 4
			str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
			filtro = 4;
	end

	% Get file name
	[FileName,PathName,handles] = put_or_get_file(handles, str1,'Select earhquakes file', 'get');
	if isequal(FileName,0),		return,		end
	fname = [PathName,FileName];

	try
		set(handles.figure1,'Pointer','watch')
		if (filtro == 1)        % Read a ISF formated catalog
			if (handles.is_projected)               % Image is projected, we need this
				opt_R = sprintf('-R%f/%f/%f/%f',handles.lims_geogs(1),handles.lims_geogs(2), ...
						handles.lims_geogs(3),handles.lims_geogs(4));
			else
				opt_R = sprintf('-R%f/%f/%f/%f',handles.x_min,handles.x_max,handles.y_min,handles.y_max);
			end
			[out_d,out_i] = read_isf(fname,opt_R);
			if (isempty(out_d)),    return;     end     % Nothing inside region
			lon = out_d(1,:)';      lat = out_d(2,:)';
			depth = out_d(3,:)';    mag = out_d(4,:)';      clear out_d;
			year = double(out_i(1,:)');     mo = out_i(2,:)';
			day = out_i(3,:)';      hh = out_i(4,:)';       clear out_i;
			year_dec = dec_year(year,double(mo),double(day),double(hh));
		elseif (filtro == 2 || filtro == 3 || filtro == 4)      % Read a lon,lat,dep,mag,yy,mm,dd file (3,4) or posit file (2)
			fid = fopen(fname,'r');
			if (fid < 0)
				errordlg(['Could not open file: ' fname],'Error'),		return
			end
			todos = fread(fid,'*char');
			if (filtro == 2)                    % posit file
				try
					[year julio d_h d_m d1 lat lon d1 d1 d1 mag d1] = strread(todos,'%d %d %d %d %d %f %f %f %f %f %f %f');
					d_h = d_h + d_m / 60;
				catch
					[tempo d1 d1 lat lon d1 d1 d1 mag d1] = strread(todos,'%s %d %s %f %f %f %f %f %f %f');
					d1 = cell2mat(tempo);
					year = str2num(d1(:,1:4));
					julio = str2num(d1(:,5:7));
					d_h = str2num(d1(:,8:9));
					d_m = str2num(d1(:,10:11));
				end
				clear d_m d1;
				year_dec = year + (julio - 1 + d_h / 24) ./ (365 + isleapyear(year));      % decimal year up to minuts.
				[mo,day] = jd2monday(julio,year);
				depth = zeros(length(year),1);
			elseif (filtro == 3)                % lon,lat,mag,dep,yy,mm,dd,hh,mm,ss
				[lon lat mag depth year mo day hh mm ss] = strread(todos,'%f %f %f %f %d %d %d %d %d %d');
				year_dec = dec_year(year,mo,day,hh,mm,ss);  clear hh mm ss;
			else                                % lon,lat,dep,mag,yy,mm,dd file
				[lon lat depth mag year mo day] = strread(todos,'%f %f %f %f %d %d %d');
				year_dec = dec_year(year,mo,day);
			end
			fclose(fid);	clear todos

			if (handles.is_projected && handles.defCoordsIn > 0)        % Image is projected, we need to use this
				x_min = handles.lims_geogs(1);      x_max = handles.lims_geogs(2);
				y_min = handles.lims_geogs(3);      y_max = handles.lims_geogs(4);
			else
				x_min = handles.x_min;      x_max = handles.x_max;
				y_min = handles.y_min;      y_max = handles.y_max;
			end

			% Get rid of events that are outside the map limits
			handMir = guidata(handles.hMirFig);
			[lon,lat,indx,indy] = aux_funs('in_map_region', handMir, lon, lat, 0, [handles.x_min handles.x_max handles.y_min handles.y_max]);
			year(indx) = [];	%mo(indx) = [];		day(indx) = [];
			depth(indx) = [];	mag(indx) = [];		year_dec(indx) = [];
			year(indy) = [];	%mo(indy) = [];		day(indy) = [];
			depth(indy) = [];	mag(indy) = [];		year_dec(indy) = [];
			if (isempty(lon))			% Nothing inside region
				set(handles.figure1,'Pointer','arrow'),		return
			end
		end
		handles.got_userFile = 1;

		mag(mag < 0) = 0;       % Take care of mags < 0 .They are very likely false
		handles.external_dat  = [lon lat depth mag];
		%handles.external_date = [day mo year year_dec];
		handles.external_date = [year year_dec];

		handles.usr_DepthMin   = min(depth);    handles.usr_DepthMax = max(depth);
		handles.usr_MagMin     = min(mag);      handles.usr_MagMax = max(mag);
		handles.usr_StartYear  = double(min(year));     handles.usr_EndYear = double(max(year));

		% Find start and end Month & Day
		tmp = min(year_dec);    tmp = tmp - fix(tmp);
		jd0 = fix(tmp * (365 + isleapyear(handles.usr_StartYear))) + 1;
		tmp = max(year_dec);    tmp = tmp - fix(tmp);
		jd1 = fix(tmp * (365 + isleapyear(handles.usr_EndYear))) + 1;
		[handles.usr_StartMonth,handles.usr_StartDay] = jd2monday(jd0,handles.usr_StartYear);
		[handles.usr_EndMonth,handles.usr_EndDay] = jd2monday(jd1,handles.usr_EndYear);

		if (filtro == 3 || filtro == 4)
			handles.have_mag_nans = any(isnan(mag));
			handles.have_dep_nans = any(isnan(depth));
		else			% On the ISF catalogs I replaced no data values by 0
			handles.have_mag_nans = 0;
			handles.have_dep_nans = 0;
		end
		set(handles.figure1,'Pointer','arrow')
	catch				% In case of error, set the pointer back to "normal" 
		set(handles.figure1,'Pointer','arrow')
		msg{1} = 'An error occured while reading file. The error message was:';
		msg{2} = '';
		msg{3} = lasterr;
		warndlg(msg,'Warning')
		return
	end

	set(handles.edit_StartYear,'String',int2str_m(handles.usr_StartYear))
	set(handles.edit_EndYear,'String',int2str_m(handles.usr_EndYear))
	set(handles.edit_StartMonth,'String',int2str_m(handles.usr_StartMonth))
	set(handles.edit_EndMonth,'String',int2str_m(handles.usr_EndMonth))
	set(handles.edit_StartDay,'String',int2str_m(handles.usr_StartDay))
	set(handles.edit_EndDay,'String',int2str_m(handles.usr_EndDay))
	set(handles.edit_MagMin,'String',num2str(handles.usr_MagMin))
	set(handles.edit_MagMax,'String',num2str(handles.usr_MagMax))
	set(handles.edit_DepthMin,'String',num2str(handles.usr_DepthMin))
	set(handles.edit_DepthMax,'String',num2str(handles.usr_DepthMax))
	if (handles.have_mag_nans),  set(handles.check_allMagnitudes,'Enable','on');  end
	if (handles.have_dep_nans),  set(handles.check_allDepths,'Enable','on');  end

	guidata(hObject,handles)

% -----------------------------------------------------------------------------
function check_allMagnitudes_CB(hObject, handles)
	if (get(hObject,'Value') && ~get(handles.check_magSlices,'Value'))
		set(hObject,'Value',0)
	end

% -----------------------------------------------------------------------------
function check_allDepths_CB(hObject, handles)
	if (get(hObject,'Value') && ~get(handles.check_depSlices,'Value'))
		set(hObject,'Value',0)
	end

% -----------------------------------------------------------------------------
function check_magSlices_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.popup_mag03,'Enable','on');    set(handles.popup_mag35,'Enable','on')
		set(handles.popup_mag56,'Enable','on');    set(handles.popup_mag67,'Enable','on')
		set(handles.popup_mag78,'Enable','on');    set(handles.popup_mag8, 'Enable','on')
	else
		set(handles.popup_mag03,'Enable','off');   set(handles.popup_mag35,'Enable','off')
		set(handles.popup_mag56,'Enable','off');   set(handles.popup_mag67,'Enable','off')
		set(handles.popup_mag78,'Enable','off');   set(handles.popup_mag8, 'Enable','off')
	end

% -----------------------------------------------------------------------------
function check_depSlices_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.popup_dep0_33,'Enable','on');       set(handles.popup_dep33_70,'Enable','on')
		set(handles.popup_dep70_150,'Enable','on');     set(handles.popup_dep150_300,'Enable','on')
		set(handles.popup_dep300,'Enable','on')
	else
		set(handles.popup_dep0_33,'Enable','off');      set(handles.popup_dep33_70,'Enable','off')
		set(handles.popup_dep70_150,'Enable','off');    set(handles.popup_dep150_300,'Enable','off')
		set(handles.popup_dep300,'Enable','off')
	end

% -----------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles = guidata(hObject);
        delete(handles.figure1);
	end

% -----------------------------------------------------------------------------
function [month, day] = jd2monday(jday,year)
%JD2MONDAY Julian day number to Julian calendar date.
%
%   [MONTH, DAY] = JD2MONDAY(JDAY,YEAR) returns the
%   Julian calendar date (month, day)
%   corresponding to the Julian day JDAY.

%   Author:      Peter J. Acklam
%   Hacked to work with the "fake" JD that start at 1 at 1fst January of each year.

	t = ( ~rem(year, 4) & rem(year, 100) ) | ~rem(year, 400);       % Check for leap-years
	tt = (~t & jday > 59);
	jday(tt) = jday(tt) + 1;            % Trick to make this algo work also for not leap-years

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
function s = int2str_m(x)
%INT2STR Convert integer to string.
%   Round the scalar X to integer and converts the result into a string.
%   Return NaN and Inf elements as strings 'NaN' and 'Inf', respectively.

	if (~isa(x,'double'))
		x = double(x);
	end

	x = round(real(x));
	% handle special case of single infinite or NaN element
	s = sprintf('%.1f',x);
	if (~strcmp(s, '-Inf') && ~strcmp(s, 'Inf') && ~strcmp(s, 'NaN'))
		s(end-1:end) = [];
	end

% --- Creates and returns a handle to the GUI figure. 
function earthquakes_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Plot seismicity',...
'NumberTitle','off',...
'Position',[520 288 391 512],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[10 49 371 141],'Style','frame');
uicontrol('Parent',h1,'Position',[10 199 371 141],'Style','frame');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[63 441 251 61],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_readFilter');

uicontrol('Parent',h1,...
'Call',@earthquakes_uiCB,...
'FontWeight','bold',...
'Position',[312 461 23 23],...
'Tooltip','Browse for wanted file',...
'Tag','push_externalFile');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@earthquakes_uiCB,...
'Position',[63 405 47 21],...
'Style','edit','Tag','edit_StartYear');

uicontrol('Parent',h1,'Position',[23 395 36 33],'String',{'Start'; 'year'},'Style','text','Tag','text1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@earthquakes_uiCB,...
'Position',[187 405 47 21],...
'Style','edit','Tag','edit_StartMonth');

uicontrol('Parent',h1,'Position',[143 399 41 30],'String',{'Start'; 'month'},'Style','text','Tag','text2');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@earthquakes_uiCB,...
'Position',[313 405 47 21],...
'Style','edit','Tag','edit_StartDay');

uicontrol('Parent',h1,'Position',[277 399 36 30],'String',{'Start'; 'day'},'Style','text','Tag','text3');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@earthquakes_uiCB,...
'Position',[63 370 47 21],...
'Style','edit','Tag','edit_EndYear');

uicontrol('Parent',h1,...
'Position',[23 361 36 30],...
'String',{'End'; 'year'},...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@earthquakes_uiCB,...
'Position',[187 370 47 21],...
'Style','edit','Tag','edit_EndMonth');

uicontrol('Parent',h1,...
'Position',[143 361 41 30],...
'String',{'End'; 'month'},...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@earthquakes_uiCB,...
'Position',[313 370 47 21],...
'Style','edit','Tag','edit_EndDay');

uicontrol('Parent',h1,...
'Position',[283 365 27 30],...
'String',{'End'; 'day'},...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@earthquakes_uiCB,...
'Position',[71 307 47 21],...
'Style','edit','Tag','edit_MagMin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@earthquakes_uiCB,...
'Position',[195 307 47 21],...
'Style','edit','Tag','edit_MagMax');

uicontrol('Parent',h1, 'Position',[17 300 52 30],...
'String',{'Minimum'; 'magnitude'},...
'Style','text');

uicontrol('Parent',h1, 'Position',[130 301 60 30],...
'String',{'Maximum'; 'magnitude'},...
'Style','text','HorizontalAlignment','right');

uicontrol('Parent',h1, 'Position',[268 310 110 15],...
'Call',@earthquakes_uiCB,...
'Enable','off',...
'String','All magnitudes',...
'Style','checkbox',...
'Tooltip','Use all mgnitudes - Including unknow magnitudes',...
'Tag','check_allMagnitudes');

uicontrol('Parent',h1, 'Position',[18 264 260 15],...
'Call',@earthquakes_uiCB,...
'String','Use different sizes for magnitude intervals',...
'Style','checkbox',...
'Tooltip','Destinguish the seismic magnitudes by size',...
'Tag','check_magSlices');

uicontrol('Parent',h1, 'Position',[10 145 51 30],...
'String',{'Minimum'; 'depth' },...
'Style','text','HorizontalAlignment','right');

uicontrol('Parent',h1, 'Position',[128 144 55 30],...
'String',{'Maximum'; 'depth'},...
'Style','text','HorizontalAlignment','right');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[18 207 52 22],...
'String',{'3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; '11'; '12'; '13'; '14'; '15' },...
'Style','popupmenu',...
'Tooltip','Symbol size for this interval',...
'Value',2,...
'Tag','popup_mag03');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[78 207 52 22],...
'String',{'3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; '11'; '12'; '13'; '14'; '15' },...
'Style','popupmenu',...
'Tooltip','Symbol size for this interval',...
'Value',4,...
'Tag','popup_mag35');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[138 207 52 22],...
'String',{'3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; '11'; '12'; '13'; '14'; '15' },...
'Style','popupmenu',...
'Tooltip','Symbol size for this interval',...
'Value',6,...
'Tag','popup_mag56');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[198 207 52 22],...
'String',{'3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; '11'; '12'; '13'; '14'; '15' },...
'Style','popupmenu',...
'Tooltip','Symbol size for this interval',...
'Value',8,...
'Tag','popup_mag67');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[258 207 52 22],...
'String',{'3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; '11'; '12'; '13'; '14'; '15' },...
'Style','popupmenu',...
'Tooltip','Symbol size for this interval',...
'Value',10,...
'Tag','popup_mag78');

uicontrol('Parent',h1,'FontSize',10,'Position',[20 230 31 16], 'String','0-3','Style','text','Tag','text11');
uicontrol('Parent',h1,'FontSize',10,'Position',[81 230 31 16], 'String','3-5','Style','text','Tag','text12');
uicontrol('Parent',h1,'FontSize',10,'Position',[139 230 31 16],'String','5-6','Style','text','Tag','text13');
uicontrol('Parent',h1,'FontSize',10,'Position',[200 230 31 16],'String','6-7','Style','text','Tag','text14');
uicontrol('Parent',h1,'FontSize',10,'Position',[260 230 31 16],'String','7-8','Style','text','Tag','text15');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[318 207 52 22],...
'String',{'3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; '11'; '12'; '13'; '14'; '15' },...
'Style','popupmenu',...
'Tooltip','Symbol size for this interval',...
'Value',13,...
'Tag','popup_mag8');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@earthquakes_uiCB,...
'Position',[63 150 47 21],...
'String','0',...
'Style','edit',...
'Tag','edit_DepthMin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@earthquakes_uiCB,...
'Position',[187 150 47 21],...
'Style','edit',...
'Tag','edit_DepthMax');

uicontrol('Parent',h1,...
'Call',@earthquakes_uiCB,...
'Enable','off',...
'Position',[268 153 88 15],...
'String','All depths',...
'Style','checkbox',...
'Tooltip','Use all depths - Including unknow depths',...
'Tag','check_allDepths');

uicontrol('Parent',h1,'FontSize',10,'Position',[320 230 31 16],'String','> 8','Style','text','Tag','text16');

uicontrol('Parent',h1,...
'Call',@earthquakes_uiCB,...
'Position',[19 115 230 15],...
'String','Use different colors for depth intervals',...
'Style','checkbox',...
'Tooltip','Destinguish the epicenter depths by color',...
'Tag','check_depSlices');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[18 58 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'Tooltip','Symbol color for this depth interval',...
'Value',1,...
'Tag','popup_dep0_33');

uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[21 81 47 16],...
'String','0-33 km',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[91 58 62 22],...
'String',{'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'Tooltip','Symbol color for this depth interval',...
'Value',2,...
'Tag','popup_dep33_70');

uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[94 81 54 16],...
'String','33-70 km',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[165 58 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'Tooltip','Symbol color for this depth interval',...
'Value',3,...
'Tag','popup_dep70_150');

uicontrol('Parent',h1,'FontSize',10,'Position',[166 81 61 16],'String','70-150 km','Style','text','Tag','text19');
uicontrol('Parent',h1,'FontSize',10,'Position',[312 81 55 16],'String','> 300 km','Style','text','Tag','text20');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[238 58 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'Tooltip','Symbol color for this depth interval',...
'Value',4,...
'Tag','popup_dep150_300');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[310 58 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'Tooltip','Symbol color for this depth interval',...
'Value',5,...
'Tag','popup_dep300');

uicontrol('Parent',h1,'FontSize',10,'Position',[235 81 68 16],'String','150-300 km','Style','text','Tag','text21');

uicontrol('Parent',h1,...
'Call',@earthquakes_uiCB,...
'Position',[315 10 66 21],...
'String','OK',...
'Tag','push_OK');

function earthquakes_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

function varargout = coordinate_system(varargin)
% Helper window to select a projection string for GMT's -J

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
 
	hObject = figure('Vis','off');
	coordinate_system_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center');

% Declare the currently existing groups
if (numel(varargin) == 3)
	handles.proj_groups = {'Longitude / Latitude';
            'Linear Projection';...
            'All GMT Projections'};
else
	handles.proj_groups = {'Longitude / Latitude';
            'Universal Transversal Mercator';...
            'World Projections';...
            'All GMT Projections';...
            'Portuguese Mess'};
end

handles.bg_color = [.764,.603,.603];
handles.lon0_txt = 'Lon_of_Central_Meridian';
handles.lat0_txt = 'Lat_of_True_Scale';
handles.latOrig_txt = 'Lat_of_Origin';
handles.lon1_txt = 'Lon_of_2_Point';
handles.lat1_txt = 'Lat_of_2_Point';
handles.lonp_txt = 'Lon_of_Pole';
handles.latp_txt = 'Lat_of_Pole';
handles.azimuth_txt = 'Azimuth';
handles.proj_name = [];
handles.sys_name = [];
handles.fig_hand = hObject;
handles.map_scale_factor = [];  % for example = 1 for TM, = 0.9996 for UTM & Polar
handles.system_FE_FN = [];      % For false eastings and false northings
handles.all_datums = varargin{2};

% Build the UTM zones list
% First Northern hemisphere
u = cell(120,1);
for (i = 1:30) 
	u{i}=sprintf(['Zone ' sprintf('%.2d',i) 'N (' sprintf('%d',180 -(i-1)*6) 'W to ', sprintf('%d',180 -i*6) 'W)']);
end
for (i = 1:30)
	u{i+30}=sprintf(['Zone ' sprintf('%.2d',i+30) 'N (' sprintf('%d',(i-1)*6) 'E to ', sprintf('%d',i*6) 'E)']);
end
% And now the Southern hemisphere
for (i = 1:30) u{i+60}=sprintf(['Zone ' sprintf('%.2d',i) 'S (' sprintf('%d',180 -(i-1)*6) 'W to ', num2str(180 -i*6) 'W)']);
end
for (i = 1:30)
	u{i+90}=sprintf(['Zone ' sprintf('%.2d',i+30) 'S (' sprintf('%d',(i-1)*6) 'E to ', sprintf('%d',i*6) 'E)']);
end
handles.UTMzones = u;

set(handles.popup_Group,'String',handles.proj_groups)       % Set group
set(handles.popup_Group,'Value',varargin{1}.group_val)      % Set group old value
handles.group_val = varargin{1}.group_val;

set(handles.popup_System,'Value',varargin{1}.system_val)                        % Set old system value
set(handles.popup_CilindricalProjections,'Value',varargin{1}.cilindrical_val)   % Set old cilind value
set(handles.popup_AzimuthalProjections,'Value',varargin{1}.azimuthal_val)       % Set old azim value
set(handles.popup_ConicalProjections,'Value',varargin{1}.conic_val)             % Set old conic value
set(handles.popup_MiscelaneousProjections,'Value',varargin{1}.miscelaneous_val) % Set old misc value
set(handles.popup_Datum,'String',varargin{2}(:,1))                              % Set Datum
set(handles.popup_Datum,'Value',varargin{1}.datum_val)                          % Set Datum old value

switch handles.proj_groups{varargin{1}.group_val}           % Set apropriate system for the current group
    case 'Longitude / Latitude'
        set(handles.popup_System,'String','Longitude / Latitude')
        set_enable(handles, 0);        set_enable(handles, 50);
        handles.nParameters = 0;       handles.projection{1} = [];
        handles.proj_name = 'Longitude / Latitude';
    case 'Universal Transversal Mercator'
        set(handles.popup_System,'String',handles.UTMzones)
        set_enable(handles, 0);        set_enable(handles, 50);
        handles.proj_name = 'UTM';     handles.nParameters = 0;
        if (varargin{1}.system_val <= 60)
            handles.projection{1} = ['-Ju' num2str(varargin{1}.system_val) '/1'];
        else
            zone = varargin{1}.system_val - 60;
            handles.projection{1} = ['-Ju-' num2str(zone) '/1'];
        end
        handles.map_scale_factor = 0.9996;
    case 'World Projections'
        caso.type = 'input';  caso.val = varargin{1}.system_val;
        handles = set_WorldProj(handles,caso);
    case 'All GMT Projections'
        set_enable(handles, 51);    handles.nParameters = [];
        set(handles.popup_System,'Value',1,'String',{'Use the right side boxes'})
        handles.proj_name = varargin{1}.ProjName;
        if (~isempty(varargin{1}.projection))
            %tmp = varargin{1}.projection(1:4);
            tmp = '-Ju29';
            if (strcmp(tmp(4),'-'))
                str_proj = varargin{1}.projection(1:3);
            elseif (strcmp(tmp(3:4),'m0'))             % Case of -Jm0/0/1
                str_proj = varargin{1}.projection(1:4);
            elseif (strmatch(tmp(4),{'0';'1';'2';'3';'4';'5';'6';'7';'8';'9'}))
                str_proj = varargin{1}.projection(1:3);
            else
                str_proj = varargin{1}.projection(1:4);
            end
            handles.projection{1} = str_proj;
        end
    case 'Portuguese Mess'
        caso.type = 'input';  caso.val = varargin{1}.system_val;
        handles = set_PtMess(handles,caso);
        handles.proj_name = varargin{1}.ProjName;
    case 'Linear Projection'        % Used only in write_script GUI
        set(handles.popup_System,'String','Linear Projection')
        set_enable(handles, 0);        set_enable(handles, 50);
        handles.nParameters = 0;       handles.projection{1} = '-Ju';
        handles.proj_name = 'Linear Projection';
end

	% If old projection parameters exist, apply them
	try
		erro = varargin{1}.ProjParameterValue;     erro = 0;
	catch
		erro = 1;
	end
	if (~erro && ~isempty(varargin{1}.ProjParameterValue))
		n_params = length(fieldnames(varargin{1}.ProjParameterValue));
		set_enable(handles, n_params);
		for (i=1:n_params)
			switch i
				case 1
					handles.nParameters = 1;
					str_par = strtok(varargin{1}.ProjParameterValue.p1,'/');
					set(handles.edit_ProjParameterValue_1,'String',str_par)
					handles.projection{2} = varargin{1}.ProjParameterValue.p1;
				case 2
					handles.nParameters = 2;
					str_par = strtok(varargin{1}.ProjParameterValue.p2,'/');
					set(handles.edit_ProjParameterValue_2,'String',str_par)
					handles.projection{3} = varargin{1}.ProjParameterValue.p2;
				case 3
					handles.nParameters = 3;
					str_par = strtok(varargin{1}.ProjParameterValue.p3,'/');
					set(handles.edit_ProjParameterValue_3,'String',str_par)
					handles.projection{4} = varargin{1}.ProjParameterValue.p3;
				case 4
					handles.nParameters = 4;
					str_par = strtok(varargin{1}.ProjParameterValue.p4,'/');
					set(handles.edit_ProjParameterValue_4,'String',str_par)
					handles.projection{5} = varargin{1}.ProjParameterValue.p4;
			end
		end
	else
		guidata(hObject, handles);
		set_enable(handles, 0); 
	end

	% Not very elegant but it will have to do until a better solution is found
	if (strcmp(handles.proj_groups{varargin{1}.group_val},'Portuguese Mess'))
		set_enable(handles, 0);
	end

	% If an old projection string was transmited
	if (isempty(varargin{1}.projection))
		for (i=1:5),    handles.projection{i} = [];     end
	end

	% Make sure projection name has something
	if (isempty(handles.proj_name))
		handles.proj_name = 'FDS Sei la';
	end

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) ------------------------------

	%------------ Write info text about datum parameters -------------------
	num = varargin{1}.datum_val;
	info.DX = varargin{2}{num,4};		info.DY = varargin{2}{num,5};
	info.DZ = varargin{2}{num,6};		info.ellipsoide = varargin{2}{num,2};
	datum_info(handles, info)

	% Choose default command line output for coordinate_system_export
	handles.output = hObject;
	guidata(hObject, handles);

	set(hObject,'Visible','on');
	% UIWAIT makes coordinate_system_export wait for user response (see UIRESUME)
	uiwait(handles.figure1);

	handles = guidata(hObject);
	varargout{1} = handles.output;
	delete(hObject)

%-------------------------------------------------------------------------------------
function popup_Group_CB(hObject, handles)
contents = get(hObject,'String');   val = get(hObject,'Value');   group = contents{val};
switch group
    case 'Longitude / Latitude'
        set(handles.popup_System,'Value',1)
        set_enable(handles,0);      set_enable(handles, 50);      handles.nParameters = 0;
        handles.projection{1} = [];
        set(handles.popup_System,'String','Longitude / Latitude')
        handles.proj_name = 'Longitude / Latitude';
        handles.sys_name = 'Geographic';
        set(handles.popup_Datum,'Value',221);       info.DX = handles.all_datums{221,4};
        info.DY = handles.all_datums{221,5};        info.DZ = handles.all_datums{221,6};
        info.ellipsoide = handles.all_datums{221,2};        datum_info(handles, info)
    case 'Universal Transversal Mercator'
        set_enable(handles,0);      set_enable(handles, 50);      handles.nParameters = 0;
        set(handles.popup_System,'String',handles.UTMzones)
        handles.projection{1} = '-Ju1/1';     % For never being empty, but I must do better
        handles.proj_name = 'UTM';
        handles.sys_name = 'Zone 01N (180W to 174W)';
        handles.map_scale_factor = 0.9996;
        set(handles.popup_Datum,'Value',221);       info.DX = handles.all_datums{221,4};
        info.DY = handles.all_datums{221,5};        info.DZ = handles.all_datums{221,6};
        info.ellipsoide = handles.all_datums{221,2};        datum_info(handles, info)
    case 'World Projections'
        set(handles.popup_System,'Value',1,'String',{'Mercator';'Miller Cylindrical';...
                'Mollweide';'Robinson';'Sinusoidal'})
        set_enable(handles, 0);     set_enable(handles, 50);      handles.nParameters = 0;
        handles.projection{1} = '-Jm0/0/1';
        set(handles.popup_Datum,'Value',221);       info.DX = handles.all_datums{221,4};
        info.DY = handles.all_datums{221,5};        info.DZ = handles.all_datums{221,6};
        info.ellipsoide = handles.all_datums{221,2};        datum_info(handles, info)
    case 'All GMT Projections'
        set(handles.popup_System,'Value',1,'String',{'Use the right side boxes'})
        set_enable(handles, 0);     set_enable(handles, 51);      handles.nParameters = [];
    case 'Portuguese Mess'
        caso.type = 'popup';  caso.val = 1;
        handles = set_PtMess(handles,caso);
    case 'Linear Projection'        % Used only in write_script GUI
        set(handles.popup_System,'String','Linear Projection')
        set_enable(handles, 0);        set_enable(handles, 50);
        handles.nParameters = 0;       handles.projection{1} = '-Jx';
        handles.proj_name = 'Linear Projection';
end
handles.group_val = val;
guidata(hObject,handles)

%-------------------------------------------------------------------------------------
function popup_System_CB(hObject, handles)
contents = get(hObject,'String');   n_value = get(hObject,'Value');
if (iscell(contents))
    system = contents{n_value};
else    % For example the 'Longitude / Latitude' group has only one field in Systems
    system = contents;
end
switch handles.proj_groups{handles.group_val}
    case 'Longitude / Latitude'
        set_enable(handles,0);      handles.nParameters = 0;
        handles.projection{1} = []; handles.proj_name = 'Longitude / Latitude';
        handles.sys_name = 'Geographic';
    case 'Universal Transversal Mercator'
        if (n_value <= 60)
            handles.projection{1} = ['-Ju' system(6:7) '/1'];
        else            % Southern hemisphere
            handles.projection{1} = ['-Ju-' system(6:7) '/1'];
        end
        set_enable(handles,0);      handles.nParameters = 0;
        handles.proj_name = 'UTM';  handles.map_scale_factor = 0.9996;
        handles.sys_name = system;
    case 'World Projections'
        caso.type = 'system';  caso.val = system;
        handles = set_WorldProj(handles,caso);
        handles.sys_name = system;
    case 'All GMT Projections'
        set(hObject,'String','Use the right side boxes')
        set_enable(handles,51)
    case 'Portuguese Mess'
        caso.type = 'system';  caso.val = system;
        handles = set_PtMess(handles,caso);
        handles.sys_name = system;
end
guidata(hObject,handles)

%-------------------------------------------------------------------------------------
function popup_Datum_CB(hObject, handles)
	% Only updates the info text about datum parameters
	val = get(hObject,'Value');
	info.DX = handles.all_datums{val,4};       info.DY = handles.all_datums{val,5};
	info.DZ = handles.all_datums{val,6};       info.ellipsoide = handles.all_datums{val,2};
	datum_info(handles, info)

%-------------------------------------------------------------------------------------
function popup_CilindricalProjections_CB(hObject, handles)
% This function only sets the apropiate -J and opens the parameters boxes. The OK button reads them
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case 'CYLINDRICAL'                  % Reset parameter names to void 
        set_enable(handles,0);      handles.nParameters = [];
    case 'Cassini Cylindrical'
        set_enable(handles,2);      handles.nParameters = 2;
        set_params_str(handles, {handles.lon0_txt; handles.lat0_txt}, 2, [2 3 4])
        handles.projection{1} = '-Jc';
    case 'Miller Cylindrical'
        set_enable(handles,1);      handles.nParameters = 1;
        set_params_str(handles, {handles.lon0_txt}, 1, [2 3 4])
        handles.projection{1} = '-Jj';
    case 'Mercator - Greenwich and Equator as origin'
        set_enable(handles,0);      handles.nParameters = 0;
        handles.projection{1} = '-Jm0/0/1';
    case 'Mercator - Give meridian and standard parallel'
        set_enable(handles,2);      handles.nParameters = 2;
        set_params_str(handles, {handles.lon0_txt; handles.lat0_txt}, 2, [2 3 4])
        handles.projection{1} = '-Jm';
    case 'Oblique Mercator - point & azimuth'
        set_enable(handles,3);      handles.nParameters = 3;
        set_params_str(handles, {handles.lon0_txt; handles.lat0_txt;'Azimuth'}, 3, [2 3 4])
        handles.projection{1} = '-Joa';
    case 'Oblique Mercator - two points'
        set_enable(handles,4);      handles.nParameters = 4;
        set_params_str(handles, {'Lon_of_First_Point';'Lat_of_First_Point';'Lon_of_Second_Point';'Lat_of_Second_Point'}, 4, [2 3 4])
        handles.projection{1} = '-Job';
    case 'Oblique Mercator - point & pole'
        set_enable(handles,4);      handles.nParameters = 4;
        set_params_str(handles, {handles.lon0_txt;handles.lat0_txt;'Lon_of_the_Pole';'Lat_of_the_Pole'}, 4, [2 3 4])
        handles.projection{1} = '-Joc';
    case 'Cylindrical Equidistant'
        set_enable(handles,1);      handles.nParameters = 1;
        set_params_str(handles, {handles.lon0_txt}, 1, [2 3 4])
        handles.projection{1} = '-Jq';
    case 'TM - Transverse Mercator, with Equator as y=0'
        set_enable(handles,1);      handles.nParameters = 1;
        set_params_str(handles, {handles.lon0_txt}, 1, [2 3 4])
        handles.projection{1} = '-Jt';
        handles.map_scale_factor = 1;
    case 'TM - Transverse Mercator, set origin'
        set_enable(handles,2);      handles.nParameters = 2;
        set_params_str(handles, {handles.lon0_txt; handles.latOrig_txt}, 2, [2 3 4])
        handles.projection{1} = '-Jt';
        handles.map_scale_factor = 1;
    case 'UTM - Universal Transverse Mercator'
        set(handles.popup_System,'String',handles.UTMzones)
        handles.group_val = 2;      % Need this - special case
        set_enable(handles,0);      handles.nParameters = 0;
        handles.projection{1} = '-Ju';
        handles.map_scale_factor = 0.9996;
    case 'General Cylindrical'
        set_enable(handles,0);      handles.nParameters = 2;
        set_params_str(handles, {handles.lon0_txt; handles.lat0_txt}, 2, [2 3 4])
        handles.projection{1} = '-Jy';
end
handles.proj_name = str{val};
guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function popup_AzimuthalProjections_CB(hObject, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case 'AZIMUTHAL'                  % Reset parameter names to void 
        set_enable(handles,0);      handles.nParameters = [];
    case 'Lambert Azimuthal Equal-Area'
        set_enable(handles,2);      handles.nParameters = 2;
        set_params_str(handles, {handles.lon0_txt; handles.lat0_txt}, 2, [1 3 4])
        handles.projection{1} = '-JA';
    case 'Azimuthal Equidistant'
        set_enable(handles,2);      handles.nParameters = 2;
        set_params_str(handles, {handles.lon0_txt; handles.lat0_txt}, 2, [1 3 4])
        handles.projection{1} = '-JE';
    case 'Gnomonic'
        set_enable(handles,3);      handles.nParameters = 3;
        set_params_str(handles, {handles.lon0_txt; handles.lat0_txt;...
                'number of degrees from the center to the edge'}, 3, [1 3 4])
        handles.projection{1} = '-JF';
    case 'Orthographic'
        set_enable(handles,2);      handles.nParameters = 2;
        set_params_str(handles, {handles.lon0_txt; handles.lat0_txt}, 2, [1 3 4])
        handles.projection{1} = '-JG';
    case 'Stereographic Equal-Angle'
        set_enable(handles,2);      handles.nParameters = 2;
        set_params_str(handles, {handles.lon0_txt; handles.lat0_txt}, 2, [1 3 4])
        handles.projection{1} = '-JS';
end
handles.proj_name = str{val};
guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function popup_ConicalProjections_CB(hObject, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
	case 'CONIC'                  % Reset parameter names to void 
	    set_enable(handles,0);      handles.nParameters = [];
	case 'Albers Conic Equal-Area'
    	set_enable(handles,4);      handles.nParameters = 4;
    	set_params_str(handles, {handles.lon0_txt; handles.lat0_txt;'Lat_of_Southern_Parallel';'Lat-of_Northern_Parallel'}, 4, [1 2 4])
	    handles.projection{1} = '-Jb';
	case 'Equidistant Conic'
    	set_enable(handles,4);      handles.nParameters = 4;
	    set_params_str(handles, {handles.lon0_txt; handles.lat0_txt;'Lat_of_Southern_Parallel';'Lat_of_Northern_Parallel'}, 4, [1 2 4])
	    handles.projection{1} = '-Jd';
	case 'Lambert Conic Conformal'
	    set_enable(handles,4);      handles.nParameters = 4;
	    set_params_str(handles, {handles.lon0_txt; handles.lat0_txt;'Lat_of_Southern_Parallel';'Lat_of_Northern_Parallel'}, 4, [1 2 4])
	    handles.projection{1} = '-Jl';
end
handles.proj_name = str{val};
guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function popup_MiscelaneousProjections_CB(hObject, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
switch str{val};
    case 'MISCELLANEOUS'                  % Reset parameter names to void
        set_enable(handles,0);      handles.nParameters = [];
    case 'Hammer'
        set_enable(handles,1);      handles.nParameters = 1;
        set_params_str(handles, {handles.lon0_txt}, 1, [1 2 3])
        handles.projection{1} = '-Jh';
    case 'Sinusoidal'
        set_enable(handles,1);      handles.nParameters = 1;
        set_params_str(handles, {handles.lon0_txt}, 1, [1 2 3])
        handles.projection{1} = '-Ji';
    case 'Eckert IV'
        set_enable(handles,1);      handles.nParameters = 1;
        set_params_str(handles, {handles.lon0_txt}, 1, [1 2 3])
        handles.projection{1} = '-Jkf';
    case 'Eckert VI'
        set_enable(handles,1);      handles.nParameters = 1;
        set_params_str(handles, {handles.lon0_txt}, 1, [1 2 3])
        handles.projection{1} = '-Jks';
    case 'Robinson'
        set_enable(handles,1);      handles.nParameters = 1;
        set_params_str(handles, {handles.lon0_txt}, 1, [1 2 3])
        handles.projection{1} = '-Jn';
    case 'Winkel Tripel'
        set_enable(handles,1);      handles.nParameters = 1;
        set_params_str(handles, {handles.lon0_txt}, 1, [1 2 3])
        handles.projection{1} = '-Jr';
    case 'Van der Grinten'
        set_enable(handles,1);      handles.nParameters = 1;
        set_params_str(handles, {handles.lon0_txt}, 1, [1 2 3])
        handles.projection{1} = '-Jv';
    case 'Mollweide'
        set_enable(handles,1);      handles.nParameters = 1;
        set_params_str(handles, {handles.lon0_txt}, 1, [1 2 3])
        handles.projection{1} = '-Jw';
end
handles.proj_name = str{val};
guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function edit_ProjParameterValue_1_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || isempty(xx)),  set(hObject, 'String', '');  return;    end    % Just a stupid user error
	xx = get(hObject,'String');     val = test_dms(xx);     x = 0;
	if (str2double(val{1}) > 0)
		for i=1:length(val),  x = x + str2double(val{i}) / (60^(i-1));    end
	else
		for i=1:length(val),  x = x - abs(str2double(val{i})) / (60^(i-1));   end
	end
	x = num2str(x);
	if (handles.nParameters == 1),  handles.projection{2} = [x '/1'];
	else                            handles.projection{2} = [x '/'];    end
	guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function edit_ProjParameterValue_2_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || isempty(xx)),  set(hObject, 'String', '');  return;    end    % Just a stupid user error
	xx = get(hObject,'String');     val = test_dms(xx);     x = 0;
	if (str2double(val{1}) > 0)
        for i=1:length(val),  x = x + str2double(val{i}) / (60^(i-1));    end
	else
        for i=1:length(val),  x = x - abs(str2double(val{i})) / (60^(i-1));   end
	end
	x = num2str(x);
	if (handles.nParameters == 2),  handles.projection{3} = [x '/1'];
	else                            handles.projection{3} = [x '/'];    end
	guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function edit_ProjParameterValue_3_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || isempty(xx)),  set(hObject, 'String', '');  return;    end    % Just a stupid user error
	xx = get(hObject,'String');     val = test_dms(xx);     x = 0;
	if (str2double(val{1}) > 0)
        for i=1:length(val),  x = x + str2double(val{i}) / (60^(i-1));    end
	else
        for i=1:length(val),  x = x - abs(str2double(val{i})) / (60^(i-1));   end
	end
	x = num2str(x);
	if (handles.nParameters == 3),  handles.projection{4} = [x '/1'];
	else                            handles.projection{4} = [x '/'];    end
	guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function edit_ProjParameterValue_4_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || isempty(xx)),  set(hObject, 'String', '');  return;    end    % Just a stupid user error
	xx = get(hObject,'String');     val = test_dms(xx);     x = 0;
	if (str2double(val{1}) > 0)
		for i=1:length(val),  x = x + str2double(val{i}) / (60^(i-1));    end
	else
		for i=1:length(val),  x = x - abs(str2double(val{i})) / (60^(i-1));   end
	end
	x = num2str(x);
	handles.projection{5} = [x '/1'];
	guidata(hObject, handles);

%-------------------------------------------------------------------------------------
function pushbutton_OK_CB(hObject, handles)
if (isempty(handles.nParameters))
	errordlg('You did not select anything. So what are you expecting?','Chico Clever');    return
end
erro = 0;
for (i=1:handles.nParameters)
	if (isempty(handles.projection{i+1}))
		errordlg(['Missing projection parameter number ' num2str(i)],'Error');    erro = 1;
	end
end
if (erro),     return;     end
Out.group_val = get(handles.popup_Group,'Value');
Out.system_val = get(handles.popup_System,'Value');
Out.datum_val = get(handles.popup_Datum,'Value');
Out.cilindrical_val = get(handles.popup_CilindricalProjections,'Value');
Out.azimuthal_val = get(handles.popup_AzimuthalProjections,'Value');
Out.conic_val = get(handles.popup_ConicalProjections,'Value');
Out.miscelaneous_val = get(handles.popup_MiscelaneousProjections,'Value');
Out.ProjName = handles.proj_name;
Out.SysName = handles.sys_name;
Out.map_scale_factor = handles.map_scale_factor;
Out.system_FE_FN = handles.system_FE_FN;

if (handles.nParameters == 0 && ~strcmp(handles.projection{1},'-Ju'))
    Out.projection = handles.projection{1};
    Out.ProjParameterValue = [];
elseif (handles.nParameters == 0 && strcmp(handles.projection{1},'-Ju'))
    n_value = get(handles.popup_System,'Value');
    UTMzone_txt = get(handles.popup_System,'String');
    system = UTMzone_txt{n_value};
    if (n_value <= 60)  % Northern hemisphere
        Out.projection = ['-Ju' system(6:7) '/1'];
    else                % Southern hemisphere
        Out.projection = ['-Ju-' system(6:7) '/1'];
    end    
    Out.ProjParameterValue = [];
elseif (handles.nParameters == 0 && strcmp(handles.projection{1},'-Jx'))
    Out.projection = handles.projection{1};
    Out.ProjParameterValue = [];
else
    Out.projection = cat(2,handles.projection{:});
	for (i=1:handles.nParameters)
        switch i
            case 1,     Out.ProjParameterValue.p1 = handles.projection{i+1};
            case 2,     Out.ProjParameterValue.p2 = handles.projection{i+1};
            case 3,     Out.ProjParameterValue.p3 = handles.projection{i+1};
            case 4,     Out.ProjParameterValue.p4 = handles.projection{i+1};
        end
	end
end

handles.output = Out;           guidata(hObject,handles)
uiresume(handles.figure1);

%-------------------------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	if (exist('OCTAVE_VERSION','builtin'))		% To know if we are running under Octave
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	else
		do_uiresume = strcmp(get(handles.figure1, 'waitstatus'), 'waiting');
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		handles.output = [];		% User gave up, return nothing
		guidata(handles.figure1, handles);	uiresume(handles.figure1);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

%-------------------------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata)
handles = guidata(hObject);
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = [];        % User said no by hitting escape
    guidata(hObject, handles);    uiresume(handles.figure1);
end

%-------------------------------------------------------------------------------------
function set_enable(handles, opt)
cor = handles.bg_color;     str = 'NOTHING HERE';
switch opt
	case 0
		set(handles.edit_ProjParameterValue_1,'String',str,'Enable','off','Backgroundcolor',cor)
		set(handles.edit_ProjParameterValue_2,'String',str,'Enable','off','Backgroundcolor',cor)
		set(handles.edit_ProjParameterValue_3,'String',str,'Enable','off','Backgroundcolor',cor)
		set(handles.edit_ProjParameterValue_4,'String',str,'Enable','off','Backgroundcolor',cor)
	case 1
		set(handles.edit_ProjParameterValue_1,'String','','Enable','on','Backgroundcolor',[1 1 1])
		set(handles.edit_ProjParameterValue_2,'String',str,'Enable','off','Backgroundcolor',cor)
		set(handles.edit_ProjParameterValue_3,'String',str,'Enable','off','Backgroundcolor',cor)
		set(handles.edit_ProjParameterValue_4,'String',str,'Enable','off','Backgroundcolor',cor)
	case 2
		set(handles.edit_ProjParameterValue_1,'String','','Enable','on','Backgroundcolor',[1 1 1])
		set(handles.edit_ProjParameterValue_2,'String','','Enable','on','Backgroundcolor',[1 1 1])
		set(handles.edit_ProjParameterValue_3,'String',str,'Enable','off','Backgroundcolor',cor)
		set(handles.edit_ProjParameterValue_4,'String',str,'Enable','off','Backgroundcolor',cor)
	case 3
		set(handles.edit_ProjParameterValue_1,'String','','Enable','on','Backgroundcolor',[1 1 1])
		set(handles.edit_ProjParameterValue_2,'String','','Enable','on','Backgroundcolor',[1 1 1])
		set(handles.edit_ProjParameterValue_3,'String','','Enable','on','Backgroundcolor',[1 1 1])
		set(handles.edit_ProjParameterValue_4,'String',str,'Enable','off','Backgroundcolor',cor)
	case 4
		set(handles.edit_ProjParameterValue_1,'String','','Enable','on','Backgroundcolor',[1 1 1])
		set(handles.edit_ProjParameterValue_2,'String','','Enable','on','Backgroundcolor',[1 1 1])
		set(handles.edit_ProjParameterValue_3,'String','','Enable','on','Backgroundcolor',[1 1 1])
		set(handles.edit_ProjParameterValue_4,'String','','Enable','on','Backgroundcolor',[1 1 1])
	case 50     % Desable the "All GMT projections popups"
		set(handles.popup_CilindricalProjections,'Value',1,'Enable','off')
		set(handles.popup_AzimuthalProjections,'Value',1,'Enable','off')
		set(handles.popup_ConicalProjections,'Value',1,'Enable','off')
		set(handles.popup_MiscelaneousProjections,'Value',1,'Enable','off')
	case 51     % Enable the "All GMT projections popups"
		set(handles.popup_CilindricalProjections,'Value',1,'Enable','on')
		set(handles.popup_AzimuthalProjections,'Value',1,'Enable','on')
		set(handles.popup_ConicalProjections,'Value',1,'Enable','on')
		set(handles.popup_MiscelaneousProjections,'Value',1,'Enable','on')
end

% Clean eventual remaining projection parameters from a previous use
if (opt <= 4)
    for (i=opt+2:5),    handles.projection{i} = [];     end
    guidata(handles.fig_hand, handles);
end

%-------------------------------------------------------------------------------------
function set_params_str(handles, str, num, top)
if nargin == 3, top = [0 0 0];    end
switch num
    case 1
        set(handles.edit_ProjParameterValue_1, 'String', str{1});
    case 2
        set(handles.edit_ProjParameterValue_1, 'String', str{1});
        set(handles.edit_ProjParameterValue_2, 'String', str{2});
    case 3
        set(handles.edit_ProjParameterValue_1, 'String', str{1});
        set(handles.edit_ProjParameterValue_2, 'String', str{2});
        set(handles.edit_ProjParameterValue_3, 'String', str{3});
    case 4
        set(handles.edit_ProjParameterValue_1, 'String', str{1});
        set(handles.edit_ProjParameterValue_2, 'String', str{2});
        set(handles.edit_ProjParameterValue_3, 'String', str{3});
        set(handles.edit_ProjParameterValue_4, 'String', str{4});
end

% Set the other groups to their name group
for (i=1:3)
	if (top(i) == 1)
        set(handles.popup_CilindricalProjections,'Value',1);
	elseif (top(i) == 2)
        set(handles.popup_AzimuthalProjections,'Value',1);
	elseif (top(i) == 3)
        set(handles.popup_ConicalProjections,'Value',1);
	elseif (top(i) == 4)
        set(handles.popup_MiscelaneousProjections,'Value',1);
	end
end

%-------------------------------------------------------------------------------------
function handles = set_PtMess(handles, caso)
% Handle the Pt case
val = caso.val;     set_enable(handles, 50);
if (strcmp(caso.type,'popup') || strcmp(caso.type,'input'))
    set(handles.popup_System,'Value',val,'String',{'UTM (ED50)';'Gauss (D73)';'Gauss-Militar (DLx)';...
            'Gauss (WGS84)';'Base SW';'S. Bras';'Observatorio';'Porto Santo'})
end
switch val
    case {1, 'UTM (ED50)'}
        id = 69;        handles.nParameters = 0;
        handles.projection{1} = '-Ju29/1';
        handles.proj_name = 'UTM - Universal Transverse Mercator';
        handles.sys_name = 'UTM';
        handles.map_scale_factor = 0.9996;
    case {2, 'Gauss (D73)'}
        id = 49;        handles.nParameters = 2;
        handles.projection{1} = '-Jt';
        handles.projection{2} = '-8:07:54.862/';
        handles.projection{3} = '39:40:00/1';
        set_params_str(handles, {'-8:07:54.862'; '39:40:00'}, 2, [1 2 3 4])
        handles.proj_name = 'TM - Transverse Mercator, set origin';
        handles.sys_name = 'Gauss (D73)';
        handles.map_scale_factor = 1;
        handles.system_FE_FN = [180.598 -86.990];
    case {3, 'Gauss-Militar (DLx)'}
        id = 50;        handles.nParameters = 2;
        handles.projection{1} = '-Jt';
        handles.projection{2} = '-8:07:54.862/';
        handles.projection{3} = '39:40:00/1';
        set_params_str(handles, {'-8:07:54.862'; '39:40:00'}, 2, [1 2 3 4])
        handles.proj_name = 'TM - Transverse Mercator, set origin';
        handles.sys_name = 'Gauss-Militar (DLx)';
        handles.map_scale_factor = 1;
        handles.system_FE_FN = [200000 300000];
    case {4, 'Gauss (WGS84)'}
        id = 221;        handles.nParameters = 2;
        handles.projection{1} = '-Jt';
        handles.projection{2} = '-8:07:59.191/';
        handles.projection{3} = '39:40:5.73/1';
        set_params_str(handles, {'-8:07:59.191'; '39:40:5.73'}, 2, [1 2 3 4])
        handles.proj_name = 'TM - Transverse Mercator, set origin';
        handles.sys_name = 'Gauss (WGS84)';
        handles.map_scale_factor = 1;
    case {5, 'Base SW'}
        id = 75;        handles.nParameters = 0;
        handles.projection{1} = '-Ju26/1';
        handles.proj_name = 'UTM - Universal Transverse Mercator';
        handles.sys_name = 'UTM';
        handles.map_scale_factor = 0.9996;
    case {6, 'S. Bras'}
        id = 189;        handles.nParameters = 0;
        handles.projection{1} = '-Ju26/1';
        handles.proj_name = 'UTM - Universal Transverse Mercator';
        handles.sys_name = 'UTM';
        handles.map_scale_factor = 0.9996;
    case {7, 'Observatorio'}
        id = 146;        handles.nParameters = 0;
        handles.projection{1} = '-Ju25/1';
        handles.proj_name = 'UTM - Universal Transverse Mercator';
        handles.sys_name = 'UTM';
        handles.map_scale_factor = 0.9996;
    case {8, 'Porto Santo'}
        id = 163;        handles.nParameters = 0;
        handles.projection{1} = '-Ju28/1';
        handles.proj_name = 'UTM - Universal Transverse Mercator';
        handles.sys_name = 'UTM';
        handles.map_scale_factor = 0.9996;
end
set(handles.popup_Datum,'Value',id);
set_enable(handles, 0);
info.DX = handles.all_datums{id,4};       info.DY = handles.all_datums{id,5};
info.DZ = handles.all_datums{id,6};       info.ellipsoide = handles.all_datums{id,2};
datum_info(handles, info)

%-------------------------------------------------------------------------------------
function handles = set_WorldProj(handles, caso)
% Handle the World Projections case
if (strcmp(caso.type,'input'))
    set(handles.popup_System,'String',{'Mercator';'Miller Cylindrical';'Mollweide';'Robinson';'Sinusoidal'})
end
set_enable(handles,50);
val = caso.val;             id = 221;   % World Projs use wgs by default
switch val
    case {1, 'Mercator'}
        set_enable(handles,0);     handles.nParameters = 0;
        handles.projection{1} = '-Jm0/0/1'; handles.proj_name = 'Mecator (Origin in equator)';
    case {2, 'Miller Cylindrical'}
        set_enable(handles,1);     handles.nParameters = 1;
        handles.projection{1} = '-Jj';  handles.proj_name = 'Miller Cylindrical';
        set_params_str(handles, {handles.lon0_txt}, 1, [2 3 4])
    case {3, 'Mollweide'}
        set_enable(handles,1);     handles.nParameters = 1;
        handles.projection{1} = '-Jw';  handles.proj_name = 'Mollweide';
        set_params_str(handles, {handles.lon0_txt}, 1, [2 3 4])
    case {4, 'Robinson'}
        set_enable(handles,1);     handles.nParameters = 1;
        handles.projection{1} = '-Jn';  handles.proj_name = 'Robinson';
        set_params_str(handles, {handles.lon0_txt}, 1, [2 3 4])
    case {5, 'Sinusoidal'}
        set_enable(handles,1);     handles.nParameters = 1;
        handles.projection{1} = '-Ji';  handles.proj_name = 'Sinusoidal';
        set_params_str(handles, {handles.lon0_txt}, 1, [2 3 4])
end
for (i=2:handles.nParameters+1)     % Make sure to reset the projection
    handles.projection{i} = '';
end
set(handles.popup_Datum,'Value',id);
info.DX = handles.all_datums{id,4};       info.DY = handles.all_datums{id,5};
info.DZ = handles.all_datums{id,6};       info.ellipsoide = handles.all_datums{id,2};
datum_info(handles, info)

%------------ Write info text about datum parameters -------------------
function datum_info(handles, info)
	h_txt = findobj(handles.fig_hand,'Tag','datum_desc');
	pos = get(h_txt,'Position');
	sft = ['DX = ' num2str(info.DX) '    DY = ' num2str(info.DY) '    DZ = ' num2str(info.DZ)];
	string = {['Ellipsoise -> ' info.ellipsoide]; ['Shift ->   ' sft]};
	[outstring,newpos] = textwrap(h_txt,string);
	pos(4) = newpos(4);
	set(h_txt,'String',outstring,'Position',[pos(1),pos(2),pos(3),pos(4)])


% --- Creates and returns a handle to the GUI figure. 
function coordinate_system_LayoutFcn(h1)

set(h1, 'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','coordinate_system',...
'NumberTitle','off',...
'Position',[520 580 530 220],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[10 78 251 139],'Style','frame');
uicontrol('Parent',h1,'Position',[270 78 251 139],'Style','frame');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@coordinate_system_uiCB,...
'Position',[51 186 201 22],...
'String',{'Longitude / Latitude'; 'Universal Transversal Mercator'; 'World Projections'; 'All GMT Projections'; 'Portuguese Mess' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_Group');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[14 190 35 15],...
'String','Group',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@coordinate_system_uiCB,...
'Position',[51 155 201 22],...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_System');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[12 159 38 15],...
'String','System',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@coordinate_system_uiCB,...
'Position',[51 47 321 22],...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_Datum');

uicontrol('Parent',h1,...
'HorizontalAlignment','right',...
'Position',[13 51 35 15],...
'String','Datum',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@coordinate_system_uiCB,...
'Position',[280 186 231 22],...
'String',{'CYLINDRICAL'; 'Cassini Cylindrical'; 'Miller Cylindrical'; 'Mercator - Greenwich and Equator as origin'; ...
	'Mercator - Give meridian and standard parallel'; 'Oblique Mercator - point & azimuth'; 'Oblique Mercator - two points';...
	'Oblique Mercator - point & pole'; 'Cylindrical Equidistant'; 'TM - Transverse Mercator, with Equator as y=0'; ...
	'TM - Transverse Mercator, set origin'; 'UTM - Universal Transverse Mercator'; 'General Cylindrical'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_CilindricalProjections');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@coordinate_system_uiCB,...
'Position',[280 155 231 22],...
'String',{'AZIMUTHAL'; 'Lambert Azimuthal Equal-Area'; 'Azimuthal Equidistant'; 'Gnomonic'; 'Orthographic'; 'Stereographic Equal-Angle' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_AzimuthalProjections');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@coordinate_system_uiCB,...
'Position',[280 126 231 22],...
'String',{'CONIC'; 'Albers Conic Equal-Area'; 'Equidistant Conic'; 'Lambert Conic Conformal'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_ConicalProjections');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@coordinate_system_uiCB,...
'Position',[280 95 231 22],...
'String',{'MISCELLANEOUS'; 'Hammer'; 'Sinusoidal'; 'Eckert IV'; 'Eckert VI'; 'Robinson'; 'Winkel Tripel'; 'Van der Grinten'; 'Mollweide' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_MiscelaneousProjections');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@coordinate_system_uiCB,...
'HorizontalAlignment','center',...
'Position',[20 115 111 21],...
'Style','edit',...
'Tooltip','Enter value of projection parameter 1',...
'Tag','edit_ProjParameterValue_1');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[75 137 135 14],...
'String','Map Projection Parameters',...
'Style','text');

uicontrol('Parent',h1,...
'Call',@coordinate_system_uiCB,...
'Position',[450 6 66 21],...
'String','OK',...
'Tag','pushbutton_OK');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@coordinate_system_uiCB,...
'HorizontalAlignment','center',...
'Position',[140 115 111 21],...
'Style','edit',...
'Tooltip','Enter value of projection parameter 2',...
'Tag','edit_ProjParameterValue_2');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@coordinate_system_uiCB,...
'HorizontalAlignment','center',...
'Position',[20 90 111 21],...
'Style','edit',...
'Tooltip','Enter value of projection parameter 3',...
'Tag','edit_ProjParameterValue_3');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@coordinate_system_uiCB,...
'HorizontalAlignment','center',...
'Position',[140 90 111 21],...
'Style','edit',...
'Tooltip','Enter value of projection parameter 4',...
'Tag','edit_ProjParameterValue_4');

uicontrol('Parent',h1,...
'Enable','inactive',...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[51 8 351 31],...
'String','map',...
'Style','text',...
'Tag','datum_desc');

function coordinate_system_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

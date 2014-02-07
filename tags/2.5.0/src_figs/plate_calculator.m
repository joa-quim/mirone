function varargout = plate_calculator(varargin)
% Calculate plate velocities

%	Copyright (c) 2004-2013 by J. Luis
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

% $Id: $

	hObject = figure('Tag','figure1','Visible','off');
	plate_calculator_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	mir_dirs = getappdata(0,'MIRONE_DIRS');
	if (~isempty(mir_dirs))
		handles.home_dir = mir_dirs.home_dir;		% Start in values
		handles.work_dir = mir_dirs.work_dir;
		handles.last_dir = mir_dirs.last_dir;
	else
		handles.home_dir = cd;		handles.work_dir = cd;		handles.last_dir = cd;
	end
	handles.path_data = [handles.home_dir filesep 'data' filesep];

	handles.first_NNR = true;
	handles.first_PB = true;
	handles.first_AKIM2000 = true;
	handles.first_DEOS2K = true;
	handles.first_REVEL = true;
	handles.first_MORVEL = true;
	handles.first_GEODVEL = true;
	handles.absolute_motion = 0;		% when == 1, it signals an absolute motion model
	handles.abs2rel = 0;				% when == 1, flags that an absolute model was turned relative

	set(handles.checkbox_Abs2Rel,'Visible','off')

	% Read the Nuvel-1A poles file as they are the default
	fid = fopen([handles.path_data 'Nuvel1A_poles.dat'],'r');
	[abbrev name lat lon omega] = strread(fread(fid,'*char'),'%s %s %f %f %f');
	fclose(fid);

	% Save the poles parameters in the handles structure
	handles.Nuvel1A_abbrev = abbrev;
	handles.Nuvel1A_name = name;
	handles.Nuvel1A_lat = lat;
	handles.Nuvel1A_lon = lon;
	handles.Nuvel1A_omega = omega;

	% Fill the popupmenus with the Plate's names
	set([handles.popup_FixedPlate handles.popup_MovingPlate],'String',name)

	set_Nuvel1Aplate_model(hObject,handles)
	setappdata(hObject,'current_model','Nuvel1A')
	handles.Nuvel1A_comb = do_plate_comb('Nuvel1A');

	% Need to change the ButtonDownFcn call arguments (I didn't set it directly for a question of generality)
	h_patch = findobj(hObject,'Type','patch');
	set(h_patch,'ButtonDownFcn',{@bdn_plate,handles,'Nuvel1A'})

	% Initialize (but out of the map) two symbols for ploting the pole position
	line(-500,-500,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',7,'Tag','pole_out');
	line(-500,-500,'Marker','+','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,'Tag','pole_in1');
	line(-500,-500,'Marker','o','MarkerEdgeColor','k','MarkerSize',8,'Tag','pole_in2');

	%----------- Give a Pro look (3D) to the frame box  ---------
	new_frame3D(hObject, [handles.txt_Abs handles.txt_Rel])
	%------------- END Pro look (3D) ----------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

%--------------------------------------------------------------------------------------------------
function popup_PickPlate_CB(hObject, handles, qual)
	D2R = pi/180;
	if (strcmp(qual, 'fixed'))
		ind_fix = get(hObject,'Value');
		ind_mov = get(handles.popup_MovingPlate,'Value');
	else
		ind_mov = get(hObject,'Value');
		ind_fix = get(handles.popup_FixedPlate,'Value');
	end
	model = getappdata(handles.figure1, 'current_model');

	switch model
		case 'Nuvel1A'
			lat2 = handles.Nuvel1A_lat(ind_mov);		lon2 = handles.Nuvel1A_lon(ind_mov);
			omega2 = handles.Nuvel1A_omega(ind_mov);	handles.abb_mov = handles.Nuvel1A_abbrev{ind_mov};
		case 'NNR'
			lat2 = handles.Nuvel1A_NNR_lat(ind_mov);	lon2 = handles.Nuvel1A_NNR_lon(ind_mov);
			omega2 = handles.Nuvel1A_NNR_omega(ind_mov);handles.abb_mov = handles.Nuvel1A_NNR_abbrev{ind_mov};
		case 'PB'
			lat2 = handles.PB_lat(ind_mov);				lon2 = handles.PB_lon(ind_mov);
			omega2 = handles.PB_omega(ind_mov);			handles.abb_mov = handles.PB_abbrev{ind_mov};
		case 'AKIM2000'
			lat2 = handles.AKIM2000_lat(ind_mov);		lon2 = handles.AKIM2000_lon(ind_mov);
			omega2 = handles.AKIM2000_omega(ind_mov);	handles.abb_mov = handles.AKIM2000_abbrev{ind_mov};
		case 'REVEL'
			lat2 = handles.REVEL_lat(ind_mov);			lon2 = handles.REVEL_lon(ind_mov);
			omega2 = handles.REVEL_omega(ind_mov);		handles.abb_mov = handles.REVEL_abbrev{ind_mov};
		case 'DEOS2K'
			lat2 = handles.DEOS2K_lat(ind_mov);			lon2 = handles.DEOS2K_lon(ind_mov);
			omega2 = handles.DEOS2K_omega(ind_mov);		handles.abb_mov = handles.DEOS2K_abbrev{ind_mov};
		case 'MORVEL'
			lat2 = handles.MORVEL_lat(ind_mov);			lon2 = handles.MORVEL_lon(ind_mov);
			omega2 = handles.MORVEL_omega(ind_mov);		handles.abb_mov = handles.MORVEL_abbrev{ind_mov};
		case 'GEODVEL'
			lat2 = handles.GEODVEL_lat(ind_mov);		lon2 = handles.GEODVEL_lon(ind_mov);
			omega2 = handles.GEODVEL_omega(ind_mov);	handles.abb_mov = handles.GEODVEL_abbrev{ind_mov};
	end

	if ~(handles.absolute_motion)       % That is, if relative motion
		switch model
			case 'Nuvel1A'
				lat1 = handles.Nuvel1A_lat(ind_fix);		lon1 = handles.Nuvel1A_lon(ind_fix);
				omega1 = handles.Nuvel1A_omega(ind_fix);	handles.abb_fix = handles.Nuvel1A_abbrev{ind_fix};
			case 'NNR'
				lat1 = handles.Nuvel1A_NNR_lat(ind_fix);	lon1 = handles.Nuvel1A_NNR_lon(ind_fix);
				omega1 = handles.Nuvel1A_NNR_omega(ind_fix);handles.abb_fix = handles.Nuvel1A_NNR_abbrev{ind_fix};
			case 'PB'
				lat1 = handles.PB_lat(ind_fix);				lon1 = handles.PB_lon(ind_fix);
				omega1 = handles.PB_omega(ind_fix);			handles.abb_fix = handles.PB_abbrev{ind_fix};
			case 'AKIM2000'
				lat1 = handles.AKIM2000_lat(ind_fix);		lon1 = handles.AKIM2000_lon(ind_fix);
				omega1 = handles.AKIM2000_omega(ind_fix);	handles.abb_fix = handles.AKIM2000_abbrev{ind_fix};       
			case 'REVEL'
				lat1 = handles.REVEL_lat(ind_fix);			lon1 = handles.REVEL_lon(ind_fix);
				omega1 = handles.REVEL_omega(ind_fix);		handles.abb_fix = handles.REVEL_abbrev{ind_fix};
			case 'DEOS2K'
				lat1 = handles.DEOS2K_lat(ind_fix);			lon1 = handles.DEOS2K_lon(ind_fix);
				omega1 = handles.DEOS2K_omega(ind_fix);		handles.abb_fix = handles.DEOS2K_abbrev{ind_fix};       
			case 'MORVEL'
				lat1 = handles.MORVEL_lat(ind_fix);			lon1 = handles.MORVEL_lon(ind_fix);
				omega1 = handles.MORVEL_omega(ind_fix);		handles.abb_fix = handles.MORVEL_abbrev{ind_fix};
			case 'GEODVEL'
				lat2 = handles.GEODVEL_lat(ind_fix);		lon2 = handles.GEODVEL_lon(ind_fix);
				omega2 = handles.GEODVEL_omega(ind_fix);	handles.abb_mov = handles.GEODVEL_abbrev{ind_fix};
		end
		[lon,lat,omega] = calculate_pole(lon1,lat1,omega1,lon2,lat2,omega2);
	else                                % Absolute motion
		lon = lon2;     lat = lat2;     omega = omega2;
		handles.abb_fix = 'absolute';
	end

	if (omega == 0)     % This works as a test for when the same plate is selected as Fixed and Moving
		set([handles.edit_PoleLon handles.edit_PoleLat handles.edit_PoleRate],'String','')
		return
	end

	if ~(handles.absolute_motion)		% They were computed and are still in radians
		lon = lon / D2R;		lat = lat / D2R;
	end
	set(handles.edit_PoleLon,'String',sprintf('%3.2f',lon))
	set(handles.edit_PoleLat,'String',sprintf('%2.2f',lat))
	set(handles.edit_PoleRate,'String',sprintf('%1.4f',omega))
	push_Calculate_CB(hObject,handles,'nada')
	guidata(hObject, handles);

%--------------------------------------------------------------------------------------------------
function radio_Nuvel1A_CB(hObject, handles, tipo)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end

	if (nargin == 2)	tipo = 'Nuvel1A';	end
	set([handles.radio_Nuvel1A_NNR handles.radio_PBird handles.radio_DEOS2K ...
		handles.radio_REVEL handles.radio_GEODVEL],'Value',0)

	if (tipo(1) == 'N')		set(handles.radio_MORVEL,'Value',0)
	else					set(handles.radio_Nuvel1A,'Value',0)
	end
	set(handles.checkbox_Abs2Rel,'Visible','off')

	set(handles.popup_FixedPlate,'Enable','on')
	handles.absolute_motion = 0;				% The Nuvel1A is a relative motion model

	if (handles.first_MORVEL && tipo(1) == 'M')	% Load and read poles deffinition
		fid = fopen([handles.path_data 'MORVEL_poles.dat'],'r');
		[abbrev name lat lon omega] = strread(fread(fid,'*char'),'%s %s %f %f %f');
		fclose(fid);
		% Save the poles parameters in the handles structure
		handles.MORVEL_abbrev = abbrev;
		handles.MORVEL_name = name;
		handles.MORVEL_lat = lat;
		handles.MORVEL_lon = lon;
		handles.MORVEL_omega = omega;
		handles.MORVEL_comb = do_plate_comb('MORVEL');			% Based on PB because is the closest one
		handles.first_MORVEL = 0;
		setappdata(handles.figure1,'current_model','MORVEL')	% Need because Nuvel1A could be still active
	end

	% Fill the popupmenus with the Plate's names
	set([handles.popup_FixedPlate handles.popup_MovingPlate],'Value',1)
	if (tipo(1) == 'N')
		set([handles.popup_FixedPlate handles.popup_MovingPlate],'String',handles.Nuvel1A_name)
		handles.Nuvel1A_comb = do_plate_comb('Nuvel1A');
	else
		set([handles.popup_FixedPlate handles.popup_MovingPlate],'String',handles.MORVEL_name)
		handles.MORVEL_comb = do_plate_comb('MORVEL');			% Need change when we build the MORVEL boundaries
	end

	if (strncmp(tipo,'Nuvel1A',7))							% For both Nuvel1A and Nuvel1A_NNR
		set_Nuvel1Aplate_model(hObject,handles)
	else													% MORVEL
		set_PBplate_model(hObject,handles)					% Not mistake
	end

	% Clear the pole edit boxes fields
	set([handles.edit_PoleLon handles.edit_PoleLat handles.edit_PoleRate],'String','')

	% Remove info about the previously calculated velocity results
	set(handles.text_Speed,'String','Speed   = ');      set(handles.text_Azim,'String','Azimuth = ')

	% Flag in appdata which model is currently loaded
	setappdata(handles.figure1,'current_model',tipo)

	% Need to change the ButtonDownFcn call arguments
	h_patch = findobj('Type','patch');
	set(h_patch,'ButtonDownFcn',{@bdn_plate,handles,tipo})
	guidata(hObject, handles);

%--------------------------------------------------------------------------------------------------
function radio_Nuvel1A_NNR_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end

	D2R = pi/180;
	set([handles.radio_Nuvel1A handles.radio_PBird handles.radio_DEOS2K ...
		handles.radio_REVEL handles.radio_GEODVEL],'Value',0)
	set(handles.checkbox_Abs2Rel,'Visible','on')

	model = getappdata(handles.figure1,'current_model');
	if ~any(strcmp(model,{'Nuvel1A','Nuvel1A_NNR'}))    % Another plate model was loaded
		set_Nuvel1Aplate_model(hObject,handles)
	end

	if (handles.first_NNR)      % Load and read poles deffinition
		fid = fopen([handles.path_data 'Nuvel1A_NNR_poles.dat'],'r');
		[abbrev name lat lon omega] = strread(fread(fid,'*char'),'%s %s %f %f %f');
		fclose(fid);
		% Save the poles parameters in the handles structure
		handles.Nuvel1A_NNR_abbrev = abbrev;
		handles.Nuvel1A_NNR_name = name;
		handles.Nuvel1A_NNR_lat = lat;
		handles.Nuvel1A_NNR_lon = lon;
		handles.Nuvel1A_NNR_omega = omega;
		handles.Nuvel1A_NNR_comb = do_plate_comb('Nuvel1A_NNR');
		handles.first_NNR = 0;
	end

	% Fill the Moving plate popupmenu with the plate's names (we have to this in every case)
	set(handles.popup_MovingPlate,'Value',1,'String',handles.Nuvel1A_NNR_name)

	if (handles.abs2rel)                        % We are in "relativized absolute" motion mode
		set(handles.popup_FixedPlate,'Value',1)
		set(handles.popup_FixedPlate,'String',handles.Nuvel1A_NNR_name)
		set(handles.popup_FixedPlate,'Enable','on')
		lon1 = handles.Nuvel1A_NNR_lon(1);      lat1 = handles.Nuvel1A_NNR_lat(1);
		omega1 = handles.Nuvel1A_NNR_omega(1);
		ind = get(handles.popup_MovingPlate,'Value');
		lon2 = handles.Nuvel1A_NNR_lon(ind);    lat2 = handles.Nuvel1A_NNR_lat(ind);
		omega2 = handles.Nuvel1A_NNR_omega(ind);
		[lon,lat,omega] = calculate_pole(lon1,lat1,omega1,lon2,lat2,omega2);
		lon = lon/D2R;     lat = lat/D2R;
	else                                        % On the original absolute motion mode
		handles.absolute_motion = 1;
		set(handles.popup_FixedPlate,'Enable','off')
		lon = handles.Nuvel1A_NNR_lon(1);       lat = handles.Nuvel1A_NNR_lat(1);
		omega = handles.Nuvel1A_NNR_omega(1);
	end

	try         % Delete any existing pole representation
		delete(findobj('Tag','pole_out'));     delete(findobj('Tag','pole_in1'));
		delete(findobj('Tag','pole_in2'));
	end

	% Actualize the pole edit boxes fields and plot the pole
	if (omega ~= 0)         % That is, if the pole exists
		set(handles.edit_PoleLon,'String',num2str(lon,'%3.2f'))
		set(handles.edit_PoleLat,'String',num2str(lat,'%2.2f'))
		set(handles.edit_PoleRate,'String',num2str(omega,'%1.4f'))

		tmp = lon;
		if (tmp > 180),     tmp = tmp - 360; end
		line(tmp,lat,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k',...
			'MarkerSize',7,'Tag','pole_out');
		tmp = lon+180;
		if (tmp > 180),     tmp = tmp - 360;   end
		line(tmp,-lat,'Marker','+','MarkerFaceColor','k','MarkerEdgeColor','k',...
			'MarkerSize',8,'Tag','pole_in1');
		line(tmp,-lat,'Marker','o','MarkerEdgeColor','k','MarkerSize',8,'Tag','pole_in2');
	else
		set([handles.edit_PoleLon handles.edit_PoleLat handles.edit_PoleRate],'String','')
	end

	% Remove info about the previously calculated velocity results
	set(handles.text_Speed,'String','Speed   = ');      set(handles.text_Azim,'String','Azimuth = ')

	% Flag in appdata which model is currently loaded
	setappdata(handles.figure1,'current_model','NNR')

	% Need to change the ButtonDownFcn call arguments
	h_patch = findobj('Type','patch');
	set(h_patch,'ButtonDownFcn',{@bdn_plate,handles,'NNR'})
	guidata(hObject, handles);

%--------------------------------------------------------------------------------------------------
function radio_PBird_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end

	set([handles.radio_Nuvel1A handles.radio_Nuvel1A_NNR handles.radio_MORVEL ...
		handles.radio_DEOS2K handles.radio_REVEL handles.radio_GEODVEL],'Value',0)
	set(handles.checkbox_Abs2Rel,'Visible','off')

	set(handles.popup_FixedPlate,'Enable','on')
	handles.absolute_motion = 0;            % The PB is a relative motion model

	if (handles.first_PB)      % Load and read poles deffinition
		fid = fopen([handles.path_data 'PB_poles.dat'],'r');
		[abbrev name lat lon omega] = strread(fread(fid,'*char'),'%s %s %f %f %f');
		fclose(fid);
		% Save the poles parameters in the handles structure
		handles.PB_abbrev = abbrev;
		handles.PB_name = name;
		handles.PB_lat = lat;
		handles.PB_lon = lon;
		handles.PB_omega = omega;
		handles.PB_comb = do_plate_comb('PB');
		handles.first_PB = 0;
	end

	% Fill the popupmenus with the Plate's names
	set([handles.popup_FixedPlate handles.popup_MovingPlate],'Value',1)
	set([handles.popup_FixedPlate handles.popup_MovingPlate],'String',handles.PB_name)

	set_PBplate_model(hObject,handles)

	% Clear the pole edit boxes fields
	set([handles.edit_PoleLon handles.edit_PoleLat handles.edit_PoleRate],'String','')

	% Remove info about the previously calculated velocity results
	set(handles.text_Speed,'String','Speed   = ');      set(handles.text_Azim,'String','Azimuth = ')

	% Flag in appdata which model is currently loaded
	setappdata(handles.figure1,'current_model','PB')

	% Need to change the ButtonDownFcn call arguments
	h_patch = findobj('Type','patch');
	set(h_patch,'ButtonDownFcn',{@bdn_plate,handles,'PB'})
	guidata(hObject, handles);

%--------------------------------------------------------------------------------------------------
function radio_GEODVEL_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end

	set([handles.radio_Nuvel1A handles.radio_Nuvel1A_NNR handles.radio_MORVEL ...
		handles.radio_DEOS2K handles.radio_REVEL handles.radio_PBird],'Value',0)
	set(handles.checkbox_Abs2Rel,'Visible','off')

	set(handles.popup_FixedPlate,'Enable','on')
	handles.absolute_motion = 0;            % The PB is a relative motion model

	if (handles.first_GEODVEL)		% Load and read poles deffinition
		fid = fopen([handles.path_data 'GEODVEL_poles.dat'],'r');
		[abbrev name lat lon omega] = strread(fread(fid,'*char'),'%s %s %f %f %f');
		fclose(fid);
		% Save the poles parameters in the handles structure
		handles.GEODVEL_abbrev = abbrev;
		handles.GEODVEL_name = name;
		handles.GEODVEL_lat = lat;
		handles.GEODVEL_lon = lon;
		handles.GEODVEL_omega = omega;
		handles.GEODVEL_comb = do_plate_comb('GEODVEL');
		handles.first_GEODVEL = 0;
	end

	% Fill the popupmenus with the Plate's names
	set([handles.popup_FixedPlate handles.popup_MovingPlate],'Value',1)
	set([handles.popup_FixedPlate handles.popup_MovingPlate],'String',handles.GEODVEL_name)

	set_GEODVELplate_model(hObject,handles)

	% Clear the pole edit boxes fields
	set([handles.edit_PoleLon handles.edit_PoleLat handles.edit_PoleRate],'String','')

	% Remove info about the previously calculated velocity results
	set(handles.text_Speed,'String','Speed   = ');
	set(handles.text_Azim,'String','Azimuth = ')

	% Flag in appdata which model is currently loaded
	setappdata(handles.figure1,'current_model','GEODVEL')

	% Need to change the ButtonDownFcn call arguments
	h_patch = findobj('Type','patch');
	set(h_patch,'ButtonDownFcn',{@bdn_plate,handles,'GEODVEL'})
	guidata(hObject, handles);

%--------------------------------------------------------------------------------------------------
function radio_AKIM2000_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end

	D2R = pi/180;
	set([handles.radio_Nuvel1A handles.radio_Nuvel1A_NNR handles.radio_PBird ...
		handles.radio_REVEL handles.radio_GEODVEL],'Value',0)
	set(handles.checkbox_Abs2Rel,'Visible','on')

	set_AKIM2000plate_model(hObject,handles)

	if (handles.first_AKIM2000)      % Load and read poles deffinition
		fid = fopen([handles.path_data 'AKIM2000_poles.dat'],'r');
		[abbrev name lat lon omega] = strread(fread(fid,'*char'),'%s %s %f %f %f');
		fclose(fid);
		% Save the poles parameters in the handles structure
		handles.AKIM2000_abbrev = abbrev;
		handles.AKIM2000_name = name;
		handles.AKIM2000_lat = lat;
		handles.AKIM2000_lon = lon;
		handles.AKIM2000_omega = omega;
		handles.AKIM2000_comb = do_plate_comb('AKIM2000');
		handles.first_AKIM2000 = 0;
	end

	% Fill the Moving plate popupmenu with the plate's names (we have to this in every case)
	set(handles.popup_MovingPlate,'Value',1,'String',handles.AKIM2000_name)

	if (handles.abs2rel)                        % We are in "relativized absolute" motion mode
		set(handles.popup_FixedPlate,'Value',1)
		set(handles.popup_FixedPlate,'String',handles.AKIM2000_name)
		set(handles.popup_FixedPlate,'Enable','on')
		lon1 = handles.AKIM2000_lon(1);      lat1 = handles.AKIM2000_lat(1);
		omega1 = handles.AKIM2000_omega(1);
		ind = get(handles.popup_MovingPlate,'Value');
		lon2 = handles.AKIM2000_lon(ind);    lat2 = handles.AKIM2000_lat(ind);
		omega2 = handles.AKIM2000_omega(ind);
		[lon,lat,omega] = calculate_pole(lon1,lat1,omega1,lon2,lat2,omega2);
		lon = lon/D2R;     lat = lat/D2R;
	else                                        % On the original absolute motion mode
		handles.absolute_motion = 1;
		set(handles.popup_FixedPlate,'Enable','off')
		lon = handles.AKIM2000_lon(1);       lat = handles.AKIM2000_lat(1);
		omega = handles.AKIM2000_omega(1);
	end

	try         % Delete any existing pole representation
		delete(findobj('Tag','pole_out'));     delete(findobj('Tag','pole_in1'));
		delete(findobj('Tag','pole_in2'));
	end

	% Actualize the pole edit boxes fields and plot the pole
	if (omega ~= 0)         % That is, if the pole exists
		set(handles.edit_PoleLon,'String',num2str(lon,'%3.2f'))
		set(handles.edit_PoleLat,'String',num2str(lat,'%2.2f'))
		set(handles.edit_PoleRate,'String',num2str(omega,'%1.4f'))

		tmp = lon;
		if (tmp > 180),     tmp = tmp - 360; end
		line(tmp,lat,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k',...
			'MarkerSize',7,'Tag','pole_out');
		tmp = lon+180;
		if (tmp > 180),     tmp = tmp - 360;   end
		line(tmp,-lat,'Marker','+','MarkerFaceColor','k','MarkerEdgeColor','k',...
			'MarkerSize',8,'Tag','pole_in1');
		line(tmp,-lat,'Marker','o','MarkerEdgeColor','k','MarkerSize',8,'Tag','pole_in2');
	else
		set([handles.edit_PoleLon handles.edit_PoleLat handles.edit_PoleRate],'String','')
	end

	% Remove info about the previously calculated velocity results
	set(handles.text_Speed,'String','Speed   = ');      set(handles.text_Azim,'String','Azimuth = ')

	% Flag in appdata which model is currently loaded
	setappdata(handles.figure1,'current_model','AKIM2000')

	% Need to change the ButtonDownFcn call arguments
	h_patch = findobj('Type','patch');
	set(h_patch,'ButtonDownFcn',{@bdn_plate,handles,'AKIM2000'})
	guidata(hObject, handles);

%--------------------------------------------------------------------------------------------------
function radio_REVEL_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end

	D2R = pi/180;
	set([handles.radio_Nuvel1A handles.radio_Nuvel1A_NNR handles.radio_PBird ...
		handles.radio_MORVEL handles.radio_DEOS2K handles.radio_GEODVEL],'Value',0)
	set(handles.checkbox_Abs2Rel,'Visible','on')

	set_REVELplate_model(hObject,handles)

	if (handles.first_REVEL)      % Load and read poles deffinition
		fid = fopen([handles.path_data 'REVEL_poles.dat'],'r');
		[abbrev name lat lon omega] = strread(fread(fid,'*char'),'%s %s %f %f %f');
		fclose(fid);
		% Save the poles parameters in the handles structure
		handles.REVEL_abbrev = abbrev;
		handles.REVEL_name = name;
		handles.REVEL_lat = lat;
		handles.REVEL_lon = lon;
		handles.REVEL_omega = omega;
		handles.REVEL_comb = do_plate_comb('REVEL');
		handles.first_REVEL = 0;
	end

	% Fill the Moving plate popupmenu with the plate's names (we have to this in every case)
	set(handles.popup_MovingPlate,'Value',1,'String',handles.REVEL_name)

	if (handles.abs2rel)                        % We are in "relativized absolute" motion mode
		set(handles.popup_FixedPlate,'Value',1)
		set(handles.popup_FixedPlate,'String',handles.REVEL_name)
		set(handles.popup_FixedPlate,'Enable','on')
		lon1 = handles.REVEL_lon(1);      lat1 = handles.REVEL_lat(1);
		omega1 = handles.REVEL_omega(1);
		ind = get(handles.popup_MovingPlate,'Value');
		lon2 = handles.REVEL_lon(ind);    lat2 = handles.REVEL_lat(ind);
		omega2 = handles.REVEL_omega(ind);
		[lon,lat,omega] = calculate_pole(lon1,lat1,omega1,lon2,lat2,omega2);
		lon = lon/D2R;     lat = lat/D2R;
	else                                        % On the original absolute motion mode
		handles.absolute_motion = 1;
		set(handles.popup_FixedPlate,'Enable','off')
		lon = handles.REVEL_lon(1);       lat = handles.REVEL_lat(1);
		omega = handles.REVEL_omega(1);
	end

	try         % Delete any existing pole representation
		delete(findobj('Tag','pole_out'));     delete(findobj('Tag','pole_in1'));
		delete(findobj('Tag','pole_in2'));
	end

	% Actualize the pole edit boxes fields and plot the pole
	if (omega ~= 0)         % That is, if the pole exists
		set(handles.edit_PoleLon,'String',num2str(lon,'%3.2f'))
		set(handles.edit_PoleLat,'String',num2str(lat,'%2.2f'))
		set(handles.edit_PoleRate,'String',num2str(omega,'%1.4f'))

		tmp = lon;
		if (tmp > 180),     tmp = tmp - 360; end
		line(tmp,lat,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k',...
			'MarkerSize',7,'Tag','pole_out');
		tmp = lon+180;
		if (tmp > 180),     tmp = tmp - 360;   end
		line(tmp,-lat,'Marker','+','MarkerFaceColor','k','MarkerEdgeColor','k',...
			'MarkerSize',8,'Tag','pole_in1');
		line(tmp,-lat,'Marker','o','MarkerEdgeColor','k','MarkerSize',8,'Tag','pole_in2');
	else
		set([handles.edit_PoleLon handles.edit_PoleLat handles.edit_PoleRate],'String','')
	end

	% Remove info about the previously calculated velocity results
	set(handles.text_Speed,'String','Speed   = ');      set(handles.text_Azim,'String','Azimuth = ')

	% Flag in appdata which model is currently loaded
	setappdata(handles.figure1,'current_model','REVEL')

	% Need to change the ButtonDownFcn call arguments
	h_patch = findobj('Type','patch');
	set(h_patch,'ButtonDownFcn',{@bdn_plate,handles,'REVEL'})
	guidata(hObject, handles);

%--------------------------------------------------------------------------------------------------
function radio_DEOS2K_CB(hObject, handles)
	if ( ~get(hObject,'Val') ),		set(hObject,'Val',1),	return,		end

	D2R = pi/180;
	set([handles.radio_Nuvel1A handles.radio_Nuvel1A_NNR handles.radio_PBird ...
		handles.radio_MORVEL handles.radio_REVEL handles.radio_GEODVEL],'Value',0)
	set(handles.checkbox_Abs2Rel,'Visible','on')

	set_DEOS2Kplate_model(hObject,handles)

	if (handles.first_DEOS2K)      % Load and read poles deffinition
		fid = fopen([handles.path_data 'DEOS2K_poles.dat'],'r');
		[abbrev name lat lon omega] = strread(fread(fid,'*char'),'%s %s %f %f %f');
		fclose(fid);
		% Save the poles parameters in the handles structure
		handles.DEOS2K_abbrev = abbrev;
		handles.DEOS2K_name = name;
		handles.DEOS2K_lat = lat;
		handles.DEOS2K_lon = lon;
		handles.DEOS2K_omega = omega;
		handles.DEOS2K_comb = do_plate_comb('DEOS2K');
		handles.first_DEOS2K = 0;
	end

	% Fill the Moving plate popupmenu with the plate's names (we have to this in every case)
	set(handles.popup_MovingPlate,'Value',1,'String',handles.DEOS2K_name)

	if (handles.abs2rel)                        % We are in "relativized absolute" motion mode
		set(handles.popup_FixedPlate,'Value',1)
		set(handles.popup_FixedPlate,'String',handles.DEOS2K_name)
		set(handles.popup_FixedPlate,'Enable','on')
		lon1 = handles.DEOS2K_lon(1);      lat1 = handles.DEOS2K_lat(1);
		omega1 = handles.DEOS2K_omega(1);
		ind = get(handles.popup_MovingPlate,'Value');
		lon2 = handles.DEOS2K_lon(ind);    lat2 = handles.DEOS2K_lat(ind);
		omega2 = handles.DEOS2K_omega(ind);
		[lon,lat,omega] = calculate_pole(lon1,lat1,omega1,lon2,lat2,omega2);
		lon = lon/D2R;     lat = lat/D2R;
	else                                        % On the original absolute motion mode
		handles.absolute_motion = 1;
		set(handles.popup_FixedPlate,'Enable','off')
		lon = handles.DEOS2K_lon(1);       lat = handles.DEOS2K_lat(1);
		omega = handles.DEOS2K_omega(1);
	end

	try         % Delete any existing pole representation
		delete(findobj('Tag','pole_out'));     delete(findobj('Tag','pole_in1'));
		delete(findobj('Tag','pole_in2'));
	end

	% Actualize the pole edit boxes fields and plot the pole
	if (omega ~= 0)         % That is, if the pole exists
		set(handles.edit_PoleLon,'String',num2str(lon,'%3.2f'))
		set(handles.edit_PoleLat,'String',num2str(lat,'%2.2f'))
		set(handles.edit_PoleRate,'String',num2str(omega,'%1.4f'))

		tmp = lon;
		if (tmp > 180),     tmp = tmp - 360; end
		line(tmp,lat,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k',...
			'MarkerSize',7,'Tag','pole_out');
		tmp = lon+180;
		if (tmp > 180),     tmp = tmp - 360;   end
		line(tmp,-lat,'Marker','+','MarkerFaceColor','k','MarkerEdgeColor','k',...
			'MarkerSize',8,'Tag','pole_in1');
		line(tmp,-lat,'Marker','o','MarkerEdgeColor','k','MarkerSize',8,'Tag','pole_in2');
	else
		set([handles.edit_PoleLon handles.edit_PoleLat handles.edit_PoleRate],'String','')
	end

	% Remove info about the previously calculated velocity results
	set(handles.text_Speed,'String','Speed   = ');      set(handles.text_Azim,'String','Azimuth = ')

	% Flag in appdata which model is currently loaded
	setappdata(handles.figure1,'current_model','DEOS2K')

	% Need to change the ButtonDownFcn call arguments
	h_patch = findobj('Type','patch');
	set(h_patch,'ButtonDownFcn',{@bdn_plate,handles,'DEOS2K'})
	guidata(hObject, handles);

%--------------------------------------------------------------------------------------------------
function edit_PtLon_CB(hObject, handles)
	xx = get(hObject,'String');
	if (isempty(xx) || str2double(xx) < -180 || str2double(xx) > 360)
		set(hObject,'String','0')
	end

%--------------------------------------------------------------------------------------------------
function edit_PtLat_CB(hObject, handles)
	xx = get(hObject,'String');
	if (isempty(xx) || str2double(xx) < -90 || str2double(xx) > 90)
		set(hObject,'String','0')
	end

%--------------------------------------------------------------------------------------------------
function push_Calculate_CB(hObject, handles, opt)
% Calculate the Euler velocity based on the Euler pole read on respective edit boxes.

	if (nargin == 3),   opt = [];   end
	D2R = pi/180;
	earth_rad = 6371e3;    % Earth radius in km

	plon = str2double(get(handles.edit_PoleLon,'String'))*D2R;
	plat = str2double(get(handles.edit_PoleLat,'String'))*D2R;
	omega = str2double(get(handles.edit_PoleRate,'String'));
	if isnan(plon) || isnan(plat) || isnan(omega)
		errordlg('Euler Pole parameters are wrong.','Error');
		return
	end

	% Plot the pole using the arrow extremities analogy
	% A brief note. I tried to have these markers drawn only once and than changing their position
	% using set. However, no matter how many times I used uistack, after deleting the patches and
	% drawing new ones (following a model change), the markers never sohwn up again. So, the drastic
	% solution is to delete them and star over again (that is, drawing new ones).
	h_pole_out = findobj('Tag','pole_out');     h_pole_in1 = findobj('Tag','pole_in1');
	h_pole_in2 = findobj('Tag','pole_in2');
	delete(h_pole_out);     delete(h_pole_in1);     delete(h_pole_in2);
	tmp = plon;
	if (tmp > pi),		tmp = tmp - 2*pi;	end
	line(tmp/D2R,plat/D2R,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',7,'Tag','pole_out');
	tmp = plon+pi;
	if (tmp > pi),		tmp = tmp - 2*pi;	end
	line(tmp/D2R,-plat/D2R,'Marker','+','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,'Tag','pole_in1');
	line(tmp/D2R,-plat/D2R,'Marker','o','MarkerEdgeColor','k','MarkerSize',8,'Tag','pole_in2');

	% Get interest point coordinates
	alon = str2double(get(handles.edit_PtLon,'String'))*D2R;
	alat = str2double(get(handles.edit_PtLat,'String'))*D2R;
	if isnan(alon) || isnan(alat)
		if isempty(opt)         % User hit "Calculate" but no point where calculate is provided
			errordlg('Calculate velocity where?','Chico Clever');
		end
		return
	end
 	ecc = 0.0818191908426215;		% WGS84
 	alat = atan2( (1-ecc^2)*sin(alat), cos(alat) );			% Convert to geocentric

	x = cos(plat)*sin(plon)*sin(alat) - cos(alat)*sin(alon)*sin(plat);    % East vel
	y = cos(alat)*cos(alon)*sin(plat) - cos(plat)*cos(plon)*sin(alat);    % North vel
	z = cos(plat)*cos(alat)*sin(alon-plon);
	vlon = -sin(alon)*x + cos(alon)*y;
	vlat = -sin(alat)*cos(alon)*x-sin(alat)*sin(alon)*y + cos(alat)*z;
	azim = 90 - atan2(vlat,vlon) / D2R;

	if (azim < 0)       % Give allways the result in the 0-360 range
		azim = azim + 360;
	end

	x = sin(alat)*sin(plat) + cos(alat)*cos(plat)*cos(plon-alon);
	delta = acos(x);
	vel = omega*D2R/1e+3 * earth_rad * sin(delta);		% to give velocity in cm/Ma

	% Get the position of the text objects that will contain the velocity results
	set(handles.text_Speed,'String',['Speed   = ' sprintf('%.1f',vel) ' mm/yr'])
	set(handles.text_Azim,'String',['Azimuth = ' sprintf('%3.1f',azim) ' degree (cw from N)'])

% -----------------------------------------------------------------------------------------
function [vel,azim] = compute_velocity(alat,alon,plat,plon,omega)
% alat & alon are the point coords. plat, plon & omega are the pole parameters
	D2R = pi/180;
	earth_rad = 6371e3;    % Earth radius in km

	x = cos(plat)*sin(plon)*sin(alat) - cos(alat)*sin(alon)*sin(plat);    % East vel
	y = cos(alat)*cos(alon)*sin(plat) - cos(plat)*cos(plon)*sin(alat);    % North vel
	z = cos(plat)*cos(alat)*sin(alon-plon);
	vlon = -sin(alon)*x + cos(alon)*y;
	vlat = -sin(alat)*cos(alon)*x-sin(alat)*sin(alon)*y + cos(alat)*z;
	azim = 90 - atan2(vlat,vlon) / D2R;

	x = sin(alat)*sin(plat) + cos(alat)*cos(plat)*cos(plon-alon);
	delta = acos(x);
	vel = omega*D2R/1e+4 * earth_rad * sin(delta);      % to give velocity in cm/Ma


% -----------------------------------------------------------------------------------------
function bdn_plate(obj,eventdata,handles,opt)
% This is the ButtonDownFcn function that finds the selected point and from it, guesses
% the plate pairs involved in the movement.

	if (nargin == 3),   opt = [];   end

	D2R = pi/180;
	tag = get(gcbo,'Tag');
	pt = get(handles.axes1, 'CurrentPoint');
	set(handles.edit_PtLon,'String',pt(1,1))
	set(handles.edit_PtLat,'String',pt(1,2))

	if strcmp(opt,'Nuvel1A')
		mod_abb = handles.Nuvel1A_abbrev;
		mod_comb = handles.Nuvel1A_comb;
	elseif strcmp(opt,'NNR')
		mod_abb = handles.Nuvel1A_NNR_abbrev;
		mod_comb = handles.Nuvel1A_NNR_comb;
	elseif strcmp(opt,'PB')
		mod_abb = handles.PB_abbrev;
		mod_comb = handles.PB_comb;
	elseif strcmp(opt,'AKIM2000')
		mod_abb = handles.AKIM2000_abbrev;
		mod_comb = handles.AKIM2000_comb;
	elseif strcmp(opt,'REVEL')
		mod_abb = handles.REVEL_abbrev;
		mod_comb = handles.REVEL_comb;
	elseif strcmp(opt,'DEOS2K')
		mod_abb = handles.DEOS2K_abbrev;
		mod_comb = handles.DEOS2K_comb;
	elseif strcmp(opt,'MORVEL')
		if (strcmp(tag, 'AF'))		tag = 'NB';		end
		mod_abb = handles.MORVEL_abbrev;
		mod_comb = handles.MORVEL_comb;
	elseif strcmp(opt,'GEODVEL')
		mod_abb = handles.GEODVEL_abbrev;
		mod_comb = handles.GEODVEL_comb;
	end

	% Find (and set it on the popup) the moving plate name
	ind = find(strcmp(tag, mod_abb));

	if (isempty(ind))
		msg = {'Sorry, I don''t know a pole for this place.' ...
			'You will have to enter parameters manualy'};
		h(1) = text(-170,0,msg,'FontSize',15, 'FontName','bold');
		pause(3)
		delete(h);
		return
	end
	set(handles.popup_MovingPlate,'Value',ind)

	% If it is a absolute motion, call the computing function and return.
	handles = guidata(handles.figure1);		% Use updated version (in case something important changed)
	if (handles.absolute_motion)
		switch opt
			case 'NNR'
				lat = handles.Nuvel1A_NNR_lat(ind);     lon = handles.Nuvel1A_NNR_lon(ind);
				omega = handles.Nuvel1A_NNR_omega(ind);
				set(handles.edit_PoleLon,'String',sprintf('%3.2f', lon))
				set(handles.edit_PoleLat,'String',sprintf('%2.2f', lat))
				set(handles.edit_PoleRate,'String',sprintf('%1.4f', omega))
				push_Calculate_CB(obj,handles,'NNR')
				return
			case 'AKIM2000'
				lat = handles.AKIM2000_lat(ind);     lon = handles.AKIM2000_lon(ind);
				omega = handles.AKIM2000_omega(ind);
				set(handles.edit_PoleLon,'String',sprintf('%3.2f', lon))
				set(handles.edit_PoleLat,'String',sprintf('%2.2f', lat))
				set(handles.edit_PoleRate,'String',sprintf('%1.4f', omega))
				push_Calculate_CB(obj,handles,'AKIM2000')
				return
			case 'REVEL'
				lat = handles.REVEL_lat(ind);     lon = handles.REVEL_lon(ind);
				omega = handles.REVEL_omega(ind);
				set(handles.edit_PoleLon,'String',sprintf('%3.2f', lon))
				set(handles.edit_PoleLat,'String',sprintf('%2.2f', lat))
				set(handles.edit_PoleRate,'String',sprintf('%1.4f', omega))
				push_Calculate_CB(obj,handles,'REVEL')
				return
			case 'DEOS2K'
				lat = handles.DEOS2K_lat(ind);     lon = handles.DEOS2K_lon(ind);
				omega = handles.DEOS2K_omega(ind);
				set(handles.edit_PoleLon,'String',sprintf('%3.2f', lon))
				set(handles.edit_PoleLat,'String',sprintf('%2.2f', lat))
				set(handles.edit_PoleRate,'String',sprintf('%1.4f', omega))
				push_Calculate_CB(obj,handles,'DEOS2K')
				return
		end
	end

	% Find (and set it on the popup) the closest plate. This will be assume as the fix plate
	neigh = mod_comb{ind}(3:end);			% It starts at 3 because the first 2 have the abbrev of the moving plate
	dash = strfind(neigh,'-');
	c_tet = cos((90-pt(1,2))*D2R);			% cosinus of current point co-lat
	s_tet = sin((90-pt(1,2))*D2R);			% sinus of current point co-lat
	d_min = zeros(1,numel(dash));
	for (i = 1:numel(dash))					% Loop on all neighbour plates
		abb = neigh(dash(i)+1:dash(i)+2);
		h_pol = findobj(handles.axes1,'Tag',abb);
		x = get(h_pol,'XData');			y = get(h_pol,'YData');
		if (iscell(x))						% This occurs when a polygon wsa splited in two (Date-line jump)
			x = [x{1}; x{2}];
			y = [y{1}; y{2}];
		end
		dist = c_tet*cos((90-y)*D2R) + s_tet*sin((90-y)*D2R).*cos((pt(1,1)-x)*D2R);
		try
			d_min(i) = min(abs(acos(dist)));
		end
	end

	[c,i] = min(d_min);         % i will contain the index to closest plate in "neigh" string
	abb = neigh(dash(i)+1:dash(i)+2);
	ind = find(strcmp(abb, mod_abb));
	set(handles.popup_FixedPlate,'Value',ind)

	% If we reach here, it means we are dealing with a relative motion
	ind_mov = get(handles.popup_MovingPlate,'Value');
	ind_fix = get(handles.popup_FixedPlate,'Value');
	model = getappdata(handles.figure1,'current_model');
	switch model
		case 'Nuvel1A'
			lat2 = handles.Nuvel1A_lat(ind_mov);		lon2 = handles.Nuvel1A_lon(ind_mov);
			omega2 = handles.Nuvel1A_omega(ind_mov);	handles.abb_mov = handles.Nuvel1A_abbrev{ind_mov};
			lat1 = handles.Nuvel1A_lat(ind_fix);		lon1 = handles.Nuvel1A_lon(ind_fix);
			omega1 = handles.Nuvel1A_omega(ind_fix);	handles.abb_fix = handles.Nuvel1A_abbrev{ind_fix};
		case 'NNR'
			lat2 = handles.Nuvel1A_NNR_lat(ind_mov);	lon2 = handles.Nuvel1A_NNR_lon(ind_mov);
			omega2 = handles.Nuvel1A_NNR_omega(ind_mov);handles.abb_mov = handles.Nuvel1A_NNR_abbrev{ind_mov};
			lat1 = handles.Nuvel1A_NNR_lat(ind_fix);	lon1 = handles.Nuvel1A_NNR_lon(ind_fix);
			omega1 = handles.Nuvel1A_NNR_omega(ind_fix);handles.abb_fix = handles.Nuvel1A_NNR_abbrev{ind_fix};
		case 'PB'
			lat2 = handles.PB_lat(ind_mov);				lon2 = handles.PB_lon(ind_mov);
			omega2 = handles.PB_omega(ind_mov);			handles.abb_mov = handles.PB_abbrev{ind_mov};
			lat1 = handles.PB_lat(ind_fix);				lon1 = handles.PB_lon(ind_fix);
			omega1 = handles.PB_omega(ind_fix);			handles.abb_fix = handles.PB_abbrev{ind_fix};
		case 'AKIM2000'
			lat2 = handles.AKIM2000_lat(ind_mov);		lon2 = handles.AKIM2000_lon(ind_mov);
			omega2 = handles.AKIM2000_omega(ind_mov);	handles.abb_mov = handles.AKIM2000_abbrev{ind_mov};        
			lat1 = handles.AKIM2000_lat(ind_fix);		lon1 = handles.AKIM2000_lon(ind_fix);
			omega1 = handles.AKIM2000_omega(ind_fix);	handles.abb_fix = handles.AKIM2000_abbrev{ind_fix};        
		case 'REVEL'
			lat2 = handles.REVEL_lat(ind_mov);			lon2 = handles.REVEL_lon(ind_mov);
			omega2 = handles.REVEL_omega(ind_mov);		handles.abb_mov = handles.REVEL_abbrev{ind_mov};
			lat1 = handles.REVEL_lat(ind_fix);			lon1 = handles.REVEL_lon(ind_fix);
			omega1 = handles.REVEL_omega(ind_fix);		handles.abb_fix = handles.REVEL_abbrev{ind_fix};
		case 'DEOS2K'
			lat2 = handles.DEOS2K_lat(ind_mov);			lon2 = handles.DEOS2K_lon(ind_mov);
			omega2 = handles.DEOS2K_omega(ind_mov);		handles.abb_mov = handles.DEOS2K_abbrev{ind_mov};        
			lat1 = handles.DEOS2K_lat(ind_fix);			lon1 = handles.DEOS2K_lon(ind_fix);
			omega1 = handles.DEOS2K_omega(ind_fix);		handles.abb_fix = handles.DEOS2K_abbrev{ind_fix};        
		case 'MORVEL'
			try
				lat2 = handles.MORVEL_lat(ind_mov);		lon2 = handles.MORVEL_lon(ind_mov);
				omega2 = handles.MORVEL_omega(ind_mov);	handles.abb_mov = handles.MORVEL_abbrev{ind_mov};
				lat1 = handles.MORVEL_lat(ind_fix);		lon1 = handles.MORVEL_lon(ind_fix);
				omega1 = handles.MORVEL_omega(ind_fix);	handles.abb_fix = handles.MORVEL_abbrev{ind_fix};
			catch
				msg = {'Sorry, I don''t know a pole for this place.' ...
						'You will have to enter parameters manualy'};
				h = text(-170,0,msg,'FontSize',15, 'FontName','bold');
				pause(2),		delete(h);
				return
			end
		case 'GEODVEL'
			try
				lat2 = handles.GEODVEL_lat(ind_mov);	lon2 = handles.GEODVEL_lon(ind_mov);
				omega2 = handles.GEODVEL_omega(ind_mov);handles.abb_mov = handles.GEODVEL_abbrev{ind_mov};
				lat1 = handles.GEODVEL_lat(ind_fix);	lon1 = handles.GEODVEL_lon(ind_fix);
				omega1 = handles.GEODVEL_omega(ind_fix);handles.abb_fix = handles.GEODVEL_abbrev{ind_fix};
			catch
				msg = {'Sorry, I don''t know a pole for this place.' ...
						'You will have to enter parameters manualy'};
				h = text(-170,0,msg,'FontSize',15, 'FontName','bold');
				pause(2),		delete(h);
				return
			end
	end
	[lon,lat,omega] = calculate_pole(lon1,lat1,omega1,lon2,lat2,omega2);
	lon = lon/D2R;     lat = lat/D2R;

	if (omega == 0)     % This should not happen, but just in case
		set([handles.edit_PoleLon handles.edit_PoleLat handles.edit_PoleRate],'String','')
		return
	end

	set(handles.edit_PoleLon,'String',sprintf('%3.2f', lon))
	set(handles.edit_PoleLat,'String',sprintf('%2.2f', lat))
	set(handles.edit_PoleRate,'String',sprintf('%1.4f', omega))
	push_Calculate_CB(obj,handles,'nada')

% -----------------------------------------------------------------------------------------
function set_Nuvel1Aplate_model(hObject,handles)
% Draw the Nuvel1A model plates (as patches)
% The order by which the plates are drawn seams to be very important to give correct results
% I tested quite a lot, but I cannot guarantie that everything works correctly.
load([handles.path_data 'nuvel_polyg.mat'])

% Find and kill the patches of "other" plate model. I could just hide them, but I'm
% afraid that would noticebly slow down the program speed (well, everything does anyway).
h_patch = findobj('Type','patch');
if ~isempty(h_patch)    % h_patch is not empty when the user had previously selected another model
    delete(h_patch)
end

% The Nuvel1A was the first implemented model. The troubles found in drawing North_American and
% Antarctic plates where solved in a pedestrian way (by adding some lat = +/- points in the original
% polygon file). A more clever solution was followed in subsequent models (e.g. PB model), where
% those need points were added inside the program code.

ind = (North_Am(:,1) > 180);
North_Am(ind,1) = North_Am(ind,1) - 360;
patch(North_Am(:,1),North_Am(:,2),'b','FaceAlpha',0.5,'Tag','NA','ButtonDownFcn',{@bdn_plate,handles})

ind = (Antarctic(:,1) > 180);
Antarctic(ind,1) = Antarctic(ind,1) - 360;
patch(Antarctic(:,1),Antarctic(:,2),'y','FaceAlpha',0.5,'Tag','AN','ButtonDownFcn',{@bdn_plate,handles})

patch(Nazca(:,1)-360,Nazca(:,2),'c','FaceAlpha',0.5,'Tag','NZ','ButtonDownFcn',{@bdn_plate,handles})
patch(Juan(:,1)-360,Juan(:,2),'r','FaceAlpha',0.5,'Tag','JF','ButtonDownFcn',{@bdn_plate,handles})
patch(Cocos(:,1)-360,Cocos(:,2),'r','FaceAlpha',0.5,'Tag','CO','ButtonDownFcn',{@bdn_plate,handles})
patch(Scotia(:,1)-360,Scotia(:,2),'g','FaceAlpha',0.5,'Tag','SA','ButtonDownFcn',{@bdn_plate,handles})
patch(Caribbe(:,1)-360,Caribbe(:,2),'m','FaceAlpha',0.5,'Tag','CA','ButtonDownFcn',{@bdn_plate,handles})

% Those cross the "Date-line" and need special care
ind = (Australia(:,1) <= 180);
patch(Australia(ind,1),Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})
ind = (Australia(:,1) > 180);
patch(Australia(ind,1)-360,Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})

ind = (Africa(:,1) > 180);
Africa(ind,1) = Africa(ind,1) - 360;
patch(Africa(:,1),Africa(:,2),'r','FaceAlpha',0.5,'Tag','AF','ButtonDownFcn',{@bdn_plate,handles})

ind = (Eurasia(:,1) > 180);
Eurasia(ind,1) = Eurasia(ind,1) - 360;
patch(Eurasia(:,1),Eurasia(:,2),'g','FaceAlpha',0.5,'Tag','EU','ButtonDownFcn',{@bdn_plate,handles})

ind = (Pacific(:,1) <= 180);
patch(Pacific(ind,1),Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})
ind = (Pacific(:,1) > 180);
patch(Pacific(ind,1)-360,Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})

% These come here because there seams to exist problems with their neighbouring patches.
% That is revealed by the fact that a click on them says that we are at another plate.
patch(South_Am(:,1)-360,South_Am(:,2),'g','FaceAlpha',0.5,'Tag','SA','ButtonDownFcn',{@bdn_plate,handles});
patch(Philippine(:,1),Philippine(:,2),'r','FaceAlpha',0.5,'Tag','PS','ButtonDownFcn',{@bdn_plate,handles})
patch(Arabia(:,1),Arabia(:,2),'y','FaceAlpha',0.5,'Tag','AR','ButtonDownFcn',{@bdn_plate,handles})
patch(India(:,1),India(:,2),'c','FaceAlpha',0.5,'Tag','IN','ButtonDownFcn',{@bdn_plate,handles})

% -----------------------------------------------------------------------------------------
function set_PBplate_model(hObject,handles)
% Draw the PB model plates  (as patches)
load([handles.path_data 'PB_polyg.mat'])

% Find and kill the patches of "other" plate model. I could just hide them, but I'm
% afraid that would noticebly slow down the program operations.
h_patch = findobj('Type','patch');
delete(h_patch)

% I have to put them here otherwise there will be many failures. That is, clicking on
% a plate (tipicaly on a small plate) near to those gives a wrong result. These clearly
% smells a patch-gcbo bug.
Antarctic = [179.999 -90; Antarctic; -179.999 -90];
patch(Antarctic(:,1),Antarctic(:,2),'y','FaceAlpha',0.5,'Tag','AN','ButtonDownFcn',{@bdn_plate,handles})

North_Am = [-179.999 90; North_Am; 179.999 90];
patch(North_Am(:,1),North_Am(:,2),'b','FaceAlpha',0.5,'Tag','NA','ButtonDownFcn',{@bdn_plate,handles})

ind = (Australia(:,1) <= 180 & Australia(:,1) >= 0);
patch(Australia(ind,1),Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})
ind = (Australia(:,1) > -180 & Australia(:,1) < 0);
patch(Australia(ind,1),Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})

ind = (Pacific(:,1) <= 180 & Pacific(:,1) >= 0);
patch(Pacific(ind,1),Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})
ind = (Pacific(:,1) > -180 & Pacific(:,1) < 0);
patch(Pacific(ind,1),Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})

patch(Africa(:,1),Africa(:,2),'r','FaceAlpha',0.5,'Tag','AF','ButtonDownFcn',{@bdn_plate,handles})
patch(Amur(:,1),Amur(:,2),'y','FaceAlpha',0.5,'Tag','AM','ButtonDownFcn',{@bdn_plate,handles})
% patch(Arabia(:,1),Arabia(:,2),'y','FaceAlpha',0.5,'Tag','AR','ButtonDownFcn',{@bdn_plate,handles})
% patch(Agean_Sea(:,1),Agean_Sea(:,2),'b','FaceAlpha',0.5,'Tag','AS','ButtonDownFcn',{@bdn_plate,handles})
% patch(Anatolia(:,1),Anatolia(:,2),'m','FaceAlpha',0.5,'Tag','AT','ButtonDownFcn',{@bdn_plate,handles})
patch(Birds_Head(:,1),Birds_Head(:,2),'m','FaceAlpha',0.5,'Tag','BH','ButtonDownFcn',{@bdn_plate,handles})
patch(Banda_Sea(:,1),Banda_Sea(:,2),'m','FaceAlpha',0.5,'Tag','BS','ButtonDownFcn',{@bdn_plate,handles})
patch(Burma(:,1),Burma(:,2),'m','FaceAlpha',0.5,'Tag','BU','ButtonDownFcn',{@bdn_plate,handles})
patch(Caribbe(:,1),Caribbe(:,2),'m','FaceAlpha',0.5,'Tag','CA','ButtonDownFcn',{@bdn_plate,handles})
patch(Caroline(:,1),Caroline(:,2),'g','FaceAlpha',0.5,'Tag','CL','ButtonDownFcn',{@bdn_plate,handles})
patch(Cocos(:,1),Cocos(:,2),'r','FaceAlpha',0.5,'Tag','CO','ButtonDownFcn',{@bdn_plate,handles})
patch(Conway_Reef(:,1),Conway_Reef(:,2),'r','FaceAlpha',0.5,'Tag','CR','ButtonDownFcn',{@bdn_plate,handles})
patch(Easter(:,1),Easter(:,2),'r','FaceAlpha',0.5,'Tag','EA','ButtonDownFcn',{@bdn_plate,handles})
patch(Eurasia(:,1),Eurasia(:,2),'g','FaceAlpha',0.5,'Tag','EU','ButtonDownFcn',{@bdn_plate,handles})
% Here is a clear example on the patch problem (bug I would say). If the declaration of these
% patches is donne above (I let them there as comented lines) instead of here, clicking on these
% plates (patches) gives wrong results. However, if I let them here, they work correctly.
patch(Arabia(:,1),Arabia(:,2),'y','FaceAlpha',0.5,'Tag','AR','ButtonDownFcn',{@bdn_plate,handles})
patch(Agean_Sea(:,1),Agean_Sea(:,2),'b','FaceAlpha',0.5,'Tag','AS','ButtonDownFcn',{@bdn_plate,handles})
patch(Anatolia(:,1),Anatolia(:,2),'m','FaceAlpha',0.5,'Tag','AT','ButtonDownFcn',{@bdn_plate,handles})

patch(Futuna(:,1),Futuna(:,2),'g','FaceAlpha',0.5,'Tag','FT','ButtonDownFcn',{@bdn_plate,handles})
patch(Galapagos(:,1),Galapagos(:,2),'r','FaceAlpha',0.5,'Tag','GP','ButtonDownFcn',{@bdn_plate,handles})
patch(Juan(:,1),Juan(:,2),'r','FaceAlpha',0.5,'Tag','JF','ButtonDownFcn',{@bdn_plate,handles})
patch(Juan_Fernandez(:,1),Juan_Fernandez(:,2),'r','FaceAlpha',0.5,'Tag','JZ','ButtonDownFcn',{@bdn_plate,handles})
patch(India(:,1),India(:,2),'c','FaceAlpha',0.5,'Tag','IN','ButtonDownFcn',{@bdn_plate,handles})
patch(Mariana(:,1),Mariana(:,2),'c','FaceAlpha',0.5,'Tag','MA','ButtonDownFcn',{@bdn_plate,handles})
patch(Manus(:,1),Manus(:,2),'c','FaceAlpha',0.5,'Tag','MN','ButtonDownFcn',{@bdn_plate,handles})
patch(Maoke(:,1),Maoke(:,2),'c','FaceAlpha',0.5,'Tag','MO','ButtonDownFcn',{@bdn_plate,handles})
patch(Molluca_Sea(:,1),Molluca_Sea(:,2),'c','FaceAlpha',0.5,'Tag','MS','ButtonDownFcn',{@bdn_plate,handles})
patch(Nazca(:,1),Nazca(:,2),'c','FaceAlpha',0.5,'Tag','NZ','ButtonDownFcn',{@bdn_plate,handles})
patch(Niaofoou(:,1),Niaofoou(:,2),'c','FaceAlpha',0.5,'Tag','NI','ButtonDownFcn',{@bdn_plate,handles})
patch(North_Andes(:,1),North_Andes(:,2),'b','FaceAlpha',0.5,'Tag','ND','ButtonDownFcn',{@bdn_plate,handles})
patch(North_Bismarck(:,1),North_Bismarck(:,2),'b','FaceAlpha',0.5,'Tag','NB','ButtonDownFcn',{@bdn_plate,handles})
patch(New_Hebrides(:,1),New_Hebrides(:,2),'g','FaceAlpha',0.5,'Tag','NH','ButtonDownFcn',{@bdn_plate,handles})
patch(Okhotsk(:,1),Okhotsk(:,2),'c','FaceAlpha',0.5,'Tag','OK','ButtonDownFcn',{@bdn_plate,handles})
patch(Okinawa(:,1),Okinawa(:,2),'c','FaceAlpha',0.5,'Tag','ON','ButtonDownFcn',{@bdn_plate,handles})
patch(Panama(:,1),Panama(:,2),'r','FaceAlpha',0.5,'Tag','PM','ButtonDownFcn',{@bdn_plate,handles})
patch(Philippine(:,1),Philippine(:,2),'r','FaceAlpha',0.5,'Tag','PS','ButtonDownFcn',{@bdn_plate,handles})
patch(Rivera(:,1),Rivera(:,2),'r','FaceAlpha',0.5,'Tag','RI','ButtonDownFcn',{@bdn_plate,handles})
patch(South_Am(:,1),South_Am(:,2),'g','FaceAlpha',0.5,'Tag','SA','ButtonDownFcn',{@bdn_plate,handles});
patch(Altiplano(:,1),Altiplano(:,2),'r','FaceAlpha',0.5,'Tag','AP','ButtonDownFcn',{@bdn_plate,handles})
patch(South_Bismarck(:,1),South_Bismarck(:,2),'g','FaceAlpha',0.5,'Tag','SB','ButtonDownFcn',{@bdn_plate,handles});
patch(Scotia(:,1),Scotia(:,2),'b','FaceAlpha',0.5,'Tag','SC','ButtonDownFcn',{@bdn_plate,handles})
patch(Shetland(:,1),Shetland(:,2),'r','FaceAlpha',0.5,'Tag','SL','ButtonDownFcn',{@bdn_plate,handles})
patch(Solomon_Sea(:,1),Solomon_Sea(:,2),'m','FaceAlpha',0.5,'Tag','SS','ButtonDownFcn',{@bdn_plate,handles});
patch(Somalia(:,1),Somalia(:,2),'m','FaceAlpha',0.5,'Tag','SO','ButtonDownFcn',{@bdn_plate,handles});
patch(Sunda(:,1),Sunda(:,2),'y','FaceAlpha',0.5,'Tag','SU','ButtonDownFcn',{@bdn_plate,handles});
patch(Sandwich(:,1),Sandwich(:,2),'r','FaceAlpha',0.5,'Tag','SW','ButtonDownFcn',{@bdn_plate,handles})
patch(Timor(:,1),Timor(:,2),'g','FaceAlpha',0.5,'Tag','TI','ButtonDownFcn',{@bdn_plate,handles})
patch(Tonga(:,1),Tonga(:,2),'g','FaceAlpha',0.5,'Tag','TO','ButtonDownFcn',{@bdn_plate,handles})
patch(Woodlark(:,1),Woodlark(:,2),'r','FaceAlpha',0.5,'Tag','WL','ButtonDownFcn',{@bdn_plate,handles})
patch(Yangtze(:,1),Yangtze(:,2),'b','FaceAlpha',0.5,'Tag','YZ','ButtonDownFcn',{@bdn_plate,handles})

ind = (Balmoral_Reef(:,1) <= 180 & Balmoral_Reef(:,1) > 0);
patch(Balmoral_Reef(ind,1),Balmoral_Reef(ind,2),'c','FaceAlpha',0.5,'Tag','BR','ButtonDownFcn',{@bdn_plate,handles})
ind = (Balmoral_Reef(:,1) > -180 & Balmoral_Reef(:,1) < 0);
patch(Balmoral_Reef(ind,1),Balmoral_Reef(ind,2),'c','FaceAlpha',0.5,'Tag','BR','ButtonDownFcn',{@bdn_plate,handles})

ind = (Kermadec(:,1) <= 180 & Kermadec(:,1) >= 0);
patch(Kermadec(ind,1),Kermadec(ind,2),'r','FaceAlpha',0.5,'Tag','KE','ButtonDownFcn',{@bdn_plate,handles})
ind = (Kermadec(:,1) > -180 & Kermadec(:,1) < 0);
patch(Kermadec(ind,1),Kermadec(ind,2),'r','FaceAlpha',0.5,'Tag','KE','ButtonDownFcn',{@bdn_plate,handles})

% -----------------------------------------------------------------------------------------
function set_GEODVELplate_model(hObject,handles)
% Draw the DEOS2K model plates (as patches)
	load([handles.path_data 'PB_polyg.mat'])     % use PB polygons because ... it has them all

	% Find and kill the patches of "other" plate model.
	delete(findobj(handles.axes1,'Type','patch'))

	patch(Eurasia(:,1),Eurasia(:,2),'g','FaceAlpha',0.5,'Tag','EU','ButtonDownFcn',{@bdn_plate,handles})

	Antarctic = [179.999 -90; Antarctic; -179.999 -90];
	patch(Antarctic(:,1),Antarctic(:,2),'y','FaceAlpha',0.5,'Tag','AN','ButtonDownFcn',{@bdn_plate,handles})

	North_Am = [-179.999 90; North_Am; 179.999 90];
	patch(North_Am(:,1),North_Am(:,2),'b','FaceAlpha',0.5,'Tag','NA','ButtonDownFcn',{@bdn_plate,handles})

	ind = (Australia(:,1) <= 180 & Australia(:,1) >= 0);
	patch(Australia(ind,1),Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})
	ind = (Australia(:,1) > -180 & Australia(:,1) < 0);
	patch(Australia(ind,1),Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})

	ind = (Pacific(:,1) <= 180 & Pacific(:,1) >= 0);
	patch(Pacific(ind,1),Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})
	ind = (Pacific(:,1) > -180 & Pacific(:,1) < 0);
	patch(Pacific(ind,1),Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})

	patch(Africa(:,1),Africa(:,2),'r','FaceAlpha',0.5,'Tag','NU','ButtonDownFcn',{@bdn_plate,handles})
	patch(Arabia(:,1),Arabia(:,2),'y','FaceAlpha',0.5,'Tag','AR','ButtonDownFcn',{@bdn_plate,handles})
	patch(India(:,1),India(:,2),'c','FaceAlpha',0.5,'Tag','IN','ButtonDownFcn',{@bdn_plate,handles})
	patch(Nazca(:,1),Nazca(:,2),'c','FaceAlpha',0.5,'Tag','NZ','ButtonDownFcn',{@bdn_plate,handles})
	patch(South_Am(:,1),South_Am(:,2),'g','FaceAlpha',0.5,'Tag','SA','ButtonDownFcn',{@bdn_plate,handles});
	patch(Somalia(:,1),Somalia(:,2),'m','FaceAlpha',0.5,'Tag','SO','ButtonDownFcn',{@bdn_plate,handles});
	
% -----------------------------------------------------------------------------------------
function set_AKIM2000plate_model(hObject,handles)
% Draw the APKIM2000 model plates (as patches)
load([handles.path_data 'AKIM2000_polyg.mat'])

% Find and kill the patches of "other" plate model.
h_patch = findobj('Type','patch');
delete(h_patch)

ind=find(North_Am(:,1)>180);
North_Am(ind,1) = North_Am(ind,1) - 360;
patch(North_Am(:,1),North_Am(:,2),'b','FaceAlpha',0.5,'Tag','NA','ButtonDownFcn',{@bdn_plate,handles})

ind=find(Antarctic(:,1)>180);
Antarctic(ind,1) = Antarctic(ind,1) - 360;
patch(Antarctic(:,1),Antarctic(:,2),'y','FaceAlpha',0.5,'Tag','AN','ButtonDownFcn',{@bdn_plate,handles})

patch(Nazca(:,1)-360,Nazca(:,2),'c','FaceAlpha',0.5,'Tag','NZ','ButtonDownFcn',{@bdn_plate,handles})
patch(Caribbe(:,1)-360,Caribbe(:,2),'m','FaceAlpha',0.5,'Tag','CA','ButtonDownFcn',{@bdn_plate,handles})

% Those cross the "Date-line" and need special care
ind=find(Australia(:,1)<=180);
patch(Australia(ind,1),Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})
ind=find(Australia(:,1)>180);
patch(Australia(ind,1)-360,Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})

patch(Africa(:,1),Africa(:,2),'r','FaceAlpha',0.5,'Tag','AF','ButtonDownFcn',{@bdn_plate,handles})
patch(Somalia(:,1),Somalia(:,2),'m','FaceAlpha',0.5,'Tag','SO','ButtonDownFcn',{@bdn_plate,handles})

ind=find(Eurasia(:,1)>180);
Eurasia(ind,1) = Eurasia(ind,1) - 360;
patch(Eurasia(:,1),Eurasia(:,2),'g','FaceAlpha',0.5,'Tag','EU','ButtonDownFcn',{@bdn_plate,handles})

ind=find(Pacific(:,1)<=180);
patch(Pacific(ind,1),Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})
ind=find(Pacific(:,1)>180);
patch(Pacific(ind,1)-360,Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})

patch(South_Am(:,1)-360,South_Am(:,2),'g','FaceAlpha',0.5,'Tag','SA','ButtonDownFcn',{@bdn_plate,handles});
patch(Arabia(:,1),Arabia(:,2),'y','FaceAlpha',0.5,'Tag','AR','ButtonDownFcn',{@bdn_plate,handles})

% -----------------------------------------------------------------------------------------
function set_REVELplate_model(hObject,handles)
% Draw the REVEL model plates (as patches)
load([handles.path_data 'REVEL_polyg.mat'])

% Find and kill the patches of "other" plate model.
h_patch = findobj('Type','patch');
delete(h_patch)

Antarctic = [179.999 -90; Antarctic; -179.999 -90];
patch(Antarctic(:,1),Antarctic(:,2),'y','FaceAlpha',0.5,'Tag','AN','ButtonDownFcn',{@bdn_plate,handles})

North_Am = [-179.999 90; North_Am; 179.999 90];
patch(North_Am(:,1),North_Am(:,2),'b','FaceAlpha',0.5,'Tag','NA','ButtonDownFcn',{@bdn_plate,handles})

ind = (Australia(:,1) <= 180 & Australia(:,1) >= 0);
patch(Australia(ind,1),Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})
ind = (Australia(:,1) > -180 & Australia(:,1) < 0);
patch(Australia(ind,1),Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})

ind = (Pacific(:,1) <= 180 & Pacific(:,1) >= 0);
patch(Pacific(ind,1),Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})
ind = (Pacific(:,1) > -180 & Pacific(:,1) < 0);
patch(Pacific(ind,1),Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})

patch(Nubia(:,1),Nubia(:,2),'r','FaceAlpha',0.5,'Tag','NU','ButtonDownFcn',{@bdn_plate,handles})
patch(Amuria(:,1),Amuria(:,2),'y','FaceAlpha',0.5,'Tag','AM','ButtonDownFcn',{@bdn_plate,handles})
patch(Caribbe(:,1),Caribbe(:,2),'m','FaceAlpha',0.5,'Tag','CA','ButtonDownFcn',{@bdn_plate,handles})
patch(Eurasia(:,1),Eurasia(:,2),'g','FaceAlpha',0.5,'Tag','EU','ButtonDownFcn',{@bdn_plate,handles})
patch(Arabia(:,1),Arabia(:,2),'y','FaceAlpha',0.5,'Tag','AR','ButtonDownFcn',{@bdn_plate,handles})
patch(Anatolia(:,1),Anatolia(:,2),'m','FaceAlpha',0.5,'Tag','AT','ButtonDownFcn',{@bdn_plate,handles})
patch(India(:,1),India(:,2),'c','FaceAlpha',0.5,'Tag','IN','ButtonDownFcn',{@bdn_plate,handles})
patch(Nazca(:,1),Nazca(:,2),'c','FaceAlpha',0.5,'Tag','NZ','ButtonDownFcn',{@bdn_plate,handles})
patch(Okhotsk(:,1),Okhotsk(:,2),'c','FaceAlpha',0.5,'Tag','OK','ButtonDownFcn',{@bdn_plate,handles})
patch(Philippine(:,1),Philippine(:,2),'r','FaceAlpha',0.5,'Tag','PS','ButtonDownFcn',{@bdn_plate,handles})
patch(South_Am(:,1)-360,South_Am(:,2),'g','FaceAlpha',0.5,'Tag','SA','ButtonDownFcn',{@bdn_plate,handles});
patch(Somalia(:,1),Somalia(:,2),'m','FaceAlpha',0.5,'Tag','SO','ButtonDownFcn',{@bdn_plate,handles});
patch(Sunda(:,1),Sunda(:,2),'y','FaceAlpha',0.5,'Tag','SU','ButtonDownFcn',{@bdn_plate,handles});

% -----------------------------------------------------------------------------------------
function set_DEOS2Kplate_model(hObject,handles)
% Draw the DEOS2K model plates (as patches)
	load([handles.path_data 'REVEL_polyg.mat'])     % use REVEL polygons because I don't know beter

	% Find and kill the patches of "other" plate model.
	h_patch = findobj('Type','patch');
	delete(h_patch)

	North_Am = [-179.999 90; North_Am; 179.999 90];
	patch(North_Am(:,1),North_Am(:,2),'b','FaceAlpha',0.5,'Tag','NA','ButtonDownFcn',{@bdn_plate,handles})

	Antarctic = [179.999 -90; Antarctic; -179.999 -90];
	patch(Antarctic(:,1),Antarctic(:,2),'y','FaceAlpha',0.5,'Tag','AN','ButtonDownFcn',{@bdn_plate,handles})

	patch(Nubia(:,1),Nubia(:,2),'r','FaceAlpha',0.5,'Tag','NU','ButtonDownFcn',{@bdn_plate,handles})
	patch(Somalia(:,1),Somalia(:,2),'m','FaceAlpha',0.5,'Tag','SO','ButtonDownFcn',{@bdn_plate,handles})
	patch(South_Am(:,1)-360,South_Am(:,2),'g','FaceAlpha',0.5,'Tag','SA','ButtonDownFcn',{@bdn_plate,handles});

	ind = (Eurasia(:,1) > 180);
	Eurasia(ind,1) = Eurasia(ind,1) - 360;
	patch(Eurasia(:,1),Eurasia(:,2),'g','FaceAlpha',0.5,'Tag','EU','ButtonDownFcn',{@bdn_plate,handles})

	% These cross the "Date-line" and need special care
	ind = (Pacific(:,1) <= 180 & Pacific(:,1) >= 0);
	patch(Pacific(ind,1),Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})
	ind = (Pacific(:,1) > -180 & Pacific(:,1) <0 );
	patch(Pacific(ind,1),Pacific(ind,2),'m','FaceAlpha',0.5,'Tag','PA','ButtonDownFcn',{@bdn_plate,handles})

	ind = (Australia(:,1) <= 180 & Australia(:,1) >= 0);
	patch(Australia(ind,1),Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})
	ind = (Australia(:,1) > -180 & Australia(:,1) < 0);
	patch(Australia(ind,1),Australia(ind,2),'b','FaceAlpha',0.5,'Tag','AU','ButtonDownFcn',{@bdn_plate,handles})

% -----------------------------------------------------------------------------------------
function out = do_plate_comb(which)
% These function creates strings describing the neighbour plates to the first 
% plate abbreviation in the char string.
switch which
	case 'Nuvel1A'
		out(1) = {'AF-AN-AR-EU-NA-SA'};%-AN-AR-AS-AT-EU-NA-SA-SO
		out(2) = {'AN-AF-AU-NZ-PA'};
		out(3) = {'AR-AF-EU-IN'};
		out(4) = {'AU-AF-AN-EU-IN-PA'};
		out(5) = {'CA-CO-NA-NZ-SA'};
		out(6) = {'CO-CA-NA-NZ-PA'};
		out(7) = {'EU-AF-AR-AU-IN-NA-PA-PS'};
		out(8) = {'IN-AF-AU-AR-EU'};
		out(9) = {'JF-NA-PA'};
		out(10) = {'NA-AF-CA-CO-EU-JF-PA-SA'};
		out(11) = {'NZ-AN-CA-CO-PA-SA'};
		out(12) = {'PA-AN-AU-CO-EU-JF-NA-NZ-PS'};
		out(13) = {'PS-AU-EU-PA'};
		out(14) = {'SA-AF-AN-CA-NA-NZ'};
	case 'Nuvel1A_NNR'
		out(1) = {'AF-AN-AR-EU-NA-SA'};
		out(2) = {'AN-AF-AU-NZ-PA-SC'};
		out(3) = {'AR-AF-EU-IN'};
		out(4) = {'AU-AF-AN-EU-IN-PA'};
		out(5) = {'CA-CO-NA-NZ-SA'};
		out(6) = {'CO-CA-NA-NZ-PA'};
		out(7) = {'EU-AF-AR-AU-IN-NA-PA-PS'};
		out(8) = {'IN-AF-AU-AR-EU'};
		out(9) = {'JF-NA-PA'};
		out(10) = {'NA-AF-CA-CO-EU-JF-PA-SA'};
		out(11) = {'NZ-AN-CA-CO-PA-SA'};
		out(12) = {'PA-AN-AU-CO-EU-JF-NA-NZ-PS'};
		out(13) = {'PS-AU-EU-PA'};
		out(14) = {'SA-AF-AN-CA-NA-NZ-SC'};
		out(15) = {'RI-'};
		out(16) = {'SC-AN-SA'};
	case 'MORVEL'					% FAKE THING BASED ON PB MODEL (TEMPORARY)
		out(1) = {'AM-EU-PS-YZ'};
		out(2) = {'AN-NB-AU-JZ-NZ-PA-SA-SC-SO-SW'};
		out(3) = {'AR-NB-EU-IN-SO'};
		out(4) = {'AU-AN-IN-PA-SO-SU'};
		out(5) = {'CA-CO-NA-SA'};
		out(6) = {'CO-CA-NA-NZ-PA-RI'};
		out(7) = {'CP-'};		% We still don't have this
		out(8) = {'EU-NB-AM-AR-IN-NA-SU-YZ'};
		out(9) = {'IN-AR-AU-EU-SO'};
		out(10) = {'JF-NA-PA'};
		out(11) = {'LW-'};		% We still don't have this
		out(12) = {'MQ-'};		% We still don't have this
		out(13) = {'NA-NB-CA-CO-EU-JF-PA-RI-SA'};
		out(14) = {'NB-AN-AR-EU-NA-SA-SO'};
		out(15) = {'NZ-AN-CO-JZ-PA-SA'};
		out(16) = {'PA-AN-AU-CO-JF-JZ-NA-NZ-PS-RI'};
		out(17) = {'PS-AM-PA-SU-YZ'};
		out(18) = {'RI-CO-NA-PA'};
		out(19) = {'SA-NB-AN-AP-CA-NA-NZ-SC-SW'};
		out(20) = {'SC-AN-SA-SW'};
		out(21) = {'SO-NB-AN-AR-AU-IN'};
		out(22) = {'SR-'};		% We still don't have this
		out(23) = {'SU-AU-EU-PS'};
		out(24) = {'SW-AN-SA-SC'};
		out(25) = {'YZ-AM-EU-PS'};
	case 'AKIM2000'
		out(1) = {'AF-AN-AR-EU-NA-SA-SO'};
		out(2) = {'AN-AF-AU-NZ-PA'};
		out(3) = {'AR-AF-EU'};
		out(4) = {'AU-AF-AN-EU-PA'};
		out(5) = {'CA-NA-NZ-SA'};
		out(6) = {'EU-AF-AR-AU-NA-PA'};
		out(7) = {'NA-AF-CA-EU-PA-SA'};
		out(8) = {'NZ-AN-CA-PA-SA'};
		out(9) = {'PA-AN-AU-EU-NA-NZ'};
		out(10) = {'SA-AF-AN-CA-NA-NZ'};
		out(11) = {'AS-'};
		out(12) = {'SO-AF-AN-AU-AR'};
	case 'PB'
		out(1) = {'AF-AN-AR-AS-AT-EU-NA-SA-SO'};
		out(2) = {'AM-EU-OK-ON-PS-YZ'};
		out(3) = {'AN-AF-AU-JZ-NZ-PA-SA-SC-SL-SO-SW'};
		out(4) = {'AP-NZ-SA'};
		out(5) = {'AR-AF-AT-EU-IN-SO'};
		out(6) = {'AS-AF-AT-EU'};
		out(7) = {'AT-AF-AR-AS-EU'};
		out(8) = {'AU-AN-BH-BR-BS-BU-CR-FT-IN-KE-MO-NH-NI-PA-SO-SU-TI-TO-WL'};
		out(9) = {'BH-AU-BS-CL-MO-MS-PS-SU'};
		out(10) = {'BR-AU-CR-NH-PA'};
		out(11) = {'BS-AU-BH-MS-SU-TI'};
		out(12) = {'BU-AU-EU-IN-SU'};
		out(13) = {'CA-CO-NA-ND-PM-SA'};
		out(14) = {'CL-BH-NB-PA-PS-WL'};
		out(15) = {'CO-CA-GP-NA-NZ-PA-PM-RI'};
		out(16) = {'CR-AU-BR-NH'};
		out(17) = {'EA-NZ-PA'};
		out(18) = {'EU-AF-AM-AR-AS-AT-BU-IN-NA-OK-SU-YZ'};
		out(19) = {'FT-AU-NI-PA'};
		out(20) = {'GP-CO-NZ-PA'};
		out(21) = {'IN-AR-AU-BU-EU-SO'};
		out(22) = {'JF-NA-PA'};
		out(23) = {'JZ-AN-NZ-PA'};
		out(24) = {'KE-AU-PA-TO'};
		out(25) = {'MA-PA-PS'};
		out(26) = {'MN-NB-SB'};
		out(27) = {'MO-AU-BH-WL'};
		out(28) = {'MS-BH-BS-SU'};
		out(29) = {'NA-AF-CA-CO-EU-JF-OK-PA-RI-SA'};
		out(30) = {'NB-CL-MN-PA-SB-SS-WL'};
		out(31) = {'ND-CA-NZ-PM-SA'};
		out(32) = {'NH-AU-BR-CR-PA'};
		out(33) = {'NI-AU-FT-PA-TO'};
		out(34) = {'NZ-AN-AP-CO-EA-GP-JZ-ND-PA-PM-SA'};
		out(35) = {'OK-AM-EU-NA-PA-PS'};
		out(36) = {'ON-AM-PS-YZ'};
		out(37) = {'PA-AN-AU-BR-CL-CO-EA-FT-GP-JF-JZ-KE-MA-NA-NB-NH-NI-NZ-OK-PS-RI-TO-WL'};
		out(38) = {'PM-CA-CO-ND-NZ'};
		out(39) = {'PS-AM-BH-CL-MA-OK-ON-PA-SU-YZ'};
		out(40) = {'RI-CO-NA-PA'};
		out(41) = {'SA-AF-AN-AP-CA-NA-ND-NZ-SC-SW'};
		out(42) = {'SB-MN-NB-SS-WL'};
		out(43) = {'SC-AN-SA-SL-SW'};
		out(44) = {'SL-AN-SC'};
		out(45) = {'SO-AF-AN-AR-AU-IN'};
		out(46) = {'SS-NB-SB-WL'};
		out(47) = {'SU-AU-BH-BS-BU-EU-MS-PS-TI'};
		out(48) = {'SW-AN-SA-SC'};
		out(49) = {'TI-AU-BS-SU'};
		out(50) = {'TO-AU-KE-NI-PA'};
		out(51) = {'WL-AU-CL-MO-NB-PA-SB-SS'};
		out(52) = {'YZ-AM-EU-ON-PS'};
	case 'REVEL'
		out(1) = {'AM-EU-OK-PS'};
		out(2) = {'AN-AU-NU-NZ-PA-SA-SO'};
		out(3) = {'AR-AT-EU-IN-NU-SO'};
		out(4) = {'AT-AR-EU-NU'};
		out(5) = {'AU-AN-IN-PA-SO-SU'};
		out(6) = {'CA-NA-SA'};
		out(7) = {'CS-'};
		out(8) = {'EU-AM-AR-AT-IN-NA-NU-OK-SU'};
		out(9) = {'IN-AR-AU-NU-EU-SO'};
		out(10) = {'NA-CA-EU-NU-OK-PA-SA'};
		out(11) = {'NU-AN-AR-AT-EU-NA-SA-SO'};
		out(12) = {'NZ-AN-PA-SA'};
		out(13) = {'OK-AM-EU-NA-PA-PS'};
		out(14) = {'PA-AN-AU-NA-NZ-OK-PS'};
		out(15) = {'PS-AM-OK-PA-SU'};
		out(16) = {'SA-AN-CA-NA-NU-NZ'};
		out(17) = {'SO-AN-AR-AU-IN-NU'};
		out(18) = {'SR-'};
		out(19) = {'SU-AU-EU-PS'};
	case 'DEOS2K'
		out(1) = {'AN-AU-PA-SA-NU-SO'};
		out(2) = {'AU-AN-EU-PA-SO'};
		out(3) = {'EU-AU-NA-NU-PA-SO'};
		out(4) = {'NA-EU-PA-NU-SA'};
		out(5) = {'NU-AN-EU-NA-SA-SO'};
		out(6) = {'PA-AN-AU-EU-NA'};
		out(7) = {'SA-AN-NA-NU'};
		out(8) = {'SO-AN-AU-NU'};
	case 'GEODVEL'
		out(1) = {'AN-AU-PA-SA-NU-SO'};
		out(2) = {'AU-AN-EU-IN-PA-SO'};
		out(3) = {'AR-EU-IN-NU-SO'};
		out(4) = {'EU-AU-AR-IN-NA-NU-PA-SO'};
		out(5) = {'IN-AR-AU-NU-EU-SO'};
		out(6) = {'NZ-AN-PA-SA'};
		out(7) = {'NA-EU-NU-PA-SA'};
		out(8) = {'NU-AN-AR-EU-NA-SA-SO'};
		out(9) = {'PA-AN-AU-EU-NA'};
		out(10) = {'SA-AN-NA-NU-NZ'};
		out(11) = {'SO-AN-AU-AR-NU'};
end

% --- Executes on button press in push_Readme.
function push_Readme_CB(hObject, handles)
message = {'There is not much to say about the program use (specially if you know what'
	'is a plate calculator). You just click arround and watch the results. A'
	'little trial will also show that you can change the point where velocity'
	'is beeing computed, which plates (in the relative motion mode, or what'
	'plate in the absolute motion mode) are  used, or even the rotation pole'
	'itself. The black dot that shows on the map indicates the position of the'
	'current Euler pole when seen from outside the Earth. The cross inside a'
	'circle represents the position where the pole enters in Earth.'
	' '
	'However, regarding the DEOS2K and the REVEL models a couple of issues'
	'must be adressed. This program works by finding the plate currently'
	'clicked and than either computing the the velocity if it is an absolute'
	'motion, or by finding first the nearest neighbour plate (if it is a'
	'relative motion), compute the relative motion Euler pole and finally'
	'compute the velocity vector. It happens that for all of this to work, we'
	'must know the plate borders for each model and the guessing of the'
	'currently clicked plate must be correct. Regarding the first, those plate'
	'borders were not available to me, so I had to improvise using the Nuvel1A'
	'and Peter Bird''s plate models. That means that not all model poles are'
	'selectable through a mouse click. And that is also why some areas show up'
	'as white regions (particularly in the REVEL model that was based on P.'
	'Bird''s plate model). If you realy want to calculate velocities in one of'
	'these areas, manualy enter the point coordinates, select the desired plate'
	'in the popupmenu and hit "Calculate".'
	' '
	'Concearning the correct guessing of the clicked plate, you may find that'
	'the program does sometimes fail (when small plates are clicked). Well,'
	'that''s not at all my fault. Complain to the Mathworks Inc for its still'
	'somewhat buggy product.'
	' '
	'The "Relativize" checkbox appears when any of the absolute models is'
	'active. When checked, the absolute model is used to compute relative'
	'motion (poles) from the guessed (or selected) neighbour plates.'};
message_win('create',message,'figname','Help on Plate Calculator');

%--------------------------------------------------------------------------------------------------
function [plon,plat,omega] = calculate_pole(lon1,lat1,omega1,lon2,lat2,omega2)
% To calculate the relative motion, we have first to calculate relative Euler
% pole. This is because the pole list is relative to the Pacific plate. So, anyother
% plate combination that does not include the Pacific plate, has to be computed.
%
% In the following let aWb denote the rotation of the (moving) plate b relative to the (fixed) plate a 
% Given that all poles are relative to the Pacific plate (p), the closing circuit implies:
% pWa + aWb + bWp = 0
% and the desired pole (aWb) is then equal to
% aWb = -pWa - bWp = aWp + pWb
% Note that from the poles list we know pWa (= -aWp) and pWb. So:
% aWb = -pWa + pWb
%
	D2R = pi/180;
	if(lon1 == lon2 && lat1 == lat2)     % The two poles are equal
		plon = 0;   plat = 0;   omega = 0;
		return
	end

	pWa_x = omega1 * cos(lat1*D2R) * cos(lon1*D2R);
	pWa_y = omega1 * cos(lat1*D2R) * sin(lon1*D2R);
	pWa_z = omega1 * sin(lat1*D2R);

	pWb_x = omega2 * cos(lat2*D2R) * cos(lon2*D2R);
	pWb_y = omega2 * cos(lat2*D2R) * sin(lon2*D2R);
	pWb_z = omega2 * sin(lat2*D2R);

	aWb_x = -pWa_x + pWb_x;
	aWb_y = -pWa_y + pWb_y;
	aWb_z = -pWa_z + pWb_z;

	% Convert cartesian pole coordinates back to spherical coordinates
	plat = atan(aWb_z/sqrt(aWb_x*aWb_x + aWb_y*aWb_y));
	plon = atan2(aWb_y,aWb_x);
	omega = sqrt(aWb_x*aWb_x + aWb_y*aWb_y + aWb_z*aWb_z);

%--------------------------------------------------------------------------------------------------
function checkbox_Abs2Rel_CB(hObject, handles)
% Use absolute models to compute relative relative motions
if ~get(hObject,'Value')        % If we turn back to absolute motion
	handles.abs2rel = 0;
	handles.absolute_motion = 1;
	set(handles.popup_FixedPlate,'Enable','off')
	set([handles.edit_PoleLon handles.edit_PoleLat handles.edit_PoleRate],'String','')
	guidata(hObject, handles);
	return
end

D2R = pi/180;
model = getappdata(handles.figure1,'current_model');

% Fill the fixed plate popupmenus with the current model plate names
switch model
	case 'NNR'
		set(handles.popup_FixedPlate,'String',handles.Nuvel1A_NNR_name)
	case 'AKIM2000'
		set(handles.popup_FixedPlate,'String',handles.AKIM2000_name)
	case 'REVEL'
		set(handles.popup_FixedPlate,'String',handles.REVEL_name)
	case 'DEOS2K'
		set(handles.popup_FixedPlate,'String',handles.DEOS2K_name)
end
set(handles.popup_FixedPlate,'Value',1)

ind_mov = get(handles.popup_MovingPlate,'Value');
ind_fix = get(handles.popup_FixedPlate,'Value');
handles.abs2rel = 1;                % Flag that an absolute model was turned relative
handles.absolute_motion = 0;
set(handles.popup_FixedPlate,'Enable','on')
guidata(hObject, handles);

switch model
	case 'NNR'
		lat2 = handles.Nuvel1A_NNR_lat(ind_mov);     lon2 = handles.Nuvel1A_NNR_lon(ind_mov);
		omega2 = handles.Nuvel1A_NNR_omega(ind_mov); handles.abb_mov = handles.Nuvel1A_NNR_abbrev{ind_mov};
	case 'AKIM2000'
		lat2 = handles.AKIM2000_lat(ind_mov);       lon2 = handles.AKIM2000_lon(ind_mov);
		omega2 = handles.AKIM2000_omega(ind_mov);   handles.abb_mov = handles.AKIM2000_abbrev{ind_mov};        
	case 'REVEL'
		lat2 = handles.REVEL_lat(ind_mov);          lon2 = handles.REVEL_lon(ind_mov);
		omega2 = handles.REVEL_omega(ind_mov);      handles.abb_mov = handles.REVEL_abbrev{ind_mov};
	case 'DEOS2K'
		lat2 = handles.DEOS2K_lat(ind_mov);         lon2 = handles.DEOS2K_lon(ind_mov);
		omega2 = handles.DEOS2K_omega(ind_mov);     handles.abb_mov = handles.DEOS2K_abbrev{ind_mov};        
end

switch model
    case 'NNR'
        lat1 = handles.Nuvel1A_NNR_lat(ind_fix);    lon1 = handles.Nuvel1A_NNR_lon(ind_fix);
        omega1 = handles.Nuvel1A_NNR_omega(ind_fix);handles.abb_fix = handles.Nuvel1A_NNR_abbrev{ind_fix};
    case 'AKIM2000'
        lat1 = handles.AKIM2000_lat(ind_fix);       lon1 = handles.AKIM2000_lon(ind_fix);
        omega1 = handles.AKIM2000_omega(ind_fix);   handles.abb_fix = handles.AKIM2000_abbrev{ind_fix};        
    case 'REVEL'
        lat1 = handles.REVEL_lat(ind_fix);          lon1 = handles.REVEL_lon(ind_fix);
        omega1 = handles.REVEL_omega(ind_fix);      handles.abb_fix = handles.REVEL_abbrev{ind_fix};
    case 'DEOS2K'
        lat1 = handles.DEOS2K_lat(ind_fix);         lon1 = handles.DEOS2K_lon(ind_fix);
        omega1 = handles.DEOS2K_omega(ind_fix);     handles.abb_fix = handles.DEOS2K_abbrev{ind_fix};        
end
[lon,lat,omega] = calculate_pole(lon1,lat1,omega1,lon2,lat2,omega2);
lon = lon/D2R;     lat = lat/D2R;

if (omega == 0)     % This works as a test for when the same plate is selected as Fixed and Moving
	set([handles.edit_PoleLon handles.edit_PoleLat handles.edit_PoleRate],'String','')
    return
end

set(handles.edit_PoleLon,'String',num2str(lon,'%3.2f'))
set(handles.edit_PoleLat,'String',num2str(lat,'%2.2f'))
set(handles.edit_PoleRate,'String',num2str(omega,'%1.4f'))
push_Calculate_CB(hObject,handles,'nada')

% --- Creates and returns a handle to the GUI figure. 
function plate_calculator_LayoutFcn(h1)
set(h1,'PaperUnits','centimeters',...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Plate Calculator',...
'NumberTitle','off',...
'Position',[520 483 709 317],...
'RendererMode','manual',...
'doublebuffer','on',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[101 221 101 91], 'Style','frame');
uicontrol('Parent',h1, 'Position',[10 221 81 91], 'Style','frame');

uicontrol('Parent',h1, 'Position',[20 303 51 15],...
'String','Relative',...
'Style','text',...
'Tag','txt_Rel');

uicontrol('Parent',h1, 'Position',[115 303 51 15],...
'String','Absolute',...
'Style','text',...
'Tag','txt_Abs');

uicontrol('Parent',h1, 'Position',[14 285 75 15],...
'Call',{@plate_calculator_uiCB,h1,'radio_Nuvel1A_CB','Nuvel1A'},...
'String','Nuvel-1A',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_Nuvel1A');

uicontrol('Parent',h1, 'Position',[14 265 75 15],...
'Call',{@plate_calculator_uiCB,h1,'radio_Nuvel1A_CB','MORVEL'},...
'String','MORVEL',...
'Style','radiobutton',...
'Value',0,...
'Tag','radio_MORVEL');

uicontrol('Parent',h1, 'Position',[14 245 71 15],...
'Call',{@plate_calculator_uiCB,h1,'radio_PBird_CB'},...
'String','P. Bird',...
'Style','radiobutton',...
'Tag','radio_PBird');

uicontrol('Parent',h1, 'Position',[14 225 71 15],...
'Call',{@plate_calculator_uiCB,h1,'radio_GEODVEL_CB'},...
'String','GEODVEL',...
'Style','radiobutton',...
'Tag','radio_GEODVEL');

uicontrol('Parent',h1, 'Position',[108 285 100 15],...
'Call',{@plate_calculator_uiCB,h1,'radio_Nuvel1A_NNR_CB'},...
'String','Nuvel-1A NNR',...
'Style','radiobutton',...
'Tag','radio_Nuvel1A_NNR');

uicontrol('Parent',h1, 'Position',[108 265 87 15],...
'Call',{@plate_calculator_uiCB,h1,'radio_DEOS2K_CB'},...
'String','DEOS2K',...
'Style','radiobutton',...
'Tag','radio_DEOS2K');

uicontrol('Parent',h1, 'Position',[108 245 87 15],...
'Call',{@plate_calculator_uiCB,h1,'radio_REVEL_CB'},...
'String','REVEL',...
'Style','radiobutton',...
'Tag','radio_REVEL');

uicontrol('Parent',h1, 'Position',[10 174 121 22],...
'BackgroundColor',[1 1 1],...
'Call',{@plate_calculator_uiCB,h1,'popup_PickPlate_CB', 'fixed'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_FixedPlate');

uicontrol('Parent',h1, 'Position',[150 173 121 22],...
'BackgroundColor',[1 1 1],...
'Call',{@plate_calculator_uiCB,h1,'popup_PickPlate_CB', 'mobile'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_MovingPlate');

uicontrol('Parent',h1, 'Position',[10 196 52 15],...
'String','Fixed Plate',...
'Style','text');

uicontrol('Parent',h1, 'Position',[150 196 71 15],...
'HorizontalAlignment','left',...
'String','Moving Plate',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 122 71 21],...
'BackgroundColor',[1 1 1],...
'Call',{@plate_calculator_uiCB,h1,'edit_PtLon_CB'},...
'Style','edit',...
'Tag','edit_PtLon');

axes('Parent',h1, 'Units','pixels', 'Position',[300 46 401 260],...
'Color',get(0,'defaultaxesColor'),...
'xlim',[-180 180],...
'ylim',[-90 90],...
'Tag','axes1');

uicontrol('Parent',h1, 'Position',[110 122 71 21],...
'BackgroundColor',[1 1 1],...
'Call',{@plate_calculator_uiCB,h1,'edit_PtLat_CB'},...
'Style','edit',...
'Tag','edit_PtLat');

uicontrol('Parent',h1, 'Position',[205 120 66 21],...
'Call',{@plate_calculator_uiCB,h1,'push_Calculate_CB'},...
'FontSize',9,...
'FontWeight','bold',...
'String','Calculate',...
'Tag','push_Calculate');

uicontrol('Parent',h1, 'Position',[10 144 72 15],...
'String','Lon (-180:180)',...
'Style','text');

uicontrol('Parent',h1, 'Position',[110 143 72 15],...
'String','Lat (-90:90)',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 73 71 21],...
'BackgroundColor',[1 1 1],...
'Style','edit',...
'Tag','edit_PoleLon');

uicontrol('Parent',h1, 'Position',[10 95 72 15],...
'String','Pole Longitude',...
'Style','text');

uicontrol('Parent',h1, 'Position',[100 73 71 21],...
'BackgroundColor',[1 1 1],...
'Style','edit',...
'Tag','edit_PoleLat');

uicontrol('Parent',h1, 'Position',[100 95 72 15],...
'String','Pole Latitude',...
'Style','text');

uicontrol('Parent',h1, 'Position',[200 73 71 21],...
'BackgroundColor',[1 1 1],...
'Style','edit',...
'Tag','edit_PoleRate');

uicontrol('Parent',h1, 'Position',[200 95 72 15],...
'String','Rate (deg/Ma)',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 5 241 51], 'Style','frame');

uicontrol('Parent',h1, 'Position',[20 34 160 17],...
'FontSize',10,...
'HorizontalAlignment','left',...
'String','Speed  =',...
'Style','text',...
'Tag','text_Speed');

uicontrol('Parent',h1, 'Position',[20 11 220 17],...
'FontSize',10,...
'HorizontalAlignment','left',...
'String','Azimuth =',...
'Style','text',...
'Tag','text_Azim');

uicontrol('Parent',h1, 'Position',[210 285 66 21],...
'Call',{@plate_calculator_uiCB,h1,'push_Readme_CB'},...
'FontSize',9,...
'FontWeight','demi',...
'ForegroundColor',[0 0 1],...
'String','Readme',...
'Tag','push_Readme');

uicontrol('Parent',h1, 'Position',[210 232 70 15],...
'Call',{@plate_calculator_uiCB,h1,'checkbox_Abs2Rel_CB'},...
'String','Relativize',...
'Style','checkbox',...
'Tooltip','Compute relative motion from asolute model',...
'Tag','checkbox_Abs2Rel');

function plate_calculator_uiCB(hObject, eventdata, h1, callback_name, opt)
% This function is executed by the callback and than the handles is allways updated.
	if (nargin == 4)
		feval(callback_name,hObject,guidata(h1));
	else
		feval(callback_name,hObject,guidata(h1),opt);
	end

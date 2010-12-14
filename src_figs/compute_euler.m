function varargout = compute_euler(varargin)
% Helper window to calculate Euler poles

%	Copyright (c) 2004-2010 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

	if ( numel(varargin) < 1 || (numel(varargin) >= 2 && numel(varargin) <= 4) )
        errordlg('COMPUTE EULER: wrong number of input arguments.','Error'),	return
	elseif (ishandle(varargin{1}))
		isGUI = true;
	else
		isGUI = false;
	end

	if (~isGUI)			% Command line running. Doesn't create a figure.
		[handles, msg] = parse_noGUI(varargin{:});
		if (~isempty(msg) && msg(1) == 'E'),	error(msg);
		elseif (~isempty(msg)),					disp(msg),		return
		end

		handles.isoca1 = le_fiche(varargin{1});
		handles.isoca2 = le_fiche(varargin{2});
		handles.pLon_ini = varargin{3};
		handles.pLat_ini = varargin{4};
		handles.pAng_ini = varargin{5};
		[pLon, pLat, pAng, resid] = calca_pEuler(handles, true, false);
		if (nargout),	varargout{1} = [pLon, pLat, pAng, resid];	end
		return
	end

	hObject = figure('Tag','figure1','Visible','off');
	compute_euler_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	% Initialize those
	handles.hCallingFig = [];
	handles.LonRange = 30;
	handles.LatRange = 30;
	handles.AngRange = 4;
	handles.nInt_lon = 21;
	handles.nInt_lat = 21;
	handles.nInt_ang = 21;
	handles.isoca1 = [];
	handles.isoca2 = [];
	handles.pLon_ini = [];
	handles.pLat_ini = [];
	handles.pAng_ini = [];
	handles.do_graphic = false;
	handles.DP_tol = 0.05;
	handles.residGrdName = [];
	set(handles.slider_wait,'Max',handles.nInt_lon * handles.nInt_lat)

	handles.hCallingFig = varargin{1};        % This is the Mirone's fig handle

	handMir = guidata(handles.hCallingFig);
	if (handMir.no_file)
		errordlg('You didn''t even load a file. What are you expecting then?','ERROR')
		delete(hObject);    return
	end
	if (~handMir.geog)
		errordlg('This operation is currently possible only for geographic type data','ERROR')
		delete(hObject);    return
	end
	handles.path_continent = [handMir.home_dir filesep 'continents' filesep];
	handles.IamCompiled = handMir.IamCompiled;		% Need to know due to crazy issue of nc_funs

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, [handles.txtSP handles.txtDS handles.txtCS])
	%------------- END Pro look (3D) -----------------------------------------------------

	set(hObject,'Visible','on');
	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end

% -------------------------------------------------------------------------------------
function edit_first_file_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),		return,		end
	% Let the push_first_file_CB do all the work
	push_first_file_CB(hObject,handles,fname)

% -------------------------------------------------------------------------------------
function push_first_file_CB(hObject, handles,opt)
	if (nargin == 3),	fname = opt;
	else				opt = [];
	end

	if (isempty(opt))    % Otherwise we already know fname from the 4th input argument
		handMir = guidata(handles.hCallingFig);
		[FileName,PathName] = put_or_get_file(handMir,{'*.dat;*.DAT', 'Mag file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select file','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	end
	handles.isoca1 = le_fiche(fname);      % A maluca
	%handles.isoca1(:,2) = geog2auth(handles.isoca1(:,2));   % Convert to authalic lats
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_second_file_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),		return,		end
	% Let the push_first_file_CB do all the work
	push_second_file_CB(hObject,handles,fname)

% -------------------------------------------------------------------------------------
function push_second_file_CB(hObject, handles, opt)
	if (nargin == 3),	fname = opt;
	else				opt = [];
	end

	if (isempty(opt))    % Otherwise we already know fname from the 4th input argument
		handMir = guidata(handles.hCallingFig);
		[FileName,PathName] = put_or_get_file(handMir,{'*.dat;*.DAT', 'Mag file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select file','get');
		if isequal(FileName,0),		return,		end	   
		fname = [PathName FileName];
	end
	handles.isoca2 = le_fiche(fname);      % The fixed line
	%handles.isoca2(:,2) = geog2auth(handles.isoca2(:,2));   % Convert to authalic lats
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_pLon_ini_CB(hObject, handles)
    handles.pLon_ini = str2double(get(hObject,'String'));
    guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_pLat_ini_CB(hObject, handles)
    handles.pLat_ini = str2double(get(hObject,'String'));
    guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_pAng_ini_CB(hObject, handles)
    handles.pAng_ini = str2double(get(hObject,'String'));
    guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function push_polesList_CB(hObject, handles)
fid = fopen([handles.path_continent 'lista_polos.dat'],'rt');
c = fread(fid,'*char').';
fclose(fid);
s = strread(c,'%s','delimiter','\n');

[s,v] = choosebox('Name','One Euler list',...
                    'PromptString','List of poles:',...
                    'SelectString','Selected poles:',...
                    'ListSize',[380 300],'ListString',s);

if (v == 1)         % Finite pole
    handles.pLon_ini = s(1);
    handles.pLat_ini = s(2);
    handles.pAng_ini = s(3);
    set(handles.edit_pLon_ini, 'String', num2str(s(1)))
    set(handles.edit_pLat_ini, 'String', num2str(s(2)))
    set(handles.edit_pAng_ini, 'String', num2str(s(3)))
    guidata(hObject,handles)
elseif (v == 2)     % Stage poles
    %set(handles.edit_polesFile,'String',s)
end

% -------------------------------------------------------------------------------------
function edit_LonRange_CB(hObject, handles)
	handles.LonRange = str2double(get(hObject,'String'));
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_nInt_CB(hObject, handles)
	if (~get(handles.check_hellinger,'Val'))
		nInt = abs(sscanf(get(hObject,'Str'),'%d'));	% We want an odd number
		if (~rem(nInt,2)),		nInt = nInt + 1;	end
		tag = get(hObject,'Tag');
		if (strcmp(tag(end-2:end),'lon')),			handles.nInt_lon = nInt;
		elseif (strcmp(tag(end-2:end),'lat')),		handles.nInt_lat = nInt;
		else										handles.nInt_ang = nInt;
		end
		set(handles.slider_wait,'Max',handles.nInt_lon * handles.nInt_lat)
	else
		handles.DP_tol = str2double(get(hObject,'Str'));
	end
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_LatRange_CB(hObject, handles)
	handles.LatRange = str2double(get(hObject,'String'));
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_AngRange_CB(hObject, handles)
	handles.AngRange = str2double(get(hObject,'String'));
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function toggle_pickLines_CB(hObject, handles)
if (get(hObject,'Value'))
    % Test if we have potential target lines and their type
    h_mir_lines = findobj(handles.hCallingFig,'Type','line');     % Fish all objects of type line in Mirone figure
    if (isempty(h_mir_lines))                                       % We don't have any lines
        str = ['If you hited this button on purpose, than you deserve the following insult.',...
                'You #!|"*!%!?~^)--$&.',... 
                'THERE ARE NO LINES IN THAT FIGURE.'];
        errordlg(str,'Chico Clever');   set(hObject,'Value',0);     return;
    end
    if (length(h_mir_lines) == 1)                                    % We don't have at least two lines
        str = ['If you hited this button on purpose, than you deserve the following insult.',...
                'You -$&#!*!%!?~^)-.|"/',... 
                'THERE IS ONLY ONE LINE IN THAT FIGURE.'];
        errordlg(str,'Chico Clever');   set(hObject,'Value',0);     return;
    end
    
    % The above test is not enough. For exemple, coastlines are not eligible neither,
    % but is very cumbersome to test all the possibilities of pure non-eligible lines.
    set(handles.hCallingFig,'pointer','crosshair')
    h_line = get_polygon(handles.hCallingFig);		% Get first line handle
    if (~isempty(h_line))
        x = get(h_line,'XData');		y = get(h_line,'YData');
        %y = geog2auth(y);							% Convert to authalic latitudes
        handles.isoca1 = [x(:) y(:)];
        set(handles.edit_first_file,'String','Got left line','FontAngle','italic')
    else
        handles.isoca1 = [];
        set(handles.edit_first_file,'String','')
    end
    h_line = get_polygon(handles.hCallingFig);		% Get second line handle
    if (~isempty(h_line))
        x = get(h_line,'XData');		y = get(h_line,'YData');
        %y = geog2auth(y);							% Convert to authalic latitudes
        handles.isoca2 = [x(:) y(:)];
        set(handles.edit_second_file,'String','Got right line','FontAngle','italic')
    else
        handles.isoca2 = [];
        set(handles.edit_second_file,'String','')
    end
    set(handles.hCallingFig,'pointer','arrow')
    if (isempty(handles.isoca1) || isempty(handles.isoca2))
        set(hObject,'Value',0)
        handles.do_graphic = false;
    else
        handles.do_graphic = true;
    end
    set(hObject,'Value',0)
    figure(handles.figure1)         % Bring this figure to front again
else        % What should I do?
    %handles.do_graphic = 0;
end
guidata(hObject, handles);

% -----------------------------------------------------------------------------
function check_hellinger_CB(hObject, handles)
if (get(hObject,'Value'))
	set([handles.edit_LonRange handles.edit_LatRange handles.edit_AngRange],'Enable','off')
	set([handles.edit_nInt_lat handles.edit_nInt_ang],'Vis','off')
	set(handles.edit_nInt_lon,'String', handles.DP_tol)
	set(handles.textNint,'String','DP tolerance')
	set(handles.edit_nInt_lon,'Tooltip', sprintf(['Tolerance used to break up the isochron into\n' ...
			'linear chunks (the Heillinger segments).\n' ...
			'The units of the tolerance are degrees\n', ...
			'of arc on the surface of a sphere']))
else
	set([handles.edit_LonRange handles.edit_LatRange handles.edit_AngRange],'Enable','on')
	set([handles.edit_nInt_lat handles.edit_nInt_ang],'Vis','on')
	set(handles.textNint,'String','N Intervals')
	set(handles.edit_nInt_lon,'Str', handles.nInt_lon,'Tooltip','The range parameters are divided into this number of intervals steps')
end

% -----------------------------------------------------------------------------
function edit_err_file_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname)
		set([handles.radio_netcd handles.radio_VTK],'Val',0)
		handles.residGrdName = [];
		guidata(handles.figure1, handles);
		return
	end
	% Let the push_err_file_CB do all the work
	push_err_file_CB(hObject, handles, fname)

% -----------------------------------------------------------------------------
function push_err_file_CB(hObject, handles, opt)
	if (nargin == 3),	fname = opt;
	else				opt = [];
	end

	if (isempty(opt))    % Otherwise we already know fname from the 3th input argument
		handMir = guidata(handles.hCallingFig);
		[FileName,PathName] = put_or_get_file(handMir,{'*.nc;*.grd;*.vtk', 'Error file (*.nc,*.grd,*.vtk)';'*.*', 'All Files (*.*)'},'Select file','put');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	end
	ind = strfind(fname,'.');
	if (~isempty(ind))			% Make a guess based on extension
		ext = fname(ind(end)+1:end);
		if (strcmpi(ext,'nc') || strcmpi(ext,'grd'))
			set(handles.radio_netcdf,'Val',1),		set(handles.radio_VTK,'Val',0)
		elseif (strcmpi(ext,'vtk'))
			set(handles.radio_VTK,'Val',1),			set(handles.radio_netcdf,'Val',0)
		end
	end
	set(handles.edit_err_file,'Str',fname);
	handles.residGrdName = fname;
	guidata(handles.figure1, handles);

% -----------------------------------------------------------------------------
function radio_netcdf_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_VTK, 'Val',0)

% -----------------------------------------------------------------------------
function radio_VTK_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_netcdf, 'Val',0)

% -------------------------------------------------------------------------------
function push_stop_CB(hObject, handles)
    set(handles.slider_wait,'Value',0)  % We have to do this first
    set(handles.slider_wait,'Max',1)    % This will signal the fit_pEuler function to stop

% -------------------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
% OK. See if we have all the information needed to compute the Euler pole

	set(handles.edit_BFresidue,'String','');		set(handles.edit_InitialResidue,'String','')
	set(handles.edit_pLon_fim,'String','');			set(handles.edit_pLat_fim,'String','')
	set(handles.edit_pAng_fim,'String','')

	if (isempty(handles.isoca1) || isempty(handles.isoca2))
		errordlg('Compute Euler pole with what? It would help if you provide me TWO lines.','Chico Clever')
		return
	end
	if (isempty(handles.pLon_ini) || isempty(handles.pLat_ini) || isempty(handles.pAng_ini))
		errordlg(['I need a first guess of the Euler pole you are seeking for.' ...
			'Pay attention to the "Starting Pole Section"'],'Error')
		return
	end

	if (~get(handles.check_hellinger,'Val'))		% Our method
		do_weighted = true;
		calca_pEuler(handles, do_weighted, true);
	else											% Try with Hellinger's (pfiu)
		[pLon,pLat,pAng] = hellinger(handles.pLon_ini,handles.pLat_ini,handles.pAng_ini, handles.isoca1, handles.isoca2, handles.DP_tol);
		set(handles.edit_pLon_fim,'String',pLon);			set(handles.edit_pLat_fim,'String',pLat)
		set(handles.edit_pAng_fim,'String',pAng)
		if (handles.do_graphic)     % Create a empty line handle
			[rlon,rlat] = rot_euler(handles.isoca1(:,1),handles.isoca1(:,2),pLon, pLat, pAng,-1);
			h_line = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',rlon,'YData',rlat, ...
			'LineStyle','-.','LineWidth',2,'Tag','Fitted Line','Userdata',1);
			draw_funs(h_line,'isochron',{'Fitted Line'})
		end
	end

% -------------------------------------------------------------------------------
function numeric_data = le_fiche(fname)
	[bin,n_column,multi_seg,n_headers] = guess_file(fname);
	% If error in reading file
	if (isempty(bin))
		errordlg(['Error reading file ' fname],'Error');
		numeric_data = [];
		return
	end
	if (isempty(n_headers)),     n_headers = NaN;    end
	if (multi_seg)
		numeric_data = text_read(fname,NaN,n_headers,'>');
	else
		numeric_data = text_read(fname,NaN,n_headers);
	end

% -------------------------------------------------------------------------------
function [polLon, polLat, polAng, area_f] = calca_pEuler(handles, do_weighted, isGUI)

	D2R = pi / 180;
	if (handles.do_graphic)     % Create a empty line handle
		h_line = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',[],'YData',[], ...
			'LineStyle','-.','LineWidth',2,'Tag','Fitted Line','Userdata',1);
	else
		h_line = [];
	end

	% Compute distances between vertices of the moving isoc
	xd = diff( (handles.isoca1(:,1) .* cos(handles.isoca1(:,2) * D2R) ) * D2R * 6371 );
	yd = diff( handles.isoca1(:,2) * D2R * 6371 );
	lenRot1 = sqrt(xd.*xd + yd.*yd);

	[rlon,rlat] = rot_euler(handles.isoca1(:,1),handles.isoca1(:,2),handles.pLon_ini,handles.pLat_ini,handles.pAng_ini,-1);
	[dist1, segLen] = distmin(handles.isoca2(:,1)*D2R, handles.isoca2(:,2)*D2R, rlon*D2R, rlat*D2R, lenRot1, 1e20, do_weighted);
	sum1 = weightedSum(dist1, segLen, do_weighted);

	xd = diff( (handles.isoca2(:,1) .* cos(handles.isoca2(:,2) * D2R) ) * D2R * 6371 );
	yd = diff( handles.isoca2(:,2) * D2R * 6371 );
	lenRot2 = sqrt(xd.*xd + yd.*yd);

	[dist2, segLen] = distmin(rlon*D2R, rlat*D2R, handles.isoca2(:,1)*D2R, handles.isoca2(:,2)*D2R, lenRot2, 1e20, do_weighted);
	sum2 = weightedSum(dist2, segLen, do_weighted);
	area0 = (sum1 + sum2) / 2;

	if (handles.do_graphic)
		set(handles.edit_InitialResidue,'String',sprintf('%.3f', area0));    pause(0.01)
		set(h_line,'XData',rlon,'YData',rlat)
	end

	% Now comes the semi-brute force aproach to compute the pole
	dLon = handles.LonRange / 2;
	dLat = handles.LatRange / 2;
	dAng = handles.AngRange / 2;
	p_lon = (handles.pLon_ini + linspace(-dLon,dLon, handles.nInt_lon)) * D2R;
	p_lat = (handles.pLat_ini + linspace(-dLat,dLat, handles.nInt_lat)) * D2R;
	p_omeg = (handles.pAng_ini + linspace(-dAng,dAng,handles.nInt_ang)) * D2R;

	% Sanitize p_lat so that it does not go out of N/S poles
	ind = ( (p_lat > pi/2) | (p_lat < -pi/2) );
	if (any(ind))
		p_lat(ind) = [];
		handles.nInt_lat = numel(p_lat);
		if (handles.do_graphic)
			set(handles.slider_wait,'Max',handles.nInt_lon * handles.nInt_lat)
		end
	end

	[polLon, polLat, polAng, area_f, resid] = ...
		fit_pEuler(handles, p_lon, p_lat, p_omeg, area0, h_line, lenRot1, lenRot2, do_weighted, isGUI);

	if (handles.do_graphic),	draw_funs(h_line,'isochron',{'Fitted Line'}),	end

	if (~isempty(resid) && isGUI)
		if (get(handles.radio_netcdf,'Val'))
			write_netcdf(handles, p_lon/D2R, p_lat/D2R, p_omeg/D2R, resid)
		else
			write_vtk(handles, p_lon/D2R, p_lat/D2R, p_omeg/D2R, resid)
		end
	elseif (~isempty(resid))	% Command line runn.
		write_netcdf(handles, p_lon/D2R, p_lat/D2R, p_omeg/D2R, resid)
	end

% -----------------------------------------------------------------------------------------
function write_netcdf(handles, lon, lat, ang, resid)
% Write the 3D matrix in netCDF

	nz = numel(ang);

	% No problem in changig handles here as it will live only locally to this function
	handles.head = [lon(1) lon(end) lat(1) lat(end) 0 0 0 diff(lon(1:2)) diff(lat(1:2))];
	handles.geog = 1;		handles.was_int16 = false;
	Z = resid(:,:,1);
	handles.head(5:6) = [min(Z(:)) max(Z(:))];			Z = single(Z);
	if (~handles.IamCompiled)
		nc_io(handles.residGrdName, sprintf('w-%s/angle',ang(1)), handles, reshape(Z,[1 size(Z)]))
	else
		handles.levelVec = ang;
		nc_io(handles.residGrdName,sprintf('w%d/angle',nz), handles, reshape(Z,[1 size(Z)]))
	end
	for (k = 2:nz)
		Z = resid(:,:,k);
		handles.head(5:6) = [min(Z(:)) max(Z(:))];		Z = single(Z);
		if (~handles.IamCompiled)
			nc_io(handles.residGrdName, sprintf('w%d\\%s', k-1, ang(k)), handles, Z)
		else
			nc_io(handles.residGrdName, sprintf('w%d', k-1), handles, Z)
		end
	end

% -----------------------------------------------------------------------------------------
function write_vtk(handles, lon, lat, ang, resid)
% Write in the VTK format

	nx = numel(lon);	ny = numel(lat);	nz = numel(ang);
	fid = fopen(handles.residGrdName, 'wb','b');
	fprintf(fid, '# vtk DataFile Version 2.0\n');
	fprintf(fid, 'converted from A B\n');
	fprintf(fid, 'BINARY\n');
	fprintf(fid, 'DATASET RECTILINEAR_GRID\n');
	fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
	fprintf(fid, 'X_COORDINATES %d float\n', nx);
	fwrite(fid, lon, 'real*4');
	fprintf(fid, 'Y_COORDINATES %d float\n', ny);
	fwrite(fid, lat, 'real*4');
	fprintf(fid, 'Z_COORDINATES %d float\n', nz);
	fwrite(fid, ang, 'real*4');
	fprintf(fid, 'POINT_DATA %d\n', nx * ny * nz);
	fprintf(fid, 'SCALARS dono float 1\n');
	fprintf(fid, 'LOOKUP_TABLE default\n');

	for (k = 1:nz)
		Z = single(resid(:,:,k))';
		fwrite(fid, Z(:), 'real*4');
	end
	fclose(fid);

% -------------------------------------------------------------------------------
function [lon_bf, lat_bf, omega_bf, area_f, resid] = ...
		fit_pEuler(handles, p_lon, p_lat, p_omeg, area0, h_line, lenRot1,  lenRot2, do_weighted, isGUI)
% This is the function that does the real work of computing the best Euler pole
% that fits the lines "isoca1" "isoca2"
	lon_bf = [];		lat_bf = [];	omega_bf = [];
	D2R = pi / 180;
	nLon = numel(p_lon);	nLat = numel(p_lat);	nAng = numel(p_omeg);
	isoca1 = handles.isoca1 * D2R;		% A maluca
	isoca2 = handles.isoca2 * D2R;		% A fixa

	if (~isempty(handles.residGrdName))	% Store all residues in a 3D array
		save_resid = true;		resid = zeros(nLat,nLon,nAng) * NaN;	testResidue = 1e20;
	else
		save_resid = false;		resid = [];		testResidue = area0;
	end

	ecc = 0.0818191908426215;	% WGS84

	i_m = 1;	j_m = 1;	k_m = 1;	ij = 0;
	for (i = 1:nLon)			% Loop over n lon intervals
		if (isGUI && get(handles.slider_wait,'Max') == 1)			% The STOP button was pushed. So, stop
			set(handles.slider_wait,'Max',handles.nInt_lon * handles.nInt_lat)     % Reset it for the next run
			area_f = -1;		resid = [];
			return
		end
		for (j = 1:nLat)		% Loop over n lat intervals
			ij = ij + 1;
			if (isGUI && get(handles.slider_wait,'Max') == 1)		% The STOP button was pushed. So, stop
				set(handles.slider_wait,'Max',handles.nInt_lon * handles.nInt_lat)     % Reset it for the next run
				area_f = -1;	resid = [];
				return
			end
			if (isGUI),		set(handles.slider_wait,'Value',ij);     pause(0.01);   end	% Otherwise slider may stall
			for (k = 1:nAng)	% Loop over n omega intervals
				%[rlon,rlat] = rot_euler(isoca1(:,1),isoca1(:,2),p_lon(i),p_lat(j),p_omeg(k),'radians',-1);
                
                lat = atan2( (1-ecc^2)*sin(isoca1(:,2)), cos(isoca1(:,2)) );
				p_sin_lat = sin(p_lat(j));				p_cos_lat = cos(p_lat(j));
				s_lat = sin(lat);						c_lat = cos(lat);
				s_lon = sin(isoca1(:,1) - p_lon(i));	c_lon = cos(isoca1(:,1) - p_lon(i));
				cc = c_lat .* c_lon;
			
				tlon = atan2(c_lat .* s_lon, p_sin_lat * cc - p_cos_lat * s_lat);
				s_lat = p_sin_lat * s_lat + p_cos_lat * cc;
				c_lat = sqrt(1 - s_lat .* s_lat);
			
				s_lon = sin(tlon + p_omeg(k));		c_lon = cos(tlon + p_omeg(k));
				cc = c_lat .* c_lon;
			
				rlat = asin(p_sin_lat * s_lat - p_cos_lat * cc);
				rlon = p_lon(i) + atan2(c_lat .* s_lon, p_sin_lat * cc + p_cos_lat * s_lat);
			
				rlat = atan2( sin(rlat), (1-ecc^2)*cos(rlat) );
				ind = (rlon > pi);					rlon(ind) = rlon(ind) - 2*pi;

				
				[dist1, segLen] = distmin(isoca2(:,1), isoca2(:,2), rlon,rlat, lenRot1, testResidue, do_weighted);
				sum1 = weightedSum(dist1, segLen, do_weighted);
				[dist2, segLen] = distmin(rlon,rlat, isoca2(:,1), isoca2(:,2), lenRot2, testResidue, do_weighted);
				sum2 = weightedSum(dist2, segLen, do_weighted);

				area = (sum1 + sum2) / 2;
				if (area < area0)
					area0 = area;
					i_m = i;	j_m = j;	k_m = k;
					if (handles.do_graphic)
						set(h_line,'XData',rlon/D2R,'YData',rlat/D2R)
						pause(0.01)
					end
					if (isGUI)
						set(handles.edit_pLon_fim,'String',sprintf('%.3f', p_lon(i) / D2R))
						set(handles.edit_pLat_fim,'String',sprintf('%.3f', p_lat(j) / D2R))
						set(handles.edit_pAng_fim,'String',sprintf('%.3f', p_omeg(k) / D2R))
						set(handles.edit_BFresidue,'String',sprintf('%.3f', area))
					end
				end

				if (save_resid),	resid(j, i, k) = area;		end

			end
		end
		if (~isGUI),	fprintf('%d\t Lons rots out of %d\r', i, nLon),	end
	end

	if (isGUI),		set(handles.slider_wait,'Value',0),		end      % Reset it for the next run
	lon_bf = p_lon(i_m) / D2R;
	lat_bf = p_lat(j_m) / D2R;
	omega_bf = p_omeg(k_m) / D2R;
	area_f = area0;

% -------------------------------------------------------------------------------
function [dist, segLen] = distmin(lon, lat, r_lon, r_lat, lengthsRot, lastResidue, do_weighted)
% Compute the shortest distance between each point in (lon,lat) and the polyline (r_lon,r_lat)
% SEGLEN holds the segment length of polyline (r_lon,r_lat) corresponding to elements of DIST
% At each quater of the total number of points in (lon.lat) we check if the current residue is
% already larger than last better residue. If yes, we stop since for sure this rotation is worst.
% The above test is short-circuited (not carried out) when saving the residues to a grid.
% Angles are already in radians.

	r_lon = r_lon .* cos(r_lat) * 6371;		% VERY strange. This should improve the fit
	lon = lon .* cos(lat) * 6371;			% But, on the contrary, it degrades it ??
	r_lon = r_lon(:)';		r_lat = r_lat(:)' * 6371;      % Make sure they are row vectors
	lat = lat * 6371;
	%eps1 = 1e-1;			eps2 = eps1 * eps1;
	n_pt = numel(lon);		n_pt_rot = numel(r_lon);
	dist = zeros(n_pt,1);
	segLen = ones(n_pt,1);
	outliners = false(n_pt,1);      % To remove vertex where lines do not intersect

	chunks = round([(n_pt * 0.25) (n_pt / 2) (n_pt * 0.75) (n_pt * 0.9) (n_pt+1)]);	% Checkpoints where current res is checked against min res
	for (k = 1:n_pt)				% Loop over vertices of fixed isoc
		Dsts = sqrt((lon(k)-r_lon).^2 + (lat(k)-r_lat).^2);		% distances from pt k in (lon,lat) to all pts in (r_lon,r_lat)
		[D,ind] = min(Dsts);		% Compute the min dist from pt k to vertex of line and its location
		if (ind == 1 || ind == n_pt_rot || ind >= n_pt && n_pt > 4)		% This point is outside the lines intersection. Flag it to die.
			outliners(k) = true;						% Actually we waste all points that are closest to the rotated line end points
			continue
		end

 		P  = [lon(k) lat(k)];
		%d = abs(det([Q2-Q1,P-Q1]))/norm(Q2-Q1); % for col. vectors
		Q1 = [r_lon(ind-1) r_lat(ind-1)];		Q2 = [r_lon(ind) r_lat(ind)];
		D1 = abs(det([Q2-Q1; P-Q1])) / norm(Q2-Q1); % for row vectors.
		Q1 = [r_lon(ind) r_lat(ind)];			Q2 = [r_lon(ind+1) r_lat(ind+1)];
		D2 = abs(det([Q2-Q1; P-Q1])) / norm(Q2-Q1);

		[dist(k),i] = min([D1 D2]);				% Store shortest dist from point k to polyline (r_lon,r_lat)
		if (i == 1),		segLen(k) = lengthsRot(ind-1);	% Store the segment length that is closest to pt k
		else				segLen(k) = lengthsRot(ind);
		end

		if ( lastResidue < 1e20 && (k == chunks(1)) )	% At 25, 50, and 75% of the points check if residue is already larger than minimum
			res = weightedSum(dist, segLen, do_weighted);
			if (res > lastResidue)
				chunks = chunks(2:end);
				break
			end
		end
		
% 		dist0 = 1e20;
% 		d1 = 0;		d2 = 0;		xm = 0;		ym = 0;
% 		for (j = ind-1:ind)
% 			x1 = r_lon(j);      x2 = r_lon(j+1);
% 			y1 = r_lat(j);      y2 = r_lat(j+1);
% 			dd = 1e20;
% 			while (dd > eps2)
% 				d1 = sqrt((lon(k)-x1)^2 + (lat(k)-y1)^2);
% 				d2 = sqrt((lon(k)-x2)^2 + (lat(k)-y2)^2);
% 				dd = (x1-x2)^2 + (y1-y2)^2;
% 				xm = (x1+x2) * 0.5;
% 				ym = (y1+y2) * 0.5;
% 				if (d1 < d2)
% 					x2 = xm;    y2 = ym;
% 				else
% 					x1 = xm;    y1 = ym;
% 				end
% 			end
% 			dist1 = (d1+d2) * 0.5;
% 			if(dist1 < dist0)		% When true, the closest point on the other line was found
%      			segLen(k) = lengthsRot(j);
% 				dist0 = dist1;
% 			end 
% 		end
% 		% 'dist' has the least distance of each vertex of line (lon,lat) to rotated line (r_lon,r_lat)
% 		dist(k) = dist0;
	end

	dist(outliners) = [];		% Delete points not belonging to the lines intersection zone
	segLen(outliners) = [];		% and the same for the weights
	if (isempty(dist)),     dist = 0;   segLen = 1;		end     % Don't let it go empty

% --------------------------------------------------------------------------------
function soma = weightedSum(dists, segLen, do_weighted)
% Convert the segment lengths along an isochron into a series of weights
% Lengths < 50 km weight 1. In the [50 80] interval weight 0.25. Longer weight zero
	if (do_weighted)
		weights = ones(numel(segLen), 1);
		ind = (segLen > 50 & segLen < 80);
		weights(ind) = 0.25;
		weights(segLen >= 80) = 0;
		soma = sum( dists .* weights ) / sum(weights);		% weighted sum
	else
		soma = sum( dists ) / numel(dists);
	end


% --------------------------------------------------------------------------------
function lat = geog2auth(lat0)
% Convert geodetic (geographic) latitudes to authalic latitudes. 
% Degrees are assumed in input
	lat0 = lat0 * pi / 180;
	flat = 1/298.257222101;         % GRS80 flatness
	ecc = sqrt(2*flat - flat.^2);   % exccentricity
	f1 = ecc^2 /3 + 31*ecc^4 / 180 + 59*ecc^6 / 560;
	f2 = 17*ecc^4 / 360 + 61*ecc^6 / 1260;
	f3 = 383*ecc^6 / 45360;

	lat = lat0 - f1*sin(2*lat0) + f2*sin(4*lat0) - f3*sin(6*lat0);
	lat = lat * 180 / pi;  %  Convert back to degrees

% --------------------------------------------------------------------------------
function lat = auth2geog(lat0)
% Convert authalic  latitudes to geodetic (geographic) latitudes. 
% Degrees are assumed in input
	lat0 = lat0 * pi / 180;
	flat = 1/298.257222101;         % GRS80 flatness
	ecc = sqrt(2*flat - flat.^2);   % exccentricity
	f1 = ecc^2 / 3 + 31*ecc^4 / 180 + 517*ecc^6 / 5040;
	f2 = 23*ecc^4 / 360 + 251*ecc^6 / 3780;
	f3 = 761*ecc^6 / 45360;

	lat = lat0 + f1*sin(2*lat0) + f2*sin(4*lat0) + f3*sin(6*lat0);
	lat = lat * 180 / pi;  %  Convert back to degrees

% ---------------------------------------------------------------------------------
function [handles, msg] = parse_noGUI(varargin)
%	
	handles.LonRange = 30;
	handles.LatRange = 30;
	handles.AngRange = 4;
	handles.nInt_lon = 21;
	handles.nInt_lat = 21;
	handles.nInt_ang = 21;
	handles.do_graphic = false;
	handles.residGrdName = [];
	handles.IamCompiled = false;	% Até ver
	msg = [];	verbose = false;

	if (numel(varargin) < 6),	return,		end		% No option provided. Just use the defaults

	for (k = 6:numel(varargin))
		switch (varargin{k}(2))
			case 'H'
				msg = '[-DLonRange/LatRange/AngRange] [-ILonInt/LatInt/AngInt or -Iint] [-Eresid_fname] [-V]';
				return
			case 'D'
				ind = strfind(varargin{k},'/');
				if (numel(ind) ~= 2)
					msg = 'Error usin option -D [-DLonRange/LatRange/AngRange]';		return
				end
				handles.LonRange = sscanf(varargin{k}(3:ind(1)-1), '%f');
				handles.LatRange = sscanf(varargin{k}(ind(1)+1:ind(2)-1), '%f');
				handles.AngRange = sscanf(varargin{k}(ind(2)+1:end), '%f');

			case 'I'
				ind = strfind(varargin{k},'/');
				if (~isempty(ind) && numel(ind) ~= 2)
					msg = 'Error usin option -I: -ILonInt/LatInt/AngInt or -Iint';		return
				end
				if (isempty(ind))
					handles.nInt_lon = sscanf(varargin{k}(3:ind(1)-1), '%f');
					handles.nInt_lat = handles.nInt_lon;	handles.nInt_ang = handles.nInt_lon;
					if (~rem(handles.nInt_lon,2)),			handles.nInt_lon = handles.nInt_lon + 1;	end
					if (~rem(handles.nInt_lat,2)),			handles.nInt_lat = handles.nInt_lat + 1;	end
					if (~rem(handles.nInt_ang,2)),			handles.nInt_ang = handles.nInt_ang + 1;	end
				else
					handles.nInt_lon = sscanf(varargin{k}(3:ind(1)-1), '%f');
					handles.nInt_lat = sscanf(varargin{k}(ind(1)+1:ind(2)-1), '%f');
					handles.nInt_ang = sscanf(varargin{k}(ind(2)+1:end), '%f');					
				end
				
			case 'E'
				handles.residGrdName = varargin{k}(3:end);

			case 'V'
				verbose = true;
		end
	end

	if (verbose)
		fprintf('LonRange = %d\tLatRange = %d\tAngRange = %d\n', handles.LonRange,handles.LatRange,handles.AngRange);
		fprintf('LonInt = %d\tLatInt = %d\tAngInt = %d\n', handles.nInt_lon,handles.nInt_lat,handles.nInt_ang);
		if (~isempty(handles.residGrdName))
			fprintf('Save residues to file: %s\n', handles.residGrdName);
		end
	end
     
% --- Creates and returns a handle to the GUI figure. 
function compute_euler_LayoutFcn(h1)
set(h1, 'Pos',[520 379 520 421],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Compute Euler pole',...
'NumberTitle','off',...
'Resize','off',...
'Tag','figure1')

uicontrol('Parent',h1,'Pos',[10 50 501 157],'String',{''},'Style','frame')
uicontrol('Parent',h1,'Pos',[10 230 501 64],'String',{''},'Style','frame');
uicontrol('Parent',h1,'Pos',[10 313 501 101],'String',{''},'Style','frame');

uicontrol('Parent',h1, 'Pos',[20 363 211 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_first_file_CB'},...
'Style','edit','Tag','edit_first_file');

uicontrol('Parent',h1, 'Pos',[231 363 21 21],...
'Call',{@compute_euler_uiCB,h1,'push_first_file_CB'},...
'FontSize',10,...
'FontWeight','bold',...
'String','...','Tag','push_first_file');

uicontrol('Parent',h1, 'Pos',[267 363 211 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_second_file_CB'},...
'Style','edit','Tag','edit_second_file');

uicontrol('Parent',h1, 'Pos',[478 363 21 21],...
'Call',{@compute_euler_uiCB,h1,'push_second_file_CB'},...
'FontSize',10,...
'FontWeight','bold',...
'String','...','Tag','push_second_file');

uicontrol('Parent',h1, 'Pos',[186 326 141 21],...
'Call',{@compute_euler_uiCB,h1,'toggle_pickLines_CB'},...
'String','Pick lines from Figure',...
'Style','togglebutton',...
'Tooltip','Allows you to mouse select the two lines from a Mirone figure',...
'Tag','toggle_pickLines');

uicontrol('Parent',h1, 'Pos',[20 249 61 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_pLon_ini_CB'},...
'Style','edit',...
'Tooltip','Start Euler pole Longitude',...
'Tag','edit_pLon_ini');

uicontrol('Parent',h1, 'Pos',[110 249 61 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_pLat_ini_CB'},...
'Style','edit',...
'Tooltip','Start Euler pole Latitude',...
'Tag','edit_pLat_ini');

uicontrol('Parent',h1, 'Pos',[200 249 61 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_pAng_ini_CB'},...
'Style','edit',...
'Tooltip','Start Euler pole Angle',...
'Tag','edit_pAng_ini');

uicontrol('Parent',h1, 'Pos',[289 248 121 21],...
'Call',{@compute_euler_uiCB,h1,'push_polesList_CB'},...
'String','Poles selector',...
'Tooltip','Select a pole from the default list',...
'Tag','push_polesList');

uicontrol('Parent',h1, 'Pos',[20 168 84 16],'HorizontalAlignment','left','Str','Longitude Range','Style','text');

uicontrol('Parent',h1, 'Pos',[104 165 41 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_LonRange_CB'},...
'String','30',...
'Style','edit',...
'Tooltip','The pole will be searched arround it''s starting longitude +/- half this range',...
'Tag','edit_LonRange');

uicontrol('Parent',h1, 'Pos',[158 186 60 15],'String','N Intervals','Style','text','Tag','textNint');
uicontrol('Parent',h1, 'Pos',[20 139 75 18],'HorizontalAlignment','left','Str','Latitude Range','Style','text');

uicontrol('Parent',h1, 'Pos',[104 138 41 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_LatRange_CB'},...
'String','30',...
'Style','edit',...
'Tooltip','The pole will be searched arround it''s starting latitude +/- half this range',...
'Tag','edit_LatRange');

uicontrol('Parent',h1, 'Pos',[104 110 41 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_AngRange_CB'},...
'String','4',...
'Style','edit',...
'Tooltip','The pole will be searched arround it''s starting angle +/- half this range',...
'Tag','edit_AngRange');

uicontrol('Parent',h1, 'Pos',[171 165 35 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_nInt_CB'},...
'String','21',...
'Style','edit',...
'Tooltip','The range parameters are divided into this number of intervals steps',...
'Tag','edit_nInt_lon');

uicontrol('Parent',h1, 'Pos',[171 138 35 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_nInt_CB'},...
'String','21',...
'Style','edit',...
'Tooltip','The range parameters are divided into this number of intervals steps',...
'Tag','edit_nInt_lat');

uicontrol('Parent',h1, 'Pos',[171 110 35 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_nInt_CB'},...
'String','21',...
'Style','edit',...
'Tooltip','The range parameters are divided into this number of intervals steps',...
'Tag','edit_nInt_ang');

uicontrol('Parent',h1,'Pos',[20  114 72  15],'HorizontalAlignment','left','Str','Angular Range','Style','text');
uicontrol('Parent',h1,'Pos',[136 329 41  16],'FontSize',10,'FontWeight','bold','String','OR','Style','text');
uicontrol('Parent',h1,'Pos',[97  386 81  16],'FontSize',10,'String','First Line','Style','text');
uicontrol('Parent',h1,'Pos',[348 386 81  16],'FontSize',10,'String','Second Line','Style','text');
uicontrol('Parent',h1,'Pos',[248 286 154 16],'FontSize',10,'String','Starting Pole Section','Style','text','Tag','txtSP');
uicontrol('Parent',h1,'Pos',[247 405 113 16],'FontSize',10,'String','Data Section','Style','text','Tag','txtDS');
uicontrol('Parent',h1,'Pos',[24  273 51  15],'String','Longitude','Style','text');
uicontrol('Parent',h1,'Pos',[116 273 51  15],'String','Latitude','Style','text');
uicontrol('Parent',h1,'Pos',[205 273 51  15],'String','Angle','Style','text');
uicontrol('Parent',h1,'Pos',[248 199 122 16],'FontSize',10,'Str','Compute Section','Style','text','Tag','txtCS');

uicontrol('Parent',h1, 'Pos',[352 173 61 21],...
'BackgroundColor',[1 1 1],...
'Style','edit',...
'Tooltip','Computed Euler pole Longitude',...
'Tag','edit_pLon_fim');

uicontrol('Parent',h1, 'Pos',[352 145 61 21],...
'BackgroundColor',[1 1 1],...
'Style','edit',...
'Tooltip','Computed Euler pole Latitude',...
'Tag','edit_pLat_fim');

uicontrol('Parent',h1, 'Pos',[352 117 61 21],...
'BackgroundColor',[1 1 1],...
'Style','edit',...
'Tooltip','Computed Euler pole Angle',...
'Tag','edit_pAng_fim');

uicontrol('Parent',h1,'Pos',[298 176 51 15],'HorizontalAlignment','left','String','Longitude','Style','text');
uicontrol('Parent',h1,'Pos',[298 147 47 17],'HorizontalAlignment','left','String','Latitude','Style','text');
uicontrol('Parent',h1,'Pos',[298 121 46 15],'HorizontalAlignment','left','String','Angle','Style','text');

uicontrol('Parent',h1,'Pos',[431 117 61 21],...
'BackgroundColor',[1 1 1],...
'Style','edit','Tooltip','Residue of the cost function',...
'Tag','edit_BFresidue');

uicontrol('Parent',h1,'Pos',[433 140 57 15],'String','BF Residue','Style','text');

uicontrol('Parent',h1, 'Pos',[430 161 61 21],...
'BackgroundColor',[1 1 1],...
'Style','edit',...
'Tooltip','Starting residue (starting pole) of the cost function',...
'Tag','edit_InitialResidue');

uicontrol('Parent',h1, 'Pos',[432 184 60 15],'String','St Residue','Style','text');

uicontrol('Parent',h1, 'Pos',[20 90 110 15],...
'Call',{@compute_euler_uiCB,h1,'check_hellinger_CB'},...
'String','Hellinger method',...
'Style','checkbox',...
'Tooltip','Use the Hellinger method',...
'Tag','check_hellinger');


uicontrol('Parent',h1, 'Pos',[170 61 225 21],...
'BackgroundColor',[1 1 1],...
'Call',{@compute_euler_uiCB,h1,'edit_err_file_CB'},...
'String','',...
'Style','edit',...
'HorizontalAlignment','left',...
'Tooltip','Name of residues grid. Leave blank if not wanted',...
'Tag','edit_err_file');

uicontrol('Parent',h1, 'Pos',[394 61 21 21],...
'Call',{@compute_euler_uiCB,h1,'push_err_file_CB'},...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'Tag','push_err_file')

uicontrol('Parent',h1, 'Pos',[220 84 155 18],...
'FontSize',10,...
'String','Residues grid (optional)',...
'Style','text')

uicontrol('Parent',h1, 'Pos',[431 79 75 21],...
'Call',{@compute_euler_uiCB,h1,'radio_netcdf_CB'},...
'String','3d netCDF',...
'Style','radiobutton',...
'Tooltip','Save residues grid as a 3D netCDF file',...
'Tag','radio_netcdf')

uicontrol('Parent',h1, 'Pos',[431 58 60 21],...
'Call',{@compute_euler_uiCB,h1,'radio_VTK_CB'},...
'String','3d VTK',...
'Style','radiobutton',...
'Tooltip','Save residues grid as a 3D VTK file',...
'Tag','radio_VTK')

uicontrol('Parent',h1, 'Pos',[10 14 231 16],...
'BackgroundColor',[0.9 0.9 0.9],...
'Enable','inactive',...
'Style','slider',...
'Tag','slider_wait');

uicontrol('Parent',h1, 'Pos',[241 12 37 19],...
'Call',{@compute_euler_uiCB,h1,'push_stop_CB'},...
'String','STOP','Tag','push_stop');

uicontrol('Parent',h1,'Pos',[435 10 76 21],...
'Call',{@compute_euler_uiCB,h1,'push_compute_CB'},...
'String','Compute','Tag','push_compute');

function compute_euler_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

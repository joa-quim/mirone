function varargout = compute_euler(varargin)
% Helper window to calculate Euler poles but that works also as a command line tool
%
% To work from the command line VARARGIN must have (and this is not a complete list)
% (isoc1, isoc2, plon, plat, pang, [[-DLonRange/LatRange/AngRange], [-ILonInt/LatInt/AngInt or -I<int>], ...
%	[-E<resid_fname>], [-Rx_min/x_max/y_min/y_max/ang_min/ang_max], [-V]])
% 
% Where isoc1 & isoc2 can either be a a nx2 matrix with the isochrons coordinates or files names
% PLON, PLAT, PANG are the parameters of the initial Euler pole.
% See function parse_noGUI for the definition of the remsining optional arguments

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

% $Id: compute_euler.m 4566 2014-09-23 16:49:31Z j $

	if (isempty(varargin) || (numel(varargin) >= 2 && numel(varargin) <= 5))
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

		if (ischar(varargin{1})),		handles.isoca1 = le_fiche(varargin{1});
		else							handles.isoca1 = varargin{1};
		end
		
		if (~handles.IamRegion)				% Otherwise those were already set inside parse_noGUI()
			handles.pLon_ini = varargin{3};
			handles.pLat_ini = varargin{4};
			handles.pAng_ini = varargin{5};
		end

		if (isempty(handles.noise))			% "Regular" mode
			if (ischar(varargin{2})),		handles.isoca2 = le_fiche(varargin{2});
			else							handles.isoca2 = varargin{2};
			end
			% The distmin MEX algo relies on the assumption that both lines have vertex growing in
			% the same sense. That is, Lat is increasing or decreasing for both (we test only Lat)
			Dlat1 = handles.isoca1(end,2) - handles.isoca1(1,2);
			Dlat2 = handles.isoca2(end,2) - handles.isoca2(1,2);
			if (sign(Dlat1) ~= sign(Dlat2))		% Lines have oposite orientation. Reverse one of them.
				handles.isoca2 = handles.isoca2(end:-1:1,:);
			end

			[pLon, pLat, pAng, resid] = calca_pEuler(handles, true, false);
		else
			[rlon, rlat] = rot_euler(handles.isoca1(:,1),handles.isoca1(:,2),handles.pLon_ini,handles.pLat_ini,handles.pAng_ini,-1);
			%alea = (-1 + 2 * rand(numel(rlon), 2)) * handles.noise;
			alea = randn(numel(rlon), 2) * handles.noise;
			handles.isoca2 = [(rlon + (alea(:, 1)))	(rlat + (alea(:, 2)))];
			[pLon, pLat, pAng, resid] = calca_pEuler(handles, true, false);
		end
		if (nargout),	varargout{1} = [pLon, pLat, pAng, resid];	end
		return
	end

	hObject = figure('Vis','off');
	compute_euler_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	% Initialize those
	handles.hCallingFig = [];
	handles.LonRange = 30;
	handles.LatRange = 20;
	handles.AngRange = 1;
	handles.nInt_lon = 61;
	handles.nInt_lat = 41;
	handles.nInt_ang = 21;
	handles.isoc_age = 0;
	handles.isoca1   = [];
	handles.isoca2   = [];
	handles.pLon_ini = [];
	handles.pLat_ini = [];
	handles.pAng_ini = [];
	handles.do_graphic = false;
	handles.DP_tol = 0.05;
	handles.residGrdName = [];
	handles.is_spheric = true;
	set(handles.slider_wait,'Max',handles.nInt_lon)

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
	handles.home_dir = handMir.home_dir;
	handles.path_continent  = [handMir.home_dir filesep 'continents' filesep];
	handles.deflation_level = handMir.deflation_level;
	handles.IamCompiled = handMir.IamCompiled;		% Need to know due to crazy issue of nc_funs

	str = sprintf(['The range interval is divided into this number of equally spaced points\n' ...
		'Alternatively use the form N*Delta (e.g. 100*0.1) to set up both range and resolution\n' ...
		'Actual point spacing is = 0.5']);
	set(handles.edit_nInt, 'Tooltip', str)

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, [handles.txtSP handles.txtDS handles.txtCS])
	%------------- END Pro look (3D) -----------------------------------------------------

	% Add this figure handle to the carra?as list
	plugedWin = getappdata(handMir.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handMir.figure1,'dependentFigs',plugedWin);

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
                    'ListSize',[450 300],'ListString',s);

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
function push_reciclePole_CB(hObject, handles)
% Get the parameters of computed pole and put them on starting pole edit boxes
	set(handles.edit_pLon_fim, 'Str', get(handles.edit_pLon_ini, 'Str'))
	set(handles.edit_pLat_fim, 'Str', get(handles.edit_pLat_ini, 'Str'))
	set(handles.edit_pAng_fim, 'Str', get(handles.edit_pAng_ini, 'Str'))

% -------------------------------------------------------------------------------------
function edit_LonRange_CB(hObject, handles)
	handles.LonRange = str2double(get(hObject,'String'));
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_nInt_CB(hObject, handles)
	if (~get(handles.check_hellinger,'Val'))
		str = get(hObject,'Str');
		ind = strfind(str, '*');
		tag = get(hObject,'UserData');
		tooltip = dataread('string',get(hObject,'Tooltip'),'%s','delimiter','\n');
		tooltip = sprintf('%s\n%s',tooltip{1},tooltip{2});		% Retain only first 2 lines, third will be updated here
		if (strcmp(tag,'lon')),			rang = get(handles.edit_LonRange, 'Str');
		elseif (strcmp(tag,'lat')),		rang = get(handles.edit_LatRange, 'Str');
		else							rang = get(handles.edit_AngRange, 'Str');
		end
		if (~isempty(ind))
			nInt = sscanf(str(1:ind(1)-1),'%d');
			d = sscanf(str(ind(1)+1:end),'%f');
			if (isnan(nInt) || isnan(d)),	set(hObject,'Str','3'),		return,		end
			if (rem(nInt,2)),	nInt = nInt + 1;	end		% We want an EVEN number of intervals
			nPts = nInt + 1;		% But this one must be ODD
			rang = sprintf('%.6g', nInt * d);
		else
			nPts = abs(sscanf(str,'%d'));
			if (~rem(nPts,2)),		nPts = nPts + 1;	end	% We want an ODD number of points
			d = str2double(rang) / (nPts - 1);
		end
		set(hObject,'Tooltip', sprintf('%s\nActual point spacing is = %.6g', tooltip, d))
		if (strcmp(tag,'lon')),			handles.nInt_lon = nPts;	set(handles.edit_LonRange, 'Str', rang)
		elseif (strcmp(tag,'lat')),		handles.nInt_lat = nPts;	set(handles.edit_LatRange, 'Str', rang)
		else							handles.nInt_ang = nPts;	set(handles.edit_AngRange, 'Str', rang)
		end
		set(handles.slider_wait,'Max',handles.nInt_lon)
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
		h_mir_lines = findobj(handles.hCallingFig,'Type','line');	% Fish all objects of type line in Mirone figure
		if (isempty(h_mir_lines))									% We don't have any lines
			str = ['If you hited this button on purpose, than you deserve the following insult.',...
					'You #!|"*!%!?~^)--$&.',... 
					'THERE ARE NO LINES IN THAT FIGURE.'];
			errordlg(str,'Chico Clever');   set(hObject,'Value',0);     return;
		end
		if (length(h_mir_lines) == 1)								% We don't have at least two lines
			str = ['If you hited this button on purpose, than you deserve the following insult.',...
					'You -$&#!*!%!?~^)-.|"/',... 
					'THERE IS ONLY ONE LINE IN THAT FIGURE.'];
			errordlg(str,'Chico Clever');   set(hObject,'Value',0);     return;
		end

		% The above test is not enough. For exemple, coastlines are not eligible neither,
		% but is very cumbersome to test all the possibilities of pure non-eligible lines.
		set(handles.hCallingFig,'pointer','crosshair')
		h_line1 = get_polygon(handles.hCallingFig);		% Get first line handle
		if (~isempty(h_line1))
			x = get(h_line1,'XData');		y = get(h_line1,'YData');
			handles.isoca1 = [x(:) y(:)];
			set(handles.edit_first_file,'String','Got left line','FontAngle','italic')
		else
			handles.isoca1 = [];
			set(handles.edit_first_file,'String','')
		end
		h_line2 = get_polygon(handles.hCallingFig);		% Get second line handle
		if (~isempty(h_line2))
			x = get(h_line2,'XData');		y = get(h_line2,'YData');
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
			handles.hLines = [NaN NaN];
		else
			handles.hLines = [h_line1 h_line2];
			handles.do_graphic = true;
		end
		set(hObject,'Value',0)
		figure(handles.figure1)			% Bring this figure to front again
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
		set(handles.edit_nInt,'Tooltip', sprintf(['Tolerance used to break up the isochron into\n' ...
				'linear chunks (the Heillinger segments).\n' ...
				'The units of the tolerance are degrees\n', ...
				'of arc on the surface of a sphere']))
	else
		set([handles.edit_LonRange handles.edit_LatRange handles.edit_AngRange],'Enable','on')
		set([handles.edit_nInt_lat handles.edit_nInt_ang],'Vis','on')
		set(handles.textNint,'String','N Intervals')
		str = sprintf(['The range interval is divided into this number of equally spaced points\n' ...
			'Alternatively use the form N*Delta (e.g. 100*0.1) to set up both range and resolution']);
		set(handles.edit_nInt,'Str', handles.nInt,'Tooltip', str)		% and so we loose the resolution info
	end

% -----------------------------------------------------------------------------
function edit_err_file_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname)
		set([handles.radio_netcdf handles.radio_VTK],'Val',0)
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

% -------------------------------------------------------------------------------
function out = outward_displacement_treta(handles)
% See if the OPTcontrol.txt file has an outward displacement file entry
	opt_file = [handles.home_dir filesep 'data' filesep 'OPTcontrol.txt'];
	out = [];
	if ( exist(opt_file, 'file') == 2 )
		fid = fopen(opt_file, 'r');
		c = (fread(fid,'*char'))';      fclose(fid);
		lines = strread(c,'%s','delimiter','\n');   clear c fid;
		m = numel(lines);
		odFile = [];
		for (k = 1:m)
			if (~strncmp(lines{k},'MIR_OD',6)),	continue,	end
			odFile = ddewhite(lines{k}(8:end));
			if (exist(odFile,'file') ~= 2)
				errordlg(['Outward displacement file ' odFile ' does not exist. Ignoring OD request'],'Error')
				odFile = [];
			end
			break
		end
		if (~isempty(odFile))
			out = text_read(odFile);
			if (isempty(out) || size(out,2) == 1)
				errordlg(['File ' odFile ' doesn''t have any recognized nymeric data, or one column only.'],'Error');
			end
		end
	end

% -------------------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
% OK. See if we have all the information needed to compute the Euler pole

	set(handles.edit_BFresidue,'String','');		set(handles.edit_InitialResidue,'String','')
	set(handles.edit_pLon_fim,'String','');			set(handles.edit_pLat_fim,'String','')
	set(handles.edit_pAng_fim,'String','')

	if (isempty(handles.isoca1) || isempty(handles.isoca2))
		errordlg('Compute Euler pole with what? It would help if you provide me TWO lines.','Chico Clever')
		return
	else
		% Fish the polyline coordinates again so that eventual line edits are taken into account right away
		x = get(handles.hLines(1),'XData');		y = get(handles.hLines(1),'YData');
		handles.isoca1 = [x(:) y(:)];
		x = get(handles.hLines(2),'XData');		y = get(handles.hLines(2),'YData');
		handles.isoca2 = [x(:) y(:)];		
	end
	if (isempty(handles.pLon_ini) || isempty(handles.pLat_ini) || isempty(handles.pAng_ini))
		errordlg(['I need a first guess of the Euler pole you are seeking for.' ...
			'Pay attention to the "Starting Pole Section"'],'Error')
		return
	end

	% Check if we have an Outward Displacement request placed in OPTcontrol.txt. If not 'out' is []
	OD = outward_displacement_treta(handles);	

	% The distmin MEX algo relies on the assumption that both lines have vertex growing in
	% the same sense. That is, Lat is increasing or decreasing for both (we test only Lat)
	Dlat1 = handles.isoca1(end,2) - handles.isoca1(1,2);
	Dlat2 = handles.isoca2(end,2) - handles.isoca2(1,2);
	if (sign(Dlat1) ~= sign(Dlat2))			% Lines have oposite orientation. Reverse one of them.
		handles.isoca2 = handles.isoca2(end:-1:1,:);
	end

	if (~get(handles.check_hellinger,'Val'))		% Our method
		do_weighted = true;
		calca_pEuler(handles, do_weighted, true, OD);
	else											% Try with Hellinger's (pfiu)
		[pLon,pLat,pAng] = hellinger(handles.pLon_ini,handles.pLat_ini,handles.pAng_ini, handles.isoca1, handles.isoca2, handles.DP_tol);
		set(handles.edit_pLon_fim,'String',pLon);			set(handles.edit_pLat_fim,'String',pLat)
		set(handles.edit_pAng_fim,'String',pAng)
		if (handles.do_graphic)     % Create a empty line handle
			[rlon,rlat] = rot_euler(handles.isoca1(:,1),handles.isoca1(:,2),pLon, pLat, pAng,-1);
			hLine = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',rlon,'YData',rlat, ...
			'LineStyle','-.','LineWidth',2,'Tag','Fitted Line','Userdata',1);
			setappdata(hLine,'LineInfo','Fitted Line');
			draw_funs(hLine,'isochron',{'Fitted Line'})
		end
	end

	if (~isempty(get(handles.edit_pLon_fim,'Str')))	% If a new pole was computed set this button visible (driven by lazyness)
		set(handles.push_reciclePole,'Vis', 'on')
	else
		set(handles.push_reciclePole,'Vis', 'off')
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
	if (isa(numeric_data, 'cell'))
		numeric_data = numeric_data{1};
	end

% -------------------------------------------------------------------------------
function [polLon, polLat, polAng, area_f] = calca_pEuler(handles, do_weighted, isGUI, OD)
% ...
	D2R = pi / 180;		h_line = [];
	ecc = 0.0818191908426215;			% WGS84
	if (handles.do_graphic)				% Create a empty line handle
		h_line = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',[],'YData',[], ...
			'LineStyle','-.','LineWidth',2,'Tag','Fitted Line','Userdata',1);
	end

	isoca1 = handles.isoca1 * D2R;		% A maluca
	isoca2 = handles.isoca2 * D2R;		% A fixa

	% Check if we have an Outward Displacement (in)correction request
	if (nargin == 4 && ~isempty(OD))	% If yes OD has ods in km
		R2D = 180 / pi;
	 	[rlon1,rlat1] = rot_euler(isoca1(:,1), isoca1(:,2),handles.pLon_ini*D2R,handles.pLat_ini*D2R,handles.pAng_ini*D2R,'rad',-1);
	 	[rlon2,rlat2] = rot_euler(isoca2(:,1), isoca2(:,2),handles.pLon_ini*D2R,handles.pLat_ini*D2R,-handles.pAng_ini*D2R,'rad',-1);
		shifts_1 = interp1(OD(:,1)*D2R, OD(:,2), rlat1, 'linear', 0);		% Shifts in km
		shifts_2 = interp1(OD(:,1)*D2R, OD(:,2), rlat2, 'linear', 0);

		isoca1 = isoca1 * R2D;			isoca2 = isoca2 * R2D;
		rlat1  = rlat1 * R2D;			rlon1  = rlon1 * R2D;
		rlat2  = rlat2 * R2D;			rlon2  = rlon2 * R2D;
		[s,a12] = vdist(isoca1(:,2),isoca1(:,1),rlat1,rlon1);	% coords must be in degrees here
		for (k = 1:numel(rlon1))
			[rlat1(k), rlon1(k)] = vreckon(rlat1(k), rlon1(k), shifts_1(k)*1000, a12(k), 1);
		end
		[s,a12] = vdist(isoca2(:,2),isoca2(:,1),rlat2,rlon2);
		for (k = 1:numel(rlon2))
			[rlat2(k), rlon2(k)] = vreckon(rlat2(k), rlon2(k), shifts_2(k)*1000, a12(k), 1);
		end

		isoca1 = isoca1 * D2R;			isoca2 = isoca2 * D2R;	% Put them back in radians
		rlat1  = rlat1 * D2R;			rlon1  = rlon1 * D2R;
		rlat2  = rlat2 * D2R;			rlon2  = rlon2 * D2R;
		% Now invert rotations to get the "original" corrected positions
	 	[isoca1(:,1),isoca1(:,2)] = rot_euler(rlon1, rlat1, handles.pLon_ini*D2R,handles.pLat_ini*D2R,-handles.pAng_ini*D2R,'rad',-1);
	 	[isoca2(:,1),isoca2(:,2)] = rot_euler(rlon2, rlat2, handles.pLon_ini*D2R,handles.pLat_ini*D2R,handles.pAng_ini*D2R,'rad',-1);
	end

	isoca1(:,2) = atan2( (1-ecc^2)*sin(isoca1(:,2)), cos(isoca1(:,2)) );	% Lat da isoca1 geocentrica
	isoca2(:,2) = atan2( (1-ecc^2)*sin(isoca2(:,2)), cos(isoca2(:,2)) );	% Lat da isoca2 geocentrica

	% Compute distances between vertices of the moving isoc
	X = isoca2(:,1) .* cos(isoca2(:,2)) * 6371;
	Y = isoca2(:,2) * 6371;
	xd = diff(X);			yd = diff(Y);
	lenRot2 = sqrt(xd .* xd + yd .* yd);

 	[rlon,rlat] = rot_euler(isoca1(:,1), isoca1(:,2),handles.pLon_ini*D2R,handles.pLat_ini*D2R,handles.pAng_ini*D2R,'rad',0);
	X_rot = rlon .* cos(rlat) * 6371;		Y_rot = rlat * 6371;
	xd = diff(X_rot);		yd = diff(Y_rot);
	lenRot1 = sqrt(xd .* xd + yd .* yd);

	if (handles.is_spheric)
	 	area0 = distmin(isoca2(:,1), isoca2(:,2), lenRot2, rlon, rlat, lenRot1);
		%[dist1, segLen1, sum1] = distmin_o(isoca2(:,1), isoca2(:,2), rlon, rlat, lenRot1, do_weighted, 1e20);
		%[dist2, segLen2, sum2] = distmin_o(rlon, rlat, isoca2(:,1), isoca2(:,2), lenRot2, do_weighted, 1e20);
		%area0 = (sum1 + sum2) / 2;
	else
		errordlg('This branch - NON-SPHERICAL - is not working','Error'),	return
		%[dist1, segLen1, sum1] = distmin(X, Y, X_rot, Y_rot, lenRot1, do_weighted, 1e20);
		%[dist2, segLen2, sum2] = distmin(X_rot, Y_rot, X, Y, lenRot2, do_weighted, 1e20);
		%area0 = (sum1 + sum2) / 2;
	end
	%sum22 = weightedSum(dist2, segLen2, do_weighted);

	% See if there is a request to plot only the signed residues of the starting pole.
	if (get(handles.check_plotRes,'Val'))
		resid_along_isoca(handles, isoca2, lenRot2, rlon, rlat, lenRot1)
		return
	end

	if (handles.do_graphic)
		set(handles.edit_InitialResidue,'String',sprintf('%.4f', area0));    pause(0.01)
		set(h_line,'XData',rlon,'YData',rlat)
	end

	% Now comes the semi-brute force aproach to compute the pole
	dLon = handles.LonRange / 2;
	dLat = handles.LatRange / 2;
	dAng = handles.AngRange / 2;
	if (handles.nInt_lon > 1)
		p_lon = (handles.pLon_ini + linspace(-dLon,dLon, handles.nInt_lon)) * D2R;
	else
		p_lon = handles.pLon_ini * D2R;
	end
	if (handles.nInt_lat > 1)
		p_lat = (handles.pLat_ini + linspace(-dLat,dLat, handles.nInt_lat)) * D2R;
	else
		p_lat = handles.pLat_ini * D2R;
	end
	if (handles.nInt_ang > 1)
		p_omeg = (handles.pAng_ini + linspace(-dAng,dAng,handles.nInt_ang)) * D2R;
	else
		p_omeg = handles.pAng_ini * D2R;
	end

	% Sanitize p_lat so that it does not go out of N/S poles
	ind = ((p_lat > pi/2) | (p_lat < -pi/2));
	if (any(ind))
		p_lat(ind) = [];
		handles.nInt_lat = numel(p_lat);
		if (handles.do_graphic)
			set(handles.slider_wait,'Max',handles.nInt_lon)
		end
	end

	[polLon, polLat, polAng, area_f, resid] = ...
		fit_pEuler(handles, isoca1, isoca2, p_lon, p_lat, p_omeg, area0, h_line, lenRot1, lenRot2, do_weighted, isGUI);

	if (handles.do_graphic)
		setappdata(h_line,'LineInfo','Fitted Line');
		draw_funs(h_line,'isochron',{'Fitted Line'})
	end

	if (~isempty(resid) && isGUI)
		if (get(handles.radio_netcdf, 'Val'))
			write_netcdf(handles, p_lon/D2R, p_lat/D2R, p_omeg/D2R, resid)
		else
			write_vtk(handles, p_lon/D2R, p_lat/D2R, p_omeg/D2R, resid)
		end
	elseif (~isempty(resid))	% Command line runn.
		if (handles.isoc_age),	p_omeg = p_omeg / handles.isoc_age;		end		% Store angular velocity instead
		write_netcdf(handles, p_lon/D2R, p_lat/D2R, p_omeg/D2R, resid)
	end

% -------------------------------------------------------------------------------
function [lon_bf, lat_bf, omega_bf, area_f, resid] = ...
		fit_pEuler(handles, isoca1, isoca2, p_lon, p_lat, p_omeg, area0, h_line, lenRot1,  lenRot2, do_weighted, isGUI)
% This is the function that does the real work of computing the best Euler pole
% that fits the lines "isoca1" "isoca2"
% ISOCA1 -> A maluca	- Normally they are already converted to geocentrics
% ISOCA2 -> Fix one
	lon_bf = [];		lat_bf = [];	omega_bf = [];
	D2R = pi / 180;
	ecc = 0.0818191908426215;			% WGS84
	nLon = numel(p_lon);	nLat = numel(p_lat);	nAng = numel(p_omeg);

	if (~isempty(handles.residGrdName))	% Store all residues in a 3D array
		save_resid = true;		resid = zeros(nLat,nLon,nAng) * NaN;	testResidue = 1e20;
	else
		save_resid = false;		resid = [];		testResidue = area0;
	end

	if (~isGUI)
		strAdv = fprintf('\n%000\t Lons rots out of %d', nLon);		% For the Text progressbar		
	end

	s_lat = sin(isoca1(:,2));			c_lat = cos(isoca1(:,2));

	X = isoca2(:,1) .* cos(isoca2(:,2)) * 6371;		Y = isoca2(:,2) * 6371;
	i_m = 1;	j_m = 1;	k_m = 0;
	for i = 1:nLon			% Loop over n lon intervals
		if (isGUI && get(handles.slider_wait,'Max') == 1)			% The STOP button was pushed. So, stop
			set(handles.slider_wait,'Max',handles.nInt_lon)			% Reset it for the next run
			area_f = -1;		resid = [];
			return
		end
		s_lon = sin(isoca1(:,1) - p_lon(i));	c_lon = cos(isoca1(:,1) - p_lon(i));
		cc = c_lat .* c_lon;
		for j = 1:nLat		% Loop over n lat intervals
			if (isGUI && get(handles.slider_wait,'Max') == 1)		% The STOP button was pushed. So, stop
				set(handles.slider_wait,'Max',handles.nInt_lon)		% Reset it for the next run
				area_f = -1;	resid = [];
				return
			end
			p_sin_lat = sin(p_lat(j));			p_cos_lat = cos(p_lat(j));
			tlon = atan2(c_lat .* s_lon, p_sin_lat * cc - p_cos_lat * s_lat);
			s_lat_ = p_sin_lat * s_lat + p_cos_lat * cc;
			c_lat_ = sqrt(1 - s_lat_ .* s_lat_);
			if (isGUI && ~rem(j,10)),		set(handles.slider_wait,'Value',i),		drawnow;	end	% Otherwise slider may stall
			for k = 1:nAng	% Loop over n omega intervals
				%[rlon,rlat] = rot_euler(isoca1(:,1),isoca1(:,2),p_lon(i),p_lat(j),p_omeg(k),'radians',-1);

                %lat = atan2( (1-ecc^2)*sin(isoca1(:,2)), cos(isoca1(:,2)) );
				%s_lat = sin(lat);						c_lat = cos(lat);
				%s_lon = sin(isoca1(:,1) - p_lon(i));	c_lon = cos(isoca1(:,1) - p_lon(i));
				%cc = c_lat .* c_lon;
				%p_sin_lat = sin(p_lat(j));				p_cos_lat = cos(p_lat(j));

				%tlon = atan2(c_lat .* s_lon, p_sin_lat * cc - p_cos_lat * s_lat);
				%s_lat_ = p_sin_lat * s_lat + p_cos_lat * cc;
				%c_lat_ = sqrt(1 - s_lat_ .* s_lat_);

				s_lon_ = sin(tlon + p_omeg(k));			c_lon_ = cos(tlon + p_omeg(k));
				cc_ = c_lat_ .* c_lon_;

				rlat = asin(p_sin_lat * s_lat_ - p_cos_lat * cc_);
				rlon = p_lon(i) + atan2(c_lat_ .* s_lon_, p_sin_lat * cc_ + p_cos_lat * s_lat_);

				%rlat = atan2( sin(rlat), (1-ecc^2)*cos(rlat) );		% Convert back to geodetic latitudes
				ind = (rlon > pi);					
				if (any(ind)),		rlon(ind) = rlon(ind) - 2*pi;	end

				if (handles.is_spheric)
				 	[area, lix, dists, weights] = distmin(isoca2(:,1), isoca2(:,2), lenRot2, rlon, rlat, lenRot1);
					%[dist1, segLen1, sum1] = distmin_o(isoca2(:,1), isoca2(:,2), rlon, rlat, lenRot1, do_weighted, 1e20);
					%[dist2, segLen2, sum2] = distmin_o(rlon, rlat, isoca2(:,1), isoca2(:,2), lenRot2, do_weighted, 1e20);
					%area = (sum1 + sum2) / 2;
				else
					X_rot = rlon .* cos(rlat) * 6371;		Y_rot = rlat * 6371;
					[dist1, segLen1, sum1] = distmin(X, Y, X_rot, Y_rot, lenRot1, do_weighted, testResidue);
					[dist2, segLen2, sum2] = distmin(X_rot, Y_rot, X, Y, lenRot2, do_weighted, testResidue);
					area = (sum1 + sum2) / 2;
				end

				if (area < area0)
					area0 = area;
					i_m = i;	j_m = j;	k_m = k;
					if (handles.do_graphic)
						rlat = atan2(sin(rlat), (1-ecc^2)*cos(rlat));		% Convert back to geodetic latitudes
						set(h_line,'XData',rlon/D2R,'YData',rlat/D2R)
						pause(0.01)
					end
					if (isGUI)
						set(handles.edit_pLon_fim,'String',sprintf('%.3f', p_lon(i) / D2R))
						set(handles.edit_pLat_fim,'String',sprintf('%.3f', p_lat(j) / D2R))
						set(handles.edit_pAng_fim,'String',sprintf('%.3f', p_omeg(k) / D2R))
						set(handles.edit_BFresidue,'String',sprintf('%.4f', area))
						% Compute also the Weighted Standard Deviation
						n_non_zero = numel(find(weights));
						wstd = sqrt(sum((weights .* (dists - area)).^2) / (sum(weights) * (n_non_zero - 1) / n_non_zero));
						set(handles.edit_BFresidue,'Tooltip',sprintf('Average residue in km\nWeighted STD = %.3f',wstd))
					end
				end

				if (save_resid),	resid(j, i, k) = area;		end

			end
		end
		if (~isGUI)			% Make a Text progressbar
			fprintf(repmat('\b',1,numel(strAdv)))		% Come back to the begining of line
			strAdv = sprintf('%0.3d\t Lons rots out of %d', i, nLon);
			fprintf(strAdv)
		end
	end

	if (~isGUI),	fprintf('\n');		end
	if (isGUI),		set(handles.slider_wait,'Value',0),		end      % Reset it for the next run
	lon_bf = p_lon(i_m) / D2R;
	lat_bf = p_lat(j_m) / D2R;
	if (k_m > 0)
		omega_bf = p_omeg(k_m) / D2R;
	else			% Nothing better than initial angle was found. Return that exact value
		omega_bf = handles.pAng_ini;
	end
	area_f = area0;

% -------------------------------------------------------------------------------
function resid_along_isoca(handles, isoca2, lenRot2, rlon, rlat, lenRot1)
% ...
	R2D = 180 / pi;		ecc = 0.0818191908426215;			% WGS84
	[area0, xy, distsA, weights] = distmin(isoca2(:,1), isoca2(:,2), lenRot2, rlon, rlat, lenRot1);
	ind = (weights < 0.25);		% Ad-hoc criterium
	distsA(ind) = [];	xy(ind,:) = [];		isoBak = isoca2;	isoca2(ind,:) = [];

	% Compute the azims of points and their the ones that are closest on the otehr line, next
	% compare them with the azimuth given by the Euler pole and thus decide the ones on left and right
	azPtsNear = azimuth_geo(xy(:,2),xy(:,1),isoca2(:,2),isoca2(:,1),'rad') * R2D;
	azEuler   = compute_EulerAzim(xy(:,2),xy(:,1), handles.pLat_ini/R2D,handles.pLon_ini/R2D) * R2D;
	difa = abs(azPtsNear - azEuler);
	inds = zeros(numel(azPtsNear),1);
	inds(difa < 60) = 1;	inds(difa > 120) = -1;
	distsA = distsA .* inds;
	distsA(~inds) = [];		xy(~inds,:) = [];		% Those should be points edging Fracture Zones
	latA = atan2( sin(xy(:,2)), (1-ecc^2)*cos(xy(:,2)) ) * R2D;		% Convert back to geodetic latitudes
	
	% Now again but swaping order of lines
	[area0, xy, distsB, weights] = distmin(rlon, rlat, lenRot1, isoBak(:,1), isoBak(:,2), lenRot2);
	ind = (weights < 0.25);
	distsB(ind) = [];	xy(ind,:) = [];		rlon(ind,:) = [];	rlat(ind,:) = [];
	azPtsNear = azimuth_geo(xy(:,2),xy(:,1),rlat,rlon,'rad') * R2D;
	azEuler   = compute_EulerAzim(xy(:,2),xy(:,1), handles.pLat_ini/R2D,handles.pLon_ini/R2D) * R2D;
	difa = abs(azPtsNear - azEuler);
	inds = zeros(numel(azPtsNear),1);
	inds(difa < 60) = -1;	inds(difa > 120) = 1;	% NOTE: swapp signs because we reversed order of lines
	distsB = distsB .* inds;
	distsB(~inds) = [];		xy(~inds,:) = [];
	latB = atan2( sin(xy(:,2)), (1-ecc^2)*cos(xy(:,2)) ) * R2D;

	% The correct thing to do would be to interpolate along ridge to get a common reference
	% but we'll do the crude aproximation of using Lat as the X axis (good enough for the Atlantic)
	lat = [latA; latB];		dists = ([distsA; distsB]);
	[lat, ind] = sort(lat);
	dists = dists(ind);

	ecran(lat, dists, 'Residues along isoc')

% -----------------------------------------------------------------------------------------
function azim = compute_EulerAzim(alat,alon,plat,plon)
% alat & alon are the point coords. plat, plon are the pole parameters (All in radians)

	x = cos(plat).*sin(plon).*sin(alat) - cos(alat).*sin(alon).*sin(plat);    % East vel
	y = cos(alat).*cos(alon).*sin(plat) - cos(plat).*cos(plon).*sin(alat);    % North vel
	z = cos(plat).*cos(alat).*sin(alon-plon);
	vlon = -sin(alon).*x + cos(alon).*y;
	vlat = -sin(alat).*cos(alon).*x-sin(alat).*sin(alon).*y + cos(alat).*z;
	azim = pi/2 - atan2(vlat,vlon);
	if (azim < 0),		azim = azim + 2*pi;		end		% Give allways the result in the 0-360 range

% -------------------------------------------------------------------------------
function [dist, segLen] = distmin_(lon, lat, r_lon, r_lat, lengthsRot, do_weighted, lastResidue)
% Compute the shortest distance between each point in (lon,lat) and the polyline (r_lon,r_lat)
% SEGLEN holds the segment length of polyline (r_lon,r_lat) corresponding to elements of DIST
% At each quater of the total number of points in (lon.lat) we check if the current residue is
% already larger than last better residue. If yes, we stop since for sure this rotation is worst.
% The above test is short-circuited (not carried out) when saving the residues to a grid.
% Angles are already in radians.
	
	r_lon = r_lon(:)';		r_lat = r_lat(:)';      % Make sure they are row vectors
	
	%eps1 = 1e-1;			eps2 = eps1 * eps1;
	n_pt = numel(lon);		n_pt_rot = numel(r_lon);
	dist = zeros(n_pt,1);
	segLen = ones(n_pt,1);
	outliners = false(n_pt,1);      % To remove vertex where lines do not intersect

	chunks = round([(n_pt * 0.25) (n_pt / 2) (n_pt * 0.75) (n_pt * 0.9) (n_pt+1)]);	% Checkpoints where current res is checked against min res
	for (k = 1:n_pt)				% Loop over vertices of fixed isoc
		Dsts = ((lon(k)-r_lon).^2 + (lat(k)-r_lat).^2);		% distances^2 from pt k in (lon,lat) to all pts in (r_lon,r_lat)
		[D,ind] = min(Dsts);		% Compute the min dist from pt k to vertex of line and its location
		if (ind == 1 || ind == n_pt_rot && n_pt > 4)		% This point is outside the lines intersection. Flag it to die.
			outliners(k) = true;						% Actually we waste all points that are closest to the rotated line end points
			continue
		end

%  		P  = [lon(k) lat(k)];
% 		Q1 = [r_lon(ind-1) r_lat(ind-1)];		Q2 = [r_lon(ind) r_lat(ind)];
% 		D1 = abs(det([Q2-Q1; P-Q1])) / norm(Q2-Q1); % for row vectors.
% 		Q1 = [r_lon(ind) r_lat(ind)];			Q2 = [r_lon(ind+1) r_lat(ind+1)];
% 		D2 = abs(det([Q2-Q1; P-Q1])) / norm(Q2-Q1);

		% This is the same as above but runs at least twice fast
 		r_lon_t = r_lon(ind-1:ind+1);
 		r_lat_t = r_lat(ind-1:ind+1);
 		Q1(1) = r_lon_t(1);			Q1(2) = r_lat_t(1);
		Q2(1) = r_lon_t(2);			Q2(2) = r_lat_t(2);
		DQ(1) = Q2(1) - Q1(1);		DQ(2) = Q2(2) - Q1(2);
		D1 = abs(DQ(1)*(lat(k)-Q1(2)) - DQ(2)*(lon(k)-Q1(1))) / sqrt(DQ(1)*DQ(1) + DQ(2)*DQ(2));
 		Q1(1) = r_lon_t(2);			Q1(2) = r_lat_t(2);
		Q2(1) = r_lon_t(3);			Q2(2) = r_lat_t(3);
		DQ(1) = Q2(1) - Q1(1);		DQ(2) = Q2(2) - Q1(2);
		D2 = abs(DQ(1)*(lat(k)-Q1(2)) - DQ(2)*(lon(k)-Q1(1))) / sqrt(DQ(1)*DQ(1) + DQ(2)*DQ(2));

		if (D1 < D2)				% Store shortest dist from point k to polyline (r_lon,r_lat)
			dist(k) = D1;	segLen(k) = lengthsRot(ind-1);	% Store the segment length that is closest to pt k
		else
			dist(k) = D2;	segLen(k) = lengthsRot(ind);
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
	ind = (dists ~= 0);
	dists = dists(ind);		segLen = segLen(ind);
	if (do_weighted)
		weights = ones(numel(segLen), 1);
		ind = (segLen > 50 & segLen < 80);
		weights(ind) = 0.25;
		weights(segLen >= 80) = 0;
		soma = sum( dists .* weights ) / sum(weights);		% weighted sum
	else
		soma = sum( dists ) / numel(dists);
	end

% -----------------------------------------------------------------------------------------
function write_netcdf(handles, lon, lat, ang, resid)
% Write the 3D matrix in netCDF

	nz = numel(ang);

	% No problem in changing handles here as it will live only locally to this function
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
	fprintf(fid, 'BINARY\nDATASET RECTILINEAR_GRID\n');
	fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
	fprintf(fid, 'X_COORDINATES %d float\n', nx);
	fwrite(fid, lon, 'real*4');
	fprintf(fid, 'Y_COORDINATES %d float\n', ny);
	fwrite(fid, lat, 'real*4');
	fprintf(fid, 'Z_COORDINATES %d float\n', nz);
	fwrite(fid, ang, 'real*4');
	fprintf(fid, 'POINT_DATA %d\n', nx * ny * nz);
	fprintf(fid, 'SCALARS dono float 1\nLOOKUP_TABLE default\n');

	for (k = 1:nz)
		Z = single(resid(:,:,k));
		fwrite(fid, Z(:), 'real*4');
	end
	fclose(fid);

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
% Use -D or -R in alternative. The first will center the searching region about the provided pole position
% -R on the other hand impose a fixed searching region and a fake starting pole will be computed at its center.
%	This option is usefull when we want to compute residues on the exact same region.
% -I can be given as the numbert of intervals (old way), or as increments.

	handles.LonRange = 30;
	handles.LatRange = 20;
	handles.AngRange = 1;
	handles.nInt_lon = 61;
	handles.nInt_lat = 41;
	handles.nInt_ang = 21;
	handles.do_graphic = false;
	handles.is_spheric = true;
	handles.residGrdName = [];
	handles.noise = [];
	handles.check_plotRes = [];
	handles.isoc_age = 0;
	handles.IamCompiled = false;	% Até ver
	handles.IamRegion = false;		% Will be set to true if -R is provided	
	handles.deflation_level = 5;	% Since we don't have access to Mirone handles we use this as default

	msg = [];	verbose = false;

	if (numel(varargin) < 6),	return,		end		% No option provided. Just use the defaults

	for (k = 6:numel(varargin))
		switch (varargin{k}(2))
			case 'H'
				msg = sprintf(['[-DLonRange/LatRange/AngRange] [-ILonInt/LatInt/AngInt or -Iint] [-E<resid_fname>]\n',...
						'[-A<age>] [-Rx_min/x_max/y_min/y_max/ang_min/ang_max] [-V]']);
				return
			case 'A'
				handles.isoc_age = sscanf(varargin{k}(3:end), '%f');
				if (handles.isoc_age <= 0 || isnan(handles.isoc_age))
					msg = 'Error using option -A [-A<isoc_age>]';		return
				end
			case 'C'
				handles.is_spheric = false;			% Not used now
			case 'D'
				ind = strfind(varargin{k},'/');
				if (numel(ind) ~= 2)
					msg = 'Error using option -D [-DLonRange/LatRange/AngRange]';		return
				end
				handles.LonRange = sscanf(varargin{k}(3:ind(1)-1), '%f');
				handles.LatRange = sscanf(varargin{k}(ind(1)+1:ind(2)-1), '%f');
				handles.AngRange = sscanf(varargin{k}(ind(2)+1:end), '%f');

			case 'I'
				ind = strfind(varargin{k},'/');
				if (~isempty(ind) && numel(ind) ~= 2)
					msg = 'Error using option -I: -ILonInt/LatInt/AngInt or -Iint';		return
				end
				if (isempty(ind))
					handles.nInt_lon = sscanf(varargin{k}(3:end), '%f');
					handles.nInt_lat = handles.nInt_lon;	handles.nInt_ang = handles.nInt_lon;
					if (handles.nInt_lon > 1 && ~rem(handles.nInt_lon,2)),	handles.nInt_lon = handles.nInt_lon + 1;	end
					if (handles.nInt_lat > 1 && ~rem(handles.nInt_lat,2)),	handles.nInt_lat = handles.nInt_lat + 1;	end
					if (handles.nInt_ang > 1 && ~rem(handles.nInt_ang,2)),	handles.nInt_ang = handles.nInt_ang + 1;	end
				else
					handles.nInt_lon = sscanf(varargin{k}(3:ind(1)-1), '%f');
					handles.nInt_lat = sscanf(varargin{k}(ind(1)+1:ind(2)-1), '%f');
					handles.nInt_ang = sscanf(varargin{k}(ind(2)+1:end), '%f');					
				end

			case 'N'
				handles.noise = sscanf(varargin{k}(3:end), '%f');

			case 'E'
				handles.residGrdName = varargin{k}(3:end);

			case 'R'
				ind = strfind(varargin{k},'/');
				if (numel(ind) ~= 5)
					msg = 'Error using option -R [-Rx_min/x_max/y_min/y_max/ang_min/ang_max]';		return
				end
				% OK, here -R was used and we must compute a fake -D.
				x_min = sscanf(varargin{k}(3:ind(1)-1),'%f');
				x_max = sscanf(varargin{k}(ind(1)+1:ind(2)-1),'%f');
				y_min = sscanf(varargin{k}(ind(2)+1:ind(3)-1),'%f');
				y_max = sscanf(varargin{k}(ind(3)+1:ind(4)-1),'%f');
				a_min = sscanf(varargin{k}(ind(4)+1:ind(5)-1),'%f');
				a_max = sscanf(varargin{k}(ind(5)+1:end),'%f');
				handles.LonRange = x_max - x_min;
				handles.LatRange = y_max - y_min;
				handles.AngRange = a_max - a_min;
				handles.pLon_ini = x_min + handles.LonRange / 2;
				handles.pLat_ini = y_min + handles.LatRange / 2;
				handles.pAng_ini = a_min + handles.AngRange / 2;
				handles.IamRegion = true;

			case 'V'
				verbose = true;
		end
	end

	if (handles.nInt_lon <= 1)		% -I was given as an increment instead of n of points
		handles.nInt_lon = round(handles.LonRange / handles.nInt_lon) + 1;
		handles.nInt_lat = round(handles.LatRange / handles.nInt_lat) + 1;
		handles.nInt_ang = round(handles.AngRange / handles.nInt_ang) + 1;
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

uicontrol('Parent',h1,'Pos',[10 50  501 157],'Style','frame')
uicontrol('Parent',h1,'Pos',[10 230 501  64],'Style','frame');
uicontrol('Parent',h1,'Pos',[10 313 501 101],'Style','frame');

uicontrol('Parent',h1, 'Pos',[20 363 211 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'Style','edit','Tag','edit_first_file');

uicontrol('Parent',h1, 'Pos',[231 363 21 21],...
'Call',@compute_euler_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','...','Tag','push_first_file');

uicontrol('Parent',h1, 'Pos',[267 363 211 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'Style','edit','Tag','edit_second_file');

uicontrol('Parent',h1, 'Pos',[478 363 21 21],...
'Call',@compute_euler_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','...','Tag','push_second_file');

uicontrol('Parent',h1, 'Pos',[186 326 141 21],...
'Call',@compute_euler_uiCB,...
'String','Pick lines from Figure',...
'Style','togglebutton',...
'Tooltip','Allows you to mouse select the two lines from a Mirone figure',...
'Tag','toggle_pickLines');

uicontrol('Parent',h1, 'Pos',[20 249 61 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'Style','edit',...
'Tooltip','Start Euler pole Longitude',...
'Tag','edit_pLon_ini');

uicontrol('Parent',h1, 'Pos',[110 249 61 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'Style','edit',...
'Tooltip','Start Euler pole Latitude',...
'Tag','edit_pLat_ini');

uicontrol('Parent',h1, 'Pos',[200 249 61 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'Style','edit',...
'Tooltip','Start Euler pole Angle',...
'Tag','edit_pAng_ini');

uicontrol('Parent',h1, 'Pos',[289 248 121 21],...
'Call',@compute_euler_uiCB,...
'String','Poles selector',...
'Tooltip','Select a pole from the default list',...
'Tag','push_polesList');

uicontrol('Parent',h1, 'Pos',[440 248 21 21],...
'Call',@compute_euler_uiCB,...
'String','^',...
'Tooltip','Re-initialize with computed pole',...
'Vis', 'off',...
'Tag','push_reciclePole');

uicontrol('Parent',h1, 'Pos',[20 168 84 16],'HorizontalAlignment','left','Str','Longitude Range','Style','text');

uicontrol('Parent',h1, 'Pos',[104 165 41 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'String','30',...
'Style','edit',...
'Tooltip','The pole will be searched arround it''s starting longitude +/- half this range',...
'Tag','edit_LonRange');

uicontrol('Parent',h1, 'Pos',[161 186 60 15],'String','N Points','Style','text','Tag','textNint');
uicontrol('Parent',h1, 'Pos',[20 139 75 18],'HorizontalAlignment','left','Str','Latitude Range','Style','text');

uicontrol('Parent',h1, 'Pos',[104 138 41 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'String','20',...
'Style','edit',...
'Tooltip','The pole will be searched arround it''s starting latitude +/- half this range',...
'Tag','edit_LatRange');

uicontrol('Parent',h1, 'Pos',[104 110 41 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'String','1',...
'Style','edit',...
'Tooltip','The pole will be searched arround it''s starting angle +/- half this range',...
'Tag','edit_AngRange');

uicontrol('Parent',h1, 'Pos',[160 165 61 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'String','61',...
'Style','edit',...
'UserData','lon',...
'Tag','edit_nInt');

uicontrol('Parent',h1, 'Pos',[160 138 61 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'String','41',...
'Style','edit',...
'UserData','lat',...
'Tag','edit_nInt');

uicontrol('Parent',h1, 'Pos',[160 110 61 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'String','21',...
'Style','edit',...
'UserData','ang',...
'Tag','edit_nInt');

uicontrol('Parent',h1,'Pos',[20  114 72  15],'HorizontalAlignment','left','Str','Angular Range','Style','text');
uicontrol('Parent',h1,'Pos',[136 329 41  16],'FontSize',10,'FontWeight','bold','String','OR','Style','text','Tag','txtOR');
uicontrol('Parent',h1,'Pos',[97  386 81  16],'FontSize',10,'String','First Line','Style','text','Tag','txtFirstLine');
uicontrol('Parent',h1,'Pos',[348 386 81  16],'FontSize',10,'String','Second Line','Style','text','Tag','txtSecondLine');
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
'Call',@compute_euler_uiCB,...
'String','Hellinger method',...
'Style','checkbox',...
'Tooltip','Use the Hellinger method',...
'Tag','check_hellinger');

uicontrol('Parent',h1, 'Pos',[20 60 120 15],...
'String','Plot residues only',...
'Style','checkbox',...
'Tooltip','Plot the signed residues between isochrn and rotated isochron. No pole is computed',...
'Tag','check_plotRes');

uicontrol('Parent',h1, 'Pos',[170 61 225 21],...
'BackgroundColor',[1 1 1],...
'Call',@compute_euler_uiCB,...
'String','',...
'Style','edit',...
'HorizontalAlignment','left',...
'Tooltip','Name of residues grid. Leave blank if not wanted',...
'Tag','edit_err_file');

uicontrol('Parent',h1, 'Pos',[394 61 21 21],...
'Call',@compute_euler_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'Tag','push_err_file')

uicontrol('Parent',h1, 'Pos',[220 84 155 18],...
'FontSize',10,...
'String','Residues grid (optional)',...
'Style','text')

uicontrol('Parent',h1, 'Pos',[431 79 75 21],...
'Call',@compute_euler_uiCB,...
'String','3d netCDF',...
'Style','radiobutton',...
'Tooltip','Save residues grid as a 3D netCDF file',...
'Tag','radio_netcdf')

uicontrol('Parent',h1, 'Pos',[431 58 60 21],...
'Call',@compute_euler_uiCB,...
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
'Call',@compute_euler_uiCB,...
'String','STOP','Tag','push_stop');

uicontrol('Parent',h1,'Pos',[435 10 76 21],...
'Call',@compute_euler_uiCB,...
'String','Compute','Tag','push_compute');

function compute_euler_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

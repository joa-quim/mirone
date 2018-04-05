function varargout = TMD_OSU(varargin)
% Helper window to compute tides with the OSU grids from within Mirone.
%
% This is a rework of the Tidal Model Driver (TMD)
% COPYRIGHT: OREGON STATE UNIVERSITY, 2012
%
% NONCOMMERCIAL USE:
% OSU grants to you (hereafter, User) a royalty-free, nonexclusive right to execute, copy,
% modify and distribute both the binary and source code solely for academic, research and
% other similar noncommercial uses, subject to the following conditions:
% See entire license at http://volkov.oce.orst.edu/tides/COPYRIGHT.pdf

% --------------------------------------------------------------------
% Mironified by J. Luis
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id: TMD_OSU.m 10360 2018-04-05 15:44:12Z j $

	if (nargin > 1 && ischar(varargin{1}))
		gui_CB = str2func(varargin{1});
		[varargout{1:nargout}] = feval(gui_CB,varargin{2:end});
	else
		[varargout{1:nargout}] = TMD_OSU_OF(varargin{:});
	end

% ---------------------------------------------------------------------------------
function varargout = TMD_OSU_OF(varargin)

	hObject = figure('Vis','off');
	TMD_OSU_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right');

	handles.harmonic = ones(1, 14);
	handles.selected_comp = [{handles.check_M2} {1}];	% The by default selected component
	handles.comp_names = {'M2' 'S2' 'N2' 'K2' 'K1' 'O1' 'P1' 'Q1' 'MF' 'MM' 'M4' 'MS4', 'MN4' 'ALL'};

	% See if we have a TMP dir set by a ENV var that will take precedence over the default one
	path_tmp = '';
	t = set_gmt('MIRONE_TMP', 'whatever');				% Inquire if MIRONE_TMP exists
	if (~isempty(t))		% Now check that the dir exists
		if (exist(t, 'dir') == 7),		path_tmp = [t '/'];		end
	end
	if (isempty(path_tmp))		% If the above failed
		home_dir = fileparts(mfilename('fullpath'));	% Get the Mirone home dir and set path
		path_tmp = [home_dir '/../tmp/'];
	end

	try
		x = load([path_tmp 'used_pref.mat'], 'TMDmodel');
		if (isfield(x, 'TMDmodel')),	TMDmodel = x.TMDmodel;
		else,							TMDmodel = '';
		end
		[gfile, hfile, ufile, msg] = find_model_files(TMDmodel);
		if (~isempty(msg))		% A second chance to correct has been given inside find_model_files()
			errordlg({'And again I am warning you that the Model file has fake data.'; '';msg;''; 'Suiciding.'}, 'Error')
			delete(handles.figure1),		return
		end
	catch
		[gfile, hfile, ufile, msg, TMDmodel] = find_model_files([]);
		if (~isempty(msg))		% A second chance to correct has been given inside find_model_files()
			w = {'It looks that you do not have the model files. They are not shipped'
				'with Mirone because they are big, so you will have to get them first'
				'and come back here when you have them. Please see: '
				''
				'   http://volkov.oce.orst.edu/tides/global.html'
				''
				'and possibly the file you want is (attention: binary, NOT nc format)'
				''
				'   ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo9.tar.gz'};
			warndlg(w,'Warning','modal')
			delete(handles.figure1),		return
		end
	end
	if (isempty(gfile))
		delete(handles.figure1),		return	% Silently suiciding
	end

	handles.gfile = gfile;
	handles.hfile = hfile;
	handles.ufile = ufile;

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, [])
	%------------- END Pro look (3D) ------------------------------

	try
		% Plota a Bat
		[xy_lims,H,mz,iob,dt] = grd_in(handles.gfile);
		[n,m] = size(H);
		H = single(H)';
		H(H == 0) = NaN;
		dx = (xy_lims(2) - xy_lims(1)) / n;
		dy = (xy_lims(4) - xy_lims(3)) / m;
		x = xy_lims(1)+dx/2:dx:xy_lims(2)-dx/2;
		y = xy_lims(3)+dy/2:dy:xy_lims(4)-dy/2;
		handles.head = double([x(1) x(end) y(1) y(end) min(H(:)) max(H(:)) 0 dx dy]);
		handles.X = x;		handles.Y = y;		handles.Z = H;
		tmp.X = x;			tmp.Y = y;			tmp.head = handles.head;	tmp.name = 'Bathy for Maregs';
		pos = get(hObject, 'Pos');
		tmp.screenSize = get(0, 'ScreenSize');	tmp.screenSize(3) = tmp.screenSize(3) - pos(3) - 10;
		handles.hMirFig = mirone(H, tmp);
		move2side(hObject, handles.hMirFig,'left');
	catch
		errordlg(sprintf('Died because of this:\n\n%s', lasterr), 'Error')
		delete(handles.figure1),		return
	end

	handles.startTimeStr = datestr(now);
	set(handles.edit_startDate, 'Str', handles.startTimeStr)
	set(handles.edit_dataFiles, 'Str', TMDmodel)

	% Save Model's file name in used_pref.mat file for future lazziness.
	store_in_pref(handles, TMDmodel)

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

	if (nargin > 0),	external_drive(handles, 'TDM_OSU', varargin{:}),	end

% ----------------------------------------------------------------------------------------
function BDNfcn(hObject, evt, comp)
	handles = guidata(hObject);
	set(hObject, 'val', handles.harmonic(comp))	% Make sure the right-click did not unselect
	set(hObject, 'BackgroundColor', 'g')
	set(handles.selected_comp{1}, 'BackgroundColor', get(0,'factoryUicontrolBackgroundColor'))
	handles.selected_comp = [{hObject} {comp}];
	plot_fun(handles, comp)
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------------------
function worker(handles, hObject, comp)
% Do the tasks common to all components checkboxes
% HOBJECT is the handle of the clicked checkbox
	handles.harmonic(comp) = get(hObject, 'val');
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------------------
function plot_fun(handles, comp)
% ...
	if (get(handles.radio_plotBat, 'VAL'))
		return
	end

	[Z, name] = pick_component(handles, comp);
	Z = single(abs(Z))';
	setappdata(handles.hMirFig,'dem_z',Z);		% Update grid so that cursor display correct values
	img = scaleto8(Z);
	handMir = guidata(handles.hMirFig);
	set(handMir.hImg, 'CData', img)
	set(handles.hMirFig, 'Name', name)
	% Since u,v are on different grids we must set the coords too
	if (get(handles.radio_comp_u, 'Val') || get(handles.radio_comp_ut, 'Val'))
		X = handles.X - handles.head(8) / 2;	Y = handles.Y;
	elseif (get(handles.radio_comp_v, 'Val') || get(handles.radio_comp_vt, 'Val'))
		Y = handles.Y - handles.head(9) / 2;	X = handles.X;
	else
		X = handles.X;		Y = handles.Y;
	end
	setappdata(handles.hMirFig,'dem_x',X);
	setappdata(handles.hMirFig,'dem_y',Y);

% ----------------------------------------------------------------------------------------
function [Z, name] = pick_component(handles, comp)
% Pick (compute) the component that is select in GUI
	if (get(handles.radio_comp_z, 'Val'))
		if (comp > 13)			% Do them all
			Z = h_in(handles.hfile, 1);
			for (k = 2:13)
				Z = Z + h_in(handles.hfile, k);
			end
		else
			Z = h_in(handles.hfile, comp);		% reads in elevation for constituent COMP
		end
		Z(Z == 0) = NaN;		% A CONFIRMAR NECESSIDADE
		name = 'Elevation ';
		if (get(handles.radio_plotPhase, 'Val'))
			Z = atan2(-imag(Z), real(Z)) / pi*180;
			name = 'Phase ';
		end
	elseif (get(handles.radio_comp_u, 'Val') || get(handles.radio_comp_ut, 'Val'))
		if (comp > 13)			% Do them all
			Z = u_in(handles.hfile, 1, 'u');
			for (k = 2:13)
				Z = Z + u_in(handles.hfile, k, 'u');
			end
		else
			Z = u_in(handles.ufile, comp, 'u');
		end
		Z(Z == 0) = NaN;
		if (get(handles.radio_plotAmp, 'Val'))
			name = 'East transport ';
			if (get(handles.radio_comp_u, 'Val'))	% Velocity needs to be divided by
				hv = Huv(double(handles.Z'), 'u');
				Z = abs(Z) ./ hv * 100;
				name = 'East velocity ';
			end
		else
			Z = atan2(-imag(Z), real(Z)) / pi*180;
			name = 'Phase ';
		end
	elseif (get(handles.radio_comp_v, 'Val') || get(handles.radio_comp_vt, 'Val'))
		if (comp > 13)			% Do them all
			Z = u_in(handles.hfile, 1, 'v');
			for (k = 2:13)
				Z = Z + u_in(handles.hfile, k, 'v');
			end
		else
			Z = u_in(handles.ufile, comp, 'v');
		end
		Z(Z == 0) = NaN;
		if (get(handles.radio_plotAmp, 'Val'))
			name = 'North transport ';
			if (get(handles.radio_comp_v, 'Val'))	% Velocity needs to be divided by
				hu = Huv(double(handles.Z'), 'v');
				Z = abs(Z) ./ hu * 100;
				name = 'North velocity ';
			end
		else
			Z = atan2(-imag(Z), real(Z)) / pi*180;
			name = 'Phase ';
		end
	end
	name = [name handles.comp_names{comp}];

% ----------------------------------------------------------------------------------------
function check_M2_CB(hObject, handles)
	worker(handles, hObject, 1)

% ----------------------------------------------------------------------------------------
function check_S2_CB(hObject, handles)
	worker(handles, hObject, 2)

% ----------------------------------------------------------------------------------------
function check_N2_CB(hObject, handles)
	worker(handles, hObject, 3)

% ----------------------------------------------------------------------------------------
function check_K2_CB(hObject, handles)
	worker(handles, hObject, 4)

% ----------------------------------------------------------------------------------------
function check_K1_CB(hObject, handles)
	worker(handles, hObject, 5)

% ----------------------------------------------------------------------------------------
function check_O1_CB(hObject, handles)
	worker(handles, hObject, 6)

% ----------------------------------------------------------------------------------------
function check_P1_CB(hObject, handles)
	worker(handles, hObject, 7)

% ----------------------------------------------------------------------------------------
function check_Q1_CB(hObject, handles)
	worker(handles, hObject, 8)

% ----------------------------------------------------------------------------------------
function check_MF_CB(hObject, handles)
	worker(handles, hObject, 9)

% ----------------------------------------------------------------------------------------
function check_MM_CB(hObject, handles)
	worker(handles, hObject, 10)

% ----------------------------------------------------------------------------------------
function check_M4_CB(hObject, handles)
	worker(handles, hObject, 11)

% ----------------------------------------------------------------------------------------
function check_MS4_CB(hObject, handles)
	worker(handles, hObject, 12)

% ----------------------------------------------------------------------------------------
function check_MN4_CB(hObject, handles)
	worker(handles, hObject, 13)

% ----------------------------------------------------------------------------------------
function check_All_CB(hObject, handles)
	v = [handles.check_K2 handles.check_K1 handles.check_O1 handles.check_P1 handles.check_Q1 ...
		handles.check_MF handles.check_MM handles.check_M4 handles.check_MS4 handles.check_MN4];
	if (get(hObject,'Value'))
		set(v, 'Val', 1)
		handles.harmonic(4:13) = 1;
	else
		set(v, 'Val', 0)
		handles.harmonic(4:13) = 0;
	end
	set([handles.check_M2 handles.check_S2 handles.check_M2], 'Val', 1)		% These are always set on
	handles.harmonic(1:3) = 1;
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------------------
function radio_plotBat_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_plotAmp handles.radio_plotPhase],'Val',0)
	handMir = guidata(handles.hMirFig);
	set(handMir.hImg, 'CData', scaleto8(handles.Z))
	set(handles.hMirFig, 'Name', 'Bathymetry')

% ----------------------------------------------------------------------------------------
function radio_plotAmp_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_plotBat handles.radio_plotPhase],'Val',0)

% ----------------------------------------------------------------------------------------
function radio_plotPhase_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_plotAmp handles.radio_plotBat],'Val',0)

% ----------------------------------------------------------------------------------------
function radio_comp_z_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_comp_u handles.radio_comp_v handles.radio_comp_ut handles.radio_comp_vt],'Val',0)
	plot_fun(handles, handles.selected_comp{2})

% ----------------------------------------------------------------------------------------
function radio_comp_u_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_comp_z handles.radio_comp_v handles.radio_comp_ut handles.radio_comp_vt],'Val',0)
	plot_fun(handles, handles.selected_comp{2})

% ----------------------------------------------------------------------------------------
function radio_comp_v_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_comp_z handles.radio_comp_u handles.radio_comp_vt handles.radio_comp_ut],'Val',0)
	plot_fun(handles, handles.selected_comp{2})

% ----------------------------------------------------------------------------------------
function radio_comp_vt_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_comp_z handles.radio_comp_u handles.radio_comp_v handles.radio_comp_ut],'Val',0)
	plot_fun(handles, handles.selected_comp{2})

% ----------------------------------------------------------------------------------------
function radio_comp_ut_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_comp_z handles.radio_comp_u handles.radio_comp_v handles.radio_comp_vt],'Val',0)
	plot_fun(handles, handles.selected_comp{2})

% ----------------------------------------------------------------------------------------
function push_dataFiles_CB(hObject, handles)
	edit_dataFiles_CB(handles.edit_dataFiles, handles)	% Let it do the work

% ----------------------------------------------------------------------------------------
function edit_dataFiles_CB(hObject, handles)
% ...
	[gfile, hfile, ufile, msg, TMDmodel] = find_model_files([]);
	if (isempty(gfile)),	return,		end		% User gave up
	if (~isempty(msg)),		errordlg(msg, 'Error'),		end
	set(handles.edit_dataFiles, 'Str', TMDmodel)
	store_in_pref(handles, TMDmodel)

% ----------------------------------------------------------------------------------------
function edit_startDate_CB(hObject, handles)
	if (~strcmp(get(hObject, 'Str'), handles.startTimeStr))
		str = uisetdate;
		if (isempty(str)),		return,		end
		set(handles.edit_startDate, 'Str', str)
		handles.startTimeStr = str;
		guidata(handles.figure1, handles)
	end

% ----------------------------------------------------------------------------------------
function push_startDate_CB(hObject, handles)
	str = uisetdate;
	if (isempty(str)),		return,		end
	set(handles.edit_startDate, 'Str', str)
	handles.startTimeStr = str;
	guidata(handles.figure1, handles)

% ----------------------------------------------------------------------------------------
function edit_duration_CB(hObject, handles)
	xx = str2double(get(hObject, 'Str'));
	if (isnan(xx) || xx <= 0)
		set(hObject, 'Str', '48')
	end

% ----------------------------------------------------------------------------------------
function edit_step_CB(hObject, handles)
	xx = str2double(get(hObject, 'Str'));
	if (isnan(xx) || xx <= 0)
		set(hObject, 'Str', '60')
	end

% ----------------------------------------------------------------------------------------
function push_plotTidePt_CB(hObject, handles)
% Call the Draw_CB() Mirone function to plot one point (Tide gauge)
	if (~ishandle(handles.hMirFig))
		errordlg('You killed it (the Mir window). Bye Bye','Error')
		delete(handles.figure1),	return
	end
	handMir = guidata(handles.hMirFig);
	figure(handles.hMirFig)
	h = mirone('Draw_CB', handMir, 'Symbol', 'o');
	set(h, 'Tag', 'TMDpt')

% ----------------------------------------------------------------------------------------
function push_loadTidesFile_CB(hObject, handles)
% Load a file with maregraph positions
	if (~ishandle(handles.hMirFig))
		errordlg('You killed it (the Mir window). Bye Bye','Error')
		delete(handles.figure1),	return
	end
	handMir = guidata(handles.hMirFig);
	figure(handles.hMirFig)
	load_xyz(handMir, [], 'AsPoint');
	h = findobj(handMir.axes1, 'Tag','Pointpolyline');
	set(h, 'Tag', 'TMDpt')

% ----------------------------------------------------------------------------------------
function push_plotTide_CB(hObject, handles)
% ...
	[SerialDay, TimeSeries] = push_compute(handles);
	hf = ecran(SerialDay, TimeSeries);
	hAxes = findobj(hf, 'Tag', 'axes1');
	setappdata(hAxes, 'LabelFormatType', 'Date');		% Tell pixval_stsbar to display XX coords as Date-Time
	h = findobj(hf,'Tag','add_uictx');
	cb = get(h, 'Callback');
	feval(cb, h, guidata(hf))			% Call the ecran's add_uictx_CB function

% ----------------------------------------------------------------------------------------
function push_saveTide_CB(hObject, handles)
% Compute and save the Tide series in a file

	handMir = guidata(handles.hMirFig);
	h = findobj(handMir.axes1, 'Tag', 'TMDpt');
	if (isempty(h))
		warndlg('Select point(s) where to compute tides first. Fair isn''t it?', 'Warning')
		return
	else
		lon = get(h, 'XData');		lat = get(h, 'YData');
		if (isa(lon, 'cell'))
			lon = cell2mat(lon);	lat = cell2mat(lat);
		end
	end

	[SerialDay, TimeSeries] = push_compute(handles);
	
	str1 = {'*.dat;*.DAT', 'Tide file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select Tide File name','put', name);
	if isequal(FileName,0),			return,		end

	fid = fopen([PathName FileName],'w');
	fprintf(fid,'# DATENUM dd-mmm-yyyy HH:MM:SS\n');
	fprintf(fid,'# Constituents included: ');
	fprintf(fid,'%s  ', handles.comp_names{handles.harmonic});
	fprintf(fid,'\n');

	oname = 'z(m)';
	if (get(handles.radio_comp_ut, 'Val') || get(handles.radio_comp_vt, 'Val'))
		oname = 'U (m^2/s)';
		if (get(handles.radio_comp_vt, 'Val')),		oname(1) = 'V';		end
	elseif (get(handles.radio_comp_u, 'Val') || get(handles.radio_comp_v, 'Val'))
		oname = 'u (cm/s)';
		if (get(handles.radio_comp_vt, 'Val')),		oname(1) = 'v';		end
	end
	fprintf(fid,'# Parameter: %s\n', oname);

	frmt = [repmat('%12.5f\t', 1, numel(lon)-1) '%12.5f\n'];
	fprintf(fid, ['# Lon:\t' frmt], lon');
	fprintf(fid, ['# Lat:\t' frmt], lat');

	frmt = ['%s ' repmat('%10.4f\t', 1, size(TimeSeries,2)-1) '%10.4f\n'];
	for (k = 1:size(SerialDay,1))
		fprintf(fid, frmt, datestr(SerialDay(k,1)), TimeSeries(k,:));
	end
	fclose(fid);

% ----------------------------------------------------------------------------------------
function [SerialDay, TimeSeries] = push_compute(handles)
% ...
	handMir = guidata(handles.hMirFig);
	h = findobj(handMir.axes1, 'Tag', 'TMDpt');
	if (isempty(h))
		lon = 347;	lat = 36;		% Por enquanto nao Erra
	else
		lon = get(h, 'XData');		lat = get(h, 'YData');
		if (isa(lon, 'cell'))
			lon = cell2mat(lon);	lat = cell2mat(lat);
		end
	end
	n_pts = numel(lon);
	nz = min(sum(handles.harmonic), 13);
	amp = zeros(nz, n_pts);		pha = amp;

	conList = rd_con(handles.hfile);
	nz = 0;
	nc = min(numel(handles.harmonic), 13);
	Cid = repmat('    ', min(sum(handles.harmonic),13), 1);
	for (k = 1:nc)		% Constroi um vector com os nomes das componentes
		if (handles.harmonic(k))
			nz = nz + 1;
			Cid(nz,:) = conList(k,:);
		end
	end

	ic2 = 1;
	for (k = 1:13)
		if (~handles.harmonic(k)),	continue,	end
		Z = pick_component(handles, k);
		if (get(handles.radio_comp_z, 'Val'))
			zi_r = grdtrack_m(single(real(Z))', handles.head, [lon(:) lat(:)], '-Z');
			zi_i = grdtrack_m(single(imag(Z))', handles.head, [lon(:) lat(:)], '-Z');
			zi = zi_r + 1i*zi_i;
		elseif (get(handles.radio_comp_u, 'Val'))
		end
		amp(ic2,:) = abs(zi)';
		pha(ic2,:) = atan2(-imag(zi), real(zi))';
		ic2 = ic2 + 1;
	end

	[yy, mm, dd, hh, mi, ss] = datevec(get(handles.edit_startDate, 'Str'));

	duration = str2double(get(handles.edit_duration, 'Str')) / 24;		% In decimal days
	dt = str2double(get(handles.edit_step, 'Str')) / (24 * 60);			% Step decimal days
	d1 = datenum(yy, mm, dd, hh, 0, 0);									% Start at round hours
	SerialDay = (d1:dt:d1+duration)';
	
	d0 = datenum(1992,1,1);		% corresponds to 48622mjd
	TimeSeries = zeros(numel(SerialDay), n_pts);
	for (k = 1:n_pts)
		cph = -1i*pha(:,k);% * pi/180;		Already in rads
		cam = amp(:,k) .* exp(cph);
		dh   = InferMinor(cam, conList, SerialDay);
		hhat = harp1(SerialDay-d0, cam, Cid);
		TimeSeries(:, k) = hhat + dh;
	end

% ----------------------------------------------------------------------------------------
function dh = InferMinor(zmaj,cid,SDtime)
% Lana Erofeeva, re-make for matlab OCT 2004
% usage:
% [dh]=InferMinor(zmaj,cid,SDtime);
%
% Based on Richard Ray's code perth2
% Return correction for the 16 minor constituents
% zeros, if not enough input constituents for inference
% Input:
% cid(ncon,4)   - GIVEN constituents
% zmaj(ncon,... - Complex HC for GIVEN constituents/points
% SDtime(...    - time expressed in Serial Days (see help datenum)
% Modes:
%      Time series: zmaj(ncon,1),  SDtime(nt,1)  ->Output: dh(nt,1)
%      Drift Track: zmaj(ncon,nt), SDtime(nt,1)  ->Output: dh(nt,1)
%      Map:         zmaj(ncon,N,M),SDtime(1,1)   ->Output: dh(N,M)

rad=pi/180;		PP=282.8;	dh=NaN;
% check size
[ncon,dum] = size(cid);
if (dum ~= 4)
	fprintf('InferMinor call: Wrong constituents\n')
	return
end
[n1,n2,n3] = size(zmaj);
if (n1~=ncon)
	zmaj=conj(zmaj');
	[n1,n2,n3]=size(zmaj);
end
if (n1~=ncon)
	fprintf('InferMinor call: Wrong zmaj size\n')
	return
end
nt=length(SDtime);
%[n1 n2 n3 nt]
if n2==1 && n3==1
	N=nt;M=1;
elseif n3==1 && n2>1 && nt==n2
	N=nt;M=1;
elseif n2>1 && nt==1
	N=n2;M=n3;
else
	fprintf('InferMinor call: Wrong zmaj/SDtime size\n');
	return
end
%
dh=zeros(N,M);
d0=datenum(1992,1,1); % corresponds to 48622mjd
d1=SDtime;
time=d1-d0; % time (days) relatively Jan 1 1992 GMT
time_mjd=48622+time;
[ncon,dum]=size(cid);
cid8=['q1  ';'o1  ';'p1  ';'k1  ';'n2  ';'m2  ';'s2  ';'k2  '];
% re-order zmaj to correspond to cid8
z8=zeros(8,N,M);
ni=0;
for i1=1:8
	for j1=1:ncon
		if cid(j1,:)==cid8(i1,:)
			z8(i1,:,:)=zmaj(j1,:,:);
			if (i1~=3 && i1~=8), ni=ni+1; end
		end
	end
end
if (ni < 6)				 % not enough constituents for inference
	return
end
%
zmin=zeros(18,N,M);   
zmin(1,:,:)  = 0.263 *z8(1,:,:) - 0.0252*z8(2,:,:);   %2Q1
zmin(2,:,:)  = 0.297 *z8(1,:,:) - 0.0264*z8(2,:,:);   %sigma1
zmin(3,:,:)  = 0.164 *z8(1,:,:) + 0.0048*z8(2,:,:);   %rho1 +
zmin(4,:,:)  = 0.0140*z8(2,:,:) + 0.0101*z8(4,:,:);   %M1
zmin(5,:,:)  = 0.0389*z8(2,:,:) + 0.0282*z8(4,:,:);   %M1
zmin(6,:,:)  = 0.0064*z8(2,:,:) + 0.0060*z8(4,:,:);   %chi1
zmin(7,:,:)  = 0.0030*z8(2,:,:) + 0.0171*z8(4,:,:);   %pi1
zmin(8,:,:)  =-0.0015*z8(2,:,:) + 0.0152*z8(4,:,:);   %phi1
zmin(9,:,:)  =-0.0065*z8(2,:,:) + 0.0155*z8(4,:,:);   %theta1
zmin(10,:,:) =-0.0389*z8(2,:,:) + 0.0836*z8(4,:,:);   %J1 +
zmin(11,:,:) =-0.0431*z8(2,:,:) + 0.0613*z8(4,:,:);   %OO1 +
zmin(12,:,:) = 0.264 *z8(5,:,:) - 0.0253*z8(6,:,:);   %2N2 +
zmin(13,:,:) = 0.298 *z8(5,:,:) - 0.0264*z8(6,:,:);   %mu2 +
zmin(14,:,:) = 0.165 *z8(5,:,:) + 0.00487*z8(6,:,:);  %nu2 +
zmin(15,:,:) = 0.0040*z8(6,:,:) + 0.0074*z8(7,:,:);   %lambda2
zmin(16,:,:) = 0.0131*z8(6,:,:) + 0.0326*z8(7,:,:);   %L2 +
zmin(17,:,:) = 0.0033*z8(6,:,:) + 0.0082*z8(7,:,:);   %L2 +
zmin(18,:,:) = 0.0585*z8(7,:,:);                      %t2 + 
%
hour = (time - floor(time))*24.D0;
t1 = 15.*hour;
t2 = 30.*hour;
[S,H,P,omega]=astrol(time_mjd); % nsites x nt
%
arg=zeros(18,N,M);
arg(1,:,:) = t1 - 4.*S + H + 2.*P - 90.;     % 2Q1
arg(2,:,:) = t1 - 4.*S + 3.*H - 90.;         % sigma1
arg(3,:,:) = t1 - 3.*S + 3.*H - P - 90.;     % rho1
arg(4,:,:) = t1 - S + H - P + 90.;           % M1
arg(5,:,:) = t1 - S + H + P + 90.;           % M1
arg(6,:,:) = t1 - S + 3.*H - P + 90.;        % chi1
arg(7,:,:) = t1 - 2.*H + PP - 90.;           % pi1
arg(8,:,:) = t1 + 3.*H + 90.;                % phi1
arg(9,:,:) = t1 + S - H + P + 90.;           % theta1
arg(10,:,:) = t1 + S + H - P + 90.;          % J1
arg(11,:,:) = t1 + 2.*S + H + 90.;           % OO1
arg(12,:,:) = t2 - 4.*S + 2.*H + 2.*P;       % 2N2
arg(13,:,:) = t2 - 4.*S + 4.*H;              % mu2
arg(14,:,:) = t2 - 3.*S + 4.*H - P;          % nu2
arg(15,:,:) = t2 - S + P + 180.D0;           % lambda2
arg(16,:,:) = t2 - S + 2.*H - P + 180.D0;    % L2
arg(17,:,:) = t2 - S + 2.*H + P;             % L2
arg(18,:,:) = t2 - H + PP;                   % t2
%
%     determine nodal corrections f and u
sinn = sin(omega*rad);
cosn = cos(omega*rad);
sin2n = sin(2.*omega*rad);
cos2n = cos(2.*omega*rad);
f = ones(18,N,M);
f(1,:,:) = sqrt((1.0 + 0.189*cosn - 0.0058*cos2n).^2 + (0.189*sinn - 0.0058*sin2n).^2);
f(2,:,:) = f(1,:,:);
f(3,:,:) = f(1);
f(4,:,:) = sqrt((1.0 + 0.185*cosn).^2 + (0.185*sinn).^2);
f(5,:,:) = sqrt((1.0 + 0.201*cosn).^2 + (0.201*sinn).^2);
f(6,:,:) = sqrt((1.0 + 0.221*cosn).^2 + (0.221*sinn).^2);
f(10,:,:) = sqrt((1.0 + 0.198*cosn).^2 + (0.198*sinn).^2);
f(11,:,:) = sqrt((1.0 + 0.640*cosn + 0.134*cos2n).^2 + (0.640*sinn + 0.134*sin2n).^2 );
f(12,:,:) = sqrt((1.0 - 0.0373*cosn).^2 + (0.0373*sinn).^2);
f(13,:,:) = f(12,:,:);
f(14,:,:) = f(12,:,:);
f(16,:,:) = f(12,:,:);
f(17,:,:) = sqrt((1.0 + 0.441*cosn).^2 + (0.441*sinn).^2);
%
u = zeros(18,N,M);
u(1,:,:) = atan2(0.189*sinn - 0.0058*sin2n, 1.0 + 0.189*cosn - 0.0058*sin2n)/rad;
u(2,:,:) = u(1,:,:);
u(3,:,:) = u(1,:,:);
u(4,:,:) = atan2( 0.185*sinn, 1.0 + 0.185*cosn)/rad;
u(5,:,:) = atan2(-0.201*sinn, 1.0 + 0.201*cosn)/rad;
u(6,:,:) = atan2(-0.221*sinn, 1.0 + 0.221*cosn)/rad;
u(10,:,:) = atan2(-0.198*sinn, 1.0 + 0.198*cosn)/rad;
u(11,:,:) = atan2(-0.640*sinn - 0.134*sin2n, 1.0 + 0.640*cosn + 0.134*cos2n)/rad;
u(12,:,:) = atan2(-0.0373*sinn, 1.0 - 0.0373*cosn)/rad;
u(13,:,:) = u(12,:,:);
u(14,:,:) = u(12,:,:);
u(16,:,:) = u(12,:,:);
u(17,:,:) = atan2(-0.441*sinn, 1.0 + 0.441*cosn)/rad;
%     sum over all tides
%     ------------------
for k=1:18
  tmp=squeeze(real(zmin(k,:,:)).*f(k,:,:).*...
                      cos((arg(k,:,:) + u(k,:,:))*rad)-...
                      imag(zmin(k,:,:)).*f(k,:,:).*...
                      sin((arg(k,:,:)+u(k,:,:))*rad));
  [n1,n2]=size(tmp);if n1==1,tmp=tmp';end
  dh(:,:) = dh(:,:) + tmp;
end

% ----------------------------------------------------------------------------------------
function  [s,h,p,N] = astrol(time)
%  Computes the basic astronomical mean longitudes  s, h, p, N.
%  Note N is not N', i.e. N is decreasing with time.
%  These formulae are for the period 1990 - 2010, and were derived
%  by David Cartwright (personal comm., Nov. 1990).
%  time is UTC in decimal MJD.
%  All longitudes returned in degrees.
%  R. D. Ray    Dec. 1990
%  Non-vectorized version. Re-make for matlab by Lana Erofeeva, 2003
% usage: [s,h,p,N]=astrol(time)
%        time, MJD
	circle=360;
	T = time - 51544.4993;
	% mean longitude of moon
	% ----------------------
	s = 218.3164 + 13.17639648 * T;
	% mean longitude of sun
	% ---------------------
	h = 280.4661 +  0.98564736 * T;
	% mean longitude of lunar perigee
	% -------------------------------
	p =  83.3535 +  0.11140353 * T;
	% mean longitude of ascending lunar node
	% --------------------------------------
	N = 125.0445D0 -  0.05295377D0 * T;
	%
	s = mod(s,circle);
	h = mod(h,circle);
	p = mod(p,circle);
	N = mod(N,circle);
	%
	if s<0, s = s + circle; end
	if h<0, h = h + circle; end
	if p<0, p = p + circle; end
	if N<0, N = N + circle; end

% ----------------------------------------------------------------------------------------
function [gfile, hfile, ufile, msg, TMDmodel] = find_model_files(TMDmodel, count)
% If called with COUNT = 0 no second chance via a call to uigetfile is provided.
	if (nargin == 1),	count = 1;	end			% To not let recursivity calls go beyond level 2
	msg = '';	gfile = '';		hfile = '';		ufile = '';
	if (~isempty(TMDmodel))
		pato = fileparts(TMDmodel);
		fid  = fopen(TMDmodel);
		if (fid < 1)
			msg = ['Cannot open the Master model file (' TMDmodel ')'];		return
		end
		hfile = fgetl(fid);		ufile = fgetl(fid);		gfile = fgetl(fid);
		fclose(fid);
		[p, fn, ext] = fileparts(hfile);	hfile = [pato '/' fn ext];
		[p, fn, ext] = fileparts(ufile);	ufile = [pato '/' fn ext];
		[p, fn, ext] = fileparts(gfile);	gfile = [pato '/' fn ext];
		if (~exist(hfile,'file') || ~exist(ufile,'file') || ~exist(gfile,'file'))
			msg = ['One or more of the model files listed in the Master model file (' TMDmodel ') does not exist'];
		end
	end
	if ((isempty(TMDmodel) || ~isempty(msg)) && count == 1)		% This branch is called directly when TMDmodel = ''
		if (~isempty(msg)),	warndlg(msg, 'Warning'),	end
		[fname, pato] = uigetfile('Model*', 'Open MODEL file');
		if (fname ~= 0)
			TMDmodel = [pato fname];
			[gfile, hfile, ufile, msg] = find_model_files(TMDmodel, 2);		% Try again
		end
	end

% ----------------------------------------------------------------------------------------
function store_in_pref(handles, TMDmodel)
% Save Model file name in used_pref.mat file for future lazziness.
	handMir = guidata(handles.hMirFig);
	V7 = handMir.version7;
	try			% Try if file already exist and if yes, append
		x.x = load([handMir.path_tmp 'used_pref.mat']);
		if (~V7),		save([handMir.path_tmp 'used_pref.mat'],'TMDmodel','-append')		% R <= 13
		else,			save([handMir.path_tmp 'used_pref.mat'],'TMDmodel','-append', '-v6')
		end
	catch		% First time. Create a new one
		if (~V7),		save([handMir.path_tmp 'used_pref.mat'],'TMDmodel','V7')
		else,			save([handMir.path_tmp 'used_pref.mat'],'TMDmodel','V7','-v6')
		end
	end

% ----------------------------------------------------------------------------------------
function conList = rd_con(fname)
% read constituents from h*out or u*out file
	fid = fopen(fname,'r','b');
	fread(fid,1,'long');
	nm = fread(fid,3,'long');
	nc = nm(3);
	fread(fid,2,'float');		% th_lim
	fread(fid,2,'float');		% ph_lim
	C = fread(fid,nc*4,'uchar');
	C = reshape(C,4,nc);	C=C';
	conList = char(C);
	fclose(fid);

% ----------------------------------------------------------------------------------------
function hhat = harp(time, hc, con)
% function to predict tidal elevation fields at "time" using harmonic constants
%
% INPUT: time (days) relatively Jan 1, 1992 (48622mjd)
%        con(nc,4) - char*4 tidal constituent IDs 
%        hc(n,m,nc) - harmonic constant vector  (complex)
% OUTPUT:hhat - tidal field at time
%  
%        Nodal corrections included
%
% usage: [hhat]=harp(time,hc,con);
%
% see harp1 for making time series at ONE location

	[n, m, nc] = size(hc);	hhat = zeros(n, m);
	nc1 = size(con, 1);
	if (nc1 ~= nc)
		fprintf('Imcompatible hc and con: check dimensions\n');
		return
	end
	if (length(time) > 1),	time = time(1);		end
	ispec = zeros(1,nc);	amp = zeros(1,nc);	ph = zeros(1,nc);	omega = zeros(1,nc);	alpha = zeros(1,nc);
	for (k = 1:nc)
		[ispec(k),amp(k),ph(k),omega(k),alpha(k)] = constit(con(k,:));
	end
	%
	igood = find(ispec ~= -1);
	con1  = con(igood,:);
	%
	[pu1,pf1] = nodal(time+48622,con1);
	pu = zeros(1,nc);		pf = ones(1,nc);
	pu(igood) = pu1;		pf(igood) = pf1;
	%
	x = zeros(n,m,2*nc);
	x(:,:,1:2:end) = real(hc);
	x(:,:,2:2:end) = imag(hc);
	for k = 1:nc
		arg = pf(k)*x(:,:,2*k-1)*cos(omega(k)*time*86400+ph(k)+pu(k))...
			- pf(k)*x(:,:,2*k)  *sin(omega(k)*time*86400+ph(k)+pu(k));
		hhat = hhat + arg;
	end

% ----------------------------------------------------------------------------------------
function hhat = harp1(time,hc,con)
% function to predict tidal time series
% using harmonic constants
% INPUT: time (days) relatively Jan 1, 1992 (48622mjd)
%        con(nc,4) - char*4 tidal constituent IDs 
%        hc(nc) - harmonic constant vector  (complex)
% OUTPUT:hhat - time series reconstructed using HC
%  
%        Nodal corrections included
%
% usage: [hhat]=harp1(time,hc,con);
%
	L = length(time);
	[n1,n2] = size(time);
	if (n1 == 1),	time = time';	end
	[nc,dum] = size(con);

	ispec = zeros(1,nc);	amp = zeros(1,nc);	ph = zeros(1,nc);	omega = zeros(1,nc);	alpha = zeros(1,nc);
	for (k = 1:nc)
		[ispec(k),amp(k),ph(k),omega(k),alpha(k),cNum]=constit(con(k,:));
	end
	%
	igood = find(ispec ~= -1);
	con1  = con(igood,:);
	%
	[pu1,pf1] = nodal(time+48622,con1);
	pu = zeros(L,nc);		pf=ones(L,nc);
	pu(:,igood) = pu1;pf(:,igood) = pf1;
	%
	hhat = zeros(size(time));
	x = zeros(2*nc,1);
	x(1:2:end)=real(hc);
	x(2:2:end)=imag(hc);
	for k=1:nc
		arg = pf(:,k).*x(2*k-1).*cos(omega(k)*time*86400+ph(k)+pu(:,k))...
			- pf(:,k).*x(2*k)  .*sin(omega(k)*time*86400+ph(k)+pu(:,k));
		hhat=hhat+arg;
	end

% ----------------------------------------------------------------------------------------
function [ispec,amp,ph,omega,alpha,constitNum] = constit(c)
%  constit returns amplitude, phase, frequency,alpha, species for
%  a tidal constituent identified by a 4 character string
%   This version returns zeros for all phases
%  Usage: [ispec,amp,ph,omega,alpha,constitNum] = constit(c);
	amp=0;	ph=0;	omega=0;	alpha=0;	constitNum=0;
	%  c : character string ID of constituents that the program knows about
	%   all other variables are in the same order

	c_data = [	'm2  ';'s2  ';'k1  ';'o1  '; ...
				'n2  ';'p1  ';'k2  ';'q1  '; ...
				'2n2 ';'mu2 ';'nu2 ';'l2  '; ...
				't2  ';'j1  ';'no1 ';'oo1 '; ...
				'rho1';'mf  ';'mm  ';'ssa ';'m4  ';'ms4 ';'mn4 ' ];

	constitNum_data = [1,2,5,6,3,7,4,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
	%   ispec : species type (spherical harmonic dependence of quadropole potential)
	ispec_data =   [2;2;1;1; ...
					2;1;2;1; ...
					2;2;2;2; ...
					2;1;1;1; ...
					1;0;0;0;0;0;0];
	%   alpha : loading love number ... frequncy dependance here is suspect
	alpha_data = [0.693;0.693;0.736;0.695; ...
				0.693;0.706;0.693;0.695; ...
				0.693;0.693;0.693;0.693; ...
				0.693;0.693;0.693;0.693; ...
				0.693;0.693;0.693;0.693;0.693;0.693;0.693];

	%  omega : frequencies
	omega_data = [1.405189e-04;1.454441e-04;7.292117e-05;6.759774e-05; ...
				1.378797e-04;7.252295e-05;1.458423e-04;6.495854e-05; ...
				1.352405e-04;1.355937e-04;1.382329e-04;1.431581e-04; ...
				1.452450e-04;7.556036e-05;7.028195e-05;7.824458e-05; ...
				6.531174e-05;0.053234e-04;0.026392e-04;0.003982e-04; ...
				2.810377e-04;2.859630e-04;2.783984e-04];

	%  phase : Astronomical arguments (relative to t0 = 1 Jan 0:00 1992]
	% Richrad Ray subs: "arguments" and "astrol"
	phase_data = [ 1.731557546;0.000000000;0.173003674;1.558553872;...
				6.050721243;6.110181633;3.487600001;5.877717569;
				4.086699633;3.463115091;5.427136701;0.553986502;
				0.052841931;2.137025284;2.436575100;1.929046130;
				5.254133027;1.756042456;1.964021610;3.487600001;
				3.463115091;1.731557546;1.499093481];
	%  amp   :   amplitudes
	%amp_data = [0.242334;0.112743;0.141565;0.100661; ...
	amp_data = [0.2441  ;0.112743;0.141565;0.100661; ...
				0.046397;0.046848;0.030684;0.019273; ...
				0.006141;0.007408;0.008811;0.006931; ...
				0.006608;0.007915;0.007915;0.004338; ...
				0.003661;0.042041;0.022191;0.019567;0;0;0];

	nspecies = length(ispec_data);
	nl = length(c);
	kk = 0;
	for k = 1:nspecies
		if (lower(c) == c_data(k,1:nl))
			kk = k;
			break
		end
	end
	if (kk == 0)
		ispec = -1;
		return
	else
		constitNum = constitNum_data(kk);
		ispec = ispec_data(kk);
		amp   = amp_data(kk);
		alpha = alpha_data(kk);
		omega = omega_data(kk);
		ph = phase_data(kk);
	end

% ----------------------------------------------------------------------------------------
% ARGUMENTS and ASTROL FORTRAN subroutines SUPPLIED by RICHARD RAY, March 1999
% This is matlab remake of ARGUMENTS by Lana Erofeeva, Jan 2003
% NOTE - "no1" in constit.h corresponds to "M1" in arguments
% usage: [pu,pf]=nodal(time,cid);
% time - mjd, cid(nc,4) - tidal constituents array char*4
% pu(:,nc),pf(:,nc) - nodal corrections for the constituents
% 
function [pu,pf] = nodal(time,cid)
%
	cid0 = ['m2  ';'s2  ';'k1  ';'o1  '; ...
			'n2  ';'p1  ';'k2  ';'q1  '; ...
			'2n2 ';'mu2 ';'nu2 ';'l2  '; ...
			't2  ';'j1  ';'no1 ';'oo1 '; ...
			'rho1';'mf  ';'mm  ';'ssa ';'m4  ';...
			'ms4 ';'mn4 '];
	index = [30,35,19,12,27,17,37,10,25,26,28,33,34,23,14,24,11,5,3,2,45,46,44];
	%     Determine equilibrium arguments
	%     -------------------------------
% 	pp  = 282.94; % solar perigee at epoch 2000
	rad = pi/180;
	n1 = size(time, 1);
	if (n1 == 1),	time = time';	end
	[s,h,p,omega] = astrol(time);
	hour = (time - floor(time))*24.;
	t1   = 15*hour;		%t2 = 30.*hour;
	nT   = length(t1);
% 	arg  = zeros(nT,53);
% 	% arg is not needed for now, but let it be here
% 	% arg is kept in constit.m (phase_data) for time= Jan 1, 1992, 00:00 GMT
% 	arg(:, 1) = h - pp;                        % Sa
% 	arg(:, 2) = 2*h;                           % Ssa
% 	arg(:, 3) = s - p;                         % Mm
% 	arg(:, 4) = 2*s - 2*h;                     % MSf
% 	arg(:, 5) = 2*s;                           % Mf
% 	arg(:, 6) = 3*s - p;                       % Mt
% 	arg(:, 7) = t1 - 5*s + 3*h + p - 90;       % alpha1
% 	arg(:, 8) = t1 - 4*s + h + 2*p - 90;       % 2Q1
% 	arg(:, 9) = t1 - 4*s + 3*h - 90;           % sigma1
% 	arg(:,10) = t1 - 3*s + h + p - 90;         % q1
% 	arg(:,11) = t1 - 3*s + 3*h - p - 90;       % rho1
% 	arg(:,12) = t1 - 2*s + h - 90;             % o1
% 	arg(:,13) = t1 - 2*s + 3*h + 90;           % tau1
% 	arg(:,14) = t1 - s + h + 90;               % M1
% 	arg(:,15) = t1 - s + 3*h - p + 90;         % chi1
% 	arg(:,16) = t1 - 2*h + pp - 90;            % pi1
% 	arg(:,17) = t1 - h - 90;                   % p1
% 	arg(:,18) = t1 + 90;                       % s1
% 	arg(:,19) = t1 + h + 90;                   % k1
% 	arg(:,20) = t1 + 2*h - pp + 90;            % psi1
% 	arg(:,21) = t1 + 3*h + 90;                 % phi1
% 	arg(:,22) = t1 + s - h + p + 90;           % theta1
% 	arg(:,23) = t1 + s + h - p + 90;           % J1
% 	arg(:,24) = t1 + 2*s + h + 90;             % OO1
% 	arg(:,25) = t2 - 4*s + 2*h + 2*p;          % 2N2
% 	arg(:,26) = t2 - 4*s + 4*h;                % mu2
% 	arg(:,27) = t2 - 3*s + 2*h + p;            % n2
% 	arg(:,28) = t2 - 3*s + 4*h - p;            % nu2
% 	arg(:,29) = t2 - 2*s + h + pp;             % M2a
% 	arg(:,30) = t2 - 2*s + 2*h;                % M2
% 	arg(:,31) = t2 - 2*s + 3*h - pp;           % M2b
% 	arg(:,32) = t2 - s + p + 180.;             % lambda2
% 	arg(:,33) = t2 - s + 2*h - p + 180.;       % L2
% 	arg(:,34) = t2 - h + pp;                   % T2
% 	arg(:,35) = t2;                            % S2
% 	arg(:,36) = t2 + h - pp + 180;             % R2
% 	arg(:,37) = t2 + 2*h;                      % K2
% 	arg(:,38) = t2 + s + 2*h - pp;             % eta2
% 	arg(:,39) = t2 - 5*s + 4.0*h + p;          % MNS2
% 	arg(:,40) = t2 + 2*s - 2*h;                % 2SM2
% 	arg(:,41) = 1.5*arg(:,30);                 % M3
% 	arg(:,42) = arg(:,19) + arg(:,30);         % MK3
% 	arg(:,43) = 3*t1;                          % S3
% 	arg(:,44) = arg(:,27) + arg(:,30);         % MN4
% 	arg(:,45) = 2*arg(:,30);                   % M4
% 	arg(:,46) = arg(:,30) + arg(:,35);         % MS4
% 	arg(:,47) = arg(:,30) + arg(:,37);         % MK4
% 	arg(:,48) = 4*t1;                          % S4
% 	arg(:,49) = 5*t1;                          % S5
% 	arg(:,50) = 3*arg(:,30);                   % M6
% 	arg(:,51) = 3*t2;                          % S6
% 	arg(:,52) = 7.0*t1;                        % S7
% 	arg(:,53) = 4*t2;                          % S8
%
	%     determine nodal corrections f and u 
	%     -----------------------------------
	sinn = sin(omega*rad);
	cosn = cos(omega*rad);
	sin2n = sin(2*omega*rad);
	cos2n = cos(2*omega*rad);
	sin3n = sin(3*omega*rad);

	f = zeros(nT,53);
	f(:,1) = 1;                                     % Sa
	f(:,2) = 1;                                     % Ssa
	f(:,3) = 1 - 0.130*cosn;                        % Mm
	f(:,4) = 1;                                     % MSf
	f(:,5) = 1.043 + 0.414*cosn;                    % Mf
	f(:,6) = sqrt((1+.203*cosn+.040*cos2n).^2 + (.203*sinn+.040*sin2n).^2);        % Mt

	f(:,7) = 1;                                     % alpha1
	f(:,8) = sqrt((1.+.188*cosn).^2+(.188*sinn).^2);% 2Q1
	f(:,9) = f(:,8);                                % sigma1
	f(:,10) = f(:,8);                               % q1
	f(:,11) = f(:,8);                               % rho1
	f(:,12) = sqrt((1.0+0.189*cosn-0.0058*cos2n).^2 + (0.189*sinn-0.0058*sin2n).^2);% O1
	f(:, 13) = 1;                                   % tau1
	tmp1  = 1.36*cos(p*rad)+.267*cos((p-omega)*rad);% Ray's
	tmp2  = 0.64*sin(p*rad)+.135*sin((p-omega)*rad);
	f(:,14) = sqrt(tmp1.^2 + tmp2.^2);                % M1
	f(:,15) = sqrt((1.+.221*cosn).^2+(.221*sinn).^2);% chi1
	f(:,16) = 1;                                    % pi1
	f(:,17) = 1;                                    % P1
	f(:,18) = 1;                                    % S1
	f(:,19) = sqrt((1.+.1158*cosn-.0029*cos2n).^2 + (.1554*sinn-.0029*sin2n).^2);  % K1
	f(:,20) = 1;                                    % psi1
	f(:,21) = 1;                                    % phi1
	f(:,22) = 1;                                    % theta1
	f(:,23) = sqrt((1.+.169*cosn).^2+(.227*sinn).^2); % J1
	f(:,24) = sqrt((1.0+0.640*cosn+0.134*cos2n).^2 + (0.640*sinn+0.134*sin2n).^2 ); % OO1
	f(:,25) = sqrt((1.-.03731*cosn+.00052*cos2n).^2 + (.03731*sinn-.00052*sin2n).^2);% 2N2
	f(:,26) = f(:,25);                                % mu2
	f(:,27) = f(:,25);                                % N2
	f(:,28) = f(:,25);                                % nu2
	f(:,29) = 1;                                    % M2a
	f(:,30) = f(:,25);                                % M2
	f(:,31) = 1;                                    % M2b
	f(:,32) = 1;                                    % lambda2
	temp1 = 1.-0.25*cos(2*p*rad)-0.11*cos((2*p-omega)*rad)-0.04*cosn;
	temp2 = 0.25*sin(2*p*rad)+0.11*sin((2*p-omega)*rad)+ 0.04*sinn;
	f(:,33) = sqrt(temp1.^2 + temp2.^2);              % L2
	f(:,34) = 1;                                    % T2
	f(:,35) = 1;                                    % S2
	f(:,36) = 1;                                    % R2
	f(:,37) = sqrt((1.+.2852*cosn+.0324*cos2n).^2 + (.3108*sinn+.0324*sin2n).^2);  % K2
	f(:,38) = sqrt((1.+.436*cosn).^2+(.436*sinn).^2); % eta2
	f(:,39) = f(:,30).^2;                            % MNS2
	f(:,40) = f(:,30);                              % 2SM2
	f(:,41) = 1;   % wrong                          % M3
	f(:,42) = f(:,19).*f(:,30);                     % MK3
	f(:,43) = 1;                                    % S3
	f(:,44) = f(:,30).^2;                           % MN4
	f(:,45) = f(:,44);                              % M4
	f(:,46) = f(:,44);                              % MS4
	f(:,47) = f(:,30).*f(:,37);                     % MK4
	f(:,48) = 1;                                    % S4
	f(:,49) = 1;                                    % S5
	f(:,50) = f(:,30).^3;                           % M6
	f(:,51) = 1;                                    % S6
	f(:,52) = 1;                                    % S7
	f(:,53) = 1;                                    % S8

	u = zeros(nT,53);
	u(:, 1) = 0;                                       % Sa
	u(:, 2) = 0;                                       % Ssa
	u(:, 3) = 0;                                       % Mm
	u(:, 4) = 0;                                       % MSf
	u(:, 5) = -23.7*sinn + 2.7*sin2n - 0.4*sin3n;      % Mf
	u(:, 6) = atan(-(.203*sinn+.040*sin2n) ./ (1+.203*cosn+.040*cos2n))/rad;        % Mt
	u(:, 7) = 0;                                       % alpha1
	u(:, 8) = atan(.189*sinn./(1.+.189*cosn))/rad;     % 2Q1
	u(:, 9) = u(:,8);                                  % sigma1
	u(:,10) = u(:,8);                                  % q1
	u(:,11) = u(:,8);                                  % rho1
	u(:,12) = 10.8*sinn - 1.3*sin2n + 0.2*sin3n;       % O1
	u(:,13) = 0;                                       % tau1
	u(:,14) = atan2(tmp2,tmp1)/rad;                    % M1
	u(:,15) = atan(-.221*sinn./(1.+.221*cosn))/rad;    % chi1
	u(:,16) = 0;                                       % pi1
	u(:,17) = 0;                                       % P1
	u(:,18) = 0;                                       % S1
	u(:,19) = atan((-.1554*sinn+.0029*sin2n) ./ (1.+.1158*cosn-.0029*cos2n))/rad;       % K1
	u(:,20) = 0;                                       % psi1
	u(:,21) = 0;                                       % phi1
	u(:,22) = 0;                                       % theta1
	u(:,23) = atan(-.227*sinn./(1.+.169*cosn))/rad;    % J1
	u(:,24) = atan(-(.640*sinn+.134*sin2n) ./ (1.+.640*cosn+.134*cos2n))/rad;         % OO1
	u(:,25) = atan((-.03731*sinn+.00052*sin2n) ./ (1.-.03731*cosn+.00052*cos2n))/rad;     % 2N2
	u(:,26) = u(:,25);                                 % mu2
	u(:,27) = u(:,25);                                 % N2
	u(:,28) = u(:,25);                                 % nu2
	u(:,29) = 0;                                       % M2a
	u(:,30) = u(:,25);                                   % M2
	u(:,31) = 0;                                       % M2b
	u(:,32) = 0;                                       % lambda2
	u(:,33) = atan(-temp2./temp1)/rad ;                % L2
	u(:,34) = 0;                                       % T2
	u(:,35) = 0;                                       % S2
	u(:,36) = 0;                                       % R2
	u(:,37) = atan(-(.3108*sinn+.0324*sin2n) ./ (1.+.2852*cosn+.0324*cos2n))/rad;     % K2
	u(:,38) = atan(-.436*sinn./(1.+.436*cosn))/rad;    % eta2
	u(:,39) = u(:,30)*2;                               % MNS2
	u(:,40) = u(:,30);                                 % 2SM2
	u(:,41) = 1.5d0*u(:,30);                           % M3
	u(:,42) = u(:,30) + u(:,19);                       % MK3
	u(:,43) = 0;                                       % S3
	u(:,44) = u(:,30)*2;                               % MN4
	u(:,45) = u(:,44);                                 % M4
	u(:,46) = u(:,30);                                 % MS4
	u(:,47) = u(:,30)+u(:,37);                         % MK4
	u(:,48) = 0;                                       % S4
	u(:,49) = 0;                                       % S5
	u(:,50) = u(:,30)*3;                               % M6
	u(:,51) = 0;                                       % S6
	u(:,52) = 0;                                       % S7
	u(:,53) = 0;                                       % S8
	% set correspondence between given constituents and supported in OTIS
	ncmx = size(cid0, 1);
	PU = zeros(size(u,1), ncmx);
	PF = zeros(size(u,1), ncmx);
	for (ic = 1:ncmx)
		PU(:,ic) = u(:,index(ic)) * rad;
		PF(:,ic) = f(:,index(ic));
	end
	% take pu,pf for the set of given cid only
	nc = size(cid, 1);
	pu = [];	pf = [];
	for (ic = 1:nc)
		ic0 = 0;
		for (k = 1:ncmx)
			if (cid(ic,:) == cid0(k,:)),	ic0 = k;	break,	end
		end
		if (ic0 > 0)
			pu = [pu, PU(:,ic0)];
			pf = [pf, PF(:,ic0)];
		end
	end

% ----------------------------------------------------------------------------------------
function [h,th_lim,ph_lim] = h_in(cfile,ic)
% reads in elevation for constituent # ic in file cfile
	fid = fopen(cfile,'r','b');
	ll = fread(fid,1,'long');
	nm = fread(fid,3,'long');
	n = nm(1);
	m = nm(2);
	th_lim = fread(fid,2,'float');
	ph_lim = fread(fid,2,'float');
	%nskip = (ic-1)*(nm(1)*nm(2)*8+8) + 8 + 4*nc;
	nskip = (ic-1)*(nm(1)*nm(2)*8+8) + 8 + ll - 28;
	fseek(fid,nskip,'cof');
	htemp = fread(fid,[2*n,m],'float');
	h = htemp(1:2:2*n-1,:) + 1i * htemp(2:2:2*n,:);
	fclose(fid);

% ----------------------------------------------------------------------------------------
function [ll_lims,hz,mz,iob,dt] = grd_in(cfile)
%  reads a grid file in matlab
	fid = fopen(cfile,'r','b');
	fseek(fid,4,'bof');
	n = fread(fid,1,'long');
	m = fread(fid,1,'long');
	lats = fread(fid,2,'float');
	lons = fread(fid,2,'float');
	dt = fread(fid,1,'float');
	if(lons(1) < 0) && (lons(2) < 0 ) && dt>0
		lons = lons + 360;
	end
	ll_lims = [lons; lats];
	nob = fread(fid,1,'long');
	if nob == 0
		fseek(fid,20,'cof');
		iob = [];
	else
		fseek(fid,8,'cof');
		iob = fread(fid,[2,nob],'long');
		fseek(fid,8,'cof');
	end
	hz = fread(fid,[n,m],'*float');
	fseek(fid,8,'cof');
	mz = fread(fid,[n,m],'long');
	fclose(fid);

% ----------------------------------------------------------------------------------------
function [u,v,th_lim,ph_lim] = u_in(cfile,ic, opt)
% reads in transports for constituent # ic in file cfile
% OPT is only used when only one argout is requested. Than it must be == 'u' or 'v'
	fid = fopen(cfile,'r','b');
	ll = fread(fid,1,'long');
	nm = fread(fid,3,'long');
	n = nm(1);
	m = nm(2);
	th_lim = fread(fid,2,'float');
	ph_lim = fread(fid,2,'float');
	%nskip = (ic-1)*(nm(1)*nm(2)*16+8) + 8 + 4*nc;
	nskip = (ic-1)*(nm(1)*nm(2)*16+8) + 8 + ll-28;
	fseek(fid,nskip,'cof');
	htemp = fread(fid,[4*n,m],'float');
	fclose(fid);
	if (nargout == 1)
		if (opt == 'u')
			ur = htemp(1:4:4*n-3,:);		ui = htemp(2:4:4*n-2,:);
			clear htemp;
			u = ur + 1i * ui;
		else
			vr = htemp(3:4:4*n-1,:);		vi = htemp(4:4:4*n,:);
			clear htemp;
			u = vr + 1i * vi;
		end
	else
		ur = htemp(1:4:4*n-3,:);		ui = htemp(2:4:4*n-2,:);
		vr = htemp(3:4:4*n-1,:);		vi = htemp(4:4:4*n,:);
		clear htemp;
		u = ur + 1i * ui;
		clear ur ui; 
		v = vr + 1i * vi;
	end

% ----------------------------------------------------------------------------------------
function [hu,hv] = Huv(hz, opt)
	[n,m] = size(hz);
	indxm = [n,1:n-1];		indym = [m,1:m-1];
	if (nargout == 1)
		if (opt == 'u')
			mu = Muv(hz, opt);
			hu = mu .* (hz + hz(indxm,:)) / 2;
		else
			mv = Muv(hz, opt);
			hu = mv .* (hz + hz(:,indym)) / 2;
		end
	else
		[mu,mv] = Muv(hz);
		hu = mu .* (hz + hz(indxm,:)) / 2;
		hv = mv .* (hz + hz(:,indym)) / 2;
	end

% ----------------------------------------------------------------------------------------
function [mz,mu,mv] = Muv(hz, opt)
%  Given a rectangular bathymetry grid, construct masks for zeta, u and v nodes on a C-grid
	[n,m] = size(hz);
	mz = (hz > 0);
	if (nargout == 1)
		if (opt == 'u')
			return
		else
			indx = [2:n 1];
			mu(indx,:) = mz .* mz(indx,:);
			mz = mu;
			return
		end
	else
		indx = [2:n 1];
		mu(indx,:) = mz .* mz(indx,:);
	end
	if (nargout == 3)
		indy = [2:m 1];
		mv(:,indy) = mz .* mz(:,indy);
	end

% ----------------------------------------------------------------------------------------
function h1 = TMD_OSU_LayoutFcn(h1)

	set(h1, 'Position',[520 405 191 395],...
		'Color',get(0,'factoryUicontrolBackgroundColor'),...
		'MenuBar','none',...
		'Name','Tidal Prediction (OSU)',...
		'NumberTitle','off',...
		'Resize','off',...
		'HandleVisibility','callback',...
		'Tag','figure1');

uicontrol('Parent',h1, 'Position',[90 154 50 14],...
'HorizontalAlignment','left',...
'String','Duration',...
'Style','text');

uicontrol('Parent',h1, 'Position',[90 194 70 14],...
'HorizontalAlignment','left',...
'String','Starting time',...
'Style','text');

uicontrol('Parent',h1, 'Position',[90 303 91 81], 'Style','frame');
uicontrol('Parent',h1, 'Position',[10 83 71 301], 'Style','frame');
uicontrol('Parent',h1, 'Position',[90 213 91 71], 'Style','frame');

uicontrol('Parent',h1, 'Position',[20 359 50 15],...
'Callback',@TMD_uiCB,...
'String','M2',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 1}, ...
'Tag','check_M2');

uicontrol('Parent',h1, 'Position',[20 339 50 15],...
'Callback',@TMD_uiCB,...
'String','S2',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 2}, ...
'Tag','check_S2');

uicontrol('Parent',h1, 'Position',[20 319 50 15],...
'Callback',@TMD_uiCB,...
'String','N2',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 3}, ...
'Tag','check_N2');

uicontrol('Parent',h1, 'Position',[20 299 50 15],...
'Callback',@TMD_uiCB,...
'String','K2',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 4}, ...
'Tag','check_K2');

uicontrol('Parent',h1, 'Position',[20 279 50 15],...
'Callback',@TMD_uiCB,...
'String','K1',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 5}, ...
'Tag','check_K1');

uicontrol('Parent',h1, 'Position',[20 259 50 15],...
'Callback',@TMD_uiCB,...
'String','O1',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 6}, ...
'Tag','check_O1');

uicontrol('Parent',h1, 'Position',[20 239 50 15],...
'Callback',@TMD_uiCB,...
'String','P1',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 7}, ...
'Tag','check_P1');

uicontrol('Parent',h1, 'Position',[20 219 50 15],...
'Callback',@TMD_uiCB,...
'String','Q1',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 8}, ...
'Tag','check_Q1');

uicontrol('Parent',h1, 'Position',[20 199 50 15],...
'Callback',@TMD_uiCB,...
'String','MF',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 9}, ...
'Tag','check_MF');

uicontrol('Parent',h1, 'Position',[20 179 50 15],...
'Callback',@TMD_uiCB,...
'String','MM',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 10}, ...
'Tag','check_MM');

uicontrol('Parent',h1, 'Position',[20 159 50 15],...
'Callback',@TMD_uiCB,...
'String','M4',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 11}, ...
'Tag','check_M4');

uicontrol('Parent',h1, 'Position',[20 139 50 15],...
'Callback',@TMD_uiCB,...
'String','MS4',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 12}, ...
'Tag','check_MS4');

uicontrol('Parent',h1, 'Position',[20 119 50 15],...
'Callback',@TMD_uiCB,...
'String','MN4',...
'Style','checkbox',...
'Tooltip','Left click to select/unselect; right click to plot',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 13}, ...
'Tag','check_MN4');

uicontrol('Parent',h1, 'Position',[20 89 50 15],...
'Callback',@TMD_uiCB,...
'String','All',...
'Style','checkbox',...
'Tooltip','Select all tidal components',...
'Value',1,...
'ButtonDownFcn', {@BDNfcn, 14}, ...
'Tag','check_All');

uicontrol('Parent',h1, 'Position',[100 349 75 15],...
'Callback',@TMD_uiCB,...
'String','Bathymetry',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_plotBat');

uicontrol('Parent',h1, 'Position',[100 329 79 15],...
'Callback',@TMD_uiCB,...
'String','Amplitude',...
'Style','radiobutton',...
'Tag','radio_plotAmp');

uicontrol('Parent',h1, 'Position',[100 310 79 15],...
'Callback',@TMD_uiCB,...
'String','Phase',...
'Style','radiobutton',...
'Tag','radio_plotPhase');

uicontrol('Parent',h1, 'Position',[113 375 41 16],...
'FontSize',10, 'String','Plot',...
'Style','text');

uicontrol('Parent',h1, 'Position',[100 258 40 15],...
'Callback',@TMD_uiCB,...
'String','z',...
'Style','radiobutton',...
'Tooltip','Water elevation',...
'Value',1,...
'Tag','radio_comp_z');

uicontrol('Parent',h1, 'Position',[100 239 30 15],...
'Callback',@TMD_uiCB,...
'String','u',...
'Style','radiobutton',...
'Tooltip','East velocity',...
'Tag','radio_comp_u');

uicontrol('Parent',h1, 'Position',[140 239 30 15],...
'Callback',@TMD_uiCB,...
'String','v',...
'Style','radiobutton',...
'Tooltip','North velocity',...
'Tag','radio_comp_v');

uicontrol('Parent',h1, 'Position',[96 274 80 18],...
'FontSize',10,...
'String','Component',...
'Style','text');

uicontrol('Parent',h1, 'Position',[100 219 30 15],...
'Callback',@TMD_uiCB,...
'String','U',...
'Style','radiobutton',...
'Tooltip','East transport',...
'Tag','radio_comp_ut');

uicontrol('Parent',h1, 'Position',[140 219 30 15],...
'Callback',@TMD_uiCB,...
'String','V',...
'Style','radiobutton',...
'Tooltip','North transport',...
'Tag','radio_comp_vt');

uicontrol('Parent',h1, 'Position',[90 173 71 21],...
'BackgroundColor',[1 1 1],...
'Callback',@TMD_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tooltip','Start date for tide computation',...
'Tag','edit_startDate');

uicontrol('Parent',h1, 'Position',[160 173 23 21],...
'Callback',@TMD_uiCB,...
'String','',...
'Tooltip','Callendar tool',...
'Tag','push_startDate');

uicontrol('Parent',h1, 'Position',[150 154 30 14],...
'HorizontalAlignment','left',...
'String','Step', 'Style','text');

uicontrol('Parent',h1, 'Position',[90 133 31 21],...
'BackgroundColor',[1 1 1],...
'Callback',@TMD_uiCB,...
'String','48',...
'Style','edit',...
'Tooltip','Durantion in hours',...
'Tag','edit_duration');

uicontrol('Parent',h1, 'Position',[150 133 31 21],...
'BackgroundColor',[1 1 1],...
'Callback',@TMD_uiCB,...
'String','15',...
'Style','edit',...
'Tooltip','Step for simulation in minutes',...
'Tag','edit_step');

uicontrol('Parent',h1, 'Position',[90 92 41 23],...
'Callback',@TMD_uiCB,...
'String','Point',...
'Tooltip','Plot a point for a tide station',...
'Tag','push_plotTidePt');

uicontrol('Parent',h1, 'Position',[140 92 41 23],...
'Callback',@TMD_uiCB,...
'String','File',...
'Tooltip','Import a file with tide station locations',...
'Tag','push_loadTidesFile');

uicontrol('Parent',h1, 'Position',[10 43 151 21],...
'BackgroundColor',[1 1 1],...
'Callback',@TMD_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tooltip','Directory with coeficient files',...
'Tag','edit_dataFiles');

uicontrol('Parent',h1, 'Position',[160 43 23 21],...
'Callback',@TMD_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'String','...',...
'Tooltip','Browse for data dir',...
'Tag','push_dataFiles');

uicontrol('Parent',h1, 'Position',[12 64 100 15],...
'HorizontalAlignment','left',...
'String','Data files location',...
'Style','text');

uicontrol('Parent',h1, 'Position',[15 377 62 16],...
'FontSize',10,...
'String','Harmonic',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 11 80 23],...
'Callback',@TMD_uiCB,...
'String','Save Tide',...
'Tooltip','Compute tide at selected location(s) and save them',...
'Tag','push_saveTide');

uicontrol('Parent',h1, 'Position',[103 11 80 23],...
'Callback',@TMD_uiCB,...
'String','Plot Tide',...
'Tooltip','Compute tide at selected location(s) and plot them',...
'Tag','push_plotTide');


function TMD_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
